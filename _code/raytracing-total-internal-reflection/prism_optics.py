"""프리즘 반사/굴절 계산 (다각형 경계 + 전반사).

MATLAB 광선추적기의 MakePrism.m / CalculateRayPath_Prism.m을 이 글에 필요한
만큼(정삼각형 프리즘, 여러 면을 오가며 전반사할 수 있는 경우)만 뽑아 Python으로
재구성했다.

렌즈(sphere_optics.py)와 다른 점은 두 가지다.
1. 경계가 원호가 아니라 선분들의 모음(다각형)이라, 교차점과 법선을
   arccos/arcsin이 아니라 직선 교차와 벡터 연산으로 구한다.
2. 한 번 진입·퇴장으로 끝나지 않고, 더 이상 어떤 경계와도 만나지 않을 때까지
   while 루프를 돈다 — 그 안에서 전반사(TIR)가 몇 번이든 일어날 수 있다.

실행: python prism_optics.py
"""
import numpy as np


def load_nk_file(path):
    """[파장(A), n, k] 3열 텍스트 파일을 읽어 [파장(nm), n, k]로 변환."""
    data = np.loadtxt(path)
    data[:, 0] = data[:, 0] / 10.0  # Angstrom -> nm
    return data


def refractive_index(nk_data, wavelength_nm):
    """목표 파장 주변 2점으로 Cauchy 모델(n = A + B/lambda^2)을 국소 피팅한다."""
    idx = np.searchsorted(nk_data[:, 0], wavelength_nm)
    idx = np.clip(idx, 1, len(nk_data) - 1)
    lam = nk_data[idx - 1 : idx + 1, 0]
    n_vals = nk_data[idx - 1 : idx + 1, 1]
    A = np.column_stack([np.ones(2), 1.0 / lam**2])
    coeff_n = np.linalg.solve(A, n_vals)
    return coeff_n[0] + coeff_n[1] / wavelength_nm**2


def make_prism(ex, ey, cx=0.0, cy=0.0, tx=0.0, ty=0.0, theta=0.0):
    """다각형 꼭짓점 (ex, ey)로 프리즘을 정의한다.

    (cx, cy) 기준으로 theta만큼 회전한 뒤 (tx, ty)만큼 평행이동한다.
    BOUNDARY는 인접한 꼭짓점을 잇는 선분들의 리스트(마지막-첫 꼭짓점도 닫는다).
    """
    ex, ey = np.asarray(ex, dtype=float), np.asarray(ey, dtype=float)
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    XY = rot @ np.vstack([ex - cx, ey - cy]) + np.array([[tx], [ty]])
    verts = XY.T  # (N, 2)

    boundary = [verts[[i, i + 1]] for i in range(len(verts) - 1)]
    boundary.append(verts[[-1, 0]])
    return dict(verts=verts, BOUNDARY=boundary)


def trace_ray_through_prism(ray, prism, n_lens, n_air, max_bounces=20, verbose=False):
    """프리즘 경계와 더 이상 만나지 않을 때까지 광선을 추적한다.

    반환값은 (N, 4) 배열 -- 각 행은 충돌점 [x, y, vx, vy] (그 지점에서 나가는
    방향). 원본 CalculateRayPath_Prism.m의 TRANSIS와 동일한 형태다.
    """
    xs, ys, vx, vy = map(float, ray)
    boundary = prism["BOUNDARY"]
    n_seg = len(boundary)
    air2glass = True
    path = [[xs, ys, vx, vy]]

    for _ in range(max_bounces):
        a_r, b_r = -vy, vx
        c_r = a_r * xs + b_r * ys

        best = None
        for x1, y1, x2, y2 in ((*seg[0], *seg[1]) for seg in boundary):
            a_e, b_e = y2 - y1, x1 - x2
            c_e = a_e * x1 + b_e * y1
            M = np.array([[a_r, b_r], [a_e, b_e]])
            if abs(np.linalg.det(M)) < 1e-8:
                continue
            J = np.linalg.solve(M, [c_r, c_e])
            seg_len_sq = (x1 - x2) ** 2 + (y1 - y2) ** 2
            on_segment = (
                np.sum((J - [x1, y1]) ** 2) <= seg_len_sq
                and np.sum((J - [x2, y2]) ** 2) <= seg_len_sq
            )
            forward = vx * (J[0] - xs) + vy * (J[1] - ys) > 0
            if on_segment and forward:
                dist = np.hypot(xs - J[0], ys - J[1])
                if best is None or dist < best[0]:
                    best = (dist, J, (x1, y1, x2, y2))

        if best is None:
            break
        _, (xc, yc), (x1, y1, x2, y2) = best

        # 면의 법선 -- 선분에 수직인 방향. 광선이 들어오는 쪽을 향하도록 부호를 맞춘다
        # (반사 공식 v'=v-2(v.n)n은 법선 부호에 무관하지만, 부호를 맞춰두면 이후
        # 계산이 더 읽기 쉽다).
        nx, ny = -(y2 - y1), (x2 - x1)
        nnorm = np.hypot(nx, ny)
        nx, ny = nx / nnorm, ny / nnorm
        if vx * nx + vy * ny < 0:
            nx, ny = -nx, -ny

        theta_in = np.arccos(np.clip((vx * nx + vy * ny) / np.hypot(vx, vy), -1, 1))
        n1, n2 = (n_air, n_lens) if air2glass else (n_lens, n_air)
        sin_out = n1 * np.sin(theta_in) / n2

        if abs(sin_out) <= 1:
            theta_out = np.arcsin(sin_out)
            rot_sign = np.sign(vx * ny - vy * nx)
            rot = np.array(
                [
                    [np.cos(theta_in - theta_out), -rot_sign * np.sin(theta_in - theta_out)],
                    [rot_sign * np.sin(theta_in - theta_out), np.cos(theta_in - theta_out)],
                ]
            )
            vout = rot @ np.array([vx, vy])
            air2glass = not air2glass
            if verbose:
                print(f"  face=({x1:.1f},{y1:.1f})-({x2:.1f},{y2:.1f})  "
                      f"theta_in={np.rad2deg(theta_in):.3f}deg  투과  "
                      f"theta_out={np.rad2deg(theta_out):.3f}deg")
        else:
            vout = np.array([vx, vy]) - 2 * (vx * nx + vy * ny) * np.array([nx, ny])
            if verbose:
                print(f"  face=({x1:.1f},{y1:.1f})-({x2:.1f},{y2:.1f})  "
                      f"theta_in={np.rad2deg(theta_in):.3f}deg  전반사(TIR)")

        vout = vout / np.linalg.norm(vout)
        vx, vy = vout
        # 다음 탐색 시작점을 진행방향으로 살짝 밀어둔다 -- 안 그러면 부동소수점
        # 오차 때문에 같은 면을 자기 자신과 다시 충돌한 것으로 잘못 검출한다.
        xs, ys = xc + 0.5 * vx, yc + 0.5 * vy
        path.append([xc, yc, vx, vy])

    return np.array(path)


if __name__ == "__main__":
    nk = load_nk_file("N-BK7.nk")
    n_lens = refractive_index(nk, 750.0)
    print(f"n(N-BK7, 750nm) = {n_lens:.6f}")

    prism = make_prism([0, 300 * np.sqrt(3), 0], [300, 0, -300])
    ray = np.array([-150.0, 100.0, 1.0, 0.0])
    path = trace_ray_through_prism(ray, prism, n_lens, 1.0, verbose=True)
    print(path)

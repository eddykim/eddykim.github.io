"""구면 렌즈 굴절 계산 — 최종(수정) 버전.

MATLAB 광선추적기의 MakeSphericalLens.m / CalculateRayPath_SphericalLens.m /
FindNearleastCollision.m을 이 글에 필요한 만큼만(단일 구면렌즈, 광선 하나가
두 면을 통과하는 경우) 뽑아 Python으로 재구성했다. 원본은 20개 넘는 부품
타입(미러·프리즘·빔스플리터·조리개...)을 한 함수에서 다 처리하는 큰 코드라,
여기서는 그 중 구면렌즈 굴절만 잘라냈다.

실행: python sphere_optics.py
"""
import numpy as np


# ---------------------------------------------------------------------------
# 굴절률 모델 — NK 데이터 파일 로딩 + Cauchy 분산식 국소 피팅
# ---------------------------------------------------------------------------

def load_nk_file(path):
    """[파장(Å), n, k] 3열 텍스트 파일을 읽어 [파장(nm), n, k]로 변환."""
    data = np.loadtxt(path)
    data[:, 0] = data[:, 0] / 10.0  # Å -> nm
    return data


def refractive_index(nk_data, wavelength_nm):
    """목표 파장 주변 2점으로 Cauchy 모델(n = A + B/lambda^2)을 국소 피팅한다.

    전체 파장 대역을 하나의 Cauchy 식으로 맞추면 오차가 크지만, 목표 파장
    근방 2점만 쓰면 국소적으로는 잘 맞는다 (원본 NK_Cauch_Model.m과 동일한 방식).
    """
    idx = np.searchsorted(nk_data[:, 0], wavelength_nm)
    idx = np.clip(idx, 1, len(nk_data) - 1)
    lam = nk_data[idx - 1 : idx + 1, 0]
    n_vals = nk_data[idx - 1 : idx + 1, 1]

    A = np.column_stack([np.ones(2), 1.0 / lam**2])
    coeff_n = np.linalg.solve(A, n_vals)
    n = coeff_n[0] + coeff_n[1] / wavelength_nm**2
    return n


# ---------------------------------------------------------------------------
# 부품 생성
# ---------------------------------------------------------------------------

def make_spherical_lens(R1, R2, t, D, x=0.0, y=0.0, theta=0.0):
    """두 구면(R1, R2)으로 이루어진 렌즈를 정의한다.

    R1: 입사 쪽(양수 = 볼록), R2: 출사 쪽(음수 = 볼록), t: 중심 두께, D: 지름.
    반환하는 dict의 BOUNDARY는 각 면을 20개 점으로 근사한 아크(arc), BD_OFFSET은
    렌즈 가장자리(D/2)에서 구면과 그 접평면 사이의 sag(처짐량) — 부동소수점
    오차로 렌즈 경계 안쪽 점이 "경계 밖"으로 잘못 판정되는 걸 막는 여유값이다.
    """
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    cx1, cy1 = rot @ np.array([R1 - t / 2, 0]) + np.array([x, y])
    cx2, cy2 = rot @ np.array([R2 + t / 2, 0]) + np.array([x, y])

    w1 = np.linspace(np.arcsin(D / (2 * R1)), -np.arcsin(D / (2 * R1)), 21)
    X1, Y1 = -R1 * np.cos(w1) + R1 - t / 2, R1 * np.sin(w1)
    w2 = np.linspace(-np.arcsin(D / (2 * R2)), np.arcsin(D / (2 * R2)), 21)
    X2, Y2 = -R2 * np.cos(w2) + R2 + t / 2, R2 * np.sin(w2)

    boundary1 = (rot @ np.vstack([X1, Y1]) + np.array([[x], [y]])).T
    boundary2 = (rot @ np.vstack([X2, Y2]) + np.array([[x], [y]])).T

    bd_offset = 0.01 * np.array(
        [np.abs(R1) - np.sqrt(R1**2 - (D / 2) ** 2), np.abs(R2) - np.sqrt(R2**2 - (D / 2) ** 2)]
    )

    return dict(
        R1=R1, R2=R2, t=t, D=D, x=x, y=y, theta=theta,
        cx1=cx1, cy1=cy1, cx2=cx2, cy2=cy2,
        BOUNDARY=[boundary1, boundary2], BD_OFFSET=bd_offset,
    )


def make_parallel_beam(x, y, vx, vy, D, num):
    """(x,y)를 중심으로 (vx,vy) 방향 진행벡터에 수직으로 폭 D에 걸쳐 num개
    평행광선을 만든다. 각 광선은 [x, y, vx, vy] 4원소 배열이다."""
    norm = np.hypot(vx, vy)
    unx, uny = -vy / norm, vx / norm
    xs = np.linspace(x - D / 2 * unx, x + D / 2 * unx, num)
    ys = np.linspace(y - D / 2 * uny, y + D / 2 * uny, num)
    return [np.array([xs[i], ys[i], vx, vy]) for i in range(num)]


# ---------------------------------------------------------------------------
# 광선-렌즈 교차 판정
# ---------------------------------------------------------------------------

def find_entry_surface(ray, lens):
    """광선이 렌즈의 두 면 중 어느 쪽에 먼저 닿는지 찾는다.

    광선을 ax+by=c 직선으로, 각 면의 경계를 선분들의 모음으로 보고 교차를
    검사한다 (원본 FindNearleastCollision.m의 단일 부품 버전). 반환값은
    (면 번호(1 또는 2), 교차까지 거리) 또는 못 맞았으면 None.
    """
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y

    best = None
    for surf_idx, boundary in enumerate(lens["BOUNDARY"], start=1):
        for i in range(len(boundary) - 1):
            x1, y1 = boundary[i]
            x2, y2 = boundary[i + 1]
            a_e, b_e = y2 - y1, x1 - x2
            c_e = a_e * x1 + b_e * y1

            M = np.array([[a_r, b_r], [a_e, b_e]])
            if abs(np.linalg.det(M)) < 1e-8:
                continue
            J = np.linalg.solve(M, np.array([c_r, c_e]))

            seg_len_sq = (x1 - x2) ** 2 + (y1 - y2) ** 2
            on_segment = (
                np.sum((J - [x1, y1]) ** 2) <= seg_len_sq
                and np.sum((J - [x2, y2]) ** 2) <= seg_len_sq
            )
            forward = vx * (J[0] - x) + vy * (J[1] - y) > 0
            if not (on_segment and forward):
                continue

            dist = np.hypot(x - J[0], y - J[1])
            if best is None or dist < best[1]:
                best = (surf_idx, dist)
    return best


# ---------------------------------------------------------------------------
# 굴절 계산 (Snell's law) — 최종 버전
# ---------------------------------------------------------------------------

def refract_through_lens(ray, entry_surface, lens, n_lens, n_air):
    """구면렌즈의 두 면을 통과하는 굴절을 계산해 [입사점, 출사점] 2점 경로를
    반환한다. entry_surface가 1이면 R1 쪽에서, 아니면 R2 쪽에서 들어온다."""
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y
    phi_c, phi_s = np.arctan2(b_r, a_r), np.arctan2(a_r, b_r)

    if entry_surface == 1:
        cx1, cy1, R1 = lens["cx1"], lens["cy1"], lens["R1"]
        boundary1, off1 = lens["BOUNDARY"][0], lens["BD_OFFSET"][0]

        w1_c = np.arccos((c_r - a_r * cx1 - b_r * cy1) / (R1 * np.hypot(a_r, b_r))) + phi_c
        w1_s = np.arcsin((c_r - a_r * cx1 - b_r * cy1) / (R1 * np.hypot(a_r, b_r))) - phi_s
        x1_c, y1_c = R1 * np.cos(w1_c) + cx1, R1 * np.sin(w1_c) + cy1
        x1_s, y1_s = R1 * np.cos(w1_s) + cx1, R1 * np.sin(w1_s) + cy1
        w1, x1, y1 = _pick_candidate(boundary1, off1, w1_c, x1_c, y1_c, w1_s, x1_s, y1_s, x, y)

        theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
        vx1, vy1 = np.cos(theta_t1 - (np.pi - w1)), np.sin(theta_t1 - (np.pi - w1))

        cx2, cy2, R2 = lens["cx2"], lens["cy2"], lens["R2"]
        boundary2, off2 = lens["BOUNDARY"][1], lens["BD_OFFSET"][1]
        w2_c = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R2 * np.hypot(vy1, vx1))) + np.arctan2(vx1, -vy1)
        w2_s = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R2 * np.hypot(vy1, vx1))) - np.arctan2(-vy1, vx1)
        x2_c, y2_c = R2 * np.cos(w2_c) + cx2, R2 * np.sin(w2_c) + cy2
        x2_s, y2_s = R2 * np.cos(w2_s) + cx2, R2 * np.sin(w2_s) + cy2
        w2, x2, y2 = _pick_candidate(boundary2, off2, w2_c, x2_c, y2_c, w2_s, x2_s, y2_s, x1, y1)

        phi2 = np.arctan2(vy1, vx1)
        theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)
        vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
    else:
        cx1, cy1, R2 = lens["cx2"], lens["cy2"], lens["R2"]
        boundary1, off1 = lens["BOUNDARY"][1], lens["BD_OFFSET"][1]

        w1_c = np.arccos((c_r - a_r * cx1 - b_r * cy1) / (R2 * np.hypot(a_r, b_r))) + phi_c
        w1_s = np.arcsin((c_r - a_r * cx1 - b_r * cy1) / (R2 * np.hypot(a_r, b_r))) - phi_s
        x1_c, y1_c = R2 * np.cos(w1_c) + cx1, R2 * np.sin(w1_c) + cy1
        x1_s, y1_s = R2 * np.cos(w1_s) + cx1, R2 * np.sin(w1_s) + cy1
        w1, x1, y1 = _pick_candidate(boundary1, off1, w1_c, x1_c, y1_c, w1_s, x1_s, y1_s, x, y)

        theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
        vx1, vy1 = np.cos(-theta_t1 - (np.pi - w1)), np.sin(-theta_t1 - (np.pi - w1))

        cx2, cy2, R1 = lens["cx1"], lens["cy1"], lens["R1"]
        boundary2, off2 = lens["BOUNDARY"][0], lens["BD_OFFSET"][0]
        w2_c = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R1 * np.hypot(vy1, vx1))) + np.arctan2(vx1, -vy1)
        w2_s = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R1 * np.hypot(vy1, vx1))) - np.arctan2(-vy1, vx1)
        x2_c, y2_c = R1 * np.cos(w2_c) + cx2, R1 * np.sin(w2_c) + cy2
        x2_s, y2_s = R1 * np.cos(w2_s) + cx2, R1 * np.sin(w2_s) + cy2
        w2, x2, y2 = _pick_candidate(boundary2, off2, w2_c, x2_c, y2_c, w2_s, x2_s, y2_s, x1, y1)

        phi2 = np.arctan2(vy1, vx1)
        theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)
        vx2, vy2 = np.cos(theta_t2 + w2), np.sin(theta_t2 + w2)

    return np.array([[x1, y1, vx1, vy1], [x2, y2, vx2, vy2]])


def _pick_candidate(boundary, offset, wc, xc, yc, ws, xs, ys, ref_x, ref_y):
    """arccos/arcsin 두 후보해 중 렌즈 물리 경계 안에 있는 쪽을 고른다.

    두 후보 다 경계 안이면(보통의 경우) 기준점에 더 가까운 쪽을 쓴다.
    """
    def inside(px, py):
        return (
            px <= np.max(boundary[:, 0]) + offset and px >= np.min(boundary[:, 0]) - offset
            and py <= np.max(boundary[:, 1]) + offset and py >= np.min(boundary[:, 1]) - offset
        )

    in_c, in_s = inside(xc, yc), inside(xs, ys)
    if in_c and in_s:
        if np.hypot(xc - ref_x, yc - ref_y) < np.hypot(xs - ref_x, ys - ref_y):
            return wc, xc, yc
        return ws, xs, ys
    elif in_c:
        return wc, xc, yc
    else:
        return ws, xs, ys


# ---------------------------------------------------------------------------
# 데모: 평행광이 양볼록 렌즈를 통과하는 전체 경로
# ---------------------------------------------------------------------------

def trace_parallel_beam_through_lens(lens, source_rays, wavelength_nm, n_air, nk_data, roi):
    """각 광선에 대해 [입사 전, 렌즈 통과, ROI 경계까지] 전체 경로(N,4)를 계산."""
    n_lens = refractive_index(nk_data, wavelength_nm)
    paths = []
    for ray in source_rays:
        entry = find_entry_surface(ray, lens)
        if entry is None:
            paths.append(np.array([ray]))
            continue
        surf_idx, _ = entry
        inside = refract_through_lens(ray, surf_idx, lens, n_lens, n_air)
        end = _extend_to_roi(inside[-1], roi)
        paths.append(np.vstack([ray, inside, end]))
    return paths


def _extend_to_roi(ray, roi):
    """마지막 광선 위치·방향에서 ROI(x_min,x_max,y_min,y_max) 경계까지 직선으로 연장."""
    x, y, vx, vy = ray
    x_min, x_max, y_min, y_max = roi
    candidates = []
    if vx > 0:
        candidates.append(((x_max - x) / vx, x_max, y + vy * (x_max - x) / vx))
    elif vx < 0:
        candidates.append(((x_min - x) / vx, x_min, y + vy * (x_min - x) / vx))
    if vy > 0:
        candidates.append(((y_max - y) / vy, x + vx * (y_max - y) / vy, y_max))
    elif vy < 0:
        candidates.append(((y_min - y) / vy, x + vx * (y_min - y) / vy, y_min))
    t, xe, ye = min((c for c in candidates if c[0] > 0), key=lambda c: c[0])
    return np.array([xe, ye, 0.0, 0.0])


if __name__ == "__main__":
    nk = load_nk_file("N-BK7.nk")
    lens = make_spherical_lens(R1=1000, R2=-1000, t=100, D=500)
    source = make_parallel_beam(x=1000, y=0, vx=-1, vy=0, D=200, num=11)
    paths = trace_parallel_beam_through_lens(lens, source, 750, 1.0002778, nk, roi=(-1500, 1500, -500, 500))
    for i, p in enumerate(paths):
        print(f"ray {i}: {p.shape[0]} points, end={p[-1][:2]}")

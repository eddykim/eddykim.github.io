"""구면 렌즈 굴절 + 출사면 전반사(TIR) 처리 -- 최종(수정) 버전.

1편(`raytracing-spherical-lens-refraction`)의 sphere_optics.py를 이 글에
필요한 만큼(단일 렌즈, 출사면에서 임계각을 넘을 수 있는 경우) 다시 가져와,
prism_optics.py와 같은 패턴(법선 벡터 기준 반사, v'=v-2(v.n)n)으로 출사면 TIR
처리를 추가했다. 실제 저장소(code_python_raytrace_250921/calculate.py)의
CalculateRayPath_SphericalLens에 적용한 수정과 같은 내용이다.

실행: python lens_tir.py
"""
import numpy as np


def load_nk_file(path):
    data = np.loadtxt(path)
    data[:, 0] = data[:, 0] / 10.0  # Angstrom -> nm
    return data


def refractive_index(nk_data, wavelength_nm):
    idx = np.searchsorted(nk_data[:, 0], wavelength_nm)
    idx = np.clip(idx, 1, len(nk_data) - 1)
    lam = nk_data[idx - 1 : idx + 1, 0]
    n_vals = nk_data[idx - 1 : idx + 1, 1]
    A = np.column_stack([np.ones(2), 1.0 / lam**2])
    coeff_n = np.linalg.solve(A, n_vals)
    return coeff_n[0] + coeff_n[1] / wavelength_nm**2


def make_spherical_lens(R1, R2, t, D, x=0.0, y=0.0, theta=0.0):
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
    norm = np.hypot(vx, vy)
    unx, uny = -vy / norm, vx / norm
    xs = np.linspace(x - D / 2 * unx, x + D / 2 * unx, num)
    ys = np.linspace(y - D / 2 * uny, y + D / 2 * uny, num)
    return [np.array([xs[i], ys[i], vx, vy]) for i in range(num)]


def _pick_candidate(boundary, offset, wc, xc, yc, ws, xs, ys, ref_x, ref_y):
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


def find_entry_surface(ray, lens):
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


def refract_through_lens(ray, entry_surface, lens, n_lens, n_air):
    """구면렌즈의 두 면을 통과하는 굴절을 계산한다.

    출사면(두 번째 면)에서 임계각을 넘으면 -- Snell's law의 arcsin이 정의역을
    벗어나면 -- 굴절 대신 그 면의 법선을 기준으로 반사시킨다(전반사, TIR).
    반사 공식은 prism_optics.py와 완전히 같다: v' = v - 2(v.n)n.
    """
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
        sin_t2 = n_lens * np.sin(w2 - phi2) / n_air
        if abs(sin_t2) <= 1:
            theta_t2 = np.arcsin(sin_t2)
            vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
        else:
            n2x, n2y = np.cos(w2), np.sin(w2)
            dot2 = vx1 * n2x + vy1 * n2y
            vx2, vy2 = vx1 - 2 * dot2 * n2x, vy1 - 2 * dot2 * n2y
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
        sin_t2 = n_lens * np.sin(w2 - phi2) / n_air
        if abs(sin_t2) <= 1:
            theta_t2 = np.arcsin(sin_t2)
            vx2, vy2 = np.cos(theta_t2 + w2), np.sin(theta_t2 + w2)
        else:
            n2x, n2y = np.cos(w2), np.sin(w2)
            dot2 = vx1 * n2x + vy1 * n2y
            vx2, vy2 = vx1 - 2 * dot2 * n2x, vy1 - 2 * dot2 * n2y

    return np.array([[x1, y1, vx1, vy1], [x2, y2, vx2, vy2]])


def exit_surface_is_tir(path_row0, entry_surface, lens, n_lens, n_air):
    """refract_through_lens가 반환한 첫 행([x1,y1,vx1,vy1], 출사면으로 향하는
    내부 광선)만 가지고, 그 출사면에서 전반사가 일어났는지를 다시 판정한다.

    vx2 부호로는 판별할 수 없다 -- 출사면 법선이 충분히 기울어 있으면 반사
    후에도 진행방향의 x성분이 양수로 남을 수 있다(이 글의 평철렌즈가 그런
    경우다). 그래서 실제로 arcsin 정의역을 다시 검사한다. entry_surface는
    find_entry_surface가 돌려준 값을 그대로 넘긴다 -- 어느 면으로 들어왔는지에
    따라 출사면(R1/R2 중 나머지 하나)이 달라진다.
    """
    x1, y1, vx1, vy1 = path_row0
    if entry_surface == 1:
        cx2, cy2, R2 = lens["cx2"], lens["cy2"], lens["R2"]
        boundary2, off2 = lens["BOUNDARY"][1], lens["BD_OFFSET"][1]
    else:
        cx2, cy2, R2 = lens["cx1"], lens["cy1"], lens["R1"]
        boundary2, off2 = lens["BOUNDARY"][0], lens["BD_OFFSET"][0]
    w2_c = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R2 * np.hypot(vy1, vx1))) + np.arctan2(vx1, -vy1)
    w2_s = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R2 * np.hypot(vy1, vx1))) - np.arctan2(-vy1, vx1)
    x2_c, y2_c = R2 * np.cos(w2_c) + cx2, R2 * np.sin(w2_c) + cy2
    x2_s, y2_s = R2 * np.cos(w2_s) + cx2, R2 * np.sin(w2_s) + cy2
    w2, x2, y2 = _pick_candidate(boundary2, off2, w2_c, x2_c, y2_c, w2_s, x2_s, y2_s, x1, y1)
    phi2 = np.arctan2(vy1, vx1)
    sin_t2 = n_lens * np.sin(w2 - phi2) / n_air
    return abs(sin_t2) > 1


if __name__ == "__main__":
    nk = load_nk_file("N-BK7.nk")
    n_lens = refractive_index(nk, 750.0)
    lens = make_spherical_lens(R1=1e6, R2=-100.0, t=20.0, D=180.0)
    ray = np.array([-1000.0, -88.2, 1.0, 0.0])
    out = refract_through_lens(ray, 1, lens, n_lens, 1.0)
    print(f"n(N-BK7,750nm)={n_lens:.6f}")
    print(out)

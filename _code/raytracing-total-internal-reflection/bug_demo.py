"""출사면 TIR 미처리 재현과 수정.

1편 마지막에 남긴 숙제: refract_through_lens/CalculateRayPath_SphericalLens는
출사면에서 arcsin 정의역을 확인하지 않는다. 여기서는 그걸 실제로 임계각을
넘는 렌즈(평철렌즈, R1=1e6(평면 근사), R2=-50, D=78, y0=36)로 재현하고,
수정 전/후를 나란히 보여준다.

이 함수는 광선이 실제로 어느 면으로 먼저 들어오는지를 스스로 판단하지
않는다 -- 호출하는 쪽이 entry_surface를 지정해야 한다. R1/R2는 렌즈를
정의할 때의 매개변수 순서일 뿐 광원에서 봤을 때 어느 면이 먼저인지와는
무관하므로(1편 참고), 항상 find_entry_surface로 먼저 판정한 값을 넘겨야
한다. 여기서도 그렇게 한다.

실행: python bug_demo.py
"""
import warnings

import numpy as np

from lens_tir import find_entry_surface, load_nk_file, make_spherical_lens, refractive_index
from lens_tir import refract_through_lens as refract_through_lens_fixed


def refract_through_lens_original(ray, entry_surface, lens, n_lens, n_air):
    """1편에서 쓰던, 출사면 TIR을 처리하지 않는 원래 버전 (entry_surface=1만)."""
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y
    phi_c, phi_s = np.arctan2(b_r, a_r), np.arctan2(a_r, b_r)

    cx1, cy1, R1 = lens["cx1"], lens["cy1"], lens["R1"]
    boundary1, off1 = lens["BOUNDARY"][0], lens["BD_OFFSET"][0]
    w1_c = np.arccos((c_r - a_r * cx1 - b_r * cy1) / (R1 * np.hypot(a_r, b_r))) + phi_c
    w1_s = np.arcsin((c_r - a_r * cx1 - b_r * cy1) / (R1 * np.hypot(a_r, b_r))) - phi_s
    x1_c, y1_c = R1 * np.cos(w1_c) + cx1, R1 * np.sin(w1_c) + cy1
    x1_s, y1_s = R1 * np.cos(w1_s) + cx1, R1 * np.sin(w1_s) + cy1

    def inside(px, py, boundary, offset):
        return (
            px <= np.max(boundary[:, 0]) + offset and px >= np.min(boundary[:, 0]) - offset
            and py <= np.max(boundary[:, 1]) + offset and py >= np.min(boundary[:, 1]) - offset
        )

    in_c, in_s = inside(x1_c, y1_c, boundary1, off1), inside(x1_s, y1_s, boundary1, off1)
    if in_c and in_s:
        w1, x1, y1 = (
            (w1_c, x1_c, y1_c)
            if np.hypot(x1_c - x, y1_c - y) < np.hypot(x1_s - x, y1_s - y)
            else (w1_s, x1_s, y1_s)
        )
    else:
        w1, x1, y1 = (w1_c, x1_c, y1_c) if in_c else (w1_s, x1_s, y1_s)

    theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
    vx1, vy1 = np.cos(theta_t1 - (np.pi - w1)), np.sin(theta_t1 - (np.pi - w1))

    cx2, cy2, R2 = lens["cx2"], lens["cy2"], lens["R2"]
    boundary2, off2 = lens["BOUNDARY"][1], lens["BD_OFFSET"][1]
    w2_c = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R2 * np.hypot(vy1, vx1))) + np.arctan2(vx1, -vy1)
    w2_s = np.arccos((-vy1 * x1 + vx1 * y1 + vy1 * cx2 - vx1 * cy2) / (R2 * np.hypot(vy1, vx1))) - np.arctan2(-vy1, vx1)
    x2_c, y2_c = R2 * np.cos(w2_c) + cx2, R2 * np.sin(w2_c) + cy2
    x2_s, y2_s = R2 * np.cos(w2_s) + cx2, R2 * np.sin(w2_s) + cy2
    in_c2, in_s2 = inside(x2_c, y2_c, boundary2, off2), inside(x2_s, y2_s, boundary2, off2)
    if in_c2 and in_s2:
        w2, x2, y2 = (
            (w2_c, x2_c, y2_c)
            if np.hypot(x2_c - x1, y2_c - y1) < np.hypot(x2_s - x1, y2_s - y1)
            else (w2_s, x2_s, y2_s)
        )
    else:
        w2, x2, y2 = (w2_c, x2_c, y2_c) if in_c2 else (w2_s, x2_s, y2_s)

    phi2 = np.arctan2(vy1, vx1)
    theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)  # <- 여기서 도메인 체크가 없다
    vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
    return np.array([[x1, y1, vx1, vy1], [x2, y2, vx2, vy2]])


nk = load_nk_file("N-BK7.nk")
n_lens = refractive_index(nk, 750.0)
print(f"n(N-BK7, 750nm) = {n_lens:.6f}")

lens = make_spherical_lens(R1=1e6, R2=-50.0, t=20.0, D=78.0)
ray = np.array([-1000.0, 36.0, 1.0, 0.0])
entry = find_entry_surface(ray, lens)
print(f"find_entry_surface: {entry[0]}번 면")

with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    out_before = refract_through_lens_original(ray, entry[0], lens, n_lens, 1.0)
print("\n수정 전:")
print(out_before)
if w:
    print(f"-> RuntimeWarning: {w[0].message}")

out_after = refract_through_lens_fixed(ray, entry[0], lens, n_lens, 1.0)
print("\n수정 후:")
print(out_after)

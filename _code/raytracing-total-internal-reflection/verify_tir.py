"""2편에서 인용하는 세 가지 검증 수치를 실제로 계산한다.

1. 프리즘 임계 입사각(theta1) 이분법 스윕 vs 이론값
2. 렌즈 출사면 TIR 반사의 반사법칙(입사각=반사각) 확인
3. 렌즈 출사면 TIR 발생 경계(y0) 이분법 스윕과 그 지점의 내부입사각 vs 임계각

실행: python verify_tir.py
"""
import numpy as np

from lens_tir import (
    _pick_candidate,
    exit_surface_is_tir,
    find_entry_surface,
    load_nk_file as load_nk_lens,
    make_spherical_lens,
    refract_through_lens,
    refractive_index as refractive_index_lens,
)
from prism_optics import (
    load_nk_file as load_nk_prism,
    make_prism,
    refractive_index as refractive_index_prism,
    trace_ray_through_prism,
)

WAVELENGTH_NM = 750.0
TOL = 1e-9  # 이분법 종료 조건 (동일 단위: 각도는 deg, 위치는 mm)


def bisect(is_true_at, lo, hi, tol=TOL):
    """is_true_at(lo)=True, is_true_at(hi)=False라고 가정하고 경계를 이분법으로 좁힌다."""
    assert is_true_at(lo) and not is_true_at(hi)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if is_true_at(mid):
            lo = mid
        else:
            hi = mid
    return lo


# ============================================================
# 1. 프리즘 임계 입사각 -- 이분법 스윕 vs 이론값
# ============================================================
nk_prism = load_nk_prism("N-BK7.nk")
n_glass = refractive_index_prism(nk_prism, WAVELENGTH_NM)
prism = make_prism([0, 300 * np.sqrt(3), 0], [300, 0, -300])


def fire(theta1_deg, y_hit=100.0, x0=-150.0):
    theta1 = np.deg2rad(theta1_deg)
    vx, vy = np.cos(theta1), np.sin(theta1)
    y0 = y_hit - (0 - x0) * np.tan(theta1)
    ray = np.array([x0, y0, vx, vy])
    return trace_ray_through_prism(ray, prism, n_glass, 1.0)


def prism_is_tir(theta1_deg):
    return len(fire(theta1_deg)) >= 4  # TIR 1회 이상이면 경로 점이 4개 이상


theta1_boundary = bisect(prism_is_tir, 0.0, 60.0)

theta_c = np.arcsin(1.0 / n_glass)
A_apex = np.deg2rad(60.0)
theta1_theory = np.rad2deg(np.arcsin(n_glass * np.sin(A_apex - theta_c)))

print("1. 프리즘 임계 입사각")
print(f"   이분법 스윕: {theta1_boundary:.6f}deg")
print(f"   이론값:      {theta1_theory:.6f}deg")
print(f"   오차:        {abs(theta1_boundary - theta1_theory):.1e}deg")

# ============================================================
# 2. 렌즈 출사면 TIR 반사의 반사법칙 확인
# ============================================================
nk_lens = load_nk_lens("N-BK7.nk")
n_lens = refractive_index_lens(nk_lens, WAVELENGTH_NM)
lens = make_spherical_lens(R1=1e6, R2=-50.0, t=20.0, D=78.0)


def exit_normal(path_row0, entry_surface):
    """refract_through_lens가 반환한 첫 행(x1,y1,vx1,vy1)으로부터 출사면
    법선(cos w2, sin w2)을 다시 구한다 -- exit_surface_is_tir 내부와 동일한 절차."""
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
    w2, _, _ = _pick_candidate(boundary2, off2, w2_c, x2_c, y2_c, w2_s, x2_s, y2_s, x1, y1)
    return np.array([np.cos(w2), np.sin(w2)])


ray = np.array([-1000.0, 36.0, 1.0, 0.0])
entry = find_entry_surface(ray, lens)
out = refract_through_lens(ray, entry[0], lens, n_lens, 1.0)
n_hat = exit_normal(out[0], entry[0])
vin, vout = out[0, 2:4], out[1, 2:4]
ang_in = np.rad2deg(np.arccos(abs(np.dot(vin, n_hat))))
ang_out = np.rad2deg(np.arccos(abs(np.dot(vout, n_hat))))

print("\n2. 렌즈 출사면 TIR 반사법칙 (y0=36mm)")
print(f"   입사각: {ang_in:.8f}deg")
print(f"   반사각: {ang_out:.8f}deg")
print(f"   차이:   {ang_in - ang_out:.1e}deg")

# ============================================================
# 3. 렌즈 출사면 TIR 발생 경계(y0) -- 이분법 스윕 vs 임계각
# ============================================================


def lens_tir_at(y0):
    ray = np.array([-1000.0, y0, 1.0, 0.0])
    entry = find_entry_surface(ray, lens)
    out = refract_through_lens(ray, entry[0], lens, n_lens, 1.0)
    return exit_surface_is_tir(out[0], entry[0], lens, n_lens, 1.0)


y0_boundary = bisect(lambda y0: not lens_tir_at(y0), 0.0, 39.0)
# not lens_tir_at(y0): y0=0(투과)에서 True, y0=39(TIR)에서 False가 되도록
# lens_tir_at 자체를 뒤집어 bisect에 넣었으므로, 경계는 "TIR이 막 시작되는 y0"다.

ray = np.array([-1000.0, y0_boundary, 1.0, 0.0])
entry = find_entry_surface(ray, lens)
out = refract_through_lens(ray, entry[0], lens, n_lens, 1.0)
n_hat = exit_normal(out[0], entry[0])
vin = out[0, 2:4]
internal_angle = np.rad2deg(np.arccos(np.clip(abs(np.dot(vin, n_hat)), -1, 1)))
theta_c_lens = np.rad2deg(np.arcsin(1.0 / n_lens))

print("\n3. 렌즈 출사면 TIR 발생 경계 (y0 스윕)")
print(f"   경계 y0:       {y0_boundary:.6f}mm")
print(f"   그 지점 내부입사각: {internal_angle:.6f}deg")
print(f"   이론 임계각:       {theta_c_lens:.6f}deg")
print(f"   오차:             {abs(internal_angle - theta_c_lens):.1e}deg")

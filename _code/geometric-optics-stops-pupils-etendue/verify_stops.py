"""stops.py 를 독립된 경로들과 교차검증한다.

실행: python verify_stops.py

  A. Hecht 예제 5.6             ↔  책의 답 (구멍의 상 −40 cm·지름 200 mm, 구멍이 조리개, 상점 20 cm), 조리개가 바뀌는 거리
  B. 시험계의 입사동·출사동       ↔  Optiland paraxial (EPL·EPD·XPL·XPD), opticore entrance_pupil
  C. 라그랑주 불변량             ↔  면마다 보존, 물체 높이 × 물체 쪽 개구 = 상 높이 × 상 쪽 개구
  D. cos⁴ 법칙                  ↔  동공 원판 위 조도 적분 (작은 동공), Born & Wolf 4.8.3 식 (28)
  E. 비네팅 비율 (격자)           ↔  두 원이 겹친 넓이의 닫힌 식 (한 구경만 자를 때)
  F. 텔레센트릭                  ↔  조리개를 뒤 초점면에 두면 크기가 초점 어긋남에 무관, 일반 배치는 y·L/(L−δ)
  G. 밝기 정리                   ↔  돋보기 초점의 조도 = 태양 복사휘도 × π sin²θ′ ≤ 태양 표면의 방출도
"""
import json
import os
import subprocess

import numpy as np

from generate_figures import N_GLASS, Y_OBJ, Z_OBJ, image_plane, measured_size, test_system
from stops import (Element, Pupil, find_aperture_stop, image_of_aperture, lagrange_invariant,
                   marginal_and_chief, pupil_irradiance, trace, vignetting_fraction)

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "repo_opticore")


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<60s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def verify_hecht():
    print("A. Hecht 예제 5.6")
    els = [Element(0.0, "thin", (100.0,), 70.0), Element(80.0, "stop", (), 20.0)]
    stop, table = find_aperture_stop(els, -200.0)
    ep = table[1][1]
    ok = check("구멍이 조리개다", stop, 1, 0.5)
    ok &= check("구멍의 상 위치: 렌즈 오른쪽 400 mm (si = −40 cm)", ep.z, 400.0, 1e-6)
    ok &= check("구멍의 상 지름 200 mm", 2 * ep.radius, 200.0, 1e-6)
    pts = trace(els, -200.0, 0.0, 0.1)
    z_img = pts[-1, 0] - pts[-1, 1] / pts[-1, 2]
    ok &= check("상점 P 는 렌즈 뒤 200 mm", z_img, 200.0, 1e-6)
    print(f"      물점에서 본 각: 렌즈 테두리 {np.degrees(table[0][2]):.2f}°, 입사동 {np.degrees(table[1][2]):.2f}°")
    # 물점을 멀리 옮기면 렌즈 테두리가 조리개가 된다: 70/s = 100/(400+s) 에서 s = 2800/3 mm
    s_switch = 2800.0 / 3.0
    ok &= check("물점 900 mm: 구멍이 조리개", find_aperture_stop(els, -900.0)[0], 1, 0.5)
    ok &= check("물점 1000 mm: 렌즈가 조리개", find_aperture_stop(els, -1000.0)[0], 0, 0.5)
    print(f"      조리개가 바뀌는 물체 거리 {s_switch:.1f} mm")
    return ok


OPTILAND = r"""
import json, sys, numpy as np
from optiland import optic
from optiland.materials import IdealMaterial
n = float(sys.argv[1]); t_img = float(sys.argv[2])
s = lambda v: float(np.asarray(v).reshape(-1)[0])
glass, air = IdealMaterial(n=n), IdealMaterial(n=1.0)
L = optic.Optic()
L.surfaces.add(index=0, radius=float("inf"), thickness=150.0, material=air)
L.surfaces.add(index=1, radius=60.0, thickness=8.0, material=glass)
L.surfaces.add(index=2, radius=-60.0, thickness=25.0, material=air)
L.surfaces.add(index=3, radius=float("inf"), thickness=22.0, material=air, is_stop=True)
L.surfaces.add(index=4, radius=80.0, thickness=6.0, material=glass)
L.surfaces.add(index=5, radius=-80.0, thickness=t_img, material=air)
L.surfaces.add(index=6)
L.set_aperture(aperture_type="float_by_stop_size", value=10.0)
L.fields.set_type("object_height"); L.fields.add(y=0.0); L.fields.add(y=10.0)
L.wavelengths.add(value=0.5876, is_primary=True)
p = L.paraxial
print(json.dumps(dict(EPL=s(p.EPL()), EPD=s(p.EPD()), XPL=s(p.XPL()), XPD=s(p.XPD()))))
"""

OPTICORE = r"""
import json, sys
from opticore.geometry import Plane, Sphere
from opticore.interactions import Refractive
from opticore.materials import ConstantMaterial
from opticore.surface import Surface
from opticore.system import OpticalSystem
from opticore.pupil import entrance_pupil
n = float(sys.argv[1])
g, a = ConstantMaterial(n), ConstantMaterial(1.0)
s = OpticalSystem(ambient=a)
s.add(Surface(Sphere(60.0), g, Refractive(), position=(0, 0, 0.0)))
s.add(Surface(Sphere(-60.0), a, Refractive(), position=(0, 0, 8.0)))
s.add(Surface(Plane(), a, Refractive(), position=(0, 0, 33.0), stop_radius=5.0))
s.add(Surface(Sphere(80.0), g, Refractive(), position=(0, 0, 55.0)))
s.add(Surface(Sphere(-80.0), a, Refractive(), position=(0, 0, 61.0)))
p = entrance_pupil(s, 0.5876)
print(json.dumps(dict(z=float(p.z), r=float(p.radius_mm))))
"""


def run(py, code, *args):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    out = subprocess.run([py, "-c", code, *map(str, args)], capture_output=True, text=True,
                         env=env, check=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


def verify_pupils():
    print("B. 시험계의 입사동·출사동")
    els = test_system()
    z_img = image_plane(els)
    ep = image_of_aperture(els, 2, "object")
    xp = image_of_aperture(els, 2, "image")
    print(f"      입사동 z {ep.z:.4f} 반지름 {ep.radius:.4f} | 출사동 z {xp.z:.4f} 반지름 {xp.radius:.4f} | 상면 z {z_img:.4f}")
    ok = True
    py = os.path.join(REPO, "reference_libs", ".venv-oracle", "bin", "python")
    if os.path.exists(py):
        o = run(py, OPTILAND, N_GLASS, z_img - 61.0)
        ok &= check("입사동 위치 (첫 면 기준) ↔ Optiland EPL", ep.z, o["EPL"], 1e-6)
        ok &= check("입사동 지름 ↔ Optiland EPD", 2 * ep.radius, o["EPD"], 1e-6)
        ok &= check("출사동 위치 (상면 기준) ↔ Optiland XPL", xp.z - z_img, o["XPL"], 1e-6)
        ok &= check("출사동 지름 ↔ Optiland XPD", 2 * xp.radius, o["XPD"], 1e-6)
    else:
        print(f"  [SKIP] Optiland 오라클 환경이 없다 ({py})")
    py = os.path.join(REPO, "opticore", ".venv", "bin", "python")
    if os.path.exists(py):
        o = run(py, OPTICORE, N_GLASS)
        # opticore 는 광축 근처 실광선으로 근축값을 잰다 (탐침 기울기 1e-4)
        ok &= check("입사동 위치 ↔ opticore entrance_pupil", ep.z, o["z"], 1e-4)
        ok &= check("입사동 반지름 ↔ opticore entrance_pupil", ep.radius, o["r"], 1e-4)
    else:
        print(f"  [SKIP] opticore 환경이 없다 ({py})")
    return ok


def verify_lagrange():
    print("C. 라그랑주 불변량")
    els = test_system()
    z_img = image_plane(els)
    _, ep, marg, chief = marginal_and_chief(els, Z_OBJ, Y_OBJ, z_end=z_img)
    H = lagrange_invariant(marg, chief)
    ok = check("면마다 H 가 같다", H, H[0], 1e-12)
    ok &= check("물체 높이 × 물체 쪽 개구 = 상 높이 × 상 쪽 개구 (크기)",
                abs(Y_OBJ * marg[0, 2]), abs(chief[-1, 1] * marg[-1, 2]), 1e-12)
    print(f"      H = {H[0]:.5f} mm·rad (물체 높이 {Y_OBJ} mm × 개구 {marg[0, 2]:.5f})")
    return ok


def verify_cos4():
    print("D. cos⁴ 법칙")
    phis = np.deg2rad([10, 20, 30, 40])
    small = [pupil_irradiance(0.01, 1.0, np.tan(p)) / pupil_irradiance(0.01, 1.0, 0.0) for p in phis]
    ok = check("작은 동공 (a/L = 0.01) ↔ cos⁴ φ", small, np.cos(phis) ** 4, 2e-4)
    big = [pupil_irradiance(0.5, 1.0, np.tan(p)) / pupil_irradiance(0.5, 1.0, 0.0) for p in phis]
    print("      큰 동공 (a/L = 0.5) / cos⁴: " + ", ".join(f"{b / c:.3f}" for b, c in zip(big, np.cos(phis) ** 4)))
    # 축상 조도의 절대값: B&W 식 (24) E = π B sin²θ (B = 1)
    a = 0.3
    ok &= check("축상 조도 ↔ π sin²θ (B&W 4.8.3 식 24)", pupil_irradiance(a, 1.0, 0.0),
                np.pi * a**2 / (1 + a**2), 1e-6)
    return ok


def circle_overlap(r1, r2, d):
    if d >= r1 + r2:
        return 0.0
    if d <= abs(r1 - r2):
        return np.pi * min(r1, r2) ** 2
    a1 = r1**2 * np.arccos((d**2 + r1**2 - r2**2) / (2 * d * r1))
    a2 = r2**2 * np.arccos((d**2 + r2**2 - r1**2) / (2 * d * r2))
    k = 0.5 * np.sqrt((-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2))
    return a1 + a2 - k


def verify_vignetting():
    print("E. 비네팅")
    els = test_system()
    ep = image_of_aperture(els, 2, "object")
    y = 35.0
    # 첫 면의 구경을 물점에서 입사동 평면으로 투영한 원 (이 높이에서 이것만 자른다)
    s = (ep.z - Z_OBJ) / (0.0 - Z_OBJ)
    cen, rad = y * (1 - s), els[0].radius * s
    exact = circle_overlap(ep.radius, rad, abs(cen)) / (np.pi * ep.radius**2)
    grid = vignetting_fraction(els, Z_OBJ, y, ep, n_grid=1201)
    ok = check(f"물체 높이 {y:.0f} mm: 격자 ↔ 두 원 겹침 닫힌 식", grid, exact, 2e-3)
    print(f"      통과 비율 {exact:.4f}")
    return ok


def verify_telecentric():
    print("F. 텔레센트릭")
    f = 50.0
    d = np.linspace(-5, 5, 11)
    ref = measured_size(0.0, 0.0)
    tele = [measured_size(f * (1 - 1e-12), x) / ref for x in d]
    conv = [measured_size(0.0, x) / ref for x in d]
    ok = check("조리개가 뒤 초점면: 크기가 δ 에 무관", tele, 1.0, 1e-9)
    ok &= check("조리개가 렌즈: 크기 ∝ 2f/(2f − δ)", conv, 2 * f / (2 * f - d), 1e-12)
    return ok


def verify_brightness():
    print("G. 밝기 정리 (돋보기로 햇빛 모으기)")
    alpha = 4.65e-3                 # 태양의 각반지름 [rad]
    B = 1.0                         # 태양 복사휘도 (상대값)
    ok = True
    for D, f in ((50.0, 100.0), (50.0, 50.0)):
        flux = B * np.pi * (D / 2) ** 2 * np.pi * alpha**2   # 렌즈 넓이 × 태양이 차지하는 입체각
        E = flux / (np.pi * (f * alpha) ** 2)                 # 상의 넓이로 나눈다
        ok &= check(f"D/f = {D / f:.1f}: 상의 조도 ↔ π B sin²θ′ (근축)", E, np.pi * B * (D / (2 * f)) ** 2, 1e-12)
        print(f"      상의 조도 / 태양 표면 방출도(πB) = {E / (np.pi * B):.3f}")
    return ok


def main():
    results = [verify_hecht(), verify_pupils(), verify_lagrange(), verify_cos4(), verify_vignetting(),
               verify_telecentric(), verify_brightness()]
    print("\n전체:", "통과" if all(results) else "실패 있음")


if __name__ == "__main__":
    main()

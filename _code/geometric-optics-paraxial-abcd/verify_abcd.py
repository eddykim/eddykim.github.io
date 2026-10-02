"""abcd.py 를 독립된 경로들과 교차검증한다.

실행: python verify_abcd.py

  A. 행렬로 구한 초점·주평면  ↔  1편 실광선 추적의 극한, Born & Wolf 4.4.3 두꺼운 렌즈 닫힌 식(25)·(28)
  B. 행렬                     ↔  Optiland paraxial (오라클 환경이 있을 때만)
  C. 행렬식                   ↔  감소각 규약 1, 보통 각 규약 n_in/n_out (Saleh & Teich)
  D. 라그랑주(스미스-헬름홀츠) 불변량 ↔ 면마다 보존, 켤레면에서 횡배율 × 각배율 = 1 (Goodman B-16)
  E. B = 0 결상               ↔  주평면 기준 가우스 식, 초점 기준 뉴턴 식, 종배율 = 횡배율²
  F. 초점면 → 초점면 행렬      ↔  A = D = 0 (Goodman B.4)
  G. 연속 GRIN 행렬           ↔  얇은 판 곱, 1편 광선 방정식(작은 각), Saleh & Teich GRIN 렌즈 식
  H. 포물선형 GRIN 의 실제 반주기 ↔ (π/α) cos θ0
  I. 주기 렌즈열의 안정 조건   ↔  2000 단 반복 적용, 주기 6·4 (Saleh & Teich 예제 1.4-1)
  J. 망원 조건                 ↔  Born & Wolf 식 (34) t = 2nr/(n-1)
  K. 두 얇은 렌즈의 합성 초점   ↔  Born & Wolf 식 (39), 망원 렌즈의 주평면이 렌즈 앞으로 나간다
  L. 위상공간 넓이             ↔  광선 다발의 테두리를 행렬로 옮긴 다각형 넓이
"""
import json
import os
import subprocess

import numpy as np

from abcd import (cardinal_points, grin, grin_sliced, image_distance, iterate, periodic_cell,
                  prop, stability, surface, system, thin_lens, to_angle_convention)
from rays import Singlet, n2_parabolic, sample, trace_ray

np.seterr(invalid="ignore")


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<60s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


LENS = Singlet()
N, NA, T, R1, R2 = LENS.n_glass, LENS.n_air, LENS.thickness, LENS.r1, LENS.r2
M = system(surface(NA, N, R1), prop(T, N), surface(N, NA, R2))
CP = cardinal_points(M, NA, NA)


def verify_thick_lens():
    print("A. 주요점 ↔ 실광선 극한, Born & Wolf 닫힌 식")
    ok = check("뒤 초점 [mm] ↔ 광축 근처 실광선(h = 1e-4 mm)",
               T + CP.z_back_focus, LENS.axis_crossing([1e-4])[0], 1e-6)
    # Born & Wolf 4.4.3: n0 = n2 = NA, n1 = N
    D = (N - NA) * (NA - N) * T - N * ((NA - N) * R1 + (N - NA) * R2)
    # 식 (25): f = -n0 n1 r1 r2 / D. B&W 는 데카르트 부호 규약이라 같은 매질에서 f' = -f 다
    f_bw = -NA * N * R1 * R2 / D
    d_bw = -NA * (NA - N) * R1 * T / D          # 식 (28) 의 d  (앞 꼭짓점 → H)
    dp_bw = NA * (N - NA) * R2 * T / D          # 식 (28) 의 d' (뒤 꼭짓점 → H')
    ok &= check("초점거리 [mm] ↔ B&W (25) f", CP.f_obj, f_bw, 1e-9)
    ok &= check("앞 주평면 H [mm] ↔ B&W (28) d", CP.z_front_principal, d_bw, 1e-9)
    ok &= check("뒤 주평면 H' [mm] ↔ B&W (28) d'", CP.z_back_principal, dp_bw, 1e-9)
    print(f"      EFL {CP.f_obj:.4f} mm, 앞 초점 {CP.z_front_focus:.4f}, H {CP.z_front_principal:.4f}, "
          f"H' {CP.z_back_principal:.4f}, 뒤 초점 {CP.z_back_focus:.4f} (각 꼭짓점 기준)")
    return ok


OPTILAND_SNIPPET = r"""
import json, sys, numpy as np
from optiland import optic
from optiland.materials import IdealMaterial
s = lambda v: float(np.asarray(v).reshape(-1)[0])
out = {}
for na in (1.0002778, 1.0):
    air = IdealMaterial(n=na)
    L = optic.Optic()
    L.surfaces.add(index=0, radius=float("inf"), thickness=float("inf"), material=air)
    L.surfaces.add(index=1, radius=1000.0, thickness=100.0, material="N-BK7", is_stop=True)
    L.surfaces.add(index=2, radius=-1000.0, thickness=1000.0, material=air)
    L.surfaces.add(index=3)
    L.set_aperture(aperture_type="EPD", value=20.0); L.fields.set_type("angle"); L.fields.add(y=0.0)
    L.wavelengths.add(value=0.75, is_primary=True)
    p = L.paraxial
    out[str(na)] = dict(f1=s(p.f1()), f2=s(p.f2()), F1=s(p.F1()), F2=s(p.F2()), P1=s(p.P1()), P2=s(p.P2()))
print(json.dumps(out))
"""

ORACLE_PYTHON = os.environ.get("OPTILAND_PYTHON", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
    "repo_opticore", "reference_libs", ".venv-oracle", "bin", "python"))


def verify_optiland():
    print("B. 행렬 ↔ Optiland paraxial")
    py = os.path.abspath(ORACLE_PYTHON)
    if not os.path.exists(py):
        print(f"  [SKIP] Optiland 오라클 환경이 없다 ({py})")
        return True
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    res = json.loads(subprocess.run([py, "-c", OPTILAND_SNIPPET], capture_output=True, text=True,
                                    env=env, check=True).stdout.strip().splitlines()[-1])
    o = res["1.0002778"]
    ok = check("물체 쪽 초점거리 f1 [mm]", -o["f1"], CP.f_obj, 1e-6)
    ok &= check("앞 초점 F1 (앞 꼭짓점 기준) [mm]", o["F1"], CP.z_front_focus, 1e-6)
    ok &= check("앞 주점 P1 (앞 꼭짓점 기준) [mm]", o["P1"], CP.z_front_principal, 1e-6)
    # 상 쪽: Optiland 의 f2 는 1/Φ 로, 상 공간 굴절률 n' 을 곱하지 않는다
    ok &= check("Optiland f2 × n_air ↔ 상 쪽 초점거리 n'/Φ", o["f2"] * NA, CP.f_img, 1e-6)
    # 공기를 정확히 1 로 두면 두 규약이 같아진다
    m1 = system(surface(1.0, N, R1), prop(T, N), surface(N, 1.0, R2))
    c1 = cardinal_points(m1, 1.0, 1.0)
    o1 = res["1.0"]
    ok &= check("n_air = 1: f2 ↔ 행렬", o1["f2"], c1.f_img, 1e-6)
    # Optiland 의 F2, P2 는 상면(마지막 면에서 1000 mm 뒤) 기준이다
    ok &= check("n_air = 1: F2 (상면 기준) ↔ 행렬", o1["F2"], c1.z_back_focus - 1000.0, 1e-6)
    ok &= check("n_air = 1: P2 (상면 기준) ↔ 행렬", o1["P2"], c1.z_back_principal - 1000.0, 1e-6)
    print(f"      n_air = 1.0002778 에서 Optiland f2 = {o['f2']:.6f} (행렬 n'/Φ = {CP.f_img:.6f})")
    return ok


def verify_determinants():
    print("C. 행렬식")
    mats = [prop(123.0, N), surface(NA, N, R1), surface(N, NA, R2), thin_lens(0.01), M,
            grin(1.5, 0.5, 2.7)]
    ok = check("감소각 규약: 모든 행렬의 det = 1", [np.linalg.det(m) for m in mats], 1.0, 1e-12)
    ms = to_angle_convention(surface(NA, N, R1), NA, N)
    ok &= check("보통 각 규약: 구면 하나의 det = n1/n2", np.linalg.det(ms), NA / N, 1e-12)
    # 보통 각 규약의 구면 행렬을 Saleh & Teich 식 (1.4-6) 꼴로 직접 적어 비교한다
    saleh = np.array([[1.0, 0.0], [-(N - NA) / (N * R1), NA / N]])
    ok &= check("보통 각 규약 ↔ Saleh & Teich 식 (1.4-6) 꼴", ms, saleh, 1e-15)
    return ok


def verify_lagrange():
    print("D. 라그랑주 불변량")
    r1, r2 = np.array([1.0, 0.0]), np.array([0.0, 1e-3])
    elems = [surface(NA, N, R1), prop(T, N), surface(N, NA, R2), prop(500.0, NA)]
    inv = []
    for e in [np.eye(2)] + elems:
        r1, r2 = e @ r1, e @ r2
        inv.append(r1[0] * r2[1] - r2[0] * r1[1])
    ok = check("두 광선의 y1 v2 - y2 v1 이 면마다 같다", inv, inv[0], 1e-15)
    s_img, tot = image_distance(M, 3000.0, NA, NA)
    a, b, c, d = tot.ravel()
    ok &= check("켤레면에서 B = 0", b, 0.0, 1e-9)
    ok &= check("횡배율 A × 각배율 D = 1", a * d, 1.0, 1e-12)
    print(f"      물체 3000 mm → 상 {s_img:.3f} mm (뒤 꼭짓점 기준), 횡배율 {a:.5f}")
    return ok


def verify_imaging():
    print("E. 결상식")
    ok = True
    for s in (1500.0, 3000.0, 10000.0):
        s_img, _ = image_distance(M, s, NA, NA)
        # 주평면에서 잰 거리: 물체는 H 앞, 상은 H' 뒤
        so, si = s + CP.z_front_principal, s_img - CP.z_back_principal
        ok &= check(f"가우스 식 1/so + 1/si = 1/f (물체 {s:.0f} mm)", 1 / so + 1 / si, 1 / CP.f_obj, 1e-12)
        zo, zi = s + CP.z_front_focus, s_img - CP.z_back_focus
        ok &= check(f"뉴턴 식 z z' = f f' (물체 {s:.0f} mm)", zo * zi, CP.f_obj * CP.f_img, 1e-6)
    # 종배율: 물체를 광축 방향으로 ds 움직일 때 상이 움직이는 거리. 같은 매질이면 횡배율의 제곱
    s0, ds = 3000.0, 1e-3
    s_a, tot = image_distance(M, s0, NA, NA)
    s_b, _ = image_distance(M, s0 + ds, NA, NA)
    m_long = -(s_b - s_a) / ds          # 물체가 멀어지면(ds>0) 상은 렌즈 쪽으로 온다
    ok &= check("종배율 = 횡배율²", m_long, tot[0, 0] ** 2, 1e-6)
    print(f"      물체 3000 mm: 횡배율 {tot[0, 0]:.4f}, 종배율 {m_long:.4f}")
    return ok


def verify_focal_planes():
    print("F. 초점면 → 초점면")
    ff = system(prop(-CP.z_front_focus, NA), M, prop(CP.z_back_focus, NA))
    ok = check("A = D = 0", [ff[0, 0], ff[1, 1]], 0.0, 1e-12)
    ok &= check("B·C = -1", ff[0, 1] * ff[1, 0], -1.0, 1e-12)
    return ok


def verify_grin():
    print("G. 연속 GRIN 행렬")
    n0, al, length = 1.5, 0.5, 3.0
    g = grin(n0, al, length)
    ok = check("↔ 얇은 판 4000 장의 곱", g, grin_sliced(n0, al, length, 4000), 1e-7)
    # 근축 행렬과 광선 방정식의 차이는 근축 근사의 오차이므로 진폭의 제곱으로 줄어야 한다
    rel = []
    for scale in (1.0, 0.1):
        th, y0 = np.deg2rad(0.5 * scale), 0.1 * scale
        sol = trace_ray(n2_parabolic(n0, al), (0, y0), th, 3 * length, stop=lambda x, y: x - length)
        y_ode = sol.y_events[0][0][1]
        y_mat = (g @ np.array([y0, n0 * np.sin(th)]))[0]
        rel.append(abs(y_ode - y_mat) / abs(y_mat))
    print(f"      광선 방정식과의 상대 차이: 진폭 1 → {rel[0]:.2e}, 진폭 0.1 → {rel[1]:.2e}")
    ok &= check("↔ 1편 광선 방정식: 진폭 1/10 에 차이 1/100 (근축 오차)", rel[0] / rel[1], 100.0, 5.0)
    # Saleh & Teich 연습 1.3-1: 길이 d 의 GRIN 판은 초점거리 1/(n0 α sin αd) 인 렌즈
    c = cardinal_points(g, 1.0, 1.0)
    ok &= check("초점거리 ↔ 1/(n0 α sin αd)", c.f_img, 1 / (n0 * al * np.sin(al * length)), 1e-12)
    ok &= check("주점 ↔ 끝면에서 tan(αd/2)/(n0 α) 안쪽", -c.z_back_principal,
                np.tan(al * length / 2) / (n0 * al), 1e-12)
    return ok


def verify_grin_half_pitch():
    print("H. 포물선형 GRIN 의 실제 반주기")
    medium = n2_parabolic(1.5, 0.5)
    degs = np.array([1, 10, 30, 45, 60.0])
    got = []
    for deg in degs:
        sol = trace_ray(medium, (0, 0), np.deg2rad(deg), 40, stop=lambda x, y: y if x > 0.5 else 1.0)
        got.append(sol.y[0][-1])
    ok = check("↔ (π/α) cos θ0", got, np.pi / 0.5 * np.cos(np.deg2rad(degs)), 1e-8)
    for d, g in zip(degs, got):
        print(f"      θ0 = {d:4.0f}°  반주기 {g:.4f}  (근축 {np.pi / 0.5:.4f}, {100 * (g / (np.pi / 0.5) - 1):+.1f}%)")
    return ok


def verify_periodic():
    print("I. 주기 렌즈열")
    ok = True
    for r in np.linspace(0.05, 5.0, 100):
        b, _ = stability(periodic_cell(1.0, r))
        traj = iterate(periodic_cell(1.0, r), [1.0, 0.3], 2000)
        bounded = np.max(np.abs(traj[:, 0])) < 1e3
        stable_pred = abs(b) <= 1 and not np.isclose(abs(b), 1)
        if bounded != stable_pred and not np.isclose(r, 4.0, atol=0.06):
            ok = False
            print(f"  [FAIL] d/f = {r:.3f}: 예측 {stable_pred}, 반복 {bounded}")
    print(f"  [{'OK  ' if ok else 'FAIL'}] 0 < d/f < 4 에서만 2000 단 동안 유계 (100 개 간격)")
    ok &= check("d = f: M^6 = I", np.linalg.matrix_power(periodic_cell(1.0, 1.0), 6), np.eye(2), 1e-12)
    ok &= check("d = 2f: M^4 = I", np.linalg.matrix_power(periodic_cell(1.0, 2.0), 4), np.eye(2), 1e-12)
    return ok


def verify_telescopic():
    print("J. 망원 조건")
    nr = N / NA
    t_tel = 2 * nr * R1 / (nr - 1)
    m = system(surface(NA, N, R1), prop(t_tel, N), surface(N, NA, R2))
    ok = check("t = 2nr/(n-1) 에서 C = 0", m[1, 0], 0.0, 1e-15)
    print(f"      이 렌즈를 두께 {t_tel:.1f} mm 로 만들면 초점이 무한대로 간다")
    return ok


def verify_telephoto():
    print("K. 두 얇은 렌즈 (망원 렌즈)")
    f1, f2, d = 100.0, -50.0, 60.0
    m = system(thin_lens(1 / f1), prop(d), thin_lens(1 / f2))
    c = cardinal_points(m)
    f_bw = 1 / (1 / f1 + 1 / f2 - d / (f1 * f2))
    ok = check("합성 초점거리 ↔ B&W (39)", c.f_img, f_bw, 1e-9)
    length = d + c.z_back_focus               # 첫 렌즈에서 초점까지 실제 길이
    hp_from_first = d + c.z_back_principal     # 첫 렌즈 기준 H' 위치
    ok &= check("H' 가 첫 렌즈 앞에 있다 (부호)", np.sign(hp_from_first), -1.0, 0.5)
    print(f"      f = {c.f_img:.1f} mm, 첫 렌즈→초점 {length:.1f} mm (길이/초점거리 {length / c.f_img:.2f}), "
          f"H' 는 첫 렌즈 앞 {-hp_from_first:.1f} mm")
    return ok


def verify_phase_space_area():
    print("L. 위상공간 넓이")
    side = np.linspace(-1, 1, 201)
    edge = np.r_[np.c_[side, -np.ones_like(side)], np.c_[np.ones_like(side), side],
                 np.c_[side[::-1], np.ones_like(side)], np.c_[-np.ones_like(side), side[::-1]]]
    box = edge * np.array([50.0, 0.01])         # 높이 ±50 mm, 감소각 ±0.01
    full = system(prop(300.0, NA), M, prop(700.0, NA))

    def area(p):
        x, y = p[:, 0], p[:, 1]
        return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

    return check("1편 렌즈 앞뒤로 광선 다발의 넓이 보존 (상대)", area(box @ full.T) / area(box), 1.0, 1e-12)


def main():
    results = [verify_thick_lens(), verify_optiland(), verify_determinants(), verify_lagrange(),
               verify_imaging(), verify_focal_planes(), verify_grin(), verify_grin_half_pitch(),
               verify_periodic(), verify_telescopic(), verify_telephoto(), verify_phase_space_area()]
    print("\n전체:", "통과" if all(results) else "실패 있음")


if __name__ == "__main__":
    main()

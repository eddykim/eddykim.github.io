"""aberr.py 를 독립된 경로들과 교차검증한다.

실행: python verify_aberr.py

  A. 실광선 추적기             ↔  1편의 2차원 추적기 rays.py (종구면수차 15.845 mm)
  B. 자이델 합 S_I…S_V          ↔  Optiland SeidelAberrations (평볼록 렌즈, 쿡 삼중렌즈)
  C. 자이델 계수 → 파면수차      ↔  같은 렌즈의 실광선 OPD 에 다항식을 맞춘 값 (이 모듈의 추적기)
  D. 실광선 OPD 지도             ↔  opticore sample_pupil_opd (다른 사람이 짠 추적기·기준 구면)
  E. 광선수차 = 파면의 기울기     ↔  실광선 횡수차와 −(1/n′NA) dW/dρ
  F. 렌즈 휘기                   ↔  얇은 렌즈의 닫힌 식: 코마 0 은 B&W 5.6 (n = 1.5 에서 r1 = 5f/9, r2 = −5f),
                                   구면수차 최소는 q = 2(n² − 1)/(n + 2)
  G. 페츠발 곡률과 3:1           ↔  Hecht 식 (6.49) Δx = y²/2 · Σ 1/(n f), 실광선 T·S 초점
  H. 구면수차의 초점             ↔  최소 착란원 3/4, 회절 초점 1/2 (B&W 9.3)
  I. 허용치 표                   ↔  B&W 표 9.3 (0.94λ, 0.60λ, 0.35λ 에서 스트렐 ≥ 0.8), 마레샬 λ/14
  J. 색지움 렌즈                 ↔  얇은 렌즈의 2차 스펙트럼 f·(P1 − P2)/(V1 − V2)
  K. 포물면 거울의 코마          ↔  3차 식: 구결 코마 θf/(16N²), 자오 코마 3θf/(16N²)
  L. 커버글라스                  ↔  평판의 3차 구면수차 W040 = t (n² − 1) NA⁴ / (8 n³) (낮은 NA)
"""
import json
import os
import subprocess

import numpy as np

from aberr import (AIR, LINE_C, LINE_D, LINE_F, Surface, abbe_number, axis_crossing, fit_seidel_terms, opd,
                   paraxial_trace, pupil_grid, seidel_sums, setup, spot, strehl_center, wave_coefficients)
from rays import Singlet
from systems import (BK7, F2, WL, achromat, center_stop_z, bent_lens, cooke_triplet, plano_convex, TRIPLET_RX, TRIPLET_STOP)

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "repo_opticore")
LAM = WL * 1e-3


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<62s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def run(py, code, *args):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
    out = subprocess.run([py, "-c", code, *map(str, args)], capture_output=True, text=True,
                         env=env, check=True)
    return json.loads(out.stdout.strip().splitlines()[-1])


# ------------------------------------------------------------------ A

def verify_tracer():
    print("A. 3차원 실광선 추적기 ↔ 1편 rays.py")
    s = Singlet()
    na = s.n_air
    lens = [Surface(0.0, s.r1, lambda wl: s.n_glass), Surface(s.thickness, s.r2, lambda wl: na)]
    hs = np.array([1e-3, 25.0, 50.0, 75.0, 100.0])
    z_new = axis_crossing(lens, s.wavelength_um, hs, z_start=-300.0, n0=na)
    z_old = s.axis_crossing(hs)
    ok = check("광축을 지나는 z (높이 0.001~100 mm)", z_new, z_old, 1e-8)
    ok &= check("종구면수차 = 15.845 mm (1편)", z_new[0] - z_new[-1], 15.845, 5e-4)
    return ok


# ------------------------------------------------------------------ B

OPTILAND = r"""
import json, sys, numpy as np
from optiland import optic
from optiland.materials import IdealMaterial
from optiland.aberrations.seidel import SeidelAberrations
spec = json.loads(sys.argv[1])
L = optic.Optic()
L.surfaces.add(index=0, radius=float("inf"), thickness=float("inf"))
for i, (R, t, n, stop) in enumerate(spec["surfaces"]):
    kw = dict(index=i + 1, radius=float(R), thickness=float(t), is_stop=bool(stop))
    if n != 1.0:
        kw["material"] = IdealMaterial(n=n)
    L.surfaces.add(**kw)
L.surfaces.add(index=len(spec["surfaces"]) + 1)
L.set_aperture(aperture_type="EPD", value=spec["epd"])
L.fields.set_type("angle"); L.fields.add(y=0); L.fields.add(y=spec["field"])
L.wavelengths.add(value=spec["wl"], is_primary=True)
print(json.dumps([float(v) for v in np.asarray(SeidelAberrations(L).seidels()).ravel()]))
"""


def optiland_spec(surfs, epd, field_deg, wl):
    """면 목록 → Optiland 처방 (반지름, 다음 면까지 두께, 뒤 매질 굴절률, 조리개)."""
    st = setup(surfs, wl, 0.0)
    rows = []
    for i, s in enumerate(surfs):
        nxt = surfs[i + 1].z if i + 1 < len(surfs) else st.z_img
        rows.append([s.R if np.isfinite(s.R) else 1e30, nxt - s.z, s.medium(wl), s.stop])
    return json.dumps(dict(surfaces=rows, epd=epd, field=field_deg, wl=wl)).replace("Infinity", "1e30")


def verify_seidel_optiland():
    print("B. 자이델 합 ↔ Optiland SeidelAberrations")
    py = os.path.join(REPO, "reference_libs", ".venv-oracle", "bin", "python")
    if not os.path.exists(py):
        print(f"  [SKIP] Optiland 오라클 환경이 없다 ({py})")
        return True
    ok = True
    cases = [("평볼록 F/10, 시야 5°", plano_convex(D=10.0), 10.0, 5.0),
             ("뒤집은 평볼록 + 곡률 중심 조리개, 시야 10°", plano_convex(D=10.0, stop_z=center_stop_z(), flip=True), None, 10.0),
             ("쿡 삼중렌즈 F/5, 시야 20°", cooke_triplet(), 10.0, 20.0)]
    for name, surfs, epd, fdeg in cases:
        st = setup(surfs, WL, np.radians(fdeg))
        epd = 2 * st.r_ep if epd is None else epd
        S = seidel_sums(surfs, WL, st)
        o = np.array(run(py, OPTILAND, optiland_spec(surfs, epd, fdeg, WL)))
        # Optiland 은 전체 부호가 반대인 규약을 쓴다 (S_I > 0 이 덜 보정된 양렌즈인 Welford 규약과 반대)
        ok &= check(f"{name}: S_I…S_V (상대)", (S + o) / np.abs(S).max(), 0.0, 1e-9)
    return ok


# ------------------------------------------------------------------ C

def verify_seidel_vs_real():
    print("C. 자이델 계수 ↔ 실광선 OPD 다항식 맞춤 (작은 개구·시야에서 3차가 지배)")
    lens = plano_convex(D=4.0)
    st = setup(lens, WL, np.radians(2.0))
    W3 = wave_coefficients(seidel_sums(lens, WL, st))
    X, Y, M = pupil_grid(41)
    W, okm = opd(lens, WL, st, X[M], Y[M], z_image=st.z_img)
    f = fit_seidel_terms(W, X[M], Y[M])
    ok = True
    for k_fit, k3 in [("W040", "W040"), ("W131", "W131"), ("W222", "W222"), ("defocus", "W220")]:
        ok &= check(f"{k3}: 3차 {W3[k3] / LAM:+.5f} λ ↔ 실광선 {f[k_fit] / LAM:+.5f} λ (상대)",
                    (f[k_fit] - W3[k3]) / W3[k3], 0.0, 0.02)
    return ok


# ------------------------------------------------------------------ D

OPTICORE = r"""
import json, sys, numpy as np
from opticore.geometry import Plane, Sphere
from opticore.interactions import Refractive
from opticore.materials import ConstantMaterial
from opticore.surface import Surface
from opticore.system import OpticalSystem
from opticore.pupil import Field
from opticore.wavefront.pupil import sample_pupil_opd
n, R, t, wl, D, fdeg, zimg = map(float, sys.argv[1:8])
g, a = ConstantMaterial(n), ConstantMaterial(1.0)
s = OpticalSystem(ambient=a)
s.add(Surface(Sphere(R), g, Refractive(), position=(0, 0, 0.0)))
s.add(Surface(Plane(), a, Refractive(), position=(0, 0, t)))
f = sample_pupil_opd(s, wl, entrance_z=-1.0, pupil_radius_mm=D / 2, focus_z=zimg, n_side=41,
                     field=Field(y=fdeg))
m = f.mask
print(json.dumps(dict(px=f.px[m].tolist(), py=f.py[m].tolist(), opd=f.opd_um[m].tolist())))
"""


def verify_opticore_opd():
    print("D. 실광선 OPD 지도 ↔ opticore sample_pupil_opd")
    py = os.path.join(REPO, "opticore", ".venv", "bin", "python")
    if not os.path.exists(py):
        print(f"  [SKIP] opticore 환경이 없다 ({py})")
        return True
    ok = True
    for D, fdeg in [(20.0, 0.0), (10.0, 5.0)]:
        lens = plano_convex(D=D, stop_z=-1.0)       # opticore 의 입사동 평면(z = −1)과 같은 자리에 조리개
        st = setup(lens, WL, np.radians(fdeg))
        o = run(py, OPTICORE, BK7(WL), lens[1].R, lens[2].z, WL, D, fdeg, st.z_img)
        px, py_, W_oc = np.array(o["px"]), np.array(o["py"]), np.array(o["opd"]) * 1e-3
        W, _ = opd(lens, WL, st, px, py_, z_image=st.z_img)
        f_me, f_oc = fit_seidel_terms(W, px, py_), fit_seidel_terms(W_oc, px, py_)
        keys = ["W040", "W131", "W222", "defocus"] if fdeg else ["W040", "defocus"]
        for k in keys:
            ok &= check(f"F/{100 / D:.0f}, 시야 {fdeg}°: {k} {f_me[k] / LAM:+.4f} λ ↔ opticore {f_oc[k] / LAM:+.4f} λ",
                        (f_me[k] - f_oc[k]) / LAM, 0.0, 2e-3 if fdeg == 0 else 2e-2)
        if fdeg == 0:
            # 기준 구면이 지나는 점이 다르다 (opticore: 마지막 면 위의 주광선 점, 여기: 출사동). 같은 중심의 두 구면
            # 사이를 광선이 비스듬히 지나는 만큼만 어긋나므로 수차의 제곱 정도로 작다.
            ok &= check("축상 OPD 지도 전체 (λ, W 최대 4.6λ)", (W - W_oc) / LAM, 0.0, 5e-3)
    return ok


# ------------------------------------------------------------------ E

def verify_ray_wave():
    print("E. 광선수차 = 파면의 기울기")
    lens = plano_convex(D=20.0)
    st = setup(lens, WL, 0.0)
    rho = np.linspace(-1, 1, 401)
    W, _ = opd(lens, WL, st, np.zeros_like(rho), rho)
    x, y, _ = spot(lens, WL, st, np.zeros_like(rho), rho, st.z_img)
    eps = -np.gradient(W, rho) / (st.n_img * st.sin_u_img)
    inner = slice(5, -5)
    return check("F/5 축상: 횡수차 ↔ −(1/n′NA) dW/dρ (최대 횡수차 대비)",
                 (y[inner] - eps[inner]) / np.abs(y).max(), 0.0, 0.02)


# ------------------------------------------------------------------ F

def verify_bending():
    print("F. 렌즈 휘기 ↔ 얇은 렌즈의 닫힌 식")
    ok = True
    for n in [1.5, float(BK7(WL))]:
        g = lambda wl, n=n: n
        f = 100.0

        def lens(q):
            P = 1 / ((n - 1) * f)
            c1, c2 = P * (q + 1) / 2, P * (q - 1) / 2
            R1 = np.inf if abs(c1) < 1e-15 else 1 / c1
            R2 = np.inf if abs(c2) < 1e-15 else 1 / c2
            return [Surface(0.0, np.inf, AIR, 1.0, True), Surface(0.0, R1, g, 50), Surface(1e-4, R2, AIR, 50)]

        qs = np.linspace(0.4, 1.2, 1601)
        w040 = [wave_coefficients(seidel_sums(lens(q), WL, setup(lens(q), WL, np.radians(1.0))))["W040"] for q in qs]
        w131 = [wave_coefficients(seidel_sums(lens(q), WL, setup(lens(q), WL, np.radians(1.0))))["W131"] for q in qs]
        q_min = qs[np.argmin(w040)]
        q_c0 = np.interp(0.0, w131, qs)
        ok &= check(f"n = {n:.4f}: 구면수차 최소 q ↔ 2(n²−1)/(n+2) = {2 * (n**2 - 1) / (n + 2):.4f}",
                    q_min, 2 * (n**2 - 1) / (n + 2), 1e-3)
        ok &= check(f"n = {n:.4f}: 코마 0 인 q ↔ (2n+1)(n−1)/(n+1) = {(2 * n + 1) * (n - 1) / (n + 1):.4f}",
                    q_c0, (2 * n + 1) * (n - 1) / (n + 1), 1e-3)
        if n == 1.5:
            P = 1 / ((n - 1) * f)
            c1, c2 = P * (q_c0 + 1) / 2, P * (q_c0 - 1) / 2
            ok &= check("n = 1.5: 코마 0 인 렌즈 r1 = 5f/9, r2 = −5f (B&W 5.6 식 19)", [1 / c1, 1 / c2],
                        [5 * f / 9, -5 * f], 0.05)
    return ok


# ------------------------------------------------------------------ G

def verify_petzval():
    print("G. 페츠발 곡률과 3:1")
    from generate_figures import focal_tangential_sagittal
    ok = True
    lens = plano_convex(D=2.0)
    n = BK7(WL)
    st0 = setup(lens, WL, 0.0)
    m = paraxial_trace(lens, *st0.marg, st0.z0, WL)
    f = st0.marg[0] / -m[-1, 2]
    for deg in [1.0, 2.0]:
        st = setup(lens, WL, np.radians(deg))
        Wc = wave_coefficients(seidel_sums(lens, WL, st))
        u2 = st.n_img * st.sin_u_img**2
        P = -2 * Wc["W220P"] / u2
        y = -np.tan(np.radians(deg)) * f          # 근축 상 높이의 크기
        ok &= check(f"{deg}°: 페츠발 면 위치 ↔ Hecht (6.49) −y²/(2nf) (상대)", P / (-(y**2) / (2 * n * f)) - 1, 0.0, 0.03)
        (zt, zs), _ = focal_tangential_sagittal(lens, np.radians(deg), dp=0.002)
        ok &= check(f"{deg}°: 실광선 (T − P)/(S − P) = 3", (zt - st0.z_img - P) / (zs - st0.z_img - P), 3.0, 0.05)
    # 조리개를 곡면의 곡률 중심(겉보기 자리)에 두면 주광선이 곡면에 수직: 코마·비점수차가 사라진다 (Hecht 5.7, 슈미트 카메라)
    lens = plano_convex(D=10.0, stop_z=center_stop_z(), flip=True)
    st = setup(lens, WL, np.radians(10.0))
    Wc = wave_coefficients(seidel_sums(lens, WL, st))
    ok &= check("곡률 중심 조리개: W131, W222 = 0 (λ)", [Wc["W131"] / LAM, Wc["W222"] / LAM], 0.0, 1e-9)
    (zt, zs), _ = focal_tangential_sagittal(lens, np.radians(10.0), dp=0.002)
    zP = setup(lens, WL, 0.0).z_img - 2 * Wc["W220P"] / (st.n_img * st.sin_u_img**2)
    ok &= check("곡률 중심 조리개, 10°: 실광선 T·S 초점이 페츠발 면 위 (mm)", [zt - zP, zs - zP], 0.0, 0.03)
    return ok


# ------------------------------------------------------------------ H

def verify_focus_positions():
    print("H. 구면수차의 초점: 최소 착란원 3/4, 회절 초점 1/2")
    ok = True
    rho = np.linspace(0, 1, 20001)
    fr = np.linspace(0.0, 1.2, 12001)
    blur = [np.max(np.abs(4 * rho**3 - 4 * f_ * rho)) for f_ in fr]   # W040 = 1, W020 = −2·frac
    ok &= check("순수 3차: 번짐 최소 위치 = 3/4", fr[np.argmin(blur)], 0.75, 2e-4)
    X, Y, M = pupil_grid(201)
    r2 = X**2 + Y**2
    fr2 = np.linspace(0.3, 0.7, 401)
    S = [strehl_center(0.3 * r2**2 - 0.6 * f_ * r2, M) for f_ in fr2]
    ok &= check("순수 3차 (W040 = 0.3λ): 스트렐 최대 위치 = 1/2", fr2[np.argmax(S)], 0.5, 2e-3)
    lens = plano_convex(D=40.0)
    zp = setup(lens, WL, 0.0).z_img
    zm = axis_crossing(lens, WL, [20.0])[0]
    hs = np.linspace(0, 20, 4001)[1:]
    from aberr import trace_real, to_plane
    P = np.stack([np.zeros_like(hs), hs, np.full_like(hs, -10.0)], 1)
    D = np.tile([0.0, 0.0, 1.0], (len(hs), 1))
    Pf, Df, _, _, _ = trace_real(lens, P, D, WL)
    zs = np.linspace(zm, zp, 4001)
    radius = [np.max(np.abs(to_plane(Pf, Df, z)[:, 1])) for z in zs]
    z_lc = zs[np.argmin(radius)]
    ok &= check("F/2.5 실광선: 최소 착란원 (근축 초점에서, LSA 비) ≈ 3/4", (zp - z_lc) / (zp - zm), 0.75, 0.02)
    return ok


# ------------------------------------------------------------------ I

def verify_tolerances():
    print("I. B&W 표 9.3 허용치와 마레샬 λ/14")
    from generate_figures import best_strehl
    X, Y, M = pupil_grid(201)
    ok = True
    for name, Wf, A, rms_formula in [("구면수차 A'040 = 0.94λ", lambda x, y: (x**2 + y**2) ** 2, 0.94, 0.94 / (6 * np.sqrt(5))),
                                    ("코마 A'031 = 0.60λ", lambda x, y: y * (x**2 + y**2), 0.60, 0.60 / (6 * np.sqrt(2))),
                                    ("비점수차 A'022 = 0.35λ", lambda x, y: y**2, 0.35, 0.35 / (2 * np.sqrt(6)))]:
        S, r = best_strehl(Wf, X, Y, M, A)
        ok &= check(f"{name}: 균형 후 RMS ↔ 닫힌 식 {rms_formula:.4f}λ", r, rms_formula, 5e-4)
        ok &= check(f"{name}: RMS ≈ λ/14", r, 1 / 14, 1.5e-3)
        ok &= check(f"{name}: 스트렐 0.80~0.83", S, 0.815, 0.016)
    return ok


# ------------------------------------------------------------------ J

def verify_achromat():
    print("J. 색지움 렌즈의 2차 스펙트럼 ↔ 얇은 렌즈 식")
    lens, info = achromat()

    def bf(wl):
        m = paraxial_trace(lens, 1.0, 0.0, 0.0, wl)
        return lens[-1].z - m[-1, 0] / m[-1, 2]

    m = paraxial_trace(lens, 1.0, 0.0, 0.0, LINE_D)
    f = -1.0 / m[-1, 2]
    P = lambda g: (g(LINE_D) - g(LINE_C)) / (g(LINE_F) - g(LINE_C))
    V1, V2 = abbe_number(BK7), abbe_number(F2)
    pred = f * (P(BK7) - P(F2)) / (V1 - V2)
    ok = check("F·C 선 뒤초점 일치 (µm)", (bf(LINE_F) - bf(LINE_C)) * 1e3, 0.0, 1e-6)
    ok &= check(f"d 선 − C 선 초점 {1e3 * (bf(LINE_D) - bf(LINE_C)):+.2f} µm ↔ f(P1−P2)/(V1−V2) = {1e3 * pred:+.2f} µm (상대)",
                (bf(LINE_D) - bf(LINE_C)) / -pred - 1, 0.0, 0.1)
    print(f"      초점거리 {f:.3f} mm, 2차 스펙트럼 ≈ f/{f / abs(pred):.0f}")
    return ok


# ------------------------------------------------------------------ K

def verify_paraboloid():
    print("K. 포물면 거울의 코마 ↔ 3차 식")
    from generate_figures import paraboloid_spot
    f, D, th = 300.0, 100.0, np.radians(0.05)
    N = f / D
    spots = paraboloid_spot(f, D, th)
    rho, _, x, y = spots[-1]
    chief = -f * np.tan(th)
    ok = check("구결 코마 (가장 바깥 고리 반폭) ↔ θf/(16N²) (상대)", (np.ptp(x) / 2) / (th * f / (16 * N**2)) - 1, 0.0, 0.02)
    ok &= check("자오 코마 (꼬리 길이) ↔ 3θf/(16N²) (상대)", np.max(np.abs(y - chief)) / (3 * th * f / (16 * N**2)) - 1, 0.0, 0.02)
    return ok


# ------------------------------------------------------------------ L

def verify_coverslip():
    print("L. 커버글라스의 구면수차 ↔ 평판의 3차 식")
    from generate_figures import coverslip_strehl
    n, dt, lam = 1.523, 0.010, 0.55e-3
    ok = True
    for NA in [0.1, 0.2]:
        sig, S = coverslip_strehl(NA, dt, n=n, lam_mm=lam)
        w040 = dt * (n**2 - 1) * NA**4 / (8 * n**3)
        ok &= check(f"NA {NA}: 최선 초점 RMS ↔ W040/(6√5) (상대)", sig / (w040 / (6 * np.sqrt(5)) / lam) - 1, 0.0, 0.05)
    return ok


def main():
    results = [f() for f in (verify_tracer, verify_seidel_optiland, verify_seidel_vs_real, verify_opticore_opd,
                             verify_ray_wave, verify_bending, verify_petzval, verify_focus_positions,
                             verify_tolerances, verify_achromat, verify_paraboloid, verify_coverslip)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

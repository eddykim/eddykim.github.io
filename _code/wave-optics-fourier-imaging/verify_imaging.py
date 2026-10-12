"""imaging.py 를 독립된 경로들과 교차검증한다.

실행: python verify_imaging.py

  A. 렌즈의 푸리에 변환          ↔  물체를 앞 초점면에 두면 뒤 초점면의 장 = 물체의 푸리에 변환 (이차 위상 없음, Goodman 5-19)
  B. 결맞는 차단                 ↔  정면 조명이면 주기 λ/NA, 비스듬한 조명이면 λ/(2NA) 까지 (아베)
  C. OTF                         ↔  |PSF|² 의 푸리에 변환 = 동공의 자기상관 = Goodman 6-32 닫힌 식 (opticore 있으면 대조)
  D. 결맞는 가장자리             ↔  경계에서 세기 1/4, 넘침 (1/2 + Si(π)/π)² = 1.187. 결맞지 않으면 1/2 (Goodman 6.5)
  E. 두 점                       ↔  결맞음: 위상차 0 은 골 없음, π/2 는 0.735, π 는 0 (Goodman 그림 6.19).
                                     같은 위상에서 골 0.735 가 되는 간격 0.82 λ/NA (보른과 울프 8.6.3)
  F. 부분 결맞음                 ↔  차단 주파수 (1+σ)NA/λ. σ = 0 이면 결맞는 결상과 같다
  G. 초점 어긋남                 ↔  정사각 동공의 OTF = Λ·sinc[(8W/λ)(f/2f0)(1 − f/2f0)] (Goodman 6-41)
  H. 위상차                      ↔  약한 위상 물체: 0차를 λ/4 앞당기면 1 − 2φ, 늦추면 1 + 2φ, 0차 세기를 1/9 로 줄이면 대비 3배,
                                     명시야·암시야는 φ² 차수 (보른과 울프 8.6.3)
  I. 아베-포터                   ↔  0차를 막으면 열린 비율 m > 1/2 인 그물의 명암이 뒤집힌다: 철사 m², 구멍 (1 − m)²
"""
import os
import subprocess

import numpy as np
from scipy.special import j0, j1, sici

from imaging import (amplitude_psf, coherent_image, dark_field, filter_4f, grid, incoherent_image, lens_fourier_plane,
                     modulation, otf_by_overlap, otf_from_psf, partially_coherent_image, phase_contrast, pupil_grid)

OPTICORE_PY = os.environ.get("OPTICORE_PY", os.path.expanduser(
    "~/Library/CloudStorage/Dropbox-Yjedward.kim/Kim Youngjoon/repository/repo_opticore/opticore/.venv/bin/python"))


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<66s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def verify_lens_ft():
    print("A. 렌즈의 푸리에 변환")
    lam, f, a = 0.5, 2000.0, 20.0
    n, dx = 1024, 1.0
    x, _ = grid(n, dx)
    X, Y = np.meshgrid(x, x)
    U0 = ((np.hypot(X - 6, Y) <= a) * (1 + 0.5 * np.cos(2 * np.pi * X / 8))).astype(complex)   # 중심을 비켜 위상이 생기게
    Uf = lens_fourier_plane(U0, dx, lam, f)
    us = np.array([0.0, 12.0, 25.0, 40.0, 60.0])            # 각 0.03 rad 까지 (더 크면 근축 근사의 위상 오차가 보인다)
    got = Uf[n // 2, n // 2 + us.astype(int)]
    # 직접 푸리에 적분: (1/iλf) ∬ U0 exp(−i2π u ξ/λf)
    sel = np.abs(U0) > 0
    ref = np.array([np.sum(U0[sel] * np.exp(-2j * np.pi * u * X[sel] / (lam * f))) * dx * dx / (1j * lam * f) for u in us])
    ok = check("세기: 렌즈 뒤 초점면 ↔ 물체의 푸리에 변환", np.abs(got) ** 2 / np.abs(ref[0]) ** 2, np.abs(ref) ** 2 / np.abs(ref[0]) ** 2, 0.01)
    ratio = got / ref
    ratio /= ratio[0]
    ok &= check("위상까지: 비율이 상수 (이차 위상 없음)", np.abs(ratio - 1), 0, 0.03)
    return ok


def grating_modulation(p, NA, lam, illum=(0.0, 0.0), sigma=None):
    n, dx = 512, p / 16                       # 창에 주기가 정수 번 (512/16 = 32 주기)
    x, _ = grid(n, dx)
    U = np.tile(1 + 0.5 * np.cos(2 * np.pi * x / p), (n, 1)).astype(complex)
    if sigma is None:
        I = np.abs(coherent_image(U, dx, NA, lam, illum=illum)) ** 2
    else:
        I = partially_coherent_image(U, dx, NA, lam, sigma)
    return modulation(I[n // 2])


def verify_cutoff():
    print("B. 결맞는 차단 (아베)")
    lam, NA = 0.5, 0.5                        # λ/NA = 1 µm
    m = [grating_modulation(p, NA, lam) for p in (1.05, 0.95)]
    ok = check("정면 조명: 주기 1.05 λ/NA 는 보이고 0.95 λ/NA 는 균일", [m[0] > 0.3, m[1] < 1e-9], [1, 1], 0.5)
    fs = -0.98 * NA / lam                       # 동공 가장자리 가까이에서 비스듬히
    m = [grating_modulation(p, NA, lam, illum=(fs, 0.0)) for p in (0.52, 0.48)]
    ok &= check("비스듬한 조명: 0.52 λ/NA 는 보이고 0.48 λ/NA 는 균일", [m[0] > 0.1, m[1] < 1e-9], [1, 1], 0.5)
    return ok


def verify_otf():
    print("C. OTF ↔ 동공의 자기상관 ↔ Goodman 6-32")
    lam, NA = 0.5, 0.5
    n, dx = 1024, 0.05
    O = np.real(otf_from_psf(n, dx, NA, lam))
    _, fgrid = grid(n, dx)
    nu = np.linspace(0, 2, 21)                  # f / (NA/λ)
    f = nu * NA / lam
    route1 = np.interp(f, fgrid[:n // 2], O[0, :n // 2])
    route2 = np.real(otf_by_overlap(f, NA, lam, m=1201))
    s = nu / 2
    closed = 2 / np.pi * (np.arccos(np.clip(s, 0, 1)) - s * np.sqrt(np.clip(1 - s**2, 0, 1)))
    ok = check("|PSF|² 의 FT ↔ 닫힌 식 (차단 2NA/λ 까지)", route1, closed, 0.01)
    ok &= check("동공 겹침 넓이 ↔ 닫힌 식", route2, closed, 0.005)
    if os.path.exists(OPTICORE_PY):
        code = f"""
import numpy as np
from opticore.wavefront import diffraction_limited_mtf
fno = 1 / (2 * {NA})
print(repr([float(v) for v in diffraction_limited_mtf(np.array({list(f * 1000)}), fno, {lam})]))
"""
        try:
            r = subprocess.run([OPTICORE_PY, "-c", code], capture_output=True, text=True, timeout=120,
                               env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            vals = eval(r.stdout.strip().splitlines()[-1])
            ok &= check("동공 겹침 넓이 ↔ opticore diffraction_limited_mtf (F수 1/(2NA))", route2, vals, 0.005)
        except Exception as e:  # noqa: BLE001
            print(f"  [SKIP] opticore 대조 실패: {e}")
    else:
        print("  [SKIP] opticore 환경 없음")
    return ok


def verify_edge():
    print("D. 가장자리")
    lam, NA = 0.5, 0.5
    n, dx = 1024, 0.02
    x, _ = grid(n, dx)
    xe = x + dx / 2                              # 경계를 표본 사이 (x = 0 과 dx 사이의 dx/2) 에 둔다
    step = ((xe > 0) & (xe < n * dx / 2 - 3)).astype(float)
    U = np.tile(step, (n, 1)).astype(complex)
    Ic = np.abs(coherent_image(U, dx, NA, lam)) ** 2
    Ii = incoherent_image(np.abs(U) ** 2, dx, NA, lam)
    row = n // 2
    at_edge = lambda I: 0.5 * (I[row, n // 2 - 1] + I[row, n // 2])      # noqa: E731  경계 양옆 두 표본의 평균
    over = (0.5 + sici(np.pi)[0] / np.pi) ** 2
    sel = (x > -3) & (x < 6)
    ok = check("결맞음: 경계 세기 1/4", at_edge(Ic), 0.25, 3e-3)
    ok &= check("결맞음: 첫 봉우리 (1/2 + Si(π)/π)² = 1.187", Ic[row, sel].max(), over, 3e-3)
    ok &= check("결맞지 않음: 경계 세기 1/2", at_edge(Ii), 0.5, 3e-3)
    return ok


def two_point_image(sep, phase, lam, NA, m=601):
    """두 점을 sep 만큼 띄운 결맞는 상의 가로 단면. 동공 위 표본으로 직접 합한다:
    U(x) = Σ_f P(f) [e^{−iπ f_x sep} + e^{iφ} e^{iπ f_x sep}] e^{i2π f_x x}."""
    c = NA / lam
    g = np.linspace(-c, c, m)
    FX, FY = np.meshgrid(g, g)
    P = (np.hypot(FX, FY) <= c).astype(float)
    fx = FX[P > 0]
    spec = np.exp(-1j * np.pi * fx * sep) + np.exp(1j * phase) * np.exp(1j * np.pi * fx * sep)
    x = np.linspace(-1.2, 1.2, 241) * lam / NA
    U = np.array([np.sum(spec * np.exp(2j * np.pi * fx * xx)) for xx in x])
    return x, np.abs(U) ** 2


def verify_two_point():
    print("E. 두 점의 상")
    lam, NA = 0.5, 0.5
    s = 0.61 * lam / NA
    dips = []
    for ph in (0.0, np.pi / 2, np.pi):
        x, I = two_point_image(s, ph, lam, NA)
        dips.append(I[np.argmin(np.abs(x))] / I.max())
    ok = check("레일리 간격 0.61 λ/NA, 위상차 0, π/2, π: 골/봉우리 ↔ 1, 0.735, 0", dips, [1.0, 0.735, 0.0], 3e-3)
    lo, hi = 0.7 * lam / NA, 1.0 * lam / NA
    for _ in range(30):
        mid = (lo + hi) / 2
        x, I = two_point_image(mid, 0.0, lam, NA)
        d = I[np.argmin(np.abs(x))] / I.max()
        lo, hi = (mid, hi) if d > 0.735 else (lo, mid)
    ok &= check("같은 위상에서 골 0.735 가 되는 간격 ↔ 0.82 λ/NA", (lo + hi) / 2 * NA / lam, 0.82, 0.005)
    return ok


def verify_partial():
    print("F. 부분 결맞음")
    lam, NA = 0.5, 0.5
    m = [grating_modulation(0.8, NA, lam, sigma=sg) for sg in (0.2, 0.3)]
    ok = check("주기 0.8 λ/NA (1/p = 1.25 NA/λ): σ 0.2 는 균일, 0.3 은 보임", [m[0] < 1e-9, m[1] > 0.01], [1, 1], 0.5)
    n, dx = 256, 0.05
    x, _ = grid(n, dx)
    U = np.tile((np.abs(x) < 2).astype(float), (n, 1)).astype(complex)
    ok &= check("σ = 0 ↔ 결맞는 결상", partially_coherent_image(U, dx, NA, lam, 0.0), np.abs(coherent_image(U, dx, NA, lam)) ** 2, 1e-12)
    return ok


def verify_defocus():
    print("G. 초점 어긋남 (정사각 동공, Goodman 6-41)")
    lam, NA = 0.5, 0.5
    nu = np.linspace(0, 1.95, 40)                 # f / (2 f0)… 여기서는 f / f0 (정사각형 반폭 f0 = NA/λ)
    f = nu * NA / lam
    ok = True
    for Wm in (0.25, 0.5, 1.0):
        W = lambda rho, th, Wm=Wm: Wm * rho**2           # noqa: E731  가장자리 (모서리 아닌 변) 에서 Wm
        got = np.real(otf_by_overlap(f, NA, lam, W=W, m=1201, shape="square"))
        s = nu / 2
        closed = (1 - s) * np.sinc(8 * Wm * s * (1 - s))
        ok &= check(f"W_m = {Wm} λ: 동공 겹침 ↔ Λ·sinc[(8W/λ)(f/2f0)(1 − f/2f0)]", got, closed, 0.01)
    return ok


def verify_phase_contrast():
    print("H. 위상차")
    lam, NA = 0.5, 0.5
    n, dx = 256, 0.05
    x, f = grid(n, dx)
    phi0 = 0.02
    p = n * dx / 8
    phi = np.tile(phi0 * np.cos(2 * np.pi * x / p), (n, 1))       # 평균 0 인 약한 위상
    U = np.exp(1j * phi)
    dot = 0.5 / (n * dx)                                            # 0 주파수 한 칸만
    I_pos = np.abs(phase_contrast(U, dx, NA, lam, dot, alpha=-np.pi / 2)) ** 2
    I_neg = np.abs(phase_contrast(U, dx, NA, lam, dot, alpha=+np.pi / 2)) ** 2
    ok = check("0차를 λ/4 앞당김 ↔ 1 − 2φ (오차는 φ² 차수)", I_pos, 1 - 2 * phi, 3 * phi0**2)
    ok &= check("0차를 λ/4 늦춤 ↔ 1 + 2φ", I_neg, 1 + 2 * phi, 3 * phi0**2)
    I_abs = np.abs(phase_contrast(U, dx, NA, lam, dot, a=1 / 3, alpha=-np.pi / 2)) ** 2
    gain = (modulation(I_abs[0]) / modulation(I_pos[0]))
    ok &= check("0차 세기를 1/9 (진폭 1/3) 로 줄이면 대비 3배", gain, 3.0, 0.02)
    I_bf = np.abs(coherent_image(U, dx, NA, lam)) ** 2
    I_df = np.abs(dark_field(U, dx, NA, lam, dot)) ** 2
    ok &= check("명시야: 세기 변화 ≤ φ² 차수 (1 − 2φ 의 1/10 미만)", modulation(I_bf[0]) < 0.1 * modulation(I_pos[0]), 1, 0.5)
    ok &= check("암시야: 세기 ≈ φ² (최대 φ0²)", I_df[0].max(), phi0**2, 0.1 * phi0**2)
    return ok


def verify_abbe_porter():
    print("I. 아베-포터: 0차 막기")
    n = 256
    y, x = np.mgrid[:n, :n]
    wire = ((x % 16) < 2) | ((y % 16) < 2)                     # 철사 폭 2, 간격 16
    t = (~wire).astype(complex)
    m = t.real.mean()
    mask = np.ones((n, n))
    mask[0, 0] = 0.0
    I = np.abs(filter_4f(t, mask)) ** 2
    return check(f"열린 비율 m = {m:.3f}: 철사 자리 m², 구멍 자리 (1 − m)²", [I[0, 0], I[8, 8]], [m**2, (1 - m) ** 2], 1e-12)


def main():
    results = [f() for f in (verify_lens_ft, verify_cutoff, verify_otf, verify_edge, verify_two_point, verify_partial,
                             verify_defocus, verify_phase_contrast, verify_abbe_porter)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

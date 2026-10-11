"""coherence.py 를 독립된 경로들과 교차검증한다.

실행: python verify_coherence.py

  A. 위너-힌친 정리           ↔  무작위 파동 표본에서 잰 자기상관 = 스펙트럼의 푸리에 변환 (Saleh & Teich 식 12.1-17)
  B. 결맞음 시간과 선폭       ↔  Saleh & Teich 표 12.1-1 (사각형 1, 로런츠 1/π, 가우스 √(2ln2/π))
  C. 나트륨 이중선            ↔  무늬가 처음 사라지는 광로차 λ²/(2Δλ), 되돌아오는 주기 λ²/Δλ
  D. 푸리에 변환 분광         ↔  긴 주사로 두 선의 파장을 되찾고, 1/Δmax 가 선 간격보다 크면 갈라지지 않는다
  E. 박막 무늬의 가시도       ↔  파장별 평균으로 잰 무늬 포락선 = 두 빔 가시도 × |g(2nt)|
  F. OCT 깊이 분해능          ↔  가우스 스펙트럼에서 (2 ln2/π) λ²/Δλ
  G. 공간 결맞음              ↔  반 시테르트-제르니케: 2J1(x)/x, 첫 영점 1.22 λ/θ (Hecht 식 12.28–12.29)
  H. 베텔게우스               ↔  Hecht 12.5.1: 121 인치, 570 nm → 0.047 초각
  I. 스페클 통계              ↔  한 장의 세기는 지수분포(평균 = 표준편차), N 장 평균의 대비 1/√N
"""
import numpy as np
from scipy.special import j1

from coherence import (autocorrelation, contrast, degree_from_spectrum, gaussian_line, lorentz_line,
                       michelson_interferogram, observed_film_reflectance, oct_ascan, random_wave,
                       spectrum_from_interferogram, speckle_pattern, visibility_extended_source, wavelength_band)

NA_D1, NA_D2 = 0.588995, 0.589592


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<64s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def verify_wiener_khinchin():
    print("A. 위너-힌친: 표본의 자기상관 ↔ 스펙트럼의 푸리에 변환")
    rng = np.random.default_rng(3)
    nu = np.linspace(0, 2, 4001)
    S = gaussian_line(nu, 1.0, 0.1)
    n_t, dt, lags = 2**17, 0.05, 400
    g_meas = np.zeros(lags, complex)
    for _ in range(8):
        g_meas += autocorrelation(random_wave(S, nu, n_t, dt, rng), lags)
    g_meas /= 8
    g_th = degree_from_spectrum(nu, S, np.arange(lags) * dt)
    return check("|g(τ)| 8 표본 평균 ↔ 스펙트럼의 푸리에 변환 (통계 오차)", np.abs(g_meas), np.abs(g_th), 0.03)


def tau_c(sigma, S, dmax, n=8001):
    """Saleh 식 (12.1-10): τc = ∫ |g|² dτ (양쪽). 광로차 단위로 쓴다."""
    D = np.linspace(-dmax, dmax, n)
    g = degree_from_spectrum(sigma, S, D)
    return np.trapezoid(np.abs(g) ** 2, D)


def verify_coherence_time():
    print("B. 결맞음 시간 × 선폭 (Saleh & Teich 표 12.1-1)")
    s0, w = 2.0, 0.01
    sig = np.linspace(s0 - 60 * w, s0 + 60 * w, 24001)
    ok = True
    sig_r = np.linspace(s0 - w / 2, s0 + w / 2, 20001)   # 사각형은 받침 안만 촘촘히 (긴 지연에서 겹침이 없게)
    rect = np.ones_like(sig_r)
    # 사각형 스펙트럼의 |g| 는 sinc 라 꼬리가 1/τ² 로만 줄어든다: 적분 구간을 넓게 잡는다
    ok &= check("사각형: τc·Δν = 1", tau_c(sig_r, rect, 200 / w, n=16001) * w, 1.0, 5e-3)
    ok &= check("로런츠: τc·Δν = 1/π", tau_c(sig, lorentz_line(sig, s0, w), 60 / w) * w, 1 / np.pi, 5e-3)
    ok &= check("가우스: τc·Δν = √(2ln2/π)", tau_c(sig, gaussian_line(sig, s0, w), 20 / w) * w, np.sqrt(2 * np.log(2) / np.pi), 5e-3)
    return ok


def verify_sodium():
    print("C. 나트륨 D 이중선의 맥놀이")
    sig = np.linspace(1 / 0.5900, 1 / 0.5884, 30001)
    S = gaussian_line(sig, 1 / NA_D1, 2e-6) + gaussian_line(sig, 1 / NA_D2, 2e-6)
    D = np.linspace(0, 1200, 12001)
    g = np.abs(degree_from_spectrum(sig, S, D))
    lam = 2 / (1 / NA_D1 + 1 / NA_D2)
    beat = 1 / abs(1 / NA_D1 - 1 / NA_D2)          # 정확한 값: 파수 차의 역수
    ok = check("첫 영점 광로차 ↔ 1/(2Δσ) ≈ λ²/(2Δλ) (µm)", D[np.argmin(g[:6000])], beat / 2, 0.1)
    i2 = 4000 + np.argmax(g[4000:9000])
    ok &= check("되돌아오는 주기 ↔ λ²/Δλ (µm, 상대)", D[i2] / (lam**2 / (NA_D2 - NA_D1)) - 1, 0.0, 2e-3)
    print(f"      광로차 주기 {beat:.1f} µm = 거울 이동 {beat / 2:.1f} µm")
    return ok


def verify_fts():
    print("D. 푸리에 변환 분광")
    sig = np.linspace(1 / 0.5900, 1 / 0.5884, 30001)
    S = gaussian_line(sig, 1 / NA_D1, 2e-6) + gaussian_line(sig, 1 / NA_D2, 2e-6)
    out = np.linspace(1 / 0.5900, 1 / 0.5884, 4001)
    ok = True
    for Dmax, resolved in [(1800.0, True), (300.0, False)]:
        D = np.linspace(0, Dmax, int(Dmax * 100) + 1)
        I = michelson_interferogram(sig, S, D)
        Ds = np.concatenate([-D[::-1], D[1:]])
        Is = np.concatenate([I[::-1], I[1:]])
        spec = spectrum_from_interferogram(Ds, Is, out)
        peaks = [i for i in range(1, len(spec) - 1) if spec[i] > spec[i - 1] and spec[i] > spec[i + 1] and spec[i] > 0.5 * spec.max()]
        lams = np.sort(1 / out[peaks])
        if resolved:
            ok &= check(f"Δmax = {Dmax / 1000} mm: 두 봉우리 파장 (nm)", lams * 1e3, [NA_D1 * 1e3, NA_D2 * 1e3], 0.02)
        else:
            ok &= check(f"Δmax = {Dmax / 1000} mm: 1/Δmax > Δσ 이면 봉우리 하나", len(lams), 1, 0.5)
    print(f"      선 간격 Δσ = {abs(1 / NA_D1 - 1 / NA_D2) * 1e4:.2f} cm⁻¹, 1/Δmax = {1e4 / 1800:.2f}, {1e4 / 300:.2f} cm⁻¹")
    return ok


def verify_film():
    print("E. 박막 무늬의 가시도 ↔ 두 빔 가시도 × |g(2nt)|")
    n = 1.33
    sig = np.linspace(1 / 0.75, 1 / 0.40, 3001)
    s0, ds = wavelength_band(0.55, 0.15)
    S = gaussian_line(sig, s0, ds)
    r = (1 - n) / (1 + n)
    I1, I2 = r**2, ((1 - r**2) * r) ** 2
    V0 = 2 * np.sqrt(I1 * I2) / (I1 + I2)
    ok = True
    # 포락선이 무늬 한 주기 안에서 거의 변하지 않는 두께에서 잰다 (얇을수록 포락선이 가파르다)
    for tc in [0.6, 1.0, 1.5]:
        tt = np.linspace(tc - 0.25, tc + 0.25, 2001)
        R = observed_film_reflectance(n, tt, sig, S)
        win = np.abs(tt - tc) < 0.5 / (2 * n * s0) * 0.6        # 무늬 한 주기 남짓
        Vm = (R[win].max() - R[win].min()) / (R[win].max() + R[win].min())
        g = np.abs(degree_from_spectrum(sig, S, [2 * n * tc]))[0]
        ok &= check(f"두께 {tc} µm: 잰 가시도 ↔ V0·|g(2nt)|", Vm, V0 * g, 0.03)
    return ok


def verify_oct():
    print("F. OCT 깊이 분해능 ↔ (2 ln2/π) λ²/Δλ")
    from scipy.signal import hilbert
    ok = True
    for lam0, dlam in [(0.84, 0.05), (1.31, 0.1)]:
        s0, ds = wavelength_band(lam0, dlam)
        sig = np.linspace(s0 - 6 * ds, s0 + 6 * ds, 3001)
        S = gaussian_line(sig, s0, ds)
        z = np.linspace(60, 140, 16001)
        I = oct_ascan(sig, S, z, [(100.0, 0.05)])
        env = np.abs(hilbert(I - np.median(I)))
        above = z[env > env.max() / 2]
        pred = 2 * np.log(2) / np.pi / ds            # λ²/Δλ 의 정확한 꼴은 파수 반치폭의 역수 1/Δσ
        ok &= check(f"{lam0 * 1e3:.0f} nm, Δλ {dlam * 1e3:.0f} nm: 포락선 반치폭 (µm, 상대)", (above.max() - above.min()) / pred - 1, 0.0, 0.01)
        print(f"      반치폭 {above.max() - above.min():.2f} µm, λ²/Δλ 꼴 {2 * np.log(2) / np.pi * lam0**2 / dlam:.2f} µm")
    return ok


def verify_spatial():
    print("G. 반 시테르트-제르니케 ↔ 2J1(x)/x")
    lam, th = 0.55, np.radians(0.53)
    d = np.array([5.0, 20.0, 40.0, 60.0, 100.0, 150.0])
    V = np.array([visibility_extended_source(dd, lam, th) for dd in d])
    x = np.pi * d * th / lam
    ok = check("가시도 ↔ |2J1(x)/x|", V, np.abs(2 * j1(x) / x), 2e-3)
    ds_ = np.linspace(60, 85, 2501)
    Vz = np.array([visibility_extended_source(dd, lam, th, n_pts=2001) for dd in ds_[::25]])
    ok &= check("첫 영점 ↔ 1.22 λ/θ (µm)", ds_[::25][np.argmin(Vz)], 1.22 * lam / th, 0.3)
    print(f"      햇빛(0.53°, 550 nm) 결맞음 거리 1.22λ/θ = {1.22 * lam / th:.1f} µm")
    return ok


def verify_betelgeuse():
    print("H. 베텔게우스 (Hecht 12.5.1)")
    h = 121 * 0.0254
    theta = 1.22 * 570e-9 / h
    return check("각 지름 (초각)", np.degrees(theta) * 3600, 0.047, 5e-4)


def verify_speckle():
    print("I. 스페클 통계")
    rng = np.random.default_rng(5)
    I = speckle_pattern(512, 0.12, rng)
    ok = check("한 장: 대비 = 1 (지수분포)", contrast(I), 1.0, 0.05)
    ok &= check("한 장: P(I > 2<I>) = e⁻²", np.mean(I > 2), np.exp(-2), 0.01)
    for N in [4, 16]:
        Iavg = sum(speckle_pattern(512, 0.12, rng) for _ in range(N)) / N
        ok &= check(f"{N} 장 평균: 대비 = 1/√N", contrast(Iavg), 1 / np.sqrt(N), 0.03)
    return ok


def main():
    results = [f() for f in (verify_wiener_khinchin, verify_coherence_time, verify_sodium, verify_fts, verify_film,
                             verify_oct, verify_spatial, verify_betelgeuse, verify_speckle)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

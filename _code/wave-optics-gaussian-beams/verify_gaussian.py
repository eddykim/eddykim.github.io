"""gaussian.py 의 두 경로(q 파라미터 해석 / 격자 위 파동 전파)와 교재 수치를 대조한다.

실행: python verify_gaussian.py

  A. 자유 전파                ↔  각스펙트럼으로 전파한 장의 폭·곡률·축 위 위상 = W0√(1+(z/z0)²), z(1+(z0/z)²), tan⁻¹(z/z0)
  B. 원 안의 에너지           ↔  반지름 W 안에 1 − e⁻² = 86.5 % (Saleh & Teich 식 3.1-17 의 86 %)
  C. 렌즈                     ↔  렌즈 위상을 곱해 전파한 장의 허리 위치·크기 = ABCD 법칙 q2 = (Aq1+B)/(Cq1+D),
                                  평행 빔의 초점 2W0' ≈ (4/π) λ F# (Saleh & Teich 식 3.2-17)
  D. 공진기 안정성            ↔  광선을 왕복시켜 유계로 남는 조건 = 0 ≤ g1 g2 ≤ 1 (g = 1 − L/R)
  E. 공진기 모드              ↔  폭스-리 반복으로 수렴한 장의 폭 = 왕복 ABCD 의 자기 일치 q. 공초점 He-Ne(L = 28 cm):
                                  w0 = 0.168 mm, 전각 발산 2.40 mrad (Hecht 예제 13.6)
  F. 고차 모드의 구이 위상    ↔  HG_lm 은 (l+m+1)ζ, LG_p^l 은 (2p+|l|+1)ζ (Saleh & Teich 식 3.3-10, 3.4-1)
  G. 빔 품질                  ↔  HG00 과 HG10 을 같은 세기로 결맞지 않게 섞은 빔: 2차 모멘트 폭이 쌍곡선으로 자라고 M² = 2
"""
import numpy as np

from gaussian import (M_free, M_lens, W_R_from_q, abcd, chain, encircled_power, fox_li, gaussian_field, grid,
                      hermite_gaussian, is_stable_rays, laguerre_gaussian, q_at_waist, rayleigh_range,
                      resonator_mode_q, second_moment_width, symmetric_mode_waist, angular_spectrum, thin_lens_phase, waist_from_q)


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<66s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def fit_curvature(U, X, Y, W):
    """장의 위상을 r² 에 대해 직선으로 맞춰 곡률 반지름을 읽는다 (exp(ik r²/2R))."""
    r2 = (X**2 + Y**2).ravel()
    ph = np.unwrap(np.angle(U[U.shape[0] // 2]))
    x = X[U.shape[0] // 2]
    sel = np.abs(x) < 0.8 * W
    c = np.polyfit(x[sel] ** 2, ph[sel], 1)
    return np.inf if abs(c[0]) < 1e-15 else np.pi / (lam_g * c[0])


lam_g = 0.5


def verify_free():
    print("A. 자유 전파: 각스펙트럼 ↔ 가우시안 빔 공식")
    lam, W0 = lam_g, 5.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 512, 0.5
    x, (X, Y) = grid(n, dx)
    U0 = gaussian_field(X, Y, W0, lam)
    zs = np.array([0.5, 1.0, 2.0, 3.0]) * z0
    W, R, G = [], [], []
    for z in zs:
        U = angular_spectrum(U0, dx, lam, z)
        Wz = second_moment_width(np.abs(U) ** 2, x)
        W.append(Wz)
        R.append(fit_curvature(U, X, Y, Wz))
        G.append(-np.angle(U[n // 2, n // 2]))                    # 축 위 위상의 늦음
    ok = check("폭 W(z) ↔ W0 √(1 + (z/z0)²)", np.array(W) / (W0 * np.sqrt(1 + (zs / z0) ** 2)), 1, 2e-3)
    ok &= check("곡률 R(z) ↔ z (1 + (z0/z)²)  (상대)", np.array(R) / (zs * (1 + (z0 / zs) ** 2)), 1, 0.01)
    ok &= check("축 위 위상의 늦음 ↔ 구이 위상 tan⁻¹(z/z0)", G, np.arctan(zs / z0), 2e-3)
    Wq, Rq = W_R_from_q(zs + q_at_waist(W0, lam), lam)
    ok &= check("q 파라미터 경로 ↔ 같은 공식 (폭)", Wq / (W0 * np.sqrt(1 + (zs / z0) ** 2)), 1, 1e-12)
    return ok


def verify_power():
    print("B. 반지름 W 안의 에너지")
    lam, W0 = lam_g, 5.0
    x, (X, Y) = grid(512, 0.1)
    I = np.abs(gaussian_field(X, Y, W0, lam)) ** 2
    return check("반지름 W0 안 ↔ 1 − e⁻² = 0.865", encircled_power(I, X, Y, W0), 1 - np.exp(-2), 3e-3)


def numeric_waist(U_after_lens, dx, lam, x, z_lo, z_hi, n_z=81):
    zz = np.linspace(z_lo, z_hi, n_z)
    Ws = np.array([second_moment_width(np.abs(angular_spectrum(U_after_lens, dx, lam, z)) ** 2, x) for z in zz])
    i = np.argmin(Ws)
    lo, hi = max(i - 3, 0), min(i + 4, n_z)
    c = np.polyfit(zz[lo:hi], Ws[lo:hi] ** 2, 2)                  # W² 는 허리 근처에서 z 의 이차식
    zmin = -c[1] / (2 * c[0])
    return zmin, np.sqrt(np.polyval(c, zmin))


def verify_lens():
    print("C. 렌즈: 파동 전파 ↔ ABCD 법칙")
    lam, W0, f = lam_g, 20.0, 2000.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 1024, 0.5
    x, (X, Y) = grid(n, dx)
    ok = True
    for d1 in (0.0, z0, 3 * z0):
        U = angular_spectrum(gaussian_field(X, Y, W0, lam), dx, lam, d1) * thin_lens_phase(X, Y, f, lam)
        q2 = abcd(M_lens(f), d1 + q_at_waist(W0, lam))
        zw, Ww = waist_from_q(q2, lam)                             # 허리는 렌즈 뒤 −zw 에 있다
        zn, Wn = numeric_waist(U, dx, lam, x, 0.6 * (-zw), 1.4 * (-zw) if -zw > 0 else 1.0)
        ok &= check(f"렌즈 앞 {d1 / z0:.0f} z0: 허리 위치·크기 (µm)", [zn, Wn], [-zw, Ww], 2.0)
    # 평행 빔(허리가 렌즈에 있고 z0 ≫ f) 의 초점: 2W0' ≈ (4/π) λ F#, F# = f / (2W0)
    W0b, fb = 2000.0, 50e3
    q2 = abcd(M_lens(fb), q_at_waist(W0b, lam))
    _, Wf = waist_from_q(q2, lam)
    ok &= check("평행 빔의 초점 2W0' ↔ (4/π) λ F#  (상대)", 2 * Wf / (4 / np.pi * lam * fb / (2 * W0b)), 1, 1e-3)
    return ok


def verify_stability():
    print("D. 공진기 안정성: 광선 왕복 ↔ 0 ≤ g1 g2 ≤ 1")
    L = 1.0
    rng = np.random.default_rng(0)
    agree = 0
    total = 0
    for _ in range(400):
        g1, g2 = rng.uniform(-2, 2, 2)
        prod = g1 * g2
        if abs(prod) < 0.02 or abs(prod - 1) < 0.02 or abs(g1 - 1) < 1e-3 or abs(g2 - 1) < 1e-3:
            continue                                               # 경계 근처는 빼고
        R1, R2 = L / (1 - g1), L / (1 - g2)
        total += 1
        agree += is_stable_rays(L, R1, R2) == (0 < prod < 1)
    return check(f"무작위 공진기 {total} 개: 판정 일치 비율", agree / total, 1.0, 1e-12)


def verify_resonator():
    print("E. 공진기 모드: 폭스-리 ↔ 자기 일치 q, Hecht 예제 13.6")
    lam = 0.6328
    L, R = 300e3, 600e3                                            # 30 cm, 거울 곡률 60 cm (g = 0.5)
    q = resonator_mode_q(L, R, R, lam)
    W_mirror, R_mirror = W_R_from_q(q, lam)
    U, losses, _ = fox_li(512, 8.0, lam, L, R, a=600.0, n_iter=1500)
    x, _ = grid(512, 8.0)
    W_num = second_moment_width(np.abs(U) ** 2, x)
    ok = check("폭스-리 1500 왕복 뒤 거울 위 폭 ↔ 자기 일치 q (상대)", W_num / W_mirror, 1, 0.03)
    ok &= check("자기 일치 q 의 파면 곡률 = 거울 곡률 (상대)", abs(R_mirror) / R, 1, 1e-9)
    # Hecht 예제 13.6: 공초점, L = 28 cm, He-Ne
    Lc = 280e3                                                     # 공초점은 왕복 행렬이 −I 라 거울 곡률 맞춤으로 푼다
    w0, z0c = symmetric_mode_waist(Lc, Lc, lam)
    Wm, _ = W_R_from_q(Lc / 2 + q_at_waist(w0, lam), lam)
    ok &= check("공초점 L = 28 cm: 허리 w0 (mm) ↔ 0.168", w0 / 1000, 0.168, 5e-4)
    ok &= check("  전각 발산 2λ/(π w0) (mrad) ↔ 2.40", 2 * lam / (np.pi * w0) * 1000, 2.40, 5e-3)
    ok &= check("  거울 위 폭 = √2 w0", Wm / w0, np.sqrt(2), 1e-9)
    return ok


def overlap_phase(U, T):
    return np.angle(np.sum(U * np.conj(T)))


def verify_higher_modes():
    print("F. 고차 모드의 구이 위상")
    lam, W0 = lam_g, 5.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 512, 0.4
    x, (X, Y) = grid(n, dx)
    z = z0
    Wz = W0 * np.sqrt(2)
    Rz = 2 * z0
    curv = np.exp(1j * np.pi * (X**2 + Y**2) / (lam * Rz))
    zeta = np.pi / 4
    got, want = [], []
    for (l, m) in [(0, 0), (1, 0), (1, 1), (2, 1)]:
        U = angular_spectrum(hermite_gaussian(X, Y, W0, l, m).astype(complex), dx, lam, z)
        T = hermite_gaussian(X, Y, Wz, l, m) * curv               # 같은 모양을 W(z) 로 늘리고 곡률만 붙인 틀 (구이 위상 없음)
        got.append(-overlap_phase(U, T))
        want.append((l + m + 1) * zeta)
    for (p, l) in [(0, 1), (1, 0), (1, 2)]:
        U = angular_spectrum(laguerre_gaussian(X, Y, W0, p, l), dx, lam, z)
        T = laguerre_gaussian(X, Y, Wz, p, l) * curv
        got.append(-overlap_phase(U, T))
        want.append((2 * p + abs(l) + 1) * zeta)
    return check("z = z0 에서 늦음: HG 00,10,11,21 / LG (0,1),(1,0),(1,2) ↔ (N+1)π/4", np.mod(np.array(got) + np.pi, 2 * np.pi) - np.pi,
                 np.mod(np.array(want) + np.pi, 2 * np.pi) - np.pi, 3e-3)


def verify_m2():
    print("G. 빔 품질 M²")
    lam, W0 = lam_g, 5.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 512, 0.5
    x, (X, Y) = grid(n, dx)
    modes = [hermite_gaussian(X, Y, W0, 0, 0), hermite_gaussian(X, Y, W0, 1, 0)]
    modes = [m / np.sqrt(np.sum(np.abs(m) ** 2)) for m in modes]
    zs = np.linspace(-3, 3, 13) * z0
    Wx = []
    for z in zs:
        I = sum(np.abs(angular_spectrum(m.astype(complex), dx, lam, z)) ** 2 for m in modes)   # 결맞지 않게: 세기를 더한다
        Wx.append(second_moment_width(I, x))
    Wx = np.array(Wx)
    c = np.polyfit(zs, Wx**2, 2)                                   # W² = W0m² + θm² z²  (허리 z = 0)
    W0m, thm = np.sqrt(c[2]), np.sqrt(c[0])
    M2 = np.pi * W0m * thm / lam
    ok = check("W(z)² 가 z 의 이차식 (맞춤 잔차 / W0²)", np.max(np.abs(np.polyval(c, zs) - Wx**2)) / W0**2, 0, 1e-3)
    ok &= check("M² = π W0m θm / λ ↔ (1 + 3)/2 = 2", M2, 2.0, 0.01)
    return ok


def main():
    results = [f() for f in (verify_free, verify_power, verify_lens, verify_stability, verify_resonator,
                             verify_higher_modes, verify_m2)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

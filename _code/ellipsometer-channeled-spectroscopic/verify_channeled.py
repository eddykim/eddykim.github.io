"""channeled.py 를 독립된 경로들과 교차검증한다.

실행: python verify_channeled.py

1~4편에서 참고문헌 식을 옮겨 적다 생긴 오독을 이 방식으로 매번 잡았다.
여기서도 한쪽 결과를 다른 쪽에서 유도하지 않는다.

  신호 합성   : 부품 뮬러 행렬의 곱 (channeled.py)
  대조 대상 A : Oka & Kato 1999 식 (4) — 반송 주파수 셋
  대조 대상 B : Hagen 2022 식 (1) — 헤테로다인 일곱
  대조 대상 C : 이승우 2021 식 (3.28) — 지연자 하나의 닫힌 식
  대조 대상 D : Okabe 2009 — 방해석의 투과율 비 각도 44.320 도
  대조 대상 E : Hagen 2022 실측 — 석영 3.89/7.77 mm 가 만드는 OPD 39.3/78.6 um
  대조 대상 F : 이승우 2021 식 (6.4) — 지연량 표류가 만드는 회전
"""
import numpy as np

from channeled import (BETA_QUARTZ, D, N_E_CALCITE, N_O_CALCITE, channel_positions,
                       dsigma_avg, fresnel_gamma, mueller_polarizer, mueller_retarder,
                       mueller_sample, opd, resolving_power, retardance,
                       spectrum_single, spectrum_stokes, thickness_window, to_opd)

TOL = 1e-9
SIGMA = np.linspace(1.0 / 0.80, 1.0 / 0.40, 4096)       # 400~800 nm, 파수 [1/um]
# OPD 축의 눈금 간격. 봉우리 위치는 이보다 정밀하게 잴 수 없다.
DH = 1.0 / (SIGMA[-1] - SIGMA[0])


def peak_centers(h, C, frac=0.02, hmin=None):
    """|C(h)| 의 국소 최대 위치. 봉우리는 폭이 있으므로 문턱만으로는 못 센다."""
    mag = np.abs(C)
    thr = frac * mag.max()
    sel = mag > thr
    if hmin is not None:
        sel &= h > hmin
    idx = np.nonzero(sel)[0]
    if len(idx) == 0:
        return np.array([])
    # 연속한 구간을 하나의 봉우리로 묶고 그 안의 최대점을 대표로 삼는다
    out, start = [], idx[0]
    for a, b in zip(idx[:-1], idx[1:]):
        if b != a + 1:
            seg = np.arange(start, a + 1)
            out.append(h[seg[np.argmax(mag[seg])]])
            start = b
    seg = np.arange(start, idx[-1] + 1)
    out.append(h[seg[np.argmax(mag[seg])]])
    return np.array(out)


def check(name, got, want, tol=TOL):
    err = np.max(np.abs(np.asarray(got) - np.asarray(want)))
    ok = err < tol
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<54s} 최대오차 {err:.3e}")
    return ok


def verify_single_closed_form():
    """지연자 하나: 뮬러 곱이 이승우 2021 식 3.28 과 맞는가.

    I = (Iin/4)[1 - cos2psi cos phi - sin2psi sin(Delta) sin phi]
    """
    print("이승우 식 3.28 — 지연자 하나의 닫힌 식 vs 뮬러 곱")
    ok = True
    for psi_d, dl_d in ((16.8, 87.0), (35.0, 50.0), (27.0, 2.0), (45.0, -120.0)):
        psi, dl = D(psi_d), D(dl_d)
        got = spectrum_single(psi, dl, SIGMA, 2000.0)
        phi = retardance(SIGMA, 2000.0)
        want = 0.25 * (1.0 - np.cos(2 * psi) * np.cos(phi)
                       - np.sin(2 * psi) * np.sin(dl) * np.sin(phi))
        ok &= check(f"psi={psi_d:5.1f}° Delta={dl_d:6.1f}°", got, want)
    return ok


def verify_single_channels():
    """지연자 하나면 OPD 영역에 봉우리가 0, +-L 세 자리에만 서는가."""
    print(f"채널 위치 — 지연자 하나는 봉우리 셋 (OPD 눈금 {DH:.2f} um)")
    ok = True
    for t_um in (1500.0, 3000.0):
        sig = spectrum_single(D(35.0), D(50.0), SIGMA, t_um)
        h, C = to_opd(sig, SIGMA)
        pk = peak_centers(h, C, frac=0.05)
        want = np.array([-opd(t_um), 0.0, opd(t_um)])
        ok &= check(f"t={t_um:.0f} um -> 봉우리 {len(pk)}개, L={opd(t_um):.2f} um",
                    pk, want, tol=DH)
    return ok


def verify_oka_carriers():
    """Oka 식 4 — 지연자 둘이면 반송 주파수가 L2, L1-L2, L1+L2 셋인가.

    Oka 배치는 R1(0) R2(45) A(0) 이고, S1 은 L2 를, S2/S3 는 L1-+L2 를 탄다.
    뮬러 곱으로 만든 신호의 OPD 봉우리 위치를 그 셋과 대조한다.
    """
    print("Oka 식 4 — 반송 주파수 셋 (봉우리 일곱)")
    ok = True
    for t1, t2 in ((3000.0, 1000.0), (2000.0, 4000.0)):
        S = np.tile(np.array([1.0, 0.3, -0.5, 0.6])[:, None], (1, len(SIGMA)))
        sig = spectrum_stokes(S, SIGMA, t1, t2)
        h, C = to_opd(sig, SIGMA)
        got = np.sort(peak_centers(h, C, frac=0.03, hmin=1.0))   # 양의 h 쪽만
        L1, L2 = opd(t1), opd(t2)
        want = np.array(sorted({L2, abs(L1 - L2), L1 + L2}))
        ok &= check(f"t1={t1:.0f} t2={t2:.0f} um  반송파 {len(got)}개", got, want, tol=DH)
    return ok


def verify_hagen_eq1():
    """Hagen 식 1 — 닫힌 식과 뮬러 곱이 같은 스펙트럼을 주는가.

    I = 1/2 s0 + 1/2 s1 cos d2
        + 1/4 [ (s2 sin d1 - s3 cos d1) ... ] 형태를 삼각함수로 풀어 쓴 것이
    아래 want 다. 지수 꼴을 오일러 공식으로 풀면

      1/4 s1 (e^{i d2} + e^{-i d2})               = 1/2 s1 cos d2
      1/8 [(s2-i s3)e^{iA} + (s2+i s3)e^{-iA}]    = 1/4 (s2 cos A + s3 sin A)
      1/8 [(-s2-i s3)e^{iB} + (-s2+i s3)e^{-iB}]  = 1/4 (-s2 cos B + s3 sin B)

    이고 A = d2+d1, B = d2-d1 이다. 직류항은 1/2 s0 다 — 이 배치에는 편광자가
    분석기 하나뿐이라 지연자 하나 배치(편광자 둘, 1/4)와 인자가 다르다.
    """
    print("Hagen 식 1 — 닫힌 식 vs 뮬러 곱")
    ok = True
    for s1, s2, s3 in ((0.3, -0.5, 0.6), (-0.8, 0.1, 0.2), (0.0, 0.0, 1.0)):
        S = np.tile(np.array([1.0, s1, s2, s3])[:, None], (1, len(SIGMA)))
        # Hagen 은 R(45, d2) 라고 적지만 그의 회전 규약은 이 시리즈와 손방향이 반대다.
        # 배치 각도를 훑어 확인한 결과 우리 규약의 az2 = -45 도가 그의 식과 일치한다.
        # 시리즈 규약(1~4편, 김영준 2025 식 2.44)은 그대로 두고 대조할 때만 맞춘다.
        got = spectrum_stokes(S, SIGMA, 3000.0, 1000.0, az2=-np.pi / 4)
        d1 = retardance(SIGMA, 3000.0)
        d2 = retardance(SIGMA, 1000.0)
        want = (0.5 + 0.5 * s1 * np.cos(d2)
                + 0.25 * (s2 * np.cos(d2 + d1) + s3 * np.sin(d2 + d1))
                + 0.25 * (-s2 * np.cos(d2 - d1) + s3 * np.sin(d2 - d1)))
        ok &= check(f"s=({s1:+.1f},{s2:+.1f},{s3:+.1f})", got, want, tol=1e-8)
    return ok


def verify_calcite_gamma():
    """Okabe 2009 — 방해석의 투과율 비 각도가 44.320 도인가.

    논문은 gamma = atan(Ts/Tf) 라고 적지만, 세기 투과율로 그대로 계산하면 43.64 도가
    나온다. 진폭 투과율 비로 보면 44.320 도가 정확히 나오므로 그쪽이 맞다.
    """
    print("Okabe 2009 — 방해석의 투과율 비 각도")
    g = fresnel_gamma(N_E_CALCITE, N_O_CALCITE)     # 음성 1축이라 e 축이 빠르다
    ok = check("gamma (진폭 투과율 비)", np.rad2deg(g), 44.320, tol=1e-3)
    print(f"         참고: 세기 투과율 비로 계산하면 "
          f"{np.rad2deg(np.arctan(np.tan(g) ** 2)):.3f}° 로 논문값과 다르다")
    return ok


def verify_calcite_extra_channel():
    """이색성이 있으면 L1 자리에 없던 채널이 생기고, 2:1 이면 L1-L2 와 겹치는가."""
    print("Okabe 2009 — 이색성이 만드는 여분 채널과 2:1 의 충돌")
    g = fresnel_gamma(N_E_CALCITE, N_O_CALCITE)
    S = np.tile(np.array([1.0, 0.3, -0.5, 0.6])[:, None], (1, len(SIGMA)))
    ok = True
    for tag, (t1, t2) in (("3:1", (3000.0, 1000.0)), ("2:1", (2000.0, 1000.0))):
        apod = np.hanning(len(SIGMA))      # 유한 대역의 사이드로브를 억눌러야 잰다
        ideal = spectrum_stokes(S, SIGMA, t1, t2) * apod
        real = spectrum_stokes(S, SIGMA, t1, t2, gamma1=g, gamma2=g) * apod
        h, Ci = to_opd(ideal, SIGMA)
        _, Cr = to_opd(real, SIGMA)
        L1, L2 = opd(t1), opd(t2)
        sel = np.abs(h - L1) < 0.5
        a_i = np.abs(Ci[sel]).max() / np.abs(Ci).max()
        a_r = np.abs(Cr[sel]).max() / np.abs(Cr).max()
        print(f"         {tag}: L1={L1:.1f} um 자리 세기  이상적 {a_i:.2e} -> 이색성 {a_r:.2e}"
              f"  ({a_r / a_i:.0f}배)   L1-L2={abs(L1-L2):.1f} um")
        if tag == "3:1":
            ok &= a_r / a_i > 100.0          # 없던 채널이 뚜렷하게 생긴다
            ok &= abs(abs(L1 - L2) - L2) > 1.0   # 겹치지 않는다
        else:
            ok &= abs(abs(L1 - L2) - L2) < 1e-9  # 2:1 은 L1-L2 = L2 로 정확히 겹친다
    print(f"  [{'OK  ' if ok else 'FAIL'}] 3:1 은 아홉 자리로 갈라지고 2:1 은 겹친다")
    return ok


def verify_hagen_measured_opd():
    """Hagen 2022 실측 — 석영 3.89 / 7.77 mm 가 39.3 / 78.6 um 의 OPD 를 주는가."""
    print("Hagen 2022 실측 — 두께에서 채널 위치가 나오는가")
    ok = True
    for t_mm, want_um in ((3.89, 39.3), (7.77, 78.6)):
        got = opd(t_mm * 1000.0)
        ok &= check(f"t={t_mm} mm -> OPD (상대오차 {abs(got/want_um-1):.1%})",
                    got / want_um, 1.0, tol=0.015)
    print(f"         논문의 실측 채널 위치에서 역산한 beta 는 "
          f"{39.3/3890:.5f}, {78.6/7770:.5f} 로 본문의 {BETA_QUARTZ} 보다 1.3% 크다")
    print(f"         beta={BETA_QUARTZ} (2편의 공칭 복굴절 0.009 보다 "
          f"{BETA_QUARTZ / 0.009 - 1:.0%} 크다 — 분산 항)")
    return ok


def verify_hagen_budget():
    """Hagen 2022 3절 — 분해 가능한 점 632 개, 채널당 90 개가 나오는가."""
    print("Hagen 2022 식 7 — 분해능 예산")
    # 논문은 식 2 에서 lambda_max = 1.045 um 를 쓰고 식 5 에서는 1035 nm 를 쓴다.
    # 본문의 2.44e-3 은 1035 nm 쪽과 맞으므로 재현에는 그 값을 쓴다.
    ok = check("dsigma_avg (논문 2.44e-3)", dsigma_avg(0.400, 1.035, 0.00101),
               2.44e-3, tol=1e-5)
    rp = resolving_power(0.400, 1.035, 0.00101)
    ok &= check("RP_sigma (논문 632)", rp, 632.0, tol=4.0)
    ok &= check("채널당 점 수 (논문 약 90)", rp / 7.0, 90.0, tol=1.0)
    print(f"         lambda_max 를 식 2 의 1.045 um 로 두면 "
          f"{dsigma_avg(0.400, 1.045, 0.00101):.3e} 로 논문값과 어긋난다.")
    print(f"         논문의 632 는 파수 범위는 1.045 um, 분해능은 1035 nm 로 계산한 값이다"
          f" (1.035 um 로 통일하면 {rp:.1f})")
    return ok


def verify_thickness_window():
    """두께 창 — 아래는 채널 분리, 위는 나이퀴스트. 뮬러면 창이 좁아지는가."""
    print("두께의 허용 창 — 아래는 분리, 위는 나이퀴스트")
    ok = True
    t_lo, t_hi = thickness_window(0.400, 0.800, 0.0005, ratio=(1, 2))
    print(f"         스토크스 CSP (1:2): {t_lo:8.0f} ~ {t_hi:8.0f} um")
    t_lo4, t_hi4 = thickness_window(0.400, 0.800, 0.0005, ratio=(2, 1, 7, 15))
    print(f"         뮬러   CSP (2:1:7:15): {t_lo4:8.0f} ~ {t_hi4:8.0f} um")
    ok &= t_lo < t_hi                       # 스토크스는 창이 열린다
    ok &= t_hi4 < t_lo4                      # 뮬러는 같은 분광기로 창이 비어 있다
    ok &= t_hi4 < t_hi                       # 지연자가 늘면 상한이 내려간다
    ok &= abs(t_lo - t_lo4) < 1e-9           # 하한은 지연자 수와 무관하다
    print(f"  [{'OK  ' if ok else 'FAIL'}] 지연자가 넷이면 상한이 "
          f"{t_hi / t_hi4:.1f} 배 내려가고 하한은 그대로다 "
          f"-> 이 분광기로는 뮬러 CSP 의 창이 비어 있다")
    for dl_nm in (0.5, 0.3, 0.24, 0.2):
        lo4, hi4 = thickness_window(0.400, 0.800, dl_nm / 1000.0, ratio=(2, 1, 7, 15))
        print(f"           dlambda={dl_nm:4.2f} nm -> {lo4:6.0f} ~ {hi4:6.0f} um  "
              f"{'창이 열린다' if lo4 < hi4 else '빈다'}")
    # 이승우 2021 의 유도를 beta 만 바꿔 재현한다
    L_min = 18.0 / (1.0 / 0.400 - 1.0 / 0.800)
    print(f"         이승우 식 3.26 재현: L > {L_min:.1f} um "
          f"-> t > {L_min / BETA_QUARTZ:.0f} um (논문은 dn=0.0045 로 3200 um)")
    return ok


def verify_demodulation_roundtrip():
    """채널을 잘라 Psi, Delta 를 되찾는가 — 반송파 위상까지 되돌려야 한다."""
    print("복조 왕복 — 채널을 잘라 넣은 값을 되찾는가")
    ok = True
    t = 3000.0
    phi = retardance(SIGMA, t)
    L = opd(t)
    mid = np.abs(SIGMA - SIGMA.mean()) < 0.3 * (SIGMA[-1] - SIGMA[0])
    for psi_d, dl_d in ((35.0, 50.0), (16.8, 87.0), (27.0, -120.0)):
        psi, dl = D(psi_d), D(dl_d)
        f = spectrum_single(psi, dl, SIGMA, t)
        h, C = to_opd(f, SIGMA)
        w = np.zeros_like(h)
        hw = L * 0.5
        sel = np.abs(h - L) <= hw
        x = (h[sel] - L) / hw
        w[sel] = 0.5 * (1.0 + np.cos(np.pi * x))
        z = np.fft.ifft(np.fft.ifftshift(C * w)) * len(SIGMA)
        band = 2.0 * z * np.exp(-1j * phi)
        a, b = np.real(band)[mid].mean(), -np.imag(band)[mid].mean()
        got = (-4.0 * a, -4.0 * b)
        want = (np.cos(2 * psi), np.sin(2 * psi) * np.sin(dl))
        ok &= check(f"psi={psi_d:5.1f}° Delta={dl_d:6.1f}° -> (cos2psi, sin2psi sinD)",
                    got, want, tol=2e-4)
    return ok


def verify_drift_rotation():
    """이승우 식 6.4 — 지연량 표류가 (cos2psi, sin2psi sinD) 평면의 회전인가."""
    print("이승우 식 6.4 — 지연량 표류는 회전이다")
    ok = True
    psi, dl = D(35.0), D(50.0)
    for dphi_d in (1.0, 5.0, 20.0):
        dphi = D(dphi_d)
        sig = spectrum_single(psi, dl, SIGMA, 2000.0, dphi=dphi)
        phi = retardance(SIGMA, 2000.0)
        # 표류를 모른 채 이상적 식으로 읽으면 나오는 두 계수
        a = -np.cos(2 * psi) * np.cos(dphi) - np.sin(2 * psi) * np.sin(dl) * np.sin(dphi)
        b = -np.sin(2 * psi) * np.sin(dl) * np.cos(dphi) + np.cos(2 * psi) * np.sin(dphi)
        want = 0.25 * (1.0 + a * np.cos(phi) + b * np.sin(phi))
        ok &= check(f"dphi={dphi_d:5.1f}°  회전한 계수로 재구성", sig, want)
    return ok


if __name__ == "__main__":
    results = [verify_single_closed_form(), verify_single_channels(),
               verify_oka_carriers(), verify_hagen_eq1(),
               verify_calcite_gamma(), verify_calcite_extra_channel(),
               verify_hagen_measured_opd(), verify_hagen_budget(),
               verify_thickness_window(), verify_demodulation_roundtrip(),
               verify_drift_rotation()]
    print("\n전체 결과:", "모두 통과" if all(results) else "실패 항목 있음")
    raise SystemExit(0 if all(results) else 1)

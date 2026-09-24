"""modulation.py 의 뮬러 곱 계산을 Fujiwara 의 해석식과 교차검증한다.

실행: python verify_modulation.py

modulation.py 는 부품 뮬러 행렬을 곱해서만 I(t) 를 만들고, 이 파일은 그 결과를
독립적으로 유도된 닫힌 식과 비교한다. 두 경로가 일치해야 부호 규약이 맞다고
본다. 배경이론 3편에서 TMM 부호 버그를 이 방식으로 잡았다.

대조 대상 (H. Fujiwara, Spectroscopic Ellipsometry, Wiley, 2007):
  식 4.18  RAE          I = I0[1 + S1 cos2A + S2 sin2A]
  식 4.28  RAE+보상자   I = I0[1 + S1 cos2A + (S2 cos d - S3 sin d) sin2A]
  식 4.32  RCE          I = I0[(2 + S1) - 2 S3 sin2C + S1 cos4C + S2 sin4C]  (A=0, d=90°)
  식 4.41  PME          I = I0[1 - S3 sin d + (-S1 sin2M + S2 cos2M) cos d]
  식 4.44  PME          sin d, cos d 의 야코비-안거 전개 (베셀 함수)

무편광에 이상적 편광자 두 장을 통과시키므로 비례상수는 I0 = 1/4 이다.
RCE 만 식 4.32 가 대괄호 안에 2 를 품은 꼴이라 I0 가 1/8 로 적힌다.
"""
import numpy as np
from scipy.special import jv

from modulation import (fourier_coefficients, intensity_pme, intensity_rae,
                        intensity_rae_compensator, intensity_rce,
                        normalized_coefficients, pem_retardation)

TOL = 1e-10
I0 = 0.25   # 무편광 + 이상적 편광자 두 장
# 특수각을 피한 임의의 시편 여러 개로 검사한다.
CASES = [(np.deg2rad(p), np.deg2rad(d))
         for p, d in [(22.0, 143.0), (45.0, -37.0), (63.5, 91.0), (12.0, 5.0)]]


def stokes_after_sample(psi, delta):
    """반사광의 정규화 스토크스 파라미터 (Fujiwara 식 4.18, 4.41 의 S1~S3)."""
    return -np.cos(2 * psi), np.sin(2 * psi) * np.cos(delta), -np.sin(2 * psi) * np.sin(delta)


def check(name, got, want, tol=TOL):
    err = np.max(np.abs(np.asarray(got) - np.asarray(want)))
    status = "OK  " if err < tol else "FAIL"
    print(f"  [{status}] {name:<46s} 최대오차 {err:.3e}")
    return err < tol


def verify_rae():
    """RAE: 뮬러 곱 파형 = 식 4.18, 그리고 푸리에 계수 = (S1, S2)."""
    print("RAE (PSA_R) — 식 4.18 / 4.19")
    ok = True
    A = np.linspace(0, np.pi, 720, endpoint=False)
    for psi, delta in CASES:
        S1, S2, _ = stokes_after_sample(psi, delta)
        got = intensity_rae(psi, delta, A)
        want = I0 * (1 + S1 * np.cos(2 * A) + S2 * np.sin(2 * A))
        ok &= check(f"파형  psi={np.rad2deg(psi):5.1f}° D={np.rad2deg(delta):6.1f}°", got, want)
        # 분석기 180° 회전이 1 optical cycle 이므로 A=wt 를 그대로 한 주기로 본다.
        c = normalized_coefficients(got, n_max=4)
        ok &= check("  푸리에 계수 (alpha, beta) = (S1, S2)", c[1], (S1, S2))
        ok &= check("  4w 성분 없음", c[2], (0.0, 0.0))
    return ok


def verify_rae_compensator():
    """RAE+보상자: 보상자가 Delta 를 이동만 시킨다(식 4.27~4.28)."""
    print("RAE+보상자 (PSCA_R) — 식 4.28")
    ok = True
    A = np.linspace(0, np.pi, 720, endpoint=False)
    for psi, delta in CASES:
        for d_deg in (0.0, 45.0, 90.0, 132.0):
            d = np.deg2rad(d_deg)
            S1, S2, S3 = stokes_after_sample(psi, delta)
            got = intensity_rae_compensator(psi, delta, A, d)
            want = I0 * (1 + S1 * np.cos(2 * A)
                         + (S2 * np.cos(d) - S3 * np.sin(d)) * np.sin(2 * A))
            ok &= check(f"파형  D={np.rad2deg(delta):6.1f}° d={d_deg:5.1f}°", got, want)
        # Delta' = Delta - d 로 "이동만" 하는지: 지연량 d 인 측정은 Delta-d 인 시편의 RAE 와 같다.
        d = np.deg2rad(73.0)
        ok &= check("  Delta' = Delta - d (RAE 와 동치)",
                    intensity_rae_compensator(psi, delta, A, d),
                    intensity_rae(psi, delta - d, A))
    return ok


def verify_rce():
    """RCE: A=0°, d=90° 에서 식 4.32 의 2w·4w 구조가 나오는가."""
    print("RCE (PSC_R A) — 식 4.32")
    ok = True
    C = np.linspace(0, np.pi, 720, endpoint=False)
    for psi, delta in CASES:
        S1, S2, S3 = stokes_after_sample(psi, delta)
        got = intensity_rce(psi, delta, C)
        want = 0.5 * I0 * ((2 + S1) - 2 * S3 * np.sin(2 * C)
                           + S1 * np.cos(4 * C) + S2 * np.sin(4 * C))
        ok &= check(f"파형  psi={np.rad2deg(psi):5.1f}° D={np.rad2deg(delta):6.1f}°",
                    got, want)
        # dc = (2 + S1) 이므로 정규화 계수는 그것으로 나눈 값이 된다.
        c = normalized_coefficients(got, n_max=4)
        ok &= check("  2w 계수: cos 쪽 0, sin 쪽 -2S3/(2+S1)",
                    c[1], (0.0, -2 * S3 / (2 + S1)))
        ok &= check("  4w 계수 = (S1, S2)/(2+S1)", c[2], (S1 / (2 + S1), S2 / (2 + S1)))
    return ok


def verify_pme():
    """PME: 식 4.41, 그리고 지연량을 F sin(wt) 로 흔들 때의 베셀 전개(식 4.44)."""
    print("PME (PSMA) — 식 4.41 / 4.44 / 4.46")
    ok = True
    d = np.deg2rad(np.linspace(-170, 170, 401))
    for psi, delta in CASES:
        S1, S2, S3 = stokes_after_sample(psi, delta)
        for M_deg in (0.0, 45.0, 30.0):
            M = np.deg2rad(M_deg)
            got = intensity_pme(psi, delta, d, mod_angle=M)
            want = I0 * (1 - S3 * np.sin(d)
                         + (-S1 * np.sin(2 * M) + S2 * np.cos(2 * M)) * np.cos(d))
            ok &= check(f"지연량 의존성  D={np.rad2deg(delta):6.1f}° M={M_deg:4.1f}°", got, want)

    # 지연량 진폭 F 를 얼마로 두어야 J0(F)=0 인가 (식 4.45~4.46)
    F0 = 2.404825557695773  # J0 의 첫 영점
    F = np.deg2rad(138.0)   # 실무에서 쓰는 반올림값
    print(f"  J0 의 첫 영점 F = {np.rad2deg(F0):.2f}°,  실무값 138° 에서 J0 = {jv(0, F):+.5f}")
    ok &= check("  J0(137.79°) = 0", jv(0, F0), 0.0, tol=1e-12)
    print(f"  F=138°: 2*J1 = {2*jv(1, F):.3f}   2*J2 = {2*jv(2, F):.3f}"
          f"   (Fujiwara 표기값 1.04, 0.86)")
    ok &= check("  2*J1(138°) = 1.04", 2 * jv(1, F), 1.04, tol=5e-3)
    ok &= check("  2*J2(138°) = 0.86", 2 * jv(2, F), 0.86, tol=5e-3)

    # 야코비-안거 전개는 근사가 아니라 항등식이다. I0 로 정규화하면 정확히 맞아야 한다.
    n_pts = 4096
    t = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)  # omega = 1 로 두고 한 주기
    for psi, delta in CASES:
        got = intensity_pme(psi, delta, pem_retardation(t, F, 1.0))
        dc, c = fourier_coefficients(got, n_max=4)
        s2p = np.sin(2 * psi)
        ok &= check(f"  dc = I0[1 + sin2psi cosD J0(F)]  (D={np.rad2deg(delta):6.1f}°)",
                    dc, I0 * (1 + s2p * np.cos(delta) * jv(0, F)))
        ok &= check("  w  계수 = I0 sin2psi sinD 2J1(F)",
                    c[1][1], I0 * s2p * np.sin(delta) * 2 * jv(1, F))
        ok &= check("  2w 계수 = I0 sin2psi cosD 2J2(F)",
                    c[2][0], I0 * s2p * np.cos(delta) * 2 * jv(2, F))
        ok &= check("  w  의 cos 성분 · 2w 의 sin 성분은 0", (c[1][0], c[2][1]), (0.0, 0.0))
    return ok


def verify_delta_ambiguity():
    """RAE 는 +Delta 와 -Delta 를 구분하지 못하고, RCE 는 구분한다."""
    print("Delta 부호 모호성")
    ok = True
    A = np.linspace(0, np.pi, 720, endpoint=False)
    for psi, delta in CASES:
        if abs(np.sin(delta)) < 1e-6:
            continue
        ok &= check(f"RAE: I(+D) = I(-D)  D={np.rad2deg(delta):6.1f}°",
                    intensity_rae(psi, delta, A), intensity_rae(psi, -delta, A))
        gap = np.max(np.abs(intensity_rce(psi, delta, A) - intensity_rce(psi, -delta, A)))
        status = "OK  " if gap > 1e-3 else "FAIL"
        print(f"  [{status}] RCE: I(+D) != I(-D)                            차이   {gap:.3e}")
        ok &= gap > 1e-3
    return ok


if __name__ == "__main__":
    results = [verify_rae(), verify_rae_compensator(), verify_rce(),
               verify_pme(), verify_delta_ambiguity()]
    print("\n전체 결과:", "모두 통과" if all(results) else "실패 항목 있음")
    raise SystemExit(0 if all(results) else 1)

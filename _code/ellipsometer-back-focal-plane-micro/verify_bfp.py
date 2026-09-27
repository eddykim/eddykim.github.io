"""bfp.py 를 독립된 경로들과 교차검증한다.

실행: python verify_bfp.py

1~5편에서 참고문헌 식을 옮겨 적다 생긴 오독을 이 방식으로 매번 잡았다.

  신호 합성   : 부품 뮬러 행렬의 곱 (bfp.py)
  대조 대상 A : Ye 2007 의 수치 — NA 0.9 -> 64.16 도, 유리 n=1.5 -> 56.3 도
  대조 대상 B : 아베 사인조건 — 반경과 입사각의 왕복
  대조 대상 C : 뮬러 곱 전개 — 고리 신호가 2phi 와 4phi 성분으로 갈리는가
  대조 대상 D : Munro & Torok 2008 — 존스 평균은 편광을 유지하고 뮬러 평균은 소멸시킨다
  대조 대상 E : 김영준 2025 식 4.2 — 브루스터 반경으로 세운 눈금이 참 입사각을 주는가
"""
import numpy as np

from bfp import (D, NA_DEFAULT, N_FUSED_SILICA, N_SI, N_SIO2, angle_to_radius,
                 annulus, averaged_mueller, bfp_intensity, brewster_angle,
                 calibrate_radius, depolarization_index, fourier_alpha,
                 invert_annulus, jones_sample, jones_to_mueller, max_angle,
                 mueller_sample, psi_delta_bare, psi_delta_film, radius_to_angle,
                 tan2_psi_profile)

TOL = 1e-9


def check(name, got, want, tol=TOL):
    err = np.max(np.abs(np.asarray(got) - np.asarray(want)))
    ok = err < tol
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<52s} 최대오차 {err:.3e}")
    return ok


def verify_ye_numbers():
    """Ye 2007 이 본문에 적은 수치를 재현하는가."""
    print("Ye 2007 — 논문이 적은 수치")
    ok = check("NA 0.9 의 최대 입사각 (논문 64.16°)",
               np.rad2deg(max_angle(0.9)), 64.16, tol=5e-3)
    ok &= check("유리 n=1.5 의 브루스터각 (논문 56.3°)",
                np.rad2deg(brewster_angle(1.5)), 56.3, tol=1e-2)
    # 논문의 금 유사 브루스터각 67.9° 는 재현되지 않는다. 정의가 다른 듯하다
    th_au = np.rad2deg(brewster_angle(0.35 - 2.45j))
    print(f"         금 N=0.35-2.45j 의 |rp| 최소각 {th_au:.1f}° "
          f"(논문은 67.9° — 정의가 달라 보여 본문에 쓰지 않는다)")
    return ok


def verify_abbe_roundtrip():
    """아베 사인조건 — 반경과 입사각이 서로의 역인가."""
    print("아베 사인조건 — 반경 <-> 입사각 왕복")
    ok = True
    for na in (0.55, 0.75, 0.9, 0.95):
        r = np.linspace(0.0, 1.0, 51)
        ok &= check(f"NA={na}: r -> theta -> r",
                    angle_to_radius(radius_to_angle(r, na), na), r)
    # 가장자리가 최대 입사각인가
    for na in (0.55, 0.9):
        ok &= check(f"NA={na}: r=1 이 arcsin(NA) 인가",
                    radius_to_angle(1.0, na), max_angle(na))
    return ok


def verify_annulus_harmonics():
    """고리 신호를 전개하면 2phi 와 4phi 성분만 남는가. ★★

    편광자와 분석기가 실험실 좌표계에 고정돼 있으므로 시편 좌표계에서는 둘이 함께
    돈다. 3편이 분류한 "편광자·분석기 동시 회전"이고, 1편의 RAE(2omega 뿐)와 다르다.
    뮬러 곱을 전개하면

      I ∝ (3 + u)/2 - 2 c2 cos2phi + (1 - u)/2 cos4phi,
      c2 = cos2Psi,  u = sin2Psi cosDelta

    이 되고 사인 성분은 항등적으로 0 이다.
    """
    print("고리 신호의 고조파 — 2phi 와 4phi 만 남는가 (사인 성분은 0)")
    ok = True
    phi = np.linspace(0.0, 2 * np.pi, 2880, endpoint=False)
    B = np.array([np.ones_like(phi), np.cos(2 * phi), np.sin(2 * phi),
                  np.cos(4 * phi), np.sin(4 * phi), np.cos(6 * phi)]).T
    for psi_d, dl_d in ((16.8, 87.0), (35.0, 50.0), (27.0, 140.0), (45.0, 90.0)):
        psi, dl = D(psi_d), D(dl_d)
        I = bfp_intensity(psi, dl, phi, pol=0.0, ana=0.0)
        got, *_ = np.linalg.lstsq(B, I, rcond=None)
        c2 = np.cos(2 * psi)
        u = np.sin(2 * psi) * np.cos(dl)
        want = np.array([(3 + u) / 2, -2 * c2, 0.0, (1 - u) / 2, 0.0, 0.0]) / 4
        ok &= check(f"psi={psi_d:5.1f}° Delta={dl_d:6.1f}° -> 여섯 계수",
                    got, want, tol=1e-9)
    return ok


def verify_annulus_roundtrip():
    """고리에서 Psi, Delta 를 되찾는가 — 실제 막 시편으로."""
    print("왕복 — 고리에서 읽은 Psi, Delta 가 참값인가")
    ok = True
    for d_nm in (0.0, 32.9, 53.6):
        for r in (0.4, 0.7, 0.95):
            phi, I, (psi, dl) = annulus(d_nm, r, na=NA_DEFAULT, n_phi=720,
                                        pol=0.0, ana=0.0)
            a2, a4 = fourier_alpha(I, phi)
            psi_r, dl_r = invert_annulus(a2, a4)
            # arccos 이라 Delta 의 부호를 잃는다 (1편의 RAE 와 같은 한계)
            ok &= check(f"막 {d_nm:5.1f} nm, r={r:.2f} "
                        f"(theta={np.rad2deg(radius_to_angle(r)):.1f}°)",
                        (np.rad2deg(psi_r), np.rad2deg(dl_r)),
                        (np.rad2deg(psi), np.rad2deg(abs(dl))), tol=1e-7)
    return ok


def verify_jones_vs_mueller_average():
    """Munro & Torok 2008 — 무엇을 평균하느냐가 편광 소멸을 가른다. ★★

    존스 행렬 평균(공초점)은 편광을 유지하고, 뮬러 행렬 평균(통상)은 소멸시킨다.
    """
    print("Munro & Torok 2008 — 존스 평균 vs 뮬러 평균")
    ok = True
    th_c = D(45.0)
    print(f"{'각도 폭':>10} {'존스 평균 지수':>14} {'뮬러 평균 지수':>14}")
    prev = None
    for hw_deg in (0.5, 2.0, 5.0, 10.0):
        hw = D(hw_deg)
        Mj = averaged_mueller(th_c, hw, 32.9, coherent=True)
        Mm = averaged_mueller(th_c, hw, 32.9, coherent=False)
        dj, dm = depolarization_index(Mj), depolarization_index(Mm)
        print(f"   ±{hw_deg:5.1f}° {dj:14.9f} {dm:14.9f}")
        ok &= abs(dj - 1.0) < 1e-9          # 존스 평균은 순수 상태를 유지한다
        ok &= dm < 1.0                       # 뮬러 평균은 편광을 소멸시킨다
        if prev is not None:
            ok &= dm < prev                  # 폭이 넓어질수록 더 소멸한다
        prev = dm
    print(f"  [{'OK  ' if ok else 'FAIL'}] 존스 평균은 지수 1 을 지키고 "
          f"뮬러 평균은 폭이 넓어질수록 떨어진다")
    return ok


def verify_focusing_blur():
    """집속이 tan(Psi) 를 뭉갠다 — Cohn/Erman 의 합성곱 효과."""
    print("집속 평균 — 각도 폭이 넓어질수록 Psi 가 참값에서 멀어진다")
    th_c = D(45.0)
    psi0, _ = psi_delta_film(np.array([th_c]), 32.9)
    errs = []
    print(f"{'각도 폭':>10} {'Psi 오차(deg)':>16}")
    for hw_deg in (0.5, 2.0, 5.0, 10.0):
        M = averaged_mueller(th_c, D(hw_deg), 32.9, coherent=False)
        psi_eff = 0.5 * np.arccos(np.clip(-M[0, 1], -1.0, 1.0))
        e = abs(np.rad2deg(psi_eff - psi0[0]))
        errs.append(e)
        print(f"   ±{hw_deg:5.1f}° {e:15.5f}")
    ok = all(b > a for a, b in zip(errs, errs[1:]))
    print(f"  [{'OK  ' if ok else 'FAIL'}] 오차가 단조 증가한다 "
          f"({errs[0]:.4f}° -> {errs[-1]:.4f}°)")
    return ok


def verify_brewster_calibration():
    """김영준 2025 식 4.2 — 브루스터 반경으로 세운 눈금이 참 입사각을 주는가. ★"""
    print("브루스터각 보정 — 반경 눈금이 참 입사각을 주는가")
    ok = True
    for na in (0.85, 0.9, 0.95):
        r = np.linspace(0.05, 1.0, 40)
        r_b, th_b, th_cal = calibrate_radius(r, na, N_FUSED_SILICA)
        th_true = radius_to_angle(r, na)
        ok &= check(f"NA={na}: 보정 눈금 vs 참 입사각 "
                    f"(r_b={r_b:.3f}, theta_b={np.rad2deg(th_b):.2f}°)",
                    np.rad2deg(th_cal), np.rad2deg(th_true), tol=2e-2)
    print(f"         용융 실리카 n={N_FUSED_SILICA} -> 브루스터각 "
          f"{np.rad2deg(brewster_angle(N_FUSED_SILICA)):.2f}° (논문 55°)")
    return ok


def verify_brewster_needs_na():
    """브루스터각에 닿으려면 NA 가 얼마여야 하는가 — 재료마다 다르다."""
    print("브루스터각에 닿는 데 필요한 NA")
    ok = True
    rows = [("용융 실리카", N_FUSED_SILICA), ("유리 n=1.5", 1.5),
            ("실리콘", N_SI)]
    for name, n in rows:
        tb = brewster_angle(n)
        na_need = np.sin(tb)
        reach = "닿는다" if na_need <= 0.95 else "닿지 않는다"
        print(f"         {name:12s} 브루스터 {np.rad2deg(tb):5.2f}° "
              f"-> NA {na_need:.3f} 필요  (NA 0.95 로 {reach})")
    ok &= np.sin(brewster_angle(N_FUSED_SILICA)) < 0.95
    ok &= np.sin(brewster_angle(N_SI)) > 0.95      # 실리콘은 못 닿는다
    print(f"  [{'OK  ' if ok else 'FAIL'}] 유전체는 닿고 실리콘은 못 닿는다")
    return ok


def verify_tan2psi_minimum():
    """tan^2(Psi) 의 최소점이 정말 브루스터각인가."""
    print("tan^2(Psi) 의 최소점 = 브루스터각")
    ok = True
    for n in (1.458, 1.5, 2.0):
        rr = np.linspace(1e-4, 1.0, 200001)
        prof = tan2_psi_profile(rr, 0.95, n)
        th_min = radius_to_angle(rr[int(np.argmin(prof))], 0.95)
        ok &= check(f"n={n}: 최소점 각도 vs arctan(n)",
                    np.rad2deg(th_min), np.rad2deg(np.arctan(n)), tol=2e-3)
    return ok


if __name__ == "__main__":
    results = [verify_ye_numbers(), verify_abbe_roundtrip(),
               verify_annulus_harmonics(), verify_annulus_roundtrip(),
               verify_jones_vs_mueller_average(), verify_focusing_blur(),
               verify_tan2psi_minimum(), verify_brewster_calibration(),
               verify_brewster_needs_na()]
    print("\n전체 결과:", "모두 통과" if all(results) else "실패 항목 있음")
    raise SystemExit(0 if all(results) else 1)

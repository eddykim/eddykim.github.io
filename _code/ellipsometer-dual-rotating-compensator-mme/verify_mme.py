"""mme.py 의 뮬러 곱 계산을 Collins & Koh 1999 의 닫힌 식과 교차검증한다.

실행: python verify_mme.py

1편에서 참고문헌 식을 옮겨 적다 생긴 오독을, 2편에서 규약 불일치를 이 방식으로 잡았다.
여기서도 한쪽 결과를 다른 쪽에서 유도하지 않는다.

  경로 A: 부품 뮬러 행렬의 곱 -> I(C) -> FFT
  경로 B: Collins & Koh 표 2 의 닫힌 식

대조 대상 (R. W. Collins and J. Koh, JOSA A 16, 1997 (1999)):
  표 2  5:3 배치의 비영 푸리에 계수. c_j = cos^2(d_j/2), s_j = sin^2(d_j/2)
  소멸  n = 9, 12, 14, 15 의 계수 8개
구조 검사:
  delta1 -> 0 이면 첫 보상자가 사라져 1편의 회전보상자 배치(RCE)로 환원되어야 한다
"""
import numpy as np

from mme import (RATIO, condition_number, fourier_coefficients, highest_harmonic,
                 intensity, mueller_polarizer, mueller_retarder, optical_cycle,
                 reduction_matrix)

TOL = 1e-9
N_PTS = 4096
rng = np.random.default_rng(0)


def check(name, got, want, tol=TOL):
    err = np.max(np.abs(np.asarray(got) - np.asarray(want)))
    ok = err < tol
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<50s} 최대오차 {err:.3e}")
    return ok


def random_sample():
    """16요소가 모두 비영인 시편 행렬. 물리적 조건은 걸지 않는다(선형성만 본다)."""
    M = rng.uniform(-0.4, 0.4, (4, 4))
    M[0, 0] = 1.0
    return M


def coeffs_of(M, d1, d2, P=0.0, A=0.0, n_max=20):
    """I(C) 를 만들고 dc 로 정규화한 푸리에 계수를 돌려준다."""
    I = intensity(M, optical_cycle(N_PTS), d1, d2, RATIO, pol=P, ana=A)
    dc, co = fourier_coefficients(I, n_max)
    return dc, {n: (a / dc, b / dc) for n, (a, b) in co.items()}


def verify_vanishing():
    """5:3 에서 n = 9, 12, 14, 15 의 계수가 소멸하는가."""
    print("소멸 고조파 — Collins & Koh 가 지목한 n = 9, 12, 14, 15")
    ok = True
    for d_deg in (90.0, 128.0, 60.0):
        d = np.deg2rad(d_deg)
        _, co = coeffs_of(random_sample(), d, d)
        dead = [np.hypot(*co[n]) for n in (9, 12, 14, 15)]
        alive = [np.hypot(*co[n]) for n in (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 16)]
        ok &= check(f"delta={d_deg:5.1f}°  n=9,12,14,15 가 0", dead, [0] * 4)
        good = min(alive) > 1e-4
        print(f"  [{'OK  ' if good else 'FAIL'}] {'  나머지 12개는 살아 있다':<48s} "
              f"최소크기 {min(alive):.3e}")
        ok &= good
    # 32C 를 넘는 고조파는 없어야 한다
    ok &= check("최고차 고조파가 n=16 (32C)", highest_harmonic(RATIO), 16, tol=0.5)
    return ok


def verify_table2():
    """표 2 의 닫힌 식 몇 개를 뮬러 곱 결과와 직접 대조한다.

    c_j = cos^2(d_j/2), s_j = sin^2(d_j/2) 이고 P', A' 는 참각이다.
    여기서는 교정량을 0 으로 두었으므로 P' = P, A' = A 다.
    """
    print("Collins & Koh 표 2 — 닫힌 식과의 대조")
    ok = True
    for d1_deg, d2_deg, P_deg, A_deg in [(90, 90, 0, 0), (128, 128, 12, -7),
                                         (110, 75, 25, 40)]:
        d1, d2 = np.deg2rad(d1_deg), np.deg2rad(d2_deg)
        P, A = np.deg2rad(P_deg), np.deg2rad(A_deg)
        M = random_sample()
        _, co = coeffs_of(M, d1, d2, P, A)
        c1, s1 = np.cos(d1 / 2) ** 2, np.sin(d1 / 2) ** 2
        c2, s2 = np.cos(d2 / 2) ** 2, np.sin(d2 / 2) ** 2
        pm, pp = 2 * (P - A), 2 * (P + A)
        tag = f"d=({d1_deg},{d2_deg})° P={P_deg}° A={A_deg}°"
        # 표 2 의 계수는 dc 로 정규화된 값이다. dc 자리의 식으로 나눠야 맞는다.
        dc_expr = (1
                   + M[0, 1] * c1 * np.cos(2 * P) + M[0, 2] * c1 * np.sin(2 * P)
                   + M[1, 0] * c2 * np.cos(2 * A) + M[2, 0] * c2 * np.sin(2 * A)
                   + M[1, 1] * c1 * c2 * np.cos(2 * P) * np.cos(2 * A)
                   + M[1, 2] * c1 * c2 * np.sin(2 * P) * np.cos(2 * A)
                   + M[2, 1] * c1 * c2 * np.cos(2 * P) * np.sin(2 * A)
                   + M[2, 2] * c1 * c2 * np.sin(2 * P) * np.sin(2 * A))
        def nz(pair):
            return (pair[0] / dc_expr, pair[1] / dc_expr)

        # 4C: M44 만 실린다
        want = (-0.5 * M[3, 3] * np.sin(d1) * np.sin(d2) * np.cos(pm),
                -0.5 * M[3, 3] * np.sin(d1) * np.sin(d2) * np.sin(pm))
        ok &= check(f"{tag}  4C  = M44 sin d1 sin d2", co[2], nz(want))

        # 16C: 같은 M44 인데 P'+A' 쪽이다
        want = (0.5 * M[3, 3] * np.sin(d1) * np.sin(d2) * np.cos(pp),
                0.5 * M[3, 3] * np.sin(d1) * np.sin(d2) * np.sin(pp))
        ok &= check(f"{tag}  16C = M44 (P'+A' 쪽)", co[8], nz(want))

        # 2C: M24, M34 가 sin d1 * s2 를 타고 실린다
        want = (-0.5 * M[1, 3] * np.sin(d1) * s2 * np.sin(pm)
                + 0.5 * M[2, 3] * np.sin(d1) * s2 * np.cos(pm),
                -0.5 * M[1, 3] * np.sin(d1) * s2 * np.cos(pm)
                - 0.5 * M[2, 3] * np.sin(d1) * s2 * np.sin(pm))
        ok &= check(f"{tag}  2C  = M24, M34", co[1], nz(want))

        # 6C: 4행이 sin d2 를 타고 실린다
        want = (M[3, 0] * np.sin(d2) * np.sin(2 * A)
                + M[3, 1] * c1 * np.sin(d2) * np.cos(2 * P) * np.sin(2 * A)
                + M[3, 2] * c1 * np.sin(d2) * np.sin(2 * P) * np.sin(2 * A),
                -M[3, 0] * np.sin(d2) * np.cos(2 * A)
                - M[3, 1] * c1 * np.sin(d2) * np.cos(2 * P) * np.cos(2 * A)
                - M[3, 2] * c1 * np.sin(d2) * np.sin(2 * P) * np.cos(2 * A))
        ok &= check(f"{tag}  6C  = M41, M42, M43", co[3], nz(want))

        # 10C: 4열이 sin d1 을 타고 실린다
        want = (-M[0, 3] * np.sin(d1) * np.sin(2 * P)
                - M[1, 3] * np.sin(d1) * c2 * np.sin(2 * P) * np.cos(2 * A)
                - M[2, 3] * np.sin(d1) * c2 * np.sin(2 * P) * np.sin(2 * A),
                M[0, 3] * np.sin(d1) * np.cos(2 * P)
                + M[1, 3] * np.sin(d1) * c2 * np.cos(2 * P) * np.cos(2 * A)
                + M[2, 3] * np.sin(d1) * c2 * np.cos(2 * P) * np.sin(2 * A))
        ok &= check(f"{tag}  10C = M14, M24, M34", co[5], nz(want))

        # 32C: 중심 2x2 가 s1*s2 를 타고 실린다
        want = (0.5 * s1 * s2 * (M[1, 1] * np.cos(pp) - M[1, 2] * np.sin(pp)
                                 - M[2, 1] * np.sin(pp) - M[2, 2] * np.cos(pp)),
                0.5 * s1 * s2 * (M[1, 1] * np.sin(pp) + M[1, 2] * np.cos(pp)
                                 + M[2, 1] * np.cos(pp) - M[2, 2] * np.sin(pp)))
        ok &= check(f"{tag}  32C = 중심 2x2", co[16], nz(want))
    return ok


def verify_fourth_column_needs_sin_delta():
    """4열·4행이 정말 sin(delta) 에만 걸리는가 — 2편의 논지를 수치로 확인."""
    print("지연량이 0 이나 180° 면 4행·4열이 사라지는가")
    ok = True
    for d_deg in (0.0, 180.0):
        d = np.deg2rad(d_deg)
        M = random_sample()
        # 4열만 남긴 시편과 4열을 0 으로 만든 시편의 신호가 같아야 한다
        M4 = M.copy()
        M4[:, 3] = 0.0
        I_all = intensity(M, optical_cycle(512), d, np.pi / 2, RATIO)
        I_no4 = intensity(M4, optical_cycle(512), d, np.pi / 2, RATIO)
        ok &= check(f"delta1={d_deg:5.1f}°  4열이 신호에 기여하지 않는다", I_all, I_no4)
        M3 = M.copy()
        M3[3, :] = 0.0
        I_no3 = intensity(M3, optical_cycle(512), np.pi / 2, d, RATIO)
        I_all2 = intensity(M, optical_cycle(512), np.pi / 2, d, RATIO)
        ok &= check(f"delta2={d_deg:5.1f}°  4행이 신호에 기여하지 않는다", I_all2, I_no3)
    return ok


def verify_reduces_to_rce():
    """delta1 = 0 이면 첫 보상자가 사라져 1편의 회전보상자 배치가 되어야 한다."""
    print("구조 검사 — delta1 = 0 이면 1편의 RCE 로 환원되는가")
    M = random_sample()
    C = optical_cycle(512)
    got = intensity(M, C, 0.0, np.pi / 2, RATIO, pol=np.pi / 4, ana=0.0)
    # 보상자 하나짜리를 직접 만들어 비교한다
    mp, ma = mueller_polarizer(np.pi / 4), mueller_polarizer(0.0)
    want = np.array([(ma @ (mueller_retarder(3 * c, np.pi / 2) @ (M @ (mp @ np.array(
        [1.0, 0, 0, 0])))))[0] for c in C])
    return check("보상자 하나짜리 배치와 일치", got, want)


def report_design():
    """설계 수치를 표로 찍는다. 본문 표와 대조하는 용도다."""
    print("회전비별 최고차 고조파 (Collins & Koh 표 1 과 대조)")
    for ratio, expect in ((5, 1), 24), ((5, 2), 28), ((5, 3), 32), ((5, 4), 36):
        got = highest_harmonic(ratio)
        ok = 2 * got == expect
        print(f"  [{'OK  ' if ok else 'FAIL'}] {ratio[0]}:{ratio[1]} -> "
              f"{2*got}C (논문 {expect}C), 최소 적분 {2*got+1}회")
    print("조건수")
    for ratio in ((5, 3), (5, 1)):
        c90 = condition_number(np.pi / 2, np.pi / 2, ratio)
        c128 = condition_number(np.deg2rad(128), np.deg2rad(128), ratio)
        print(f"  {ratio[0]}:{ratio[1]}  delta=90°: {c90:.2f}   "
              f"delta=128°: {c128:.2f}   개선 {c90/c128:.2f}배")
    same = np.allclose(condition_number(np.deg2rad(128), np.deg2rad(128), (5, 3)),
                       condition_number(np.deg2rad(128), np.deg2rad(128), (5, 1)))
    print(f"  [{'OK  ' if same else 'FAIL'}] 5:3 과 5:1 의 조건수가 같다")
    return same


if __name__ == "__main__":
    results = [verify_vanishing(), verify_table2(),
               verify_fourth_column_needs_sin_delta(), verify_reduces_to_rce(),
               report_design()]
    print("\n전체 결과:", "모두 통과" if all(results) else "실패 항목 있음")
    raise SystemExit(0 if all(results) else 1)

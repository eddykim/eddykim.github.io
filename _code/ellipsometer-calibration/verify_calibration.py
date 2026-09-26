"""calibration.py 를 독립된 경로들과 교차검증한다.

실행: python verify_calibration.py

1~3편에서 참고문헌 식을 옮겨 적다 생긴 오독을 이 방식으로 잡았다.
여기서도 한쪽 결과를 다른 쪽에서 유도하지 않는다.

  신호 합성   : 부품 뮬러 행렬의 곱
  대조 대상 A : Fujiwara 식 4.58-4.60 (As 회전 + eta 배율)
  대조 대상 B : Johs 식 2-5 (model_ab)
  왕복       : 넣은 Ps, As, eta 를 잔차 보정으로 되찾는가
  두 보정    : 잔차 보정 결과와 회귀 보정(LM) 결과가 일치하는가
  환원       : 오차를 0 으로 두면 1편의 이상적 RAE 가 되는가
  역산 왕복  : 이상적 역산이 참 Psi, Delta 를 돌려주는가

잔차 보정의 왕복은 완벽하게 닫히지 않는다. R(P) 가 포물선이 아니기 때문이다.
이것을 오차로 보지 않고, 훑는 범위를 좁힐 때 오차가 함께 줄어드는지로 확인한다 —
줄어들면 원인이 포물선 근사임이 확정되고, 그것이 회귀 보정이 필요한 이유가 된다.
"""
import numpy as np

from calibration import (invert_ideal, measure, model_ab,
                         regression_calibration, residual, residual_calibration)

TOL = 1e-9
D = np.deg2rad
N_TRIAL = 40           # 맹점 시험의 잡음 반복 횟수

# Delta 가 90도 부근이라 잔차 보정이 잘 듣는 시편과, 0도 부근이라 안 듣는 시편
GOOD = (D(16.8), D(87.0))      # Johs 표 1 의 SiO2/Si 와 비슷한 자리
BAD = (D(27.0), D(2.0))        # Johs 표 1 의 7059 유리와 비슷한 자리


def check(name, got, want, tol=TOL):
    err = np.max(np.abs(np.asarray(got) - np.asarray(want)))
    ok = err < tol
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<52s} 최대오차 {err:.3e}")
    return ok


def verify_ideal_reduces_to_post1():
    """오차가 없으면 1편의 RAE 식 I = I0[1 + S1 cos2A + S2 sin2A] 가 되는가."""
    print("환원 — 오차 0 이면 1편의 이상적 RAE")
    ok = True
    for psi, delta in (GOOD, BAD, (D(45), D(-37))):
        a, b = measure(psi, delta, D(45.0))          # 편광자 45도
        S1, S2 = -np.cos(2 * psi), np.sin(2 * psi) * np.cos(delta)
        ok &= check(f"psi={np.rad2deg(psi):5.1f}° D={np.rad2deg(delta):6.1f}° "
                    f"-> (alpha,beta)=(S1,S2)", (a, b), (S1, S2))
    return ok


def verify_fujiwara_rotation():
    """Fujiwara 식 4.59 — 측정값이 참값에 As 회전과 eta 배율을 씌운 꼴인가."""
    print("Fujiwara 식 4.59 — As 회전 + eta 배율")
    ok = True
    for a_s_deg, eta in ((0.0, 1.0), (12.0, 0.97), (-30.0, 0.85), (64.3, 0.9687)):
        a_s = D(a_s_deg)
        for psi, delta in (GOOD, BAD):
            a0, b0 = measure(psi, delta, D(45.0))                      # 참값
            a1, b1 = measure(psi, delta, D(45.0), a_s=a_s, eta=eta)    # 측정값
            want = (eta * (a0 * np.cos(2 * a_s) - b0 * np.sin(2 * a_s)),
                    eta * (a0 * np.sin(2 * a_s) + b0 * np.cos(2 * a_s)))
            ok &= check(f"As={a_s_deg:6.1f}° eta={eta:.4f}  D={np.rad2deg(delta):5.1f}°",
                        (a1, b1), want)
    return ok


def verify_inverse_matrix():
    """식 4.60 이 식 4.59 의 역인가 — 되돌리면 참값이 나오는가."""
    print("Fujiwara 식 4.60 — 4.59 의 역행렬")
    ok = True
    for a_s_deg, eta in ((25.0, 0.93), (-40.0, 0.88)):
        a_s = D(a_s_deg)
        psi, delta = GOOD
        a0, b0 = measure(psi, delta, D(45.0))
        a1, b1 = measure(psi, delta, D(45.0), a_s=a_s, eta=eta)
        back = ((a1 * np.cos(2 * a_s) + b1 * np.sin(2 * a_s)) / eta,
                (-a1 * np.sin(2 * a_s) + b1 * np.cos(2 * a_s)) / eta)
        ok &= check(f"As={a_s_deg:6.1f}° eta={eta:.2f}  되돌리면 참값", back, (a0, b0))
    return ok


def verify_johs_model():
    """Johs 식 2-5 가 뮬러 곱으로 만든 신호와 맞는가 — 서로 다른 경로다."""
    print("Johs 식 2-5 — 닫힌 식 vs 뮬러 곱")
    ok = True
    for ps_deg, a_s_deg, eta in ((0, 0, 1.0), (42.3, 64.3, 0.9687), (-12, 25, 0.9)):
        ps, a_s = D(ps_deg), D(a_s_deg)
        for psi, delta in (GOOD, BAD):
            P = D(np.array([10.0, 30.0, 55.0, 80.0])) + ps
            got = np.array([measure(psi, delta, p, ps=ps, a_s=a_s, eta=eta) for p in P]).T
            want = model_ab(P, psi, delta, ps, a_s, eta)
            ok &= check(f"Ps={ps_deg:5.1f}° As={a_s_deg:5.1f}° D={np.rad2deg(delta):5.1f}°",
                        got, want, tol=1e-8)
    return ok


def verify_inversion_roundtrip():
    """이상적 역산(식 4.16-4.17)이 참 Psi, Delta 를 돌려주는가.

    이 식은 그림 4 에서 '보정하지 않은 눈금값을 그대로 넣는' 용도로 쓴다.
    그 전에 오차가 없을 때 정확히 닫히는지 확인해 둔다.
    """
    print("역산 왕복 — 오차 0 이면 참 Psi, Delta 가 돌아온다")
    ok = True
    for psi_d, dl_d in ((16.8, 87.0), (27.0, 2.0), (45.0, 120.0), (8.0, 160.0)):
        got = []
        for P_d in (30.0, 45.0, 70.0):
            a, b = measure(D(psi_d), D(dl_d), D(P_d))
            psi, delta = invert_ideal(a, b, D(P_d))
            got.append((np.rad2deg(psi), np.rad2deg(delta)))
        ok &= check(f"psi={psi_d:5.1f}° D={dl_d:6.1f}°  편광자 3 자리에서 동일",
                    got, [(psi_d, dl_d)] * 3, tol=1e-8)
    return ok


def verify_locus_ellipse():
    """편광자를 훑을 때 (alpha', beta') 가 그리는 타원의 반축이 1 과 |cos Delta| 인가.

    본문 2절이 이 기하로 세 오차를 가른다. u = tan(P-Ps) 를 u = tan(Psi) tan(th)
    로 두면 alpha0 = cos 2th, beta0 = cos(Delta) sin 2th 가 되어 Psi 가 사라진다.
    """
    print("자리의 기하 — (alpha', beta') 타원의 반축이 (1, |cos Delta|) 인가")
    ok = True
    P = np.linspace(-np.pi / 2, np.pi / 2, 2001)
    for psi_d, dl_d in ((35.0, 50.0), (16.8, 87.0), (27.0, 2.0), (45.0, 120.0)):
        ab = np.array([measure(D(psi_d), D(dl_d), p) for p in P])
        got = (np.abs(ab[:, 0]).max(), np.abs(ab[:, 1]).max())
        ok &= check(f"psi={psi_d:5.1f}° D={dl_d:6.1f}°  반축 (1, |cos D|)",
                    got, (1.0, abs(np.cos(D(dl_d)))), tol=1e-6)
    return ok


def verify_residual_roundtrip():
    """넣은 Ps, As, eta 를 잔차 보정이 되찾는가.

    Ps 와 As 는 R(P) 의 최소 '위치'에서 나오므로 포물선 근사에 둔감하다.
    eta 는 최소 '값'에서 나오므로 민감하다 — 아래에서 따로 다룬다.
    """
    print("왕복 — 잔차 보정이 주입값을 되찾는가 (훑는 범위 ±2°)")
    ok = True
    psi, delta = GOOD
    for ps_deg, a_s_deg, eta in ((0.0, 0.0, 1.0), (0.26, 30.0, 0.97),
                                 (-1.5, 64.3, 0.9687)):
        r = residual_calibration(psi, delta, np.arange(-2, 2.01, 0.2) + ps_deg,
                                 ps=D(ps_deg), a_s=D(a_s_deg), eta=eta)
        ok &= check(f"Ps  주입 {ps_deg:6.2f}° -> 회수 {np.rad2deg(r['ps']):6.2f}°",
                    np.rad2deg(r['ps']), ps_deg, tol=5e-3)
        ok &= check(f"eta 주입 {eta:.4f}  -> 회수 {r['eta']:.4f}", r['eta'], eta, tol=2e-4)
        ok &= check(f"As  주입 {a_s_deg:6.2f}° -> 회수 {np.rad2deg(r['a_s']):6.2f}°",
                    np.rad2deg(r['a_s']), a_s_deg, tol=5e-3)
    return ok


def verify_parabolic_approximation():
    """남는 오차가 포물선 근사 탓인가 — 범위를 좁히면 함께 줄어야 한다."""
    print("포물선 근사 — 훑는 범위를 좁히면 eta 오차가 줄어드는가")
    psi, delta = GOOD
    ps, a_s, eta = D(0.26), D(64.3), 0.9687
    print(f"{'범위':>8} {'eta 회수':>10} {'eta 오차':>11} {'Ps 오차':>11}")
    errs = []
    for half in (5.0, 3.0, 2.0, 1.0, 0.5):
        r = residual_calibration(psi, delta, np.arange(-half, half + 1e-9, half / 10) +
                                 np.rad2deg(ps), ps=ps, a_s=a_s, eta=eta)
        e = abs(r['eta'] - eta)
        errs.append(e)
        print(f"   ±{half:4.1f}° {r['eta']:10.6f} {e:11.3e} "
              f"{abs(np.rad2deg(r['ps']) - np.rad2deg(ps)):10.5f}°")
    # 단조 감소 여부로 판정한다
    ok = all(b < a for a, b in zip(errs, errs[1:]))
    print(f"  [{'OK  ' if ok else 'FAIL'}] 오차가 단조 감소한다 "
          f"({errs[0]:.1e} -> {errs[-1]:.1e}, {errs[0] / errs[-1]:.0f} 배)")
    # 구간이 최소점에 대칭이 아니면 짝수차 항이 상쇄되지 않아 Ps 도 흔들린다
    sym = residual_calibration(psi, delta, np.arange(-5, 5.01, 0.25) + np.rad2deg(ps),
                               ps=ps, a_s=a_s, eta=eta)
    off = residual_calibration(psi, delta, np.arange(-5, 5.01, 0.25),
                               ps=ps, a_s=a_s, eta=eta)
    e_sym = abs(np.rad2deg(sym['ps']) - np.rad2deg(ps))
    e_off = abs(np.rad2deg(off['ps']) - np.rad2deg(ps))
    print(f"         ±5° 에서 Ps 오차: 최소점 대칭 {e_sym:.6f}°, 눈금 0 중심 {e_off:.6f}°")
    ok &= e_sym < 1e-6 < e_off
    # R(P) 자체가 포물선이 아님을 4차항으로 확인한다
    P = D(np.arange(-5, 5.01, 0.5)) + ps
    R = np.array([residual(psi, delta, p, ps=ps, a_s=a_s, eta=eta) for p in P])
    r2 = np.max(np.abs(np.polyval(np.polyfit(P, R, 2), P) - R))
    r4 = np.max(np.abs(np.polyval(np.polyfit(P, R, 4), P) - R))
    print(f"         ±5° 에서 R(P) 맞춤 잔차: 2차 {r2:.3e}, 4차 {r4:.3e}")
    ok &= r4 < r2 / 10
    return ok


def verify_regression_roundtrip():
    """회귀 보정이 다섯 값을 되찾는가 — 잔차가 듣지 않는 시편에서도."""
    print("왕복 — 회귀 보정이 다섯 값을 되찾는가")
    ok = True
    for tag, (psi, delta) in (("Delta 좋음", GOOD), ("Delta 나쁨", BAD)):
        ps, a_s, eta = D(0.26), D(64.3), 0.9687
        P = D(np.linspace(5.0, 85.0, 90)) + ps
        a_exp, b_exp = np.array([measure(psi, delta, p, ps=ps, a_s=a_s, eta=eta)
                                 for p in P]).T
        res = regression_calibration(P, a_exp, b_exp,
                                     guess=(D(20), D(60), 0.0, D(60), 0.95))
        got = res['params']
        want = (psi, delta, ps, a_s, eta)
        ok &= check(f"{tag}  다섯 파라미터 회수", got, want, tol=1e-6)
        print(f"         MSE = {res['mse']:.3e}  (잡음이 없으므로 0 에 가까워야 한다)")
    return ok


def verify_blind_spot():
    """Delta 가 0·180° 근처면 잔차 보정이 실제로 무너지는가.

    잡음이 없으면 곡률이 아무리 작아도 최소점은 정확히 찾힌다. 맹점은 곡률이
    잡음에 묻힐 때 드러나므로, 검출 세기에 가우스 잡음을 넣고 반복 시행한다.
    """
    print("맹점 — Delta 가 0·180° 근처면 잔차 보정이 무너진다 (잡음 sigma=2e-4)")
    psi, ps = D(27.0), D(0.26)
    print(f"{'Delta':>8} {'R(P) 곡률':>13} {'Ps 오차 중앙값':>16} {'맞춤 실패':>10}")
    ok = True
    rows = []
    for d_deg in (90.0, 60.0, 30.0, 10.0, 3.0, 1.0):
        errs, fail = [], 0
        for k in range(N_TRIAL):
            r = residual_calibration(psi, D(d_deg), np.arange(-5, 5.01, 0.5) + 0.26,
                                     ps=ps, noise=2e-4,
                                     rng=np.random.default_rng(1000 * int(d_deg) + k))
            if r['ok']:
                errs.append(abs(np.rad2deg(r['ps']) - 0.26))
            else:
                fail += 1
        clean = residual_calibration(psi, D(d_deg), np.arange(-5, 5.01, 0.5) + 0.26,
                                     ps=ps)
        med = float(np.median(errs)) if errs else np.nan
        rows.append((d_deg, clean['curvature'], med, fail))
        print(f"{d_deg:7.1f}° {clean['curvature']:13.3e} {med:15.4f}° "
              f"{f'{fail}/{N_TRIAL}':>10}")
    # 곡률이 무너지면 Ps 오차가 커져야 한다
    ok &= rows[0][2] < 0.02 and rows[-1][2] > 20 * rows[0][2]
    print(f"  [{'OK  ' if ok else 'FAIL'}] Delta 90° -> 1° 에서 Ps 오차가 "
          f"{rows[-1][2] / rows[0][2]:.0f} 배로 커진다")
    # 회귀 보정은 같은 자리에서도 버텨야 한다
    a_s, eta = D(64.3), 0.9687
    P = D(np.linspace(5.0, 85.0, 90)) + ps
    a_exp, b_exp = np.array([measure(psi, D(1.0), p, ps=ps, a_s=a_s, eta=eta)
                             for p in P]).T
    res = regression_calibration(P, a_exp, b_exp, guess=(D(20), D(10), 0.0, D(60), 0.95))
    ok &= check("  Delta=1° 에서도 회귀 보정은 Ps 를 되찾는다",
                res['params'][2], ps, tol=1e-6)
    return ok


def verify_model_mismatch():
    """모형에 없는 효과가 eta 를 오염시키는가 — Johs 표 2 의 재현."""
    print("모형 불일치 — 넣지 않은 PDS 가 eta 를 오염시킨다")
    psi, delta = D(16.8), D(87.0)
    ps, a_s, eta = D(0.26), D(64.3), 1.0
    P = D(np.linspace(5.0, 85.0, 90)) + ps
    a_exp, b_exp = np.array([measure(psi, delta, p, ps=ps, a_s=a_s, eta=eta, pds=0.25)
                             for p in P]).T
    res = regression_calibration(P, a_exp, b_exp, guess=(D(20), D(80), 0.0, D(60), 0.95))
    eta_fit, psi_fit = res['params'][4], np.rad2deg(res['params'][0])
    print(f"  참값     : eta = {eta:.4f}, psi = {np.rad2deg(psi):.3f}°")
    print(f"  PDS 무시 : eta = {eta_fit:.4f}, psi = {psi_fit:.3f}°, MSE = {res['mse']:.3e}")
    drift = abs(eta_fit - eta)
    ok = drift > 1e-3
    print(f"  [{'OK  ' if ok else 'FAIL'}] eta 가 참값에서 {drift:.4f} 만큼 밀렸다")
    return ok


if __name__ == "__main__":
    results = [verify_ideal_reduces_to_post1(), verify_fujiwara_rotation(),
               verify_inverse_matrix(), verify_johs_model(),
               verify_inversion_roundtrip(), verify_locus_ellipse(),
               verify_residual_roundtrip(), verify_parabolic_approximation(),
               verify_regression_roundtrip(),
               verify_blind_spot(), verify_model_mismatch()]
    print("\n전체 결과:", "모두 통과" if all(results) else "실패 항목 있음")
    raise SystemExit(0 if all(results) else 1)

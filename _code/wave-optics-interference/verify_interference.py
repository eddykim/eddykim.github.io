"""interference.py 를 독립된 경로들과 교차검증한다.

실행: python verify_interference.py

  A. 두 파동의 간섭           ↔  B&W 7.2 식 (15)–(17), 공간 평균 = I1 + I2 (에너지 보존), 가시도 2√(I1I2)/(I1+I2)
  B. 영의 무늬 간격           ↔  λL/d (Hecht 9.3, Saleh & Teich 연습문제 2.5-2)
  C. 박막 반사율 (급수 합)    ↔  전달행렬법 ↔ opticore multilayer_rt (다른 사람이 짠 Abelès 행렬)
  D. 뉴턴 링                  ↔  어두운 고리 반지름 √(mλR) (Hecht 9.4)
  E. 빔스플리터와 마이컬슨     ↔  유니타리, 두 출력의 합 = 입력, 반사 위상 π/2 (Saleh & Teich 2.5)
  F. M 개 빔                  ↔  정점 M², 평균 M, 첫 영점 2π/M (Saleh & Teich 식 2.5-12)
  G. 파브리-페로              ↔  에어리 함수 (1−R)²/(1+R²−2R cosφ), 피네스 π√R/(1−R) (큰 R), 판의 전달행렬
  H. 수차의 간섭무늬          ↔  어두운 고리 수 = 파면수차의 파장 수
  I. LIGO 수치                ↔  Saleh & Teich 예제 2.5-1 (2F/π ≈ 286, Δd ≈ 2 am)
"""
import json
import os
import subprocess

import numpy as np

from interference import (beamsplitter, etalon_transmission, fwhm_of_peak, interferogram, michelson_outputs,
                          multibeam_sum, newton_gap, superpose, thin_layer_reflectance, tmm_normal, visibility_of,
                          young_screen)

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "repo_opticore")


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<62s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def verify_two_wave():
    print("A. 두 파동의 간섭")
    d = np.linspace(0, 2 * np.pi, 4001)[:-1]
    ok = True
    for I1, I2 in [(1.0, 1.0), (1.0, 0.25), (2.0, 0.01)]:
        I = superpose(I1, I2, d)
        ok &= check(f"I1 = {I1}, I2 = {I2}: I1 + I2 + 2√(I1I2) cosδ", I, I1 + I2 + 2 * np.sqrt(I1 * I2) * np.cos(d), 1e-12)
        ok &= check("  위상 평균 = I1 + I2 (에너지는 자리만 옮긴다)", I.mean(), I1 + I2, 1e-12)
        ok &= check("  가시도 = 2√(I1I2)/(I1+I2)", visibility_of(I), 2 * np.sqrt(I1 * I2) / (I1 + I2), 1e-6)
    return ok


def verify_young():
    print("B. 영의 무늬 간격 ↔ λL/d")
    ok = True
    for dd, LL, lam in [(500.0, 1e6, 0.55), (100.0, 5e5, 0.633), (2000.0, 2e6, 0.45)]:
        p = lam * LL / dd
        x = np.linspace(-3 * p, 3 * p, 600001)
        I = young_screen(x, dd, LL, lam)
        mins = x[1:-1][(I[1:-1] < I[:-2]) & (I[1:-1] < I[2:])]
        # λL/d 는 각이 작을 때의 근사다. 가운데 두 어두운 줄 사이로 잰다 (바깥으로 갈수록 각의 제곱만큼 벌어진다)
        centre = np.sort(mins[np.argsort(np.abs(mins))[:2]])
        ok &= check(f"d = {dd / 1000} mm, L = {LL / 1e6} m, λ = {lam} µm: 가운데 무늬 간격 (상대)", np.diff(centre)[0] / p - 1, 0.0, 1e-5)
    return ok


OPTICORE = r"""
import json, sys, numpy as np
from opticore.polarization.thinfilm import multilayer_rt
spec = json.loads(sys.argv[1])
out = []
for n0, n1, n2, t, lam in spec:
    rs, rp = multilayer_rt(complex(n0), complex(n2), [(complex(n1), t)], np.array([0.0]), lam, "reflect")
    out.append(float(abs(rs[0]) ** 2))
print(json.dumps(out))
"""


def verify_thin_film():
    print("C. 박막 반사율: 급수 합 ↔ 전달행렬 ↔ opticore")
    cases = [(1.0, 1.33, 1.0, t, 0.5893) for t in (0.05, 0.11, 0.3, 0.77)] + \
            [(1.52, 1.0, 1.52, t, 0.5893) for t in (0.1, 0.4)] + [(1.0, 1.38, 1.52, 0.1068, 0.5893)]
    R_series = np.array([thin_layer_reflectance(*c) for c in cases])
    R_tmm = np.array([abs(tmm_normal([c[0], c[1], c[2]], [c[3]], c[4])[0]) ** 2 for c in cases])
    ok = check("다중반사 급수 ↔ 전달행렬", R_series, R_tmm, 1e-12)
    py = os.path.join(REPO, "opticore", ".venv", "bin", "python")
    if os.path.exists(py):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        out = subprocess.run([py, "-c", OPTICORE, json.dumps(cases)], capture_output=True, text=True, env=env, check=True)
        R_oc = np.array(json.loads(out.stdout.strip().splitlines()[-1]))
        ok &= check("다중반사 급수 ↔ opticore multilayer_rt", R_series, R_oc, 1e-10)
    else:
        print(f"  [SKIP] opticore 환경이 없다 ({py})")
    # MgF2 (n = 1.38) λ/4 무반사막: 2 빔 근사로 최소 반사율 위치 = λ/(4n)
    t = np.linspace(0.05, 0.16, 11001)
    Rm = thin_layer_reflectance(1.0, 1.38, 1.52, t, 0.5893)
    ok &= check("MgF₂ 무반사막: 최소 반사율 두께 = λ/(4n)", t[np.argmin(Rm)], 0.5893 / (4 * 1.38), 2e-5)
    return ok


def verify_newton():
    print("D. 뉴턴 링 ↔ √(mλR)")
    lam, R = 0.5893, 1.0e6
    r = np.linspace(0, 4000, 800001)
    Rr = thin_layer_reflectance(1.52, 1.0, 1.52, newton_gap(r, R), lam)
    mins = r[1:-1][(Rr[1:-1] < Rr[:-2]) & (Rr[1:-1] < Rr[2:])]
    m = np.arange(1, 9)
    return check("어두운 고리 m = 1…8 반지름 (상대)", mins[:8] / np.sqrt(m * lam * R) - 1, 0.0, 1e-4)


def verify_michelson():
    print("E. 빔스플리터와 마이컬슨")
    ok = True
    for tt in [1 / np.sqrt(2), 0.6, 0.9]:
        B = beamsplitter(tt)
        ok &= check(f"t = {tt:.3f}: 산란 행렬이 유니타리", B.conj().T @ B, np.eye(2), 1e-12)
        a, b = michelson_outputs(np.linspace(0, 4 * np.pi, 101), tt)
        ok &= check(f"t = {tt:.3f}: 두 출력의 합 = 1", a + b, 1.0, 1e-12)
    B = beamsplitter()
    ok &= check("반사와 투과의 위상차 = π/2", np.angle(B[0, 1] / B[0, 0]), np.pi / 2, 1e-12)
    a, b = michelson_outputs(np.linspace(0, 2 * np.pi, 1001))
    ok &= check("50:50 일 때 한쪽 출력은 완전히 0 까지 (가시도 1)", visibility_of(a), 1.0, 1e-9)
    return ok


def verify_multibeam():
    print("F. M 개 빔")
    ok = True
    phi = np.linspace(-np.pi, np.pi, 200001)
    for M in [2, 5, 20]:
        I = multibeam_sum(phi, 1.0, M)
        ok &= check(f"M = {M}: 정점 = M²", I.max(), M**2, 1e-6)
        ok &= check(f"M = {M}: 위상 평균 = M", I[:-1].mean(), M, 1e-6)
        ok &= check(f"M = {M}: 첫 영점 = 2π/M", multibeam_sum(2 * np.pi / M, 1.0, M), 0.0, 1e-20)
    return ok


def verify_fabry_perot():
    print("G. 파브리-페로")
    ok = True
    phi = np.linspace(-np.pi, np.pi, 400001)
    for R in [0.04, 0.5, 0.9, 0.98]:
        T = etalon_transmission(phi, R)
        ok &= check(f"R = {R}: 빔 합 ↔ 에어리 (1−R)²/(1+R²−2R cosφ)", T, (1 - R) ** 2 / (1 + R**2 - 2 * R * np.cos(phi)), 1e-9)
    for R in [0.9, 0.98]:
        T = etalon_transmission(phi, R)
        F_num = 2 * np.pi / fwhm_of_peak(phi, T)
        ok &= check(f"R = {R}: 피네스 {F_num:.2f} ↔ π√R/(1−R) (상대)", F_num / (np.pi * np.sqrt(R) / (1 - R)) - 1, 0.0, 2e-3)
    # 실리콘 판 (n = 3.5) 의 투과율: 전달행렬 ↔ 에어리 + 프레넬 R
    n, lam = 3.5, 1.55
    R = ((n - 1) / (n + 1)) ** 2
    ts = np.linspace(100.0, 100.5, 501)
    T_tmm = np.array([abs(tmm_normal([1.0, n, 1.0], [t], lam)[1]) ** 2 for t in ts])
    T_airy = (1 - R) ** 2 / (1 + R**2 - 2 * R * np.cos(4 * np.pi * n * ts / lam))
    ok &= check(f"실리콘 판 (R = {R:.3f}): 전달행렬 ↔ 에어리", T_tmm, T_airy, 1e-10)
    return ok


def verify_interferograms():
    print("H. 수차의 간섭무늬: 어두운 고리 수 = 파면수차 (파장)")
    rho = np.linspace(0, 1, 200001)
    ok = True
    for name, W, n_expect in [("초점 2λ", 2 * rho**2, 2), ("구면수차 3λ", 3 * rho**4, 3), ("초점 5.0λ", 5.0 * rho**2, 5)]:
        I = interferogram(W)
        dark = I < 0.05
        n_dark = int(np.sum(dark[1:] & ~dark[:-1]))      # 어두운 띠로 들어가는 횟수
        ok &= check(f"{name}: 어두운 고리 {n_dark}", n_dark, n_expect, 0.5)
    return ok


def verify_ligo():
    print("I. LIGO 수치 (Saleh & Teich 예제 2.5-1)")
    F = 450.0
    ok = check("팔의 파브리-페로가 위상을 키우는 배율 2F/π ≈ 286", 2 * F / np.pi, 286, 1.0)
    dd = 5e-22 * 4e3
    ok &= check("변형 5×10⁻²² × 4 km = 2 am", dd, 2e-18, 1e-20)
    ok &= check("양성자 반지름(0.84 fm)보다 약 400배 작다", 0.84e-15 / dd, 420, 30)
    return ok


def main():
    results = [f() for f in (verify_two_wave, verify_young, verify_thin_film, verify_newton, verify_michelson,
                             verify_multibeam, verify_fabry_perot, verify_interferograms, verify_ligo)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

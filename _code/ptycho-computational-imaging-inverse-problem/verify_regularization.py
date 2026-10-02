"""1편 계산의 교차검증. 한쪽 결과를 다른 쪽에서 유도하지 않는 두 경로를 대조한다.

  1. 티호노프     정규방정식 (선형계)  <->  SVD 필터 인자
  2. 위너          순환 행렬 정규방정식  <->  푸리에 영역 나눗셈
  3. 랜드웨버      반복 계산  <->  SVD 필터 인자 1 - (1 - w s^2)^k 의 닫힌 해
  4. TV-ADMM      ADMM  <->  같은 문제를 범용 제약 최적화기(trust-constr)로 푼 목적함수 값
  5. 피카르 조건   잡음의 특잇벡터 성분 평균 |u_i . e| = sigma sqrt(2/pi) (Hansen 식 2.16 의 가우시안 판)
  6. 부호 커널     전달함수에 영점이 없고 상자에는 있다 (Raskar 2006 의 주장)
  7. 불일치 원리   고른 lam 의 잔차 = sigma sqrt(m), 그 오차가 최적 lam 오차의 10 % 이내 (여러 난수)
"""
import numpy as np
from scipy.optimize import minimize, LinearConstraint

from forward_models import (box_kernel, coded_kernel, smear_matrix, circulant_matrix,
                            transfer_function, test_signal, measure)
from regularize import (tikhonov, tikhonov_svd, wiener_circulant, landweber, landweber_svd,
                        tv_admm, tv_objective, difference_matrix,
                        discrepancy_lambda, best_lambda, rel_error)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
x = test_signal()
n = x.size
Ab = smear_matrix(box_kernel(), n)
y = measure(Ab, x, 0.01, rng)

# 1. 티호노프 두 경로
for lam in (1e-3, 0.05, 0.5):
    d = np.abs(tikhonov(Ab, y, lam) - tikhonov_svd(Ab, y, lam)).max()
    check(f"티호노프 정규방정식 = SVD (lam={lam})", d < 1e-8, f"최대 차이 {d:.1e}")

# 2. 위너 = 순환 티호노프
N = 256
C = circulant_matrix(box_kernel(), N)
yc = C @ rng.normal(size=N) + 0.01 * rng.normal(size=N)
d = np.abs(tikhonov(C, yc, 0.05) - wiener_circulant(box_kernel(), yc, 0.05)).max()
check("위너 필터 = 순환 합성곱 티호노프", d < 1e-10, f"최대 차이 {d:.1e}")

# 3. 랜드웨버 반복 = 닫힌 해
for k in (10, 300):
    d = np.abs(landweber(Ab, y, k)[0] - landweber_svd(Ab, y, k)).max()
    check(f"랜드웨버 {k}회 반복 = SVD 닫힌 해", d < 1e-9, f"최대 차이 {d:.1e}")

# 4. TV: ADMM 목적함수 값을 제약 최적화기의 값과 비교 (작은 문제)
m = 60
xs = np.zeros(m); xs[10:25] = 1; xs[35:45] = 0.5
As = smear_matrix(box_kernel(8), m)
ys = measure(As, xs, 0.01, rng)
lam = 0.02
x_admm = tv_admm(As, ys, lam, n_iter=5000)
D = difference_matrix(m)
# 변수 [x, t], 목적 1/2||A x - y||^2 + lam sum t, 조건 -t <= D x <= t (t 가 |Dx| 의 상한)
def f(v):
    r = As @ v[:m] - ys
    return 0.5 * r @ r + lam * v[m:].sum()
def g(v):
    return np.r_[As.T @ (As @ v[:m] - ys), lam * np.ones(m - 1)]
I = np.eye(m - 1)
cons = LinearConstraint(np.block([[D, -I], [-D, -I]]), -np.inf, 0)
v0 = np.r_[np.linalg.lstsq(As, ys, rcond=None)[0], np.ones(m - 1)]
res = minimize(f, v0, jac=g, constraints=[cons], method="trust-constr",
               options={"maxiter": 5000, "gtol": 1e-12, "xtol": 1e-14})
fa, fq = tv_objective(As, ys, x_admm, lam), tv_objective(As, ys, res.x[:m], lam)
check("TV-ADMM 목적함수 = 제약 최적화기", abs(fa - fq) / fq < 1e-4,
      f"ADMM {fa:.8f}, trust-constr {fq:.8f}, 해 차이 {np.abs(x_admm - res.x[:m]).max():.1e}")

# 5. 피카르: 잡음의 특잇벡터 성분 크기 평균
U, s, Vt = np.linalg.svd(Ab, full_matrices=False)
sig = 0.01
coef = np.abs(U.T @ (sig * rng.normal(size=(Ab.shape[0], 400)))).mean()
check("잡음 바닥 |u.e| 평균 = sigma sqrt(2/pi)", abs(coef / (sig * np.sqrt(2 / np.pi)) - 1) < 0.02,
      f"{coef:.5f} vs {sig * np.sqrt(2 / np.pi):.5f}")

# 6. 전달함수 영점
_, Hb = transfer_function(box_kernel())
_, Hc = transfer_function(coded_kernel())
check("상자 커널은 영점, 부호 커널은 영점 없음", Hb.min() < 1e-12 and Hc.min() > 0.02,
      f"min|H| 상자 {Hb.min():.1e}, 부호 {Hc.min():.3f}")
sb = np.linalg.svd(Ab, compute_uv=False)
sc = np.linalg.svd(smear_matrix(coded_kernel(), n), compute_uv=False)
print(f"       조건수 상자 {sb[0] / sb[-1]:.0f}, 부호 {sc[0] / sc[-1]:.1f}")

# 7. 불일치 원리
worst = 0.0
for seed in range(5):
    for sg in (1e-3, 1e-2, 3e-2):
        yy = measure(Ab, x, sg, np.random.default_rng(100 + seed))
        ld = discrepancy_lambda(Ab, yy, sg)
        r = np.linalg.norm(Ab @ tikhonov(Ab, yy, ld) - yy)
        assert abs(r / (sg * np.sqrt(Ab.shape[0])) - 1) < 1e-6
        eb = best_lambda(lambda A, v, l: tikhonov(A, v, l), Ab, yy, x, np.logspace(-5, 0, 80))[1]
        worst = max(worst, rel_error(tikhonov(Ab, yy, ld), x) / eb - 1)
check("불일치 원리 lam 의 오차가 최적의 10 % 이내 (15 조건)", worst < 0.10, f"최악 +{100 * worst:.1f} %")

print("\n전부 통과" if ok else "\n실패 있음")

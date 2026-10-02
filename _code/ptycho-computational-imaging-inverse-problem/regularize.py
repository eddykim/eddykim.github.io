"""역문제를 푸는 방법들 — 계산 이미징과 타이코그래피 1편.

  최소자승        || A x - y || 를 최소로. 잡음을 믿는 만큼 그대로 되돌린다.
  티호노프        || A x - y ||^2 + lam^2 || L x ||^2. 해가 작거나(L = I) 매끈하다고(L = 미분) 가정한다.
  위너 필터       주기 경계에서 티호노프를 푸리에 영역으로 옮긴 것. 주파수마다 나눗셈 한 번.
  랜드웨버        최소자승 경사하강. 일찍 멈추면 티호노프와 비슷한 정규화가 된다.
  TV (ADMM)       || A x - y ||^2 / 2 + lam || D x ||_1. 기울기가 드물다(계단형)고 가정한다.

모든 함수는 번진 신호 y 에서 원 신호 x 의 추정을 돌려준다.
"""
import numpy as np


def least_squares(A, y):
    return np.linalg.lstsq(A, y, rcond=None)[0]


def tikhonov(A, y, lam, L=None):
    """정규방정식 (A^T A + lam^2 L^T L) x = A^T y 를 푼다. y 가 2D 이면 열마다 푼다."""
    n = A.shape[1]
    L = np.eye(n) if L is None else L
    return np.linalg.solve(A.T @ A + lam**2 * (L.T @ L), A.T @ y)


def tikhonov_svd(A, y, lam):
    """같은 티호노프 해를 SVD 필터 인자 f_i = s_i^2 / (s_i^2 + lam^2) 로 쓴 것.

        x = sum_i f_i (u_i . y / s_i) v_i

    최소자승 해(f_i = 1)에서 작은 특잇값 쪽 항을 누르는 것이 정규화의 전부라는 것을 보여준다.
    """
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    f = s**2 / (s**2 + lam**2)
    return Vt.T @ (f * (U.T @ y) / s)


def wiener_circulant(kernel, y, lam):
    """주기 경계 합성곱에서 티호노프(L = I) 해. conj(H) / (|H|^2 + lam^2) 를 곱한다.

    Goodman 의 위너 필터에서 잡음·신호 전력 스펙트럼 비가 주파수와 무관한 상수
    lam^2 인 경우와 같다.
    """
    n = y.size
    H = np.fft.fft(kernel, n)
    return np.real(np.fft.ifft(np.conj(H) / (np.abs(H) ** 2 + lam**2) * np.fft.fft(y)))


def difference_matrix(n):
    """(n - 1) x n 전방 차분 행렬. D x 는 이웃 화소의 차이다."""
    return (np.eye(n, k=1) - np.eye(n))[:-1]


def landweber(A, y, n_iter, record=(), step=None):
    """x_{k+1} = x_k + w A^T (y - A x_k). 최소자승 목적함수의 경사하강이다.

    w < 2 / s_max^2 이면 수렴한다. record 에 준 반복 횟수마다 추정을 저장해 돌려준다.
    """
    w = step or 1.0 / np.linalg.norm(A, 2) ** 2
    x = np.zeros(A.shape[1])
    snaps = {}
    for it in range(1, n_iter + 1):
        x = x + w * (A.T @ (y - A @ x))
        if it in record:
            snaps[it] = x.copy()
    return x, snaps


def landweber_svd(A, y, n_iter, step=None):
    """랜드웨버 k 회 반복의 닫힌 해. 필터 인자 f_i = 1 - (1 - w s_i^2)^k.

    반복을 돌리지 않고 SVD 로 같은 답을 내므로, landweber() 의 검증 경로가 된다.
    k 가 작을 때 f_i 는 작은 s_i 에서 0 에 가깝고, k 가 커질수록 전부 1 로 간다.
    """
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    w = step or 1.0 / s[0] ** 2
    f = 1 - (1 - w * s**2) ** n_iter
    return Vt.T @ (f * (U.T @ y) / s)


def tv_admm(A, y, lam, rho=None, n_iter=500):
    """1D TV 정규화를 ADMM 으로 푼다 (Boyd 등 2011, Wetzstein EE367 노트).

        min_x  1/2 || A x - y ||^2 + lam || z ||_1   조건 D x = z

    x 갱신은 선형계 하나, z 갱신은 부드러운 문턱(soft threshold), u 는 조건 위반의 누적이다.
    """
    n = A.shape[1]
    D = difference_matrix(n)
    rho = rho or 10 * lam
    M = np.linalg.inv(A.T @ A + rho * D.T @ D)
    Aty = A.T @ y
    z = np.zeros(n - 1)
    u = np.zeros(n - 1)
    for _ in range(n_iter):
        x = M @ (Aty + rho * D.T @ (z - u))
        v = D @ x + u
        z = np.sign(v) * np.maximum(np.abs(v) - lam / rho, 0)
        u = v - z
    return x


def tv_objective(A, y, x, lam):
    return 0.5 * np.sum((A @ x - y) ** 2) + lam * np.sum(np.abs(np.diff(x)))


def rel_error(est, ref):
    return np.linalg.norm(est - ref) / np.linalg.norm(ref)


def best_lambda(solver, A, y, x_true, lams):
    """정답을 알 때 오차가 가장 작은 lam. 방법끼리 공정하게 비교하려고 쓴다.

    실제 측정에는 정답이 없으므로 L곡선 같은 기준으로 골라야 한다. 여기서는 각 방법이
    낼 수 있는 가장 좋은 결과끼리 비교하려는 것이다.
    """
    errs = np.array([rel_error(solver(A, y, l), x_true) for l in lams])
    i = int(np.argmin(errs))
    return lams[i], errs[i], errs


def l_curve(A, y, lams, L=None):
    """lam 마다 (잔차 노름, 해의 반노름). 로그-로그로 그리면 L 자 모양이 된다고 알려져 있다.

    Hansen 은 이 모양이 특잇값이 서서히 0 으로 줄어드는 문제에서 나온다고 설명한다.
    상자 번짐처럼 조건수가 수백에 그치면 모서리가 뚜렷하지 않다.
    """
    L = np.eye(A.shape[1]) if L is None else L
    res, sol = [], []
    for l in lams:
        x = tikhonov(A, y, l, L)
        res.append(np.linalg.norm(A @ x - y))
        sol.append(np.linalg.norm(L @ x))
    return np.array(res), np.array(sol)


def discrepancy_lambda(A, y, sigma, lo=1e-6, hi=10.0):
    """불일치 원리(Morozov): 잔차가 잡음의 크기 sigma sqrt(m) 와 같아지는 lam.

    잡음이 섞인 측정을 잡음보다 더 정확하게 맞추려 하면 잡음까지 복원하게 된다.
    그래서 잔차가 딱 잡음만큼 남도록 lam 을 고른다. 잔차는 lam 에 대해 단조 증가하므로
    이분법으로 찾는다. 검출기의 잡음 수준을 알아야 쓸 수 있다.
    """
    from scipy.optimize import brentq
    target = sigma * np.sqrt(A.shape[0])
    f = lambda ll: np.linalg.norm(A @ tikhonov(A, y, 10**ll) - y) - target
    return 10 ** brentq(f, np.log10(lo), np.log10(hi))

"""실험 5-3 — "feature 공간을 넓힌다"가 무슨 뜻인지 가른다.

실험 5-2 에서 파장 구간을 450~750 → 450~1000 으로 넓혀도 corr(d, n1) 이
-0.9979 → -0.9967 로 거의 변하지 않았다. 수직입사 반사율은 위상 beta 를 통해
n1 과 d 의 곱에 주로 의존하므로, 같은 종류의 데이터를 더 모아도 그 축퇴 방향에
대한 정보가 늘지 않기 때문이다.

여기서는 **데이터 점 수를 고정한 채** 측정의 종류만 바꿔 비교한다.
총 300점으로 고정 → 단일각 300점 vs 2각 150점씩 vs 3각 100점씩.
"""
import warnings
warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

import numpy as np

from reflectance_model import N_SIO2
from angle_model import reflectance_angle

TRUE_D, NOISE_STD = 1490.0, 0.004
TOTAL_PTS = 300


def build_grid(angles, lo=450.0, hi=750.0, total=TOTAL_PTS):
    """각도마다 total/len(angles) 점씩 배분. 총 점수는 항상 total."""
    per = total // len(angles)
    return [(th, np.linspace(lo, hi, per)) for th in angles]


def model(p, grid):
    d, n1 = p
    return np.concatenate([reflectance_angle(d, wl, th, n1=n1) for th, wl in grid])


def measure(grid, seed=0, d=TRUE_D, n1=N_SIO2):
    y = model([d, n1], grid)
    return y + np.random.default_rng(seed).normal(0, NOISE_STD, y.size)


def jac(p, grid, h=(1e-3, 1e-6)):
    cols = []
    for i, hi in enumerate(h):
        pp = np.array(p, float); pp[i] += hi
        pm = np.array(p, float); pm[i] -= hi
        cols.append((model(pp, grid) - model(pm, grid)) / (2 * hi))
    return np.column_stack(cols)


def corr_cond(p, grid):
    A = jac(p, grid).T @ jac(p, grid)
    C = np.linalg.inv(A)
    return float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1])), float(np.linalg.cond(A))


def lm_multi(p0, grid, y, n_iter=80, tau=1e-3):
    """2파라미터 LM (Nielsen 댐핑) — 여러 각도를 이어붙인 잔차에 대해."""
    p = np.array(p0, float)
    r = model(p, grid) - y
    J = jac(p, grid); A = J.T @ J; g = J.T @ r
    mu = tau * float(np.max(np.diag(A))); nu = 2.0
    for _ in range(n_iter):
        if np.linalg.norm(g, np.inf) <= 1e-14:
            break
        try:
            h = np.linalg.solve(A + mu * np.eye(2), -g)
        except np.linalg.LinAlgError:
            mu *= nu; nu *= 2; continue
        pn = p + h
        # n1 이 물리적으로 말이 안 되는 값(<=1)으로 새면 모델이 NaN 을 낸다.
        # 스텝을 기각하고 댐핑을 키운다(비유한 잔차도 동일 처리).
        if not np.all(np.isfinite(pn)) or pn[1] <= 1.0:
            mu *= nu; nu *= 2; continue
        F0 = 0.5 * float(np.sum(r ** 2))
        rn = model(pn, grid) - y
        if not np.all(np.isfinite(rn)):
            mu *= nu; nu *= 2; continue
        F1 = 0.5 * float(np.sum(rn ** 2))
        den = 0.5 * float(h @ (mu * h - g))
        rho = (F0 - F1) / den if den > 0 else -1.0
        if rho > 0:
            p = pn; r = rn
            J = jac(p, grid); A = J.T @ J; g = J.T @ r
            mu *= max(1 / 3, 1 - (2 * rho - 1) ** 3); nu = 2.0
        else:
            mu *= nu; nu *= 2
    return p


CONFIGS = [
    ("수직입사만 450~750nm (300점)",        build_grid([0.0])),
    ("수직입사만 450~1000nm (300점)",       build_grid([0.0], hi=1000.0)),
    ("0°+60° 각 150점 (300점)",            build_grid([0.0, 60.0])),
    ("0°+45°+70° 각 100점 (300점)",        build_grid([0.0, 45.0, 70.0])),
]

if __name__ == "__main__":
    print("실험 5-3 — 총 데이터 점수 300 고정, 측정의 '종류'만 바꾼다")
    print(f"{'구성':34s}{'corr(d,n1)':>13s}{'cond(JtJ)':>12s}{'d 산포(nm)':>12s}{'n1 산포':>11s}")
    for name, grid in CONFIGS:
        c, k = corr_cond([TRUE_D, N_SIO2], grid)
        ds, ns = [], []
        for s in range(40):
            y = measure(grid, seed=300 + s)
            p = lm_multi([TRUE_D, N_SIO2], grid, y)
            ds.append(p[0]); ns.append(p[1])
        print(f"{name:34s}{c:13.6f}{k:12.3e}{np.std(ds, ddof=1):12.4f}{np.std(ns, ddof=1):11.5f}")

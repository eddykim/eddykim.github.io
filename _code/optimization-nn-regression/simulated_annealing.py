"""Simulated Annealing: 두께 d 위에서의 랜덤워크 + Metropolis 기준.

Kirkpatrick, Gelatt, Vecchi (1983)의 구조를 그대로 따른다 — 국소 최적화기를
쓰지 않고, 매 스텝 무작위로 제안한 d_prop을 Metropolis 기준
P(accept) = min(1, exp(-ΔJ/T))으로 채택/기각한 뒤 T를 기하급수적으로 냉각한다.
basin_hopping.py와 달리 로컬 옵티마이저 없이 순수하게 담금질 자체로 basin을 넘는다.
"""
import numpy as np
from levenberg_marquardt import objective


def simulated_annealing(d0, wavelength_nm, measured_R, n_iter=300, T0=0.02,
                         cooling=0.97, step_sigma=80.0, rng=None):
    """반환: d_hist, J_hist, T_hist, best_d, best_J.

    T0: 초기 온도. cooling: 스텝마다 T *= cooling으로 기하급수적 냉각.
    step_sigma: 제안 분포 N(0, step_sigma^2)의 표준편차 [nm].
    best_d/best_J는 현재 위치가 아니라 지금까지 방문한 것 중 최솟값을 추적한다
    (Metropolis는 확률적으로 나쁜 스텝도 받아들이므로 현재 위치가 최적이라는
    보장이 없다).
    """
    if rng is None:
        rng = np.random.default_rng()
    d = d0
    J = objective(d, wavelength_nm, measured_R)
    best_d, best_J = d, J
    d_hist = [d]
    J_hist = [J]
    T_hist = [T0]
    T = T0
    for _ in range(n_iter):
        d_prop = d + rng.normal(0, step_sigma)
        J_prop = objective(d_prop, wavelength_nm, measured_R)
        dJ = J_prop - J
        if dJ < 0 or rng.random() < np.exp(-dJ / T):
            d, J = d_prop, J_prop
        if J < best_J:
            best_d, best_J = d, J
        d_hist.append(d)
        J_hist.append(J)
        T *= cooling
        T_hist.append(T)
    return np.array(d_hist), np.array(J_hist), np.array(T_hist), best_d, best_J

"""Basin-Hopping: 로컬 최적화(LM) + 무작위 섭동 + Metropolis 기준.

Wales & Doye (1997)의 구조를 그대로 따른다. 매 hop마다
(1) 현재 해에 큰 무작위 섭동을 준 뒤 (2) levenberg_marquardt()로 그 지점에서
수렴할 때까지 로컬 최적화하고 (3) 두 국소 최솟값(현재 vs 새로 수렴한 곳) 사이를
Metropolis 기준으로 채택/기각한다. SA와 달리 담금질 대상이 원래 좌표 d가 아니라
"어느 basin의 바닥에 있는가"이므로 매 스텝이 항상 어떤 국소 최솟값 위에 있다.
"""
import numpy as np
from levenberg_marquardt import levenberg_marquardt


def basin_hopping(d0, wavelength_nm, measured_R, n_hops=15, T0=0.05,
                   cooling=1.0, perturb_sigma=150.0, lm_kwargs=None, rng=None):
    """반환: trace_d, trace_J, best_d, best_J, n_accept.

    perturb_sigma: 현재 해 주변 제안 분포 N(0, perturb_sigma^2)의 표준편차 [nm].
    basin 간격(~200nm, 1편 그림2)보다 커야 다른 basin으로 건너뛸 수 있다.
    cooling=1.0(기본값)은 온도를 고정한다 — basin-hopping 원 논문도 보통 고정
    온도를 쓰고, 냉각은 필요에 따라 추가하는 변형이다.
    """
    if rng is None:
        rng = np.random.default_rng()
    lm_kwargs = lm_kwargs or {}
    d_hist_lm, J_hist_lm, _ = levenberg_marquardt(d0, wavelength_nm, measured_R, **lm_kwargs)
    d, J = d_hist_lm[-1], J_hist_lm[-1]
    best_d, best_J = d, J
    trace_d = [d]
    trace_J = [J]
    T = T0
    n_accept = 0
    for _ in range(n_hops):
        d_trial0 = d + rng.normal(0, perturb_sigma)
        d_hist_trial, J_hist_trial, _ = levenberg_marquardt(d_trial0, wavelength_nm, measured_R, **lm_kwargs)
        d_trial, J_trial = d_hist_trial[-1], J_hist_trial[-1]
        dJ = J_trial - J
        if dJ < 0 or rng.random() < np.exp(-dJ / T):
            d, J = d_trial, J_trial
            n_accept += 1
        if J < best_J:
            best_d, best_J = d, J
        trace_d.append(d)
        trace_J.append(J)
        T *= cooling
    return np.array(trace_d), np.array(trace_J), best_d, best_J, n_accept

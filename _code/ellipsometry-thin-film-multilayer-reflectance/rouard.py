"""로아드 방법(Rouard's method)에 의한 다층 박막 반사계수 계산.

가장 아래(기판에 가장 가까운) 계면에서 시작해, "그 아래에 있는 모든 층"을 하나의
등가 매질로 접어 올리는 재귀다. 재귀의 각 단계는 단층 박막의 유효 반사계수 식과
형태가 완전히 같다.

    rho_{M+1} = r_{M,M+1}                                    (기판 계면에서 출발)
    rho_m     = (r_{m-1,m} + rho_{m+1} e^{-2j beta_m})
                / (1 + r_{m-1,m} rho_{m+1} e^{-2j beta_m})
    r_total   = rho_1

핵심은 재귀의 분자·분모에 들어가는 것이 바로 아래 계면의 "맨" 반사계수
r_{m,m+1} 이 아니라, 그 아래 전체를 이미 접어 올린 rho_{m+1} 이라는 점이다.
r_{m,m+1} 을 그대로 넣으면 재귀가 성립하지 않고 한 층 아래까지만 반영된다.
또 r_total 은 rho_1 그 자체이며, 여기에 r_{01} 과 beta_1 을 한 번 더 적용하면
같은 계면과 같은 층을 두 번 세게 된다.

이 구현이 TMM/SMM 과 같은 값을 주는지는 verify_methods.py 에서 확인한다.
"""
import numpy as np


def _cos_theta(n0, theta0, n):
    """복소굴절률 n 에서의 굴절각 cosine. 스넬의 법칙 N0 sin(theta0) = n sin(theta)."""
    sin_t = (n0 / n) * np.sin(theta0)
    return np.sqrt(1 - sin_t ** 2 + 0j)


def _fresnel_r(n_a, cos_a, n_b, cos_b, pol):
    """계면 반사계수 r_{ab}. 1편·smm_tensor.py 와 같은 Fresnel 부호 규약을 쓴다."""
    if pol == "s":
        return (n_a * cos_a - n_b * cos_b) / (n_a * cos_a + n_b * cos_b)
    return (n_b * cos_a - n_a * cos_b) / (n_b * cos_a + n_a * cos_b)


def rouard_r(n_list, d_list, wavelength, theta0=0.0, pol="s"):
    """다층 박막의 총반사계수 (단일 파장, 단일 입사각).

    n_list: [N0(ambient), N1, ..., NM(박막들), Nsub] 길이 M+2
    d_list: [d1, ..., dM] 각 박막 두께, 길이 M
    theta0: 입사각 [rad]
    """
    n0 = n_list[0]
    cos_list = [_cos_theta(n0, theta0, n) for n in n_list]
    M = len(d_list)

    # 맨 아래 계면(마지막 박막/기판)에서 출발한다.
    rho = _fresnel_r(n_list[M], cos_list[M], n_list[M + 1], cos_list[M + 1], pol)

    # 위로 한 층씩 접어 올린다.
    for m in range(M, 0, -1):
        beta = 2 * np.pi / wavelength * n_list[m] * d_list[m - 1] * cos_list[m]
        r_top = _fresnel_r(n_list[m - 1], cos_list[m - 1], n_list[m], cos_list[m], pol)
        phase = np.exp(-2j * beta)
        rho = (r_top + rho * phase) / (1 + r_top * rho * phase)

    return rho

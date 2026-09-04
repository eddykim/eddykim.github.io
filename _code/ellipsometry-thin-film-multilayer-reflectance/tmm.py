"""전달행렬법(TMM, Abeles characteristic matrix) 다층 박막 반사·투과계수 계산.

식 2.30~2.37 (김영준, 서울대 박사논문, 2.5절)의 특성행렬 Q_m을 그대로 구현한다.
편광에 따른 광학 어드미턴스(optical admittance) eta = n*cos(theta) (s-편광)
또는 n/cos(theta) (p-편광)를 도입하면 같은 행렬식으로 두 편광을 모두 다룰 수 있다.
"""
import numpy as np


def _cos_theta(n0, theta0, n):
    """복소굴절률 n에서의 굴절각 cosine. 스넬의 법칙 N0 sin(theta0) = n sin(theta)."""
    sin_t = (n0 / n) * np.sin(theta0)
    return np.sqrt(1 - sin_t ** 2 + 0j)


def _admittance(n, cos_t, pol):
    return n * cos_t if pol == "s" else n / cos_t


def tmm_rt(n_list, d_list, wavelength, theta0=0.0, pol="s"):
    """다층 박막 TMM 반사·투과계수 (단일 파장, 단일 입사각).

    n_list: [N0(ambient), N1, ..., Nm(박막들), Nt(substrate)] 길이 m+2
    d_list: [d1, ..., dm] 박막 두께(각 층), 길이 m. ambient/substrate는 두께 없음.
    wavelength, theta0: 스칼라. theta0 단위는 rad.
    """
    n0 = n_list[0]
    cos_list = [_cos_theta(n0, theta0, n) for n in n_list]
    eta = [_admittance(n, c, pol) for n, c in zip(n_list, cos_list)]

    D0 = np.array([[1, 1], [eta[0], -eta[0]]], dtype=complex)
    Dsub = np.array([[1, 1], [eta[-1], -eta[-1]]], dtype=complex)

    M = np.linalg.inv(D0)
    for n, d, c, e in zip(n_list[1:-1], d_list, cos_list[1:-1], eta[1:-1]):
        beta = 2 * np.pi / wavelength * n * d * c
        Q = np.array([
            [np.cos(beta), 1j / e * np.sin(beta)],
            [1j * e * np.sin(beta), np.cos(beta)],
        ])
        M = M @ Q
    M = M @ Dsub

    r_total = M[1, 0] / M[0, 0]
    t_total = 1.0 / M[0, 0]
    return r_total, t_total


def reflectance_spectrum(n_list, d_list, wavelengths, theta0=0.0, pol="s"):
    """파장 배열에 대해 반사율 R = |r_total|^2 을 계산."""
    r = np.array([tmm_rt(n_list, d_list, wl, theta0, pol)[0] for wl in wavelengths])
    return np.abs(r) ** 2

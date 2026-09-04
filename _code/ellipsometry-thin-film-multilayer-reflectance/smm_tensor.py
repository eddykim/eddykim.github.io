"""산란행렬법(SMM)을 입사각×파장 텐서로 벡터화한 구현.

tmm.py의 tmm_rt는 (파장, 입사각) 한 점씩 호출해야 하는 스칼라 버전이다.
입사각 분해(angle-resolved) 계측처럼 파장 축과 입사각 축을 동시에 반복 계산해야
하는 경우, 계면행렬(I)·층행렬(L)을 넘파이 배열 전체에 대해 한 번에 만들고
numpy의 배치 행렬곱(@)으로 층 수만큼만 반복하면 파이썬 반복문 없이 계산된다.
"""
import numpy as np


def _cos_theta_grid(n_list, theta0):
    """theta0: (n_angle,) 입사각 배열. 반환: 매질별 (n_angle, 1) 코사인 리스트."""
    n0 = n_list[0]
    sin0 = n0 * np.sin(theta0)[:, None]  # (n_angle, 1)
    return [np.sqrt(1 - (sin0 / n) ** 2 + 0j) for n in n_list]


def _interface_matrix(n_a, n_b, cos_a, cos_b, pol):
    """두 매질 경계의 계면행렬 I (식 2.38~2.39). shape: (n_angle, n_wav, 2, 2)."""
    if pol == "s":
        r = (n_a * cos_a - n_b * cos_b) / (n_a * cos_a + n_b * cos_b)
        t = 2 * n_a * cos_a / (n_a * cos_a + n_b * cos_b)
    else:
        r = (n_b * cos_a - n_a * cos_b) / (n_b * cos_a + n_a * cos_b)
        t = 2 * n_a * cos_a / (n_b * cos_a + n_a * cos_b)
    r, t = np.broadcast_arrays(r, t)
    I = np.empty((*r.shape, 2, 2), dtype=complex)
    I[..., 0, 0], I[..., 0, 1] = 1 / t, r / t
    I[..., 1, 0], I[..., 1, 1] = r / t, 1 / t
    return I


def _layer_matrix(n, d, cos_n, wavelength):
    """박막 내부 위상지연을 담는 층행렬 L (식 2.40). shape: (n_angle, n_wav, 2, 2)."""
    beta = 2 * np.pi / wavelength[None, :] * n * d * cos_n
    L = np.zeros((*beta.shape, 2, 2), dtype=complex)
    L[..., 0, 0] = np.exp(1j * beta)
    L[..., 1, 1] = np.exp(-1j * beta)
    return L


def smm_reflectance_tensor(n_list, d_list, theta0, wavelength, pol="s"):
    """다층 박막 SMM 반사계수를 (입사각, 파장) 텐서로 한 번에 계산.

    n_list: [N0, N1, ..., Nm, Nt] (상수 굴절률 리스트, 비분산 근사)
    d_list: [d1, ..., dm] 박막 두께
    theta0: (n_angle,) 입사각 배열(rad), wavelength: (n_wav,) 파장 배열
    반환: r_total, shape (n_angle, n_wav)
    """
    cos_list = _cos_theta_grid(n_list, theta0)  # 각 원소 shape (n_angle, 1), 브로드캐스트로 (n_angle, n_wav)와 결합됨
    S = _interface_matrix(n_list[0], n_list[1], cos_list[0], cos_list[1], pol)
    for j in range(1, len(d_list) + 1):
        L = _layer_matrix(n_list[j], d_list[j - 1], cos_list[j], wavelength)
        I = _interface_matrix(n_list[j], n_list[j + 1], cos_list[j], cos_list[j + 1], pol)
        S = S @ L @ I
    return S[..., 1, 0] / S[..., 0, 0]


def sample_mueller_matrix(r_p, r_s):
    """이상적(편광 소멸 없는) 시료의 뮬러 행렬. r_p, r_s: 동일 shape 텐서.

    2편 식(뮬러 행렬 형식)과 논문 식 2.44를 그대로 따른다.
    """
    rho = r_p / r_s
    psi = np.arctan(np.abs(rho))
    delta = np.angle(rho)
    m00 = 0.5 * (np.abs(r_p) ** 2 + np.abs(r_s) ** 2)
    cos2psi, sin2psi = np.cos(2 * psi), np.sin(2 * psi)
    M = np.zeros((*rho.shape, 4, 4))
    M[..., 0, 0] = M[..., 1, 1] = 1.0
    M[..., 0, 1] = M[..., 1, 0] = -cos2psi
    M[..., 2, 2] = M[..., 3, 3] = sin2psi * np.cos(delta)
    M[..., 2, 3] = sin2psi * np.sin(delta)
    M[..., 3, 2] = -sin2psi * np.sin(delta)
    return m00[..., None, None] * M, psi, delta

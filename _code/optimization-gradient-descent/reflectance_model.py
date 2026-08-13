"""단일 SiO2 박막 (Air / SiO2 / Si) 수직입사 반사율 모델.

repo_sr/src/Reflectance_models.py 의 다층 transfer-matrix 코드를
단일 박막·수직입사 특수해(Airy 공식)로 단순화해서 다시 작성했다.
SiO2/Si 굴절률은 400~800nm 구간의 대표값으로 근사한 상수다
(파장 의존성을 정확히 반영한 실측 nk 값이 아니라 시연용 근사치).
"""
import numpy as np

N_AIR = 1.0
N_SIO2 = 1.46
N_SI = 3.9 - 0.02j


def reflectance(thickness_nm, wavelength_nm, n1=N_SIO2, n2=N_SI):
    """Air/SiO2/Si 3층 구조의 수직입사 반사율 R(lambda).

    thickness_nm: SiO2 두께 [nm], 스칼라
    wavelength_nm: 파장 배열 [nm]
    """
    n0 = N_AIR

    r01 = (n0 - n1) / (n0 + n1)
    r12 = (n1 - n2) / (n1 + n2)

    beta = 2 * np.pi * n1 * thickness_nm / wavelength_nm
    phase = np.exp(-2j * beta)

    r = (r01 + r12 * phase) / (1 + r01 * r12 * phase)
    return np.abs(r) ** 2


if __name__ == "__main__":
    wavelength_nm = np.linspace(450, 750, 5)
    print(reflectance(150.0, wavelength_nm))

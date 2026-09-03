"""단일 SiO2 박막 (Air / SiO2 / Si) 수직입사 반사율 모델.

1편(_code/optimization-gradient-descent/reflectance_model.py)과 동일하다.
포스트 자산을 slug 폴더에 자기완결적으로 모아두기 위해 복사해 둔다.
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

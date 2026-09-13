"""경사입사 반사율 — 실험 5의 "다른 종류의 feature" 용.

reflectance_model.py 는 수직입사 전용이다. 수직입사에서 반사율은 위상 항
beta = 2*pi*n1*d/lambda 를 통해 n1 과 d 의 **곱**에 주로 의존하므로, 파장을 더
모아도 (d, n1) 축퇴가 풀리지 않는다(실험 5-2에서 확인됨).

경사입사는 beta 에 cos(theta1) 이 붙고 Fresnel 계수가 편광별로 갈라지므로
n1 에 대한 의존이 d 와 다른 방식으로 들어온다. Kwak & Kim 2023 이 말한
"다각도 측정이 유일성을 개선한다"를 우리 예제에서 확인하기 위한 모델이다.

원본 reflectance_model.py 는 손대지 않았다.
"""
import numpy as np

from reflectance_model import N_AIR, N_SIO2, N_SI


def reflectance_angle(thickness_nm, wavelength_nm, theta0_deg,
                      n1=N_SIO2, n2=N_SI, pol="unpol"):
    """Air/SiO2/Si 경사입사 반사율.

    theta0_deg: 공기 중 입사각 [deg]. 0 이면 수직입사 모델과 일치한다.
    pol: "s", "p", "unpol"(s,p 평균)
    """
    n0 = N_AIR
    th0 = np.deg2rad(theta0_deg)
    sin0 = n0 * np.sin(th0)

    # 복소 굴절률에서도 성립하도록 복소 제곱근을 쓴다
    c0 = np.cos(th0)
    c1 = np.sqrt(1.0 - (sin0 / n1) ** 2 + 0j)
    c2 = np.sqrt(1.0 - (sin0 / n2) ** 2 + 0j)

    # Fresnel 계수
    r01_s = (n0 * c0 - n1 * c1) / (n0 * c0 + n1 * c1)
    r12_s = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
    r01_p = (n1 * c0 - n0 * c1) / (n1 * c0 + n0 * c1)
    r12_p = (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2)

    beta = 2 * np.pi * n1 * thickness_nm * c1 / wavelength_nm
    ph = np.exp(-2j * beta)

    rs = (r01_s + r12_s * ph) / (1 + r01_s * r12_s * ph)
    rp = (r01_p + r12_p * ph) / (1 + r01_p * r12_p * ph)

    Rs, Rp = np.abs(rs) ** 2, np.abs(rp) ** 2
    if pol == "s":
        return Rs
    if pol == "p":
        return Rp
    return 0.5 * (Rs + Rp)

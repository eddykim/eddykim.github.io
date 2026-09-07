"""로아드·TMM·SMM 세 계산법이 같은 반사계수를 내놓는지 검증한다.

본문 6절의 "이상적인 시료에 대해서는 세 방법 모두 같은 반사·투과계수를 내놓는다"
는 서술의 근거다. 단층과 3층 구조를, s/p 편광과 수직·경사 입사에서 모두 비교한다.

실행: python verify_methods.py
"""
import numpy as np

from rouard import rouard_r
from smm_tensor import smm_reflectance_tensor
from tmm import tmm_rt

N_AIR, N_SI = 1.0, 3.88 - 0.02j

CASES = [
    ("단층  SiO2 300nm", [N_AIR, 1.46, N_SI], [300.0]),
    ("3층   SiO2/TiO2/SiN", [N_AIR, 1.46, 2.35, 2.02, N_SI], [120.0, 80.0, 45.0]),
]
WAVELENGTH = 633.0

print(f"파장 {WAVELENGTH:.0f}nm, 세 방법의 총반사계수 r 비교")
print()
worst = 0.0
for name, n_list, d_list in CASES:
    for theta_deg in (0.0, 65.0):
        for pol in ("s", "p"):
            th = np.radians(theta_deg)
            r_tmm = tmm_rt(n_list, d_list, WAVELENGTH, th, pol)[0]
            r_rou = rouard_r(n_list, d_list, WAVELENGTH, th, pol)
            r_smm = smm_reflectance_tensor(
                n_list, d_list, np.array([th]), np.array([WAVELENGTH]), pol
            ).ravel()[0]
            gap = max(abs(r_rou - r_tmm), abs(r_smm - r_tmm))
            worst = max(worst, gap)
            print(f"  {name}  theta={theta_deg:4.0f}deg  {pol}-pol")
            print(f"    Rouard {r_rou.real:+.8f}{r_rou.imag:+.8f}j")
            print(f"    TMM    {r_tmm.real:+.8f}{r_tmm.imag:+.8f}j")
            print(f"    SMM    {r_smm.real:+.8f}{r_smm.imag:+.8f}j   최대 차이 {gap:.2e}")
print()
print(f"세 방법의 최대 불일치: {worst:.3e}")

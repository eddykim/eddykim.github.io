"""1편 그림 4개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-gradient-descent/ 에 fig1~4 저장

측정치는 실제 spectrometer 데이터가 아니라, 실제 두께를 알고 있는 SiO2/Si
샘플(과거 Metrol_LAB_SiO2 측정 시리즈에 실존하는 1490nm 두께)을 가정하고
물리 모델에 가우시안 노이즈를 더해 합성한 것이다. 150nm 근방은 박막이
광학적으로 얇아 반사율이 두께 변화에 둔감해(파라미터 모호성이 큼) 예시로
적합하지 않아, 간섭 무늬가 뚜렷하게 나타나는 1490nm 샘플로 바꿨다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from gradient_descent import gradient_descent, objective
from reflectance_model import reflectance

# 한글 텍스트(제목/축이름/범례)에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면
# log축 눈금의 마이너스 기호(-1, -2 지수)가 한글 폰트에 없어 깨지기 때문에,
# 숫자 눈금은 기본 폰트(DejaVu Sans)를 그대로 쓴다.
# macOS 전용 설정이다. 다른 OS는 나눔고딕 등 설치된 한글 폰트 이름으로 교체할 것.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}  # legend(prop=...)는 FontProperties 인자 이름(family)을 쓴다

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-gradient-descent",
)
os.makedirs(OUT_DIR, exist_ok=True)

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0  # nm
noise_std = 0.004
initial_guess = 1540.0

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# 그림1: 모델 초기값 vs 측정치
R_initial = reflectance(initial_guess, wavelength_nm)
plt.figure(figsize=(6, 4))
plt.plot(wavelength_nm, R_measured, ".", color="gray", markersize=3, label="측정치 (합성)")
plt.plot(wavelength_nm, R_initial, "-", color="C0", label=f"초기 모델 (d={initial_guess:.0f} nm)")
plt.xlabel("파장 (nm)", **KFONT)
plt.ylabel("반사율", **KFONT)
plt.legend(prop=LFONT)
plt.title("그림1. 모델 초기값 vs 측정치", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig1-model-vs-measurement.png"), dpi=150)
plt.close()

# 그림2: 두께에 따른 목적함수 J(d) - 오차 지형(landscape)
d_scan = np.linspace(1150, 1850, 400)
J_scan = np.array([objective(d, wavelength_nm, R_measured) for d in d_scan])
plt.figure(figsize=(6, 4))
plt.plot(d_scan, J_scan, "-", color="C0")
plt.axvline(true_thickness, color="gray", linestyle="--", label="실제 두께 (global minimum)")
plt.xlabel("두께 d (nm)", **KFONT)
plt.ylabel("목적함수 J(d)", **KFONT)
plt.legend(prop=LFONT)
plt.title("그림2. 두께에 따른 목적함수 J(d)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig2-objective-landscape.png"), dpi=150)
plt.close()

# alpha=1500(너무 큼)과 alpha=300(적절함) 두 경우를 같은 조건(n_iter=20)으로 돌려서
# 같은 두 종류의 그림(두께 vs iteration, J vs iteration)으로 나란히 비교한다.
N_ITER_COMPARE = 20
d_hist_bad, J_hist_bad = gradient_descent(
    initial_guess, wavelength_nm, R_measured, alpha=1500, n_iter=N_ITER_COMPARE
)
d_hist_good, J_hist_good = gradient_descent(
    initial_guess, wavelength_nm, R_measured, alpha=300, n_iter=N_ITER_COMPARE
)

# 그림3: 두께 추정값 vs iteration - step size 비교
plt.figure(figsize=(6, 4))
plt.plot(d_hist_bad, "o-", color="C3", markersize=4, label="alpha=1500 (너무 큼)")
plt.plot(d_hist_good, "o-", color="C0", markersize=4, label="alpha=300 (적절함)")
plt.axhline(true_thickness, color="gray", linestyle="--", label="실제 두께")
plt.xlabel("iteration")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT)
plt.title("그림3. 두께 추정값 vs iteration - step size 비교", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig3-thickness-vs-iteration.png"), dpi=150)
plt.close()

# 그림4: 목적함수 J vs iteration - step size 비교
plt.figure(figsize=(6, 4))
plt.plot(J_hist_bad, "o-", color="C3", markersize=4, label="alpha=1500 (너무 큼)")
plt.plot(J_hist_good, "o-", color="C0", markersize=4, label="alpha=300 (적절함)")
plt.xlabel("iteration")
plt.ylabel("목적함수 J (error)", **KFONT)
plt.yscale("log")
plt.legend(prop=LFONT)
plt.title("그림4. 목적함수 J vs iteration - step size 비교", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4-objective-vs-iteration.png"), dpi=150)
plt.close()

print(f"[alpha=300]  최종 추정 두께: {d_hist_good[-1]:.2f} nm (실제: {true_thickness} nm)")
print(f"[alpha=1500] 최종 추정 두께: {d_hist_bad[-1]:.2f} nm (발산)")

"""1편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-gradient-descent/      (한국어)
      ../../assets/img/posts/optimization-gradient-descent/en/   (영문)

같은 seed 를 쓰므로 두 언어의 데이터는 완전히 동일하고 라벨만 다르다.

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

# 영문판은 한글 폰트가 필요 없으므로 matplotlib 기본 폰트를 그대로 쓴다.
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-gradient-descent",
)

# 언어별 라벨. 데이터 생성 코드는 공유하고 문자열만 갈아 끼운다.
LABELS = {
    "ko": {
        "dir": BASE_DIR,
        "font": KFONT, "legend": LFONT,
        "measured": "측정치 (합성)", "initial": "초기 모델 (d={:.0f} nm)",
        "wavelength": "파장 (nm)", "reflectance": "반사율",
        "fig1": "그림1. 모델 초기값 vs 측정치",
        "true_thickness_gm": "실제 두께 (global minimum)",
        "thickness": "두께 d (nm)", "objective": "목적함수 J(d)",
        "fig2": "그림2. 두께에 따른 목적함수 J(d)",
        "too_large": "alpha=1500 (너무 큼)", "just_right": "alpha=300 (적절함)",
        "true_thickness": "실제 두께", "thickness_est": "두께 추정값 (nm)",
        "fig3": "그림3. 두께 추정값 vs iteration - step size 비교",
        "objective_err": "목적함수 J (error)",
        "fig4": "그림4. 목적함수 J vs iteration - step size 비교",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"),
        "font": EFONT, "legend": ELFONT,
        "measured": "Measured (synthetic)", "initial": "Initial model (d={:.0f} nm)",
        "wavelength": "Wavelength (nm)", "reflectance": "Reflectance",
        "fig1": "Fig 1. Initial model vs measurement",
        "true_thickness_gm": "True thickness (global minimum)",
        "thickness": "Thickness d (nm)", "objective": "Objective J(d)",
        "fig2": "Fig 2. Objective J(d) versus thickness",
        "too_large": "alpha=1500 (too large)", "just_right": "alpha=300 (appropriate)",
        "true_thickness": "True thickness", "thickness_est": "Estimated thickness (nm)",
        "fig3": "Fig 3. Thickness estimate vs iteration - step size",
        "objective_err": "Objective J (error)",
        "fig4": "Fig 4. Objective J vs iteration - step size",
    },
}

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0  # nm
noise_std = 0.004
initial_guess = 1540.0

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)
R_initial = reflectance(initial_guess, wavelength_nm)

# 목적함수 지형 스캔
d_scan = np.linspace(1150, 1850, 400)
J_scan = np.array([objective(d, wavelength_nm, R_measured) for d in d_scan])

# alpha=1500(너무 큼)과 alpha=300(적절함)을 같은 조건(n_iter=20)으로 돌려 나란히 비교한다.
N_ITER_COMPARE = 20
d_hist_bad, J_hist_bad = gradient_descent(
    initial_guess, wavelength_nm, R_measured, alpha=1500, n_iter=N_ITER_COMPARE
)
d_hist_good, J_hist_good = gradient_descent(
    initial_guess, wavelength_nm, R_measured, alpha=300, n_iter=N_ITER_COMPARE
)


def render(L):
    """주어진 라벨 묶음으로 그림 4개를 그린다. 데이터는 위에서 한 번만 계산한다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림1: 모델 초기값 vs 측정치
    plt.figure(figsize=(6, 4))
    plt.plot(wavelength_nm, R_measured, ".", color="gray", markersize=3, label=L["measured"])
    plt.plot(wavelength_nm, R_initial, "-", color="C0", label=L["initial"].format(initial_guess))
    plt.xlabel(L["wavelength"], **F)
    plt.ylabel(L["reflectance"], **F)
    plt.legend(prop=LF)
    plt.title(L["fig1"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig1-model-vs-measurement.png"), dpi=150)
    plt.close()

    # 그림2: 두께에 따른 목적함수 J(d) - 오차 지형(landscape)
    plt.figure(figsize=(6, 4))
    plt.plot(d_scan, J_scan, "-", color="C0")
    plt.axvline(true_thickness, color="gray", linestyle="--", label=L["true_thickness_gm"])
    plt.xlabel(L["thickness"], **F)
    plt.ylabel(L["objective"], **F)
    plt.legend(prop=LF)
    plt.title(L["fig2"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig2-objective-landscape.png"), dpi=150)
    plt.close()

    # 그림3: 두께 추정값 vs iteration
    plt.figure(figsize=(6, 4))
    plt.plot(d_hist_bad, "o-", color="C3", markersize=4, label=L["too_large"])
    plt.plot(d_hist_good, "o-", color="C0", markersize=4, label=L["just_right"])
    plt.axhline(true_thickness, color="gray", linestyle="--", label=L["true_thickness"])
    plt.xlabel("iteration")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF)
    plt.title(L["fig3"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig3-thickness-vs-iteration.png"), dpi=150)
    plt.close()

    # 그림4: 목적함수 J vs iteration
    plt.figure(figsize=(6, 4))
    plt.plot(J_hist_bad, "o-", color="C3", markersize=4, label=L["too_large"])
    plt.plot(J_hist_good, "o-", color="C0", markersize=4, label=L["just_right"])
    plt.xlabel("iteration")
    plt.ylabel(L["objective_err"], **F)
    plt.yscale("log")
    plt.legend(prop=LF)
    plt.title(L["fig4"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig4-objective-vs-iteration.png"), dpi=150)
    plt.close()


for lang, labels in LABELS.items():
    render(labels)
    print(f"[{lang}] 그림 4개 저장: {os.path.normpath(labels['dir'])}")

print(f"[alpha=300]  최종 추정 두께: {d_hist_good[-1]:.2f} nm (실제: {true_thickness} nm)")
print(f"[alpha=1500] 최종 추정 두께: {d_hist_bad[-1]:.2f} nm (발산)")

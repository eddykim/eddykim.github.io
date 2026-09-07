"""최적화 3편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-levenberg-marquardt/      (한국어)
      ../../assets/img/posts/optimization-levenberg-marquardt/en/   (영문)

측정 데이터는 1·2편과 동일(SiO2/Si 단층, 두께 1490nm, seed=0 노이즈).
그림4의 "실패 사례"만 예외적으로, 피팅 모델의 SiO2 굴절률을 실제(1.46)와
다르게(1.02) 잘못 가정해 간섭 콘트라스트가 거의 없는 상황을 인위적으로 만든다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from gradient_descent import gradient_descent
from newton import newton
from gauss_newton import gauss_newton
from levenberg_marquardt import levenberg_marquardt, objective as lm_objective
from reflectance_model import reflectance

# 한글 텍스트에만 한글 폰트를 지정한다. 시리즈의 다른 편과 같은 AppleGothic 을 쓴다.
# macOS 전용 설정이다. 다른 OS는 나눔고딕 등 설치된 한글 폰트 이름으로 교체할 것.
# 그리스 문자(tau, mu)는 폰트에 글리프가 없어 그대로 넣으면 두부(tofu)가 되므로
# 반드시 $\tau$ 처럼 mathtext 로 쓴다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-levenberg-marquardt",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "mu_small": r"$\mu_0$ 작게 ($\tau$=1e-6, GN에 가까움)",
        "mu_large": r"$\mu_0$ 크게 ($\tau$=1e6, steepest descent에 가까움)",
        "mu_small_short": r"$\mu_0$ 작게 ($\tau$=1e-6)",
        "mu_large_short": r"$\mu_0$ 크게 ($\tau$=1e6)",
        "true_thickness": "실제 두께", "thickness_est": "두께 추정값 (nm)",
        "fig1": r"그림1. $\mu_0$ 크기에 따른 LM 수렴 궤적 (d0=1540nm)",
        "damping": r"damping parameter $\mu$",
        "fig2": r"그림2. gain ratio 기반 $\mu$ 자동조정 추이",
        "gd": "Gradient Descent (alpha=300)", "newton": "Newton (발산)",
        "gn": "Gauss-Newton", "lm": r"Levenberg-Marquardt ($\tau$=1e-6)",
        "fig3": "그림3. 네 방법 비교 (d0=1540nm)",
        "gn_cycle": "Gauss-Newton (영원히 진동)", "lm_plain": "Levenberg-Marquardt",
        "wrong_floor": "이 모델의 실제 최솟값",
        "fig4": "그림4. 굴절률을 잘못 가정한 모델에서: GN 실패 vs LM 성공",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "mu_small": r"small $\mu_0$ ($\tau$=1e-6, close to GN)",
        "mu_large": r"large $\mu_0$ ($\tau$=1e6, close to steepest descent)",
        "mu_small_short": r"small $\mu_0$ ($\tau$=1e-6)",
        "mu_large_short": r"large $\mu_0$ ($\tau$=1e6)",
        "true_thickness": "true thickness", "thickness_est": "Estimated thickness (nm)",
        "fig1": r"Fig 1. LM trajectories by size of $\mu_0$ (d0=1540nm)",
        "damping": r"damping parameter $\mu$",
        "fig2": r"Fig 2. Adaptation of $\mu$ driven by the gain ratio",
        "gd": "Gradient Descent (alpha=300)", "newton": "Newton (diverging)",
        "gn": "Gauss-Newton", "lm": r"Levenberg-Marquardt ($\tau$=1e-6)",
        "fig3": "Fig 3. The four methods compared (d0=1540nm)",
        "gn_cycle": "Gauss-Newton (cycles forever)", "lm_plain": "Levenberg-Marquardt",
        "wrong_floor": "true minimum of this model",
        "fig4": "Fig 4. With a wrong refractive index: GN fails, LM rescues",
    },
}

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0
noise_std = 0.004
initial_guess = 1540.0

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

# 그림1·2: mu0을 작게(tau=1e-6, GN에 가까움) vs 크게(tau=1e6, steepest descent에
# 가까움) 잡았을 때 LM이 어떻게 다르게 움직이는지 비교.
N_ITER_MU = 20
d_small, J_small, mu_small = levenberg_marquardt(
    initial_guess, wavelength_nm, R_measured, tau=1e-6, n_iter=N_ITER_MU)
d_large, J_large, mu_large = levenberg_marquardt(
    initial_guess, wavelength_nm, R_measured, tau=1e6, n_iter=N_ITER_MU)

# 그림3: GD / Newton / Gauss-Newton / LM(tau=1e-6) 네 방법 비교, d0=1540, 10스텝.
N_ITER = 10
d_gd, J_gd = gradient_descent(initial_guess, wavelength_nm, R_measured, alpha=300, n_iter=N_ITER)
d_nt, J_nt = newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)
d_gn, J_gn = gauss_newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)
d_lm, J_lm, _ = levenberg_marquardt(initial_guess, wavelength_nm, R_measured, tau=1e-6, n_iter=N_ITER)

# 그림4: 실패 사례 실증. 피팅 모델의 SiO2 굴절률을 1.02로 잘못 가정하면
# 간섭 콘트라스트가 거의 사라져 Jacobian이 거의 특이(near-singular)해진다.
# 이 상황에서 Gauss-Newton은 세 값 사이를 영원히 순환하며 수렴하지 않는다.
WRONG_N1 = 1.02


def residual_wrong(d):
    return reflectance(d, wavelength_nm, n1=WRONG_N1) - R_measured


def objective_wrong(d):
    r = residual_wrong(d)
    return 0.5 * np.sum(r ** 2)


def jacobian_wrong(d, h=1e-3):
    return (residual_wrong(d + h) - residual_wrong(d - h)) / (2 * h)


def gauss_newton_wrong(d0, n_iter):
    d = d0
    d_hist = [d]
    for _ in range(n_iter):
        r = residual_wrong(d)
        Jr = jacobian_wrong(d)
        step = -np.sum(Jr * r) / np.sum(Jr ** 2)
        d = d + step
        d_hist.append(d)
    return np.array(d_hist)


def levenberg_marquardt_wrong(d0, n_iter, tau=1e-3, eps1=1e-12, eps2=1e-12):
    d = d0
    r = residual_wrong(d)
    Jr = jacobian_wrong(d)
    A = np.sum(Jr ** 2)
    g = np.sum(Jr * r)
    mu = tau * A
    nu = 2.0
    d_hist = [d]
    found = abs(g) <= eps1
    for _ in range(n_iter):
        if found:
            break
        h_lm = -g / (A + mu)
        d_new = d + h_lm
        F_old = objective_wrong(d)
        F_new = objective_wrong(d_new)
        L0_minus_Lh = 0.5 * h_lm * (mu * h_lm - g)
        rho = (F_old - F_new) / L0_minus_Lh if L0_minus_Lh > 0 else -1.0
        if rho > 0:
            d = d_new
            r = residual_wrong(d)
            Jr = jacobian_wrong(d)
            A = np.sum(Jr ** 2)
            g = np.sum(Jr * r)
            found = (abs(g) <= eps1) or (abs(h_lm) <= eps2 * (abs(d) + eps2))
            mu = mu * max(1.0 / 3.0, 1.0 - (2 * rho - 1) ** 3)
            nu = 2.0
        else:
            mu = mu * nu
            nu = 2.0 * nu
        d_hist.append(d)
    return np.array(d_hist)


N_ITER_FAIL = 20
d_gn_wrong = gauss_newton_wrong(initial_guess, N_ITER_FAIL)
d_lm_wrong = levenberg_marquardt_wrong(initial_guess, N_ITER_FAIL, tau=1e-3)

# 이 (잘못된) 모델의 실제 최솟값 위치를 촘촘한 스캔으로 찾는다.
_ds = np.linspace(1400, 1700, 3000)
_Js = np.array([objective_wrong(d) for d in _ds])
wrong_floor_d = float(_ds[np.argmin(_Js)])
wrong_floor_J = float(_Js.min())


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L):
    """주어진 라벨 묶음으로 그림 4개를 그린다. 계산은 위에서 한 번만 했다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림1: mu0 크기에 따른 수렴 궤적
    plt.figure(figsize=(6, 4))
    plt.plot(d_small, "o-", color="C0", markersize=4, label=L["mu_small"])
    plt.plot(d_large, "o-", color="C3", markersize=4, label=L["mu_large"])
    plt.axhline(true_thickness, color="gray", linestyle="--", label=L["true_thickness"])
    plt.xlabel("iteration")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig1"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig1-mu0-comparison.png"), dpi=150)
    plt.close()

    # 그림2: gain ratio 기반 mu 자동조정
    plt.figure(figsize=(6, 4))
    plt.plot(range(len(mu_small)), mu_small, "o-", color="C0", markersize=4,
             label=L["mu_small_short"])
    plt.plot(range(len(mu_large)), mu_large, "o-", color="C3", markersize=4,
             label=L["mu_large_short"])
    plt.xlabel("iteration")
    plt.ylabel(L["damping"], **F)
    plt.yscale("log")
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig2"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig2-mu-adaptation.png"), dpi=150)
    plt.close()

    # 그림3: 네 방법 비교
    plt.figure(figsize=(6, 4))
    plt.plot(d_gd, "o-", color="C0", markersize=4, label=L["gd"])
    plt.plot(d_nt, "o-", color="C3", markersize=4, label=L["newton"])
    plt.plot(d_gn, "o-", color="C2", markersize=4, label=L["gn"])
    plt.plot(d_lm, "x--", color="C4", markersize=6, label=L["lm"])
    plt.axhline(true_thickness, color="gray", linestyle="--", label=L["true_thickness"])
    plt.xlabel("iteration")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig3"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig3-four-methods-comparison.png"), dpi=150)
    plt.close()

    # 그림4: GN 실패 vs LM 성공
    plt.figure(figsize=(6, 4))
    plt.plot(d_gn_wrong, "o-", color="C3", markersize=4, label=L["gn_cycle"])
    plt.plot(d_lm_wrong, "o-", color="C4", markersize=4, label=L["lm_plain"])
    plt.axhline(wrong_floor_d, color="gray", linestyle="--", label=L["wrong_floor"])
    plt.xlabel("iteration")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig4"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig4-gn-failure-lm-rescue.png"), dpi=150)
    plt.close()


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

# ---------------------------------------------------------------------------
# 본문에 인용할 수치 출력 (언어와 무관한 수치이므로 한 번만 찍는다)
# ---------------------------------------------------------------------------
print()
print("=== mu0 comparison (d0=1540) ===")
print(f"tau=1e-6: mu0={mu_small[0]:.4g}, final d={d_small[-1]:.3f}")
print(f"tau=1e6 : mu0={mu_large[0]:.4g}, final d={d_large[-1]:.3f}")
print(f"A = sum(Jr^2) @ d0=1540: {mu_large[0] / 1e6:.6f}")
print()
print("=== four-method comparison (d0=1540, 10 steps) ===")
print(f"GD(alpha=300): final d={d_gd[-1]:.2f}")
print(f"Newton: final d={d_nt[-1]:.2f} (발산, 개별 반복값은 재현되지 않는다)")
print(f"Gauss-Newton: final d={d_gn[-1]:.2f}")
print(f"LM(tau=1e-6): final d={d_lm[-1]:.2f}")
print()
print("=== basin dependence: LM vs GN at d0=1300, 1690 ===")
for d0 in [1300.0, 1690.0]:
    dg, _ = gauss_newton(d0, wavelength_nm, R_measured, n_iter=30)
    row = [f"GN={dg[-1]:.2f}"]
    for tau in [1e-3, 1.0, 100.0]:
        db, _, _ = levenberg_marquardt(d0, wavelength_nm, R_measured, tau=tau, n_iter=30)
        row.append(f"LM(tau={tau:g})={db[-1]:.2f}")
    print(f"d0={d0:.0f}: " + "  ".join(row))
print()
print("=== failure case (wrong n1=1.02) ===")
Jr_ok = (reflectance(1490.0 + 1e-3, wavelength_nm) - reflectance(1490.0 - 1e-3, wavelength_nm)) / 2e-3
Jr_bad = jacobian_wrong(1490.0)
print(f"sum(Jr^2) @ d=1490  correct model: {np.sum(Jr_ok ** 2):.4e}")
print(f"sum(Jr^2) @ d=1490  wrong model  : {np.sum(Jr_bad ** 2):.4e}")
print(f"ratio: {np.sum(Jr_ok ** 2) / np.sum(Jr_bad ** 2):.0f}x")
print(f"objective_wrong floor: d={wrong_floor_d:.3f}, J={wrong_floor_J:.5f}")
print(f"GN cycles among: {np.round(d_gn_wrong[-3:], 1)}")
print(f"LM converges to: {d_lm_wrong[-1]:.3f}")

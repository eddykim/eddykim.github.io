"""3편 그림 4개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-levenberg-marquardt/ 에 fig1~4 저장

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

KFONT = {"fontfamily": "NanumGothic"}
LFONT = {"family": "NanumGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-levenberg-marquardt",
)
os.makedirs(OUT_DIR, exist_ok=True)

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0
noise_std = 0.004
initial_guess = 1540.0

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# ---------------------------------------------------------------------------
# 그림1·2 공통: mu0을 작게(tau=1e-6, GN에 가까움) vs 크게(tau=1e6, steepest
# descent에 가까움) 잡았을 때 LM이 어떻게 다르게 움직이는지 비교.
# ---------------------------------------------------------------------------
N_ITER_MU = 20
d_small, J_small, mu_small = levenberg_marquardt(
    initial_guess, wavelength_nm, R_measured, tau=1e-6, n_iter=N_ITER_MU)
d_large, J_large, mu_large = levenberg_marquardt(
    initial_guess, wavelength_nm, R_measured, tau=1e6, n_iter=N_ITER_MU)

plt.figure(figsize=(6, 4))
plt.plot(d_small, "o-", color="C0", markersize=4, label=r"$\mu_0$ 작게 (τ=1e-6, GN에 가까움)")
plt.plot(d_large, "o-", color="C3", markersize=4, label=r"$\mu_0$ 크게 (τ=1e6, steepest descent에 가까움)")
plt.axhline(true_thickness, color="gray", linestyle="--", label="실제 두께")
plt.xlabel("iteration")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title(r"그림1. $\mu_0$ 크기에 따른 LM 수렴 궤적 (d0=1540nm)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig1-mu0-comparison.png"), dpi=150)
plt.close()

plt.figure(figsize=(6, 4))
plt.plot(range(len(mu_small)), mu_small, "o-", color="C0", markersize=4, label=r"$\mu_0$ 작게 (τ=1e-6)")
plt.plot(range(len(mu_large)), mu_large, "o-", color="C3", markersize=4, label=r"$\mu_0$ 크게 (τ=1e6)")
plt.xlabel("iteration")
plt.ylabel(r"damping parameter $\mu$", **KFONT)
plt.yscale("log")
plt.legend(prop=LFONT, fontsize=8)
plt.title(r"그림2. gain ratio 기반 $\mu$ 자동조정 추이", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig2-mu-adaptation.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림3: GD / Newton / Gauss-Newton / LM(tau=1e-6) 네 방법 비교, d0=1540, 10스텝.
# ---------------------------------------------------------------------------
N_ITER = 10
d_gd, J_gd = gradient_descent(initial_guess, wavelength_nm, R_measured, alpha=300, n_iter=N_ITER)
d_nt, J_nt = newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)
d_gn, J_gn = gauss_newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)
d_lm, J_lm, _ = levenberg_marquardt(initial_guess, wavelength_nm, R_measured, tau=1e-6, n_iter=N_ITER)

plt.figure(figsize=(6, 4))
plt.plot(d_gd, "o-", color="C0", markersize=4, label="Gradient Descent (alpha=300)")
plt.plot(d_nt, "o-", color="C3", markersize=4, label="Newton (발산)")
plt.plot(d_gn, "o-", color="C2", markersize=4, label="Gauss-Newton")
plt.plot(d_lm, "x--", color="C4", markersize=6, label=r"Levenberg-Marquardt ($\tau$=1e-6)")
plt.axhline(true_thickness, color="gray", linestyle="--", label="실제 두께")
plt.xlabel("iteration")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림3. 네 방법 비교 (d0=1540nm)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig3-four-methods-comparison.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림4: 실패 사례 실증. 피팅 모델의 SiO2 굴절률을 1.02로 잘못 가정하면
# 간섭 콘트라스트가 거의 사라져 Jacobian이 거의 특이(near-singular)해진다.
# 이 상황에서 Gauss-Newton은 세 값 사이를 영원히 순환하며 수렴하지 않는다.
# ---------------------------------------------------------------------------
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
wrong_floor_d = 1486.63  # objective_wrong 최솟값 지점, 아래 print로 확인

plt.figure(figsize=(6, 4))
plt.plot(d_gn_wrong, "o-", color="C3", markersize=4, label="Gauss-Newton (영원히 진동)")
plt.plot(d_lm_wrong, "o-", color="C4", markersize=4, label="Levenberg-Marquardt")
plt.axhline(wrong_floor_d, color="gray", linestyle="--", label="이 모델의 실제 최솟값")
plt.xlabel("iteration")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림4. 굴절률을 잘못 가정한 모델에서: GN 실패 vs LM 성공", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4-gn-failure-lm-rescue.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 본문에 인용할 수치 출력
# ---------------------------------------------------------------------------
print("=== mu0 comparison (d0=1540) ===")
print(f"tau=1e-6: mu0={mu_small[0]:.4g}, final d={d_small[-1]:.3f}, steps to converge~{len(d_small)-1}")
print(f"tau=1e6 : mu0={mu_large[0]:.4g}, final d={d_large[-1]:.3f}, steps to converge~{len(d_large)-1}")
print()
print("=== four-method comparison (d0=1540, 10 steps) ===")
print(f"GD(alpha=300): final d={d_gd[-1]:.2f}")
print(f"Newton: final d={d_nt[-1]:.2f} (발산)")
print(f"Gauss-Newton: final d={d_gn[-1]:.2f}")
print(f"LM(tau=1e-6): final d={d_lm[-1]:.2f}")
print()
print("=== basin dependence: LM vs GN at d0=1300, 1690 ===")
for d0 in [1300.0, 1690.0]:
    d_hist_b, J_hist_b, _ = levenberg_marquardt(d0, wavelength_nm, R_measured, tau=1.0, n_iter=30)
    print(f"d0={d0}: LM final d={d_hist_b[-1]:.2f}, J={J_hist_b[-1]:.4g}")
print()
print("=== failure case (wrong n1=1.02) ===")
ds = np.linspace(1400, 1700, 3000)
Js = np.array([objective_wrong(d) for d in ds])
idx = np.argmin(Js)
print(f"objective_wrong floor: d={ds[idx]:.3f}, J={Js[idx]:.5f}")
print(f"GN cycles among: {np.round(d_gn_wrong[-4:], 1)}")
print(f"LM converges to: {d_lm_wrong[-1]:.3f}, J={objective_wrong(d_lm_wrong[-1]):.5f}")
Jr_wrong_at_1490 = jacobian_wrong(1490.0)
print(f"sum(Jr^2) wrong model @ d=1490: {np.sum(Jr_wrong_at_1490**2):.6g}")

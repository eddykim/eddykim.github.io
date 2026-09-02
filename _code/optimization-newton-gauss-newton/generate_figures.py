"""2편 그림 4개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-newton-gauss-newton/ 에 fig1~4 저장

측정 데이터는 1편과 동일한 방식(SiO2/Si 단층, 두께 1490nm, 가우시안 노이즈,
seed=0)으로 합성한다. 초기값 d0=1540nm도 1편과 동일하게 맞춰서, "1편에서 gradient
descent가 성공했던 바로 그 시작점에서 Newton법이 실제로는 실패한다"는 대비를
그대로 보여준다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from gradient_descent import gradient_descent
from gauss_newton import gauss_newton
from newton import newton, numerical_grad_hess, objective
from reflectance_model import reflectance

# Linux(나눔고딕)용 한글 폰트 설정. macOS라면 1편처럼 AppleGothic 등으로 교체할 것.
KFONT = {"fontfamily": "NanumGothic"}
LFONT = {"family": "NanumGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-newton-gauss-newton",
)
os.makedirs(OUT_DIR, exist_ok=True)

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0  # nm, 1편과 동일
noise_std = 0.004
initial_guess = 1540.0  # nm, 1편과 동일한 초기값

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# ---------------------------------------------------------------------------
# 그림1: 개념도 - 목적함수의 국소 2차(포물선) 근사
# 곡률이 정상적으로 양(+)인 지점(d=1500)에서 Newton의 국소 모델을 보여준다.
# ---------------------------------------------------------------------------
d_anchor = 1500.0
grad_a, hess_a = numerical_grad_hess(d_anchor, wavelength_nm, R_measured)
d_next = d_anchor - grad_a / hess_a

d_local = np.linspace(1480, 1520, 200)
J_local = np.array([objective(d, wavelength_nm, R_measured) for d in d_local])
J_quad = (
    objective(d_anchor, wavelength_nm, R_measured)
    + grad_a * (d_local - d_anchor)
    + 0.5 * hess_a * (d_local - d_anchor) ** 2
)

plt.figure(figsize=(6, 4))
plt.plot(d_local, J_local, "-", color="C0", label="실제 J(d)")
plt.plot(d_local, J_quad, "--", color="C3", label="d=1500에서의 2차 근사")
plt.axvline(d_anchor, color="gray", linestyle=":", label=f"현재점 d_k={d_anchor:.0f}")
plt.axvline(d_next, color="C3", linestyle=":", alpha=0.6, label=f"Newton 스텝 후 d={d_next:.1f}")
plt.xlabel("두께 d (nm)", **KFONT)
plt.ylabel("목적함수 J(d)", **KFONT)
plt.ylim(bottom=0)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림1. 목적함수의 국소 2차 근사와 Newton 스텝", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig1-quadratic-approximation.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 실험 공통: d0=1540에서 GD / Newton / Gauss-Newton 세 방법을 나란히 돌린다.
# ---------------------------------------------------------------------------
N_ITER = 10
d_hist_gd, J_hist_gd = gradient_descent(
    initial_guess, wavelength_nm, R_measured, alpha=300, n_iter=N_ITER
)
d_hist_nt, J_hist_nt = newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)
d_hist_gn, J_hist_gn = gauss_newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)

# ---------------------------------------------------------------------------
# 그림2: 두께 추정값 vs iteration - GD vs Newton vs Gauss-Newton (d0=1540)
# ---------------------------------------------------------------------------
plt.figure(figsize=(6, 4))
plt.plot(d_hist_gd, "o-", color="C0", markersize=4, label="Gradient Descent (alpha=300)")
plt.plot(d_hist_nt, "o-", color="C3", markersize=4, label="Newton (발산)")
plt.plot(d_hist_gn, "o-", color="C2", markersize=4, label="Gauss-Newton")
plt.axhline(true_thickness, color="gray", linestyle="--", label="실제 두께")
plt.xlabel("iteration")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림2. 두께 추정값 vs iteration (d0=1540nm)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig2-thickness-vs-iteration.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림3: 목적함수 J vs iteration (log) - GD vs Newton vs Gauss-Newton (d0=1540)
# ---------------------------------------------------------------------------
plt.figure(figsize=(6, 4))
plt.plot(J_hist_gd, "o-", color="C0", markersize=4, label="Gradient Descent (alpha=300)")
plt.plot(J_hist_nt, "o-", color="C3", markersize=4, label="Newton (발산)")
plt.plot(J_hist_gn, "o-", color="C2", markersize=4, label="Gauss-Newton")
plt.xlabel("iteration")
plt.ylabel("목적함수 J (error)", **KFONT)
plt.yscale("log")
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림3. 목적함수 J vs iteration (d0=1540nm)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig3-objective-vs-iteration.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림4: Newton의 두 얼굴 - 초기값에 따라 성공(d0=1500)/실패(d0=1540)가 갈린다.
# ---------------------------------------------------------------------------
d0_success = 1500.0
d_hist_success, J_hist_success = newton(d0_success, wavelength_nm, R_measured, n_iter=N_ITER)

plt.figure(figsize=(6, 4))
plt.plot(J_hist_success, "o-", color="C2", markersize=4, label=f"d0={d0_success:.0f} (성공, 2차수렴)")
plt.plot(J_hist_nt, "o-", color="C3", markersize=4, label=f"d0={initial_guess:.0f} (실패, 발산)")
plt.xlabel("iteration")
plt.ylabel("목적함수 J (error)", **KFONT)
plt.yscale("log")
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림4. Newton법의 두 얼굴 - 초기값에 따른 성공/실패", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4-newton-success-vs-failure.png"), dpi=150)
plt.close()

hess_bad = numerical_grad_hess(initial_guess, wavelength_nm, R_measured)[1]
hess_good = numerical_grad_hess(d0_success, wavelength_nm, R_measured)[1]

print(f"[d0={initial_guess}] Newton 최종 두께: {d_hist_nt[-1]:.2f} nm (실패)")
print(f"[d0={initial_guess}] Gauss-Newton 최종 두께: {d_hist_gn[-1]:.2f} nm")
print(f"[d0={initial_guess}] Gradient Descent 최종 두께: {d_hist_gd[-1]:.2f} nm")
print(f"[d0={d0_success}] Newton 최종 두께: {d_hist_success[-1]:.2f} nm (성공)")
print(f"Newton Hessian @ d0={initial_guess}: {hess_bad:.6f} (음수 -> 비볼록)")
print(f"Newton Hessian @ d0={d0_success}: {hess_good:.6f} (양수 -> 정상)")

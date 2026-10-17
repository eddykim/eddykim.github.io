"""최적화 2편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-newton-gauss-newton/      (한국어)
      ../../assets/img/posts/optimization-newton-gauss-newton/en/   (영문)

측정 데이터는 1편과 동일한 방식(SiO2/Si 단층, 두께 1490nm, 가우시안 노이즈,
seed=0)으로 합성한다. 초기값 d0=1540nm도 1편과 동일하게 맞춰서, "1편에서 gradient
descent가 성공했던 바로 그 시작점에서 Newton법이 실제로는 실패한다"는 대비를
그대로 보여준다.

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은 데이터가
완전히 동일하고 표기만 다르다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from gradient_descent import gradient_descent
from gauss_newton import gauss_newton
from newton import newton, numerical_grad_hess, objective
from reflectance_model import reflectance

# 한글 텍스트에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면 log축 눈금의 마이너스
# 기호가 한글 폰트에 없어 깨진다. 시리즈의 다른 편과 같은 AppleGothic 을 쓴다.
# macOS 전용 설정이다. 다른 OS는 나눔고딕 등 설치된 한글 폰트 이름으로 교체할 것.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-newton-gauss-newton",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR,
        "font": KFONT, "legend": LFONT,
        "thickness": "두께 d (nm)", "objective": "목적함수 J(d)",
        "actual_J": "실제 J(d)", "quad": "d={:.0f}에서의 2차 근사",
        "current": "현재점 d_k={:.0f}", "after_step": "Newton 스텝 후 d={:.1f}",
        "gd": "Gradient Descent (alpha=300)", "newton": "Newton (발산)",
        "gn": "Gauss-Newton", "true_thickness": "실제 두께",
        "thickness_est": "두께 추정값 (nm)",
        "objective_err": "목적함수 J (error)",
        "success": "d0={:.0f} (성공)", "failure": "d0={:.0f} (실패, 발산)",
        "p_nt": "[d0={:.0f}] Newton 최종 두께: {:.2f} nm (실패)",
        "p_gn": "[d0={:.0f}] Gauss-Newton 최종 두께: {:.2f} nm",
        "p_gd": "[d0={:.0f}] Gradient Descent 최종 두께: {:.2f} nm",
        "p_ok": "[d0={:.0f}] Newton 최종 두께: {:.2f} nm (성공)",
        "p_hess_bad": "Newton Hessian @ d0={:.0f}: {:.6f} (음수 -> 비볼록)",
        "p_hess_good": "Newton Hessian @ d0={:.0f}: {:.6f} (양수 -> 정상)",
        "p_chaos_head": "발산 궤적의 초기값 민감도 (10스텝째 두께):",
        "p_chaos_row": "  d0={:.12f} -> {:.2f} nm",
        "concave_quad": "d={:.0f}에서의 2차 근사 (위로 볼록)",
        "concave_top": "포물선 꼭대기 = Newton 스텝 후 d={:.1f}",
        "p_concave": "Newton 스텝 @ d0={:.0f}: {:+.2f} nm -> d={:.2f}, 포물선 꼭대기 J={:.2f}",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"),
        "font": EFONT, "legend": ELFONT,
        "thickness": "Thickness d (nm)", "objective": "Objective J(d)",
        "actual_J": "actual J(d)", "quad": "quadratic approximation at d={:.0f}",
        "current": "current point d_k={:.0f}", "after_step": "after Newton step, d={:.1f}",
        "gd": "Gradient Descent (alpha=300)", "newton": "Newton (diverging)",
        "gn": "Gauss-Newton", "true_thickness": "true thickness",
        "thickness_est": "Estimated thickness (nm)",
        "objective_err": "Objective J (error)",
        "success": "d0={:.0f} (converges)", "failure": "d0={:.0f} (fails, diverges)",
        "p_nt": "[d0={:.0f}] Newton final thickness: {:.2f} nm (failed)",
        "p_gn": "[d0={:.0f}] Gauss-Newton final thickness: {:.2f} nm",
        "p_gd": "[d0={:.0f}] Gradient Descent final thickness: {:.2f} nm",
        "p_ok": "[d0={:.0f}] Newton final thickness: {:.2f} nm (succeeded)",
        "p_hess_bad": "Newton Hessian @ d0={:.0f}: {:.6f} (negative -> non-convex)",
        "p_hess_good": "Newton Hessian @ d0={:.0f}: {:.6f} (positive -> well behaved)",
        "p_chaos_head": "sensitivity of the diverging trajectory to the initial value (10th-step thickness):",
        "p_chaos_row": "  d0={:.12f} -> {:.2f} nm",
        "concave_quad": "quadratic approximation at d={:.0f} (opens downward)",
        "concave_top": "top of the parabola = after Newton step, d={:.1f}",
        "p_concave": "Newton step @ d0={:.0f}: {:+.2f} nm -> d={:.2f}, top of parabola J={:.2f}",
    },
}

# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0  # nm, 1편과 동일
noise_std = 0.004
initial_guess = 1540.0  # nm, 1편과 동일한 초기값

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# 그림1용: 곡률이 정상적으로 양(+)인 지점(d=1500)에서 Newton의 국소 모델
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

# d0=1540에서 GD / Newton / Gauss-Newton 세 방법을 나란히 돌린다.
N_ITER = 10
d_hist_gd, J_hist_gd = gradient_descent(
    initial_guess, wavelength_nm, R_measured, alpha=300, n_iter=N_ITER
)
d_hist_nt, J_hist_nt = newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)
d_hist_gn, J_hist_gn = gauss_newton(initial_guess, wavelength_nm, R_measured, n_iter=N_ITER)

# 초기값에 따라 성공(d0=1500)/실패(d0=1540)가 갈린다.
d0_success = 1500.0
d_hist_success, J_hist_success = newton(d0_success, wavelength_nm, R_measured, n_iter=N_ITER)

hess_bad = numerical_grad_hess(initial_guess, wavelength_nm, R_measured)[1]
hess_good = numerical_grad_hess(d0_success, wavelength_nm, R_measured)[1]

# 발산 궤적은 Hessian 이 0 에 가까운 구간을 지나며 반올림 오차를 증폭시킨다.
# 초기값을 1e-12 nm 만 흔들어도 10스텝째 값이 수백 nm 달라진다는 것을 보여둔다.
# 개별 반복값을 인용해서는 안 된다는 근거다.
CHAOS_EPS = [0.0, 1e-12, 1e-10]
chaos = [
    (initial_guess + eps, newton(initial_guess + eps, wavelength_nm, R_measured, n_iter=N_ITER)[0][-1])
    for eps in CHAOS_EPS
]

# 그림5용: 곡률이 음(-)인 지점(d=1540)에서 Newton의 국소 모델.
# 포물선이 위로 볼록하므로 기울기가 0인 점은 바닥이 아니라 꼭대기이고,
# Newton 스텝은 그 꼭대기로 뛴다(약 +336nm).
grad_c, hess_c = numerical_grad_hess(initial_guess, wavelength_nm, R_measured)
J_c = objective(initial_guess, wavelength_nm, R_measured)
d_top = initial_guess - grad_c / hess_c
J_top = J_c + grad_c * (d_top - initial_guess) + 0.5 * hess_c * (d_top - initial_guess) ** 2
d_wide = np.linspace(1440, 1940, 501)
J_wide = np.array([objective(d, wavelength_nm, R_measured) for d in d_wide])
J_quad_c = J_c + grad_c * (d_wide - initial_guess) + 0.5 * hess_c * (d_wide - initial_guess) ** 2


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L):
    """주어진 라벨 묶음으로 그림 5개를 그리고 검증 수치를 출력한다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림1: 목적함수의 국소 2차(포물선) 근사
    plt.figure(figsize=(6, 4))
    plt.plot(d_local, J_local, "-", color="C0", label=L["actual_J"])
    plt.plot(d_local, J_quad, "--", color="C3", label=L["quad"].format(d_anchor))
    plt.axvline(d_anchor, ymax=0.6, color="gray", linestyle=":", label=L["current"].format(d_anchor))
    plt.axvline(d_next, ymax=0.6, color="C3", linestyle=":", alpha=0.6, label=L["after_step"].format(d_next))
    plt.xlabel(L["thickness"], **F)
    plt.ylabel(L["objective"], **F)
    plt.ylim(bottom=0)
    plt.legend(prop=LF, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig1-quadratic-approximation.png"), dpi=150)
    plt.close()

    # 그림2: 두께 추정값 vs iteration
    plt.figure(figsize=(6, 4))
    plt.plot(d_hist_gd, "o-", color="C0", markersize=4, label=L["gd"])
    plt.plot(d_hist_nt, "o-", color="C3", markersize=4, label=L["newton"])
    plt.plot(d_hist_gn, "o-", color="C2", markersize=4, label=L["gn"])
    plt.axhline(true_thickness, color="gray", linestyle="--", label=L["true_thickness"])
    plt.xlabel("iteration")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig2-thickness-vs-iteration.png"), dpi=150)
    plt.close()

    # 그림3: 목적함수 J vs iteration (log)
    plt.figure(figsize=(6, 4))
    plt.plot(J_hist_gd, "o-", color="C0", markersize=4, label=L["gd"])
    plt.plot(J_hist_nt, "o-", color="C3", markersize=4, label=L["newton"])
    plt.plot(J_hist_gn, "o-", color="C2", markersize=4, label=L["gn"])
    plt.xlabel("iteration")
    plt.ylabel(L["objective_err"], **F)
    plt.yscale("log")
    plt.legend(prop=LF, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig3-objective-vs-iteration.png"), dpi=150)
    plt.close()

    # 그림4: Newton의 두 얼굴
    plt.figure(figsize=(6, 4))
    plt.plot(J_hist_success, "o-", color="C2", markersize=4, label=L["success"].format(d0_success))
    plt.plot(J_hist_nt, "o-", color="C3", markersize=4, label=L["failure"].format(initial_guess))
    plt.xlabel("iteration")
    plt.ylabel(L["objective_err"], **F)
    plt.yscale("log")
    plt.legend(prop=LF, fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig4-newton-success-vs-failure.png"), dpi=150)
    plt.close()

    # 그림5: 곡률이 음수인 지점의 2차 근사 (위로 볼록한 포물선)
    plt.figure(figsize=(6, 4))
    plt.plot(d_wide, J_wide, "-", color="C0", label=L["actual_J"])
    plt.plot(d_wide, J_quad_c, "--", color="C3", label=L["concave_quad"].format(initial_guess))
    plt.axvline(initial_guess, ymax=0.62, color="gray", linestyle=":", label=L["current"].format(initial_guess))
    plt.axvline(d_top, ymax=0.62, color="C3", linestyle=":", alpha=0.6, label=L["concave_top"].format(d_top))
    plt.plot([d_top], [J_top], "o", color="C3", markersize=5)
    plt.xlabel(L["thickness"], **F)
    plt.ylabel(L["objective"], **F)
    plt.ylim(0, 30)  # 범례가 포물선 꼭대기를 가리지 않도록 위에 여유를 둔다
    plt.legend(prop=LF, fontsize=8, loc="upper left")
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig5-concave-quadratic-approximation.png"), dpi=150)
    plt.close()

    print(L["p_nt"].format(initial_guess, d_hist_nt[-1]))
    print(L["p_gn"].format(initial_guess, d_hist_gn[-1]))
    print(L["p_gd"].format(initial_guess, d_hist_gd[-1]))
    print(L["p_ok"].format(d0_success, d_hist_success[-1]))
    print(L["p_hess_bad"].format(initial_guess, hess_bad))
    print(L["p_hess_good"].format(d0_success, hess_good))
    print(L["p_chaos_head"])
    for d0, last in chaos:
        print(L["p_chaos_row"].format(d0, last))
    print(L["p_concave"].format(initial_guess, d_top - initial_guess, d_top, J_top))


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

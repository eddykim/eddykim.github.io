"""XPS 기초 5편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/xps-thickness-arxps-pitfalls/      (한국어)
      ../../assets/img/posts/xps-thickness-arxps-pitfalls/en/   (영문)

계산은 overlayer.py에 있고 여기서는 그리기만 한다. 그림 4의 스퍼터 프로파일은
파라미터를 임의로 둔 개념 모사다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

import overlayer as O

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "xps-thickness-arxps-pitfalls",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "vac": "진공", "ox": "산화막 (SiO$_2$)", "sub": "기판 (Si)", "normal": "법선",
        "path_o": "산화막 신호", "path_m": "기판 신호: 산화막을 지나며 감쇠",
        "fig1a": "(a) 균일 오버레이어 모델", "fig1b": r"(b) 두께에 따른 강도비 $I_o/I_m$",
        "d": "두께 $d$ (nm)", "R": r"강도비 $I_o / I_m$", "sub5": r"기판 신호 5% ($3\lambda$)",
        "fig1": "오버레이어 두께와 강도비 (SiO$_2$/Si, $\\lambda$ = 3 nm)",
        "fig2a": "(a) 통계 오차: 두께의 표준편차", "fig2b": "(b) 계통 오차: 역산 두께 대 참 두께",
        "err": "두께 오차 (nm)", "true_d": "참 두께 (nm)", "est_d": "역산 두께 (nm)",
        "r5": "강도비 5% 오차", "c4": "계수 잡음 (기판 10$^4$ 계수)", "c5": "계수 잡음 (기판 10$^5$ 계수)",
        "ideal": "참값", "lam_hi": r"$\lambda$ +10%", "lam_lo": r"$\lambda$ $-$10%",
        "r0": "$R_0$ 계산값 0.53 사용 (실측 0.88)",
        "fig2": "두께 역산의 정밀도와 정확도 (합성 계산)",
        "depth": "깊이 $z$ (nm)", "conc": "산화물 비율 $c(z)$",
        "step": "급격한 계면", "graded": "확산 계면, 폭 {w} nm",
        "ang": r"검출 각도 $\theta$ (°, 법선 기준)", "rel": "급격한 계면 대비 강도비 차이 (%)\n($R_0$ 불확도 안에서 상수 보정 뒤)",
        "band": "불확도 띠 (각도 ±0.5° + 계수 잡음)",
        "fig3a": "(a) 서로 다른 깊이 분포", "fig3b": "(b) 각도별 강도비의 차이",
        "fig3": "ARXPS: 각도 데이터로 계면 구조를 가를 수 있는가 (합성 계산)",
        "xo": "산소 원자 분율", "sd": "스퍼터 깊이 (nm)",
        "s_true": "참 분포", "s_mix": "원자 혼합 뒤", "s_obs": "+ 정보깊이 가중 (관측)",
        "s_red": "+ 산소 우선 제거",
        "fig4": "스퍼터 깊이 프로파일이 계면을 뭉개는 과정 (개념 모사, 파라미터 임의)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "vac": "vacuum", "ox": "oxide (SiO$_2$)", "sub": "substrate (Si)", "normal": "normal",
        "path_o": "oxide signal", "path_m": "substrate signal: attenuated by the oxide",
        "fig1a": "(a) Uniform overlayer model", "fig1b": r"(b) Intensity ratio $I_o/I_m$ vs thickness",
        "d": "Thickness $d$ (nm)", "R": r"Intensity ratio $I_o / I_m$",
        "sub5": r"substrate signal 5% ($3\lambda$)",
        "fig1": "Overlayer thickness and intensity ratio (SiO$_2$/Si, $\\lambda$ = 3 nm)",
        "fig2a": "(a) Statistical error: s.d. of thickness",
        "fig2b": "(b) Systematic error: inferred vs true thickness",
        "err": "Thickness error (nm)", "true_d": "True thickness (nm)",
        "est_d": "Inferred thickness (nm)",
        "r5": "5% error in ratio", "c4": "counting noise (10$^4$ substrate counts)",
        "c5": "counting noise (10$^5$ substrate counts)",
        "ideal": "true", "lam_hi": r"$\lambda$ +10%", "lam_lo": r"$\lambda$ $-$10%",
        "r0": "using calculated $R_0$ = 0.53 (measured 0.88)",
        "fig2": "Precision and accuracy of thickness inversion (synthetic calculation)",
        "depth": "Depth $z$ (nm)", "conc": "Oxide fraction $c(z)$",
        "step": "abrupt interface", "graded": "graded interface, width {w} nm",
        "ang": r"Emission angle $\theta$ (°, from normal)",
        "rel": "Ratio difference from abrupt interface (%)\n(after scaling within $R_0$ uncertainty)",
        "band": "uncertainty band (±0.5° angle + counting)",
        "fig3a": "(a) Different depth profiles", "fig3b": "(b) Difference in intensity ratio by angle",
        "fig3": "ARXPS: can angle data distinguish interface structure? (synthetic calculation)",
        "xo": "Oxygen atomic fraction", "sd": "Sputter depth (nm)",
        "s_true": "true profile", "s_mix": "after atomic mixing",
        "s_obs": "+ information-depth weighting (observed)", "s_red": "+ preferential O removal",
        "fig4": "How a sputter depth profile smears an interface (schematic, arbitrary parameters)",
    },
}

# ---------------------------------------------------------------------------
# 계산 (한 번만)
# ---------------------------------------------------------------------------
d_grid = np.linspace(0.05, 12, 400)
R0_ang = {t: O.ratio(d_grid, t) for t in (0, 60)}
dd_r5 = O.ratio_error(d_grid)
dd_c4 = O.counting_sigma(d_grid, 1e4)
dd_c5 = O.counting_sigma(d_grid, 1e5)
d_lam_hi = O.thickness(O.ratio(d_grid), lam=1.1 * O.L)
d_lam_lo = O.thickness(O.ratio(d_grid), lam=0.9 * O.L)
d_r0 = O.thickness(O.ratio(d_grid), r0=O.R0_CALC)

angles = np.arange(0, 71, 2.5)
D_STEP = 2.0
WIDTHS = (0.3, 0.5, 0.8, 1.2)
step_prof = O.profile_step(D_STEP)
step_ratio = O.arxps_ratio(step_prof, angles)
graded = {}
for w in WIDTHS:
    c0, sc, rel, _ = O.best_graded_match(D_STEP, w, angles)
    graded[w] = (c0, O.profile_graded(c0, w), 100 * rel, sc)
band = 100 * O.total_band(D_STEP, angles)

sd = np.linspace(0, 14, 1401)
s_true, s_mix, s_obs, s_red = O.sputter_profile(sd)


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L_):
    out, F, LF = L_["dir"], L_["font"], L_["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. 오버레이어 기하와 강도비
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.6, 4.3), gridspec_kw={"width_ratios": [1, 1.2]})
    a1.fill_between([0, 4], 0, 1.0, color="tab:orange", alpha=0.25)
    a1.fill_between([0, 4], -2.2, 0, color="0.8")
    a1.axhline(1.0, color="black", lw=1.2)
    a1.axhline(0, color="black", lw=0.8, ls="--")
    a1.text(0.1, 1.15, L_["vac"], fontsize=9, **F)
    a1.text(2.4, 0.45, L_["ox"], fontsize=9, **F)
    a1.text(0.1, -1.9, L_["sub"], fontsize=9, **F)
    a1.annotate("", xy=(3.75, 1.0), xytext=(3.75, 0),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1))
    a1.text(3.82, 0.5, "$d$", fontsize=12, va="center")
    th = np.radians(35)
    u = np.array([np.sin(th), np.cos(th)])
    exits = {}
    for key, start, col in (("o", np.array([0.9, 0.5]), "tab:orange"),
                            ("m", np.array([1.2, -1.0]), "tab:blue")):
        end = start + u * ((1.0 - start[1]) / u[1])
        exits[key] = end
        a1.plot(*zip(start, end), color=col, lw=2.2)
        a1.annotate("", xy=end + 0.6 * u, xytext=end, arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6))
        a1.plot(*start, "o", color=col)
    a1.text(*(exits["o"] + 0.7 * u + np.array([-0.15, 0.08])), L_["path_o"], color="tab:orange",
            fontsize=8.5, ha="right", **F)
    a1.text(*(exits["m"] + 0.7 * u + np.array([0.05, 0.08])), L_["path_m"], color="tab:blue",
            fontsize=8.5, ha="left", **F)
    xm = exits["m"][0]
    a1.plot([xm, xm], [1.0, 2.1], color="gray", lw=0.8, ls=":")
    a1.text(xm - 0.05, 2.15, L_["normal"], fontsize=8, color="gray", ha="center", **F)
    arc = np.linspace(np.pi / 2 - th, np.pi / 2, 20)
    a1.plot(xm + 0.45 * np.cos(arc), 1.0 + 0.45 * np.sin(arc), color="gray", lw=1)
    a1.text(xm + 0.12, 1.55, r"$\theta$", fontsize=11, color="gray")
    a1.set_xlim(0, 4.3)
    a1.set_ylim(-2.2, 3.0)
    a1.axis("off")
    a1.set_title(L_["fig1a"], fontsize=11, **F)
    for t, ls in ((0, "-"), (60, "--")):
        a2.semilogy(d_grid, R0_ang[t], color="black", ls=ls, lw=1.8, label=rf"$\theta$ = {t}°")
    a2.axvline(3 * O.L, color="tab:red", lw=1, ls=":")
    a2.text(3 * O.L - 0.15, 0.12, L_["sub5"], color="tab:red", fontsize=8.5, ha="right", **F)
    a2.set_xlim(0, 12)
    a2.set_ylim(0.05, 1e3)
    a2.set_xlabel(L_["d"], fontsize=11, **F)
    a2.set_ylabel(L_["R"], fontsize=11, **F)
    a2.legend(prop={**LF, "size": 9}, loc="upper left")
    a2.grid(alpha=0.3, which="both")
    a2.set_title(L_["fig1b"], fontsize=11, **F)
    fig.suptitle(L_["fig1"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-overlayer-geometry.png"), dpi=150)
    plt.close(fig)

    # 그림 2. 정밀도와 정확도
    fig2, (b1, b2) = plt.subplots(1, 2, figsize=(10.6, 4.3))
    b1.plot(d_grid, dd_r5, color="tab:purple", lw=1.8, label=L_["r5"])
    b1.plot(d_grid, dd_c4, color="tab:blue", lw=1.8, ls="--", label=L_["c4"])
    b1.plot(d_grid, dd_c5, color="tab:blue", lw=1.8, label=L_["c5"])
    b1.set_xlim(0, 12)
    b1.set_ylim(0, 1.2)
    b1.set_xlabel(L_["true_d"], fontsize=11, **F)
    b1.set_ylabel(L_["err"], fontsize=11, **F)
    b1.legend(prop={**LF, "size": 8.5}, loc="upper left")
    b1.grid(alpha=0.3)
    b1.set_title(L_["fig2a"], fontsize=11, **F)
    b2.plot(d_grid, d_grid, color="black", lw=1.2, label=L_["ideal"])
    b2.fill_between(d_grid, d_lam_lo, d_lam_hi, color="tab:green", alpha=0.25)
    b2.plot(d_grid, d_lam_hi, color="tab:green", lw=1.2, label=L_["lam_hi"])
    b2.plot(d_grid, d_lam_lo, color="tab:green", lw=1.2, ls="--", label=L_["lam_lo"])
    b2.plot(d_grid, d_r0, color="tab:red", lw=1.8, label=L_["r0"])
    b2.set_xlim(0, 12)
    b2.set_ylim(0, 14)
    b2.set_xlabel(L_["true_d"], fontsize=11, **F)
    b2.set_ylabel(L_["est_d"], fontsize=11, **F)
    b2.legend(prop={**LF, "size": 8.5}, loc="upper left")
    b2.grid(alpha=0.3)
    b2.set_title(L_["fig2b"], fontsize=11, **F)
    fig2.suptitle(L_["fig2"], fontsize=12, **F)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "fig2-thickness-sensitivity.png"), dpi=150)
    plt.close(fig2)

    # 그림 3. ARXPS
    fig3, (c1, c2) = plt.subplots(1, 2, figsize=(10.6, 4.3))
    colors = {0.3: "tab:green", 0.5: "tab:blue", 0.8: "tab:orange", 1.2: "tab:red"}
    c1.plot(O.Z, step_prof, color="black", lw=2, label=L_["step"])
    for w in WIDTHS:
        c1.plot(O.Z, graded[w][1], color=colors[w], lw=1.5, label=L_["graded"].format(w=w))
    c1.set_xlim(0, 6)
    c1.set_ylim(-0.05, 1.1)
    c1.set_xlabel(L_["depth"], fontsize=11, **F)
    c1.set_ylabel(L_["conc"], fontsize=11, **F)
    c1.legend(prop={**LF, "size": 8.5}, loc="upper right")
    c1.grid(alpha=0.3)
    c1.set_title(L_["fig3a"], fontsize=11, **F)
    c2.fill_between(angles, -band, band, color="0.85", label=L_["band"])
    c2.axhline(0, color="black", lw=1)
    for w in WIDTHS:
        c2.plot(angles, graded[w][2], color=colors[w], lw=1.8, marker="o", ms=3,
                label=L_["graded"].format(w=w))
    c2.set_xlim(0, 70)
    c2.set_xlabel(L_["ang"], fontsize=11, **F)
    c2.set_ylabel(L_["rel"], fontsize=11, **F)
    c2.set_ylim(-35, 35)
    c2.legend(prop={**LF, "size": 8}, loc="upper left")
    c2.grid(alpha=0.3)
    c2.set_title(L_["fig3b"], fontsize=11, **F)
    fig3.suptitle(L_["fig3"], fontsize=12, **F)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, "fig3-arxps.png"), dpi=150)
    plt.close(fig3)

    # 그림 4. 스퍼터 인공물
    fig4, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.plot(sd, s_true, color="black", lw=2, label=L_["s_true"])
    ax.plot(sd, s_mix, color="tab:blue", lw=1.5, ls="--", label=L_["s_mix"])
    ax.plot(sd, s_obs, color="tab:blue", lw=2, label=L_["s_obs"])
    ax.plot(sd, s_red, color="tab:red", lw=2, label=L_["s_red"])
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel(L_["sd"], fontsize=11, **F)
    ax.set_ylabel(L_["xo"], fontsize=11, **F)
    ax.legend(prop={**LF, "size": 9}, loc="upper right")
    ax.grid(alpha=0.3)
    ax.set_title(L_["fig4"], fontsize=11, **F)
    fig4.tight_layout()
    fig4.savefig(os.path.join(out, "fig4-sputter-artifact.png"), dpi=150)
    plt.close(fig4)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

print("\n그림 3 각도별 (폭 0.3 / 0.5 / 0.8 / 1.2 nm 차이, 띠)")
for i in range(0, angles.size, 4):
    print(f"  {angles[i]:5.1f}°: " + " / ".join(f"{graded[w][2][i]:+.2f}%" for w in WIDTHS)
          + f"  | 띠 ±{band[i]:.2f}%")
for w in WIDTHS:
    over = angles[np.abs(graded[w][2]) > band]
    print(f"  폭 {w} nm: 중심 {graded[w][0]:.2f} nm, 상수 {100 * (graded[w][3] - 1):+.2f}%, 띠를 넘는 각도 수 {over.size}/{angles.size}")

"""전자현미경 배경이론 5편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/electron-microscopy-tem-diffraction-contrast/      (한국어)
      ../../assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/   (영문)

역공간 길이는 nm^-1, 실공간 길이는 nm 로 계산한다.
그림1은 Ewald 작도, 그림2는 BF/DF 광선도(개념도), 그림3은 두빔 근사 계산이다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "electron-microscopy-tem-diffraction-contrast",
)

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREEN = "#27795b"
ORANGE = "#d68910"
GRAY = "#7f8c8d"
DARK = "#2c3e50"

# ── 물리 상수와 조건 ──────────────────────────────────────
H, M0, QE, C = 6.62607015e-34, 9.1093837015e-31, 1.602176634e-19, 2.99792458e8
D_HKL = 0.2               # nm, 격자면 간격 (일반적인 금속·반도체 수준)
LAMBDA_XRAY = 0.15406     # nm, Cu K-alpha
XI_G = 60.0               # nm, 소광거리. Si(111) 200 kV 수준의 대표값
THICKNESS = 90.0          # nm, 그림3(b) 의 시편 두께. 1.5 xi_g 라 s=0 에서 극대다


def wavelength(volts):
    """상대론 보정을 넣은 전자 파장 [nm]. 1편과 같은 식이다."""
    v = float(volts)
    p2 = 2.0 * M0 * QE * v * (1.0 + QE * v / (2.0 * M0 * C**2))
    return H / np.sqrt(p2) * 1e9


LAMBDA_E = wavelength(200e3)


def bragg_angle_deg(lam, d=D_HKL):
    """1차 브래그 각 [도]. 2 d sin(theta_B) = lambda."""
    return np.degrees(np.arcsin(lam / (2.0 * d)))


def ewald_arc(lam, x):
    """원점에서 출발하는 Ewald 구의 표면 높이 [nm^-1].

    반지름 1/lambda 인 구가 역격자 원점을 지난다. 원점 근처에서
    y = R - sqrt(R^2 - x^2) 이고, 파장이 짧을수록(R 이 클수록) 평평해진다.
    """
    r = 1.0 / lam
    inside = r**2 - np.asarray(x, dtype=float) ** 2
    return np.where(inside > 0, r - np.sqrt(np.clip(inside, 0, None)), np.nan)


def two_beam_intensity(t, s, xi_g=XI_G):
    """두빔 근사에서 회절파의 세기.

        I_g = sin^2(pi t s_eff) / (xi_g s_eff)^2,   s_eff = sqrt(s^2 + 1/xi_g^2)

    s = 0 이면 주기가 정확히 xi_g 인 진동이 되고, 이것이 두께 프린지다.
    """
    s_eff = np.sqrt(np.asarray(s, dtype=float) ** 2 + 1.0 / xi_g**2)
    return np.sin(np.pi * np.asarray(t, dtype=float) * s_eff) ** 2 / (xi_g * s_eff) ** 2


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "fig1": "Ewald 구가 평평해지면 반사가 한꺼번에 켜진다",
        "f1x": "역공간 $g_x$ (nm$^{-1}$)", "f1y": "역공간 $g_y$ (nm$^{-1}$)",
        "f1xray": "X선 (Cu Kα, λ = {:.3f} nm)\n구 반지름 1/λ = {:.1f} nm$^{{-1}}$",
        "f1elec": "전자 (200 kV, λ = {:.2f} pm)\n구 반지름 1/λ = {:.0f} nm$^{{-1}}$",
        "f1lat": "역격자점", "f1sphere": "Ewald 구",
        "f1bragg": "브래그 각 $θ_B$ = {:.2f}°",
        "f1hit": "구 위에 놓인 반사", "f1none": "정확히 놓이는 반사가 없다",
        "fig2": "대물 조리개를 어디에 두느냐가 명암을 뒤집는다",
        "f2bf": "명시야 (BF) — 투과파만 통과", "f2df": "암시야 (DF) — 회절파만 통과",
        "f2spec": "시편", "f2lens": "대물렌즈", "f2bfp": "후초점면",
        "f2ap": "대물 조리개", "f2img": "상면",
        "f2t": "투과파 000", "f2g": "회절파 g",
        "f2imgbf": "회절이 센 곳이 어둡다", "f2imgdf": "회절이 센 곳이 밝다",
        "fig3": "편위 파라미터가 두께 프린지와 회절 대비를 만든다",
        "f3a": "두께에 따른 회절파 세기", "f3ax": "시편 두께 t (nm)", "f3ay": "회절파 세기",
        "f3b": "편위 파라미터에 따른 세기 (t = {:.0f} nm)",
        "f3bx": "편위 파라미터 s (nm$^{-1}$)", "f3by": "회절파 세기",
        "f3s0": "s = 0 (브래그 정확히 맞춤)", "f3s1": "s = 0.01 nm$^{-1}$",
        "f3s2": "s = 0.02 nm$^{-1}$",
        "f3xi": "주기 = 소광거리 $ξ_g$ = {:.0f} nm",
        "f3peak": "s 가 조금만 틀어져도\n세기가 급격히 떨어진다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "A flattened Ewald sphere switches many reflections on at once",
        "f1x": "Reciprocal space $g_x$ (nm$^{-1}$)",
        "f1y": "Reciprocal space $g_y$ (nm$^{-1}$)",
        "f1xray": "X-ray (Cu Kα, λ = {:.3f} nm)\nsphere radius 1/λ = {:.1f} nm$^{{-1}}$",
        "f1elec": "Electron (200 kV, λ = {:.2f} pm)\nsphere radius 1/λ = {:.0f} nm$^{{-1}}$",
        "f1lat": "Reciprocal lattice point", "f1sphere": "Ewald sphere",
        "f1bragg": "Bragg angle $θ_B$ = {:.2f}°",
        "f1hit": "Reflections on the sphere", "f1none": "No reflection lies on it",
        "fig2": "Where the objective aperture sits inverts the contrast",
        "f2bf": "Bright field (BF) — transmitted beam only",
        "f2df": "Dark field (DF) — diffracted beam only",
        "f2spec": "Specimen", "f2lens": "Objective lens", "f2bfp": "Back focal plane",
        "f2ap": "Objective aperture", "f2img": "Image plane",
        "f2t": "Transmitted 000", "f2g": "Diffracted g",
        "f2imgbf": "Strongly diffracting regions go dark",
        "f2imgdf": "Strongly diffracting regions go bright",
        "fig3": "The excitation error makes thickness fringes and diffraction contrast",
        "f3a": "Diffracted intensity against thickness",
        "f3ax": "Specimen thickness t (nm)", "f3ay": "Diffracted intensity",
        "f3b": "Intensity against excitation error (t = {:.0f} nm)",
        "f3bx": "Excitation error s (nm$^{-1}$)", "f3by": "Diffracted intensity",
        "f3s0": "s = 0 (exact Bragg)", "f3s1": "s = 0.01 nm$^{-1}$",
        "f3s2": "s = 0.02 nm$^{-1}$",
        "f3xi": "Period = extinction distance $ξ_g$ = {:.0f} nm",
        "f3peak": "A small change in s\ndrops the intensity sharply",
    },
}


def _ewald_panel(ax, L, lam, key, lattice_g, xlim, ylim, color):
    """한 패널에 역격자와 Ewald 구를 그린다."""
    gx = np.arange(-4, 5) * lattice_g
    gy = np.arange(0, 5) * lattice_g
    xs, ys = np.meshgrid(gx, gy)
    ax.plot(xs.ravel(), ys.ravel(), "o", color=GRAY, ms=4.5, zorder=3)

    x = np.linspace(xlim[0], xlim[1], 3000)
    ax.plot(x, ewald_arc(lam, x), color=color, lw=2.4, zorder=4)

    # 구 위에 (거의) 놓인 역격자점을 표시한다.
    on = [(px, py) for px in gx for py in gy
          if not (px == 0 and py == 0)           # 000 은 투과파이지 반사가 아니다
          and abs(ewald_arc(lam, px) - py) < 0.35]
    if on:
        ax.plot([p[0] for p in on], [p[1] for p in on], "o", color=ACCENT,
                ms=9, zorder=5)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.set_xlabel(L["f1x"], **L["font"])
    ax.set_ylabel(L["f1y"], **L["font"])
    ax.grid(alpha=0.25)
    # 전자 패널만 파장을 pm 로 표시한다.
    shown = lam * 1000 if key == "f1elec" else lam
    ax.set_title(L[key].format(shown, 1.0 / lam), fontsize=10, **L["font"])
    return len(on)


def fig1_ewald(L):
    g0 = 1.0 / D_HKL                      # nm^-1
    xlim, ylim = (-16, 16), (-1.5, 12)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.2))

    n_x = _ewald_panel(ax1, L, LAMBDA_XRAY, "f1xray", g0, xlim, ylim, BLUE)
    n_e = _ewald_panel(ax2, L, LAMBDA_E, "f1elec", g0, xlim, ylim, ACCENT)

    ax1.text(0.5, 0.93, L["f1bragg"].format(bragg_angle_deg(LAMBDA_XRAY)),
             transform=ax1.transAxes, ha="center", fontsize=9.5, color=BLUE,
             **L["font"])
    ax2.text(0.5, 0.93, L["f1bragg"].format(bragg_angle_deg(LAMBDA_E)),
             transform=ax2.transAxes, ha="center", fontsize=9.5, color=ACCENT,
             **L["font"])
    ax1.text(0.5, 0.72, L["f1none"] if n_x == 0 else f"{L['f1hit']}: {n_x}",
             transform=ax1.transAxes, ha="center", fontsize=9.5, color=DARK,
             **L["font"])
    ax2.text(0.5, 0.72, f"{L['f1hit']}: {n_e}", transform=ax2.transAxes,
             ha="center", fontsize=9.5, color=DARK, **L["font"])

    handles = [plt.Line2D([], [], color=GRAY, marker="o", ls="", label=L["f1lat"]),
               plt.Line2D([], [], color=DARK, lw=2.0, label=L["f1sphere"])]
    fig.legend(handles=handles, loc="lower center", ncol=2,
               prop=L["legend"] or None, fontsize=9, frameon=False)
    fig.suptitle(L["fig1"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    return fig


def _bf_df_panel(ax, L, title, aperture_on_g, img_key):
    """BF 또는 DF 광선도 한 장."""
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-0.4, 8.4)
    ax.axis("off")
    ax.set_title(title, fontsize=10.5, **L["font"])

    Y_SPEC, Y_LENS, Y_BFP, Y_IMG = 7.5, 5.6, 3.6, 0.6
    dx = 0.85                                  # 후초점면에서 회절파가 벌어진 거리

    for y, key, col in [(Y_SPEC, "f2spec", DARK), (Y_LENS, "f2lens", DARK),
                        (Y_BFP, "f2bfp", GRAY), (Y_IMG, "f2img", DARK)]:
        ax.plot([-2.1, 2.1], [y, y], color=col, lw=1.6 if key != "f2bfp" else 1.0,
                ls="-" if key != "f2bfp" else "--")
        ax.text(-2.25, y, L[key], ha="right", va="center", fontsize=8.5,
                color=col, **L["font"])

    # 투과파는 축을 따라 곧장 내려가고, 회절파는 시편에서 각도를 갖고 갈라진다.
    # 렌즈가 한 번 꺾어 후초점면의 다른 높이를 지나게 하고, 상면에서 다시 만난다.
    x_lens = dx * (Y_LENS - Y_IMG) / (Y_BFP - Y_IMG)
    beams = [(0.0, 0.0, BLUE, "f2t"), (x_lens, dx, ACCENT, "f2g")]
    blocked_x = dx if not aperture_on_g else 0.0
    for xl, xb, col, key in beams:
        ax.plot([0, xl], [Y_SPEC, Y_LENS], color=col, lw=1.6)      # 시편 → 렌즈
        ax.plot([xl, xb], [Y_LENS, Y_BFP], color=col, lw=1.6)      # 렌즈 → 후초점면
        ax.text(xb + 0.12, Y_BFP - 0.42, L[key], fontsize=8.5, color=col,
                ha="left", **L["font"])
        if abs(xb - blocked_x) > 1e-9:
            ax.plot([xb, 0], [Y_BFP, Y_IMG], color=col, lw=1.6)    # 조리개 통과
        else:
            ax.plot([xb], [Y_BFP], "x", color=col, ms=11, mew=2.6, zorder=6)

    # 대물 조리개: 통과시킬 빔의 자리만 열어 둔다
    open_x = dx if aperture_on_g else 0.0
    for x0, w in [(-2.1, 2.1 + open_x - 0.22), (open_x + 0.22, 2.1 - open_x - 0.22)]:
        if w > 0:
            ax.add_patch(Rectangle((x0, Y_BFP - 0.13), w, 0.26,
                                   facecolor="#d5d8dc", edgecolor=DARK, lw=1.0))
    ax.text(-1.15, Y_BFP + 0.34, L["f2ap"], ha="center", fontsize=8.5, color=DARK,
            **L["font"])

    ax.text(0, Y_IMG - 0.35, L[img_key], ha="center", va="top", fontsize=9,
            color=DARK, **L["font"])


def fig2_bf_df(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 6.2))
    _bf_df_panel(ax1, L, L["f2bf"], False, "f2imgbf")
    _bf_df_panel(ax2, L, L["f2df"], True, "f2imgdf")
    fig.suptitle(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def fig3_excitation_error(L):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))

    t = np.linspace(0, 300, 2000)
    for s, color, key, ls in [(0.0, ACCENT, "f3s0", "-"),
                              (0.01, BLUE, "f3s1", "--"),
                              (0.02, GREEN, "f3s2", ":")]:
        ax.plot(t, two_beam_intensity(t, s), color=color, lw=2.0, ls=ls, label=L[key])
    ax.set_xlabel(L["f3ax"], **L["font"])
    ax.set_ylabel(L["f3ay"], **L["font"])
    ax.set_title(L["f3a"], **L["font"])
    ax.grid(alpha=0.3)
    ax.legend(prop=L["legend"] or None, fontsize=9, loc="upper right",
              framealpha=0.95)
    # s=0 일 때 주기가 소광거리와 같다는 것을 화살표로 표시한다.
    ax.annotate("", xy=(XI_G / 2, 1.08), xytext=(XI_G * 1.5, 1.08),
                arrowprops=dict(arrowstyle="<->", color=ACCENT, lw=1.4))
    ax.text(XI_G, 1.13, L["f3xi"].format(XI_G), ha="center", fontsize=9,
            color=ACCENT, **L["font"])
    ax.set_ylim(0, 1.55)

    s = np.linspace(-0.06, 0.06, 2000)
    ax2.plot(s, two_beam_intensity(THICKNESS, s), color=DARK, lw=2.2)
    ax2.set_xlabel(L["f3bx"], **L["font"])
    ax2.set_ylabel(L["f3by"], **L["font"])
    ax2.set_title(L["f3b"].format(THICKNESS), **L["font"])
    ax2.grid(alpha=0.3)
    peak = two_beam_intensity(THICKNESS, 0.0)
    ax2.annotate(L["f3peak"], xy=(0.008, two_beam_intensity(THICKNESS, 0.008)),
                 xytext=(0.030, peak * 0.86), fontsize=9, color=DARK, **L["font"],
                 arrowprops=dict(arrowstyle="->", color=DARK))

    fig.suptitle(L["fig3"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def main():
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for name, builder in [("fig1-ewald-sphere", fig1_ewald),
                              ("fig2-bright-dark-field", fig2_bf_df),
                              ("fig3-excitation-error", fig3_excitation_error)]:
            fig = builder(L)
            fig.savefig(os.path.join(L["dir"], name + ".png"), dpi=150,
                        bbox_inches="tight")
            plt.close(fig)
        print(f"[{lang}] 3 figures written to {L['dir']}")

    print()
    print(f"200 kV 전자 : lambda = {LAMBDA_E*1000:.3f} pm, 1/lambda = {1/LAMBDA_E:.1f} nm^-1")
    print(f"Cu K-alpha  : lambda = {LAMBDA_XRAY:.5f} nm, 1/lambda = {1/LAMBDA_XRAY:.2f} nm^-1")
    print(f"d = {D_HKL} nm 에서 브래그 각 : 전자 {bragg_angle_deg(LAMBDA_E):.3f}°, "
          f"X선 {bragg_angle_deg(LAMBDA_XRAY):.2f}°")
    print(f"Ewald 구 반지름 / 역격자 간격 : 전자 {(1/LAMBDA_E)/(1/D_HKL):.1f}배, "
          f"X선 {(1/LAMBDA_XRAY)/(1/D_HKL):.2f}배")
    for x in (5.0, 10.0, 15.0):
        print(f"  g_x={x:4.0f} nm^-1 에서 구의 높이 : 전자 {ewald_arc(LAMBDA_E, x):6.3f}, "
              f"X선 {ewald_arc(LAMBDA_XRAY, x):7.3f} nm^-1")
    print(f"두께 프린지 주기 (s=0) = 소광거리 {XI_G:.0f} nm")


if __name__ == "__main__":
    main()

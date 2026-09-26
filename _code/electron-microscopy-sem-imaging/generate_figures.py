"""전자현미경 배경이론 4편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/electron-microscopy-sem-imaging/      (한국어)
      ../../assets/img/posts/electron-microscopy-sem-imaging/en/   (영문)

그림1은 검출기 배치 개념도이고, 그림2와 그림3은 계산 결과다.
3편과 달리 몬테카를로가 필요 없다 — 이차전자 수율의 기울기 의존과 전체
수율 곡선은 모두 해석식으로 충분하다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "electron-microscopy-sem-imaging",
)

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREEN = "#27795b"
ORANGE = "#d68910"
GRAY = "#7f8c8d"
DARK = "#2c3e50"

# ── 그림2: 표면 형상 ──────────────────────────────────────
FEATURE_H = 200.0      # nm, 구조물 높이
FEATURE_W = 300.0      # nm, 평탄부 반폭
FEATURE_S = 45.0       # nm, 측벽이 서는 정도 (작을수록 가파르다)
SEC_CAP = 6.0          # sec(theta) 를 자르는 상한. 실제 표면에서 발산하지 않는다.

# ── 그림3: 전체 수율 곡선 ─────────────────────────────────
# (이름, delta_max, E_max[keV], eta) — 교재에 실린 대표값 수준이다.
# 둘 다 절연체다. 도체는 전하가 빠져나가므로 수율과 무관하게 대전하지 않는다.
YIELD_MATERIALS = [
    ("sio2", 2.4, 0.40, 0.15),
    ("pmma", 2.1, 0.30, 0.10),
]


def surface_profile(x):
    """가로 위치 x [nm] 에서의 표면 높이 [nm]. 측벽이 완만한 사다리꼴."""
    return FEATURE_H * 0.5 * (np.tanh((x + FEATURE_W) / FEATURE_S)
                              - np.tanh((x - FEATURE_W) / FEATURE_S))


def se_yield_vs_tilt(theta):
    """기울기 theta 에서의 이차전자 수율. 평탄면 대비 비율로 돌려준다.

    빔이 기울어진 면에 들어가면 이차전자 탈출 깊이(3편에서 본 ~5 nm) 안에
    머무는 경로 길이가 1/cos(theta) 배로 늘어난다. 그만큼 탈출할 수 있는
    이차전자가 많아지므로 수율이 sec(theta) 를 따른다.
    """
    return np.minimum(1.0 / np.cos(theta), SEC_CAP)


def total_yield(e0_kev, delta_max, e_max_kev, eta):
    """전체 전자 수율 sigma = delta + eta.

    delta 는 Joy 의 보편 수율 곡선을 쓴다.
        delta/delta_max = 1.28 (E/E_max)^-0.67 [1 - exp(-1.614 (E/E_max)^1.67)]
    eta 는 이 에너지 범위에서 거의 일정하다고 본다.
    """
    u = np.asarray(e0_kev, dtype=float) / e_max_kev
    delta = delta_max * 1.28 * u**-0.67 * (1.0 - np.exp(-1.614 * u**1.67))
    return delta + eta


def crossover_energies(delta_max, e_max_kev, eta):
    """전체 수율이 1 이 되는 두 지점 E1, E2 [keV] 를 수치로 찾는다."""
    e = np.logspace(-2.3, 1.7, 40000)
    s = total_yield(e, delta_max, e_max_kev, eta) - 1.0
    sign_change = np.where(np.diff(np.sign(s)))[0]
    return [float(e[i]) for i in sign_change]


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림1
        "fig1": "검출기가 어디에 놓이느냐가 무엇을 볼지 정한다",
        "f1beam": "입사 빔", "f1lens": "대물렌즈\n폴피스", "f1spec": "시편",
        "f1et": "ET 검출기\n(+200 V 그리드)", "f1inlens": "인렌즈 검출기",
        "f1bsed": "환형 BSE 검출기",
        "f1se": "이차전자 — 에너지가 낮아\n전기장에 휘어 끌려온다",
        "f1bse": "후방산란전자 — 에너지가 커서\n직진한다. 시선이 닿아야 잡힌다",
        # 그림2
        "fig2": "기울어진 면이 밝은 이유, 그리고 엣지 효과",
        "f2a": "기울기에 따른 이차전자 수율",
        "f2ax": "표면 기울기 θ (도)", "f2ay": "수율 (평탄면 = 1)",
        "f2sec": "$\\sec θ$", "f2cap": "실제 표면에서는 이만큼 오르지 않는다",
        "f2b": "사다리꼴 구조를 가로지른 라인 스캔",
        "f2bx": "가로 위치 (nm)", "f2by1": "표면 높이 (nm)", "f2by2": "SE 신호 (임의 단위)",
        "f2prof": "표면 형상", "f2sig": "SE 신호",
        "f2edge": "측벽에서 신호가 솟는다\n— 엣지 효과",
        # 그림3
        "fig3": "전체 수율이 1을 지나는 곳에서 대전이 사라진다",
        "f3x": "가속전압 (kV)", "f3y": "전체 수율 $\\sigma_T$ = δ + η",
        "f3sio2": "SiO$_2$ (절연체)", "f3pmma": "PMMA (고분자, 절연체)",
        "f3pos": "$\\sigma_T$ > 1: 나가는 전자가 더 많다 → 양(+)으로 대전",
        "f3neg": "$\\sigma_T$ < 1: 들어오는 전자가 더 많다 → 음(-)으로 대전",
        "f3e1": "$E_1$ = {:.0f} eV", "f3e2": "$E_2$ = {:.2f} kV",
        "f3op": "보통의 30 kV\n(강한 음 대전)",
        "f3note": "저전압 SEM 은 $E_2$ 근처에서 찍어\n코팅 없이 절연체를 본다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Where a detector sits decides what it sees",
        "f1beam": "Incident beam", "f1lens": "Objective\npole piece", "f1spec": "Specimen",
        "f1et": "ET detector\n(+200 V grid)", "f1inlens": "In-lens detector",
        "f1bsed": "Annular BSE detector",
        "f1se": "Secondary electrons — low energy,\nbent in and collected by the field",
        "f1bse": "Backscattered electrons — high energy,\ntravel straight. Line of sight required",
        "fig2": "Why tilted surfaces look bright, and the edge effect",
        "f2a": "Secondary electron yield against tilt",
        "f2ax": "Surface tilt θ (degrees)", "f2ay": "Yield (flat surface = 1)",
        "f2sec": "$\\sec θ$", "f2cap": "Real surfaces do not rise this far",
        "f2b": "Line scan across a trapezoidal feature",
        "f2bx": "Lateral position (nm)", "f2by1": "Surface height (nm)",
        "f2by2": "SE signal (arb. units)",
        "f2prof": "Surface profile", "f2sig": "SE signal",
        "f2edge": "Signal peaks on the sidewalls\n— the edge effect",
        "fig3": "Charging vanishes where the total yield passes one",
        "f3x": "Accelerating voltage (kV)", "f3y": "Total yield $\\sigma_T$ = δ + η",
        "f3sio2": "SiO$_2$ (insulator)", "f3pmma": "PMMA (polymer, insulator)",
        "f3pos": "$\\sigma_T$ > 1: more electrons leave than arrive → charges positive",
        "f3neg": "$\\sigma_T$ < 1: more electrons arrive than leave → charges negative",
        "f3e1": "$E_1$ = {:.0f} eV", "f3e2": "$E_2$ = {:.2f} kV",
        "f3op": "Typical 30 kV\n(strong negative charging)",
        "f3note": "Low-voltage SEM works near $E_2$,\nimaging insulators without a coating",
    },
}


def fig1_detector_layout(L):
    """검출기 배치 개념도.

    좌표는 임의 단위다. 광축은 x=0, 시편 표면은 y=0 에 둔다.
    BSE 화살표가 폴피스를 뚫고 지나가지 않도록 도달 지점을 직접 계산한다.
    """
    fig, ax = plt.subplots(figsize=(10.5, 7))
    ax.set_xlim(-6.6, 8.6)
    ax.set_ylim(-1.2, 6.4)
    ax.axis("off")

    BORE = 0.9          # 폴피스 사이 구멍의 반폭
    LENS_Y0, LENS_Y1 = 2.4, 4.3
    BSE_Y0, BSE_Y1 = 2.02, 2.34
    INLENS_Y = 4.55

    # 대물렌즈 폴피스
    for x0 in (-2.6, BORE):
        ax.add_patch(Rectangle((x0, LENS_Y0), 2.6 - BORE, LENS_Y1 - LENS_Y0,
                               facecolor="#d6e4f0", edgecolor=DARK, lw=1.4))
    ax.annotate(L["f1lens"], xy=(-2.6, 3.4), xytext=(-3.3, 3.4),
                ha="right", va="center", fontsize=9.5, color=DARK, **L["font"],
                arrowprops=dict(arrowstyle="-", color=DARK, lw=0.9))

    # 환형 BSE 검출기 — 렌즈 바로 아래, 광축을 둘러싼다
    for x0 in (-2.25, 0.7):
        ax.add_patch(Rectangle((x0, BSE_Y0), 1.55, BSE_Y1 - BSE_Y0,
                               facecolor="#f5d5d0", edgecolor=ACCENT, lw=1.3))
    ax.annotate(L["f1bsed"], xy=(-2.25, 2.18), xytext=(-3.3, 1.95),
                ha="right", va="center", fontsize=9.5, color=ACCENT, **L["font"],
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=0.9))

    # 인렌즈 검출기 — 구멍 위쪽, 축을 따라 올라온 전자만 받는다
    for x0 in (-0.82, 0.3):
        ax.add_patch(Rectangle((x0, INLENS_Y), 0.52, 0.3,
                               facecolor="#d5efe4", edgecolor=GREEN, lw=1.3))
    ax.annotate(L["f1inlens"], xy=(0.82, INLENS_Y + 0.15), xytext=(1.6, INLENS_Y + 0.5),
                ha="left", va="center", fontsize=9.5, color=GREEN, **L["font"],
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.9))

    # 시편
    ax.add_patch(Rectangle((-5.0, -0.85), 10.4, 0.85, facecolor="#e8e8e8",
                           edgecolor=DARK, lw=1.4))
    ax.text(-4.7, -0.42, L["f1spec"], ha="left", va="center", fontsize=10,
            color=DARK, **L["font"])

    # 입사 빔
    ax.add_patch(FancyArrowPatch((0, 6.2), (0, 0.02), arrowstyle="-|>",
                                 mutation_scale=17, color=DARK, lw=2.2))
    ax.text(0.18, 5.75, L["f1beam"], ha="left", fontsize=10, color=DARK, **L["font"])

    # ET 검출기 (측면)
    ax.add_patch(Rectangle((4.3, 0.55), 1.7, 1.15, facecolor="#fdf0d5",
                           edgecolor=ORANGE, lw=1.4))
    ax.text(5.15, 1.85, L["f1et"], ha="center", va="bottom", fontsize=9.5,
            color=ORANGE, **L["font"])

    # 후방산란전자 — 직진한다. 폴피스에 막히므로 도달 지점까지만 그린다.
    for ang in (58, 66, 114, 122):          # 환형 BSE 검출기로 가는 것
        a = np.radians(ang)
        t = BSE_Y0 / np.sin(a)
        ax.add_patch(FancyArrowPatch((0, 0.02), (t * np.cos(a), BSE_Y0),
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=ACCENT, lw=1.3, alpha=0.85))
    for ang in (82, 98):                     # 구멍을 지나 인렌즈로 가는 것
        a = np.radians(ang)
        t = INLENS_Y / np.sin(a)
        ax.add_patch(FancyArrowPatch((0, 0.02), (t * np.cos(a), INLENS_Y),
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=ACCENT, lw=1.3, alpha=0.85))
    ax.annotate(L["f1bse"], xy=(-1.05, 1.35), xytext=(-3.3, 0.85),
                ha="right", va="center", fontsize=9, color=ACCENT, **L["font"],
                arrowprops=dict(arrowstyle="-", color=ACCENT, lw=0.9))

    # 이차전자 — 에너지가 낮아 ET 의 전기장에 휘어 들어간다
    for y0 in (0.72, 0.98, 1.24):
        ax.add_patch(FancyArrowPatch((0.05, 0.02), (4.25, y0),
                                     connectionstyle="arc3,rad=-0.38",
                                     arrowstyle="-|>", mutation_scale=11,
                                     color=ORANGE, lw=1.4, alpha=0.95))
    # 영문 라벨이 한국어보다 길어 폴피스를 덮으므로 왼쪽 정렬로 오른쪽에 편다.
    ax.annotate(L["f1se"], xy=(3.3, 1.28), xytext=(2.85, 3.2),
                ha="left", va="bottom", fontsize=9, color=ORANGE, **L["font"],
                arrowprops=dict(arrowstyle="-", color=ORANGE, lw=0.9))

    ax.set_title(L["fig1"], fontsize=13, **L["font"])
    fig.tight_layout()
    return fig


def fig2_edge_effect(L):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))

    # (a) sec(theta) 법칙
    deg = np.linspace(0, 82, 400)
    th = np.radians(deg)
    raw = 1.0 / np.cos(th)
    ax.plot(deg, raw, color=GRAY, lw=1.4, ls=":", label=L["f2sec"])
    ax.plot(deg, se_yield_vs_tilt(th), color=BLUE, lw=2.3)
    ax.axhline(SEC_CAP, color=GRAY, ls="--", lw=1.0)
    ax.text(31, SEC_CAP + 0.28, L["f2cap"], fontsize=8.5, color=GRAY, **L["font"])
    for d in (0, 45, 60, 75):
        y = min(1 / np.cos(np.radians(d)), SEC_CAP)
        ax.plot([d], [y], "o", color=BLUE, ms=6, zorder=5)
        ax.annotate(f"{d}° : {y:.2f}", xy=(d, y), xytext=(6, -12),
                    textcoords="offset points", fontsize=8.5, color=DARK)
    ax.set_xlabel(L["f2ax"], **L["font"])
    ax.set_ylabel(L["f2ay"], **L["font"])
    ax.set_title(L["f2a"], **L["font"])
    ax.set_ylim(0, SEC_CAP + 1.3)
    ax.grid(alpha=0.3)
    ax.legend(prop=L["legend"] or None, loc="upper left", fontsize=9)

    # (b) 사다리꼴 구조 라인 스캔
    x = np.linspace(-700, 700, 2000)
    z = surface_profile(x)
    slope = np.gradient(z, x)
    signal = se_yield_vs_tilt(np.arctan(np.abs(slope)))

    ax2.plot(x, z, color=DARK, lw=2.0, label=L["f2prof"])
    ax2.set_xlabel(L["f2bx"], **L["font"])
    ax2.set_ylabel(L["f2by1"], **L["font"])
    ax2.set_ylim(-30, FEATURE_H * 1.35)
    ax2.grid(alpha=0.3)

    ax3 = ax2.twinx()
    ax3.plot(x, signal, color=BLUE, lw=2.2, label=L["f2sig"])
    ax3.set_ylabel(L["f2by2"], color=BLUE, **L["font"])
    ax3.tick_params(axis="y", labelcolor=BLUE)
    ax3.set_ylim(0.8, signal.max() * 1.45)

    peak = int(np.argmax(signal))
    ax3.annotate(L["f2edge"], xy=(x[peak], signal[peak]),
                 xytext=(x[peak] + 130, signal[peak] * 1.27),
                 fontsize=9, color=BLUE, ha="left", **L["font"],
                 arrowprops=dict(arrowstyle="->", color=BLUE))

    handles = [plt.Line2D([], [], color=DARK, lw=2.0, label=L["f2prof"]),
               plt.Line2D([], [], color=BLUE, lw=2.2, label=L["f2sig"])]
    ax2.legend(handles=handles, prop=L["legend"] or None, loc="upper left",
               fontsize=9)
    ax2.set_title(L["f2b"], **L["font"])

    fig.suptitle(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def fig3_charging(L):
    e = np.logspace(-2.2, 1.6, 2000)
    fig, ax = plt.subplots(figsize=(10, 5.2))

    styles = {"sio2": (ACCENT, "-", "f3sio2"), "pmma": (BLUE, "--", "f3pmma")}
    crossings = {}
    for name, dmax, emax, eta in YIELD_MATERIALS:
        color, ls, key = styles[name]
        ax.plot(e, total_yield(e, dmax, emax, eta), color=color, lw=2.3, ls=ls,
                label=L[key])
        crossings[name] = crossover_energies(dmax, emax, eta)

    ax.axhline(1.0, color=DARK, lw=1.3)
    ax.set_xscale("log")
    ax.set_xlim(e[0], e[-1])
    ax.set_ylim(0, 3.0)

    # 절연체(SiO2)의 두 교차점
    e1, e2 = crossings["sio2"][0], crossings["sio2"][-1]
    for ev, key, fmt_val in [(e1, "f3e1", e1 * 1000), (e2, "f3e2", e2)]:
        ax.plot([ev], [1.0], "o", color=ACCENT, ms=8, zorder=6)
        ax.annotate(L[key].format(fmt_val), xy=(ev, 1.0), xytext=(0, -26),
                    textcoords="offset points", ha="center", fontsize=9.5,
                    color=ACCENT, **L["font"])

    ax.axvspan(e1, e2, color=ORANGE, alpha=0.12, zorder=0)
    ax.text(np.sqrt(e1 * e2), 2.72, L["f3pos"], ha="center", fontsize=9,
            color="#9c6500", **L["font"])
    ax.text(e[-1] * 0.92, 0.07, L["f3neg"], ha="right", va="bottom", fontsize=9,
            color=DARK, **L["font"])

    ax.plot([30], [total_yield(30, *YIELD_MATERIALS[0][1:])], "v",
            color=DARK, ms=9, zorder=6)
    ax.annotate(L["f3op"], xy=(30, total_yield(30, *YIELD_MATERIALS[0][1:])),
                xytext=(0, 22), textcoords="offset points", ha="center",
                fontsize=9, color=DARK, **L["font"],
                arrowprops=dict(arrowstyle="->", color=DARK))

    ax.annotate(L["f3note"], xy=(e2, 1.0), xytext=(e2 * 2.1, 1.85),
                fontsize=9, color=ACCENT, **L["font"],
                arrowprops=dict(arrowstyle="->", color=ACCENT))

    ax.set_xlabel(L["f3x"], **L["font"])
    ax.set_ylabel(L["f3y"], **L["font"])
    ax.set_title(L["fig3"], fontsize=13, **L["font"])
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop=L["legend"] or None, loc="upper right")
    fig.tight_layout()
    return fig


def main():
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for name, builder in [("fig1-detector-layout", fig1_detector_layout),
                              ("fig2-edge-effect", fig2_edge_effect),
                              ("fig3-charging-crossover", fig3_charging)]:
            fig = builder(L)
            fig.savefig(os.path.join(L["dir"], name + ".png"), dpi=150,
                        bbox_inches="tight")
            plt.close(fig)
        print(f"[{lang}] 3 figures written to {L['dir']}")

    # 본문에 인용한 수치를 검산할 수 있게 출력한다.
    print()
    for d in (0, 30, 45, 60, 75):
        print(f"기울기 {d:2d}° : 수율 x{1/np.cos(np.radians(d)):.2f}")
    print()
    for name, dmax, emax, eta in YIELD_MATERIALS:
        cr = crossover_energies(dmax, emax, eta)
        txt = ", ".join(f"{c*1000:.0f} eV" if c < 1 else f"{c:.2f} kV" for c in cr)
        print(f"{name:>5} (delta_max={dmax}, E_max={emax} keV, eta={eta}) "
              f": 교차점 {txt}, 30 kV 에서 sigma={total_yield(30, dmax, emax, eta):.3f}")


if __name__ == "__main__":
    main()

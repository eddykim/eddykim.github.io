"""전자현미경 배경이론 7편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/electron-microscopy-stem-analytical/      (한국어)
      ../../assets/img/posts/electron-microscopy-stem-analytical/en/   (영문)

그림2 의 Z 지수는 3편에서 쓴 차폐 러더퍼드 단면적을 환형 검출기 범위로
적분해서 얻는다. 새 가정을 넣지 않고 3편의 차폐 상수를 그대로 쓴다.
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
    "assets", "img", "posts", "electron-microscopy-stem-analytical",
)

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREEN = "#27795b"
ORANGE = "#d68910"
PURPLE = "#7d3c98"
GRAY = "#7f8c8d"
DARK = "#2c3e50"

E_KEV = 200.0            # 가속전압 [keV]
THETA_MAX = 200e-3       # rad, 검출기 바깥 반각
ELEMENTS = [("C", 6), ("Si", 14), ("Cu", 29), ("Ag", 47), ("W", 74), ("Au", 79)]

# 검출기 구역 (내각 범위 [mrad], 색, 라벨 키)
DETECTOR_BANDS = [
    (0, 10, BLUE, "f2bf"),
    (10, 20, GREEN, "f2abf"),
    (20, 50, ORANGE, "f2adf"),
    (50, 200, ACCENT, "f2haadf"),
]


def screening(z, e_kev=E_KEV):
    """3편에서 쓴 차폐 상수 alpha_s = 3.4e-3 Z^0.67 / E."""
    return 3.4e-3 * np.asarray(z, dtype=float) ** 0.67 / e_kev


def annular_cross_section(z, theta_1, theta_2=THETA_MAX):
    """내각 theta_1, 외각 theta_2 인 환형 검출기가 받는 산란 단면적 (상대값).

    차폐 러더퍼드 미분 단면적을 입체각으로 적분하면 닫힌 형태가 나온다.
        dsigma/dOmega ∝ Z^2 / (sin^2(theta/2) + alpha_s)^2
        sigma ∝ Z^2 [ 1/(sin^2(t1/2)+alpha_s) - 1/(sin^2(t2/2)+alpha_s) ]
    """
    z = np.asarray(z, dtype=float)
    a = screening(z)
    u1, u2 = np.sin(theta_1 / 2) ** 2, np.sin(theta_2 / 2) ** 2
    return z**2 * (1.0 / (u1 + a) - 1.0 / (u2 + a))


def z_exponent(theta_1):
    """sigma ∝ Z^n 의 유효 지수 n. 여섯 원소에 로그-로그 직선을 맞춘다."""
    zs = np.array([z for _, z in ELEMENTS], dtype=float)
    s = annular_cross_section(zs, theta_1)
    return float(np.polyfit(np.log(zs), np.log(s), 1)[0])


def screening_angle_mrad(z):
    """차폐가 듣기 시작하는 각 theta_0 = 2 sqrt(alpha_s) [mrad]."""
    return 2.0 * np.sqrt(screening(z)) * 1e3


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림1
        "fig1": "광원과 검출기를 맞바꾸면 같은 상이 된다",
        "f1tem": "TEM", "f1stem": "STEM (관례대로 뒤집어 그림)",
        "f1src": "광원", "f1det": "검출기", "f1spec": "시편", "f1lens": "렌즈",
        "f1angles_tem": "위쪽 각 = 조명 각\n아래쪽 각 = 수집 각",
        "f1angles_stem": "위쪽 각 = 수집 각\n아래쪽 각 = 조명 각",
        "f1note": "광선 경로는 같고 광원과 검출기의 자리만 바뀐다 — 상반성 원리",
        # 그림2
        "fig2": "검출기 내각이 원자번호 의존성을 정한다",
        "f2a": "내각에 따른 유효 지수", "f2ax": "검출기 내각 (mrad)",
        "f2ay": "유효 지수 n  ($\\sigma \\propto Z^n$)",
        "f2bf": "BF", "f2abf": "ABF", "f2adf": "ADF", "f2haadf": "HAADF",
        "f2scr": "차폐각 $θ_0$\nSi {:.0f} / Au {:.0f} mrad",
        "f2mark": "HAADF 상용 조건\n{:.0f} mrad → n = {:.2f}",
        "f2b": "원소별 단면적 (내각 두 조건)",
        "f2bx": "원자번호 Z", "f2by": "환형 검출기 단면적 (임의 단위)",
        "f2low": "내각 20 mrad (n = {:.2f})", "f2high": "내각 100 mrad (n = {:.2f})",
        # 그림3
        "fig3": "4D-STEM은 화소마다 회절 패턴 한 장을 저장한다",
        "f3scan": "주사 위치 (x, y)", "f3diff": "회절 패턴 ($k_x$, $k_y$)",
        "f3cube": "네 축짜리 데이터",
        "f3bf": "가운데만 합\n→ 명시야",
        "f3adf": "고리만 합\n→ 환형 암시야",
        "f3dpc": "무게중심 이동\n→ 전기장 지도",
        "f3ptycho": "패턴 전체로 위상 역산\n→ 타이코그래피",
        "f3note": "같은 한 번의 측정에서 뽑아내는 상들이다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Swap source and detector and the image is the same",
        "f1tem": "TEM", "f1stem": "STEM (drawn inverted, as is conventional)",
        "f1src": "Source", "f1det": "Detector", "f1spec": "Specimen", "f1lens": "Lens",
        "f1angles_tem": "Upper angle = illumination\nLower angle = collection",
        "f1angles_stem": "Upper angle = collection\nLower angle = illumination",
        "f1note": "The ray paths are identical; only source and detector swap places — reciprocity",
        "fig2": "The detector inner angle sets the atomic number dependence",
        "f2a": "Effective exponent against inner angle",
        "f2ax": "Detector inner angle (mrad)",
        "f2ay": "Effective exponent n  ($\\sigma \\propto Z^n$)",
        "f2bf": "BF", "f2abf": "ABF", "f2adf": "ADF", "f2haadf": "HAADF",
        "f2scr": "Screening angle $θ_0$\nSi {:.0f} / Au {:.0f} mrad",
        "f2mark": "Typical HAADF\n{:.0f} mrad → n = {:.2f}",
        "f2b": "Cross section by element, two inner angles",
        "f2bx": "Atomic number Z", "f2by": "Annular cross section (arb. units)",
        "f2low": "Inner angle 20 mrad (n = {:.2f})",
        "f2high": "Inner angle 100 mrad (n = {:.2f})",
        "fig3": "4D-STEM stores a whole diffraction pattern per pixel",
        "f3scan": "Scan position (x, y)", "f3diff": "Diffraction pattern ($k_x$, $k_y$)",
        "f3cube": "A four-axis dataset",
        "f3bf": "Sum the centre\n→ bright field",
        "f3adf": "Sum the ring\n→ annular dark field",
        "f3dpc": "Centre-of-mass shift\n→ field map",
        "f3ptycho": "Invert phase from the pattern\n→ ptychography",
        "f3note": "All of these come out of one and the same measurement",
    },
}


def _reciprocity_panel(ax, L, reverse):
    """TEM 또는 STEM 광선도 한 장. reverse=True 면 STEM 이다."""
    ax.set_xlim(-2.4, 2.4)
    ax.set_ylim(-1.4, 8.2)
    ax.axis("off")

    Y_TOP, Y_LENS1, Y_SPEC, Y_LENS2, Y_BOT = 7.6, 5.9, 4.0, 2.1, 0.4
    top_key, bot_key = ("f1det", "f1src") if reverse else ("f1src", "f1det")

    for y, key in [(Y_TOP, top_key), (Y_SPEC, "f1spec"), (Y_BOT, bot_key)]:
        ax.plot([-1.9, 1.9], [y, y], color=DARK, lw=1.7)
        ax.text(-2.05, y, L[key], ha="right", va="center", fontsize=9.5,
                color=DARK, **L["font"])
    for y in (Y_LENS1, Y_LENS2):
        ax.plot([-1.4, 1.4], [y, y], color=BLUE, lw=1.5)
        ax.text(1.55, y, L["f1lens"], ha="left", va="center", fontsize=8.5,
                color=BLUE, **L["font"])

    # 상반성의 요점은 두 광선도의 경로가 **같다**는 것이다. 기하를 바꾸면
    # 오히려 대응이 보이지 않으므로, 광선은 한 벌만 그리고 광원과 검출기의
    # 자리만 맞바꾼다. 그래서 STEM 은 관례대로 뒤집어 그린다.
    for off in (-1.0, -0.4, 0.4, 1.0):
        pts = [(off * 0.12, Y_TOP), (off * 0.12, Y_LENS1), (off * 0.5, Y_SPEC),
               (off, Y_LENS2), (off, Y_BOT)]
        for (x0, y0), (x1, y1) in zip(pts[:-1], pts[1:]):
            ax.plot([x0, x1], [y0, y1], color=ACCENT, lw=1.3, alpha=0.85)

    # 진행 방향 화살표 — 두 패널에서 반대로 그린다
    ay0, ay1 = (Y_TOP - 0.35, Y_TOP - 1.1) if not reverse else (Y_BOT + 1.1, Y_BOT + 0.35)
    ax.add_patch(FancyArrowPatch((1.75, ay0), (1.75, ay1), arrowstyle="-|>",
                                 mutation_scale=15, color=DARK, lw=1.8))

    ax.text(0, -0.15, L["f1angles_stem"] if reverse else L["f1angles_tem"],
            ha="center", va="top", fontsize=9, color=GRAY, **L["font"])

    ax.set_title(L["f1stem"] if reverse else L["f1tem"], fontsize=11.5,
                 **L["font"])


def fig1_reciprocity(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 6.4))
    _reciprocity_panel(ax1, L, reverse=False)
    _reciprocity_panel(ax2, L, reverse=True)
    fig.suptitle(L["fig1"], fontsize=13, **L["font"])
    fig.text(0.5, 0.03, L["f1note"], ha="center", fontsize=9.5, color=DARK,
             **L["font"])
    fig.subplots_adjust(left=0.06, right=0.97, top=0.88, bottom=0.09, wspace=0.28)
    return fig


def fig2_z_exponent(L):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8))

    t1 = np.linspace(5, 180, 500)
    n = np.array([z_exponent(t * 1e-3) for t in t1])
    for lo, hi, color, key in DETECTOR_BANDS:
        ax.axvspan(lo, min(hi, 180), color=color, alpha=0.12, zorder=0)
        ax.text((lo + min(hi, 180)) / 2, 2.02, L[key], ha="center", fontsize=9,
                color=color, **L["font"])
    ax.plot(t1, n, color=DARK, lw=2.4, zorder=3)

    th_si, th_au = screening_angle_mrad(14), screening_angle_mrad(79)
    ax.axvline(th_si, color=GRAY, ls=":", lw=1.2)
    ax.axvline(th_au, color=GRAY, ls=":", lw=1.2)
    ax.text(th_au + 3, 1.42, L["f2scr"].format(th_si, th_au), fontsize=8.5,
            color=GRAY, **L["font"])

    n80 = z_exponent(80e-3)
    ax.plot([80], [n80], "o", color=ACCENT, ms=8, zorder=5)
    ax.annotate(L["f2mark"].format(80, n80), xy=(80, n80), xytext=(105, 1.72),
                fontsize=9, color=ACCENT, **L["font"],
                arrowprops=dict(arrowstyle="->", color=ACCENT))

    ax.set_xlabel(L["f2ax"], **L["font"])
    ax.set_ylabel(L["f2ay"], **L["font"])
    ax.set_title(L["f2a"], **L["font"])
    ax.set_xlim(0, 180)
    ax.set_ylim(1.3, 2.12)
    ax.grid(alpha=0.3)

    zs = np.array([z for _, z in ELEMENTS], dtype=float)
    for t, color, key in [(20e-3, ORANGE, "f2low"), (100e-3, ACCENT, "f2high")]:
        s = annular_cross_section(zs, t)
        s = s / s[0]
        ax2.loglog(zs, s, "o-", color=color, lw=1.9, ms=7,
                   label=L[key].format(z_exponent(t)))
    # W(74) 와 Au(79) 는 로그축에서 붙으므로 라벨을 좌우로 갈라 놓는다.
    offsets = {"W": (-13, 6), "Au": (11, 6)}
    ref = annular_cross_section(zs, 100e-3)
    for (name, z), s in zip(ELEMENTS, ref / ref[0]):
        dx, dy = offsets.get(name, (0, 9))
        ax2.annotate(name, xy=(z, s), xytext=(dx, dy), textcoords="offset points",
                     ha="center", fontsize=9, color=DARK)
    ax2.set_xlabel(L["f2bx"], **L["font"])
    ax2.set_ylabel(L["f2by"], **L["font"])
    ax2.set_title(L["f2b"], **L["font"])
    ax2.grid(alpha=0.3, which="both")
    ax2.legend(prop=L["legend"] or None, fontsize=9, loc="upper left")

    fig.suptitle(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def fig3_4dstem(L):
    fig, ax = plt.subplots(figsize=(11, 6.0))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7.4)
    ax.axis("off")

    # 왼쪽: 주사 격자. 각 화소가 회절 패턴 한 장을 갖는다.
    for i in range(4):
        for j in range(4):
            ax.add_patch(Rectangle((0.5 + i * 0.62, 4.3 + j * 0.62), 0.58, 0.58,
                                   facecolor="#dfe6ec", edgecolor=DARK, lw=0.8))
    ax.add_patch(Rectangle((0.5 + 2 * 0.62, 4.3 + 1 * 0.62), 0.58, 0.58,
                           facecolor="#f5d5d0", edgecolor=ACCENT, lw=1.8))
    ax.text(1.74, 3.95, L["f3scan"], ha="center", va="top", fontsize=9.5,
            color=DARK, **L["font"])

    # 오른쪽: 그 화소 하나의 회절 패턴
    cx, cy, r = 6.0, 5.2, 1.25
    ax.add_patch(Rectangle((cx - r, cy - r), 2 * r, 2 * r,
                           facecolor="#111111", edgecolor=DARK, lw=1.2))
    ax.add_patch(plt.Circle((cx, cy), 0.30, color="#f7f7f7", zorder=4))
    for rad, alpha in [(0.62, 0.55), (0.92, 0.30)]:
        ax.add_patch(plt.Circle((cx, cy), rad, fill=False, color="#f0a868",
                                lw=6, alpha=alpha, zorder=3))
    ax.text(cx, cy - r - 0.35, L["f3diff"], ha="center", va="top", fontsize=9.5,
            color=DARK, **L["font"])
    ax.add_patch(FancyArrowPatch((3.15, 5.2), (cx - r - 0.15, 5.2),
                                 arrowstyle="-|>", mutation_scale=16,
                                 color=ACCENT, lw=2.0))
    ax.text((3.15 + cx - r) / 2, 5.45, L["f3cube"], ha="center", fontsize=9,
            color=ACCENT, **L["font"])

    # 아래: 같은 데이터에서 뽑아내는 상들
    outs = [(1.5, "f3bf", BLUE), (4.3, "f3adf", ORANGE),
            (7.1, "f3dpc", GREEN), (9.9, "f3ptycho", PURPLE)]
    for x, key, color in outs:
        ax.add_patch(FancyArrowPatch((cx, cy - r - 0.75), (x + 0.75, 2.05),
                                     arrowstyle="-|>", mutation_scale=12,
                                     color=color, lw=1.4, alpha=0.8,
                                     connectionstyle="arc3,rad=0.12"))
        ax.add_patch(Rectangle((x - 0.15, 1.15), 1.8, 1.0, facecolor="white",
                               edgecolor=color, lw=1.5))
        ax.text(x + 0.75, 1.65, L[key], ha="center", va="center", fontsize=8.5,
                color=color, **L["font"])

    ax.text(6.0, 0.72, L["f3note"], ha="center", fontsize=9.5, color=DARK,
            **L["font"])
    ax.set_title(L["fig3"], fontsize=13, **L["font"])
    fig.tight_layout()
    return fig


def main():
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for name, builder in [("fig1-reciprocity", fig1_reciprocity),
                              ("fig2-z-exponent", fig2_z_exponent),
                              ("fig3-4dstem", fig3_4dstem)]:
            fig = builder(L)
            fig.savefig(os.path.join(L["dir"], name + ".png"), dpi=150,
                        bbox_inches="tight")
            plt.close(fig)
        print(f"[{lang}] 3 figures written to {L['dir']}")

    print()
    print("차폐각 theta_0 = 2 sqrt(alpha_s):")
    for name, z in [("Si", 14), ("Au", 79)]:
        print(f"  {name} (Z={z}) : {screening_angle_mrad(z):.1f} mrad")
    print("\n내각에 따른 유효 지수 n (sigma ∝ Z^n):")
    for t in (20, 30, 50, 80, 100, 150):
        print(f"  {t:4d} mrad : n = {z_exponent(t * 1e-3):.3f}")


if __name__ == "__main__":
    main()

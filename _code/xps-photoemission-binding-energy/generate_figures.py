"""XPS 기초 1편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/xps-photoemission-binding-energy/      (한국어)
      ../../assets/img/posts/xps-photoemission-binding-energy/en/   (영문)

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은 데이터가
완전히 동일하고 표기만 다르다. 스펙트럼은 전부 합성(모사)이며 실측이 아니다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from math import erf

# 한글 텍스트에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면 눈금의 마이너스
# 기호가 한글 폰트에 없어 깨진다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "xps-photoemission-binding-energy",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "sample": "시료 (도체)", "spectrometer": "분광기",
        "core": "코어 준위", "valence": "가전자대",
        "fermi": r"페르미 준위 $E_F$ (공통)",
        "vac_s": r"진공 준위 (시료)", "vac_sp": r"진공 준위 (분광기)",
        "contact": "전기적 접촉\n→ $E_F$ 정렬",
        "fig1": "광전 방출의 에너지 준위와 페르미 준위 정렬",
        "xlabel_be": "결합에너지 (eV)", "ylabel": "계수율 (임의 단위)",
        "fig2": "합성 서베이 스펙트럼 (가상 SiO$_2$/Si 시료, Al K$\\alpha$)",
        "auger": "오제",
        "fig3": "화학이동: 에틸 트리플루오로아세테이트형 C 1s (모사)",
        "xlabel_rel": r"상대 결합에너지 (CH$_3$ 기준, eV)",
        "observed": "관측 스펙트럼", "envelope": "성분 합",
        "more_neg": "주변 원자의 전기음성도 ↑  →  결합에너지 ↑",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "sample": "Sample (conductor)", "spectrometer": "Spectrometer",
        "core": "core level", "valence": "valence band",
        "fermi": r"Fermi level $E_F$ (common)",
        "vac_s": "vacuum level (sample)", "vac_sp": "spectrometer\nvacuum level",
        "contact": "electrical contact\n→ $E_F$ aligned",
        "fig1": "Energy levels in photoemission and Fermi-level alignment",
        "xlabel_be": "Binding energy (eV)", "ylabel": "Count rate (arb. units)",
        "fig2": "Synthetic survey spectrum (hypothetical SiO$_2$/Si, Al K$\\alpha$)",
        "auger": "Auger",
        "fig3": "Chemical shift: ethyl-trifluoroacetate-like C 1s (simulated)",
        "xlabel_rel": r"Relative binding energy (vs. CH$_3$, eV)",
        "observed": "observed", "envelope": "sum of components",
        "more_neg": "more electronegative neighbors  →  higher binding energy",
    },
}

HV_AL = 1486.6  # Al Kα 광자 에너지 (eV)


def gaussian(x, center, fwhm, area):
    """면적이 area인 Gaussian 피크."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return area / (sigma * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def peak_with_step(be, center, fwhm, area, step_frac):
    """광전자 피크 + 높은 결합에너지 쪽으로 이어지는 계단형 배경.

    비탄성 산란으로 에너지를 잃은 전자는 운동에너지가 낮아지므로,
    결합에너지 축에서는 피크의 왼쪽(높은 BE 쪽)에 쌓인다.
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    step = 0.5 * (1 + np.vectorize(erf)((be - center) / (np.sqrt(2) * sigma)))
    return gaussian(be, center, fwhm, area) + step_frac * area * step


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------
rng = np.random.default_rng(0)

# 그림2: 가상 SiO2(얇은 산화막)/Si 시료 + 표면 오염 탄소. 위치는 대략값이다.
be = np.linspace(0, 1350, 6750)
survey_peaks = [
    # (라벨, BE, FWHM, 면적, 계단 비율, 종류)
    ("O 1s", 533.0, 2.0, 900.0, 0.020, "pe"),
    ("C 1s", 284.8, 2.0, 150.0, 0.015, "pe"),
    ("Si 2s", 154.0, 2.5, 180.0, 0.015, "pe"),
    ("Si 2p", 103.3, 2.0, 200.0, 0.012, "pe"),
    ("", 99.4, 2.0, 60.0, 0.010, "pe"),          # 산화막 아래 Si(0)
    ("O 2s", 25.0, 3.0, 60.0, 0.010, "pe"),
    ("O KLL", HV_AL - 508.0, 8.0, 350.0, 0.020, "auger"),
    ("C KLL", HV_AL - 263.0, 8.0, 40.0, 0.010, "auger"),
]
survey = np.full_like(be, 8.0) + 0.004 * be       # 완만한 이차전자 배경
for _, c, w, a, s, _k in survey_peaks:
    survey += peak_with_step(be, c, w, a, s)
survey_noisy = rng.poisson(survey * 20) / 20

# 그림3: C 1s 네 성분. 인접 간격 1.7~3.1 eV, 전체 폭 약 8 eV
# (Travnikova et al. 2012)에 맞춘 모사값이다.
rel = np.linspace(-3, 11, 1400)
c1s_components = [
    (r"CH$_3$", 0.0), (r"O–CH$_2$", 1.7), (r"O–C=O", 4.8), (r"CF$_3$", 7.9),
]
C1S_FWHM, C1S_AREA = 1.0, 42.6   # 피크 높이 약 40
c1s_parts = [gaussian(rel, pos, C1S_FWHM, C1S_AREA) for _, pos in c1s_components]
c1s_sum = np.sum(c1s_parts, axis=0) + 2.0
c1s_obs = rng.poisson(c1s_sum * 30) / 30


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L_):
    """주어진 라벨 묶음으로 그림 3개를 그린다."""
    out, F, LF = L_["dir"], L_["font"], L_["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. 에너지 준위도: 시료와 분광기가 페르미 준위를 공유한다
    fig, ax = plt.subplots(figsize=(7.6, 5.8))
    y_core, y_vb_top, y_f = 0.0, 5.6, 6.0
    phi_s, phi_sp = 2.0, 1.3               # 도식용 크기 (시료 일함수 > 분광기 일함수)
    y_vs, y_vsp = y_f + phi_s, y_f + phi_sp
    sx, px = (0.0, 3.0), (5.0, 8.0)        # 시료·분광기 기둥의 x 범위

    ax.fill_between(sx, 4.4, y_vb_top, color="tab:gray", alpha=0.25)
    ax.text(sx[0] + 0.6, 4.95, L_["valence"], fontsize=9, **F)
    ax.plot(sx, [y_core] * 2, color="black", lw=2.5)
    ax.text(sx[0] + 0.1, y_core - 0.45, L_["core"], fontsize=9, **F)
    ax.plot([sx[0], px[1]], [y_f] * 2, color="tab:blue", lw=1.5, ls="--")
    ax.text(px[1], y_f - 0.35, L_["fermi"], fontsize=9, color="tab:blue",
            ha="right", **F)
    ax.plot(sx, [y_vs] * 2, color="black", lw=1.2)
    ax.text(sx[0] + 0.6, y_vs + 0.12, L_["vac_s"], fontsize=9, **F)
    ax.plot(px, [y_vsp] * 2, color="black", lw=1.2)
    ax.text(px[0] + 0.1, y_vsp + 0.12, L_["vac_sp"], fontsize=9, ha="left", **F)

    ax.text(1.5, 11.6, L_["sample"], ha="center", fontsize=11, weight="bold", **F)
    ax.text(6.5, 11.6, L_["spectrometer"], ha="center", fontsize=11, weight="bold", **F)
    ax.text(4.0, y_f + 0.25, L_["contact"], ha="center", fontsize=8.5,
            color="tab:blue", **F)

    y_top = y_core + 10.8                  # 광전자의 최종 에너지 준위 (hν 위)

    def vline(x, y0, y1, text, color, side="right", style="<->", ty=None):
        ax.annotate("", xy=(x, y1), xytext=(x, y0),
                    arrowprops=dict(arrowstyle=style, color=color, lw=1.5))
        dx = 0.12 if side == "right" else -0.12
        ax.text(x + dx, (y0 + y1) / 2 if ty is None else ty, text, color=color, fontsize=11,
                ha="left" if side == "right" else "right", va="center")

    vline(0.4, y_core, y_top, r"$h\nu$", "tab:purple", side="left", style="-|>", ty=9.4)
    vline(1.9, y_core, y_f, r"$E_B$", "black")
    vline(2.6, y_f, y_vs, r"$\phi_s$", "tab:gray")
    vline(5.4, y_f, y_vsp, r"$\phi_{sp}$", "tab:gray")
    vline(7.3, y_vsp, y_top, r"$E_K$", "tab:red")
    ax.plot([0.4, 7.3], [y_top] * 2, color="tab:red", lw=1, ls=":")
    ax.text(4.1, y_top - 0.9, r"$E_K = h\nu - E_B - \phi_{sp}$",
            fontsize=12, color="tab:red", ha="center", va="top")

    ax.set_xlim(-0.6, 8.6)
    ax.set_ylim(-1.0, 12.3)
    ax.axis("off")
    ax.set_title(L_["fig1"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-energy-diagram.png"), dpi=150)
    plt.close(fig)

    # 그림 2. survey 스펙트럼 (결합에너지 축을 오른쪽→왼쪽으로 증가)
    fig2, ax2 = plt.subplots(figsize=(8.0, 4.4))
    ax2.plot(be, survey_noisy, color="black", lw=0.7)
    for label, c, _w, _a, _s, kind in survey_peaks:
        if not label:
            continue
        y = survey[np.argmin(np.abs(be - c))]
        text = f"{label}\n({L_['auger']})" if kind == "auger" else label
        ax2.annotate(text, xy=(c, y), xytext=(c, y + 45), ha="center",
                     fontsize=9, color="tab:red" if kind == "auger" else "tab:blue",
                     arrowprops=dict(arrowstyle="-", color="gray", lw=0.6), **F)
    ax2.set_xlim(1350, 0)
    ax2.set_ylim(0, survey.max() * 1.35)
    ax2.set_xlabel(L_["xlabel_be"], fontsize=11, **F)
    ax2.set_ylabel(L_["ylabel"], fontsize=11, **F)
    ax2.set_yticks([])
    ax2.set_title(L_["fig2"], fontsize=12, **F)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "fig2-survey-schematic.png"), dpi=150)
    plt.close(fig2)

    # 그림 3. C 1s 화학이동 성분 분해
    fig3, ax3 = plt.subplots(figsize=(7.4, 4.6))
    colors = ["tab:green", "tab:orange", "tab:purple", "tab:red"]
    ax3.plot(rel, c1s_obs, "o", ms=2, color="gray", alpha=0.6, label=L_["observed"])
    for (name, pos), part, col in zip(c1s_components, c1s_parts, colors):
        ax3.fill_between(rel, 2.0, part + 2.0, color=col, alpha=0.35)
        ax3.text(pos, part.max() + 4, name, ha="center", fontsize=10, color=col)
    ax3.plot(rel, c1s_sum, color="black", lw=1.4, label=L_["envelope"])
    ax3.annotate("", xy=(9.8, 55), xytext=(-1.8, 55),
                 arrowprops=dict(arrowstyle="-|>", color="black", lw=1))
    ax3.text(4.0, 57.5, L_["more_neg"], ha="center", fontsize=9, **F)
    ax3.set_xlim(11, -3)
    ax3.set_ylim(0, 74)
    ax3.set_xlabel(L_["xlabel_rel"], fontsize=11, **F)
    ax3.set_ylabel(L_["ylabel"], fontsize=11, **F)
    ax3.set_yticks([])
    ax3.legend(prop=LF, fontsize=9, loc="upper center", ncol=2, frameon=False)
    ax3.set_title(L_["fig3"], fontsize=12, **F)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, "fig3-chemical-shift.png"), dpi=150)
    plt.close(fig3)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

print()
for label, c, *_ in survey_peaks:
    if label.endswith("KLL"):
        print(f"{label}: BE {c:.1f} eV  (KE {HV_AL - c:.1f} eV)")
print("C 1s 성분 간격:", np.diff([p for _, p in c1s_components]))

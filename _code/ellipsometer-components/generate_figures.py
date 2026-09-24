"""타원계측기 2편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-components/      (한국어)
      ../../assets/img/posts/ellipsometer-components/en/   (영문)

계산에 무작위 요소가 없으므로 두 언어의 그림은 데이터가 완전히 동일하고 표기만 다르다.
지연량 계산은 retarder.py 에서 오고, 그 구현은 verify_retarder.py 가
존스 경로와 뮬러 경로를 대조해 검증한 것이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon

from retarder import (DN_QUARTZ, critical_angle, effective_fast_axis_jones,
                      effective_retardance_jones, jones_stack, plate_thickness,
                      retardance, tir_phase_difference)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-components",
)

# 분광 엘립소미터가 흔히 다루는 대역
SE_BAND_UM = (0.20, 1.0)
# 그림 3·4 의 파장 축
WL = np.linspace(200e-9, 1000e-9, 801)
# 두 소자 보상자 설계값 (본문 참고: Johs 등의 착상을 따르되 수치는 이 글의 설계다)
TWO_EL = dict(l1=300e-9, l2=550e-9, offset=np.deg2rad(67.5))

C_SIMPLE, C_TWO, C_MULTI, C_RHOMB = "tab:blue", "tab:red", "tab:orange", "tab:green"
DANGER = "#c0392b"


def two_element_curves():
    """두 소자 보상자의 유효 지연량과 유효 빠른 축을 파장별로 계산."""
    d1 = plate_thickness(0.25, TWO_EL["l1"])
    d2 = plate_thickness(0.25, TWO_EL["l2"])
    de, ax = [], []
    for lam in WL:
        J = jones_stack([retardance(lam, d1), retardance(lam, d2)],
                        [0.0, TWO_EL["offset"]])
        de.append(np.rad2deg(effective_retardance_jones(J)))
        ax.append(np.rad2deg(effective_fast_axis_jones(J)))
    return np.array(de), np.array(ax)


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림 1
        "fig1": "부품이 쓸 수 있는 파장 영역 — 대역이 부품을 고른다",
        "band": "분광 엘립소미터가 흔히 쓰는 대역",
        "wl_um": "파장 (μm)",
        "items": [
            ("방해석 (글랜-테일러 편광자)", 0.21, 5.0, False),
            ("MgF$_2$ (로숑 편광자·보상자)", 0.12, 7.0, True),
            ("운모 (보상자)", 0.29, 7.0, True),
            ("와이어그리드 편광자", 2.5, 100.0, False),
        ],
        "openend": "위쪽 한계는 출처에 명시되지 않음",
        # 그림 2
        "fig2": "편광자 세 가지 — 무엇이 자외선에서 탈락시키는가",
        "gt": "글랜-테일러", "gth": "글랜-톰슨", "roch": "로숑",
        "gt_note": "사이가 공기층\n자외선까지 통과",
        "gth_note": "사이가 접착제\n접착제가 자외선을 흡수한다",
        "roch_note": "광로가 꺾이지 않는다\n회전시켜도 빔이 안 흔들린다",
        "air": "공기", "glue": "접착제", "calcite": "방해석", "quartz": "MgF$_2$",
        "reject": "제거되는 성분", "pass": "투과축 성분",
        # 그림 3
        "fig3": "파장이 넓어지면 단순 파장판은 무너진다",
        "fig3a": "(a) 유효 지연량",
        "fig3b": "(b) 대가 — 유효 빠른 축이 파장 따라 흔들린다",
        "wl_nm": "파장 (nm)",
        "ret_deg": "유효 지연량 (deg)",
        "axis_deg": "유효 빠른 축 (deg)",
        "simple": "단순 영차 $\\lambda/4$판 (550 nm 설계)",
        "two": "두 소자 보상자 (300 nm + 550 nm, 축차 67.5°)",
        "fatal": "0°·180°·360° — 정보가 사라지는 지연량",
        "cross": "여기서 180°를 가로지른다\n자외선에서 못 쓴다",
        "stay": "0°·180° 에 닿지 않는다",
        "axis_swing": "축이 {:.0f}° 움직인다\n보정이 복잡해지는 이유",
        # 그림 4
        "fig4": "영차와 다중차 — 같은 $\\lambda/4$판인데 파장을 넓히면 다르다",
        "zero": "영차 ($N=1/4$, 두께 {:.1f} μm)",
        "multi": "다중차 ($N=2.25$, 두께 {:.1f} μm)",
        "design": "설계 파장 550 nm",
        "fig4a": "(a) 지연량 자체 — 세로축이 로그다",
        "fig4b": "(b) 360°로 나눈 나머지 — 실제로 편광에 작용하는 값",
        "ret_mod": "지연량 mod 360° (deg)",
        "multi_note": "다중차는 이 대역에서\n치명적인 값을 {:d}번 지나간다",
        # 그림 5
        "fig5": "프레넬 롬 — 복굴절 없이 전반사만으로 지연을 만든다",
        "aoi": "롬 내부 입사각 (deg)",
        "ps_deg": "전반사 한 번의 p–s 위상차 크기 (deg)",
        "nlabel": "$n$ = {:.2f}",
        "target": "22.5° — 네 번 반사하면 $\\lambda/4$",
        "crit": "임계각",
        "achro": "굴절률이 1.50→1.52로 변해도\n곡선이 거의 겹친다 (무채색)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Usable spectral range of each component — the band picks the part",
        "band": "band commonly used in spectroscopic ellipsometry",
        "wl_um": "wavelength (μm)",
        "items": [
            ("calcite (Glan-Taylor polarizer)", 0.21, 5.0, False),
            ("MgF$_2$ (Rochon polarizer, compensator)", 0.12, 7.0, True),
            ("mica (compensator)", 0.29, 7.0, True),
            ("wire-grid polarizer", 2.5, 100.0, False),
        ],
        "openend": "upper limit not stated in the source",
        "fig2": "Three polarizers — what rules one out in the ultraviolet",
        "gt": "Glan-Taylor", "gth": "Glan-Thompson", "roch": "Rochon",
        "gt_note": "air gap between prisms\npasses the ultraviolet",
        "gth_note": "cemented with glue\nthe glue absorbs the ultraviolet",
        "roch_note": "the beam is not deviated\nrotating it does not shake the beam",
        "air": "air", "glue": "glue", "calcite": "calcite", "quartz": "MgF$_2$",
        "reject": "rejected", "pass": "transmitted",
        "fig3": "A simple wave plate breaks down once the band gets wide",
        "fig3a": "(a) effective retardance",
        "fig3b": "(b) the cost — the effective fast axis swings with wavelength",
        "wl_nm": "wavelength (nm)",
        "ret_deg": "effective retardance (deg)",
        "axis_deg": "effective fast axis (deg)",
        "simple": "simple zero-order $\\lambda/4$ plate (designed at 550 nm)",
        "two": "two-element compensator (300 + 550 nm, axes 67.5° apart)",
        "fatal": "0°, 180°, 360° — retardances that destroy the information",
        "cross": "crosses 180° here\nunusable in the ultraviolet",
        "stay": "never reaches 0° or 180°",
        "axis_swing": "the axis moves {:.0f}°\nwhich is why calibration gets harder",
        "fig4": "Zero-order vs multi-order — same $\\lambda/4$ plate, different over a band",
        "zero": "zero order ($N=1/4$, {:.1f} μm thick)",
        "multi": "multi order ($N=2.25$, {:.1f} μm thick)",
        "design": "design wavelength 550 nm",
        "fig4a": "(a) retardance itself — note the logarithmic axis",
        "fig4b": "(b) retardance mod 360° — what actually acts on the polarization",
        "ret_mod": "retardance mod 360° (deg)",
        "multi_note": "the multi-order plate sweeps past\nthe fatal values {:d} times here",
        "fig5": "The Fresnel rhomb — retardance from total internal reflection alone",
        "aoi": "angle of incidence inside the rhomb (deg)",
        "ps_deg": "p–s phase difference per reflection (deg)",
        "nlabel": "$n$ = {:.2f}",
        "target": "22.5° — four reflections give $\\lambda/4$",
        "crit": "critical angle",
        "achro": "the curves nearly coincide as $n$ goes\n1.50 → 1.52 (achromatic)",
    },
}


# ---------------------------------------------------------------------------
# 그림 1 — 파장 영역 막대
# ---------------------------------------------------------------------------

def figure1(L):
    fig, ax = plt.subplots(figsize=(9.4, 4.0))
    items = L["items"]
    colors = ["#5b8ff9", "#61ddaa", "#9270ca", "#f6bd16"]
    for i, ((name, lo, hi, open_end), col) in enumerate(zip(items, colors)):
        y = len(items) - 1 - i
        ax.barh(y, hi - lo, left=lo, height=0.5, color=col, ec="#444", lw=0.9)
        if open_end:
            ax.add_patch(FancyArrowPatch((hi, y), (hi * 1.5, y), arrowstyle="-|>",
                                         mutation_scale=13, color=col, lw=2.0))
        ax.text(lo * 0.92, y, name, ha="right", va="center", fontsize=9.5, **L["font"])
    ax.axvspan(*SE_BAND_UM, color=DANGER, alpha=0.10, lw=0, zorder=0)
    ax.text(np.sqrt(SE_BAND_UM[0] * SE_BAND_UM[1]), len(items) - 0.32, L["band"],
            ha="center", fontsize=9, color=DANGER, **L["font"])
    ax.set_xscale("log")
    ax.set_xlim(0.05, 200)
    ax.set_ylim(-0.7, len(items) - 0.1)
    ax.set_yticks([])
    ax.set_xlabel(L["wl_um"], fontsize=10.5, **L["font"])
    ax.grid(axis="x", alpha=0.3, which="both")
    ax.text(0.99, 0.02, L["openend"], transform=ax.transAxes, ha="right",
            fontsize=8, color="#666", **L["font"])
    ax.set_title(L["fig1"], fontsize=12.5, pad=10, **L["font"])
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig1-wavelength-coverage.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 편광자 세 가지
# ---------------------------------------------------------------------------

def figure2(L):
    """글랜형 프리즘은 직육면체를 대각으로 잘라 사이에 얇은 층을 둔 구조다."""
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.1))
    specs = [(L["gt"], L["gt_note"], L["air"], "#eaf4ff", L["calcite"], "#2a9d4a"),
             (L["gth"], L["gth_note"], L["glue"], "#fff0f0", L["calcite"], DANGER),
             (L["roch"], L["roch_note"], None, "#eefaf0", L["quartz"], None)]
    # 대각면 양끝
    a, b = np.array([0.38, 0.30]), np.array([0.70, 0.74])
    for ax, (title, note, gap, fc, mat, gapcol) in zip(axes, specs):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(Polygon([(0.16, 0.30), a, b, (0.16, 0.74)],
                             fc=fc, ec="#444", lw=1.2, zorder=2))
        ax.add_patch(Polygon([a + [0.03, 0], (0.86, 0.30), (0.86, 0.74), b + [0.03, 0]],
                             fc=fc, ec="#444", lw=1.2, zorder=2))
        ax.text(0.22, 0.36, mat, fontsize=8.5, zorder=4, **L["font"])
        if gap is not None:
            ax.plot([a[0] + 0.015, b[0] + 0.015], [a[1], b[1]], color=gapcol,
                    lw=3.0, solid_capstyle="butt", zorder=3)
            ax.text(b[0] + 0.02, 0.80, gap, ha="center", fontsize=9.5,
                    color=gapcol, zorder=4, **L["font"])
        # 광선: y=0.52 에서 대각면과 만난다
        t = (0.52 - a[1]) / (b[1] - a[1])
        hit = a + t * (b - a)
        ax.add_patch(FancyArrowPatch((0.02, 0.52), (0.16, 0.52), arrowstyle="-|>",
                                     mutation_scale=13, color="tab:blue", lw=2.0, zorder=5))
        ax.add_patch(FancyArrowPatch((0.86, 0.52), (0.98, 0.52), arrowstyle="-|>",
                                     mutation_scale=13, color="tab:blue", lw=2.0, zorder=5))
        ax.plot([0.16, 0.86], [0.52, 0.52], color="tab:blue",
                lw=1.8, ls=":" if title == L["roch"] else "-", zorder=5)
        if title != L["roch"]:
            ax.add_patch(FancyArrowPatch(tuple(hit), (hit[0] - 0.13, 0.26),
                                         arrowstyle="-|>", mutation_scale=11,
                                         color="#999", lw=1.5, zorder=5))
            ax.text(hit[0] - 0.15, 0.19, L["reject"], fontsize=8.2, color="#888",
                    ha="center", zorder=5, **L["font"])
        ax.set_title(title, fontsize=11.5, pad=8, **L["font"])
        ax.text(0.5, 0.015, note, ha="center", va="bottom", fontsize=8.8,
                color="#444", linespacing=1.45, **L["font"])
    fig.suptitle(L["fig2"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(L["dir"], "fig2-polarizer-types.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — 지연량의 파장 의존성
# ---------------------------------------------------------------------------

def figure3(L):
    nm = WL * 1e9
    simple = np.rad2deg(retardance(WL, plate_thickness(0.25, 550e-9)))
    two, axis = two_element_curves()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.4, 7.0), sharex=True,
                                   gridspec_kw={"height_ratios": [2.0, 1.0]})
    for v in (0, 180, 360):
        ax1.axhline(v, color=DANGER, lw=1.1, ls="--", alpha=0.75, zorder=1)
    ax1.axhspan(170, 190, color=DANGER, alpha=0.10, lw=0, zorder=0)
    ax1.plot(nm, simple, color=C_SIMPLE, lw=2.2, label=L["simple"], zorder=3)
    ax1.plot(nm, two, color=C_TWO, lw=2.2, label=L["two"], zorder=3)

    i180 = int(np.argmin(np.abs(simple - 180)))
    ax1.annotate(L["cross"], xy=(nm[i180], 180), xytext=(nm[i180] + 150, 268),
                 fontsize=9, color=C_SIMPLE, ha="left",
                 arrowprops=dict(arrowstyle="->", color=C_SIMPLE, lw=1.3), **L["font"])
    j = int(np.argmin(np.abs(nm - 820)))
    ax1.annotate(L["stay"], xy=(nm[j], two[j]), xytext=(600, 140), fontsize=9.2,
                 color=C_TWO, ha="left",
                 arrowprops=dict(arrowstyle="->", color=C_TWO, lw=1.3), **L["font"])
    ax1.set_ylim(0, 380)
    ax1.set_yticks(np.arange(0, 361, 90))
    ax1.set_ylabel(L["ret_deg"], fontsize=10.5, **L["font"])
    ax1.set_title(L["fig3a"], fontsize=11.5, loc="left", **L["font"])
    ax1.plot([], [], color=DANGER, lw=1.1, ls="--", label=L["fatal"])
    ax1.legend(prop=L["legend"], fontsize=9.0, loc="lower left", framealpha=0.95)
    ax1.grid(alpha=0.25)

    ax2.plot(nm, axis, color=C_TWO, lw=2.2)
    ax2.set_ylabel(L["axis_deg"], fontsize=10.5, **L["font"])
    ax2.set_xlabel(L["wl_nm"], fontsize=10.5, **L["font"])
    ax2.set_title(L["fig3b"], fontsize=11.5, loc="left", **L["font"])
    ax2.text(0.98, 0.88, L["axis_swing"].format(axis.max() - axis.min()),
             transform=ax2.transAxes, ha="right", va="top", fontsize=9,
             color=C_TWO, linespacing=1.4, **L["font"])
    ax2.grid(alpha=0.25)
    ax2.set_xlim(nm[0], nm[-1])

    fig.suptitle(L["fig3"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig3-retardance-vs-wavelength.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — 영차 vs 다중차
# ---------------------------------------------------------------------------

def figure4(L):
    """지연량은 2pi 주기이므로 절대값과 360도 나머지를 함께 본다."""
    nm = WL * 1e9
    d0, dm = plate_thickness(0.25, 550e-9), plate_thickness(2.25, 550e-9)
    r0, rm = np.rad2deg(retardance(WL, d0)), np.rad2deg(retardance(WL, dm))

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9.4, 7.2), sharex=True,
                                   gridspec_kw={"height_ratios": [1.0, 1.25]})
    ax1.plot(nm, r0, color=C_SIMPLE, lw=2.4, label=L["zero"].format(d0 * 1e6))
    ax1.plot(nm, rm, color=C_MULTI, lw=2.2, label=L["multi"].format(dm * 1e6))
    ax1.set_yscale("log")
    ax1.set_ylabel(L["ret_deg"], fontsize=10.5, **L["font"])
    ax1.set_title(L["fig4a"], fontsize=11.5, loc="left", **L["font"])
    ax1.legend(prop=L["legend"], fontsize=9.2, loc="upper right")
    ax1.grid(alpha=0.25, which="both")

    for v in (0, 180, 360):
        ax2.axhline(v, color=DANGER, lw=1.1, ls="--", alpha=0.8, zorder=1)
    for r, col, lw in ((r0, C_SIMPLE, 2.4), (rm, C_MULTI, 1.8)):
        m = np.mod(r, 360.0)
        m[np.abs(np.diff(m, prepend=m[0])) > 180] = np.nan   # 되돌아오는 선 끊기
        ax2.plot(nm, m, color=col, lw=lw, zorder=3)
    n_cross = int(np.sum(np.abs(np.diff(np.mod(rm, 360.0))) > 180))
    ax2.text(0.985, 0.94, L["multi_note"].format(n_cross), transform=ax2.transAxes,
             ha="right", va="top", fontsize=9.2, color=C_MULTI,
             linespacing=1.45, **L["font"])
    ax2.axvline(550, color="#666", lw=1.0, ls=":")
    ax2.text(556, 340, L["design"], fontsize=9, color="#666", **L["font"])
    ax2.set_ylim(-15, 375)
    ax2.set_yticks([0, 90, 180, 270, 360])
    ax2.set_xlim(nm[0], nm[-1])
    ax2.set_xlabel(L["wl_nm"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["ret_mod"], fontsize=10.5, **L["font"])
    ax2.set_title(L["fig4b"], fontsize=11.5, loc="left", **L["font"])
    ax2.grid(alpha=0.25)

    fig.suptitle(L["fig4"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig4-zero-vs-multi-order.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 5 — 프레넬 롬
# ---------------------------------------------------------------------------

def figure5(L):
    fig, ax = plt.subplots(figsize=(9.4, 4.6))
    for n, style in zip((1.50, 1.51, 1.52), ("-", "-", "-")):
        tc = critical_angle(n)
        th = np.linspace(tc + 1e-9, np.pi / 2 - 1e-9, 4000)
        d = np.abs(np.rad2deg(tir_phase_difference(th, n)))
        ax.plot(np.rad2deg(th), d, style, lw=2.0, alpha=0.85,
                label=L["nlabel"].format(n))
        if n == 1.51:
            ax.axvline(np.rad2deg(tc), color="#888", lw=1.0, ls=":")
            ax.text(np.rad2deg(tc) + 0.5, 3, L["crit"], fontsize=9, color="#666",
                    **L["font"])
            hits = np.rad2deg(th[np.where(np.diff(np.sign(d - 22.5)))[0]])
            ax.plot(hits, np.full_like(hits, 22.5), "o", color=DANGER, ms=7, zorder=5)
            for h in hits:
                ax.annotate(f"{h:.1f}°", xy=(h, 22.5), xytext=(h, 10.5),
                            fontsize=9, ha="center", color=DANGER,
                            arrowprops=dict(arrowstyle="->", color=DANGER, lw=1.0))
    ax.axhline(22.5, color=DANGER, lw=1.3, ls="--")
    ax.text(46, 24.0, L["target"], ha="left", fontsize=9.2, color=DANGER, **L["font"])
    ax.text(0.985, 0.66, L["achro"], transform=ax.transAxes, ha="right", va="top",
            fontsize=9.2, color="#333", linespacing=1.4, **L["font"])
    ax.set_xlim(38, 90)
    ax.set_ylim(0, 50)
    ax.set_xlabel(L["aoi"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["ps_deg"], fontsize=10.5, **L["font"])
    ax.set_title(L["fig5"], fontsize=12.5, pad=10, **L["font"])
    ax.legend(prop=L["legend"], fontsize=9.5, loc="upper right")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig5-fresnel-rhomb.png"), dpi=150)
    plt.close(fig)


def main():
    two, axis = two_element_curves()
    d0 = plate_thickness(0.25, 550e-9) * 1e6
    dm = plate_thickness(2.25, 550e-9) * 1e6
    simple = np.rad2deg(retardance(WL, plate_thickness(0.25, 550e-9)))
    print(f"석영 dn = {DN_QUARTZ}")
    print(f"영차 lambda/4 @550nm 두께 {d0:.2f} um, 다중차(N=2.25) {dm:.2f} um")
    print(f"단순 영차판 지연량 범위 (200-1000nm): {simple.min():.1f} ~ {simple.max():.1f}°"
          f"  -> 180° 통과 파장 {np.interp(180, simple[::-1], (WL*1e9)[::-1]):.0f} nm")
    print(f"두 소자 보상자 지연량 {two.min():.1f} ~ {two.max():.1f}°, "
          f"유효 빠른 축 변동 {axis.max()-axis.min():.1f}°")
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        figure1(L); figure2(L); figure3(L); figure4(L); figure5(L)
        print(f"  [{lang}] 그림 5개 저장 -> {os.path.normpath(L['dir'])}")


if __name__ == "__main__":
    main()

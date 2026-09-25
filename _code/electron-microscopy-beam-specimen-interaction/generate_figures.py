"""전자현미경 배경이론 3편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py   (몬테카를로 때문에 1~2분 걸린다)
출력: ../../assets/img/posts/electron-microscopy-beam-specimen-interaction/      (한국어)
      ../../assets/img/posts/electron-microscopy-beam-specimen-interaction/en/   (영문)

산란 계산은 monte_carlo.py 에 있다. 이 파일은 그 결과를 그리기만 한다.
난수 씨앗을 고정했으므로 두 언어의 그림은 같은 데이터를 쓰고 라벨만 다르다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import monte_carlo as mc

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "electron-microscopy-beam-specimen-interaction",
)

# 그림1: 물질 2종 × 가속전압 2종
PANELS = [(mc.SI, 5.0), (mc.SI, 30.0), (mc.W, 5.0), (mc.W, 30.0)]
N_TRAJ = 220
N_ETA = 3000       # eta 통계용 표본 (궤적은 그리지 않는다)

# 그림2: 깊이별 분포를 볼 조건
FIG2_MAT, FIG2_E0, FIG2_N = mc.SI, 15.0, 1200
SI_K_EDGE = 1.839            # keV, Si K 각질 이온화 에너지
SE_ESCAPE_NM = 5.0           # nm, 이차전자 탈출 깊이의 대표값

# 그림3: 후방산란계수를 잴 원소
FIG3_MATS = [mc.C, mc.SI, mc.CU, mc.AG, mc.W, mc.AU]
FIG3_N = 3000
FIG3_E0 = 20.0

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREEN = "#27795b"
GRAY = "#7f8c8d"

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "fig1": "가속전압과 원자번호가 바꾸는 상호작용 부피",
        "f1title": "{} · {:.0f} kV",
        "f1x": "가로 위치 ({})", "f1y": "깊이 ({})",
        "f1abs": "흡수", "f1bse": "후방산란", "f1ko": "K-O 침투깊이 {:.3g} {}",
        "f1eta": "후방산란계수 η = {:.2f}",
        "f1note": "네 패널의 축 눈금이 서로 다르다 — 각 패널 안의 숫자로 크기를 비교할 것",
        "fig2": "같은 빔이 만드는 신호들이 서로 다른 깊이에서 나온다",
        "f2a": "깊이별 에너지 침착과 신호 발생",
        "f2b": "표면 20 nm 확대 — 이차전자가 빠져나오는 층",
        "f2x": "깊이 ($\\mu$m)", "f2x2": "깊이 (nm)", "f2y": "상대 빈도",
        "f2dep": "에너지 침착 (이차전자 생성에 비례)",
        "f2xray": "특성 X선 발생 (E > 1.84 keV)",
        "f2bse": "후방산란 전자가 도달했던 최대 깊이",
        "f2se": "이차전자 탈출 깊이 ~5 nm",
        "f2ratio": "K-O 침투깊이 {rng:.2f} $\\mu$m 의 {pct:.3f} %",
        "fig3": "후방산란계수는 원자번호를 따라간다",
        "f3x": "원자번호 Z", "f3y": "후방산란계수 η",
        "f3mc": "몬테카를로 ({:.0f} kV, 전자 {:,}개)", "f3reuter": "Reuter 경험식",
        "f3note": "고Z 에서 몬테카를로가 계통적으로 높다 —\n차폐 러더퍼드 단면적의 한계다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Interaction volume reshaped by accelerating voltage and atomic number",
        "f1title": "{} · {:.0f} kV",
        "f1x": "Lateral position ({})", "f1y": "Depth ({})",
        "f1abs": "Absorbed", "f1bse": "Backscattered", "f1ko": "K-O range {:.3g} {}",
        "f1eta": "Backscatter coefficient η = {:.2f}",
        "f1note": "Axis scales differ between panels — compare using the numbers inside each",
        "fig2": "One beam, several signals, each from a different depth",
        "f2a": "Energy deposition and signal generation with depth",
        "f2b": "Top 20 nm — the layer secondary electrons escape from",
        "f2x": "Depth ($\\mu$m)", "f2x2": "Depth (nm)", "f2y": "Relative frequency",
        "f2dep": "Energy deposition (proportional to SE generation)",
        "f2xray": "Characteristic X-ray generation (E > 1.84 keV)",
        "f2bse": "Maximum depth reached by backscattered electrons",
        "f2se": "SE escape depth ~5 nm",
        "f2ratio": "{pct:.3f} % of the {rng:.2f} $\\mu$m K-O range",
        "fig3": "The backscatter coefficient tracks atomic number",
        "f3x": "Atomic number Z", "f3y": "Backscatter coefficient η",
        "f3mc": "Monte Carlo ({:.0f} kV, {:,} electrons)", "f3reuter": "Reuter empirical fit",
        "f3note": "Monte Carlo runs systematically high at high Z —\na limit of the screened Rutherford cross section",
    },
}


def _scaled(value_cm):
    """cm 값을 보기 좋은 단위로 바꾼다. (배율, 단위이름)"""
    return (1e4, "$\\mu$m") if value_cm >= 1e-4 else (1e7, "nm")


def fig1_interaction_volume(L, cache):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
    for ax, (mat, e0) in zip(axes.ravel(), PANELS):
        res = mc.simulate(mat, e0, N_TRAJ, seed=7)
        r_ko = mc.kanaya_okayama_range(mat, e0)
        k, unit = _scaled(r_ko)

        for r in res:
            p = r["path"]
            color = ACCENT if r["backscattered"] else BLUE
            ax.plot(p[:, 0] * k, p[:, 2] * k, color=color, lw=0.35,
                    alpha=0.55 if r["backscattered"] else 0.25, zorder=3)

        lim = r_ko * k * 1.15
        theta = np.linspace(0, np.pi, 200)
        ax.plot(r_ko * k * np.cos(theta), r_ko * k * np.sin(theta),
                color="#2c3e50", ls="--", lw=1.3, zorder=4)
        # y 축이 뒤집혀 있으므로 va="top" 이라야 호 아래(화면 기준)로 간다.
        ax.text(0, r_ko * k * 1.04, L["f1ko"].format(r_ko * k, unit),
                ha="center", va="top", fontsize=8.5, color="#2c3e50", **L["font"])

        eta = cache["eta"][(mat.name, e0)]
        ax.text(0.03, 0.05, L["f1eta"].format(eta), transform=ax.transAxes,
                fontsize=9, color=ACCENT, **L["font"])

        ax.axhline(0, color="#2c3e50", lw=1.2)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(lim * 1.12, -lim * 0.12)          # 깊이는 아래로 증가
        ax.set_aspect("equal")
        ax.set_title(L["f1title"].format(mat.name, e0), **L["font"])
        ax.set_xlabel(L["f1x"].format(unit), **L["font"])
        ax.set_ylabel(L["f1y"].format(unit), **L["font"])
        ax.grid(alpha=0.25)

    handles = [plt.Line2D([], [], color=BLUE, lw=1.6, label=L["f1abs"]),
               plt.Line2D([], [], color=ACCENT, lw=1.6, label=L["f1bse"])]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 0.975),
               prop=L["legend"] or None, fontsize=9)
    fig.suptitle(L["fig1"], fontsize=13, **L["font"])
    fig.text(0.5, 0.005, L["f1note"], ha="center", fontsize=8.5, color=GRAY, **L["font"])
    fig.tight_layout(rect=(0, 0.022, 1, 0.955))
    return fig


def fig2_depth_distribution(L, cache):
    res = cache["fig2"]
    r_ko_um = mc.kanaya_okayama_range(FIG2_MAT, FIG2_E0) * 1e4

    dep_z = np.concatenate([r["depths"] for r in res]) * 1e4      # µm
    dep_w = np.concatenate([r["losses"] for r in res])
    ene = np.concatenate([r["energies"] for r in res])
    xray = dep_z[ene > SI_K_EDGE]
    xray_w = dep_w[ene > SI_K_EDGE]
    bse_depth = np.array([r["max_depth"] for r in res if r["backscattered"]]) * 1e4

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.6))

    fine = np.linspace(0, r_ko_um * 1.05, 90)
    coarse = np.linspace(0, r_ko_um * 1.05, 32)     # 표본이 적은 후방산란용
    for data, weights, color, key, ls, bins in [
            (dep_z, dep_w, BLUE, "f2dep", "-", fine),
            (xray, xray_w, GREEN, "f2xray", "--", fine),
            (bse_depth, None, ACCENT, "f2bse", ":", coarse)]:
        h, edges = np.histogram(data, bins=bins, weights=weights)
        if h.max() > 0:
            h = h / h.max()
        ax.plot(0.5 * (edges[1:] + edges[:-1]), h, color=color, lw=1.9, ls=ls,
                label=L[key])

    se_um = SE_ESCAPE_NM * 1e-3
    ax.axvspan(0, se_um, color="#f39c12", alpha=0.45, zorder=5)
    ax.annotate(L["f2se"], xy=(se_um, 0.10), xytext=(r_ko_um * 0.33, 0.09),
                fontsize=9, color="#b9770e", **L["font"],
                arrowprops=dict(arrowstyle="->", color="#b9770e"))
    ax.set_xlabel(L["f2x"], **L["font"])
    ax.set_ylabel(L["f2y"], **L["font"])
    ax.set_title(L["f2a"], **L["font"])
    ax.set_xlim(0, r_ko_um * 1.05)
    ax.set_ylim(0, 1.15)
    ax.grid(alpha=0.3)
    ax.legend(prop=L["legend"] or None, fontsize=8.5, loc="upper right",
              framealpha=0.95)

    # 오른쪽: 표면 20 nm 확대
    bins2 = np.linspace(0, 20e-3, 11)
    h, edges = np.histogram(dep_z, bins=bins2, weights=dep_w)
    centers = 0.5 * (edges[1:] + edges[:-1]) * 1e3            # nm
    ax2.bar(centers, h / h.max(), width=(centers[1] - centers[0]) * 0.9,
            color=BLUE, alpha=0.75, zorder=3)
    ax2.axvspan(0, SE_ESCAPE_NM, color="#f39c12", alpha=0.45, zorder=4)
    ax2.set_ylim(0, 1.35)
    ax2.text(SE_ESCAPE_NM * 1.3, 1.24, L["f2se"], fontsize=9.5, color="#b9770e",
             **L["font"])
    ax2.text(SE_ESCAPE_NM * 1.3, 1.12,
             L["f2ratio"].format(rng=r_ko_um,
                                 pct=SE_ESCAPE_NM * 1e-3 / r_ko_um * 100),
             fontsize=9, color=GRAY, **L["font"])
    ax2.set_xlabel(L["f2x2"], **L["font"])
    ax2.set_ylabel(L["f2y"], **L["font"])
    ax2.set_title(L["f2b"], **L["font"])
    ax2.set_xlim(0, 20)
    ax2.grid(alpha=0.3, axis="y")

    fig.suptitle(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def fig3_backscatter_vs_z(L, cache):
    zs = np.array([m.Z for m in FIG3_MATS])
    etas = cache["fig3"]

    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    zline = np.linspace(4, 82, 300)
    ax.plot(zline, mc.reuter_eta(zline), color=GRAY, lw=2.0, ls="--",
            label=L["f3reuter"])
    ax.plot(zs, etas, "o", color=ACCENT, ms=8, zorder=5,
            label=L["f3mc"].format(FIG3_E0, FIG3_N))
    for m, z, e in zip(FIG3_MATS, zs, etas):
        ax.annotate(m.name, xy=(z, e), xytext=(0, 11), textcoords="offset points",
                    ha="center", fontsize=9.5, color="#2c3e50")

    ax.set_xlabel(L["f3x"], **L["font"])
    ax.set_ylabel(L["f3y"], **L["font"])
    ax.set_title(L["fig3"], fontsize=13, **L["font"])
    ax.set_xlim(0, 86)
    ax.set_ylim(0, 0.62)
    ax.grid(alpha=0.3)
    ax.legend(prop=L["legend"] or None, loc="upper left")
    ax.text(0.97, 0.06, L["f3note"], transform=ax.transAxes, fontsize=8.5,
            color=GRAY, ha="right", **L["font"])
    fig.tight_layout()
    return fig


def main():
    # 무거운 계산은 한 번만 하고 두 언어가 나눠 쓴다.
    print("몬테카를로 계산 중...")
    cache = {"eta": {(m.name, e): mc.backscatter_coefficient(m, e, N_ETA, seed=5)
                     for m, e in PANELS}}
    cache["fig2"] = (mc.simulate(FIG2_MAT, FIG2_E0, FIG2_N, seed=11, keep_path=False))
    cache["fig3"] = [mc.backscatter_coefficient(m, FIG3_E0, FIG3_N, seed=3)
                     for m in FIG3_MATS]

    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for name, builder in [
                ("fig1-interaction-volume", lambda l: fig1_interaction_volume(l, cache)),
                ("fig2-depth-distribution", lambda l: fig2_depth_distribution(l, cache)),
                ("fig3-backscatter-vs-z", lambda l: fig3_backscatter_vs_z(l, cache))]:
            fig = builder(L)
            fig.savefig(os.path.join(L["dir"], name + ".png"), dpi=150,
                        bbox_inches="tight")
            plt.close(fig)
        print(f"[{lang}] 3 figures written to {L['dir']}")

    # 본문에 인용한 수치를 검산할 수 있게 출력한다.
    print()
    for mat, e0 in PANELS:
        print(f"{mat.name:>3} {e0:4.0f} kV : K-O = "
              f"{mc.kanaya_okayama_range(mat, e0)*1e4:8.4f} µm")
    print()
    for m, e in zip(FIG3_MATS, cache["fig3"]):
        print(f"{m.name:>3} (Z={m.Z:.0f}) : MC eta = {e:.3f}, "
              f"Reuter = {mc.reuter_eta(m.Z):.3f}")
    r_ko = mc.kanaya_okayama_range(FIG2_MAT, FIG2_E0) * 1e4
    print(f"\nSi {FIG2_E0:.0f} kV : K-O = {r_ko:.3f} µm, "
          f"SE 탈출깊이 {SE_ESCAPE_NM:.0f} nm = 그 {SE_ESCAPE_NM*1e-3/r_ko*100:.3f} %")


if __name__ == "__main__":
    main()

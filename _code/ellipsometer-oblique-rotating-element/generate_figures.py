"""타원계측기 1편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-oblique-rotating-element/      (한국어)
      ../../assets/img/posts/ellipsometer-oblique-rotating-element/en/   (영문)

계산에 무작위 요소가 없으므로 두 언어의 그림은 데이터가 완전히 동일하고 표기만 다르다.
그림에 쓰는 I(t) 는 전부 modulation.py 의 뮬러 곱 계산에서 나오며,
그 계산은 verify_modulation.py 가 Fujiwara 해석식과 대조해 검증한 것이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Polygon
from scipy.special import jv

from modulation import (amplitude_spectrum, intensity_pme, intensity_rae,
                        intensity_rce, normalized_coefficients, pem_retardation)
from smm_tensor import smm_reflectance_tensor

# 한글 텍스트에만 한글 폰트를 지정한다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-oblique-rotating-element",
)

# 그림에 쓰는 기준 시편: SiO2 100 nm / Si, 입사각 75°, 632.8 nm
N_AIR, N_SIO2, N_SI = 1.0, 1.457, 3.882 - 0.019j
FILM_NM, AOI_DEG, WAVELENGTH_NM = 100.0, 75.0, 632.8
PEM_F = np.deg2rad(138.0)

C_RAE, C_RCE, C_PME = "tab:blue", "tab:red", "tab:green"


def sample_psi_delta():
    """기준 시편의 (Psi, Delta) 를 3편의 SMM 구현으로 계산한다(rad)."""
    wl, aoi = np.array([WAVELENGTH_NM]), np.deg2rad(np.array([AOI_DEG]))
    layers = [N_AIR, N_SIO2, N_SI]
    rp = smm_reflectance_tensor(layers, [FILM_NM], aoi, wl, pol="p")
    rs = smm_reflectance_tensor(layers, [FILM_NM], aoi, wl, pol="s")
    rho = (rp / rs).ravel()[0]
    return np.arctan(np.abs(rho)), np.angle(rho)


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림 1
        "fig1": "경사입반사 구조와 PSG–시편–PSA 사슬",
        "source": "광원", "detector": "검출기",
        "psg": "편광 생성단 (PSG)", "psa": "편광 분석단 (PSA)",
        "psg_sub": "편광자 (+ 보상자)", "psa_sub": "(보상자 +) 분석기",
        "sample": "시편", "known": "아는 것", "unknown": "모르는 것",
        "known_txt": "들여보낸 편광 상태", "meas_txt": "검출 세기 $I$ 하나뿐",
        "sample_txt": "3편의 $M_{sample}(\\Psi, \\Delta)$",
        "plane": "입사면", "film": "박막", "substrate": "기판",
        # 그림 2
        "fig2": "네 가지 광학 배치 — 아래첨자 R 이 변조를 담당하는 소자다",
        "pol": "편광자\nP", "ana": "분석기\nA", "comp": "보상자\nC", "pem": "위상변조기\nM",
        "smp": "시편\nS",
        "rae": "RAE  $PSA_R$", "raec": "RAE+보상자  $PSCA_R$",
        "rce": "RCE  $PSC_RA$", "pme": "PME  $PSMA$",
        "rae_note": "분석기를 돌린다\n$S_3$ 를 못 잰다",
        "raec_note": "보상자는 고정, 분석기가 회전\n$\\delta$ 를 바꿔 두 번 재야 한다",
        "rce_note": "보상자를 돌린다\n단일 측정으로 $S_1$–$S_3$",
        "pme_note": "회전 부품이 없다\n지연량을 전기로 흔든다",
        "rot": "회전", "mod": "전기 변조",
        # 그림 3
        "fig3": "RAE 는 $\\Delta$ 의 부호를 구분하지 못하고, RCE 는 구분한다",
        "fig3a": "(a) RAE — 두 곡선이 완전히 겹친다",
        "fig3b": "(b) RCE — $2\\omega$ 성분의 부호가 갈린다",
        "ana_angle": "분석기 방위각 $A$ (deg)",
        "comp_angle": "보상자 방위각 $C$ (deg)",
        "norm_int": "정규화 검출 세기 $I/I_0$",
        "dpos": "$\\Delta = +{:.1f}°$", "dneg": "$\\Delta = -{:.1f}°$",
        "identical": "두 파형의 최대 차이\n{:.1e} (수치오차 수준)",
        "beta2": "$\\beta_2$ 부호 반전",
        # 그림 4
        "fig4": "같은 시편을 세 방식으로 측정했을 때의 검출 세기와 고조파 구성",
        "fig4sub": "SiO$_2$ {:.0f} nm / Si, 입사각 {:.0f}°, $\\Psi$={:.1f}°, $\\Delta$={:.1f}°",
        "cycle": "변조 1주기  $\\omega t$ (deg)",
        "harmonic": "고조파 차수 $n$ (변조 주파수 $\\omega$ 의 배수)",
        "coef": "정규화 푸리에 계수 크기",
        "rae_t": "RAE — 분석기 회전", "rce_t": "RCE — 보상자 회전",
        "pme_t": "PME — 광탄성 변조",
        "only2w": "$2\\omega$ 하나", "two24": "$2\\omega$ 와 $4\\omega$",
        "onetwo": "$\\omega$ 와 $2\\omega$",
        # 그림 5
        "fig5": "$\\Delta$ 측정 오차의 증폭률과 보상자를 이용한 회피",
        "fig5a": "(a) RAE — $\\Delta \\cong 0°, 180°$ 에서 발산한다",
        "fig5b": "(b) 보상자로 작업점을 옮겨 취약 구간을 메운다",
        "delta_ax": "시편의 $\\Delta$ (deg)",
        "amp_ax": "오차 증폭률  $1/|\\sin\\Delta'|$",
        "rae_curve": "보상자 없음 ($\\delta = 0°$)",
        "dlabel": "$\\delta = {:.0f}°$",
        "envelope": "네 측정 중 가장 좋은 것만 고르면",
        "danger": "취약 구간",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Oblique-incidence geometry and the PSG–sample–PSA chain",
        "source": "source", "detector": "detector",
        "psg": "Polarization State Generator", "psa": "Polarization State Analyzer",
        "psg_sub": "polarizer (+ compensator)", "psa_sub": "(compensator +) analyzer",
        "sample": "sample", "known": "known", "unknown": "unknown",
        "known_txt": "the launched polarization state", "meas_txt": "one intensity $I$, nothing else",
        "sample_txt": "$M_{sample}(\\Psi, \\Delta)$ from Part 3",
        "plane": "plane of incidence", "film": "film", "substrate": "substrate",
        "fig2": "Four optical configurations — subscript R marks the modulating element",
        "pol": "polarizer\nP", "ana": "analyzer\nA", "comp": "compensator\nC",
        "pem": "modulator\nM", "smp": "sample\nS",
        "rae": "RAE  $PSA_R$", "raec": "RAE+compensator  $PSCA_R$",
        "rce": "RCE  $PSC_RA$", "pme": "PME  $PSMA$",
        "rae_note": "rotate the analyzer\n$S_3$ is not measurable",
        "raec_note": "fixed compensator, rotating analyzer\ntwo runs at different $\\delta$",
        "rce_note": "rotate the compensator\n$S_1$–$S_3$ in a single run",
        "pme_note": "nothing rotates\nretardation driven electrically",
        "rot": "rotation", "mod": "electrical",
        "fig3": "RAE cannot tell the sign of $\\Delta$; RCE can",
        "fig3a": "(a) RAE — the two curves coincide exactly",
        "fig3b": "(b) RCE — the $2\\omega$ component flips sign",
        "ana_angle": "analyzer azimuth $A$ (deg)",
        "comp_angle": "compensator azimuth $C$ (deg)",
        "norm_int": "normalized intensity $I/I_0$",
        "dpos": "$\\Delta = +{:.1f}°$", "dneg": "$\\Delta = -{:.1f}°$",
        "identical": "max difference between\nthe two: {:.1e}",
        "beta2": "$\\beta_2$ flips sign",
        "fig4": "Detected intensity and harmonic content for one sample, three methods",
        "fig4sub": "SiO$_2$ {:.0f} nm / Si, AOI {:.0f}°, $\\Psi$={:.1f}°, $\\Delta$={:.1f}°",
        "cycle": "one modulation cycle  $\\omega t$ (deg)",
        "harmonic": "harmonic order $n$ (multiple of $\\omega$)",
        "coef": "normalized Fourier coefficient",
        "rae_t": "RAE — rotating analyzer", "rce_t": "RCE — rotating compensator",
        "pme_t": "PME — photoelastic modulation",
        "only2w": "$2\\omega$ only", "two24": "$2\\omega$ and $4\\omega$",
        "onetwo": "$\\omega$ and $2\\omega$",
        "fig5": "Error amplification in $\\Delta$ and how a compensator removes it",
        "fig5a": "(a) RAE — diverges at $\\Delta \\cong 0°, 180°$",
        "fig5b": "(b) shifting the operating point with a compensator",
        "delta_ax": "sample $\\Delta$ (deg)",
        "amp_ax": "error amplification  $1/|\\sin\\Delta'|$",
        "rae_curve": "no compensator ($\\delta = 0°$)",
        "dlabel": "$\\delta = {:.0f}°$",
        "envelope": "best of the four measurements",
        "danger": "weak region",
    },
}


# ---------------------------------------------------------------------------
# 그림 1 — 경사입반사 구조
# ---------------------------------------------------------------------------

def figure1(L):
    """입사각이 실제 각도로 보이도록 aspect=equal 로 그린다."""
    fig, ax = plt.subplots(figsize=(9.2, 5.2))
    ax.set_aspect("equal")
    ax.set_xlim(-0.02, 2.22)
    ax.set_ylim(-0.50, 0.80)
    ax.axis("off")

    hit = np.array([1.10, 0.0])
    rise, theta_deg = 0.45, 60.0
    run = rise * np.tan(np.deg2rad(theta_deg))
    p_in = hit + np.array([-run, rise])
    p_out = hit + np.array([run, rise])

    # 시편(박막 + 기판)
    ax.add_patch(plt.Rectangle((0.52, -0.07), 1.16, 0.07, fc="#cfe3f7", ec="#3b6ea5", lw=1.2))
    ax.add_patch(plt.Rectangle((0.52, -0.28), 1.16, 0.21, fc="#dcdcdc", ec="#777", lw=1.2))
    ax.text(1.72, -0.035, L["film"], fontsize=9.5, va="center", **L["font"])
    ax.text(1.72, -0.175, L["substrate"], fontsize=9.5, va="center", **L["font"])

    # 법선과 광선
    ax.plot([hit[0], hit[0]], [0, 0.66], ls=":", color="#888", lw=1.0, zorder=1)
    for a, b, col in [(p_in, hit, "tab:blue"), (hit, p_out, "tab:red")]:
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=16,
                                     color=col, lw=2.0, shrinkA=0, shrinkB=0, zorder=3))

    # 입사각
    arc = np.deg2rad(np.linspace(90, 90 + theta_deg, 60))
    ax.plot(hit[0] + 0.20 * np.cos(arc), 0.20 * np.sin(arc), color="#555", lw=1.0, zorder=2)
    ax.text(hit[0] - 0.135, 0.195, r"$\theta$", fontsize=13, color="#333", zorder=4)

    # p / s 표기 (광선 위 40% 지점)
    for frac, end_pt, col in [(0.45, p_in, "tab:blue"), (0.45, p_out, "tab:red")]:
        q = hit + frac * (end_pt - hit)
        ax.text(q[0], q[1] + 0.045, "p", fontsize=10.5, color=col, weight="bold",
                ha="center", zorder=4)
        ax.text(q[0] + 0.075, q[1] + 0.045, "s", fontsize=10.5, color=col, weight="bold",
                alpha=0.5, ha="center", zorder=4)

    def box(cx, y, w, h, title, sub, fc):
        ax.add_patch(FancyBboxPatch((cx - w / 2, y), w, h, boxstyle="round,pad=0.015",
                                    fc=fc, ec="#444", lw=1.1, zorder=5))
        ax.text(cx, y + h * 0.63, title, ha="center", fontsize=10.5, weight="bold",
                zorder=6, **L["font"])
        ax.text(cx, y + h * 0.22, sub, ha="center", fontsize=8.8, color="#444",
                zorder=6, **L["font"])

    box(p_in[0], rise, 0.64, 0.19, L["psg"], L["psg_sub"], "#eaf4ff")
    box(p_out[0], rise, 0.64, 0.19, L["psa"], L["psa_sub"], "#ffeceb")
    ax.text(p_in[0], 0.70, f"{L['known']}: {L['known_txt']}", ha="center", fontsize=9,
            color="#2a8a4a", **L["font"])
    ax.text(p_out[0], 0.70, f"{L['unknown']}: {L['meas_txt']}", ha="center", fontsize=9,
            color="#c0392b", **L["font"])

    # 사슬 한가운데 앉는 시편 행렬 — 기판 아래로 빼서 광로를 가리지 않게 한다
    ax.add_patch(FancyBboxPatch((0.70, -0.46), 0.80, 0.135, boxstyle="round,pad=0.015",
                                fc="#fff3d6", ec="#b8860b", lw=1.2, zorder=5))
    ax.text(1.10, -0.392, L["sample_txt"], ha="center", fontsize=10, zorder=6, **L["font"])
    ax.add_patch(FancyArrowPatch((1.10, -0.325), (1.10, -0.275), arrowstyle="-|>",
                                 mutation_scale=13, color="#b8860b", lw=1.5, zorder=4))

    ax.set_title(L["fig1"], fontsize=12.5, pad=10, **L["font"])
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig1-oblique-configuration.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 네 가지 광학 배치
# ---------------------------------------------------------------------------

def figure2(L):
    fig, axes = plt.subplots(4, 1, figsize=(9.2, 6.3))
    rows = [
        (L["rae"], L["rae_note"], [("pol", 0), ("smp", 0), ("ana", 1)]),
        (L["raec"], L["raec_note"], [("pol", 0), ("smp", 0), ("comp", 0), ("ana", 1)]),
        (L["rce"], L["rce_note"], [("pol", 0), ("smp", 0), ("comp", 1), ("ana", 0)]),
        (L["pme"], L["pme_note"], [("pol", 0), ("smp", 0), ("pem", 2), ("ana", 0)]),
    ]
    colors = {"pol": "#eaf4ff", "ana": "#ffeceb", "comp": "#e9f7e9",
              "pem": "#f3e9ff", "smp": "#fff3d6"}

    for ax, (title, note, elems) in zip(axes, rows):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        n = len(elems)
        w, gap = 0.112, 0.045
        total = n * w + (n - 1) * gap
        x0 = 0.335
        # 광축
        ax.plot([x0 - 0.045, x0 + total + 0.045], [0.45, 0.45], color="#bbb", lw=1.2, zorder=0)
        for i, (kind, moving) in enumerate(elems):
            x = x0 + i * (w + gap)
            ax.add_patch(FancyBboxPatch((x, 0.24), w, 0.42, boxstyle="round,pad=0.012",
                                        fc=colors[kind], ec="#444", lw=1.1, zorder=2))
            ax.text(x + w / 2, 0.45, L[kind], ha="center", va="center",
                    fontsize=9, zorder=3, **L["font"])
            if moving == 1:      # 기계적 회전
                ax.add_patch(FancyArrowPatch((x + w * 0.22, 0.73), (x + w * 0.78, 0.73),
                                             connectionstyle="arc3,rad=-0.55",
                                             arrowstyle="-|>", mutation_scale=11,
                                             color="#c0392b", lw=1.6))
                ax.text(x + w / 2, 0.86, L["rot"], ha="center", fontsize=8,
                        color="#c0392b", **L["font"])
            elif moving == 2:    # 전기적 변조
                tt = np.linspace(0, 2 * np.pi, 120)
                ax.plot(x + w * (0.18 + 0.64 * tt / (2 * np.pi)), 0.76 + 0.05 * np.sin(2 * tt),
                        color="#7d3c98", lw=1.5)
                ax.text(x + w / 2, 0.88, L["mod"], ha="center", fontsize=8,
                        color="#7d3c98", **L["font"])
        ax.text(0.012, 0.66, title, fontsize=10.8, weight="bold", **L["font"])
        ax.text(0.012, 0.47, note, fontsize=8.3, color="#555", va="top",
                linespacing=1.5, **L["font"])

    axes[0].set_title(L["fig2"], fontsize=12.5, pad=12, **L["font"])
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig2-instrument-configurations.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — Delta 부호 모호성
# ---------------------------------------------------------------------------

def figure3(L, psi, delta):
    ang = np.linspace(0, np.pi, 721)
    deg = np.rad2deg(ang)
    d_deg = abs(np.rad2deg(delta))

    rae_p, rae_m = intensity_rae(psi, delta, ang), intensity_rae(psi, -delta, ang)
    rce_p, rce_m = intensity_rce(psi, delta, ang), intensity_rce(psi, -delta, ang)
    gap = np.max(np.abs(rae_p - rae_m))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.3))
    for ax, (yp, ym), xlabel, title in [
            (ax1, (rae_p, rae_m), L["ana_angle"], L["fig3a"]),
            (ax2, (rce_p, rce_m), L["comp_angle"], L["fig3b"])]:
        norm = yp.mean()
        ax.plot(deg, yp / norm, color="tab:blue", lw=3.2, alpha=0.85,
                label=L["dpos"].format(d_deg))
        ax.plot(deg, ym / norm, color="tab:red", lw=1.6, ls="--",
                label=L["dneg"].format(d_deg))
        ax.set_xlim(0, 180)
        ax.set_xticks(np.arange(0, 181, 45))
        ax.set_xlabel(xlabel, fontsize=10.5, **L["font"])
        ax.set_ylabel(L["norm_int"], fontsize=10.5, **L["font"])
        ax.set_title(title, fontsize=11.5, **L["font"])
        ax.grid(alpha=0.25)
        ax.legend(prop=L["legend"], fontsize=9.5, loc="upper right")

    ax1.text(0.04, 0.06, L["identical"].format(gap), transform=ax1.transAxes,
             fontsize=9, color="#c0392b", va="bottom", **L["font"])
    c_p = normalized_coefficients(intensity_rce(psi, delta, ang[:-1]))[1][1]
    c_m = normalized_coefficients(intensity_rce(psi, -delta, ang[:-1]))[1][1]
    ax2.text(0.04, 0.06, f"{L['beta2']}\n{c_p:+.3f}  /  {c_m:+.3f}",
             transform=ax2.transAxes, fontsize=9, color="#c0392b", va="bottom", **L["font"])

    fig.suptitle(L["fig3"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(L["dir"], "fig3-rae-delta-ambiguity.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 5 — 세 방식의 파형과 고조파
# ---------------------------------------------------------------------------

def figure_waveforms(L, psi, delta):
    n_pts, n_max = 2048, 6
    wt = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    waves = [
        (L["rae_t"], C_RAE, intensity_rae(psi, delta, wt), L["only2w"]),
        (L["rce_t"], C_RCE, intensity_rce(psi, delta, wt), L["two24"]),
        (L["pme_t"], C_PME, intensity_pme(psi, delta, pem_retardation(wt, PEM_F, 1.0)),
         L["onetwo"]),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(10.6, 7.6),
                             gridspec_kw={"width_ratios": [2.1, 1.0]})
    deg = np.rad2deg(wt)
    for (name, col, sig, tag), (axw, axs) in zip(waves, axes):
        axw.plot(deg, sig / sig.mean(), color=col, lw=1.9)
        axw.set_xlim(0, 360)
        axw.set_xticks(np.arange(0, 361, 90))
        axw.set_ylabel(L["norm_int"], fontsize=10, **L["font"])
        axw.set_title(name, fontsize=11, loc="left", **L["font"])
        axw.grid(alpha=0.25)

        spec = amplitude_spectrum(sig, n_max=n_max)
        idx = np.arange(1, n_max + 1)
        axs.bar(idx, spec, color=col, alpha=0.85, width=0.62)
        axs.set_xticks(idx)
        axs.set_ylim(0, max(1.18 * spec.max(), 1e-3))
        axs.set_ylabel(L["coef"], fontsize=9, **L["font"])
        axs.set_title(tag, fontsize=10.5, loc="left", color=col, **L["font"])
        axs.grid(alpha=0.25, axis="y")
        for i, v in zip(idx, spec):
            if v > 0.02 * spec.max():
                axs.text(i, v, f"{v:.2f}", ha="center", va="bottom", fontsize=7.5)

    axes[2][0].set_xlabel(L["cycle"], fontsize=10.5, **L["font"])
    axes[2][1].set_xlabel(L["harmonic"], fontsize=9.5, **L["font"])

    sub = L["fig4sub"].format(FILM_NM, AOI_DEG, np.rad2deg(psi), np.rad2deg(delta))
    fig.suptitle(f"{L['fig4']}\n{sub}", fontsize=12, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(os.path.join(L["dir"], "fig5-waveform-comparison.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — Delta 민감도와 보상자
# ---------------------------------------------------------------------------

def figure_sensitivity(L):
    """RAE 는 beta = sin2Psi cos(Delta') 를 재므로 오차 증폭률이 1/|sin Delta'| 이다."""
    d = np.linspace(-180, 180, 1441)
    clip = 12.0

    def amp(shift_deg):
        s = np.abs(np.sin(np.deg2rad(d - shift_deg)))
        with np.errstate(divide="ignore"):
            return np.where(s > 1 / clip, 1 / np.maximum(s, 1e-12), np.nan)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.3), sharey=True)

    ax1.plot(d, amp(0), color=C_RAE, lw=2.2, label=L["rae_curve"])
    for x in (-180, 0, 180):
        ax1.axvspan(x - 12, x + 12, color="#c0392b", alpha=0.12, lw=0)
    ax1.annotate(L["danger"], xy=(6, 7.0), xytext=(72, 9.4), fontsize=9.5,
                 color="#c0392b", ha="left", va="center",
                 arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.2), **L["font"])
    ax1.set_title(L["fig5a"], fontsize=11.5, **L["font"])
    ax1.legend(prop=L["legend"], fontsize=9.5, loc="upper right")

    shifts = [0, 45, 90, 135]
    curves = []
    for sh, c in zip(shifts, ["#1f77b4", "#2ca02c", "#ff7f0e", "#9467bd"]):
        a = amp(sh)
        curves.append(a)
        ax2.plot(d, a, color=c, lw=1.4, alpha=0.75, label=L["dlabel"].format(sh))
    best = np.nanmin(np.vstack(curves), axis=0)
    ax2.plot(d, best, color="k", lw=2.6, label=L["envelope"])
    ax2.set_title(L["fig5b"], fontsize=11.5, **L["font"])
    ax2.legend(prop=L["legend"], fontsize=8.8, loc="upper right", ncol=2)

    for ax in (ax1, ax2):
        ax.set_xlim(-180, 180)
        ax.set_ylim(0, clip)
        ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_xlabel(L["delta_ax"], fontsize=10.5, **L["font"])
        ax.grid(alpha=0.25)
    ax1.set_ylabel(L["amp_ax"], fontsize=10.5, **L["font"])

    fig.suptitle(L["fig5"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(L["dir"], "fig4-rae-sensitivity.png"), dpi=150)
    plt.close(fig)


def main():
    psi, delta = sample_psi_delta()
    print(f"기준 시편: SiO2 {FILM_NM:.0f} nm / Si, 입사각 {AOI_DEG:.0f}°, "
          f"{WAVELENGTH_NM} nm  ->  Psi = {np.rad2deg(psi):.2f}°, "
          f"Delta = {np.rad2deg(delta):.2f}°")
    print(f"PEM: F = 138°, 2J1 = {2*jv(1, PEM_F):.3f}, 2J2 = {2*jv(2, PEM_F):.3f}, "
          f"2J3 = {2*jv(3, PEM_F):.3f}, 2J4 = {2*jv(4, PEM_F):.3f}")
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        figure1(L)
        figure2(L)
        figure3(L, psi, delta)
        figure_sensitivity(L)      # 본문 그림4
        figure_waveforms(L, psi, delta)  # 본문 그림5
        print(f"  [{lang}] 그림 5개 저장 -> {os.path.normpath(L['dir'])}")


if __name__ == "__main__":
    main()

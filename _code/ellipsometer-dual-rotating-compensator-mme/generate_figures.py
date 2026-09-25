"""타원계측기 3편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-dual-rotating-compensator-mme/      (한국어)
      ../../assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/   (영문)

계산은 mme.py 의 뮬러 곱에서 나오고, 그 구현은 verify_mme.py 가
Collins & Koh 표 1·2 의 닫힌 식과 대조해 검증한 것이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from mme import (RATIO, condition_number, harmonic_magnitudes, highest_harmonic,
                 intensity, optical_cycle)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-dual-rotating-compensator-mme",
)

N_PTS = 4096
DEG = np.arange(60.0, 175.1, 1.0)
C_MME, C_RCE, C_DEAD = "tab:red", "tab:blue", "#c0392b"

# 그림에 쓰는 기준 시편 — 16요소가 모두 비영인 임의 행렬
_rng = np.random.default_rng(0)
SAMPLE = _rng.uniform(-0.4, 0.4, (4, 4))
SAMPLE[0, 0] = 1.0

# Collins & Koh 표 1 의 회전비별 설계값
RATIOS = [(5, 1), (5, 2), (5, 3), (5, 4)]

_SWEEP: dict = {}


def cond_sweep(ratio):
    """지연량에 따른 조건수. 그림 4·5 와 두 언어가 같은 결과를 쓰므로 캐시한다."""
    if ratio not in _SWEEP:
        _SWEEP[ratio] = np.array(
            [condition_number(np.deg2rad(d), np.deg2rad(d), ratio) for d in DEG])
    return _SWEEP[ratio]


def harmonic_owner():
    """각 고조파에 실리는 뮬러 요소와 세기 인자 (Collins & Koh 표 2)."""
    return {
        1:  ("M24 M34", "sin d1 · s2"),
        2:  ("M44", "sin d1 · sin d2"),
        3:  ("M41 M42 M43", "sin d2"),
        4:  ("2x2", "s1 · s2"),
        5:  ("M14 M24 M34", "sin d1"),
        6:  ("M21 M31 2x2", "s2"),
        7:  ("M42 M43", "s1 · sin d2"),
        8:  ("M44", "sin d1 · sin d2"),
        10: ("M12 M13 2x2", "s1"),
        11: ("M24 M34", "sin d1 · s2"),
        13: ("M42 M43", "s1 · sin d2"),
        16: ("2x2", "s1 · s2"),
    }


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림 1
        "fig1": "회전형 배치별로 닿을 수 있는 뮬러 요소",
        "cfg": ["편광자·분석기 동시 회전", "시편 뒤에 보상자 하나",
                "시편 앞에 보상자 하나", "보상자 둘을 다른 속도로"],
        "cfg_note": ["4행·4열을 못 잰다", "처음 세 열", "처음 세 행", "16개 전부"],
        "got": "측정되는 요소", "miss": "닿지 않는 요소",
        # 그림 2
        "fig2": "같은 시편을 RCE와 DRC-MME로 볼 때의 파형과 고조파",
        "cycle": "기본 회전각 $C$ (deg)", "norm_int": "정규화 검출 세기",
        "order": "고조파 차수 $n$ (주파수 $2nC$)", "mag": "정규화 계수 크기",
        "rce_t": "1편의 RCE — 보상자 하나",
        "mme_t": "DRC-MME — 보상자 둘 (5:3)",
        "rce_note": "$2\\omega$ 와 $4\\omega$ 둘뿐",
        "mme_note": "$2C$ 부터 $32C$ 까지 24개",
        "dead": "$n$ = 9, 12, 14, 15 는 소멸한다",
        # 그림 3
        "fig3": "어느 뮬러 요소가 어느 고조파에 실리는가",
        "elem": "실리는 요소", "factor": "세기 인자",
        "fig3note": "4행·4열은 오직 $\\sin\\delta$ 를 탄다 — 지연량이 0이나 180°면 통째로 사라진다",
        # 그림 4
        "fig4": "회전비가 정하는 것은 정밀도가 아니라 읽어야 하는 횟수다",
        "ratio_ax": "회전비 $m_1 : m_2$",
        "high_ax": "최고차 고조파", "nint_ax": "주기당 최소 적분 횟수",
        "fig4a": "(a) 비율을 키우면 고조파가 올라간다",
        "fig4b": "(b) 그런데 조건수는 전부 같다",
        "cond_ax": "데이터 환산 행렬의 조건수",
        "ret_ax": "보상자 지연량 (deg)",
        # 그림 5
        "fig5": "정밀도를 정하는 것은 지연량이다",
        "opt": "최소 {:.0f}°  (조건수 {:.2f})",
        "at90": "90°에서 {:.2f}",
        "gain": "90° 대신 {:.0f}°를 쓰면\n잡음 증폭이 {:.1f}배 줄어든다",
        "overlap": "5:3 과 5:1 이 완전히 겹친다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Mueller elements each rotating-element configuration can reach",
        "cfg": ["polarizer and analyzer rotate together", "one compensator after the sample",
                "one compensator before the sample", "two compensators, different speeds"],
        "cfg_note": ["4th row and column missing", "first three columns",
                     "first three rows", "all sixteen"],
        "got": "measured", "miss": "out of reach",
        "fig2": "One sample seen by RCE and by DRC-MME: waveform and harmonics",
        "cycle": "base rotation angle $C$ (deg)", "norm_int": "normalized intensity",
        "order": "harmonic order $n$ (frequency $2nC$)", "mag": "normalized coefficient",
        "rce_t": "RCE from post 1 — one compensator",
        "mme_t": "DRC-MME — two compensators (5:3)",
        "rce_note": "only $2\\omega$ and $4\\omega$",
        "mme_note": "24 of them, $2C$ through $32C$",
        "dead": "$n$ = 9, 12, 14, 15 vanish",
        "fig3": "Which Mueller element rides on which harmonic",
        "elem": "elements carried", "factor": "strength factor",
        "fig3note": "the 4th row and column ride on $\\sin\\delta$ alone — at 0° or 180° they disappear",
        "fig4": "The ratio sets how often you must read, not how precise you are",
        "ratio_ax": "rotation ratio $m_1 : m_2$",
        "high_ax": "highest harmonic", "nint_ax": "minimum integrations per cycle",
        "fig4a": "(a) a larger ratio pushes the harmonics up",
        "fig4b": "(b) yet the condition numbers are identical",
        "cond_ax": "condition number of the reduction matrix",
        "ret_ax": "compensator retardance (deg)",
        "fig5": "It is the retardance that sets the precision",
        "opt": "minimum at {:.0f}°  (condition number {:.2f})",
        "at90": "{:.2f} at 90°",
        "gain": "choosing {:.0f}° over 90° cuts the\nnoise amplification by {:.1f}x",
        "overlap": "5:3 and 5:1 coincide exactly",
    },
}


# ---------------------------------------------------------------------------
# 그림 1 — 배치별 도달 범위
# ---------------------------------------------------------------------------

def figure1(L):
    reach = [
        np.array([[1, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0], [0, 0, 0, 0]]),
        np.array([[1, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0], [1, 1, 1, 0]]),
        np.array([[1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1], [0, 0, 0, 0]]),
        np.ones((4, 4), dtype=int),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(10.6, 3.5))
    for ax, R, title, note in zip(axes, reach, L["cfg"], L["cfg_note"]):
        for j in range(4):
            for k in range(4):
                on = R[j, k]
                ax.plot(k, 3 - j, "o", ms=17 if on else 8,
                        color=C_MME if on else "#ccc",
                        mec="#444" if on else "#aaa", mew=0.8)
        ax.set_xlim(-0.7, 3.7)
        ax.set_ylim(-0.7, 3.7)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
        for sp in ax.spines.values():
            sp.set_color("#bbb")
        ax.set_title(title, fontsize=9.6, pad=8, **L["font"])
        ax.text(1.5, -1.15, note, ha="center", fontsize=9.6,
                color=C_MME if R.sum() == 16 else "#555",
                weight="bold" if R.sum() == 16 else "normal", **L["font"])
    fig.suptitle(L["fig1"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fig.savefig(os.path.join(L["dir"], "fig1-mme-configurations.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 파형과 고조파
# ---------------------------------------------------------------------------

def figure2(L):
    C = optical_cycle(N_PTS)
    d = np.pi / 2
    I_mme = intensity(SAMPLE, C, d, d, RATIO)
    I_rce = intensity(SAMPLE, C, 0.0, d, RATIO, pol=np.pi / 4)   # 보상자 하나
    n_max = 18
    m_mme = harmonic_magnitudes(I_mme, n_max)
    m_rce = harmonic_magnitudes(I_rce, n_max)

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 6.6),
                             gridspec_kw={"width_ratios": [1.5, 1.0]})
    deg = np.rad2deg(C)
    for row, (I, mag, col, title, note) in enumerate(
            [(I_rce, m_rce, C_RCE, L["rce_t"], L["rce_note"]),
             (I_mme, m_mme, C_MME, L["mme_t"], L["mme_note"])]):
        aw, ah = axes[row]
        aw.plot(deg, I / I.mean(), color=col, lw=1.5)
        aw.set_xlim(0, 180)
        aw.set_xticks(np.arange(0, 181, 45))
        aw.set_ylabel(L["norm_int"], fontsize=10, **L["font"])
        aw.set_title(title, fontsize=11, loc="left", **L["font"])
        aw.grid(alpha=0.25)

        idx = np.arange(1, n_max + 1)
        ah.bar(idx, mag, color=col, alpha=0.9, width=0.7)
        if row == 1:
            for n in (9, 12, 14, 15):
                ah.plot(n, 0, "v", color=C_DEAD, ms=7, clip_on=False)
            ah.text(0.97, 0.93, L["dead"], transform=ah.transAxes, ha="right",
                    va="top", fontsize=9, color=C_DEAD, **L["font"])
        ah.set_xlim(0.3, n_max + 0.7)
        ah.set_xticks([1, 4, 8, 12, 16])
        ah.set_ylabel(L["mag"], fontsize=9.5, **L["font"])
        ah.set_title(note, fontsize=10.5, loc="left", color=col, **L["font"])
        ah.grid(alpha=0.25, axis="y")

    axes[1][0].set_xlabel(L["cycle"], fontsize=10.5, **L["font"])
    axes[1][1].set_xlabel(L["order"], fontsize=10, **L["font"])
    fig.suptitle(L["fig2"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(L["dir"], "fig2-waveform-spectrum.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — 계수 지도
# ---------------------------------------------------------------------------

def pretty(text, core):
    """표에 쓸 표기로 다듬는다. Mjk -> $M_{jk}$, s1/d1 -> 첨자, 2x2 -> 중심 2x2."""
    out = []
    for tok in text.split():
        if tok == "2x2":
            out.append(core)
        elif tok.startswith("M") and tok[1:].isdigit():
            out.append(f"$M_{{{tok[1:]}}}$")
        elif tok == "sin":
            out.append("sin")
        elif tok in ("d1", "d2"):
            out.append(f"$\\delta_{tok[1]}$")
        elif tok in ("s1", "s2"):
            out.append(f"$s_{tok[1]}$")
        else:
            out.append(tok)
    return " ".join(out)


def figure3(L):
    owner = harmonic_owner()
    ns = sorted(owner)
    fig, ax = plt.subplots(figsize=(9.6, 5.0))
    for i, n in enumerate(ns):
        elem, fac = owner[n]
        y = len(ns) - 1 - i
        sin_driven = fac.startswith("sin") or "sin" in fac
        col = C_MME if "sin d" in fac else "#5b8ff9"
        ax.barh(y, 1.0, left=0, height=0.62, color=col, alpha=0.20, lw=0)
        ax.text(0.02, y, f"{2*n}C", va="center", fontsize=10, weight="bold",
                color=col, **L["font"])
        core = "중심 $2\\times2$" if L is LABELS["ko"] else "central $2\\times2$"
        ax.text(0.17, y, pretty(elem, core), va="center", fontsize=10, **L["font"])
        ax.text(0.62, y, pretty(fac, core), va="center", fontsize=10,
                color="#444", **L["font"])
    ax.text(0.17, len(ns) - 0.25, L["elem"], fontsize=9.5, color="#888", **L["font"])
    ax.text(0.62, len(ns) - 0.25, L["factor"], fontsize=9.5, color="#888", **L["font"])
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-1.1, len(ns) + 0.1)
    ax.axis("off")
    ax.text(0.5, -0.85, L["fig3note"], ha="center", fontsize=10, color=C_MME,
            **L["font"])
    ax.set_title(L["fig3"], fontsize=12.5, pad=12, **L["font"])
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig3-coefficient-map.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — 회전비가 정하는 것
# ---------------------------------------------------------------------------

def figure4(L):
    highs = [2 * highest_harmonic(r) for r in RATIOS]
    names = [f"{a}:{b}" for a, b in RATIOS]
    x = np.arange(len(RATIOS))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.3))
    bars = ax1.bar(x, highs, color=["#5b8ff9", "#61ddaa", C_MME, "#f6bd16"],
                   width=0.6, ec="#444", lw=0.8)
    for xi, h in zip(x, highs):
        ax1.text(xi, h + 0.8, f"{h}C\n({h+1}회)" if L is LABELS["ko"] else f"{h}C\n({h+1})",
                 ha="center", fontsize=9.5, **L["font"])
    ax1.set_xticks(x)
    ax1.set_xticklabels(names)
    ax1.set_ylim(0, max(highs) * 1.28)
    ax1.set_xlabel(L["ratio_ax"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["high_ax"], fontsize=10.5, **L["font"])
    ax1.set_title(L["fig4a"], fontsize=11.5, loc="left", **L["font"])
    ax1.grid(alpha=0.25, axis="y")

    for r, col, ls in zip([(5, 3), (5, 1)], [C_MME, "#5b8ff9"], ["-", "--"]):
        cond = cond_sweep(r)
        ax2.plot(DEG, cond, ls, color=col, lw=2.6 if ls == "-" else 1.6,
                 label=f"{r[0]}:{r[1]}")
    ax2.set_yscale("log")
    ax2.set_xlabel(L["ret_ax"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["cond_ax"], fontsize=10.5, **L["font"])
    ax2.set_title(L["fig4b"], fontsize=11.5, loc="left", **L["font"])
    ax2.legend(prop=L["legend"], fontsize=10)
    ax2.text(0.5, 0.55, L["overlap"], transform=ax2.transAxes, ha="center",
             fontsize=10, color=C_MME, **L["font"])
    ax2.grid(alpha=0.25, which="both")

    fig.suptitle(L["fig4"], fontsize=12.5, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(L["dir"], "fig4-ratio-tradeoff.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 5 — 조건수와 최적 지연량
# ---------------------------------------------------------------------------

def figure5(L):
    cond = cond_sweep(RATIO)
    i = int(np.argmin(cond))
    d_opt, c_opt = DEG[i], cond[i]
    c90 = float(np.interp(90.0, DEG, cond))

    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    ax.plot(DEG, cond, color=C_MME, lw=2.6)
    ax.plot([d_opt], [c_opt], "o", color=C_MME, ms=9, zorder=5)
    ax.plot([90.0], [c90], "o", color="#5b8ff9", ms=9, zorder=5)
    ax.axvline(90.0, color="#5b8ff9", lw=1.0, ls=":")
    ax.axvline(d_opt, color=C_MME, lw=1.0, ls=":")
    ax.annotate(L["opt"].format(d_opt, c_opt), xy=(d_opt, c_opt),
                xytext=(d_opt + 4, c_opt + 13.0), fontsize=10, color=C_MME,
                arrowprops=dict(arrowstyle="->", color=C_MME, lw=1.2), **L["font"])
    ax.annotate(L["at90"].format(c90), xy=(90.0, c90), xytext=(97, c90 + 9.0),
                fontsize=10, color="#3a6ea5",
                arrowprops=dict(arrowstyle="->", color="#3a6ea5", lw=1.2), **L["font"])
    ax.text(0.62, 0.93, L["gain"].format(d_opt, c90 / c_opt),
            transform=ax.transAxes, ha="center", va="top", fontsize=10.5,
            color="#333", linespacing=1.45, **L["font"])
    ax.set_xlim(DEG[0], DEG[-1])
    ax.set_ylim(0, min(cond.max(), 34))
    ax.set_xlabel(L["ret_ax"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["cond_ax"], fontsize=10.5, **L["font"])
    ax.set_title(L["fig5"], fontsize=12.5, pad=10, **L["font"])
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig5-condition-number.png"), dpi=150)
    plt.close(fig)


def main():
    cond = cond_sweep(RATIO)
    i = int(np.argmin(cond))
    print(f"회전비 {RATIO[0]}:{RATIO[1]}, 최고차 고조파 {2*highest_harmonic(RATIO)}C, "
          f"최소 적분 {2*highest_harmonic(RATIO)+1}회")
    print(f"조건수 최소 {DEG[i]:.0f}° ({cond[i]:.2f}), 90° 에서 "
          f"{np.interp(90.0,DEG,cond):.2f} -> {np.interp(90.0,DEG,cond)/cond[i]:.2f}배 개선")
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        figure1(L); figure2(L); figure3(L); figure4(L); figure5(L)
        print(f"  [{lang}] 그림 5개 저장 -> {os.path.normpath(L['dir'])}")


if __name__ == "__main__":
    main()

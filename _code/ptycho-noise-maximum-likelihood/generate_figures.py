"""계산 이미징과 타이코그래피 6편 그림 6개 생성 (한국어판·영문판). 수치는 studies6.py 의 캐시를 읽는다."""
import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import studies6 as S

SLUG = "ptycho-noise-maximum-likelihood"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
COL = {"rPIE": "#c8553d", "rPIE_decay": "#e09f3e", "DM": "#555555", "amp": "#2f6690", "poisson": "#6a994e",
       "int": "#9a6fb0", "rPIE_nopc": "#8c2f39"}

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "ph": "패턴당 광자", "le1": "0~1개인 화소", "count": "광자 수",
        "lam": "평균 광자 수", "var": "분산", "v_n": "n 의 분산 (= 평균)", "v_s": r"$\sqrt{n}$ 의 분산", "quarter": "1/4",
        "photons": "가장 밝은 패턴의 광자 수", "err": "시편 오차", "slope": r"$1/\sqrt{N}$ 기울기",
        "rPIE": "rPIE", "rPIE_decay": "rPIE, 걸음을 줄임", "DM": "DM", "amp": "ML (진폭)", "poisson": "ML (푸아송)",
        "int": "세기 최소자승", "rPIE_nopc": "rPIE, 세기 제어 없음",
        "truth": "참 시편 (위상)", "iters": "반복 횟수", "unlit": "조명이 닿지 않는 화소의 진폭",
        "npos": "위치 수 (총광자 수 고정)",
        "fig1": "그림1. 같은 회절 패턴을 세 선량으로 뽑았다(로그 색). 광자가 적으면 바깥쪽 화소는 대부분 0이나 1이다",
        "fig2": "그림2. 광자 수 n의 분산은 평균과 같지만, √n의 분산은 평균이 몇 개만 넘으면 1/4 근처로 고정된다",
        "fig3": "그림3. 프로브를 알 때 선량별 최종 오차(300회, 노이즈 표본 3개의 중앙값과 범위). 진폭과 푸아송 ML이 가장 낮다",
        "fig4": "그림4. 패턴당 $10^5$ 광자에서 프로브를 알 때의 복원 위상. 같은 데이터에서 비용 함수와 갱신 방식만 다르다",
        "fig5": "그림5. 패턴당 $10^5$ 광자, 프로브를 모를 때. 계속 돌린 rPIE는 척도 방향으로 표류해 무너지고, 50회부터 ML로 이어 정련하면 곧바로 내려가 그대로 머문다",
        "fig6": "그림6. 총광자 수를 $1.44\\times10^8$으로 고정하고 위치 수를 바꿨다. 오차는 위치 수보다 총광자 수가 정한다",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "ph": "photons per pattern", "le1": "pixels with 0–1", "count": "Photon count",
        "lam": "Mean photon count", "var": "Variance", "v_n": "Variance of n (= mean)", "v_s": r"Variance of $\sqrt{n}$",
        "quarter": "1/4",
        "photons": "Photons in the brightest pattern", "err": "Specimen error", "slope": r"$1/\sqrt{N}$ slope",
        "rPIE": "rPIE", "rPIE_decay": "rPIE, decaying step", "DM": "DM", "amp": "ML (amplitude)", "poisson": "ML (Poisson)",
        "int": "Intensity least squares", "rPIE_nopc": "rPIE, no power control",
        "truth": "True specimen (phase)", "iters": "Iterations", "unlit": "Amplitude of an unilluminated pixel",
        "npos": "Number of positions (fixed total photons)",
        "fig1": "Fig 1. The same diffraction pattern drawn at three doses (log colour). At low dose most outer pixels hold 0 or 1 photon",
        "fig2": "Fig 2. The variance of the photon count n equals its mean, but the variance of √n settles near 1/4 once the mean exceeds a few",
        "fig3": "Fig 3. Final error versus dose with the probe known (300 iterations, median and range of 3 noise draws). Amplitude and Poisson ML are lowest",
        "fig4": "Fig 4. Reconstructed phase at $10^5$ photons per pattern with the probe known. Same data; only the cost and update differ",
        "fig5": "Fig 5. $10^5$ photons per pattern, probe unknown. Continued rPIE drifts along the scale direction and collapses; ML refinement from iteration 50 drops at once and stays down",
        "fig6": "Fig 6. Total photons fixed at $1.44\\times10^8$ while the number of positions changes. The total, not the number of positions, sets the error",
    },
}

R = S.load()


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.07, 1, 1)):
        t = fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **Fk)
        fig.canvas.draw()
        if t.get_window_extent().width > 0.98 * fig.bbox.width:  # 영문 캡션이 폭을 넘으면 줄바꿈
            text = textwrap.fill(L[key], width=int(len(L[key]) * 0.95 * fig.bbox.width / t.get_window_extent().width))
            t.set_text(text); t.set_va("bottom"); t.set_y(0.012)
            rect = (rect[0], rect[1] + 0.03 * text.count("\n"), rect[2], rect[3])
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림1: 세 선량의 회절 패턴
    Pt = R["patterns"]
    fig, ax = plt.subplots(1, 3, figsize=(12, 4.2))
    for a, ph, lab in zip(ax, (1e4, 1e6, 1e8), (r"$10^4$", r"$10^6$", r"$10^8$")):
        v = Pt[ph]
        im = a.imshow(v["pattern"] + 0.5, norm=LogNorm(vmin=0.5, vmax=max(v["max"], 1)), cmap="magma")
        a.set_title(f"{lab} {L['ph']}  ({L['le1']} {v['frac01']:.1%})", fontsize=9.5, **Fk); a.set_axis_off()
        cb = fig.colorbar(im, ax=a, fraction=0.046, pad=0.02); cb.set_label(L["count"], **Fk)
    cap(fig, "fig1"); save(fig, "fig1-photon-patterns.png")

    # 그림2: 분산 안정화
    V = R["varstab"]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.loglog(V["lam"], V["var_n"], color=COL["int"], lw=1.8, label=L["v_n"])
    ax.loglog(V["lam"], V["var_sqrt"], color=COL["amp"], lw=1.8, label=L["v_s"])
    ax.axhline(0.25, color="k", ls="--", lw=0.9, label=L["quarter"])
    ax.set_xlabel(L["lam"], **Fk); ax.set_ylabel(L["var"], **Fk); ax.grid(alpha=0.3, which="both")
    ax.legend(prop=LG, fontsize=9)
    cap(fig, "fig2", rect=(0, 0.08, 1, 1)); save(fig, "fig2-variance-stabilization.png")

    # 그림3: 선량별 오차
    D = R["dose"]
    doses = np.array(S.DOSES)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for me in S.METHODS:
        e = np.array([[D[(ph, me, sd)][300] for sd in S.SEEDS] for ph in S.DOSES])
        ax.loglog(doses, np.median(e, 1), "o-", color=COL[me], ms=4.5, label=L[me])
        ax.fill_between(doses, e.min(1), e.max(1), color=COL[me], alpha=0.2)
    ref = np.median([D[(1e4, "poisson", sd)][300] for sd in S.SEEDS]) * np.sqrt(1e4 / doses) * 0.7
    ax.loglog(doses, ref, "k--", lw=0.9, label=L["slope"])
    ax.set_xlabel(L["photons"], **Fk); ax.set_ylabel(L["err"], **Fk); ax.grid(alpha=0.3, which="both")
    ax.legend(prop=LG, fontsize=8.5, ncol=2, loc="lower left")
    cap(fig, "fig3", rect=(0, 0.06, 1, 1)); save(fig, "fig3-error-vs-dose.png")

    # 그림4: 10^5 광자의 복원 위상
    Im = R["images"]
    obj, m = Im["obj"], Im["mask"]
    ys, xs = np.where(m); bb = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.9))
    ph0 = np.angle(obj)
    vmin, vmax = ph0[m].min(), ph0[m].max()
    panels = [(obj, L["truth"])] + [(Im["rec"][k][0], f"{L[k]}: {Im['rec'][k][1]:.2g}") for k in ("rPIE", "poisson", "int")]
    for a, (img, t) in zip(ax, panels):
        a.imshow(np.where(m, np.angle(img), np.nan)[bb], cmap="gray", vmin=vmin, vmax=vmax)
        a.set_title(t, fontsize=9.5, **Fk); a.set_axis_off()
    cap(fig, "fig4", rect=(0, 0.06, 1, 1)); save(fig, "fig4-reconstructions.png")

    # 그림5: 낮은 선량, 프로브를 모를 때
    Lw = R["lowdose"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for me in ("rPIE", "rPIE_nopc", "amp", "poisson"):
        for sd in S.SEEDS:
            r = Lw[(1e5, me, sd)]["rows"]
            ax[0].semilogy(r[:, 0], r[:, 1], color=COL[me], lw=1.2, alpha=0.85, label=L[me] if sd == 0 else None)
            if me != "rPIE_nopc":
                ax[1].semilogy(r[:, 0], r[:, 3], color=COL[me], lw=1.2, alpha=0.85, label=L[me] if sd == 0 else None)
    for a in ax:
        a.axvline(50, color="#999999", ls=":", lw=0.9)
        a.set_xlabel(L["iters"], **Fk); a.grid(alpha=0.3, which="both")
    ax[0].set_ylabel(L["err"], **Fk); ax[1].set_ylabel(L["unlit"], **Fk)
    ax[0].legend(prop=LG, fontsize=8.5); ax[1].legend(prop=LG, fontsize=8.5)
    cap(fig, "fig5"); save(fig, "fig5-low-dose-drift.png")

    # 그림6: 광자 나누기
    F = R["frac"]
    st = sorted(F)
    npos = np.array([F[s]["npos"] for s in st])
    e = np.array([F[s]["err"] for s in st])
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    ax.plot(npos, np.median(e, 1), "o-", color=COL["amp"])
    ax.fill_between(npos, e.min(1), e.max(1), color=COL["amp"], alpha=0.2)
    ax.set_xscale("log"); ax.set_ylim(0, 0.03)
    ax.set_xticks([50, 100, 200, 400], ["50", "100", "200", "400"]); ax.minorticks_off()
    ax.set_xlabel(L["npos"], **Fk); ax.set_ylabel(L["err"], **Fk); ax.grid(alpha=0.3, which="both")
    cap(fig, "fig6", rect=(0, 0.08, 1, 1)); save(fig, "fig6-dose-fractionation.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")

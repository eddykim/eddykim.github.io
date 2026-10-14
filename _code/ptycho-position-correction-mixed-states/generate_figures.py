"""계산 이미징과 타이코그래피 7편 그림 6개 생성 (한국어판·영문판). 수치는 studies7.py 의 캐시를 읽는다."""
import os

import matplotlib.pyplot as plt
import numpy as np

import studies7 as S

SLUG = "ptycho-position-correction-mixed-states"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
COL = {None: "#c8553d", "anneal": "#2f6690", "grad": "#6a994e", 0.012: "#2f6690", 0.02: "#c8553d"}

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "truth": "참 시편 (위상)", None: "보정 없음", "anneal": "어닐링", "grad": "기울기",
        "sig": "위치 오차 (rms, 화소)", "err": "시편 오차", "iters": "반복 횟수", "prms": "남은 위치 오차 (rms, 화소)",
        "before": "명목 위치의 오차", "after": "어닐링 뒤 남은 오차 (같은 축척)",
        "coh": "결맞은 광원", "src": "광원 폭", "mode": "모드", "K": "모드 수 K", "resid": "데이터 잔차",
        "rpie": "혼합 상태 rPIE", "ml": "ML 정련", "trunc": "참 상태를 K개로 자름",
        "k1": "모드 하나", "k3": "모드 셋", "true_mode": "참 모드", "rec_mode": "복원한 모드",
        "fig1": "위치를 rms 1.4화소 틀린 채 복원한 위상(같은 데이터). 보정 없이는 무너지고, 위치를 함께 고치면 되살아난다",
        "fig2": "왼쪽: 위치 오차에 따른 최종 시편 오차(오차 표본 3개). 오른쪽: rms 1.4화소에서 남은 위치 오차의 수렴",
        "fig3": "위치 144개의 오차를 화살표로 그렸다(5배 확대). 어닐링 뒤에는 거의 점으로 줄어든다",
        "fig4": "넓은 광원이 만든 프로브의 결맞음 모드(진폭)와 점유율. 광원이 넓을수록 첫 모드의 몫이 줄어든다",
        "fig5": "모드 수에 따른 시편 오차와 데이터 잔차(순서 3개). 점선은 참 상태를 K개 모드로 잘랐을 때의 잔차다",
        "fig6": "지배 모드 72 %인 데이터에서 모드 하나와 셋으로 복원한 위상, 그리고 참 모드와 복원한 모드",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "truth": "True specimen (phase)", None: "No correction", "anneal": "Annealing", "grad": "Gradient",
        "sig": "Position error (rms, px)", "err": "Specimen error", "iters": "Iterations", "prms": "Remaining position error (rms, px)",
        "before": "Error of the nominal positions", "after": "Error left after annealing (same scale)",
        "coh": "Coherent source", "src": "source width", "mode": "mode", "K": "Number of modes K", "resid": "Data residual",
        "rpie": "Mixed-state rPIE", "ml": "ML refinement", "trunc": "True state truncated to K modes",
        "k1": "One mode", "k3": "Three modes", "true_mode": "True mode", "rec_mode": "Recovered mode",
        "fig1": "Phase reconstructed from the same data with positions wrong by 1.4 px rms. Uncorrected it collapses; correcting the positions recovers it",
        "fig2": "Left: final specimen error versus position error (3 error draws). Right: convergence of the remaining position error at 1.4 px rms",
        "fig3": "Errors of the 144 positions drawn as arrows (magnified 5×). After annealing they shrink almost to dots",
        "fig4": "Coherent modes (amplitude) and occupancies of the probe made by an extended source. The wider the source, the smaller the first mode's share",
        "fig5": "Specimen error and data residual versus number of modes (3 orders). The dashed line is the residual of the true state truncated to K modes",
        "fig6": "Phase reconstructed with one and three modes from data with a 72 % dominant mode, and the true and recovered modes",
    },
}

R = S.load()


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.07, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **Fk)
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림1: 위치 오차 1 화소의 복원
    Pi = R["posimg"]
    obj, m = Pi["obj"], Pi["mask"]
    ys, xs = np.where(m); bb = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    ph = np.angle(obj); vmin, vmax = ph[m].min(), ph[m].max()
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.9))
    panels = [(obj, L["truth"])] + [(Pi["runs"][c]["rec"], f"{L[c]}: {Pi['runs'][c]['err']:.2g}") for c in (None, "anneal", "grad")]
    for a, (img, t) in zip(ax, panels):
        a.imshow(np.where(m, np.angle(img), np.nan)[bb], cmap="gray", vmin=vmin, vmax=vmax)
        a.set_title(t, fontsize=9.5, **Fk); a.set_axis_off()
    cap(fig, "fig1", rect=(0, 0.06, 1, 1)); save(fig, "fig1-position-error-recon.png")

    # 그림2: 위치 오차에 따른 오차, 위치 수렴
    P = R["poserr"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    for c in (None, "anneal", "grad"):
        sig = [s for s in S.SIGMAS if (s, c, 0) in P]
        x = np.array([np.median([P[(s, c, sd)]["rms0"] for sd in S.SEEDS]) for s in sig])
        e = np.array([[P[(s, c, sd)]["trace"][300][0] for sd in S.SEEDS] for s in sig])
        ax[0].semilogy(x, np.median(e, 1), "o-", color=COL[c], label=L[c])
        ax[0].fill_between(x, e.min(1), e.max(1), color=COL[c], alpha=0.2)
    base = np.median([P[(0.0, None, sd)]["trace"][300][0] for sd in S.SEEDS])
    ax[0].axhline(base, color="k", ls=":", lw=0.9)
    ax[0].set_xlabel(L["sig"], **Fk); ax[0].set_ylabel(L["err"], **Fk); ax[0].legend(prop=LG, fontsize=9); ax[0].grid(alpha=0.3, which="both")
    for c in ("anneal", "grad"):
        for sd in S.SEEDS:
            tr = P[(1.0, c, sd)]["trace"]
            ax[1].semilogy([0] + list(S.CK), [P[(1.0, c, sd)]["rms0"]] + [tr[k][1] for k in S.CK], "-", color=COL[c], alpha=0.8,
                           label=L[c] if sd == 0 else None)
    ax[1].set_xlabel(L["iters"], **Fk); ax[1].set_ylabel(L["prms"], **Fk); ax[1].legend(prop=LG, fontsize=9); ax[1].grid(alpha=0.3, which="both")
    cap(fig, "fig2"); save(fig, "fig2-position-correction.png")

    # 그림3: 위치 오차 화살표
    pt, pn, pa = Pi["pos_t"], Pi["pos_n"], Pi["runs"]["anneal"]["pos"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 5.4))
    for a, d, key in ((ax[0], pn - pt, "before"), (ax[1], pa - pt, "after")):
        d = d - d.mean(0)
        a.quiver(pt[:, 1], -pt[:, 0], 5 * d[:, 1], -5 * d[:, 0], angles="xy", scale_units="xy", scale=1, width=0.004, color="#2f6690")
        a.plot(pt[:, 1], -pt[:, 0], "k.", ms=2)
        a.set_aspect("equal"); a.set_axis_off(); a.set_title(L[key], fontsize=10, **Fk)
    cap(fig, "fig3", rect=(0, 0.05, 1, 1)); save(fig, "fig3-position-arrows.png")

    # 그림4: 결맞음 모드
    MI = R["modes_info"]
    fig, ax = plt.subplots(3, 4, figsize=(11, 8.6))
    for i, st in enumerate((0.0,) + S.SOURCES):
        info = MI[st]
        for k in range(4):
            a = ax[i, k]
            if k < len(info["modes"]):
                a.imshow(np.abs(info["modes"][k]), cmap="gray")
                a.set_title(f"{L['mode']} {k + 1}: {info['occ'][k]:.0%}", fontsize=9, **Fk)
            a.set_axis_off()
        ax[i, 0].text(-0.15, 0.5, L["coh"] if st == 0 else f"{L['src']} {st}", transform=ax[i, 0].transAxes,
                      rotation=90, va="center", ha="center", fontsize=9.5, **Fk)
    cap(fig, "fig4", rect=(0.02, 0.04, 1, 1)); save(fig, "fig4-coherent-modes.png", dpi=120)

    # 그림5: 모드 수
    Mo, Tr = R["modes"], R["trunc"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    for st in S.SOURCES:
        Ks = np.array(S.KS)
        for i, key in enumerate(("err", "resid")):
            v = np.array([[Mo[(st, K, od)]["rpie"][key] for od in S.SEEDS] for K in S.KS])
            ax[i].semilogy(Ks, np.median(v, 1), "o-", color=COL[st], label=f"{L['rpie']}, {L['src']} {st}")
            ax[i].fill_between(Ks, v.min(1), v.max(1), color=COL[st], alpha=0.2)
            refined = [K for K in S.KS if "ml" in Mo[(st, K, 0)]]
            for K in refined:
                vv = [Mo[(st, K, od)]["ml"][key] for od in S.SEEDS]
                ax[i].semilogy([K] * len(vv), vv, "s", mfc="none", color=COL[st], label=L["ml"] if (K == refined[0] and st == 0.02) else None)
        tk = sorted(Tr[st])
        ax[1].semilogy(tk, [Tr[st][k] for k in tk], "--", color=COL[st], lw=1, label=L["trunc"] if st == 0.012 else None)
    for a, key in zip(ax, ("err", "resid")):
        a.set_xlabel(L["K"], **Fk); a.set_ylabel(L[key], **Fk); a.grid(alpha=0.3, which="both")
    ax[0].legend(prop=LG, fontsize=8); ax[1].legend(prop=LG, fontsize=8)
    cap(fig, "fig5"); save(fig, "fig5-modes-vs-k.png")

    # 그림6: 모드 하나와 셋
    Mi = R["modeimg"]
    obj, m = Mi["obj"], Mi["mask"]
    ys, xs = np.where(m); bb = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    ph = np.angle(obj); vmin, vmax = ph[m].min(), ph[m].max()
    fig = plt.figure(figsize=(12, 6.3))
    gs = fig.add_gridspec(2, 6)
    for j, (img, t) in enumerate([(obj, L["truth"]), (Mi["runs"][1]["rec"], f"{L['k1']}: {Mi['runs'][1]['err']:.2g}"),
                                  (Mi["runs"][3]["rec"], f"{L['k3']}: {Mi['runs'][3]['err']:.2g}")]):
        a = fig.add_subplot(gs[0, 2 * j:2 * j + 2])
        a.imshow(np.where(m, np.angle(img), np.nan)[bb], cmap="gray", vmin=vmin, vmax=vmax)
        a.set_title(t, fontsize=9.5, **Fk); a.set_axis_off()
    rec = Mi["runs"][3]["modes"]
    pw = np.sum(np.abs(rec) ** 2, (1, 2)); pw = pw / pw.sum()
    for k in range(3):
        a = fig.add_subplot(gs[1, k]); a.imshow(np.abs(Mi["modes_t"][k]), cmap="gray"); a.set_axis_off()
        a.set_title(f"{L['true_mode']} {k + 1}: {Mi['occ'][k]:.0%}", fontsize=8.5, **Fk)
        a = fig.add_subplot(gs[1, 3 + k]); a.imshow(np.abs(rec[k]), cmap="gray"); a.set_axis_off()
        a.set_title(f"{L['rec_mode']} {k + 1}: {pw[k]:.0%}", fontsize=8.5, **Fk)
    cap(fig, "fig6", rect=(0, 0.04, 1, 1)); save(fig, "fig6-mode-reconstruction.png", dpi=120)


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")

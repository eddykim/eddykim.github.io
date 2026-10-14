"""계산 이미징과 타이코그래피 5편 그림 6개 생성 (한국어판·영문판). 수치는 studies5.py 의 캐시를 읽는다."""
import os

import matplotlib.pyplot as plt
import numpy as np

import studies5 as S

SLUG = "ptycho-momentum-difference-map"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
COL = {"rPIE": "#c8553d", "rPIE_warm": "#2f6690", "mPIE": "#6a994e", "mPIE_ptylab": "#bc8a2c", "DM": "#555555"}

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "gd": "경사하강", "mom6": r"모멘텀 $\eta$ = 0.6", "mom9": r"모멘텀 $\eta$ = 0.9", "step": "걸음 수", "dist": "최솟값까지 거리", "x": "x (완만한 방향)", "y": "y (가파른 방향)",
        "rPIE": r"rPIE ($\gamma_{obj}$ = 1)", "rPIE_warm": r"rPIE, 처음 10회만 $\gamma_{obj}$ = 0.25", "mPIE": "mPIE (Maiden 값)",
        "mPIE_ptylab": "mPIE (PtyLab 기본값)", "DM": "DM",
        "iters": "반복 횟수", "err": "시편 오차", "resid": "데이터 잔차",
        "pinhole": "핀홀 프로브", "defocus": "초점이 어긋난 집속 프로브", "diffuser": "확산판 프로브",
        "power": "프로브 세기 / 목표 세기", "nopc": "mPIE, 세기 제어 없음", "pc": "세기 제어 있음",
        "unlit": "조명이 닿지 않는 화소의 진폭", "d_mPIE": "mPIE (모멘텀 둘 다)", "d_mom_obj": "모멘텀을 시편에만",
        "d_mom_prb": "모멘텀을 프로브에만", "d_rPIE": "rPIE (모멘텀 없음)",
        "eta": r"모멘텀 $\eta$ (마찰 = 피드백)", "gamma": r"시편 걸음 $\gamma_{obj}$",
        "t_p": "참 프로브", "found": "복원한 프로브", "err_plain": "시편 오차 (프로브 복원)", "err_known": "시편 오차 (프로브를 앎)",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "gd": "Gradient descent", "mom6": r"Momentum $\eta$ = 0.6", "mom9": r"Momentum $\eta$ = 0.9", "step": "Step", "dist": "Distance to minimum",
        "x": "x (shallow direction)", "y": "y (steep direction)",
        "rPIE": r"rPIE ($\gamma_{obj}$ = 1)", "rPIE_warm": r"rPIE, $\gamma_{obj}$ = 0.25 for the first 10", "mPIE": "mPIE (Maiden values)",
        "mPIE_ptylab": "mPIE (PtyLab defaults)", "DM": "DM",
        "iters": "Iterations", "err": "Specimen error", "resid": "Data residual",
        "pinhole": "Pinhole probe", "defocus": "Defocused focused probe", "diffuser": "Diffuser probe",
        "power": "Probe power / target", "nopc": "mPIE, no power control", "pc": "With power control",
        "unlit": "Amplitude of an unilluminated pixel", "d_mPIE": "mPIE (momentum on both)", "d_mom_obj": "Momentum on specimen only",
        "d_mom_prb": "Momentum on probe only", "d_rPIE": "rPIE (no momentum)",
        "eta": r"Momentum $\eta$ (friction = feedback)", "gamma": r"Specimen step $\gamma_{obj}$",
        "t_p": "True probe", "found": "Recovered probe", "err_plain": "Specimen error (probe recovered)",
        "err_known": "Specimen error (probe known)",
    },
}

R = S.load()
CK = np.array(S.CK)


def traces(case, alg, i):
    rs = [v for (c, a, s), v in R["compare"].items() if c == case and a == alg]
    return np.array([[r["trace"][k][i] for k in S.CK] for r in rs])


def band(ax, case, alg, i, label):
    v = traces(case, alg, i)
    with np.errstate(all="ignore"):
        ax.loglog(CK, np.median(v, 0), "o-", ms=3.5, color=COL[alg], label=label)
        if len(v) > 1:
            ax.fill_between(CK, v.min(0), v.max(0), color=COL[alg], alpha=0.18)


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0, 1, 1)):
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림1: 좁은 골짜기
    T = R["toy"]
    TC = ((0.0, COL["rPIE"], "gd"), (0.6, COL["rPIE_warm"], "mom6"), (0.9, COL["mPIE"], "mom9"))
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4), gridspec_kw={"width_ratios": [1.6, 1]})
    xx, yy = np.meshgrid(np.linspace(-10, 6, 320), np.linspace(-1.6, 1.6, 200))
    ax[0].contour(xx, yy, 0.5 * (xx**2 + 25 * yy**2), levels=np.geomspace(0.05, 80, 12), colors="#bbbbbb", linewidths=0.7)
    for e, cc, k in TC:
        p = T[e][:45]
        ax[0].plot(p[:, 0], p[:, 1], "o-", ms=2.2, lw=0.9, color=cc, label=L[k])
    ax[0].plot(0, 0, "k*", ms=9)
    ax[0].set_xlabel(L["x"], **Fk); ax[0].set_ylabel(L["y"], **Fk); ax[0].legend(prop=LG, fontsize=8.5, loc="lower right")
    for e, cc, k in TC:
        ax[1].semilogy(np.linalg.norm(T[e], axis=1), color=cc, label=L[k])
    ax[1].set_ylim(1e-3, 30)
    ax[1].set_xlabel(L["step"], **Fk); ax[1].set_ylabel(L["dist"], **Fk); ax[1].grid(alpha=0.3, which="both")
    cap(fig, "fig1"); save(fig, "fig1-momentum-valley.png")

    # 그림2: 확산판
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4))
    for i, key in enumerate(("err", "resid")):
        for alg in ("rPIE", "rPIE_warm", "mPIE"):
            band(ax[i], "diffuser", alg, i, L[alg])
        ax[i].set_xlabel(L["iters"], **Fk); ax[i].set_ylabel(L[key], **Fk); ax[i].grid(alpha=0.3, which="both")
    ax[0].legend(prop=LG, fontsize=8.5)
    cap(fig, "fig2"); save(fig, "fig2-diffuser-step.png")

    # 그림3: 표류와 그 누적
    D = R["drift"]
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4))
    for s_ in S.ORDERS:
        r = D[("mPIE", False, s_)]["rows"]
        ax[0].semilogy(r[:, 0], r[:, 1] / D[("mPIE", False, s_)]["power"], color=COL["mPIE"], alpha=0.85, lw=1.2)
        for alg, cc, ls in (("mPIE", COL["mPIE"], "-"), ("mom_obj", COL["rPIE"], "-"), ("mom_prb", COL["rPIE_warm"], "-"),
                            ("rPIE", COL["DM"], "--")):
            r = D[(alg, True, s_)]["rows"]
            ax[1].semilogy(r[:, 0], r[:, 6], ls, color=cc, alpha=0.85, lw=1.2, label=L["d_" + alg] if s_ == 0 else None)
    ax[0].set_title(L["nopc"], fontsize=10, **Fk); ax[0].set_ylabel(L["power"], **Fk)
    ax[1].set_title(L["pc"], fontsize=10, **Fk); ax[1].set_ylabel(L["unlit"], **Fk)
    ax[1].legend(prop=LG, fontsize=8, loc="lower left")
    for a_ in ax:
        a_.set_xlabel(L["iters"], **Fk); a_.grid(alpha=0.3, which="both")
    cap(fig, "fig3"); save(fig, "fig3-drift-divergence.png")

    # 그림4: 매개변수 지도
    PM = R["param_map"]
    etas, gammas = PM["etas"], PM["gammas"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.1), sharey=True)
    for a, case in zip(ax, S.PROBES):
        Z = np.array([[np.median([t[0] for t in PM["res"][(case, e, g)]]) for g in gammas] for e in etas])
        ok = np.array([[sum(t[0] < 0.05 for t in PM["res"][(case, e, g)]) for g in gammas] for e in etas])
        im = a.imshow(np.log10(np.clip(Z, 1e-5, 1)), cmap="viridis_r", vmin=-4.5, vmax=0, origin="lower", aspect="auto")
        for i in range(len(etas)):
            for j in range(len(gammas)):
                a.text(j, i, str(ok[i, j]), ha="center", va="center", fontsize=9, color="w" if Z[i, j] > 0.05 else "k")
        a.add_patch(plt.Rectangle((gammas.index(0.2) - 0.5, etas.index(0.9) - 0.5), 1, 1, fill=False, ec="r", lw=2))
        a.set_xticks(range(len(gammas)), [str(g) for g in gammas]); a.set_yticks(range(len(etas)), [str(e) for e in etas])
        a.set_title(L[case], fontsize=10, **Fk); a.set_xlabel(L["gamma"], **Fk)
    ax[0].set_ylabel(L["eta"], **Fk)
    fig.subplots_adjust(left=0.05, right=0.9, bottom=0.16, top=0.9, wspace=0.08)
    cb = fig.colorbar(im, cax=fig.add_axes([0.915, 0.16, 0.012, 0.74])); cb.set_label("log10 " + L["err"], **Fk)
    save(fig, "fig4-parameter-map.png")

    # 그림5: 세 프로브, 다섯 방법
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2), sharey=True)
    for a, case in zip(ax, S.PROBES):
        for alg in ("rPIE", "rPIE_warm", "mPIE", "mPIE_ptylab", "DM"):
            band(a, case, alg, 0, L[alg])
        a.set_title(L[case], fontsize=10, **Fk); a.set_xlabel(L["iters"], **Fk); a.grid(alpha=0.3, which="both")
    ax[0].set_ylabel(L["err"], **Fk); ax[0].legend(prop=LG, fontsize=8, loc="lower left")
    cap(fig, "fig5"); save(fig, "fig5-three-probes.png")

    # 그림6: 확산판 정체
    st = R["stall"]
    plain, known = st[("plain", 0)], st[("known", 0)]
    P, p = plain["P"], plain["probe"]
    p = p * np.vdot(p, P) / np.vdot(p, p)
    m = plain["mask"]; ys, xs = np.where(m); bb = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    fig, ax = plt.subplots(1, 4, figsize=(13, 3.7))
    for a, img, t in ((ax[0], np.abs(P), L["t_p"]), (ax[1], np.abs(p), L["found"])):
        a.imshow(img, cmap="gray")
        cy, cx = S.centroid(P)
        a.plot(cx, cy, "r+", ms=10, mew=1.5)
        if t == L["found"]:
            cy2, cx2 = S.centroid(p); a.plot(cx2, cy2, "c+", ms=10, mew=1.5)
            a.contour(np.abs(P), levels=[0.5 * np.abs(P).max()], colors="r", linewidths=0.8)
        a.set_title(t, fontsize=9.5, **Fk); a.set_axis_off()
    for a, r, key in ((ax[2], plain, "err_plain"), (ax[3], known, "err_known")):
        a.imshow(np.where(m, np.abs(r["rec"] - r["obj"]), np.nan)[bb], cmap="magma", vmin=0, vmax=0.1)
        a.set_title(f"{L[key]}: {r['err']:.2g}", fontsize=9.5, **Fk); a.set_axis_off()
    cap(fig, "fig6", rect=(0, 0, 1, 0.97)); save(fig, "fig6-diffuser-stall.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")

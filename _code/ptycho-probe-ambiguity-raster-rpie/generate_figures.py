"""계산 이미징과 타이코그래피 4편 그림 6개 생성 (한국어판·영문판)."""
import os

import matplotlib.pyplot as plt
import numpy as np

from ptycho4_core import N, pinhole_probe, make_object, scan, object_shape, intensities, probe_d90
import studies

SLUG = "ptycho-probe-ambiguity-raster-rpie"
BASE = os.path.join(os.path.dirname(__file__), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
C_TRUE, C_A, C_B, C_C = "#555555", "#c8553d", "#2f6690", "#6a994e"

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "v_true": "원래 (O, P)", "v_scale": "척도 교환 (cO, P/c)", "v_ramp": "위상 기울기",
        "v_lat": "격자 주기 함수 f (규칙 격자에서만)",
        "o_ph": "시편 위상", "p_amp": "프로브 진폭", "p_ph": "프로브 위상",
        "fig1": "그림1. 넷 다 모든 회절 패턴이 똑같다. 마지막 것은 주사가 규칙 격자일 때만 그렇다",
        "regular": "규칙 격자", "jitter": "±1화소 흔든 격자", "fermat": "페르마 나선",
        "rec": "복원 진폭", "errmap": "오차 |O - O_true|", "spec": "오차 스펙트럼 (로그)",
        "fig2": "그림2. 규칙 격자에서만 오차가 격자 모양으로 남고, 그 스펙트럼은 격자 주파수에 몰린다",
        "pn": "프로브 세기 |P| / |P|max", "w": "시편에 반영하는 비율 w",
        "l_epie": "ePIE (rPIE $\\alpha = 1$)", "l_pie": "PIE ($\\alpha = 10^{-4}$)",
        "l_r25": "rPIE $\\alpha = 0.25$", "l_r05": "rPIE $\\alpha = 0.05$",
        "fig3": "그림3. ePIE는 밝은 곳 말고는 거의 고치지 않는다. rPIE는 α로 그 곡선을 끌어올린다",
        "iters": "반복 횟수", "err": "상대 오차",
        "pinhole": "핀홀 프로브", "defocus": "초점이 어긋난 집속 프로브", "diffuser": "확산판 프로브",
        "fig4": "그림4. rPIE는 ePIE보다 열 배 이상 빨리 수렴한다. 확산판 프로브는 둘 다 풀지 못한다",
        "t_p": "참 프로브", "init": "초기 추정 (초점에 맞춤)", "e_found": "ePIE 200회", "r_found": "rPIE 200회",
        "amp": "진폭", "phase": "위상",
        "fig5": "그림5. 초점에 맞춘 추정에서 출발해 ePIE와 rPIE 모두 휜 파면을 찾는다. rPIE가 훨씬 빠르다",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "v_true": "Original (O, P)", "v_scale": "Scale exchange (cO, P/c)", "v_ramp": "Phase ramp",
        "v_lat": "Lattice-periodic f (regular grid only)",
        "o_ph": "object phase", "p_amp": "probe amplitude", "p_ph": "probe phase",
        "fig1": "Fig 1. All four give identical diffraction patterns; the last only on a regular grid",
        "regular": "Regular grid", "jitter": "Grid jittered by ±1 px", "fermat": "Fermat spiral",
        "rec": "Reconstructed amplitude", "errmap": "Error |O - O_true|", "spec": "Error spectrum (log)",
        "fig2": "Fig 2. Only on the regular grid does the error form a grid pattern, concentrated at grid frequencies",
        "pn": "Probe brightness |P| / |P|max", "w": "Fraction applied to the object, w",
        "l_epie": "ePIE (rPIE $\\alpha = 1$)", "l_pie": "PIE ($\\alpha = 10^{-4}$)",
        "l_r25": "rPIE $\\alpha = 0.25$", "l_r05": "rPIE $\\alpha = 0.05$",
        "fig3": "Fig 3. ePIE barely corrects anything but the bright region; rPIE lifts the curve through α",
        "iters": "Iterations", "err": "Relative error",
        "pinhole": "Pinhole probe", "defocus": "Defocused focused probe", "diffuser": "Diffuser probe",
        "fig4": "Fig 4. rPIE converges more than ten times faster than ePIE; neither solves the diffuser probe",
        "t_p": "True probe", "init": "Initial guess (in focus)", "e_found": "ePIE, 200 it.", "r_found": "rPIE, 200 it.",
        "amp": "amplitude", "phase": "phase",
        "fig5": "Fig 5. From an in-focus guess, both ePIE and rPIE recover the curved wavefront; rPIE is far faster",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
SC = studies.scans()
AL = studies.algorithms()

# 그림1 의 네 쌍
P = pinhole_probe()
pos_reg = scan("regular")
shp = object_shape(pos_reg)
OBJ = make_object(shp)
Y, X = np.indices(shp)
yy, xx = np.indices((N, N))
y0, x0 = pos_reg[0]
c = 0.6 * np.exp(1.0j)
g = (0.10, -0.06)
f_o = (1 + 0.3 * np.cos(2 * np.pi * (X - x0) / 8)) * np.exp(0.6j * np.sin(2 * np.pi * (Y - y0) / 8))
f_p = (1 + 0.3 * np.cos(2 * np.pi * xx / 8)) * np.exp(0.6j * np.sin(2 * np.pi * yy / 8))
VARIANTS = [("v_true", OBJ, P), ("v_scale", OBJ * c, P / c),
            ("v_ramp", OBJ * np.exp(1j * (g[0] * X + g[1] * Y)), P * np.exp(-1j * (g[0] * xx + g[1] * yy))),
            ("v_lat", OBJ * f_o, P / f_p)]
I0 = intensities(OBJ, P, pos_reg)
SAME = [float(np.abs(intensities(o, p, pos_reg) - I0).max() / I0.max()) for _, o, p in VARIANTS]


def med_range(case, name, its):
    vals = np.array([[AL[(case, name, s)]["trace"][k] for k in its] for s in (0, 1, 2)])
    return np.median(vals, 0), vals.min(0), vals.max(0)


SUMMARY = {"same": SAME, "D90": probe_d90(P),
           "scans": {k: dict({a: b for a, b in v.items() if isinstance(b, float)}, runs=v["runs"]) for k, v in SC.items()},
           "algos": {k: (v["trace"], v["probe_err"]) for k, v in AL.items() if k[1] != "_truth"}}


# ── 그리기 ───────────────────────────────────────────────────
def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.06, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **Fk)
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    def img(ax, a, title, cmap="gray", **kw):
        ax.imshow(a, cmap=cmap, **kw); ax.set_title(title, fontsize=8.5, **Fk); ax.set_axis_off()

    # 그림1: 모호성
    fig, ax = plt.subplots(3, 4, figsize=(11, 8.4))
    crop = (slice(20, 100), slice(20, 100))
    for j, (key, o, p) in enumerate(VARIANTS):
        img(ax[0, j], np.angle(o[crop]), f"{L[key]}\n{L['o_ph']}", cmap="twilight", vmin=-np.pi, vmax=np.pi)
        img(ax[1, j], np.abs(p), L["p_amp"], vmin=0)
        img(ax[2, j], np.where(np.abs(P) > 0.05 * np.abs(P).max(), np.angle(p), np.nan), L["p_ph"],
            cmap="twilight", vmin=-np.pi, vmax=np.pi)
    cap(fig, "fig1", rect=(0, 0.04, 1, 1)); save(fig, "fig1-ambiguities.png", dpi=120)

    # 그림2: 규칙 격자
    fig, ax = plt.subplots(3, 4, figsize=(12.5, 9.2), gridspec_kw={"width_ratios": [0.9, 1, 1, 1]})
    for i, kind in enumerate(("regular", "jitter", "fermat")):
        r = SC[kind]
        pos = r["pos"]
        ax[i, 0].plot(pos[:, 1], -pos[:, 0], ".", ms=3, color=C_B)
        ax[i, 0].set_aspect("equal"); ax[i, 0].set_axis_off()
        ax[i, 0].set_title(L[kind], fontsize=9.5, **Fk)
        m = r["mask"]; yy_, xx_ = np.where(m)
        bb = (slice(yy_.min(), yy_.max() + 1), slice(xx_.min(), xx_.max() + 1))
        img(ax[i, 1], np.abs(np.where(m, r["rec"], np.nan)[bb]), f"{L['rec']} ({r['err']:.2g})", vmin=0, vmax=1)
        img(ax[i, 2], np.abs(r["err_crop"]), L["errmap"], cmap="magma", vmin=0, vmax=0.12)
        spec = np.fft.fftshift(np.abs(np.fft.fft2(r["err_crop"])) ** 2)
        img(ax[i, 3], np.log10(spec + 1e-12 * spec.max()), f"{L['spec']}: {r['grid_frac']:.0%}", cmap="viridis",
            vmin=np.log10(spec.max()) - 6, vmax=np.log10(spec.max()))
    cap(fig, "fig2", rect=(0, 0.03, 1, 1)); save(fig, "fig2-raster-grid.png", dpi=120)

    # 그림3: 가중 곡선
    x = np.linspace(1e-4, 1, 400)
    fig, ax = plt.subplots(figsize=(7.5, 4))
    ax.plot(x, x**2, color=C_A, lw=1.6, label=L["l_epie"])
    ax.plot(x, x * x**2 / (x**2 + 1e-4), "--", color=C_TRUE, lw=1.3, label=L["l_pie"])
    for a, cc, key in ((0.25, C_C, "l_r25"), (0.05, C_B, "l_r05")):
        ax.plot(x, x**2 / ((1 - a) * x**2 + a), color=cc, lw=1.6, label=L[key])
    ax.set_xlabel(L["pn"], **Fk); ax.set_ylabel(L["w"], **Fk); ax.legend(prop=LG, fontsize=8.5)
    ax.grid(alpha=0.3)
    cap(fig, "fig3", rect=(0, 0.07, 1, 1)); save(fig, "fig3-update-weights.png")

    # 그림4: 수렴
    its = studies.REC
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for a_, case in zip(ax, ("pinhole", "defocus", "diffuser")):
        for name, cc in (("ePIE", C_A), ("rPIE", C_B)):
            med, lo, hi = med_range(case, name, its)
            a_.loglog(its, med, "o-", color=cc, ms=4, label=name)
            a_.fill_between(its, lo, hi, color=cc, alpha=0.2)
        a_.set_title(L[case], fontsize=10, **Fk); a_.set_xlabel(L["iters"], **Fk); a_.grid(alpha=0.3, which="both")
    ax[0].set_ylabel(L["err"], **Fk); ax[0].legend(prop=LG, fontsize=8.5)
    cap(fig, "fig4", rect=(0, 0.07, 1, 1)); save(fig, "fig4-convergence.png")

    # 그림5: 초점 어긋난 프로브
    tr = AL[("defocus", "_truth")]
    Pt, P0 = tr["probe"], tr["probe0"]
    from ptycho4_core import align
    pe = align(AL[("defocus", "ePIE", 0)]["probe"], Pt)
    pr = align(AL[("defocus", "rPIE", 0)]["probe"], Pt)
    fig, ax = plt.subplots(2, 4, figsize=(11, 6.4))
    msk = np.abs(Pt) > 0.05 * np.abs(Pt).max()
    for j, (pp, key) in enumerate(((Pt, "t_p"), (P0, "init"), (pe, "e_found"), (pr, "r_found"))):
        img(ax[0, j], np.abs(pp), f"{L[key]}: {L['amp']}")
        img(ax[1, j], np.where(msk, np.angle(pp), np.nan), L["phase"], cmap="twilight", vmin=-np.pi, vmax=np.pi)
    cap(fig, "fig5", rect=(0, 0.05, 1, 1)); fig.subplots_adjust(hspace=0.18); save(fig, "fig5-defocus-probe.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 5개 저장: {L['dir']}")
    print("\n본문 수치")
    for k, v in SUMMARY.items():
        print(f"  {k}: {v}")

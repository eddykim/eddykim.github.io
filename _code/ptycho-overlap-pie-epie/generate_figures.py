"""계산 이미징과 타이코그래피 3편 그림 7개 생성 (한국어판·영문판)."""
import os

import matplotlib.pyplot as plt
import numpy as np

from ptycho_core import fft2c, pie_weight, align, intensities
import studies

SLUG = "ptycho-overlap-pie-epie"
BASE = os.path.join(os.path.dirname(__file__), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
C_TRUE, C_A, C_B, C_C, C_D = "#555555", "#c8553d", "#2f6690", "#6a994e", "#e09f3e"

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "obj_t": "시편 진폭과 주사 위치 (원 = 프로브 90 % 지름)", "pat": "위치 {} 의 회절 패턴",
        "fig1": "그림1. 프로브를 겹치게 옮겨 가며 위치마다 회절 패턴을 한 장씩 찍는다",
        "pamp": "프로브 진폭", "pph": "프로브 위상", "pos_ax": "프로브 중심에서의 거리 (화소)",
        "naive": "그냥 나눌 때의 증폭 1/|P| (최댓값 1로 맞춤)", "frac": "PIE가 반영하는 비율",
        "pnorm": "|P| / |P|max",
        "fig2": "그림2. PIE는 출사파를 프로브로 나누는 대신, 프로브가 밝은 만큼만 고친다",
        "iters": "반복 횟수", "err": "상대 오차",
        "l_pie": "PIE, 정확한 프로브", "l_wrong": "PIE, 틀린 프로브", "l_epie": "ePIE, 틀린 프로브에서 출발",
        "l_probe": "ePIE의 프로브 오차",
        "t_amp": "참 진폭", "t_ph": "참 위상", "r_amp": "PIE 진폭 (100회)", "r_ph": "PIE 위상 (100회)",
        "fig3": "그림3. 프로브를 알면 PIE는 수십 번 만에 수렴하고, 프로브가 틀리면 멈춘다",
        "ov_ax": "중첩률 (프로브 90 % 지름 기준)", "l_epie_o": "ePIE, 시편", "l_epie_p": "ePIE, 프로브",
        "rec_lo": "PIE, 중첩 {:.0%}", "rec_ep": "ePIE, 중첩 {:.0%}",
        "fig4": "그림4. 프로브를 알면 절반만 겹쳐도 되지만, 프로브까지 찾으려면 더 많이 겹쳐야 한다",
        "init": "초기 프로브 (원판)", "found": "ePIE가 찾은 프로브", "true": "참 프로브",
        "amp": "진폭", "phase": "위상",
        "fig5": "그림5. ePIE는 크기만 대충 맞는 원판에서 출발해 프로브의 무늬와 위상까지 찾는다",
        "dose_ax": "조명된 화소당 광자 수", "l_cdi": "단일 패턴 CDI (2편의 HIO)",
        "c_cdi": "단일 패턴 CDI", "c_pie": "PIE", "photons": "화소당 {} 광자",
        "fig6": "그림6. 같은 선량이면 겹쳐 찍은 쪽의 오차가 작다 (PIE 기준 2~4배)",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "obj_t": "Object amplitude and scan positions (circles = 90 % probe diameter)", "pat": "Pattern at position {}",
        "fig1": "Fig 1. The probe is stepped with overlap, and one diffraction pattern is recorded at each position",
        "pamp": "Probe amplitude", "pph": "Probe phase", "pos_ax": "Distance from probe centre (pixels)",
        "naive": "Gain of plain division 1/|P| (scaled to max 1)", "frac": "Fraction applied by PIE",
        "pnorm": "|P| / |P|max",
        "fig2": "Fig 2. Instead of dividing by the probe, PIE corrects only as much as the probe is bright",
        "iters": "Iterations", "err": "Relative error",
        "l_pie": "PIE, correct probe", "l_wrong": "PIE, wrong probe", "l_epie": "ePIE, from the wrong probe",
        "l_probe": "ePIE probe error",
        "t_amp": "True amplitude", "t_ph": "True phase", "r_amp": "PIE amplitude (100 it.)", "r_ph": "PIE phase (100 it.)",
        "fig3": "Fig 3. With a known probe PIE converges in tens of iterations; with a wrong probe it stalls",
        "ov_ax": "Overlap (by 90 % probe diameter)", "l_epie_o": "ePIE, object", "l_epie_p": "ePIE, probe",
        "rec_lo": "PIE, {:.0%} overlap", "rec_ep": "ePIE, {:.0%} overlap",
        "fig4": "Fig 4. With a known probe half overlap suffices; recovering the probe as well needs more",
        "init": "Initial probe (disc)", "found": "Probe found by ePIE", "true": "True probe",
        "amp": "amplitude", "phase": "phase",
        "fig5": "Fig 5. Starting from a disc of roughly the right size, ePIE recovers the probe's fringes and phase",
        "dose_ax": "Photons per illuminated pixel", "l_cdi": "Single-pattern CDI (HIO, Post 2)",
        "c_cdi": "Single-pattern CDI", "c_pie": "PIE", "photons": "{} photons/pixel",
        "fig6": "Fig 6. At equal dose, overlapping measurements give a smaller error (2-4 times, for PIE)",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
P = studies.P
CV = studies.convergence()
OV = studies.overlap()
NZ = studies.noise()
POS, OBJ, MASK = CV["pos"], CV["obj"], CV["mask"]
I_DEMO = intensities(OBJ, P, POS)
ys, xs = np.where(MASK)
BB = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))


def show_ptycho(rec, obj, mask):
    a = align(rec, obj, mask)
    return np.where(mask, a, np.nan)[BB]


def show_cdi(rec, box):
    F = lambda x: np.fft.fft2(x, norm="ortho")
    iF = lambda x: np.fft.ifft2(x, norm="ortho")
    best = None
    for cand in (rec, np.conj(np.roll(np.flip(rec), 1, axis=(0, 1)))):
        xc = iF(F(box) * np.conj(F(cand)))
        sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
        c2 = np.roll(cand, sh, axis=(0, 1)); a = np.vdot(c2, box) / np.vdot(c2, c2)
        e = np.linalg.norm(a * c2 - box)
        if best is None or e < best[0]:
            best = (e, a * c2)
    return best[1][40:88, 40:88]


SUMMARY = {"pie": CV["pie"][0], "pie_wrong": CV["pie_wrong"][0], "epie": CV["epie"][0], "epie_probe": CV["epie"][2],
           "overlap": {s: {k: v for k, v in r.items() if not k.startswith("rec") and k not in ("obj", "mask")}
                       for s, r in OV.items()},
           "noise": {d: {k: v for k, v in r.items() if k in ("pie", "epie", "cdi", "per_pattern")} for d, r in NZ.items()},
           "D90": studies.D90, "n_pos_demo": len(POS)}


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
        ax.imshow(a, cmap=cmap, **kw); ax.set_title(title, fontsize=9, **Fk); ax.set_axis_off()

    n = P.shape[0]
    # 그림1: 배치
    fig = plt.figure(figsize=(12, 4.4))
    gs = fig.add_gridspec(2, 4, width_ratios=[2, 1, 1, 1])
    a0 = fig.add_subplot(gs[:, 0])
    a0.imshow(np.abs(OBJ), cmap="gray")
    for j, (y, x) in enumerate(POS):
        hl = j in (0, 1, 13)
        a0.add_patch(plt.Circle((x + n / 2, y + n / 2), studies.D90 / 2, fill=False,
                                ec=C_D if hl else C_B, lw=1.6 if hl else 0.5, alpha=1 if hl else 0.6))
    a0.set_title(L["obj_t"], fontsize=9, **Fk); a0.set_axis_off()
    for k, j in enumerate((0, 1, 13)):
        a = fig.add_subplot(gs[:, k + 1])
        img(a, np.log10(I_DEMO[j] + 1e-6 * I_DEMO[j].max()), L["pat"].format(j + 1), cmap="magma")
    cap(fig, "fig1"); save(fig, "fig1-setup.png")

    # 그림2: 프로브와 PIE 가중
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8), gridspec_kw={"width_ratios": [1, 1, 1.8]})
    img(ax[0], np.abs(P), L["pamp"])
    img(ax[1], np.where(np.abs(P) > 0.05 * np.abs(P).max(), np.angle(P), np.nan), L["pph"], cmap="twilight",
        vmin=-np.pi, vmax=np.pi)
    row = P[n // 2, n // 2:]
    r = np.arange(row.size)
    pn = np.abs(row) / np.abs(P).max()
    ax[2].semilogy(r, pn, color=C_TRUE, lw=1.4, label=L["pnorm"])
    ax[2].semilogy(r, np.abs(pie_weight(P)[n // 2, n // 2:] * row), "--", color=C_B, lw=1.6, label=L["frac"])
    gain = 1 / np.maximum(np.abs(row), 1e-30)
    ax[2].semilogy(r, gain / gain.min(), color=C_A, lw=1.2, label=L["naive"])
    ax[2].set_ylim(1e-3, 3e2); ax[2].set_xlabel(L["pos_ax"], **Fk)
    ax[2].legend(prop=LG, fontsize=7.5, loc="lower left"); ax[2].grid(alpha=0.3, which="both")
    cap(fig, "fig2", rect=(0, 0.07, 1, 0.94)); save(fig, "fig2-pie-weight.png")

    # 그림3: 수렴
    fig = plt.figure(figsize=(12.5, 4.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.7, 1, 1])
    a0 = fig.add_subplot(gs[:, 0])
    for key, c, ls, lab in (("pie", C_B, "o-", "l_pie"), ("pie_wrong", C_A, "s--", "l_wrong"),
                            ("epie", C_C, "^-", "l_epie")):
        tr = CV[key][0]
        a0.loglog(list(tr), list(tr.values()), ls, color=c, label=L[lab], ms=5)
    trp = CV["epie"][2]
    a0.loglog(list(trp), list(trp.values()), ":", color=C_C, label=L["l_probe"])
    a0.set_xlabel(L["iters"], **Fk); a0.set_ylabel(L["err"], **Fk); a0.legend(prop=LG, fontsize=8)
    a0.grid(alpha=0.3, which="both")
    rec = show_ptycho(CV["pie"][1], OBJ, MASK)
    tru = np.where(MASK, OBJ, np.nan)[BB]
    img(fig.add_subplot(gs[0, 1]), np.abs(tru), L["t_amp"], vmin=0, vmax=1)
    img(fig.add_subplot(gs[0, 2]), np.abs(rec), L["r_amp"], vmin=0, vmax=1)
    img(fig.add_subplot(gs[1, 1]), np.angle(tru), L["t_ph"], cmap="twilight", vmin=-np.pi, vmax=np.pi)
    img(fig.add_subplot(gs[1, 2]), np.angle(rec), L["r_ph"], cmap="twilight", vmin=-np.pi, vmax=np.pi)
    cap(fig, "fig3", rect=(0, 0.06, 1, 1)); save(fig, "fig3-pie-convergence.png")

    # 그림4: 중첩률
    fig = plt.figure(figsize=(12.5, 4.3))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.8, 1, 1, 1])
    a0 = fig.add_subplot(gs[0])
    steps = sorted(OV, reverse=True)
    ov = [OV[s]["overlap"] for s in steps]
    a0.semilogy(ov, [OV[s]["pie"] for s in steps], "o-", color=C_B, label="PIE")
    a0.semilogy(ov, [OV[s]["epie"] for s in steps], "^-", color=C_C, label=L["l_epie_o"])
    a0.semilogy(ov, [OV[s]["probe"] for s in steps], "^:", color=C_C, mfc="none", label=L["l_epie_p"])
    a0.set_xlabel(L["ov_ax"], **Fk); a0.set_ylabel(L["err"], **Fk); a0.legend(prop=LG, fontsize=8)
    a0.grid(alpha=0.3, which="both")
    for k, (s, key, lab) in enumerate(((24, "rec_pie", "rec_lo"), (12, "rec_epie", "rec_ep"), (6, "rec_epie", "rec_ep"))):
        r_ = OV[s]
        m = r_["mask"]; yy, xx = np.where(m)
        bb = (slice(yy.min(), yy.max() + 1), slice(xx.min(), xx.max() + 1))
        a_ = align(r_[key], r_["obj"], m)
        img(fig.add_subplot(gs[k + 1]), np.abs(np.where(m, a_, np.nan)[bb]), L[lab].format(r_["overlap"]),
            vmin=0, vmax=1)
    cap(fig, "fig4", rect=(0, 0.06, 1, 1)); save(fig, "fig4-overlap.png")

    # 그림5: 프로브 복원
    fig, ax = plt.subplots(2, 3, figsize=(9, 6.2))
    found = align(CV["epie"][3], P)
    for k, (pp, key) in enumerate(((studies.WRONG_PROBE, "init"), (found, "found"), (P, "true"))):
        img(ax[0, k], np.abs(pp), f"{L[key]}: {L['amp']}")
        img(ax[1, k], np.where(np.abs(pp) > 0.05 * np.abs(pp).max(), np.angle(pp), np.nan), L["phase"],
            cmap="twilight", vmin=-np.pi, vmax=np.pi)
    cap(fig, "fig5", rect=(0, 0.05, 1, 1)); save(fig, "fig5-epie-probe.png")

    # 그림6: 선량
    fig = plt.figure(figsize=(12, 4.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.6, 1, 1])
    a0 = fig.add_subplot(gs[0])
    ds = sorted(NZ)
    a0.loglog(ds, [NZ[d]["cdi"] for d in ds], "s--", color=C_A, label=L["l_cdi"])
    a0.loglog(ds, [NZ[d]["pie"] for d in ds], "o-", color=C_B, label="PIE")
    a0.loglog(ds, [NZ[d]["epie"] for d in ds], "^-", color=C_C, label="ePIE")
    a0.set_xlabel(L["dose_ax"], **Fk); a0.set_ylabel(L["err"], **Fk); a0.legend(prop=LG, fontsize=8)
    a0.grid(alpha=0.3, which="both")
    d = 1255
    r_ = NZ[d]
    img(fig.add_subplot(gs[1]), np.abs(show_cdi(r_["rec_cdi"], r_["box"])),
        f"{L['c_cdi']} ({r_['cdi']:.2f})\n{L['photons'].format(d)}", vmin=0, vmax=1)
    c = r_["obj"].shape[0] // 2
    sub = align(r_["rec_pie"], r_["obj"], r_["mask"])[c - 24:c + 24, c - 24:c + 24]
    img(fig.add_subplot(gs[2]), np.abs(sub), f"{L['c_pie']} ({r_['pie']:.2f})\n{L['photons'].format(d)}",
        vmin=0, vmax=1)
    cap(fig, "fig6", rect=(0, 0.06, 1, 0.92)); save(fig, "fig6-dose.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")
    print("\n본문 수치")
    for k, v in SUMMARY.items():
        print(f"  {k}: {v}")

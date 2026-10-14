"""계산 이미징과 타이코그래피 10편 계산 그림 6개 생성 (한국어판·영문판). 수치는 studies10.py 의 캐시를 읽는다.

그림1 은 Tanksalvala et al. 2021 의 그림(CC BY 4.0, ext-tanksalvala2021-setup.jpg)이라 여기서 만들지 않는다.
위상은 본문의 규약(시간 인자 exp(+j w t), n = 1 - delta - j beta)에 맞춰 코드 값의 켤레로 그린다.
"""
import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse

import ptycho10_core as C
import studies10 as S

SLUG = "ptycho-euv-imaging-reflectometry"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
COL = {"amp": "#c8553d", "phase": "#2f6690", "both": "#333333", "ideal": "#6a994e", "corrected": "#2f6690",
       "uncorrected": "#c8553d", 0: "#2f6690", 1: "#6a994e", 2: "#c8553d"}
CURVE_COL = ["#555555", "#9fc5e8", "#6a994e", "#9a6fb0"]

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "sample": "시편", "det": "검출기", "inc": "입사빔", "spec": "정반사", "foot": "조명 자국은 입사 방향으로 1/sin θ배 길다",
        "qx": r"$q_x / k$ (입사 방향)", "qy": r"$q_y / k$", "grid_det": "검출기 화소가 받는 q", "grid_uni": "복원에 쓰는 균일 격자",
        "g": "입사각 (표면에서 잰 각, 도)", "R": "반사율 (%)", "phi": "반사 위상 (도)", "band": "Tanksalvala 2021\n측정 구간",
        "curves": ["Si (자연 산화막 1.5 nm)", "Si", "SiO$_2$", "Si$_3$N$_4$"],
        "true": "참 시편", "corr": "보정함", "uncorr": "보정하지 않음 (선형 사상)", "na": "검출기 수치구경",
        "err": "복원 오차", "ideal": "균일 격자에서 측정 (기준)",
        "os": "검출기 과표본 배수", "bias": "영역 위상 차의 편향 (rad)", "phase_img": "위상 (rad)",
        "ratio_amp": "반사율 비 |r구조 / r기판|", "ratio_ph": "위상 차 (rad)",
        "model": "모형 (참 매개변수)", "meas": "복원 (노이즈 표본 3개)", "deg": "입사각 (도)",
        "img21": "21도, 기판 기준 위상", "img25": "25도, 기판 기준 위상",
        "narrow": "좁은 각도 (21–25도, 5개)", "wide": "넓은 각도 (15–35도, 5개)",
        "dN": "질화막 두께 (nm)", "rho": "질화막 밀도 배율", "amp": "반사율 비만", "phase": "위상 차만", "both": "둘 다",
        "fits": "맞춤 결과", "truth": "참값",
        "mat": ["Si$_3$N$_4$", "SiO$_2$", "Al"], "vac": "진공", "thick": "두께 (nm)",
        "la": r"$-\ln|T|$", "ph": r"위상 $\varphi$ (rad)", "cal": "진공 기준으로 맞춤", "raw": "기준 없이 (임의의 복소 상수)",
        "slope": "기울기 = δ/β (표값)", "layout": "물질 배치와 두께",
        "fig3": "그림3. 반사 타이코그래피의 기하(왼쪽)와 22도에서 검출기 화소가 받는 표면 공간주파수(오른쪽). 회색 사각형은 복원에 쓰는 균일 격자의 범위다",
        "fig2": "그림2. 29.3 nm에서 계산한 s 편광 반사율과 반사 위상. 21–25도 구간에서 물질마다 값이 크게 갈린다",
        "fig4": "그림4. 22도 반사 타이코그래피의 복원 위상(위). 보정하지 않으면 무너진다. 아래 왼쪽은 보정의 세기 보간이 영역 위상 차에 남기는 편향으로, 검출기를 촘촘히 둘수록 줄어든다. 아래 오른쪽은 수치구경에 따른 복원 오차다",
        "fig5": "그림5. 각도별 복원에서 기판에 대한 구조 영역의 반사율 비와 위상 차. 점은 노이즈 표본 3개, 선은 참 매개변수의 모형이다",
        "fig6": "그림6. 질화막 두께와 밀도의 1σ 신뢰 타원(정보 행렬)과 실제 맞춤 결과. 좁은 각도에서 반사율 비만 쓰면 두 값이 묶인다",
        "fig7": "그림7. 13.5 nm 투과 복원의 화소별 위상과 로그 진폭. 진공으로 맞추면 물질마다 원점을 지나는 직선에 놓이고(가운데), 기준이 없으면 흩어진다(오른쪽)",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "sample": "Sample", "det": "Detector", "inc": "Incident", "spec": "Specular", "foot": "Footprint stretched by 1/sin θ",
        "qx": r"$q_x / k$ (incidence direction)", "qy": r"$q_y / k$", "grid_det": "q received by detector pixels",
        "grid_uni": "Uniform grid used in reconstruction",
        "g": "Grazing angle (deg)", "R": "Reflectivity (%)", "phi": "Reflection phase (deg)", "band": "Tanksalvala 2021\nrange",
        "curves": ["Si (1.5 nm native oxide)", "Si", "SiO$_2$", "Si$_3$N$_4$"],
        "true": "True sample", "corr": "Corrected", "uncorr": "Uncorrected (linear map)", "na": "Detector NA",
        "err": "Reconstruction error", "ideal": "Measured on uniform grid (reference)",
        "os": "Detector oversampling factor", "bias": "Bias of region phase step (rad)", "phase_img": "Phase (rad)",
        "ratio_amp": "Reflectance ratio |r_struct / r_sub|", "ratio_ph": "Phase step (rad)",
        "model": "Model (true parameters)", "meas": "Reconstructed (3 noise draws)", "deg": "Grazing angle (deg)",
        "img21": "21°, phase relative to substrate", "img25": "25°, phase relative to substrate",
        "narrow": "Narrow range (21–25°, 5 angles)", "wide": "Wide range (15–35°, 5 angles)",
        "dN": "Nitride thickness (nm)", "rho": "Nitride density factor", "amp": "Amplitude ratio only",
        "phase": "Phase step only", "both": "Both", "fits": "Fits", "truth": "Truth",
        "mat": ["Si$_3$N$_4$", "SiO$_2$", "Al"], "vac": "Vacuum", "thick": "Thickness (nm)",
        "la": r"$-\ln|T|$", "ph": r"Phase $\varphi$ (rad)", "cal": "Calibrated to vacuum",
        "raw": "No reference (arbitrary complex constant)", "slope": "Slope = δ/β (tabulated)", "layout": "Materials and thickness",
        "fig3": "Fig 3. Geometry of reflection ptychography (left) and the surface spatial frequencies received by the detector pixels at 22° (right). The grey square is the extent of the uniform grid used in reconstruction",
        "fig2": "Fig 2. Computed s-polarized reflectivity and reflection phase at 29.3 nm. Between 21° and 25° the values separate strongly by material",
        "fig4": "Fig 4. Reconstructed phase from reflection ptychography at 22° (top); without correction it fails. Bottom left: bias that intensity interpolation in the correction leaves in region phase steps, shrinking as the detector is sampled more finely. Bottom right: reconstruction error versus NA",
        "fig5": "Fig 5. Reflectance ratio and phase step of the structure relative to the substrate, from reconstructions at each angle. Dots are 3 noise draws; lines are the model with the true parameters",
        "fig6": "Fig 6. 1σ confidence ellipses (information matrix) for nitride thickness and density, with the actual fits. Over a narrow angle range the amplitude ratio alone ties the two together",
        "fig7": "Fig 7. Pixel-wise phase and log amplitude of a 13.5 nm transmission reconstruction. Calibrated to vacuum, each material lies on a line through the origin (centre); without a reference the points scatter (right)",
    },
}

R = S.load()


def ellipse(ax, cov, centre, color, label, ls="-"):
    w, v = np.linalg.eigh(cov)
    ang = np.degrees(np.arctan2(v[1, -1], v[0, -1]))
    ax.add_patch(Ellipse(centre, 2 * np.sqrt(w[-1]), 2 * np.sqrt(w[0]), angle=ang, fill=False, color=color, lw=1.8,
                         ls=ls, label=label))


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.07, 1, 1)):
        t = fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **Fk)
        fig.canvas.draw()
        if t.get_window_extent().width > 0.98 * fig.bbox.width:  # 캡션이 폭을 넘으면 줄바꿈
            text = textwrap.fill(L[key], width=int(len(L[key]) * 0.95 * fig.bbox.width / t.get_window_extent().width))
            t.set_text(text); t.set_va("bottom"); t.set_y(0.012)
            rect = (rect[0], rect[1] + 0.03 * text.count("\n"), rect[2], rect[3])
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림3: 기하와 q 사상
    Gm = R["geometry"]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2), gridspec_kw=dict(width_ratios=[1.3, 1]))
    a = ax[0]
    g = np.deg2rad(22)
    a.plot([-3, 3], [0, 0], color="k", lw=3)
    a.fill_between([-3, 3], -0.5, 0, color="#dddddd")
    a.text(2.9, -0.32, L["sample"], ha="right", **Fk)
    xi = -2.6
    a.annotate("", xy=(0, 0), xytext=(xi, -xi * np.tan(g)), arrowprops=dict(arrowstyle="->", color="#c8553d", lw=2))
    a.text(xi, -xi * np.tan(g) + 0.12, L["inc"], color="#c8553d", **Fk)
    xe = 2.6
    a.annotate("", xy=(xe, xe * np.tan(g)), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="#2f6690", lw=2))
    a.text(xe - 0.9, xe * np.tan(g) + 0.15, L["spec"], color="#2f6690", **Fk)
    d = np.array([np.cos(g), np.sin(g)]) * 3.0
    e = np.array([-np.sin(g), np.cos(g)]) * 0.9
    a.plot([d[0] - e[0], d[0] + e[0]], [d[1] - e[1], d[1] + e[1]], color="#333333", lw=4)
    a.text(d[0] + e[0] + 0.05, d[1] + e[1], L["det"], **Fk)
    th = np.linspace(0, g, 30)
    a.plot(0.9 * np.cos(th), 0.9 * np.sin(th), color="k", lw=1)
    a.text(1.0, 0.12, "θ", fontsize=12)
    a.plot(-0.9 * np.cos(th), 0.9 * np.sin(th), color="k", lw=1)
    a.text(-1.15, 0.12, "θ", fontsize=12)
    a.annotate("", xy=(0.75, -0.18), xytext=(-0.75, -0.18), arrowprops=dict(arrowstyle="<->", color="#c8553d"))
    a.text(0, -0.45, L["foot"], ha="center", fontsize=9, color="#c8553d", **Fk)
    a.set_xlim(-3.1, 4.4); a.set_ylim(-0.6, 2.4); a.set_aspect("equal"); a.axis("off")
    a = ax[1]
    qx, qy = Gm["qe"][0] / Gm["k"], Gm["qe"][1] / Gm["k"]
    for i in range(0, C.N, 4):
        a.plot(qx[i, :], qy[i, :], color="#2f6690", lw=0.7)
        a.plot(qx[:, i], qy[:, i], color="#2f6690", lw=0.7)
    a.plot([], [], color="#2f6690", label=L["grid_det"])
    ux, uy = Gm["qu"][0] / Gm["k"], Gm["qu"][1] / Gm["k"]
    a.add_patch(plt.Rectangle((ux.min(), uy.min()), ux.max() - ux.min(), uy.max() - uy.min(), fill=True,
                              color="#999999", alpha=0.25, label=L["grid_uni"]))
    a.set_xlabel(L["qx"], **Fk); a.set_ylabel(L["qy"], **Fk)
    a.set_ylim(-0.23, 0.29)
    a.legend(prop=LG, fontsize=9, loc="upper center", ncol=1); a.grid(alpha=0.3)
    cap(fig, "fig3"); save(fig, "fig3-geometry.png")

    # 그림2: 반사율과 반사 위상
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for (name, r), col, lab in zip(Gm["curves"].items(), CURVE_COL, L["curves"]):
        rc = np.conj(r)
        ax[0].plot(Gm["gs"], 100 * np.abs(rc) ** 2, color=col, lw=2, label=lab)
        ax[1].plot(Gm["gs"], np.degrees(np.angle(rc)) % 360, color=col, lw=2, label=lab)
    for a in ax:
        a.axvspan(21, 25, color="#cccccc", alpha=0.5)
        a.set_xlabel(L["g"], **Fk); a.grid(alpha=0.3); a.set_xlim(0, 60)
    ax[0].text(23, 92, L["band"], ha="center", fontsize=8, **Fk)
    ax[0].set_ylabel(L["R"], **Fk); ax[1].set_ylabel(L["phi"], **Fk)
    ax[0].legend(prop=LG, fontsize=9)
    cap(fig, "fig2"); save(fig, "fig2-reflectivity.png")

    # 그림4: 기울어진 평면 보정
    T = R["tilt"]
    t20 = {r["mode"]: r for r in T if r["na"] == 0.20}
    fig = plt.figure(figsize=(13.5, 7.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[0.8, 1])
    asp = t20["ideal"]["dx"] / t20["ideal"]["dy"]
    m = t20["ideal"]["m"]
    ys, xs = np.nonzero(m)
    sl = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
    ref = t20["ideal"]["true"]
    sub = np.median(ref[m])
    for j, (img, title) in enumerate([(ref, L["true"]), (t20["corrected"]["obj"], L["corr"]),
                                      (t20["uncorrected"]["obj"], L["uncorr"])]):
        a = fig.add_subplot(gs[0, j])
        a.imshow(np.angle(np.conj(img / sub))[sl].T, cmap="twilight", vmin=-np.pi, vmax=np.pi, aspect=1 / asp)
        e = "" if j == 0 else f"  ({L['err']} {t20['corrected' if j == 1 else 'uncorrected']['err']:.2f})"
        a.set_title(title + e, fontsize=10, **Fk); a.set_xticks([]); a.set_yticks([])
    a = fig.add_subplot(gs[1, 1:])
    for mode, lab in [("ideal", L["ideal"]), ("corrected", L["corr"]), ("uncorrected", L["uncorr"])]:
        rr = sorted([r for r in T if r["mode"] == mode], key=lambda r: r["na"])
        a.plot([r["na"] for r in rr], [r["err"] for r in rr], "o-", color=COL[mode], label=lab)
    a.set_yscale("log"); a.set_xlabel(L["na"], **Fk); a.set_ylabel(L["err"], **Fk); a.grid(alpha=0.3, which="both")
    a.legend(prop=LG, fontsize=9)
    a = fig.add_subplot(gs[1, 0])
    Bi = R["bias"]
    for gg, col in zip((22, 30, 40), ("#2f6690", "#6a994e", "#c8553d")):
        rr = sorted([r for r in Bi if r["g"] == gg], key=lambda r: r["os"])
        a.plot([r["os"] for r in rr], [max(abs(r["dphB"]), abs(r["dphC"])) for r in rr], "o-", color=col,
               label=f"{gg}°")
    a.set_xscale("log", base=2); a.set_yscale("log"); a.set_xticks([1, 2, 4]); a.set_xticklabels(["1", "2", "4"])
    a.set_xlabel(L["os"], **Fk); a.set_ylabel(L["bias"], **Fk); a.grid(alpha=0.3, which="both")
    a.legend(prop=LG, fontsize=9)
    cap(fig, "fig4", rect=(0, 0.05, 1, 1)); save(fig, "fig4-tilted-plane.png")

    # 그림5: 영역 사이 비
    P = R["pipeline"]
    fig = plt.figure(figsize=(12.5, 7.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.75, 1])
    for j, gg in enumerate((21, 25)):
        r = next(x for x in P if x["g"] == gg and x["seed"] == 0)
        mm = r["A"] | r["B"] | r["Cp"]
        ys, xs = np.nonzero(mm)
        sl = (slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1))
        a = fig.add_subplot(gs[0, j])
        im = a.imshow(np.angle(np.conj(r["img"]))[sl].T, cmap="twilight", vmin=-np.pi, vmax=np.pi,
                      aspect=r["dy"] / r["dx"])
        a.set_title(L["img21" if gg == 21 else "img25"], fontsize=10, **Fk); a.set_xticks([]); a.set_yticks([])
        if j == 1:
            plt.colorbar(im, ax=a, fraction=0.04).set_label(L["phase_img"], **Fk)
    gfine = np.linspace(14, 46, 129)
    ptrue = (S.TRUE["dN"], S.TRUE["dO"], S.TRUE["rho"])
    mod = np.conj(S.model_ratio(ptrue, gfine))
    gg = np.array([x["g"] for x in P])
    zz = np.conj(np.array([x["rB"] for x in P]))
    mod_ph = np.unwrap(np.angle(mod))
    # 측정 위상은 같은 각도의 모형 위상 근처로 2π 를 맞춰 그린다
    meas_ph = np.interp(gg, gfine, mod_ph) + np.angle(zz / np.conj(S.model_ratio(ptrue, gg)))
    for j, (ym, yv, lab) in enumerate([(np.abs(mod), np.abs(zz), L["ratio_amp"]), (mod_ph, meas_ph, L["ratio_ph"])]):
        a = fig.add_subplot(gs[1, j])
        a.plot(gfine, ym, color="#333333", lw=1.5, label=L["model"])
        a.plot(gg, yv, "o", color="#c8553d", ms=4, label=L["meas"])
        a.axvspan(21, 25, color="#cccccc", alpha=0.5)
        a.set_xlabel(L["deg"], **Fk); a.set_ylabel(lab, **Fk); a.grid(alpha=0.3)
        a.legend(prop=LG, fontsize=9)
    cap(fig, "fig5", rect=(0, 0.05, 1, 1)); save(fig, "fig5-phase-steps.png")

    # 그림6: 신뢰 타원
    Fi = R["fit"]
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))
    for a, setname in zip(ax, ("narrow", "wide")):
        f = Fi[setname]
        ctr = (S.TRUE["dN"], S.TRUE["rho"])
        for use in ("amp", "phase", "both"):
            cov = f["crlb"][use]["cov"][np.ix_([0, 2], [0, 2])]
            ellipse(a, cov, ctr, COL[use], L[use])
            pts = np.array([x["x"] for x in f["fits"][use]])
            a.plot(pts[:, 0], pts[:, 2], "o" if use == "both" else "x", color=COL[use], ms=6)
        a.plot(*ctr, "k+", ms=12, mew=2, label=L["truth"])
        sd = f["crlb"]["amp"]["sd"]                 # 반사율 비만 쓴 타원이 다 보이도록. 위상 차만의 타원은 넘칠 수 있다
        a.set_xlim(ctr[0] - 2.6 * sd[0], ctr[0] + 2.6 * sd[0])
        a.set_ylim(ctr[1] - 2.6 * sd[2], ctr[1] + 2.6 * sd[2])
        a.set_title(L[setname], fontsize=10, **Fk); a.set_xlabel(L["dN"], **Fk); a.set_ylabel(L["rho"], **Fk)
        a.grid(alpha=0.3); a.legend(prop=LG, fontsize=8, loc="upper right")
    cap(fig, "fig6"); save(fig, "fig6-fit.png")

    # 그림7: 산란 지수
    Tr = R["transmission"]
    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.6))
    a = ax[0]
    show = np.where(Tr["mat"] >= 0, Tr["thick"], np.nan)
    im = a.imshow(np.ma.masked_invalid(show), cmap="Greys")
    for i in range(3):
        ys, xs = np.nonzero(Tr["mat"] == i)
        a.text(xs.mean(), 8, L["mat"][i], ha="center", color=COL[i], fontsize=10, fontweight="bold")
    a.text(np.nonzero(Tr["mat"] == -1)[1].mean(), 8, L["vac"], ha="center", fontsize=10, **Fk)
    a.contour(Tr["m"], [0.5], colors="#c8553d", linewidths=1)
    plt.colorbar(im, ax=a, fraction=0.04).set_label(L["thick"], **Fk)
    a.set_title(L["layout"], fontsize=10, **Fk); a.set_xticks([]); a.set_yticks([])
    for a, key in zip(ax[1:], ("cal", "raw")):
        d = Tr[key]
        for i in range(3):
            s = d["mat"] == i
            a.plot(-d["logamp"][s], -d["phase"][s], ".", ms=1.5, color=COL[i], alpha=0.4)
            a.plot([], [], "o", color=COL[i], label=L["mat"][i])
        if key == "cal":
            xx = np.linspace(0, max(-Tr["cal"]["logamp"].min(), 1.0), 10)
            for i, (mname, _, _) in enumerate(Tr["bands"]):
                a.plot(xx, Tr["fq_true"][mname] * xx, "--", color=COL[i], lw=1)
            a.plot([], [], "k--", lw=1, label=L["slope"])
        a.axhline(0, color="k", lw=0.5); a.axvline(0, color="k", lw=0.5)
        a.set_xlabel(L["la"], **Fk); a.set_ylabel(L["ph"], **Fk); a.set_title(L[key], fontsize=10, **Fk)
        a.grid(alpha=0.3); a.legend(prop=LG, fontsize=8, markerscale=1)
    cap(fig, "fig7"); save(fig, "fig7-scattering-quotient.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")

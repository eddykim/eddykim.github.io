"""계산 이미징과 타이코그래피 11편 그림 6개 생성 (한국어판·영문판). 수치는 studies11.py 의 캐시를 읽는다.
캡션은 본문에 둔다. 그림 안에는 넣지 않는다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

import fpm_core as C
import studies11 as S

SLUG = "ptycho-experimental-fpm"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
BLUE, RED, GREY, GREEN = "#2f6690", "#c8553d", "#555555", "#6a994e"

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "raw_bf": "가운데 LED 영상 (명시야)", "raw_df": "바깥 LED 영상 (암시야, 밝기 따로 맞춤)", "ledmap": "LED 225개의 조명 NA",
        "na_x": "조명 NA (x)", "na_y": "조명 NA (y)", "pupil": "대물렌즈 NA 0.1", "bf": "명시야", "df": "암시야",
        "na_ill": "조명 NA", "mean": "영상 평균 밝기 (가운데 = 1)", "usaf": "USAF 표적", "blood": "혈액 도말",
        "spec": "가운데 영상의 세기 스펙트럼 (log)", "circ2": "반지름 2r (2NA)", "radius": "반지름 (푸리에 화소)",
        "power": "평균 세기 스펙트럼 (log)", "nf": "노이즈 바닥 (진폭 상대)",
        "steps": ["기본", "LED 밝기", "+ 동공", "+ LED 위치"], "resid": "데이터 잔차", "floor": "노이즈 바닥 (명시야)",
        "iters": "반복 횟수", "median": "밝기 배율의 중앙값", "unbound": "묶지 않음", "bound": "중앙값 1로 묶음",
        "diff": "모형 - 측정 (가운데 영상)", "na_obj": "가정한 대물렌즈 NA",
        "z": "초점 어긋남 z ($\\mu$m)", "zfit": "보정 동공에 가장 가까운 z", "given": "제공된 보정 동공 (위상)",
        "epry": "평평한 동공에서 EPRY", "zepry": "z 를 맞춘 뒤 EPRY",
        "rawz": "원본 (가운데 LED)", "fpm": "FPM 복원", "freq": "공간주파수 (주기/$\\mu$m)", "frc": "FRC (LED 반씩 두 복원)",
        "half": "문턱 0.5", "axial": "축 방향 LED 덮개의 끝", "corner": "대각선 LED 덮개의 끝", "nodata": "측정 없음",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "raw_bf": "Centre LED image (bright field)", "raw_df": "Outer LED image (dark field, own scale)", "ledmap": "Illumination NA of 225 LEDs",
        "na_x": "Illumination NA (x)", "na_y": "Illumination NA (y)", "pupil": "Objective NA 0.1", "bf": "Bright field", "df": "Dark field",
        "na_ill": "Illumination NA", "mean": "Mean image brightness (centre = 1)", "usaf": "USAF target", "blood": "Blood smear",
        "spec": "Intensity spectrum of centre image (log)", "circ2": "Radius 2r (2NA)", "radius": "Radius (Fourier pixels)",
        "power": "Mean intensity spectrum (log)", "nf": "Noise floor (relative amplitude)",
        "steps": ["Baseline", "LED brightness", "+ pupil", "+ LED positions"], "resid": "Data residual", "floor": "Noise floor (bright field)",
        "iters": "Iterations", "median": "Median of brightness factors", "unbound": "Unbound", "bound": "Median tied to 1",
        "diff": "Model - measured (centre image)", "na_obj": "Assumed objective NA",
        "z": "Defocus z (µm)", "zfit": "z closest to calibrated pupil", "given": "Supplied calibrated pupil (phase)",
        "epry": "EPRY from flat pupil", "zepry": "EPRY after fitting z",
        "rawz": "Raw (centre LED)", "fpm": "FPM reconstruction", "freq": "Spatial frequency (cycles/µm)", "frc": "FRC (two half-LED reconstructions)",
        "half": "Threshold 0.5", "axial": "Edge of axial LED coverage", "corner": "Edge of diagonal LED coverage", "nodata": "No data",
    },
}

R = S.load()


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def save(fig, name, dpi=150):
        fig.tight_layout()
        fig.savefig(os.path.join(out, name), dpi=dpi, bbox_inches="tight", pad_inches=0.1); plt.close(fig)

    d = C.load_uconn("USAF_red.mat")
    bf = d["na_ill"] < d["NA"]

    # 그림1: 원시 영상과 LED 지도
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.7))
    j_df = int(np.argsort(d["na_ill"])[60])
    ax[0].imshow(np.sqrt(d["I"][0]), cmap="gray"); ax[0].set_title(L["raw_bf"], fontsize=10, **Fk)
    ax[1].imshow(np.sqrt(d["I"][j_df]), cmap="gray", vmax=np.percentile(np.sqrt(d["I"][j_df]), 99.5))
    ax[1].set_title(L["raw_df"], fontsize=10, **Fk)
    for a in ax[:2]:
        a.set_axis_off()
    na_k = d["k"] * d["wl"] / d["fov"]
    ax[2].scatter(na_k[bf, 1], na_k[bf, 0], s=14, color=BLUE, label=L["bf"])
    ax[2].scatter(na_k[~bf, 1], na_k[~bf, 0], s=14, color=RED, label=L["df"])
    ax[2].add_patch(Circle((0, 0), d["NA"], fill=False, ls="--", color=GREY, label=L["pupil"]))
    ax[2].set_aspect("equal"); ax[2].set_xlabel(L["na_x"], **Fk); ax[2].set_ylabel(L["na_y"], **Fk)
    ax[2].set_title(L["ledmap"], fontsize=10, **Fk); ax[2].legend(prop=LG, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False)
    save(fig, "fig1-data-geometry.png")

    # 그림2: 명·암시야 경계
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    for f, key, col in (("USAF_red.mat", "usaf", BLUE), ("bloodsmear_green.mat", "blood", RED)):
        b = R["boundary"][f]
        ax.semilogy(b["na"], b["mean"] / b["mean"][0], "o", ms=3.5, color=col, label=L[key])
    ax.axvline(0.1, color=GREY, ls="--", lw=1, label=L["pupil"])
    ax.set_xlabel(L["na_ill"], **Fk); ax.set_ylabel(L["mean"], **Fk); ax.grid(alpha=0.3, which="both")
    ax.legend(prop=LG, fontsize=9)
    save(fig, "fig2-brightfield-boundary.png")

    # 그림3: 노이즈 바닥
    N = R["noise"]
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    ax[0].imshow(N["spec"], cmap="magma"); c = d["nl"] // 2
    ax[0].add_patch(Circle((c, c), 2 * N["r"], fill=False, color="w", ls="--", label=L["circ2"]))
    ax[0].set_title(L["spec"], fontsize=10, **Fk); ax[0].set_axis_off(); ax[0].legend(prop=LG, fontsize=8, loc="lower right")
    ax[1].semilogy(N["prof"][: d["nl"] // 2 * 141 // 100], color=BLUE)
    ax[1].axvline(2 * N["r"], color=GREY, ls="--", label=L["circ2"])
    ax[1].set_xlabel(L["radius"], **Fk); ax[1].set_ylabel(L["power"], **Fk); ax[1].grid(alpha=0.3, which="both")
    ax[1].legend(prop=LG, fontsize=9)
    nb = N["na"] < N["NA"]
    ax[2].semilogy(N["na"][nb], N["nf"][nb], "o", ms=3.5, color=BLUE, label=L["bf"])
    ax[2].semilogy(N["na"][~nb], N["nf"][~nb], "o", ms=3.5, color=RED, label=L["df"])
    ax[2].axvline(0.1, color=GREY, ls="--", lw=1)
    ax[2].set_xlabel(L["na_ill"], **Fk); ax[2].set_ylabel(L["nf"], **Fk); ax[2].grid(alpha=0.3, which="both")
    ax[2].legend(prop=LG, fontsize=9)
    save(fig, "fig3-noise-floor.png")

    # 그림4: 보정 사다리
    Ld = R["ladder"]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6), gridspec_kw=dict(width_ratios=[1.15, 1, 1]))
    vals = [Ld["기본"], Ld["LED 밝기"], Ld["밝기 + 동공"], Ld["+ LED 위치"]]
    ax[0].bar(range(4), vals, color=[GREY, BLUE, BLUE, BLUE])
    for i, v in enumerate(vals):
        ax[0].text(i, v + 0.002, f"{v:.3f}", ha="center", fontsize=9)
    ax[0].axhline(Ld["노이즈 바닥 명시야"], color=RED, ls="--", label=L["floor"])
    ax[0].set_xticks(range(4)); ax[0].set_xticklabels(L["steps"], fontsize=9, **Fk); ax[0].set_ylabel(L["resid"], **Fk)
    ax[0].legend(prop=LG, fontsize=9)
    tr = Ld["묶지 않은 배율 중앙값"]
    ax[1].semilogy(range(1, len(tr) + 1), tr, "o-", ms=3.5, color=RED, label=L["unbound"])
    ax[1].semilogy([1, len(tr)], [1, 1], "-", color=BLUE, label=L["bound"])
    ax[1].set_xlabel(L["iters"], **Fk); ax[1].set_ylabel(L["median"], **Fk); ax[1].grid(alpha=0.3, which="both")
    ax[1].legend(prop=LG, fontsize=9)
    a, m, _ = Ld["차이"][min(Ld["차이"], key=lambda j: Ld["차이"][j][2])]
    lim = 0.4 * a.max()
    ax[2].imshow(m - a, cmap="RdBu_r", vmin=-lim, vmax=lim); ax[2].set_title(L["diff"], fontsize=10, **Fk); ax[2].set_axis_off()
    save(fig, "fig4-calibration-ladder.png")

    # 그림5: 수차 (혈액 데이터)
    B = R["blood"]
    supp = B["지지"]
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.4), gridspec_kw=dict(width_ratios=[1.4, 1, 1, 1]))
    zs = sorted(B["z 훑기"])
    ax[0].plot(np.array(zs) * 1e6, [B["z 훑기"][z] for z in zs], "o-", ms=3.5, color=BLUE)
    ax[0].axvline(B["보정 동공에 맞춘 z"][0] * 1e6, color=RED, ls="--", label=L["zfit"])
    ax[0].set_xlabel(L["z"], **Fk); ax[0].set_ylabel(L["resid"], **Fk); ax[0].grid(alpha=0.3); ax[0].legend(prop=LG, fontsize=9)
    ref = B["보정 동공 배열"]
    for a_, P, key in ((ax[1], ref, "given"), (ax[2], B["이상 동공 + EPRY"]["P"], "epry"), (ax[3], B["z 맞춘 뒤 EPRY"]["P"], "zepry")):
        z = P * np.exp(-1j * np.angle(np.sum(P[supp] * np.conj(ref[supp]))))
        a_.imshow(np.where(supp, np.angle(z), np.nan), cmap="twilight", vmin=-np.pi, vmax=np.pi)
        ys, xs = np.where(supp); a_.set_xlim(xs.min() - 2, xs.max() + 2); a_.set_ylim(ys.max() + 2, ys.min() - 2)
        a_.set_title(L[key], fontsize=10, **Fk); a_.set_axis_off()
    save(fig, "fig5-aberration-defocus.png")

    # 그림6: 분해능
    Q = R["resol"]
    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.9), gridspec_kw=dict(width_ratios=[1, 1, 1.2]))
    up = S.UP
    y0, y1, x0, x1 = 50, 90, 65, 105                      # 원본 화소 (8·9군 영역)
    ax[0].imshow(Q["raw"][y0:y1, x0:x1], cmap="gray", interpolation="nearest"); ax[0].set_title(L["rawz"], fontsize=10, **Fk)
    ax[1].imshow(np.abs(Q["o"])[y0 * up:y1 * up, x0 * up:x1 * up], cmap="gray"); ax[1].set_title(L["fpm"], fontsize=10, **Fk)
    for a_ in ax[:2]:
        a_.set_axis_off()
    fr = Q["f"] / Q["dxh"] * 1e-6
    ax[2].plot(fr[1:], Q["v"][1:], color=BLUE, label=L["frc"])
    ax[2].axhline(0.5, color=GREY, ls="--", lw=1, label=L["half"])
    na_ax = float(np.abs(d["k"] * d["wl"] / d["fov"]).max())
    f_ax, f_co = (Q["NA"] + na_ax) / (Q["wl"] * 1e6), (Q["NA"] + Q["na_max"]) / (Q["wl"] * 1e6)
    ax[2].axvline(f_ax, color=RED, ls=":", label=L["axial"])
    ax[2].axvline(f_co, color=GREEN, ls=":", label=L["corner"])
    ax[2].axvspan(f_co, fr.max(), color=GREY, alpha=0.15, label=L["nodata"])
    ax[2].set_xlabel(L["freq"], **Fk); ax[2].set_ylabel("FRC"); ax[2].set_ylim(0, 1.02); ax[2].grid(alpha=0.3)
    ax[2].legend(prop=LG, fontsize=9)
    save(fig, "fig6-resolution.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")

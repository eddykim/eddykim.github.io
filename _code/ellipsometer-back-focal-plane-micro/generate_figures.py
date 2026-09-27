"""타원계측기 6편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-back-focal-plane-micro/      (한국어)
      ../../assets/img/posts/ellipsometer-back-focal-plane-micro/en/   (영문)

계산은 bfp.py 의 뮬러 곱에서 나오고, 그 구현은 verify_bfp.py 가 Ye 2007 의 수치,
아베 사인조건, Munro & Torok 의 존스/뮬러 평균, 김영준 식 4.2 와 대조해 검증한 것이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyArrow

from bfp import (D, NA_DEFAULT, N_FUSED_SILICA, N_SI, N_SIO2, annulus,
                 averaged_mueller, bfp_map, brewster_angle, calibrate_radius,
                 depolarization_index, fourier_alpha, invert_annulus, max_angle,
                 psi_delta_film, radius_to_angle, tan2_psi_profile)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-back-focal-plane-micro",
)

D_FILM = 32.9           # nm, Ye 2007 의 시편
NA_LIST = (0.55, 0.75, 0.90, 0.95)
C_P, C_S, C_NA, C_BR = "tab:red", "tab:blue", "#2c3e50", "tab:green"
C_J, C_M, C_BAD = "tab:green", "#c0392b", "#c0392b"

_CACHE: dict = {}


def cached(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


# ---------------------------------------------------------------------------
# 그림 1 — 반경이 입사각이다
# ---------------------------------------------------------------------------


def figure1(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    r = np.linspace(0.0, 1.0, 400)
    for na, ls in zip(NA_LIST, ("-", "--", "-.", ":")):
        ax1.plot(r, np.rad2deg(radius_to_angle(r, na)), ls, lw=1.8,
                 label=f"NA = {na:.2f}  ($\\theta_{{max}}$ = "
                       f"{np.rad2deg(max_angle(na)):.1f}°)")
    for n, name, c in ((N_FUSED_SILICA, L["fused"], C_BR), (N_SI, L["si"], C_BAD)):
        tb = np.rad2deg(brewster_angle(n))
        ax1.axhline(tb, color=c, lw=1.2, ls="--")
        ax1.text(0.985, tb + 1.4, L["brew"].format(name, tb, np.sin(D(tb))),
                 fontsize=8.6, color=c, ha="right", **L["font"])
    ax1.set_xlabel(L["r_norm"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["aoi"], fontsize=10.5, **L["font"])
    ax1.set_title(L["f1a"], fontsize=11.5, pad=8, **L["font"])
    ax1.set_xlim(0, 1)
    ax1.set_ylim(0, 85)
    ax1.legend(fontsize=8.4, prop=L["legend"], loc="lower right", framealpha=0.92)
    ax1.grid(alpha=0.25)

    # 후초점면 원판에 등입사각 고리를 그린다
    ax2.add_patch(Circle((0, 0), 1.0, fc="#f2f4f7", ec=C_NA, lw=1.6))
    for t_deg in (10, 20, 30, 40, 50, 60):
        t = D(t_deg)
        rr = np.sin(t) / NA_DEFAULT
        if rr <= 1.0:
            ax2.add_patch(Circle((0, 0), rr, fc="none", ec="gray", lw=0.9, ls=":"))
            if t_deg in (20, 40, 60):
                ax2.text(0.0, rr, f"{t_deg}°", fontsize=8.4, color="gray",
                         ha="center", va="bottom",
                         bbox=dict(boxstyle="square,pad=0.05", fc="#f2f4f7",
                                   ec="none"))
    tb = brewster_angle(N_FUSED_SILICA)
    ax2.add_patch(Circle((0, 0), np.sin(tb) / NA_DEFAULT, fc="none", ec=C_BR,
                         lw=2.0))
    ax2.text(0, -1.30, L["f1note"].format(NA_DEFAULT, np.rad2deg(max_angle(NA_DEFAULT))),
             fontsize=9.0, ha="center", **L["font"])
    ax2.set_xlim(-1.25, 1.25)
    ax2.set_ylim(-1.52, 1.18)
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.set_title(L["f1b"], fontsize=11.5, pad=8, **L["font"])
    fig.suptitle(L["fig1"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(L["dir"], "fig1-bfp-mapping.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 방위각이 편광축을 돌린다
# ---------------------------------------------------------------------------


def figure2(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    ax1.add_patch(Circle((0, 0), 1.0, fc="#f2f4f7", ec=C_NA, lw=1.6))
    for phi_deg in range(0, 360, 30):
        p = D(phi_deg)
        x, y = 0.72 * np.cos(p), 0.72 * np.sin(p)
        ax1.add_patch(FancyArrow(x - 0.13 * np.cos(p), y - 0.13 * np.sin(p),
                                 0.26 * np.cos(p), 0.26 * np.sin(p),
                                 width=0.012, head_width=0.06, head_length=0.07,
                                 color=C_P, length_includes_head=True))
        ax1.add_patch(FancyArrow(x + 0.13 * np.sin(p), y - 0.13 * np.cos(p),
                                 -0.26 * np.sin(p), 0.26 * np.cos(p),
                                 width=0.008, head_width=0.05, head_length=0.06,
                                 color=C_S, length_includes_head=True, alpha=0.8))
    ax1.arrow(-1.18, -1.28, 0.55, 0.0, width=0.02, head_width=0.09,
              head_length=0.09, color="k", length_includes_head=True)
    ax1.text(-0.55, -1.28, L["f2_pol"], fontsize=9.2, va="center", **L["font"])
    ax1.text(0, 1.14, L["f2_p"], fontsize=9.4, color=C_P, ha="center", **L["font"])
    ax1.text(0, -1.10, L["f2_s"], fontsize=9.4, color=C_S, ha="center", **L["font"])
    ax1.set_xlim(-1.3, 1.3)
    ax1.set_ylim(-1.45, 1.3)
    ax1.set_aspect("equal")
    ax1.axis("off")
    ax1.set_title(L["f2a"], fontsize=11.5, pad=8, **L["font"])

    phi = np.linspace(0, 360, 721)
    ax2.plot(phi, np.cos(D(phi)) ** 2, "-", color=C_P, lw=2.0, label=L["f2_fp"])
    ax2.plot(phi, np.sin(D(phi)) ** 2, "--", color=C_S, lw=2.0, label=L["f2_fs"])
    for x, lab, c in ((0, L["f2_allp"], C_P), (90, L["f2_alls"], C_S)):
        ax2.axvline(x, color=c, lw=1.0, ls=":")
        ax2.text(x + 9, 0.52, lab, fontsize=8.8, color=c, rotation=90,
                 va="center", ha="left", **L["font"])
    ax2.set_xlabel(L["azim"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["f2_y"], fontsize=10.5, **L["font"])
    ax2.set_title(L["f2b"], fontsize=11.5, pad=8, **L["font"])
    ax2.set_xlim(0, 360)
    ax2.set_xticks(range(0, 361, 90))
    ax2.legend(fontsize=9.0, prop=L["legend"], loc="upper right", framealpha=0.92)
    ax2.grid(alpha=0.25)
    fig.suptitle(L["fig2"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(L["dir"], "fig2-azimuth-rotation.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — 고리 하나에서 읽는다
# ---------------------------------------------------------------------------

R_RING = 0.80


def figure3(L):
    fig, axes = plt.subplots(1, 3, figsize=(13.4, 4.4))
    X, Y, I = cached("map", lambda: bfp_map(D_FILM, NA_DEFAULT, n_pix=241,
                                            pol=0.0, ana=0.0))
    im = axes[0].pcolormesh(X, Y, I, cmap="magma", shading="auto")
    axes[0].add_patch(Circle((0, 0), R_RING, fc="none", ec="w", lw=1.8, ls="--"))
    axes[0].set_aspect("equal")
    axes[0].axis("off")
    axes[0].set_title(L["f3a"].format(D_FILM), fontsize=11.0, pad=8, **L["font"])
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.03)

    phi, Iring, (psi, dl) = cached("ring", lambda: annulus(
        D_FILM, R_RING, NA_DEFAULT, n_phi=720, pol=0.0, ana=0.0))
    a2, a4 = fourier_alpha(Iring, phi)
    dc = Iring.mean()
    pd = np.rad2deg(phi)
    axes[1].plot(pd, Iring / dc, "-", color=C_NA, lw=1.8, label=L["f3_meas"])
    axes[1].plot(pd, 1 + a2 * np.cos(2 * phi), "--", color=C_P, lw=1.3,
                 label=L["f3_2phi"].format(a2))
    axes[1].plot(pd, 1 + a4 * np.cos(4 * phi), ":", color=C_S, lw=1.6,
                 label=L["f3_4phi"].format(a4))
    axes[1].set_xlabel(L["azim"], fontsize=10.5, **L["font"])
    axes[1].set_ylabel(L["f3_y"], fontsize=10.5, **L["font"])
    axes[1].set_title(L["f3b"].format(np.rad2deg(radius_to_angle(R_RING))),
                      fontsize=11.0, pad=8, **L["font"])
    axes[1].set_xlim(0, 360)
    axes[1].set_xticks(range(0, 361, 90))
    axes[1].set_ylim(-0.12, 3.35)
    axes[1].legend(fontsize=8.2, prop=L["legend"], loc="upper center", ncol=3,
                   framealpha=0.92, columnspacing=0.8, handlelength=1.8)
    axes[1].grid(alpha=0.25)

    rs = np.linspace(0.12, 1.0, 45)
    got, true = [], []
    for r in rs:
        ph, Ir, (p0, d0) = annulus(D_FILM, r, NA_DEFAULT, n_phi=720,
                                   pol=0.0, ana=0.0)
        b2, b4 = fourier_alpha(Ir, ph)
        pr, dr = invert_annulus(b2, b4)
        got.append((np.rad2deg(pr), np.rad2deg(dr)))
        true.append((np.rad2deg(p0), np.rad2deg(abs(d0))))
    got, true = np.array(got), np.array(true)
    th = np.rad2deg(radius_to_angle(rs, NA_DEFAULT))
    axes[2].plot(th, true[:, 0], "-", color=C_P, lw=2.4, alpha=0.35)
    axes[2].plot(th, got[:, 0], "o", color=C_P, ms=3.4, label=r"$\Psi$")
    axes[2].plot(th, true[:, 1], "-", color=C_S, lw=2.4, alpha=0.35)
    axes[2].plot(th, got[:, 1], "s", color=C_S, ms=3.4, label=r"$\Delta$")
    axes[2].set_xlabel(L["aoi"], fontsize=10.5, **L["font"])
    axes[2].set_ylabel(L["f3_deg"], fontsize=10.5, **L["font"])
    axes[2].set_title(L["f3c"], fontsize=11.0, pad=8, **L["font"])
    axes[2].legend(fontsize=9.0, loc="center left", framealpha=0.92)
    axes[2].text(0.5, 0.03, L["f3note"], transform=axes[2].transAxes,
                 fontsize=8.8, ha="center", va="bottom", **L["font"],
                 bbox=dict(boxstyle="round,pad=0.3", fc="#fffbe6", ec="#d9c77a",
                           lw=0.8))
    axes[2].grid(alpha=0.25)
    fig.suptitle(L["fig3"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(L["dir"], "fig3-annular-readout.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — 집속은 평균이다
# ---------------------------------------------------------------------------

HALF = np.linspace(0.2, 12.0, 30)


def sweep_average():
    th_c = D(45.0)
    psi0, _ = psi_delta_film(np.array([th_c]), D_FILM)
    out = []
    for hw in HALF:
        Mj = averaged_mueller(th_c, D(hw), D_FILM, coherent=True)
        Mm = averaged_mueller(th_c, D(hw), D_FILM, coherent=False)
        pe = abs(np.rad2deg(0.5 * np.arccos(np.clip(-Mm[0, 1], -1, 1)) - psi0[0]))
        out.append((depolarization_index(Mj), depolarization_index(Mm), pe))
    return np.array(out)


def figure4(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    res = cached("avg", sweep_average)
    ax1.plot(HALF, res[:, 0], "-", color=C_J, lw=2.0, label=L["f4_j"])
    ax1.plot(HALF, res[:, 1], "--", color=C_M, lw=2.0, label=L["f4_m"])
    ax1.set_xlabel(L["f4_x"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["f4_y1"], fontsize=10.5, **L["font"])
    ax1.set_title(L["f4a"], fontsize=11.5, pad=8, **L["font"])
    ax1.legend(fontsize=9.0, prop=L["legend"], loc="lower left", framealpha=0.92)
    ax1.grid(alpha=0.25)
    ax1.text(0.97, 0.95, L["f4note"], transform=ax1.transAxes, fontsize=9.0,
             ha="right", va="top", **L["font"],
             bbox=dict(boxstyle="round,pad=0.3", fc="#fffbe6", ec="#d9c77a", lw=0.8))

    ax2.semilogy(HALF, np.maximum(res[:, 2], 1e-6), "-", color=C_M, lw=2.0)
    ax2.set_xlabel(L["f4_x"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["f4_y2"], fontsize=10.5, **L["font"])
    ax2.set_title(L["f4b"], fontsize=11.5, pad=8, **L["font"])
    i = int(np.argmin(np.abs(HALF - 5.0)))
    ax2.plot([HALF[i]], [res[i, 2]], "o", color=C_M, ms=7)
    ax2.annotate(L["f4_ann"].format(HALF[i], res[i, 2]),
                 xy=(HALF[i], res[i, 2]), xytext=(-12, 30),
                 textcoords="offset points", fontsize=9.0, ha="right", **L["font"],
                 arrowprops=dict(arrowstyle="->", color=C_M, lw=1.0))
    ax2.grid(alpha=0.25, which="both")
    fig.suptitle(L["fig4"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(L["dir"], "fig4-focusing-average.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 5 — 입사각 눈금을 브루스터각으로 맞춘다
# ---------------------------------------------------------------------------

NA_TRUE, NA_NOMINAL = 0.90, 0.95


def figure5(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.6))
    rr = np.linspace(0.02, 1.0, 600)
    prof = tan2_psi_profile(rr, NA_TRUE, N_FUSED_SILICA)
    r_b, th_b, _ = calibrate_radius(rr, NA_TRUE, N_FUSED_SILICA)
    ax1.semilogy(rr, prof, "-", color=C_NA, lw=1.9)
    ax1.set_ylim(prof.min() * 0.25, 3.0)
    ax1.axvline(r_b, color=C_BR, lw=1.4, ls="--")
    ax1.plot([r_b], [tan2_psi_profile(np.array([r_b]), NA_TRUE, N_FUSED_SILICA)[0]],
             "*", color=C_BR, ms=16, zorder=5)
    ax1.annotate(L["f5_min"].format(r_b, np.rad2deg(th_b)),
                 xy=(r_b, prof.min()), xytext=(0.12, prof.min() * 60),
                 fontsize=9.2, color=C_BR, **L["font"],
                 arrowprops=dict(arrowstyle="->", color=C_BR, lw=1.0))
    ax1.set_xlabel(L["r_norm"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(r"$\tan^2\Psi = R_p/R_s$", fontsize=10.5)
    ax1.set_title(L["f5a"], fontsize=11.5, pad=8, **L["font"])
    ax1.grid(alpha=0.25, which="both")

    th_true = np.rad2deg(radius_to_angle(rr, NA_TRUE))
    th_nom = np.rad2deg(radius_to_angle(rr, NA_NOMINAL))
    _, _, th_cal = calibrate_radius(rr, NA_TRUE, N_FUSED_SILICA)
    ax2.plot(rr, th_true, "-", color=C_NA, lw=2.4, alpha=0.35, label=L["f5_true"])
    ax2.plot(rr, th_nom, "--", color=C_BAD, lw=1.8,
             label=L["f5_nom"].format(NA_NOMINAL))
    ax2.plot(rr[::14], np.rad2deg(th_cal)[::14], "o", color=C_BR, ms=4.2,
             label=L["f5_cal"])
    ax2.set_xlabel(L["r_norm"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["aoi"], fontsize=10.5, **L["font"])
    ax2.set_title(L["f5b"], fontsize=11.5, pad=8, **L["font"])
    ax2.legend(fontsize=8.8, prop=L["legend"], loc="upper left", framealpha=0.92)
    err = np.max(np.abs(th_nom - th_true))
    ax2.text(0.97, 0.05, L["f5_err"].format(err), transform=ax2.transAxes,
             fontsize=9.0, ha="right", va="bottom", color=C_BAD, **L["font"])
    ax2.grid(alpha=0.25)
    fig.suptitle(L["fig5"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(os.path.join(L["dir"], "fig5-brewster-calibration.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 라벨
# ---------------------------------------------------------------------------

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "r_norm": "정규화 반경 $r$", "aoi": "입사각 (deg)",
        "azim": "방위각 $\\phi$ (deg)",
        "fused": "용융 실리카", "si": "실리콘",
        "brew": "{} 브루스터각 {:.1f}° (NA {:.2f} 필요)",
        # 그림 1
        "fig1": "후초점면의 반경이 입사각이다",
        "f1a": "(a) 아베 사인조건 $\\theta = \\sin^{-1}(r\\,\\mathrm{NA})$",
        "f1b": "(b) 후초점면의 등입사각 고리",
        "f1note": "NA = {:.2f}, 가장자리가 {:.1f}°\n초록 고리가 용융 실리카의 브루스터각",
        # 그림 2
        "fig2": "방위각이 편광축을 돌린다",
        "f2a": "(a) 후초점면의 국소 p·s 방향",
        "f2b": "(b) 고정된 편광자가 p와 s로 갈리는 비율",
        "f2_p": "p 방향 (방사형)", "f2_s": "s 방향 (접선형)",
        "f2_pol": "고정된 편광자",
        "f2_fp": "p 성분 $\\cos^2\\phi$", "f2_fs": "s 성분 $\\sin^2\\phi$",
        "f2_allp": "전부 p", "f2_alls": "전부 s",
        "f2_y": "세기 비율",
        # 그림 3
        "fig3": "고리 하나에서 읽는다 — 편광자와 분석기가 함께 돈다",
        "f3a": "(a) 합성한 후초점면 세기 (막 {:.1f} nm)",
        "f3b": "(b) 흰 고리의 방위각별 세기 (입사각 {:.1f}°)",
        "f3c": "(c) 반경을 바꿔 얻은 입사각별 값",
        "f3_meas": "합성 신호",
        "f3_2phi": "$2\\phi$ 성분 ($\\alpha_2$ = {:.3f})",
        "f3_4phi": "$4\\phi$ 성분 ($\\alpha_4$ = {:.3f})",
        "f3_y": "정규화 세기", "f3_deg": "각도 (deg)",
        "f3note": "옅은 선이 참값, 점이 고리에서 읽은 값",
        # 그림 4
        "fig4": "집속은 평균이다 — 무엇을 평균하는지가 갈린다",
        "f4a": "(a) 편광 소멸 지수",
        "f4b": "(b) 뮬러 평균이 남기는 $\\Psi$ 오차",
        "f4_x": "입사각 반폭 (deg)",
        "f4_y1": "편광 소멸 지수", "f4_y2": "$\\Psi$ 오차 (deg)",
        "f4_j": "존스 행렬 평균 (공초점 검출)",
        "f4_m": "뮬러 행렬 평균 (통상 검출)",
        "f4note": "존스 평균은 1을 지키고\n뮬러 평균만 편광을 소멸시킨다",
        "f4_ann": "반폭 {:.0f}°에서\n$\\Psi$ 오차 {:.4f}°",
        # 그림 5
        "fig5": "입사각 눈금을 브루스터각으로 맞춘다",
        "f5a": "(a) $\\tan^2\\Psi$의 최소점이 브루스터각이다",
        "f5b": "(b) 제조사 NA를 믿으면 생기는 오차",
        "f5_min": "최소점 $r_b$ = {:.3f}\n$\\to$ 이 자리가 {:.2f}°",
        "f5_true": "참 눈금 (NA 0.90)",
        "f5_nom": "제조사 NA {:.2f}를 그대로 쓸 때",
        "f5_cal": "브루스터각으로 보정한 눈금",
        "f5_err": "최대 {:.1f}° 어긋난다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "r_norm": "normalized radius $r$", "aoi": "angle of incidence (deg)",
        "azim": "azimuth $\\phi$ (deg)",
        "fused": "fused silica", "si": "silicon",
        "brew": "{} Brewster {:.1f}° (needs NA {:.2f})",
        "fig1": "Radius in the back focal plane is angle of incidence",
        "f1a": "(a) Abbe sine condition $\\theta = \\sin^{-1}(r\\,\\mathrm{NA})$",
        "f1b": "(b) iso-angle rings in the back focal plane",
        "f1note": "NA = {:.2f}, the rim is {:.1f}°\nthe green ring is fused silica's Brewster angle",
        "fig2": "Azimuth rotates the polarization axis",
        "f2a": "(a) local p and s directions across the pupil",
        "f2b": "(b) how a fixed polarizer splits into p and s",
        "f2_p": "p direction (radial)", "f2_s": "s direction (tangential)",
        "f2_pol": "fixed polarizer",
        "f2_fp": "p share $\\cos^2\\phi$", "f2_fs": "s share $\\sin^2\\phi$",
        "f2_allp": "all p", "f2_alls": "all s",
        "f2_y": "intensity share",
        "fig3": "Reading one ring: polarizer and analyzer turn together",
        "f3a": "(a) pupil intensity, {:.1f} nm film",
        "f3b": "(b) along the white ring ({:.1f}° incidence)",
        "f3c": "(c) recovered at each radius",
        "f3_meas": "signal",
        "f3_2phi": "$2\\phi$ ($\\alpha_2$={:.3f})",
        "f3_4phi": "$4\\phi$ ($\\alpha_4$={:.3f})",
        "f3_y": "normalized intensity", "f3_deg": "degrees",
        "f3note": "pale lines: truth; markers: read from the ring",
        "fig4": "Focusing is an average — what gets averaged matters",
        "f4a": "(a) depolarization index",
        "f4b": "(b) error in $\\Psi$ left by the Mueller average",
        "f4_x": "half-width in angle of incidence (deg)",
        "f4_y1": "depolarization index", "f4_y2": "error in $\\Psi$ (deg)",
        "f4_j": "Jones-matrix average (confocal detection)",
        "f4_m": "Mueller-matrix average (conventional detection)",
        "f4note": "the Jones average stays at 1;\nonly the Mueller average depolarizes",
        "f4_ann": "at {:.0f}° half-width,\n$\\Psi$ error {:.4f}°",
        "fig5": "Calibrating the angle scale with the Brewster angle",
        "f5a": "(a) the minimum of $\\tan^2\\Psi$ is the Brewster angle",
        "f5b": "(b) error from trusting the nominal NA",
        "f5_min": "minimum at $r_b$ = {:.3f}\n$\\to$ this radius is {:.2f}°",
        "f5_true": "true scale (NA 0.90)",
        "f5_nom": "using the nominal NA {:.2f}",
        "f5_cal": "calibrated with the Brewster angle",
        "f5_err": "off by up to {:.1f}°",
    },
}


def main():
    print(f"NA 별 최대 입사각: " + ", ".join(
        f"{na}->{np.rad2deg(max_angle(na)):.1f}°" for na in NA_LIST))
    for n, name in ((N_FUSED_SILICA, "용융실리카"), (1.5, "유리"), (N_SI, "실리콘")):
        tb = brewster_angle(n)
        print(f"  {name}: 브루스터 {np.rad2deg(tb):.2f}° -> NA {np.sin(tb):.3f} 필요")
    phi, Ir, (psi, dl) = cached("ring", lambda: annulus(
        D_FILM, R_RING, NA_DEFAULT, n_phi=720, pol=0.0, ana=0.0))
    a2, a4 = fourier_alpha(Ir, phi)
    pr, dr = invert_annulus(a2, a4)
    print(f"고리 r={R_RING} (입사각 {np.rad2deg(radius_to_angle(R_RING)):.1f}°): "
          f"a2={a2:.4f}, a4={a4:.4f} -> Psi {np.rad2deg(pr):.3f}° "
          f"(참 {np.rad2deg(psi):.3f}°), Delta {np.rad2deg(dr):.3f}° "
          f"(참 {np.rad2deg(abs(dl)):.3f}°)")
    res = cached("avg", sweep_average)
    for i in (0, len(HALF) // 2, len(HALF) - 1):
        print(f"  반폭 {HALF[i]:5.1f}°: 존스 지수 {res[i,0]:.9f}, "
              f"뮬러 지수 {res[i,1]:.6f}, Psi 오차 {res[i,2]:.5f}°")
    rr = np.linspace(0.02, 1.0, 600)
    err = np.max(np.abs(np.rad2deg(radius_to_angle(rr, NA_NOMINAL))
                        - np.rad2deg(radius_to_angle(rr, NA_TRUE))))
    print(f"제조사 NA {NA_NOMINAL} 를 믿으면 최대 {err:.2f}° 어긋난다 "
          f"(실제 NA {NA_TRUE})")
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        figure1(L); figure2(L); figure3(L); figure4(L); figure5(L)
        print(f"  [{lang}] 그림 5개 저장 -> {os.path.normpath(L['dir'])}")


if __name__ == "__main__":
    main()

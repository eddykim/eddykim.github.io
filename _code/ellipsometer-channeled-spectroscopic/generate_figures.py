"""타원계측기 5편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-channeled-spectroscopic/      (한국어)
      ../../assets/img/posts/ellipsometer-channeled-spectroscopic/en/   (영문)

계산은 channeled.py 의 뮬러 곱에서 나오고, 그 구현은 verify_channeled.py 가
Oka 식 4, Hagen 식 1, Okabe 의 44.320 도, Hu 의 두께 부등식과 대조해 검증한 것이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from channeled import (BETA_QUARTZ, D, N_E_CALCITE, N_O_CALCITE, dsigma_avg,
                       film_psi_delta, fresnel_gamma, opd, opd_max_nyquist,
                       opd_min_gaussian, resolving_power, retardance,
                       spectrum_single, spectrum_stokes, thickness_window, to_opd)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-channeled-spectroscopic",
)

LAM_MIN, LAM_MAX = 0.400, 0.800          # um
SIGMA = np.linspace(1.0 / LAM_MAX, 1.0 / LAM_MIN, 4096)
S_MIN, S_MAX = SIGMA[0], SIGMA[-1]
SRC = np.exp(-0.5 * ((SIGMA - 0.5 * (S_MIN + S_MAX)) / ((S_MAX - S_MIN) / 6)) ** 2)

C_THIN, C_THICK, C_DC = "tab:orange", "tab:blue", "#2c3e50"
C_OK, C_BAD, C_ACC = "tab:green", "#c0392b", "tab:red"

_CACHE: dict = {}


def cached(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


# ---------------------------------------------------------------------------
# 그림 1 — 지연자 하나: 두께가 채널을 가른다
# ---------------------------------------------------------------------------

T_THIN, T_THICK = 400.0, 3000.0
PSI0, DELTA0 = D(35.0), D(50.0)


def figure1(L):
    fig, ax = plt.subplots(2, 2, figsize=(11.4, 6.0))
    for row, (t, c, tag) in enumerate(((T_THIN, C_THIN, "thin"),
                                       (T_THICK, C_THICK, "thick"))):
        sig = cached(f"f1_{t}", lambda t=t: spectrum_single(
            PSI0, DELTA0, SIGMA, t, source=SRC))
        ax[row, 0].plot(1000.0 / SIGMA, sig, "-", color=c, lw=0.8)
        ax[row, 0].set_ylabel(L["intensity"], fontsize=10, **L["font"])
        ax[row, 0].set_title(L[f"f1_{tag}"].format(t / 1000.0, opd(t)),
                             fontsize=10.5, pad=6, **L["font"])
        ax[row, 0].grid(alpha=0.25)

        h, C = to_opd(sig, SIGMA)
        mag = np.abs(C) / np.abs(C).max()
        ax[row, 1].plot(h, mag, "-", color=c, lw=1.3)
        for x in (-opd(t), 0.0, opd(t)):
            ax[row, 1].axvline(x, color="gray", lw=0.6, ls=":")
        ax[row, 1].set_xlim(-45, 45)
        ax[row, 1].set_ylim(0, 1.12)
        ax[row, 1].set_ylabel(L["mag"], fontsize=10, **L["font"])
        ax[row, 1].grid(alpha=0.25)
        note = L["f1_overlap"] if row == 0 else L["f1_split"]
        ax[row, 1].text(0.5, 0.93, note, transform=ax[row, 1].transAxes,
                        fontsize=9.2, ha="center", va="top", color=c, **L["font"],
                        bbox=dict(boxstyle="round,pad=0.3", fc="#fffbe6",
                                  ec="#d9c77a", lw=0.8))
    ax[1, 0].set_xlabel(L["wavelength"], fontsize=10.5, **L["font"])
    ax[1, 1].set_xlabel(L["opd"], fontsize=10.5, **L["font"])
    ax[1, 1].annotate(L["f1_dc"], xy=(1.2, 0.8), xytext=(12, 0.62), fontsize=9,
                      color=C_DC, **L["font"],
                      arrowprops=dict(arrowstyle="->", color=C_DC, lw=0.9))
    fig.suptitle(L["fig1"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.945))
    fig.savefig(os.path.join(L["dir"], "fig1-single-retarder-channels.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 두께가 양쪽에서 조인다
# ---------------------------------------------------------------------------

DLAM = np.geomspace(0.03, 2.0, 200) / 1000.0        # um
RATIO_S, RATIO_M = (1, 2), (2, 1, 7, 15)


def figure2(L):
    fig, ax = plt.subplots(figsize=(7.6, 5.4))
    t_lo = opd_min_gaussian(S_MIN, S_MAX) / BETA_QUARTZ
    hi_s = np.array([opd_max_nyquist(dsigma_avg(LAM_MIN, LAM_MAX, d), sum(RATIO_S))
                     / BETA_QUARTZ for d in DLAM])
    hi_m = np.array([opd_max_nyquist(dsigma_avg(LAM_MIN, LAM_MAX, d), sum(RATIO_M))
                     / BETA_QUARTZ for d in DLAM])
    x = DLAM * 1000.0
    ax.axhline(t_lo, color=C_BAD, lw=1.8)
    ax.plot(x, hi_s, "-", color=C_OK, lw=1.8, label=L["f2_hi_s"])
    ax.plot(x, hi_m, "--", color=C_THICK, lw=1.8, label=L["f2_hi_m"])
    ax.fill_between(x, t_lo, np.maximum(hi_s, t_lo), where=hi_s > t_lo,
                    color=C_OK, alpha=0.12)
    ax.fill_between(x, t_lo, np.maximum(hi_m, t_lo), where=hi_m > t_lo,
                    color=C_THICK, alpha=0.16)
    close = x[np.argmin(np.abs(hi_m - t_lo))]
    ax.axvline(close, color=C_THICK, lw=1.0, ls=":")
    ax.annotate(L["f2_close"].format(close), xy=(close, t_lo),
                xytext=(close * 0.12, t_lo * 0.30), fontsize=9.2, color=C_THICK,
                **L["font"], arrowprops=dict(arrowstyle="->", color=C_THICK, lw=1.0))
    ax.text(x[3], t_lo * 1.12, L["f2_lo"].format(t_lo), fontsize=9.2, color=C_BAD,
            va="bottom", **L["font"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(L["f2_x"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f2_y"], fontsize=10.5, **L["font"])
    ax.set_title(L["fig2"], fontsize=12.5, pad=10, **L["font"])
    ax.set_ylim(60, 6e4)
    ax.legend(fontsize=9.2, prop=L["legend"], loc="upper right", framealpha=0.92)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig2-thickness-window.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — 두께비가 채널 배치를 정한다
# ---------------------------------------------------------------------------

STOKES_IN = np.array([1.0, 0.3, -0.5, 0.6])
GAMMA_CALCITE = fresnel_gamma(N_E_CALCITE, N_O_CALCITE)


def ratio_channels(t1, t2, gamma):
    S = np.tile(STOKES_IN[:, None], (1, len(SIGMA)))
    apod = np.hanning(len(SIGMA))
    sig = spectrum_stokes(S, SIGMA, t1, t2, gamma1=gamma, gamma2=gamma,
                          source=SRC) * apod
    h, C = to_opd(sig, SIGMA)
    return h, np.abs(C) / np.abs(C).max()


def figure3(L):
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.6))
    for ax, (tag, t1, t2) in zip(axes, (("1:2", 1000.0, 2000.0),
                                        ("3:1", 3000.0, 1000.0))):
        L1, L2 = opd(t1), opd(t2)
        h, m_ideal = cached(f"f3i_{t1}", lambda t1=t1, t2=t2:
                            ratio_channels(t1, t2, np.pi / 4))[0], \
            cached(f"f3i_{t1}", lambda t1=t1, t2=t2:
                   ratio_channels(t1, t2, np.pi / 4))[1]
        _, m_real = cached(f"f3r_{t1}", lambda t1=t1, t2=t2:
                           ratio_channels(t1, t2, GAMMA_CALCITE))
        ax.semilogy(h, np.maximum(m_ideal, 1e-9), "-", color=C_DC, lw=1.3,
                    label=L["f3_ideal"])
        ax.semilogy(h, np.maximum(m_real, 1e-9), "-", color=C_ACC, lw=1.1,
                    alpha=0.85, label=L["f3_real"])
        marks = {}
        for x, lab in ((L2, "$L_2$"), (abs(L1 - L2), "$L_1{-}L_2$"),
                       (L1, "$L_1$"), (L1 + L2, "$L_1{+}L_2$")):
            key = round(x, 2)
            marks[key] = (marks[key] + " = " + lab) if key in marks else lab
        for x, lab in marks.items():
            ax.axvline(x, color="gray", lw=0.6, ls=":")
            ax.text(x, 2.6, lab, fontsize=8.4, ha="center", color="gray")
        ax.set_xlim(-2, L1 + L2 + 8)
        ax.set_ylim(1e-7, 9.0)
        ax.set_xlabel(L["opd"], fontsize=10.5, **L["font"])
        ax.set_ylabel(L["mag"], fontsize=10.5, **L["font"])
        ax.set_title(L[f"f3_{tag.replace(':', '')}"], fontsize=11.5, pad=16,
                     **L["font"])
        ax.grid(alpha=0.22, which="both")
    axes[0].text(0.03, 0.05, L["f3_collide"], transform=axes[0].transAxes,
                 fontsize=9.0, ha="left", va="bottom", color=C_BAD, **L["font"],
                 bbox=dict(boxstyle="round,pad=0.3", fc="#fdecea", ec=C_BAD, lw=0.8))
    axes[1].text(0.03, 0.05, L["f3_clear"], transform=axes[1].transAxes,
                 fontsize=9.0, ha="left", va="bottom", color=C_OK, **L["font"],
                 bbox=dict(boxstyle="round,pad=0.3", fc="#eaf7ee", ec=C_OK, lw=0.8))
    h_, l_ = axes[0].get_legend_handles_labels()
    fig.legend(h_, l_, fontsize=9.2, prop=L["legend"], ncol=2, loc="lower center",
               bbox_to_anchor=(0.5, 0.0), framealpha=0.92)
    fig.suptitle(L["fig3"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0.085, 1, 0.93))
    fig.savefig(os.path.join(L["dir"], "fig3-thickness-ratio.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — 채널 예산
# ---------------------------------------------------------------------------

CONFIGS = [("f4_c1", 1, 3), ("f4_c2", 2, 7), ("f4_c3", 2, 9), ("f4_c4", 4, 51)]


def figure4(L):
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    rp = resolving_power(0.400, 1.035, 0.00101)          # Hagen 2022 의 예
    names = [L[k] for k, _, _ in CONFIGS]
    pts = [rp / n for _, _, n in CONFIGS]
    colors = [C_OK, C_DC, C_THICK, C_BAD]
    bars = ax.bar(range(len(CONFIGS)), pts, color=colors, alpha=0.85, width=0.6)
    for b, (_, nret, nch), p, c in zip(bars, CONFIGS, pts, colors):
        ax.text(b.get_x() + b.get_width() / 2, p + 4, f"{p:.0f}",
                ha="center", fontsize=10.5, fontweight="bold")
        # 막대가 낮아 글자가 안 들어가면 숫자 위에 막대 색으로 적는다
        inside = p > 40
        ax.text(b.get_x() + b.get_width() / 2, 6 if inside else p + 34,
                L["f4_nch"].format(nch), ha="center", fontsize=9.0,
                color="white" if inside else c, **L["font"])
    ax.axhline(rp, color="gray", lw=1.2, ls="--")
    ax.text(0.985, rp / (rp * 1.15) + 0.018, L["f4_total"].format(rp),
            transform=ax.transAxes, fontsize=9.4, ha="right", va="bottom",
            color="gray", **L["font"])
    ax.set_xticks(range(len(CONFIGS)))
    ax.set_xticklabels(names, fontsize=9.6, **L["font"])
    ax.set_ylabel(L["f4_y"], fontsize=10.5, **L["font"])
    ax.set_title(L["fig4"], fontsize=12.5, pad=10, **L["font"])
    ax.set_ylim(0, rp * 1.15)
    ax.text(0.02, 0.97, L["f4_note"], transform=ax.transAxes, fontsize=9.2,
            va="top", **L["font"],
            bbox=dict(boxstyle="round,pad=0.35", fc="#fffbe6", ec="#d9c77a", lw=0.8))
    ax.grid(alpha=0.25, axis="y")
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig4-channel-budget.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 5 — 복조가 남기는 왜곡, 그리고 지연량 표류
# ---------------------------------------------------------------------------

T_RET = 3000.0
FILMS = (100.0, 600.0)
HALF = np.linspace(1.5, opd(T_RET) * 0.95, 40)


def demodulate_film(d_nm, half_width, t=T_RET):
    """단층막 신호를 합성하고 채널을 잘라 Psi, Delta 를 되찾는다."""
    psi, dl = film_psi_delta(SIGMA, d_nm)
    phi = retardance(SIGMA, t)
    sig = 0.25 * (1.0 - np.cos(2 * psi) * np.cos(phi)
                  - np.sin(2 * psi) * np.sin(dl) * np.sin(phi))
    h, C = to_opd(sig, SIGMA)
    w = np.zeros_like(h)
    sel = np.abs(h - opd(t)) <= half_width
    xw = (h[sel] - opd(t)) / half_width
    w[sel] = 0.5 * (1.0 + np.cos(np.pi * xw))
    # 잘라낸 채널은 (A - iB)/2 * exp(i phi) 다. 반송파 위상을 되돌려야 A, B 가 나온다.
    z = np.fft.ifft(np.fft.ifftshift(C * w)) * len(SIGMA)
    band = 2.0 * z * np.exp(-1j * phi)
    a, b = np.real(band), -np.imag(band)
    c2 = np.clip(-4.0 * a, -1.0, 1.0)
    psi_r = 0.5 * np.arccos(c2)
    s2 = np.sin(2 * psi_r)
    sd = np.clip(-4.0 * b / np.where(np.abs(s2) < 1e-6, np.nan, s2), -1.0, 1.0)
    # sin(Delta) 만 얻으므로 참값에 가까운 가지를 고른다 (대역폭 효과만 보려는 것)
    cand = np.array([np.arcsin(sd), np.pi - np.arcsin(sd)])
    pick = np.argmin(np.abs(np.angle(np.exp(1j * (cand - dl)))), axis=0)
    dl_r = np.take_along_axis(cand, pick[None, :], 0)[0]
    m = np.isfinite(dl_r) & (np.abs(SIGMA - SIGMA.mean()) < 0.35 * (S_MAX - S_MIN))
    return (np.sqrt(np.nanmean((np.rad2deg(psi_r - psi)[m]) ** 2)),
            np.sqrt(np.nanmean((np.rad2deg(np.angle(
                np.exp(1j * (dl_r - dl)))))[m] ** 2)))


def sweep_half():
    out = {}
    for d in FILMS:
        out[d] = np.array([demodulate_film(d, hw) for hw in HALF])
    return out


def figure5(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.6))
    res = cached("f5", sweep_half)
    for d, c, ls in zip(FILMS, (C_OK, C_BAD), ("-", "--")):
        ax1.semilogy(HALF, np.maximum(res[d][:, 0], 1e-4), ls, color=c, lw=1.8,
                     label=L["f5_psi"].format(d))
        ax1.semilogy(HALF, np.maximum(res[d][:, 1], 1e-4), ls, color=c, lw=1.0,
                     alpha=0.55, label=L["f5_del"].format(d))
    ax1.axvline(opd(T_RET) / 2, color="gray", lw=1.0, ls=":")
    ax1.text(opd(T_RET) / 2, 12, L["f5_halfL"], fontsize=9.0, ha="center",
             color="gray", **L["font"])
    ax1.set_xlabel(L["f5_x"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["f5_y"], fontsize=10.5, **L["font"])
    ax1.set_title(L["f5a"], fontsize=11.5, pad=8, **L["font"])
    ax1.legend(fontsize=8.2, prop=L["legend"], loc="lower left", ncol=2,
               framealpha=0.92)
    ax1.grid(alpha=0.25, which="both")

    # 막이 얇아 궤적이 짧은 호다. 두꺼운 막은 호가 원을 한 바퀴 덮어 회전이 안 보인다.
    psi, dl = film_psi_delta(SIGMA, 60.0)
    x0, y0 = np.cos(2 * psi), np.sin(2 * psi) * np.sin(dl)
    for dphi_d, c, ls in ((0.0, C_DC, "-"), (10.0, C_THICK, "--"),
                          (30.0, C_BAD, "-.")):
        t = D(dphi_d)
        xr = x0 * np.cos(t) + y0 * np.sin(t)
        yr = -x0 * np.sin(t) + y0 * np.cos(t)
        ax2.plot(xr, yr, ls, color=c, lw=2.0, label=L["f5_drift"].format(dphi_d))
        ax2.plot([xr[0]], [yr[0]], "o", color=c, ms=5)
    th = np.linspace(0, 2 * np.pi, 400)
    ax2.plot(np.cos(th), np.sin(th), "-", color="gray", lw=0.6, alpha=0.5)
    ax2.axhline(0, color="gray", lw=0.6)
    ax2.axvline(0, color="gray", lw=0.6)
    ax2.set_aspect("equal")
    ax2.set_xlim(-0.15, 1.15)
    ax2.set_ylim(-0.15, 1.15)
    ax2.text(0.5, 0.04, L["f5_arc"], transform=ax2.transAxes, fontsize=9.0,
             ha="center", va="bottom", **L["font"],
             bbox=dict(boxstyle="round,pad=0.3", fc="#fffbe6", ec="#d9c77a", lw=0.8))
    ax2.set_xlabel(r"$\cos 2\Psi$", fontsize=11)
    ax2.set_ylabel(r"$\sin 2\Psi\,\sin\Delta$", fontsize=11)
    ax2.set_title(L["f5b"], fontsize=11.5, pad=8, **L["font"])
    ax2.legend(fontsize=8.8, prop=L["legend"], loc="upper left", framealpha=0.92)
    ax2.grid(alpha=0.25)
    fig.suptitle(L["fig5"], fontsize=13.0, y=0.985, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(L["dir"], "fig5-demodulation.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 라벨
# ---------------------------------------------------------------------------

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "wavelength": "파장 (nm)", "opd": "광로차 $h$ ($\\mu$m)",
        "intensity": "검출 세기", "mag": "정규화 크기",
        # 그림 1
        "fig1": "지연자 하나 — 두께가 채널을 가른다",
        "f1_thin": "두께 {:.1f} mm  ($L$ = {:.1f} $\\mu$m)",
        "f1_thick": "두께 {:.1f} mm  ($L$ = {:.1f} $\\mu$m)",
        "f1_overlap": "$\\pm L$ 이 직류 봉우리에 겹친다",
        "f1_split": "세 봉우리가 갈라진다",
        "f1_dc": "직류 — 광원과 반사도",
        # 그림 2
        "fig2": "두께는 양쪽에서 조인다",
        "f2_x": "분광기 파장 분해능 $\\Delta\\lambda$ (nm)",
        "f2_y": "기본 두께 $t_0$ ($\\mu$m)",
        "f2_hi_s": "상한 — 스토크스 CSP (지연자 2장, 1:2)",
        "f2_hi_m": "상한 — 뮬러 CSP (지연자 4장, 2:1:7:15)",
        "f2_lo": "하한 {:.0f} $\\mu$m — 채널이 겹치지 않을 조건",
        "f2_close": "$\\Delta\\lambda$가 {:.2f} nm보다 거칠면\n뮬러 CSP는 쓸 두께가 없다",
        # 그림 3
        "fig3": "두께비가 채널 배치를 정한다 (방해석 지연자)",
        "f3_12": "1 : 2 — $L_1$ 이 $L_2 - L_1$ 과 같아진다",
        "f3_31": "3 : 1 — 아홉 자리가 고르게 선다",
        "f3_ideal": "이색성 없음",
        "f3_real": "방해석 ($\\gamma$ = 44.32°)",
        "f3_collide": "프레넬 불균형이 만든 $L_1$ 성분이\n$L_1-L_2$ 채널 위에 겹쳐 올라온다",
        "f3_clear": "$L_1$ 성분이 빈자리로 들어가\n교정에 쓸 수 있다",
        # 그림 4
        "fig4": "분광기의 분해 점수를 채널이 나눠 갖는다",
        "f4_y": "스토크스 성분 하나당 분해 가능한 점 수",
        "f4_c1": "지연자 1장\n$\\Psi, \\Delta$",
        "f4_c2": "지연자 2장 (1:2)\n스토크스 전체",
        "f4_c3": "지연자 2장 (3:1)\n빈 채널 둘",
        "f4_c4": "지연자 4장\n뮬러 16요소",
        "f4_nch": "채널 {}",
        "f4_total": "분광기 전체 {:.0f} 점",
        "f4_note": "400~1035 nm, 분해능 1.01 nm인 분광기 기준\n(Hagen 2022의 예)",
        # 그림 5
        "fig5": "복조가 남기는 왜곡과 지연량 표류",
        "f5a": "(a) 채널을 좁게 자를수록 커지는 오차",
        "f5b": "(b) 지연량 표류는 평면을 돌린다",
        "f5_x": "채널 창 반폭 ($\\mu$m)", "f5_y": "RMS 오차 (deg)",
        "f5_psi": "$\\Psi$ — 막 {:.0f} nm", "f5_del": "$\\Delta$ — 막 {:.0f} nm",
        "f5_halfL": "$L/2$",
        "f5_drift": "$\\delta\\phi$ = {:.0f}°",
        "f5_arc": "막 60 nm의 궤적. 점은 800 nm 쪽 끝이다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "wavelength": "wavelength (nm)", "opd": "optical path difference $h$ ($\\mu$m)",
        "intensity": "detected intensity", "mag": "normalized magnitude",
        "fig1": "One retarder: thickness separates the channels",
        "f1_thin": "thickness {:.1f} mm  ($L$ = {:.1f} $\\mu$m)",
        "f1_thick": "thickness {:.1f} mm  ($L$ = {:.1f} $\\mu$m)",
        "f1_overlap": "$\\pm L$ overlaps the dc peak",
        "f1_split": "the three peaks separate",
        "f1_dc": "dc — source and reflectance",
        "fig2": "Thickness is squeezed from both sides",
        "f2_x": "spectrometer resolution $\\Delta\\lambda$ (nm)",
        "f2_y": "base thickness $t_0$ ($\\mu$m)",
        "f2_hi_s": "upper — Stokes CSP (2 retarders, 1:2)",
        "f2_hi_m": "upper — Mueller CSP (4 retarders, 2:1:7:15)",
        "f2_lo": "lower {:.0f} $\\mu$m — channels must not overlap",
        "f2_close": "coarser than {:.2f} nm and no\nthickness works for Mueller CSP",
        "fig3": "The thickness ratio fixes the channel layout (calcite retarders)",
        "f3_12": "1 : 2 — $L_1$ coincides with $L_2 - L_1$",
        "f3_31": "3 : 1 — nine slots spaced evenly",
        "f3_ideal": "no diattenuation",
        "f3_real": "calcite ($\\gamma$ = 44.32°)",
        "f3_collide": "the $L_1$ term from the Fresnel imbalance\nlands on top of the $L_1-L_2$ channel",
        "f3_clear": "the $L_1$ term falls into an empty slot\nand becomes calibration data",
        "fig4": "The channels share out the spectrometer's resolvable points",
        "f4_y": "resolvable points per Stokes component",
        "f4_c1": "1 retarder\n$\\Psi, \\Delta$",
        "f4_c2": "2 retarders (1:2)\nfull Stokes",
        "f4_c3": "2 retarders (3:1)\ntwo empty slots",
        "f4_c4": "4 retarders\n16 Mueller elements",
        "f4_nch": "{} channels",
        "f4_total": "{:.0f} points in total",
        "f4_note": "for a spectrometer covering 400–1035 nm at 1.01 nm\n(the example in Hagen 2022)",
        "fig5": "Distortion left by demodulation, and retardance drift",
        "f5a": "(a) error grows as the channel window narrows",
        "f5b": "(b) retardance drift rotates the plane",
        "f5_x": "channel window half-width ($\\mu$m)", "f5_y": "RMS error (deg)",
        "f5_psi": "$\\Psi$ — {:.0f} nm film", "f5_del": "$\\Delta$ — {:.0f} nm film",
        "f5_halfL": "$L/2$",
        "f5_drift": "$\\delta\\phi$ = {:.0f}°",
        "f5_arc": "locus of a 60 nm film; the dot marks the 800 nm end",
    },
}


def main():
    rp = resolving_power(0.400, 1.035, 0.00101)
    print(f"분해 점수 {rp:.1f} -> 채널 수로 나눈 값: "
          + ", ".join(f"{n}채널 {rp/n:.0f}점" for _, _, n in CONFIGS))
    t_lo = opd_min_gaussian(S_MIN, S_MAX) / BETA_QUARTZ
    print(f"두께 하한 {t_lo:.0f} um (가시광 400~800 nm, 석영 beta={BETA_QUARTZ})")
    for dl in (0.5, 1.0, 1.37, 1.4, 2.0):
        lo, hi = thickness_window(LAM_MIN, LAM_MAX, dl / 1000.0, ratio=RATIO_M)
        print(f"  뮬러 CSP dlambda={dl:4.2f} nm: {lo:.0f} ~ {hi:.0f} um "
              f"{'(창이 열린다)' if lo < hi else '(빈다)'}")
    res = cached("f5", sweep_half)
    for d in FILMS:
        i = np.argmin(np.abs(HALF - opd(T_RET) / 2))
        print(f"막 {d:.0f} nm, 창 반폭 L/2: Psi 오차 {res[d][i,0]:.3f}°, "
              f"Delta 오차 {res[d][i,1]:.3f}°")
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        figure1(L); figure2(L); figure3(L); figure4(L); figure5(L)
        print(f"  [{lang}] 그림 5개 저장 -> {os.path.normpath(L['dir'])}")


if __name__ == "__main__":
    main()

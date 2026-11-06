"""타원계측기 7편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-snapshot-acse-system/      (한국어)
      ../../assets/img/posts/ellipsometer-snapshot-acse-system/en/   (영문)

계산은 acse.py 의 뮬러 곱에서 나오고, 그 구현은 verify_acse.py 가 저널 식 5,
a = cos2Psi, b = sin2Psi sinDelta 의 닫힌 식, 5편의 부호와 대조해 검증한 것이다.
그림은 모두 계산이며 실측 데이터를 쓰지 않는다. 그림 4 의 부품 오차는 가상값이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

import acse as A

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-snapshot-acse-system",
)

SIG = A.sigma_grid(800)
LAM = 1e3 / SIG                      # nm
EDGE = int(0.05 * len(SIG))          # 대역 양 끝 5 % 를 "끝"으로 본다
C_EEM, C_FT, C_TRUE, C_HT = "tab:red", "tab:blue", "black", "tab:green"
BLUR = dict(A.IDEAL, lss_blur=A.DLAM_LSS)

_CACHE: dict = {}


def cached(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


def row(theta_deg, d_nm, err):
    """한 입사각 줄의 검출 신호와 참 기저 신호."""
    t = np.array([np.deg2rad(theta_deg)])
    y = A.detect(SIG, t, d_nm, err)[0]
    y_pi = A.detect(SIG, t, d_nm, dict(err, phi0=err["phi0"] + np.pi))[0]
    return y, 0.5 * (y + y_pi)


def mid(x):
    return x[EDGE:-EDGE]


def ends(x):
    return np.r_[x[:EDGE], x[-EDGE:]]


def shade_edges(ax):
    for lo, hi in ((LAM[EDGE], LAM[0]), (LAM[-1], LAM[-EDGE])):
        ax.axvspan(min(lo, hi), max(lo, hi), color="0.85", alpha=0.6, lw=0)


# ---------------------------------------------------------------------------
# 그림 1 — (파수, 입사각) 결합 신호 영상
# ---------------------------------------------------------------------------

D1 = 300.0


def fig1_data():
    th = A.theta_grid(241)
    img = A.detect(SIG, th, D1, A.IDEAL)
    img_pi = A.detect(SIG, th, D1, dict(A.IDEAL, phi0=np.pi))
    dc = 0.5 * (img + img_pi)
    bs = dict(A.IDEAL, bs_diatt=0.04)
    y0_bs = A.detect(SIG, np.array([0.0]), D1, bs)[0]
    return th, img, img / dc - 1.0, y0_bs


def figure1(L):
    th, img, ac, y0_bs = cached("f1", fig1_data)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.6),
                                   gridspec_kw=dict(width_ratios=(1.05, 1)))
    pc = ax1.pcolormesh(LAM, np.rad2deg(th), ac, cmap="RdBu_r", vmin=-1, vmax=1,
                        shading="auto", rasterized=True)
    cb = fig.colorbar(pc, ax=ax1, pad=0.02)
    cb.set_label(L["f1_cb"], fontsize=9.5, **L["font"])
    ax1.axhline(0, color="k", lw=0.8, ls=":")
    ax1.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["aoi"], fontsize=10.5, **L["font"])
    ax1.set_title(L["f1a"].format(D1), fontsize=11.0, pad=8, **L["font"])

    scale = np.max(img)
    for tdeg, c in ((0, "black"), (25, "tab:blue"), (55, "tab:red")):
        i = int(np.argmin(np.abs(np.rad2deg(th) - tdeg)))
        ax2.plot(LAM, img[i] / scale, color=c, lw=1.1,
                 label=L["f1_row"].format(tdeg))
    ax2.plot(LAM, y0_bs / scale, color="0.45", lw=1.0, ls="--",
             label=L["f1_bs"])
    ax2.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["f1_y"], fontsize=10.5, **L["font"])
    ax2.set_title(L["f1b"], fontsize=11.0, pad=8, **L["font"])
    ax2.set_xlim(LAM.min(), LAM.max())
    ax2.legend(prop={**L["legend"], "size": 8.2}, loc="upper right", framealpha=0.92)
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig1-signal-image.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 기저 신호: 극값 포락선(EEM) 대 푸리에 저역통과(FT)
# ---------------------------------------------------------------------------

TH2, D2 = 50.0, 300.0
T_SCAN = np.array([1000, 1500, 2000, 2500, 3000, 3500, 4092.3, 4600, 5100, 5500])


def dc_errors(theta_deg, d_nm, err, thickness=A.T_RETARDER):
    e = dict(err, t_ret=thickness)
    y, dct = row(theta_deg, d_nm, e)
    L = A.opd(thickness)
    dce = A.eem_baseline(SIG, y, L)[0]
    dcf = A.ft_baseline(SIG, y, L)
    return (dce - dct) / dct, (dcf - dct) / dct


def fig2_scan():
    out = []
    for t in T_SCAN:
        re, rf = dc_errors(TH2, D2, BLUR, t)
        out.append((A.fringe_count(t),
                    np.sqrt(np.mean(mid(re) ** 2)), np.max(np.abs(ends(re))),
                    np.sqrt(np.mean(mid(rf) ** 2)), np.max(np.abs(ends(rf)))))
    return np.array(out)


T2_SECOND = A.T_RETARDER / 2      # 둘째 지연자 (두께 1/2)


def two_retarder_row(theta_deg, d_nm, off1=0.0, off2=0.0):
    """편광 생성단에 지연자를 둘 둔 배치: P(0) -> R(45, phi1) -> R(0, phi2) -> 시편 -> A(-45).

    학위논문 그림 1.9 가 실측으로 보인 맥놀이 무늬를 뮬러 곱으로 다시 만든다.
    """
    th = np.deg2rad(theta_deg)
    n = len(SIG)
    psi, dl, R = A.film_psi_delta(SIG, np.full(n, th), d_nm)
    om = np.full(n, np.deg2rad(90.0))
    p1 = A.retardance(SIG) + off1
    p2 = A.retardance(SIG, T2_SECOND) + off2
    M = A._chain(A.mueller_polarizer(np.full(n, np.deg2rad(-45.0))),
                 A.mueller_rotation(-om), A.mueller_sample(psi, dl), A.mueller_rotation(om),
                 A.mueller_retarder(np.zeros(n), p2),
                 A.mueller_retarder(np.full(n, np.deg2rad(45.0)), p1),
                 A.mueller_polarizer(np.zeros(n)))
    return A.led_spectrum(SIG) * R * M[..., 0, 0]


def fig2_two():
    y = two_retarder_row(TH2, D2)
    q = (0.0, np.pi / 2, np.pi, 3 * np.pi / 2)
    dct = np.mean([two_retarder_row(TH2, D2, o1, o2) for o1 in q for o2 in q], axis=0)
    dce, up, lo = A.eem_baseline(SIG, y)
    return y, dct, dce, up, lo


def figure2(L):
    y, dct = cached("f2row", lambda: row(TH2, D2, BLUR))
    dce, up, lo = A.eem_baseline(SIG, y)
    dcf = A.ft_baseline(SIG, y)
    scan = cached("f2scan", fig2_scan)

    fig, axs = plt.subplots(2, 2, figsize=(12.6, 8.8))
    axs = axs.ravel()
    ax = axs[0]
    s = np.max(y)
    ax.plot(LAM, y / s, color="0.55", lw=0.8, label=L["f2_sig"])
    ax.plot(LAM, up / s, color=C_EEM, lw=0.9, ls=":", label=L["f2_env"])
    ax.plot(LAM, lo / s, color=C_EEM, lw=0.9, ls=":")
    ax.plot(LAM, dct / s, color=C_TRUE, lw=1.6, label=L["f2_true"])
    ax.plot(LAM, dce / s, color=C_EEM, lw=1.2, ls="--", label=L["f2_eem"])
    ax.plot(LAM, dcf / s, color=C_FT, lw=1.2, ls="-.", label=L["f2_ft"])
    ax.set_xlim(LAM.min(), LAM.max())
    ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f1_y"], fontsize=10.5, **L["font"])
    ax.set_title(L["f2a"].format(TH2, D2), fontsize=11.0, pad=8, **L["font"])
    ax.legend(prop={**L["legend"], "size": 7.8}, loc="upper right", framealpha=0.92)
    ax.grid(alpha=0.25)

    ax = axs[1]
    shade_edges(ax)
    ax.semilogy(LAM, np.abs((dce - dct) / dct), color=C_EEM, lw=1.0, label="EEM")
    ax.semilogy(LAM, np.abs((dcf - dct) / dct), color=C_FT, lw=1.0, label="FT")
    ax.set_xlim(LAM.min(), LAM.max())
    ax.set_ylim(1e-6, 3)
    ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f2_rel"], fontsize=10.5, **L["font"])
    ax.set_title(L["f2b"], fontsize=11.0, pad=8, **L["font"])
    ax.text(0.5, 0.80, L["f2_edge"], transform=ax.transAxes, ha="center",
            fontsize=8.4, color="0.35", **L["font"])
    ax.legend(prop={"size": 8.5}, loc="upper center", ncol=2, framealpha=0.92)
    ax.grid(alpha=0.25, which="both")

    ax = axs[2]
    N = scan[:, 0]
    ax.semilogy(N, scan[:, 1], "o-", color=C_EEM, ms=4, lw=1.3, label=L["f2_mid"].format("EEM"))
    ax.semilogy(N, scan[:, 2], "o--", color=C_EEM, ms=4, lw=1.0, mfc="white",
                label=L["f2_end"].format("EEM"))
    ax.semilogy(N, scan[:, 3], "s-", color=C_FT, ms=4, lw=1.3, label=L["f2_mid"].format("FT"))
    ax.semilogy(N, scan[:, 4], "s--", color=C_FT, ms=4, lw=1.0, mfc="white",
                label=L["f2_end"].format("FT"))
    n0 = A.fringe_count()
    ax.axvline(n0, color="0.4", lw=0.9, ls=":")
    ax.text(n0 + 0.5, 1.6e-2, L["f2_this"].format(A.T_RETARDER / 1e3, n0),
            fontsize=8.2, color="0.3", ha="left", **L["font"])
    ax.set_xlabel(L["f2_N"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f2_rel"], fontsize=10.5, **L["font"])
    ax.set_title(L["f2c"], fontsize=11.0, pad=8, **L["font"])
    ax.set_ylim(3e-4, 5)
    ax.legend(prop={**L["legend"], "size": 7.6}, loc="upper right", ncol=2, framealpha=0.92)
    ax.grid(alpha=0.25, which="both")

    ax = axs[3]
    y2, dct2, dce2, up2, lo2 = cached("f2two", fig2_two)
    s2 = np.max(y2)
    ax.plot(LAM, y2 / s2, color="0.55", lw=0.7, label=L["f2_sig"])
    ax.plot(LAM, up2 / s2, color=C_EEM, lw=0.9, ls=":", label=L["f2_env"])
    ax.plot(LAM, lo2 / s2, color=C_EEM, lw=0.9, ls=":")
    ax.plot(LAM, dct2 / s2, color=C_TRUE, lw=1.6, label=L["f2_true"])
    ax.plot(LAM, dce2 / s2, color=C_EEM, lw=1.2, ls="--", label=L["f2_eem2"])
    ax.set_xlim(LAM.min(), LAM.max())
    ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f1_y"], fontsize=10.5, **L["font"])
    ax.set_title(L["f2d"], fontsize=11.0, pad=8, **L["font"])
    ax.legend(prop={**L["legend"], "size": 7.8}, loc="upper right", framealpha=0.92)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig2-baseline-eem-ft.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — 힐베르트 복조와 부호, 경계 오차
# ---------------------------------------------------------------------------

TH3, D3 = 40.0, 1000.0


def fig3_data():
    y, dct = row(TH3, D3, A.IDEAL)
    t = np.array([np.deg2rad(TH3)])
    a0, b0 = [v[0] for v in A.ideal_ab(SIG, t, D3)]
    phin = A.retardance(SIG)
    dce = A.eem_baseline(SIG, y)[0]
    ac = y / dce - 1.0
    from scipy.signal import hilbert
    z = hilbert(ac)
    zc = z * np.exp(-1j * phin)
    af, bf = A.ft_demod(SIG, y, phin)
    return dict(ac=ac, z=z, zc=zc, a0=a0, b0=b0, af=af, bf=bf)


def figure3(L):
    d = cached("f3", fig3_data)
    fig, axs = plt.subplots(1, 3, figsize=(15.2, 4.5))
    ax = axs[0]
    ax.plot(LAM, d["ac"], color="0.45", lw=0.8, label=L["f3_ac"])
    ax.plot(LAM, np.abs(d["z"]), color=C_HT, lw=1.4, label=L["f3_env"])
    ax.plot(LAM, -np.abs(d["z"]), color=C_HT, lw=1.4)
    ax.set_xlim(LAM.min(), LAM.max())
    ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f3_y"], fontsize=10.5, **L["font"])
    ax.set_title(L["f3a"].format(TH3, D3), fontsize=11.0, pad=8, **L["font"])
    ax.legend(prop={**L["legend"], "size": 8.2}, loc="lower right", framealpha=0.92)
    ax.grid(alpha=0.25)

    ax = axs[1]
    k = slice(None, None, 14)
    ax.plot(LAM, d["a0"], color="tab:blue", lw=2.2, alpha=0.35, label=L["f3_a0"])
    ax.plot(LAM, d["b0"], color="tab:red", lw=2.2, alpha=0.35, label=L["f3_b0"])
    ax.plot(LAM[k], np.real(d["zc"])[k], "o", color="tab:blue", ms=3.4, label=L["f3_re"])
    ax.plot(LAM[k], -np.imag(d["zc"])[k], "o", color="tab:red", ms=3.4, label=L["f3_mim"])
    ax.plot(LAM, np.imag(d["zc"]), color="tab:red", lw=0.9, ls="--", label=L["f3_pim"])
    ax.set_xlim(LAM.min(), LAM.max())
    ax.set_ylim(-1.15, 1.15)
    ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f3_ab"], fontsize=10.5, **L["font"])
    ax.set_title(L["f3b"], fontsize=11.0, pad=8, **L["font"])
    ax.legend(prop={**L["legend"], "size": 7.6}, loc="lower left", ncol=2, framealpha=0.92)
    ax.grid(alpha=0.25)

    ax = axs[2]
    shade_edges(ax)
    eh = np.hypot(np.real(d["zc"]) - d["a0"], -np.imag(d["zc"]) - d["b0"])
    ef = np.hypot(d["af"] - d["a0"], d["bf"] - d["b0"])
    ax.semilogy(LAM, eh, color=C_EEM, lw=1.0, label="EEM + HT")
    ax.semilogy(LAM, ef, color=C_FT, lw=1.0, label=L["f3_ft"])
    ax.set_xlim(LAM.min(), LAM.max())
    ax.set_ylim(1e-5, 2)
    ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f3_err"], fontsize=10.5, **L["font"])
    ax.set_title(L["f3c"], fontsize=11.0, pad=8, **L["font"])
    ax.legend(prop={**L["legend"], "size": 8.2}, loc="upper center", framealpha=0.92)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig3-hilbert-demod.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — 진폭·위상 오차 지도와 보정 전후
# ---------------------------------------------------------------------------

SIG4 = A.sigma_grid(600)
TH4 = A.theta_grid(61)
D_REF, D_TEST = 0.0, 300.0
E_ALL = A.VIRTUAL_ERRORS
E_PV = dict(A.IDEAL, phi0=E_ALL["phi0"], kappa=E_ALL["kappa"],
            depol=E_ALL["depol"], lss_blur=E_ALL["lss_blur"])   # 위상·가시도 오차만


ERRS = {"all": E_ALL, "pv": E_PV}


def measure(d_nm, key):
    return cached(("m", d_nm, key), lambda: A.analyze_image(
        SIG4, A.detect(SIG4, TH4, d_nm, ERRS[key])))


def calibrate(key):
    def run():
        am, bm = measure(D_REF, key)
        ar, br = A.ideal_ab(SIG4, TH4, D_REF)
        return A.error_maps(am, bm, ar, br)
    return cached(("cal", key), run)


def fig4_data():
    out = {}
    for key in ERRS:
        ea, ec = calibrate(key)
        a, b = measure(D_TEST, key)
        a0, b0 = A.ideal_ab(SIG4, TH4, D_TEST)
        ac, bc = A.apply_correction(a, b, ea, ec)
        e4 = int(0.05 * len(SIG4))
        rms = lambda r: np.sqrt(np.nanmean(r[:, e4:-e4] ** 2, axis=1))
        out[key] = dict(ea=ea, ec=ec, before=rms(np.hypot(a - a0, b - b0)),
                        after=rms(np.hypot(ac - a0, bc - b0)))
    return out


def figure4(L):
    d = cached("f4", fig4_data)
    lam4 = 1e3 / SIG4
    thd = np.rad2deg(TH4)
    fig, axs = plt.subplots(1, 3, figsize=(15.6, 4.6),
                            gridspec_kw=dict(width_ratios=(1, 1, 1.05)))
    ea = np.where(np.abs(thd)[:, None] < 3.0, np.nan, d["all"]["ea"])
    pc = axs[0].pcolormesh(lam4, thd, ea, cmap="magma",
                           norm=LogNorm(vmin=0.5, vmax=20), shading="auto",
                           rasterized=True)
    cb = fig.colorbar(pc, ax=axs[0], pad=0.02)
    cb.set_label("$\\varepsilon_A$", fontsize=10.5)
    axs[0].set_title(L["f4a"], fontsize=11.0, pad=8, **L["font"])

    ec = np.where(np.abs(thd)[:, None] < 3.0, np.nan, np.rad2deg(d["all"]["ec"]))
    pc = axs[1].pcolormesh(lam4, thd, ec, cmap="coolwarm", vmin=-60, vmax=60,
                           shading="auto", rasterized=True)
    cb = fig.colorbar(pc, ax=axs[1], pad=0.02)
    cb.set_label("$\\varepsilon_\\chi$ (deg)", fontsize=10.5)
    axs[1].set_title(L["f4b"], fontsize=11.0, pad=8, **L["font"])
    for ax in axs[:2]:
        ax.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
        ax.set_ylabel(L["aoi"], fontsize=10.5, **L["font"])
        ax.text(0.5, 0.5, L["f4_mask"], transform=ax.transAxes, ha="center",
                va="center", fontsize=7.8, color="0.25", **L["font"])

    ax = axs[2]
    sel = np.abs(thd) >= 3.0
    pos = sel & (thd > 0)
    ax.semilogy(thd[pos], d["all"]["before"][pos], "o-", color="0.35", ms=3.5,
                lw=1.2, label=L["f4_before"])
    ax.semilogy(thd[pos], d["all"]["after"][pos], "o-", color=C_EEM, ms=3.5,
                lw=1.4, label=L["f4_after"])
    ax.semilogy(thd[pos], d["pv"]["after"][pos], "s--", color=C_HT, ms=3.2,
                lw=1.1, label=L["f4_pv"])
    ax.set_xlabel(L["aoi"], fontsize=10.5, **L["font"])
    ax.set_ylabel(L["f4_y"], fontsize=10.5, **L["font"])
    ax.set_title(L["f4c"].format(D_TEST), fontsize=11.0, pad=8, **L["font"])
    ax.set_xlim(0, 61)
    ax.legend(prop={**L["legend"], "size": 7.8}, loc="lower right", framealpha=0.92)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig4-error-maps.png"), dpi=150)
    plt.close(fig)


def thickness_table():
    """보정 전후 두께 회귀. 기준 시편은 맨 실리콘, 측정 시편은 SiO2 세 두께.

    회귀는 입사각 8° 이상, 대역 양 끝 5 % 를 뺀 영역에서 격자를 솎아(입사각 둘에 하나,
    파수 셋에 하나) 한다. 계산 시간을 줄이려는 것이다.
    """
    e4 = int(0.05 * len(SIG4))
    core = np.abs(TH4)[:, None] >= np.deg2rad(8.0)
    core = core & np.ones((len(TH4), len(SIG4)), bool)
    core[:, :e4] = False
    core[:, -e4:] = False
    ti, si = slice(None, None, 2), slice(None, None, 3)
    sub = lambda x: x[ti, si]
    rows = []
    for key in ERRS:
        ea, ec = calibrate(key)
        for dt in (100.0, 300.0, 1000.0):
            a, b = measure(dt, key)
            ac, bc = A.apply_correction(a, b, ea, ec)
            fb = A.fit_thickness(SIG4[si], TH4[ti], sub(a), sub(b), dt, mask=sub(core))[0]
            m = sub(core & np.isfinite(ac) & np.isfinite(bc))
            fa = A.fit_thickness(SIG4[si], TH4[ti], sub(ac), sub(bc), dt, mask=m)[0]
            rows.append((key, dt, fb, fa))
    return rows


# ---------------------------------------------------------------------------
# 그림 5 — LSS 분해능과 지연자 두께
# ---------------------------------------------------------------------------


def figure5(L):
    t_min, t_max = A.thickness_window()
    t = np.linspace(200, 8000, 800)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.5))
    ax1.axvspan(t_min / 1e3, t_max / 1e3, color="tab:green", alpha=0.10, lw=0)
    for lam, c in ((0.475, "tab:blue"), (0.575, "tab:green"), (0.675, "tab:red")):
        ax1.plot(t / 1e3, A.visibility_box(lam, t), color=c, lw=1.5,
                 label=L["f5_v"].format(lam * 1e3))
    for x in (t_min, t_max):
        ax1.axvline(x / 1e3, color="tab:green", lw=1.0, ls="--")
    ax1.axvline(A.T_RETARDER / 1e3, color="k", lw=1.0, ls=":")
    ax1.text(A.T_RETARDER / 1e3 + 0.08, 0.06, L["f5_this"].format(A.T_RETARDER / 1e3),
             fontsize=8.4, **L["font"])
    ax1.text(t_min / 1e3 + 0.08, 0.30, L["f5_lo"].format(t_min / 1e3), fontsize=8.2,
             color="tab:green", **L["font"])
    ax1.text(t_max / 1e3 + 0.08, 0.30, L["f5_hi"].format(t_max / 1e3), fontsize=8.2,
             color="tab:green", **L["font"])
    ax1.set_xlim(0.2, 8.0)
    ax1.set_ylim(0, 1.05)
    ax1.set_xlabel(L["f5_t"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["f5_vy"], fontsize=10.5, **L["font"])
    ax1.set_title(L["f5a"], fontsize=11.0, pad=22, **L["font"])
    ax1.legend(prop={**L["legend"], "size": 8.2}, loc="lower left", framealpha=0.92)
    ax1.grid(alpha=0.25)
    top = ax1.secondary_xaxis(
        "top", functions=(lambda x: A.fringe_count(x * 1e3),
                          lambda n: n / A.fringe_count(1e3)))
    top.set_xlabel(L["f5_N"], fontsize=9.5, **L["font"])

    lam = np.linspace(475, 675, 300)
    period = (lam * 1e-3) ** 2 / A.opd() * 1e3
    ax2.plot(lam, period, color="k", lw=1.6, label=L["f5_per"])
    ax2.axhline(A.DLAM_LSS * 1e3, color="tab:purple", lw=1.4, ls="--",
                label=L["f5_res"].format(A.DLAM_LSS * 1e3))
    ax2.set_xlim(475, 675)
    ax2.set_ylim(0, 14)
    ax2.set_xlabel(L["lam"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["f5_py"], fontsize=10.5, **L["font"])
    ax2b = ax2.twinx()
    ax2b.plot(lam, A.visibility_box(lam * 1e-3, A.T_RETARDER), color="tab:orange",
              lw=1.4, label=L["f5_vis"])
    ax2b.set_ylim(0.6, 1.0)
    ax2b.set_ylabel(L["f5_vy"], fontsize=10.5, color="tab:orange", **L["font"])
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2b.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, prop={**L["legend"], "size": 8.2}, loc="upper left",
               framealpha=0.92)
    ax2.set_title(L["f5b"].format(A.T_RETARDER / 1e3), fontsize=11.0, pad=22, **L["font"])
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(os.path.join(L["dir"], "fig5-resolution-window.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 라벨
# ---------------------------------------------------------------------------

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "lam": "파장 (nm)", "aoi": "입사각 (deg)",
        # 그림 1
        "f1a": "(a) 정규화한 무늬 (SiO$_2$ {:.0f} nm / Si)",
        "f1_cb": "기저 신호로 나눈 무늬",
        "f1b": "(b) 입사각별 검출 신호",
        "f1_row": "입사각 {}°", "f1_bs": "0°, 빔 스플리터 이색성 0.04",
        "f1_y": "검출 세기 (최댓값 = 1)",
        # 그림 2
        "f2a": "(a) 포락선과 기저 신호 (입사각 {:.0f}°, 막 {:.0f} nm)",
        "f2_sig": "검출 신호", "f2_env": "위·아래 포락선",
        "f2_true": "참 기저 신호", "f2_eem": "EEM", "f2_ft": "FT 저역통과",
        "f2b": "(b) 기저 신호의 상대 오차",
        "f2_rel": "상대 오차", "f2_edge": "회색 띠: 대역 양 끝 5 %",
        "f2c": "(c) 무늬 수에 따른 오차",
        "f2_N": "대역 안의 무늬 수 $N$",
        "f2_mid": "{} 가운데 (RMS)", "f2_end": "{} 양 끝 (최대)",
        "f2_this": "이 장비\n{:.2f} mm, N = {:.1f}",
        "f2d": "(d) 지연자 둘: 맥놀이 무늬의 포락선",
        "f2_eem2": "포락선 평균",
        # 그림 3
        "f3a": "(a) 무늬와 해석 신호의 크기 ({:.0f}°, {:.0f} nm)",
        "f3_ac": "정규화 무늬 $AC'$", "f3_env": "$\\pm|Z| = \\pm A'$",
        "f3_y": "기저 신호로 나눈 무늬",
        "f3b": "(b) 반송파를 뺀 $Z$에서 읽은 $a$, $b$",
        "f3_a0": "참 $a$", "f3_b0": "참 $b$",
        "f3_re": "실수부 $\\to a$", "f3_mim": "허수부 부호 반전 $\\to b$",
        "f3_pim": "허수부 그대로 ($= -b$)",
        "f3_ab": "계수 값",
        "f3c": "(c) 복조 오차", "f3_ft": "FT 채널 복조 (5편)",
        "f3_err": "$|(\\Delta a, \\Delta b)|$",
        # 그림 4
        "f4a": "(a) 진폭 오차 지도", "f4b": "(b) 위상 오차 지도",
        "f4_mask": "",
        "f4c": "(c) SiO$_2$ {:.0f} nm 막의 잔차",
        "f4_before": "보정 전", "f4_after": "보정 후",
        "f4_pv": "보정 후 (위상·가시도 오차만 있을 때)",
        "f4_y": "$(a, b)$ 잔차 RMS (대역 가운데)",
        # 그림 5
        "f5a": "(a) 두께 창과 무늬 가시도",
        "f5_v": "{:.0f} nm", "f5_t": "지연자 두께 $t$ (mm)",
        "f5_vy": "무늬 가시도", "f5_N": "대역 안의 무늬 수 $N$",
        "f5_this": "{:.2f} mm", "f5_lo": "하한\n{:.2f} mm", "f5_hi": "상한\n{:.2f} mm",
        "f5b": "(b) 무늬 한 주기와 분해능 ($t$ = {:.2f} mm)",
        "f5_per": "무늬 한 주기", "f5_res": "LSS 분해능 {:.2f} nm",
        "f5_py": "파장 폭 (nm)", "f5_vis": "무늬 가시도",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": {}, "legend": {},
        "lam": "wavelength (nm)", "aoi": "angle of incidence (deg)",
        "f1a": "(a) normalized fringes (SiO$_2$ {:.0f} nm / Si)",
        "f1_cb": "fringe divided by the baseline",
        "f1b": "(b) detected signal at each angle",
        "f1_row": "{}° incidence", "f1_bs": "0°, beam-splitter diattenuation 0.04",
        "f1_y": "detected intensity (max = 1)",
        "f2a": "(a) envelopes and baseline ({:.0f}°, {:.0f} nm film)",
        "f2_sig": "detected signal", "f2_env": "upper and lower envelopes",
        "f2_true": "true baseline", "f2_eem": "EEM", "f2_ft": "FT low-pass",
        "f2b": "(b) relative error of the baseline",
        "f2_rel": "relative error", "f2_edge": "grey bands: outer 5 % of the band",
        "f2c": "(c) error versus fringe count",
        "f2_N": "fringes in the band, $N$",
        "f2_mid": "{} centre (RMS)", "f2_end": "{} ends (max)",
        "f2_this": "this instrument\n{:.2f} mm, N = {:.1f}",
        "f2d": "(d) two retarders: envelopes of a beat pattern",
        "f2_eem2": "mean of the envelopes",
        "f3a": "(a) fringe and analytic-signal magnitude ({:.0f}°, {:.0f} nm)",
        "f3_ac": "normalized fringe $AC'$", "f3_env": "$\\pm|Z| = \\pm A'$",
        "f3_y": "fringe divided by the baseline",
        "f3b": "(b) $a$, $b$ read from $Z$ with the carrier removed",
        "f3_a0": "true $a$", "f3_b0": "true $b$",
        "f3_re": "real part $\\to a$", "f3_mim": "imaginary part, sign flipped $\\to b$",
        "f3_pim": "imaginary part as is ($= -b$)",
        "f3_ab": "coefficient",
        "f3c": "(c) demodulation error", "f3_ft": "FT channel demodulation (Part 5)",
        "f3_err": "$|(\\Delta a, \\Delta b)|$",
        "f4a": "(a) amplitude error map", "f4b": "(b) phase error map",
        "f4_mask": "",
        "f4c": "(c) residual for a {:.0f} nm SiO$_2$ film",
        "f4_before": "before correction", "f4_after": "after correction",
        "f4_pv": "after correction (phase and visibility errors only)",
        "f4_y": "RMS residual of $(a, b)$ (band centre)",
        "f5a": "(a) thickness window and fringe visibility",
        "f5_v": "{:.0f} nm", "f5_t": "retarder thickness $t$ (mm)",
        "f5_vy": "fringe visibility", "f5_N": "fringes in the band, $N$",
        "f5_this": "{:.2f} mm", "f5_lo": "lower limit\n{:.2f} mm",
        "f5_hi": "upper limit\n{:.2f} mm",
        "f5b": "(b) one fringe period and the resolution ($t$ = {:.2f} mm)",
        "f5_per": "one fringe period", "f5_res": "LSS resolution {:.2f} nm",
        "f5_py": "wavelength width (nm)", "f5_vis": "fringe visibility",
    },
}


def ab_error_blur(theta_deg=50.0, d_nm=100.0):
    """분해능 유무와 기저 신호(EEM / 참값)에 따른 (a, b) 오차. 8절 본문의 수치."""
    t = np.array([np.deg2rad(theta_deg)])
    a0, b0 = [v[0] for v in A.ideal_ab(SIG, t, d_nm)]
    phin = A.retardance(SIG)
    out = {}
    for name, err in (("분해능 없음", A.IDEAL), ("분해능 1.92 nm", BLUR)):
        y, dct = row(theta_deg, d_nm, err)
        for dname, dc in (("EEM", A.eem_baseline(SIG, y)[0]), ("참 기저", dct)):
            a, b, _ = A.ht_demod(SIG, y, dc, phin)
            out[f"{name}/{dname}"] = round(float(np.sqrt(np.mean(mid(np.hypot(a - a0, b - b0)) ** 2))), 4)
    return out


def phase_map_parts():
    """위상 지도 가운데 값을 부품을 하나씩 빼며 본다(입사각 15°, 50°). 7절 본문의 수치."""
    th = np.deg2rad(np.array([15.0, 50.0]))
    ar, br = A.ideal_ab(SIG4, th, D_REF)
    c = len(SIG4) // 2
    out = {}
    for name, err in (("전부", E_ALL),
                      ("빔 스플리터 뺌", dict(E_ALL, bs_diatt=0.0, bs_ret=0.0)),
                      ("빔 스플리터·대물렌즈 뺌", dict(E_ALL, bs_diatt=0.0, bs_ret=0.0,
                                                   obj_ret=0.0, obj_diatt=0.0))):
        am, bm = A.analyze_image(SIG4, A.detect(SIG4, th, D_REF, err))
        _, ec = A.error_maps(am, bm, ar, br)
        out[name] = [round(float(np.rad2deg(ec[i, c])), 1) for i in range(2)]
    return out


def report():
    """본문에 쓰는 수치를 출력한다."""
    print(f"광로차 L = {A.opd():.2f} um, 무늬 수 N = {A.fringe_count():.1f}")
    t_min, t_max = A.thickness_window()
    print(f"두께 창 {t_min/1e3:.2f} ~ {t_max/1e3:.2f} mm, 이 장비 {A.T_RETARDER/1e3:.4f} mm")
    for lam in (0.475, 0.575, 0.675):
        print(f"  {lam*1e3:.0f} nm: 무늬 주기 {lam**2/A.opd()*1e3:.2f} nm, "
              f"가시도 {A.visibility_box(lam, A.T_RETARDER):.3f}")
    print(f"  상한 두께에서 475 nm 가시도 {A.visibility_box(0.475, t_max):.3f}")

    th, img, ac, y0_bs = cached("f1", fig1_data)
    i0 = int(np.argmin(np.abs(th)))
    print(f"그림1 0° 줄 무늬 진폭(이상) {np.ptp(ac[i0])/2:.2e}, "
          f"빔 스플리터 이색성 0.04 일 때 {np.ptp(y0_bs/np.mean(y0_bs))/2:.3f}(정규화 반진폭, 광원 모양 포함)")
    sym = np.max(np.abs(ac - ac[::-1]))
    print(f"  +-theta 대칭 최대 차 {sym:.1e}")

    y, dct = cached("f2row", lambda: row(TH2, D2, BLUR))
    dce = A.eem_baseline(SIG, y)[0]
    dcf = A.ft_baseline(SIG, y)
    re, rf = (dce - dct) / dct, (dcf - dct) / dct
    print(f"그림2 ({TH2:.0f}°, {D2:.0f} nm): EEM 가운데 RMS {np.sqrt(np.mean(mid(re)**2)):.1e}, "
          f"끝 최대 {np.max(np.abs(ends(re))):.1e} / FT {np.sqrt(np.mean(mid(rf)**2)):.1e}, "
          f"{np.max(np.abs(ends(rf))):.1e}")
    y2, dct2, dce2, _, _ = cached("f2two", fig2_two)
    r2 = (dce2 - dct2) / dct2
    print(f"  지연자 둘(두께 1:1/2): 포락선 평균의 상대 오차 가운데 RMS {np.sqrt(np.mean(mid(r2)**2)):.1e}, "
          f"최대 {np.max(np.abs(mid(r2))):.2f}")
    for r in cached("f2scan", fig2_scan):
        print(f"  N={r[0]:5.1f}: EEM {r[1]:.1e}/{r[2]:.1e}  FT {r[3]:.1e}/{r[4]:.1e}")
    print("  입사각별 (막 300 nm, 분해능 포함): EEM 가운데 / FT 가운데")
    for tdeg in (2, 4, 6, 8, 10, 12, 15, 20, 30, 60):
        re, rf = dc_errors(tdeg, D2, BLUR)
        print(f"   {tdeg:3d}°: {np.sqrt(np.mean(mid(re)**2)):.1e} / {np.sqrt(np.mean(mid(rf)**2)):.1e}")

    d = cached("f3", fig3_data)
    eh = np.hypot(np.real(d["zc"]) - d["a0"], -np.imag(d["zc"]) - d["b0"])
    ef = np.hypot(d["af"] - d["a0"], d["bf"] - d["b0"])
    print(f"그림3 ({TH3:.0f}°, {D3:.0f} nm): EEM+HT 가운데 RMS {np.sqrt(np.mean(mid(eh)**2)):.1e}, "
          f"끝 최대 {np.max(ends(eh)):.2f} / FT {np.sqrt(np.mean(mid(ef)**2)):.1e}, {np.max(ends(ef)):.2f}")
    print(f"  Re-a 최대 {np.max(np.abs(mid(np.real(d['zc'])-d['a0']))):.1e}, "
          f"(+Im)+b 최대 {np.max(np.abs(mid(np.imag(d['zc'])+d['b0']))):.1e}")

    print(f"8절 가시도의 몫 (50°, 100 nm, a·b 오차 가운데 RMS): {ab_error_blur()}")
    print(f"그림4(b) 위상 지도의 부품별 몫: {phase_map_parts()}")
    f4 = cached("f4", fig4_data)
    thd = np.rad2deg(TH4)
    for key in ("all", "pv"):
        sel = (np.abs(thd) >= 8)
        print(f"그림4 [{key}] 잔차 RMS 보정 전 {np.sqrt(np.mean(f4[key]['before'][sel]**2)):.3f}, "
              f"후 {np.sqrt(np.mean(f4[key]['after'][sel]**2)):.4f} (|theta|>=8°)")
    ea = f4["all"]["ea"]
    for tdeg in (5, 10, 20, 40, 60):
        i = int(np.argmin(np.abs(thd - tdeg)))
        print(f"   {tdeg}°: eps_A 중앙값 {np.nanmedian(ea[i]):.2f}, "
              f"eps_chi 범위 {np.rad2deg(np.nanmin(f4['all']['ec'][i])):.1f} ~ "
              f"{np.rad2deg(np.nanmax(f4['all']['ec'][i])):.1f}°")
    print("두께 회귀 (기준: 맨 Si)")
    for key, dt, fb, fa in thickness_table():
        print(f"  [{key}] {dt:6.0f} nm: 보정 전 {fb:8.2f} ({fb-dt:+.2f}), 후 {fa:8.2f} ({fa-dt:+.2f})")


def main():
    import sys
    figs = {"1": figure1, "2": figure2, "3": figure3, "4": figure4, "5": figure5}
    sel = [a for a in sys.argv[1:] if a in figs] or list(figs)
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for k in sel:
            figs[k](L)
        print(f"  [{lang}] 그림 {', '.join(sel)} 저장 -> {os.path.normpath(L['dir'])}")
    if "--no-report" not in sys.argv:
        report()


if __name__ == "__main__":
    main()

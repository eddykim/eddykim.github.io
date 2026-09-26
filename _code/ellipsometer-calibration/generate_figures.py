"""타원계측기 4편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometer-calibration/      (한국어)
      ../../assets/img/posts/ellipsometer-calibration/en/   (영문)

계산은 calibration.py 의 뮬러 곱에서 나오고, 그 구현은 verify_calibration.py 가
Fujiwara 식 4.59-4.61 과 Johs 식 2-5 의 닫힌 식과 대조해 검증한 것이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from calibration import (invert_ideal, measure, regression_calibration,
                         residual, residual_calibration, waveform)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometer-calibration",
)

D = np.deg2rad
C_IDEAL, C_PS, C_AS, C_ETA = "#2c3e50", "tab:blue", "tab:orange", "tab:green"
C_FIT, C_BAD = "tab:red", "#c0392b"

# 기준 시편 — Johs 표 1 의 SiO2/Si 와 비슷한 자리
PSI0, DELTA0 = D(16.8), D(87.0)
# 그림 1 전용 시편. (alpha', beta') 자리는 반축이 1 과 |cos Delta| 인 타원이므로,
# 기준 시편(|cos 87°| = 0.05)으로 그리면 납작해서 기하가 보이지 않는다.
PSI1, DELTA1 = D(35.0), D(50.0)
NOISE = 2e-4            # 검출 세기의 가우스 잡음 표준편차 (dc 는 0.25)
N_TRIAL = 40

_CACHE: dict = {}


def cached(key, fn):
    """두 언어가 같은 수치를 쓰므로 계산은 한 번만 한다."""
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


# ---------------------------------------------------------------------------
# 그림 1 — 세 오차가 신호에 하는 일
# ---------------------------------------------------------------------------


def locus(ps=0.0, a_s=0.0, eta=1.0, n=181):
    """편광자를 한 바퀴 돌릴 때 (alpha', beta') 가 그리는 자리."""
    P = np.linspace(-np.pi / 2, np.pi / 2, n)
    ab = np.array([measure(PSI1, DELTA1, p, ps=ps, a_s=a_s, eta=eta) for p in P])
    return ab[:, 0], ab[:, 1]


def figure1(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.7))

    # (a) 분석기 한 바퀴 동안의 파형
    A, I0 = cached("wf_ideal", lambda: waveform(PSI1, DELTA1, D(45.0)))
    sets = [(L["ideal"], {}, C_IDEAL, "-", 2.0),
            (L["e_ps"], dict(ps=D(4.0)), C_PS, "--", 1.6),
            (L["e_as"], dict(a_s=D(8.0)), C_AS, "-.", 1.6),
            (L["e_eta"], dict(eta=0.7), C_ETA, ":", 1.9)]
    for lab, kw, c, ls, lw in sets:
        I = cached(f"wf_{kw}", lambda kw=kw: waveform(PSI1, DELTA1, D(45.0), **kw))[1]
        ax1.plot(np.rad2deg(A), I / I0.mean(), ls, color=c, lw=lw, label=lab)
    ax1.set_xlabel(L["ana_ax"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["norm_int"], fontsize=10.5, **L["font"])
    ax1.set_title(L["fig1a"], fontsize=11.5, pad=8, **L["font"])
    ax1.set_xlim(0, 180)
    ax1.set_ylim(0.0, 2.55)
    ax1.set_xticks(np.arange(0, 181, 45))
    ax1.legend(fontsize=8.0, prop=L["legend"], loc="upper center", ncol=2,
               framealpha=0.9, columnspacing=0.8, handlelength=2.0,
               borderaxespad=0.3, handletextpad=0.5)
    ax1.grid(alpha=0.25)

    # (b) (alpha', beta') 평면에서의 자리
    for lab, kw, c, ls, lw in sets:
        if "ps" in kw:          # Ps 는 타원을 바꾸지 않는다. 점의 위치만 민다
            continue
        a, b = cached(f"loc_{kw}", lambda kw=kw: locus(**kw))
        ax2.plot(a, b, ls, color=c, lw=lw, label=lab)
    ax2.plot(*measure(PSI1, DELTA1, D(45.0)), "o", color=C_IDEAL, ms=6, zorder=5)
    ax2.plot(*measure(PSI1, DELTA1, D(45.0), ps=D(4.0)), "o", color=C_PS, ms=6,
             zorder=5)
    ax2.axhline(0, color="gray", lw=0.6)
    ax2.axvline(0, color="gray", lw=0.6)
    ax2.set_aspect("equal")
    ax2.set_xlim(-1.15, 1.15)
    ax2.set_ylim(-1.62, 1.15)
    ax2.set_xlabel(r"$\alpha'$", fontsize=11.5)
    ax2.set_ylabel(r"$\beta'$", fontsize=11.5)
    ax2.set_title(L["fig1b"], fontsize=11.5, pad=8, **L["font"])
    ax2.annotate(L["slide"], xy=measure(PSI1, DELTA1, D(45.0)),
                 xytext=(-1.10, 1.10), ha="left", va="top", fontsize=8.4,
                 color=C_PS, **L["font"], linespacing=1.3,
                 arrowprops=dict(arrowstyle="->", color=C_PS, lw=1.0))
    ax2.text(0.5, 0.02, L["fig1note"], transform=ax2.transAxes, fontsize=7.9,
             ha="center", va="bottom", **L["font"], linespacing=1.35,
             bbox=dict(boxstyle="round,pad=0.35", fc="#fffbe6", ec="#d9c77a", lw=0.8))
    ax2.grid(alpha=0.25)

    fig.suptitle(L["fig1"], fontsize=13.0, y=0.99, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig1-error-anatomy.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 2 — 잔차 보정과 그 포물선 근사
# ---------------------------------------------------------------------------

PS_TRUE, AS_TRUE, ETA_TRUE = D(0.26), D(64.3), 0.9687
HALVES = (5.0, 3.0, 2.0, 1.0, 0.5)


def eta_convergence():
    out = []
    for half in HALVES:
        r = residual_calibration(PSI0, DELTA0,
                                 np.arange(-half, half + 1e-9, half / 20)
                                 + np.rad2deg(PS_TRUE),
                                 ps=PS_TRUE, a_s=AS_TRUE, eta=ETA_TRUE)
        out.append((half, r["eta"], abs(r["eta"] - ETA_TRUE)))
    return np.array(out)


def figure2(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))

    # (a) R(P) 와 포물선 맞춤
    scan = np.arange(-5, 5.01, 0.25) + np.rad2deg(PS_TRUE)
    r = cached("rescal", lambda: residual_calibration(
        PSI0, DELTA0, scan, ps=PS_TRUE, a_s=AS_TRUE, eta=ETA_TRUE))
    Pd = np.rad2deg(r["P"])
    ax1.plot(Pd, r["R"], "o", color=C_IDEAL, ms=3.0, label=L["r_meas"])
    c = np.polyfit(r["P"], r["R"], 2)
    fine = np.linspace(r["P"][0], r["P"][-1], 400)
    ax1.plot(np.rad2deg(fine), np.polyval(c, fine), "-", color=C_FIT, lw=1.6,
             label=L["r_fit"])
    ps_d, r_min = np.rad2deg(r["ps"]), 1.0 - r["eta"] ** 2
    ax1.axvline(ps_d, color=C_PS, lw=1.2, ls="--")
    ax1.axhline(r_min, color=C_ETA, lw=1.2, ls="--")
    ax1.plot([ps_d], [r_min], "*", color=C_FIT, ms=14, zorder=6)
    ax1.annotate(L["a_ps"].format(ps_d), xy=(ps_d, r["R"].max() * 0.86),
                 xytext=(ps_d - 5.0, r["R"].max() * 0.99), fontsize=9.2,
                 color=C_PS, ha="left", **L["font"],
                 arrowprops=dict(arrowstyle="->", color=C_PS, lw=1.0))
    ax1.annotate(L["a_eta"].format(r["eta"]), xy=(Pd[-1], r_min),
                 xytext=(Pd[-1] - 4.2, r_min + 0.012 * r["R"].max() / 0.06),
                 fontsize=9.2, color=C_ETA, **L["font"],
                 arrowprops=dict(arrowstyle="->", color=C_ETA, lw=1.0))
    ax1.set_xlabel(L["pol_ax"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["res_ax"], fontsize=10.5, **L["font"])
    ax1.set_title(L["fig2a"], fontsize=11.5, pad=8, **L["font"])
    ax1.legend(fontsize=9.0, prop=L["legend"], loc="lower right", framealpha=0.9)
    ax1.grid(alpha=0.25)

    # (b) 포물선 근사가 남기는 eta 오차
    conv = cached("conv", eta_convergence)
    ax2.loglog(conv[:, 0], conv[:, 2], "o-", color=C_FIT, lw=1.8, ms=6)
    for half, _, err in conv:
        ax2.annotate(f"{err:.0e}".replace("e-0", "e-"), xy=(half, err),
                     xytext=(4, 6), textcoords="offset points", fontsize=8.2)
    ax2.set_xlabel(L["half_ax"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["eta_err_ax"], fontsize=10.5, **L["font"])
    ax2.set_title(L["fig2b"], fontsize=11.5, pad=8, **L["font"])
    ax2.invert_xaxis()
    ax2.set_xticks(list(HALVES))
    ax2.set_xticklabels([f"{h:g}" for h in HALVES])
    ax2.minorticks_off()
    ax2.text(0.04, 0.06, L["fig2note"], transform=ax2.transAxes, fontsize=9.0,
             va="bottom", **L["font"],
             bbox=dict(boxstyle="round,pad=0.35", fc="#fffbe6", ec="#d9c77a", lw=0.8))
    ax2.grid(alpha=0.25, which="both")

    fig.suptitle(L["fig2"], fontsize=13.0, y=0.99, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig3-residual-calibration.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 3 — Delta 가 0·180° 근처일 때의 맹점
# ---------------------------------------------------------------------------

PSI_BAD = D(27.0)
DELTAS = (90.0, 60.0, 30.0, 10.0, 3.0, 1.0)


def blind_spot_data():
    scan = np.arange(-5, 5.01, 0.25) + 0.26
    curves, rows = [], []
    for d_deg in DELTAS:
        R = np.array([residual(PSI_BAD, D(d_deg), p, ps=D(0.26))
                      for p in np.deg2rad(scan)])
        curves.append((d_deg, scan, R - R.min()))
        errs, fail = [], 0
        for k in range(N_TRIAL):
            r = residual_calibration(PSI_BAD, D(d_deg),
                                     np.arange(-5, 5.01, 0.5) + 0.26, ps=D(0.26),
                                     noise=NOISE,
                                     rng=np.random.default_rng(1000 * int(d_deg) + k))
            errs.append(abs(np.rad2deg(r["ps"]) - 0.26)) if r["ok"] else None
            fail += 0 if r["ok"] else 1
        rows.append((d_deg, float(np.median(errs)) if errs else np.nan, fail))
    # R 에 실리는 잡음 크기를 같은 조건에서 경험적으로 잰다
    ref = np.array([residual(PSI_BAD, D(90.0), 0.0, noise=NOISE,
                             rng=np.random.default_rng(k)) for k in range(400)])
    return curves, np.array(rows), float(ref.std())


def figure3(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))
    curves, rows, sigma_R = cached("blind", blind_spot_data)

    cmap = plt.get_cmap("viridis")
    for i, (d_deg, scan, dR) in enumerate(curves):
        ax1.semilogy(scan, np.maximum(dR, 1e-9), "-", lw=1.7,
                     color=cmap(i / (len(curves) - 1)),
                     label=f"$\\Delta$ = {d_deg:g}°")
    ax1.axhspan(1e-9, sigma_R, color=C_BAD, alpha=0.13)
    ax1.text(0.5, 0.035, L["noise_band"].format(sigma_R),
             transform=ax1.transAxes, fontsize=8.4, ha="center", va="bottom",
             color=C_BAD, **L["font"])
    ax1.set_xlabel(L["pol_ax"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["depth_ax"], fontsize=10.5, **L["font"])
    ax1.set_title(L["fig3a"], fontsize=11.5, pad=8, **L["font"])
    ax1.set_ylim(1e-8, 8.0)
    ax1.legend(fontsize=8.5, prop=L["legend"], loc="upper center", ncol=3,
               framealpha=0.9)
    ax1.grid(alpha=0.25, which="both")

    ax2.loglog(rows[:, 0], rows[:, 1], "o-", color=C_BAD, lw=1.8, ms=6)
    for d_deg, med, fail in rows:
        if fail:
            ax2.annotate(L["fail"].format(int(fail), N_TRIAL), xy=(d_deg, med),
                         xytext=(-10, 10), textcoords="offset points", fontsize=8.4,
                         ha="right", va="bottom", color=C_BAD, **L["font"])
    ax2.axhline(0.01, color=C_IDEAL, lw=1.1, ls="--")
    ax2.text(0.98, 0.335, L["target"], transform=ax2.transAxes, fontsize=9.0,
             ha="right", va="bottom", color=C_IDEAL, **L["font"])
    ax2.set_xlabel(L["delta_ax"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["ps_err_ax"], fontsize=10.5, **L["font"])
    ax2.set_title(L["fig3b"], fontsize=11.5, pad=8, **L["font"])
    ax2.invert_xaxis()
    ax2.set_xticks(list(DELTAS))
    ax2.set_xticklabels([f"{d:g}" for d in DELTAS])
    ax2.minorticks_off()
    ax2.grid(alpha=0.25, which="both")

    fig.suptitle(L["fig3"], fontsize=13.0, y=0.99, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig4-blind-spot.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 4 — 보정하지 않으면 Psi, Delta 가 얼마나 틀리는가
# ---------------------------------------------------------------------------

DELTA_SWEEP = np.linspace(3.0, 177.0, 175)
CASES = (("e_ps", dict(ps=D(0.1)), C_PS, "-"),
         ("e_as", dict(a_s=D(0.1)), C_AS, "--"),
         ("e_eta", dict(eta=0.999), C_ETA, "-."))


def uncalibrated_error():
    out = {}
    for key, kw, _, _ in CASES:
        dpsi, ddel = [], []
        for d_deg in DELTA_SWEEP:
            a, b = measure(PSI0, D(d_deg), D(45.0), **kw)
            psi, delta = invert_ideal(a, b, D(45.0))
            dpsi.append(np.rad2deg(psi) - np.rad2deg(PSI0))
            ddel.append(np.rad2deg(delta) - d_deg)
        out[key] = (np.array(dpsi), np.array(ddel))
    return out


def figure4(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))
    data = cached("uncal", uncalibrated_error)
    FLOOR = 1e-5          # 로그축이 정확한 0 을 그릴 수 없으므로 바닥을 둔다
    for key, kw, c, ls in CASES:
        dpsi, ddel = data[key]
        ax1.semilogy(DELTA_SWEEP, np.maximum(np.abs(dpsi), FLOOR), ls, color=c,
                     lw=1.7, label=L[key + "4"])
        ax2.semilogy(DELTA_SWEEP, np.maximum(np.abs(ddel), FLOOR), ls, color=c,
                     lw=1.7, label=L[key + "4"])
    # Ps 는 Delta 에 전혀 번지지 않는다 — cos(Delta) 식에서 편광자 각이 약분된다
    ax2.text(90, FLOOR * 1.45, L["zero_note"], fontsize=9.0, ha="center",
             va="bottom", color=C_PS, **L["font"])
    ax2.set_ylim(FLOOR * 0.7, None)
    for ax, ylab, title in ((ax1, L["dpsi_ax"], L["fig4a"]),
                            (ax2, L["ddel_ax"], L["fig4b"])):
        ax.axhline(0.01, color=C_IDEAL, lw=1.1, ls=":")
        ax.set_xlabel(L["delta_true_ax"], fontsize=10.5, **L["font"])
        ax.set_ylabel(ylab, fontsize=10.5, **L["font"])
        ax.set_title(title, fontsize=11.5, pad=8, **L["font"])
        ax.set_xlim(0, 180)
        ax.set_xticks(np.arange(0, 181, 45))
        ax.grid(alpha=0.25, which="both")
    h, lab = ax1.get_legend_handles_labels()
    fig.legend(h, lab, fontsize=9.2, prop=L["legend"], ncol=3, loc="lower center",
               framealpha=0.9, bbox_to_anchor=(0.5, 0.0))

    fig.suptitle(L["fig4"], fontsize=13.0, y=0.99, **L["font"])
    fig.tight_layout(rect=(0, 0.075, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig2-uncalibrated-error.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 그림 5 — 모형에 없는 효과는 eta 로 숨고, MSE 가 그것을 알린다
# ---------------------------------------------------------------------------

# 첫 점은 0 (모형이 완전히 맞는 경우), 나머지는 MSE 가 멱함수로 오르므로 로그 간격
# 0.25 는 본문과 verify_calibration.py 가 인용하는 값이므로 정확히 포함시킨다
PDS_SWEEP = np.unique(np.concatenate([[0.0, 0.25], np.geomspace(0.003, 0.30, 30)]))


def model_mismatch():
    P = D(np.linspace(5.0, 85.0, 90)) + PS_TRUE
    eta, psi, mse = [], [], []
    for pds in PDS_SWEEP:
        a_exp, b_exp = np.array([measure(PSI0, DELTA0, p, ps=PS_TRUE, a_s=AS_TRUE,
                                         eta=1.0, pds=pds) for p in P]).T
        res = regression_calibration(P, a_exp, b_exp,
                                     guess=(D(20), D(80), 0.0, D(60), 0.95))
        eta.append(res["params"][4])
        psi.append(np.rad2deg(res["params"][0]))
        mse.append(res["mse"])
    return np.array(eta), np.array(psi), np.array(mse)


def figure5(L):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 4.5))
    eta, psi, mse = cached("mismatch", model_mismatch)

    ax1.plot(PDS_SWEEP, eta, "-", color=C_ETA, lw=1.9, label=L["eta_fit"])
    ax1.axhline(1.0, color=C_IDEAL, lw=1.1, ls="--")
    ax1.text(0.42, 0.012, L["eta_true"], transform=ax1.transAxes, fontsize=9.0,
             ha="center", va="bottom", color=C_IDEAL, **L["font"])
    axr = ax1.twinx()
    axr.plot(PDS_SWEEP, psi - np.rad2deg(PSI0), "--", color=C_BAD, lw=1.7,
             label=L["psi_shift"])
    axr.set_ylabel(L["dpsi_ax"], fontsize=10.5, color=C_BAD, **L["font"])
    axr.tick_params(axis="y", colors=C_BAD)
    ax1.set_xlabel(L["pds_ax"], fontsize=10.5, **L["font"])
    ax1.set_ylabel(L["eta_fit_ax"], fontsize=10.5, color=C_ETA, **L["font"])
    ax1.tick_params(axis="y", colors=C_ETA)
    ax1.set_title(L["fig5a"], fontsize=11.5, pad=8, **L["font"])
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = axr.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, fontsize=8.8, prop=L["legend"], loc="center left",
               framealpha=0.9)
    ax1.grid(alpha=0.25)

    ax2.loglog(PDS_SWEEP[1:], mse[1:], "-", color=C_FIT, lw=1.9)
    ax2.set_xlabel(L["pds_ax"], fontsize=10.5, **L["font"])
    ax2.set_ylabel(L["mse_ax"], fontsize=10.5, **L["font"])
    ax2.set_title(L["fig5b"], fontsize=11.5, pad=8, **L["font"])
    i = int(np.argmin(np.abs(PDS_SWEEP - 0.25)))
    ax2.plot([PDS_SWEEP[i]], [mse[i]], "o", color=C_FIT, ms=7, zorder=5)
    ax2.annotate(L["a_mse"].format(PDS_SWEEP[i], mse[i]), xy=(PDS_SWEEP[i], mse[i]),
                 xytext=(-16, -40), textcoords="offset points", fontsize=9.0,
                 ha="right", **L["font"],
                 arrowprops=dict(arrowstyle="->", color=C_FIT, lw=1.0))
    ax2.text(0.97, 0.06, L["zero_mse"].format(mse[0]), transform=ax2.transAxes,
             fontsize=9.0, ha="right", va="bottom", color=C_IDEAL, **L["font"])
    ax2.text(0.04, 0.96, L["fig5note"], transform=ax2.transAxes, fontsize=9.0,
             va="top", **L["font"],
             bbox=dict(boxstyle="round,pad=0.35", fc="#fffbe6", ec="#d9c77a", lw=0.8))
    ax2.grid(alpha=0.25, which="both")

    fig.suptitle(L["fig5"], fontsize=13.0, y=0.99, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(L["dir"], "fig5-model-mismatch.png"), dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# 라벨
# ---------------------------------------------------------------------------

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 공통
        "pol_ax": "편광자 눈금 $P$ (deg)",
        "delta_ax": "시편의 $\\Delta$ (deg)",
        "delta_true_ax": "시편의 참 $\\Delta$ (deg)",
        "dpsi_ax": "$\\Psi$ 오차 (deg)",
        "ddel_ax": "$\\Delta$ 오차 (deg)",
        # 그림 1
        "fig1": "세 가지 오차가 측정 신호에 하는 일",
        "fig1a": "(a) 분석기 한 바퀴 동안의 파형",
        "fig1b": "(b) 편광자를 훑을 때 $(\\alpha', \\beta')$ 가 그리는 타원",
        "ana_ax": "분석기 눈금 $A$ (deg)",
        "norm_int": "정규화 검출 세기",
        "ideal": "오차 없음",
        "e_ps": "$P_s$ = 4° (편광자 눈금)",
        "e_as": "$A_s$ = 8° (분석기 눈금)",
        "e_eta": "$\\eta$ = 0.7 (검출기 감쇠)",
        "slide": "$P_s$ 는 타원을 바꾸지 않고\n점만 그 위에서 민다",
        "fig1note": "$P_s$ 는 점을 타원 위에서 밀고,  $A_s$ 는 타원을 돌리고,\n"
                    "$\\eta$ 는 타원을 줄인다 — 셋은 서로 다른 일을 한다",
        # 그림 2
        "fig2": "잔차 보정 — 최소점의 위치와 깊이에서 세 값을 읽는다",
        "fig2a": "(a) 잔차 $R(P)$ 와 포물선 맞춤",
        "fig2b": "(b) 포물선 근사가 남기는 $\\eta$ 오차",
        "res_ax": "잔차 $R = 1 - (\\alpha'^2 + \\beta'^2)$",
        "r_meas": "합성한 측정값",
        "r_fit": "2차 맞춤",
        "a_ps": "최소점의 위치 → $P_s$ = {:.2f}°",
        "a_eta": "최솟값의 깊이 → $\\eta$ = {:.4f}  (참값 0.9687)",
        "half_ax": "훑는 범위 (± deg, 오른쪽으로 갈수록 좁다)",
        "eta_err_ax": "$\\eta$ 회수 오차",
        "fig2note": "범위를 좁히면 오차가 함께 줄어든다.\n원인은 잡음이 아니라 포물선 근사다.",
        # 그림 3
        "fig3": "$\\Delta$ 가 0·180°에 가까우면 잔차 보정은 최소점을 잃는다",
        "fig3a": "(a) 최소점의 깊이 $R(P) - R_{\\min}$",
        "fig3b": "(b) 잡음이 있을 때 $P_s$ 회수 오차",
        "depth_ax": "최소점으로부터의 깊이",
        "ps_err_ax": "$P_s$ 회수 오차 중앙값 (deg)",
        "noise_band": "잡음이 $R$ 에 남기는 흔들림 ({:.1e})",
        "fail": "{}/{} 회 맞춤 실패",
        "target": "목표 0.01°",
        # 그림 4
        "fig4": "보정하지 않은 눈금값을 그대로 넣으면 — 세 오차는 서로 다른 곳으로 번진다",
        "fig4a": "(a) $\\Psi$ 에 번지는 오차",
        "fig4b": "(b) $\\Delta$ 에 번지는 오차",
        "e_ps4": "$P_s$ = 0.1° 만 있을 때",
        "e_as4": "$A_s$ = 0.1° 만 있을 때",
        "e_eta4": "$\\eta$ = 0.999 만 있을 때",
        "fig4note": "$\\eta$ 의 오차는 $\\Delta$ 가 0·180° 로 갈수록\n훨씬 크게 번진다",
        "zero_note": "$P_s$ 는 $\\Delta$ 에 전혀 번지지 않는다 (정확히 0)",
        # 그림 5
        "fig5": "모형에 없는 효과는 $\\eta$ 로 숨고, MSE 가 그것을 알린다",
        "fig5a": "(a) 넣지 않은 편광 의존 감도가 밀어내는 값",
        "fig5b": "(b) 같은 맞춤의 MSE",
        "pds_ax": "실제로 있는 편광 의존 감도 크기",
        "eta_fit_ax": "맞춰 나온 $\\eta$",
        "eta_fit": "맞춰 나온 $\\eta$",
        "eta_true": "참값 1",
        "psi_shift": "$\\Psi$ 오차 (오른쪽 축)",
        "mse_ax": "MSE",
        "a_mse": "크기 {:.2f} 에서\nMSE = {:.1e}",
        "zero_mse": "편광 의존 감도가 아예 없으면 MSE = {:.0e}",
        "fig5note": "맞춤은 수렴한다. 수렴했다는 것이\n맞았다는 뜻은 아니다.",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "pol_ax": "polarizer reading $P$ (deg)",
        "delta_ax": "sample $\\Delta$ (deg)",
        "delta_true_ax": "true sample $\\Delta$ (deg)",
        "dpsi_ax": "error in $\\Psi$ (deg)",
        "ddel_ax": "error in $\\Delta$ (deg)",
        "fig1": "What the three errors do to the measured signal",
        "fig1a": "(a) waveform over one analyzer revolution",
        "fig1b": "(b) ellipse traced by $(\\alpha', \\beta')$ as the polarizer is scanned",
        "ana_ax": "analyzer reading $A$ (deg)",
        "norm_int": "normalized intensity",
        "ideal": "no error",
        "e_ps": "$P_s$ = 4° (polarizer offset)",
        "e_as": "$A_s$ = 8° (analyzer offset)",
        "e_eta": "$\\eta$ = 0.7 (detector gain)",
        "slide": "$P_s$ slides the point,\nleaving the ellipse alone",
        "fig1note": "$P_s$ slides the point along the ellipse, $A_s$ rotates it,\n"
                    "$\\eta$ shrinks it \u2014 the three act differently",
        "fig2": "Residual calibration: three numbers from the position and depth of one minimum",
        "fig2a": "(a) residual $R(P)$ and its parabolic fit",
        "fig2b": "(b) error in $\\eta$ left by the parabolic approximation",
        "res_ax": "residual $R = 1 - (\\alpha'^2 + \\beta'^2)$",
        "r_meas": "synthesized data",
        "r_fit": "quadratic fit",
        "a_ps": "position of minimum $\\rightarrow$ $P_s$ = {:.2f}°",
        "a_eta": "depth of minimum $\\rightarrow$ $\\eta$ = {:.4f}  (true 0.9687)",
        "half_ax": "scan range (± deg, narrowing to the right)",
        "eta_err_ax": "error in recovered $\\eta$",
        "fig2note": "Narrowing the scan shrinks the error \u2014 the\ncause is the approximation, not noise.",
        "fig3": "Near $\\Delta$ = 0° or 180°, residual calibration loses its minimum",
        "fig3a": "(a) depth of the minimum, $R(P) - R_{\\min}$",
        "fig3b": "(b) error in recovered $P_s$ with noise present",
        "depth_ax": "depth above the minimum",
        "ps_err_ax": "median error in $P_s$ (deg)",
        "noise_band": "scatter noise leaves in $R$ ({:.1e})",
        "fail": "{}/{} fits failed",
        "target": "target 0.01°",
        "fig4": "Uncalibrated readings fed straight into the ideal formulas: each error goes its own way",
        "fig4a": "(a) error propagated into $\\Psi$",
        "fig4b": "(b) error propagated into $\\Delta$",
        "e_ps4": "$P_s$ = 0.1° alone",
        "e_as4": "$A_s$ = 0.1° alone",
        "e_eta4": "$\\eta$ = 0.999 alone",
        "fig4note": "An error in $\\eta$ propagates far more strongly\nas $\\Delta$ approaches 0° or 180°",
        "zero_note": "$P_s$ does not propagate into $\\Delta$ at all (exactly zero)",
        "fig5": "An unmodeled effect hides in $\\eta$, and the MSE reports it",
        "fig5a": "(a) values pushed by a polarization-dependent sensitivity left out of the model",
        "fig5b": "(b) MSE of the same fits",
        "pds_ax": "magnitude of the polarization-dependent sensitivity present",
        "eta_fit_ax": "fitted $\\eta$",
        "eta_fit": "fitted $\\eta$",
        "eta_true": "true value 1",
        "psi_shift": "error in $\\Psi$ (right axis)",
        "mse_ax": "MSE",
        "a_mse": "at magnitude {:.2f},\nMSE = {:.1e}",
        "zero_mse": "with no such sensitivity at all, MSE = {:.0e}",
        "fig5note": "The fit converges. Converging is not\nthe same as being right.",
    },
}


def main():
    conv = cached("conv", eta_convergence)
    print(f"포물선 근사 eta 오차: ±{conv[0,0]:.1f}° 에서 {conv[0,2]:.2e} -> "
          f"±{conv[-1,0]:.1f}° 에서 {conv[-1,2]:.2e} ({conv[0,2]/conv[-1,2]:.0f} 배)")
    _, rows, sigma_R = cached("blind", blind_spot_data)
    print(f"맹점: Delta 90° 에서 Ps 오차 {rows[0,1]:.4f}°, 1° 에서 {rows[-1,1]:.4f}° "
          f"({rows[-1,1]/rows[0,1]:.0f} 배), 잡음이 R 에 남기는 흔들림 {sigma_R:.2e}")
    data = cached("uncal", uncalibrated_error)
    for key, _, _, _ in CASES:
        dpsi, ddel = data[key]
        print(f"  {key}: Delta=90° 에서 dPsi {np.abs(dpsi[87]):.4f}° "
              f"dDelta {np.abs(ddel[87]):.4f}° / "
              f"Delta=3° 에서 dDelta {np.abs(ddel[0]):.4f}°")
    eta, psi, mse = cached("mismatch", model_mismatch)
    i = int(np.argmin(np.abs(PDS_SWEEP - 0.25)))
    print(f"모형 불일치: 크기 0.25 에서 eta {eta[i]:.4f}, Psi {psi[i]:.3f}°, "
          f"MSE {mse[i]:.2e} (크기 0 에서 MSE {mse[0]:.1e})")
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        figure1(L); figure2(L); figure3(L); figure4(L); figure5(L)
        print(f"  [{lang}] 그림 5개 저장 -> {os.path.normpath(L['dir'])}")


if __name__ == "__main__":
    main()

"""XPS 기초 4편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py      (그림 4의 스윕 때문에 2분 남짓 걸린다)
출력: ../../assets/img/posts/xps-quantification-peak-fitting/      (한국어)
      ../../assets/img/posts/xps-quantification-peak-fitting/en/   (영문)

모든 스펙트럼은 합성 데이터다. 계산은 xps_fit.py에 있고 여기서는 그리기만 한다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import voigt_profile

import xps_fit as X

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "xps-quantification-peak-fitting",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "be": "결합에너지 (eV)", "cps": "계수", "x": "에너지 (eV)", "int": "세기 (면적 1로 정규화)",
        "data": "데이터", "true_bg": "참 배경", "bg": {"linear": "선형", "Shirley": "Shirley", "Tougaard": "Tougaard"},
        "area_err": "면적 오차", "fig1": "배경 모델에 따른 피크 면적 (합성 Si 2p, 참 배경은 Tougaard형)",
        "gauss": "Gaussian", "lor": "Lorentzian", "voigt": "Voigt (G ⊗ L)",
        "fig2a": "(a) 선형 축", "fig2b": "(b) 로그 축: 꼬리의 차이",
        "fig2": "Gaussian, Lorentzian, 그리고 합성곱 Voigt",
        "fit": "피팅 합", "shirley": "Shirley 배경", "resid": "정규화 잔차",
        "si0": "Si$^0$ 2p$_{3/2}$ / 2p$_{1/2}$", "si4": "Si$^{4+}$ 2p$_{3/2}$ / 2p$_{1/2}$",
        "truth": "참값", "est": "추정",
        "area": "면적", "frac": "Si$^{4+}$ 몫",
        "fig3": "제약 피팅 결과 (합성 Si 2p, Shirley 배경)",
        "ncomp": "제약 없는 단일선 개수 n", "chi2": r"$\chi^2$ (잡음 실현 20회 평균)",
        "constrained": "제약 모델(이중선 2개)", "fake": "가짜 중간 산화물(Si$^{1+}$–Si$^{3+}$) 몫",
        "ox": "Si$^{4+}$ 몫", "share": "보고되는 면적 몫 (%)",
        "true_ox": "참 Si$^{4+}$ 몫", "true_fake": "참 중간 산화물 = 0",
        "fig4a": r"(a) 성분을 늘리면 $\chi^2$는 계속 준다", "fig4b": "(b) 그런데 없는 화학종이 보고된다",
        "fig4": "성분 수를 늘릴 때 (합성 Si 2p)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "be": "Binding energy (eV)", "cps": "Counts", "x": "Energy (eV)",
        "int": "Intensity (unit area)",
        "data": "data", "true_bg": "true background",
        "bg": {"linear": "linear", "Shirley": "Shirley", "Tougaard": "Tougaard"},
        "area_err": "area error",
        "fig1": "Peak area by background model (synthetic Si 2p, Tougaard-type true background)",
        "gauss": "Gaussian", "lor": "Lorentzian", "voigt": "Voigt (G ⊗ L)",
        "fig2a": "(a) Linear axis", "fig2b": "(b) Log axis: difference in the tails",
        "fig2": "Gaussian, Lorentzian and their convolution, the Voigt",
        "fit": "fit sum", "shirley": "Shirley background", "resid": "normalized residual",
        "si0": "Si$^0$ 2p$_{3/2}$ / 2p$_{1/2}$", "si4": "Si$^{4+}$ 2p$_{3/2}$ / 2p$_{1/2}$",
        "truth": "true", "est": "fit",
        "area": "area", "frac": "Si$^{4+}$ share",
        "fig3": "Constrained fit (synthetic Si 2p, Shirley background)",
        "ncomp": "Number of unconstrained singlets n",
        "chi2": r"$\chi^2$ (mean of 20 noise realizations)",
        "constrained": "constrained model (two doublets)",
        "fake": "spurious intermediate oxide (Si$^{1+}$–Si$^{3+}$)",
        "ox": "Si$^{4+}$ share", "share": "Reported area share (%)",
        "true_ox": "true Si$^{4+}$ share", "true_fake": "true intermediate oxide = 0",
        "fig4a": r"(a) More components, lower $\chi^2$", "fig4b": "(b) But nonexistent species are reported",
        "fig4": "Adding components (synthetic Si 2p)",
    },
}

# ---------------------------------------------------------------------------
# 계산 (한 번만)
# ---------------------------------------------------------------------------
E = X.E
y1, bgs, areas, a_true = X.background_comparison()
bg_true = X.true_background()

xg = np.linspace(-6, 6, 2401)
FW = 1.0
g = voigt_profile(xg, FW * X.FWHM_TO_SIGMA, 0.0)
lz = voigt_profile(xg, 1e-9, FW / 2)
vo = voigt_profile(xg, FW * X.FWHM_TO_SIGMA, FW / 2)


def fwhm(x, y):
    half = y.max() / 2
    idx = np.where(y >= half)[0]
    return x[idx[-1]] - x[idx[0]]


fw_v = fwhm(xg, vo)

y3, bg3, n_iter, p3, e3, c2r3, _ = X.constrained_demo()
resid3 = (y3 - bg3 - X.doublet(E, *p3[:3]) - X.doublet(E, *p3[3:])) / np.sqrt(np.clip(y3, 1, None))

print("스윕 계산 중 (2분 남짓)...")
sweep_out, sweep_con = X.sweep()
ns = sorted(sweep_out)


def comps(p):
    """이중선을 2p3/2, 2p1/2 두 성분으로 나눠 돌려준다."""
    c, a, w = p
    a32 = a / (1 + X.SO_RATIO)
    return X.voigt(E, c, a32, w), X.voigt(E, c + X.SO_SPLIT, a32 * X.SO_RATIO, w)


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L_):
    out, F, LF = L_["dir"], L_["font"], L_["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. 배경 모델 비교
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.plot(E, y1, ".", ms=2, color="gray", label=L_["data"])
    ax.plot(E, bg_true, color="black", lw=1.2, ls=":", label=L_["true_bg"])
    colors = {"linear": "tab:green", "Shirley": "tab:blue", "Tougaard": "tab:red"}
    for k, b in bgs.items():
        err = 100 * (areas[k] / a_true - 1)
        ax.plot(E, b, color=colors[k], lw=1.6,
                label=f"{L_['bg'][k]}  ({L_['area_err']} {err:+.1f}%)")
    ax.set_xlim(E[-1], E[0])
    ax.set_ylim(0, 300)                 # 배경을 보이려고 세로축을 자른다
    ax.set_xlabel(L_["be"], fontsize=11, **F)
    ax.set_ylabel(L_["cps"], fontsize=11, **F)
    ax.legend(prop={**LF, "size": 9}, loc="upper left")
    ax.set_title(L_["fig1"], fontsize=11.5, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-backgrounds.png"), dpi=150)
    plt.close(fig)

    # 그림 2. Gaussian / Lorentzian / Voigt
    fig2, (a1, a2) = plt.subplots(1, 2, figsize=(10.2, 4.0))
    for ax in (a1, a2):
        ax.plot(xg, g, color="tab:blue", lw=1.8, label=f"{L_['gauss']} (FWHM {FW:.2f})")
        ax.plot(xg, lz, color="tab:orange", lw=1.8, label=f"{L_['lor']} (FWHM {FW:.2f})")
        ax.plot(xg, vo, color="black", lw=1.8, ls="--", label=f"{L_['voigt']} (FWHM {fw_v:.2f})")
        ax.set_xlabel(L_["x"], fontsize=11, **F)
        ax.grid(alpha=0.3)
    for y, col, w in ((g, "tab:blue", FW), (vo, "black", fw_v)):
        a1.annotate("", xy=(w / 2, y.max() / 2), xytext=(-w / 2, y.max() / 2),
                    arrowprops=dict(arrowstyle="<->", color=col, lw=1))
    a1.set_xlim(-4, 4)
    a1.set_ylabel(L_["int"], fontsize=11, **F)
    a1.legend(prop={**LF, "size": 8.5}, loc="upper left")
    a1.set_title(L_["fig2a"], fontsize=11, **F)
    a2.set_yscale("log")
    a2.set_ylim(1e-4, 2)
    a2.set_xlim(-6, 6)
    a2.set_title(L_["fig2b"], fontsize=11, **F)
    fig2.suptitle(L_["fig2"], fontsize=12, **F)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "fig2-voigt-composition.png"), dpi=150)
    plt.close(fig2)

    # 그림 3. 제약 피팅 결과
    fig3, (t, b) = plt.subplots(2, 1, figsize=(8.4, 6.0), sharex=True,
                                gridspec_kw={"height_ratios": [3, 1]})
    t.plot(E, y3, ".", ms=2.5, color="gray", label=L_["data"])
    t.plot(E, bg3, color="tab:green", lw=1.2, label=L_["shirley"])
    for p, col, lab in ((p3[:3], "tab:blue", L_["si0"]), (p3[3:], "tab:red", L_["si4"])):
        c32, c12 = comps(p)
        t.fill_between(E, bg3, bg3 + c32, color=col, alpha=0.35, label=lab)
        t.fill_between(E, bg3, bg3 + c12, color=col, alpha=0.18)
    t.plot(E, bg3 + X.doublet(E, *p3[:3]) + X.doublet(E, *p3[3:]), color="black", lw=1.2,
           label=L_["fit"])
    tr, fr = X.TRUE, p3
    rows = [
        [f"Si$^0$ {L_['area']}", f"{tr['Si0'][1]:.0f}", f"{fr[1]:.0f} ± {e3[1]:.0f}"],
        [f"Si$^{{4+}}$ {L_['area']}", f"{tr['Si4+'][1]:.0f}", f"{fr[4]:.0f} ± {e3[4]:.0f}"],
        ["Si$^{4+}$ FWHM", f"{tr['Si4+'][2]:.2f}", f"{fr[5]:.2f} ± {e3[5]:.2f}"],
        [L_["frac"], f"{100 * X.TRUE_OX_FRACTION:.1f}%", f"{100 * fr[4] / (fr[1] + fr[4]):.1f}%"],
        [r"$\chi^2_{red}$", "1", f"{c2r3:.2f}"],
    ]
    tab = t.table(cellText=rows, colLabels=["", L_["truth"], L_["est"]], cellLoc="center",
                  bbox=[0.02, 0.50, 0.40, 0.46])
    tab.auto_set_font_size(False)
    for cell in tab.get_celld().values():
        cell.set_fontsize(8.5)
        cell.set_edgecolor("0.8")
        if F:
            cell.get_text().set_fontfamily(F["fontfamily"])
    t.set_ylabel(L_["cps"], fontsize=11, **F)
    t.legend(prop={**LF, "size": 8.5}, loc="upper right")
    t.set_title(L_["fig3"], fontsize=11.5, **F)
    b.plot(E, resid3, color="tab:purple", lw=0.8)
    b.axhline(0, color="black", lw=0.6)
    b.set_ylabel(L_["resid"], fontsize=9.5, **F)
    b.set_xlim(E[-1], E[0])
    b.set_xlabel(L_["be"], fontsize=11, **F)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, "fig3-fit-result.png"), dpi=150)
    plt.close(fig3)

    # 그림 4. 과적합
    fig4, (c1, c2) = plt.subplots(1, 2, figsize=(10.4, 4.3))
    chi = np.array([np.mean(sweep_out[n]["chi2"]) for n in ns])
    chis = np.array([np.std(sweep_out[n]["chi2"]) for n in ns])
    c1.errorbar(ns, chi, yerr=chis, fmt="o-", color="black", capsize=3)
    c1.axhline(np.mean(sweep_con["chi2"]), color="tab:blue", ls="--", lw=1.2, label=L_["constrained"])
    c1.set_xlabel(L_["ncomp"], fontsize=10.5, **F)
    c1.set_ylabel(L_["chi2"], fontsize=10.5, **F)
    c1.set_title(L_["fig4a"], fontsize=11, **F)
    c1.legend(prop={**LF, "size": 9})
    c1.grid(alpha=0.3)
    fake = np.array([100 * np.mean(sweep_out[n]["fake"]) for n in ns])
    fakes = np.array([100 * np.std(sweep_out[n]["fake"]) for n in ns])
    ox = np.array([100 * np.mean(sweep_out[n]["ox"]) for n in ns])
    oxs = np.array([100 * np.std(sweep_out[n]["ox"]) for n in ns])
    c2.errorbar(ns, fake, yerr=fakes, fmt="o-", color="tab:red", capsize=3, label=L_["fake"])
    c2.errorbar(ns, ox, yerr=oxs, fmt="s-", color="tab:blue", capsize=3, label=L_["ox"])
    c2.axhline(100 * X.TRUE_OX_FRACTION, color="tab:blue", ls=":", lw=1, label=L_["true_ox"])
    c2.axhline(0, color="tab:red", ls=":", lw=1, label=L_["true_fake"])
    c2.set_ylim(-3, 66)
    c2.set_xlabel(L_["ncomp"], fontsize=10.5, **F)
    c2.set_ylabel(L_["share"], fontsize=10.5, **F)
    c2.set_title(L_["fig4b"], fontsize=11, **F)
    c2.legend(prop={**LF, "size": 8.5}, loc="upper left")
    c2.grid(alpha=0.3)
    fig4.suptitle(L_["fig4"], fontsize=12, **F)
    fig4.tight_layout()
    fig4.savefig(os.path.join(out, "fig4-overfitting.png"), dpi=150)
    plt.close(fig4)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

print(f"\nVoigt FWHM (G {FW}, L {FW}) = {fw_v:.3f} eV")
print("배경 면적 오차:", {k: f"{100 * (v / a_true - 1):+.1f}%" for k, v in areas.items()})

"""XPS 기초 3편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/xps-spectrum-structure-satellites/      (한국어)
      ../../assets/img/posts/xps-spectrum-structure-satellites/en/   (영문)

그림 1~3의 스펙트럼은 전부 합성(모사)이며 실측이 아니다. 피크 위치는 아래 출처의
문헌값을 따랐다.
  Ti 2p (TiO2): Thermo Fisher XPS Periodic Table — 458.5 eV, 분리 5.7 eV
  Cu 2p3/2: M. C. Biesinger et al., Appl. Surf. Sci. 257, 887 (2010) 표 8
그림 4는 같은 논문 표 7·10의 문헌 평균값(결합에너지, 수정 오제 파라미터)을 그대로 찍었다.
Doniach–Šunjić 식은 CasaXPS "Peak Fitting in XPS"(2006)의 형태를 따랐다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import voigt_profile

# 한글 텍스트에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면 눈금의 마이너스
# 기호가 한글 폰트에 없어 깨진다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "xps-spectrum-structure-satellites",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "be": "결합에너지 (eV)", "rel_be": "상대 결합에너지 (eV)", "cps": "계수율 (임의 단위)",
        "fig1a": "(a) TiO$_2$형 Ti 2p: p 준위, 면적비 1:2", "fig1b": "(b) d 준위 이중선: 면적비 2:3",
        "sum": "합", "area": "설정 면적비", "window": "창 안 적분", "split": "분리",
        "fig1": "스핀-궤도 이중선 (합성 스펙트럼)",
        "fig2a": "(a) Cu(0)형: 위성 없음", "fig2b": "(b) Cu(II)형: 셰이크업 위성",
        "main": "주 피크", "sat": "셰이크업 위성",
        "fig2": "Cu 2p$_{3/2}$ 영역의 위성 유무 (합성 스펙트럼)",
        "fig3a": "(a) 같은 위치·같은 면적", "fig3b": "(b) 대칭 Voigt 하나로 피팅",
        "fig3c": "(c) 대칭 Voigt 둘로 피팅",
        "ds": "비대칭 (DS ⊗ Gauss)", "voigt": "대칭 (Voigt)", "data": "데이터",
        "fit": "피팅", "resid": "잔차 (×3, 아래로 이동)", "fake": "가짜 성분",
        "fig3": "금속의 비대칭 선형을 대칭 함수로 맞추면 (합성 스펙트럼)",
        "wag_x": r"광전자 결합에너지 $E_B$ (eV)", "wag_y": r"오제 운동에너지 $E_K$ (eV)",
        "fig4a": "(a) Cu: 2p$_{3/2}$ 와 L$_3$M$_{45}$M$_{45}$",
        "fig4b": "(b) Zn: 2p$_{3/2}$ 와 L$_3$M$_{45}$M$_{45}$",
        "names": {"Cu0": "Cu 금속", "Cu2O": "Cu$_2$O", "CuO": "CuO", "Zn0": "Zn 금속", "ZnO": "ZnO"},
        "fig4": "Wagner 도표 (Biesinger et al. 2010의 문헌 평균값)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "be": "Binding energy (eV)", "rel_be": "Relative binding energy (eV)",
        "cps": "Count rate (arb. units)",
        "fig1a": "(a) TiO$_2$-like Ti 2p: p level, area ratio 1:2",
        "fig1b": "(b) d-level doublet: area ratio 2:3",
        "sum": "sum", "area": "set area ratio", "window": "window:", "split": "splitting",
        "fig1": "Spin-orbit doublets (synthetic spectra)",
        "fig2a": "(a) Cu(0)-like: no satellite", "fig2b": "(b) Cu(II)-like: shake-up satellite",
        "main": "main peak", "sat": "shake-up satellite",
        "fig2": "Presence of satellites in the Cu 2p$_{3/2}$ region (synthetic spectra)",
        "fig3a": "(a) Same position, same area", "fig3b": "(b) Fit with one symmetric Voigt",
        "fig3c": "(c) Fit with two symmetric Voigts",
        "ds": "asymmetric (DS ⊗ Gauss)", "voigt": "symmetric (Voigt)", "data": "data",
        "fit": "fit", "resid": "residual (×3, offset)", "fake": "spurious component",
        "fig3": "Fitting a metal's asymmetric line with symmetric functions (synthetic spectra)",
        "wag_x": r"Photoelectron binding energy $E_B$ (eV)",
        "wag_y": r"Auger kinetic energy $E_K$ (eV)",
        "fig4a": "(a) Cu: 2p$_{3/2}$ and L$_3$M$_{45}$M$_{45}$",
        "fig4b": "(b) Zn: 2p$_{3/2}$ and L$_3$M$_{45}$M$_{45}$",
        "names": {"Cu0": "Cu metal", "Cu2O": "Cu$_2$O", "CuO": "CuO", "Zn0": "Zn metal", "ZnO": "ZnO"},
        "fig4": "Wagner plot (literature averages from Biesinger et al. 2010)",
    },
}

FWHM_TO_SIGMA = 1 / (2 * np.sqrt(2 * np.log(2)))
rng = np.random.default_rng(0)


def voigt(x, center, area, sigma, gamma):
    """면적이 area인 대칭 Voigt. sigma: Gaussian 표준편차, gamma: Lorentzian 반폭(HWHM)."""
    return area * voigt_profile(x - center, sigma, gamma)


def doublet(x, e_main, split, ratio, sigma, gamma_main, gamma_minor, area_main=1.0):
    """스핀-궤도 이중선. e_main: 큰 쪽(j = l+1/2) 성분의 위치,
    split: 분리 폭(작은 성분이 높은 결합에너지 쪽), ratio: 작은/큰 면적비."""
    major = voigt(x, e_main, area_main, sigma, gamma_main)
    minor = voigt(x, e_main + split, area_main * ratio, sigma, gamma_minor)
    return major, minor


def doniach_sunjic(be, e0, alpha, f):
    """Doniach–Šunjić 선형 (Γ(1-α) 상수 인자 생략). be: 결합에너지 축.
    원래 식은 운동에너지 x에 대해 쓰여 꼬리가 x가 작은 쪽으로 늘어진다.
    결합에너지 축에서 꼬리가 높은 쪽에 오도록 u = -(be - e0)로 뒤집어 넣는다."""
    u = -(be - e0)
    return (np.cos(np.pi * alpha / 2 + (1 - alpha) * np.arctan(u / f))
            / (f**2 + u**2) ** ((1 - alpha) / 2))


def gauss_convolve(x, y, fwhm):
    """등간격 격자 위에서 Gaussian(FWHM)과 합성곱."""
    dx = x[1] - x[0]
    k = np.arange(-5 * fwhm, 5 * fwhm + dx, dx)
    g = np.exp(-0.5 * (k / (fwhm * FWHM_TO_SIGMA)) ** 2)
    return np.convolve(y, g / g.sum(), mode="same")


def area(x, y):
    return float(np.sum(y) * (x[1] - x[0]))


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------
# 그림1 (a) TiO2형 Ti 2p: 2p3/2 458.5 eV, 분리 5.7 eV, 면적비 1:2.
#          2p1/2는 Coster–Kronig 과정으로 수명이 짧아 Lorentzian 폭을 더 크게 준다.
x_ti = np.linspace(452, 470, 1800)
ti_major, ti_minor = doublet(x_ti, 458.5, 5.7, 0.5, sigma=0.9 * FWHM_TO_SIGMA,
                             gamma_main=0.15, gamma_minor=0.6)
# 그림1 (b) d 준위: j=5/2 : j=3/2 = 3 : 2, 분리 폭은 임의(6 eV)
x_d = np.linspace(-4, 10, 1400)
d_major, d_minor = doublet(x_d, 0.0, 6.0, 2 / 3, sigma=0.8 * FWHM_TO_SIGMA,
                           gamma_main=0.2, gamma_minor=0.2)
ratio_ti = area(x_ti, ti_minor) / area(x_ti, ti_major)
ratio_d = area(x_d, d_minor) / area(x_d, d_major)

# 그림2 Cu 2p3/2 영역. Cu(0) 932.63 eV(좁음), Cu(II) 산화물 933.76 eV(FWHM 3.0 eV, Biesinger 표 8).
# 위성 묶음은 940~945 eV 부근(Thermo Fisher: ~943 eV). 위성의 모양·세기는 모사다.
x_cu = np.linspace(926, 950, 2400)
cu0_main = voigt(x_cu, 932.63, 1.0, 0.7 * FWHM_TO_SIGMA, 0.15)
cu2_main = voigt(x_cu, 933.76, 1.0, 2.6 * FWHM_TO_SIGMA, 0.3)
cu2_sat = (voigt(x_cu, 941.2, 0.22, 1.4 * FWHM_TO_SIGMA, 0.4)
           + voigt(x_cu, 943.4, 0.30, 1.6 * FWHM_TO_SIGMA, 0.4))

# 그림3 금속 피크: DS(α = 0.12, f = 0.15 eV) ⊗ Gauss(FWHM 0.8 eV), 위치 0, 면적 1
x_m = np.linspace(-6, 12, 3601)
DS_ALPHA, DS_F, G_FWHM = 0.12, 0.15, 0.8
ds_raw = gauss_convolve(x_m, doniach_sunjic(x_m, 0.0, DS_ALPHA, DS_F), G_FWHM)
metal = ds_raw / area(x_m, ds_raw)
sym_same = voigt(x_m, 0.0, 1.0, G_FWHM * FWHM_TO_SIGMA, DS_F)       # 같은 위치·면적의 대칭형
tail_high = area(x_m[x_m > 0], metal[x_m > 0])                     # 높은 BE 쪽 몫
peak_pos = x_m[np.argmax(metal)]
noise_level = 0.004
metal_noisy = metal + rng.normal(0, noise_level, x_m.size)


def one_voigt(x, c, a, s, g):
    return voigt(x, c, a, abs(s), abs(g))


def two_voigt(x, c1, a1, s1, g1, c2, a2, s2, g2):
    return voigt(x, c1, a1, abs(s1), abs(g1)) + voigt(x, c2, a2, abs(s2), abs(g2))


p1, _ = curve_fit(one_voigt, x_m, metal_noisy, p0=[0.2, 1.0, 0.4, 0.2])
fit1 = one_voigt(x_m, *p1)
p2, _ = curve_fit(two_voigt, x_m, metal_noisy,
                  p0=[0.0, 0.8, 0.35, 0.1, 1.0, 0.2, 0.6, 0.3], maxfev=20000)
fit2 = two_voigt(x_m, *p2)
comp_a = voigt(x_m, p2[0], p2[1], abs(p2[2]), abs(p2[3]))
comp_b = voigt(x_m, p2[4], p2[5], abs(p2[6]), abs(p2[7]))
if p2[0] > p2[4]:                                                  # a가 주 성분이 되게 정렬
    comp_a, comp_b, p2 = comp_b, comp_a, np.r_[p2[4:], p2[:4]]
rms1 = float(np.sqrt(np.mean((metal_noisy - fit1) ** 2)))
rms2 = float(np.sqrt(np.mean((metal_noisy - fit2) ** 2)))
fake_share = abs(p2[5]) / (abs(p2[1]) + abs(p2[5]))

# 그림4 Wagner 도표: Biesinger 2010 표 7(Cu)·표 10(Zn)의 문헌 평균 (E_B, α')
WAGNER = {
    "Cu": {"Cu0": (932.61, 1851.23), "Cu2O": (932.43, 1849.19), "CuO": (933.57, 1851.49)},
    "Zn": {"Zn0": (1021.62, 2013.85), "ZnO": (1021.96, 2010.14)},
}


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L_):
    """주어진 라벨 묶음으로 그림 4개를 그린다."""
    out, F, LF = L_["dir"], L_["font"], L_["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. 스핀-궤도 이중선
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.4, 4.0))
    for ax, x, maj, mnr, lab_maj, lab_mnr, r, r_set, title, split in (
        (a1, x_ti, ti_major, ti_minor, r"2p$_{3/2}$", r"2p$_{1/2}$", ratio_ti, 0.5, L_["fig1a"], 5.7),
        (a2, x_d, d_major, d_minor, r"d$_{5/2}$", r"d$_{3/2}$", ratio_d, 2 / 3, L_["fig1b"], 6.0),
    ):
        ax.fill_between(x, maj, color="tab:blue", alpha=0.3, label=lab_maj)
        ax.fill_between(x, mnr, color="tab:orange", alpha=0.3, label=lab_mnr)
        ax.plot(x, maj + mnr, color="black", lw=1.4, label=L_["sum"])
        xm, xn = x[np.argmax(maj)], x[np.argmax(mnr)]
        ytop = max(maj.max(), mnr.max()) * 1.12
        ax.annotate("", xy=(xm, ytop), xytext=(xn, ytop),
                    arrowprops=dict(arrowstyle="<->", color="gray", lw=1))
        ax.text((xm + xn) / 2, ytop * 1.03, f"{L_['split']} {split} eV", ha="center",
                fontsize=9, color="gray", **F)
        ax.text(0.03, 0.95, f"{L_['area']} {r_set:.2f}\n({L_['window']} {r:.2f})", transform=ax.transAxes,
                fontsize=9.5, va="top", **F)
        ax.set_xlim(x[-1], x[0])
        ax.set_ylim(0, ytop * 1.25)
        ax.set_yticks([])
        ax.set_title(title, fontsize=10.5, **F)
        ax.legend(prop={**LF, "size": 9}, loc="upper right")
    a1.set_xlabel(L_["be"], fontsize=11, **F)
    a2.set_xlabel(L_["rel_be"], fontsize=11, **F)
    a1.set_ylabel(L_["cps"], fontsize=11, **F)
    fig.suptitle(L_["fig1"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-spin-orbit-doublet.png"), dpi=150)
    plt.close(fig)

    # 그림 2. 위성 유무
    fig2, (b1, b2) = plt.subplots(1, 2, figsize=(10.4, 3.9), sharey=True)
    b1.plot(x_cu, cu0_main, color="black", lw=1.5)
    b1.fill_between(x_cu, cu0_main, color="tab:blue", alpha=0.3, label=L_["main"])
    b2.plot(x_cu, cu2_main + cu2_sat, color="black", lw=1.5)
    b2.fill_between(x_cu, cu2_main, color="tab:blue", alpha=0.3, label=L_["main"])
    b2.fill_between(x_cu, cu2_sat, color="tab:red", alpha=0.35, label=L_["sat"])
    for ax, title in ((b1, L_["fig2a"]), (b2, L_["fig2b"])):
        ax.set_xlim(x_cu[-1], x_cu[0])
        ax.set_xlabel(L_["be"], fontsize=11, **F)
        ax.set_yticks([])
        ax.set_title(title, fontsize=10.5, **F)
        ax.legend(prop={**LF, "size": 9}, loc="upper left")
    b1.set_ylabel(L_["cps"], fontsize=11, **F)
    fig2.suptitle(L_["fig2"], fontsize=12, **F)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "fig2-satellites.png"), dpi=150)
    plt.close(fig2)

    # 그림 3. 비대칭 선형과 대칭 피팅
    fig3, (c1, c2, c3) = plt.subplots(1, 3, figsize=(12.6, 4.2), sharey=True)
    c1.plot(x_m, metal, color="tab:red", lw=2, label=L_["ds"])
    c1.plot(x_m, sym_same, color="tab:blue", lw=1.6, ls="--", label=L_["voigt"])
    off = -0.12
    for ax, fit, title in ((c2, fit1, L_["fig3b"]), (c3, fit2, L_["fig3c"])):
        ax.plot(x_m, metal_noisy, ".", ms=1.5, color="gray", label=L_["data"])
        ax.plot(x_m, fit, color="black", lw=1.4, label=L_["fit"])
        ax.plot(x_m, off + 3 * (metal_noisy - fit), color="tab:purple", lw=0.8,
                label=L_["resid"])
        ax.axhline(off, color="tab:purple", lw=0.5, ls=":")
        ax.set_title(title, fontsize=10.5, **F)
    c3.fill_between(x_m, comp_a, color="tab:blue", alpha=0.25)
    c3.fill_between(x_m, comp_b, color="tab:orange", alpha=0.45, label=L_["fake"])
    c3.annotate(f"{L_['fake']}\n{100 * fake_share:.0f}%", xy=(p2[4], comp_b.max()),
                xytext=(p2[4] + 3.5, comp_b.max() + 0.25), fontsize=9, color="tab:orange",
                arrowprops=dict(arrowstyle="->", color="tab:orange"), **F)
    c2.text(0.97, 0.55, f"RMS {rms1:.4f}", transform=c2.transAxes, ha="right", fontsize=9)
    c3.text(0.97, 0.55, f"RMS {rms2:.4f}", transform=c3.transAxes, ha="right", fontsize=9)
    for ax in (c1, c2, c3):
        ax.set_xlim(6, -3)
        ax.set_xlabel(L_["rel_be"], fontsize=11, **F)
        ax.legend(prop={**LF, "size": 8.5}, loc="upper left")
    c1.set_ylim(-0.3, max(metal.max(), sym_same.max()) * 1.3)
    c1.set_yticks([])
    c1.set_ylabel(L_["cps"], fontsize=11, **F)
    c1.set_title(L_["fig3a"], fontsize=10.5, **F)
    fig3.suptitle(L_["fig3"], fontsize=12, **F)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, "fig3-asymmetric-lineshape.png"), dpi=150)
    plt.close(fig3)

    # 그림 4. Wagner 도표 (결합에너지 축은 왼쪽으로 증가, α' 등고선은 화면상 기울기 +1)
    fig4, axes = plt.subplots(1, 2, figsize=(10.4, 5.4))
    SPAN = 6.0                                     # 두 축의 범위를 같게 해 α' 선이 45°로 보이게 한다
    for ax, elem, title, pad in ((axes[0], "Cu", L_["fig4a"], 1.2), (axes[1], "Zn", L_["fig4b"], 1.2)):
        pts = WAGNER[elem]
        bes = [b for b, _ in pts.values()]
        kes = [a - b for b, a in pts.values()]
        xc, yc = (min(bes) + max(bes)) / 2, (min(kes) + max(kes)) / 2
        xlo, xhi = xc - SPAN / 2, xc + SPAN / 2
        ylo, yhi = yc - SPAN / 2, yc + SPAN / 2
        a_lo, a_hi = int(np.floor(xlo + ylo)), int(np.ceil(xhi + yhi))
        xx = np.linspace(xlo, xhi, 50)
        for a in range(a_lo, a_hi + 1):
            ax.plot(xx, a - xx, color="0.8", lw=0.8, zorder=0)
            xt = xhi - 0.15
            yt = a - xt
            if ylo + 0.3 < yt < yhi - 0.3:
                ax.text(xt, yt, rf"$\alpha'$ = {a}", fontsize=7.5, color="0.5", ha="left", va="bottom",
                        rotation=45, rotation_mode="anchor")
        for key, (b, a) in pts.items():
            ax.plot(b, a - b, "o", ms=9, color="tab:red" if "O" in key else "tab:blue", zorder=3)
            ax.annotate(L_["names"][key] + "\n" + rf"$\alpha'$ = {a:.2f}", xy=(b, a - b), xytext=(8, 8),
                        textcoords="offset points", fontsize=8.5, **F)
        ax.set_xlim(xhi, xlo)
        ax.set_ylim(ylo, yhi)
        ax.set_aspect("equal")
        ax.set_xlabel(L_["wag_x"], fontsize=10.5, **F)
        ax.set_ylabel(L_["wag_y"], fontsize=10.5, **F)
        ax.set_title(title, fontsize=10.5, **F)
        ax.grid(alpha=0.2)
    fig4.suptitle(L_["fig4"], fontsize=12, **F)
    fig4.tight_layout()
    fig4.savefig(os.path.join(out, "fig4-wagner-plot.png"), dpi=150)
    plt.close(fig4)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

print()
print(f"그림1 면적비: Ti 2p1/2 : 2p3/2 = {ratio_ti:.3f}, d3/2 : d5/2 = {ratio_d:.3f}")
print(f"그림3 DS: 최대점 {peak_pos:+.2f} eV, 높은 BE 쪽 면적 몫 {100 * tail_high:.0f}% "
      f"(대칭이면 50%)")
print(f"  Voigt 1개: 중심 {p1[0]:+.2f} eV, RMS {rms1:.4f}")
print(f"  Voigt 2개: 주 {p2[0]:+.2f} eV, 가짜 {p2[4]:+.2f} eV, 가짜 성분 면적 몫 "
      f"{100 * fake_share:.0f}%, RMS {rms2:.4f} (잡음 σ = {noise_level})")
for elem, pts in WAGNER.items():
    for key, (b, a) in pts.items():
        print(f"  {key}: E_B {b:.2f}, α' {a:.2f}, E_K(Auger) = {a - b:.2f} eV")
print(f"  Zn: ΔE_B = {WAGNER['Zn']['ZnO'][0] - WAGNER['Zn']['Zn0'][0]:+.2f} eV, "
      f"Δα' = {WAGNER['Zn']['ZnO'][1] - WAGNER['Zn']['Zn0'][1]:+.2f} eV")
print(f"  Cu(0)→Cu2O: ΔE_B = {WAGNER['Cu']['Cu2O'][0] - WAGNER['Cu']['Cu0'][0]:+.2f} eV, "
      f"Δα' = {WAGNER['Cu']['Cu2O'][1] - WAGNER['Cu']['Cu0'][1]:+.2f} eV")
print(f"Al Kα ↔ Mg Kα: 오제 피크의 겉보기 결합에너지 이동 = {1486.6 - 1253.6:.1f} eV")

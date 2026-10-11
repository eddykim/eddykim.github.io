"""파동광학 배경이론 2편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/wave-optics-coherence/      (한국어)
      ../../assets/img/posts/wave-optics-coherence/en/   (영문)

그림 1·7 은 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import j1

from coherence import (autocorrelation, contrast, degree_from_spectrum, gaussian_line, michelson_interferogram,
                       observed_film_reflectance, oct_ascan, random_wave, spectrum_from_interferogram,
                       speckle_pattern, visibility_extended_source, wavelength_band, young_extended_source)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "wave-optics-coherence")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"
NA_D1, NA_D2 = 0.588995, 0.589592            # 나트륨 D 선 (µm)


def t(L, ko, en):
    s = ko if L["lang"] == "ko" else en
    return s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$").replace("θ", r"$\theta$").replace("τ", r"$\tau$")


def font(L):
    return {"fontfamily": ["AppleGothic", "DejaVu Sans"]} if L["lang"] == "ko" else {}


def setfonts(L):
    plt.rcParams["font.family"] = ["AppleGothic", "DejaVu Sans"] if L["lang"] == "ko" else ["DejaVu Sans"]


def legend(ax, L, **kw):
    ax.legend(prop=dict(family=plt.rcParams["font.family"], size=kw.pop("size", 8.5)), **kw)


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    fig.savefig(os.path.join(L["dir"], name), dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("  저장:", os.path.join(L["dir"], name))


# ------------------------------------------------------------------ 그림 2: 무작위 빛과 결맞음 시간

def figure2(L):
    rng = np.random.default_rng(7)
    n_t, dt = 2**16, 0.05                    # 시간 단위: 반송파 주기 = 1
    nu = np.linspace(0, 2, 4001)
    fig, axes = plt.subplots(2, 3, figsize=(13, 5.6), gridspec_kw=dict(width_ratios=[1.6, 1, 1]))
    for row, (width, col, name_ko, name_en) in enumerate([(0.02, C_RAY, "좁은 스펙트럼 (선폭 2%)", "narrow spectrum (2% linewidth)"),
                                                           (0.2, C_RED, "넓은 스펙트럼 (선폭 20%)", "broad spectrum (20% linewidth)")]):
        S = gaussian_line(nu, 1.0, width)
        U = random_wave(S, nu, n_t, dt, rng)
        U /= np.sqrt(np.mean(np.abs(U) ** 2))
        tt = np.arange(n_t) * dt
        a, b, c = axes[row]
        sel = tt < 120
        a.plot(tt[sel], U.real[sel], color=col, lw=0.7)
        a.plot(tt[sel], np.abs(U)[sel], color="k", lw=0.8, alpha=0.6)
        a.set_ylabel(t(L, "전기장", "field"), **font(L))
        a.set_title(t(L, f"(a) {name_ko}: 전기장(색)과 포락선(검정)", f"(a) {name_en}: field (colour) and envelope (black)") if row == 0
                    else t(L, f"(d) {name_ko}", f"(d) {name_en}"), fontsize=10.5, **font(L))
        lags = int(60 / dt)
        g = autocorrelation(U, lags)
        tau = np.arange(lags) * dt
        b.plot(tau, np.abs(g), color=col, lw=2, label=t(L, "무작위 파동에서 잰 값", "measured on the random wave"))
        g_th = degree_from_spectrum(nu, S, tau)        # 같은 식을 시간 단위로 (σ ↔ ν, Δ ↔ τ)
        b.plot(tau[::25], np.abs(g_th)[::25], "o", mfc="none", color="k", ms=4, label=t(L, "스펙트럼의 푸리에 변환", "Fourier transform of spectrum"))
        b.set_ylim(0, 1.08)
        b.set_ylabel("|g(τ)|".replace("τ", r"$\tau$"))
        legend(b, L, loc="upper right", size=8)
        b.set_title(t(L, "(b) 결맞음도: 지연에 따른 상관" if row == 0 else "(e)", "(b) Degree of coherence vs delay" if row == 0 else "(e)"),
                    fontsize=10.5, **font(L))
        c.plot(nu, S / S.max(), color=col, lw=2)
        c.set_xlim(0.4, 1.6)
        c.set_ylabel(t(L, "스펙트럼 (상대)", "spectrum (relative)"), **font(L))
        c.set_title(t(L, "(c) 파워 스펙트럼" if row == 0 else "(f)", "(c) Power spectrum" if row == 0 else "(f)"), fontsize=10.5, **font(L))
    axes[1, 0].set_xlabel(t(L, "시간 (반송파 주기 단위)", "time (carrier periods)"), **font(L))
    axes[1, 1].set_xlabel(t(L, "지연 τ (반송파 주기 단위)", "delay τ (carrier periods)"), **font(L))
    axes[1, 2].set_xlabel(t(L, "주파수 (반송파 = 1)", "frequency (carrier = 1)"), **font(L))
    fig.tight_layout()
    save(fig, L, "fig2-random-light-coherence-time.png")


# ------------------------------------------------------------------ 그림 3: 마이컬슨 간섭무늬와 푸리에 변환 분광

def figure3(L):
    fig = plt.figure(figsize=(13, 7.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.55, 1], hspace=0.55)
    # (a) LED: 짧은 무늬 묶음
    s0, ds = wavelength_band(0.85, 0.05)
    sig = np.linspace(s0 - 6 * ds, s0 + 6 * ds, 3001)
    S = gaussian_line(sig, s0, ds)
    D = np.linspace(-40, 40, 8001)
    I = michelson_interferogram(sig, S, D)
    a = fig.add_subplot(gs[0, 0])
    a.plot(D, I, color=C_RED, lw=0.8)
    a.set_xlim(-40, 40)
    a.set_ylabel(t(L, "세기", "intensity"), **font(L))
    a.set_title(t(L, "(a) LED (850 nm, 반치폭 50 nm): 광로차 ±10 µm 남짓에서만 무늬", "(a) LED (850 nm, 50 nm FWHM): fringes only within about ±10 µm"),
                fontsize=10.5, **font(L))
    a.set_xlabel(t(L, r"광로차 $\Delta$ ($\mu$m)", r"path difference $\Delta$ ($\mu$m)"), **font(L))
    # (b) 나트륨 D 이중선: 맥놀이
    sigN = np.linspace(1 / 0.5900, 1 / 0.5884, 30001)
    SN = gaussian_line(sigN, 1 / NA_D1, 2e-6) + gaussian_line(sigN, 1 / NA_D2, 2e-6)
    DN = np.linspace(0, 1800, 180001)
    IN = michelson_interferogram(sigN, SN, DN)
    gN = np.abs(degree_from_spectrum(sigN, SN, DN[::50]))
    b = fig.add_subplot(gs[1, 0])
    b.plot(DN / 1000, IN, color=C_KEY, lw=0.3, alpha=0.7)
    b.plot(DN[::50] / 1000, 0.5 * (1 + gN), color="k", lw=1.2)
    b.plot(DN[::50] / 1000, 0.5 * (1 - gN), color="k", lw=1.2)
    beat = 0.5893**2 / (NA_D2 - NA_D1)
    for k in range(1, 4):
        b.axvline(k * beat / 1000, color=C_GREY, ls=":", lw=1)
    b.set_xlabel(t(L, r"광로차 $\Delta$ (mm)", r"path difference $\Delta$ (mm)"), **font(L))
    b.set_ylabel(t(L, "세기", "intensity"), **font(L))
    b.set_title(t(L, f"(b) 나트륨 D 이중선: 광로차 {beat / 1000:.2f} mm 마다 무늬가 사라졌다 돌아온다 (검정: 가시도 포락선)",
                  f"(b) Sodium D doublet: fringes vanish and return every {beat / 1000:.2f} mm of path difference (black: visibility envelope)"),
                fontsize=10.5, **font(L))
    # (c) 백색광: 한 줄의 흰 무늬
    s0w, dsw = wavelength_band(0.57, 0.25)
    sigW = np.linspace(1 / 0.85, 1 / 0.38, 3001)
    SW = gaussian_line(sigW, s0w, dsw)
    DW = np.linspace(-3, 3, 3001)
    c = fig.add_subplot(gs[2, 0])
    c.plot(DW, michelson_interferogram(sigW, SW, DW), color=C_GRN, lw=1.2)
    c.set_xlabel(t(L, r"광로차 $\Delta$ ($\mu$m)", r"path difference $\Delta$ ($\mu$m)"), **font(L))
    c.set_ylabel(t(L, "세기", "intensity"), **font(L))
    c.set_title(t(L, "(c) 백색광 (중심 570 nm, 반치폭 250 nm): 몇 줄뿐", "(c) White light (570 nm centre, 250 nm FWHM): only a few fringes"),
                fontsize=10.5, **font(L))
    # (d) 푸리에 변환 분광: (b) 의 간섭무늬에서 이중선을 되찾는다
    d_ax = fig.add_subplot(gs[:, 1])
    sig_out = np.linspace(1 / 0.5900, 1 / 0.5884, 1601)
    for Dmax, col, lab_ko, lab_en in [(1800.0, C_RAY, "광로차 1.8 mm 까지 잼", "scan to 1.8 mm"),
                                       (300.0, C_RED, "광로차 0.3 mm 까지만", "scan to 0.3 mm only")]:
        sel = DN <= Dmax
        Dsym = np.concatenate([-DN[sel][::-1], DN[sel][1:]])
        Isym = np.concatenate([IN[sel][::-1], IN[sel][1:]])
        spec = spectrum_from_interferogram(Dsym, Isym, sig_out)
        d_ax.plot(1e3 / sig_out, spec / spec.max(), color=col, lw=1.8, label=t(L, lab_ko, lab_en))
    for lam in (NA_D1, NA_D2):
        d_ax.axvline(lam * 1e3, color="k", ls=":", lw=0.8)
    d_ax.set_xlim(588.5, 590.1)
    d_ax.set_xticks([588.5, 589.0, 589.5, 590.0])
    d_ax.set_xlabel(t(L, "파장 (nm)", "wavelength (nm)"), **font(L))
    d_ax.set_ylabel(t(L, "되찾은 스펙트럼 (상대)", "recovered spectrum (relative)"), **font(L))
    legend(d_ax, L, loc="upper right")
    d_ax.set_title(t(L, "(d) (b)의 간섭무늬를 푸리에 변환하면\n두 선이 갈라진다 (점선: 589.0, 589.6 nm)",
                     "(d) Fourier transforming the interferogram of (b)\nseparates the two lines (dotted: 589.0, 589.6 nm)"),
                   fontsize=10.5, **font(L))
    save(fig, L, "fig3-michelson-interferograms.png")
    print(f"    그림3: 나트륨 맥놀이 광로차 {beat:.1f} µm (거울 {beat / 2:.1f} µm)")


# ------------------------------------------------------------------ 그림 4: 박막 무늬는 언제 보이는가

def figure4(L):
    n = 1.33
    t_um = np.linspace(0.0, 30.0, 6001)
    sig = np.linspace(1 / 0.75, 1 / 0.40, 3001)
    s0, ds = wavelength_band(0.55, 0.15)
    S_eye = gaussian_line(sig, s0, ds)                   # 눈·카메라가 받는 넓은 띠
    s1, ds1 = wavelength_band(0.55, 0.002)
    sig1 = np.linspace(s1 - 6 * ds1, s1 + 6 * ds1, 1201)
    S_spec = gaussian_line(sig1, s1, ds1)                # 분광기 한 화소의 띠 (2 nm)
    R_eye = observed_film_reflectance(n, t_um, sig, S_eye)
    R_spec = observed_film_reflectance(n, t_um, sig1, S_spec)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.0), gridspec_kw=dict(width_ratios=[1.5, 1]))
    a.plot(t_um, R_spec * 100, color=C_RAY, lw=0.6, label=t(L, "분광기 한 화소 (반치폭 2 nm)", "one spectrometer pixel (2 nm FWHM)"))
    a.plot(t_um, R_eye * 100, color=C_RED, lw=1.6, label=t(L, "눈 (반치폭 150 nm)", "the eye (150 nm FWHM)"))
    a.set_xlabel(t(L, "물막 두께 (µm)", "water film thickness (µm)"), **font(L))
    a.set_ylabel(t(L, "반사율 (%)", "reflectance (%)"), **font(L))
    a.set_ylim(0, 9.5)
    legend(a, L, loc="upper right", ncol=2)
    a.set_title(t(L, "(a) 같은 막을 눈과 분광기로 볼 때 (중심 550 nm)", "(a) The same film seen by the eye and by a spectrometer (550 nm centre)"),
                fontsize=10.5, **font(L))
    # (b) 두꺼운 막의 반사 스펙트럼: 분광기에는 무늬가 보인다
    lam = np.linspace(0.50, 0.60, 4001)
    for tt, col in [(5.0, C_GRN), (30.0, C_PUR)]:
        # 분광기 화소마다: 그 화소의 2 nm 띠를 파장 l 로 옮겨 평균한 반사율
        Rlam = np.array([observed_film_reflectance(n, [tt], sig1 - s1 + 1 / l, S_spec)[0] for l in lam])
        b.plot(lam * 1e3, Rlam * 100, color=col, lw=0.8, label=t(L, f"두께 {tt:g} µm", f"thickness {tt:g} µm"))
    b.set_xlabel(t(L, "파장 (nm)", "wavelength (nm)"), **font(L))
    b.set_ylabel(t(L, "반사율 (%)", "reflectance (%)"), **font(L))
    b.set_ylim(0, 9.5)
    legend(b, L, loc="upper right")
    b.set_title(t(L, "(b) 분광기로 본 반사 스펙트럼: 무늬 간격이 두께를 알려 준다", "(b) Reflection spectra: the fringe spacing gives the thickness"),
                fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig4-thin-film-coherence.png")


# ------------------------------------------------------------------ 그림 5: OCT

def figure5(L):
    s0, ds = wavelength_band(0.84, 0.05)
    sig = np.linspace(s0 - 6 * ds, s0 + 6 * ds, 2001)
    S = gaussian_line(sig, s0, ds)
    layers = [(100.0, 0.05), (130.0, 0.03), (138.0, 0.03), (200.0, 0.02)]
    z = np.linspace(80, 220, 28001)
    I = oct_ascan(sig, S, z, layers)
    ac = I - np.median(I)
    fig, (a, b) = plt.subplots(2, 1, figsize=(12, 5.2), sharex=True, gridspec_kw=dict(height_ratios=[1, 0.6]))
    a.plot(z, ac, color=C_RAY, lw=0.4)
    for zi, ri in layers:
        a.axvline(zi, color=C_GREY, ls=":", lw=1)
    a.set_ylabel(t(L, "간섭 신호", "interference signal"), **font(L))
    a.set_title(t(L, "(a) 시간 영역 OCT: 기준 거울을 움직이며 잰 간섭 신호 (840 nm, 반치폭 50 nm)",
                  "(a) Time-domain OCT: interference signal as the reference mirror moves (840 nm, 50 nm FWHM)"),
                fontsize=10.5, **font(L))
    # (b) 포락선: 힐베르트 변환으로
    from scipy.signal import hilbert
    envelope = np.abs(hilbert(ac))
    b.plot(z, envelope, color=C_RED, lw=1.5)
    for zi, ri in layers:
        b.axvline(zi, color=C_GREY, ls=":", lw=1)
    b.set_xlabel(t(L, "깊이 (µm, 공기 환산)", "depth (µm, in air)"), **font(L))
    b.set_ylabel(t(L, "포락선", "envelope"), **font(L))
    b.set_title(t(L, "(b) 포락선: 경계마다 봉우리 하나. 8 µm 떨어진 두 경계도 갈라진다",
                  "(b) Envelope: one peak per boundary; two boundaries 8 µm apart are resolved"), fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig5-oct-ascan.png")
    # 반치폭
    sel = (z > 90) & (z < 110)
    e = envelope[sel]; zz = z[sel]
    half = e.max() / 2
    above = zz[e > half]
    print(f"    그림5: 첫 경계 포락선 반치폭 {above.max() - above.min():.2f} µm")


# ------------------------------------------------------------------ 그림 6: 공간 결맞음

def figure6(L):
    lam = 0.55
    th_sun = np.radians(0.53)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.3), gridspec_kw=dict(width_ratios=[1, 1.15]))
    # (a) 광원 크기에 따른 무늬
    d = 50.0
    x = np.linspace(-2, 2, 801) * lam * 1e6 / d       # 두 주기
    for th, col, lab in [(0.0, "k", "0"), (0.3 * th_sun, C_RAY, "0.16°"), (0.6 * th_sun, C_KEY, "0.32°"),
                         (th_sun, C_RED, "0.53°")]:
        I = young_extended_source(x, d, 1e6, lam, max(th, 1e-12))
        a.plot(x / 1000, I / 2, color=col, lw=1.6, label=t(L, f"광원 각 지름 {lab}", f"source diameter {lab}"))
    a.set_xlabel(t(L, "1 m 뒤 스크린 위 위치 (mm)", "position on a screen 1 m away (mm)"), **font(L))
    a.set_ylabel(t(L, "세기 (평균 = 0.5)", "intensity (mean = 0.5)"), **font(L))
    legend(a, L, loc="lower right", size=8)
    a.set_ylim(0, 1.25)
    a.set_title(t(L, "(a) 간격 50 µm 두 구멍, 550 nm: 광원이 클수록 무늬가 흐려진다",
                  "(a) Holes 50 µm apart, 550 nm: larger sources wash out the fringes"), fontsize=10.5, **font(L))
    # (b) 가시도 vs 구멍 간격 (태양)
    ds_ = np.linspace(1, 200, 120)
    V_sim = np.array([visibility_extended_source(dd, lam, th_sun, n_pts=2001) for dd in ds_])
    xx = np.pi * ds_ * th_sun / lam
    b.plot(ds_, np.abs(2 * j1(xx) / xx), color="k", lw=1.2, ls="--", label=t(L, r"$|2J_1(x)/x|$, $x = \pi d\theta_s/\lambda$", r"$|2J_1(x)/x|$, $x = \pi d\theta_s/\lambda$"))
    b.plot(ds_[::3], V_sim[::3], "o", mfc="none", color=C_RED, ms=5, label=t(L, "독립 점광원들의 무늬를 더한 값", "sum of fringes from independent points"))
    zero = 1.22 * lam / th_sun
    b.axvline(zero, color=C_GREY, ls=":", lw=1)
    b.text(zero + 3, 0.85, f"{zero:.0f} µm", fontsize=9.5)
    b.set_xlabel(t(L, "두 구멍 사이 거리 d (µm)", "hole separation d (µm)"), **font(L))
    b.set_ylabel(t(L, "가시도 = 공간 결맞음도", "visibility = degree of spatial coherence"), **font(L))
    legend(b, L, loc="upper right")
    b.set_title(t(L, "(b) 햇빛(각 지름 0.53°, 550 nm)의 공간 결맞음", "(b) Spatial coherence of sunlight (0.53°, 550 nm)"),
                fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig6-spatial-coherence.png")
    print(f"    그림6: 햇빛 결맞음 첫 영점 {zero:.1f} µm")


# ------------------------------------------------------------------ 그림 8: 스페클

def figure8(L):
    rng = np.random.default_rng(11)
    n = 512
    I1 = speckle_pattern(n, 0.12, rng)
    I16 = sum(speckle_pattern(n, 0.12, rng) for _ in range(16)) / 16
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.15]))
    for ax, I, ttl_ko, ttl_en in [(axes[0], I1, "(a) 결맞는 빛 한 번", "(a) Coherent light, one realization"),
                                  (axes[1], I16, "(b) 독립인 스페클 16장의 평균", "(b) Average of 16 independent speckle patterns")]:
        ax.imshow(I[:256, :256], cmap="gray", vmin=0, vmax=3)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(t(L, f"{ttl_ko}\n대비 {contrast(I):.2f}", f"{ttl_en}\ncontrast {contrast(I):.2f}"), fontsize=10.5, **font(L))
    c = axes[2]
    bins = np.linspace(0, 6, 61)
    c.hist(I1.ravel(), bins=bins, density=True, color=C_RAY, alpha=0.6, label=t(L, "한 번", "one"))
    c.hist(I16.ravel(), bins=bins, density=True, color=C_RED, alpha=0.6, label=t(L, "16장 평균", "average of 16"))
    xs = np.linspace(0, 6, 200)
    c.plot(xs, np.exp(-xs), "k--", lw=1.2, label=r"$e^{-I/\langle I\rangle}$")
    c.set_yscale("log")
    c.set_ylim(1e-3, 3)
    c.set_xlabel(t(L, "세기 / 평균 세기", "intensity / mean"), **font(L))
    c.set_ylabel(t(L, "확률 밀도", "probability density"), **font(L))
    legend(c, L, loc="upper right")
    c.set_title(t(L, "(c) 세기의 분포", "(c) Intensity distribution"), fontsize=10.5, **font(L))
    c.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    fig.tight_layout()
    save(fig, L, "fig8-speckle.png")
    print(f"    그림8: 대비 {contrast(I1):.3f} → {contrast(I16):.3f}")


FIGS = {2: figure2, 3: figure3, 4: figure4, 5: figure5, 6: figure6, 8: figure8}


def main():
    which = [int(a) for a in sys.argv[1:]] or sorted(FIGS)
    for lang, d in [("ko", BASE_DIR), ("en", os.path.join(BASE_DIR, "en"))]:
        L = {"lang": lang, "dir": d}
        setfonts(L)
        for k in which:
            print(f"[{lang}] 그림 {k}")
            FIGS[k](L)


if __name__ == "__main__":
    main()

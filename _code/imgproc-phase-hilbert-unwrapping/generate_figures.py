"""이미지처리 6편 그림 6개 생성 (한국어판·영문판)."""
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import hilbert

from phase_extract import (phase_error, phase_fourier, phase_hilbert,
                           phase_shift_5_hariharan, phase_shift_n, residues,
                           unwrap_itoh, wrap)
from phase_signal import background, fringe, fringe_2d, modulation, true_phase

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
BASE = os.path.join(os.path.dirname(__file__), "..", "..",
                    "assets", "img", "posts", "imgproc-phase-hilbert-unwrapping")

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "pos": "위치 (화소)", "inten": "세기", "shift": "천이량 $\\delta$ (rad)",
        "fringe": "간섭 무늬 $I(x)$", "bg": "배경 $a(x)$", "env": "포락선 $a \\pm b$",
        "onepix": "한 화소에서 천이량에 따른 세기", "fitted": "맞춘 정현파",
        "phi_here": "이 정현파의 위상이 $\\phi$",
        "fig1": "그림1. 한 장에는 미지수 셋이 겹쳐 있고 식은 하나뿐이다",
        "ac": "배경을 뺀 신호", "hil": "힐베르트 변환", "amp": "순간 진폭 $b(x)$",
        "true_amp": "참 변조 $b(x)$", "witherr": "배경을 두고 계산한 위상 오차",
        "clean": "배경을 빼고 계산한 위상 오차", "perr": "위상 오차 (rad)",
        "fig2": "그림2. 해석 신호의 크기가 진폭, 편각이 위상이다",
        "noise_ax": "프레임 노이즈", "rms": "위상 오차 RMS (rad)",
        "m_h": "힐베르트 (1장)", "m_f": "푸리에 (1장)", "m_p4": "위상천이 (4장)",
        "band_ax": "측대역 폭 ($f_0$ 단위)", "e_noise": "노이즈 하 오차 (노이즈 0.05)",
        "e_shape": "무잡음, 형상이 급한 구간의 오차",
        "fig3": "그림3. 푸리에법의 노이즈 내성은 대역통과에서 오고, 그 폭이 또 하나의 교환이다",
        "eps": "천이량 오차 $\\varepsilon$", "s3": "3단계", "s4": "4단계",
        "s5": "5단계 최소자승", "sh": "Hariharan 5단계",
        "slope1": "1차 기울기", "slope2": "2차 기울기",
        "fig4": "그림4. Hariharan 5단계만 천이량 오차에 2차로 반응한다",
        "carrier_ax": "화소당 위상 증가 (rad)", "uerr": "언래핑 최대 오차 (rad)",
        "limit": "화소당 $\\pi$", "jump": "$2\\pi$ 도약 발생률 (%)",
        "fig5": "그림5. 위상차가 $\\pi$ 를 넘는 순간 오차가 $10^{-13}$ 에서 수백 rad 으로 뛴다",
        "wrapped": "래핑된 위상", "res_u": "잔여점 — 언샘플링 기인",
        "res_n": "잔여점 — 노이즈 기인", "pathA": "경로 A", "pathB": "경로 B",
        "fig6": "그림6. 잔여점이 있으면 언래핑 결과가 경로에 따라 달라진다",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "pos": "Position (pixel)", "inten": "Intensity", "shift": "Shift $\\delta$ (rad)",
        "fringe": "Fringe $I(x)$", "bg": "Background $a(x)$", "env": "Envelope $a \\pm b$",
        "onepix": "Intensity vs shift at one pixel", "fitted": "Fitted sinusoid",
        "phi_here": "its phase is $\\phi$",
        "fig1": "Fig 1. One frame holds three unknowns and offers one equation",
        "ac": "Background removed", "hil": "Hilbert transform", "amp": "Instantaneous amplitude",
        "true_amp": "True modulation $b(x)$", "witherr": "Phase error, background left in",
        "clean": "Phase error, background removed", "perr": "Phase error (rad)",
        "fig2": "Fig 2. The analytic signal gives amplitude as magnitude and phase as argument",
        "noise_ax": "Frame noise", "rms": "Phase error RMS (rad)",
        "m_h": "Hilbert (1 frame)", "m_f": "Fourier (1 frame)", "m_p4": "Phase shifting (4 frames)",
        "band_ax": "Sideband width (units of $f_0$)", "e_noise": "Error under noise 0.05",
        "e_shape": "Error on steep shape, noiseless",
        "fig3": "Fig 3. The Fourier method's noise immunity comes from its bandpass, whose width is another trade",
        "eps": "Shift error $\\varepsilon$", "s3": "3-step", "s4": "4-step",
        "s5": "5-step least squares", "sh": "Hariharan 5-step",
        "slope1": "first order", "slope2": "second order",
        "fig4": "Fig 4. Only the Hariharan 5-step responds quadratically to shift error",
        "carrier_ax": "Phase increment per pixel (rad)", "uerr": "Max unwrap error (rad)",
        "limit": "$\\pi$ per pixel", "jump": "$2\\pi$ jump rate (%)",
        "fig5": "Fig 5. Past $\\pi$ the error jumps from $10^{-13}$ to hundreds of radians",
        "wrapped": "Wrapped phase", "res_u": "Residues from undersampling",
        "res_n": "Residues from noise", "pathA": "Path A", "pathB": "Path B",
        "fig6": "Fig 6. With residues present, unwrapping depends on the path",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
X = np.arange(512.0)
SH4 = np.arange(4) * np.pi / 2
I0, A0, B0, PHI = fringe(X)
PIX = 300

# 그림1: 한 화소에서 천이량에 따른 세기
dfine = np.linspace(0, 2 * np.pi, 200)
Ipix = A0[PIX] + B0[PIX] * np.cos(PHI[PIX] + dfine)
Ipix4 = A0[PIX] + B0[PIX] * np.cos(PHI[PIX] + SH4)

# 그림2: 해석 신호
from scipy.ndimage import gaussian_filter1d
AC = I0 - gaussian_filter1d(I0, 15)
Z = hilbert(AC)
ph_clean, amp_clean = phase_hilbert(I0, detrend_sigma=15)
ph_raw, _ = phase_hilbert(I0, detrend_sigma=None)
err_clean = phase_error(ph_clean, PHI)
err_raw = phase_error(ph_raw, PHI)

# 그림3: 노이즈에 따른 세 방법
NOISE = np.array([0.002, 0.005, 0.01, 0.02, 0.05, 0.10])
TRIALS = 60
rng3 = np.random.default_rng(2)
res3 = {"m_h": [], "m_f": [], "m_p4": []}
for nz in NOISE:
    acc = {k: [] for k in res3}
    for _ in range(TRIALS):
        one = fringe(X, noise=nz, rng=rng3)[0]
        acc["m_h"].append(np.sqrt(np.mean(phase_error(
            phase_hilbert(one, detrend_sigma=15)[0], PHI)[30:-30] ** 2)))
        acc["m_f"].append(np.sqrt(np.mean(phase_error(
            phase_fourier(one)[0], PHI)[30:-30] ** 2)))
        fr = [fringe(X, noise=nz, rng=rng3, phase_shift=d)[0] for d in SH4]
        acc["m_p4"].append(np.sqrt(np.mean(phase_error(
            phase_shift_n(fr, SH4)[0], PHI)[30:-30] ** 2)))
    for k in res3:
        res3[k].append(np.mean(acc[k]))
res3 = {k: np.array(v) for k, v in res3.items()}

# 푸리에 측대역 폭의 교환. 좁히면 노이즈는 막지만 급한 형상이 뭉개진다.
_clean = fringe(X)[0]
_F = np.fft.fft(_clean - _clean.mean())
_fr = np.fft.fftfreq(len(X))
F0 = _fr[np.argmax(np.abs(_F) * (_fr > 0))]
_d2 = np.abs(np.gradient(np.gradient(PHI)))
STEEP = _d2 > np.percentile(_d2, 85)
STEEP[:40] = False
STEEP[-40:] = False
BANDS = np.array([0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9])
rng3b = np.random.default_rng(6)
band_noise, band_shape = [], []
for w in BANDS:
    acc = []
    for _ in range(40):
        noisy = fringe(X, noise=0.05, rng=rng3b)[0]
        acc.append(np.sqrt(np.mean(phase_error(
            phase_fourier(noisy, band=(F0, w * F0))[0], PHI)[30:-30] ** 2)))
    band_noise.append(np.mean(acc))
    band_shape.append(np.sqrt(np.mean(phase_error(
        phase_fourier(_clean, band=(F0, w * F0))[0], PHI)[STEEP] ** 2)))
band_noise, band_shape = np.array(band_noise), np.array(band_shape)

# 그림4: 디튜닝
EPS = np.array([0.005, 0.01, 0.02, 0.05, 0.10, 0.20])
det = {"s3": [], "s4": [], "s5": [], "sh": []}
for e in EPS:
    for key, N in (("s3", 3), ("s4", 4), ("s5", 5)):
        nom = np.arange(N) * 2 * np.pi / N if N == 3 else np.arange(N) * np.pi / 2
        fr = [fringe(X, phase_shift=d)[0] for d in nom * (1 + e)]
        det[key].append(np.sqrt(np.mean(phase_error(phase_shift_n(fr, nom)[0], PHI) ** 2)))
    fr5 = [fringe(X, phase_shift=d)[0] for d in np.arange(5) * np.pi / 2 * (1 + e)]
    det["sh"].append(np.sqrt(np.mean(phase_error(phase_shift_5_hariharan(fr5), PHI) ** 2)))
det = {k: np.array(v) for k, v in det.items()}

# 그림5: 언래핑 한계
CARRIERS = np.array([0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 2.9, 3.0, 3.05, 3.1, 3.2, 3.5, 4.0])
uerr = []
for c in CARRIERS:
    _, _, _, p = fringe(X, carrier=c)
    fr = [fringe(X, carrier=c, phase_shift=d)[0] for d in SH4]
    w, _, _ = phase_shift_n(fr, SH4)
    u = unwrap_itoh(w)
    e = u - p
    uerr.append(np.abs(e - e[0]).max())
uerr = np.array(uerr)
step_max = np.array([np.abs(np.diff(fringe(X, carrier=c)[3])).max() for c in CARRIERS])

NZ5 = np.array([0.0, 0.02, 0.05, 0.08, 0.10, 0.14, 0.20, 0.28])
rng5 = np.random.default_rng(0)
jump = []
for nz in NZ5:
    f = 0
    for _ in range(200):
        fr = [fringe(X, noise=nz, rng=rng5, phase_shift=d)[0] for d in SH4]
        w, _, _ = phase_shift_n(fr, SH4)
        e = unwrap_itoh(w) - PHI
        f += abs(np.ptp(e - e[0])) > np.pi
    jump.append(100.0 * f / 200)
jump = np.array(jump)

# 그림6: 2D 잔여점
def measured2d(amp, nz, seed=1):
    fr = [fringe_2d(amp=amp, noise=nz, seed=seed + k, phase_shift=d)[0]
          for k, d in enumerate(SH4)]
    W, _, _ = phase_shift_n([f.ravel() for f in fr], SH4)
    return W.reshape(fr[0].shape)


W_ok = measured2d(14, 0.0)
W_under = measured2d(140, 0.0)
W_noisy = measured2d(14, 0.35)
R_under, R_noisy = residues(W_under), residues(W_noisy)

# 경로 의존성: 잔여점 하나를 감싸는 두 경로로 적분한 값의 차이
def integrate_path(W, path):
    vals = [W[path[0]]]
    for (y0, x0), (y1, x1) in zip(path[:-1], path[1:]):
        vals.append(vals[-1] + wrap(W[y1, x1] - W[y0, x0]))
    return np.array(vals)


ry, rx = np.argwhere(R_under != 0)[len(np.argwhere(R_under != 0)) // 2]
y0, y1 = max(ry - 18, 0), min(ry + 18, W_under.shape[0] - 1)
x0, x1 = max(rx - 18, 0), min(rx + 18, W_under.shape[1] - 1)
pathA = ([(y0, x) for x in range(x0, x1 + 1)] + [(y, x1) for y in range(y0 + 1, y1 + 1)])
pathB = ([(y, x0) for y in range(y0, y1 + 1)] + [(y1, x) for x in range(x0 + 1, x1 + 1)])
endA = integrate_path(W_under, pathA)[-1]
endB = integrate_path(W_under, pathB)[-1]


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    F, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.06, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **F)
        fig.tight_layout(rect=rect)

    # 그림1
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.9))
    sl = slice(220, 380)
    ax[0].plot(X[sl], I0[sl], "-", lw=1.4, color="#555", label=L["fringe"])
    ax[0].plot(X[sl], A0[sl], "--", lw=1.3, color="#2b5c9b", label=L["bg"])
    ax[0].plot(X[sl], (A0 + B0)[sl], ":", lw=1.2, color="#c0392b", label=L["env"])
    ax[0].plot(X[sl], (A0 - B0)[sl], ":", lw=1.2, color="#c0392b")
    ax[0].axvline(PIX, color="k", lw=0.8, ls="-.")
    ax[0].set_xlabel(L["pos"], **F); ax[0].set_ylabel(L["inten"], **F)
    ax[0].legend(prop=LG, fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(dfine, Ipix, "-", lw=1.5, color="#1e8449", label=L["fitted"])
    ax[1].plot(SH4, Ipix4, "o", ms=7, color="#c0392b")
    ax[1].axvline((-PHI[PIX]) % (2 * np.pi), ls=":", lw=1.2, color="k")
    ax[1].set_xlabel(L["shift"], **F); ax[1].set_ylabel(L["inten"], **F)
    ax[1].set_title(L["onepix"], fontsize=9.5, **F)
    ax[1].legend(prop=LG, fontsize=8.5); ax[1].grid(alpha=0.3)
    cap(fig, "fig1", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig1-where-is-phase.png"), dpi=150)
    plt.close(fig)

    # 그림2
    fig, ax = plt.subplots(1, 3, figsize=(11.4, 3.7))
    ax[0].plot(X[sl], AC[sl], "-", lw=1.3, color="#555", label=L["ac"])
    ax[0].plot(X[sl], np.imag(Z)[sl], "-", lw=1.1, color="#e67e22", label=L["hil"])
    ax[0].set_xlabel(L["pos"], **F); ax[0].legend(prop=LG, fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(X, amp_clean, "-", lw=1.4, color="#1e8449", label=L["amp"])
    ax[1].plot(X, B0, "--", lw=1.2, color="k", label=L["true_amp"])
    ax[1].set_xlabel(L["pos"], **F); ax[1].legend(prop=LG, fontsize=8); ax[1].grid(alpha=0.3)
    ax[2].plot(X, err_raw, "-", lw=1.1, color="#c0392b", label=L["witherr"])
    ax[2].plot(X, err_clean, "-", lw=1.3, color="#2b5c9b", label=L["clean"])
    ax[2].axhline(0, color="k", lw=0.6)
    ax[2].set_xlabel(L["pos"], **F); ax[2].set_ylabel(L["perr"], **F)
    ax[2].legend(prop=LG, fontsize=8); ax[2].grid(alpha=0.3)
    cap(fig, "fig2", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig2-analytic-signal.png"), dpi=150)
    plt.close(fig)

    # 그림3
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.0))
    for k, c, m in (("m_h", "#c0392b", "o"), ("m_f", "#e67e22", "v"), ("m_p4", "#2b5c9b", "s")):
        ax[0].loglog(NOISE, res3[k], m + "-", ms=4, color=c, label=L[k])
    ax[0].set_xlabel(L["noise_ax"], **F); ax[0].set_ylabel(L["rms"], **F)
    ax[0].legend(prop=LG, fontsize=9); ax[0].grid(alpha=0.3, which="both")
    ax[1].semilogy(BANDS, band_noise, "o-", ms=4, color="#e67e22", label=L["e_noise"])
    ax[1].semilogy(BANDS, band_shape, "s-", ms=4, color="#1e8449", label=L["e_shape"])
    ax[1].set_xlabel(L["band_ax"], **F); ax[1].set_ylabel(L["rms"], **F)
    ax[1].legend(prop=LG, fontsize=8.5); ax[1].grid(alpha=0.3, which="both")
    cap(fig, "fig3", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig3-three-methods.png"), dpi=150)
    plt.close(fig)

    # 그림4
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for k, c, m in (("s3", "#c0392b", "o"), ("s4", "#e67e22", "v"),
                    ("s5", "#1e8449", "s"), ("sh", "#2b5c9b", "^")):
        ax.loglog(EPS, det[k], m + "-", ms=4, color=c, label=L[k])
    ref1 = det["s4"][0] * (EPS / EPS[0])
    ref2 = det["sh"][0] * (EPS / EPS[0]) ** 2
    ax.loglog(EPS, ref1, ":", lw=1.1, color="#888", label=L["slope1"])
    ax.loglog(EPS, ref2, "--", lw=1.1, color="#888", label=L["slope2"])
    ax.set_xlabel(L["eps"], **F); ax.set_ylabel(L["rms"], **F)
    ax.legend(prop=LG, fontsize=8.5, ncol=2); ax.grid(alpha=0.3, which="both")
    cap(fig, "fig4")
    fig.savefig(os.path.join(out, "fig4-detuning.png"), dpi=150)
    plt.close(fig)

    # 그림5
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.9))
    ax[0].semilogy(step_max, np.maximum(uerr, 1e-12), "o-", ms=4, color="#2b5c9b")
    ax[0].axvline(np.pi, ls="--", lw=1.3, color="#c0392b")
    ax[0].text(np.pi * 0.99, uerr.max() * 0.02, L["limit"], fontsize=9,
               color="#c0392b", ha="right", **F)
    ax[0].set_xlabel(L["carrier_ax"], **F); ax[0].set_ylabel(L["uerr"], **F)
    ax[0].grid(alpha=0.3, which="both")
    ax[1].plot(NZ5, jump, "s-", ms=5, color="#c0392b")
    ax[1].set_xlabel(L["noise_ax"], **F); ax[1].set_ylabel(L["jump"], **F)
    ax[1].grid(alpha=0.3)
    cap(fig, "fig5", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig5-unwrap-limit.png"), dpi=150)
    plt.close(fig)

    # 그림6
    fig, ax = plt.subplots(1, 3, figsize=(11.4, 4.3))
    ax[0].imshow(W_ok, cmap="twilight"); ax[0].set_title(L["wrapped"], fontsize=9.5, **F)
    for a, (W, R, t) in zip(ax[1:], [(W_under, R_under, L["res_u"]),
                                     (W_noisy, R_noisy, L["res_n"])]):
        a.imshow(W, cmap="twilight")
        yy, xx = np.nonzero(R > 0); a.plot(xx, yy, ".", ms=2.5, color="#e8f000")
        yy, xx = np.nonzero(R < 0); a.plot(xx, yy, ".", ms=2.5, color="#00e5ff")
        a.set_title(f"{t} ({np.abs(R).sum()})", fontsize=9.5, **F)
    pa = np.array(pathA); pb = np.array(pathB)
    ax[1].plot(pa[:, 1], pa[:, 0], "-", lw=1.4, color="#ff2d55", label=L["pathA"])
    ax[1].plot(pb[:, 1], pb[:, 0], "-", lw=1.4, color="#ffffff", label=L["pathB"])
    ax[1].legend(prop=LG, fontsize=7.5, loc="lower left")
    for a in ax:
        a.set_xticks([]); a.set_yticks([])
    cap(fig, "fig6", rect=(0, 0.06, 1, 0.94))
    fig.savefig(os.path.join(out, "fig6-residues-2d.png"), dpi=110)
    plt.close(fig)


for lang, L in LABELS.items():
    draw(L)
    print(f"[{lang}] 그림 6개 저장: {L['dir']}")

print()
print("=== 그림2: 배경 제거의 효과 ===")
print(f"  배경을 두고: 위상 오차 RMS {np.sqrt(np.mean(err_raw**2)):.4f} rad, "
      f"최대 {np.abs(err_raw).max():.4f}")
print(f"  배경을 빼고: 위상 오차 RMS {np.sqrt(np.mean(err_clean**2)):.4f} rad, "
      f"최대 {np.abs(err_clean).max():.4f}")
print(f"  순간 진폭 복원 오차 RMS {np.sqrt(np.mean((amp_clean-B0)[40:-40]**2)):.4f}")
print()
print("=== 그림3: 노이즈에 따른 위상 오차 RMS ===")
for i, nz in enumerate(NOISE):
    print(f"  노이즈 {nz:5.3f}: 힐베르트 {res3['m_h'][i]:.4f}  "
          f"푸리에 {res3['m_f'][i]:.4f}  위상천이4 {res3['m_p4'][i]:.4f}")
print()
print("=== 그림3(b): 푸리에 측대역 폭의 교환 ===")
for w, bn, bs in zip(BANDS, band_noise, band_shape):
    print(f"  폭 {w:4.2f}*f0: 노이즈 하 {bn:.4f}, 무잡음 급한 구간 {bs:.5f}")
print(f"  무잡음 위상천이 오차: {np.sqrt(np.mean(phase_error(phase_shift_n([fringe(X,phase_shift=d)[0] for d in SH4],SH4)[0],PHI)**2)):.2e}")
print()
print("=== 그림4: 디튜닝 ===")
for i, e in enumerate(EPS):
    print(f"  eps {e:5.3f}: 3단계 {det['s3'][i]:.5f}  4단계 {det['s4'][i]:.5f}  "
          f"5단계LSQ {det['s5'][i]:.5f}  Hariharan {det['sh'][i]:.6f}")
print(f"  eps 0.01 -> 0.10 배율: 4단계 {det['s4'][1]/det['s4'][1]*det['s4'][4]/det['s4'][1]:.0f}배, "
      f"Hariharan {det['sh'][4]/det['sh'][1]:.0f}배")
print()
print("=== 그림5: 언래핑 ===")
for c, sm, ue in zip(CARRIERS, step_max, uerr):
    print(f"  carrier {c:4.2f}: 최대 위상차 {sm:5.3f} rad/화소, 언래핑 오차 {ue:.4g}")
print("  노이즈별 2pi 도약률: " + ", ".join(f"{n:.2f}->{j:.0f}%" for n, j in zip(NZ5, jump)))
print()
print("=== 그림6: 2D 잔여점 ===")
print(f"  정상(amp=14, 무잡음): {np.abs(residues(W_ok)).sum()}개")
print(f"  언샘플링(amp=140, 무잡음): {np.abs(R_under).sum()}개, "
      f"최대 기울기 {max(np.abs(np.diff(fringe_2d(amp=140)[3],axis=1)).max(), np.abs(np.diff(fringe_2d(amp=140)[3],axis=0)).max()):.2f} rad/화소")
print(f"  노이즈(amp=14, 0.35): {np.abs(R_noisy).sum()}개, 최대 기울기 0.77 rad/화소")
print(f"  잔여점 하나를 감싼 두 경로의 적분 차이: {endB-endA:+.3f} rad "
      f"({(endB-endA)/2/np.pi:+.2f} x 2pi)")

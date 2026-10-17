"""계측 영상처리 3편 그림 6개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/imgproc-frequency-domain-filtering/      (한국어)
      ../../assets/img/posts/imgproc-frequency-domain-filtering/en/   (영문)

계산은 한 번만 하고 라벨만 갈아 끼운다. 시리즈의 다른 편과 같은 방식이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle
from scipy.ndimage import maximum_filter
from skimage.data import camera

from baseline_extraction import eem_baseline, fourier_lpf_baseline
from frequency_tools import (hann2d, impulse_response, lowpass, notch_filter,
                             spectrum)
from spectral_signal import RETARDER_DND, make_signal

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                        "assets", "img", "posts", "imgproc-frequency-domain-filtering")

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "grating": "합성 격자", "natural": "실제 영상", "spec": "진폭 스펙트럼 (로그)",
        "spec_zoom": "진폭 스펙트럼 (로그, 가운데 ±64빈)",
        "contaminated": "주기성 노이즈가 섞인 영상", "peaks": "스펙트럼 (원: 노이즈 봉우리)",
        "restored_on": "노치 후: 빈에 맞는 주파수",
        "restored_off": "노치 후: 빈 사이에 낀 주파수",
        "ideal": "이상적 저역통과", "gauss": "가우시안 저역통과",
        "filt": "필터 $H$", "psf_cut": "공간 영역 단면", "result": "필터 통과 결과",
        "raw_spec": "창 없이 구한 스펙트럼", "win_spec": "한(Hann) 창을 곱한 뒤",
        "prof_ax": "스펙트럼 세로선 단면 (직류 성분으로 정규화)", "row_ax": "세로 주파수 (빈)",
        "prof_raw": "창 없음", "prof_win": "한(Hann) 창",
        "wavenumber": r"파수 $\sigma$ (1/$\mu$m)", "intensity": "세기 (임의 단위)",
        "signal": "분광 간섭 신호", "upper": "상부 포락선", "lower": "하부 포락선",
        "eem_bl": "EEM 기저 신호", "true_bl": "참 기저 신호",
        "zoom": "극값 하나를 확대 (점선: 국소 2차 피팅)",
        "err": "기저 신호 추출 오차", "lpf_lab": "FT+LPF", "eem_lab": "EEM",
        "noise_ax": "노이즈 표준편차 (임의 단위)", "rms_ax": "기저 신호 RMS 오차",
        "cross": "교차점",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": {}, "legend": {},
        "grating": "Synthetic grating", "natural": "Natural image",
        "spec": "Amplitude spectrum (log)",
        "spec_zoom": "Amplitude spectrum (log, central ±64 bins)",
        "contaminated": "Image with periodic noise", "peaks": "Spectrum (circles: noise peaks)",
        "restored_on": "Notched: frequency on a bin",
        "restored_off": "Notched: frequency between bins",
        "ideal": "Ideal low-pass", "gauss": "Gaussian low-pass",
        "filt": "filter $H$", "psf_cut": "Spatial-domain cut", "result": "Filtered result",
        "raw_spec": "Spectrum without a window", "win_spec": "After a Hann window",
        "prof_ax": "Vertical cut, normalized to DC",
        "row_ax": "Vertical frequency (bin)",
        "prof_raw": "No window", "prof_win": "Hann window",
        "wavenumber": r"Wavenumber $\sigma$ (1/$\mu$m)", "intensity": "Intensity (a.u.)",
        "signal": "Spectral interference signal", "upper": "Upper envelope",
        "lower": "Lower envelope", "eem_bl": "EEM baseline", "true_bl": "True baseline",
        "zoom": "One extremum, zoomed (dashed: local quadratic fit)",
        "err": "Baseline extraction error", "lpf_lab": "FT+LPF", "eem_lab": "EEM",
        "noise_ax": "Noise standard deviation (a.u.)", "rms_ax": "Baseline RMS error",
        "cross": "crossover",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
rng = np.random.default_rng(5)
cam = camera().astype(float)
N = cam.shape[0]

# 그림1: 방향성 격자와 실제 영상의 스펙트럼
yy, xx = np.mgrid[0:N, 0:N].astype(float)
# 주파수를 푸리에 빈에 맞춘다(31, 13, 51 빈). 빈 사이에 끼면 4절의 누설 십자선이
# 1절에서 미리 그려져 점 네 개가 선의 격자에 묻힌다.
G_TILT = (31.0 / N, 13.0 / N)          # (fx, fy) 기울어진 줄무늬
G_HORIZ = 51.0 / N                     # fy 수평 줄무늬
grating = (128 + 60 * np.cos(2 * np.pi * (G_TILT[0] * xx + G_TILT[1] * yy))
           + 40 * np.cos(2 * np.pi * G_HORIZ * yy))
spec_grating, spec_cam = spectrum(grating), spectrum(cam)
# 점이 한 빈짜리라 그대로 보이면 안 보인다. 가운데 ±64빈만 잘라 보이고,
# 보기 쉽게 3x3 최댓값 필터로 점을 키운다(그림 표시용, 계산과 무관).
ZOOM = 64
spec_grating_zoom = maximum_filter(
    spec_grating[N // 2 - ZOOM:N // 2 + ZOOM + 1, N // 2 - ZOOM:N // 2 + ZOOM + 1], size=3)

# 그림2: 주기성 노이즈. 조명 리플이나 격자 반사가 실제로 이런 모양으로 들어온다.
# 주파수가 푸리에 빈에 정확히 떨어질 때와 아닐 때를 함께 본다. 노치가 잘 듣는지는
# 노이즈 세기가 아니라 주파수가 격자에 맞는지가 정한다.
NOTCH_R = 0.012
F_ON = (44.0 / N, 72.0 / N)            # 빈에 정확히 맞는 주파수
F_OFF = (0.085, 0.140)                 # 빈 사이에 낀 주파수
AMP_NOISE = 38.0
stripes_on = AMP_NOISE * np.cos(2 * np.pi * (F_ON[1] * xx + F_ON[0] * yy))
stripes_off = AMP_NOISE * np.cos(2 * np.pi * (F_OFF[1] * xx + F_OFF[0] * yy))
contaminated = cam + stripes_on
spec_contaminated = spectrum(contaminated)
restored_on, notch_mask = notch_filter(contaminated, [F_ON], radius=NOTCH_R)
restored_off, _ = notch_filter(cam + stripes_off, [F_OFF], radius=NOTCH_R)
collateral, _ = notch_filter(cam, [F_ON], radius=NOTCH_R)   # 깨끗한 영상에 같은 노치
rms_before = np.sqrt(np.mean(stripes_on ** 2))
rms_on = np.sqrt(np.mean((restored_on - cam) ** 2))
rms_off = np.sqrt(np.mean((restored_off - cam) ** 2))
rms_collateral = np.sqrt(np.mean((collateral - cam) ** 2))

# 그림3: 이상적 필터와 가우시안 필터
CUT = 0.06
ideal_img, H_ideal = lowpass(cam, CUT, "ideal")
gauss_img, H_gauss = lowpass(cam, CUT, "gaussian")
h_ideal, h_gauss = impulse_response(H_ideal), impulse_response(H_gauss)
cut_i, cut_g = h_ideal[N // 2], h_gauss[N // 2]
undershoot = abs(cut_i.min() / cut_i.max()) * 100.0

# 그림4: 주기 확장이 만드는 십자. 창을 곱하면 눌린다.
# 로그 눈금 그림만으로는 차이가 잘 안 보이므로, 스펙트럼 세로선의 단면을 함께 잰다.
# 두 스펙트럼의 전체 밝기가 다르므로 각자의 DC 성분으로 나눠 견준다.
win = hann2d(cam.shape)
spec_win = spectrum(cam, window=win)
CEN = N // 2


def _column_profile(img):
    F = np.abs(np.fft.fftshift(np.fft.fft2(img)))
    F = F / F[CEN, CEN]
    on_axis = F[:, CEN]                      # 십자선 위
    off_axis = 0.5 * (F[:, CEN - 3] + F[:, CEN + 3])   # 바로 옆
    return on_axis, off_axis


prof_on_raw, prof_off_raw = _column_profile(cam)
prof_on_win, prof_off_win = _column_profile(cam * win)
_sel = np.r_[10:CEN - 5, CEN + 5:N - 10]
cross_raw = float(np.median(prof_on_raw[_sel] / prof_off_raw[_sel]))
cross_win = float(np.median(prof_on_win[_sel] / prof_off_win[_sel]))
edge_gap = float(np.sqrt(np.mean((cam[0] - cam[-1]) ** 2)))

# 그림5~6: 분광 간섭 신호
DETECT_SIGMA, HALF = 6.0, 4
sig_x, sig_y, sig_bls, _ = make_signal()
eem_bl, (env_up, env_lo), maxima, minima = eem_baseline(
    sig_x, sig_y, half=HALF, detect_sigma=DETECT_SIGMA)
valid = ~np.isnan(eem_bl)
LPF_CUT = 17
lpf_bl = fourier_lpf_baseline(sig_y, LPF_CUT)
err_lpf, err_eem = lpf_bl - sig_bls, eem_bl - sig_bls

# 확대할 극값 하나와 그 국소 적합
zi = len(maxima) // 2
zx0 = maxima[zi, 0]
zmask = np.abs(sig_x - zx0) < 0.028
zc = int(np.argmin(np.abs(sig_x - zx0)))
zlo, zhi = zc - HALF, zc + HALF + 1
zpoly = np.polyfit(sig_x[zlo:zhi] - sig_x[zlo:zhi].mean(), sig_y[zlo:zhi], 2)
zfit_x = np.linspace(sig_x[zlo] - 0.0035, sig_x[zhi - 1] + 0.0035, 80)
zfit_y = np.polyval(zpoly, zfit_x - sig_x[zlo:zhi].mean())

# 노이즈를 올려가며 두 방법을 겨룬다
noise_levels = np.array([0.0, 0.002, 0.005, 0.010, 0.020, 0.035, 0.050, 0.070])
rms_lpf, rms_eem = [], []
for nz in noise_levels:
    x, y, b, _ = make_signal(noise=nz, seed=1)
    e, _, _, _ = eem_baseline(x, y, half=HALF, detect_sigma=DETECT_SIGMA)
    v = ~np.isnan(e)
    rms_eem.append(np.sqrt(np.mean((e[v] - b[v]) ** 2)))
    rms_lpf.append(min(np.sqrt(np.mean((fourier_lpf_baseline(y, c)[v] - b[v]) ** 2))
                       for c in range(8, 23)))
rms_lpf, rms_eem = np.array(rms_lpf), np.array(rms_eem)
cross_i = int(np.argmax(rms_eem > rms_lpf))


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    F, LG = L["font"], L["legend"]

    def caption(fig):
        # 그림 번호와 설명은 본문 캡션이 맡는다. 이미지 안에는 넣지 않는다.
        fig.tight_layout()

    # 그림1
    fig, ax = plt.subplots(2, 2, figsize=(7.4, 6.7))
    for r, (im, sp, name, st) in enumerate([
            (grating, spec_grating_zoom, L["grating"], L["spec_zoom"]),
            (cam, spec_cam, L["natural"], L["spec"])]):
        ax[r, 0].imshow(im, cmap="gray"); ax[r, 0].set_title(name, fontsize=10, **F)
        ax[r, 1].imshow(sp, cmap="inferno", interpolation="nearest")
        ax[r, 1].set_title(st, fontsize=10, **F)
        for a in ax[r]:
            a.set_xticks([]); a.set_yticks([])
    caption(fig)
    fig.savefig(os.path.join(out, "fig1-reading-the-spectrum.png"), dpi=150)
    plt.close(fig)

    # 그림2
    fig, ax = plt.subplots(1, 4, figsize=(11.2, 3.0))
    for a, (im, t, cm) in zip(ax, [(contaminated, L["contaminated"], "gray"),
                                   (spec_contaminated, L["peaks"], "inferno"),
                                   (restored_on, L["restored_on"], "gray"),
                                   (restored_off, L["restored_off"], "gray")]):
        a.imshow(im, cmap=cm); a.set_title(t, fontsize=9.5, **F)
        a.set_xticks([]); a.set_yticks([])
    # 노이즈 봉우리는 한 빈짜리 점이라 그대로는 안 보인다. 두 자리에 원을 그린다.
    for sgn in (1, -1):
        ax[1].add_patch(Circle((N // 2 + sgn * F_ON[1] * N, N // 2 + sgn * F_ON[0] * N),
                               14, fill=False, lw=1.4, color="#3ff0ff"))
    caption(fig)
    fig.savefig(os.path.join(out, "fig2-notch-filter.png"), dpi=150)
    plt.close(fig)

    # 그림3
    fig, ax = plt.subplots(2, 3, figsize=(9.6, 6.1))
    for r, (H, h, cut, img, name) in enumerate([
            (H_ideal, h_ideal, cut_i, ideal_img, L["ideal"]),
            (H_gauss, h_gauss, cut_g, gauss_img, L["gauss"])]):
        ax[r, 0].imshow(H, cmap="gray"); ax[r, 0].set_title(f"{name} {L['filt']}",
                                                            fontsize=9.5, **F)
        ax[r, 0].set_xticks([]); ax[r, 0].set_yticks([])
        c = N // 2
        ax[r, 1].plot(np.arange(-40, 41), cut[c - 40:c + 41], lw=1.4, color="#2b5c9b")
        ax[r, 1].axhline(0, color="k", lw=0.6)
        ax[r, 1].set_title(L["psf_cut"], fontsize=9.5, **F); ax[r, 1].grid(alpha=0.3)
        ax[r, 2].imshow(img, cmap="gray"); ax[r, 2].set_title(L["result"], fontsize=9.5, **F)
        ax[r, 2].set_xticks([]); ax[r, 2].set_yticks([])
    caption(fig)
    fig.savefig(os.path.join(out, "fig3-ideal-vs-gaussian.png"), dpi=150)
    plt.close(fig)

    # 그림4
    fig, ax = plt.subplots(1, 3, figsize=(10.2, 3.4))
    vlo, vhi = spec_cam.min(), spec_cam.max()
    for a, (im, t) in zip(ax[:2], [(spec_cam, L["raw_spec"]), (spec_win, L["win_spec"])]):
        a.imshow(im, cmap="inferno", vmin=vlo, vmax=vhi)
        a.set_title(t, fontsize=9.5, **F); a.set_xticks([]); a.set_yticks([])
    rows = np.arange(N) - CEN
    ax[2].semilogy(rows, prof_on_raw, lw=0.9, color="#c0392b", label=L["prof_raw"])
    ax[2].semilogy(rows, prof_on_win, lw=0.9, color="#2b5c9b", label=L["prof_win"])
    ax[2].set_xlabel(L["row_ax"], **F)
    ax[2].set_title(L["prof_ax"], fontsize=8.5, **F)
    ax[2].set_ylim(top=1e2)   # 위쪽 여백에 범례를 둬 DC 봉우리를 가리지 않는다
    ax[2].legend(prop=LG, fontsize=8.5, loc="upper right"); ax[2].grid(alpha=0.3, which="both")
    caption(fig)
    fig.savefig(os.path.join(out, "fig4-periodic-extension.png"), dpi=150)
    plt.close(fig)

    # 그림5
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.8),
                           gridspec_kw={"width_ratios": [2.1, 1]})
    ax[0].plot(sig_x, sig_y, lw=0.7, color="#888", label=L["signal"])
    ax[0].plot(sig_x, env_up, lw=1.3, color="#c0392b", label=L["upper"])
    ax[0].plot(sig_x, env_lo, lw=1.3, color="#2b5c9b", label=L["lower"])
    ax[0].plot(sig_x, eem_bl, lw=1.8, color="#1e8449", label=L["eem_bl"])
    ax[0].plot(sig_x, sig_bls, "--", lw=1.1, color="k", label=L["true_bl"])
    ax[0].plot(maxima[:, 0], maxima[:, 1], "x", ms=4, color="#c0392b")
    ax[0].plot(minima[:, 0], minima[:, 1], "x", ms=4, color="#2b5c9b")
    ax[0].set_xlabel(L["wavenumber"], **F); ax[0].set_ylabel(L["intensity"], **F)
    ax[0].set_ylim(top=1.6)   # 범례가 들어갈 위쪽 여백
    ax[0].legend(prop=LG, fontsize=8, ncol=2, loc="upper left"); ax[0].grid(alpha=0.3)
    ax[1].plot(sig_x[zmask], sig_y[zmask], "o", ms=3.5, color="#888")
    ax[1].plot(zfit_x, zfit_y, "--", lw=2.0, color="#c0392b", zorder=3)
    ax[1].plot(maxima[zi, 0], maxima[zi, 1], "x", ms=9, mew=2, color="#c0392b")
    ax[1].set_title(L["zoom"], fontsize=8.5, **F)
    ax[1].set_xlabel(L["wavenumber"], **F); ax[1].grid(alpha=0.3)
    caption(fig)
    fig.savefig(os.path.join(out, "fig5-extrema-envelope.png"), dpi=150)
    plt.close(fig)

    # 그림6
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 3.7))
    ax[0].plot(sig_x, err_lpf, lw=1.2, color="#c0392b", label=L["lpf_lab"])
    ax[0].plot(sig_x, err_eem, lw=1.2, color="#1e8449", label=L["eem_lab"])
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_xlabel(L["wavenumber"], **F); ax[0].set_ylabel(L["err"], **F)
    ax[0].legend(prop=LG, fontsize=9); ax[0].grid(alpha=0.3)
    ax[1].loglog(noise_levels[1:], rms_lpf[1:], "o-", ms=4, color="#c0392b",
                 label=L["lpf_lab"])
    ax[1].loglog(noise_levels[1:], rms_eem[1:], "s-", ms=4, color="#1e8449",
                 label=L["eem_lab"])
    ax[1].axvline(noise_levels[cross_i], ls=":", lw=1.2, color="#555")
    ax[1].text(noise_levels[cross_i] * 1.08, rms_eem[1] * 1.5, L["cross"],
               fontsize=8.5, color="#555", **F)
    ax[1].set_xlabel(L["noise_ax"], **F); ax[1].set_ylabel(L["rms_ax"], **F)
    ax[1].legend(prop=LG, fontsize=9); ax[1].grid(alpha=0.3, which="both")
    # 로그 눈금에서 10^-2 하나만 읽히지 않도록 1-2-5 눈금에 숫자를 단다.
    ax[1].set_xticks([0.002, 0.005, 0.01, 0.02, 0.05])
    ax[1].set_xticklabels(["0.002", "0.005", "0.01", "0.02", "0.05"])
    ax[1].xaxis.set_minor_formatter(plt.NullFormatter())
    caption(fig)
    fig.savefig(os.path.join(out, "fig6-lpf-vs-eem.png"), dpi=150)
    plt.close(fig)


for lang, L in LABELS.items():
    draw(L)
    print(f"[{lang}] 그림 6개 저장: {L['dir']}")

print()
print("=== 그림1: 합성 격자 ===")
print(f"  격자 주파수 (빈): 기울어진 줄무늬 ({G_TILT[0]*N:.0f}, {G_TILT[1]*N:.0f}), 수평 줄무늬 (0, {G_HORIZ*N:.0f})")
_Fg = np.abs(np.fft.fft2(grating)); _Fg[0, 0] = 0
_top = np.sort(_Fg.ravel())[::-1]
print(f"  교류 에너지 중 봉우리 4개의 몫 {(_top[:4]**2).sum()/(_Fg**2).sum():.3f}")
print()
print("=== 그림2: 노치 필터 ===")
print(f"  줄무늬 진폭 {AMP_NOISE:.0f} DN, 오염 RMS {rms_before:.2f}, 노치 반경 {NOTCH_R}")
print(f"  빈에 맞는 주파수 {F_ON[0]*N:.0f}/{N}, {F_ON[1]*N:.0f}/{N}")
print(f"    노치 후 RMS {rms_on:.2f}  ({rms_before/rms_on:.0f}배 개선)")
print(f"    같은 노치를 깨끗한 영상에 걸었을 때 손상 {rms_collateral:.2f} "
      f"-> 줄무늬는 사실상 완전히 사라졌다")
print(f"  빈 사이에 낀 주파수 {F_OFF[0]*N:.2f}/{N}, {F_OFF[1]*N:.2f}/{N}")
print(f"    노치 후 RMS {rms_off:.2f}  ({rms_before/rms_off:.1f}배 개선) "
      f"-> 누설 때문에 잔여가 손상량의 {rms_off/rms_collateral:.0f}배로 남는다")
print()
print("=== 그림3: 이상적 vs 가우시안 ===")
print(f"  차단주파수 {CUT} (나이키스트 0.5 기준)")
print(f"  이상적 임펄스 응답: 최대 {cut_i.max():.5f}, 첫 음의 골 {cut_i.min():.5f} "
      f"= {undershoot:.1f} %")
print(f"  가우시안 임펄스 응답 최소 {cut_g.min():.2e} (음수 골 없음)")
print(f"  결과 범위: 이상적 [{ideal_img.min():.1f}, {ideal_img.max():.1f}], "
      f"가우시안 [{gauss_img.min():.1f}, {gauss_img.max():.1f}], 원본 [0, 255]")
print()
print("=== 그림4: 주기 확장 ===")
print(f"  상하 끝 화소값 차이 RMS = {edge_gap:.1f} DN (이 불연속이 십자를 만든다)")
print(f"  십자선 대 주변 밝기 비: 창 없음 {cross_raw:.2f} -> 한 창 {cross_win:.2f}")
print()
print("=== 그림5: EEM ===")
print(f"  줄무늬 {RETARDER_DND*(sig_x.max()-sig_x.min()):.0f}개, "
      f"극대 {len(maxima)}개 극소 {len(minima)}개")
print(f"  유효 구간 {valid.sum()}/{len(sig_x)} 점")
print(f"  검출용 평활화 sigma = {DETECT_SIGMA} 화소 (줄무늬 한 주기 "
      f"{len(sig_x)/(RETARDER_DND*(sig_x.max()-sig_x.min())):.0f} 화소)")
print()
print("=== 그림6: 무잡음 비교 ===")
n = len(sig_x); edge = np.zeros(n, bool); edge[:n//20] = True; edge[-n//20:] = True
for name, e in (("FT+LPF", err_lpf), ("EEM", err_eem)):
    m = valid & np.isfinite(e)
    print(f"  {name:7s} 전구간 RMS {np.sqrt(np.mean(e[m]**2)):.5f}  "
          f"가장자리 5% {np.sqrt(np.mean(e[edge&m]**2)):.5f}  "
          f"중앙 {np.sqrt(np.mean(e[(~edge)&m]**2)):.5f}")
print()
print("=== 그림6: 노이즈 의존성 ===")
for nz, a, b in zip(noise_levels, rms_lpf, rms_eem):
    mark = "  <- 역전" if b > a else ""
    print(f"  노이즈 {nz:5.3f}: FT+LPF {a:.5f}  EEM {b:.5f}{mark}")

# ── 본문 보조 수치 (그림과 무관) ─────────────────────────────
print()
print("=== 4절: 노치 반경을 바꿔 보기 (빈 사이에 낀 줄무늬) ===")
stripes_off_only = cam + stripes_off
for r in (0.006, 0.012, 0.024, 0.048):
    rest, _ = notch_filter(stripes_off_only, [F_OFF], radius=r)
    coll, _ = notch_filter(cam, [F_OFF], radius=r)
    print(f"  반경 {r:.3f}: 잔여 RMS {np.sqrt(np.mean((rest - cam) ** 2)):.2f} DN, "
          f"노치 자신의 손상 {np.sqrt(np.mean((coll - cam) ** 2)):.2f} DN")

print()
print("=== 6절: 극값 정밀화를 빼면 ===")
import baseline_extraction as _be
_orig = _be._refine_extremum
_be._refine_extremum = lambda x_, y_, i, half=3: (x_[i], y_[i])   # 표본점을 그대로 극값으로
raw_bl, _, _, _ = _be.eem_baseline(sig_x, sig_y, half=HALF, detect_sigma=DETECT_SIGMA)
_be._refine_extremum = _orig
_m = valid & ~np.isnan(raw_bl)
print(f"  정밀화 없음 RMS {np.sqrt(np.mean((raw_bl[_m] - sig_bls[_m]) ** 2)):.5f}, "
      f"정밀화 RMS {np.sqrt(np.mean(err_eem[valid] ** 2)):.5f}")

print()
print("=== 7절: 교차점과 변조 진폭 ===")
_amp = make_signal()[3]
_d = np.log(rms_eem[1:]) - np.log(rms_lpf[1:])
_k = int(np.argmax(_d > 0))
_lx = np.log(noise_levels[1:])
_xc = np.exp(_lx[_k - 1] + (_lx[_k] - _lx[_k - 1]) * (-_d[_k - 1]) / (_d[_k] - _d[_k - 1]))
print(f"  교차 노이즈(로그 보간) {_xc:.4f}, 평균 진폭 {_amp.mean():.4f} -> {_xc/_amp.mean()*100:.1f} %")
print(f"  검출 평활화 sigma {DETECT_SIGMA} 화소가 56화소 주기 줄무늬 진폭에 곱하는 값 "
      f"{np.exp(-2 * np.pi**2 * DETECT_SIGMA**2 / 56.0**2):.2f}")

print()
print("=== 6절: 극값에서 I = I_BLS ± A 근사의 오차 ===")
import spectral_signal as _ss
_x = np.linspace(1.30, 2.50, 200001)            # 극값을 정확히 잡으려고 촘촘히 계산
_b, _a = _ss.baseline(_x), _ss.amplitude(_x)
_y = _b + _a * np.cos(2 * np.pi * _ss.RETARDER_DND * _x - _ss.phase_offset(_x))
_dy = np.gradient(_y)
_ix = np.flatnonzero(np.signbit(_dy[:-1]) != np.signbit(_dy[1:]))
_kind = np.where(_dy[_ix] > 0, 1, -1)
_apx = np.abs(_y[_ix] - (_b[_ix] + _kind * _a[_ix]))
print(f"  극값 {len(_ix)}개, 최대 오차 {_apx.max():.2e} (평균 진폭의 {_apx.max()/_a.mean()*100:.2f} %)")

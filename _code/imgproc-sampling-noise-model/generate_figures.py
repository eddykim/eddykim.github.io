"""이미지처리 1편 그림 6개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/imgproc-sampling-noise-model/      (한국어)
      ../../assets/img/posts/imgproc-sampling-noise-model/en/   (영문)

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은 데이터가
완전히 동일하고 표기만 다르다. 시리즈의 다른 편과 같은 방식이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy.ndimage import gaussian_filter
from skimage.data import camera

from sensor_model import Sensor, fit_ptc, photon_transfer_curve

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "imgproc-sampling-noise-model",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림1
        "chain": ["광자", "신호 전자", "전압", "DN"],
        "chain_ops": ["양자효율\n$\\eta$", "전하-전압\n변환", "ADC\n이득 $K$"],
        "noise_shot": "샷 노이즈\n(신호에 비례)",
        "noise_dark": "암전류\n+ 고정패턴",
        "noise_read": "읽기 노이즈\n(신호와 무관)",
        "noise_quant": "양자화\n노이즈",
        "fig1": "그림1. 광자에서 픽셀 값까지의 사슬과 각 단계에서 끼어드는 노이즈",
        # 그림2
        "mean_dn": "평균 신호 $\\mu_y - y_0$ (DN)",
        "var_dn": "분산 $\\sigma_y^2$ (DN$^2$)",
        "ptc_data": "측정 (두 장 차이로 구한 분산)",
        "ptc_fit": "선형 구간 직선 맞춤",
        "ptc_sat": "포화 시작",
        "ptc_zoom": "저신호 구간 확대 — 절편이 분해되는 영역",
        "fig2": "그림2. 광전달곡선 — 기울기가 이득 $K$, 절편이 어두운 영역 노이즈를 알려준다",
        # 그림3
        "bit_depth": "비트 심도 (bit)",
        "noise_e": "노이즈 (e$^-$ rms)",
        "quant_noise": "양자화 노이즈 $\\Delta/\\sqrt{12}$",
        "read_noise": "읽기 노이즈 (3 e$^-$)",
        "total_noise": "합성 노이즈",
        "bit_cross": "여기서부터\n양자화가 묻힌다",
        "fig3": "그림3. 비트 심도를 늘려도 읽기 노이즈 아래로는 내려가지 않는다",
        # 그림4
        "n_frames": "평균한 프레임 수 $N$",
        "resid_dn": "남은 공간적 산포 (DN)",
        "with_fpn": "고정패턴 있음",
        "without_fpn": "고정패턴 보정 후",
        "sqrt_law": r"$1/\sqrt{N}$ 기준선",
        "fpn_floor": "고정패턴 바닥",
        "fig4": "그림4. 프레임 평균은 고정패턴 노이즈 앞에서 멈춘다",
        # 그림5
        "zp": "존 플레이트", "cam": "실제 영상",
        "orig": "원본 (고해상)",
        "naive": "4배 나이브 샘플링",
        "prefilt": "저역통과 후 4배 샘플링",
        "fig5": "그림6. 나이키스트를 넘는 성분은 사라지지 않고 낮은 주파수로 접혀 들어온다",
        # 그림6
        "f_true": "실제 공간주파수 (나이키스트 단위)",
        "f_apparent": "보이는 공간주파수 (나이키스트 단위)",
        "ideal": "접힘이 없다면",
        "folded": "실제로 보이는 주파수",
        "nyq": "나이키스트",
        "fig6": "그림5. 겉보기 주파수의 접힘 — 계측에서 위험한 이유",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": {}, "legend": {},
        "chain": ["Photons", "Signal electrons", "Voltage", "DN"],
        "chain_ops": ["Quantum\nefficiency $\\eta$", "Charge-to-\nvoltage", "ADC\ngain $K$"],
        "noise_shot": "Shot noise\n(scales with signal)",
        "noise_dark": "Dark current\n+ fixed pattern",
        "noise_read": "Read noise\n(signal independent)",
        "noise_quant": "Quantization\nnoise",
        "fig1": "Fig 1. From photons to pixel values, and where each noise source enters",
        "mean_dn": "Mean signal $\\mu_y - y_0$ (DN)",
        "var_dn": "Variance $\\sigma_y^2$ (DN$^2$)",
        "ptc_data": "Measured (variance from frame difference)",
        "ptc_fit": "Linear-range fit",
        "ptc_sat": "Saturation onset",
        "ptc_zoom": "Low-signal zoom — where the intercept resolves",
        "fig2": "Fig 2. Photon transfer curve — slope gives gain $K$, intercept gives dark noise",
        "bit_depth": "Bit depth (bit)",
        "noise_e": "Noise (e$^-$ rms)",
        "quant_noise": "Quantization noise $\\Delta/\\sqrt{12}$",
        "read_noise": "Read noise (3 e$^-$)",
        "total_noise": "Combined noise",
        "bit_cross": "quantization is\nburied past here",
        "fig3": "Fig 3. More bits never take you below the read noise",
        "n_frames": "Number of averaged frames $N$",
        "resid_dn": "Residual spatial spread (DN)",
        "with_fpn": "With fixed pattern",
        "without_fpn": "After fixed-pattern correction",
        "sqrt_law": r"$1/\sqrt{N}$ reference",
        "fpn_floor": "Fixed-pattern floor",
        "fig4": "Fig 4. Frame averaging stops at the fixed-pattern floor",
        "zp": "Zone plate", "cam": "Natural image",
        "orig": "Original (high resolution)",
        "naive": "Naive 4x sampling",
        "prefilt": "Low-pass filtered, then 4x",
        "fig5": "Fig 6. Content above Nyquist does not vanish — it folds down in frequency",
        "f_true": "True spatial frequency (Nyquist units)",
        "f_apparent": "Apparent spatial frequency (Nyquist units)",
        "ideal": "If there were no folding",
        "folded": "What is actually seen",
        "nyq": "Nyquist",
        "fig6": "Fig 5. Frequency folding — why aliasing is dangerous in metrology",
    },
}

# ─────────────────────────────────────────────────────────────
# 계산 — 언어와 무관하게 한 번만
# ─────────────────────────────────────────────────────────────
rng = np.random.default_rng(1)
sensor = Sensor(shape=(192, 192), seed=7)

# 그림2: 광전달곡선. 포화 직전까지 노출을 훑는다.
photon_levels = np.concatenate([
    np.linspace(0, 1200, 16), np.linspace(1600, 30000, 20),
])
EXPOSURE_S = 0.01
# 분산 자체가 통계량이라 추정값이 흔들린다. 쌍을 여러 번 찍어 평균해야 절편이
# 분해된다 — 실제 카메라 특성 측정에서도 레벨당 수십 장을 찍는 이유다.
ptc_mean, ptc_var = photon_transfer_curve(
    sensor, photon_levels, EXPOSURE_S, rng, n_pairs=24)
ptc_signal = ptc_mean - sensor.offset_dn
FIT_LO, FIT_HI = 0.0, 200.0
SIGMA_Q_SQ = 1.0 / 12.0   # 양자화 노이즈 분산 [DN^2], 균등 분포 가정
K_hat, sigma_read_hat, slope, intercept = fit_ptc(
    ptc_signal, ptc_var, (FIT_LO, FIT_HI), sigma_q_sq=SIGMA_Q_SQ)
K_raw, sigma_read_raw, _, _ = fit_ptc(ptc_signal, ptc_var, (FIT_LO, FIT_HI))
sat_idx = int(np.argmax(ptc_var))
sat_signal = ptc_signal[sat_idx]

# 그림3: 비트 심도. 양자화 간격은 최대 전하 용량을 2^B 로 나눈 값이다.
bits = np.arange(6, 17)
quant_step_e = sensor.full_well_e / (2.0 ** bits)
quant_noise_e = quant_step_e / np.sqrt(12.0)
read_e = sensor.read_noise_e
total_noise_e = np.sqrt(read_e ** 2 + quant_noise_e ** 2)
bit_cross = np.log2(sensor.full_well_e / (np.sqrt(12.0) * read_e))

# 그림4: 프레임 평균. 밝기를 고정하고 N 장 평균의 공간적 산포를 본다.
PHOTONS = 8000.0
n_list = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256, 512])
avg_rng = np.random.default_rng(11)
resid_fpn, resid_clean = [], []
frames_fpn = [sensor.capture(PHOTONS, EXPOSURE_S, avg_rng, fpn=True)
              for _ in range(int(n_list[-1]))]
frames_clean = [sensor.capture(PHOTONS, EXPOSURE_S, avg_rng, fpn=False)
                for _ in range(int(n_list[-1]))]
for n in n_list:
    resid_fpn.append(np.mean(frames_fpn[:n], axis=0).std())
    resid_clean.append(np.mean(frames_clean[:n], axis=0).std())
resid_fpn = np.array(resid_fpn)
resid_clean = np.array(resid_clean)
sqrt_line = resid_clean[0] / np.sqrt(n_list)
fpn_floor = resid_fpn[-1]

# 그림5: 에일리어싱. 존 플레이트는 중심에서 멀어질수록 주파수가 올라가므로
# 나이키스트를 넘는 지점이 눈에 그대로 보인다.
N_HI, FACTOR = 512, 4
yy, xx = np.mgrid[0:N_HI, 0:N_HI].astype(float)
r2 = ((xx - N_HI / 2) ** 2 + (yy - N_HI / 2) ** 2) / (N_HI / 2) ** 2
zone_plate = 0.5 + 0.5 * np.cos(np.pi * 44.0 * r2)
cam_hi = camera().astype(float) / 255.0

def naive_sample(img, factor):
    """픽셀을 건너뛰며 고르는 것 — 실제로 픽셀이 커진 상황과 같다."""
    return img[::factor, ::factor]

def prefiltered_sample(img, factor):
    """나이키스트 위를 먼저 지우고 고른다. 잃을 것은 잃되 접히지는 않는다."""
    return naive_sample(gaussian_filter(img, sigma=0.5 * factor), factor)

zp_naive, zp_pre = naive_sample(zone_plate, FACTOR), prefiltered_sample(zone_plate, FACTOR)
cam_naive, cam_pre = naive_sample(cam_hi, FACTOR), prefiltered_sample(cam_hi, FACTOR)

# 그림6: 겉보기 주파수 접힘. 나이키스트를 넘으면 삼각파 모양으로 되접힌다.
f_true = np.linspace(0, 4, 800)
f_apparent = np.abs(f_true - 2.0 * np.round(f_true / 2.0))


# ─────────────────────────────────────────────────────────────
# 그리기
# ─────────────────────────────────────────────────────────────
def draw(labels):
    out = labels["dir"]
    os.makedirs(out, exist_ok=True)
    F, L = labels["font"], labels["legend"]

    def caption(ax_or_plt, text):
        ax_or_plt.suptitle("") if False else None

    # ── 그림1: 사슬 개념도
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    ax.set_xlim(0, 10.6); ax.set_ylim(0, 4.2); ax.axis("off")
    box_x = [0.3, 3.0, 5.7, 8.4]
    for x, name in zip(box_x, labels["chain"]):
        ax.add_patch(FancyBboxPatch((x, 1.7), 1.9, 0.9, boxstyle="round,pad=0.08",
                                    fc="#eef3fa", ec="#2b5c9b", lw=1.4))
        ax.text(x + 0.95, 2.15, name, ha="center", va="center", fontsize=11, **F)
    for i, op in enumerate(labels["chain_ops"]):
        x0, x1 = box_x[i] + 1.9, box_x[i + 1]
        ax.add_patch(FancyArrowPatch((x0, 2.15), (x1, 2.15), arrowstyle="-|>",
                                     mutation_scale=15, lw=1.3, color="#2b5c9b"))
        ax.text((x0 + x1) / 2, 2.62, op, ha="center", va="bottom", fontsize=8.5, **F)
    noises = [(1.25, labels["noise_shot"], "#c0392b"),
              (1.25, labels["noise_dark"], "#c0392b"),
              (6.65, labels["noise_read"], "#7d3c98"),
              (9.35, labels["noise_quant"], "#7d3c98")]
    for (x, txt, col), y0 in zip(noises, [0.35, 0.35, 0.35, 0.35]):
        pass
    for x, txt, col, dx in [(1.25, labels["noise_shot"], "#c0392b", -0.9),
                            (1.25, labels["noise_dark"], "#c0392b", 0.95),
                            (6.65, labels["noise_read"], "#7d3c98", 0.0),
                            (9.35, labels["noise_quant"], "#7d3c98", 0.0)]:
        ax.add_patch(FancyArrowPatch((x + dx, 1.05), (x + dx * 0.25, 1.65),
                                     arrowstyle="-|>", mutation_scale=12,
                                     lw=1.1, color=col, ls="--"))
        ax.text(x + dx, 0.95, txt, ha="center", va="top", fontsize=8, color=col, **F)
    fig.text(0.5, 0.015, labels["fig1"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig1-photon-to-dn-chain.png"), dpi=150)
    plt.close(fig)

    # ── 그림2: 광전달곡선
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(ptc_signal, ptc_var, "o", ms=4, color="#2b5c9b", label=labels["ptc_data"])
    xs = np.linspace(0, sat_signal * 1.02, 50)
    ax.plot(xs, slope * xs + intercept, "-", lw=1.6,
            color="#c0392b", label=labels["ptc_fit"])
    ax.axvline(sat_signal, ls=":", lw=1.2, color="#555")
    ax.text(sat_signal, ax.get_ylim()[1] * 0.55, "  " + labels["ptc_sat"],
            fontsize=8.5, color="#555", **F)
    ax.set_xlabel(labels["mean_dn"], **F)
    ax.set_ylabel(labels["var_dn"], **F)
    ax.legend(prop=L, fontsize=9, loc="upper left")
    ax.grid(alpha=0.3)
    # 절편은 기울기 항보다 세 자릿수 작다. 확대해야 비로소 보인다.
    axi = ax.inset_axes([0.55, 0.13, 0.42, 0.34])
    mz = ptc_signal <= 260
    axi.plot(ptc_signal[mz], ptc_var[mz], "o", ms=3, color="#2b5c9b")
    xz = np.linspace(0, 260, 30)
    axi.plot(xz, slope * xz + intercept, "-", lw=1.3, color="#c0392b")
    axi.axhline(intercept, ls=":", lw=1.0, color="#555")
    axi.set_title(labels["ptc_zoom"], fontsize=8, **F)
    axi.tick_params(labelsize=7)
    axi.grid(alpha=0.3)
    fig.text(0.5, 0.015, labels["fig2"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig2-photon-transfer-curve.png"), dpi=150)
    plt.close(fig)

    # ── 그림3: 비트 심도
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogy(bits, quant_noise_e, "o-", ms=4, color="#2b5c9b",
                label=labels["quant_noise"])
    ax.axhline(read_e, ls="--", lw=1.4, color="#c0392b", label=labels["read_noise"])
    ax.semilogy(bits, total_noise_e, "s-", ms=4, color="#1e8449",
                label=labels["total_noise"])
    ax.axvline(bit_cross, ls=":", lw=1.2, color="#555")
    ax.annotate(labels["bit_cross"], xy=(bit_cross, read_e * 4),
                xytext=(bit_cross + 1.2, read_e * 9), fontsize=8.5, color="#555",
                arrowprops=dict(arrowstyle="->", color="#555", lw=1.0), **F)
    ax.set_xlabel(labels["bit_depth"], **F)
    ax.set_ylabel(labels["noise_e"], **F)
    ax.legend(prop=L, fontsize=9)
    ax.grid(alpha=0.3, which="both")
    fig.text(0.5, 0.015, labels["fig3"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig3-bit-depth-vs-noise.png"), dpi=150)
    plt.close(fig)

    # ── 그림4: 프레임 평균
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.loglog(n_list, resid_fpn, "o-", ms=4, color="#c0392b", label=labels["with_fpn"])
    ax.loglog(n_list, resid_clean, "s-", ms=4, color="#2b5c9b",
              label=labels["without_fpn"])
    ax.loglog(n_list, sqrt_line, ":", lw=1.3, color="#555", label=labels["sqrt_law"])
    ax.axhline(fpn_floor, ls="--", lw=1.0, color="#c0392b", alpha=0.5)
    ax.text(n_list[-1], fpn_floor * 1.06, labels["fpn_floor"], fontsize=8.5,
            color="#c0392b", ha="right", va="bottom", **F)
    ax.set_xlabel(labels["n_frames"], **F)
    ax.set_ylabel(labels["resid_dn"], **F)
    ax.legend(prop=L, fontsize=9)
    ax.grid(alpha=0.3, which="both")
    fig.text(0.5, 0.015, labels["fig4"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig4-frame-averaging-fpn-floor.png"), dpi=150)
    plt.close(fig)

    # ── 그림5: 에일리어싱
    fig, axes = plt.subplots(2, 3, figsize=(8.4, 6.3))
    panels = [
        (zone_plate, labels["orig"], labels["zp"]),
        (zp_naive, labels["naive"], None),
        (zp_pre, labels["prefilt"], None),
        (cam_hi, labels["orig"], labels["cam"]),
        (cam_naive, labels["naive"], None),
        (cam_pre, labels["prefilt"], None),
    ]
    for ax, (img, title, rowlab) in zip(axes.ravel(), panels):
        ax.imshow(img, cmap="gray", interpolation="nearest")
        ax.set_title(title, fontsize=9.5, **F)
        ax.set_xticks([]); ax.set_yticks([])
        if rowlab:
            ax.set_ylabel(rowlab, fontsize=10, **F)
    fig.text(0.5, 0.015, labels["fig5"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.045, 1, 1), h_pad=2.2)
    fig.savefig(os.path.join(out, "fig6-aliasing.png"), dpi=150)
    plt.close(fig)

    # ── 그림6: 주파수 접힘
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(f_true, f_true, ":", lw=1.3, color="#555", label=labels["ideal"])
    ax.plot(f_true, f_apparent, "-", lw=2.0, color="#c0392b", label=labels["folded"])
    ax.axvline(1.0, ls="--", lw=1.1, color="#2b5c9b")
    ax.text(1.04, 3.4, labels["nyq"], fontsize=9, color="#2b5c9b", **F)
    ax.set_xlabel(labels["f_true"], **F)
    ax.set_ylabel(labels["f_apparent"], **F)
    ax.set_ylim(0, 4.1)
    ax.legend(prop=L, fontsize=9, loc="upper right")
    ax.grid(alpha=0.3)
    fig.text(0.5, 0.015, labels["fig6"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig5-frequency-folding.png"), dpi=150)
    plt.close(fig)


for lang, labels in LABELS.items():
    draw(labels)
    print(f"[{lang}] 그림 6개 저장: {labels['dir']}")

# ─────────────────────────────────────────────────────────────
# 본문에 인용할 수치
# ─────────────────────────────────────────────────────────────
print()
print("=== 그림2: 광전달곡선으로 되찾은 센서 상수 ===")
print(f"  참값   K = {sensor.gain_K:.4f} DN/e-, 읽기 노이즈 = {sensor.read_noise_e:.2f} e-")
print(f"  맞춤 절편 = {intercept:.4f} DN^2 "
      f"(이론: K^2*sigma_d^2 = {sensor.gain_K**2*(sensor.read_noise_e**2 + sensor.dark_rate_e*EXPOSURE_S):.4f} "
      f"+ sigma_q^2 = {SIGMA_Q_SQ:.4f})")
print(f"  추정값 K = {K_hat:.4f} DN/e-")
print(f"    절편을 그대로 환산  = {sigma_read_raw:.2f} e-")
print(f"    양자화분을 빼고 환산 = {sigma_read_hat:.2f} e-")
print(f"  (암전류 기여 포함 이론값 = "
      f"{np.sqrt(sensor.read_noise_e**2 + sensor.dark_rate_e*EXPOSURE_S):.2f} e-)")
print(f"  포화 시작 신호 = {sat_signal:.0f} DN "
      f"(= {sat_signal/sensor.gain_K:.0f} e-, 최대 전하 용량 {sensor.full_well_e} e-)")
print()
print("=== 그림3: 비트 심도 ===")
for b, q, t in zip(bits, quant_noise_e, total_noise_e):
    print(f"  {b:2d} bit: 양자화 {q:8.3f} e-, 합성 {t:6.3f} e-")
print(f"  양자화 = 읽기 노이즈가 되는 지점: {bit_cross:.1f} bit")
print()
print("=== 그림4: 프레임 평균 ===")
for n, a, b in zip(n_list, resid_fpn, resid_clean):
    print(f"  N={n:4d}: 고정패턴 있음 {a:7.3f} DN, 보정 후 {b:7.3f} DN")
print(f"  고정패턴 바닥 = {fpn_floor:.3f} DN "
      f"(신호 {sensor.gain_K*sensor.qe*PHOTONS:.0f} DN 의 "
      f"{100*fpn_floor/(sensor.gain_K*sensor.qe*PHOTONS):.2f} %)")
print()
print("=== 그림5: 에일리어싱 ===")
print(f"  존 플레이트 {N_HI}x{N_HI} -> {zp_naive.shape[0]}x{zp_naive.shape[1]} ({FACTOR}배)")
print(f"  저역통과 가우시안 sigma = {0.5*FACTOR:.1f} 픽셀")

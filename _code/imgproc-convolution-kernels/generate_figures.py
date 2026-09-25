"""이미지처리 2편 그림 6개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/imgproc-convolution-kernels/      (한국어)
      ../../assets/img/posts/imgproc-convolution-kernels/en/   (영문)

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼운다. 시리즈의 다른 편과 같은 방식이다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, median_filter, uniform_filter
from skimage.data import camera

from kernels import (CENTRAL_DIFF, SCHARR_SMOOTH, SOBEL_SMOOTH, box_kernel,
                     convolve1d_edge, frequency_response, gaussian_kernel,
                     noise_gain, separable_cost)
from edge_experiment import make_edge, sweep

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                        "assets", "img", "posts", "imgproc-convolution-kernels")

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "tap": "커널 탭 위치 (픽셀)", "weight": "가중치",
        "freq": "공간주파수 (나이키스트 = 1)", "resp": "응답 크기",
        "box": "박스 (폭 9)", "gauss": r"가우시안 ($\sigma$ = 2)",
        "shape": "커널 모양", "response": "주파수 응답",
        "fig1": "그림2. 박스 커널은 고주파를 깔끔하게 자르지 못하고 부엽을 남긴다",

        "sigma_ax": r"평활화 폭 $\sigma$ (픽셀)", "ratio": r"$\sigma_{out}/\sigma_{in}$",
        "measured": "측정", "theory": r"이론 $\sum w^2$ (2차원)",
        "gauss_m": "가우시안 (측정)", "box_m": "같은 폭 박스 (측정)",
        "fig2": "그림1. 평활화가 노이즈를 줄이는 정도는 커널 제곱합이 정한다",

        "ideal_d": "이상적 미분", "central": "중심차분", "sobel": "Sobel 평활 결합",
        "dog": r"DoG ($\sigma$ = 2)",
        "fig3": "그림3. 미분 커널의 응답은 주파수와 함께 커진다 — 평활화와 묶어야 하는 이유",

        "pos_std": "엣지 위치 표준편차 (픽셀)", "width_ax": "엣지 폭 10-90% (픽셀)",
        "std_lab": "위치 산포", "width_lab": "엣지 폭", "opt": "산포 최솟값",
        "fig4": "그림4. 평활화의 교환비 — 위치는 안정되고 엣지는 뭉개진다",

        "pos": "위치 (픽셀)", "err": "참값과의 차이 (DN)",
        "zero": "0으로 채움", "replicate": "끝값 반복", "reflect": "거울 반사",
        "fig5": "그림5. 경계 처리 방식이 가장자리에 남기는 오차 (오른쪽은 세로축 확대)",

        "orig_sp": "스파이크 섞인 원본", "gauss_sp": r"가우시안 ($\sigma$ = 1.5)",
        "med_sp": "메디안 (3x3)",
        "fig6": "그림6. 선형 필터는 스파이크를 지우지 못하고 번지게 한다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": {}, "legend": {},
        "tap": "Kernel tap position (pixel)", "weight": "Weight",
        "freq": "Spatial frequency (Nyquist = 1)", "resp": "Response magnitude",
        "box": "Box (width 9)", "gauss": r"Gaussian ($\sigma$ = 2)",
        "shape": "Kernel shape", "response": "Frequency response",
        "fig1": "Fig 2. A box kernel cuts high frequencies poorly, leaving sidelobes",

        "sigma_ax": r"Smoothing width $\sigma$ (pixel)", "ratio": r"$\sigma_{out}/\sigma_{in}$",
        "measured": "Measured", "theory": r"Theory $\sum w^2$ (2-D)",
        "gauss_m": "Gaussian (measured)", "box_m": "Box of equal width (measured)",
        "fig2": "Fig 1. Noise reduction is set by the sum of squared kernel weights",

        "ideal_d": "Ideal derivative", "central": "Central difference",
        "sobel": "Sobel (with smoothing)", "dog": r"DoG ($\sigma$ = 2)",
        "fig3": "Fig 3. Derivative response grows with frequency — hence the pairing with smoothing",

        "pos_std": "Std of edge position (pixel)", "width_ax": "Edge width 10-90% (pixel)",
        "std_lab": "Position scatter", "width_lab": "Edge width", "opt": "Scatter minimum",
        "fig4": "Fig 4. The smoothing trade — position stabilizes, the edge blurs",

        "pos": "Position (pixel)", "err": "Deviation from truth (DN)",
        "zero": "Zero padding", "replicate": "Replicate", "reflect": "Reflect",
        "fig5": "Fig 5. Error left at the borders by each padding mode (right panel rescaled)",

        "orig_sp": "Original with spikes", "gauss_sp": r"Gaussian ($\sigma$ = 1.5)",
        "med_sp": "Median (3x3)",
        "fig6": "Fig 6. A linear filter smears a spike instead of removing it",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
rng = np.random.default_rng(3)

# 그림1: 커널 모양과 주파수 응답
BOX_K, GAUSS_S = 9, 2.0
w_box, w_gauss = box_kernel(BOX_K), gaussian_kernel(GAUSS_S)
f_ax, H_box = frequency_response(w_box)
_, H_gauss = frequency_response(w_gauss)

# 그림2: 노이즈 이득. 백색 노이즈를 실제로 걸러 이론값과 대조한다.
# 2차원 가우시안은 1차원 커널을 두 방향에 거는 것과 같으므로, 이득도 1차원 이득의
# 제곱이다. 1차원 이득이 sqrt(sum w^2) 이니 2차원 이득은 sum w^2 이 된다.
sig_list = np.array([0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
white = rng.standard_normal((768, 768))
meas_g = np.array([gaussian_filter(white, s)[24:-24, 24:-24].std() for s in sig_list])
theo_g = np.array([noise_gain(gaussian_kernel(s)) for s in sig_list])
# 같은 분산을 갖는 박스 폭은 k = sqrt(12)*sigma 다. 홀수 폭으로 반올림되므로
# k 가 5 이상으로 커지는 구간만 비교한다.
box_widths = np.array([int(round(np.sqrt(12.0) * s)) | 1 for s in sig_list])
box_mask = box_widths >= 5
meas_b = np.array([uniform_filter(white, int(k))[24:-24, 24:-24].std()
                   for k in box_widths[box_mask]])

# 그림3: 미분 커널의 주파수 응답
f3, H_cd = frequency_response(CENTRAL_DIFF)
H_ideal = np.pi * f3                                  # 이상적 미분의 크기
_, H_sobel = frequency_response(np.convolve(CENTRAL_DIFF, SOBEL_SMOOTH))
w_dog = np.gradient(gaussian_kernel(GAUSS_S))
_, H_dog = frequency_response(w_dog)

# 그림4: 핵심 실험
sweep_sigmas = np.array([0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5,
                         3.0, 4.0, 5.0, 6.0, 8.0, 10.0])
sw_bias, sw_std, sw_width = sweep(sweep_sigmas, 400, np.random.default_rng(0))
opt_i = int(np.argmin(sw_std))

# 그림5: 경계 처리. 밝은 평지 위의 완만한 경사 — 가장자리 오차가 드러난다.
# 신호는 배열 바깥으로도 실제로 이어진다. 긴 신호를 먼저 만들고 가운데 창만
# "촬영된 배열"로 삼아야, 각 경계 처리가 바깥을 얼마나 잘못 짐작하는지 드러난다.
n5, PAD5 = 96, 64
t_long = np.arange(-PAD5, n5 + PAD5)
def scene(t):
    return 800.0 + 300.0 * np.sin(2.0 * np.pi * t / 70.0) + 1.6 * t
base_long = scene(t_long)
base = scene(np.arange(n5))          # 센서가 실제로 담아낸 구간
w5 = gaussian_kernel(3.0)
edge_modes = {m: convolve1d_edge(base, w5, m) for m in ("zero", "replicate", "reflect")}
truth = np.convolve(base_long, w5[::-1], mode="same")[PAD5:PAD5 + n5]

# 그림6: 스파이크 픽셀
img = camera().astype(float)[120:280, 200:360]
spiky = img.copy()
ys, xs = rng.integers(0, img.shape[0], 260), rng.integers(0, img.shape[1], 260)
spiky[ys[:130], xs[:130]] = 255.0
spiky[ys[130:], xs[130:]] = 0.0
sp_gauss = gaussian_filter(spiky, 1.5)
sp_med = median_filter(spiky, size=3)


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    F, LG = L["font"], L["legend"]

    # 그림1
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.6))
    xb = np.arange(len(w_box)) - len(w_box) // 2
    xg = np.arange(len(w_gauss)) - len(w_gauss) // 2
    ax[0].stem(xb, w_box, linefmt="C0-", markerfmt="C0o", basefmt=" ", label=L["box"])
    ax[0].stem(xg, w_gauss, linefmt="C3-", markerfmt="C3s", basefmt=" ", label=L["gauss"])
    ax[0].set_xlabel(L["tap"], **F); ax[0].set_ylabel(L["weight"], **F)
    ax[0].set_title(L["shape"], fontsize=10, **F)
    ax[0].legend(prop=LG, fontsize=8.5); ax[0].grid(alpha=0.3)
    ax[1].plot(f_ax, np.abs(H_box), "C0-", lw=1.6, label=L["box"])
    ax[1].plot(f_ax, np.abs(H_gauss), "C3-", lw=1.6, label=L["gauss"])
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel(L["freq"], **F); ax[1].set_ylabel(L["resp"], **F)
    ax[1].set_title(L["response"], fontsize=10, **F)
    ax[1].legend(prop=LG, fontsize=8.5); ax[1].grid(alpha=0.3)
    fig.text(0.5, 0.015, L["fig1"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig2-kernel-shape-response.png"), dpi=150)
    plt.close(fig)

    # 그림2
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(sig_list, theo_g, "-", lw=1.6, color="#c0392b", label=L["theory"])
    ax.plot(sig_list, meas_g, "o", ms=5, color="#2b5c9b", label=L["gauss_m"])
    ax.plot(sig_list[box_mask], meas_b, "s", ms=5, color="#1e8449", label=L["box_m"])
    ax.set_xlabel(L["sigma_ax"], **F); ax.set_ylabel(L["ratio"], **F)
    ax.set_yscale("log"); ax.legend(prop=LG, fontsize=9); ax.grid(alpha=0.3, which="both")
    fig.text(0.5, 0.015, L["fig2"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig1-noise-gain.png"), dpi=150)
    plt.close(fig)

    # 그림3
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(f3, H_ideal, ":", lw=1.4, color="#555", label=L["ideal_d"])
    ax.plot(f3, np.abs(H_cd), "-", lw=1.7, color="#c0392b", label=L["central"])
    ax.plot(f3, np.abs(H_sobel), "-", lw=1.5, color="#1e8449", label=L["sobel"])
    ax.plot(f3, np.abs(H_dog), "-", lw=1.5, color="#2b5c9b", label=L["dog"])
    ax.set_xlabel(L["freq"], **F); ax.set_ylabel(L["resp"], **F)
    ax.legend(prop=LG, fontsize=9); ax.grid(alpha=0.3)
    fig.text(0.5, 0.015, L["fig3"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig3-derivative-kernels.png"), dpi=150)
    plt.close(fig)

    # 그림4
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.plot(sweep_sigmas, sw_std, "o-", ms=4, color="#c0392b", label=L["std_lab"])
    ax.plot(sweep_sigmas[opt_i], sw_std[opt_i], "*", ms=15, color="#c0392b")
    ax.annotate(L["opt"], xy=(sweep_sigmas[opt_i], sw_std[opt_i]),
                xytext=(sweep_sigmas[opt_i] + 0.9, sw_std[opt_i] - 0.0035),
                fontsize=8.5, color="#c0392b", va="top",
                arrowprops=dict(arrowstyle="->", color="#c0392b", lw=1.0), **F)
    ax.set_ylim(sw_std.min() - 0.008, sw_std.max() * 1.06)
    ax.set_xlabel(L["sigma_ax"], **F)
    ax.set_ylabel(L["pos_std"], color="#c0392b", **F)
    ax.tick_params(axis="y", colors="#c0392b")
    ax.grid(alpha=0.3)
    ax2 = ax.twinx()
    ax2.plot(sweep_sigmas, sw_width, "s--", ms=4, color="#2b5c9b", label=L["width_lab"])
    ax2.set_ylabel(L["width_ax"], color="#2b5c9b", **F)
    ax2.tick_params(axis="y", colors="#2b5c9b")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, prop=LG, fontsize=9, loc="upper right")
    fig.text(0.5, 0.015, L["fig4"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig4-edge-tradeoff.png"), dpi=150)
    plt.close(fig)

    # 그림5
    fig, axp = plt.subplots(1, 2, figsize=(9.0, 3.8))
    cols = {"zero": "#c0392b", "replicate": "#e67e22", "reflect": "#2b5c9b"}
    for m in ("zero", "replicate", "reflect"):
        axp[0].plot(edge_modes[m] - truth, lw=1.5, color=cols[m], label=L[m])
    for m in ("replicate", "reflect"):
        axp[1].plot(edge_modes[m] - truth, lw=1.5, color=cols[m], label=L[m])
    for a in axp:
        a.axhline(0, color="k", lw=0.6)
        a.set_xlabel(L["pos"], **F); a.set_ylabel(L["err"], **F)
        a.legend(prop=LG, fontsize=9); a.grid(alpha=0.3)
    ax = axp[0]
    fig.text(0.5, 0.015, L["fig5"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig5-boundary-modes.png"), dpi=150)
    plt.close(fig)

    # 그림6
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 3.2))
    for a, (im, t) in zip(axes, [(spiky, L["orig_sp"]), (sp_gauss, L["gauss_sp"]),
                                 (sp_med, L["med_sp"])]):
        a.imshow(im, cmap="gray", vmin=0, vmax=255, interpolation="nearest")
        a.set_title(t, fontsize=9.5, **F); a.set_xticks([]); a.set_yticks([])
    fig.text(0.5, 0.015, L["fig6"], ha="center", fontsize=9, **F)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig6-median-vs-gaussian.png"), dpi=150)
    plt.close(fig)


for lang, L in LABELS.items():
    draw(L)
    print(f"[{lang}] 그림 6개 저장: {L['dir']}")

print()
print("=== 그림1: 커널 ===")
print(f"  박스(폭 {BOX_K}) 최대 부엽 = {np.abs(H_box)[np.argmax(f_ax > 0.25):].max():.3f}")
print(f"  가우시안(sigma={GAUSS_S}) 나이키스트 응답 = {np.abs(H_gauss)[-1]:.2e}")
print()
print("=== 그림2: 노이즈 이득 ===")
for s, t, g in zip(sig_list, theo_g, meas_g):
    print(f"  sigma={s:4.2f}: 이론 {t:.4f}  가우시안(측정) {g:.4f}")
for k, b in zip(box_widths[box_mask], meas_b):
    print(f"  같은 분산 박스 폭 {k:2d}: {b:.4f}")
print()
print("=== 그림3: 미분 커널 ===")
print(f"  중심차분: 나이키스트 응답 {np.abs(H_cd)[-1]:.3f} (이상적 {H_ideal[-1]:.3f})")
print(f"  중심차분 노이즈 이득 sqrt = {np.sqrt(noise_gain(CENTRAL_DIFF)):.4f}")
print(f"  DoG(sigma=2) 노이즈 이득 sqrt = {np.sqrt(noise_gain(w_dog)):.4f}")
print(f"  응답 최대 지점: 중심차분 f={f3[np.argmax(np.abs(H_cd))]:.2f}, "
      f"DoG f={f3[np.argmax(np.abs(H_dog))]:.2f}")
print()
print("=== 그림4: 교환비 ===")
for s, b, sd, w in zip(sweep_sigmas, sw_bias, sw_std, sw_width):
    print(f"  sigma={s:5.2f}: 편향 {b:+.4f}  산포 {sd:.4f}  엣지폭 {w:6.2f}")
print(f"  산포 최솟값: sigma={sweep_sigmas[opt_i]:.2f}, {sw_std[opt_i]:.4f} px, "
      f"엣지폭 {sw_width[opt_i]:.2f} px (평활화 없을 때 {sw_std[0]:.4f} px)")
print(f"  개선 배수 = {sw_std[0]/sw_std[opt_i]:.2f}")
print()
print("=== 그림5: 경계 ===")
for m in ("zero", "replicate", "reflect"):
    e = edge_modes[m] - truth
    print(f"  {m:10s}: 최대 오차 {np.abs(e).max():7.2f} DN, "
          f"가장자리 5화소 RMS {np.sqrt((e[:5]**2).mean()):6.2f} DN")
print()
print("=== 그림6: 스파이크 ===")
print(f"  스파이크 화소 비율 = {260/img.size*100:.2f} %")
for name, out in (("스파이크 원본", spiky), ("가우시안", sp_gauss), ("메디안", sp_med)):
    print(f"  {name:12s}: 원본 대비 RMS {np.sqrt(((out-img)**2).mean()):6.2f} DN, "
          f"최대 편차 {np.abs(out-img).max():6.1f} DN")
print()
print("=== 분리가능성 비용 ===")
for k in (5, 9, 15, 25):
    a, b = separable_cost(k)
    print(f"  {k}x{k}: 통째 {a:4d} 곱셈/화소, 분리 {b:3d} -> {a/b:.1f}배")

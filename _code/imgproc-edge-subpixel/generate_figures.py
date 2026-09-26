"""이미지처리 4편 그림 6개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/imgproc-edge-subpixel/      (한국어)
      ../../assets/img/posts/imgproc-edge-subpixel/en/   (영문)
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from edge_model import blurred_step, edge_width_10_90, sample_edge
from subpixel import ESTIMATORS, _gradient, partial_area

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                        "assets", "img", "posts", "imgproc-edge-subpixel")

NAMES_KO = {"parabola": "포물선 보간", "moment": "모멘트", "area": "면적"}
NAMES_EN = {"parabola": "Parabola", "moment": "Moment", "area": "Partial area"}
COLORS = {"parabola": "#c0392b", "moment": "#1e8449", "area": "#2b5c9b"}
MARKS = {"parabola": "o", "moment": "s", "area": "^"}

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT, "names": NAMES_KO,
        "pos": "위치 (화소)", "dn": "밝기 (DN)",
        "ideal": "이상적 계단", "blurred": "PSF 통과 후", "sampled": "화소가 담아낸 값",
        "fig1": "그림1. 엣지는 광학계에서 한 번, 화소 적분에서 또 한 번 뭉개진다",
        "grad": "기울기 크기", "para_fit": "포물선 적합", "true": "참 위치",
        "profile": "단면", "recovered": "모멘트가 되찾은 계단",
        "norm_area": r"정규화 단면 $(A-v)/(A-B)$", "cum": "누적합 = 엣지까지 거리",
        "fig2": "그림2. 세 계열은 같은 단면에서 서로 다른 것을 읽는다",
        "phase": "참 위치의 소수부 (화소)", "bias": "편향 (화소)",
        "fig3": "그림3. 포물선 보간의 편향은 화소 격자를 따라 되풀이된다",
        "blur_ax": r"PSF 흐림 폭 $w$ (화소)", "maxbias": "최대 절대 편향 (화소)",
        "w10": "10-90% 엣지 폭 (화소)",
        "fig4": "그림4. 포물선은 봉우리, 모멘트는 U자, 면적은 단조 증가 — 최적 흐림이 다 다르다",
        "noise_ax": "노이즈 (DN)", "scatter": "추정 위치 표준편차 (화소)",
        "winlen": "면적법의 창 길이 (화소)", "wb": "편향 (창을 자르지 않았을 때)",
        "ws": r"산포 (노이즈 5 DN)", "sqrtn": r"$\sqrt{n}$ 기준선",
        "fig5": "그림5. 편향이 작은 방법이 산포는 크고, 면적법의 창 길이는 그 안의 또 다른 교환이다",
        "rmse": "총 오차 RMSE (화소)", "winner": "가장 작은 총 오차를 내는 방법",
        "fig6": "그림6. 흐림 폭과 노이즈가 승자를 바꾼다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": {}, "legend": {}, "names": NAMES_EN,
        "pos": "Position (pixel)", "dn": "Intensity (DN)",
        "ideal": "Ideal step", "blurred": "After the PSF", "sampled": "Sampled by pixels",
        "fig1": "Fig 1. An edge is blurred once by the optics and again by pixel integration",
        "grad": "Gradient magnitude", "para_fit": "Parabola fit", "true": "True position",
        "profile": "Profile", "recovered": "Step recovered by moments",
        "norm_area": r"Normalized profile $(A-v)/(A-B)$", "cum": "Cumulative sum = distance to edge",
        "fig2": "Fig 2. The three families read different things from the same profile",
        "phase": "Fractional part of the true position (pixel)", "bias": "Bias (pixel)",
        "fig3": "Fig 3. The bias of parabolic interpolation repeats with the pixel grid",
        "blur_ax": r"PSF blur width $w$ (pixel)", "maxbias": "Max absolute bias (pixel)",
        "w10": "10-90% edge width (pixel)",
        "fig4": "Fig 4. A peak, a U and a rise — each method has its own optimal blur",
        "noise_ax": "Noise (DN)", "scatter": "Std of estimated position (pixel)",
        "winlen": "Window length of the area method (pixel)", "wb": "Bias (uncropped window)",
        "ws": r"Scatter (noise 5 DN)", "sqrtn": r"$\sqrt{n}$ reference",
        "fig5": "Fig 5. Smaller bias comes with larger scatter, and the area window is a trade of its own",
        "rmse": "Total error RMSE (pixel)", "winner": "Method with the smallest total error",
        "fig6": "Fig 6. Blur width and noise change the winner",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
N = 24
X0_REF = 11.37
PHASES = np.linspace(0.0, 1.0, 121)[:-1]

# 그림1: 두 번의 뭉갬
x_fine = np.linspace(6.0, 17.0, 1200)
step_ideal = blurred_step(x_fine, X0_REF, 0.0)
step_blur = blurred_step(x_fine, X0_REF, 1.0)
sampled = sample_edge(N, X0_REF, 1.0)

# 그림2: 세 계열의 원리
prof_demo = sample_edge(N, X0_REF, 1.0)
g_demo = np.abs(_gradient(prof_demo))
gi = int(np.argmax(g_demo[1:-1])) + 1
ga, gb, gc = g_demo[gi - 1], g_demo[gi], g_demo[gi + 1]
gpoly = np.polyfit([gi - 1, gi, gi + 1], [ga, gb, gc], 2)
gfit_x = np.linspace(gi - 2.2, gi + 2.2, 80)
lo_d, hi_d = prof_demo[:3].mean(), prof_demo[-3:].mean()
norm_demo = (hi_d - prof_demo) / (hi_d - lo_d)
cum_demo = np.cumsum(norm_demo)
mom_p = ESTIMATORS["moment"](prof_demo)

# 그림3: 픽셀 로킹
BLURS_LOCK = (0.4, 0.6, 1.0, 2.0)
lock = {b: {k: np.array([ESTIMATORS[k](sample_edge(N, 11.0 + ph, b)) - (11.0 + ph)
                         for ph in PHASES]) for k in ESTIMATORS} for b in BLURS_LOCK}

# 그림4: 편향 대 흐림 폭
BLUR_GRID = np.array([0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0])
maxbias = {k: [] for k in ESTIMATORS}
for b in BLUR_GRID:
    for k in ESTIMATORS:
        maxbias[k].append(max(abs(ESTIMATORS[k](sample_edge(N, 11.0 + ph, b)) - (11.0 + ph))
                              for ph in PHASES))
maxbias = {k: np.array(v) for k, v in maxbias.items()}

# 그림5: 노이즈 대 산포
NOISE_GRID = np.array([1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
N_TRIAL = 600
rng = np.random.default_rng(0)
scatter = {k: [] for k in ESTIMATORS}
for nz in NOISE_GRID:
    est = {k: np.empty(N_TRIAL) for k in ESTIMATORS}
    for t in range(N_TRIAL):
        p = sample_edge(N, X0_REF, 1.0) + nz * rng.standard_normal(N)
        for k in ESTIMATORS:
            est[k][t] = ESTIMATORS[k](p)
    for k in ESTIMATORS:
        scatter[k].append(est[k].std())
scatter = {k: np.array(v) for k, v in scatter.items()}

# 면적법의 창 길이 의존성. 산포는 창 길이의 제곱근을 따라 커지고, 창이 너무 짧으면
# 평탄부를 잘못 잡아 편향이 폭발한다. 두 곡선이 만나는 곳이 실용적인 창 길이다.
WIN_LENS = np.array([6, 8, 10, 12, 16, 20, 24, 32, 40, 52, 64])
rngw = np.random.default_rng(7)
win_bias, win_std = [], []
for n in WIN_LENS:
    true = n / 2.0 + 0.37
    clean = sample_edge(int(n), true, 1.0)
    win_bias.append(abs(partial_area(clean, window=None) - true))
    est = np.array([partial_area(sample_edge(int(n), true, 1.0)
                                 + 5.0 * rngw.standard_normal(int(n)), window=None)
                    for _ in range(400)])
    win_std.append(est.std())
win_bias, win_std = np.array(win_bias), np.array(win_std)

# 그림6: 총 오차. 편향은 무잡음 스윕의 제곱평균, 산포는 시뮬레이션으로 따로 구해 합친다.
BLUR6 = np.array([0.3, 0.5, 0.8, 1.2, 2.0])
NOISE6 = np.array([1.0, 3.0, 10.0, 30.0])
rng6 = np.random.default_rng(1)
rmse = {k: np.zeros((len(BLUR6), len(NOISE6))) for k in ESTIMATORS}
for i, b in enumerate(BLUR6):
    bias_rms = {k: np.sqrt(np.mean([(ESTIMATORS[k](sample_edge(N, 11.0 + ph, b)) - (11.0 + ph)) ** 2
                                    for ph in PHASES])) for k in ESTIMATORS}
    for j, nz in enumerate(NOISE6):
        est = {k: np.empty(300) for k in ESTIMATORS}
        for t in range(300):
            p = sample_edge(N, X0_REF, b) + nz * rng6.standard_normal(N)
            for k in ESTIMATORS:
                est[k][t] = ESTIMATORS[k](p)
        for k in ESTIMATORS:
            rmse[k][i, j] = np.hypot(bias_rms[k], est[k].std())
winner = np.empty((len(BLUR6), len(NOISE6)), dtype=object)
for i in range(len(BLUR6)):
    for j in range(len(NOISE6)):
        winner[i, j] = min(ESTIMATORS, key=lambda k: rmse[k][i, j])


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    F, LG, NM = L["font"], L["legend"], L["names"]

    def cap(fig, key, rect=(0, 0.06, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **F)
        fig.tight_layout(rect=rect)

    # 그림1
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(x_fine, step_ideal, ":", lw=1.3, color="#888", label=L["ideal"])
    ax.plot(x_fine, step_blur, "-", lw=1.6, color="#2b5c9b", label=L["blurred"])
    ax.plot(np.arange(N), sampled, "o", ms=5, color="#c0392b", label=L["sampled"])
    ax.axvline(X0_REF, ls="--", lw=1.0, color="k")
    ax.text(X0_REF + 0.15, 250, L["true"], fontsize=9, **F)
    ax.set_xlim(6, 17); ax.set_xlabel(L["pos"], **F); ax.set_ylabel(L["dn"], **F)
    ax.legend(prop=LG, fontsize=9); ax.grid(alpha=0.3)
    cap(fig, "fig1")
    fig.savefig(os.path.join(out, "fig1-edge-formation.png"), dpi=150)
    plt.close(fig)

    # 그림2
    fig, ax = plt.subplots(1, 3, figsize=(11.0, 3.6))
    ax[0].plot(np.arange(N), g_demo, "o-", ms=4, lw=1.0, color="#888", label=L["grad"])
    ax[0].plot(gfit_x, np.polyval(gpoly, gfit_x), "--", lw=1.8,
               color=COLORS["parabola"], label=L["para_fit"])
    ax[0].axvline(X0_REF, ls=":", lw=1.1, color="k")
    ax[0].set_xlim(gi - 4, gi + 4)
    ax[0].set_ylim(-0.05 * g_demo.max(), 1.15 * g_demo.max())
    ax[0].set_title(NM["parabola"], fontsize=10, **F)
    ax[0].legend(prop=LG, fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(np.arange(N), prof_demo, "o-", ms=4, lw=1.0, color="#888", label=L["profile"])
    ax[1].step(np.arange(N), np.where(np.arange(N) < mom_p, lo_d, hi_d), where="mid",
               lw=1.8, color=COLORS["moment"], label=L["recovered"])
    ax[1].axvline(X0_REF, ls=":", lw=1.1, color="k")
    ax[1].set_title(NM["moment"], fontsize=10, **F)
    ax[1].legend(prop=LG, fontsize=8); ax[1].grid(alpha=0.3)
    ax[2].bar(np.arange(N), norm_demo, color="#bcd", label=L["norm_area"])
    ax[2].plot(np.arange(N), cum_demo, "-", lw=1.8, color=COLORS["area"], label=L["cum"])
    ax[2].axvline(X0_REF, ls=":", lw=1.1, color="k")
    ax[2].set_title(NM["area"], fontsize=10, **F)
    ax[2].legend(prop=LG, fontsize=8, loc="center right"); ax[2].grid(alpha=0.3)
    for a in ax:
        a.set_xlabel(L["pos"], **F)
    cap(fig, "fig2", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig2-three-families.png"), dpi=150)
    plt.close(fig)

    # 그림3
    fig, ax = plt.subplots(2, 2, figsize=(9.2, 6.0), sharex=True, sharey=True)
    for a, b in zip(ax.ravel(), BLURS_LOCK):
        for k in ESTIMATORS:
            a.plot(PHASES, lock[b][k], lw=1.5, color=COLORS[k], label=NM[k])
        a.axhline(0, color="k", lw=0.6)
        a.set_title(f"$w$ = {b}", fontsize=10)
        a.grid(alpha=0.3)
    ax[0, 0].legend(prop=LG, fontsize=8)
    for a in ax[1]:
        a.set_xlabel(L["phase"], **F)
    for a in ax[:, 0]:
        a.set_ylabel(L["bias"], **F)
    cap(fig, "fig3", rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(out, "fig3-pixel-locking.png"), dpi=150)
    plt.close(fig)

    # 그림4
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for k in ESTIMATORS:
        ax.semilogy(BLUR_GRID, np.maximum(maxbias[k], 1e-6), MARKS[k] + "-", ms=4,
                    color=COLORS[k], label=NM[k])
    ax.set_xlabel(L["blur_ax"], **F); ax.set_ylabel(L["maxbias"], **F)
    ax.legend(prop=LG, fontsize=9); ax.grid(alpha=0.3, which="both")
    sec = ax.secondary_xaxis("top", functions=(edge_width_10_90, lambda v: v / 2.5631))
    sec.set_xlabel(L["w10"], **F)
    cap(fig, "fig4", rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(out, "fig4-bias-vs-blur.png"), dpi=150)
    plt.close(fig)

    # 그림5
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.0))
    for k in ESTIMATORS:
        ax[0].loglog(NOISE_GRID, scatter[k], MARKS[k] + "-", ms=4,
                     color=COLORS[k], label=NM[k])
    ax[0].set_xlabel(L["noise_ax"], **F); ax[0].set_ylabel(L["scatter"], **F)
    ax[0].legend(prop=LG, fontsize=9); ax[0].grid(alpha=0.3, which="both")
    ax[1].loglog(WIN_LENS, np.maximum(win_bias, 5e-5), "o-", ms=4,
                 color="#c0392b", label=L["wb"])
    ax[1].loglog(WIN_LENS, np.maximum(win_std, 5e-5), "^-", ms=4,
                 color=COLORS["area"], label=L["ws"])
    ref = win_std[4] * np.sqrt(WIN_LENS / WIN_LENS[4])
    ax[1].loglog(WIN_LENS, ref, ":", lw=1.2, color="#555", label=L["sqrtn"])
    ax[1].set_ylim(3e-5, 3.0)
    ax[1].set_xlabel(L["winlen"], **F)
    ax[1].legend(prop=LG, fontsize=8.5); ax[1].grid(alpha=0.3, which="both")
    cap(fig, "fig5", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig5-noise-scatter.png"), dpi=150)
    plt.close(fig)

    # 그림6
    fig, ax = plt.subplots(1, len(NOISE6), figsize=(11.6, 3.4), sharey=True)
    for j, (a, nz) in enumerate(zip(ax, NOISE6)):
        for k in ESTIMATORS:
            a.semilogy(BLUR6, rmse[k][:, j], MARKS[k] + "-", ms=4,
                       color=COLORS[k], label=NM[k])
        a.set_title(f"{L['noise_ax']} = {nz:.0f}", fontsize=9.5, **F)
        a.set_xlabel(L["blur_ax"], **F); a.grid(alpha=0.3, which="both")
    ax[0].set_ylabel(L["rmse"], **F); ax[0].legend(prop=LG, fontsize=8)
    cap(fig, "fig6", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig6-total-error.png"), dpi=150)
    plt.close(fig)


for lang, L in LABELS.items():
    draw(L)
    print(f"[{lang}] 그림 6개 저장: {L['dir']}")

print()
print("=== 그림3~4: 무잡음 편향 ===")
print(f"{'blur':>6s} {'10-90폭':>8s} " + " ".join(f"{k:>10s}" for k in ESTIMATORS))
for i, b in enumerate(BLUR_GRID):
    print(f"{b:6.2f} {edge_width_10_90(b):8.2f} "
          + " ".join(f"{maxbias[k][i]:10.4f}" for k in ESTIMATORS))
for k in ESTIMATORS:
    i = int(np.argmin(maxbias[k]))
    j = int(np.argmax(maxbias[k]))
    print(f"  {k:9s} 최소 편향 {maxbias[k][i]:.4f} @ w={BLUR_GRID[i]}, "
          f"최대 {maxbias[k][j]:.4f} @ w={BLUR_GRID[j]}")
print()
print("=== 그림5: 노이즈 대 산포 (w=1.0) ===")
for i, nz in enumerate(NOISE_GRID):
    print(f"  노이즈 {nz:4.0f} DN: "
          + "  ".join(f"{k} {scatter[k][i]:.4f}" for k in ESTIMATORS))
print(f"  면적/모멘트 산포비 = {np.mean(scatter['area']/scatter['moment']):.1f}배")
print()
print("=== 그림5(b): 면적법 창 길이 (w=1.0, 노이즈 5 DN) ===")
for n, b, sd in zip(WIN_LENS, win_bias, win_std):
    print(f"  창 {n:3d}: 편향 {b:8.4f}  산포 {sd:.4f}")
print(f"  자동 창을 썼을 때 (w=1.0): 창 길이 12, 산포 {scatter['area'][2]:.4f}")
print()
print("=== 그림6: 총 오차 승자 ===")
print(f"{'blur':>6s} " + " ".join(f"{'노이즈 '+str(int(n)):>14s}" for n in NOISE6))
for i, b in enumerate(BLUR6):
    print(f"{b:6.2f} " + " ".join(f"{winner[i,j]:>14s}" for j in range(len(NOISE6))))

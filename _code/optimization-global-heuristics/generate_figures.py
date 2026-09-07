"""최적화 4편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-global-heuristics/      (한국어)
      ../../assets/img/posts/optimization-global-heuristics/en/   (영문)

측정 데이터는 1~3편과 동일(SiO2/Si 단층, 두께 1490nm, seed=0 노이즈).
"어려운" 초기값 d0=1300은 3편에서 LM이 1293nm 국소최솟값에 갇히던 지점을 그대로 쓴다.

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은 데이터가
완전히 동일하고 표기만 다르다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from reflectance_model import reflectance
from levenberg_marquardt import objective, levenberg_marquardt
from simulated_annealing import simulated_annealing
from basin_hopping import basin_hopping

# 한글 텍스트에만 한글 폰트를 지정한다. 시리즈의 다른 편과 같은 AppleGothic 을 쓴다.
# macOS 전용 설정이다. 다른 OS는 나눔고딕 등 설치된 한글 폰트 이름으로 교체할 것.
# 그리스 문자(T_0 의 아래첨자 등)는 폰트에 글리프가 없을 수 있으므로 mathtext 로 쓴다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-global-heuristics",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT, "tickfont": "AppleGothic",
        "thickness": "두께 d (nm)", "J_log": "J(d)  (log scale)",
        "fig1": "그림1. 목적함수 지형 — 세 개의 basin (점선: basin 경계)",
        "sa_low": r"$T_0$ 매우 낮음 (0.0005)", "sa_mid": r"$T_0$ 적당 (0.02)",
        "sa_high": r"$T_0$ 매우 높음 (5.0)",
        "thickness_est": "두께 추정값 (nm)",
        "fig2": "그림2. SA 온도 스케줄에 따른 탐색 궤적 (d0=1300nm)",
        "bh_low": r"$T_0$ 매우 낮음 (0.0005)", "bh_mid": r"$T_0$ 적당 (0.05)",
        "bh_high": r"$T_0$ 매우 높음 (5.0)",
        "fig3": "그림3. Basin-hopping 온도 스케줄에 따른 궤적 (d0=1300nm)",
        "bars": ["단일 LM\n(1회 시도)", "Simulated\nAnnealing", "Basin-\nHopping"],
        "success_rate": "전역 최솟값 도달 성공률 (%)",
        "fig4": "그림4. 무작위 시작 {}회 성공률 비교",
        "start": "시작 (d0=1300 → LM 수렴)", "last": "마지막 위치",
        "fig5": "그림5. 지형 위에 겹친 basin-hopping 궤적",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "tickfont": None,
        "thickness": "Thickness d (nm)", "J_log": "J(d)  (log scale)",
        "fig1": "Fig 1. The objective landscape — three basins (dashed: basin boundaries)",
        "sa_low": r"$T_0$ very low (0.0005)", "sa_mid": r"$T_0$ moderate (0.02)",
        "sa_high": r"$T_0$ very high (5.0)",
        "thickness_est": "Estimated thickness (nm)",
        "fig2": "Fig 2. SA search trajectories by temperature schedule (d0=1300nm)",
        "bh_low": r"$T_0$ very low (0.0005)", "bh_mid": r"$T_0$ moderate (0.05)",
        "bh_high": r"$T_0$ very high (5.0)",
        "fig3": "Fig 3. Basin-hopping trajectories by temperature schedule (d0=1300nm)",
        "bars": ["single LM\n(one attempt)", "Simulated\nAnnealing", "Basin-\nHopping"],
        "success_rate": "Rate of reaching the global minimum (%)",
        "fig4": "Fig 4. Success rate over {} random starts",
        "start": "start (d0=1300, after LM)", "last": "final position",
        "fig5": "Fig 5. A basin-hopping trajectory over the landscape",
    },
}

# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0
noise_std = 0.004
initial_guess_hard = 1300.0  # 3편에서 LM이 1293nm에 갇히던 시작점

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# 그림1: 목적함수 지형 J(d), 세 basin과 그 경계.
ds = np.linspace(1200, 1800, 6001)
Js = np.array([objective(d, wavelength_nm, R_measured) for d in ds])
minima = [(ds[i], Js[i]) for i in range(1, len(ds) - 1) if Js[i] < Js[i-1] and Js[i] < Js[i+1]]
basin_left, basin_right = 1391.3, 1589.2

# 그림2: SA 온도 스케줄 비교 (d0=1300, n_iter=300, step_sigma=80).
SA_ITER = 300
sa_runs = []
for T0, key, color in [(0.0005, "sa_low", "C0"), (0.02, "sa_mid", "C2"), (5.0, "sa_high", "C3")]:
    d_hist, J_hist, T_hist, best_d, best_J = simulated_annealing(
        initial_guess_hard, wavelength_nm, R_measured,
        n_iter=SA_ITER, T0=T0, cooling=0.97, step_sigma=80.0,
        rng=np.random.default_rng(1))
    # basin 경계를 처음 넘은 스텝과, 전역 최솟값 5nm 이내에 처음 든 스텝
    crossed = next((i for i, d in enumerate(d_hist) if d > basin_left), None)
    hit = next((i for i, d in enumerate(d_hist) if abs(d - true_thickness) < 5), None)
    sa_runs.append((T0, key, color, d_hist, best_d, best_J, crossed, hit))

# 그림3: Basin-hopping 온도 스케줄 비교 (d0=1300, n_hops=25, perturb_sigma=150).
BH_HOPS = 25
bh_runs = []
for T0, key, color in [(0.0005, "bh_low", "C0"), (0.05, "bh_mid", "C2"), (5.0, "bh_high", "C3")]:
    trace_d, trace_J, best_d, best_J, n_accept = basin_hopping(
        initial_guess_hard, wavelength_nm, R_measured, n_hops=BH_HOPS, T0=T0,
        perturb_sigma=150.0, lm_kwargs=dict(tau=1e-3, n_iter=30),
        rng=np.random.default_rng(1))
    hit = next((i for i, d in enumerate(trace_d) if abs(d - true_thickness) < 5), None)
    bh_runs.append((T0, key, color, trace_d, best_d, best_J, n_accept, hit))

# 그림4: 성공률 비교 (무작위 시작 50회, 단일 LM vs SA vs Basin-hopping).
N_SEEDS = 50
success_ms = success_bh = success_sa = 0
hops_bh, iters_sa = [], []
for seed in range(N_SEEDS):
    d0_random = np.random.default_rng(seed).uniform(1200, 1800)

    d_hist_ms, J_hist_ms, _ = levenberg_marquardt(
        d0_random, wavelength_nm, R_measured, tau=1e-3, n_iter=30)
    success_ms += abs(d_hist_ms[-1] - true_thickness) < 5

    tr_d, tr_J, bd, bJ, n_acc = basin_hopping(
        d0_random, wavelength_nm, R_measured, n_hops=15, T0=0.05, perturb_sigma=150.0,
        lm_kwargs=dict(tau=1e-3, n_iter=30), rng=np.random.default_rng(seed + 1000))
    ok = abs(bd - true_thickness) < 5
    success_bh += ok
    if ok:
        hops_bh.append(next(i for i, d in enumerate(tr_d) if abs(d - true_thickness) < 5))

    dh_sa, Jh_sa, Th_sa, bd2, bJ2 = simulated_annealing(
        d0_random, wavelength_nm, R_measured, n_iter=300, T0=0.02, cooling=0.97,
        step_sigma=80.0, rng=np.random.default_rng(seed + 2000))
    ok2 = abs(bd2 - true_thickness) < 5
    success_sa += ok2
    if ok2:
        iters_sa.append(next(i for i, d in enumerate(dh_sa) if abs(d - true_thickness) < 5))

rates = [success_ms / N_SEEDS * 100, success_sa / N_SEEDS * 100, success_bh / N_SEEDS * 100]

# 그림5: 지형 위에 basin-hopping 실제 궤적 오버레이 (d0=1300, T0=0.05).
# 그림3의 T0=0.05 실행과 같은 seed 를 쓰므로 궤적이 같다.
overlay_d, overlay_J, overlay_best_d, overlay_best_J, overlay_accept = basin_hopping(
    initial_guess_hard, wavelength_nm, R_measured, n_hops=BH_HOPS, T0=0.05,
    perturb_sigma=150.0, lm_kwargs=dict(tau=1e-3, n_iter=30),
    rng=np.random.default_rng(1))


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L):
    """주어진 라벨 묶음으로 그림 5개를 그린다. 계산은 위에서 한 번만 했다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림1: 목적함수 지형
    plt.figure(figsize=(6, 4))
    plt.plot(ds, Js, color="C0", linewidth=1.2)
    for d, J in minima:
        plt.plot(d, J, "o", color="C3", markersize=5)
    plt.axvline(basin_left, color="gray", linestyle="--", linewidth=1)
    plt.axvline(basin_right, color="gray", linestyle="--", linewidth=1)
    plt.yscale("log")
    plt.xlabel(L["thickness"], **F)
    plt.ylabel(L["J_log"], **F)
    plt.title(L["fig1"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig1-objective-landscape.png"), dpi=150)
    plt.close()

    # 그림2: SA 온도 스케줄 비교
    plt.figure(figsize=(6, 4))
    for _, key, color, d_hist, *_ in sa_runs:
        plt.plot(d_hist, "-", color=color, linewidth=1, label=L[key])
    plt.axhline(true_thickness, color="gray", linestyle="--", linewidth=1)
    plt.xlabel("iteration")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig2"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig2-sa-temperature-schedule.png"), dpi=150)
    plt.close()

    # 그림3: basin-hopping 온도 스케줄 비교
    plt.figure(figsize=(6, 4))
    for _, key, color, trace_d, *_ in bh_runs:
        plt.plot(trace_d, "o-", color=color, markersize=4, linewidth=1, label=L[key])
    plt.axhline(true_thickness, color="gray", linestyle="--", linewidth=1)
    plt.xlabel("hop")
    plt.ylabel(L["thickness_est"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig3"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig3-basinhopping-temperature-schedule.png"), dpi=150)
    plt.close()

    # 그림4: 성공률 비교
    plt.figure(figsize=(6, 4))
    bars = plt.bar(L["bars"], rates, color=["C3", "C2", "C0"])
    if L["tickfont"]:
        plt.xticks(fontfamily=L["tickfont"])
    for bar, rate in zip(bars, rates):
        plt.text(bar.get_x() + bar.get_width() / 2, rate + 2, f"{rate:.0f}%",
                 ha="center", **F)
    plt.ylim(0, 110)
    plt.ylabel(L["success_rate"], **F)
    plt.title(L["fig4"].format(N_SEEDS), **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig4-success-rate-comparison.png"), dpi=150)
    plt.close()

    # 그림5: 지형 위에 겹친 basin-hopping 궤적
    plt.figure(figsize=(6, 4))
    plt.plot(ds, Js, color="lightgray", linewidth=1.2, zorder=1)
    plt.plot(overlay_d, overlay_J, "o-", color="C3", markersize=5, linewidth=1, zorder=2)
    plt.plot(overlay_d[0], overlay_J[0], "s", color="C0", markersize=9, zorder=3, label=L["start"])
    plt.plot(overlay_d[-1], overlay_J[-1], "*", color="C2", markersize=14, zorder=3, label=L["last"])
    plt.yscale("log")
    plt.xlabel(L["thickness"], **F)
    plt.ylabel(L["J_log"], **F)
    plt.legend(prop=LF, fontsize=8)
    plt.title(L["fig5"], **F)
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig5-trajectory-on-landscape.png"), dpi=150)
    plt.close()


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

# ---------------------------------------------------------------------------
# 본문에 인용할 수치 출력 (언어와 무관하므로 한 번만 찍는다)
# ---------------------------------------------------------------------------
print()
print("=== 그림1: local minima ===")
for d, J in minima:
    print(f"  d={d:.2f} J={J:.5f}")
print(f"basin 경계: {basin_left}, {basin_right} (폭 {basin_right - basin_left:.1f}nm)")
print()
print("=== 그림2: SA temperature sweep (d0=1300) ===")
for T0, key, _, d_hist, best_d, best_J, crossed, hit in sa_runs:
    print(f"  T0={T0:<7g} d_range=[{d_hist.min():.1f},{d_hist.max():.1f}] "
          f"basin 경계 통과 step={crossed} 전역최솟값 5nm 이내 첫 도달 step={hit} "
          f"best_J={best_J:.5f}")
print()
print("=== 그림3: basin-hopping temperature sweep (d0=1300) ===")
for T0, key, _, trace_d, best_d, best_J, n_accept, hit in bh_runs:
    print(f"  T0={T0:<7g} accept={n_accept}/{BH_HOPS} best_d={best_d:.2f} "
          f"best_J={best_J:.5f} 전역최솟값 도달 hop={hit}")
print(f"  T0=5.0 궤적: {np.round(np.unique(bh_runs[2][3]), 1)}")
print(f"  T0=0.0005 와 T0=0.05 궤적이 같은가: "
      f"{np.allclose(bh_runs[0][3], bh_runs[1][3])}")
print()
print("=== 그림4: success rate (N=50 random start) ===")
print(f"  single LM: {success_ms}/{N_SEEDS} = {rates[0]:.0f}%")
print(f"  SA: {success_sa}/{N_SEEDS} = {rates[1]:.0f}%, avg iters-to-hit={np.mean(iters_sa):.1f}")
print(f"  basin-hopping: {success_bh}/{N_SEEDS} = {rates[2]:.0f}%, "
      f"avg hops-to-hit={np.mean(hops_bh):.1f}")
print()
print("=== 그림5: 오버레이 궤적 (T0=0.05) ===")
print(f"  첫 바닥={overlay_d[0]:.2f}, 마지막={overlay_d[-1]:.2f}, "
      f"전역최솟값 도달 hop={next(i for i, d in enumerate(overlay_d) if abs(d - true_thickness) < 5)}")

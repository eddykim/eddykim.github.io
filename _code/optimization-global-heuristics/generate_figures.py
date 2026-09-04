"""4편 그림 5개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-global-heuristics/ 에 fig1~5 저장

측정 데이터는 1~3편과 동일(SiO2/Si 단층, 두께 1490nm, seed=0 노이즈).
"어려운" 초기값 d0=1300은 3편에서 LM이 1293nm 국소최솟값에 갇히던 지점을 그대로 쓴다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from reflectance_model import reflectance
from levenberg_marquardt import objective, levenberg_marquardt
from simulated_annealing import simulated_annealing
from basin_hopping import basin_hopping

KFONT = {"fontfamily": "NanumGothic"}
LFONT = {"family": "NanumGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "optimization-global-heuristics",
)
os.makedirs(OUT_DIR, exist_ok=True)

np.random.seed(0)

wavelength_nm = np.linspace(450, 750, 300)
true_thickness = 1490.0
noise_std = 0.004
initial_guess_hard = 1300.0  # 3편에서 LM이 1293nm에 갇히던 시작점

R_true = reflectance(true_thickness, wavelength_nm)
R_measured = R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

# ---------------------------------------------------------------------------
# 그림1: 목적함수 지형 J(d), 세 basin과 그 경계.
# ---------------------------------------------------------------------------
ds = np.linspace(1200, 1800, 6001)
Js = np.array([objective(d, wavelength_nm, R_measured) for d in ds])
minima = [(ds[i], Js[i]) for i in range(1, len(ds) - 1) if Js[i] < Js[i-1] and Js[i] < Js[i+1]]
basin_left, basin_right = 1391.3, 1589.2

plt.figure(figsize=(6, 4))
plt.plot(ds, Js, color="C0", linewidth=1.2)
for d, J in minima:
    plt.plot(d, J, "o", color="C3", markersize=5)
plt.axvline(basin_left, color="gray", linestyle="--", linewidth=1)
plt.axvline(basin_right, color="gray", linestyle="--", linewidth=1)
plt.yscale("log")
plt.xlabel("두께 d (nm)", **KFONT)
plt.ylabel("J(d)  (log scale)", **KFONT)
plt.title("그림1. 목적함수 지형 — 세 개의 basin (점선: basin 경계)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig1-objective-landscape.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림2: SA 온도 스케줄 비교 (d0=1300, n_iter=300, step_sigma=80).
# ---------------------------------------------------------------------------
SA_ITER = 300
sa_runs = {}
for T0, label, color in [
    (0.0005, r"$T_0$ 매우 낮음 (0.0005)", "C0"),
    (0.02, r"$T_0$ 적당 (0.02)", "C2"),
    (5.0, r"$T_0$ 매우 높음 (5.0)", "C3"),
]:
    rng = np.random.default_rng(1)
    d_hist, J_hist, T_hist, best_d, best_J = simulated_annealing(
        initial_guess_hard, wavelength_nm, R_measured,
        n_iter=SA_ITER, T0=T0, cooling=0.97, step_sigma=80.0, rng=rng)
    sa_runs[label] = (d_hist, best_d, best_J, color)

plt.figure(figsize=(6, 4))
for label, (d_hist, best_d, best_J, color) in sa_runs.items():
    plt.plot(d_hist, color=color, linewidth=1, alpha=0.8, label=label)
plt.axhline(true_thickness, color="gray", linestyle="--", label="실제 두께")
plt.xlabel("iteration")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림2. SA 온도 스케줄에 따른 탐색 궤적 (d0=1300nm)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig2-sa-temperature-schedule.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림3: Basin-hopping 온도 스케줄 비교 (d0=1300, n_hops=25, perturb_sigma=150).
# ---------------------------------------------------------------------------
BH_HOPS = 25
bh_runs = {}
for T0, label, color in [
    (0.0005, r"$T_0$ 매우 낮음 (0.0005)", "C0"),
    (0.05, r"$T_0$ 적당 (0.05)", "C2"),
    (5.0, r"$T_0$ 매우 높음 (5.0)", "C3"),
]:
    rng = np.random.default_rng(1)
    trace_d, trace_J, best_d, best_J, n_accept = basin_hopping(
        initial_guess_hard, wavelength_nm, R_measured, n_hops=BH_HOPS, T0=T0,
        perturb_sigma=150.0, lm_kwargs=dict(tau=1e-3, n_iter=30), rng=rng)
    bh_runs[label] = (trace_d, best_d, best_J, n_accept, color)

plt.figure(figsize=(6, 4))
for label, (trace_d, best_d, best_J, n_accept, color) in bh_runs.items():
    plt.plot(trace_d, "o-", color=color, markersize=4, linewidth=1, label=label)
plt.axhline(true_thickness, color="gray", linestyle="--", label="실제 두께")
plt.xlabel("hop")
plt.ylabel("두께 추정값 (nm)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림3. Basin-hopping 온도 스케줄에 따른 궤적 (d0=1300nm)", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig3-basinhopping-temperature-schedule.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림4: 성공률 비교 (무작위 시작 50회, 단일 LM vs SA vs Basin-hopping).
# ---------------------------------------------------------------------------
N_SEEDS = 50
success_ms = success_bh = success_sa = 0
hops_bh, iters_sa = [], []
for seed in range(N_SEEDS):
    rng = np.random.default_rng(seed)
    d0_random = rng.uniform(1200, 1800)

    d_hist_ms, J_hist_ms, _ = levenberg_marquardt(d0_random, wavelength_nm, R_measured, tau=1e-3, n_iter=30)
    success_ms += abs(d_hist_ms[-1] - true_thickness) < 5

    rng_bh = np.random.default_rng(seed + 1000)
    trace_d, trace_J, best_d, best_J, n_accept = basin_hopping(
        d0_random, wavelength_nm, R_measured, n_hops=15, T0=0.05, perturb_sigma=150.0,
        lm_kwargs=dict(tau=1e-3, n_iter=30), rng=rng_bh)
    ok = abs(best_d - true_thickness) < 5
    success_bh += ok
    if ok:
        hops_bh.append(next(i for i, d in enumerate(trace_d) if abs(d - true_thickness) < 5))

    rng_sa = np.random.default_rng(seed + 2000)
    d_hist_sa, J_hist_sa, T_hist_sa, best_d2, best_J2 = simulated_annealing(
        d0_random, wavelength_nm, R_measured, n_iter=300, T0=0.02, cooling=0.97,
        step_sigma=80.0, rng=rng_sa)
    ok2 = abs(best_d2 - true_thickness) < 5
    success_sa += ok2
    if ok2:
        iters_sa.append(next(i for i, d in enumerate(d_hist_sa) if abs(d - true_thickness) < 5))

rates = [success_ms / N_SEEDS * 100, success_sa / N_SEEDS * 100, success_bh / N_SEEDS * 100]
plt.figure(figsize=(6, 4))
bars = plt.bar(["단일 LM\n(1회 시도)", "Simulated\nAnnealing", "Basin-\nHopping"], rates,
                color=["C3", "C2", "C0"])
plt.xticks(fontfamily="NanumGothic")
for bar, rate in zip(bars, rates):
    plt.text(bar.get_x() + bar.get_width() / 2, rate + 2, f"{rate:.0f}%",
              ha="center", **KFONT)
plt.ylim(0, 110)
plt.ylabel("전역 최솟값 도달 성공률 (%)", **KFONT)
plt.title(f"그림4. 무작위 시작 {N_SEEDS}회 성공률 비교", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig4-success-rate-comparison.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 그림5: 지형 위에 basin-hopping 실제 궤적 오버레이 (d0=1300, T0=0.05).
# ---------------------------------------------------------------------------
rng = np.random.default_rng(1)
trace_d, trace_J, best_d, best_J, n_accept = basin_hopping(
    initial_guess_hard, wavelength_nm, R_measured, n_hops=BH_HOPS, T0=0.05,
    perturb_sigma=150.0, lm_kwargs=dict(tau=1e-3, n_iter=30), rng=rng)

plt.figure(figsize=(6, 4))
plt.plot(ds, Js, color="lightgray", linewidth=1.2, zorder=1)
plt.plot(trace_d, trace_J, "o-", color="C3", markersize=5, linewidth=1, zorder=2)
plt.plot(trace_d[0], trace_J[0], "s", color="C0", markersize=9, zorder=3, label="시작 (d0=1300 → LM 수렴)")
plt.plot(trace_d[-1], trace_J[-1], "*", color="C2", markersize=14, zorder=3, label="마지막 위치")
plt.yscale("log")
plt.xlabel("두께 d (nm)", **KFONT)
plt.ylabel("J(d)  (log scale)", **KFONT)
plt.legend(prop=LFONT, fontsize=8)
plt.title("그림5. 지형 위에 겹친 basin-hopping 궤적", **KFONT)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "fig5-trajectory-on-landscape.png"), dpi=150)
plt.close()

# ---------------------------------------------------------------------------
# 본문에 인용할 수치 출력
# ---------------------------------------------------------------------------
print("=== 그림1: local minima ===")
for d, J in minima:
    print(f"  d={d:.2f} J={J:.5f}")
print(f"basin 경계: {basin_left}, {basin_right}")
print()
print("=== 그림2: SA temperature sweep (d0=1300) ===")
for label, (d_hist, best_d, best_J, color) in sa_runs.items():
    print(f"{label}: final_d={d_hist[-1]:.2f} best_d={best_d:.2f} best_J={best_J:.5f} "
          f"d_range=[{d_hist.min():.1f},{d_hist.max():.1f}]")
print()
print("=== 그림3: basin-hopping temperature sweep (d0=1300) ===")
for label, (trace_d, best_d, best_J, n_accept, color) in bh_runs.items():
    print(f"{label}: accept={n_accept}/{BH_HOPS} final_d={trace_d[-1]:.2f} best_d={best_d:.2f} best_J={best_J:.5f}")
print()
print("=== 그림4: success rate (N=50 random start) ===")
print(f"single LM: {success_ms}/{N_SEEDS} = {rates[0]:.0f}%")
print(f"SA: {success_sa}/{N_SEEDS} = {rates[1]:.0f}%, avg iters-to-hit={np.mean(iters_sa):.1f}")
print(f"basin-hopping: {success_bh}/{N_SEEDS} = {rates[2]:.0f}%, avg hops-to-hit={np.mean(hops_bh):.1f}")
print()
print("=== 그림5: trajectory ===")
print(f"trace_d: {np.round(trace_d, 2)}")

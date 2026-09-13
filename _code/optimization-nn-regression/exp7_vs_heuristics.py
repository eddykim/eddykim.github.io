"""실험 7 — 4편의 전역 탐색(SA / Basin-Hopping)과 같은 잣대로 비교.

4편은 무작위 시작 50회에서 단일 LM 36%, SA 100%, BH 100% 를 보고했다.
회귀망→LM 도 100% 이므로 성공률로는 구분되지 않는다. 차이는 **비용**에서 난다.

공정한 비용 지표로 '전방 모델(Airy 공식) 평가 횟수'를 센다. 벽시계 시간은 구현에
좌우되지만 전방 모델 호출 수는 알고리즘의 성질이다. SA·BH·LM 모두
levenberg_marquardt 모듈의 reflectance 를 거치므로 그 이름을 감싸 계수한다.
"""
import warnings
warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

import time

import numpy as np

import levenberg_marquardt as lmmod
import train_data as tdmod
from common import (WL, TRUE_D, NOISE_STD, D_RANGE, success, lm_fit,
                    Scaler, train_net, predict_d, rmse_nm)
from train_data import make_dataset, measured_reference
from simulated_annealing import simulated_annealing
from basin_hopping import basin_hopping

# --- 전방 모델 호출 계수기 ------------------------------------------------
# reflectance 를 이름으로 가져다 쓴 모듈이 둘 있으므로 양쪽 다 감싼다.
# (levenberg_marquardt: LM·SA·BH 경로 / train_data: 학습데이터 생성 경로)
_orig = lmmod.reflectance
CNT = {"n": 0}


def _counting(*a, **k):
    CNT["n"] += 1
    return _orig(*a, **k)


lmmod.reflectance = _counting
tdmod.reflectance = _counting


def measure(fn):
    """fn() 실행 중의 전방 모델 호출 수와 벽시계 시간을 함께 돌려준다."""
    CNT["n"] = 0
    t0 = time.perf_counter()
    out = fn()
    return out, CNT["n"], time.perf_counter() - t0


# --- 회귀망 준비 (일회성 비용) --------------------------------------------
print("회귀망 학습 (일회성 비용)")
rng = np.random.default_rng(0)
CNT["n"] = 0
t0 = time.perf_counter()
X, d = make_dataset(8000, D_RANGE, WL, rng)
gen_calls, gen_s = CNT["n"], time.perf_counter() - t0
Xv, dv = make_dataset(2000, D_RANGE, WL, np.random.default_rng(99))
sc = Scaler(*D_RANGE).fit_inputs(X)
net, hist, train_s = train_net(X, d, sc, [300, 256, 1], np.random.default_rng(1),
                               epochs=600, lr=1e-3, X_val=Xv, d_val=dv)
print(f"  학습데이터 8000개 생성: 전방모델 {gen_calls}회, {gen_s:.1f}s")
print(f"  Adam 600 epoch 학습   : 전방모델 0회,   {train_s:.0f}s")
print(f"  검증 RMSE {rmse_nm(net, sc, Xv, dv):.2f} nm")

# --- 시료당 비교 ----------------------------------------------------------
N = 100
methods = {}
start_rng = np.random.default_rng(2024)
starts = start_rng.uniform(*D_RANGE, size=N)

for name in ("단일 LM (무작위 시작)", "Simulated Annealing", "Basin-Hopping", "회귀망 → LM"):
    methods[name] = {"succ": 0, "calls": [], "time": []}

for i in range(N):
    R = measured_reference(WL, seed=i)
    d0 = starts[i]

    (res, c, t) = measure(lambda: lm_fit(d0, WL, R))
    methods["단일 LM (무작위 시작)"]["succ"] += success(res[0])
    methods["단일 LM (무작위 시작)"]["calls"].append(c)
    methods["단일 LM (무작위 시작)"]["time"].append(t)

    (res, c, t) = measure(lambda: simulated_annealing(
        d0, WL, R, n_iter=300, T0=0.02, cooling=0.97, step_sigma=80.0,
        rng=np.random.default_rng(2000 + i)))
    methods["Simulated Annealing"]["succ"] += success(res[3])
    methods["Simulated Annealing"]["calls"].append(c)
    methods["Simulated Annealing"]["time"].append(t)

    (res, c, t) = measure(lambda: basin_hopping(
        d0, WL, R, n_hops=15, T0=0.05, perturb_sigma=150.0,
        lm_kwargs=dict(tau=1e-3, n_iter=30), rng=np.random.default_rng(1000 + i)))
    methods["Basin-Hopping"]["succ"] += success(res[2])
    methods["Basin-Hopping"]["calls"].append(c)
    methods["Basin-Hopping"]["time"].append(t)

    def nn_then_lm():
        dn = float(predict_d(net, sc, R)[0])
        return lm_fit(dn, WL, R)
    (res, c, t) = measure(nn_then_lm)
    methods["회귀망 → LM"]["succ"] += success(res[0])
    methods["회귀망 → LM"]["calls"].append(c)
    methods["회귀망 → LM"]["time"].append(t)

print(f"\n[B] 시료마다 측정이 다른 경우 — 노이즈 seed 0~{N-1}, 무작위 시작 {N}회")
print(f"{'방법':26s}{'성공률(%)':>10s}{'전방모델 호출/시료':>18s}{'시간/시료(ms)':>15s}")
for k, v in methods.items():
    print(f"{k:26s}{v['succ']:10d}{np.mean(v['calls']):18.0f}{np.mean(v['time'])*1e3:15.2f}")

# --- [A] 4편 설정 재현: 측정 고정(seed 0), 초기값만 무작위 ------------------
R_fix = measured_reference(WL, seed=0)
mA = {k: {"succ": 0, "calls": []} for k in methods}
for i in range(N):
    d0 = starts[i]
    (res, c, _) = measure(lambda: lm_fit(d0, WL, R_fix))
    mA["단일 LM (무작위 시작)"]["succ"] += success(res[0]); mA["단일 LM (무작위 시작)"]["calls"].append(c)
    (res, c, _) = measure(lambda: simulated_annealing(
        d0, WL, R_fix, n_iter=300, T0=0.02, cooling=0.97, step_sigma=80.0,
        rng=np.random.default_rng(2000 + i)))
    mA["Simulated Annealing"]["succ"] += success(res[3]); mA["Simulated Annealing"]["calls"].append(c)
    (res, c, _) = measure(lambda: basin_hopping(
        d0, WL, R_fix, n_hops=15, T0=0.05, perturb_sigma=150.0,
        lm_kwargs=dict(tau=1e-3, n_iter=30), rng=np.random.default_rng(1000 + i)))
    mA["Basin-Hopping"]["succ"] += success(res[2]); mA["Basin-Hopping"]["calls"].append(c)
    (res, c, _) = measure(lambda: lm_fit(float(predict_d(net, sc, R_fix)[0]), WL, R_fix))
    mA["회귀망 → LM"]["succ"] += success(res[0]); mA["회귀망 → LM"]["calls"].append(c)

print(f"\n[A] 4편 설정 재현 — 측정 고정(seed 0), 초기값만 무작위 {N}회")
print(f"  (4편 보고치: 단일 LM 36%, SA 100%, BH 100% — 무작위 시작 50회)")
print(f"{'방법':26s}{'성공률(%)':>10s}{'전방모델 호출/시료':>18s}")
for k, v in mA.items():
    print(f"{k:26s}{v['succ']:10d}{np.mean(v['calls']):18.0f}")

# --- 손익분기 -------------------------------------------------------------
print("\n일회성 비용 대 시료당 절감 (전방모델 호출 기준)")
nn = np.mean(methods["회귀망 → LM"]["calls"])
for k in ("Simulated Annealing", "Basin-Hopping"):
    per = np.mean(methods[k]["calls"])
    save = per - nn
    if save > 0:
        print(f"  vs {k:22s} 시료당 {save:6.0f}회 절감 → 손익분기 {gen_calls/save:6.1f} 시료")
print(f"  (일회성: 학습데이터 생성 {gen_calls}회 + Adam 학습 {train_s:.0f}s)")

print("\n벽시계 기준 손익분기 (Adam 학습시간까지 포함)")
nn_t = np.mean(methods["회귀망 → LM"]["time"])
one_time_ms = (gen_s + train_s) * 1e3
for k in ("Simulated Annealing", "Basin-Hopping"):
    per_t = np.mean(methods[k]["time"])
    save_t = (per_t - nn_t) * 1e3
    if save_t > 0:
        print(f"  vs {k:22s} 시료당 {save_t:6.2f} ms 절감 → 손익분기 {one_time_ms/save_t:8.0f} 시료")
print("  → 전방모델이 싼 이 예제에서는 학습시간이 지배한다. 전방모델이 비싼 문제"
      "(다층 RCWA/TMM 등)에서는 위의 '호출 수' 기준이 지배하며 훨씬 빨리 회수된다.")

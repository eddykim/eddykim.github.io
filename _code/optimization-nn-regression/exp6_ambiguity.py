"""실험 6(재설계) — 모호한 입력에서 MSE 회귀망이 조건부 평균을 내는가.

1차 시도(550~650nm, 1.22 프린지)는 모호성이 부족해 편차가 회귀망 자체 산포에
묻혔다. 여기서는 밴드를 훨씬 좁혀 차수 모호성을 실제로 만든다.

수직입사에서 반사율은 두께에 대해 주기 lambda/(2*n1) ~ 205nm 로 거의 주기적이다.
밴드가 좁아 프린지가 1개도 안 들어오면 두께는 그 주기의 배수만큼 구분되지 않는다.
학습 구간 [1200,1800] 에는 약 3주기가 들어가므로 3중 모호성이 생긴다.

Bishop 1994 의 예측: MSE 최소화는 조건부 평균 <t|x> 로 수렴하므로,
모호성이 비대칭인 가장자리에서는 출력이 학습 구간 중심(1500nm)쪽으로 끌린다.
"""
import warnings
warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

import numpy as np

from common import TRUE_D, NOISE_STD, D_RANGE, Scaler, train_net, predict_d, rmse_nm
from train_data import make_dataset, make_spectra
from reflectance_model import N_SIO2

LAYERS, EPOCHS = [300, 256, 1], 600
PERIOD = 600.0 / (2 * N_SIO2)   # 중심파장 600nm 에서의 차수 주기 ~205nm


def run(lo, hi, noise, tag, quiet=False):
    wl = np.linspace(lo, hi, 300)
    cycles = 4 * np.pi * N_SIO2 * TRUE_D * (1 / lo - 1 / hi) / (2 * np.pi)
    rng = np.random.default_rng(0)
    X, d = make_dataset(8000, D_RANGE, wl, rng, noise_std=noise)
    Xv, dv = make_dataset(2000, D_RANGE, wl, np.random.default_rng(99), noise_std=noise)
    sc = Scaler(*D_RANGE).fit_inputs(X)
    net, hist, tt = train_net(X, d, sc, LAYERS, np.random.default_rng(1),
                              epochs=EPOCHS, lr=1e-3, X_val=Xv, d_val=dv)
    vr = rmse_nm(net, sc, Xv, dv)
    if not quiet:
        print(f"\n{tag}: {lo:.0f}~{hi:.0f}nm, 노이즈 {noise}, 밴드 내 {cycles:.2f} 프린지"
              f"  → val RMSE {vr:.1f} nm")
        print(f"  {'참 두께':>8s}{'예측 평균':>11s}{'편차':>9s}{'중심(1500)쪽':>13s}")
    curve_t, curve_p = [], []
    n_pull, n_tot = 0, 0
    for dt in (1250, 1350, 1450, 1500, 1550, 1650, 1750):
        ps = [float(predict_d(net, sc, make_spectra([float(dt)], wl)[0] +
                    np.random.default_rng(900 + s).normal(0, noise, wl.size))[0])
              for s in range(30)]
        m = float(np.mean(ps))
        curve_t.append(float(dt)); curve_p.append(m)
        if dt != 1500:
            n_tot += 1
            pulled = (m - dt) * (1500 - dt) > 0
            n_pull += pulled
            mark = "예" if pulled else "아니오"
        else:
            mark = "-"
        if not quiet:
            print(f"  {dt:8d}{m:11.1f}{m-dt:9.1f}{mark:>13s}")
    if not quiet:
        print(f"  → 중심 쪽으로 끌린 경우: {n_pull}/{n_tot}")
    return vr, np.array(curve_t), np.array(curve_p)


if __name__ == "__main__":
    print(f"두께 차수 주기(중심 600nm) ≈ {PERIOD:.0f} nm — 학습구간 600nm 에 약 3주기")
    run(450, 750, NOISE_STD, "[기준] 넓은 밴드")
    run(550, 650, NOISE_STD, "[1차 시도] 좁은 밴드")
    run(595, 605, NOISE_STD, "[재설계] 아주 좁은 밴드")
    run(595, 605, 0.02,      "[재설계+노이즈] 아주 좁은 밴드 + 5배 노이즈")

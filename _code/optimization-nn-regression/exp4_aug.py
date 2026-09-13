"""실험 4(재설계) — augmentation 이 언제 도움이 되는가.

1차 시도의 결함 두 가지를 고쳤다.
 (1) 증강군이 서로 다른 두께를 1,000개만 봤다(기준선은 8,000개). 두께 다양성과
     증강 효과가 섞여 있었다 → 동일 두께 개수로 맞춘다.
 (2) 시험 조건에 파장축·광량 오프셋이 없었다. 증강이 대비하는 오차 모드가
     시험에 없으니 이득이 나올 수 없다 → 그 오차 모드를 시험에 넣는다.

두 체제를 비교한다.
 A. 데이터 풍부 : 두께 8,000개 (전방 모델이 공짜인 우리 상황)
 B. 데이터 희소 : 두께 250개  (Kwak et al. 2021 의 실제 상황 — 양산 라인에서
                  정상 시료 125개뿐이라 시료당 40개씩 증강했다)
"""
import warnings
warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

import numpy as np

from common import WL, TRUE_D, NOISE_STD, D_RANGE, Scaler, train_net, predict_d, rmse_nm
from train_data import make_spectra, augment, measured_reference
from reflectance_model import reflectance

LAYERS, EPOCHS = [300, 256, 1], 600
N_TOTAL = 8000
AUG_KW = dict(noise_std_range=(0.0, 0.02), vertical_offset=0.04,
              vertical_scale=0.05, lateral_shift_nm=6.0)


def make_train(n_distinct, n_per, augmented, seed=0):
    rng = np.random.default_rng(seed)
    base = rng.uniform(*D_RANGE, size=n_distinct)
    if augmented:
        return augment(base, WL, rng, n_per_sample=n_per, **AUG_KW)
    X = np.repeat(make_spectra(base, WL), n_per, axis=0)
    d = np.repeat(base, n_per)
    return X, d


def build(n_distinct, n_per, augmented, seed=0):
    X, d = make_train(n_distinct, n_per, augmented, seed)
    vr_rng = np.random.default_rng(99)
    vbase = vr_rng.uniform(*D_RANGE, size=500)
    # 검증셋은 항상 '실측을 닮은' 분포로 둔다 — early stopping 기준을 통일한다.
    Xv, dv = augment(vbase, WL, vr_rng, n_per_sample=2, **AUG_KW)
    sc = Scaler(*D_RANGE).fit_inputs(X)
    net, hist, tt = train_net(X, d, sc, LAYERS, np.random.default_rng(seed + 1),
                              epochs=EPOCHS, lr=1e-3, X_val=Xv, d_val=dv)
    return net, sc, rmse_nm(net, sc, Xv, dv)


# ---- 시험 조건: 조건마다 무작위 두께 200개로 RMSE 를 낸다 ------------------
# 회귀망은 시료 하나에서의 편차가 크게 흔들리므로(실험 1-4) 단일 시료 평가는 신뢰할 수 없다.
N_EVAL = 200


def test_sets(seed=4242):
    """조건 이름 -> (X (N,300), d_true (N,))"""
    rng = np.random.default_rng(seed)
    d = rng.uniform(1250, 1750, size=N_EVAL)
    g = lambda sd, k: np.random.default_rng(k).normal(0, sd, size=(N_EVAL, WL.size))
    T = {}
    T["정상"] = (make_spectra(d, WL) + g(NOISE_STD, 1), d)
    T["노이즈 5배"] = (make_spectra(d, WL) + g(0.02, 2), d)
    T["모델오차 n1=1.50"] = (make_spectra(d, WL, n1=1.50) + g(NOISE_STD, 3), d)
    # 파장축이 +4nm 어긋난 분광기
    T["파장축 +4nm"] = (np.stack([reflectance(x, WL + 4.0) for x in d]) + g(NOISE_STD, 4), d)
    # 광량 스케일 1.03, 베이스라인 +0.02
    T["광량 x1.03 +0.02"] = (1.03 * make_spectra(d, WL) + 0.02 + g(NOISE_STD, 5), d)
    return T


def evaluate(net, sc, T):
    """조건별 RMSE(nm)"""
    return [float(np.sqrt(np.mean((predict_d(net, sc, X) - dt) ** 2)))
            for X, dt in T.values()]


if __name__ == "__main__":
    T = test_sets()
    names = list(T)
    regimes = [
        ("A. 풍부 8000개, 증강 없음", 8000, 1, False),
        ("A. 풍부 8000개, 증강",      8000, 1, True),
        ("B. 희소  250개, 증강 없음",  250, 32, False),
        ("B. 희소  250개, 증강",       250, 32, True),
    ]
    print("표의 수치 = 조건별 RMSE (nm), 무작위 두께 200개 평균. 두께 개수와 총 표본수(8000)를 맞췄다.\n")
    print(f"{'훈련':26s}" + "".join(f"{n:>17s}" for n in names) + f"{'valRMSE':>10s}")
    rows = {}
    for label, nd, npr, aug in regimes:
        net, sc, vr = build(nd, npr, aug)
        devs = evaluate(net, sc, T)
        rows[label] = devs
        print(f"{label:26s}" + "".join(f"{v:17.2f}" for v in devs) + f"{vr:10.2f}", flush=True)

    print("\n증강 효과 (RMSE 차이, 음수면 증강이 개선):")
    print(f"{'체제':16s}" + "".join(f"{n:>17s}" for n in names))
    for reg, a, b in (("A. 풍부", "A. 풍부 8000개, 증강 없음", "A. 풍부 8000개, 증강"),
                      ("B. 희소", "B. 희소  250개, 증강 없음", "B. 희소  250개, 증강")):
        diff = [y - x for x, y in zip(rows[a], rows[b])]
        print(f"{reg:16s}" + "".join(f"{v:17.2f}" for v in diff))

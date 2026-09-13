"""최적화 5편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/optimization-nn-regression/      (한국어)
      ../../assets/img/posts/optimization-nn-regression/en/   (영문)

측정 데이터는 1~4편과 동일(SiO2/Si 단층, 두께 1490nm, seed=0 노이즈).
계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은 데이터가
완전히 동일하고 표기만 다르다. 수치는 run_experiments.py / exp4_aug.py /
exp5_angles.py / exp6_ambiguity.py / exp7_vs_heuristics.py 의 결과와 같다.
캐시: 계산 결과를 figdata.npz 에 저장한다. 그림만 손볼 때는 캐시를 재사용해
즉시 다시 그린다. 강제 재계산은 `python generate_figures.py --recompute`.
"""
import os
import sys
import time
import warnings

warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

import matplotlib.pyplot as plt
import numpy as np

import levenberg_marquardt as lmmod
import train_data as tdmod
from common import (WL, TRUE_D, NOISE_STD, D_RANGE, success, lm_fit,
                    Scaler, train_net, predict_d, rmse_nm, objective)
from train_data import make_dataset, make_spectra, augment, measured_reference
from reflectance_model import N_SIO2, reflectance
from simulated_annealing import simulated_annealing
from basin_hopping import basin_hopping
from exp5_angles import build_grid, corr_cond, lm_multi, measure as measure5
from exp4_aug import build as build_aug, test_sets as tests4, evaluate as eval4, AUG_KW
from exp6_ambiguity import run as run_ambiguity

# 한글 텍스트에만 한글 폰트를 지정한다. 시리즈의 다른 편과 같은 AppleGothic 을 쓴다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                        "assets", "img", "posts", "optimization-nn-regression")

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT, "tickfont": "AppleGothic",
        "f1_title": "그림1. 성공률과 시료당 비용 (무작위 시작 100회)",
        "f1_methods": ["단일 LM", "Simulated\nAnnealing", "Basin-\nHopping", "회귀망\n→ LM"],
        "f1_y1": "전역 최솟값 도달 성공률 (%)", "f1_y2": "전방모델 호출 / 시료",
        "f2_title": "그림2. 잔차 재계산이 알려주는 것과 못 알려주는 것",
        "f2_cases": ["정상\n(기준선)", "학습범위 밖\n(2000nm)", "노이즈 5배\n(std=0.02)", "모델 오차\n(n1=1.50)"],
        "f2_y": "기준선 대비 잔차 J 배수 (log)", "f2_dev": "두께 편차",
        "f3_title": "그림3. 증강이 고치는 오차와 못 고치는 오차",
        "f3_conds": ["정상", "노이즈\n5배", "모델오차\nn1=1.50", "파장축\n+4nm", "광량\nx1.03+0.02"],
        "f3_y": "예측 RMSE (nm)", "f3_no": "증강 없음", "f3_yes": "증강",
        "f3_rich": "학습 두께 8000개", "f3_poor": "학습 두께 250개",
        "f4_title": "그림4. 모호한 입력에서 회귀망은 조건부 평균으로 간다",
        "f4_x": "참 두께 (nm)", "f4_y": "회귀망 예측 (nm)",
        "f4_ideal": "이상적 (예측 = 참값)", "f4_wide": "450~750nm (프린지 3.9개)",
        "f4_narrow": "595~605nm + 5배 노이즈 (프린지 0.12개)", "f4_center": "학습구간 중심 1500nm",
        "f5_title": "그림5. 같은 300점, 측정의 종류만 바꿨을 때 (d, n1) 추정값 산포",
        "f5_x": "두께 d 편차 (nm)", "f5_y": "굴절률 n1 편차",
        "f5_a": "수직입사만 450~750nm", "f5_b": "수직입사만 450~1000nm",
        "f5_c": "0°+60° 각 150점", "f5_corr": "상관",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "tickfont": None,
        "f1_title": "Fig 1. Success rate and per-sample cost (100 random starts)",
        "f1_methods": ["single LM", "Simulated\nAnnealing", "Basin-\nHopping", "network\n→ LM"],
        "f1_y1": "Rate of reaching the global minimum (%)", "f1_y2": "Forward-model calls / sample",
        "f2_title": "Fig 2. What recomputing the residual does and does not reveal",
        "f2_cases": ["normal\n(baseline)", "outside training\n(2000nm)", "5x noise\n(std=0.02)", "model error\n(n1=1.50)"],
        "f2_y": "Residual J relative to baseline (log)", "f2_dev": "thickness error",
        "f3_title": "Fig 3. Errors augmentation fixes, and errors it does not",
        "f3_conds": ["normal", "5x\nnoise", "model err\nn1=1.50", "wavelength\n+4nm", "photometric\nx1.03+0.02"],
        "f3_y": "Prediction RMSE (nm)", "f3_no": "no augmentation", "f3_yes": "augmented",
        "f3_rich": "8000 training thicknesses", "f3_poor": "250 training thicknesses",
        "f4_title": "Fig 4. On ambiguous input the network returns the conditional mean",
        "f4_x": "True thickness (nm)", "f4_y": "Network prediction (nm)",
        "f4_ideal": "ideal (prediction = truth)", "f4_wide": "450-750nm (3.9 fringes)",
        "f4_narrow": "595-605nm + 5x noise (0.12 fringes)", "f4_center": "training-range centre 1500nm",
        "f5_title": "Fig 5. (d, n1) estimate scatter - same 300 points, different kinds of measurement",
        "f5_x": "Thickness error d (nm)", "f5_y": "Refractive index error n1",
        "f5_a": "normal incidence 450-750nm", "f5_b": "normal incidence 450-1000nm",
        "f5_c": "0deg+60deg, 150 points each", "f5_corr": "corr",
    },
}

# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------
CACHE = os.path.join(os.path.dirname(__file__), "figdata.npz")
RECOMPUTE = "--recompute" in sys.argv


def compute():
    print("계산 시작 — 신경망 6개를 학습하므로 몇 분 걸린다", flush=True)
    T0 = time.perf_counter()

    LAYERS, EPOCHS = [300, 256, 1], 600


    def build(wl=WL, noise=0.0, aug=None, n_distinct=8000, n_per=1, seed=0):
        rng = np.random.default_rng(seed)
        base = rng.uniform(*D_RANGE, size=n_distinct)
        if aug is None:
            X = np.repeat(make_spectra(base, wl), n_per, axis=0)
            d = np.repeat(base, n_per)
            if noise > 0:
                X = X + rng.normal(0, noise, X.shape)
        else:
            X, d = augment(base, wl, rng, n_per_sample=n_per, **aug)
        vr = np.random.default_rng(99)
        vb = vr.uniform(*D_RANGE, size=500)
        if aug is None:
            Xv, dv = make_dataset(2000, D_RANGE, wl, vr, noise_std=noise)
        else:
            Xv, dv = augment(vb, wl, vr, n_per_sample=2, **aug)
        sc = Scaler(*D_RANGE).fit_inputs(X)
        net, _, _ = train_net(X, d, sc, LAYERS, np.random.default_rng(seed + 1),
                              epochs=EPOCHS, lr=1e-3, X_val=Xv, d_val=dv)
        return net, sc


    NET, SC = build()
    print(f"  기준 회귀망 학습 완료 ({time.perf_counter()-T0:.0f}s)", flush=True)

    # --- 그림1: 성공률과 비용 (실험 7) -----------------------------------------
    _orig = lmmod.reflectance
    CNT = {"n": 0}


    def _count(*a, **k):
        CNT["n"] += 1
        return _orig(*a, **k)


    lmmod.reflectance = _count


    def timed(fn):
        CNT["n"] = 0
        out = fn()
        return out, CNT["n"]


    N = 100
    starts = np.random.default_rng(2024).uniform(*D_RANGE, size=N)
    f1_succ, f1_calls = [0, 0, 0, 0], [[], [], [], []]
    for i in range(N):
        R = measured_reference(WL, seed=i)
        d0 = starts[i]
        (r, c) = timed(lambda: lm_fit(d0, WL, R));                      f1_succ[0] += success(r[0]); f1_calls[0].append(c)
        (r, c) = timed(lambda: simulated_annealing(d0, WL, R, n_iter=300, T0=0.02,
                       cooling=0.97, step_sigma=80.0, rng=np.random.default_rng(2000 + i)))
        f1_succ[1] += success(r[3]); f1_calls[1].append(c)
        (r, c) = timed(lambda: basin_hopping(d0, WL, R, n_hops=15, T0=0.05, perturb_sigma=150.0,
                       lm_kwargs=dict(tau=1e-3, n_iter=30), rng=np.random.default_rng(1000 + i)))
        f1_succ[2] += success(r[2]); f1_calls[2].append(c)
        (r, c) = timed(lambda: lm_fit(float(predict_d(NET, SC, R)[0]), WL, R))
        f1_succ[3] += success(r[0]); f1_calls[3].append(c)
    lmmod.reflectance = _orig
    f1_calls = [float(np.mean(c)) for c in f1_calls]
    print(f"  그림1 계산 완료: 성공률 {f1_succ}, 호출 {[round(c) for c in f1_calls]}", flush=True)

    # --- 그림2: 잔차 재계산 (실험 3) -------------------------------------------
    R_base = measured_reference(WL)
    d_base = float(predict_d(NET, SC, R_base)[0])
    J_base = objective(d_base, WL, R_base)
    f2_cases = [
        (R_base, TRUE_D),
        (make_spectra([2000.0], WL)[0] + np.random.default_rng(11).normal(0, NOISE_STD, WL.size), 2000.0),
        (make_spectra([TRUE_D], WL)[0] + np.random.default_rng(12).normal(0, 0.02, WL.size), TRUE_D),
        (make_spectra([TRUE_D], WL, n1=1.50)[0] + np.random.default_rng(13).normal(0, NOISE_STD, WL.size), TRUE_D),
    ]
    f2_ratio, f2_dev = [], []
    for R, dt in f2_cases:
        dn = float(predict_d(NET, SC, R)[0])
        f2_ratio.append(objective(dn, WL, R) / J_base)
        f2_dev.append(dn - dt)
    print(f"  그림2 계산 완료: 배수 {[round(x,1) for x in f2_ratio]}, 편차 {[round(x,1) for x in f2_dev]}", flush=True)

    # --- 그림3: 증강 (실험 4) ---------------------------------------------------
    # exp4_aug.py 의 build/test_sets/evaluate 를 그대로 쓴다. 여기서 학습 파이프라인을
    # 다시 구현하면 early stopping 검증셋 같은 조건이 어긋나 본문 표와 다른 수치가 나온다.
    TESTS4 = tests4()
    f3 = {}
    for tag, nd, npr, aug in (("rich_no", 8000, 1, False), ("rich_aug", 8000, 1, True),
                              ("poor_no", 250, 32, False), ("poor_aug", 250, 32, True)):
        net, sc, _ = build_aug(nd, npr, aug)
        f3[tag] = eval4(net, sc, TESTS4)
        print(f"  그림3 {tag}: {[round(v,2) for v in f3[tag]]}", flush=True)

    # --- 그림4: 조건부 평균 (실험 6) --------------------------------------------
    # exp6_ambiguity.run() 을 그대로 호출한다(본문 표와 같은 수치를 보장).
    _, f4_true, f4_wide = run_ambiguity(450, 750, NOISE_STD, "", quiet=True)
    _, _, f4_narrow = run_ambiguity(595, 605, 0.02, "", quiet=True)
    print(f"  그림4 계산 완료: 좁은밴드 {[round(v,1) for v in f4_narrow]}", flush=True)

    # --- 그림5: (d, n1) 신뢰타원 (실험 5-3) -------------------------------------
    f5 = []
    for grid in (build_grid([0.0]), build_grid([0.0], hi=1000.0), build_grid([0.0, 60.0])):
        corr, _ = corr_cond([TRUE_D, N_SIO2], grid)
        ds, ns = [], []
        for s in range(120):
            y = measure5(grid, seed=300 + s)
            p = lm_multi([TRUE_D, N_SIO2], grid, y)
            ds.append(p[0] - TRUE_D); ns.append(p[1] - N_SIO2)
        f5.append((np.array(ds), np.array(ns), corr))
    print(f"  그림5 계산 완료: 상관 {[round(c,4) for _,_,c in f5]}", flush=True)
    print(f"계산 총 {time.perf_counter()-T0:.0f}s", flush=True)
    np.savez(CACHE,
             f1_succ=np.array(f1_succ), f1_calls=np.array(f1_calls),
             f2_ratio=np.array(f2_ratio), f2_dev=np.array(f2_dev),
             f3_keys=np.array(list(f3.keys())),
             f3_vals=np.array([f3[k] for k in f3]),
             f4_true=f4_true, f4_wide=np.array(f4_wide), f4_narrow=np.array(f4_narrow),
             f5_ds=np.array([a for a, _, _ in f5]),
             f5_ns=np.array([b for _, b, _ in f5]),
             f5_corr=np.array([c for _, _, c in f5]))
    return (f1_succ, f1_calls, f2_ratio, f2_dev, f3, f4_true, f4_wide, f4_narrow, f5)


if RECOMPUTE or not os.path.exists(CACHE):
    (f1_succ, f1_calls, f2_ratio, f2_dev, f3,
     f4_true, f4_wide, f4_narrow, f5) = compute()
else:
    print(f"캐시 사용: {os.path.basename(CACHE)}  (재계산은 --recompute)")
    z = np.load(CACHE, allow_pickle=False)
    f1_succ, f1_calls = list(z["f1_succ"]), list(z["f1_calls"])
    f2_ratio, f2_dev = list(z["f2_ratio"]), list(z["f2_dev"])
    f3 = {k: list(v) for k, v in zip(z["f3_keys"], z["f3_vals"])}
    f4_true, f4_wide, f4_narrow = z["f4_true"], z["f4_wide"], z["f4_narrow"]
    f5 = [(z["f5_ds"][i], z["f5_ns"][i], float(z["f5_corr"][i])) for i in range(3)]



# ---------------------------------------------------------------------------
def render(L):
    out = os.path.normpath(L["dir"])
    os.makedirs(out, exist_ok=True)
    tf = {"fontfamily": L["tickfont"]} if L["tickfont"] else {}

    # 그림1 — 성공률(막대) + 비용(로그 점)
    fig, ax1 = plt.subplots(figsize=(6, 4))
    x = np.arange(4)
    ax1.bar(x, f1_succ, color=["#c0c0c0", "#8fb8de", "#5b8fc9", "#e08a3c"], width=0.6)
    ax1.set_ylim(0, 112)
    ax1.set_ylabel(L["f1_y1"], **L["font"])
    ax1.set_xticks(x)
    ax1.set_xticklabels(L["f1_methods"], **tf)
    for i, v in enumerate(f1_succ):
        ax1.text(i, v - 6, f"{v}%", ha="center", fontsize=10, color="white",
                 fontweight="bold")
    ax2 = ax1.twinx()
    ax2.plot(x, f1_calls, "ko--", ms=6, lw=1.2)
    ax2.set_yscale("log")
    ax2.set_ylabel(L["f1_y2"], **L["font"])
    for i, v in enumerate(f1_calls):
        ax2.annotate(f"{v:.0f}", (i, v), textcoords="offset points",
                     xytext=(0, 9), ha="center", fontsize=9,
                     bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    ax1.set_title(L["f1_title"], **L["font"])
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig1-success-vs-cost.png"), dpi=150)
    plt.close()

    # 그림2 — 잔차 배수(로그 막대) + 편차 주석
    plt.figure(figsize=(6, 4))
    cols = ["#c0c0c0", "#d95f5f", "#8fb8de", "#e08a3c"]
    plt.bar(np.arange(4), f2_ratio, color=cols, width=0.6)
    plt.yscale("log")
    plt.ylim(0.35, 4e3)
    plt.axhline(1.0, color="k", lw=0.8, ls=":")
    plt.ylabel(L["f2_y"], **L["font"])
    plt.xticks(np.arange(4), L["f2_cases"], **tf)
    for i, (r, d) in enumerate(zip(f2_ratio, f2_dev)):
        plt.annotate(f"x{r:,.0f}" if r >= 10 else f"x{r:.1f}",
                     (i, r), textcoords="offset points", xytext=(0, 6),
                     ha="center", fontsize=9)
        plt.annotate(f"{L['f2_dev']} {d:+.1f} nm", (i, r), textcoords="offset points",
                     xytext=(0, -14), ha="center", fontsize=8, color="#444", **L["font"])
    plt.title(L["f2_title"], **L["font"])
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig2-residual-alarm.png"), dpi=150)
    plt.close()

    # 그림3 — 증강 유무 × 풍부/희소
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), sharey=True)
    x = np.arange(5)
    for ax, (no, yes, sub) in zip(axes, (("rich_no", "rich_aug", L["f3_rich"]),
                                          ("poor_no", "poor_aug", L["f3_poor"]))):
        ax.bar(x - 0.2, f3[no], width=0.4, color="#c0c0c0", label=L["f3_no"])
        ax.bar(x + 0.2, f3[yes], width=0.4, color="#e08a3c", label=L["f3_yes"])
        ax.set_xticks(x)
        ax.set_xticklabels(L["f3_conds"], fontsize=8, **tf)
        ax.set_title(sub, fontsize=10, **L["font"])
        ax.set_yscale("log")
    axes[0].set_ylabel(L["f3_y"], **L["font"])
    axes[0].legend(prop=L["legend"], fontsize=9)
    fig.suptitle(L["f3_title"], **L["font"])
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig3-augmentation-rmse.png"), dpi=150)
    plt.close()

    # 그림4 — 예측 vs 참값
    plt.figure(figsize=(6, 4))
    lim = [1200, 1800]
    plt.plot(lim, lim, "k:", lw=1, label=L["f4_ideal"])
    plt.axhline(1500, color="#999", lw=0.8, ls="--")
    plt.plot(f4_true, f4_wide, "o-", color="#5b8fc9", ms=5, label=L["f4_wide"])
    plt.plot(f4_true, f4_narrow, "s-", color="#d95f5f", ms=5, label=L["f4_narrow"])
    plt.annotate(L["f4_center"], (1210, 1505), fontsize=8, color="#666", **L["font"])
    plt.xlabel(L["f4_x"], **L["font"])
    plt.ylabel(L["f4_y"], **L["font"])
    plt.xlim(lim); plt.ylim(lim)
    plt.legend(prop=L["legend"], fontsize=8, loc="lower right")
    plt.title(L["f4_title"], fontsize=10, **L["font"])
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig4-conditional-mean.png"), dpi=150)
    plt.close()

    # 그림5 — (d, n1) 산점
    plt.figure(figsize=(6, 4))
    styles = [("#c0c0c0", "o", L["f5_a"]), ("#5b8fc9", "^", L["f5_b"]), ("#e08a3c", "s", L["f5_c"])]
    for (ds, ns, corr), (c, m, lab) in zip(f5, styles):
        plt.scatter(ds, ns, s=16, c=c, marker=m, alpha=0.75,
                    label=f"{lab}  ({L['f5_corr']} {corr:.3f})")
    plt.axhline(0, color="k", lw=0.6); plt.axvline(0, color="k", lw=0.6)
    plt.xlabel(L["f5_x"], **L["font"])
    plt.ylabel(L["f5_y"], **L["font"])
    plt.legend(prop=L["legend"], fontsize=8)
    plt.title(L["f5_title"], fontsize=10, **L["font"])
    plt.tight_layout()
    plt.savefig(os.path.join(out, "fig5-correlation-ellipse.png"), dpi=150)
    plt.close()


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)
print("그림 5개 x 2언어 생성 완료")

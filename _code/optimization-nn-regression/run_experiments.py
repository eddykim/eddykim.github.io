"""최적화 5편 실험 1~6 실행.

실행: python run_experiments.py [1 2 3 4 5 6]   (인자 없으면 전부)
결과는 표준출력에 표로 찍고, results.json 에 수치를 저장한다.

측정 데이터 규약은 1~4편과 동일(WL 450~750nm 300점, d=1490nm, noise_std=0.004, seed=0).
학습 구성은 사전 스윕으로 고른 값이다 — 8000표본 / [300,256,1] / 600 epoch / Adam lr=1e-3 /
검증 오차 최소 epoch 으로 되돌리는 early stopping(Fried & Masa 1994 의 권고).
데이터·epoch 을 2배로 늘려도 val RMSE 가 거의 변하지 않아 고원으로 판단했다.
"""
import json
import sys
import time

import numpy as np

from common import (WL, TRUE_D, NOISE_STD, D_RANGE, HARD_D0, SUCCESS_TOL,
                    success, lm_fit, train_net, predict_d, rmse_nm,
                    MLP, Scaler, make_dataset, make_spectra, augment,
                    measured_reference, objective, reflectance, N_SIO2)
from lm2 import lm2, objective2, correlation_and_condition

N_TRAIN, N_VAL = 8000, 2000
LAYERS = [300, 256, 1]
EPOCHS, LR = 600, 1e-3
RES = {}


def hdr(t):
    print("\n" + "=" * 78 + f"\n{t}\n" + "=" * 78, flush=True)


def build_net(noise_std=0.0, layers=LAYERS, epochs=EPOCHS, wl=WL, seed=0,
              augment_kw=None, n_train=N_TRAIN, d_range=D_RANGE, n1=N_SIO2):
    """학습 데이터 생성 → 스케일러 → 학습. 반환: (net, scaler, hist, 학습시간, val RMSE)"""
    rng = np.random.default_rng(seed)
    if augment_kw is None:
        X, d = make_dataset(n_train, d_range, wl, rng, noise_std=noise_std, n1=n1)
    else:
        base = rng.uniform(d_range[0], d_range[1], size=n_train // augment_kw.get("n_per_sample", 1))
        X, d = augment(base, wl, rng, n1=n1, **augment_kw)
    if augment_kw is None:
        Xv, dv = make_dataset(N_VAL, d_range, wl, np.random.default_rng(99),
                              noise_std=noise_std, n1=n1)
    else:
        # 증강 훈련이면 검증셋도 같은 증강 분포에서 뽑는다.
        # 깨끗한 이론 스펙트럼으로 검증하면 early stopping 이 "깨끗한 데이터에
        # 가장 좋은" 체크포인트를 골라 증강 효과를 상쇄한다.
        vrng = np.random.default_rng(99)
        nper = augment_kw.get("n_per_sample", 1)
        vbase = vrng.uniform(d_range[0], d_range[1], size=max(1, N_VAL // nper))
        Xv, dv = augment(vbase, wl, vrng, n1=n1, **augment_kw)
    sc = Scaler(*d_range).fit_inputs(X)
    net, hist, tt = train_net(X, d, sc, layers, np.random.default_rng(seed + 1),
                              epochs=epochs, lr=LR, X_val=Xv, d_val=dv)
    return net, sc, hist, tt, rmse_nm(net, sc, Xv, dv)


# ---------------------------------------------------------------- 실험 1
def exp1():
    hdr("실험 1 — 이론 데이터로 훈련한 회귀망 vs LM")
    out = {}

    print("\n[1-1] 구조 비교 (무노이즈 이론 스펙트럼 훈련, 600 epoch)")
    print(f"{'구조':24s}{'val RMSE(nm)':>14s}{'train RMSE(nm)':>16s}{'학습(s)':>10s}")
    arch_rows = []
    for layers in ([300, 128, 1], [300, 256, 1], [300, 256, 128, 1]):
        net, sc, hist, tt, vr = build_net(layers=layers)
        Xt, dt = make_dataset(2000, D_RANGE, WL, np.random.default_rng(7))
        tr = rmse_nm(net, sc, Xt, dt)
        print(f"{str(layers):24s}{vr:14.2f}{tr:16.2f}{tt:10.0f}")
        arch_rows.append(dict(layers=layers, val_rmse=vr, test_rmse=tr, train_s=tt))
    out["arch"] = arch_rows

    print("\n[1-2] 훈련 노이즈 유무 (구조 %s)" % LAYERS)
    nets = {}
    for tag, ns in (("무노이즈(이론)", 0.0), ("노이즈 0.004", NOISE_STD)):
        net, sc, hist, tt, vr = build_net(noise_std=ns)
        nets[tag] = (net, sc)
        # 과훈련 점검: 훈련/검증 곡선의 최저점
        best_ep = hist["best_epoch"]
        print(f"  {tag:16s} valRMSE={vr:6.2f} nm  학습 {tt:4.0f}s  "
              f"val최저 epoch={best_ep}/{EPOCHS-1}  "
              f"train={hist['train'][-1]:.2e} val={hist['val'][-1]:.2e}")
        out[f"regime_{tag}"] = dict(val_rmse=vr, train_s=tt, best_epoch=best_ep,
                                    final_train=hist["train"][-1], final_val=hist["val"][-1])
    RES["_nets"] = nets

    print("\n[1-3] 1490nm 측정 1회(seed=0)에 대한 예측")
    Rm = measured_reference(WL)
    t0 = time.perf_counter(); _ = predict_d(nets["무노이즈(이론)"][0], nets["무노이즈(이론)"][1], Rm); t_inf = time.perf_counter() - t0
    rows = []
    for tag, (net, sc) in nets.items():
        d_nn = float(predict_d(net, sc, Rm)[0])
        J_nn = objective(d_nn, WL, Rm)
        rows.append((f"회귀망 {tag}", d_nn, d_nn - TRUE_D, J_nn, "-"))
    for tag, d0 in (("LM d0=1480", 1480.0), ("LM d0=1300(어려움)", HARD_D0)):
        d_lm, it, J_lm = lm_fit(d0, WL, Rm)
        rows.append((tag, d_lm, d_lm - TRUE_D, J_lm, it))
    print(f"{'방법':24s}{'추정 d(nm)':>13s}{'편차(nm)':>11s}{'최종 J':>12s}{'iter':>6s}")
    for r in rows:
        print(f"{r[0]:24s}{r[1]:13.3f}{r[2]:11.3f}{r[3]:12.3e}{str(r[4]):>6s}")
    out["single"] = [dict(method=r[0], d=r[1], dev=r[2], J=r[3], it=str(r[4])) for r in rows]

    print("\n[1-4] 정밀도 — 같은 1490nm 시료, 노이즈 seed 0~49")
    preds = {k: [] for k in list(nets) + ["LM d0=1480"]}
    t_lm = []
    for s in range(50):
        R = measured_reference(WL, seed=s)
        for tag, (net, sc) in nets.items():
            preds[tag].append(float(predict_d(net, sc, R)[0]))
        t1 = time.perf_counter(); d_lm, _, _ = lm_fit(1480.0, WL, R); t_lm.append(time.perf_counter() - t1)
        preds["LM d0=1480"].append(d_lm)
    print(f"{'방법':24s}{'평균 d':>12s}{'정확도(편차)':>14s}{'정밀도(std)':>13s}")
    for tag, v in preds.items():
        v = np.array(v)
        print(f"{tag:24s}{v.mean():12.3f}{v.mean()-TRUE_D:14.3f}{v.std(ddof=1):13.4f}")
        out.setdefault("precision", {})[tag] = dict(mean=float(v.mean()),
                                                    bias=float(v.mean() - TRUE_D),
                                                    std=float(v.std(ddof=1)))
    out["timing"] = dict(nn_inference_s=t_inf, lm_mean_s=float(np.mean(t_lm)))
    print(f"\n  회귀망 추론 1회 {t_inf*1e3:.3f} ms   |   LM 1회 평균 {np.mean(t_lm)*1e3:.1f} ms"
          f"   (회귀망 학습 1회는 위 [1-2] 참고)")
    RES["exp1"] = out


# ---------------------------------------------------------------- 실험 2
def exp2():
    hdr("실험 2 — 회귀망 출력을 LM 초기값으로 (Fried & Masa 1994 Table I 구조)")
    net, sc = RES["_nets"]["무노이즈(이론)"]
    out = {}

    def four_conditions(R, rng):
        d_nn = float(predict_d(net, sc, R)[0])
        d_b, it_b, J_b = lm_fit(d_nn, WL, R)
        d0r = rng.uniform(*D_RANGE)
        d_c, it_c, J_c = lm_fit(d0r, WL, R)
        d_d, it_d, J_d = lm_fit(HARD_D0, WL, R)
        return dict(nn=(d_nn, 0, objective(d_nn, WL, R)),
                    nn_lm=(d_b, it_b, J_b),
                    rand_lm=(d_c, it_c, J_c),
                    hard_lm=(d_d, it_d, J_d))

    for label, gen in (("[2-A] 참 두께 1490nm 고정, 노이즈 seed 0~99",
                        lambda i: (measured_reference(WL, seed=i), TRUE_D)),
                       ("[2-B] 참 두께 무작위 [1250,1750], 100회",
                        None)):
        print("\n" + label)
        rng = np.random.default_rng(2024)
        acc = {k: {"succ": 0, "dev": [], "it": []} for k in ("nn", "nn_lm", "rand_lm", "hard_lm")}
        for i in range(100):
            if gen is not None:
                R, dtrue = gen(i)
            else:
                dtrue = rng.uniform(1250, 1750)
                R = make_spectra([dtrue], WL)[0] + np.random.default_rng(5000 + i).normal(0, NOISE_STD, WL.size)
            res = four_conditions(R, rng)
            for k, (d, it, J) in res.items():
                acc[k]["succ"] += abs(d - dtrue) < SUCCESS_TOL
                acc[k]["dev"].append(abs(d - dtrue))
                acc[k]["it"].append(it)
        names = {"nn": "(a) 회귀망 단독", "nn_lm": "(b) 회귀망 → LM",
                 "rand_lm": "(c) 무작위 초기값 → LM", "hard_lm": "(d) d0=1300 → LM"}
        print(f"{'조건':26s}{'성공률(%)':>11s}{'중앙 |편차|(nm)':>17s}{'평균 iter':>11s}")
        for k in ("nn", "nn_lm", "rand_lm", "hard_lm"):
            a = acc[k]
            print(f"{names[k]:26s}{a['succ']:11d}{np.median(a['dev']):17.3f}{np.mean(a['it']):11.1f}")
            out.setdefault(label, {})[k] = dict(success=a["succ"],
                                                median_dev=float(np.median(a["dev"])),
                                                mean_iter=float(np.mean(a["it"])))
    RES["exp2"] = out


# ---------------------------------------------------------------- 실험 3
def exp3():
    hdr("실험 3 — 훈련 분포와 측정 분포를 어긋나게 하기")
    net, sc = RES["_nets"]["무노이즈(이론)"]
    out = {}
    cases = [
        ("(a) 학습범위 밖 d=2000nm", make_spectra([2000.0], WL)[0] +
         np.random.default_rng(11).normal(0, NOISE_STD, WL.size), 2000.0, dict()),
        ("(b) 노이즈 불일치 std=0.02", make_spectra([TRUE_D], WL)[0] +
         np.random.default_rng(12).normal(0, 0.02, WL.size), TRUE_D, dict()),
        ("(c) 모델오차 n1=1.50", make_spectra([TRUE_D], WL, n1=1.50)[0] +
         np.random.default_rng(13).normal(0, NOISE_STD, WL.size), TRUE_D, dict(n1=1.50)),
    ]
    print(f"{'경우':26s}{'회귀망 d':>11s}{'편차':>9s}{'회귀망 J':>12s}"
          f"{'LM d':>11s}{'LM 편차':>9s}{'LM J':>12s}")
    for name, R, dtrue, _ in cases:
        d_nn = float(predict_d(net, sc, R)[0])
        J_nn = objective(d_nn, WL, R)          # 전방 모델로 잔차 재계산 = 자기 점검
        d_lm, it, J_lm = lm_fit(d_nn, WL, R)   # 회귀망 예측에서 출발한 LM
        print(f"{name:26s}{d_nn:11.2f}{d_nn-dtrue:9.2f}{J_nn:12.3e}"
              f"{d_lm:11.2f}{d_lm-dtrue:9.2f}{J_lm:12.3e}")
        out[name] = dict(d_nn=d_nn, dev_nn=d_nn - dtrue, J_nn=J_nn,
                         d_lm=d_lm, dev_lm=d_lm - dtrue, J_lm=J_lm, d_true=dtrue)
    # 정상 조건의 J 를 기준선으로
    Rm = measured_reference(WL)
    d0 = float(predict_d(net, sc, Rm)[0])
    out["기준선(정상)"] = dict(d_nn=d0, J_nn=objective(d0, WL, Rm))
    print(f"\n  기준선(정상 측정): 회귀망 d={d0:.2f}, J={objective(d0, WL, Rm):.3e}")
    print("  → J 가 기준선보다 크게 뜨면 '전방 모델로 잔차를 재계산하는' 자기 점검이 작동한 것이다.")
    RES["exp3"] = out


# ---------------------------------------------------------------- 실험 4
def exp4():
    hdr("실험 4 — 실측을 흉내 낸 augmentation (Kwak et al. 2021 의 교란 축)")
    out = {}
    base_kw = dict(n_per_sample=8, noise_std_range=(0.0, 0.02),
                   vertical_offset=0.04, vertical_scale=0.05, lateral_shift_nm=6.0)
    variants = {
        "증강 없음(이론)": None,
        "노이즈만": dict(base_kw, use_vertical=False, use_lateral=False),
        "세로만": dict(base_kw, use_noise=False, use_lateral=False),
        "가로만": dict(base_kw, use_noise=False, use_vertical=False),
        "전부": dict(base_kw),
    }
    tests = {
        "정상(seed0)": (measured_reference(WL), TRUE_D),
        "(b) std=0.02": (make_spectra([TRUE_D], WL)[0] +
                         np.random.default_rng(12).normal(0, 0.02, WL.size), TRUE_D),
        "(c) n1=1.50": (make_spectra([TRUE_D], WL, n1=1.50)[0] +
                        np.random.default_rng(13).normal(0, NOISE_STD, WL.size), TRUE_D),
    }
    print(f"{'훈련':16s}" + "".join(f"{k:>16s}" for k in tests) + f"{'val RMSE':>11s}")
    for tag, kw in variants.items():
        net, sc, hist, tt, vr = build_net(augment_kw=kw)
        devs = []
        for k, (R, dt) in tests.items():
            devs.append(float(predict_d(net, sc, R)[0]) - dt)
        print(f"{tag:16s}" + "".join(f"{v:16.2f}" for v in devs) + f"{vr:11.2f}")
        out[tag] = dict(devs={k: d for k, d in zip(tests, devs)}, val_rmse=vr)
    print("\n  (표의 수치는 참값 대비 편차 nm. val RMSE 는 각 훈련 분포의 검증셋 기준이라 직접 비교 대상이 아니다.)")
    RES["exp4"] = out


# ---------------------------------------------------------------- 실험 5
def exp5():
    hdr("실험 5 — feature 공간을 넓히면 모호성이 줄어드는가")
    out = {}
    WL_WIDE = np.linspace(450, 1000, 300)

    print("\n[5-1] 미지수 1개(d): 파장 구간 450~750 vs 450~1000 (둘 다 300점)")
    print(f"{'구간':16s}{'국소최솟값 수':>14s}{'LM 무작위시작 성공률(%)':>24s}{'회귀망 val RMSE(nm)':>21s}")
    for tag, wl in (("450~750nm", WL), ("450~1000nm", WL_WIDE)):
        R = make_spectra([TRUE_D], wl)[0] + np.random.default_rng(0).normal(0, NOISE_STD, wl.size)
        ds = np.linspace(*D_RANGE, 6001)
        Js = np.array([objective(x, wl, R) for x in ds])
        nloc = sum(1 for i in range(1, len(ds) - 1) if Js[i] < Js[i - 1] and Js[i] < Js[i + 1])
        rng = np.random.default_rng(3)
        succ = sum(success(lm_fit(rng.uniform(*D_RANGE), wl, R)[0]) for _ in range(100))
        net, sc, hist, tt, vr = build_net(wl=wl)
        print(f"{tag:16s}{nloc:14d}{succ:24d}{vr:21.2f}")
        out[f"5-1 {tag}"] = dict(n_local_min=nloc, lm_success=succ, nn_val_rmse=vr)

    print("\n[5-2] 미지수 2개(d, n1): 파라미터 상호상관과 조건수")
    print(f"{'구간':16s}{'corr(d,n1)':>14s}{'cond(JtJ)':>13s}{'LM2 성공률(%)':>15s}"
          f"{'d 산포(nm)':>12s}{'n1 산포':>11s}")
    for tag, wl in (("450~750nm", WL), ("450~1000nm", WL_WIDE)):
        R = make_spectra([TRUE_D], wl, n1=N_SIO2)[0] + np.random.default_rng(0).normal(0, NOISE_STD, wl.size)
        corr, cond = correlation_and_condition([TRUE_D, N_SIO2], wl, R)
        rng = np.random.default_rng(4)
        ok, dd, nn_ = 0, [], []
        for _ in range(60):
            p0 = [TRUE_D + rng.uniform(-40, 40), N_SIO2 + rng.uniform(-0.04, 0.04)]
            p, it, J = lm2(p0, wl, R)
            if abs(p[0] - TRUE_D) < SUCCESS_TOL and abs(p[1] - N_SIO2) < 0.02:
                ok += 1
            dd.append(p[0]); nn_.append(p[1])
        # 노이즈 실현마다 다시 풀어 산포를 본다
        d_s, n_s = [], []
        for s in range(30):
            Rs = make_spectra([TRUE_D], wl)[0] + np.random.default_rng(700 + s).normal(0, NOISE_STD, wl.size)
            p, _, _ = lm2([TRUE_D, N_SIO2], wl, Rs)
            d_s.append(p[0]); n_s.append(p[1])
        print(f"{tag:16s}{corr:14.6f}{cond:13.3e}{ok*100//60:15d}"
              f"{np.std(d_s, ddof=1):12.4f}{np.std(n_s, ddof=1):11.5f}")
        out[f"5-2 {tag}"] = dict(corr=corr, cond=cond, lm2_success_pct=ok * 100 // 60,
                                 d_std=float(np.std(d_s, ddof=1)), n1_std=float(np.std(n_s, ddof=1)))
    RES["exp5"] = out


# ---------------------------------------------------------------- 실험 6
def exp6():
    hdr("실험 6 — 모호한 입력에서 회귀망이 조건부 평균으로 가는가 (Bishop 1994)")
    out = {}
    WL_NARROW = np.linspace(550, 650, 300)   # 프린지 ~1.2개만 들어온다
    phase_cycles = 4 * np.pi * N_SIO2 * TRUE_D * (1 / 550 - 1 / 650) / (2 * np.pi)
    print(f"  좁은 구간 550~650nm: d=1490nm 에서 밴드 내 위상 변화 ≈ {phase_cycles:.2f} 프린지")
    print("  → 차수 모호성이 커진다. 학습 구간 [1200,1800] 의 중심은 1500nm.")

    net, sc, hist, tt, vr = build_net(wl=WL_NARROW, noise_std=NOISE_STD)
    print(f"  좁은 구간 학습 val RMSE = {vr:.2f} nm  (넓은 구간 대비 크게 나빠야 정상)")

    print(f"\n{'참 두께(nm)':>12s}{'회귀망 예측':>13s}{'편차':>9s}{'중심(1500)쪽으로?':>18s}")
    rows = []
    for dtrue in (1250, 1300, 1400, 1500, 1600, 1700, 1750):
        ps = []
        for s in range(20):
            R = make_spectra([float(dtrue)], WL_NARROW)[0] + \
                np.random.default_rng(900 + s).normal(0, NOISE_STD, WL_NARROW.size)
            ps.append(float(predict_d(net, sc, R)[0]))
        m = float(np.mean(ps))
        pulled = "예" if (m - dtrue) * (1500 - dtrue) > 0 else ("-" if dtrue == 1500 else "아니오")
        print(f"{dtrue:12d}{m:13.2f}{m-dtrue:9.2f}{pulled:>18s}")
        rows.append(dict(d_true=dtrue, pred=m, dev=m - dtrue, pulled_to_center=pulled))
    out["curve"] = rows
    out["val_rmse_narrow"] = vr
    print("\n  Bishop 1994 의 예측: MSE 회귀망 출력은 조건부 평균 <t|x> 다.")
    print("  모호성이 대칭인 중심부에서는 참값과 일치하고, 가장자리에서는 중심 쪽으로 끌린다.")
    RES["exp6"] = out


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3", "4", "5", "6"]
    t0 = time.perf_counter()
    if "1" in which or any(w in which for w in ("2", "3")):
        exp1()
    for w, fn in (("2", exp2), ("3", exp3), ("4", exp4), ("5", exp5), ("6", exp6)):
        if w in which:
            fn()
    RES.pop("_nets", None)
    try:
        with open("results.json") as f:
            prev = json.load(f)
    except Exception:
        prev = {}
    prev.update(RES)
    with open("results.json", "w") as f:
        json.dump(prev, f, indent=2, default=float)
    print(f"\n총 소요 {time.perf_counter()-t0:.0f}s  → results.json 저장")

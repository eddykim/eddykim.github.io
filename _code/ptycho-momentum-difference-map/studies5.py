"""5편의 표와 그림에 들어가는 수치 실험. 모두 페르마 나선 144 위치, 초기 시편은 빈 공간(1).

  compare()   핀홀·초점 어긋남·확산판 x rPIE·작은 걸음으로 시작한 rPIE·mPIE(Maiden)·mPIE(PtyLab)·DM,
              패턴 순서 5 가지, 오차와 데이터 잔차의 궤적
  factors()   확산판에서 모멘텀 없이 시편 걸음·프로브 걸음·프로브 갱신 형태를 하나씩 바꾼 2x2x2
  warmup()    확산판에서 처음 K 회만 작은 시편 걸음을 쓰고 나머지는 1 로 돌렸을 때
  drift()     핀홀 mPIE(Maiden): 세기 제어 없이 프로브 세기가 어디로 가는지, 세기 제어가 있을 때 시편 최댓값이 어디서 커지는지
  param_map() 모멘텀 eta(마찰 = 피드백) x 시편 걸음 gamma_obj, 프로브 셋, 순서 셋의 최종 오차
  stall()     확산판의 0.02 정체: 프로브 무게중심의 밀림, 무게중심 고정, 프로브를 알 때
  toy()       좁은 골짜기 2차 함수에서 경사하강과 모멘텀의 궤적

python3 studies5.py 로 돌리면 결과를 .cache_post5.pkl 에 남긴다 (generate_figures.py 가 읽는다).
"""
import os
import pickle
import time
from functools import partial
from multiprocessing import Pool

import numpy as np

from ptycho5_core import (N, disc, pinhole_probe, diffuser_probe, focused_probe, make_object, scan, object_shape,
                          illumination, intensities, rpie, dm, error, residual, align, centroid, maiden, ptylab)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache_post5.pkl")
CK = (1, 2, 5, 10, 20, 30, 50, 75, 100, 150, 200, 300)
ORDERS = (0, 1, 2, 3, 4)
PROBES = {"pinhole": lambda: (pinhole_probe(), disc(N, 14)),
          "defocus": lambda: (focused_probe(), focused_probe(defocus=0.0)),
          "diffuser": lambda: (diffuser_probe(), disc(N, 12))}
WARM = 10                                     # 처음 몇 회를 작은 걸음으로


def warm_gamma(k, small=0.25, K=WARM):
    return small if k <= K else 1.0


ALGS = {
    "rPIE": dict(gamma_obj=1.0, probe_power="auto"),
    "rPIE_warm": dict(gamma_obj=warm_gamma, probe_power="auto"),
    "mPIE": dict(gamma_obj=0.2, momentum=maiden(), probe_power="auto"),
    "mPIE_ptylab": dict(alpha=0.1, beta=0.1, gamma_obj=0.25, gamma_prb=0.25, momentum=ptylab(), probe_power="auto"),
}


def problem(case):
    P, P0 = PROBES[case]()
    pos = scan("fermat")
    shp = object_shape(pos)
    obj = make_object(shp)
    ill = illumination(pos, P, shp)
    return P, P0, pos, shp, obj, ill, ill > 0.1 * ill.max(), intensities(obj, P, pos)


def _run(args):
    case, alg, kw, order, n_iter, ck = args
    P, P0, pos, shp, obj, ill, m, I = problem(case)
    tr = {}

    def cb(k, o, p):
        if k in ck:
            tr[k] = (error(o, obj, m), residual(o, p, pos, I))

    t = time.time()
    with np.errstate(all="ignore"):
        if alg == "DM":
            o, p = dm(I, pos, P0, shp, n_iter, callback=cb)
        else:
            o, p = rpie(I, pos, P0, shp, n_iter, rng=order, callback=cb, **kw)
    return case, alg, order, tr, error(p, P), time.time() - t


def _pool(jobs, procs=14):
    with Pool(procs) as pool:
        return pool.map(_run, jobs)


def compare(n_iter=300):
    jobs = [(c, a, kw, s, n_iter, CK) for c in PROBES for a, kw in ALGS.items() for s in ORDERS]
    jobs += [(c, "DM", None, 0, n_iter, CK) for c in PROBES]          # DM 은 순서가 없다
    return {(c, a, s): dict(trace=tr, probe_err=pe, sec=dt) for c, a, s, tr, pe, dt in _pool(jobs)}


def factors(n_iter=300, orders=(0, 1, 2)):
    jobs = []
    for g_o in (1.0, 0.25):
        for g_p in (1.0, 0.25):
            for beta in (1.0, 0.1):
                kw = dict(gamma_obj=g_o, gamma_prb=g_p, beta=beta, probe_power="auto")
                jobs += [("diffuser", (g_o, g_p, beta), kw, s, n_iter, (100, n_iter)) for s in orders]
    out = {}
    for c, a, s, tr, pe, dt in _pool(jobs):
        out.setdefault(a, []).append(tr)
    return out


def warmup(n_iter=300, orders=(0, 1, 2)):
    jobs = []
    for K in (0, 3, 10, 30):
        kw = dict(gamma_obj=partial(warm_gamma, K=K) if K else 1.0, probe_power="auto")
        jobs += [("diffuser", K, kw, s, n_iter, (K or 1, n_iter)) for s in orders]
    out = {}
    for c, a, s, tr, pe, dt in _pool(jobs):
        out.setdefault(a, []).append(tr)
    return out


DRIFT = {"mPIE": dict(gamma_obj=0.2, momentum=maiden()),
         "rPIE": dict(gamma_obj=1.0),
         "mom_obj": dict(gamma_obj=0.2, momentum=dict(T=30, fric_obj=0.9, feed_obj=0.9, fric_prb=0.0, feed_prb=0.0)),
         "mom_prb": dict(gamma_obj=0.2, momentum=dict(T=30, fric_obj=0.0, feed_obj=0.0, fric_prb=0.9, feed_prb=0.9))}


def _drift(args):
    """열: 반복, 프로브 세기, 시편 진폭 최댓값, 그 자리의 상대 조명, 오차, 잔차, 조명이 닿지 않는 화소의 진폭.

    rPIE 는 조명이 닿지 않는 화소를 건드리지 않으므로, 그 화소의 값은 세기 제어가 시편에 곱한 배율의 누적이다.
    """
    alg, pc, order, n_iter = args
    P, P0, pos, shp, obj, ill, m, I = problem("pinhole")
    dark = tuple(np.argwhere(ill <= 1e-9 * ill.max())[0])
    rows = []

    def cb(k, o, p):
        a = np.abs(o)
        i = np.unravel_index(np.argmax(a), a.shape)
        rows.append((k, float(np.sum(np.abs(p) ** 2)), float(a.max()), float(ill[i] / ill.max()),
                     error(o, obj, m), residual(o, p, pos, I), float(a[dark])))

    with np.errstate(all="ignore"):
        rpie(I, pos, P0, shp, n_iter, probe_power="auto" if pc else None, rng=order, callback=cb, **DRIFT[alg])
    return alg, pc, order, np.array(rows), float(np.max(np.sum(I, axis=(-2, -1))))


def drift(n_iter=300):
    """핀홀. mPIE(Maiden) 세기 제어 없음/있음, 비교용 rPIE, 모멘텀을 시편에만/프로브에만 건 경우."""
    jobs = [("mPIE", pc, s, n_iter) for pc in (False, True) for s in ORDERS]
    jobs += [(a, True, s, n_iter) for a in ("rPIE", "mom_obj", "mom_prb") for s in ORDERS]
    with Pool(14) as pool:
        res = pool.map(_drift, jobs)
    return {(a, pc, s): dict(rows=r, power=pw) for a, pc, s, r, pw in res}


def param_map(n_iter=300, orders=(0, 1, 2), etas=(0.0, 0.5, 0.7, 0.8, 0.9), gammas=(0.1, 0.2, 0.35, 0.5, 1.0)):
    jobs = []
    for c in PROBES:
        for e in etas:
            for g in gammas:
                kw = dict(gamma_obj=g, momentum=maiden(eta_obj=e, eta_prb=e) if e else None, probe_power="auto")
                jobs += [(c, (e, g), kw, s, n_iter, (n_iter,)) for s in orders]
    out = {}
    for c, a, s, tr, pe, dt in _pool(jobs):
        out.setdefault((c,) + a, []).append(tr[n_iter])
    return dict(etas=etas, gammas=gammas, res=out)


def _stall(args):
    kind, order, n_iter = args
    P, P0, pos, shp, obj, ill, m, I = problem("diffuser")
    with np.errstate(all="ignore"):
        if kind == "known":
            o, p = rpie(I, pos, P, shp, n_iter, probe_start=10**9, rng=order)
        else:
            o, p = rpie(I, pos, P0, shp, n_iter, gamma_obj=warm_gamma, probe_power="auto",
                        center_probe=(kind == "center"), rng=order)
    return kind, order, dict(err=error(o, obj, m), resid=residual(o, p, pos, I), probe_err=error(p, P),
                             offset=centroid(p) - centroid(P), probe=p, rec=align(o, obj, m), obj=obj, mask=m, P=P)


def stall(n_iter=500, orders=(0, 1, 2)):
    with Pool(9) as pool:
        res = pool.map(_stall, [(k, s, n_iter) for k in ("plain", "center", "known") for s in orders])
    return {(k, s): v for k, s, v in res}


def toy(n=120, etas=(0.0, 0.6, 0.9)):
    """f = (x^2 + 25 y^2)/2, 같은 걸음 1.8/25 에서 모멘텀 eta 만 바꾼다 (heavy ball).

    경사하강(eta = 0)은 가파른 y 방향 때문에 걸음을 크게 못 잡고, 완만한 x 방향으로 느리게 간다.
    알맞은 모멘텀은 x 방향으로 속도를 붙이고, 지나친 모멘텀은 최솟값을 지나쳐 출렁인다.
    """
    A = np.diag([1.0, 25.0])
    lr = 1.8 / 25.0
    out = {}
    for eta in etas:
        z = np.array([-9.0, 1.2]); v = np.zeros(2); path = [z.copy()]
        for _ in range(n):
            v = eta * v - lr * A @ z
            z = z + v
            path.append(z.copy())
        out[eta] = np.array(path)
    return out


def run_all():
    """단계마다 캐시에 덧붙인다. 이미 있는 단계는 건너뛴다 (다시 돌리려면 캐시를 지운다)."""
    t = time.time()
    R = load() if os.path.exists(CACHE) else {}
    for name, f in (("compare", compare), ("factors", factors), ("warmup", warmup), ("drift", drift),
                    ("param_map", param_map), ("stall", stall), ("toy", toy)):
        if name in R:
            continue
        R[name] = f()
        with open(CACHE, "wb") as fh:
            pickle.dump(R, fh)
        print(f"{name}: {time.time() - t:.0f} s", flush=True)
    return R


def load():
    with open(CACHE, "rb") as fh:
        return pickle.load(fh)


if __name__ == "__main__":
    run_all()

"""6편의 표와 그림에 들어가는 수치 실험. 핀홀 프로브, 페르마 나선 144 위치, 푸아송 노이즈.

  patterns()   같은 회절 패턴을 세 선량으로 뽑은 것, 그리고 광자가 0~1 개인 화소의 비율
  varstab()    푸아송 n 과 sqrt(n) 의 분산 (몬테카를로)
  dose()       프로브를 알 때 선량 x 방법 x 노이즈 표본 3개의 최종 오차, 10^6 광자의 수렴 궤적
  fromtruth()  참 해에서 출발했을 때 각 방법이 머무는 오차 (추정량 자체의 질)
  lowdose()    프로브를 모를 때: rPIE 를 계속 돌리면 / 세기 제어 없이 / ML 로 이어 정련하면
  frac()       총광자 수를 고정하고 위치 수를 바꿀 때
  images()     10^5 광자에서 rPIE 와 ML(푸아송)의 복원 영상

python3 studies6.py 로 돌리면 .cache_post6.pkl 에 단계마다 덧붙인다 (generate_figures.py 가 읽는다).
"""
import os
import pickle
import time
from functools import partial
from multiprocessing import Pool

import numpy as np

from ptycho6_core import (N, disc, pinhole_probe, make_object, scan, object_shape, illumination, intensities, rpie,
                          dm, error, align, residual, noisy, ml, intensities as _I)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache_post6.pkl")
DOSES = (1e4, 1e5, 1e6, 1e7)
SEEDS = (0, 1, 2)
METHODS = ("rPIE", "rPIE_decay", "DM", "amp", "poisson", "int")
CK = (1, 2, 5, 10, 20, 30, 50, 100, 150, 200, 300)


def decay(k, k0=30, tau=10.0):
    """처음 k0 회는 걸음 1, 그 뒤로 1/(1 + (k - k0)/tau) 로 줄인다."""
    return 1.0 if k <= k0 else 1.0 / (1.0 + (k - k0) / tau)


def warm(k):
    """5편: 처음 10 회만 시편 걸음 0.25."""
    return 0.25 if k <= 10 else 1.0


def problem(photons, seed, pos=None):
    P = pinhole_probe()
    pos = scan("fermat") if pos is None else pos
    shp = object_shape(pos)
    obj = make_object(shp)
    ill = illumination(pos, P, shp)
    I = intensities(obj, P, pos)
    n, s = noisy(I, photons, np.random.default_rng(seed))
    return P, pos, shp, obj, ill, ill > 0.1 * ill.max(), n, s


def run_known(method, n, pos, Ps, shp, n_iter, obj0=None, callback=None, rng=0):
    """프로브를 알 때(광자 단위로 맞춘 참 프로브 Ps) 한 방법을 돌린다."""
    with np.errstate(all="ignore"):
        if method == "rPIE":
            return rpie(n, pos, Ps, shp, n_iter, probe_start=10**9, obj0=obj0, rng=rng, callback=callback)[0]
        if method == "rPIE_decay":
            return rpie(n, pos, Ps, shp, n_iter, gamma_obj=decay, probe_start=10**9, obj0=obj0, rng=rng,
                        callback=callback)[0]
        if method == "DM":
            return dm(n, pos, Ps, shp, n_iter, probe_start=10**9, callback=callback)[0]
        return ml(n, pos, Ps, shp, n_iter, model=method, obj0=obj0, callback=callback)[0]


# ── 단계별 실험 ───────────────────────────────────────────────
def patterns():
    P, pos, shp, obj, ill, m, n, s = problem(1e6, 0)
    I = intensities(obj, P, pos)
    j = int(np.argmin(np.sum((pos - pos.mean(0)) ** 2, 1)))           # 가운데 위치
    out = {"I": I[j] / I.sum((-2, -1)).max()}
    for ph in (1e4, 1e6, 1e8):
        nn, _ = noisy(I, ph, np.random.default_rng(1))
        out[ph] = dict(pattern=nn[j], frac01=float(np.mean(nn <= 1)), max=float(nn[j].max()))
    return out


def varstab(trials=200000):
    lam = np.geomspace(0.05, 1e4, 40)
    x = np.random.default_rng(0).poisson(lam, size=(trials, lam.size))
    return dict(lam=lam, var_n=np.var(x, 0), var_sqrt=np.var(np.sqrt(x), 0))


def _dose(args):
    ph, method, seed = args
    P, pos, shp, obj, ill, m, n, s = problem(ph, seed)
    tr = {}
    run_known(method, n, pos, P * np.sqrt(s), shp, 300, rng=seed,
              callback=lambda k, o, p: tr.__setitem__(k, error(o, obj, m)) if k in CK else None)
    return ph, method, seed, tr


def dose():
    jobs = [(ph, me, sd) for ph in DOSES for me in METHODS for sd in SEEDS]
    with Pool(14) as pool:
        res = pool.map(_dose, jobs)
    return {(ph, me, sd): tr for ph, me, sd, tr in res}


def _fromtruth(args):
    method, n_iter = args
    P, pos, shp, obj, ill, m, n, s = problem(1e6, 0)
    o = run_known(method, n, pos, P * np.sqrt(s), shp, n_iter, obj0=obj.copy())
    o2 = run_known(method, n, pos, P * np.sqrt(s), shp, n_iter)
    return method, error(o, obj, m), error(o2, obj, m)


def fromtruth():
    with Pool(5) as pool:
        res = pool.map(_fromtruth, [("rPIE", 300), ("rPIE_decay", 300), ("amp", 300), ("poisson", 300), ("int", 1000)])
    return {me: dict(from_truth=a, from_blank=b) for me, a, b in res}


def _lowdose(args):
    ph, method, seed = args
    P, pos, shp, obj, ill, m, n, s = problem(ph, seed)
    P0 = disc(N, 14)
    dark = tuple(np.argwhere(ill <= 1e-9 * ill.max())[0])
    rows = []

    def cb(k, o, p, off=0):
        rows.append((k + off, error(o, obj, m), float(np.sum(np.abs(p) ** 2)), float(np.abs(o[dark]))))

    with np.errstate(all="ignore"):
        pp = None if method == "rPIE_nopc" else "auto"
        o, p = rpie(n, pos, P0, shp, 50, gamma_obj=warm, probe_power=pp, rng=seed, callback=cb)
        if method in ("rPIE", "rPIE_nopc"):
            rpie(n, pos, p, shp, 400, obj0=o, probe_start=1, probe_power=pp, rng=seed + 100,
                 callback=partial(cb, off=50))
        else:
            ml(n, pos, p, shp, 400, model=method, obj0=o, update_probe=True, callback=partial(cb, off=50))
    o_true = float(np.sum(np.abs(P) ** 2) * s)
    return ph, method, seed, np.array(rows), o_true


def lowdose():
    jobs = [(ph, me, sd) for ph in (1e5, 1e6) for me in ("rPIE", "rPIE_nopc", "amp", "poisson") for sd in SEEDS]
    with Pool(14) as pool:
        res = pool.map(_lowdose, jobs)
    return {(ph, me, sd): dict(rows=r, true_power=tp) for ph, me, sd, r, tp in res}


def fermat_area(step, side=96):
    """같은 넓이(side x side)를 덮되 점 사이 간격이 step 인 페르마 나선."""
    m = int(round((side / step) ** 2))
    k = np.arange(m)
    c = step / np.sqrt(np.pi)
    th = k * np.pi * (3 - np.sqrt(5))
    pos = np.rint(np.stack([c * np.sqrt(k) * np.sin(th), c * np.sqrt(k) * np.cos(th)], 1)).astype(int)
    return pos - pos.min(0) + 4


def _frac(args):
    step, total, seed = args
    pos = fermat_area(step)
    P, pos, shp, obj, ill, m, n, s = problem(total / len(pos), seed, pos)
    o = run_known("amp", n, pos, P * np.sqrt(s), shp, 300)
    return step, seed, len(pos), error(o, obj, m)


def frac(total=1.44e8):
    with Pool(14) as pool:
        res = pool.map(_frac, [(st, total, sd) for st in (5, 6, 8, 11, 14) for sd in SEEDS])
    out = {}
    for st, sd, npos, e in res:
        out.setdefault(st, dict(npos=npos, err=[]))["err"].append(e)
    return out


def _images(method):
    P, pos, shp, obj, ill, m, n, s = problem(1e5, 0)
    o = run_known(method, n, pos, P * np.sqrt(s), shp, 300)
    return method, align(o, obj, m), error(o, obj, m)


def images():
    with Pool(3) as pool:
        res = pool.map(_images, ["rPIE", "poisson", "int"])
    P, pos, shp, obj, ill, m, n, s = problem(1e5, 0)
    return dict(obj=obj, mask=m, rec={me: (a, e) for me, a, e in res})


STAGES = (("patterns", patterns), ("varstab", varstab), ("dose", dose), ("fromtruth", fromtruth),
          ("lowdose", lowdose), ("frac", frac), ("images", images))


def run_all():
    """단계마다 캐시에 덧붙인다. 이미 있는 단계는 건너뛴다 (다시 돌리려면 캐시에서 그 키를 지운다)."""
    t = time.time()
    R = load() if os.path.exists(CACHE) else {}
    for name, f in STAGES:
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

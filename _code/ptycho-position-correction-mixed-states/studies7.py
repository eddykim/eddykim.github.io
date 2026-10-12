"""7편의 표와 그림에 들어가는 수치 실험. 핀홀 프로브, 페르마 나선 144 위치, 노이즈 없음.

  poserr()     위치 오차 sigma(가우시안, 화소) x 보정 없음·어닐링·기울기 x 오차 표본 3개
  posimg()     sigma = 1 화소: 참·명목·보정 위치와 복원 영상
  modes_info() 넓은 광원 두 크기의 결맞음 모드와 점유율
  modes()      모드 수 K x 광원 크기 x 순서 3개의 혼합 상태 rPIE, 그 일부의 ML 정련
  trunc()      참 상태를 모드 K 개로 잘랐을 때의 데이터 잔차
  modeimg()    모드 하나와 셋으로 복원한 시편, 복원한 모드

python3 studies7.py 로 돌리면 .cache_post7.pkl 에 단계마다 덧붙인다 (generate_figures.py 가 읽는다).
"""
import os
import pickle
import time
from multiprocessing import Pool

import numpy as np

from ptycho7_core import (N, disc, pinhole_probe, make_object, scan, object_shape, illumination, error, align,
                          intensities_sub, object_shape_f, rpie_pos, source_probes, coherent_modes, modal_intensities,
                          residual_modes, rpie_modes, ml_modes)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache_post7.pkl")
SIGMAS = (0.0, 0.25, 0.5, 1.0, 2.0)
SEEDS = (0, 1, 2)
CK = (1, 5, 10, 20, 30, 50, 100, 150, 200, 300)
SOURCES = (0.012, 0.02)
KS = (1, 2, 3, 5, 8)


# ── 위치 오차 ─────────────────────────────────────────────────
def pos_problem(sigma, seed):
    rng = np.random.default_rng(100 + seed)
    P = pinhole_probe()
    pos_t = scan("fermat").astype(float) + rng.uniform(0, 1, (144, 2))      # 참 위치는 소수 화소
    shp = object_shape_f(pos_t)
    obj = make_object(shp)
    ill = illumination(np.rint(pos_t).astype(int), P, shp)
    err = rng.normal(0, sigma, pos_t.shape) if sigma else np.zeros_like(pos_t)
    return P, pos_t, pos_t + err, shp, obj, ill > 0.1 * ill.max(), intensities_sub(obj, P, pos_t)


def pos_rms(pos, pos_t):
    d = pos - pos_t
    d = d - d.mean(0)
    return float(np.sqrt(np.mean(np.sum(d ** 2, 1))))


def _poserr(args):
    sigma, correct, seed = args
    P, pos_t, pos_n, shp, obj, m, I = pos_problem(sigma, seed)
    tr = {}
    t = time.time()
    with np.errstate(all="ignore"):
        rpie_pos(I, pos_n, disc(N, 14), shp, 300, correct=correct, rng=seed,
                 callback=lambda k, o, p, pos: tr.__setitem__(k, (error(o, obj, m), pos_rms(pos, pos_t))) if k in CK else None)
    return sigma, correct, seed, tr, pos_rms(pos_n, pos_t), time.time() - t


def poserr():
    jobs = [(0.0, None, sd) for sd in SEEDS]
    jobs += [(s, c, sd) for s in SIGMAS[1:] for c in (None, "anneal", "grad") for sd in SEEDS]
    with Pool(14) as pool:
        res = pool.map(_poserr, jobs)
    return {(s, c, sd): dict(trace=tr, rms0=r0, sec=dt) for s, c, sd, tr, r0, dt in res}


def _posimg(correct):
    P, pos_t, pos_n, shp, obj, m, I = pos_problem(1.0, 0)
    with np.errstate(all="ignore"):
        o, p, pos = rpie_pos(I, pos_n, disc(N, 14), shp, 300, correct=correct, rng=0)
    return correct, dict(pos=pos, rec=align(o, obj, m), err=error(o, obj, m))


def posimg():
    with Pool(3) as pool:
        res = dict(pool.map(_posimg, [None, "anneal", "grad"]))
    P, pos_t, pos_n, shp, obj, m, I = pos_problem(1.0, 0)
    return dict(pos_t=pos_t, pos_n=pos_n, obj=obj, mask=m, runs=res)


# ── 혼합 상태 ─────────────────────────────────────────────────
def mode_problem(sig_t):
    Ps, ws = source_probes(sig_t)
    modes_t, occ = coherent_modes(Ps, ws)
    pos = scan("fermat")
    shp = object_shape(pos)
    obj = make_object(shp)
    ill = illumination(pos, modes_t[0], shp)
    return modes_t, occ, pos, shp, obj, ill > 0.1 * ill.max(), modal_intensities(obj, modes_t, pos)


def modes_info():
    out = {}
    for st in (0.0,) + SOURCES:
        M, occ = coherent_modes(*source_probes(st))
        out[st] = dict(modes=M[:4], occ=occ[:8])
    return out


def _modes(args):
    sig_t, K, order, refine = args
    modes_t, occ, pos, shp, obj, m, I = mode_problem(sig_t)
    out = {}
    with np.errstate(all="ignore"):
        o, M = rpie_modes(I, pos, disc(N, 14), K, shp, 300, rng=order)
        p = np.sum(np.abs(M) ** 2, (1, 2))
        out["rpie"] = dict(err=error(o, obj, m), resid=residual_modes(o, M, pos, I), occ=p / p.sum())
        if refine:
            o, M = ml_modes(I, pos, o, M, 300)
            p = np.sum(np.abs(M) ** 2, (1, 2))
            out["ml"] = dict(err=error(o, obj, m), resid=residual_modes(o, M, pos, I), occ=p / p.sum())
    return sig_t, K, order, out


def modes():
    refine = {(0.012, 3), (0.012, 5), (0.02, 8)}
    jobs = [(st, K, od, (st, K) in refine) for st in SOURCES for K in KS for od in SEEDS]
    jobs += [(0.0, 1, od, False) for od in SEEDS]
    with Pool(14) as pool:
        res = pool.map(_modes, jobs)
    return {(st, K, od): out for st, K, od, out in res}


def trunc():
    out = {}
    for st in SOURCES:
        modes_t, occ, pos, shp, obj, m, I = mode_problem(st)
        out[st] = {K: residual_modes(obj, modes_t[:K], pos, I) for K in (1, 2, 3, 5, 8, 12)}
    return out


def _modeimg(K):
    modes_t, occ, pos, shp, obj, m, I = mode_problem(0.012)
    with np.errstate(all="ignore"):
        o, M = rpie_modes(I, pos, disc(N, 14), K, shp, 300, rng=0)
    return K, dict(rec=align(o, obj, m), err=error(o, obj, m), modes=M)


def modeimg():
    with Pool(2) as pool:
        res = dict(pool.map(_modeimg, [1, 3]))
    modes_t, occ, pos, shp, obj, m, I = mode_problem(0.012)
    return dict(obj=obj, mask=m, modes_t=modes_t[:3], occ=occ[:3], runs=res)


STAGES = (("poserr", poserr), ("posimg", posimg), ("modes_info", modes_info), ("modes", modes), ("trunc", trunc),
          ("modeimg", modeimg))


def run_all():
    """단계마다 캐시에 덧붙인다. 이미 있는 단계는 건너뛴다."""
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

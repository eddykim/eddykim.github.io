"""4편 2절: 흔든 격자에서 가끔 생기는 정체(stagnation)는 무엇인가.

정체한 복원은 잔차가 크고(데이터를 설명하지 못한다), 프로브의 무게중심이 참 프로브에서 1~2 화소
비켜나 있다. 수렴한 복원은 비켜남이 0 이다. 평행이동 모호성 방향으로 프로브만 미끄러지고, 시편이
경계(조명 밖은 1 로 남아 있다) 때문에 함께 따라가지 못한 상태로 읽힌다.
"""
from multiprocessing import Pool

import numpy as np

from jitter_sweep import jittered, P, P0
from ptycho4_core import object_shape, make_object, intensities, rpie, residual, scan


def centroid(p):
    w = np.abs(p) ** 2
    y, x = np.indices(p.shape)
    return np.array([(w * y).sum(), (w * x).sum()]) / w.sum()


def run(case):
    J, seed, order = case
    pos = scan("fermat") if J == "fermat" else jittered(J, seed)
    shp = object_shape(pos)
    I = intensities(make_object(shp), P, pos)
    o, p = rpie(I, pos, P0, shp, 300, alpha=0.05, rng=order)
    return case, residual(o, p, pos, I), centroid(p) - centroid(P)


if __name__ == "__main__":
    cases = [(1, 0, 0), (1, 1, 0), (2, 0, 0), (1, 0, 1), (1, 0, 2), ("fermat", 0, 0)]
    with Pool(len(cases)) as pool:
        for case, res, d in pool.map(run, cases):
            print(case, "residual %.1e" % res, "probe centroid offset (y, x) %.2f %.2f" % tuple(d))

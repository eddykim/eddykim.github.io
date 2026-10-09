"""4편 2절: 규칙 격자를 얼마나 흔들어야 격자 무늬가 사라지는가.

같은 12 x 12 격자(간격 8 화소)의 위치마다 [-J, J] 화소의 정수 흔들림을 더하고, 프로브를 모르는
rPIE 로 300회 복원한다. 흔들림 난수 세 가지 x 패턴 순서 세 가지를 모두 돌린다.

오차만 보면 두 가지 실패가 섞인다. 모호성은 데이터를 완전히 설명하면서(잔차가 작다) 틀리고,
정체(stagnation)는 데이터를 설명하지 못한 채(잔차가 크다) 멈춘다. 잔차가 1e-2 를 넘으면 정체로 센다.
"""
from multiprocessing import Pool

import numpy as np

from ptycho4_core import N, disc, pinhole_probe, make_object, object_shape, illumination, intensities, rpie, \
    error, residual, scan

P, P0 = pinhole_probe(), disc(N, 14)
STALL = 1e-2


def jittered(J, seed, step=8, n_side=12, margin=4):
    iy, ix = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    pos = np.stack([iy.ravel(), ix.ravel()], 1) * step
    if J:
        pos = pos + np.random.default_rng(seed).integers(-J, J + 1, pos.shape)
    return pos - pos.min(0) + margin


def run(args, n_iter=300):
    J, seed, order = args
    pos = scan("fermat") if J == "fermat" else jittered(J, seed)
    shp = object_shape(pos)
    obj = make_object(shp)
    ill = illumination(pos, P, shp)
    mask = ill > 0.1 * ill.max()
    I = intensities(obj, P, pos)
    o, p = rpie(I, pos, P0, shp, n_iter, alpha=0.05, rng=order)
    return J, seed, order, error(o, obj, mask), residual(o, p, pos, I)


def sweep(Js=(0, 1, 2, 4), seeds=(0, 1, 2), orders=(0, 1, 2), procs=8):
    jobs = [(J, s, r) for J in Js for s in (seeds if J else (0,)) for r in orders]
    jobs += [("fermat", 0, r) for r in range(12)]   # 나선은 흔들림 난수가 없으니 순서를 더 많이 본다
    with Pool(procs) as pool:
        res = pool.map(run, jobs)
    out = {}
    for J, s, r, e, res_ in res:
        out.setdefault(J, []).append((s, r, e, res_))
    return out


if __name__ == "__main__":
    for J, rows in sweep().items():
        ok = [e for _, _, e, r in rows if r < STALL]
        print(J, "runs", len(rows), "stalled", len(rows) - len(ok),
              "converged median err", round(float(np.median(ok)), 5) if ok else None,
              "max", round(max(ok), 5) if ok else None, flush=True)
        for s, r, e, res_ in rows:
            print("   seed", s, "order", r, "err %.4f resid %.1e" % (e, res_), flush=True)

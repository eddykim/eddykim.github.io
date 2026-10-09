"""4편의 표와 그림에 들어가는 수치 실험.

  scans()       규칙 격자 / 흔든 격자 / 페르마 나선에서 프로브를 모를 때(rPIE)와 알 때(프로브 고정)의 오차,
                그리고 오차 전력 가운데 주사 격자 주파수에 몰린 몫
  algorithms()  핀홀, 초점이 어긋난 집속 프로브, 확산판 프로브에서 ePIE 와 rPIE 의 수렴

모멘텀(mPIE), 프로브 세기 제어, 차분 사상(DM)은 5편에서 다룬다 (ptycho4_core 에 함께 들어 있다).
"""
import numpy as np

from ptycho4_core import (N, disc, pinhole_probe, diffuser_probe, focused_probe, make_object, scan,
                          object_shape, illumination, intensities, rpie, error, align, residual, grid_power_fraction)

REC = (1, 2, 5, 10, 20, 50, 100, 150, 200)
ALGS = {"ePIE": dict(alpha=1.0, beta=1.0), "rPIE": dict(alpha=0.05, beta=1.0)}


def problem(probe, kind="fermat"):
    pos = scan(kind)
    shp = object_shape(pos)
    obj = make_object(shp)
    ill = illumination(pos, probe, shp)
    return pos, obj, ill > 0.1 * ill.max(), intensities(obj, probe, pos)


def scans(n_iter=300, orders=(0, 1, 2)):
    """패턴 순서를 바꿔 여러 번 복원하고, 오차가 중앙값인 복원을 그림에 쓴다.

    순서 하나만 보면 정체(stagnation)를 주사의 탓으로 오해한다. 정체는 데이터를 설명하지 못하므로
    잔차(residual)가 크고, 모호성은 데이터를 완전히 설명하면서도 틀린다.
    """
    P, P0 = pinhole_probe(), disc(N, 14)
    out = {}
    for kind in ("regular", "jitter", "fermat"):
        pos, obj, mask, I = problem(P, kind)
        ys, xs = np.where(mask)
        cy, cx, h = int(ys.mean()), int(xs.mean()), 40
        runs = []
        for s in orders:
            o, p = rpie(I, pos, P0, obj.shape, n_iter, rng=s, **ALGS["rPIE"])
            a = align(o, obj, mask)
            frac, share = grid_power_fraction((a - obj)[cy - h:cy + h, cx - h:cx + h], 8)
            runs.append(dict(err=error(o, obj, mask), probe=error(p, P), resid=residual(o, p, pos, I),
                             grid_frac=frac, grid_share=share, rec=a,
                             err_crop=(a - obj)[cy - h:cy + h, cx - h:cx + h]))
        o_known, _ = rpie(I, pos, P, obj.shape, n_iter, probe_start=10**9, **ALGS["rPIE"])   # 프로브 고정 = 앎
        med = sorted(runs, key=lambda r: r["err"])[len(runs) // 2]
        out[kind] = dict(med, known=error(o_known, obj, mask), obj=obj, mask=mask, pos=pos,
                         runs=[{k: v for k, v in r.items() if isinstance(v, float)} for r in runs])
    return out


def algorithms(seeds=(0, 1, 2), n_iter=200):
    cases = {"pinhole": (pinhole_probe(), disc(N, 14)),
             "defocus": (focused_probe(), focused_probe(defocus=0.0)),
             "diffuser": (diffuser_probe(), disc(N, 12))}
    out = {}
    for case, (P, P0) in cases.items():
        pos, obj, mask, I = problem(P)
        out[(case, "_truth")] = dict(obj=obj, mask=mask, probe=P, probe0=P0)
        for name, kw in ALGS.items():
            for s in seeds:
                tr = {}
                o, p = rpie(I, pos, P0, obj.shape, n_iter, rng=s, **kw,
                            callback=lambda k, o, p: tr.__setitem__(k, error(o, obj, mask)) if k in REC else None)
                out[(case, name, s)] = dict(trace=tr, probe_err=error(p, P), probe=p)
    return out


if __name__ == "__main__":
    import time
    t = time.time()
    for k, v in scans().items():
        print(k, {a: round(b, 4) for a, b in v.items() if isinstance(b, float)})
    for key, v in algorithms().items():
        if key[1] != "_truth":
            print(key, {a: round(b, 4) for a, b in v["trace"].items() if a in (10, 50, 100, 200)},
                  "probe", round(v["probe_err"], 4))
    print(f"{time.time() - t:.0f} s")

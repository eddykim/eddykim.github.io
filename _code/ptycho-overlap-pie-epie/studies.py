"""3편의 표와 그림에 들어가는 수치 실험.

공통 설정: 프로브는 반지름 12 um 핀홀을 200 um 전파한 장, 시편은 진폭 0.2-1, 위상 0-pi 인 복소 영상.
주사는 84 화소를 덮는 흔들린 격자. 오차는 조명 합이 최대의 10 % 를 넘는 곳에서
평행이동, 위상 기울기, 복소 척도를 맞춘 뒤 잰다 (ptycho_core.error).
"""
import numpy as np

from ptycho_core import (pinhole_probe, circular_aperture, make_object, raster, object_shape, illumination,
                         intensities, poisson, pie, epie, error, probe_d90, linear_overlap, cdi_hio, cdi_error)

P = pinhole_probe()
D90 = probe_d90(P)
WRONG_PROBE = circular_aperture(P.shape[0], 14)     # 핀홀 반지름 12 화소보다 2 화소 크게 (Maiden 2009 처럼)


def setup(step, seed=0):
    pos = raster(step, seed=seed)
    shape = object_shape(pos)
    obj = make_object(shape)
    mask = illumination(pos, P, shape) > 0.1 * illumination(pos, P, shape).max()
    return pos, obj, mask


def convergence(step=8, n_pie=100, n_epie=200, record=(1, 2, 5, 10, 20, 50, 100, 200)):
    """정확한 프로브의 PIE, 틀린 프로브의 PIE, 틀린 프로브에서 출발한 ePIE 의 오차 궤적과 결과."""
    pos, obj, mask = setup(step)
    I = intensities(obj, P, pos)
    out = {"pos": pos, "obj": obj, "mask": mask}
    for key, probe in (("pie", P), ("pie_wrong", WRONG_PROBE)):
        tr = {}
        o = pie(I, pos, probe, obj.shape, n_pie,
                callback=lambda k, o: tr.__setitem__(k, error(o, obj, mask)) if k in record else None)
        out[key] = (tr, o)
    tr, trp = {}, {}

    def cb(k, o, p):
        if k in record:
            tr[k] = error(o, obj, mask)
            trp[k] = error(p, P)
    o, p = epie(I, pos, WRONG_PROBE, obj.shape, n_epie, callback=cb)
    out["epie"] = (tr, o, trp, p)
    return out


def overlap(steps=(4, 6, 8, 10, 12, 16, 20, 24), n_pie=100, n_epie=200):
    out = {}
    for s in steps:
        pos, obj, mask = setup(s)
        I = intensities(obj, P, pos)
        o = pie(I, pos, P, obj.shape, n_pie)
        o2, p2 = epie(I, pos, WRONG_PROBE, obj.shape, n_epie)
        out[s] = dict(overlap=linear_overlap(s, D90), overlap_pinhole=linear_overlap(s, 24), n_pos=len(pos),
                      pie=error(o, obj, mask), epie=error(o2, obj, mask), probe=error(p2, P),
                      rec_pie=o, rec_epie=o2, obj=obj, mask=mask)
    return out


def noise(doses=(125, 1255, 12549), step=8, n_pie=100, n_epie=200):
    """조명된 화소당 광자 수가 같을 때 PIE, ePIE, 단일 패턴 CDI 의 오차.

    CDI 는 같은 시편의 가운데 48 x 48 을 128 x 128 상자에 넣어 (sigma = 7.1) 2편의 HIO 로 푼다.
    """
    pos, obj, mask = setup(step)
    I = intensities(obj, P, pos)
    c = obj.shape[0] // 2
    box = np.zeros((128, 128), complex)
    box[40:88, 40:88] = obj[c - 24:c + 24, c - 24:c + 24]
    sup = np.zeros((128, 128), bool)
    sup[40:88, 40:88] = True
    out = {}
    for dose in doses:
        per_pattern = dose * mask.sum() / len(pos)
        In = poisson(I, per_pattern, rng=1)
        o = pie(In, pos, P, obj.shape, n_pie)
        o2, _ = epie(In, pos, WRONG_PROBE, obj.shape, n_epie)
        recs = [cdi_hio(box, sup, dose * 48 * 48, seed=r) for r in range(6)]
        errs = [cdi_error(r, box) for r in recs]
        out[dose] = dict(pie=error(o, obj, mask), epie=error(o2, obj, mask), cdi=float(np.median(errs)),
                         per_pattern=per_pattern, rec_pie=o, rec_cdi=recs[int(np.argmin(errs))], box=box,
                         obj=obj, mask=mask)
    return out


if __name__ == "__main__":
    import time
    t = time.time()
    print("D90 =", D90)
    cv = convergence()
    for k in ("pie", "pie_wrong"):
        print(k, {a: round(b, 4) for a, b in cv[k][0].items()})
    print("epie", {a: round(b, 4) for a, b in cv["epie"][0].items()}, "probe", {a: round(b, 4) for a, b in cv["epie"][2].items()})
    for s, r in overlap().items():
        print(f"step {s:2d} ov(D90)={r['overlap']:.2f} ov(pinhole)={r['overlap_pinhole']:.2f} J={r['n_pos']:3d} "
              f"PIE={r['pie']:.4f} ePIE={r['epie']:.4f} probe={r['probe']:.4f}")
    for d, r in noise().items():
        print(f"dose {d}/px ({r['per_pattern']:.0f}/pattern): CDI={r['cdi']:.3f} PIE={r['pie']:.3f} ePIE={r['epie']:.3f}")
    print(f"{time.time() - t:.0f} s")

"""8편의 표와 그림에 들어가는 수치 실험. 고해상 256 x 256 시편, 카메라 64 x 64, 동공 반지름 12, 정사각 LED 격자.

  setup()    수직 조명·비스듬한 명시야·암시야 원본 영상과 스펙트럼
  synth()    LED 15 x 15 (간격 10) 복원: 저해상 영상 한 장, 복원, 덮인 대역으로 거른 참값
  overlap()  LED 간격(=동공 원의 겹침) x 순서(가운데부터 / 무작위 2개)
  pupil()    렌즈 수차 두 크기 x LED 간격 3개: 동공을 모른 채 / 함께 찾을 때
  noise()    선량 x 비용(rPIE, 진폭·푸아송·세기 ML) x 노이즈 표본 3개

python3 studies8.py 로 돌리면 .cache_post8.pkl 에 단계마다 덧붙인다 (generate_figures.py 가 읽는다).
"""
import os
import pickle
import time
from multiprocessing import Pool

import numpy as np

from ptycho8_core import (NH, NL, fft2c, ifft2c, make_object, error, fp_pupil, leds, band, fp_spectrum, fp_images,
                          coverage, lowpass, overlap_frac, fp_rpie, fp_ml)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache_post8.pkl")
R = 12
KMAX = 84                                 # 가장 바깥 LED 의 조명 파수 (푸리에 화소)
SEEDS = (0, 1, 2)
ABS = ((1.0, 0.7), (2.0, 1.5))


def obj_mask():
    obj = make_object((NH, NH), amp_range=(0.4, 1.0), phase_max=2.0)
    m = np.zeros((NH, NH), bool)
    m[32:-32, 32:-32] = True
    return obj, m


def grid(step):
    return leds(2 * int(KMAX // step) + 1, step)


def setup():
    obj, m = obj_mask()
    ks = grid(10)
    picks = {"center": (0, 0), "bright": (0, 10), "dark": (0, 30)}
    return dict(images={k: fp_images(obj, fp_pupil(R), [v])[0] for k, v in picks.items()},
                spectrum=np.abs(fp_spectrum(obj)), ks=ks, picks=picks, R=R)


def synth():
    obj, m = obj_mask()
    ks = grid(10)
    I = fp_images(obj, fp_pupil(R), ks)
    cov = coverage(R, ks)
    ref = lowpass(obj, cov)
    with np.errstate(all="ignore"):
        o, _ = fp_rpie(I, ks, fp_pupil(R), 50)
    low = lowpass(obj, coverage(R, [(0, 0)]))                  # 수직 조명 한 장이 담는 정보 (결맞은 저해상 영상)
    return dict(obj=obj, mask=m, rec=o, ref=ref, low=low, I0=I[np.argmin(np.sum(ks ** 2, 1))], err=error(o, ref, m),
                err_true=error(o, obj, m), cov_r=float(np.sqrt(cov.sum() / np.pi)), n=len(ks))


def _overlap(args):
    step, order = args
    obj, m = obj_mask()
    ks = grid(step)
    ref = lowpass(obj, coverage(R, ks))
    I = fp_images(obj, fp_pupil(R), ks)
    with np.errstate(all="ignore"):
        o, _ = fp_rpie(I, ks, fp_pupil(R), 50, order="center" if order is None else "random", rng=order or 0)
    return step, order, error(o, ref, m), len(ks)


def overlap():
    jobs = [(st, od) for st in (4, 8, 10, 12, 13, 14, 16, 20, 24) for od in (None, 1, 2)]
    with Pool(14) as pool:
        res = pool.map(_overlap, jobs)
    return {(st, od): dict(err=e, n=n, ovl=overlap_frac(R, st)) for st, od, e, n in res}


def _pupil_err(P2, Pt):
    """동공 위상 오차 (rms, rad). 상수 위상과 위상 기울기(시편 이동과 맞바뀌는 모호성)를 뺀다."""
    sup = np.abs(Pt) > 0
    y, x = np.indices(Pt.shape)
    z = np.where(sup, P2 * np.conj(Pt), 0)
    gx = np.angle(np.sum(z[:, 1:] * np.conj(z[:, :-1])))
    gy = np.angle(np.sum(z[1:] * np.conj(z[:-1])))
    z = z * np.exp(-1j * (gx * x + gy * y))
    d = np.angle(z[sup])
    d -= np.angle(np.mean(np.exp(1j * d)))
    return float(np.sqrt(np.mean(np.angle(np.exp(1j * d)) ** 2)))


def _pupil(args):
    step, ab = args
    obj, m = obj_mask()
    ks = grid(step)
    ref = lowpass(obj, coverage(R, ks))
    Pt = fp_pupil(R, ab)
    I = fp_images(obj, Pt, ks)
    with np.errstate(all="ignore"):
        o1, _ = fp_rpie(I, ks, fp_pupil(R), 150)
        o2, P2 = fp_rpie(I, ks, fp_pupil(R), 150, update_pupil=True)
    sup = np.abs(Pt) > 0
    return step, ab, dict(no_pupil=error(o1, ref, m), epry=error(o2, ref, m), pupil_rms=_pupil_err(P2, Pt),
                          ab_rms=float(np.std(np.angle(Pt[sup]))), ovl=overlap_frac(R, step), n=len(ks),
                          rec_no=o1, rec=o2, P=P2, Pt=Pt)


def pupil():
    with Pool(6) as pool:
        res = pool.map(_pupil, [(st, ab) for st in (10, 6, 4) for ab in ABS])
    return {(st, ab): v for st, ab, v in res}


def _noise(args):
    ph, model, seed = args
    obj, m = obj_mask()
    ks = grid(10)
    P = fp_pupil(R)
    ref = lowpass(obj, coverage(R, ks))
    I = fp_images(obj, P, ks)
    s = ph / I.sum((-2, -1)).max()
    n = np.random.default_rng(seed).poisson(I * s).astype(float)
    tot = n.sum((-2, -1))
    bf = np.sum(ks ** 2, 1) <= R ** 2
    with np.errstate(all="ignore"):
        if model == "rpie":
            o, _ = fp_rpie(n, ks, P, 100)
        else:
            o0, _ = fp_rpie(n, ks, P, 5)
            o = fp_ml(n, ks, P, 200, model=model, S0=fp_spectrum(o0))
    return ph, model, seed, error(o, ref, m), float(tot[bf].mean()), float(tot[~bf].mean())


def noise():
    jobs = [(ph, mo, sd) for ph in (1e5, 1e6, 1e7) for mo in ("rpie", "amp", "poisson", "int") for sd in SEEDS]
    with Pool(14) as pool:
        res = pool.map(_noise, jobs)
    return {(ph, mo, sd): dict(err=e, bf=b, df=d) for ph, mo, sd, e, b, d in res}


STAGES = (("setup", setup), ("synth", synth), ("overlap", overlap), ("pupil", pupil), ("noise", noise))


def run_all():
    """단계마다 캐시에 덧붙인다. 이미 있는 단계는 건너뛴다."""
    t = time.time()
    R_ = load() if os.path.exists(CACHE) else {}
    for name, f in STAGES:
        if name in R_:
            continue
        R_[name] = f()
        with open(CACHE, "wb") as fh:
            pickle.dump(R_, fh)
        print(f"{name}: {time.time() - t:.0f} s", flush=True)
    return R_


def load():
    with open(CACHE, "rb") as fh:
        return pickle.load(fh)


if __name__ == "__main__":
    run_all()

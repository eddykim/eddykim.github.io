"""9편의 표와 그림에 들어가는 수치 실험. 자원 규칙에 따라 프로세스 4 개, 프로세스마다 스레드 1 개로 돌린다.

  grads()      자동미분 기울기와 손 기울기(6편)의 화소별 비교
  optim()      프로브를 모를 때(10^6 광자, 노이즈 표본 3개): rPIE 50회로 시작한 뒤 손 구현 ML, Adam(걸음 3가지), L-BFGS
  blank()      빈 시편에서 출발할 때의 Adam 과 L-BFGS
  positions()  위치 오차 sigma x 오차 표본 3개: rPIE 50회로 시작한 뒤 시편·프로브·위치를 Adam 으로 함께
  multislice() 두 장 사이 간격 dz x 장 수(1, 2): 빈 시편에서 Adam 500회 (프로브 앎, 노이즈 없음)

python3 studies9.py 로 돌리면 .cache_post9.pkl 에 단계마다 덧붙인다 (generate_figures.py 가 읽는다).
"""
import os
import pickle
import time
from multiprocessing import Pool

import numpy as np
import torch
from skimage import data, transform

import ad_torch as A
from ptycho9_core import (N, disc, pinhole_probe, focused_probe, make_object, scan, object_shape, illumination,
                          intensities, error, noisy, ml, rpie, grad_obj, intensities_sub, object_shape_f, split, rpie_pos)

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache_post9.pkl")
SEEDS = (0, 1, 2)
CK = (1, 2, 5, 10, 20, 50, 100, 200, 300)
PROCS = 4


def _init():
    torch.set_num_threads(1)


def warm(k):
    return 0.25 if k <= 10 else 1.0


# ── 기울기 ────────────────────────────────────────────────────
def grads():
    rng = np.random.default_rng(0)
    P = pinhole_probe(); pos = scan("fermat")[:40]; shp = object_shape(pos); obj = make_object(shp)
    n, s = noisy(intensities(obj, P, pos), 1e5, rng)
    Ps = P * np.sqrt(s)
    O1 = obj * (1 + 0.1 * rng.normal(size=shp))
    pick = rng.choice(np.flatnonzero(illumination(pos, P, shp) > 0), 2000, replace=False)
    out = {}
    for kind in ("amp", "poisson", "int"):
        Ot = A.t(O1).requires_grad_(True)
        A.loss(A.intensities(Ot, A.t(Ps), pos), torch.as_tensor(n), kind).backward()
        g_an, _ = grad_obj(O1, Ps, pos, n, kind)
        out[kind] = dict(ad=Ot.grad.numpy().ravel()[pick], hand=g_an.ravel()[pick],
                         rel=float(np.abs(Ot.grad.numpy() - 2 * g_an).max() / np.abs(2 * g_an).max()))
    return out


# ── 옵티마이저 ────────────────────────────────────────────────
def _problem_optim(seed):
    P = pinhole_probe(); P0 = disc(N, 14); pos = scan("fermat"); shp = object_shape(pos); obj = make_object(shp)
    ill = illumination(pos, P, shp); m = ill > 0.1 * ill.max()
    n, s = noisy(intensities(obj, P, pos), 1e6, np.random.default_rng(seed))
    return P0, pos, shp, obj, m, n


def _optim(seed):
    torch.set_num_threads(1)
    P0, pos, shp, obj, m, n = _problem_optim(seed)
    with np.errstate(all="ignore"):
        o0, p0 = rpie(n, pos, P0, shp, 50, gamma_obj=warm, probe_power="auto", rng=seed)
    out = {"start": error(o0, obj, m)}
    tr = {}; t0 = time.time()

    def cb(k, o, p):
        if k in CK:
            tr[k] = (error(o, obj, m), time.time() - t0)

    ml(n, pos, p0, shp, max(CK), model="poisson", obj0=o0, update_probe=True, callback=cb)
    out["hand"] = tr
    nt = torch.as_tensor(n)
    for name, lr in (("adam_0.003", 0.003), ("adam_0.01", 0.01), ("adam_0.03", 0.03), ("lbfgs", None)):
        O = A.t(o0).requires_grad_(True); Pp = A.t(p0).requires_grad_(True)
        opt = (torch.optim.Adam([O, Pp], lr=lr) if lr else
               torch.optim.LBFGS([O, Pp], lr=1, max_iter=1, history_size=20, line_search_fn="strong_wolfe"))

        def closure():
            opt.zero_grad()
            L = A.loss(A.intensities(O, Pp, pos), nt, "poisson")
            L.backward()
            return L

        tr = {}; t0 = time.time()
        for k in range(1, max(CK) + 1):
            opt.step(closure)
            if k in CK:
                tr[k] = (error(O.detach().numpy(), obj, m), time.time() - t0)
        out[name] = tr
    return seed, out


def optim():
    with Pool(3, initializer=_init) as pool:
        return dict(pool.map(_optim, SEEDS))


def _blank(name):
    torch.set_num_threads(1)
    P0, pos, shp, obj, m, n = _problem_optim(0)
    n = n * np.sum(np.abs(P0) ** 2) / n.sum((-2, -1)).max()       # 시작 프로브 세기에 맞춘다 (파세발)
    nt = torch.as_tensor(n)
    O = torch.ones(shp, dtype=A.DT, requires_grad=True); Pp = A.t(P0).requires_grad_(True)
    opt = (torch.optim.Adam([O, Pp], lr=0.01) if name == "adam" else
           torch.optim.LBFGS([O, Pp], lr=1, max_iter=1, history_size=20, line_search_fn="strong_wolfe"))

    def closure():
        opt.zero_grad()
        L = A.loss(A.intensities(O, Pp, pos), nt, "poisson")
        L.backward()
        return L

    tr = {}
    for k in range(1, max(CK) + 1):
        opt.step(closure)
        if k in CK:
            tr[k] = error(O.detach().numpy(), obj, m)
    return name, tr


def blank():
    with Pool(2, initializer=_init) as pool:
        return dict(pool.map(_blank, ["adam", "lbfgs"]))


# ── 위치 ──────────────────────────────────────────────────────
def _positions(args):
    sig, seed, use_pos = args
    torch.set_num_threads(1)
    rng = np.random.default_rng(100 + seed)
    P = pinhole_probe(); pos_t = scan("fermat").astype(float) + rng.uniform(0, 1, (144, 2))
    shp = object_shape_f(pos_t); obj = make_object(shp)
    ill = illumination(np.rint(pos_t).astype(int), P, shp); m = ill > 0.1 * ill.max()
    I = intensities_sub(obj, P, pos_t)
    pos_n = pos_t + (rng.normal(0, sig, pos_t.shape) if sig else 0)
    with np.errstate(all="ignore"):
        o0, p0, _ = rpie_pos(I, pos_n, disc(N, 14), shp, 50, correct=None, rng=seed)
    corner, frac = split(pos_n)
    nt = torch.as_tensor(I)
    O = A.t(o0).requires_grad_(True); Pp = A.t(p0).requires_grad_(True)
    d = torch.zeros((144, 2), dtype=torch.float64, requires_grad=True)
    fr = torch.as_tensor(frac)
    lr = 0.003
    opt = torch.optim.Adam([{"params": [O], "lr": lr}, {"params": [Pp], "lr": lr * float(np.abs(p0).max())},
                            {"params": [d], "lr": 0.02}])
    for k in range(1, 301):
        opt.zero_grad()
        L = A.loss(A.intensities(O, Pp, corner, fr + d), nt, "amp")
        L.backward()
        if not use_pos:
            d.grad.zero_()
        opt.step()
        with torch.no_grad():
            d -= d.mean(0)                    # 공동 평행이동 모호성: 보정량의 평균은 0
    rms = lambda x: float(np.sqrt(np.mean(np.sum((x - x.mean(0)) ** 2, 1))))
    return sig, seed, use_pos, dict(start=error(o0, obj, m), err=error(O.detach().numpy(), obj, m),
                                    rms0=rms(pos_n - pos_t), rms=rms(corner + frac + d.detach().numpy() - pos_t),
                                    resid=float(np.sqrt(L.item() / nt.sum().item())))


def positions():
    jobs = [(0.0, sd, False) for sd in SEEDS] + [(1.0, sd, False) for sd in SEEDS]
    jobs += [(s, sd, True) for s in (0.5, 1.0, 2.0) for sd in SEEDS]
    with Pool(PROCS, initializer=_init) as pool:
        res = pool.map(_positions, jobs)
    return {(s, sd, u): v for s, sd, u, v in res}


# ── 다층 ──────────────────────────────────────────────────────
DX = 0.5e-6


def _phase(img, shp, amp=0.8):
    a = transform.resize(img.astype(float), shp, anti_aliasing=True)
    a = (a - a.min()) / (a.max() - a.min())
    return np.exp(1j * amp * a)


def ms_problem():
    P = focused_probe(na=0.15, defocus=15e-6, dx=DX)
    pos = scan("fermat", step=4); shp = object_shape(pos)
    gray = data.astronaut().astype(float) @ [0.2125, 0.7154, 0.0721]
    O1, O2 = _phase(data.camera(), shp), _phase(gray, shp)
    ill = illumination(pos, P, shp)
    return P, pos, shp, O1, O2, ill > 0.1 * ill.max()


def _multislice(args):
    dz, ns, its = args
    torch.set_num_threads(1)
    P, pos, shp, O1, O2, m = ms_problem()
    H = A.propagator(dz, DX); Pt = A.t(P)
    with torch.no_grad():
        n = A.multislice([A.t(O1), A.t(O2)], Pt, pos, H)
    R = [torch.zeros(shp + (2,), dtype=torch.float64) for _ in range(ns)]
    for r in R:
        r[..., 0] = 1.0
        r.requires_grad_(True)
    opt = torch.optim.Adam(R, lr=0.02)
    hist = {}
    for k in range(1, its + 1):
        opt.zero_grad()
        S = [torch.view_as_complex(r) for r in R]
        Im = A.multislice(S, Pt, pos, H) if ns == 2 else A.intensities(S[0], Pt, pos)
        L = A.loss(Im, n, "amp")
        L.backward()
        opt.step()
        if k in (10, 50, 100, 200, 300, 500):
            hist[k] = float(np.sqrt(L.item() / n.sum().item()))
    S = [torch.view_as_complex(r).detach().numpy() for r in R]
    out = dict(resid=hist, slices=S)
    if ns == 2:
        out.update(e1=error(S[0], O1, m), e2=error(S[1], O2, m))
    else:
        out.update(e_prod=error(S[0], O1 * O2, m))
    return dz, ns, out


def multislice():
    jobs = [(dz, ns, 500) for dz in (1e-6, 3e-6, 10e-6, 30e-6, 80e-6) for ns in (1, 2)]
    with Pool(PROCS, initializer=_init) as pool:
        res = pool.map(_multislice, jobs)
    P, pos, shp, O1, O2, m = ms_problem()
    return dict(runs={(dz, ns): v for dz, ns, v in res}, O1=O1, O2=O2, mask=m)


STAGES = (("grads", grads), ("optim", optim), ("blank", blank), ("positions", positions), ("multislice", multislice))


def run_all():
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

"""10편 계산. 결과는 .cache_post10.pkl 에 저장하고 generate_figures.py 와 verify_post10.py 가 읽는다.

  python3 studies10.py          # 전부 (노트북 CPU, 4 프로세스로 약 15분)

  geometry      검출기 화소 -> 표면 공간주파수 사상, 반사율·반사 위상 곡선, 침투 깊이, 높이 감도
  tilt          22도, 기울어진 평면 보정 유무에 따른 복원 오차 (프로브를 안다), 검출기 수치구경 0.05–0.20
                + 보정의 세기 보간이 남기는 영역 비의 편향: 검출기 과표본 1, 2, 4 배, 22·30·40도
  pipeline      좁은 각도(21–25도)와 넓은 각도(15–35도) 각 5개에서 노이즈 표본 3개: 보정 -> rPIE(프로브도 복원)
                -> 기판으로 위상 기울기 제거 -> 영역 사이 복소 반사율 비
  fit           구조 영역의 (질화막 두께, 산화막 두께, 질화막 밀도 배율) 맞춤: 진폭 비만 / 위상 차만 / 둘 다.
                크라메르-라오 하한(정보 행렬)과 실제 맞춤 결과
  dark          기준(기판)이 어두운 40·45도에서 같은 측정 흐름
  transmission  13.5 nm 투과 타이코그래피로 물질 가르기 (Eschen 2022 의 산란 지수 delta/beta)
"""
import os
import pickle
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
from scipy.ndimage import binary_erosion
from scipy.optimize import least_squares

import ptycho10_core as C

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, ".cache_post10.pkl")
SI, SIO2, SIN = (C.nk(m, C.WL) for m in ("Si", "SiO2", "Si3N4"))

# 시편: 자연 산화막 기판(기준), 구조(질화막/산화막), 산화막 패드
T_OX = 1.5                      # nm, 기판의 자연 산화막 (알고 있다고 둔다)
TRUE = dict(dN=10.0, dO=5.0, rho=1.0, dC=15.0)
NARROW = [21, 22, 23, 24, 25]
WIDE = [15, 20, 25, 30, 35]           # 좁은 쪽과 각도 수를 같게 두고 폭만 넓힌다
PHOTONS = 1e8                   # 위치당 광자 규모(프로브 세기 기준. 반사율을 곱한 만큼 검출된다)
OS = 4                          # 측정 흐름의 검출기 과표본 배수
N_ITER = 200
SEEDS = [0, 1, 2]


def stack_struct(dN, dO, rho):
    return [(1 - (1 - SIN) * rho, dN), (SIO2, dO)]


def reflectances(g, p=TRUE):
    sub = C.refl([(SIO2, T_OX)], SI, C.WL, g)
    st = C.refl(stack_struct(p["dN"], p["dO"], p["rho"]), SI, C.WL, g)
    pad = C.refl([(SIO2, p["dC"])], SI, C.WL, g)
    return sub, st, pad


def rects(Xmax, Ymax, scale=1.0):
    """표면 좌표(nm)의 시험 패턴. x = 입사 방향. 폭이 줄어드는 막대 두 묶음, L 자 막대, 산화막 패드.

    scale 은 막대 폭과 여백에 곱한다. 수치구경을 바꿀 때 화소 크기에 맞춰 늘려 패턴을 화소 기준으로 같게 둔다.
    """
    st = []
    y = 1500 * scale
    for w in [1200, 800, 500, 300, 200]:
        w *= scale
        st.append((2000 * scale, Xmax * 0.55, y, y + w))
        y += 2 * w
    x = Xmax * 0.6
    for w in [3000, 2000, 1200, 800, 500]:
        w *= scale
        st.append((x, x + w, 1500 * scale, Ymax * 0.45))
        x += 2 * w
    st.append((Xmax * 0.15, Xmax * 0.45, Ymax * 0.6, Ymax * 0.6 + 1500 * scale))
    pad = [(Xmax * 0.55, Xmax * 0.85, Ymax * 0.6, Ymax * 0.85)]
    return [st, pad]


def make_sample(g, na=C.NA):
    G = C.Geometry(g, na=na)
    pos = C.scan(10, 10)
    shape = C.object_shape(pos)
    sub, st, pad = reflectances(g)
    scale = G.dy / C.Geometry(g).dy
    obj, label = C.rasterize(shape, G, rects(shape[0] * G.dx, shape[1] * G.dy, scale), [st, pad], sub)
    return G, pos, shape, obj, label


# ── geometry ─────────────────────────────────────────────────
def study_geometry():
    G = C.Geometry(22)
    qe, ql, qu = G.q_exact(), G.q_linear(), G.q_uniform()
    _, mask = C.tilted_plane_correct(np.zeros((1, C.N, C.N)), G)
    gs = np.linspace(2, 60, 291)
    curves = {}
    for name, layers, sub in [("Si (자연 산화막 1.5 nm)", [(SIO2, T_OX)], SI), ("Si", [], SI),
                              ("SiO2", [], SIO2), ("Si3N4", [], SIN)]:
        r = np.array([C.parratt(layers, sub, C.WL, g) for g in gs])
        curves[name] = r
    pen = {m: C.WL / (4 * np.pi * C.nk(m, C.WL).imag) for m in ("Si", "SiO2", "Si3N4")}
    kz = lambda g: 2 * np.pi / C.WL * np.sin(np.deg2rad(g))
    height = {g: 2 * kz(g) * 0.2 for g in (21, 25, 45)}          # 높이 2 A 의 왕복 위상 (rad)
    aniso = {"dx": G.dx, "dy": G.dy, "ratio": G.dx / G.dy, "1/sin22": 1 / np.sin(np.deg2rad(22)),
             "paper_ratio": 172 / 64, "paper_angle": np.rad2deg(np.arcsin(64 / 172))}
    return dict(qe=qe, ql=ql, qu=qu, mask=mask, k=G.k, gs=gs, curves=curves, pen=pen, height=height, aniso=aniso,
                R_sub={g: abs(reflectances(g)[0]) ** 2 for g in range(15, 46)})


# ── tilt ─────────────────────────────────────────────────────
def _tilt_job(args):
    na, mode = args
    G, pos, shape, obj, label = make_sample(22, na)
    P = C.probe()
    m = C.illumination(pos, P, shape)
    m = m > 0.3 * m.max()
    if mode == "ideal":
        I, mk = C.measure(obj, P, pos, G, "uniform"), None
    else:
        Ie = C.measure(obj, P, pos, G, "exact")
        I, mk = C.tilted_plane_correct(Ie, G) if mode == "corrected" else (Ie, None)
    o, _ = C.rpie(I, pos, P, shape, 150, mask=mk, update_probe=False)
    out = dict(na=na, mode=mode, err=C.error(o, obj, m))
    if na == C.NA:
        out.update(obj=C.align(o, obj, m), true=obj, m=m, dx=G.dx, dy=G.dy,
                   mask_frac=float(mk.mean()) if mk is not None else 1.0)
    return out


def region_ratios(o, obj, label, m):
    """진폭 상호상관으로 정수 평행이동만 맞추고, 기판으로 위상 기울기·전역 위상·척도를 뺀 뒤 영역 평균의 비."""
    F = np.fft.fft2
    a, b = np.abs(obj) * m, np.abs(o) * m
    xc = np.fft.ifft2(F(a - a[m].mean() * m) * np.conj(F(b - b[m].mean() * m)))
    sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
    o = np.roll(o, sh, axis=(0, 1))
    A = binary_erosion(label == -1, iterations=3) & m
    B = binary_erosion(label == 0, iterations=2) & m
    Cp = binary_erosion(label == 1, iterations=2) & m
    o2 = C.remove_ramp(o, A)
    zA = np.mean(o2[A])
    return o2, A, B, Cp, np.mean(o2[B]) / zA, np.mean(o2[Cp]) / zA


def _bias_job(args):
    g, os_ = args
    G, pos, shape, obj, label = make_sample(g)
    P = C.probe()
    m = C.illumination(pos, P, shape)
    m = m > 0.3 * m.max()
    Gd = C.Geometry(g, os=os_)
    I, mk = C.tilted_plane_correct(C.measure(obj, P, pos, Gd, "exact"), Gd)
    o, _ = C.rpie(I, pos, P, shape, 200, mask=mk, update_probe=False)
    *_, rB, rC = region_ratios(o, obj, label, m)
    sub, st, pad = reflectances(g)
    return dict(g=g, os=os_, dphB=float(np.angle(rB / (st / sub))), dphC=float(np.angle(rC / (pad / sub))),
                dampB=float(abs(rB) / abs(st / sub) - 1), dampC=float(abs(rC) / abs(pad / sub) - 1))


# ── pipeline ─────────────────────────────────────────────────
def _pipe_job(args):
    g, seed = args
    G, pos, shape, obj, label = make_sample(g)
    P = C.probe()
    ill = C.illumination(pos, P, shape)
    m = ill > 0.3 * ill.max()
    Gd = C.Geometry(g, os=OS)
    Ie = C.measure(obj, P, pos, Gd, "exact")
    rng = np.random.default_rng(1000 * g + seed)
    # 촘촘한 화소 OS^2 개가 과표본 없는 화소 하나의 광자를 나눠 받는다
    s = PHOTONS / np.sum(np.abs(P) ** 2) / OS ** 2
    n = rng.poisson(Ie * s).astype(float)
    Ic, mask = C.tilted_plane_correct(n * OS ** 2, Gd)
    o, p = C.rpie(Ic, pos, C.probe(defocus=3.0), shape, N_ITER, mask=mask, rng=seed)
    err = C.error(o, obj, m)
    o2, A, B, Cp, rB, rC = region_ratios(o, obj, label, m)
    zA = np.mean(o2[A])
    out = dict(g=g, seed=seed, err=err, detected=float(n.sum(axis=(1, 2)).max()), rB=rB, rC=rC,
               nA=int(A.sum()), nB=int(B.sum()), nC=int(Cp.sum()))
    if seed == 0 and g in (21, 25):
        out.update(img=o2, A=A, B=B, Cp=Cp, true=obj / reflectances(g)[0], dx=G.dx, dy=G.dy)
    return out


# ── fit ──────────────────────────────────────────────────────
def model_ratio(p, angles):
    dN, dO, rho = p
    return np.array([C.refl(stack_struct(dN, dO, rho), SI, C.WL, g) / C.refl([(SIO2, T_OX)], SI, C.WL, g)
                     for g in angles])


def residuals(p, angles, meas, use, sig):
    """sig = 각도별 복소 비의 상대 노이즈 (ln|z| 와 arg z 에 같은 크기로 든다)."""
    z = meas / model_ratio(p, angles)
    parts = []
    if use in ("amp", "both"):
        parts.append(np.log(np.abs(z)) / sig)
    if use in ("phase", "both"):
        parts.append(np.angle(z) / sig)
    return np.concatenate(parts)


def crlb(angles, sig, p0=(TRUE["dN"], TRUE["dO"], TRUE["rho"])):
    p0 = np.array(p0, float)
    h = np.array([1e-3, 1e-3, 1e-5])
    r0 = model_ratio(p0, angles)
    J = np.array([(model_ratio(p0 + h * e, angles) - model_ratio(p0 - h * e, angles)) / (2 * h[i])
                  for i, e in enumerate(np.eye(3))]).T / r0[:, None]
    out = {}
    w = 1 / np.asarray(sig, float)[:, None]
    for use, JJ in [("amp", J.real * w), ("phase", J.imag * w), ("both", np.vstack([J.real * w, J.imag * w]))]:
        Fm = JJ.T @ JJ
        cov = np.linalg.pinv(Fm)
        out[use] = dict(cov=cov, sd=np.sqrt(np.diag(cov)), cond=float(np.linalg.cond(JJ)))
    return out


def fit_one(angles, meas, use, sig, rng, n_start=20):
    lo, hi = np.array([0.0, 0.0, 0.5]), np.array([30.0, 20.0, 1.5])
    best = None
    for s in range(n_start):
        x0 = lo + (hi - lo) * rng.uniform(size=3)
        r = least_squares(residuals, x0, bounds=(lo, hi), args=(angles, meas, use, sig), x_scale=[1, 1, 0.1])
        if best is None or r.cost < best.cost:
            best = r
    return best.x, float(best.cost)


def study_fit(pipe):
    out = {}
    rng = np.random.default_rng(7)
    for setname, angles in [("narrow", NARROW), ("wide", WIDE)]:
        rows = [r for r in pipe if r["g"] in angles]
        # 복소 비의 상대 노이즈: 같은 각도의 노이즈 표본 3개가 흩어진 정도(ln|z| 와 arg z 를 함께)
        sig = []
        for g in angles:
            z = np.log(np.array([r["rB"] for r in rows if r["g"] == g]))
            d = z - z.mean()
            sig.append(float(np.sqrt(np.sum(d.real ** 2 + d.imag ** 2) / (2 * (len(z) - 1)))))
        sig = np.array(sig)
        fits = {use: [] for use in ("amp", "phase", "both")}
        for seed in SEEDS:
            meas = np.array([next(r["rB"] for r in rows if r["g"] == g and r["seed"] == seed) for g in angles])
            for use in fits:
                x, cost = fit_one(angles, meas, use, sig, rng)
                fits[use].append(dict(x=x, cost=cost))
        out[setname] = dict(angles=angles, sig=sig, crlb=crlb(angles, sig), fits=fits)
    # 패드: 산화막 두께 하나만 미지수
    pad = {}
    for setname, angles in [("narrow", NARROW), ("wide", WIDE)]:
        rows = [r for r in pipe if r["g"] in angles]
        ests = []
        for seed in SEEDS:
            meas = np.array([next(r["rC"] for r in rows if r["g"] == g and r["seed"] == seed) for g in angles])
            def f(d, meas=meas):
                z = meas / np.array([C.refl([(SIO2, d[0])], SI, C.WL, g) / C.refl([(SIO2, T_OX)], SI, C.WL, g)
                                     for g in angles])
                return np.concatenate([np.log(np.abs(z)), np.angle(z)])
            ests.append(min((least_squares(f, [d0], bounds=([0], [40])) for d0 in (2, 10, 20, 30)),
                            key=lambda r: r.cost).x[0])
        pad[setname] = ests
    out["pad"] = pad
    return out


# ── transmission ─────────────────────────────────────────────
def study_transmission(seed=0):
    wl = 13.5
    G = C.Geometry(90, n=C.N)                                  # 기하는 쓰지 않는다. 화소 크기만 빌린다
    pos = C.scan(10, 8)
    shape = C.object_shape(pos)
    H, W = shape
    k = 2 * np.pi / wl
    # 물질 띠 셋(열 방향으로 나란히), 각 띠는 행 방향으로 두께가 늘어나는 쐐기. 왼쪽 띠는 비어 있다(진공)
    thick = np.zeros(shape)
    mat = np.full(shape, -1)
    yy, xx = np.indices(shape)
    bands = [("Si3N4", 30, 150), ("SiO2", 30, 150), ("Al", 10, 60)]
    edges = np.linspace(0.25 * W, W - 4, 4)
    for i, (mname, t0, t1) in enumerate(bands):
        sel = (xx >= edges[i]) & (xx < edges[i + 1])
        mat[sel] = i
        thick[sel] = t0 + (t1 - t0) * (yy[sel] - 4) / (H - 8)
    obj = np.ones(shape, complex)
    for i, (mname, _, _) in enumerate(bands):
        d, b = C.delta_beta(mname, wl)
        obj[mat == i] = np.exp(-1j * k * d * thick[mat == i]) * np.exp(-k * b * thick[mat == i])
    P = C.probe()
    I = C.measure(obj, P, pos, G, "uniform")
    rng = np.random.default_rng(seed)
    n = rng.poisson(I * PHOTONS / np.sum(np.abs(P) ** 2)).astype(float)
    o, _ = C.rpie(n, pos, C.probe(defocus=3.0), shape, N_ITER, rng=seed)
    ill = C.illumination(pos, P, shape)
    m = ill > 0.3 * ill.max()
    F = np.fft.fft2
    a, b = np.abs(obj) * m, np.abs(o) * m
    xc = np.fft.ifft2(F(a - a[m].mean() * m) * np.conj(F(b - b[m].mean() * m)))
    o = np.roll(o, np.unravel_index(np.argmax(np.abs(xc)), xc.shape), axis=(0, 1))
    vac = binary_erosion(mat == -1, iterations=3) & m
    cal = C.remove_ramp(o, vac)                                # 진공이 투과율 1 이 되도록: 전역 위상·척도·기울기 제거
    # 진공 기준이 없을 때: 전역 위상과 척도가 임의로 남는다(같은 복원에 임의의 복소 상수를 곱한 것과 같다)
    raw = cal * 0.8 * np.exp(1j * 0.6)
    res = dict(wl=wl, bands=bands, mat=mat, thick=thick, m=m, err=C.error(o, obj, m),
               fq_true={b[0]: C.delta_beta(b[0], wl)[0] / C.delta_beta(b[0], wl)[1] for b in bands})
    for name, z in [("cal", cal), ("raw", raw)]:
        sel = m & binary_erosion(mat >= 0, iterations=2)
        ph, la = np.angle(z[sel]), np.log(np.abs(z[sel]))
        res[name] = dict(phase=ph, logamp=la, mat=mat[sel], thick=thick[sel],
                         fq={i: float(np.median(ph[mat[sel] == i] / la[mat[sel] == i])) for i in range(3)})
    return res


def main():
    cache = {}
    if os.path.exists(CACHE):
        cache = pickle.load(open(CACHE, "rb"))
    todo = sys.argv[1:] or ["geometry", "tilt", "pipeline", "fit", "transmission"]
    with Pool(4) as pool:
        if "geometry" in todo:
            cache["geometry"] = study_geometry()
            print("geometry done", flush=True)
        if "tilt" in todo:
            jobs = [(na, mode) for na in (0.05, 0.10, 0.15, 0.20) for mode in ("ideal", "corrected", "uncorrected")]
            cache["tilt"] = list(pool.imap(_tilt_job, jobs))
            for r in cache["tilt"]:
                print("tilt", r["na"], r["mode"], round(r["err"], 4), flush=True)
            cache["bias"] = list(pool.imap(_bias_job, [(g, o) for g in (22, 30, 40) for o in (1, 2, 4)]))
            for r in cache["bias"]:
                print("bias", r, flush=True)
        if "pipeline" in todo:
            jobs = [(g, s) for g in sorted(set(NARROW + WIDE)) for s in SEEDS]
            cache["pipeline"] = list(pool.imap(_pipe_job, jobs))
            for r in cache["pipeline"]:
                print("pipe", r["g"], r["seed"], round(r["err"], 4), np.round(r["rB"], 4), flush=True)
        if "dark" in todo:
            # 기준 영역이 어두운 각도: 같은 측정 흐름을 40·45도에서 (노이즈 표본 3개)
            cache["dark"] = list(pool.imap(_pipe_job, [(g, sd) for g in (40, 45) for sd in SEEDS]))
            for r in cache["dark"]:
                print("dark", r["g"], r["seed"], round(r["err"], 4), flush=True)
        if "transmission" in todo:
            cache["transmission"] = study_transmission()
            print("transmission done", cache["transmission"]["err"], flush=True)
    if "fit" in todo:
        cache["fit"] = study_fit(cache["pipeline"])
        print("fit done", flush=True)
    pickle.dump(cache, open(CACHE, "wb"))


def load():
    return pickle.load(open(CACHE, "rb"))


if __name__ == "__main__":
    main()

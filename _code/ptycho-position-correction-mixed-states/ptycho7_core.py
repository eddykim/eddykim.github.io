"""7편: 주사 위치 보정과 혼합 상태(부분 결맞음). 6편의 모듈(전방 모델, rPIE, DM, ML, 노이즈)에 다음을 더했다.

  소수 화소 위치  위치 R = 정수 모서리 + 소수 부분. 시편은 정수 모서리에서 잘라내고, 프로브는 소수 부분만큼
                 푸리에 위상 기울기로 옮긴다:  psi_j = P(r - f_j) O(r + c_j)
  rpie_pos       rPIE 에 위치 보정을 붙였다.
                 anneal  위치마다 반경 r 안의 무작위 후보 넷과 현재 위치 가운데 푸리에 오차가 가장 작은 것을 고르고,
                         r 을 반복에 따라 줄인다 (Maiden 등 2012, Rodenburg·Maiden 핸드북 설명)
                 grad    진폭 비용을 위치로 미분해 가우스-뉴턴 걸음을 간다:  d = -Re<chi, D> / Re<D, D>,
                         D = d Psi / d(위치) = F[ -(grad P)(r - f) O ]
                 모든 위치를 함께 옮기는 것은 측정을 바꾸지 않으므로(공동 평행이동 모호성) 보정량의 평균은 0 으로 둔다.
  혼합 상태      서로 결맞지 않는 모드 P_k 의 세기를 더한다:  I_j = sum_k |F[P_k O_j]|^2  (Thibault, Menzel 2013)
                 넓은 광원: 광원 점마다 기울어진 평면파가 핀홀을 비춘다. 광원 점 프로브들을 SVD 해 결맞음 모드와
                 점유율(모드의 세기 몫)을 얻는다.
                 rpie_modes  푸리에 사영은 모든 모드에 같은 배율 sqrt(n / sum_k |Psi_k|^2) 를 곱한다.
                             반복마다 모드를 SVD 로 직교화한다.
                 ml_modes    진폭 비용의 ML 정련 (모드 판).
"""
import warnings

import numpy as np
from scipy.ndimage import gaussian_filter
from skimage import data, transform

warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

N = 64
WL = 633e-9


# ── FFT 와 전파 ──────────────────────────────────────────────
def fft2c(x):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def ifft2c(X):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(X, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def _kz(n, dx, wl=WL):
    f = (np.arange(n) - n // 2) / (n * dx)
    FX, FY = np.meshgrid(f, f)
    return FX, FY, np.sqrt(np.maximum(1 / wl**2 - FX**2 - FY**2, 0))


def propagate(u, dx, z, wl=WL):
    """각스펙트럼 전파."""
    *_, kz = _kz(u.shape[-1], dx, wl)
    return ifft2c(fft2c(u) * np.exp(2j * np.pi * z * kz))


def disc(n, radius_px):
    y, x = np.indices((n, n)) - n // 2
    return (np.hypot(x, y) <= radius_px).astype(complex)


# ── 프로브 ───────────────────────────────────────────────────
def pinhole_probe(n=N, radius_px=12, z=200e-6, dx=1e-6):
    """반지름 12 um 핀홀 뒤 200 um 의 장 (3편과 같다)."""
    return propagate(disc(n, radius_px), dx, z)


def diffuser_probe(n=N, radius_px=12, corr_px=4.0, rms_phase=2.5, z=80e-6, dx=1e-6, seed=3):
    """무작위 위상판을 붙인 핀홀 뒤의 장. 구조가 복잡한 프로브로, Maiden 2017 이 ePIE 가 고전한다고 보인 경우다.

    위상판은 상관 길이 corr_px 의 매끄러운 가우시안 잡음이고 표준편차가 rms_phase [rad] 다.
    빛이 계산 창(64 화소) 밖으로 새지 않도록 상관 길이와 전파 거리를 골랐다 (창 안 에너지 99 %).
    """
    ph = gaussian_filter(np.random.default_rng(seed).normal(size=(n, n)), corr_px)
    ph = rms_phase * ph / ph.std()
    return propagate(disc(n, radius_px) * np.exp(1j * ph), dx, z)


def focused_probe(n=N, na=0.05, defocus=100e-6, dx=0.5e-6):
    """개구수 na 렌즈의 초점에서 defocus 만큼 벗어난 면의 장. 파면이 휘어 있다."""
    FX, FY, kz = _kz(n, dx)
    pupil = (FX**2 + FY**2 <= (na / WL) ** 2) * np.exp(2j * np.pi * defocus * kz)
    p = ifft2c(pupil.astype(complex))
    return p / np.sqrt(np.sum(np.abs(p) ** 2))


def probe_d90(probe):
    n = probe.shape[-1]
    y, x = np.indices((n, n)) - n // 2
    r = np.hypot(x, y)
    I = np.abs(probe) ** 2
    rs = np.arange(1, n // 2)
    frac = np.array([I[r < q].sum() for q in rs]) / I.sum()
    return 2 * int(rs[np.searchsorted(frac, 0.9)]) if frac[-1] >= 0.9 else None


# ── 시편과 주사 ──────────────────────────────────────────────
def _norm01(a):
    a = a.astype(float)
    return (a - a.min()) / (a.max() - a.min())


def make_object(shape, amp_range=(0.2, 1.0), phase_max=np.pi):
    n = max(shape)
    amp = _norm01(transform.resize(data.camera().astype(float), (n, n), anti_aliasing=True))
    rgb = data.astronaut().astype(float)
    gray = 0.2125 * rgb[..., 0] + 0.7154 * rgb[..., 1] + 0.0721 * rgb[..., 2]
    ph = _norm01(transform.resize(gray, (n, n), anti_aliasing=True))
    lo, hi = amp_range
    return ((lo + (hi - lo) * amp) * np.exp(1j * phase_max * ph))[: shape[0], : shape[1]]


def scan(kind, step=8, n_side=12, seed=0, margin=4):
    """같은 위치 수(n_side^2)의 주사. kind = "regular" | "jitter" | "fermat".

    페르마 나선은 r = c sqrt(k), theta = k x 황금각. 점 하나가 차지하는 면적 pi c^2 를 step^2 로 맞춘다.
    """
    m = n_side * n_side
    if kind == "fermat":
        k = np.arange(m)
        c = step / np.sqrt(np.pi)
        th = k * np.pi * (3 - np.sqrt(5))
        pos = np.stack([c * np.sqrt(k) * np.sin(th), c * np.sqrt(k) * np.cos(th)], 1)
    else:
        iy, ix = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
        pos = np.stack([iy.ravel(), ix.ravel()], 1) * float(step)
        if kind == "jitter":
            pos = pos + np.random.default_rng(seed).integers(-1, 2, pos.shape)
    pos = np.rint(pos).astype(int)
    return pos - pos.min(0) + margin


def object_shape(pos, n=N, margin=4):
    return tuple(int(v) for v in pos.max(0) + n + margin)


def illumination(pos, probe, shape):
    n = probe.shape[-1]
    ill = np.zeros(shape)
    for y, x in pos:
        ill[y:y + n, x:x + n] += np.abs(probe) ** 2
    return ill


def intensities(obj, probe, pos):
    n = probe.shape[-1]
    return np.array([np.abs(fft2c(probe * obj[y:y + n, x:x + n])) ** 2 for y, x in pos])


# ── 복원 ─────────────────────────────────────────────────────
def fourier_project(psi, sqrt_I):
    return ifft2c(sqrt_I * np.exp(1j * np.angle(fft2c(psi))))


def centroid(f):
    w = np.abs(f) ** 2
    y, x = np.indices(f.shape)
    return np.array([(w * y).sum(), (w * x).sum()]) / w.sum()


def maiden(T=30, eta_obj=0.9, eta_prb=0.9):
    """Maiden 2017 의 모멘텀: 마찰과 피드백이 같은 eta."""
    return dict(T=T, fric_obj=eta_obj, feed_obj=eta_obj, fric_prb=eta_prb, feed_prb=eta_prb)


def ptylab(T=20, fric=0.7, feed=0.3):
    """PtyLab mPIE 의 기본값 (마찰 0.7, 피드백 0.3). PtyLab 은 위치마다 5 % 확률로 거는데 여기서는 20 회마다로 둔다."""
    return dict(T=T, fric_obj=fric, feed_obj=feed, fric_prb=fric, feed_prb=feed)


def rpie(I, pos, probe0, shape, n_iter, alpha=0.05, beta=1.0, gamma_obj=1.0, gamma_prb=1.0,
         momentum=None, probe_power=None, probe_start=2, obj0=None, center_probe=False, rng=0, callback=None):
    """rPIE / mPIE / ePIE (alpha = beta = 1). 반환: (obj, probe).

    gamma_obj 는 수이거나 반복 번호 k -> 걸음 의 함수다.
    momentum 은 maiden() / ptylab() 이 만드는 사전 (T, fric_*, feed_*).
    probe_power="auto" 이면 프로브 전체 세기를 가장 밝은 패턴의 총세기에 묶고 시편에 역수를 곱한다.
    (cO, P/c) 는 측정을 바꾸지 않으므로 갱신이 이 방향으로 떠내려갈 수 있는데, 모멘텀은 이 평평한
    방향의 움직임까지 가속한다. 파세발 정리로 투과율 1 인 곳의 패턴 총세기가 프로브 세기다.
    """
    rng = np.random.default_rng(rng)
    n = probe0.shape[-1]
    sI = np.sqrt(np.maximum(I, 0))
    obj = np.ones(shape, complex) if obj0 is None else obj0.astype(complex).copy()
    probe = probe0.astype(complex).copy()
    c0 = np.rint(centroid(probe0)).astype(int)
    g_o = gamma_obj if callable(gamma_obj) else (lambda k: gamma_obj)
    if probe_power == "auto":
        probe_power = float(np.max(np.sum(np.maximum(I, 0), axis=(-2, -1))))
    if momentum:
        T = momentum["T"]
        v_o, v_p = np.zeros_like(obj), np.zeros_like(probe)
        o_prev, p_prev = obj.copy(), probe.copy()

    def fix_scale():
        if probe_power is None:
            return 1.0
        s = np.sqrt(probe_power / np.sum(np.abs(probe) ** 2))
        probe[...] *= s
        obj[...] /= s
        return s

    count = 0
    for k in range(1, n_iter + 1):
        go = g_o(k)
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            view = obj[y:y + n, x:x + n]
            o_old = view.copy()
            psi = probe * o_old
            diff = fourier_project(psi, sI[j]) - psi
            a2 = np.abs(probe) ** 2
            view += go * np.conj(probe) * diff / ((1 - alpha) * a2 + alpha * a2.max())
            if k >= probe_start:
                b2 = np.abs(o_old) ** 2
                probe += gamma_prb * np.conj(o_old) * diff / ((1 - beta) * b2 + beta * b2.max())
            count += 1
            if momentum and count % T == 0:
                s = fix_scale()
                o_prev, v_o, p_prev, v_p = o_prev / s, v_o / s, p_prev * s, v_p * s
                v_o = momentum["fric_obj"] * v_o + (obj - o_prev)
                obj += momentum["feed_obj"] * v_o
                o_prev = obj.copy()
                if k >= probe_start:
                    v_p = momentum["fric_prb"] * v_p + (probe - p_prev)
                    probe += momentum["feed_prb"] * v_p
                p_prev = probe.copy()
        if not momentum:
            fix_scale()
        if center_probe:
            d = np.rint(centroid(probe)).astype(int) - c0
            if d.any():
                probe[...] = np.roll(probe, tuple(-d), axis=(0, 1))
                obj[...] = np.roll(obj, tuple(-d), axis=(0, 1))
                if momentum:
                    v_o, o_prev = np.roll(v_o, tuple(-d), axis=(0, 1)), np.roll(o_prev, tuple(-d), axis=(0, 1))
                    v_p, p_prev = np.roll(v_p, tuple(-d), axis=(0, 1)), np.roll(p_prev, tuple(-d), axis=(0, 1))
        if callback:
            callback(k, obj, probe)
    return obj, probe


def overlap_projection(psi, obj, probe, pos, update_probe=True, n_inner=3, eps=1e-8):
    """sum_j || psi_j - P O_j ||^2 를 O, P 에 대해 번갈아 최소화한다 (DM 의 중첩 사영)."""
    n = probe.shape[-1]
    for _ in range(n_inner):
        num = np.zeros_like(obj)
        den = np.zeros(obj.shape)
        for (y, x), ps in zip(pos, psi):
            num[y:y + n, x:x + n] += np.conj(probe) * ps
            den[y:y + n, x:x + n] += np.abs(probe) ** 2
        obj = num / (den + eps * den.max())
        if update_probe:
            views = np.stack([obj[y:y + n, x:x + n] for y, x in pos])
            dp = np.sum(np.abs(views) ** 2, 0)
            probe = np.sum(np.conj(views) * psi, 0) / (dp + eps * dp.max())
    return obj, probe


def dm(I, pos, probe0, shape, n_iter, probe_start=2, n_inner=3, callback=None):
    n = probe0.shape[-1]
    sI = np.sqrt(np.maximum(I, 0))
    obj = np.ones(shape, complex)
    probe = probe0.astype(complex).copy()
    psi = np.stack([probe * obj[y:y + n, x:x + n] for y, x in pos])
    for k in range(1, n_iter + 1):
        obj, probe = overlap_projection(psi, obj, probe, pos, k >= probe_start, n_inner)
        p_o = np.stack([probe * obj[y:y + n, x:x + n] for y, x in pos])
        p_f = np.stack([fourier_project(2 * a - b, s) for a, b, s in zip(p_o, psi, sI)])
        psi = psi + p_f - p_o
        if callback:
            callback(k, obj, probe)
    return overlap_projection(psi, obj, probe, pos, True, n_inner)


# ── 평가 ─────────────────────────────────────────────────────
def align(est, ref, mask=None, shift=True, ramp=True):
    """모호성을 맞춘다: 정수 평행이동, 위상 기울기, 복소 척도. 각각 끌 수 있다."""
    mask = np.ones(est.shape, bool) if mask is None else mask
    e = est
    if shift:
        F = np.fft.fft2
        xc = np.fft.ifft2(F(ref * mask) * np.conj(F(est * mask)))
        sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
        e = np.roll(est, sh, axis=(0, 1))
    if ramp:
        z = np.where(mask, np.conj(e) * ref, 0)
        gx = np.angle(np.sum(z[:, 1:] * np.conj(z[:, :-1])))
        gy = np.angle(np.sum(z[1:, :] * np.conj(z[:-1, :])))
        Y, X = np.indices(e.shape)
        e = e * np.exp(1j * (gx * X + gy * Y))
    c = np.vdot(e[mask], ref[mask]) / np.vdot(e[mask], e[mask])
    return c * e


def error(est, ref, mask=None, **kw):
    mask = np.ones(est.shape, bool) if mask is None else mask
    a = align(est, ref, mask, **kw)
    return float(np.linalg.norm((a - ref)[mask]) / np.linalg.norm(ref[mask]))


def residual(obj, probe, pos, I):
    """측정 진폭과 추정이 내는 진폭의 상대 차이. 0 이면 추정이 데이터를 완전히 설명한다."""
    n = probe.shape[-1]
    r = t = 0.0
    for (y, x), Ij in zip(pos, I):
        A = np.abs(fft2c(probe * obj[y:y + n, x:x + n]))
        r += np.sum((A - np.sqrt(np.maximum(Ij, 0))) ** 2)
        t += np.sum(np.maximum(Ij, 0))
    return float(np.sqrt(r / t))


def grid_power_fraction(err_img, step):
    """오차 영상의 전력 중 주사 격자 주파수(1/step 의 배수) 칸에 몰린 몫과, 그 칸이 전체에서 차지하는 비율."""
    E = np.abs(np.fft.fft2(err_img)) ** 2
    f = np.fft.fftfreq(err_img.shape[0])
    harm = np.array([m / step for m in range(1, step // 2 + 1)])
    on = lambda v: (np.min(np.abs(np.abs(v)[..., None] - harm), axis=-1) < 0.5 / err_img.shape[0]) | (v == 0)
    FX, FY = np.meshgrid(f, f)
    g = on(FX) & on(FY) & ~((FX == 0) & (FY == 0))
    return float(E[g].sum() / E.sum()), float(g.mean())


# ── 6편: 노이즈와 최대가능도 ─────────────────────────────────
def noisy(I, photons, rng):
    """가장 밝은 패턴의 총광자 수가 photons 가 되도록 세기를 맞추고 푸아송 표본을 뽑는다. 반환: (n, 배율)."""
    s = photons / I.sum(axis=(-2, -1)).max()
    return rng.poisson(I * s).astype(float), s


def views(obj, pos, n=N):
    return np.stack([obj[y:y + n, x:x + n] for y, x in pos])


def cost(Psi, n, model):
    I = np.abs(Psi) ** 2
    if model == "amp":
        return float(np.sum((np.sqrt(I) - np.sqrt(n)) ** 2))
    if model == "poisson":
        return float(np.sum(I - n * np.log(I + 1e-12)))
    if model == "int":
        return float(0.5 * np.sum((I - n) ** 2))
    raise ValueError(model)


def chi(Psi, n, model):
    """비용을 Psi 의 켤레로 미분한 것 (Wirtinger). 푸리에 공간에서 Psi 에 곱해지는 꼴."""
    I = np.abs(Psi) ** 2
    if model == "amp":
        return (1 - np.sqrt(n) / np.maximum(np.sqrt(I), 1e-12)) * Psi
    if model == "poisson":
        return (1 - n / np.maximum(I, 1e-12)) * Psi
    if model == "int":
        return (I - n) * Psi
    raise ValueError(model)


def total_cost(obj, probe, pos, n, model):
    return cost(fft2c(probe * views(obj, pos, probe.shape[-1])), n, model)


def grad_obj(obj, probe, pos, n, model):
    """시편에 대한 기울기와 화소마다 받은 조명 sum|P|^2. fft2c 가 유니터리라 역변환이 곧 수반 연산이다."""
    k = probe.shape[-1]
    X = ifft2c(chi(fft2c(probe * views(obj, pos, k)), n, model))
    g = np.zeros(obj.shape, complex)
    d = np.zeros(obj.shape)
    for (y, x), xj in zip(pos, X):
        g[y:y + k, x:x + k] += np.conj(probe) * xj
        d[y:y + k, x:x + k] += np.abs(probe) ** 2
    return g, d


def grad_probe(obj, probe, pos, n, model):
    V = views(obj, pos, probe.shape[-1])
    X = ifft2c(chi(fft2c(probe * V), n, model))
    return np.sum(np.conj(V) * X, 0), np.sum(np.abs(V) ** 2, 0)


def _line_search(f, x, dx, c0, a=1.0, shrink=0.5, tries=25):
    for _ in range(tries):
        x1 = x + a * dx
        c1 = f(x1)
        if c1 < c0:
            return x1, c1, a
        a *= shrink
    return x, c0, 0.0


def ml(n, pos, probe0, shape, n_iter, model="amp", obj0=None, update_probe=False, probe_start=1, callback=None):
    """최대가능도 복원 (배치 경사하강 + 선 탐색). 반환: (obj, probe).

    방향은 기울기를 화소별 조명 sum|P|^2 (프로브는 sum|O|^2) 로 나눈 것이다 (Thibault 2012 의 전처리와 같은 취지).
    선 탐색은 비용이 실제로 줄어드는 걸음만 받아들이므로, 척도 모호성처럼 비용이 평평한 방향으로는 움직이지 않는다.
    """
    obj = np.ones(shape, complex) if obj0 is None else obj0.astype(complex).copy()
    probe = probe0.astype(complex).copy()
    c = total_cost(obj, probe, pos, n, model)
    for k in range(1, n_iter + 1):
        g, d = grad_obj(obj, probe, pos, n, model)
        obj, c, _ = _line_search(lambda o: total_cost(o, probe, pos, n, model), obj, -g / (d + 1e-3 * d.max()), c)
        if update_probe and k >= probe_start:
            g, d = grad_probe(obj, probe, pos, n, model)
            probe, c, _ = _line_search(lambda p: total_cost(obj, p, pos, n, model), probe, -g / (d + 1e-3 * d.max()), c)
        if callback:
            callback(k, obj, probe)
    return obj, probe


# ── 7편: 소수 화소 위치와 위치 보정 ───────────────────────────
_F = np.fft.fftfreq(N)
_FY, _FX = np.meshgrid(_F, _F, indexing="ij")


def shift(a, d):
    """a(r - d): d = (dy, dx) 화소만큼 순환 이동 (소수 가능)."""
    return np.fft.ifft2(np.fft.fft2(a) * np.exp(-2j * np.pi * (_FY * d[0] + _FX * d[1])))


def dshift(a, d):
    """a(r - d) 를 d 의 두 성분으로 미분한 것: -(grad a)(r - d). 반환 (2, N, N)."""
    A = np.fft.fft2(a) * np.exp(-2j * np.pi * (_FY * d[0] + _FX * d[1]))
    return np.stack([np.fft.ifft2(-2j * np.pi * _FY * A), np.fft.ifft2(-2j * np.pi * _FX * A)])


def split(p):
    c = np.floor(p).astype(int)
    return c, p - c


def object_shape_f(pos, n=N, margin=6):
    return tuple(int(v) for v in np.ceil(np.max(pos, 0)) + n + margin)


def intensities_sub(obj, probe, pos):
    out = []
    for p in pos:
        (y, x), f = split(p)
        out.append(np.abs(fft2c(shift(probe, f) * obj[y:y + N, x:x + N])) ** 2)
    return np.array(out)


def fourier_err(obj, probe, p, sI):
    (y, x), f = split(p)
    return float(np.sum((np.abs(fft2c(shift(probe, f) * obj[y:y + N, x:x + N])) - sI) ** 2))


def rpie_pos(I, pos0, probe0, shape, n_iter, correct=None, alpha=0.05, r0=1.0, start=10, n_cand=4, mu=0.5,
             max_step=0.5, update_probe=True, rng=0, callback=None):
    """위치 보정을 붙인 rPIE. 반환: (obj, probe, positions). correct = None | "anneal" | "grad"."""
    rng = np.random.default_rng(rng)
    sI = np.sqrt(np.maximum(I, 0))
    obj = np.ones(shape, complex)
    probe = probe0.astype(complex).copy()
    pos = np.asarray(pos0, float).copy()
    for k in range(1, n_iter + 1):
        r = r0 * max(0.0, 1 - (k - start) / (0.7 * n_iter))
        for j in rng.permutation(len(pos)):
            if correct == "anneal" and k >= start and r > 0:
                cands = [pos[j]] + [pos[j] + rng.uniform(-r, r, 2) for _ in range(n_cand)]
                pos[j] = cands[int(np.argmin([fourier_err(obj, probe, c, sI[j]) for c in cands]))]
            (y, x), f = split(pos[j])
            Pj = shift(probe, f)
            o_old = obj[y:y + N, x:x + N].copy()
            psi = Pj * o_old
            if correct == "grad" and k >= start:
                Psi = fft2c(psi)
                D = [fft2c(dP * o_old) for dP in dshift(probe, f)]
                chi = (1 - sI[j] / np.maximum(np.abs(Psi), 1e-12)) * Psi
                g = np.array([np.real(np.vdot(Dk, chi)) for Dk in D])
                h = np.array([np.real(np.vdot(Dk, Dk)) for Dk in D])
                pos[j] = pos[j] - np.clip(mu * g / h, -max_step, max_step)
            diff = fourier_project(psi, sI[j]) - psi
            a2 = np.abs(Pj) ** 2
            obj[y:y + N, x:x + N] = o_old + np.conj(Pj) * diff / ((1 - alpha) * a2 + alpha * a2.max())
            if update_probe and k >= 2:
                probe = probe + shift(np.conj(o_old) * diff / np.max(np.abs(o_old) ** 2), -f)
        if correct:
            pos -= (pos - pos0).mean(0)          # 공동 평행이동 모호성: 보정량의 평균은 0
        if callback:
            callback(k, obj, probe, pos)
    return obj, probe, pos


# ── 7편: 혼합 상태 ─────────────────────────────────────────────
def source_probes(sig_t, n_side=9, radius_px=12, z=200e-6, dx=1e-6):
    """가우시안 넓은 광원. 광원 점(기울기 sig_t 사이클/화소의 가우시안, ±2.5 sig 를 n_side 점)마다 기울어진
    평면파가 반지름 radius_px 핀홀을 비추고 z 만큼 전파한다. 반환: (점별 프로브, 무게)."""
    y, x = np.indices((N, N)) - N // 2
    if not sig_t:
        return np.array([propagate(disc(N, radius_px), dx, z)]), np.array([1.0])
    t = np.linspace(-2.5, 2.5, n_side) * sig_t
    Ps, ws = [], []
    for ty in t:
        for tx in t:
            ws.append(np.exp(-(ty ** 2 + tx ** 2) / (2 * sig_t ** 2)))
            Ps.append(propagate(disc(N, radius_px) * np.exp(2j * np.pi * (ty * y + tx * x)), dx, z))
    ws = np.array(ws)
    return np.array(Ps), ws / ws.sum()


def coherent_modes(Ps, ws):
    """광원 점 프로브들의 SVD. 반환: (모드, 점유율). sum_k |모드_k|^2 = sum_s w_s |P_s|^2 이다."""
    M = (np.sqrt(ws)[:, None] * Ps.reshape(len(ws), -1)).T
    U, S, _ = np.linalg.svd(M, full_matrices=False)
    return (U * S).T.reshape(-1, N, N), S ** 2 / np.sum(S ** 2)


def modal_intensities(obj, modes, pos, weights=None):
    w = np.ones(len(modes)) if weights is None else weights
    return np.array([np.sum(w[:, None, None] * np.abs(fft2c(modes * obj[y:y + N, x:x + N])) ** 2, 0) for y, x in pos])


def ortho(modes):
    K = len(modes)
    U, S, _ = np.linalg.svd(modes.reshape(K, -1).T, full_matrices=False)
    return (U * S).T.reshape(K, N, N)


def modal_project(Psi, sI):
    """모드 묶음 Psi (K, N, N) 의 푸리에 사영: 모든 모드에 같은 배율 sqrt(n) / sqrt(sum_k |Psi_k|^2)."""
    return Psi * (sI / np.maximum(np.sqrt(np.sum(np.abs(Psi) ** 2, 0)), 1e-12))


def residual_modes(obj, modes, pos, I):
    r = t = 0.0
    for (y, x), Ij in zip(pos, I):
        Im = np.sum(np.abs(fft2c(modes * obj[y:y + N, x:x + N])) ** 2, 0)
        r += np.sum((np.sqrt(Im) - np.sqrt(Ij)) ** 2)
        t += np.sum(Ij)
    return float(np.sqrt(r / t))


def rpie_modes(I, pos, probe0, K, shape, n_iter, alpha=0.05, rng=0, seed_modes=7, callback=None):
    """혼합 상태 rPIE. 첫 모드는 probe0, 나머지는 probe0 의 1/10 크기에 무작위 위상. 반환: (obj, modes)."""
    rng = np.random.default_rng(rng)
    sI = np.sqrt(np.maximum(I, 0))
    obj = np.ones(shape, complex)
    modes = np.zeros((K, N, N), complex)
    modes[0] = probe0
    g = np.random.default_rng(seed_modes)
    for k in range(1, K):
        modes[k] = probe0 * 0.1 * np.exp(2j * np.pi * g.uniform(size=(N, N)))
    for it in range(1, n_iter + 1):
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            o = obj[y:y + N, x:x + N].copy()
            psi = modes * o
            diff = ifft2c(modal_project(fft2c(psi), sI[j])) - psi
            a2 = np.sum(np.abs(modes) ** 2, 0)
            obj[y:y + N, x:x + N] = o + np.sum(np.conj(modes) * diff, 0) / ((1 - alpha) * a2 + alpha * a2.max())
            if it >= 2:
                modes = modes + np.conj(o) * diff / np.max(np.abs(o) ** 2)
        modes = ortho(modes)
        if callback:
            callback(it, obj, modes)
    return obj, modes


def _mode_cost(obj, modes, pos, sI):
    V = np.stack([obj[y:y + N, x:x + N] for y, x in pos])
    Im = np.sum(np.abs(fft2c(modes[None] * V[:, None])) ** 2, 1)
    return float(np.sum((np.sqrt(Im) - sI) ** 2))


def ml_modes(I, pos, obj0, modes0, n_iter, callback=None):
    """혼합 상태의 진폭 비용 ML 정련 (배치 경사하강 + 선 탐색). 반환: (obj, 직교화한 modes)."""
    sI = np.sqrt(np.maximum(I, 0))
    obj, modes = obj0.copy(), modes0.copy()
    c = _mode_cost(obj, modes, pos, sI)
    for it in range(1, n_iter + 1):
        for part in ("obj", "modes"):
            V = np.stack([obj[y:y + N, x:x + N] for y, x in pos])
            Psi = fft2c(modes[None] * V[:, None])
            Im = np.sum(np.abs(Psi) ** 2, 1)
            X = ifft2c((1 - sI / np.maximum(np.sqrt(Im), 1e-12))[:, None] * Psi)
            if part == "obj":
                g = np.zeros(obj.shape, complex)
                d = np.zeros(obj.shape)
                for (y, x), xj in zip(pos, X):
                    g[y:y + N, x:x + N] += np.sum(np.conj(modes) * xj, 0)
                    d[y:y + N, x:x + N] += np.sum(np.abs(modes) ** 2, 0)
                obj, c, _ = _line_search(lambda o: _mode_cost(o, modes, pos, sI), obj, -g / (d + 1e-3 * d.max()), c)
            else:
                g = np.sum(np.conj(V)[:, None] * X, 0)
                d = np.sum(np.abs(V) ** 2, 0)
                modes, c, _ = _line_search(lambda m: _mode_cost(obj, m, pos, sI), modes,
                                           -g / (d + 1e-3 * d.max())[None], c)
        if callback:
            callback(it, obj, modes)
    return obj, ortho(modes)

"""6편: 노이즈와 최대가능도. 5편의 모듈(전방 모델, 프로브, 주사, rPIE/mPIE, DM, 오차)에 다음을 더했다.

  noisy      가장 밝은 패턴의 총광자 수를 정하고 푸아송 표본을 뽑는다
  비용       (Thibault, Guizar-Sicairos 2012)  Psi = F[P O_j],  I = |Psi|^2,  n = 측정 광자 수
             amp      sum (|Psi| - sqrt(n))^2               진폭 최소자승. sqrt(n) 의 분산은 약 1/4 로 고르다
             poisson  sum (I - n log I)                      푸아송 음의 로그가능도 (상수 제외)
             int      sum (I - n)^2 / 2                      세기 최소자승
             기울기는 모두 Psi 에 chi 를 곱한 꼴이다:
             amp (1 - sqrt(n)/|Psi|),  poisson (1 - n/I),  int (I - n)
  ml         시편(과 프로브)을 한꺼번에 고치는 경사하강. 방향은 sum|P|^2 로 나눈 기울기, 걸음은 비용이
             줄어들 때까지 반으로 줄이는 선 탐색으로 정한다.
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

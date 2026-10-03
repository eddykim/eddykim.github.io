"""타이코그래피의 전방 모델과 PIE·ePIE — 계산 이미징과 타이코그래피 3편.

측정 (얇은 시편, 원시야)
    psi_j(r) = P(r - R_j) O(r)          위치 j 의 출사파 = 프로브 x 시편
    I_j(k)   = | F{psi_j}(k) |^2         위치마다 회절 패턴 한 장

PIE  (Rodenburg & Faulkner 2004)  프로브를 안다.
    O_j += (|P| / |P|max) * conj(P) / (|P|^2 + alpha) * beta * (psi' - psi)
ePIE (Maiden & Rodenburg 2009)    프로브도 모른다.
    O_j += a * conj(P)   / max|P|^2   * (psi' - psi)
    P   += b * conj(O_j) / max|O_j|^2 * (psi' - psi)

psi' 는 psi 의 회절면 크기만 측정값으로 바꾼 출사파다 (2편의 모듈러스 사영).
위치는 정수 화소의 왼쪽 위 모서리다.
"""
import warnings

import numpy as np
from skimage import data, transform

warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

N_PROBE = 64
DX = 1e-6            # 시편면 화소 [m]
WL = 633e-9          # 파장 [m]


# ── FFT ──────────────────────────────────────────────────────
def fft2c(x):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def ifft2c(X):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(X, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


# ── 조명과 시편 ───────────────────────────────────────────────
def circular_aperture(n, radius):
    y, x = np.indices((n, n)) - n // 2
    return (np.hypot(x, y) <= radius).astype(complex)


def pinhole_probe(n=N_PROBE, radius=12e-6, z=200e-6, dx=DX, wl=WL):
    """반지름 radius 핀홀 뒤 z 에서의 장. 각스펙트럼 전파로 계산한다."""
    aperture = circular_aperture(n, radius / dx)
    f = (np.arange(n) - n // 2) / (n * dx)
    FX, FY = np.meshgrid(f, f)
    kz = np.sqrt(np.maximum(1 / wl**2 - FX**2 - FY**2, 0))
    return ifft2c(fft2c(aperture) * np.exp(2j * np.pi * z * kz))


def _norm01(a):
    a = a.astype(float)
    return (a - a.min()) / (a.max() - a.min())


def _gray(rgb):
    rgb = rgb.astype(float)
    return 0.2125 * rgb[..., 0] + 0.7154 * rgb[..., 1] + 0.0721 * rgb[..., 2]


def make_object(shape, amp_range=(0.2, 1.0), phase_max=np.pi):
    """진폭은 camera, 위상은 astronaut. 서로 다른 영상을 실어 누화를 드러낸다 (Maiden 2009 와 같은 범위)."""
    n = max(shape)
    amp = _norm01(transform.resize(data.camera().astype(float), (n, n), anti_aliasing=True))
    ph = _norm01(transform.resize(_gray(data.astronaut()), (n, n), anti_aliasing=True))
    lo, hi = amp_range
    return ((lo + (hi - lo) * amp) * np.exp(1j * phase_max * ph))[: shape[0], : shape[1]]


def raster(step, span=84, jitter=1, seed=0, margin=4):
    """span 화소를 덮는 정사각 격자. 규칙 격자의 인공 무늬를 피하려 +-jitter 화소 흔든다 (Maiden 2009)."""
    m = int(np.ceil(span / step)) + 1
    iy, ix = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
    pos = np.stack([iy.ravel(), ix.ravel()], 1) * step
    pos = pos + np.random.default_rng(seed).integers(-jitter, jitter + 1, pos.shape)
    return pos - pos.min(0) + margin


def object_shape(pos, n=N_PROBE, margin=4):
    return tuple(int(v) for v in pos.max(0) + n + margin)


def illumination(pos, probe, shape):
    """위치마다 받은 조명 세기의 합. 복원을 평가할 영역(충분히 조명된 곳)을 정한다."""
    n = probe.shape[-1]
    ill = np.zeros(shape)
    for y, x in pos:
        ill[y:y + n, x:x + n] += np.abs(probe) ** 2
    return ill


def intensities(obj, probe, pos):
    n = probe.shape[-1]
    return np.array([np.abs(fft2c(probe * obj[y:y + n, x:x + n])) ** 2 for y, x in pos])


def poisson(I, photons_per_pattern, rng=None):
    """패턴당 평균 광자 수에 맞춰 척도를 바꾸고 푸아송 노이즈를 넣은 뒤 원래 척도로 돌린다."""
    scale = photons_per_pattern / I.sum(axis=(1, 2)).mean()
    return np.random.default_rng(rng).poisson(I * scale) / scale


def probe_d90(probe):
    """프로브 에너지의 90 % 를 담는 원의 지름 [화소]. 중첩률을 정의하는 데 쓴다."""
    n = probe.shape[-1]
    y, x = np.indices((n, n)) - n // 2
    r = np.hypot(x, y)
    I = np.abs(probe) ** 2
    rs = np.arange(1, n // 2)
    frac = np.array([I[r < q].sum() for q in rs]) / I.sum()
    return 2 * int(rs[np.searchsorted(frac, 0.9)])


def linear_overlap(step, diameter):
    return 1 - step / diameter


# ── 복원 ─────────────────────────────────────────────────────
def fourier_project(psi, sqrt_I):
    Psi = fft2c(psi)
    return ifft2c(sqrt_I * np.exp(1j * np.angle(Psi)))


def pie_weight(probe, alpha=1e-4):
    """PIE 의 갱신 가중 (|P|/|P|max) conj(P) / (|P|^2 + alpha |P|max^2)."""
    pmax2 = np.max(np.abs(probe)) ** 2
    return np.abs(probe) / np.sqrt(pmax2) * np.conj(probe) / (np.abs(probe) ** 2 + alpha * pmax2)


def pie(I, pos, probe, shape, n_iter, alpha=1e-4, beta=1.0, obj0=None, rng=0, callback=None):
    rng = np.random.default_rng(rng)
    n = probe.shape[-1]
    sI = np.sqrt(np.maximum(I, 0))
    obj = np.ones(shape, complex) if obj0 is None else obj0.astype(complex).copy()
    w = pie_weight(probe, alpha)
    for k in range(1, n_iter + 1):
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            view = obj[y:y + n, x:x + n]
            psi = probe * view
            view += w * beta * (fourier_project(psi, sI[j]) - psi)
        if callback:
            callback(k, obj)
    return obj


def epie(I, pos, probe0, shape, n_iter, a=1.0, b=1.0, probe_start=2, obj0=None, rng=0, callback=None):
    rng = np.random.default_rng(rng)
    n = probe0.shape[-1]
    sI = np.sqrt(np.maximum(I, 0))
    obj = np.ones(shape, complex) if obj0 is None else obj0.astype(complex).copy()
    probe = probe0.astype(complex).copy()
    for k in range(1, n_iter + 1):
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            view = obj[y:y + n, x:x + n]
            o_old = view.copy()
            psi = probe * o_old
            diff = fourier_project(psi, sI[j]) - psi
            view += a * np.conj(probe) / np.max(np.abs(probe)) ** 2 * diff
            if k >= probe_start:
                probe += b * np.conj(o_old) / np.max(np.abs(o_old)) ** 2 * diff
        if callback:
            callback(k, obj, probe)
    return obj, probe


# ── 평가 ─────────────────────────────────────────────────────
def align(est, ref, mask=None):
    """타이코그래피의 모호성을 맞춘다: 정수 평행이동, 위상 기울기, 복소 척도.

    물체와 프로브를 함께 옮기거나, 물체에 위상 기울기 exp(i g.r) 를 곱하고 프로브에서 빼도
    측정은 그대로다 (위치마다 상수 위상만 달라진다). 4편에서 자세히 다룬다.
    """
    mask = np.ones(est.shape, bool) if mask is None else mask
    F = np.fft.fft2
    xc = np.fft.ifft2(F(ref * mask) * np.conj(F(est * mask)))
    sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
    e = np.roll(est, sh, axis=(0, 1))
    z = np.where(mask, np.conj(e) * ref, 0)
    gx = np.angle(np.sum(z[:, 1:] * np.conj(z[:, :-1])))
    gy = np.angle(np.sum(z[1:, :] * np.conj(z[:-1, :])))
    Y, X = np.indices(e.shape)
    e = e * np.exp(1j * (gx * X + gy * Y))
    c = np.vdot(e[mask], ref[mask]) / np.vdot(e[mask], e[mask])
    return c * e


def error(est, ref, mask=None):
    mask = np.ones(est.shape, bool) if mask is None else mask
    a = align(est, ref, mask)
    return float(np.linalg.norm((a - ref)[mask]) / np.linalg.norm(ref[mask]))


# ── 비교용: 단일 패턴 CDI (2편의 HIO + 마지막 ER) ─────────────
def cdi_hio(obj_box, sup, photons_total, n_iter=600, er_tail=50, beta=0.9, seed=0, noise_seed=7):
    F = lambda x: np.fft.fft2(x, norm="ortho")
    iF = lambda x: np.fft.ifft2(x, norm="ortho")
    I = np.abs(F(obj_box)) ** 2
    scale = photons_total / I.sum()
    mag = np.sqrt(np.random.default_rng(noise_seed).poisson(I * scale) / scale)
    x = iF(mag * np.exp(2j * np.pi * np.random.default_rng(seed).random(mag.shape)))
    for k in range(n_iter):
        pm = iF(mag * np.exp(1j * np.angle(F(x))))
        x = np.where(sup, pm, 0) if k >= n_iter - er_tail else np.where(sup, pm, x - beta * pm)
    return x


def cdi_error(est, obj_box):
    """2편과 같은 정렬: 평행이동, 쌍둥이 상(공액 반전), 복소 척도."""
    F = lambda x: np.fft.fft2(x, norm="ortho")
    iF = lambda x: np.fft.ifft2(x, norm="ortho")
    best = np.inf
    for cand in (est, np.conj(np.roll(np.flip(est), 1, axis=(0, 1)))):
        xc = iF(F(obj_box) * np.conj(F(cand)))
        sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
        c2 = np.roll(cand, sh, axis=(0, 1))
        a = np.vdot(c2, obj_box) / np.vdot(c2, c2)
        best = min(best, np.linalg.norm(a * c2 - obj_box) / np.linalg.norm(obj_box))
    return float(best)

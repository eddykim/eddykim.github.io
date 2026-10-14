"""10편: 극자외선 이미징 반사계. 이전 편의 rPIE·정렬에 다음을 더했다.

  광학 상수   CXRO Henke 표 (data/henke_*.dat), 굴절률 n = 1 - delta + i beta (진행파 exp(i(kz - wt)))
  반사        Parratt 재귀로 s 편광 반사계수. 모든 영역을 기판 윗면(z = 0)에 맞춘 공통 기준면으로 옮겨
              높이 차이가 위상으로 들어가게 한다:  r_meas = r_stack exp(-2i kz0 H)
  반사 기하   표면과 grazing 각 g 로 입사, 검출기는 정반사 방향에 수직으로 거리 L. 검출기 화소 (u, v) 가
              받는 표면 공간주파수 q = (k_out - k_in) 의 표면 성분은 u 에 대해 비선형이다.
              중심 근처에서 q_x = k sin g u / L, q_y = k v / L  -> x 방향 화소가 1/sin g 배 길다.
  전방 모델   출사파 psi(x, y) 의 푸리에 변환을 검출기 화소가 받는 정확한 q 에서 직접 합(비균일 DFT)으로 구한다.
  보정        기울어진 평면 보정: 균일한 q 격자의 각 점이 검출기의 어디에 떨어지는지 역사상으로 찾고
              측정 세기를 보간한다. 검출기 밖에 떨어지는 점은 마스크로 뺀다.
"""
import os
import warnings

import numpy as np
from scipy.ndimage import map_coordinates

warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

HERE = os.path.dirname(os.path.abspath(__file__))
N = 64                 # 검출기 화소 수 = 시편 창 화소 수
WL = 29.3              # nm, Tanksalvala 2021 의 파장
NA = 0.20              # 검출기 가장자리까지의 수치구경 (v 방향)


# ── 광학 상수 ────────────────────────────────────────────────
def nk(material, wl):
    d = np.loadtxt(os.path.join(HERE, "data", f"henke_{material}.dat"), skiprows=2)
    return 1 - np.interp(wl, d[:, 0], d[:, 1]) + 1j * np.interp(wl, d[:, 0], d[:, 2])


def delta_beta(material, wl):
    n = nk(material, wl)
    return 1 - n.real, n.imag


# ── Parratt 반사 ─────────────────────────────────────────────
def parratt(layers, sub, wl, g_deg):
    """layers = [(n, d), ...] 위에서 아래로, sub = 기판 굴절률. s 편광 반사계수(맨 윗면 기준)."""
    k = 2 * np.pi / wl
    c = np.cos(np.deg2rad(g_deg))
    ns = [1.0 + 0j] + [n for n, _ in layers] + [sub]
    kz = [k * np.sqrt(n ** 2 - c ** 2 + 0j) for n in ns]
    r = (kz[-2] - kz[-1]) / (kz[-2] + kz[-1])
    for j in range(len(layers) - 1, -1, -1):          # 층 j+1 (ns 기준) 의 두께로 위상을 쌓는다
        d = layers[j][1]
        e = np.exp(2j * kz[j + 1] * d)
        rj = (kz[j] - kz[j + 1]) / (kz[j] + kz[j + 1])
        r = (rj + r * e) / (1 + rj * r * e)
    return r


def refl(layers, sub, wl, g_deg):
    """기판 윗면(z = 0)을 공통 기준면으로 둔 반사계수. 층이 기판 위에 쌓인 높이 H 가 위상 -2 k sin g H 로 들어간다."""
    H = sum(d for _, d in layers)
    kz0 = 2 * np.pi / wl * np.sin(np.deg2rad(g_deg))
    return parratt(layers, sub, wl, g_deg) * np.exp(-2j * kz0 * H)


def transfer_matrix_refl(layers, sub, wl, g_deg):
    """검증용: 특성 행렬(Abeles) 방법으로 같은 s 편광 반사계수를 구한다(맨 윗면 기준)."""
    k = 2 * np.pi / wl
    c = np.cos(np.deg2rad(g_deg))
    q = lambda n: np.sqrt(n ** 2 - c ** 2 + 0j)
    M = np.eye(2, dtype=complex)
    for n, d in layers:
        qn = q(n)
        b = k * qn * d
        M = M @ np.array([[np.cos(b), -1j * np.sin(b) / qn], [-1j * qn * np.sin(b), np.cos(b)]])
    q0, qs = q(1.0 + 0j), q(sub)
    num = q0 * M[0, 0] + q0 * qs * M[0, 1] - M[1, 0] - qs * M[1, 1]
    den = q0 * M[0, 0] + q0 * qs * M[0, 1] + M[1, 0] + qs * M[1, 1]
    return num / den


# ── 반사 기하 ────────────────────────────────────────────────
class Geometry:
    """grazing 각 g 의 반사 타이코그래피. 길이 단위는 nm, 검출기 거리 L = 1 로 둔 무차원 좌표."""

    def __init__(self, g_deg, n=N, na=NA, wl=WL, fit_inside=False, os=1):
        """os = 검출기 과표본 배수. 검출기 화소를 os 배 촘촘히 두고(n*os 개), 복원 격자는 n 개로 둔다."""
        self.g, self.n, self.wl, self.os = g_deg, n, wl, os
        self.k = 2 * np.pi / wl
        th = np.deg2rad(g_deg)
        self.s = np.array([np.cos(th), 0.0, np.sin(th)])     # 정반사 방향
        self.eu = np.array([-np.sin(th), 0.0, np.cos(th)])   # 입사면 안, 정반사에 수직
        self.ev = np.array([0.0, 1.0, 0.0])
        self.p = np.tan(np.arcsin(na)) / (n / 2)             # 복원 격자에 대응하는 화소 간격 (L = 1)
        self.nd, self.pd = n * os, self.p / os                # 실제 검출기 화소 수와 간격
        c = (np.arange(self.nd) - self.nd // 2) * self.pd
        self.U, self.V = np.meshgrid(c, c, indexing="ij")    # 0 축 = u (입사 방향), 1 축 = v
        # 균일 격자: 중심 화소 간격에 맞춘 q 간격과 그에 대응하는 시편 화소 크기
        self.dqx = self.k * np.sin(th) * self.p
        self.dqy = self.k * self.p
        if fit_inside:
            # 균일 격자를 검출기가 빈틈없이 덮는 범위로 줄인다(측정이 없는 칸을 없앤다). 입사 방향 분해능을 내준다
            qx, qy = self.q_exact()
            half = n // 2
            self.dqx = min(-qx.min() / half, qx.max() / (half - 1)) * 0.995
            self.dqy = min(np.abs(qy[:, 0]).min(), np.abs(qy[:, -1]).min()) / half * 0.995
        self.dx = 2 * np.pi / (n * self.dqx)                 # 입사 방향 화소 (nm)
        self.dy = 2 * np.pi / (n * self.dqy)

    def q_exact(self):
        """검출기 화소 중심이 받는 표면 공간주파수 (q_x, q_y). q_x 는 중심 근처에서 +k sin g u 가 되도록 부호를 맞췄다."""
        R = self.s[None, None, :] + self.U[..., None] * self.eu + self.V[..., None] * self.ev
        kout = self.k * R / np.linalg.norm(R, axis=-1, keepdims=True)
        return -(kout[..., 0] - self.k * self.s[0]), kout[..., 1]

    def q_linear(self):
        return self.U * self.k * np.sin(np.deg2rad(self.g)), self.V * self.k

    def q_uniform(self):
        c = np.arange(self.n) - self.n // 2
        return np.meshgrid(c * self.dqx, c * self.dqy, indexing="ij")

    def q_to_detector(self, qx, qy):
        """역사상: 표면 공간주파수 -> 검출기 좌표 (u, v). 검출기에 닿지 않으면 nan."""
        kx = -qx + self.k * self.s[0]
        ky = qy
        kz2 = self.k ** 2 - kx ** 2 - ky ** 2
        ok = kz2 > 0
        kz = np.sqrt(np.where(ok, kz2, 0))
        d = np.stack([kx, ky, kz], -1) / self.k
        t = 1.0 / (d @ self.s)
        R = d * t[..., None]
        u, v = R @ self.eu, R @ self.ev
        u[~ok] = np.nan
        v[~ok] = np.nan
        return u, v


def dft_factors(geom, q):
    """출사파 (n x n, 화소 dx x dy) 의 연속 푸리에 변환을 주어진 q 점에서 구하는 두 인수 (ortho 규격).

    Psi[p] = sum_x sum_y Ex[p, x] psi[x, y] Ey[p, y]. q 점마다 (q_x, q_y) 짝이 달라 행렬 하나로 분리되지 않으므로
    psi @ Ey^T 를 먼저 구하고 Ex 와 점별로 곱해 더한다.
    """
    qx, qy = q
    c = np.arange(geom.n) - geom.n // 2
    Ex = np.exp(-1j * qx.reshape(-1, 1) * (c * geom.dx)[None, :]) / np.sqrt(geom.n)
    Ey = np.exp(-1j * qy.reshape(-1, 1) * (c * geom.dy)[None, :]) / np.sqrt(geom.n)
    return Ex, Ey


def nudft(psi, factors):
    Ex, Ey = factors
    B = psi @ Ey.T                                   # (n, 화소 수)
    return np.einsum("px,xp->p", Ex, B)


def tilted_plane_correct(I, geom, order=3):
    """검출기 세기 I (n x n) 를 균일한 q 격자로 보간한다. 반환: (보정한 세기, 마스크)."""
    qx, qy = geom.q_uniform()
    u, v = geom.q_to_detector(qx, qy)
    iu = u / geom.pd + geom.nd // 2
    iv = v / geom.pd + geom.nd // 2
    mask = np.isfinite(iu) & (iu >= 0) & (iu <= geom.nd - 1) & (iv >= 0) & (iv <= geom.nd - 1)
    out = np.zeros((len(I),) + qx.shape)
    coords = np.stack([np.where(mask, iu, 0), np.where(mask, iv, 0)])
    for j, Ij in enumerate(I):
        out[j] = np.where(mask, map_coordinates(Ij, coords, order=order, mode="nearest"), 0)
    return np.maximum(out, 0), mask


# ── 시편과 프로브 ────────────────────────────────────────────
def fft2c(x):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def ifft2c(X):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(X, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def probe(n=N, radius=0.35, defocus=6.0, seed=None):
    """화소 좌표의 집속 프로브: 원형 동공(반지름 = 나이퀴스트의 radius 배)에 초점 어긋남 위상."""
    f = (np.arange(n) - n // 2) / n
    FX, FY = np.meshgrid(f, f, indexing="ij")
    rho2 = (FX ** 2 + FY ** 2) / (0.5 * radius) ** 2
    pupil = (rho2 <= 1) * np.exp(1j * defocus * rho2)
    return ifft2c(pupil) * n / np.sqrt(np.sum(np.abs(pupil) ** 2))


def scan(n_side, step, margin=4, seed=0):
    rng = np.random.default_rng(seed)
    iy, ix = np.meshgrid(np.arange(n_side), np.arange(n_side), indexing="ij")
    pos = np.stack([iy.ravel(), ix.ravel()], 1) * float(step)
    pos = np.rint(pos + rng.uniform(-1, 1, pos.shape)).astype(int)    # 격자 무늬를 피하려고 흔든다(4편)
    return pos - pos.min(0) + margin


def object_shape(pos, n=N, margin=4):
    return tuple(int(v) for v in pos.max(0) + n + margin)


def rasterize(shape, geom, regions, values, background, sub=4):
    """표면 좌표(nm)의 사각형 영역들을 시편 격자에 면적 가중으로 그린다.

    regions = [[(x0, x1, y0, y1), ...], ...] 영역마다 사각형 목록, values = 영역마다 복소 반사율.
    반환: (복소 시편, 영역 번호 지도 (-1 = 배경, 화소가 한 영역에 완전히 들 때만)).
    """
    H, W = shape
    o = (np.arange(sub) + 0.5) / sub - 0.5
    X = (np.arange(H)[:, None] + o[None, :]).ravel() * geom.dx
    Y = (np.arange(W)[:, None] + o[None, :]).ravel() * geom.dy
    obj = np.full(shape, background, complex)
    label = np.full(shape, -1)
    for r, (rects, val) in enumerate(zip(regions, values)):
        inside = np.zeros((X.size, Y.size), bool)
        for x0, x1, y0, y1 in rects:
            inside |= ((X >= x0) & (X < x1))[:, None] & ((Y >= y0) & (Y < y1))[None, :]
        frac = inside.reshape(H, sub, W, sub).mean(axis=(1, 3))
        obj = obj * (1 - frac) + val * frac
        label[frac == 1] = r
    return obj, label


# ── 측정 ─────────────────────────────────────────────────────
def measure(obj, prb, pos, geom, mode="exact"):
    """검출기 세기. mode="exact" 는 화소가 받는 정확한 q, "uniform" 은 균일 격자(=보정이 이상적일 때)."""
    n = geom.n
    if mode == "uniform":
        return np.array([np.abs(fft2c(prb * obj[y:y + n, x:x + n])) ** 2 for y, x in pos])
    fac = dft_factors(geom, geom.q_exact())
    out = np.empty((len(pos), geom.nd, geom.nd))
    for j, (y, x) in enumerate(pos):
        out[j] = np.abs(nudft(prb * obj[y:y + n, x:x + n], fac)).reshape(geom.nd, geom.nd) ** 2
    return out


def noisy(I, photons, rng):
    """가장 밝은 패턴의 총광자 수를 photons 로 맞춘 푸아송 표본. 반환: (n, 배율)."""
    s = photons / I.sum(axis=(-2, -1)).max()
    return rng.poisson(I * s).astype(float), s


# ── 복원 (rPIE, 푸리에 마스크) ───────────────────────────────
def rpie(I, pos, probe0, shape, n_iter, mask=None, alpha=0.05, beta=1.0, gamma_obj=1.0, gamma_prb=1.0,
         update_probe=True, probe_start=2, obj0=None, rng=0):
    """rPIE. mask 가 있으면 측정이 없는 푸리에 화소는 모형 값을 그대로 둔다. 반환: (obj, probe)."""
    rng = np.random.default_rng(rng)
    n = probe0.shape[-1]
    sI = np.sqrt(np.maximum(I, 0))
    m = np.ones((n, n), bool) if mask is None else mask
    obj = np.ones(shape, complex) if obj0 is None else obj0.astype(complex).copy()
    prb = probe0.astype(complex).copy()
    power = float(np.max(np.sum(np.maximum(I, 0) * m, axis=(-2, -1))))
    for k in range(1, n_iter + 1):
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            view = obj[y:y + n, x:x + n]
            o_old = view.copy()
            psi = prb * o_old
            Psi = fft2c(psi)
            Psi = np.where(m, sI[j] * np.exp(1j * np.angle(Psi)), Psi)
            diff = ifft2c(Psi) - psi
            a2 = np.abs(prb) ** 2
            view += gamma_obj * np.conj(prb) * diff / ((1 - alpha) * a2 + alpha * a2.max())
            if update_probe and k >= probe_start:
                b2 = np.abs(o_old) ** 2
                prb += gamma_prb * np.conj(o_old) * diff / ((1 - beta) * b2 + beta * b2.max())
        if update_probe:                                     # 척도 모호성 고정(4편)
            s = np.sqrt(power / np.sum(np.abs(fft2c(prb)) ** 2 * m))
            prb *= s
            obj /= s
    return obj, prb


def illumination(pos, prb, shape):
    n = prb.shape[-1]
    ill = np.zeros(shape)
    for y, x in pos:
        ill[y:y + n, x:x + n] += np.abs(prb) ** 2
    return ill


def align(est, ref, mask, shift=True, ramp=True):
    """정수 평행이동, 위상 기울기, 복소 척도를 맞춘다.

    이전 편과 달리 평행이동은 진폭의 상호상관으로 찾는다. 프로브가 복원 중 푸리에 공간에서 옮겨 가면
    시편에 큰 위상 기울기가 남는데, 복소 상호상관은 그 기울기 때문에 엉뚱한 봉우리를 고른다.
    """
    e = est
    if shift:
        F = np.fft.fft2
        a, b = np.abs(ref) * mask, np.abs(est) * mask
        xc = np.fft.ifft2(F(a - a[mask].mean() * mask) * np.conj(F(b - b[mask].mean() * mask)))
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


def error(est, ref, mask, **kw):
    a = align(est, ref, mask, **kw)
    return float(np.linalg.norm((a - ref)[mask]) / np.linalg.norm(ref[mask]))


def remove_ramp(obj, region):
    """한 영역(평평한 기판)의 위상이 평평해지도록 위상 기울기와 전역 위상을 뺀다. 실제 측정에서 쓰는 기준면 맞춤.

    기울기가 여러 바퀴 감겨 있을 수 있으므로 위상을 직접 평면 맞춤하지 않고, 영역 안의 이웃 화소 쌍의
    위상차(감기지 않는다)로 기울기를 구한다. 남은 작은 기울기는 한 번 더 같은 방법으로 뺀다.
    """
    Y, X = np.indices(obj.shape)
    out = obj.copy()
    for _ in range(2):
        z = np.where(region, out, 0)
        py = region[1:, :] & region[:-1, :]
        px = region[:, 1:] & region[:, :-1]
        gy = np.angle(np.sum((z[1:, :] * np.conj(z[:-1, :]))[py]))
        gx = np.angle(np.sum((z[:, 1:] * np.conj(z[:, :-1]))[px]))
        out = out * np.exp(-1j * (gy * Y + gx * X))
    out = out * np.exp(-1j * np.angle(np.mean(out[region])))
    return out / np.mean(np.abs(out[region]))

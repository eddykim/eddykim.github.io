"""11편: 실측 푸리에 타이코그래피. 8편의 FP rPIE 에 실측에 필요한 것을 더했다.

  데이터      UConn Smart Imaging Lab 공개 데이터(Song & Zheng), PtyLab FPM 공개 데이터(Loetgering 등 2023).
              get_data.py 가 data/ 에 받는다. 데이터는 이 폴더에 싣지 않는다.
  조명 기하   LED 위치 -> 시편에 닿는 광선의 공기 중 sin(각) = 조명 NA -> 조명 파수(푸리에 화소).
              시편이 유리 기판 위에 있으면 광선이 유리에서 꺾여 들어오는 점이 옆으로 밀린다.
  동공        대물렌즈 NA 의 원에 초점 어긋남 위상 exp(i 2 pi z sqrt(1/lam^2 - f^2)).
  복원        fp_cal: FP rPIE + LED 밝기 배율 + 동공 복원(EPRY). 배율은 반복마다 중앙값을 1 로 묶는다
              (묶지 않으면 배율과 스펙트럼이 척도를 주고받아 배율이 0 으로 흘러간다 — 4편의 척도 모호성).
  판정        residual: 영상마다 척도 하나를 자유로 둔 진폭 잔차. noise_floor: 세기 스펙트럼의 2NA 원 밖에서 잰 노이즈.
              frc: 두 독립 복원의 푸리에 고리 상관.
"""
import os
import warnings

import numpy as np
import scipy.io as sio
from scipy.optimize import brentq

warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")

# UConn 장치 (data_description.txt)
UC = dict(NA=0.1, px=1.845e-6, H=90.88, pitch=4.0, glass_t=1.0, glass_n=1.52, n_side=15)


def fft2c(x):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def ifft2c(X):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(X, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


# ── 조명 기하 ────────────────────────────────────────────────
def spiral(n_side):
    """가운데 LED 에서 바깥으로 도는 나선의 격자 좌표 (dx, dy). 첫 걸음 +x, 다음 +y, 걸음 길이 1,1,2,2,3,3,..."""
    out = [(0, 0)]
    x = y = 0
    step, d = 1, 0
    dirs = [(1, 0), (0, 1), (-1, 0), (0, -1)]
    while len(out) < n_side ** 2:
        for _ in range(2):
            for _ in range(step):
                x, y = x + dirs[d][0], y + dirs[d][1]
                out.append((x, y))
            d = (d + 1) % 4
        step += 1
    return np.array(out[:n_side ** 2], float)


def na_through_glass(l, H, t, n):
    """축에서 l 떨어진 LED 가 높이 H 에 있고, 시편 아래에 두께 t, 굴절률 n 의 유리가 있을 때 시편에 닿는 광선의
    공기 중 sin(각). 공기 중 옆 이동 a 와 유리 속 옆 이동 t tan(theta_g) 의 합이 l 이고, sin(theta_a) = n sin(theta_g)."""
    if l == 0:
        return 0.0
    f = lambda a: a + t * np.tan(np.arcsin(np.sin(np.arctan(a / H)) / n)) - l
    return float(np.sin(np.arctan(brentq(f, 0, l) / H)))


def defocus_pupil(nl, r, z, wl, fov):
    """반지름 r (푸리에 화소) 원형 동공에 초점 어긋남 z 의 각스펙트럼 위상. fov = 영상 한 변의 시편면 길이."""
    y, x = np.indices((nl, nl)) - nl // 2
    P = ((x ** 2 + y ** 2) <= r ** 2).astype(complex)
    if z:
        f2 = (x ** 2 + y ** 2) / fov ** 2
        P *= np.exp(2j * np.pi * z * np.sqrt(np.maximum(1 / wl ** 2 - f2, 0)))
    return P


# ── 데이터 읽기 ──────────────────────────────────────────────
def load_uconn(name, glass=True, swap=False, sx=-1, sy=-1, as_intensity=False):
    """UConn .mat. imlow_HDR 는 진폭으로 저장되어 있다(as_intensity=True 면 일부러 세기로 읽는다).
    반환 dict: I (J, nl, nl), k (J, 2) 푸리에 화소 (ky, kx), r 동공 반지름, na_ill, wl, z, ab(보정 동공 또는 None)."""
    d = sio.loadmat(os.path.join(DATA, name))
    A = np.moveaxis(d["imlow_HDR"].astype(float), -1, 0)
    sq = lambda k: float(np.squeeze(d[k]))
    wl, th, xi, yi, z = sq("wlength"), sq("theta"), sq("xint"), sq("yint"), sq("z")
    ab = d["aberration"]
    nl = A.shape[-1] - A.shape[-1] % 2
    A = A[:, :nl, :nl]
    sp = spiral(UC["n_side"])
    if swap:
        sp = sp[:, ::-1]
    sp = sp * [sx, sy]
    x0, y0 = xi + sp[:, 0] * UC["pitch"], yi + sp[:, 1] * UC["pitch"]
    c, s = np.cos(np.deg2rad(th)), np.sin(np.deg2rad(th))
    x1, y1 = x0 * c - y0 * s, x0 * s + y0 * c
    l = np.hypot(x1, y1)
    na = np.array([na_through_glass(li, UC["H"], UC["glass_t"], UC["glass_n"]) if glass else li / np.hypot(li, UC["H"])
                   for li in l])
    fov = nl * UC["px"]
    k = np.stack([-na * y1, -na * x1], 1) / np.maximum(l, 1e-12)[:, None] / wl * fov
    return dict(I=A if as_intensity else A ** 2, k=k, r=UC["NA"] / wl * fov, nl=nl, wl=wl, z=z, fov=fov, na_ill=na,
                NA=UC["NA"], ab=ab.astype(complex) if ab.size > 1 else None, sp=sp)


def load_ptylab(name, bg=True):
    """PtyLab FPM hdf5. encoder = LED 의 (x, y) [m], zled = LED–시편 거리. 바닥값(1 백분위)을 뺀다."""
    import h5py
    with h5py.File(os.path.join(DATA, name)) as h:
        I = h["ptychogram"][()].astype(float)
        enc, z, wl = h["encoder"][()], float(h["zled"][()]), float(h["wavelength"][()])
        NA, dxd, M = float(h["NA"][()]), float(h["dxd"][()]), float(h["magnification"][()])
    b = np.percentile(I, 1)
    if bg:
        I = np.maximum(I - b, 0)
    nl = I.shape[-1]
    R = np.sqrt(enc[:, 0] ** 2 + enc[:, 1] ** 2 + z ** 2)
    fov = nl * dxd / M
    na = np.hypot(enc[:, 0], enc[:, 1]) / R
    k = np.stack([enc[:, 0] / R, enc[:, 1] / R], 1) / wl * fov
    return dict(I=I, k=k, r=NA / wl * fov, nl=nl, wl=wl, z=0.0, fov=fov, na_ill=na, NA=NA, ab=None, bg=b)


# ── 복원 ─────────────────────────────────────────────────────
def band(nh, nl, kv):
    c = nh // 2 + np.asarray(kv)
    h = nl // 2
    return (slice(c[0] - h, c[0] + h), slice(c[1] - h, c[1] + h))


def fp_cal(I, ks, P0, nh, n_iter, alpha=0.1, beta=1.0, pupil_from=None, bright_from=None, callback=None):
    """FP rPIE (밝은 명시야부터) + LED 밝기 배율(bright_from 번째 반복부터) + 동공 복원 EPRY(pupil_from 부터).
    반환: (고해상 스펙트럼 S, 동공 P, 영상별 밝기 배율 c). 배율은 측정 진폭에 곱해 모형 척도로 맞추는 값이다."""
    ks = np.asarray(ks, int)
    nl = P0.shape[-1]
    a = np.sqrt(np.maximum(I, 0))
    c = np.ones(len(ks))
    S = np.zeros((nh, nh), complex)
    S[band(nh, nl, (0, 0))] = fft2c(a[int(np.argmin(np.sum(ks ** 2, 1)))])
    P = P0.astype(complex).copy()
    supp = np.abs(P0) > 0
    seq = np.argsort(np.sum(ks ** 2, 1))
    for it in range(1, n_iter + 1):
        for j in seq:
            sl = band(nh, nl, ks[j])
            s_old = S[sl].copy()
            psi = P * s_old
            u = ifft2c(psi)
            if bright_from and it >= bright_from:
                c[j] = np.sum(np.abs(u) * a[j]) / np.sum(a[j] ** 2)
            diff = fft2c(c[j] * a[j] * np.exp(1j * np.angle(u))) - psi
            a2 = np.abs(P) ** 2
            S[sl] = s_old + np.conj(P) * diff / ((1 - alpha) * a2 + alpha * a2.max())
            if pupil_from and it >= pupil_from:
                b2 = np.abs(s_old) ** 2
                P = (P + np.conj(s_old) * diff / ((1 - beta) * b2 + beta * b2.max())) * supp
        if bright_from and it >= bright_from:
            c /= np.median(c)                     # 척도 모호성: 배율의 중앙값을 1 로 묶는다
        if callback:
            callback(it, S, P, c)
    return S, P, c


def model_amps(S, P, ks):
    nh, nl = S.shape[-1], P.shape[-1]
    return np.array([np.abs(ifft2c(P * S[band(nh, nl, kv)])) for kv in np.asarray(ks, int)])


def residual(S, P, ks, I, per_image=False, crop=0):
    """진폭 잔차 sqrt(sum (s m - a)^2 / sum a^2). 영상마다 척도 s 하나는 자유(LED 밝기)로 둔다."""
    m, a = model_amps(S, P, ks), np.sqrt(np.maximum(I, 0))
    if crop:
        m, a = m[:, crop:-crop, crop:-crop], a[:, crop:-crop, crop:-crop]
    s = np.sum(m * a, (1, 2)) / np.sum(m * m, (1, 2))
    r2 = np.sum((s[:, None, None] * m - a) ** 2, (1, 2))
    if per_image:
        return np.sqrt(r2 / np.sum(a ** 2, (1, 2)))
    return float(np.sqrt(r2.sum() / np.sum(a ** 2)))


def refine_k(I, ks, P0, nh, rounds=4, R=2, its=10, k_nominal=None, **kw):
    """LED 마다 조명 파수를 주변 (2R+1)^2 정수 화소에서 다시 고른다(7편 어닐링의 푸리에판). 공동 이동은 되돌린다."""
    ks = np.asarray(ks, int).copy()
    k0 = ks.copy() if k_nominal is None else np.rint(k_nominal).astype(int)
    nl = P0.shape[-1]
    a = np.sqrt(np.maximum(I, 0))
    S, P, c = fp_cal(I, ks, P0, nh, its, **kw)
    hist = [residual(S, P, ks, I)]
    for _ in range(rounds):
        for j in range(len(ks)):
            best, kb = np.inf, ks[j]
            for dy in range(-R, R + 1):
                for dx in range(-R, R + 1):
                    kv = ks[j] + [dy, dx]
                    m = np.abs(ifft2c(P * S[band(nh, nl, kv)]))
                    s = np.sum(m * a[j]) / np.sum(m * m)
                    e = np.sum((s * m - a[j]) ** 2)
                    if e < best:
                        best, kb = e, kv
            ks[j] = kb
        ks = ks - np.rint(np.median(ks - k0, 0)).astype(int)
        S, P, c = fp_cal(I, ks, P0, nh, its, **kw)
        hist.append(residual(S, P, ks, I))
    return ks, S, P, hist


# ── 판정 ─────────────────────────────────────────────────────
def noise_floor(I, r, margin=2):
    """세기 영상 |u|^2 의 스펙트럼은 동공의 자기상관이라 반지름 2r 안에만 있다. 그 밖의 에너지를 백색 노이즈로 보고
    전체로 늘린 세기 노이즈 에너지 E_n 을 진폭 상대 노이즈로 옮긴다: 진폭 노이즈 ~ 세기 노이즈 / (2 sqrt(I))."""
    nl = I.shape[-1]
    y, x = np.indices((nl, nl)) - nl // 2
    out = (x ** 2 + y ** 2) > (2 * r + margin) ** 2
    F = np.abs(fft2c(I)) ** 2
    En = F[:, out].mean(1) * nl * nl
    return np.sqrt(En / (4 * I.sum((1, 2)) * I.mean((1, 2))))


def frc(a, b):
    """푸리에 고리 상관. 반환: (고리 주파수 [주기/화소], 값)."""
    n = a.shape[-1]
    A, B = fft2c(a), fft2c(b)
    y, x = np.indices(a.shape) - n // 2
    r = np.hypot(x, y).astype(int)
    nr = n // 2
    num = np.bincount(r.ravel(), (A * np.conj(B)).real.ravel(), minlength=nr + 1)[:nr] \
        + 1j * np.bincount(r.ravel(), (A * np.conj(B)).imag.ravel(), minlength=nr + 1)[:nr]
    da = np.bincount(r.ravel(), (np.abs(A) ** 2).ravel(), minlength=nr + 1)[:nr]
    db = np.bincount(r.ravel(), (np.abs(B) ** 2).ravel(), minlength=nr + 1)[:nr]
    return (np.arange(nr) + 0.5) / n, np.abs(num) / np.sqrt(da * db + 1e-30)


def pupil_phase_rms(P, ref, supp):
    """지지 안에서 전역 위상을 뺀 두 동공의 위상 차 rms."""
    z = P[supp] * np.conj(ref[supp])
    return float(np.sqrt(np.mean(np.angle(z * np.exp(-1j * np.angle(z.sum()))) ** 2)))

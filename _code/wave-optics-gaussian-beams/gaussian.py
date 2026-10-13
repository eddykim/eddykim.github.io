"""파동광학 배경이론 5편의 계산 모듈: 가우시안 빔의 q 파라미터와 ABCD 법칙, 파동 전파, 공진기의 폭스-리 반복, 고차 모드, 빔 품질.

길이는 µm 로 쓴다. 복소 진폭은 exp(+ikz) 로 나아가는 꼴이다 (시간 인자 e^{−iωt}).
이 규약에서 가우시안 빔은
    U(r, z) = (q0 / q) · exp(i k r² / 2q) · exp(ikz),   q = z − i z0,   z0 = π W0² / λ
이고, 1/q = 1/R + i λ/(π W²) 이다.

두 갈래의 계산을 둔다.
  (가) q 파라미터와 ABCD 행렬로 빔의 폭·곡률·허리를 계산하는 해석 경로
  (나) 장을 격자에 놓고 각스펙트럼 방법(3편)과 렌즈 위상으로 직접 전파하는 파동 경로
verify_gaussian.py 는 두 경로가 같은 답을 주는지, 그리고 교재의 수치와 맞는지 본다.
"""
import numpy as np
from scipy.special import eval_genlaguerre, eval_hermite


# ---------------------------------------------------------------- (가) q 파라미터와 ABCD

def rayleigh_range(W0, lam):
    return np.pi * W0**2 / lam


def q_at_waist(W0, lam):
    """허리에서의 q = −i z0 (이 규약)."""
    return -1j * rayleigh_range(W0, lam)


def q_from_W_R(W, R, lam):
    """폭 W(1/e² 반지름), 곡률 반지름 R 에서 q. R = inf 면 평면 파면."""
    inv = (0.0 if np.isinf(R) else 1.0 / R) + 1j * lam / (np.pi * W**2)
    return 1.0 / inv


def W_R_from_q(q, lam):
    inv = 1.0 / np.asarray(q, complex)
    W = np.sqrt(lam / (np.pi * inv.imag))
    with np.errstate(divide="ignore"):
        R = np.where(np.abs(inv.real) > 0, 1.0 / inv.real, np.inf)
    return W, R


def waist_from_q(q, lam):
    """q = z − i z0 꼴로 읽어 허리까지의 거리(양수면 허리를 지나온 것)와 허리 반지름."""
    q = complex(q)
    z, z0 = q.real, -q.imag
    return z, np.sqrt(lam * z0 / np.pi)


def M_free(d):
    return np.array([[1.0, d], [0.0, 1.0]])


def M_lens(f):
    return np.array([[1.0, 0.0], [-1.0 / f, 1.0]])


def M_mirror(R):
    """곡률 반지름 R(오목이면 양수)인 거울: 초점 거리 R/2 의 렌즈와 같다."""
    return M_lens(R / 2.0)


def abcd(M, q):
    """q2 = (A q1 + B)/(C q1 + D)."""
    A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    return (A * q + B) / (C * q + D)


def chain(*mats):
    """빛이 지나는 순서대로 준 행렬들을 곱한다 (마지막에 지나는 것이 왼쪽)."""
    M = np.eye(2)
    for m in mats:
        M = m @ M
    return M


def resonator_round_trip(L, R1, R2):
    """거울 1 바로 앞에서 출발해 거울 2, 거울 1 에서 반사되어 돌아오는 한 왕복의 행렬."""
    return chain(M_free(L), M_mirror(R2), M_free(L), M_mirror(R1))


def resonator_mode_q(L, R1, R2, lam):
    """한 왕복 뒤 자기 자신으로 돌아오는 q (q = (Aq+B)/(Cq+D) 의 해 가운데 허수부가 음수, 곧 빔이 유한한 것)."""
    M = resonator_round_trip(L, R1, R2)
    A, B, C, D = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    # C q² + (D − A) q − B = 0
    roots = np.roots([C, D - A, -B]) if abs(C) > 0 else np.array([np.nan + 0j])
    good = [r for r in roots if np.imag(r) < 0]
    return good[0] if good else None


def symmetric_mode_waist(L, R, lam):
    """대칭 공진기(두 거울 곡률 R, 거리 L)의 모드 허리 W0 를 '거울 자리(허리에서 L/2)의 파면 곡률 = R' 조건으로 푼다.
    q 법칙으로 곡률을 계산해 근을 찾는다. 공초점(R = L)처럼 왕복 행렬이 −I 라 자기 일치 q 가 정해지지 않을 때 쓴다."""
    from scipy.optimize import brentq

    def mismatch(log_z0):
        z0 = np.exp(log_z0)
        _, Rz = W_R_from_q(L / 2 - 1j * z0, lam)
        return Rz - R

    z0 = np.exp(brentq(mismatch, np.log(1e-6 * L), np.log(np.sqrt(L * R))))
    return np.sqrt(lam * z0 / np.pi), z0


def is_stable_rays(L, R1, R2, n_trips=400):
    """광선 하나를 왕복시켜 유계로 남는지 직접 본다 (행렬식 1 인 왕복 행렬을 거듭 곱한다)."""
    M = resonator_round_trip(L, R1, R2)
    v = np.array([1.0, 1e-3])
    for _ in range(n_trips):
        v = M @ v
        if np.abs(v[0]) > 1e6:
            return False
    return True


# ---------------------------------------------------------------- (나) 격자 위의 장

def grid(n, dx):
    x = (np.arange(n) - n // 2) * dx
    return x, np.meshgrid(x, x)


def gaussian_field(X, Y, W0, lam, z=0.0):
    """허리에서 거리 z 인 면의 가우시안 장 (q 파라미터 꼴, exp(ikz) 는 빼고 둔다)."""
    k = 2 * np.pi / lam
    q0 = q_at_waist(W0, lam)
    q = z + q0
    return (q0 / q) * np.exp(1j * k * (X**2 + Y**2) / (2 * q))


def angular_spectrum(U, dx, lam, z, remove_carrier=True):
    """각스펙트럼 전파 (3편). remove_carrier 면 평면파 위상 exp(ikz) 를 빼서 돌려준다."""
    U = np.asarray(U, complex)
    if U.ndim == 1:
        f2 = np.fft.fftfreq(U.size, dx) ** 2
    else:
        f2 = np.fft.fftfreq(U.shape[0], dx)[:, None] ** 2 + np.fft.fftfreq(U.shape[1], dx)[None, :] ** 2
    kz = 2 * np.pi * np.sqrt((1 / lam**2 - f2).astype(complex))
    if remove_carrier:
        kz = kz - 2 * np.pi / lam
    # 뒤로 전파할 때(z < 0)도 소멸파는 줄어들어야 한다: 허수부는 |z| 로 곱한다
    H = np.exp(1j * kz.real * z - np.abs(kz.imag) * abs(z))
    return np.fft.ifftn(np.fft.fftn(U) * H)


def thin_lens_phase(X, Y, f, lam):
    return np.exp(-1j * np.pi * (X**2 + Y**2) / (lam * f))


def second_moment_width(I, x):
    """x 방향 2차 모멘트 폭 W = 2σ (가우시안이면 1/e² 반지름과 같다). I 는 2차원 세기."""
    Ix = I.sum(axis=0)
    m = np.sum(x * Ix) / Ix.sum()
    return 2 * np.sqrt(np.sum((x - m) ** 2 * Ix) / Ix.sum())


def encircled_power(I, X, Y, r):
    return I[np.hypot(X, Y) <= r].sum() / I.sum()


# ---------------------------------------------------------------- 고차 모드

def hermite_gaussian(X, Y, W0, l, m):
    """HG_lm 모드의 허리 장: H_l(√2 x/W0) H_m(√2 y/W0) exp(−r²/W0²). 전파는 파동 경로로 따로 한다."""
    u = np.sqrt(2) * X / W0
    v = np.sqrt(2) * Y / W0
    return eval_hermite(l, u) * eval_hermite(m, v) * np.exp(-(X**2 + Y**2) / W0**2)


def laguerre_gaussian(X, Y, W0, p, l):
    """LG_p^l 모드의 허리 장: (√2 r/W0)^|l| L_p^|l|(2r²/W0²) exp(−r²/W0²) exp(i l φ)."""
    r2 = X**2 + Y**2
    rho = np.sqrt(2 * r2) / W0
    phi = np.arctan2(Y, X)
    return rho ** abs(l) * eval_genlaguerre(p, abs(l), 2 * r2 / W0**2) * np.exp(-r2 / W0**2) * np.exp(1j * l * phi)


# ---------------------------------------------------------------- 폭스-리 반복

def fox_li(n, dx, lam, L, R, a, n_iter, U0=None, record=()):
    """대칭 공진기(거울 곡률 R, 거리 L, 원형 구멍 반지름 a)에서 장을 왕복시킨다 (2차원).
    한 번 = 거리 L 진행 + 거울 반사(위상 −k r²/R, 구멍 밖은 잘림). 매번 세기를 1 로 맞추고 잃은 비율을 기록한다.
    record 에 적은 횟수의 장을 함께 돌려준다."""
    x, (X, Y) = grid(n, dx)
    k = 2 * np.pi / lam
    f2 = np.fft.fftfreq(n, dx)[:, None] ** 2 + np.fft.fftfreq(n, dx)[None, :] ** 2
    H = np.exp(1j * (2 * np.pi * np.sqrt((1 / lam**2 - f2).astype(complex)) - k) * L)
    mirror = np.exp(-1j * k * (X**2 + Y**2) / R) * (np.hypot(X, Y) <= a)
    U = (np.hypot(X, Y) <= a).astype(complex) if U0 is None else U0.astype(complex)
    U = U / np.sqrt(np.sum(np.abs(U) ** 2))
    losses, snaps = [], {}
    for it in range(1, n_iter + 1):
        U = np.fft.ifft2(np.fft.fft2(U) * H) * mirror
        p = np.sum(np.abs(U) ** 2)
        losses.append(1 - p)
        U = U / np.sqrt(p)
        if it in record:
            snaps[it] = U.copy()
    return U, np.array(losses), snaps

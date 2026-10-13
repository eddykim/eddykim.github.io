"""파동광학 배경이론 6편의 계산 모듈: 고NA 무수차(aplanatic) 렌즈의 벡터 집속.

길이는 파장 λ 를 1 로 둔 단위다(매질 굴절률 n 이면 매질 속 파장 λ/n 을 쓴다). 초점은 원점이다.

두 갈래의 계산을 둔다.
  (가) 리처즈-울프 적분: 방위각 적분을 베셀 함수로 미리 해 둔 1차원 적분 I0, I1, I2 (Zhan 2009 의 식 5.21–5.23, 5.30–5.33)
  (나) 벡터 평면파의 직접 합: 기준 구면 위 각 점 (θ, φ) 에서 나온 평면파의 전기장 벡터를 그대로 더한다
       (방위각 적분을 해석적으로 하지 않는다)
verify_vector_focus.py 는 두 갈래가 같은 장을 주는지, 낮은 NA 에서 스칼라 에어리 무늬로 가는지, 문헌 수치와 맞는지 본다.

규약: 입력 빔이 렌즈를 지나면 기준 구면 위에서
  - 동공의 지름 방향 성분은 θ̂ = (cosθ cosφ, cosθ sinφ, −sinθ) 방향으로 기울고,
  - 원주 방향 성분은 φ̂ = (−sinφ, cosφ, 0) 그대로 남으며,
  - 진폭에 √cosθ 가 곱해진다 (사인 조건 + 에너지 보존).
초점의 장은 그 평면파들의 합 E(r) ∝ ∬ E∞(θ, φ) exp(i k·r) sinθ dθ dφ 이다.
"""
import numpy as np
from scipy.special import jv


# ---------------------------------------------------------------- (가) 리처즈-울프 적분

def _theta_grid(NA, n, NA_inner=0.0, N=1201):
    a_out = np.arcsin(NA / n)
    a_in = np.arcsin(NA_inner / n) if NA_inner > 0 else 0.0
    return np.linspace(a_in, a_out, N)


def _I(order, weight, rho, z, th, k):
    rho = np.atleast_1d(np.asarray(rho, float))
    out = np.empty(rho.shape, complex)
    ez = np.exp(1j * k * z * np.cos(th))
    w = weight(th) * ez
    for idx, r in np.ndenumerate(rho):
        out[idx] = np.trapezoid(w * jv(order, k * r * np.sin(th)), th)
    return out


def rw_linear(rho, phi, z, NA, n=1.0, NA_inner=0.0, N=1201):
    """x 방향 선편광 입력의 초점 장 (Ex, Ey, Ez). rho, phi 는 같은 모양의 배열."""
    k = 2 * np.pi * n
    th = _theta_grid(NA, n, NA_inner, N)
    a = lambda t: np.sqrt(np.cos(t)) * np.sin(t)                          # noqa: E731
    I0 = _I(0, lambda t: a(t) * (1 + np.cos(t)), rho, z, th, k)
    I1 = _I(1, lambda t: a(t) * np.sin(t), rho, z, th, k)
    I2 = _I(2, lambda t: a(t) * (1 - np.cos(t)), rho, z, th, k)
    phi = np.asarray(phi, float)
    Ex = -1j * (I0 + I2 * np.cos(2 * phi))
    Ey = -1j * I2 * np.sin(2 * phi)
    Ez = -2 * I1 * np.cos(phi)
    return Ex, Ey, Ez


def rw_radial(rho, z, NA, n=1.0, NA_inner=0.0, N=1201):
    """방사 편광 입력의 초점 장 (Eρ, Ez). 방위각에 무관하다."""
    k = 2 * np.pi * n
    th = _theta_grid(NA, n, NA_inner, N)
    Er = _I(1, lambda t: np.sqrt(np.cos(t)) * np.sin(2 * t), rho, z, th, k)
    Ez = 2j * _I(0, lambda t: np.sqrt(np.cos(t)) * np.sin(t) ** 2, rho, z, th, k)
    return Er, Ez


def rw_azimuthal(rho, z, NA, n=1.0, NA_inner=0.0, N=1201):
    """방위 편광 입력의 초점 장 Eφ (다른 성분은 0)."""
    k = 2 * np.pi * n
    th = _theta_grid(NA, n, NA_inner, N)
    return 2 * _I(1, lambda t: np.sqrt(np.cos(t)) * np.sin(t), rho, z, th, k)


def rw_intensity_map(kind, x, NA, n=1.0, NA_inner=0.0, z=0.0):
    """초점면(z) 위 정사각 격자 x × x 의 성분별 세기. kind = 'x', 'circular', 'radial', 'azimuthal'.
    반지름 방향 1차원 적분을 한 번 하고 보간해 2차원 지도로 편다."""
    X, Y = np.meshgrid(x, x)
    R = np.hypot(X, Y)
    P = np.arctan2(Y, X)
    rs = np.linspace(0, R.max() * 1.001, 700)
    k = 2 * np.pi * n
    th = _theta_grid(NA, n, NA_inner)
    g = lambda arr: np.interp(R, rs, arr.real) + 1j * np.interp(R, rs, arr.imag)    # noqa: E731
    if kind in ("x", "circular"):
        a = lambda t: np.sqrt(np.cos(t)) * np.sin(t)                      # noqa: E731
        I0 = g(_I(0, lambda t: a(t) * (1 + np.cos(t)), rs, z, th, k))
        I1 = g(_I(1, lambda t: a(t) * np.sin(t), rs, z, th, k))
        I2 = g(_I(2, lambda t: a(t) * (1 - np.cos(t)), rs, z, th, k))
        def lin(phi0):                                                    # 편광 방향을 φ0 만큼 돌린 선편광
            ph = P - phi0
            ex_l, ey_l, ez = -1j * (I0 + I2 * np.cos(2 * ph)), -1j * I2 * np.sin(2 * ph), -2 * I1 * np.cos(ph)
            c, s = np.cos(phi0), np.sin(phi0)
            return c * ex_l - s * ey_l, s * ex_l + c * ey_l, ez
        if kind == "x":
            Ex, Ey, Ez = lin(0.0)
        else:                                                             # (x̂ + i ŷ)/√2
            ax, ay = lin(0.0), lin(np.pi / 2)
            Ex, Ey, Ez = [(u + 1j * v) / np.sqrt(2) for u, v in zip(ax, ay)]
    elif kind == "radial":
        Er, Ez = rw_radial(rs, z, NA, n, NA_inner)
        Er, Ez = g(Er), g(Ez)
        Ex, Ey = Er * np.cos(P), Er * np.sin(P)
    elif kind == "azimuthal":
        Ep = g(rw_azimuthal(rs, z, NA, n, NA_inner))
        Ex, Ey, Ez = -Ep * np.sin(P), Ep * np.cos(P), np.zeros_like(Ep)
    else:
        raise ValueError(kind)
    return np.abs(Ex) ** 2, np.abs(Ey) ** 2, np.abs(Ez) ** 2


# ---------------------------------------------------------------- (나) 벡터 평면파의 직접 합

def pupil_field_after_lens(theta, phi, kind="x"):
    """기준 구면 위의 장 E∞ (데카르트 성분 3개). 입력 편광을 지름·원주 성분으로 나눠 각각 θ̂, φ̂ 로 옮긴다."""
    c, s = np.cos(phi), np.sin(phi)
    if kind == "x":
        e_rho, e_phi = c, -s
    elif kind == "y":
        e_rho, e_phi = s, c
    elif kind == "radial":
        e_rho, e_phi = np.ones_like(phi), np.zeros_like(phi)
    elif kind == "azimuthal":
        e_rho, e_phi = np.zeros_like(phi), np.ones_like(phi)
    else:
        raise ValueError(kind)
    ct, st = np.cos(theta), np.sin(theta)
    th_hat = np.stack([ct * c, ct * s, -st])
    ph_hat = np.stack([-s, c, np.zeros_like(c)])
    return np.sqrt(ct) * (e_rho * th_hat + e_phi * ph_hat)


def plane_wave_sum(points, NA, kind="x", n=1.0, NA_inner=0.0, n_theta=300, n_phi=360):
    """초점 근처 점들 points[(x, y, z)] 의 장을 평면파 직접 합으로 구한다 (방위각 적분도 수치로)."""
    k = 2 * np.pi * n
    a_out = np.arcsin(NA / n)
    a_in = np.arcsin(NA_inner / n) if NA_inner > 0 else 0.0
    dth = (a_out - a_in) / n_theta
    th = a_in + (np.arange(n_theta) + 0.5) * dth                          # 중점 규칙
    ph = (np.arange(n_phi) + 0.5) * 2 * np.pi / n_phi
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    E = pupil_field_after_lens(TH, PH, kind)                              # (3, nθ, nφ)
    kx, ky, kz = k * np.sin(TH) * np.cos(PH), k * np.sin(TH) * np.sin(PH), k * np.cos(TH)
    w = np.sin(TH) * dth * (2 * np.pi / n_phi)
    out = np.empty((len(points), 3), complex)
    for i, (x, y, z) in enumerate(points):
        ph_ = np.exp(1j * (kx * x + ky * y + kz * z)) * w
        out[i] = [(E[j] * ph_).sum() for j in range(3)]
    return out


# ---------------------------------------------------------------- 두 평면파의 간섭 (직관)

def two_beam_fringe_visibility(half_angle, pol, n_x=2001):
    """세기가 같은 두 평면파가 z 축에 대해 ±θ 로 엇갈릴 때 x 방향 간섭무늬의 가시도.
    s 편광: 두 전기장이 모두 y. p 편광: 각 파의 진행 방향에 수직인 xz 평면 안의 벡터. 벡터 합으로 직접 계산한다."""
    th = half_angle
    x = np.linspace(0, 1 / max(np.sin(th), 1e-6), n_x)                    # 무늬 한 주기 남짓 (λ = 1)
    k = 2 * np.pi
    p1 = np.exp(1j * k * np.sin(th) * x)
    p2 = np.exp(-1j * k * np.sin(th) * x)
    if pol == "s":
        e1 = e2 = np.array([0.0, 1.0, 0.0])
    else:
        e1 = np.array([np.cos(th), 0.0, -np.sin(th)])                     # k1 = (sinθ, 0, cosθ) 에 수직
        e2 = np.array([np.cos(th), 0.0, np.sin(th)])                      # k2 = (−sinθ, 0, cosθ) 에 수직
    E = e1[:, None] * p1[None, :] + e2[:, None] * p2[None, :]
    I = np.sum(np.abs(E) ** 2, axis=0)
    return (I.max() - I.min()) / (I.max() + I.min())


# ---------------------------------------------------------------- 공기-유리 경계와 몰타 십자

def fresnel_t(theta_i, n1, n2):
    ct = np.cos(theta_i)
    st2 = (n1 / n2 * np.sin(theta_i)) ** 2
    ctt = np.sqrt(1 - st2 + 0j)
    ts = 2 * n1 * ct / (n1 * ct + n2 * ctt)
    tp = 2 * n1 * ct / (n2 * ct + n1 * ctt)
    return tp, ts


def crossed_polarizer_leakage(theta, phi, n_glass=1.5, n_surfaces=1):
    """x 편광으로 들어온 광선이 입사각 θ, 방위 φ 로 공기→유리 경계를 n_surfaces 번 지난 뒤
    y 검광자를 지나는 세기 비율. 국소 s·p 기저에서 투과계수를 곱해 직접 계산한다."""
    tp, ts = fresnel_t(theta, 1.0, n_glass)
    tp, ts = tp**n_surfaces, ts**n_surfaces
    c, s = np.cos(phi), np.sin(phi)
    Ep, Es = c, -s                                                        # x̂ 를 p(지름)·s(원주) 로
    Ex = tp * Ep * c - ts * Es * s                                        # 경계를 지난 뒤 다시 x·y 로
    Ey = tp * Ep * s + ts * Es * c
    return np.abs(Ey) ** 2, np.abs(Ex) ** 2 + np.abs(Ey) ** 2


# ---------------------------------------------------------------- 반점의 크기

def fwhm_1d(r, I):
    """r = 0 에서 시작하는 단면의 반치 전폭 (중심이 최대라고 본다)."""
    I = I / I.max()
    i = np.argmax(I < 0.5)
    return 2 * np.interp(0.5, [I[i], I[i - 1]], [r[i], r[i - 1]])


def half_max_area(I, dx):
    return (I > I.max() / 2).sum() * dx * dx

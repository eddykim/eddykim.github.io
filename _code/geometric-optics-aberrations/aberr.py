"""기하광학 배경이론 4편의 계산 모듈: 실광선 추적, 파면수차(OPD), 자이델 합, 제르니케, 회절 PSF.

광학계는 광축(z) 위에 놓인 구면들의 목록이다. 면마다 꼭짓점 위치 z, 곡률반지름 R(평면은 inf),
그 면 뒤 매질의 굴절률 n(λ), 조리개 여부를 갖는다. 길이 단위는 mm, 파장은 µm 이다.

부호 규약
- 근축 광선은 (높이 y, 기울기 u = dy/dz) 로 적고, 굴절은 n'u' = nu − y c (n' − n) 이다.
- 파면수차 W 는 주광선이 기준 구면에 닿을 때까지 쌓은 광로길이에서 그 광선의 광로길이를 뺀 값이다.
  실제 파면이 기준 구면보다 상 쪽으로 앞서 있으면 W > 0 (Hopkins·Welford 규약). 근축보다 세게 꺾이는
  (덜 보정된) 양렌즈에서 W040 > 0 이다. Born & Wolf 5.1 의 Φ 는 이것과 부호가 반대다(Φ = −W).
- 자이델 합 S_I…S_V 는 Welford 의 형태(면마다 굴절 불변량 A = n i 와 주광선의 Ā)를 쓴다.

문헌의 닫힌 식(얇은 렌즈의 수차, B&W 허용치 표, 페츠발 3:1 등)은 여기에 쓰지 않는다.
verify_aberr.py 가 그것들과 대조할 독립 경로여야 하기 때문이다.
"""
from dataclasses import dataclass
from math import factorial
from typing import Callable

import numpy as np

# ---------------------------------------------------------------- 유리 (Schott 셀마이어 계수)

N_BK7 = ((1.03961212, 0.231792344, 1.01046945), (0.00600069867, 0.0200179144, 103.560653))
N_F2 = ((1.39757037, 0.159201403, 1.2686543), (0.00995906143, 0.0546931752, 119.248346))
N_SK16 = ((1.34317774, 0.241144399, 0.994317969), (0.00704687339, 0.0229005, 92.7508526))

LINE_D, LINE_F, LINE_C = 0.5875618, 0.4861327, 0.6562816   # 프라운호퍼 d·F·C 선 (µm)


def sellmeier(wl_um, B, C):
    w2 = np.asarray(wl_um, dtype=float) ** 2
    return np.sqrt(1.0 + sum(b * w2 / (w2 - c) for b, c in zip(B, C)))


def glass(coeffs) -> Callable:
    return lambda wl: float(sellmeier(wl, *coeffs))


def AIR(wl):
    return 1.0


def abbe_number(medium):
    nd, nf, nc = medium(LINE_D), medium(LINE_F), medium(LINE_C)
    return (nd - 1.0) / (nf - nc)


# ---------------------------------------------------------------- 광학계

@dataclass
class Surface:
    z: float                    # 꼭짓점 위치
    R: float                    # 곡률반지름 (양수면 곡률 중심이 오른쪽)
    medium: Callable            # 이 면 뒤 매질의 굴절률 n(λ)
    semi: float = np.inf        # 구경 반지름 (조리개면이면 조리개 반지름)
    stop: bool = False

    @property
    def c(self):
        return 0.0 if np.isinf(self.R) else 1.0 / self.R


def stop_index(surfs):
    return next(i for i, s in enumerate(surfs) if s.stop)


# ---------------------------------------------------------------- 근축 추적

def paraxial_trace(surfs, y0, u0, z0, wl, n0=1.0):
    """평면 z0 에서 (y0, u0) 로 출발한 근축 광선. 면마다 [y, u(입사), u'(굴절 후), n, n', c] 한 행."""
    rows = []
    y, u, z, n = float(y0), float(u0), float(z0), float(n0)
    for s in surfs:
        y = y + (s.z - z) * u
        n1 = s.medium(wl)
        u1 = (n * u - y * s.c * (n1 - n)) / n1
        rows.append((y, u, u1, n, n1, s.c))
        u, n, z = u1, n1, s.z
    return np.array(rows)


@dataclass
class Setup:
    """주변광선과 주광선, 입사동·출사동, 근축 상면. 무한 물체면 field 는 반화각(rad)."""
    z0: float                   # 두 광선이 출발하는 평면
    marg: tuple                 # (y, u) at z0
    chief: tuple
    z_ep: float
    r_ep: float
    z_xp: float
    z_img: float                # 근축 상면 (축상 물점의 상)
    n_img: float
    sin_u_img: float            # 상 쪽 주변광선 기울기의 크기 (근축)
    infinite: bool
    field: float


def setup(surfs, wl, field, z_obj=None, n0=1.0):
    """조리개면 반지름으로 주변광선을, 조리개 중심을 지나도록 주광선을 정한다."""
    i = stop_index(surfs)
    z0 = surfs[0].z if z_obj is None else z_obj
    head = surfs[: i + 1]
    a = paraxial_trace(head, 1.0, 0.0, z0, wl, n0)[-1, 0]   # 조리개면 높이 = a·y0 + b·u0
    b = paraxial_trace(head, 0.0, 1.0, z0, wl, n0)[-1, 0]
    r = surfs[i].semi
    if z_obj is None:
        t = np.tan(field)
        marg, chief = (r / a, 0.0), (-b * t / a, t)
    else:
        marg, chief = (0.0, r / b), (field, -a * field / b)
    z_ep, r_ep = z0 + b / a, r / abs(a)
    m = paraxial_trace(surfs, *marg, z0, wl, n0)
    # 출사동: 조리개 중심에서 나온 광선이 상 공간에서 광축을 지나는 곳
    tail = paraxial_trace(surfs[i + 1:], 0.0, 1.0, surfs[i].z, wl, surfs[i].medium(wl)) if i + 1 < len(surfs) else None
    if tail is None:
        z_xp = surfs[i].z
    else:
        z_xp = surfs[-1].z - tail[-1, 0] / tail[-1, 2]
    z_img = surfs[-1].z - m[-1, 0] / m[-1, 2]
    return Setup(z0, marg, chief, z_ep, r_ep, z_xp, z_img, m[-1, 4], abs(m[-1, 2]), z_obj is None, field)


def paraxial_image_height(surfs, wl, st, n0=1.0):
    c = paraxial_trace(surfs, *st.chief, st.z0, wl, n0)
    return c[-1, 0] + (st.z_img - surfs[-1].z) * c[-1, 2]


# ---------------------------------------------------------------- 자이델 합

def seidel_surface_terms(surfs, wl, st, n0=1.0):
    """면마다의 S_I…S_V (5 × 면 수). 주변광선과 주광선 두 개의 근축 데이터만 쓴다."""
    m = paraxial_trace(surfs, *st.marg, st.z0, wl, n0)
    c = paraxial_trace(surfs, *st.chief, st.z0, wl, n0)
    y, u, u1, n, n1, cc = m.T
    yb, ub = c[:, 0], c[:, 1]
    H = n0 * (st.chief[1] * st.marg[0] - st.marg[1] * st.chief[0])   # 라그랑주 불변량
    A, Ab = n * (u + y * cc), n * (ub + yb * cc)                    # 굴절 불변량 (n·입사각)
    d_un = u1 / n1 - u / n
    S1 = -A**2 * y * d_un
    S2 = -A * Ab * y * d_un
    S3 = -Ab**2 * y * d_un
    S4 = -H**2 * cc * (1.0 / n1 - 1.0 / n)
    # S_V = Σ (Ā/A)(S_III + S_IV) 이지만 A = 0 인 면(평행광이 평면에 닿는 경우 등)에서 0/0 이 된다.
    # Δ(u/n) = A·Δ(1/n²) − y c Δ(1/n) 와 Ā y − A ȳ = H 를 넣어 나눗셈을 없앤 꼴 (Welford):
    S5 = -Ab**3 * y * (1.0 / n1**2 - 1.0 / n**2) + Ab * yb * cc * (H + Ab * y) * (1.0 / n1 - 1.0 / n)
    return np.vstack([S1, S2, S3, S4, S5])


def seidel_sums(surfs, wl, st, n0=1.0):
    return seidel_surface_terms(surfs, wl, st, n0).sum(axis=1)


def wave_coefficients(S):
    """자이델 합 → 파면수차 계수 (같은 길이 단위). 정규화된 동공 ρ 와 시야 h 에 대해
    W = W040 ρ⁴ + W131 h ρ³ cosθ + W222 h² ρ² cos²θ + W220 h² ρ² + W311 h³ ρ cosθ."""
    S1, S2, S3, S4, S5 = S
    return dict(W040=S1 / 8, W131=S2 / 2, W222=S3 / 2, W220=(S3 + S4) / 4, W311=S5 / 2,
                W220P=S4 / 4)


# ---------------------------------------------------------------- 실광선 추적

def trace_real(surfs, P, D, wl, n0=1.0, opl0=None):
    """벡터 스넬 법칙으로 3차원 광선 묶음을 추적한다. P, D: (N, 3), D 는 단위벡터.
    반환: 마지막 면 위의 점, 방향, 마지막 매질 굴절률, 쌓인 광로길이, 살아남은 광선 표시."""
    P, D = np.array(P, float), np.array(D, float)
    opl = np.zeros(len(P)) if opl0 is None else np.array(opl0, float)
    ok = np.ones(len(P), bool)
    n = n0
    for s in surfs:
        if np.isinf(s.R):
            t = (s.z - P[:, 2]) / D[:, 2]
            Q = P + t[:, None] * D
            N = np.tile([0.0, 0.0, 1.0], (len(P), 1))
        else:
            C = np.array([0.0, 0.0, s.z + s.R])
            oc = P - C
            b = np.sum(oc * D, 1)
            disc = b * b - (np.sum(oc * oc, 1) - s.R**2)
            ok &= disc >= 0
            root = np.sqrt(np.maximum(disc, 0.0))
            t = -b - root if s.R > 0 else -b + root      # 꼭짓점 쪽 교점
            Q = P + t[:, None] * D
            N = (Q - C) / abs(s.R)
            N *= np.sign(N[:, 2:3])                       # 법선을 +z 쪽으로
        ok &= Q[:, 0] ** 2 + Q[:, 1] ** 2 <= s.semi**2 * (1 + 1e-12)
        opl += n * t
        n1 = s.medium(wl)
        mu = n / n1
        cosi = np.sum(D * N, 1)
        k = 1.0 - mu**2 * (1.0 - cosi**2)
        ok &= k >= 0
        D = mu * D + (np.sqrt(np.maximum(k, 0.0)) - mu * cosi)[:, None] * N
        P, n = Q, n1
    return P, D, n, opl, ok


def to_plane(P, D, z):
    t = (z - P[:, 2]) / D[:, 2]
    return P + t[:, None] * D


def launch(st, px, py, n0=1.0, back=20.0):
    """정규화 동공 좌표 (px, py) 의 광선들을 근축 입사동을 겨눠 쏜다. 시야는 +y 쪽.
    무한 물체: 기울어진 평면파. 광로길이의 시작값을 평면파의 위상면에 맞춘다."""
    px, py = np.ravel(px), np.ravel(py)
    E = np.stack([px * st.r_ep, py * st.r_ep, np.full(px.shape, st.z_ep)], 1)   # 입사동 위의 점
    if st.infinite:
        d = np.array([0.0, np.sin(st.field), np.cos(st.field)])
        z_start = min(st.z0, st.z_ep) - back
        t = (z_start - E[:, 2]) / d[2]
        P = E + t[:, None] * d
        D = np.tile(d, (len(px), 1))
        opl0 = n0 * np.sum(P * d, 1)              # 같은 평면파 위의 점은 같은 위상
    else:
        O = np.array([0.0, st.field, st.z0])
        D = E - O
        D /= np.linalg.norm(D, axis=1)[:, None]
        P = np.tile(O, (len(px), 1))
        opl0 = np.zeros(len(px))
    return P, D, opl0


def chief_image_point(surfs, wl, st, z_image, n0=1.0):
    P, D, opl0 = launch(st, np.array([0.0]), np.array([0.0]), n0)
    P, D, n, opl, ok = trace_real(surfs, P, D, wl, n0, opl0)
    return to_plane(P, D, z_image)[0]


def opd(surfs, wl, st, px, py, z_image=None, ref_point=None, n0=1.0):
    """기준 구면에 대한 파면수차 W(px, py) (mm), 앞서면 양수. 기준 구면: 중심 = 기준점(기본: 주광선이 상면에 닿는 점),
    출사동 위 주광선 점을 지난다. 반환: W, 살아남은 광선 표시."""
    z_image = st.z_img if z_image is None else z_image
    if ref_point is None:
        ref_point = chief_image_point(surfs, wl, st, z_image, n0)
    Pc, Dc, oc = launch(st, np.array([0.0]), np.array([0.0]), n0)
    Pc, Dc, n_img, oplc, _ = trace_real(surfs, Pc, Dc, wl, n0, oc)
    E = to_plane(Pc, Dc, st.z_xp)[0]
    Rr = np.linalg.norm(ref_point - E)

    def to_sphere(P, D, opl):
        oc_ = P - ref_point
        b = np.sum(oc_ * D, 1)
        disc = b * b - (np.sum(oc_ * oc_, 1) - Rr**2)
        t = -b - np.sqrt(np.maximum(disc, 0.0))
        return opl + n_img * t

    P, D, opl0 = launch(st, px, py, n0)
    P, D, n_img, opl, ok = trace_real(surfs, P, D, wl, n0, opl0)
    W = to_sphere(Pc, Dc, oplc)[0] - to_sphere(P, D, opl)
    return W.reshape(np.shape(px)), ok.reshape(np.shape(px))


def spot(surfs, wl, st, px, py, z_image, n0=1.0):
    """상면 z_image 위 광선 도착점 (x, y) 와 살아남은 광선 표시."""
    P, D, opl0 = launch(st, px, py, n0)
    P, D, n, opl, ok = trace_real(surfs, P, D, wl, n0, opl0)
    Q = to_plane(P, D, z_image)
    return Q[:, 0], Q[:, 1], ok


def axis_crossing(surfs, wl, h, z_start=None, n0=1.0):
    """광축과 나란히 높이 h 로 들어온 자오 광선이 마지막 면을 지난 뒤 광축을 지나는 z."""
    h = np.atleast_1d(np.asarray(h, float))
    z_start = surfs[0].z - 10.0 if z_start is None else z_start
    P = np.stack([np.zeros_like(h), h, np.full_like(h, z_start)], 1)
    D = np.tile([0.0, 0.0, 1.0], (len(h), 1))
    P, D, n, opl, ok = trace_real(surfs, P, D, wl, n0)
    return P[:, 2] - P[:, 1] * D[:, 2] / D[:, 1]


# ---------------------------------------------------------------- 다항식 맞춤

def fit_seidel_terms(W, px, py, h=1.0, with_tilt_defocus=True):
    """W(ρ, θ) 를 3차 항(+ 피스톤·기울기·초점) 과 5차 일부에 최소제곱으로 맞춘다.
    θ 는 +y 축에서 잰다 (시야가 +y). 반환: 각 항의 계수 dict (h 로 나눈 값)."""
    x, y = np.ravel(px), np.ravel(py)
    W = np.ravel(W)
    r2 = x**2 + y**2
    basis = {"piston": np.ones_like(x), "tilt": y, "defocus": r2,
             "W040": r2**2, "W131": y * r2, "W222": y**2,
             "W060": r2**3, "W151": y * r2**2, "W242": y**2 * r2, "W331": y**3}
    if not with_tilt_defocus:
        basis.pop("tilt"); basis.pop("defocus")
    M = np.stack(list(basis.values()), 1)
    coef, *_ = np.linalg.lstsq(M, W, rcond=None)
    return dict(zip(basis.keys(), coef))


def zernike_radial(n, m, rho):
    m = abs(m)
    out = np.zeros_like(rho, dtype=float)
    for s in range((n - m) // 2 + 1):
        out += ((-1) ** s * factorial(n - s)
                / (factorial(s) * factorial((n + m) // 2 - s) * factorial((n - m) // 2 - s))) * rho ** (n - 2 * s)
    return out


def noll_to_nm(j):
    """Noll(1976) 번호 j → (n, m). m > 0 은 cos, m < 0 은 sin 항."""
    n, rest = 0, j - 1
    while rest > n:
        n += 1
        rest -= n
    m = (-1) ** j * ((n % 2) + 2 * ((rest + (n + 1) % 2) // 2))
    return n, m


def zernike_noll(j, rho, theta):
    """정규직교(원판 위 평균제곱 = 1) 제르니케 다항식. θ 는 x 축에서 잰다."""
    n, m = noll_to_nm(j)
    R = zernike_radial(n, m, rho)
    if m == 0:
        return np.sqrt(n + 1) * R
    norm = np.sqrt(2 * (n + 1))
    return norm * R * (np.cos(m * theta) if m > 0 else np.sin(-m * theta))


# ---------------------------------------------------------------- 회절 (FFT)

def pupil_grid(n=255):
    s = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(s, s)
    return X, Y, X**2 + Y**2 <= 1.0


def psf(W_waves, mask, pad=4):
    """동공 위 파면수차(파장 단위) → 세기 PSF. 수차 없는 같은 동공의 정점이 1 이 되게 나눈다.
    화소 간격 = (λ / NA) · (n − 1) / (2 · N_fft).
    W > 0 은 파면이 앞선 것이므로 그 광선의 위상은 exp(−i2πW) 다 (Goodman 6.4 의 W 와 부호가 반대).
    이 부호여야 회절상이 광선 도착점(ε = −∇W/NA)과 같은 쪽으로 옮겨 간다."""
    n = mask.shape[0]
    N = pad * n
    field = np.zeros((N, N), complex)
    field[:n, :n] = mask * np.exp(-2j * np.pi * np.where(mask, W_waves, 0.0))
    F = np.fft.fftshift(np.fft.fft2(field))
    return np.abs(F) ** 2 / mask.sum() ** 2


def psf_pixel(n, pad):
    """psf() 한 화소의 크기 (λ/NA 단위)."""
    return (n - 1) / (2.0 * pad * n)


def strehl_center(W_waves, mask):
    """기준점(동공 중심 방향)에서의 세기 비: |평균 exp(i2πW)|²."""
    return float(np.abs(np.mean(np.exp(-2j * np.pi * W_waves[mask]))) ** 2)


def rms(W, mask, remove=("piston",), px=None, py=None):
    """동공 위 RMS. remove 에 'tilt', 'defocus' 를 주면 그 항까지 빼고 잰다."""
    w = W[mask]
    cols = [np.ones_like(w)]
    if "tilt" in remove:
        cols += [px[mask], py[mask]]
    if "defocus" in remove:
        cols += [px[mask] ** 2 + py[mask] ** 2]
    M = np.stack(cols, 1)
    coef, *_ = np.linalg.lstsq(M, w, rcond=None)
    return float(np.sqrt(np.mean((w - (M * coef).sum(1)) ** 2)))


def geometric_spot_from_W(dWdx, dWdy, na, n_img=1.0):
    """파면 기울기 → 횡광선수차 (근축 관계). W 와 같은 길이 단위로 돌려준다.
    ε = −(1 / (n′ NA)) ∂W/∂ρ  (ρ: 정규화 동공 좌표). 부호는 verify 의 실광선 대조로 확인했다."""
    return -dWdx / (n_img * na), -dWdy / (n_img * na)

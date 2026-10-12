"""파동광학 배경이론 4편의 계산 모듈: 렌즈의 푸리에 변환, 4f 공간 필터, 결맞는·결맞지 않는·부분 결맞는 결상, 전달함수, 위상차.

길이는 µm 로 쓴다. 상은 배율 1 로 두고(물체 좌표 = 상 좌표), 결상계는 물체 쪽 개구수 NA 로만 나타낸다.
복소 진폭은 exp(+ikr) 로 퍼져 나가는 꼴이다.
문헌의 닫힌 식(원형 동공 OTF, 정사각 동공의 초점 어긋남 OTF, Si(π) 넘침, 0.82λ/NA, 1 ± 2φ)은
여기서 쓰지 않는다. verify_imaging.py 가 그것들과 대조할 독립 경로여야 하기 때문이다.
"""
import numpy as np


# ---------------------------------------------------------------- 격자와 자유 공간

def grid(n, dx):
    """n×n 표본, 간격 dx 의 좌표와 공간 주파수 (FFT 순서)."""
    x = (np.arange(n) - n // 2) * dx
    f = np.fft.fftfreq(n, dx)
    return x, f


def angular_spectrum(U, dx, lam, z):
    """평면파로 쪼개 kz z 만큼 위상을 더한 뒤 다시 합친다 (3편과 같은 방법). 1차원·2차원 모두."""
    U = np.asarray(U, complex)
    if U.ndim == 1:
        f2 = np.fft.fftfreq(U.size, dx) ** 2
    else:
        f2 = np.fft.fftfreq(U.shape[0], dx)[:, None] ** 2 + np.fft.fftfreq(U.shape[1], dx)[None, :] ** 2
    kz = 2 * np.pi * np.sqrt((1 / lam**2 - f2).astype(complex))
    return np.fft.ifftn(np.fft.fftn(U) * np.exp(1j * kz * z))


def thin_lens(x, y, f, lam):
    """얇은 렌즈의 위상 변환 exp(−ik(x²+y²)/2f): 평면파를 거리 f 에 모이는 구면파로 바꾼다."""
    return np.exp(-1j * np.pi * (x**2 + y**2) / (lam * f))


def lens_fourier_plane(U, dx, lam, f):
    """물체를 렌즈 앞 f 에 두고, 렌즈를 지나 다시 f 를 진행한 장 (각스펙트럼 전파로 직접 계산)."""
    n = U.shape[0]
    x = (np.arange(n) - n // 2) * dx
    X, Y = np.meshgrid(x, x)
    U1 = angular_spectrum(U, dx, lam, f)
    return angular_spectrum(U1 * thin_lens(X, Y, f, lam), dx, lam, f)


# ---------------------------------------------------------------- 동공과 전달함수

def pupil(fx, fy, NA, lam, W=None, center=(0.0, 0.0), shape="circle"):
    """공간 주파수 평면 위의 동공: 반지름 NA/λ 의 원 (shape="square" 면 반폭 NA/λ 의 정사각형).
    W(ρ, θ) 는 파장 단위 파면수차(ρ 는 동공 반지름으로 정규화). center 는 동공 중심을 옮길 때."""
    u, v = (fx - center[0]) * lam / NA, (fy - center[1]) * lam / NA
    rho = np.hypot(u, v)
    inside = (rho <= 1.0) if shape == "circle" else ((np.abs(u) <= 1.0) & (np.abs(v) <= 1.0))
    P = inside.astype(complex)
    if W is not None:
        th = np.arctan2(fy - center[1], fx - center[0])
        P = P * np.exp(1j * 2 * np.pi * W(rho, th))
    return P


def pupil_grid(n, dx, NA, lam, W=None, shape="circle"):
    _, f = grid(n, dx)
    FX, FY = np.meshgrid(f, f)
    return pupil(FX, FY, NA, lam, W, shape=shape)


def coherent_image(U_obj, dx, NA, lam, W=None, illum=(0.0, 0.0)):
    """결맞는 결상: 물체 장의 스펙트럼에 동공을 곱한다 (진폭에 대한 필터). illum 은 조명 평면파의 공간 주파수.
    비스듬한 조명은 물체 스펙트럼을 illum 만큼 옮기는 것과 같으므로, 거꾸로 동공을 −illum 만큼 옮겨 곱한다.
    이렇게 하면 illum 이 FFT 격자 위에 있지 않아도 스펙트럼이 새지 않는다."""
    n = U_obj.shape[0]
    x, f = grid(n, dx)
    FX, FY = np.meshgrid(f, f)
    P = pupil(FX + illum[0], FY + illum[1], NA, lam, W)
    U = np.fft.ifft2(np.fft.fft2(U_obj) * P)
    if illum[0] or illum[1]:
        X, Y = np.meshgrid(x, x)
        U = U * np.exp(1j * 2 * np.pi * (illum[0] * X + illum[1] * Y))
    return U


def amplitude_psf(n, dx, NA, lam, W=None):
    """진폭 점퍼짐함수: 동공의 역 푸리에 변환 (중심이 (n/2, n/2))."""
    return np.fft.fftshift(np.fft.ifft2(pupil_grid(n, dx, NA, lam, W)))


def otf_from_psf(n, dx, NA, lam, W=None, shape="circle"):
    """경로 1: 세기 점퍼짐함수 |h|² 의 푸리에 변환을 0 주파수에서 1 로 정규화 (FFT 순서)."""
    h = np.fft.ifft2(pupil_grid(n, dx, NA, lam, W, shape))
    O = np.fft.fft2(np.abs(h) ** 2)
    return O / O[0, 0]


def otf_by_overlap(f_shift, NA, lam, W=None, m=801, shape="circle"):
    """경로 2: 동공과, f_shift 만큼 옮긴 동공의 켤레를 겹쳐 더한다 (자기상관을 직접 계산). x 방향 이동만."""
    c = NA / lam
    g = np.linspace(-c, c, m)
    FX, FY = np.meshgrid(g, g)
    P0 = pupil(FX, FY, NA, lam, W, shape=shape)
    norm = np.sum(np.abs(P0) ** 2)
    out = []
    for s in np.atleast_1d(f_shift):
        P1 = pupil(FX + s / 2, FY, NA, lam, W, shape=shape)        # P(f + s/2)
        P2 = pupil(FX - s / 2, FY, NA, lam, W, shape=shape)        # P(f − s/2)
        out.append(np.sum(P1 * np.conj(P2)) / norm)
    return np.array(out)


def incoherent_image(I_obj, dx, NA, lam, W=None):
    """결맞지 않는 결상: 물체 세기의 스펙트럼에 OTF 를 곱한다 (세기에 대한 필터)."""
    n = I_obj.shape[0]
    return np.real(np.fft.ifft2(np.fft.fft2(I_obj) * otf_from_psf(n, dx, NA, lam, W)))


def source_points(sigma, NA, lam, n_ring=12):
    """콘덴서 개구(반지름 σ·NA/λ)를 고르게 채운 광원 점들의 공간 주파수와 무게.
    고리 n_ring 개의 가운데 반지름 (i − ½)/n_ring 에 점을 놓는 중점 규칙이다. 개구 가장자리에 정확히 놓인 점이 없어
    동공 가장자리와 겹치는 경계 효과를 피한다. 무게는 고리 넓이에 비례."""
    if sigma <= 0:
        return np.array([[0.0, 0.0]]), np.array([1.0])
    R = sigma * NA / lam
    pts, w = [], []
    for i in range(1, n_ring + 1):
        r = R * (i - 0.5) / n_ring
        k = 6 * i - 3
        area = np.pi * R**2 * ((i / n_ring) ** 2 - ((i - 1) / n_ring) ** 2)
        for j in range(k):
            a = 2 * np.pi * (j + 0.5 * (i % 2)) / k
            pts.append((r * np.cos(a), r * np.sin(a)))
            w.append(area / k)
    w = np.array(w)
    return np.array(pts), w / w.sum()


def partially_coherent_image(U_obj, dx, NA, lam, sigma, W=None, n_ring=12):
    """부분 결맞음 결상 (아베의 방법): 콘덴서 개구의 각 점이 하나의 평면파 조명이다.
    점마다 결맞는 상의 세기를 구하고, 서로 결맞지 않으므로 세기로 더한다."""
    pts, w = source_points(sigma, NA, lam, n_ring)
    I = np.zeros(U_obj.shape)
    for (sx, sy), ww in zip(pts, w):
        I += ww * np.abs(coherent_image(U_obj, dx, NA, lam, W, illum=(sx, sy))) ** 2
    return I


# ---------------------------------------------------------------- 4f 공간 필터 (아베-포터)

def filter_4f(U_obj, mask):
    """4f 계: 앞 렌즈가 물체를 푸리에 변환하고, 그 면에 마스크를 두고, 뒤 렌즈가 다시 변환한다.
    상이 뒤집히는 것은 좌표를 되돌려 생략했다. mask 는 FFT 순서의 공간 주파수 배열."""
    return np.fft.ifft2(np.fft.fft2(U_obj) * mask)


def phase_contrast(U_obj, dx, NA, lam, dot_radius, a=1.0, alpha=-np.pi / 2):
    """제르니케 위상차: 동공 가운데 반지름 dot_radius(공간 주파수) 안의 빛(0차와 아주 낮은 주파수)에
    진폭 a, 위상 α 를 곱한다. exp(+ikr) 규약에서 e^{+iα} 는 빛을 늦춘다 (유리를 더 지난 것과 같다).
    α = +π/2: 0차를 1/4 파장 늦춤 → 두꺼운 곳이 밝다 (음성 위상차, 세기 ≈ 1 + 2φ).
    α = −π/2: 0차를 1/4 파장 앞당김 → 두꺼운 곳이 어둡다 (양성 위상차, 세기 ≈ 1 − 2φ)."""
    n = U_obj.shape[0]
    _, f = grid(n, dx)
    FX, FY = np.meshgrid(f, f)
    P = pupil(FX, FY, NA, lam)
    dot = np.hypot(FX, FY) <= dot_radius
    P = np.where(dot, P * a * np.exp(1j * alpha), P)
    return np.fft.ifft2(np.fft.fft2(U_obj) * P)


def dark_field(U_obj, dx, NA, lam, dot_radius):
    return phase_contrast(U_obj, dx, NA, lam, dot_radius, a=0.0, alpha=0.0)


# ---------------------------------------------------------------- 시험 물체

def modulation(I):
    return (I.max() - I.min()) / (I.max() + I.min())


def spoke_target(n, n_spokes=36):
    y, x = np.mgrid[:n, :n] - (n - 1) / 2
    th = np.arctan2(y, x)
    r = np.hypot(x, y)
    return ((np.cos(n_spokes * th) > 0) & (r < n / 2 - 2)).astype(float)

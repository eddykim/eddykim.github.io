"""파동광학 배경이론 3편의 계산 모듈: 호이겐스-프레넬 적분, 각스펙트럼 전파, 프레넬 존, 원거리 무늬, 초점 근처의 장, 격자.

길이는 µm 로 쓴다. 단색광의 복소 진폭은 exp(+ikr) 로 퍼져 나가는 꼴이다.
문헌의 닫힌 식(2J1(x)/x, 1.22λ/D, 프레넬 적분, sinc²(u/4), 축 위의 정확해, 격자 분해능 mN, 6.25 %·33.8 %)은
여기서 쓰지 않는다. verify_diffraction.py 가 그것들과 대조할 독립 경로여야 하기 때문이다.
"""
import numpy as np
from scipy.special import j0


# ---------------------------------------------------------------- 호이겐스-프레넬 (레일리-좀머펠트 제1해)

def rs_kernel(r, z, lam):
    """레일리-좀머펠트 제1해의 핵: 넓이 dS 의 구멍 조각이 거리 r 의 점에 보내는 기여 (z/r 가 기울기 인자).
    멀리서는 (1/iλ) cosθ e^{ikr}/r 로 줄어든다."""
    k = 2 * np.pi / lam
    return (z / r) * (1 / r - 1j * k) * np.exp(1j * k * r) / (2 * np.pi * r)


def rs_direct(xi, eta, U, dA, points, lam):
    """구멍 위의 표본점 (xi, eta) 에 놓인 장 U 를 관측점들 points[(x, y, z)] 로 직접 적분한다 (2차 파원의 합)."""
    out = np.empty(len(points), complex)
    for i, (x, y, z) in enumerate(points):
        r = np.sqrt((x - xi) ** 2 + (y - eta) ** 2 + z**2)
        out[i] = np.sum(U * rs_kernel(r, z, lam)) * dA
    return out


def axial_ring_sum(rho_edges, z, lam, exact=True):
    """축 위 한 점에서 본 원형 띠들의 기여 (평면파 조명, 진폭 1). rho_edges 사이의 띠마다 복소 기여를 돌려준다.
    띠 안을 잘게 나눠 핵을 적분한다. 누적하면 위상자 나선이 된다."""
    k = 2 * np.pi / lam
    out = np.empty(len(rho_edges) - 1, complex)
    for i in range(len(rho_edges) - 1):
        rho = np.linspace(rho_edges[i], rho_edges[i + 1], 401)
        r = np.sqrt(rho**2 + z**2)
        f = rs_kernel(r, z, lam) if exact else np.exp(1j * k * (z + rho**2 / (2 * z))) / (1j * lam * z)
        out[i] = np.trapezoid(f * 2 * np.pi * rho, rho)
    return out


def zone_edges(z, lam, m_max):
    """축 위 점까지의 광로가 반파장씩 늘어나는 띠(프레넬 존)의 경계 반지름 (근사 없이)."""
    m = np.arange(m_max + 1)
    return np.sqrt((z + m * lam / 2) ** 2 - z**2)


def axial_circular_aperture(a, z, lam, n=4001):
    """반지름 a 인 원형 구멍 뒤 축 위의 세기 (평면파 1 로 정규화). 구멍 전체를 반지름 방향으로 직접 적분."""
    rho = np.linspace(0, a, n)
    r = np.sqrt(rho**2 + z**2)
    return np.abs(np.trapezoid(rs_kernel(r, z, lam) * 2 * np.pi * rho, rho)) ** 2


# ---------------------------------------------------------------- 프레넬 반지름 적분 (축대칭)

def fresnel_hankel(rho, t_rho, r, z, lam):
    """축대칭 구멍 t(ρ) 를 평면파로 비춘 장을 거리 z 에서 프레넬 근사로 구한다:
    U(r) = (2π/iλz) e^{ikz} e^{iπr²/λz} ∫ t(ρ) e^{iπρ²/λz} J0(2πρr/λz) ρ dρ."""
    k = 2 * np.pi / lam
    r = np.atleast_1d(r)
    ph = np.exp(1j * np.pi * rho**2 / (lam * z)) * t_rho * rho
    out = np.array([np.trapezoid(ph * j0(2 * np.pi * rho * rr / (lam * z)), rho) for rr in r])
    return 2 * np.pi / (1j * lam * z) * np.exp(1j * k * z) * np.exp(1j * np.pi * r**2 / (lam * z)) * out


def disk_shadow(a, r, z, lam, n=20001):
    """반지름 a 인 불투명 원판 뒤의 장. 원판 = (막힘 없는 평면파) − (같은 크기의 구멍) 이므로 구멍 적분으로 구한다.
    (스칼라 회절 적분이 선형이라 장끼리는 정확히 성립한다.)"""
    k = 2 * np.pi / lam
    rho = np.linspace(0, a, n)
    return np.exp(1j * k * z) - fresnel_hankel(rho, np.ones_like(rho), r, z, lam)


def zone_plate_transmission(rho, f, lam, n_open):
    """홀수 번째 존만 열린 존 플레이트 (n_open 개). 경계는 근축 반지름 √(mλf)."""
    m = rho**2 / (lam * f)
    zone = np.floor(m).astype(int)          # 0 이 첫 존
    return ((zone % 2 == 0) & (zone < 2 * n_open)).astype(float)


def axial_zone_plate(f, lam, n_open, z, n=200001):
    rho = np.linspace(0, np.sqrt(2 * n_open * lam * f) * 1.01, n)
    t = zone_plate_transmission(rho, f, lam, n_open)
    return np.array([np.abs(fresnel_hankel(rho, t, 0.0, zz, lam)[0]) ** 2 for zz in np.atleast_1d(z)])


# ---------------------------------------------------------------- 각스펙트럼 전파

def angular_spectrum(U, dx, lam, z):
    """장을 평면파로 쪼개 각각 kz z 만큼 위상을 더한 뒤 다시 합친다 (1차원이면 x–z 평면, 2차원이면 x–y).
    공간 주파수가 1/λ 를 넘는 성분은 kz 가 허수라 지수적으로 사라진다 (소멸파)."""
    U = np.asarray(U, complex)
    if U.ndim == 1:
        f2 = np.fft.fftfreq(U.size, dx) ** 2
    else:
        fx = np.fft.fftfreq(U.shape[1], dx)
        fy = np.fft.fftfreq(U.shape[0], dx)
        f2 = fy[:, None] ** 2 + fx[None, :] ** 2
    kz = 2 * np.pi * np.sqrt((1 / lam**2 - f2).astype(complex))
    F = np.fft.fftn(U)
    return np.fft.ifftn(F * np.exp(1j * kz * z))


def fresnel_transfer(U, dx, lam, z):
    """근축(프레넬) 전달함수 exp(ikz) exp(−iπλz f²) 로 전파한다. 각스펙트럼의 제곱근을 이항 전개해 둘째 항까지 둔 것."""
    U = np.asarray(U, complex)
    if U.ndim == 1:
        f2 = np.fft.fftfreq(U.size, dx) ** 2
    else:
        f2 = np.fft.fftfreq(U.shape[0], dx)[:, None] ** 2 + np.fft.fftfreq(U.shape[1], dx)[None, :] ** 2
    return np.fft.ifftn(np.fft.fftn(U) * np.exp(1j * 2 * np.pi / lam * z) * np.exp(-1j * np.pi * lam * z * f2))


def propagate_stack(U, dx, lam, zs):
    """같은 입력을 여러 거리로 보낸 결과 (옆에서 본 지도용, 1차원 입력)."""
    U = np.asarray(U, complex)
    f2 = np.fft.fftfreq(U.size, dx) ** 2
    kz = 2 * np.pi * np.sqrt((1 / lam**2 - f2).astype(complex))
    F = np.fft.fft(U)
    return np.array([np.fft.ifft(F * np.exp(1j * kz * z)) for z in zs])


def transfer_function(f, lam, z):
    """자유 공간의 전달함수 H(f) = exp(i 2π z √(1/λ² − f²))."""
    return np.exp(1j * 2 * np.pi * z * np.sqrt((1 / lam**2 - np.asarray(f, float) ** 2).astype(complex)))


# ---------------------------------------------------------------- 원거리 (프라운호퍼) 무늬

def far_field(pupil, pad=8):
    """구멍 함수의 2차원 FFT 로 원거리 세기를 구한다. 중심이 가운데 오게 돌려준다.
    좌표: 출력 한 칸 = 1/(pad·N·dx) 의 공간 주파수 = 각도 λ/(pad·N·dx)."""
    n = pupil.shape[0]
    big = np.zeros((pad * n, pad * n), complex)
    big[:n, :n] = pupil
    F = np.fft.fftshift(np.fft.fft2(big))
    return np.abs(F) ** 2


def circle_pupil(n, radius_px, obstruction=0.0):
    y, x = np.mgrid[:n, :n] - (n - 1) / 2
    r = np.hypot(x, y)
    return ((r <= radius_px) & (r >= obstruction * radius_px)).astype(float)


def radial_profile(I, step=1):
    """중심이 (N/2, N/2) 인 세기 배열에서 반지름 방향 단면 (오른쪽 가로줄)."""
    c = I.shape[0] // 2
    return I[c, c::step]


def encircled_energy(I, radii_px):
    """세기 배열의 중심에서 반지름 radii 안에 든 에너지 비율."""
    c = I.shape[0] // 2
    y, x = np.mgrid[:I.shape[0], :I.shape[1]] - c
    r = np.hypot(x, y)
    total = I.sum()
    return np.array([I[r <= rr].sum() / total for rr in np.atleast_1d(radii_px)])


def hexagon_mask(x, y, cx, cy, flat_radius):
    """꼭짓점이 좌우를 향하는(윗변·아랫변이 수평인) 정육각형. flat_radius 는 중심에서 변까지 거리."""
    dx, dy = x - cx, y - cy
    ok = np.ones_like(x, bool)
    for ang in (90, 30, -30):
        nrm = np.radians(ang)
        ok &= np.abs(dx * np.cos(nrm) + dy * np.sin(nrm)) <= flat_radius
    return ok


def segmented_pupil(n, seg_flat_px, gap_px=1.5, strut_px=2.0):
    """웹 망원경을 흉내 낸 구멍: 육각 조각 18장(가운데 한 장 빠짐) + 부경 받침대 3개."""
    y, x = np.mgrid[:n, :n] - (n - 1) / 2
    p = np.zeros((n, n), bool)
    pitch = 2 * seg_flat_px + gap_px                 # 이웃한 조각 중심 사이 거리 (변을 맞댄 방향)
    dirs = [np.radians(a) for a in (90, 30, -30, -90, -150, 150)]
    centers = set()
    for i in range(-2, 3):
        for j in range(-2, 3):
            # 축 좌표 (i, j) → 직교 좌표. 두 기저는 변에 수직인 90° 와 30° 방향
            cx = pitch * (i * np.cos(dirs[0]) + j * np.cos(dirs[1]))
            cy = pitch * (i * np.sin(dirs[0]) + j * np.sin(dirs[1]))
            ring = max(abs(i), abs(j), abs(i + j))
            if 1 <= ring <= 2:
                centers.add((round(cx, 6), round(cy, 6)))
    for cx, cy in centers:
        p |= hexagon_mask(x, y, cx, cy, seg_flat_px)
    # 받침대: 위로 하나, 아래 양쪽으로 둘 (아래 둘은 육각형 변과 나란하다)
    r = np.hypot(x, y)
    for ang in (90, -150, -30):
        a = np.radians(ang)
        along = x * np.cos(a) + y * np.sin(a)
        across = -x * np.sin(a) + y * np.cos(a)
        p &= ~((along > 0) & (np.abs(across) < strut_px / 2) & (r < n / 2))
    return p.astype(float)


# ---------------------------------------------------------------- 초점 근처의 장

def focal_field(v, u, n=2001):
    """수차 없는 원형 동공으로 모은 빛의 초점 근처 장 (보른과 울프 8.8 의 광학 좌표):
    U(u, v) = 2 ∫₀¹ J0(vρ) exp(−i u ρ²/2) ρ dρ.  v = k NA r, u = k NA² z. 초점 중심에서 1 로 정규화."""
    rho = np.linspace(0, 1, n)
    v = np.atleast_1d(v)
    u = np.atleast_1d(u)
    out = np.empty((len(u), len(v)), complex)
    for i, uu in enumerate(u):
        ph = np.exp(-1j * uu * rho**2 / 2) * rho
        out[i] = 2 * np.trapezoid(ph[None, :] * j0(v[:, None] * rho[None, :]), rho, axis=1)
    return out


# ---------------------------------------------------------------- 격자

def slit_array_far_field(sin_theta, n_slits, period, width, lam, pts_per_slit=24):
    """틈 n_slits 개(폭 width, 간격 period)의 원거리 장을 2차 파원의 위상자 합으로 직접 구한다 (정면 입사).
    틈 하나를 pts_per_slit 개 점으로 나눈다. 모든 점이 같은 위상일 때 1 이 되게 정규화."""
    k = 2 * np.pi / lam
    xs = (np.arange(pts_per_slit) + 0.5) / pts_per_slit * width - width / 2
    pos = (np.arange(n_slits) - (n_slits - 1) / 2)[:, None] * period + xs[None, :]
    pos = pos.ravel()
    s = np.atleast_1d(sin_theta)
    out = np.empty(len(s), complex)
    for i in range(0, len(s), 512):
        blk = s[i:i + 512]
        out[i:i + 512] = np.exp(-1j * k * blk[:, None] * pos[None, :]).sum(1) / pos.size
    return out


def thin_grating_efficiencies(t_one_period, orders):
    """한 주기의 투과 함수 표본에서 각 회절 차수로 가는 세기 비율 (푸리에 계수의 제곱)."""
    c = np.fft.fft(t_one_period) / len(t_one_period)
    return np.array([np.abs(c[q]) ** 2 for q in orders])

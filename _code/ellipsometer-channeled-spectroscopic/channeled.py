"""채널 분광 타원계측기의 신호 합성과 푸리에 복조.

두꺼운 복굴절판의 지연량이 파수에 선형으로 늘어난다는 성질을 써서, 편광 정보를
파장축 간섭무늬(채널)에 싣는다. 회전형이 시간축으로 하던 변조를 파장축에서 한다.

세 배치를 다룬다.

1. 지연자 하나  P(0) R(45,phi) S A(-45)       -> 채널 3개, Psi/Delta 측정
2. 지연자 둘    R(0,d1) R(45,d2) A(0)          -> 채널 7개, 스토크스 전체 (Oka 1999)
3. 이색성 있는 지연자 둘                        -> 채널 9개 (Okabe 2009 의 방해석)

신호 합성은 부품 뮬러 행렬의 곱으로만 한다. Oka 식 4, Hagen 식 1, Hu 식 13 의
닫힌 식은 verify_channeled.py 가 대조할 독립 경로로 남겨둔다.

부호 규약은 배경이론 2·3편, 타원계측기 1~4편과 일치시킨다.
"""
import numpy as np

# 석영의 1차 유효 복굴절 계수. Hagen 2022 가 0.4~1.045 um 대역에서 쓴 값이고,
# 같은 논문의 실측 채널 위치(t=3.89 mm -> OPD 39.3 um)에서 역산한 0.0101 과 맞는다.
# 2편의 공칭 복굴절 0.009 보다 약 10% 큰데, 분산 항이 더해지기 때문이다.
BETA_QUARTZ = 0.00998          # um 단위 OPD 를 um 두께로 나눈 값 (무차원)

# 방해석 589 nm (Okabe 2009)
N_O_CALCITE, N_E_CALCITE = 1.6548, 1.4864

D = np.deg2rad


# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소 (배경이론 2편, 타원계측기 3·4편과 동일한 정의)
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])


def mueller_polarizer(theta):
    base = 0.5 * np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


def mueller_retarder(theta, delta, gamma=np.pi / 4):
    """방위각 theta, 지연량 delta 인 선형 위상지연자.

    gamma 는 투과율 비 각도다. 빠른 축과 느린 축의 진폭 투과율을 (cos g, sin g) 로
    두므로 gamma = 45 도면 이색성이 없는 이상적 지연자가 된다. 방해석처럼 입·출사면의
    프레넬 반사율이 축마다 다르면 gamma 가 45 도에서 벗어난다(Okabe 2009).
    """
    c2g, s2g = np.cos(2 * gamma), np.sin(2 * gamma)
    c, s = np.cos(delta), np.sin(delta)
    base = np.array([[1, c2g, 0, 0],
                     [c2g, 1, 0, 0],
                     [0, 0, s2g * c, s2g * s],
                     [0, 0, -s2g * s, s2g * c]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


def mueller_sample(psi, delta):
    """등방성 시편. 3·4편의 mueller_sample 과 같고 김영준 2025 식 2.44 와도 같다."""
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    cd, sd = np.cos(delta), np.sin(delta)
    return np.array([[1, -c2, 0, 0], [-c2, 1, 0, 0],
                     [0, 0, s2 * cd, s2 * sd], [0, 0, -s2 * sd, s2 * cd]])


S_UNPOLARIZED = np.array([1.0, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# 두꺼운 지연자 — 지연량이 파수에 선형으로 늘어난다
# ---------------------------------------------------------------------------


def retardance(sigma, thickness, beta=BETA_QUARTZ):
    """파수 sigma [1/um] 에서 두께 thickness [um] 인 판의 지연량 [rad].

    phi = 2 pi beta t sigma. 기울기 L = beta * t 가 곧 광로차(OPD)다.
    """
    return 2 * np.pi * beta * thickness * sigma


def opd(thickness, beta=BETA_QUARTZ):
    """지연자의 광로차 L [um]. 푸리에 영역에서 채널이 서는 자리를 정한다."""
    return beta * thickness


def fresnel_gamma(n_fast, n_slow):
    """입·출사면의 프레넬 반사가 만드는 투과율 비 각도 gamma.

    재료 자체의 이색성이 0 이어도 굴절률이 축마다 다르면 면에서의 반사율이 달라진다.
    수직 입사, 양면, 무반사 코팅 없음을 가정한다.
    """
    t = lambda n: (1.0 - ((n - 1.0) / (n + 1.0)) ** 2) ** 2   # 세기 투과율
    return np.arctan(np.sqrt(t(n_slow) / t(n_fast)))


# ---------------------------------------------------------------------------
# 1. 지연자 하나 — Psi, Delta 를 재는 최소 구성
# ---------------------------------------------------------------------------


def spectrum_single(psi, delta, sigma, thickness, beta=BETA_QUARTZ,
                    source=None, pol=0.0, ret_az=np.pi / 4, ana=-np.pi / 4,
                    gamma=np.pi / 4, dphi=0.0):
    """P(0) R(45, phi) S A(-45) 배치의 검출 스펙트럼.

    분석기 방위각의 부호가 중요하다. +45 도로 두면 sin(Delta) 항의 부호가 뒤집혀
    이승우 2021 식 3.28 과 맞지 않는다(배치 탐색으로 확인했다).

    dphi 는 온도 따위로 지연량이 표류한 양이다(이승우 2021 식 6.2).
    """
    sigma = np.asarray(sigma, dtype=float)
    src = np.ones_like(sigma) if source is None else np.asarray(source, float)
    mp, ma = mueller_polarizer(pol), mueller_polarizer(ana)
    ms = mueller_sample(psi, delta)
    phi = retardance(sigma, thickness, beta) + dphi
    out = np.empty_like(sigma)
    for i, p in enumerate(phi):
        s = ms @ (mueller_retarder(ret_az, p, gamma) @ (mp @ S_UNPOLARIZED))
        out[i] = (ma @ s)[0]
    return src * out


# ---------------------------------------------------------------------------
# 2. 지연자 둘 — 스토크스 전체 (Oka 1999 / Hagen 2022 배치)
# ---------------------------------------------------------------------------


def spectrum_stokes(stokes, sigma, t1, t2, beta=BETA_QUARTZ,
                    az1=0.0, az2=np.pi / 4, ana=0.0,
                    gamma1=np.pi / 4, gamma2=np.pi / 4, source=None):
    """R1(0, d1) R2(45, d2) A(0) 를 지난 빛의 스펙트럼.

    stokes: (4, len(sigma)) 배열. 파수마다의 입사 스토크스 벡터.
    빛은 R1 을 먼저 지나므로 행렬은 A R2 R1 순으로 곱한다.
    """
    sigma = np.asarray(sigma, dtype=float)
    S = np.asarray(stokes, dtype=float)
    src = np.ones_like(sigma) if source is None else np.asarray(source, float)
    ma = mueller_polarizer(ana)
    d1 = retardance(sigma, t1, beta)
    d2 = retardance(sigma, t2, beta)
    out = np.empty_like(sigma)
    for i in range(len(sigma)):
        m = ma @ mueller_retarder(az2, d2[i], gamma2) @ mueller_retarder(az1, d1[i], gamma1)
        out[i] = (m @ S[:, i])[0]
    return src * out


# ---------------------------------------------------------------------------
# 푸리에 복조 — 파수축 신호를 광로차(OPD) 영역으로 옮긴다
# ---------------------------------------------------------------------------


def to_opd(signal, sigma):
    """파수축 스펙트럼을 OPD 영역으로 옮긴다. 반환: (h [um], 복소 C(h))."""
    sigma = np.asarray(sigma, dtype=float)
    n = len(sigma)
    dsig = (sigma[-1] - sigma[0]) / (n - 1)
    C = np.fft.fftshift(np.fft.fft(np.asarray(signal, float))) / n
    h = np.fft.fftshift(np.fft.fftfreq(n, d=dsig))
    return h, C


def channel_positions(t1, t2=None, beta=BETA_QUARTZ):
    """이상적인 경우 채널이 서는 OPD 자리.

    지연자 하나면 0, +-L. 둘이면 0, +-L2, +-(L1-L2), +-(L1+L2) 일곱이다.
    이색성이 있으면 여기에 +-L1 이 더해져 아홉이 된다(Okabe 2009).
    """
    L1 = opd(t1, beta)
    if t2 is None:
        return np.array([-L1, 0.0, L1])
    L2 = opd(t2, beta)
    return np.array(sorted({0.0, L2, -L2, L1 - L2, -(L1 - L2), L1 + L2, -(L1 + L2)}))


def extract_channel(signal, sigma, center, half_width, window="hann"):
    """OPD 영역에서 채널 하나를 잘라 파수축으로 되돌린다.

    창을 씌우는 이 연산이 곧 복조이고, 동시에 대역폭 제한을 만든다(Okabe 2009 3.2절).
    """
    h, C = to_opd(signal, sigma)
    w = np.zeros_like(h)
    sel = np.abs(h - center) <= half_width
    if window == "hann":
        x = (h[sel] - center) / half_width
        w[sel] = 0.5 * (1.0 + np.cos(np.pi * x))
    else:
        w[sel] = 1.0
    return np.fft.ifft(np.fft.ifftshift(C * w)) * len(sigma)


# ---------------------------------------------------------------------------
# 두께의 허용 창 — 아래는 채널 분리, 위는 나이퀴스트
# ---------------------------------------------------------------------------


def opd_min_gaussian(sigma_min, sigma_max, n_std=3.0):
    """채널이 겹치지 않을 OPD 하한 (이승우 2021 식 3.26).

    광원을 파수에 대한 가우스 분포로 두고 파장 폭을 6 STD 로 잡으면 h 영역의
    표준편차가 6/(smax-smin) 이 된다. 기저 봉우리와 측대역이 각각 n_std 배만큼
    떨어져야 하므로 L > 2 * n_std * STD_h 다. n_std = 1.5 면 논문의 18/(smax-smin).
    """
    std_h = 6.0 / (sigma_max - sigma_min)
    return 2.0 * n_std * std_h


def opd_max_nyquist(dsigma, n_channel_max):
    """표본화가 허용하는 OPD 상한 (Hu 2024 식 10·11).

    최고 반송 주파수가 표본화 주파수 1/dsigma 의 절반을 넘으면 되접힌다.
    n_channel_max 는 기본 주파수 단위로 센 최고 채널 차수이고,
    스토크스 CSP 는 m+n, 뮬러 CSP 는 m+n+p+q 다.
    """
    return 1.0 / (dsigma * (2.0 * n_channel_max + 1.0))


def dsigma_avg(lam_min, lam_max, dlam):
    """분산형 분광기의 파수 분해능 기하평균 (Hu 2024 / Hagen 2022).

    파장 등간격 표본화라 파수 분해능이 균일하지 않다. 짧은 파장 쪽이 나쁘다.
    """
    return dlam / (lam_min * lam_max)


def thickness_window(lam_min, lam_max, dlam, ratio=(1, 2), beta=BETA_QUARTZ,
                     n_std=3.0):
    """기본 두께 t0 의 허용 범위 [um]. 반환: (t_min, t_max).

    ratio 는 서로소인 두께비다. 스토크스 CSP 는 (m, n), 뮬러 CSP 는 (m, n, p, q).
    """
    s_min, s_max = 1.0 / lam_max, 1.0 / lam_min
    t_min = opd_min_gaussian(s_min, s_max, n_std) / beta
    t_max = opd_max_nyquist(dsigma_avg(lam_min, lam_max, dlam), sum(ratio)) / beta
    return t_min, t_max


def resolving_power(lam_min, lam_max, dlam):
    """파수 영역의 분해 가능한 점 수 (Hagen 2022 식 7)."""
    s_min, s_max = 1.0 / lam_max, 1.0 / lam_min
    return (s_max - s_min) / dsigma_avg(lam_min, lam_max, dlam)


# ---------------------------------------------------------------------------
# 시편 — 단층막의 Psi, Delta (스펙트럼 구조를 만들기 위한 최소 모형)
# ---------------------------------------------------------------------------


def film_psi_delta(sigma, d_nm, n_film=1.46, n_sub=3.88 - 0.02j, aoi=D(70.0)):
    """기판 위 단층막의 Psi, Delta. 굴절률은 상수로 둔 최소 모형이다.

    막이 두꺼울수록 Psi, Delta 가 파수에 대해 빨리 진동한다. 채널 분광에서
    "시편의 스펙트럼 변화가 반송 변조보다 완만해야 한다"를 보이는 데 쓴다.
    """
    sigma = np.asarray(sigma, dtype=float)
    n0 = 1.0
    s0 = n0 * np.sin(aoi)
    c0 = np.cos(aoi)
    c1 = np.sqrt(1.0 - (s0 / n_film) ** 2 + 0j)
    c2 = np.sqrt(1.0 - (s0 / n_sub) ** 2 + 0j)

    def fresnel(ni, ci, nt, ct):
        rp = (nt * ci - ni * ct) / (nt * ci + ni * ct)
        rs = (ni * ci - nt * ct) / (ni * ci + nt * ct)
        return rp, rs

    r01p, r01s = fresnel(n0, c0, n_film, c1)
    r12p, r12s = fresnel(n_film, c1, n_sub, c2)
    # 위상 두께. sigma 가 1/um 이고 d_nm 이 nm 이므로 1e-3 을 곱해 um 로 맞춘다
    beta_f = 2 * np.pi * (d_nm * 1e-3) * n_film * c1 * sigma
    z = np.exp(-2j * beta_f)
    rp = (r01p + r12p * z) / (1.0 + r01p * r12p * z)
    rs = (r01s + r12s * z) / (1.0 + r01s * r12s * z)
    rho = rp / rs
    return np.arctan(np.abs(rho)), np.angle(rho)

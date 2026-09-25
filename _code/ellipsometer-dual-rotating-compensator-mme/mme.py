"""이중 회전 보상자 뮬러행렬 타원계측기(DRC-MME)의 검출 세기와 데이터 환산.

배치는 P C1r(m1*C) S C2r(m2*C) A 다. 고정 편광자·분석기 사이에서 보상자 둘이
서로 다른 속도로 돈다. C = wt 가 기본 회전각이고 pi/w 가 기본 광학주기다.

계산은 부품 뮬러 행렬의 곱으로만 한다. Collins & Koh 의 닫힌 식은 쓰지 않는다 —
verify_mme.py 가 그 식들과 대조해야 하므로 독립 경로로 남겨둔다.

부호 규약은 배경이론 2·3편, 타원계측기 1·2편과 일치시킨다.
"""
import numpy as np

# Collins & Koh 1999 가 제안한 회전비와 검출 조건
RATIO = (5, 3)
N_SCAN = 36           # 광학주기당 검출 스캔 수
DETECTOR_MIN_MS = 5.1  # 상용 포토다이오드 어레이의 최소 스캔 시간


# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소 (배경이론 2편과 동일한 정의)
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])


def mueller_polarizer(theta):
    base = 0.5 * np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


def mueller_retarder(theta, delta):
    """방위각 theta, 지연량 delta 인 선형 위상지연자."""
    c, s = np.cos(delta), np.sin(delta)
    base = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, c, s], [0, 0, -s, c]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


S_UNPOLARIZED = np.array([1.0, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# 검출 세기
# ---------------------------------------------------------------------------


def _retarder_stack(thetas, delta):
    """방위각 배열에 대한 위상지연자 뮬러 행렬 묶음. shape (len, 4, 4)."""
    c2, s2 = np.cos(2 * thetas), np.sin(2 * thetas)
    cd, sd = np.cos(delta), np.sin(delta)
    n = len(thetas)
    R = np.zeros((n, 4, 4))
    R[:, 0, 0] = 1.0
    R[:, 1, 1] = c2 ** 2 + cd * s2 ** 2
    R[:, 1, 2] = c2 * s2 * (1 - cd)
    R[:, 1, 3] = -s2 * sd
    R[:, 2, 1] = c2 * s2 * (1 - cd)
    R[:, 2, 2] = s2 ** 2 + cd * c2 ** 2
    R[:, 2, 3] = c2 * sd
    R[:, 3, 1] = s2 * sd
    R[:, 3, 2] = -c2 * sd
    R[:, 3, 3] = cd
    return R


def intensity(M, angles, delta1, delta2, ratio=RATIO, pol=0.0, ana=0.0,
              cs1=0.0, cs2=0.0):
    """P C1r S C2r A 배치의 검출 세기 I(C).

    angles: 기본 회전각 C 배열(rad). 보상자 방위각은 m1*(C-cs1), m2*(C-cs2) 다.
    cs1, cs2 는 보상자 위상 이동각으로, 교정에서 구하는 값이다(기본값 0).
    각도 축에 대해 벡터화했다 — 조건수를 파장·지연량에 대해 훑으려면 필요하다.
    """
    m1, m2 = ratio
    angles = np.asarray(angles, dtype=float)
    mp, ma = mueller_polarizer(pol), mueller_polarizer(ana)
    R1 = _retarder_stack(m1 * (angles - cs1), delta1)
    R2 = _retarder_stack(m2 * (angles - cs2), delta2)
    s = M @ (R1 @ (mp @ S_UNPOLARIZED))[..., None]      # (n,4,1)
    s = R2 @ s
    # 분석기의 첫 행과의 내적. (n,4)@(4,) 꼴 matmul 은 일부 BLAS 에서 가짜
    # 부동소수점 경고를 올리므로 원소곱 합으로 쓴다(결과는 동일하다).
    return (s[..., 0] * ma[0]).sum(axis=1)


def optical_cycle(n_pts):
    """기본 광학주기 한 바퀴. C 는 0 에서 pi 까지 돈다."""
    return np.linspace(0.0, np.pi, n_pts, endpoint=False)


# ---------------------------------------------------------------------------
# 푸리에 계수
# ---------------------------------------------------------------------------


def fourier_coefficients(signal, n_max=20):
    """I(C) 를 dc + sum_n [a_n cos(2nC) + b_n sin(2nC)] 로 분해한다.

    C 가 0~pi 를 도므로 cos(2nC) 는 이 구간에서 정확히 n 주기를 돈다.
    따라서 실수 FFT 의 n 번째 성분이 곧 고조파 차수 n 이다.
    반환: dc, {n: (a_n, b_n)}
    """
    n_pts = len(signal)
    F = np.fft.rfft(signal) / n_pts
    return F[0].real, {n: (2 * F[n].real, -2 * F[n].imag)
                       for n in range(1, n_max + 1)}


def harmonic_magnitudes(signal, n_max=20, normalize=True):
    """차수별 크기 sqrt(a^2+b^2). normalize 면 dc 로 나눈다."""
    dc, co = fourier_coefficients(signal, n_max)
    mag = np.array([np.hypot(*co[n]) for n in range(1, n_max + 1)])
    return mag / abs(dc) if normalize else mag


def highest_harmonic(ratio, delta1=np.pi / 2, delta2=np.pi / 2, n_max=40,
                     tol=1e-9, n_pts=4096, seed=0):
    """주어진 회전비에서 실제로 살아남는 최고차 고조파를 수치로 찾는다."""
    rng = np.random.default_rng(seed)
    M = rng.uniform(-0.4, 0.4, (4, 4))
    M[0, 0] = 1.0
    mag = harmonic_magnitudes(
        intensity(M, optical_cycle(n_pts), delta1, delta2, ratio), n_max)
    nz = np.nonzero(mag > tol)[0]
    return int(nz[-1] + 1) if len(nz) else 0


# ---------------------------------------------------------------------------
# 데이터 환산 행렬
# ---------------------------------------------------------------------------


def reduction_matrix(delta1, delta2, ratio=RATIO, n_max=20, n_pts=128,
                     tol=1e-9, **kw):
    """푸리에 계수 벡터 = W @ (뮬러 요소 16개) 의 W 를 만든다.

    I 가 M 에 선형이므로 단위행렬 E_jk 를 하나씩 넣으면 열이 그대로 나온다.
    항등적으로 0 인 행(소멸 고조파)은 빼고 돌려준다.

    표본 수는 n_max 차 고조파를 담을 만큼만 있으면 된다(나이퀴스트). 기본 128 점은
    n_max = 20 에 충분하다.
    """
    assert n_pts > 2 * n_max, "표본이 최고차 고조파를 담지 못한다"
    cols = []
    C = optical_cycle(n_pts)
    for j in range(4):
        for k in range(4):
            E = np.zeros((4, 4))
            E[j, k] = 1.0
            dc, co = fourier_coefficients(intensity(E, C, delta1, delta2, ratio, **kw),
                                          n_max)
            v = [dc]
            for n in range(1, n_max + 1):
                v += list(co[n])
            cols.append(v)
    W = np.array(cols).T
    return W[np.abs(W).max(axis=1) > tol]


def condition_number(delta1, delta2, ratio=RATIO, **kw):
    """데이터 환산 행렬의 조건수. 측정 잡음이 뮬러 요소로 얼마나 증폭되는지를 잰다."""
    return np.linalg.cond(reduction_matrix(delta1, delta2, ratio, **kw))


def optimal_retardance(ratio=RATIO, lo=60.0, hi=175.0, step=0.5, **kw):
    """조건수를 최소로 만드는 지연량(도)을 훑어서 찾는다. delta1 = delta2 로 둔다."""
    deg = np.arange(lo, hi + 1e-9, step)
    cond = np.array([condition_number(np.deg2rad(d), np.deg2rad(d), ratio, **kw)
                     for d in deg])
    return deg[int(np.argmin(cond))], cond.min(), deg, cond

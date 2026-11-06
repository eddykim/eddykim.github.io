"""후초점면 각도 분해와 채널 분광을 묶은 스냅샷 타원계측기의 신호 합성과 해석.

5편의 채널 분광은 편광을 파장축 간섭무늬에 실었고, 6편의 후초점면은 반경을 입사각으로
펼쳤다. 이 장비는 슬릿으로 후초점면의 지름 한 줄(방위각 Omega = 90 도)만 잘라 라인 스캔
분광기(LSS)에 넣는다. 카메라 한 축이 입사각, 다른 축이 파수를 싣고, 파수축 안의 무늬가
a = cos 2Psi, b = sin 2Psi sin Delta 를 싣는다.

  배치(빛이 지나는 순서)  P(0) -> R(45, phi) -> [빔 스플리터] -> 대물렌즈 -> 시편
                           -> 대물렌즈 -> [빔 스플리터] -> A(-45) -> 슬릿 -> LSS

신호 합성은 부품 뮬러 행렬의 곱으로만 한다. 저널 식 5 와 a, b 의 닫힌 식은
verify_acse.py 가 대조할 독립 경로로 남겨 둔다.

해석은 세 단계다.
  1. 극값 포락선 방법(EEM)으로 기저 신호를 세운다
  2. 힐베르트 변환으로 해석 신호를 만들고 반송파를 빼서 a, b 를 읽는다
  3. 기준 시편 하나로 만든 진폭·위상 오차 지도로 다른 시편의 a, b 를 바로잡는다

부호 규약은 배경이론 2·3편, 타원계측기 1~6편과 같다. 저널이 분석기를 45 도로 적은 것은
출력 쪽 좌표가 90 도 돌아간 표기이고, 이 시리즈 규약으로는 -45 도다.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize_scalar
from scipy.signal import hilbert

D = np.deg2rad

# ---------------------------------------------------------------------------
# 장비 수치 — 저널판(Kim, Lee, Pahk 2025) 3.A 절에 실린 값
# ---------------------------------------------------------------------------
LAM_MIN, LAM_MAX = 0.475, 0.675          # um
THETA_MAX = D(60.0)                       # 최대 입사각
T_RETARDER = 4092.3                       # 석영 다중차 지연자 두께 [um]
DLAM_LSS = 0.00192                        # LSS 분광 해상도 1.92 nm [um]
BETA_QUARTZ = 0.00998                     # 5편과 같은 석영 1차 유효 복굴절 계수

# 시편 굴절률은 6편과 같은 상수 최소 모형
N_SIO2, N_SI = 1.460, 3.875 - 0.0156j


# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소 — 배열을 받아 (..., 4, 4) 로 돌려준다
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    omega = np.asarray(omega, dtype=float)
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    M = np.zeros(omega.shape + (4, 4))
    M[..., 0, 0] = M[..., 3, 3] = 1.0
    M[..., 1, 1] = M[..., 2, 2] = c
    M[..., 1, 2], M[..., 2, 1] = s, -s
    return M


def mueller_polarizer(theta):
    base = 0.5 * np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return _mm(_mm(mueller_rotation(-theta), np.broadcast_to(base, np.shape(theta) + (4, 4))),
               mueller_rotation(theta))


def mueller_retarder(theta, delta, diatt=0.0):
    """방위각 theta, 지연량 delta, 이색성 diatt 인 선형 소자.

    diatt = (Tx - Ty)/(Tx + Ty). 0 이면 이상적 지연자다. 5편의 gamma 로 쓰면
    diatt = cos 2gamma 다.
    """
    delta = np.asarray(delta, dtype=float)
    diatt = np.broadcast_to(np.asarray(diatt, dtype=float), delta.shape)
    sq = np.sqrt(1.0 - diatt ** 2)
    c, s = np.cos(delta), np.sin(delta)
    B = np.zeros(delta.shape + (4, 4))
    B[..., 0, 0] = B[..., 1, 1] = 1.0
    B[..., 0, 1] = B[..., 1, 0] = diatt
    B[..., 2, 2] = B[..., 3, 3] = sq * c
    B[..., 2, 3], B[..., 3, 2] = sq * s, -sq * s
    return _mm(_mm(mueller_rotation(-theta), B), mueller_rotation(theta))


def mueller_sample(psi, delta):
    """등방성 시편. 3~6편의 mueller_sample 과 같다 (m00 = 1 로 정규화)."""
    psi, delta = np.broadcast_arrays(np.asarray(psi, float), np.asarray(delta, float))
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    cd, sd = np.cos(delta), np.sin(delta)
    M = np.zeros(psi.shape + (4, 4))
    M[..., 0, 0] = M[..., 1, 1] = 1.0
    M[..., 0, 1] = M[..., 1, 0] = -c2
    M[..., 2, 2] = M[..., 3, 3] = s2 * cd
    M[..., 2, 3], M[..., 3, 2] = s2 * sd, -s2 * sd
    return M


def mueller_depolarizer(v):
    """편광 성분을 v 배로 줄이는 등방 편광 소멸 소자."""
    v = np.asarray(v, dtype=float)
    M = np.zeros(v.shape + (4, 4))
    M[..., 0, 0] = 1.0
    M[..., 1, 1] = M[..., 2, 2] = M[..., 3, 3] = v
    return M


S_UNPOLARIZED = np.array([1.0, 0.0, 0.0, 0.0])


def _mm(X, Y):
    """쌓인 4x4 행렬의 곱. numpy 2.x 의 matmul 이 내는 가짜 경고를 피하려고 einsum 을 쓴다."""
    return np.einsum("...ij,...jk->...ik", X, Y)


def _chain(*Ms):
    """왼쪽부터 곱한다. 빛은 오른쪽 행렬부터 지난다."""
    out = Ms[0]
    for M in Ms[1:]:
        out = _mm(out, M)
    return out


# ---------------------------------------------------------------------------
# 시편 — 단층막의 반사계수 (파수, 입사각 격자 위에서)
# ---------------------------------------------------------------------------


def film_rp_rs(sigma, theta, d_nm, n_film=N_SIO2, n_sub=N_SI, n_amb=1.0):
    """파수 sigma [1/um], 입사각 theta [rad] 에서 단층막의 rp, rs.

    d_nm = 0 이면 맨 기판이다. 굴절률은 파장에 무관한 상수로 둔다(6편의 최소 모형).
    """
    sigma, theta = np.broadcast_arrays(np.asarray(sigma, float),
                                       np.abs(np.asarray(theta, float)))
    s0 = n_amb * np.sin(theta) + 0j
    c0 = np.cos(theta) + 0j
    c1 = np.sqrt(1.0 - (s0 / n_film) ** 2)
    c2 = np.sqrt(1.0 - (s0 / n_sub) ** 2)

    def fres(ni, ci, nt, ct):
        return ((nt * ci - ni * ct) / (nt * ci + ni * ct),
                (ni * ci - nt * ct) / (ni * ci + nt * ct))

    r01p, r01s = fres(n_amb, c0, n_film, c1)
    r12p, r12s = fres(n_film, c1, n_sub, c2)
    z = np.exp(-2j * 2 * np.pi * (d_nm * 1e-3) * n_film * c1 * sigma)
    rp = (r01p + r12p * z) / (1.0 + r01p * r12p * z)
    rs = (r01s + r12s * z) / (1.0 + r01s * r12s * z)
    return rp, rs


def film_psi_delta(sigma, theta, d_nm, **kw):
    """Psi, Delta 와 비편광 반사율 R = (|rp|^2 + |rs|^2)/2."""
    rp, rs = film_rp_rs(sigma, theta, d_nm, **kw)
    rho = rp / rs
    return np.arctan(np.abs(rho)), np.angle(rho), 0.5 * (np.abs(rp) ** 2 + np.abs(rs) ** 2)


# ---------------------------------------------------------------------------
# 지연자와 광원
# ---------------------------------------------------------------------------


def retardance(sigma, thickness=T_RETARDER, beta=BETA_QUARTZ):
    """5편과 같은 선형 지연량 phi = 2 pi L sigma, L = beta t."""
    return 2 * np.pi * beta * thickness * np.asarray(sigma, float)


def opd(thickness=T_RETARDER, beta=BETA_QUARTZ):
    return beta * thickness


def fringe_count(thickness=T_RETARDER, lam_min=LAM_MIN, lam_max=LAM_MAX,
                 beta=BETA_QUARTZ):
    """대역 안의 무늬 수 N = L (sigma_max - sigma_min)."""
    return opd(thickness, beta) * (1.0 / lam_min - 1.0 / lam_max)


def led_spectrum(sigma):
    """백색 LED 를 흉내 낸 광원 스펙트럼(청색 칩 + 형광체). 최댓값 1 근처.

    실제 장비의 광원 스펙트럼이 아니라, 기저 신호가 파수에 따라 크게 변하는
    상황을 만들기 위한 가상 모양이다.
    """
    lam = 1.0 / np.asarray(sigma, float) * 1e3      # nm
    return (0.55 * np.exp(-((lam - 455.0) / 14.0) ** 2)
            + 1.0 * np.exp(-((lam - 570.0) / 80.0) ** 2))


def sigma_grid(n=800, lam_min=LAM_MIN, lam_max=LAM_MAX):
    """파수 등간격 격자 [1/um]. 짧은 파수(긴 파장) 쪽부터."""
    return np.linspace(1.0 / lam_max, 1.0 / lam_min, n)


def theta_grid(n=241, theta_max=THETA_MAX):
    """후초점면 지름을 따라 본 입사각. 음수는 Omega = 270 도 쪽 반쪽이다."""
    return np.linspace(-theta_max, theta_max, n)


# ---------------------------------------------------------------------------
# 장비 — 부품 뮬러 행렬의 곱
# ---------------------------------------------------------------------------


IDEAL = dict(t_ret=T_RETARDER, phi0=0.0, kappa=0.0, bs_diatt=0.0, bs_ret=0.0,
             obj_ret=0.0, obj_diatt=0.0, depol=1.0, lss_blur=0.0, pol_c=D(45.0))

# 그림 4 의 가상 부품 오차. 크기는 문헌이 보고한 범위 안에서 고른 예시값이고,
# 실제 장비의 측정값이 아니다.
VIRTUAL_ERRORS = dict(
    t_ret=T_RETARDER,
    phi0=D(25.0),        # 지연량 상수 오프셋 (5편 1절이 보정으로 뺀다고 한 몫)
    kappa=4.0,           # 지연량의 2차 분산 [rad um^2] (Hagen 2020 의 비선형성 꼴)
    bs_diatt=0.04,       # 빔 스플리터 이색성 (실험실 좌표 0 도 축, 왕복)
    bs_ret=D(4.0),       # 빔 스플리터 지연 (같은 축, 왕복)
    obj_ret=D(6.0),      # 대물렌즈 지연, 최대 입사각에서. (sin theta / sin 60)^2 로 늘어난다
    obj_diatt=0.03,      # 대물렌즈 이색성, 같은 꼴
    depol=0.97,          # 빔 발산으로 지연량이 퍼져 무늬 가시도가 주는 몫(조명 쪽)
    lss_blur=DLAM_LSS,   # LSS 분광 해상도
    pol_c=D(45.0),
)


def true_retardance(sigma, err):
    """부품 오차가 든 실제 지연량. 해석 쪽은 retardance() 의 공칭값만 안다."""
    sigma = np.asarray(sigma, float)
    sbar = 0.5 * (1.0 / LAM_MIN + 1.0 / LAM_MAX)
    return (retardance(sigma, err["t_ret"]) + err["phi0"]
            + err["kappa"] * (sigma - sbar) ** 2)


def system_mueller(sigma, theta, d_nm, err=IDEAL, **film):
    """(sigma, theta) 격자의 계 전체 뮬러 행렬 첫 행(검출기는 세기만 본다).

    theta > 0 은 Omega = 90 도, theta < 0 은 Omega = 270 도 반쪽이다.
    """
    SIG, TH = np.meshgrid(np.asarray(sigma, float), np.asarray(theta, float))
    psi, dl, R = film_psi_delta(SIG, TH, d_nm, **film)
    omega = np.where(TH >= 0, D(90.0), D(270.0))
    phi = true_retardance(SIG, err)

    shape = SIG.shape
    zeros = np.zeros(shape)
    Mp = mueller_polarizer(zeros)
    Mc = mueller_retarder(np.full(shape, err["pol_c"]), phi)
    Ma = mueller_polarizer(np.full(shape, D(-45.0)))
    Mbs = mueller_retarder(zeros, np.full(shape, err["bs_ret"]), err["bs_diatt"])
    w = (np.sin(TH) / np.sin(THETA_MAX)) ** 2
    # 대물렌즈 오차는 시편 좌표(p, s 축)에 붙어 있다고 둔다. 왕복으로 두 번 지난다
    Mobj = mueller_retarder(zeros, err["obj_ret"] * w, err["obj_diatt"] * w)
    Ms = mueller_sample(psi, dl)
    Mdep = mueller_depolarizer(np.full(shape, err["depol"]))
    Rf, Rb = mueller_rotation(omega), mueller_rotation(-omega)

    # 가시도 저하는 조명 쪽에 둔다: 빔 발산으로 두꺼운 지연자 안의 지연량이 퍼지면
    # cos phi, sin phi 가 함께 평균되어 a, b 두 항이 같은 비율로 준다
    M = _chain(Ma, Mbs, Rb, Mobj, Ms, Mobj, Rf, Mbs, Mdep, Mc, Mp)
    return SIG, TH, R, M[..., 0, :]


def detect(sigma, theta, d_nm, err=IDEAL, source=True, n_sub_lss=15, **film):
    """검출 영상 I(theta, sigma). 행이 입사각, 열이 파수다.

    LSS 분해능은 화소마다 파장 폭 lss_blur 의 사각 창으로 적분하는 것으로 둔다.
    """
    sigma = np.asarray(sigma, float)
    blur = err["lss_blur"]
    if blur > 0:
        lam = 1.0 / sigma
        offs = (np.arange(n_sub_lss) + 0.5) / n_sub_lss - 0.5
        acc = 0.0
        for o in offs:
            s_sub = 1.0 / (lam + o * blur)
            acc = acc + _detect_once(s_sub, theta, d_nm, err, source, **film)
        return acc / n_sub_lss
    return _detect_once(sigma, theta, d_nm, err, source, **film)


def _detect_once(sigma, theta, d_nm, err, source, **film):
    SIG, TH, R, row = system_mueller(sigma, theta, d_nm, err, **film)
    I = row[..., 0]                     # 비편광 입력이면 첫 행의 첫 성분
    src = led_spectrum(SIG) if source else 1.0
    return src * R * I


def ideal_ab(sigma, theta, d_nm, **film):
    """이상적 장비가 그 시편에서 내야 할 a, b 를 뮬러 곱으로 구한다.

    지연량을 0, pi/2, pi, 3pi/2 로 바꿔 넣은 네 세기에서 푸리에 계수를 읽는다.
    a = cos 2Psi, b = sin 2Psi sin Delta 라는 닫힌 식은 쓰지 않는다(verify 가 대조).
    """
    SIG, TH = np.meshgrid(np.asarray(sigma, float), np.asarray(theta, float))
    psi, dl, _ = film_psi_delta(SIG, TH, d_nm, **film)
    omega = np.where(TH >= 0, D(90.0), D(270.0))
    zeros = np.zeros(SIG.shape)
    I = []
    for ph in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2):
        M = _chain(mueller_polarizer(np.full(SIG.shape, D(-45.0))),
                   mueller_rotation(-omega), mueller_sample(psi, dl), mueller_rotation(omega),
                   mueller_retarder(np.full(SIG.shape, D(45.0)), np.full(SIG.shape, ph)),
                   mueller_polarizer(zeros))
        I.append(M[..., 0, 0])
    dc = 0.5 * (I[0] + I[2])
    return (I[0] - I[2]) / (2 * dc), (I[1] - I[3]) / (2 * dc)


# ---------------------------------------------------------------------------
# 1. 기저 신호 — 극값 포락선 방법(EEM)과 푸리에 저역통과(FT)
# ---------------------------------------------------------------------------


def find_extrema(sig, y, half=2):
    """수치 미분의 영점 교차로 극값을 찾고 국소 2차 맞춤으로 다듬는다.

    반환: (극대 파수, 극대 세기), (극소 파수, 극소 세기)
    """
    dy = np.diff(y)
    k = np.nonzero(dy[:-1] * dy[1:] < 0)[0] + 1
    mx, mn = [], []
    for i in k:
        lo, hi = max(i - half, 0), min(i + half + 1, len(y))
        if hi - lo < 3:
            continue
        x0 = sig[i]
        c = np.polyfit(sig[lo:hi] - x0, y[lo:hi], 2)
        if c[0] == 0:
            continue
        xv = -c[1] / (2 * c[0])
        if abs(xv) > (sig[hi - 1] - sig[lo]):
            xv = 0.0
        yv = np.polyval(c, xv)
        (mx if c[0] < 0 else mn).append((x0 + xv, yv))
    return np.array(mx).T, np.array(mn).T


def eem_baseline(sig, y, L=None):
    """위·아래 포락선의 평균을 기저 신호로 삼는다(저널 식 6·7).

    포락선은 극값을 잇는 3차 스플라인이고, 마지막 극값 너머는 외삽이다.
    무늬가 거의 없어 극값이 네 개도 안 되는 줄(수직 입사 근처)은 포락선을 세울 수 없으므로
    푸리에 저역통과로 대신한다.
    """
    mx, mn = find_extrema(sig, y)
    if mx.size < 8 or mn.size < 8:
        dc = ft_baseline(sig, y, L)
        return dc, dc, dc
    (xu, yu), (xl, yl) = mx, mn
    up = CubicSpline(xu, yu, bc_type="natural", extrapolate=True)(sig)
    lo = CubicSpline(xl, yl, bc_type="natural", extrapolate=True)(sig)
    return 0.5 * (up + lo), up, lo


def _opd_axis(sig):
    n = len(sig)
    return np.fft.fftfreq(n, d=(sig[-1] - sig[0]) / (n - 1))


def ft_channel(sig, y, center, half_width):
    """5편의 extract_channel 과 같은 연산: OPD 영역에서 한 채널을 Hann 창으로 자른다."""
    h = _opd_axis(sig)
    C = np.fft.fft(y)
    x = (h - center) / half_width
    w = np.where(np.abs(x) <= 1, 0.5 * (1 + np.cos(np.pi * x)), 0.0)
    return np.fft.ifft(C * w)


def ft_baseline(sig, y, L=None):
    """0 채널만 남겨 기저 신호로 삼는다. 창 반폭은 L/2 (5편 7절과 같다)."""
    L = opd() if L is None else L
    return np.real(ft_channel(sig, y, 0.0, L / 2))


# ---------------------------------------------------------------------------
# 2. 복조 — 힐베르트 변환(HT)과 푸리에 채널 복조
# ---------------------------------------------------------------------------


def ht_demod(sig, y, dc, phi_nom):
    """정규화 무늬에 힐베르트 변환을 걸어 a, b 를 읽는다.

    AC' = y/dc - 1 = A' cos(phi - chi). 해석 신호 Z = A' exp(j(phi - chi)) 에
    exp(-j phi) 를 곱하면 A' exp(-j chi) = a - j b 다. 그러므로
    a = Re, **b = -Im** 이다. 학위논문 식 3.14 와 저널 식 11 은 이 부호를 서로 다른
    자리에서 놓쳤다.
    """
    ac = y / dc - 1.0
    z = hilbert(ac) * np.exp(-1j * phi_nom)
    return np.real(z), -np.imag(z), ac


def ft_demod(sig, y, phi_nom, L=None):
    """5편의 푸리에 복조: +L 채널과 0 채널의 비에서 a - j b 를 얻는다."""
    L = opd() if L is None else L
    c0 = np.real(ft_channel(sig, y, 0.0, L / 2))
    cp = ft_channel(sig, y, L, L / 2)
    z = 2 * cp * np.exp(-1j * phi_nom) / c0
    return np.real(z), -np.imag(z)


def analyze_image(sig, img, method="eem_ht", thickness=T_RETARDER):
    """영상 전체를 행(입사각)마다 해석한다. 반환: a, b 영상.

    해석 쪽은 지연자의 공칭 두께만 안다. 공칭 지연량으로 반송파를 뺀다.
    """
    phi_nom = retardance(sig, thickness)
    L = opd(thickness)
    A = np.empty_like(img)
    B = np.empty_like(img)
    for i, y in enumerate(img):
        if method == "eem_ht":
            dc = eem_baseline(sig, y, L)[0]
            A[i], B[i], _ = ht_demod(sig, y, dc, phi_nom)
        else:
            A[i], B[i] = ft_demod(sig, y, phi_nom, L)
    return A, B


# ---------------------------------------------------------------------------
# 3. 진폭·위상 오차 지도 — 기준 시편 하나로 만든다 (저널 식 16·17)
# ---------------------------------------------------------------------------


def amp_phase(a, b):
    """a - j b = A' exp(-j chi) 에서 진폭 A' 과 위상 오프셋 chi."""
    return np.hypot(a, b), np.arctan2(b, a)


def error_maps(a_meas, b_meas, a_ref, b_ref):
    """eps_A = A'_M / A'_ref,  eps_chi = chi_M - chi_ref.

    측정 무늬를 A'_ref eps_A cos(phi - chi_ref - eps_chi) 로 놓은 것이다.
    지연량이 공칭값에서 dphi 만큼 어긋나면 eps_chi 에 -dphi 가 실린다.
    곧 5편의 지연량 표류(상수 하나)를 (sigma, theta) 지도로 넓힌 것이다.
    """
    Am, cm = amp_phase(a_meas, b_meas)
    Ar, cr = amp_phase(a_ref, b_ref)
    eps_a = Am / np.where(Ar > 0, Ar, np.nan)
    eps_c = np.angle(np.exp(1j * (cm - cr)))
    return eps_a, eps_c


def apply_correction(a_meas, b_meas, eps_a, eps_c):
    """측정 시편의 a, b 에서 오차 지도를 되돌린다."""
    Am, cm = amp_phase(a_meas, b_meas)
    A = Am / eps_a
    c = cm - eps_c
    return A * np.cos(c), A * np.sin(c)


# ---------------------------------------------------------------------------
# 4. 두께 회귀 — a, b 전체를 한 모형으로 맞춘다 (저널 식 13 의 단순판)
# ---------------------------------------------------------------------------


def fit_thickness(sig, theta, a, b, d_guess, mask=None, span=60.0, **film):
    """막 두께 하나만 미지수로 두고 (a, b) 의 제곱합을 최소화한다."""
    if mask is None:
        mask = np.isfinite(a) & np.isfinite(b)

    def cost(d):
        at, bt = ideal_ab(sig, theta, d, **film)
        r = np.concatenate([(a - at)[mask], (b - bt)[mask]])
        return np.mean(r ** 2)

    grid = np.linspace(max(d_guess - span, 0.0), d_guess + span, 25)
    c = [cost(g) for g in grid]
    g0 = grid[int(np.argmin(c))]
    step = grid[1] - grid[0]
    res = minimize_scalar(cost, bounds=(max(g0 - step, 0.0), g0 + step),
                          method="bounded", options=dict(xatol=1e-4))
    return res.x, res.fun


# ---------------------------------------------------------------------------
# 5. LSS 분해능과 지연자 두께
# ---------------------------------------------------------------------------


def visibility_box(lam, thickness, dlam=DLAM_LSS, beta=BETA_QUARTZ):
    """파장 폭 dlam 의 사각 창이 무늬 진폭을 줄이는 비율.

    파수 폭은 dlam / lam^2 이고, 사각 창의 변환은 sinc 다.
    """
    dsig = dlam / np.asarray(lam, float) ** 2
    return np.abs(np.sinc(opd(thickness, beta) * dsig))


def thickness_window(lam_min=LAM_MIN, lam_max=LAM_MAX, dlam=DLAM_LSS,
                     beta=BETA_QUARTZ, n_std=3.0):
    """지연자 하나일 때 두께의 허용 창 [um] — 5편 thickness_window 의 단일 지연자판.

    하한: 기저 봉우리와 +-L 채널이 겹치지 않을 조건(5편 2절).
    상한: 최고 반송 주파수 L 이 표본화 한계를 넘지 않을 조건(5편 3절, 최고 차수 1).
    """
    s_min, s_max = 1.0 / lam_max, 1.0 / lam_min
    std_h = 6.0 / (2.0 * np.pi * (s_max - s_min))
    L_min = 2.0 * n_std * std_h
    dsig = dlam / (lam_min * lam_max)
    L_max = 1.0 / (dsig * 3.0)
    return L_min / beta, L_max / beta

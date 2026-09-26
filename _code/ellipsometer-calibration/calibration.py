"""회전분석기 타원계측기(RAE)의 보정을 오차 주입부터 회수까지 구현한다.

세 덩어리다.

1. 오차가 있는 신호 합성 — 편광자 오프셋 Ps, 분석기 오프셋 As, 검출기 감쇠 eta
2. 잔차 보정 — R(P) 의 포물선 최소점에서 Ps, 최솟값에서 eta, 그 자리에서 As
3. 회귀 보정 — 포물선 근사를 버리고 Psi, Delta, Ps, As, eta 다섯을 한꺼번에 피팅

신호 합성은 부품 뮬러 행렬의 곱으로만 한다. Fujiwara 식 4.58-4.60 과 Johs 식 2-5 는
verify_calibration.py 가 대조할 독립 경로로 남겨둔다.

부호 규약은 배경이론 2·3편, 타원계측기 1~3편과 일치시킨다.
"""
import numpy as np
from scipy.optimize import least_squares

N_ANGLE = 360          # 분석기 한 바퀴를 몇 점으로 표본화할지


# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소 (배경이론 2편과 동일한 정의)
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])


def mueller_polarizer(theta):
    base = 0.5 * np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


def mueller_sample(psi, delta):
    """등방성 시편. 3편 smm_tensor.sample_mueller_matrix 와 같은 꼴이다."""
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    cd, sd = np.cos(delta), np.sin(delta)
    return np.array([[1, -c2, 0, 0], [-c2, 1, 0, 0],
                     [0, 0, s2 * cd, s2 * sd], [0, 0, -s2 * sd, s2 * cd]])


S_UNPOLARIZED = np.array([1.0, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# 1. 오차가 있는 신호
# ---------------------------------------------------------------------------


def waveform(psi, delta, P, ps=0.0, a_s=0.0, eta=1.0, n=N_ANGLE, pds=0.0,
             noise=0.0, rng=None):
    """분석기 한 바퀴 동안의 검출 세기.

    P     : 편광자 눈금값. 참각은 P - ps 다.
    a_s   : 분석기 0° 위치가 입사면에서 벗어난 양. 눈금 A 의 참각은 A - a_s.
    eta   : 검출기의 ac 감쇠. dc 는 그대로 두고 ac 성분만 eta 배 한다.
    pds   : 검출기의 편광 의존 감도. 투과축 방향에 따라 감도가
            1 + pds*cos(2A) 로 변하는 가장 단순한 꼴로 넣는다(Johs 의 Xc 에 해당).
    noise : 표본마다 더할 가우스 잡음의 표준편차. dc 세기 0.25 기준이다.

    반환: (분석기 눈금각 배열, 세기 배열)
    """
    A = np.linspace(0.0, np.pi, n, endpoint=False)
    mp = mueller_polarizer(P - ps)
    ms = mueller_sample(psi, delta)
    s_after = ms @ (mp @ S_UNPOLARIZED)
    out = np.empty(n)
    for i, a in enumerate(A):
        out[i] = (mueller_polarizer(a - a_s) @ s_after)[0]
    if pds:
        out = out * (1.0 + pds * np.cos(2 * A))
    dc = out.mean()
    out = dc + eta * (out - dc)          # ac 성분만 감쇠시킨다
    if noise:
        out = out + (rng or np.random.default_rng()).normal(0.0, noise, n)
    return A, out


def fourier_ab(signal):
    """I = dc[1 + a*cos(2A) + b*sin(2A)] 의 정규화 계수 (a, b)."""
    n = len(signal)
    k = 2 * np.pi * np.arange(n) / n
    dc = signal.mean()
    return (2 * np.mean(signal * np.cos(k)) / dc,
            2 * np.mean(signal * np.sin(k)) / dc)


def measure(psi, delta, P, **kw):
    """편광자 눈금 P 에서 측정되는 (alpha', beta')."""
    return fourier_ab(waveform(psi, delta, P, **kw)[1])


# ---------------------------------------------------------------------------
# 2. 잔차 보정 (Aspnes)
# ---------------------------------------------------------------------------


def residual(psi, delta, P, **kw):
    """잔차 함수 R(P) = 1 - (alpha'^2 + beta'^2). Fujiwara 식 4.61a."""
    a, b = measure(psi, delta, P, **kw)
    return 1.0 - (a ** 2 + b ** 2)


def residual_calibration(psi, delta, scan_deg=np.arange(-5, 5.01, 0.5), **kw):
    """잔차 보정으로 Ps, eta, As 를 뽑는다.

    절차는 Fujiwara 4.3.3 그대로다.
      Ps  : R(P) 를 포물선으로 맞춘 최소점
      eta : R 의 최솟값에서  R = 1 - eta^2
      As  : P - Ps = 0 에서 (alpha', beta') = eta*(cos 2As, sin 2As)
    반환: dict. 포물선이 아래로 볼록하지 않으면 ok=False 로 돌려준다.
    """
    P = np.deg2rad(scan_deg)
    R = np.array([residual(psi, delta, p, **kw) for p in P])
    c2, c1, c0 = np.polyfit(P, R, 2)
    if c2 <= 0:                                   # 위로 볼록하면 최소점이 없다
        return dict(ok=False, ps=np.nan, eta=np.nan, a_s=np.nan,
                    curvature=c2, P=P, R=R)
    ps = -c1 / (2 * c2)
    r_min = c0 - c1 ** 2 / (4 * c2)
    eta = np.sqrt(max(1.0 - r_min, 0.0))
    a, b = measure(psi, delta, ps, **kw)          # P - Ps = 0 인 자리
    return dict(ok=True, ps=ps, eta=eta, a_s=0.5 * np.arctan2(b, a),
                curvature=c2, P=P, R=R)


# ---------------------------------------------------------------------------
# 3. 회귀 보정 (Johs)
# ---------------------------------------------------------------------------


def model_ab(P, psi, delta, ps, a_s, eta):
    """Johs 식 2-5. 이상적 RAE 응답에 As 회전과 eta 배율을 씌운 꼴이다."""
    t, u = np.tan(psi) ** 2, np.tan(P - ps)
    denom = t + u ** 2
    a0 = (t - u ** 2) / denom
    b0 = 2 * np.tan(psi) * np.cos(delta) * u / denom
    return (eta * (a0 * np.cos(2 * a_s) - b0 * np.sin(2 * a_s)),
            eta * (a0 * np.sin(2 * a_s) + b0 * np.cos(2 * a_s)))


def regression_calibration(P_scan, a_exp, b_exp, guess, sigma=1.0):
    """Psi, Delta, Ps, As, eta 다섯을 한꺼번에 맞춘다 (Levenberg-Marquardt).

    guess: (psi, delta, ps, a_s, eta) 초기값
    반환: dict(params, mse, ok)
    """
    def resid(x):
        a, b = model_ab(P_scan, *x)
        return np.concatenate([(a - a_exp), (b - b_exp)]) / sigma

    sol = least_squares(resid, np.asarray(guess, float), method="lm", max_nfev=20000)
    # MSE 는 분광 타원계측 소프트웨어의 정의를 따른다 — 축소 카이제곱의 제곱근이다.
    # 따라서 모형에 없는 효과의 크기에 제곱이 아니라 선형으로 반응한다.
    n_obs, n_par = 2 * len(P_scan), 5
    chi2 = float(np.sum(sol.fun ** 2))
    return dict(params=sol.x, mse=np.sqrt(chi2 / max(n_obs - n_par, 1)),
                ok=bool(sol.success), chi2=chi2)


# ---------------------------------------------------------------------------
# 4. 데이터 환산 — 보정하지 않은 눈금값으로 Psi, Delta 를 읽으면 어떻게 되는가
# ---------------------------------------------------------------------------


def invert_ideal(a, b, P):
    """이상적 RAE 라고 믿고 (alpha', beta') 에서 Psi, Delta 를 읽는다.

    Fujiwara 식 4.16-4.17 이다. 오차가 있는 측정값에 이 식을 그대로 쓰면
    그 오차가 Psi, Delta 로 얼마나 번지는지 볼 수 있다.

    1편에서 본 대로 arccos 이라 Delta 의 부호는 여기서 잃는다.
    """
    a = np.clip(a, -1.0 + 1e-15, 1.0 - 1e-15)
    psi = np.arctan(np.sqrt((1.0 + a) / (1.0 - a)) * np.abs(np.tan(P)))
    delta = np.arccos(np.clip(b / np.sqrt(1.0 - a ** 2), -1.0, 1.0))
    return psi, delta

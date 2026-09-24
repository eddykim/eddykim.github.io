"""엘립소미터 편광 변조 방식의 검출 세기 I(t)를 뮬러 계산법으로 구현.

RAE / RAE+보상자 / RCE / PME 네 배치를 모두 "편광 부품 뮬러 행렬의 곱"이라는
하나의 경로로 계산한다. 해석식(Fujiwara 식 4.18, 4.28, 4.32, 4.41, 4.46)은
쓰지 않는다 — verify_modulation.py 가 이 수치 결과를 그 해석식과 대조해
부호 규약이 맞는지 검증하는 독립 경로 역할을 해야 하기 때문이다.

부호 규약은 배경이론 2편(generate_figures.py 의 mueller_* 함수)과
3편(smm_tensor.py 의 sample_mueller_matrix)을 그대로 따른다.
"""
import numpy as np

# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소 (배경이론 2편과 동일한 정의)
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    """좌표 회전 행렬. omega: rad."""
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])


def mueller_polarizer(omega_p):
    """방위각 omega_p 에 놓인 이상적 선편광자."""
    base = 0.5 * np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return mueller_rotation(-omega_p) @ base @ mueller_rotation(omega_p)


def mueller_retarder(omega_c, phi):
    """방위각 omega_c, 위상 지연량 phi 인 위상지연자.

    phi > 0 은 존스 행렬 diag(exp(i*phi), 1) 에 대응한다. 즉 시편의
    sample_mueller_matrix(..., delta) 와 완전히 같은 꼴이며, 이 부품을 시편 뒤에
    두면 두 위상이 단순히 더해진다(Delta -> Delta + phi).

    Fujiwara 4.2.3(RAE+보상자)과 4.2.5(PME)는 느린 축을 p 에 두어
    Delta' = Delta - delta 가 되므로 phi = -delta 로 호출하고,
    4.2.4(RCE)는 빠른 축을 p 에 두므로 phi = +delta 로 호출한다.
    """
    c, s = np.cos(phi), np.sin(phi)
    base = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, c, s], [0, 0, -s, c]])
    return mueller_rotation(-omega_c) @ base @ mueller_rotation(omega_c)


def mueller_sample(psi, delta):
    """등방성 시편의 뮬러 행렬(전체 반사율은 1로 정규화).

    학위논문 식 1.2 및 3편 smm_tensor.sample_mueller_matrix 와 동일한 꼴이다.
    """
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    cd, sd = np.cos(delta), np.sin(delta)
    return np.array([
        [1, -c2, 0, 0],
        [-c2, 1, 0, 0],
        [0, 0, s2 * cd, s2 * sd],
        [0, 0, -s2 * sd, s2 * cd],
    ])


S_UNPOLARIZED = np.array([1.0, 0.0, 0.0, 0.0])


def _intensity(chain, s_in=S_UNPOLARIZED):
    """chain 은 빛이 지나는 순서대로 나열한 뮬러 행렬. S0 성분(검출 세기)을 반환."""
    s = s_in
    for m in chain:
        s = m @ s
    return s[0]


# ---------------------------------------------------------------------------
# 네 가지 배치의 검출 세기 파형
# ---------------------------------------------------------------------------


def intensity_rae(psi, delta, angles, pol=np.pi / 4):
    """RAE (PSA_R): 분석기만 회전. angles 는 분석기 방위각 A 배열(rad)."""
    ms, mp = mueller_sample(psi, delta), mueller_polarizer(pol)
    return np.array([_intensity([mp, ms, mueller_polarizer(a)]) for a in angles])


def intensity_rae_compensator(psi, delta, angles, comp_delta, comp_angle=0.0,
                              pol=np.pi / 4):
    """RAE+보상자 (PSCA_R): 고정 보상자 뒤에서 분석기가 회전.

    느린 축을 p 에 두는 규약이므로 위상 지연량의 부호를 뒤집어 넣는다.
    """
    ms, mp = mueller_sample(psi, delta), mueller_polarizer(pol)
    mc = mueller_retarder(comp_angle, -comp_delta)
    return np.array([_intensity([mp, ms, mc, mueller_polarizer(a)]) for a in angles])


def intensity_rce(psi, delta, angles, comp_delta=np.pi / 2, ana=0.0, pol=np.pi / 4):
    """RCE (PSC_R A): 보상자가 회전하고 분석기는 고정. angles 는 보상자 방위각 C."""
    ms, mp, ma = mueller_sample(psi, delta), mueller_polarizer(pol), mueller_polarizer(ana)
    return np.array([_intensity([mp, ms, mueller_retarder(c, comp_delta), ma])
                     for c in angles])


def intensity_pme(psi, delta, retardations, mod_angle=0.0, ana=None,
                  pol=np.pi / 4):
    """PME (PSMA): 위상변조기의 지연량이 시간에 따라 변한다.

    retardations 는 각 시각의 지연량 delta_PEM(t) 배열(rad)이며, 보통
    F*sin(omega*t) 를 넣는다. 부품은 하나도 회전하지 않는다.

    ana 를 주지 않으면 PME 의 표준 조건인 A - M = 45° 로 분석기를 놓는다.
    이 조건이 깨지면 식 4.41 의 단순한 꼴이 성립하지 않는다.
    """
    if ana is None:
        ana = mod_angle + np.pi / 4
    ms, mp, ma = mueller_sample(psi, delta), mueller_polarizer(pol), mueller_polarizer(ana)
    return np.array([_intensity([mp, ms, mueller_retarder(mod_angle, -d), ma])
                     for d in retardations])


def pem_retardation(t, amplitude, omega):
    """광탄성 변조기의 지연량 delta(t) = F sin(omega t) (Fujiwara 식 3.4)."""
    return amplitude * np.sin(omega * t)


# ---------------------------------------------------------------------------
# 푸리에 계수 추출
# ---------------------------------------------------------------------------


def fourier_coefficients(signal, n_max=4):
    """주기 신호를 I(t) = dc + sum_n [a_n cos(n w t) + b_n sin(n w t)] 로 분해한다.

    signal 은 한 주기를 균등 표본화한 배열이다(마지막 점이 첫 점과 겹치면 안 된다).
    반환: dc, {n: (a_n, b_n)} — 정규화하지 않은 날 계수다.

    RAE·RCE 의 alpha, beta 는 dc 로 나눈 값이지만, PME 는 J0(F) 항 때문에
    dc 가 I0 와 다르므로 나누는 기준을 호출부가 직접 고르게 둔다.
    """
    n_pts = len(signal)
    k = np.arange(n_pts)
    dc = signal.mean()
    coeffs = {}
    for n in range(1, n_max + 1):
        phase = 2 * np.pi * n * k / n_pts
        coeffs[n] = (2 * np.mean(signal * np.cos(phase)),
                     2 * np.mean(signal * np.sin(phase)))
    return dc, coeffs


def normalized_coefficients(signal, n_max=4):
    """dc 로 정규화한 푸리에 계수 {n: (alpha_n, beta_n)} (RAE·RCE 의 표준 정의)."""
    dc, coeffs = fourier_coefficients(signal, n_max=n_max)
    return {n: (a / dc, b / dc) for n, (a, b) in coeffs.items()}


def amplitude_spectrum(signal, n_max=6):
    """정규화 푸리에 계수의 크기 sqrt(a^2 + b^2) 를 차수별로 반환(막대 스펙트럼용)."""
    coeffs = normalized_coefficients(signal, n_max=n_max)
    return np.array([np.hypot(*coeffs[n]) for n in range(1, n_max + 1)])

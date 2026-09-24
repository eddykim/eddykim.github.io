"""위상지연자의 파장 의존성과 광대역 설계를 계산한다.

타원계측기 2편에서 쓰는 계산은 전부 여기 있다. 세 덩어리다.

1. 복굴절 파장판의 지연량        delta = 2*pi*d*dn/lambda
2. 여러 장을 겹쳤을 때의 유효 지연량과 유효 빠른 축
3. 전반사가 만드는 p-s 위상차 (프레넬 롬)

2번은 존스 경로와 뮬러 경로 두 가지로 따로 구현했다. verify_retarder.py 가
둘을 대조해야 하므로 한쪽을 다른 쪽에서 유도하지 않는다.

기호와 부호 규약은 배경이론 2·3편, 타원계측기 1편과 일치시킨다.
"""
import numpy as np

# 석영의 복굴절률. Handbook(적색광 기준)과 학위논문(400-700 nm 기준)이
# 독립적으로 같은 값을 주므로 상수로 쓴다. Handbook 이 "n_e - n_o 는 파장
# 의존성이 약하다"고 명시한 근사에 기대는 값이다.
DN_QUARTZ = 0.009


# ---------------------------------------------------------------------------
# 1. 단일 파장판
# ---------------------------------------------------------------------------


def plate_thickness(order, wavelength_design, dn=DN_QUARTZ):
    """설계 파장에서 지연량이 order 파장이 되는 판 두께.

    order = 0.25 면 영차 4분의1파장판, 2.25 면 다중차 4분의1파장판이다.
    delta = 2*pi*order 가 되도록 d = order * lambda / dn 으로 잡는다.
    """
    return order * wavelength_design / dn


def retardance(wavelength, thickness, dn=DN_QUARTZ):
    """파장판의 지연량(rad). Fujiwara 식 3.3 = Handbook 식 4.2."""
    return 2 * np.pi * thickness * dn / wavelength


# ---------------------------------------------------------------------------
# 2-a. 존스 경로
# ---------------------------------------------------------------------------


def jones_retarder(delta, theta):
    """방위각 theta, 지연량 delta 인 선형 위상지연자의 존스 행렬.

    대칭꼴 diag(exp(+i delta/2), exp(-i delta/2)) 를 쓴다. 행렬식이 1 이라
    고유값에서 지연량을 바로 읽을 수 있다(전체 위상이 섞이지 않는다).
    """
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, s], [-s, c]])
    core = np.diag([np.exp(0.5j * delta), np.exp(-0.5j * delta)])
    return rot.T @ core @ rot


def jones_stack(deltas, thetas):
    """여러 장을 겹친 존스 행렬. 리스트 순서가 빛이 지나는 순서다."""
    out = np.eye(2, dtype=complex)
    for d, t in zip(deltas, thetas):
        out = jones_retarder(d, t) @ out
    return out


def effective_retardance_jones(J):
    """존스 행렬에서 유효 지연량(rad)을 뽑는다.

    행렬식이 1 인 2x2 유니터리는 푸앵카레 구면의 어떤 축을 중심으로 한 회전이고,
    고유값이 exp(+-i*delta_eff/2) 이므로 대각합이 2*cos(delta_eff/2) 가 된다.
    """
    Jn = J / np.sqrt(np.linalg.det(J))
    half = np.clip(np.real(np.trace(Jn)) / 2.0, -1.0, 1.0)
    return 2 * np.arccos(half)


def effective_fast_axis_jones(J):
    """유효 빠른 축의 방위각(rad).

    위상이 앞서는 쪽 고유벡터의 편광 방위각을 방위각으로 삼는다.
    반환 범위는 -pi/2 ~ pi/2 (방위각은 180도 주기다).
    """
    Jn = J / np.sqrt(np.linalg.det(J))
    w, v = np.linalg.eig(Jn)
    fast = v[:, np.argmax(np.angle(w))]
    ex, ey = fast
    s1 = abs(ex) ** 2 - abs(ey) ** 2
    s2 = 2 * np.real(np.conj(ex) * ey)
    return 0.5 * np.arctan2(s2, s1)


# ---------------------------------------------------------------------------
# 2-b. 뮬러 경로 (존스에서 유도하지 않는다 — 교차검증용 독립 경로)
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    """좌표 회전 행렬 (배경이론 2편과 동일한 정의)."""
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])


def mueller_retarder(delta, theta):
    """방위각 theta, 지연량 delta 인 선형 위상지연자의 뮬러 행렬.

    배경이론 2편 mueller_retarder 와 같은 부호 규약이다.
    """
    c, s = np.cos(delta), np.sin(delta)
    base = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, c, s], [0, 0, -s, c]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


def mueller_stack(deltas, thetas):
    """여러 장을 겹친 뮬러 행렬. 리스트 순서가 빛이 지나는 순서다."""
    out = np.eye(4)
    for d, t in zip(deltas, thetas):
        out = mueller_retarder(d, t) @ out
    return out


def effective_retardance_mueller(M):
    """뮬러 행렬에서 유효 지연량(rad)을 뽑는다.

    순수 위상지연자의 뮬러 행렬은 우하단 3x3 이 푸앵카레 구면의 회전행렬이고,
    회전각 delta_eff 에 대해 대각합이 1 + 2*cos(delta_eff) 다.
    """
    R = M[1:, 1:]
    return np.arccos(np.clip((np.trace(R) - 1) / 2.0, -1.0, 1.0))


def stokes_to_jones(S):
    """완전편광 스토크스 벡터를 존스 벡터로. 위상차 순서는 phi_y - phi_x 다.

    배경이론 2편의 정의(2026-09-24 정정본)와 같은 규약이다.
    """
    S0, S1, S2, S3 = S
    ex = np.sqrt(max(S0 + S1, 0.0) / 2)
    ey = np.sqrt(max(S0 - S1, 0.0) / 2)
    return np.array([ex, ey * np.exp(1j * np.arctan2(S3, S2))])


def jones_to_stokes(E):
    """존스 벡터를 스토크스 벡터로. S2 = 2Re(Ex* Ey), S3 = 2Im(Ex* Ey)."""
    ex, ey = E
    return np.array([abs(ex) ** 2 + abs(ey) ** 2,
                     abs(ex) ** 2 - abs(ey) ** 2,
                     2 * np.real(np.conj(ex) * ey),
                     2 * np.imag(np.conj(ex) * ey)])


# ---------------------------------------------------------------------------
# 3. 전반사 위상차 (프레넬 롬)
# ---------------------------------------------------------------------------


def tir_phase_difference(theta, n_in, n_out=1.0):
    """매질 내부(n_in)에서 전반사할 때 생기는 p-s 위상차(rad).

    프레넬 반사계수를 복소수로 직접 계산해 편각 차를 취한다.
    theta 는 입사각(rad)이며 임계각보다 커야 한다.
    """
    st = n_in * np.sin(theta) / n_out
    ct = np.sqrt(1 - st ** 2 + 0j)          # 전반사 영역에서는 순허수가 된다
    rs = (n_in * np.cos(theta) - n_out * ct) / (n_in * np.cos(theta) + n_out * ct)
    rp = (n_out * np.cos(theta) - n_in * ct) / (n_out * np.cos(theta) + n_in * ct)
    return np.angle(rp) - np.angle(rs)


def tir_phase_difference_closed(theta, n_in, n_out=1.0):
    """같은 양의 닫힌 식. tan(-d/2) = cos(t) sqrt(sin^2 t - n^2) / sin^2 t.

    n = n_out/n_in 이다. 교과서는 보통 s 기준(phi_s - phi_p)으로 적지만, 이 시리즈는
    Delta = phi_p - phi_s 규약을 쓰므로 부호를 뒤집어 맞춘다. 전반사에서 이 값은 음수다.
    verify_retarder.py 가 프레넬 계수 경로와 이 식을 대조한다.
    """
    n = n_out / n_in
    st = np.sin(theta)
    return -2 * np.arctan2(np.cos(theta) * np.sqrt(np.maximum(st ** 2 - n ** 2, 0.0)),
                           st ** 2)


def critical_angle(n_in, n_out=1.0):
    """임계각(rad)."""
    return np.arcsin(n_out / n_in)

"""위상을 꺼내는 세 가지 길 — 이미지처리 6편.

  힐베르트 변환   한 장으로 푼다. 배경을 먼저 빼고 해석 신호를 만들어 편각을 읽는다.
  푸리에 변환법   한 장으로 푼다. 스펙트럼에서 +1 차 측대역만 골라 원점으로 옮긴다.
                  Takeda 가 1982 년에 제안했고, 3편의 노치 필터와 같은 도구다.
  위상천이        여러 장으로 푼다. 위상을 아는 만큼 옮겨가며 찍어 미지수를 줄인다.

앞의 둘은 위상이 반송파를 타고 천천히 변한다는 공간적 가정을 얹는다. 위상천이는
그 가정이 필요 없는 대신 여러 장을 찍어야 하고, 찍는 동안 시편이 움직이면 안 된다.
"""
import numpy as np
from scipy.signal import hilbert


def wrap(p):
    """위상을 (-pi, pi] 로 접는다."""
    return (p + np.pi) % (2.0 * np.pi) - np.pi


def phase_hilbert(I, detrend_sigma=None):
    """해석 신호의 편각으로 위상을 읽는다.

    힐베르트 변환은 신호의 모든 주파수 성분의 위상을 -pi/2 만큼 돌리는 선형 연산이다.
    원 신호에 그 결과를 허수부로 붙이면 해석 신호(analytic signal)가 되고,
    그 크기가 순간 진폭, 편각이 순간 위상이다.

        z(x) = I_ac(x) + j H[I_ac(x)] = b(x) exp(j phi(x))

    전제가 둘 있다. 배경 a 가 먼저 빠져 있어야 하고(그렇지 않으면 편각이 배경 쪽으로
    끌린다), b 와 cos(phi) 의 스펙트럼이 겹치지 않아야 한다(Bedrosian 조건).
    """
    from scipy.ndimage import gaussian_filter1d
    ac = I - (gaussian_filter1d(I, detrend_sigma) if detrend_sigma else I.mean())
    z = hilbert(ac)
    return np.angle(z), np.abs(z)


def phase_fourier(I, band=None):
    """스펙트럼에서 +1 차 측대역만 골라 원점으로 옮긴 뒤 역변환한다."""
    n = len(I)
    F = np.fft.fft(I - I.mean())
    freq = np.fft.fftfreq(n)
    if band is None:
        half = freq > 0
        peak = np.argmax(np.abs(F) * half)
        f0 = freq[peak]
        width = 0.6 * f0
    else:
        f0, width = band
    mask = np.abs(freq - f0) < width
    z = np.fft.ifft(F * mask)
    # 반송파를 빼지 않고 그대로 두면 아래 언래핑에서 처리된다.
    return np.angle(z), 2.0 * np.abs(z)


def phase_shift_n(frames, shifts):
    """일반 N 단계 위상천이. 최소자승 해라 천이량이 고르지 않아도 된다.

        I_k = a + b cos(phi + d_k) = a + (b cos phi) cos d_k - (b sin phi) sin d_k

    미지수 (a, b cos phi, b sin phi) 에 대해 선형이므로 세 장이면 풀린다.
    배경 a 가 미지수로 함께 풀려 나가는 것이 이 방법의 핵심이다.
    """
    frames = np.asarray(frames, float)
    d = np.asarray(shifts, float)
    A = np.column_stack([np.ones_like(d), np.cos(d), -np.sin(d)])
    sol, *_ = np.linalg.lstsq(A, frames, rcond=None)
    a, c, s = sol
    return np.arctan2(s, c), np.hypot(c, s), a


def phase_shift_5_hariharan(frames):
    """Schwider-Hariharan 5 단계. 천이량 오차에 1 차로 둔감하다.

        phi = atan2( 2(I4 - I2), I1 - 2 I3 + I5 )

    천이량이 pi/2 에서 벗어나도 오차의 1 차항이 상쇄되도록 계수를 고른 것이다.
    부호를 이렇게 잡아야 결과가 phi 그대로 나온다. 흔히 보이는
    atan2(2(I2-I4), 2 I3 - I5 - I1) 형태는 같은 값에 pi 가 더해진 것이다.
    """
    I1, I2, I3, I4, I5 = frames
    return np.arctan2(2.0 * (I4 - I2), I1 - 2.0 * I3 + I5)


def phase_error(estimated, reference):
    """위상 추정 오차. 절대 위상은 관측량이 아니므로 상수 차이를 먼저 없앤다.

    단순히 평균을 빼면 안 된다. 차이가 pi 근처면 접힘 때문에 +pi 와 -pi 로
    갈라져 분산이 터진다. 복소 평면에서 평균 방향(원형 평균)을 구해 빼야 한다.
    """
    d = np.exp(1j * (np.asarray(estimated) - np.asarray(reference)))
    return wrap(np.angle(d) - np.angle(d.mean()))


def unwrap_itoh(wrapped):
    """1 차원 언래핑. 이웃 간 위상차를 접어서 누적한다.

    이웃한 두 화소의 참 위상차가 pi 보다 작다는 가정 하나로 성립한다. 표본이
    성기거나 노이즈가 크면 이 가정이 깨지고, 한 번 깨지면 그 뒤로 2pi 만큼
    어긋난 채 끝까지 간다.
    """
    d = wrap(np.diff(wrapped))
    return np.concatenate([[wrapped[0]], wrapped[0] + np.cumsum(d)])


def residues(wrapped2d):
    """2 차원 위상장의 잔여점(residue).

    2x2 고리를 한 바퀴 돌며 접은 위상차를 더한다. 0 이 아니면 그 고리 안에서
    언래핑 결과가 경로에 따라 달라진다. 잔여점이 있으면 단순 적분은 실패한다.
    """
    p = wrapped2d
    d1 = wrap(p[:-1, 1:] - p[:-1, :-1])     # 위쪽 변, 오른쪽으로
    d2 = wrap(p[1:, 1:] - p[:-1, 1:])       # 오른쪽 변, 아래로
    d3 = wrap(p[1:, :-1] - p[1:, 1:])       # 아래쪽 변, 왼쪽으로
    d4 = wrap(p[:-1, :-1] - p[1:, :-1])     # 왼쪽 변, 위로
    return np.round((d1 + d2 + d3 + d4) / (2.0 * np.pi)).astype(int)

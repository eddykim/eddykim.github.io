"""2차원 주파수 영역 도구 — 스펙트럼 보기, 노치 제거, 저역통과 비교.

공간 영역에서 손댈 수 없는 것이 두 가지 있다. 하나는 주기성 배경 노이즈이고,
다른 하나는 필터의 통과 대역을 정확히 지정하는 일이다. 둘 다 주파수 영역에서는
간단해지지만, 대신 푸리에 변환이 깔고 있는 가정 — 신호가 주기적으로 되풀이된다는
가정 — 을 함께 떠안게 된다.
"""
import numpy as np


def spectrum(img, window=None):
    """중심을 0 주파수로 옮긴 진폭 스펙트럼을 로그 눈금으로 돌려준다."""
    a = img * window if window is not None else img
    return np.log1p(np.abs(np.fft.fftshift(np.fft.fft2(a))))


def hann2d(shape):
    """2차원 한(Hann) 창. 가장자리를 0 으로 눌러 주기 확장의 불연속을 없앤다."""
    wy = np.hanning(shape[0])
    wx = np.hanning(shape[1])
    return np.outer(wy, wx)


def freq_grid(shape):
    """중심이 0 인 주파수 격자. 나이키스트가 0.5 가 되도록 정규화한다."""
    fy = np.fft.fftshift(np.fft.fftfreq(shape[0]))
    fx = np.fft.fftshift(np.fft.fftfreq(shape[1]))
    return np.meshgrid(fx, fy)


def notch_filter(img, peaks, radius):
    """스펙트럼에서 지정한 좌표 주변을 지운다. 주기성 노이즈 제거용이다.

    peaks 는 (fy, fx) 정규화 주파수 쌍의 목록이다. 각 점과 그 대칭점을 함께 막는다.
    실수 영상의 스펙트럼은 원점 대칭이므로 한쪽만 지우면 허수부가 남는다.
    """
    F = np.fft.fftshift(np.fft.fft2(img))
    FX, FY = freq_grid(img.shape)
    mask = np.ones(img.shape, float)
    for fy, fx in peaks:
        for sy, sx in ((fy, fx), (-fy, -fx)):
            mask *= 1.0 - np.exp(-0.5 * (((FX - sx) ** 2 + (FY - sy) ** 2) / radius ** 2))
    return np.real(np.fft.ifft2(np.fft.ifftshift(F * mask))), mask


def lowpass(img, cutoff, kind="ideal"):
    """원형 저역통과. kind 가 ideal 이면 칼로 자르고, gaussian 이면 부드럽게 민다."""
    F = np.fft.fftshift(np.fft.fft2(img))
    FX, FY = freq_grid(img.shape)
    R = np.hypot(FX, FY)
    if kind == "ideal":
        H = (R <= cutoff).astype(float)
    elif kind == "gaussian":
        H = np.exp(-0.5 * (R / cutoff) ** 2)
    else:
        raise ValueError(kind)
    return np.real(np.fft.ifft2(np.fft.ifftshift(F * H))), H


def impulse_response(H):
    """필터의 공간 영역 모습. 이상적 필터가 왜 링을 만드는지는 여기에 다 보인다."""
    return np.real(np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(H))))

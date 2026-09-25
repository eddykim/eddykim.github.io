"""기저 신호를 뽑는 두 가지 방법 — 푸리에 저역통과와 극값 포락선.

학위논문 3.3절이 비교하는 바로 그 두 방법이다.

푸리에 저역통과(FT+LPF)는 파수축 푸리에 변환에서 변조 주파수 위쪽을 잘라낸다.
구현이 간단한 대신 두 가지를 대가로 치른다. 첫째, 기저 신호 자체에 들어 있던
차단주파수 위쪽 성분까지 함께 지운다. 둘째, 푸리에 변환이 신호를 주기 함수로
보기 때문에 양 끝값이 다르면 그 불연속이 신호 전체에 깁스(Gibbs) 진동을 퍼뜨린다.

극값 포락선(Extrema Envelope Method, EEM)은 변조 신호가 기저 신호를 중심으로
대칭으로 진동한다는 성질만 쓴다. 극대점과 극소점에서

    I_max = I_BLS + A,    I_min = I_BLS - A

이므로 두 포락선의 평균이 곧 기저 신호다. 주파수 영역을 거치지 않으므로
주기성 가정도, 차단주파수 선택도 필요 없다.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d


def fourier_lpf_baseline(signal, cutoff_bins):
    """파수축 푸리에 변환에서 cutoff_bins 위쪽을 0 으로 만들고 되돌린다.

    cutoff_bins 는 신호 전체 길이에 대한 주기 수로 센다. 변조 주파수(줄무늬 개수)
    보다 작게 잡아야 변조 신호가 지워진다.
    """
    spec = np.fft.rfft(signal)
    spec[cutoff_bins + 1:] = 0.0
    return np.fft.irfft(spec, n=len(signal))


def _refine_extremum(x, y, i, half=3):
    """i 부근에 2차 다항식을 맞춰 극값의 위치와 세기를 화소 이하로 정밀화한다.

    논문이 "국소 다항식 적합"이라 부르는 단계다. 표본점 하나를 그대로 쓰면 극값이
    화소 격자에 갇히고, 그 오차가 포락선을 타고 기저 신호로 옮겨간다.
    """
    lo, hi = max(i - half, 0), min(i + half + 1, len(x))
    xs, ys = x[lo:hi], y[lo:hi]
    x0 = xs.mean()
    a, b, c = np.polyfit(xs - x0, ys, 2)
    if a == 0.0:
        return x[i], y[i]
    xv = -b / (2.0 * a)
    if not (xs[0] - x0 <= xv <= xs[-1] - x0):
        return x[i], y[i]      # 꼭짓점이 창 밖이면 표본점을 그대로 쓴다
    return x0 + xv, a * xv ** 2 + b * xv + c


def find_extrema(x, y, half=3, detect_sigma=0.0):
    """수치 미분의 영점 교차로 극대·극소를 찾고 각각 다항식으로 정밀화한다.

    노이즈가 있으면 미분의 부호가 수없이 뒤집혀 가짜 극값이 쏟아진다. 그래서
    **검출에만** 평활화한 신호를 쓰고, 극값의 위치와 세기는 원본에서 구한다.
    평활화 폭은 줄무늬 한 주기보다 훨씬 작게 잡아야 극값 자체가 눌리지 않는다.

    극대와 극소는 반드시 번갈아 나온다. 그렇지 않은 구간은 노이즈가 만든
    가짜이므로, 연속된 같은 종류 중 가장 극단인 것만 남긴다.
    """
    yd = gaussian_filter1d(y, detect_sigma) if detect_sigma > 0 else y
    d = np.gradient(yd)
    idx = np.flatnonzero(np.signbit(d[:-1]) != np.signbit(d[1:]))
    idx = idx[(idx >= half) & (idx < len(x) - half - 1)]
    if len(idx) == 0:
        return np.empty((0, 2)), np.empty((0, 2))

    kinds = np.where(d[idx] > 0, 1, -1)          # 1 = 극대, -1 = 극소
    keep = []
    for i, k in zip(idx, kinds):
        if keep and keep[-1][1] == k:
            j = keep[-1][0]
            if (y[i] > y[j]) if k == 1 else (y[i] < y[j]):
                keep[-1] = (i, k)        # 같은 종류가 잇따르면 더 극단인 쪽만
        else:
            keep.append((i, k))

    maxima, minima = [], []
    for i, k in keep:
        (maxima if k == 1 else minima).append(_refine_extremum(x, y, i, half))
    maxima, minima = np.array(maxima), np.array(minima)
    # 정밀화가 순서를 흐트러뜨렸다면 스플라인이 세울 수 없다. 단조롭게 정리한다.
    return _monotonic(maxima), _monotonic(minima)


def _monotonic(pts):
    """x 가 증가하도록 정렬하고 중복을 없앤다."""
    if len(pts) == 0:
        return pts
    pts = pts[np.argsort(pts[:, 0])]
    keep = np.concatenate([[True], np.diff(pts[:, 0]) > 1e-12])
    return pts[keep]


def eem_baseline(x, y, half=3, detect_sigma=0.0):
    """극값 포락선 방법. 상·하 포락선을 스플라인으로 잇고 평균한다.

    포락선을 외삽하지 않는다. 두 포락선이 모두 정의되는 구간만 유효 구간으로
    돌려주고, 그 바깥은 NaN 으로 남긴다. 없는 값을 지어내지 않는 편이 낫다.
    """
    maxima, minima = find_extrema(x, y, half, detect_sigma)
    up = CubicSpline(maxima[:, 0], maxima[:, 1])
    lo = CubicSpline(minima[:, 0], minima[:, 1])
    valid = (x >= max(maxima[0, 0], minima[0, 0])) & (x <= min(maxima[-1, 0], minima[-1, 0]))
    out = np.full_like(x, np.nan)
    out[valid] = 0.5 * (up(x[valid]) + lo(x[valid]))
    env = (np.where(valid, up(x), np.nan), np.where(valid, lo(x), np.nan))
    return out, env, maxima, minima

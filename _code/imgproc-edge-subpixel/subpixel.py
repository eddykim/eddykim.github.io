"""서브픽셀 엣지 추정 세 계열 — 이미지처리 4편.

화소 격자보다 잘게 위치를 재는 방법은 크게 셋으로 나뉜다.

  보간(interpolation)  기울기의 봉우리에 포물선을 맞춰 꼭짓점을 취한다. 가장 싸고
                       가장 널리 쓰이지만, 봉우리 모양을 포물선이라고 가정한다.
  모멘트(moment)       창 안의 밝기 모멘트 세 개로 계단의 위치를 역산한다.
                       Tabatabai 와 Mitchell 이 1984년에 제안했다.
  면적(partial area)   픽셀 값이 그 픽셀이 받은 빛의 평균이라는 사실만 쓴다.
                       흐림이 대칭이기만 하면 모양을 가정하지 않는다.

세 방법 모두 "엣지 하나가 창 안에 있다"는 것 외에 서로 다른 가정을 얹는다.
그 가정이 깨지는 방식이 곧 각 방법의 편향이다.
"""
import numpy as np

from edge_model import SUBSAMPLES  # noqa: F401  (문서용)


def _gradient(profile):
    """중심차분. 2편에서 쓴 것과 같은 커널이다."""
    return np.gradient(profile)


def parabola_gradient(profile):
    """기울기 크기의 최대점에 포물선을 맞춰 꼭짓점을 돌려준다."""
    g = np.abs(_gradient(profile))
    i = int(np.argmax(g[1:-1])) + 1
    a, b, c = g[i - 1], g[i], g[i + 1]
    denom = a - 2.0 * b + c
    if denom == 0.0:
        return float(i)
    return i + 0.5 * (a - c) / denom


def moment_tabatabai(profile, window=None):
    """밝기 모멘트 세 개로 계단 위치를 역산한다.

    창 안에 이상적인 계단 하나가 있다고 보고, 낮은 쪽이 차지하는 비율 p 를 푼다.
    p 를 창 길이에 곱하면 엣지 위치가 나온다. 흐림이 없을 때 정확하다.
    """
    v = profile if window is None else profile[window[0]:window[1]]
    n = len(v)
    m1 = v.mean()
    m2 = (v ** 2).mean()
    m3 = (v ** 3).mean()
    var = m2 - m1 ** 2
    if var <= 0:
        return np.nan
    sd = np.sqrt(var)
    skew = (m3 - 3.0 * m1 * m2 + 2.0 * m1 ** 3) / sd ** 3
    p = 0.5 * (1.0 + skew * np.sqrt(1.0 / (4.0 + skew ** 2)))
    start = 0 if window is None else window[0]
    return start - 0.5 + p * n


def transition_window(profile, rel=0.05, margin=3):
    """기울기가 최댓값의 rel 배 위로 올라오는 구간에 여유를 붙여 창을 정한다.

    면적법의 산포는 창 길이의 제곱근에 비례해 커진다. 평탄부까지 통째로 넣으면
    엣지와 무관한 화소의 노이즈만 쌓이므로, 전이 구간에 바짝 맞춘 창이 필요하다.
    반대로 창이 너무 짧으면 평탄부를 잘못 잡아 큰 편향이 생긴다.
    """
    g = np.abs(_gradient(profile))
    above = np.flatnonzero(g >= rel * g.max())
    lo = max(int(above[0]) - margin, 0)
    hi = min(int(above[-1]) + margin + 1, len(profile))
    return lo, hi


def partial_area(profile, plateau=3, window="auto"):
    """픽셀 값이 받은 빛의 평균이라는 사실만 쓴다.

    낮은 쪽 평탄부 B 와 높은 쪽 평탄부 A 를 창의 양 끝에서 잡고, 정규화한 단면
    (A - v)/(A - B) 를 창 전체에 걸쳐 더한다. 그 합이 창 왼쪽 끝에서 엣지까지의
    거리가 된다. 흐림이 엣지에 대해 대칭이면 흐림 폭과 무관하게 성립한다.

    window 가 "auto" 면 전이 구간에 맞춰 창을 자른다. None 이면 단면 전체를 쓴다.
    """
    v = np.asarray(profile, float)
    if window == "auto":
        lo_i, hi_i = transition_window(v)
    elif window is None:
        lo_i, hi_i = 0, len(v)
    else:
        lo_i, hi_i = window
    seg = v[lo_i:hi_i]
    lo = seg[:plateau].mean()
    hi = seg[-plateau:].mean()
    if hi == lo:
        return np.nan
    if hi < lo:                       # 내려가는 엣지면 뒤집어서 같은 식을 쓴다
        seg, lo, hi = -seg, -hi, -lo
    return lo_i - 0.5 + np.sum((hi - seg) / (hi - lo))


ESTIMATORS = {
    "parabola": parabola_gradient,
    "moment": moment_tabatabai,
    "area": partial_area,
}

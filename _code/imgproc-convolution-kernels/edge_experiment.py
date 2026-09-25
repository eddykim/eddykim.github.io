"""엣지 위치 추정의 정밀도와 엣지 폭 — 이미지처리 2편의 핵심 실험.

평활화를 세게 걸수록 엣지 위치는 안정되지만 엣지 자체가 뭉개진다. 이 교환비를
숫자로 적는 것이 목표다. 정답 위치를 알고 있는 합성 엣지를 만들고, 평활화 폭을
바꿔가며 추정 위치의 표준편차와 엣지 폭을 함께 잰다.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.special import erf

from kernels import CENTRAL_DIFF


def make_edge(n, x0, blur_sigma, contrast, offset):
    """광학계 흐림을 거친 계단 엣지. 계단과 가우시안 PSF 의 합성곱은 erf 가 된다."""
    x = np.arange(n)
    return offset + contrast * 0.5 * (1.0 + erf((x - x0) / (np.sqrt(2.0) * blur_sigma)))


def estimate_edge(profile, smooth_sigma):
    """평활화 후 기울기 최대점을 포물선 꼭짓점으로 보간해 서브픽셀 위치를 얻는다."""
    p = gaussian_filter1d(profile, smooth_sigma) if smooth_sigma > 0 else profile
    g = np.abs(np.convolve(p, CENTRAL_DIFF[::-1], mode="same"))
    i = int(np.argmax(g[2:-2])) + 2
    a, b, c = g[i - 1], g[i], g[i + 1]
    denom = a - 2.0 * b + c
    if denom == 0.0:
        return float(i)
    return i + 0.5 * (a - c) / denom      # 포물선 꼭짓점


def edge_width_10_90(profile, smooth_sigma):
    """10%에서 90%까지 올라오는 데 걸리는 거리. 평활화가 늘리는 것이 이 값이다."""
    p = gaussian_filter1d(profile, smooth_sigma) if smooth_sigma > 0 else profile
    lo, hi = p.min(), p.max()
    t10, t90 = lo + 0.1 * (hi - lo), lo + 0.9 * (hi - lo)
    x = np.arange(len(p))
    return float(np.interp(t90, p, x) - np.interp(t10, p, x))


def sweep(smooth_sigmas, n_trial, rng, n=128, x0=63.37, blur_sigma=1.2,
          contrast=1200.0, noise_dn=17.2):
    """평활화 폭을 훑으며 위치 편향·산포와 엣지 폭을 잰다.

    기본값의 출처는 1편이다. 대비 1200 DN 은 1편 5절이 쓴 신호 세기이고,
    노이즈 17.2 DN 은 그 조건에서 한 장을 찍었을 때의 산포다.
    """
    clean = make_edge(n, x0, blur_sigma, contrast, offset=100.0)
    bias, std, width = [], [], []
    for s in smooth_sigmas:
        est = np.empty(n_trial)
        for t in range(n_trial):
            noisy = clean + noise_dn * rng.standard_normal(n)
            est[t] = estimate_edge(noisy, s)
        bias.append(est.mean() - x0)
        std.append(est.std())
        width.append(edge_width_10_90(clean, s))
    return np.array(bias), np.array(std), np.array(width)

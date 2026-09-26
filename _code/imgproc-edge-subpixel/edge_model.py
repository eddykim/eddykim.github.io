"""엣지 단면 모델 — 이미지처리 4편.

실제 센서가 담아내는 엣지는 두 번 뭉개진 뒤에 표본화된다. 먼저 광학계의 점확산함수
(PSF)가 계단을 흐리고, 그다음 픽셀이 자기 면적 안의 빛을 적분한다. 서브픽셀 추정의
정확도는 이 두 폭과 픽셀 간격의 비가 정한다.

계단을 폭 w 의 가우시안으로 흐리면 단면은 오차함수가 된다.

    I(x) = B + (A - B) * 0.5 * (1 + erf((x - x0) / (sqrt(2) w)))

픽셀 적분까지 넣으려면 이 곡선을 픽셀 폭만큼 평균해야 하는데, 가우시안 흐림과
박스 적분의 합성에는 닫힌 형태가 없다. 그래서 픽셀 안을 여러 점으로 나눠 수치
적분한다. 1편 6절에서 픽셀이 점이 아니라 면이라고 한 것을 여기서 실제로 쓴다.
"""
import numpy as np
from scipy.special import erf

SUBSAMPLES = 32          # 픽셀 하나를 몇 조각으로 나눠 적분할 것인가


def blurred_step(x, x0, blur_w, low=100.0, high=1300.0):
    """가우시안 PSF 를 지난 계단. blur_w 가 0 이면 이상적인 계단이 된다."""
    if blur_w <= 0:
        return np.where(x >= x0, high, low)
    return low + (high - low) * 0.5 * (1.0 + erf((x - x0) / (np.sqrt(2.0) * blur_w)))


def sample_edge(n, x0, blur_w, low=100.0, high=1300.0, pixel_integrate=True):
    """픽셀 중심이 0, 1, ..., n-1 인 센서가 담아낸 엣지 단면.

    pixel_integrate 가 참이면 픽셀 개구의 적분까지 넣는다. 거짓이면 픽셀을 점으로
    보고 중심 값만 취한다. 둘의 차이가 4절에서 다루는 편향의 한 원인이다.
    """
    centers = np.arange(n, dtype=float)
    if not pixel_integrate:
        return blurred_step(centers, x0, blur_w, low, high)
    offsets = (np.arange(SUBSAMPLES) + 0.5) / SUBSAMPLES - 0.5   # -0.5 ~ +0.5
    grid = centers[:, None] + offsets[None, :]
    return blurred_step(grid, x0, blur_w, low, high).mean(axis=1)


def add_noise(profile, sigma_dn, rng):
    """1편의 노이즈 모델을 단순화해 가우시안으로 얹는다."""
    return profile + sigma_dn * rng.standard_normal(profile.shape)


def edge_width_10_90(blur_w):
    """가우시안 흐림 폭 w 에 대응하는 10-90% 상승 거리. 약 2.563 w 다."""
    return 2.5631 * blur_w

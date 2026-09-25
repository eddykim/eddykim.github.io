"""컨볼루션 커널과 그 성질 — 이미지처리 2편.

1편에서 픽셀 값의 분산을 정했다. 이 파일은 그 분산이 커널을 지나면 어떻게 바뀌는지,
그리고 미분 커널이 왜 노이즈를 키우는지를 직접 계산할 수 있게 해 둔 것이다.

핵심은 한 줄로 요약된다. 입력이 화소마다 독립이고 분산이 sigma^2 이면, 가중치 w 로
컨볼루션한 출력의 분산은 sigma^2 * sum(w^2) 이다. 평활화 커널은 sum(w)=1 이면서
sum(w^2)<1 이라 분산을 줄이고, 미분 커널은 sum(w)=0 이지만 sum(w^2)>0 이라
신호는 지우고 노이즈만 남긴다.
"""
import numpy as np


def box_kernel(k):
    """폭 k 의 박스(이동평균) 커널. 가장 싸지만 주파수 응답이 나쁘다."""
    return np.ones(k) / k


def gaussian_kernel(sigma, truncate=4.0):
    """표준편차 sigma 의 가우시안 커널. 꼬리를 truncate*sigma 에서 자른다."""
    r = int(np.ceil(truncate * sigma))
    x = np.arange(-r, r + 1)
    w = np.exp(-0.5 * (x / sigma) ** 2)
    return w / w.sum()


# 미분 커널 세 가지. 전부 합이 0 이라 평탄한 영역에서는 0 을 내놓는다.
CENTRAL_DIFF = np.array([-0.5, 0.0, 0.5])
SOBEL_1D = np.array([-1.0, 0.0, 1.0]) / 2.0        # 미분 방향 성분
SOBEL_SMOOTH = np.array([1.0, 2.0, 1.0]) / 4.0     # 직교 방향 평활 성분
SCHARR_1D = np.array([-1.0, 0.0, 1.0]) / 2.0
SCHARR_SMOOTH = np.array([3.0, 10.0, 3.0]) / 16.0  # Scharr 가 고친 것은 이 쪽이다


def noise_gain(w):
    """백색 노이즈에 대한 분산 이득. 출력 분산 = 입력 분산 * 이 값."""
    return float(np.sum(np.asarray(w, float) ** 2))


def dc_gain(w):
    """평탄한 입력에 대한 이득. 평활화 커널은 1, 미분 커널은 0 이어야 한다."""
    return float(np.sum(np.asarray(w, float)))


def frequency_response(w, n_freq=512):
    """커널의 주파수 응답. 나이키스트를 1 로 정규화한 주파수축과 함께 돌려준다.

    커널이 대칭이므로 중심을 원점에 두고 계산해야 응답이 실수로 나온다.
    """
    w = np.asarray(w, float)
    c = (len(w) - 1) / 2.0
    x = np.arange(len(w)) - c
    f = np.linspace(0, 1, n_freq)              # 1 = 나이키스트
    # exp(-j*pi*f*x) 를 모아 더한다. f=1 이 화소당 반주기에 해당한다.
    # 복소수 matmul 은 일부 BLAS 에서 헛된 부동소수점 경고를 내므로 einsum 을 쓴다.
    H = np.einsum("fx,x->f", np.exp(-1j * np.pi * np.outer(f, x)), w)
    return f, H


def separable_cost(k):
    """2차원 커널을 통째로 돌릴 때와 1차원 두 번으로 나눌 때의 화소당 곱셈 수."""
    return k * k, 2 * k


def convolve1d_edge(sig, w, mode):
    """경계 처리 방식만 바꿔가며 1차원 컨볼루션을 수행한다.

    mode 는 numpy.pad 의 이름을 그대로 쓴다.
      zero      -> 바깥을 0 으로 (constant). 밝은 영상의 가장자리가 어두워진다.
      replicate -> 끝 값을 늘림 (edge). 무난하지만 기울기를 0 으로 만든다.
      reflect   -> 거울 반사. 기울기까지 이어지므로 계측에서 기본값으로 쓸 만하다.
    """
    pad = len(w) // 2
    npmode = {"zero": "constant", "replicate": "edge", "reflect": "reflect"}[mode]
    padded = np.pad(sig, pad, mode=npmode)
    return np.convolve(padded, w[::-1], mode="valid")

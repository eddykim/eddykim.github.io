"""전방 모델 — 계산 이미징과 타이코그래피 1편.

카메라가 노출 동안 수평으로 움직이면 한 점이 선분으로 번진다. 셔터를 계속 열어 두면
번짐 커널은 상자(box)이고, 셔터를 정해진 부호대로 여닫으면(flutter shutter) 커널은
그 부호가 된다. 어느 쪽이든 측정은 미지 신호에 대한 선형 합성곱이다.

    y = A x + n

A 의 각 열은 한 화소가 센서 위에 남기는 자국이다. 여기서는 합성곱의 가장자리를
잘라내지 않는 선형 합성곱을 쓴다. 번진 신호는 원 신호보다 k - 1 화소 길어진다.
"""
import warnings

import numpy as np
from scipy.linalg import toeplitz

# 이 머신(macOS Accelerate BLAS)의 행렬곱은 결과가 멀쩡해도 divide/overflow 경고를 낸다.
# 계산 결과는 verify_regularization.py 에서 독립 경로와 대조해 확인한다.
warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)

# Raskar, Agrawal, Tumblin (2006) 이 찾은 52칩 부호. 열린 칩이 26개라 빛은 상자의 절반이다.
RASKAR_CODE = np.array([int(c) for c in "1010000111000001010000110011110111010111001001100111"], float)
K = RASKAR_CODE.size


def box_kernel(k=K):
    """셔터를 k 칩 동안 계속 연 번짐. 칩 하나의 노출을 1/k 로 두어 합이 1 이다."""
    return np.ones(k) / k


def coded_kernel(code=RASKAR_CODE):
    """부호대로 여닫은 번짐. 칩 노출은 상자와 같은 1/k 라서 합은 (열린 칩 수)/k 이다."""
    return code / code.size


def smear_matrix(kernel, n):
    """길이 n 신호에 kernel 을 선형 합성곱하는 (n + k - 1) x n 행렬."""
    k = kernel.size
    return toeplitz(np.r_[kernel, np.zeros(n - 1)], np.r_[kernel[0], np.zeros(n - 1)])


def circulant_matrix(kernel, n):
    """주기 경계를 가정한 n x n 순환 합성곱 행렬. 고유벡터가 푸리에 기저다."""
    h = np.zeros(n)
    h[: kernel.size] = kernel
    return np.array([np.roll(h, i) for i in range(n)]).T


def transfer_function(kernel, n_fft=1024):
    """커널의 이산 푸리에 변환 크기. 0 에 닿는 주파수의 정보는 측정에서 사라진다."""
    return np.fft.rfftfreq(n_fft), np.abs(np.fft.rfft(kernel, n_fft))


def test_signal(n=300):
    """계단 몇 개와 매끈한 구간을 섞은 1D 신호. 바코드를 한 줄로 읽은 것쯤으로 보면 된다."""
    x = np.zeros(n)
    x[40:90] = 1.0
    x[120:130] = 0.6
    x[160:260] = 0.5 + 0.4 * np.sin(np.linspace(0, 3 * np.pi, 100))
    x[270:280] = 0.9
    return x


def measure(A, x, sigma, rng=None):
    """신호와 무관한 가우시안 잡음(표준편차 sigma)을 더한 측정."""
    rng = np.random.default_rng(rng)
    return A @ x + sigma * rng.normal(size=A.shape[0])


def measure_poisson(A, x, photons, rng=None):
    """샷 잡음 모형. 칩 노출 1/k 동안 밝기 1 인 화소가 평균 photons 개의 광자를 보낸다.

    상자 커널은 한 화소당 k 칩 전부의 빛을, 부호 커널은 열린 칩의 빛만 받는다.
    반환값은 광자 수를 다시 밝기 단위로 나눈 것이라 A x 와 같은 척도다.
    """
    rng = np.random.default_rng(rng)
    counts = rng.poisson(np.clip(A @ x, 0, None) * photons * K)
    return counts / (photons * K)


def blur_image_rows(img, kernel):
    """영상의 각 행을 kernel 로 선형 합성곱한다 (수평 모션 블러). 폭이 k - 1 늘어난다."""
    A = smear_matrix(kernel, img.shape[1])
    return img @ A.T, A

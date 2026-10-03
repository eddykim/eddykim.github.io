"""세기만 재고 위상을 되찾는 반복법 — 계산 이미징과 타이코그래피 2편.

측정은 원시야 회절 세기 |F{o}|^2 뿐이다. 알고 있는 것은 두 가지다.

  모듈러스 제약 (M)   푸리에 공간에서 크기는 측정값 m 이어야 한다.
  지지 제약 (S)       실공간에서 물체는 정해진 영역(지지) 안에만 있다. 실수 물체라면 양수이기도 하다.

두 제약을 만족하는 점, 곧 두 집합의 교점을 찾는 것이 위상 복원이다. 아래 함수들은
모두 그 교점을 찾는 방법이고, 차이는 두 사영을 어떻게 조합하느냐뿐이다.

  ER    x <- Ps Pm x                       번갈아 사영. 오차가 늘지 않지만 정체한다.
  HIO   지지 안: Pm x / 밖: x - beta Pm x   입력을 다음 출력을 미는 힘으로 본다 (Fienup 1982).
  RAAR  beta/2 (Rs Rm + I) x + (1 - beta) Pm x,  R = 2P - I   (Luke 2005)

FFT 는 norm="ortho" 로 써서 푸리에 변환이 유니터리하다. 그래야 Pm 이 실공간에서도
가장 가까운 점으로의 사영이 된다.
"""
import warnings

import numpy as np
from skimage import data, transform

# macOS Accelerate BLAS 의 행렬곱이 결과는 멀쩡한데 경고를 낸다 (1편 참고).
warnings.filterwarnings("ignore", message=".*encountered in matmul", category=RuntimeWarning)


def F(x):
    return np.fft.fft2(x, norm="ortho")


def iF(X):
    return np.fft.ifft2(X, norm="ortho")


def _gray(rgb):
    rgb = rgb.astype(float)
    return 0.2125 * rgb[..., 0] + 0.7154 * rgb[..., 1] + 0.0721 * rgb[..., 2]


def _resize(img, s):
    s0 = min(img.shape[:2])
    img = transform.resize(img[:s0, :s0], (s, s), anti_aliasing=True)
    return (img - img.min()) / (img.max() - img.min())


def make_object(n=128, s=48, kind="real"):
    """n x n 상자 가운데에 s x s 물체. kind="real" 은 양수 실수, "complex" 는 진폭·위상이 다른 영상.

    반환: (물체, 지지 마스크). 상자 대 지지의 면적비 n^2 / s^2 가 오버샘플링 비 sigma 다.
    """
    amp = _resize(data.camera().astype(float), s)
    obj = np.zeros((n, n), complex)
    sup = np.zeros((n, n), bool)
    c = n // 2 - s // 2
    if kind == "real":
        obj[c:c + s, c:c + s] = amp
    else:
        ph = _resize(_gray(data.astronaut()), s)
        obj[c:c + s, c:c + s] = (0.5 + 0.5 * amp) * np.exp(1j * 0.8 * np.pi * ph)
    sup[c:c + s, c:c + s] = True
    return obj, sup


def P_m(x, mag):
    """모듈러스 사영: 푸리에 크기를 측정값으로 바꾸고 위상은 그대로 둔다."""
    X = F(x)
    return iF(mag * np.exp(1j * np.angle(X)))


def P_s(x, sup, real=False):
    """지지 사영: 지지 밖을 0 으로. real=True 이면 실수부만 남기고 음수도 0 으로 자른다."""
    if real:
        return np.where(sup, np.maximum(x.real, 0.0), 0.0).astype(complex)
    return np.where(sup, x, 0)


def random_start(mag, rng):
    """측정 크기에 무작위 위상을 붙여 시작점을 만든다."""
    return iF(mag * np.exp(2j * np.pi * rng.random(mag.shape)))


def step(alg, x, mag, sup, real=False, beta=0.9):
    pm = P_m(x, mag)
    if alg == "ER":
        return P_s(pm, sup, real)
    if alg == "HIO":
        ok = sup & (pm.real >= 0) if real else sup
        return np.where(ok, pm, x - beta * pm)
    if alg == "RAAR":
        rm = 2 * pm - x
        rs_rm = 2 * P_s(rm, sup, real) - rm
        return 0.5 * beta * (rs_rm + x) + (1 - beta) * pm
    raise ValueError(alg)


def estimate(alg, x, mag, sup, real=False):
    """현재 반복값에서 물체 추정을 꺼낸다. HIO·RAAR 의 x 는 추정이 아니라 구동 변수라서 사영을 거친다."""
    return x if alg == "ER" else P_s(P_m(x, mag), sup, real)


def run(alg, mag, sup, n_iter, rng=None, real=False, beta=0.9, x0=None, er_tail=0, record=()):
    """n_iter 회 반복. er_tail > 0 이면 마지막 그 횟수만큼 ER 로 마무리한다 (Fienup 의 비교 관행).

    record 에 준 반복 번호마다 추정을 저장해 돌려준다.
    """
    x = random_start(mag, np.random.default_rng(rng)) if x0 is None else x0.copy()
    snaps = {}
    for k in range(1, n_iter + 1):
        a = "ER" if k > n_iter - er_tail else alg
        x = step(a, x, mag, sup, real, beta)
        if k in record:
            snaps[k] = estimate(a, x, mag, sup, real)
    return estimate("ER" if er_tail else alg, x, mag, sup, real), snaps


def twin(x):
    """공액 반전 x[n] -> conj(x[-n]). 원점은 배열의 [0, 0] (numpy FFT 규약)."""
    return np.conj(np.roll(np.flip(x), 1, axis=(0, 1)))


def align_error(est, obj):
    """자명한 모호성(평행이동, 공액 반전, 복소 척도)을 맞춘 뒤의 상대 오차.

    평행이동은 상호상관의 최댓값으로, 척도는 최소자승 닫힌 해로 맞춘다. 원래 상과 쌍둥이 상
    중 더 가까운 쪽을 쓴다. 이 정렬을 빠뜨리면 쌍둥이 상으로 수렴한 정답을 실패로 세게 된다.
    """
    best = np.inf
    Fo = F(obj)
    for cand in (est, twin(est)):
        xc = iF(Fo * np.conj(F(cand)))
        sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
        c2 = np.roll(cand, sh, axis=(0, 1))
        a = np.vdot(c2, obj) / np.vdot(c2, c2)
        best = min(best, np.linalg.norm(a * c2 - obj) / np.linalg.norm(obj))
    return best


def fourier_error(est, mag):
    return np.linalg.norm(np.abs(F(est)) - mag) / np.linalg.norm(mag)


def poisson_magnitude(obj, total_photons, rng=None):
    """총 광자 수가 total_photons 가 되도록 세기를 맞추고 푸아송 계수 노이즈를 넣은 뒤 크기로 돌린다."""
    I = np.abs(F(obj)) ** 2
    scale = total_photons / I.sum()
    counts = np.random.default_rng(rng).poisson(I * scale)
    return np.sqrt(counts / scale)

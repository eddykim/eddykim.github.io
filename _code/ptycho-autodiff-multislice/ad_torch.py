"""9편: PyTorch 자동미분으로 쓰는 타이코그래피 전방 모델.

전방 모델만 쓰고 기울기는 자동미분에 맡긴다 (Kandel 등 2019).

  views          시편에서 위치마다 N x N 을 잘라 낸다 (정수 모서리)
  shift_probe    프로브를 위치마다 소수 화소만큼 옮긴다. 이동량 d 에 대해 미분 가능하다
  intensities    I_j = |F[P_j O_j]|^2
  multislice     시편을 여러 장으로 두고 장 사이를 각스펙트럼으로 전파한다 (Maiden 등 2012):
                 psi_1 = P O_1,  psi_{k+1} = Prop_dz(psi_k) O_{k+1},  I = |F psi_last|^2
  loss           6편의 세 비용 (진폭, 푸아송, 세기)

PyTorch 의 복소 기울기는 실수 손실 L 에 대해 2 dL/dO* 다 (6편 grad_obj 의 두 배). 방향은 같다.
"""
import numpy as np
import torch

from ptycho9_core import N, _kz

DT = torch.complex128


def t(a):
    """numpy 배열을 복소 텐서로 복사한다. torch.as_tensor 는 메모리를 공유해, 옵티마이저가 텐서를 고치면 원래
    numpy 배열(시작점)까지 바뀐다. 9편 계산에서 이 때문에 옵티마이저 비교가 앞 옵티마이저의 결과에서 출발한 적이 있다."""
    return torch.tensor(np.array(a), dtype=DT)


def tfft(x):
    return torch.fft.fftshift(torch.fft.fft2(torch.fft.ifftshift(x, dim=(-2, -1)), norm="ortho"), dim=(-2, -1))


def tifft(x):
    return torch.fft.fftshift(torch.fft.ifft2(torch.fft.ifftshift(x, dim=(-2, -1)), norm="ortho"), dim=(-2, -1))


def views(O, pos, n=N):
    pos = torch.as_tensor(np.asarray(pos))
    iy = pos[:, 0, None] + torch.arange(n)[None]
    ix = pos[:, 1, None] + torch.arange(n)[None]
    return O[iy[:, :, None], ix[:, None, :]]


_F = torch.fft.fftfreq(N, dtype=torch.float64)
_FY, _FX = torch.meshgrid(_F, _F, indexing="ij")


def shift_probe(P, d):
    """P(r - d_j) 를 위치마다. d 는 (J, 2) 실수 텐서 (화소)."""
    ph = torch.exp(-2j * torch.pi * (_FY[None] * d[:, 0, None, None] + _FX[None] * d[:, 1, None, None]))
    return torch.fft.ifft2(torch.fft.fft2(P)[None] * ph)


def intensities(O, P, pos, d=None):
    Pj = P[None] if d is None else shift_probe(P, d)
    return torch.abs(tfft(Pj * views(O, pos))) ** 2


def propagator(dz, dx):
    """각스펙트럼 전파 연산자 (중심을 맞춘 주파수 격자)."""
    *_, kz = _kz(N, dx)
    return t(np.exp(2j * np.pi * dz * kz))


def multislice(slices, P, pos, H):
    psi = P[None] * views(slices[0], pos)
    for O in slices[1:]:
        psi = tifft(tfft(psi) * H) * views(O, pos)
    return torch.abs(tfft(psi)) ** 2


def loss(I, n, kind="amp"):
    if kind == "amp":
        return torch.sum((torch.sqrt(I + 1e-12) - torch.sqrt(n)) ** 2)
    if kind == "poisson":
        return torch.sum(I - n * torch.log(I + 1e-12))
    if kind == "int":
        return 0.5 * torch.sum((I - n) ** 2)
    raise ValueError(kind)

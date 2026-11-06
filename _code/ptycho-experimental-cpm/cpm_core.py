"""12편: 실공간 타이코그래피 실측(TwoLayer)의 직접 구현. numpy·scipy·h5py 만 쓴다.

데이터는 PtyLab 예제 데이터(M. Du, figshare 10.6084/m9.figshare.14915628, CC BY 4.0)의 TwoLayer_bin4.hdf5 이다
(get_data.py). 받는 곳은 data/ 이고, 환경 변수 CPM_DATA 로 바꿀 수 있다.
갱신식과 제약은 PtyLab mPIE 의 기본 설정을 따라 하나씩 다시 짰고, 손잡이마다 켜고 끌 수 있다.

  power     매 반복 Σ|P|² 를 가장 밝은 패턴의 총세기로 고정
  power_init 라운드 시작 때 Σ|P|² 를 같은 값으로 맞춤 (PtyLab 은 power 와 한 스위치)
  p0_scale  초기 프로브의 Σ|P|² 를 고정값의 이 배수로 (None 이면 진폭 1 인 원 그대로, 약 0.001 배)
  curv      초기 프로브에 곡률 위상 exp(iπr²/(2λz)) 를 곱함
  order     "ptylab"(라운드마다 처음 두 반복은 순차, 그 뒤 무작위) | "random" | "sequential"
  momentum  위치마다 5 % 확률로 모멘텀 갱신 (mPIE)
  smooth    프로브 진폭을 가우스로 조금 부드럽게 (폭 10, 섞는 비율 0.01)
  l2        짝수 라운드에서 매 반복 O, P 에 0.99 를 곱함
  com       10 반복마다 프로브 무게중심을 창 가운데로 옮김(시편·모멘텀·버퍼도 같이)
"""
import os
import time

import h5py
import numpy as np
from scipy.ndimage import fourier_gaussian

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.environ.get("CPM_DATA", os.path.join(HERE, "data"))
HDF5 = os.path.join(DATA, "TwoLayer_bin4.hdf5")

DEFAULT = dict(power=False, power_init=False, p0_scale=None, curv=True, order="ptylab",
               momentum=True, smooth=True, l2=True, com=True,
               alpha_o=0.1, alpha_p=0.1, beta_o=0.25, beta_p=0.25, friction=0.7, feedback=0.3,
               rounds=2, iters=50)
# 글이 견주는 세 설정. 밝게 = PtyLab 기본(시작 때 맞추고 매 반복 고정), 어둡게 = 원 프로브 그대로 시작해 매 반복 고정
SETTINGS = {"bright": dict(power=True, power_init=True),
            "dim": dict(power=True, power_init=False),
            "free": dict(power=False, power_init=False),
            "bright3": dict(power=False, power_init=False, p0_scale=3.0)}


def fft2c(x):
    return np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(x, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def ifft2c(x):
    return np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(x, axes=(-2, -1)), norm="ortho"), axes=(-2, -1))


def load(fname=HDF5):
    """패턴과 기하. 시편 화소 dxo = λ z / (N dxd) — 원시야(프라운호퍼) 전파의 표본 간격."""
    with h5py.File(fname, "r") as f:
        I = np.asarray(f["ptychogram"], dtype=np.float64)
        d = {k: float(np.asarray(f[k]).ravel()[0]) for k in ("dxd", "wavelength", "zo", "entrancePupilDiameter")}
        d["encoder"] = np.asarray(f["encoder"], dtype=np.float64)
    d["I"] = I
    d["N"] = I.shape[-1]
    d["dxo"] = d["wavelength"] * d["zo"] / (d["N"] * d["dxd"])
    return d


def geometry(d, margin=64):
    """화소 위치 (행, 열) — 패턴 창의 왼쪽 위 모서리. 시편 격자는 주사 범위 + 창 + 여유."""
    N, dxo = d["N"], d["dxo"]
    pix = np.round(d["encoder"] / dxo).astype(int)
    No = N + 2 * int(np.abs(pix).max()) + 2 * margin
    pos = pix + No // 2 - N // 2
    return pos, No


def initial_probe(d, curv=True):
    N, dxo = d["N"], d["dxo"]
    x = np.linspace(-N / 2, N / 2, N) * dxo        # PtyLab 의 프로브 격자 (간격이 dxo 보다 N/(N-1) 배 크다)
    X, Y = np.meshgrid(x, x)
    P = (np.hypot(X, Y) < d["entrancePupilDiameter"] / 2).astype(complex)
    if curv:
        P *= np.exp(1j * np.pi * (X ** 2 + Y ** 2) / (2 * d["wavelength"] * d["zo"]))
    return P


def smooth_amplitude(P, width=10, aleph=1e-2):
    ph = P / (np.abs(P) + 1e-5)
    F = np.fft.fft2(np.abs(P))
    for ax in (-2, -1):
        F = fourier_gaussian(F, width, axis=ax)
    return aleph * np.abs(np.fft.ifft2(F)) * ph + (1 - aleph) * P


def com_center(*arrays):
    """프로브(첫 배열)의 무게중심을 창 가운데로. 프로브만 옮기면 시편과의 짝이 어긋나므로
    시편과 모멘텀·버퍼까지 같은 양만큼 함께 굴린다(곱 P·O 는 위치마다 그대로)."""
    I = np.abs(arrays[0]) ** 2
    n = I.shape[-1]
    c = np.arange(n) - n // 2
    yc = int(np.around(np.sum(c[:, None] * I) / I.sum()))
    xc = int(np.around(np.sum(c[None, :] * I) / I.sum()))
    if xc ** 2 + yc ** 2 > 1:
        return tuple(np.roll(a, (-yc, -xc), (-2, -1)) for a in arrays)
    return arrays


def rpie_step(patch, P, I_j, o):
    """한 위치의 rPIE 갱신(mPIE 의 안쪽). 측정 진폭으로 바꾼 출사파와의 차 Δψ 를 시편과 프로브에 나눠 준다.
    분모 α·max|P|² + (1−α)|P|² 가 갱신의 크기를 정한다. 프로브 갱신은 갱신 전 시편 조각을 쓴다."""
    psi = patch * P
    Psi = fft2c(psi)
    Iest = np.abs(Psi) ** 2
    dpsi = ifft2c(Psi * np.sqrt(I_j / (Iest + 1e-10))) - psi
    aP = np.abs(P) ** 2
    new_patch = patch + o["beta_o"] * np.conj(P) * dpsi / (o["alpha_o"] * aP.max() + (1 - o["alpha_o"]) * aP)
    aO = np.abs(patch) ** 2
    new_P = P + o["beta_p"] * np.conj(patch) * dpsi / (o["alpha_p"] * aO.max() + (1 - o["alpha_p"]) * aO)
    return new_patch, new_P, Iest


def reconstruct(d, seed=0, callback=None, P0=None, pos=None, **opt):
    """mPIE. 반환: 시편, 프로브, 반복별 기록 [Σ|P|², 오차]."""
    o = dict(DEFAULT, **opt)
    rng = np.random.RandomState(seed)
    I = d["I"]
    J, N = I.shape[0], d["N"]
    if pos is None:
        pos, No = geometry(d)
    else:
        No = int(pos.max() + N + 1)
    P = initial_probe(d, o["curv"]) if P0 is None else P0.astype(complex).copy()
    O = np.ones((No, No), complex)
    mp2 = float(np.max(I.sum((-1, -2))))           # 고정값: 가장 밝은 패턴의 총세기 (정규 직교 FFT 라 Σ|ψ|² = ΣI)
    energy = I.sum((-1, -2))
    MO, MP = np.zeros_like(O), np.zeros_like(P)
    if o["p0_scale"] is not None:
        P *= np.sqrt(o["p0_scale"] * mp2 / np.sum(np.abs(P) ** 2))
    trace = []
    t0 = time.time()
    for rnd in range(1, o["rounds"] + 1):
        l2 = o["l2"] and rnd % 2 == 0
        if o["power_init"]:
            P *= np.sqrt(mp2 / np.sum(np.abs(P) ** 2))
        BO, BP = O.copy(), P.copy()
        for it in range(o["iters"]):
            if o["order"] == "sequential" or (o["order"] == "ptylab" and it < 2):
                order = np.arange(J)
            else:
                order = np.arange(J)
                rng.shuffle(order)
            err = 0.0
            for j in order:
                r, c = pos[j]
                O[r:r + N, c:c + N], P, Iest = rpie_step(O[r:r + N, c:c + N].copy(), P, I[j], o)
                err += np.sum(np.abs(I[j] - Iest)) / energy[j]
                if o["momentum"] and rng.rand() > 0.95:
                    MO = (BO - O) + o["friction"] * MO
                    O = O - o["feedback"] * MO
                    BO = O.copy()
                    MP = (BP - P) + o["friction"] * MP
                    P = P - o["feedback"] * MP
                    BP = P.copy()
            if l2:
                O *= 0.99
                P *= 0.99
            if o["power"]:
                P *= np.sqrt(mp2 / np.sum(np.abs(P) ** 2))
            if o["com"] and it % 10 == 0:
                P, MP, BP, O, MO, BO = com_center(P, MP, BP, O, MO, BO)
            if o["smooth"]:
                P = smooth_amplitude(P)
            trace.append((float(np.sum(np.abs(P) ** 2)), float(err)))
            if callback:
                callback(rnd, it, O, P, trace)
        print(f"round {rnd}: error {trace[-1][1]:.4g}, sum|P|^2/mp2 {trace[-1][0] / mp2:.3g}, "
              f"{time.time() - t0:.0f}s", flush=True)
    return O, P, np.array(trace), pos


def data_error(O, P, pos, I):
    """복원 중 오차와 같은 정의: 패턴마다 Σ|I − Î| / ΣI 를 더한 값 (202 장이면 0–404)."""
    N = P.shape[-1]
    err = 0.0
    for j, (r, c) in enumerate(pos):
        Iest = np.abs(fft2c(O[r:r + N, c:c + N] * P)) ** 2
        err += np.sum(np.abs(I[j] - Iest)) / np.sum(I[j])
    return err


def scan_crop(O, pos, N):
    """판정에 쓰는 영역: 프로브 가운데 절반이 지난 주사 영역의, 다시 가운데 절반."""
    r0, c0 = np.min(pos, 0) + N // 4
    r1, c1 = np.max(pos, 0) + 3 * N // 4
    h, w = r1 - r0, c1 - c0
    return O[r0 + h // 4:r1 - h // 4, c0 + w // 4:c1 - w // 4]


def ncc(a, b):
    """두 영상의 정규화 상호상관 최댓값. 평행 이동과 180° 회전(곱 P·O 의 켤레 쌍둥이)을 허용한다."""
    n = [max(a.shape[i], b.shape[i]) * 2 for i in (0, 1)]
    best = 0.0
    for bb in (b, b[::-1, ::-1]):
        A = (a - a.mean()) / (a.std() * np.sqrt(a.size))
        B = (bb - bb.mean()) / (bb.std() * np.sqrt(bb.size))
        c = np.fft.ifft2(np.fft.fft2(A, n) * np.conj(np.fft.fft2(B, n))).real
        best = max(best, c.max() * np.sqrt(a.size * bb.size) / max(a.size, bb.size))
    return float(best)


def fresnel_intensity(psi, dxo, wavelength, z):
    """시편면 출사파 psi(간격 dxo)를 거리 z 의 검출기로 프레넬 전파한 세기.
    검출기 간격은 λz/(N dxo). 이차 위상 exp(iπr²/λz) 를 곱한 뒤 푸리에 변환하는 꼴이다."""
    N = psi.shape[-1]
    x = (np.arange(N) - N // 2) * dxo
    X, Y = np.meshgrid(x, x)
    return np.abs(fft2c(psi * np.exp(1j * np.pi * (X ** 2 + Y ** 2) / (wavelength * z)))) ** 2

"""합성 데이터: zPIE 실험과 같은 기하(파장·검출기·주사 위치)에 정답을 아는 시편과 프로브를 놓는다.

시편은 막대 표적(세 줄 막대 묶음, 선폭 31 µm 에서 4 µm 까지, 가로·세로)과 약한 매끈한 위상이다.
프로브는 두 가지다.
  synth_zpie.npz  지름 572 µm 원을 2 mm 전파해 가장자리를 무디게 한 뒤 곡률 위상 exp(iπr²/(2λz))(완만하게 퍼지는 빔)
  synth_real.npz  실데이터를 z 35 mm 로 복원한 프로브(곡률 반경 약 5 mm, 세게 퍼지는 빔) — --probe 로 그 npz 를 준다
패턴은 원시야(프라운호퍼) 세기에 실측과 같은 총세기로 포아송 잡음을 더한다.

사용: python synth.py [출력.npz] [--probe 복원결과.npz]   (기본 data/synth_zpie.npz)
"""
import os
import sys

import numpy as np

from dz_core import DATA, asp, dxo, fft2c, grid_size, load, positions

Z_TRUE = 0.035


def bar_target(No, dx, rng, field=440):
    """막대 표적의 투과 진폭(유리 1, 크롬 0.05)과 약한 위상(±0.3 rad, 매끈함)."""
    amp = np.ones((No, No))
    widths = 31.25e-6 / 2 ** (np.arange(12) / 4)          # 31 µm → 4 µm, 네 단계마다 절반
    cells = int(np.ceil(np.sqrt(2 * len(widths))))
    cell = field // cells                              # 주사가 지나는 가운데 field 화소 안에 촘촘히
    base = No // 2 - field // 2 + cell // 2
    k = 0
    for i in range(cells):
        for j in range(cells):
            if k >= 2 * len(widths):
                break
            w = widths[k // 2] / dx                     # 선폭 [화소]
            vertical = k % 2 == 0
            r0 = base + i * cell - int(2.5 * w)
            c0 = base + j * cell - int(2.5 * w)
            L = int(5 * w)
            for b in range(3):
                s = int(round(r0 + 2 * b * w)) if not vertical else int(round(c0 + 2 * b * w))
                e = int(round(s + w))
                if vertical:
                    amp[r0:r0 + L, s:e] = 0.05
                else:
                    amp[s:e, c0:c0 + L] = 0.05
            k += 1
    ph = rng.standard_normal((No, No))
    F = np.fft.fft2(ph)
    f = np.fft.fftfreq(No)
    FX, FY = np.meshgrid(f, f)
    ph = np.real(np.fft.ifft2(F * np.exp(-(FX ** 2 + FY ** 2) / (2 * 0.01 ** 2))))
    ph *= 0.3 / np.abs(ph).max()
    return amp * np.exp(1j * ph)


def probe(d, z):
    N, dx = d["N"], dxo(d, z)
    x = (np.arange(N) - N // 2) * dx
    X, Y = np.meshgrid(x, x)
    P = (np.hypot(X, Y) < d["pupil"] / 2).astype(complex)
    P = asp(P, 2e-3, d["wavelength"], dx)
    return P * np.exp(1j * np.pi * (X ** 2 + Y ** 2) / (2 * d["wavelength"] * z))


def make(out, seed=0, probe_file=None):
    rng = np.random.RandomState(seed)
    d = load("zpie")
    total = float(np.median(d["I"].sum((-1, -2))))
    N = d["N"]
    No = grid_size(d, 0.8 * Z_TRUE)
    pos = positions(d, Z_TRUE, No)
    O = bar_target(No, dxo(d, Z_TRUE), rng)
    if probe_file:
        P = np.load(probe_file)["P"].astype(complex)
        curv_R = d["curv_R"]                       # 시작 프로브도 실데이터와 같게(7 mm)
    else:
        P = probe(d, Z_TRUE)
        curv_R = 2 * Z_TRUE
    I = np.empty((len(pos), N, N))
    for j, (r, c) in enumerate(pos):
        I[j] = np.abs(fft2c(O[r:r + N, c:c + N] * P)) ** 2
    I *= total / np.median(I.sum((-1, -2)))
    I = rng.poisson(I).astype(np.float32)
    np.savez_compressed(out, I=I, O=O, P=P, pos=pos, encoder=d["encoder"], z=Z_TRUE,
                        dxd=d["dxd"], wavelength=d["wavelength"], pupil=d["pupil"], curv_R=curv_R)
    print("저장:", out, I.shape, "총세기 중앙값", np.median(I.sum((-1, -2))))


def load_synth(fname):
    """synth 파일을 load() 와 같은 꼴의 dict 로. 정답(O, P, 위치)은 'truth' 에 둔다."""
    s = np.load(fname)
    I = s["I"].astype(np.float64)
    d = dict(I=I, N=I.shape[-1], dxd=float(s["dxd"]), wavelength=float(s["wavelength"]),
             zo=float(s["z"]), encoder=s["encoder"], pupil=float(s["pupil"]))
    d["D"] = d["N"] * d["dxd"]
    if "curv_R" in s.files:
        d["curv_R"] = float(s["curv_R"])
    d["truth"] = dict(O=s["O"], P=s["P"], pos=s["pos"], z=float(s["z"]))
    return d


if __name__ == "__main__":
    args = sys.argv[1:]
    pf = None
    if "--probe" in args:
        k = args.index("--probe")
        pf = args[k + 1]
        args = args[:k] + args[k + 2:]
    make(args[0] if args else os.path.join(DATA, "synth_zpie.npz"), probe_file=pf)

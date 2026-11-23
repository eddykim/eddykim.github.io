"""13편: 시편–검출기 거리 z 와 층 사이 간격을 데이터에서 찾는다.

12편의 mPIE(../ptycho-experimental-cpm/cpm_core.py)를 그대로 가져와 다음을 더한다.

  · z 를 바꾸면 시편 화소 dxo = λz/D 가 바뀌고 주사 격자(화소 단위)가 다시 그려진다 (zPIE 식 1)
  · 각스펙트럼 전파와 TV 선명도 S = Σ sqrt(|∂x O|² + |∂y O|² + ε)  (zPIE 식 2·3·5)
  · zPIE: 매 반복 복원한 시편을 ±1 피사계 심도 안의 11 면으로 전파해 선명도의 무게중심으로 z 를 고친다
  · 두 장 모형(3PIE): 출사파 = O₂ · A_Δz[ O₁ · P ]  (9편)

데이터는 PtyLab 예제 데이터(figshare, CC BY 4.0)이고 이 폴더에 싣지 않는다 (get_data.py).
받는 곳은 data/ 이고, 환경 변수 DZ_DATA 로 바꿀 수 있다. numpy·scipy·h5py 만 쓴다.
"""
import math
import os
import sys
import time

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "ptycho-experimental-cpm"))
from cpm_core import (DEFAULT, com_center, fft2c, ifft2c, ncc, rpie_step,  # noqa: E402
                      smooth_amplitude)

DATA = os.environ.get("DZ_DATA", os.path.join(HERE, "data"))
FILES = {"zpie": os.path.join(DATA, "ptyLab_zPIE.h5"),
         "twolayer": os.path.join(DATA, "TwoLayer_bin4.hdf5")}
# zPIE 파일에는 입사동 지름이 없다. 논문의 빔 크기(FWHM 572 µm)를 초기 원 프로브의 지름으로 쓴다.
ZPIE_BEAM = 572e-6
# zPIE 실험의 빔은 시편 바로 앞에서 모였다가 퍼진다. 검출기에 시편의 그림자가 비치고(점 투사),
# 초기 프로브의 곡률 반경을 12편처럼 2z(70 mm)로 두면 프로브가 한 점으로 무너진다.
# 반경 7·13·25 mm, 부호 반대, 곡률 없음을 20 반복씩 견주어 7 mm 만 시편이 섰다(오차 15.5, 나머지 60–83).
ZPIE_CURV_R = 7e-3


# ── 데이터와 기하 ──────────────────────────────────────────────────────────

def load(name_or_path, pupil=None):
    """패턴, 파장, 검출기 화소, 기록된 z, 주사 위치(가운데로 옮긴 [m], (J, 2) 행·열)."""
    fname = FILES.get(name_or_path, name_or_path)
    with h5py.File(fname, "r") as f:
        I = np.asarray(f["ptychogram"], dtype=np.float64)
        d = {k: float(np.asarray(f[k]).ravel()[0]) for k in ("dxd", "wavelength", "zo")}
        enc = np.asarray(f["encoder"], dtype=np.float64)
        if "entrancePupilDiameter" in f:
            d["pupil"] = float(np.asarray(f["entrancePupilDiameter"]).ravel()[0])
    if enc.shape[0] == 2 and enc.shape[1] != 2:          # zPIE 파일은 (2, J) 로 들어 있다
        enc = enc.T
        # 그리고 두 열이 TwoLayer 와 반대 순서다(열, 행). 이웃한 패턴 사이에서 그림자가 움직인 방향을
        # 주사 위치의 차와 맞춰 보면, TwoLayer 는 검출기 행 ↔ encoder 0 열인데 zPIE 는 검출기 행 ↔ 1 열이다.
        # 여덟 가지 축 뒤집기·바꾸기 가운데 이 바꿈에서만 시편이 선다.
        enc = enc[:, ::-1]
        d["curv_R"] = ZPIE_CURV_R
    d["encoder"] = enc - enc.mean(0)
    d["I"] = I
    d["N"] = I.shape[-1]
    d["D"] = d["N"] * d["dxd"]                           # 검출기 한 변의 길이
    if pupil is not None:
        d["pupil"] = pupil
    d.setdefault("pupil", ZPIE_BEAM)
    return d


def dxo(d, z):
    """시편 화소 λz/D. z 에 비례한다."""
    return d["wavelength"] * z / d["D"]


def dof(d, z):
    """피사계 심도 λ(2z/D)² = λ/NA² (zPIE 식 4)."""
    return d["wavelength"] * (2 * z / d["D"]) ** 2


def grid_size(d, zmin, margin=64):
    """z 가 zmin 까지 작아져도(화소가 작아져 격자가 부풀어도) 들어가는 시편 배열 크기(짝수)."""
    reach = np.abs(d["encoder"]).max() / dxo(d, zmin)
    No = d["N"] + 2 * math.ceil(reach) + 2 * margin
    return No + No % 2


def positions(d, z, No):
    """화소 위치(패턴 창의 왼쪽 위 모서리). p = round(r / dxo) — z 가 틀리면 격자 전체가 r 에 비례해 늘거나 준다."""
    pix = np.round(d["encoder"] / dxo(d, z)).astype(int)
    return pix + No // 2 - d["N"] // 2


def initial_probe(d, z, curv=True):
    """지름 pupil 인 원에 곡률 위상 exp(iπr²/(λR)) 를 곱한다.
    R 은 데이터에 curv_R 이 있으면 그 값(zPIE 7 mm), 없으면 12편과 같이 2z."""
    N, dx = d["N"], dxo(d, z)
    x = np.linspace(-N / 2, N / 2, N) * dx
    X, Y = np.meshgrid(x, x)
    P = (np.hypot(X, Y) < d["pupil"] / 2).astype(complex)
    if curv:
        R = d.get("curv_R", 2 * z)
        P *= np.exp(1j * np.pi * (X ** 2 + Y ** 2) / (d["wavelength"] * R))
    return P


# ── 전파와 선명도 ─────────────────────────────────────────────────────────

def asp(u, dz, wavelength, dx):
    """각스펙트럼 전파 A_dz (zPIE 식 2·3). 소멸파는 버린다."""
    n = u.shape[-1]
    f = np.fft.fftfreq(n, dx)
    FX, FY = np.meshgrid(f, f)
    arg = 1 / wavelength ** 2 - FX ** 2 - FY ** 2
    H = np.where(arg > 0, np.exp(2j * np.pi * dz * np.sqrt(np.maximum(arg, 0))), 0)
    return np.fft.ifft2(np.fft.fft2(u) * H)


def tv(u, eps=1e-2):
    """TV 선명도 Σ sqrt(|∂x u|² + |∂y u|² + ε) (zPIE 식 5, 중앙 차분)."""
    gx = np.roll(u, -1, -1) - np.roll(u, 1, -1)
    gy = np.roll(u, -1, -2) - np.roll(u, 1, -2)
    return float(np.sum(np.sqrt(np.abs(gx) ** 2 + np.abs(gy) ** 2 + eps)))


def focus_roi(O, n):
    """시편 배열 가운데 n×n 에 가장자리로 갈수록 0 이 되는 창(초가우스)을 곱한 것.
    창이 없으면 FFT 의 주기 경계가 만든 가장자리가 TV 에 섞인다 (PtyLab zPIE 와 같은 창)."""
    No = O.shape[-1]
    s = slice(No // 2 - n // 2, No // 2 + n // 2)
    k = np.arange(-n // 2, n // 2)
    X, Y = np.meshgrid(k, k)
    w = np.exp(-(np.hypot(X, Y) / (n / 2)) ** 4)
    return w * O[s, s]


def sharpness_curve(O, d, z, dzs, n=None):
    """복원한 시편을 dzs 의 각 거리로 전파한 뒤의 TV. 반환: TV 배열."""
    n = n or d["N"]
    roi = focus_roi(O, n)
    return np.array([tv(asp(roi, dz, d["wavelength"], dxo(d, z))) for dz in dzs])


# ── 복원: z 를 고정하거나(zpie=None) zPIE 로 고쳐 가며 ─────────────────────────

ZPIE = dict(step=1.0, friction=0.7, start=1, K=11, span=1.0, sign=+1, roi=None)


def reconstruct(d, z0, seed=0, zpie=None, zmin=None, P0=None, callback=None, **opt):
    """12편 mPIE 와 같은 갱신. 다른 점은 반복마다 z 에서 주사 격자를 다시 그린다는 것뿐이다.

    zpie 가 dict 이면 Algorithm 1 의 z 갱신을 한다:
        s(k) = S[A_{dz_k} O],  dz_k ∈ span·DoF·linspace(-1, 1, K)
        m ← friction·m + step·Σ s(k) dz_k / Σ s(k),   z ← z + sign·m
    반환: 시편, 프로브, 기록(반복마다 [Σ|P|², 오차, z]), 마지막 위치.
    """
    o = dict(DEFAULT, **opt)
    zp = None if zpie is None else dict(ZPIE, **zpie)
    rng = np.random.RandomState(seed)
    I = d["I"]
    J, N = I.shape[0], d["N"]
    No = grid_size(d, zmin or 0.8 * z0)
    z = float(z0)
    pos = positions(d, z, No)
    P = initial_probe(d, z, o["curv"]) if P0 is None else P0.astype(complex).copy()
    O = np.ones((No, No), complex)
    mp2 = float(np.max(I.sum((-1, -2))))
    energy = I.sum((-1, -2))
    MO, MP = np.zeros_like(O), np.zeros_like(P)
    if o["p0_scale"] is not None:
        P *= np.sqrt(o["p0_scale"] * mp2 / np.sum(np.abs(P) ** 2))
    m = 0.0
    trace, merits = [], []
    t0 = time.time()
    n_iter = 0
    for rnd in range(1, o["rounds"] + 1):
        l2 = o["l2"] and rnd % 2 == 0
        if o["power_init"]:
            P *= np.sqrt(mp2 / np.sum(np.abs(P) ** 2))
        BO, BP = O.copy(), P.copy()
        for it in range(o["iters"]):
            # ── zPIE: 이번 반복에 쓸 z ──
            if zp is not None and n_iter >= zp["start"]:
                dzs = zp["span"] * dof(d, z) * np.linspace(-1, 1, zp["K"])
                s = sharpness_curve(O, d, z, dzs, zp["roi"])
                feedback = float(np.sum(s * dzs) / np.sum(s))
                m = zp["friction"] * m + zp["step"] * feedback
                z = z + zp["sign"] * m
                merits.append(s)
                pos = positions(d, z, No)
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
            trace.append((float(np.sum(np.abs(P) ** 2)), float(err), z))
            n_iter += 1
            if callback:
                callback(rnd, it, O, P, trace)
        print(f"round {rnd}: error {trace[-1][1]:.4g}, z {z * 1e3:.4f} mm, "
              f"sum|P|^2/mp2 {trace[-1][0] / mp2:.3g}, {time.time() - t0:.0f}s", flush=True)
    return O, P, np.array(trace), pos, np.array(merits)


def _rpie(target, probe_like, dpsi, alpha, beta):
    """rPIE 꼴 갱신 한 번: target += β conj(probe)·Δψ / (α max|probe|² + (1−α)|probe|²)."""
    a = np.abs(probe_like) ** 2
    return target + beta * np.conj(probe_like) * dpsi / (alpha * a.max() + (1 - alpha) * a)


def reconstruct_two(d, gap, z=None, seed=0, P0=None, O0=None, callback=None, fresnel=True, **opt):
    """두 장 모형(3PIE). 출사파 ψ₂ = O₂ · A_gap[O₁ · P], 검출기 = F[Q ψ₂], Q = exp(iπr²/λz).

    한 장 모형은 프레넬 이차 위상 Q 를 프로브에 흡수해 원시야 꼴로 쓴다(12편 1절). 두 장에서는 그럴 수 없다.
    Q 는 뒤 장을 지난 출사파에 곱해져야 하는데, 프로브에 넣어 두면 장 사이 전파 앞에 곱해진다. 프레넬 축척 정리에
    따르면 그것은 뒤 장을 1 + gap/z 배 확대한 모형이다(gap 1.6 mm, z 38.6 mm 면 4 %). 그래서 여기서는 P 를 실제 빔으로
    두고 Q 를 검출기 쪽에 따로 곱한다(fresnel=True). P0 로 한 장 결과의 프로브를 받으면 Q 를 떼어 내고 시작한다.

    한 위치의 갱신(Maiden 2012, 갱신식은 12편과 같은 rPIE 꼴):
      ① 측정 진폭으로 ψ₂ 를 고쳐 Δ₂ = ψ₂' − ψ₂
      ② O₂ 를 φ = A_gap[O₁P] 를 프로브 삼아, φ 를 O₂ 를 시편 삼아 갱신
      ③ 고친 φ' 를 −gap 만큼 되돌려 ψ₁' = A₋gap φ',  Δ₁ = ψ₁' − O₁P
      ④ O₁ 과 P 를 한 장 때와 같이 갱신
    gap=0 이면 두 장이 한 장의 곱과 같다(대조군). O0·P0 이 있으면 한 장 복원에서 이어 시작한다
    (어둡게 시작하면 O₂ 와 φ 가 같은 보정을 둘 다 떠안아 발산한다 — 11-22 확인).
    격자 크기는 reconstruct 와 같게(zmin = 0.8 z) 잡아 한 장 결과를 그대로 넘겨받을 수 있게 한다.
    반환: O₁, O₂, P, 기록[Σ|P|², 오차], 위치.
    """
    o = dict(DEFAULT, **opt)
    rng = np.random.RandomState(seed)
    I = d["I"]
    J, N = I.shape[0], d["N"]
    z = z or d["zo"]
    No = grid_size(d, 0.8 * z) if O0 is None else O0.shape[-1]
    pos = positions(d, z, No)
    dx, wl = dxo(d, z), d["wavelength"]
    P = initial_probe(d, z, o["curv"]) if P0 is None else P0.astype(complex).copy()
    xq = (np.arange(N) - N // 2) * dx
    XQ, YQ = np.meshgrid(xq, xq)
    Q = np.exp(1j * np.pi * (XQ ** 2 + YQ ** 2) / (wl * z)) if fresnel else np.ones((N, N), complex)
    P = P * np.conj(Q)                                    # 한 장 모형의 프로브 = 실제 빔 × Q
    O1 = np.ones((No, No), complex) if O0 is None else O0.astype(complex).copy()
    O2 = np.ones((No, No), complex)
    mp2 = float(np.max(I.sum((-1, -2))))
    energy = I.sum((-1, -2))
    n = N
    f = np.fft.fftfreq(n, dx)
    FX, FY = np.meshgrid(f, f)
    arg = 1 / wl ** 2 - FX ** 2 - FY ** 2
    H = np.where(arg > 0, np.exp(2j * np.pi * gap * np.sqrt(np.maximum(arg, 0))), 0)
    Hb = np.conj(H)
    fwd = lambda u: np.fft.ifft2(np.fft.fft2(u) * H)       # noqa: E731
    bwd = lambda u: np.fft.ifft2(np.fft.fft2(u) * Hb)      # noqa: E731
    trace = []
    t0 = time.time()
    for rnd in range(1, o["rounds"] + 1):
        if o["power_init"]:
            P *= np.sqrt(mp2 / np.sum(np.abs(P) ** 2))
        for it in range(o["iters"]):
            order = np.arange(J)
            if not (o["order"] == "sequential" or (o["order"] == "ptylab" and it < 2)):
                rng.shuffle(order)
            err = 0.0
            for j in order:
                r, c = pos[j]
                s = (slice(r, r + N), slice(c, c + N))
                o1, o2 = O1[s].copy(), O2[s].copy()
                psi1 = o1 * P
                phi = fwd(psi1)
                psi2 = o2 * phi
                Psi = fft2c(Q * psi2)
                Iest = np.abs(Psi) ** 2
                d2 = np.conj(Q) * ifft2c(Psi * np.sqrt(I[j] / (Iest + 1e-10))) - psi2
                O2[s] = _rpie(o2, phi, d2, o["alpha_o"], o["beta_o"])
                phi_new = _rpie(phi, o2, d2, o["alpha_p"], o["beta_p"])
                d1 = bwd(phi_new) - psi1
                O1[s] = _rpie(o1, P, d1, o["alpha_o"], o["beta_o"])
                P = _rpie(P, o1, d1, o["alpha_p"], o["beta_p"])
                err += np.sum(np.abs(I[j] - Iest)) / energy[j]
            if o["power"]:
                P *= np.sqrt(mp2 / np.sum(np.abs(P) ** 2))
            if o["com"] and it % 10 == 0:
                P, O1, O2 = com_center(P, O1, O2)
            if o["smooth"]:
                P = smooth_amplitude(P)
            trace.append((float(np.sum(np.abs(P) ** 2)), float(err)))
            if callback:
                callback(rnd, it, O1, O2, P, trace)
        print(f"round {rnd}: gap {gap * 1e6:.0f} um, error {trace[-1][1]:.4g}, {time.time() - t0:.0f}s", flush=True)
    return O1, O2, P * Q, np.array(trace), pos         # 프로브는 한 장 모형과 같은 꼴(× Q)로 돌려준다


def data_error(O, P, pos, I):
    """복원 중 오차와 같은 정의: 패턴마다 Σ|I − Î| / ΣI 를 더한 값."""
    N = P.shape[-1]
    err = 0.0
    for j, (r, c) in enumerate(pos):
        Iest = np.abs(fft2c(O[r:r + N, c:c + N] * P)) ** 2
        err += np.sum(np.abs(I[j] - Iest)) / np.sum(I[j])
    return err


def scan_region(O, pos, N, frac=0.25):
    """프로브 가운데 절반이 지난 주사 영역(판정·그림용)."""
    r0, c0 = np.min(pos, 0) + int(N * frac)
    r1, c1 = np.max(pos, 0) + int(N * (1 - frac))
    return O[r0:r1, c0:c1]


__all__ = ["load", "dxo", "dof", "grid_size", "positions", "initial_probe", "asp", "tv",
           "focus_roi", "sharpness_curve", "reconstruct", "reconstruct_two", "data_error", "scan_region",
           "fft2c", "ifft2c", "ncc", "rpie_step", "DATA", "FILES", "ZPIE"]

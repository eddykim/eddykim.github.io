"""11편의 표와 그림에 들어가는 수치. python3 get_data.py 로 데이터를 받은 뒤 python3 studies11.py.
단계마다 .cache_post11.pkl 에 덧붙인다 (generate_figures.py 가 읽는다).

  fmt        UConn 영상을 진폭으로 읽을 때와 세기로 읽을 때의 잔차
  axes       LED 배열의 축 바꿈·부호 8가지의 잔차 (5회)
  boundary   조명 NA 와 영상 평균 밝기 (USAF, 혈액 도말) — 명·암시야 경계
  noise      영상별 노이즈 바닥, 가운데 영상의 세기 스펙트럼 반지름 분포
  ladder     보정 사다리 (기본, 유리 무시, 초점 어긋남 무시, LED 밝기, +동공, +LED 위치), 대물 NA 훑기, 차이 영상
  blood      보정 동공이 있는 혈액 데이터: 이상 동공 / 보정 동공 / EPRY / 초점 어긋남 훑기 -> EPRY
  resol      USAF 최종 복원, LED 를 바둑판으로 나눈 두 복원의 FRC
  ptylab     PtyLab USAF 데이터의 같은 사다리
  usaf       resol 의 복원에서 USAF 9군 1–4번 세로 막대의 대비, LED 덮개가 정하는 주기 한계
"""
import os
import pickle
import time

import numpy as np

import fpm_core as C

CACHE = os.path.join(C.HERE, ".cache_post11.pkl")
UP = 4                        # 고해상 격자 배율 nh = UP * nl


def recon(d, n_iter=20, ks=None, P0=None, **kw):
    ks = np.rint(d["k"]).astype(int) if ks is None else ks
    P0 = C.defocus_pupil(d["nl"], d["r"], d["z"], d["wl"], d["fov"]) if P0 is None else P0
    S, P, c = C.fp_cal(d["I"], ks, P0, UP * d["nl"], n_iter, **kw)
    return S, P, c, C.residual(S, P, ks, d["I"])


CAL = dict(bright_from=3, pupil_from=3)


def fmt():
    return {lab: recon(C.load_uconn("USAF_red.mat", as_intensity=a))[3] for lab, a in (("amp", False), ("int", True))}


def axes():
    out = {}
    for swap in (False, True):
        for sx in (1, -1):
            for sy in (1, -1):
                out[(swap, sx, sy)] = recon(C.load_uconn("USAF_red.mat", swap=swap, sx=sx, sy=sy), 5)[3]
    return out


def boundary():
    out = {}
    for f in ("USAF_red.mat", "bloodsmear_green.mat"):
        d = C.load_uconn(f)
        out[f] = dict(na=d["na_ill"], mean=d["I"].mean((1, 2)), NA=d["NA"])
    return out


def noise():
    d = C.load_uconn("USAF_red.mat")
    nf = C.noise_floor(d["I"], d["r"])
    F = np.abs(C.fft2c(d["I"][0])) ** 2
    y, x = np.indices(F.shape) - d["nl"] // 2
    rr = np.hypot(x, y).astype(int)
    prof = np.bincount(rr.ravel(), F.ravel()) / np.bincount(rr.ravel())
    return dict(nf=nf, na=d["na_ill"], NA=d["NA"], prof=prof, r=d["r"], spec=np.log10(F + 1e-12))


def ladder():
    d = C.load_uconn("USAF_red.mat")
    ks = np.rint(d["k"]).astype(int)
    out = {}
    out["기본"] = recon(d)[3]
    out["유리 무시"] = recon(C.load_uconn("USAF_red.mat", glass=False))[3]
    dz = dict(d); dz["z"] = 0
    out["초점 어긋남 무시"] = recon(dz)[3]
    out["LED 밝기"] = recon(d, bright_from=3)[3]
    S, P, c, out["밝기 + 동공"] = recon(d, **CAL)
    out["밝기 배율 범위"] = (float(c.min()), float(c.max()))
    # 배율을 묶지 않으면: 중앙값 고정을 빼고 같은 복원
    a = np.sqrt(d["I"]); c_un = np.ones(len(ks)); nl = d["nl"]; nh = UP * nl
    Sx = np.zeros((nh, nh), complex); Sx[C.band(nh, nl, (0, 0))] = C.fft2c(a[0])
    Px = C.defocus_pupil(nl, d["r"], d["z"], d["wl"], d["fov"]); trace = []
    for it in range(1, 21):
        for j in np.argsort(np.sum(ks ** 2, 1)):
            sl = C.band(nh, nl, ks[j]); so = Sx[sl].copy(); psi = Px * so; u = C.ifft2c(psi)
            if it >= 3:
                c_un[j] = np.sum(np.abs(u) * a[j]) / np.sum(a[j] ** 2)
            diff = C.fft2c(c_un[j] * a[j] * np.exp(1j * np.angle(u))) - psi
            a2 = np.abs(Px) ** 2
            Sx[sl] = so + np.conj(Px) * diff / (0.9 * a2 + 0.1 * a2.max())
        trace.append(float(np.median(c_un)))
    out["묶지 않은 배율 중앙값"] = trace
    ks2, S2, P2, hist = C.refine_k(d["I"], ks, C.defocus_pupil(nl, d["r"], d["z"], d["wl"], d["fov"]), nh, **CAL)
    out["+ LED 위치"] = hist[-1]
    out["LED 위치 이력"] = hist
    out["LED 이동 rms"] = float(np.sqrt(np.mean(np.sum((ks2 - ks) ** 2, 1))))
    out["가장자리 35 제외"] = C.residual(S2, P2, ks2, d["I"], crop=35)
    nf = C.noise_floor(d["I"], d["r"])
    out["노이즈 바닥 명시야"] = float(np.median(nf[d["na_ill"] < d["NA"]]))
    out["노이즈 바닥 암시야"] = float(np.median(nf[d["na_ill"] >= d["NA"]]))
    m = C.model_amps(S2, P2, ks2)
    s = np.sum(m * a, (1, 2)) / np.sum(m * m, (1, 2))
    order = np.argsort(d["na_ill"])
    out["차이"] = {int(j): (a[j], s[j] * m[j], float(d["na_ill"][j])) for j in order[[0, 16]]}
    na_scan = {}
    for na in (0.08, 0.09, 0.10, 0.11, 0.12):
        r = na / d["wl"] * d["fov"]
        na_scan[na] = recon(d, P0=C.defocus_pupil(nl, r, d["z"], d["wl"], d["fov"]), **CAL)[3]
    out["NA 훑기"] = na_scan
    return out


def blood():
    d = C.load_uconn("blood_aberration.mat")
    nl = d["nl"]; ab = d["ab"]; supp = np.abs(C.defocus_pupil(nl, d["r"], 0, 1, 1)) > 0
    out = {}
    ideal = C.defocus_pupil(nl, d["r"], 0, 1, 1)
    for lab, P0, kw in (("이상 동공", ideal, dict(bright_from=3)), ("이상 동공 + EPRY", ideal, CAL),
                        ("보정 동공", ab, dict(bright_from=3)), ("보정 동공 + EPRY", ab, CAL)):
        S, P, c, r = recon(d, 30, P0=P0, **kw)
        out[lab] = dict(resid=r, P=P, phase_rms=C.pupil_phase_rms(P, ab, supp), o=C.ifft2c(S))
    zs = np.arange(-150e-6, 151e-6, 2e-6)
    fit = [C.pupil_phase_rms(C.defocus_pupil(nl, d["r"], z, d["wl"], d["fov"]), ab, supp) for z in zs]
    out["보정 동공에 맞춘 z"] = (float(zs[int(np.argmin(fit))]), float(min(fit)))
    scan = {}
    for z in np.arange(-100e-6, 101e-6, 10e-6):
        scan[float(z)] = recon(d, 15, P0=C.defocus_pupil(nl, d["r"], z, d["wl"], d["fov"]), bright_from=3)[3]
    out["z 훑기"] = scan
    zb = min(scan, key=scan.get)
    S, P, c, r = recon(d, 60, P0=C.defocus_pupil(nl, d["r"], zb, d["wl"], d["fov"]), **CAL)
    out["z 맞춘 뒤 EPRY"] = dict(z=zb, resid=r, P=P, phase_rms=C.pupil_phase_rms(P, ab, supp), o=C.ifft2c(S))
    out["보정 동공 배열"] = ab
    out["지지"] = supp
    return out


def resol():
    d = C.load_uconn("USAF_red.mat")
    ks = np.rint(d["k"]).astype(int)
    S, P, c, r = recon(d, 30, **CAL)
    sp = d["sp"].astype(int)
    bf = d["na_ill"] < d["NA"]
    g1 = bf | ((sp[:, 0] + sp[:, 1]) % 2 == 0)
    g2 = bf | ((sp[:, 0] + sp[:, 1]) % 2 == 1)
    halves = []
    for g in (g1, g2):
        dd = dict(d); dd["I"], dd["k"] = d["I"][g], d["k"][g]
        halves.append(C.ifft2c(recon(dd, 30, **CAL)[0]))
    e = 40
    f, v = C.frc(np.abs(halves[0][e:-e, e:-e]), np.abs(halves[1][e:-e, e:-e]))
    dxh = C.UC["px"] / UP
    return dict(o=C.ifft2c(S), raw=np.sqrt(d["I"][0]), f=f, v=v, dxh=dxh,
                na_max=float(d["na_ill"].max()), NA=d["NA"], wl=d["wl"], resid=r)


def ptylab():
    d = C.load_ptylab("USAFTargetFPM.hdf5")
    out = dict(bg=d["bg"])
    out["기본"] = recon(d)[3]
    out["밝기 + 동공"] = recon(d, 30, **CAL)[3]
    out["노이즈 바닥 명시야"] = float(np.median(C.noise_floor(d["I"], d["r"])[d["na_ill"] < d["NA"]]))
    return out


# 9군 원소별 세로 막대 영역 (원본 화소 50:90, 65:105 를 4배 키운 창 안의 좌표)과 주기 [um] = 1000 / 2^(9 + (e-1)/6)
USAF9 = {1: ((28, 37), (110, 126)), 2: ((40, 49), (112, 128)), 3: ((54, 64), (112, 128)), 4: ((68, 77), (113, 129))}


def usaf():
    from scipy.signal import find_peaks
    Q = load()["resol"]
    sub = np.abs(Q["o"])[50 * UP:90 * UP, 65 * UP:105 * UP]
    out = {}
    for e, ((y0, y1), (x0, x1)) in USAF9.items():
        p = sub[y0:y1, x0:x1].mean(0)
        mins, _ = find_peaks(-p, distance=2)
        maxs, _ = find_peaks(p, distance=2)
        inner = [m for m in maxs if (mins < m).any() and (mins > m).any()]
        hi, lo = p[inner].mean(), p[mins].mean()
        out[e] = dict(period=1000 / 2 ** (9 + (e - 1) / 6), contrast=float((hi - lo) / (hi + lo)))
    d = C.load_uconn("USAF_red.mat")
    na_ax = float(np.abs(d["k"] * d["wl"] / d["fov"]).max())
    out["축 방향 한계 주기"] = d["wl"] / (d["NA"] + na_ax) * 1e6
    out["대각선 한계 주기"] = d["wl"] / (d["NA"] + d["na_ill"].max()) * 1e6
    out["원본 한계 주기"] = d["wl"] / d["NA"] * 1e6
    f, v = Q["f"], Q["v"]
    i = 3 + int(np.argmax(v[3:] < 0.5))
    out["FRC 0.5 주기"] = Q["dxh"] / f[i] * 1e6
    return out


STAGES = (("fmt", fmt), ("axes", axes), ("boundary", boundary), ("noise", noise), ("ladder", ladder),
          ("blood", blood), ("resol", resol), ("ptylab", ptylab), ("usaf", usaf))


def run_all():
    t = time.time()
    R = load() if os.path.exists(CACHE) else {}
    for name, f in STAGES:
        if name in R:
            continue
        R[name] = f()
        with open(CACHE, "wb") as fh:
            pickle.dump(R, fh)
        print(f"{name}: {time.time() - t:.0f} s", flush=True)
    return R


def load():
    with open(CACHE, "rb") as fh:
        return pickle.load(fh)


if __name__ == "__main__":
    run_all()

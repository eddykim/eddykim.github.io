"""diffraction.py 를 독립된 경로들과 교차검증한다.

실행: python verify_diffraction.py

  A. 각스펙트럼 전파            ↔  레일리-좀머펠트 직접 적분 (2차 파원의 합). 근축 전달함수 ↔ prysm angular_spectrum (있으면;
                                    prysm 의 이 함수는 이름과 달리 프레넬 전달함수 exp(−iπλz f²) 를 쓴다)
  B. 원형 구멍의 축 위 세기     ↔  정확해 |1 − (z/R) e^{ik(R−z)}|² 와 프레넬 근사 4 sin²(πN_F/2)
  C. 아라고 점                  ↔  원판 뒤 축 위 세기 = (z/R)² ≈ 1 (가리지 않은 세기)
  D. 직선 가장자리              ↔  프레넬 적분 (scipy): 경계에서 1/4, 첫 봉우리 1.370 (보른과 울프 8.7)
  E. 에어리 무늬                ↔  2J1(x)/x, 첫 영점 1.22 λ/D, 고리 세기·원 안 에너지 (보른과 울프 표 8.2, 식 8.5-18)
  F. 레일리 기준                ↔  두 점 사이 골의 세기 0.735 (보른과 울프 8.6.2)
  G. 초점 근처                  ↔  축 위 [sin(u/4)/(u/4)]², 첫 영점 z = 2λ/NA² (보른과 울프 8.8)
  H. 존 플레이트                ↔  초점 f = r1²/λ, f/3, f/5 에도 초점, 열린 존 20 개면 세기 1600 배 (Hecht 10.3.5).
                                    f/3, f/5 의 축 위 세기도 1600 (근축)
  I. 격자                       ↔  주봉우리 a sinθ = mλ, 나트륨 D 선을 2차에서 가르는 최소 줄 수 494 (Hecht 10.2.8)
  J. 얇은 격자 효율             ↔  정현 진폭 격자 1차 6.25 %, 정현 위상 격자 최대 33.8 % (Goodman 4.4.3–4.4.4)
  K. 틈 하나                    ↔  곁봉우리 0.0472, 0.0165 (보른과 울프 표 8.1)
"""
import os
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.special import fresnel, j1, jv

from diffraction import (rs_direct, angular_spectrum, fresnel_transfer, axial_circular_aperture, axial_ring_sum, axial_zone_plate, circle_pupil,
                         disk_shadow, encircled_energy, far_field, focal_field, radial_profile,
                         slit_array_far_field, thin_grating_efficiencies, zone_edges)

ORACLE = os.environ.get("ORACLE_PY", os.path.expanduser(
    "~/Library/CloudStorage/Dropbox-Yjedward.kim/Kim Youngjoon/repository/repo_opticore/reference_libs/.venv-oracle/bin/python"))
PRYSM_SRC = os.path.join(os.path.dirname(os.path.dirname(ORACLE)), "..", "prysm")


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<64s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def verify_propagation():
    print("A. 각스펙트럼 전파 ↔ 레일리-좀머펠트 직접 적분")
    lam, a, z = 0.5, 10.0, 60.0                     # 반지름 10 µm 구멍, 60 µm 뒤 (N_F = 3.3, 각도가 커서 프레넬 근사 밖)
    n, dx = 2048, 0.125                             # 표본화 조건 z ≤ N dx²/λ = 64 µm 안
    y, x = (np.mgrid[:n, :n] - n // 2) * dx
    U0 = (np.hypot(x, y) <= a).astype(complex)
    U = angular_spectrum(U0, dx, lam, z)
    xs = np.array([0.0, 3.0, 7.5, 12.0, 18.0])
    sel = np.hypot(x, y) <= a
    pts = [(xx, 0.0, z) for xx in xs]
    direct = rs_direct(x[sel], y[sel], U0[sel], dx * dx, pts, lam)
    got = U[n // 2, n // 2 + np.round(xs / dx).astype(int)]
    ok = check("|U|² 다섯 점: 각스펙트럼 ↔ 직접 적분", np.abs(got) ** 2, np.abs(direct) ** 2, 0.02)
    # prysm (선택): prysm 의 angular_spectrum 은 근축(프레넬) 전달함수 exp(−iπλz f²) 를 쓴다.
    # 그래서 각도가 작은 조건(반지름 50 µm, 2 mm 뒤, N_F = 2.5)에서 대조한다.
    # 전달함수를 FFT 격자에 표본화할 수 있는 거리 z ≤ N dx²/λ (= 1024 µm) 안에서 잰다.
    a2, z2, n2, dx2 = 30.0, 900.0, 2048, 0.5
    y2, x2 = (np.mgrid[:n2, :n2] - n2 // 2) * dx2
    U0b = (np.hypot(x2, y2) <= a2).astype(complex)
    xs2 = np.array([0.0, 10.0, 25.0, 40.0, 60.0])
    idx = n2 // 2 + np.round(xs2 / dx2).astype(int)
    got_as = np.abs(angular_spectrum(U0b, dx2, lam, z2)[n2 // 2, idx]) ** 2
    got2 = np.abs(fresnel_transfer(U0b, dx2, lam, z2)[n2 // 2, idx]) ** 2
    ok &= check("근축(N_F 2): 프레넬 전달함수 ↔ 각스펙트럼", got2, got_as, 3e-3)
    if os.path.exists(ORACLE):
        code = f"""
import sys, numpy as np
sys.path.insert(0, {os.path.abspath(PRYSM_SRC)!r})
from prysm.propagation import angular_spectrum
n, dx = {n2}, {dx2}
y, x = (np.mgrid[:n, :n] - n // 2) * dx
U0 = (np.hypot(x, y) <= {a2}).astype(complex)
U = angular_spectrum(U0, {lam}, dx * 1e-3, {z2} * 1e-3, Q=1)
print(repr([float(abs(v)**2) for v in U[n//2, n//2 + np.round(np.array({list(xs2)}) / dx).astype(int)]]))
"""
        try:
            r = subprocess.run([ORACLE, "-c", code], capture_output=True, text=True, timeout=300,
                               env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
            vals = eval(r.stdout.strip().splitlines()[-1])
            ok &= check("근축(N_F 2): 프레넬 전달함수 ↔ prysm angular_spectrum", got2, vals, 1e-9)
        except Exception as e:  # noqa: BLE001
            print(f"  [SKIP] prysm 대조 실패: {e}")
    else:
        print("  [SKIP] prysm 오라클 환경 없음")
    return ok


def verify_axial_aperture():
    print("B. 원형 구멍의 축 위 세기")
    lam, a = 0.5, 20.0
    k = 2 * np.pi / lam
    zs = np.array([30.0, 55.0, 100.0, 200.0, 400.0])
    got = np.array([axial_circular_aperture(a, z, lam) for z in zs])
    R = np.sqrt(zs**2 + a**2)
    exact = np.abs(1 - zs / R * np.exp(1j * k * (R - zs))) ** 2
    ok = check("직접 적분 ↔ 정확해 |1 − (z/R)e^{ik(R−z)}|²", got, exact, 1e-4)
    a2, zs2 = 500.0, np.array([1e5, 2e5, 5e5, 1e6])       # 근축 영역: 프레넬 근사와 비교
    got2 = np.array([axial_circular_aperture(a2, z, lam, n=20001) for z in zs2])
    NF = a2**2 / (lam * zs2)
    ok &= check("근축에서 ↔ 4 sin²(π N_F / 2)", got2, 4 * np.sin(np.pi * NF / 2) ** 2, 1e-3)
    # 존 경계까지 열면: 존 1 개 → 4, 2 개 → 거의 0
    z = 1e5
    e = zone_edges(z, lam, 2)
    ok &= check("첫 존만 열면 4, 둘 열면 ≈ 0", [axial_circular_aperture(e[1], z, lam), axial_circular_aperture(e[2], z, lam)], [4, 0], 1e-3)
    return ok


def verify_arago():
    print("C. 아라고 점 (원판 뒤 축 위)")
    lam = 0.55
    ok = True
    for a, z in [(2900.0, 0.83e6), (500.0, 2e5)]:
        I = np.abs(disk_shadow(a, 0.0, z, lam, n=200001)[0]) ** 2
        R = np.sqrt(z**2 + a**2)
        ok &= check(f"반지름 {a / 1000:g} mm, {z / 1e6:g} m: 중심 세기 ↔ (z/R)²", I, (z / R) ** 2, 2e-3)
    # 정확한 핵으로 띠를 직접 합해도 같다 (원판 바깥 전체 = 평면파 − 구멍)
    z, a = 50.0, 15.0
    ring = axial_ring_sum(np.linspace(0, a, 200), z, lam).sum()
    I_exact = np.abs(np.exp(1j * 2 * np.pi / lam * z) - ring) ** 2
    ok &= check("가까이(z = 50 µm, a = 15 µm), 정확한 핵 ↔ (z/R)²", I_exact, z**2 / (z**2 + a**2), 1e-3)
    return ok


def verify_edge():
    print("D. 직선 가장자리 ↔ 프레넬 적분")
    lam, z = 0.5, 2000.0
    n, dx = 2**17, 0.1
    x = (np.arange(n) - n // 2) * dx
    open_ = (x >= 0) & (x < 5500)
    taper = np.clip((5500 - x) / 2000, 0, 1)         # 먼 쪽 끝은 부드럽게 닫아 그 회절이 들어오지 않게
    U = angular_spectrum(np.where(open_, taper, 0).astype(complex), dx, lam, z)
    I = np.abs(U) ** 2
    w = x * np.sqrt(2 / (lam * z))
    sel = (w > -3) & (w < 5)
    C, S = fresnel(w[sel])
    ref = 0.5 * ((C + 0.5) ** 2 + (S + 0.5) ** 2)
    # 프레넬 적분은 근축 근사의 답이고 각스펙트럼은 근사가 없다. 둘의 차이는 z 가 짧을수록 커진다 (z = 2 mm 에서 0.4 %)
    ok = check("세기 곡선 (w ∈ [−3, 5]) ↔ ½[(C+½)² + (S+½)²]", I[sel], ref, 5e-3)
    ok &= check("경계 (x = 0) 세기 = 1/4", I[n // 2], 0.25, 2e-3)
    i = np.argmax(np.where(sel, I, 0))
    ok &= check("첫 봉우리 세기 1.370 (w = 1.217)", [I[i], w[i]], [1.3704, 1.2172], 6e-3)
    return ok


def verify_airy():
    print("E. 에어리 무늬 ↔ 2J1(x)/x")
    n, R, pad = 256, 100.0, 8
    I = far_field(circle_pupil(n, R), pad)
    I /= I.max()
    prof = radial_profile(I)
    # 한 칸 = 1/(pad n) 주기/동공칸. x = k a sinθ = 2π R · (칸/(pad n))
    xx = 2 * np.pi * R * np.arange(len(prof)) / (pad * n)
    sel = xx < 14
    ref = (2 * j1(np.maximum(xx[sel], 1e-12)) / np.maximum(xx[sel], 1e-12)) ** 2
    ref[0] = 1
    ok = check("단면 ↔ [2J1(x)/x]² (x < 14)", prof[sel], ref, 3e-3)
    # 첫 영점: 단면을 보간해 찾는다
    # 첫 영점: 대칭인 동공의 FFT 는 위상 기울기만 걷어내면 실수다. 부호가 바뀌는 곳을 선형 보간으로 찾는다
    M = pad * n
    big = np.zeros((M, M))
    big[:n, :n] = circle_pupil(n, R)
    F = np.fft.fft2(big)[0, :M // 2] * np.exp(2j * np.pi * np.arange(M // 2) * (n - 1) / 2 / M)
    Fr = F.real / F.real[0]
    xg = 2 * np.pi * R * np.arange(M // 2) / M
    j = np.where((Fr[:-1] > 0) & (Fr[1:] <= 0))[0][0]
    x0 = xg[j] + Fr[j] / (Fr[j] - Fr[j + 1]) * (xg[j + 1] - xg[j])
    ok &= check("첫 영점 x = 3.832 (= 1.22 π)", x0, 3.8317, 0.01)
    pk = []
    for lo, hi in [(4.5, 6), (7.5, 9.5), (10.5, 12.5)]:
        s = (xx > lo) & (xx < hi)
        pk.append(prof[s].max())
    ok &= check("고리 봉우리 0.0175, 0.0042, 0.0016", pk, [0.0175, 0.0042, 0.0016], 2e-4)
    zeros = np.array([3.8317, 7.0156, 10.1735])
    radii = zeros * pad * n / (2 * np.pi * R)
    ok &= check("원 안 에너지 83.8 %, 91.0 %, 93.8 %", encircled_energy(I, radii), [0.838, 0.910, 0.938], 0.006)
    return ok


def verify_rayleigh():
    print("F. 레일리 기준")
    xs = np.linspace(-6, 10, 32001)
    n, R, pad = 256, 100.0, 8
    I = far_field(circle_pupil(n, R), pad)
    prof = radial_profile(I / I.max())
    xx = 2 * np.pi * R * np.arange(len(prof)) / (pad * n)
    airy = lambda x: np.interp(np.abs(x), xx, prof)           # noqa: E731  (FFT 로 만든 단면을 보간)
    sep = 3.8317
    two = airy(xs) + airy(xs - sep)
    mid = two[np.argmin(np.abs(xs - sep / 2))]
    return check("두 봉우리 사이 골 / 봉우리 = 0.735", mid / two.max(), 0.735, 3e-3)


def verify_focus():
    print("G. 초점 근처 ↔ 보른과 울프 8.8")
    u = np.linspace(-30, 30, 241)
    I_ax = np.abs(focal_field([0.0], u)[:, 0]) ** 2
    q = u / 4
    ref = np.where(np.abs(q) < 1e-9, 1, np.sin(q) / np.where(q == 0, 1, q)) ** 2
    ok = check("축 위 세기 ↔ [sin(u/4)/(u/4)]²", I_ax, ref, 1e-5)
    v = np.linspace(0.01, 12, 300)
    I_fp = np.abs(focal_field(v, [0.0])[0]) ** 2
    ok &= check("초점면 ↔ [2J1(v)/v]²", I_fp, (2 * j1(v) / v) ** 2, 1e-5)
    # 첫 영점 u = 4π → z = 4π/(k NA²) = 2λ/NA²
    uz = minimize_scalar(lambda uu: np.abs(focal_field([0.0], [uu])[0, 0]), bounds=(11, 14), method="bounded", options=dict(xatol=1e-6)).x
    lam, NA = 0.5, 0.05
    ok &= check("축 위 첫 영점 z (µm), NA 0.05 (f/10), λ 0.5 µm ↔ 2λ/NA²", uz / (2 * np.pi / lam * NA**2), 2 * lam / NA**2, 0.05)
    return ok


def verify_zone_plate():
    print("H. 존 플레이트")
    lam, f, n_open = 0.5, 1e5, 20
    z = np.concatenate([np.linspace(0.9 * f, 1.1 * f, 201), np.linspace(0.3 * f, 0.37 * f, 141), np.linspace(0.18 * f, 0.22 * f, 81)])
    I = axial_zone_plate(f, lam, n_open, z, n=100001)
    peaks = []
    for lo, hi in [(0.9, 1.1), (0.3, 0.37), (0.18, 0.22)]:
        s = (z >= lo * f) & (z <= hi * f)
        peaks.append(z[s][np.argmax(I[s])] / f)
    ok = check("초점 위치 / f ↔ 1, 1/3, 1/5", peaks, [1, 1 / 3, 1 / 5], 3e-3)
    s = (z >= 0.9 * f) & (z <= 1.1 * f)
    ok &= check("주 초점 세기 / 평면파 ↔ (2 × 20)² = 1600 (상대)", I[s].max() / 1600 - 1, 0, 0.01)
    # f/3, f/5 의 축 위 세기도 같다: 열린 존 하나가 그 거리의 존 셋(다섯)을 덮어 둘(넷)이 서로 지우고 하나가 남는다
    hi = []
    for m in (3, 5):
        zz = np.linspace(f / m * 0.995, f / m * 1.005, 201)
        hi.append(axial_zone_plate(f, lam, n_open, zz, n=100001).max() / 1600 - 1)
    ok &= check("f/3, f/5 의 축 위 세기 ↔ 1600 (상대)", hi, [0, 0], 0.02)
    return ok


def verify_grating():
    print("I. 격자")
    lam, d, b, N = 0.5, 2.0, 0.5, 200                # 줄이 많아야 포락선의 기울기가 봉우리를 밀지 않는다
    s = np.linspace(-0.9, 0.9, 360001)
    I = np.abs(slit_array_far_field(s, N, d, b, lam)) ** 2
    pk = [s[(np.abs(s - m * lam / d) < 0.02)][np.argmax(I[np.abs(s - m * lam / d) < 0.02])] for m in (1, 2, 3)]
    ok = check("주봉우리 sinθ ↔ mλ/d (m = 1, 2, 3)", pk, [lam / d, 2 * lam / d, 3 * lam / d], 1e-4)
    # 나트륨 D 선을 2차에서: 레일리 기준으로 갈리는 줄 수를 찾는다 (두 선의 골이 0.81 이하로 떨어지는 최소 N)
    l1, l2 = 0.5889950, 0.5895924
    d2 = 1.0 / 0.6                                   # 600 줄/mm. 틈 폭을 d/2 로 두면 2차가 빠진 차수가 되므로 0.3 d
    def dip(Nn):
        th0 = np.arcsin(2 * (l1 + l2) / 2 / d2)
        ds = 2 * (l2 - l1) / d2
        ss = np.linspace(np.sin(th0) - 3 * ds, np.sin(th0) + 3 * ds, 3001)
        tot = (np.abs(slit_array_far_field(ss, Nn, d2, 0.3 * d2, l1, 4)) ** 2 +
               np.abs(slit_array_far_field(ss, Nn, d2, 0.3 * d2, l2, 4)) ** 2)
        s1, s2 = 2 * l1 / d2, 2 * l2 / d2
        mid = tot[np.argmin(np.abs(ss - (s1 + s2) / 2))]
        top = tot[np.argmin(np.abs(ss - s1))]
        return mid / top
    # 레일리 기준에서 sinc² 모양 두 봉우리의 골은 8/π² = 0.811 (보른과 울프 8.6.2)
    lo, hi = 300, 800                                # dip(N) 은 N 이 커질수록 작아진다: 정수 이분법
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (lo, mid) if dip(mid) < 8 / np.pi**2 else (mid, hi)
    Nmin = hi
    ok &= check("D 선을 2차에서 가르는 최소 줄 수 ↔ 494 (Hecht)", Nmin, 494, 6)
    return ok


def verify_thin_gratings():
    print("J. 얇은 격자의 회절 효율 (Goodman 4.4)")
    x = np.arange(4096) / 4096
    ok = check("정현 진폭 격자 (m = 1): 0차 1/4, 1차 1/16", thin_grating_efficiencies(0.5 + 0.5 * np.cos(2 * np.pi * x), [0, 1]), [0.25, 0.0625], 1e-9)
    best = minimize_scalar(lambda m: -thin_grating_efficiencies(np.exp(1j * m / 2 * np.sin(2 * np.pi * x)), [1])[0], bounds=(2, 5), method="bounded")
    ok &= check("정현 위상 격자 1차 최대 효율 33.8 %", -best.fun, 0.338, 1e-3)
    ok &= check("  그때 위상 진폭 m/2 ↔ 1.84 rad (J1 의 최대)", best.x / 2, 1.8412, 2e-3)
    ok &= check("  같은 결과 ↔ J1(m/2)²", -best.fun, jv(1, best.x / 2) ** 2, 1e-6)
    return ok


def verify_slit():
    print("K. 틈 하나")
    lam, b = 0.5, 5.0
    s = np.linspace(0.0, 0.32, 64001)
    I = np.abs(slit_array_far_field(s, 1, 0, b, lam, pts_per_slit=400)) ** 2
    pk = []
    for m in (1, 2):
        sel = (s > (m + 0.2) * lam / b) & (s < (m + 0.8) * lam / b)
        pk.append(I[sel].max())
    ok = check("곁봉우리 0.0472, 0.0165", pk, [0.04718, 0.01648], 2e-4)
    z1 = s[(s > 0.05) & (s < 0.15)][np.argmin(I[(s > 0.05) & (s < 0.15)])]
    ok &= check("첫 영점 sinθ = λ/b", z1, lam / b, 2e-4)
    return ok


def main():
    results = [f() for f in (verify_propagation, verify_axial_aperture, verify_arago, verify_edge, verify_airy, verify_rayleigh,
                             verify_focus, verify_zone_plate, verify_grating, verify_thin_gratings, verify_slit)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

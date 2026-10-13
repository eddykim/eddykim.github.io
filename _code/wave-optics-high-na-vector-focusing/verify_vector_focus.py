"""vector_focus.py 의 두 경로와 문헌 수치를 대조한다.

실행: python verify_vector_focus.py

  A. 리처즈-울프 적분         ↔  벡터 평면파의 직접 합 (선편광·방사·방위 편광, 초점면과 초점 밖의 여러 점)
  B. 낮은 NA                  ↔  스칼라 에어리 무늬 [2J1(v)/v]², v = 2π NA ρ/λ (3–4편)
  C. 세로 성분 Ez             ↔  낮은 NA 에서 사라지고(보른과 울프 8.4 의 스칼라 극한) 그 상대 세기는 NA² 에 비례
  D. 원편광                   ↔  초점의 세기가 회전 대칭
  E. 반점 넓이                ↔  NA 0.9, 고리 동공(NA 0.85–0.9): 방사 편광 0.16 λ², 선편광 0.26 λ² (Dorn, Quabis & Leuchs 2003, Zhan 2009 인용)
  F. 두 평면파의 간섭         ↔  s 편광은 가시도 1, p 편광은 |cos 2θ| (편광이 직교하면 간섭하지 않는다, 1편의 프레넬-아라고 법칙)
"""
import numpy as np
from scipy.special import j1

from vector_focus import (half_max_area, plane_wave_sum, rw_azimuthal, rw_intensity_map, rw_linear, rw_radial,
                          two_beam_fringe_visibility)


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<66s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def match_up_to_constant(a, b):
    """a ≈ c·b 인 복소 상수 c 를 최소제곱으로 맞추고, 남은 차이를 a 의 크기로 나눈 값."""
    a, b = np.ravel(a), np.ravel(b)
    c = np.vdot(b, a) / np.vdot(b, b)
    return np.max(np.abs(a - c * b)) / np.max(np.abs(a))


def verify_two_paths():
    print("A. 리처즈-울프 적분 ↔ 벡터 평면파 직접 합")
    NA = 0.9
    pts = [(0.0, 0.0, 0.0), (0.2, 0.0, 0.0), (0.0, 0.3, 0.0), (0.25, 0.25, 0.0), (0.4, -0.1, 0.3), (0.1, 0.2, -0.5)]
    rho = np.array([np.hypot(x, y) for x, y, _ in pts])
    phi = np.array([np.arctan2(y, x) for x, y, _ in pts])
    ok = True
    # 선편광
    rw = []
    for (x, y, z), r, p in zip(pts, rho, phi):
        Ex, Ey, Ez = rw_linear(np.array([r]), np.array([p]), z, NA)
        rw.append([Ex[0], Ey[0], Ez[0]])
    pw = plane_wave_sum(pts, NA, "x")
    ok &= check("선편광 (x): 세 성분 × 여섯 점, 상수배를 맞춘 뒤 남는 차이", match_up_to_constant(pw, np.array(rw)), 0, 2e-3)
    # 방사 편광
    rw = []
    for (x, y, z), r, p in zip(pts, rho, phi):
        Er, Ez = rw_radial(np.array([r]), z, NA)
        rw.append([Er[0] * np.cos(p), Er[0] * np.sin(p), Ez[0]])
    pw = plane_wave_sum(pts, NA, "radial")
    ok &= check("방사 편광", match_up_to_constant(pw, np.array(rw)), 0, 2e-3)
    rw = []
    for (x, y, z), r, p in zip(pts, rho, phi):
        Ep = rw_azimuthal(np.array([r]), z, NA)[0]
        rw.append([-Ep * np.sin(p), Ep * np.cos(p), 0.0])
    pw = plane_wave_sum(pts, NA, "azimuthal")
    ok &= check("방위 편광", match_up_to_constant(pw, np.array(rw)), 0, 2e-3)
    return ok


def verify_low_na():
    print("B. 낮은 NA ↔ 스칼라 에어리 무늬")
    NA = 0.1
    r = np.linspace(0, 1.2 / NA, 241)
    Ex, Ey, Ez = rw_linear(r, np.zeros_like(r), 0.0, NA)
    I = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
    v = 2 * np.pi * NA * np.maximum(r, 1e-9)
    airy = (2 * j1(v) / v) ** 2
    return check("NA 0.1 초점면 단면 (x 방향) ↔ [2J1(v)/v]²", I / I[0], airy, 3e-3)


def verify_ez():
    print("C. 세로 성분 Ez")
    r = np.linspace(0, 6, 1201)
    fr = []
    for NA in (0.05, 0.1, 0.2):
        rr = r / NA * 0.1
        Ex, Ey, Ez = rw_linear(rr, np.zeros_like(rr), 0.0, NA)
        fr.append((np.abs(Ez) ** 2).max() / (np.abs(Ex[0]) ** 2))
    ok = check("|Ez|² 최대 / 중심 세기: NA 0.05 에서 1 % 미만", fr[0] < 0.01, 1, 0.5)
    ok &= check("NA 를 두 배로 하면 4 배 (NA² 비례): 0.1/0.05, 0.2/0.1", [fr[1] / fr[0], fr[2] / fr[1]], [4, 4], 0.1)
    return ok


def circular_intensity(r, phis, NA):
    """원편광 (x̂ + iŷ)/√2 입력의 초점 세기를 극좌표 점에서 직접 계산한다."""
    ax = rw_linear(r, phis, 0.0, NA)                                      # x 편광 입력
    q = rw_linear(r, phis - np.pi / 2, 0.0, NA)                           # y 편광 입력 = x 편광을 90° 돌린 것
    ay = (-q[1], q[0], q[2])                                              # 성분도 90° 돌린다
    return sum(np.abs((u + 1j * v) / np.sqrt(2)) ** 2 for u, v in zip(ax, ay))


def verify_circular():
    print("D. 원편광의 초점은 회전 대칭")
    phis = np.linspace(0, 2 * np.pi, 73)[:-1]
    I = circular_intensity(np.full_like(phis, 0.3), phis, 0.9)
    ok = check("반지름 0.3λ 원 위 72 점의 세기 상대 편차", I.std() / I.mean(), 0, 1e-6)
    x = np.linspace(-1, 1, 121)                                           # 그림에 쓰는 지도 함수도 같은 값을 주는지
    Imap = sum(rw_intensity_map("circular", x, 0.9))
    I0 = circular_intensity(np.array([0.0]), np.array([0.0]), 0.9)[0]
    ok &= check("  지도 함수의 중심 세기 ↔ 같은 계산 (상대)", Imap[60, 60] / I0 - 1, 0, 1e-3)
    return ok


def verify_dorn():
    print("E. 반점 넓이 (Dorn 등 2003)")
    x = np.linspace(-1.2, 1.2, 481)
    dx = x[1] - x[0]
    areas = {}
    for kind in ("radial", "x"):
        I = sum(rw_intensity_map(kind, x, 0.9, NA_inner=0.85))
        areas[kind] = half_max_area(I, dx)
    ok = check("고리 동공 NA 0.85–0.9, 방사 편광 반치 넓이 ↔ 0.16 λ²", areas["radial"], 0.16, 0.01)
    ok &= check("  같은 조건 선편광 ↔ 0.26 λ²", areas["x"], 0.26, 0.015)
    full = {k: half_max_area(sum(rw_intensity_map(k, x, 0.9)), dx) for k in ("radial", "x")}
    print(f"      전체 동공 NA 0.9: 방사 {full['radial']:.3f} λ², 선편광 {full['x']:.3f} λ²")
    return ok


def verify_two_beam():
    print("F. 두 평면파의 간섭")
    th = np.radians([10, 30, 45, 60, 70])
    vs = [two_beam_fringe_visibility(t, "s") for t in th]
    vp = [two_beam_fringe_visibility(t, "p") for t in th]
    ok = check("s 편광 가시도 = 1", vs, 1, 1e-9)
    ok &= check("p 편광 가시도 = |cos 2θ|", vp, np.abs(np.cos(2 * th)), 1e-6)
    return ok


def main():
    results = [f() for f in (verify_two_paths, verify_low_na, verify_ez, verify_circular, verify_dorn, verify_two_beam)]
    print(f"\n{sum(results)}/{len(results)} 항목 통과")


if __name__ == "__main__":
    main()

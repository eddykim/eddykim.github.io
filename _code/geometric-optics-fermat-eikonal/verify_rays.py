"""rays.py 를 독립된 경로들과 교차검증한다.

실행: python verify_rays.py

  A. 광선 방정식 수치적분  ↔  n^2 이 선형·포물선인 매질의 닫힌 해
  B. 적분 중 보존량       ↔  층상 매질의 p_x (= n cos θ), |p| = n
  C. 스넬 법칙            ↔  광로길이를 경계면 통과점에 대해 직접 최소화
  D. 광선                ↔  경로를 잘게 나눠 광로길이를 직접 최소화 (페르마)
  E. 헤시안의 음의 고유값 수 ↔  광선족에서 따로 찾은 켤레점 수 (모스 지표)
  F. 파면과 광선의 직교   ↔  같은 광로길이 점을 이은 곡선의 접선
  G. 신기루 임계각       ↔  불변량 n cos θ 가 주는 arccos(n_지면 / n_눈)
  H. 싱글렛 2D 추적       ↔  opticore 3D 추적, Optiland 기록값, 두꺼운 렌즈 근축식
  I. 원형 거울 화선     ↔  네프로이드 닫힌 식
  J. 위상자 합의 크기     ↔  정상위상 근사
  K. 종방향 구면수차      ↔  입사 높이의 제곱에 비례 (3차 수차, Fowles 10.2 / Hecht 6.3)
  L. 최소 착란원          ↔  가장자리 광선과 화선의 교점 (Hecht 6.3)
"""
import json
import os
import subprocess

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar

from generate_figures import PH_A, PH_H, PH_LAMBDA, compute_wavefronts, phasor_sum
from rays import (Singlet, envelope_of_lines, least_confusion, n2_linear, n2_parabolic, opl_hessian,
                  opl_polyline, ring_reflection, road_air, sample, trace_ray)

np.seterr(invalid="ignore")


def check(name, got, want, tol):
    err = np.max(np.abs(np.asarray(got, dtype=float) - np.asarray(want, dtype=float)))
    ok = bool(err < tol)
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<58s} 오차 {err:.3e} (허용 {tol:.0e})")
    return ok


def verify_closed_forms():
    print("A. 광선 방정식 ↔ 닫힌 해")
    ok = True
    a, b, th = 1.0, 0.5, np.deg2rad(20)
    s, x, y, *_ = sample(trace_ray(n2_linear(a, b), (0, 0), th, 3.0))
    beta = np.sqrt(a) * np.cos(th)
    ok &= check("n^2 = a + b y : 포물선 y = x tanθ + b x^2 / 4β^2", y,
                np.tan(th) * x + b / (4 * beta**2) * x**2, 1e-8)
    n0, al, th = 1.5, 0.5, np.deg2rad(10)
    s, x, y, *_ = sample(trace_ray(n2_parabolic(n0, al), (0, 0), th, 20.0))
    k = n0 * al / (n0 * np.cos(th))
    ok &= check("n^2 = n0^2(1 - α^2 y^2) : 사인곡선", y, np.tan(th) / k * np.sin(k * x), 1e-8)
    return ok


def verify_invariants():
    print("B. 보존량")
    M = n2_parabolic(1.5, 0.5)
    s, x, y, px, py, L = sample(trace_ray(M, (0, 0.3), np.deg2rad(25), 30.0))
    ok = check("층상 매질에서 p_x 일정", px, px[0], 1e-10)
    ok &= check("|p| = n", np.hypot(px, py), M.n(x, y), 1e-10)
    return ok


def verify_snell_from_fermat():
    print("C. 스넬 법칙 ↔ 페르마 (평면 경계)")
    n1, n2 = 1.0, 1.5
    f = lambda x: n1 * np.hypot(x, 1.0) + n2 * np.hypot(2.0 - x, 1.0)
    x = minimize_scalar(f, bracket=(0.0, 1.0, 2.0), method="brent", tol=1e-12).x
    s1, s2 = x / np.hypot(x, 1.0), (2 - x) / np.hypot(2 - x, 1.0)
    return check("최소점에서 n1 sinθ1 = n2 sinθ2", n1 * s1, n2 * s2, 1e-8)


def shoot(M, X, Y, lo=-0.8, hi=0.8, n=801):
    """(0,0) 에서 출발해 (X, Y) 를 지나는 광선의 출발각을 전부 찾는다."""
    def y_at(t):
        sol = trace_ray(M, (0, 0), t, 3 * X, stop=lambda x, y: x - X)
        return sol.y_events[0][0][1] - Y
    ts = np.linspace(lo, hi, n)
    v = np.array([y_at(t) for t in ts])
    return [brentq(y_at, ts[i], ts[i + 1], xtol=1e-14)
            for i in range(n - 1) if np.sign(v[i]) != np.sign(v[i + 1])]


def ray_on_grid(M, th, X, xs):
    sol = trace_ray(M, (0, 0), th, 3 * X, stop=lambda x, y: x - X)
    s, x, y, *_ = sample(sol, 40000)
    return np.interp(xs, x, y), sol


def verify_fermat_path():
    print("D. 광선 ↔ 광로길이 직접 최소화 (포물선 GRIN, 끝점 (4, 0.5))")
    M = n2_parabolic(1.5, 0.5)
    X, Y = 4.0, 0.5
    xs = np.linspace(0, X, 161)
    th = min(shoot(M, X, Y), key=abs)
    yr, _ = ray_on_grid(M, th, X, xs)
    f = lambda v: opl_polyline(M, xs, np.r_[0.0, v, Y])

    def grad(v, h=1e-7):
        g = np.empty_like(v)
        for i in range(len(v)):
            e = np.zeros_like(v)
            e[i] = h
            g[i] = (f(v + e) - f(v - e)) / (2 * h)
        return g

    # 출발점은 두 끝점을 잇는 직선. 광선에 대한 정보는 주지 않는다
    r = minimize(f, np.linspace(0, Y, 161)[1:-1], jac=grad, method="BFGS",
                 options={"gtol": 1e-10, "maxiter": 5000})
    ok = check("최소화 경로 ↔ 광선 [높이]", r.x, yr[1:-1], 1e-4)
    ok &= check("광로길이", r.fun, opl_polyline(M, xs, yr), 1e-8)
    return ok


def conjugate_points(M, th, X, dth=1e-5):
    """출발각을 조금 바꾼 이웃 광선과의 간격 dy/dθ 가 0 이 되는 횟수."""
    xs = np.linspace(1e-3, X, 20000)
    y1, _ = ray_on_grid(M, th - dth, X, xs)
    y2, _ = ray_on_grid(M, th + dth, X, xs)
    dy = y2 - y1
    return int(np.sum(np.sign(dy[1:]) != np.sign(dy[:-1])))


def verify_morse():
    print("E. 헤시안 음의 고유값 수 ↔ 켤레점 수 (광선은 켤레점을 지나면 최소가 아니다)")
    M = n2_parabolic(1.5, 0.5)
    ok = True
    rows = []
    for X, Y in ((4.0, 0.5), (9.0, 0.5)):
        xs = np.linspace(0, X, 81)
        for th in shoot(M, X, Y):
            yr, sol = ray_on_grid(M, th, X, xs)
            neg = int(np.sum(np.linalg.eigvalsh(opl_hessian(M, xs, yr, h=1e-4)) < 0))
            cp = conjugate_points(M, th, X)
            rows.append((X, np.rad2deg(th), sol.y[4][-1], neg, cp))
            ok &= check(f"끝점 x={X:.0f}, 출발각 {np.rad2deg(th):+6.2f}°: 음의 고유값 {neg}개",
                        neg, cp, 0.5)
    for X, th, L, neg, cp in rows:
        print(f"      x={X:.0f}  출발각 {th:+7.3f}°  광로길이 {L:.5f}  켤레점 {cp}개")
    return ok


def verify_wavefront_orthogonality():
    print("F. 파면 ⟂ 광선 (그림2 의 매질)")
    medium, src, rays, fronts = compute_wavefronts()
    errs = []
    for f_idx, lv in enumerate(np.arange(0.1, 2.0, 0.1)):
        pts = fronts[f_idx]
        for i in range(1, len(rays) - 1):
            if np.any(np.isnan(pts[i - 1:i + 2])):
                continue
            tang = pts[i + 1] - pts[i - 1]
            s, x, y, px, py, L = rays[i]
            d = np.array([np.interp(lv, L, px), np.interp(lv, L, py)])
            errs.append(abs(np.dot(tang, d)) / (np.linalg.norm(tang) * np.linalg.norm(d)))
    # 이웃 광선 사이 차분이라 광선 간격(5°)에 비례하는 오차가 남는다
    return check("|cos(광선, 파면 접선)| 최대", max(errs), 0.0, 2e-2)


def verify_mirage():
    print("G. 신기루 임계각 ↔ n cosθ 불변량")
    A = road_air()
    eye = (0.0, 1.5)

    def hits(dep):
        sol = trace_ray(A, eye, -dep, 2000.0, stop=lambda x, y: y, max_step=0.5)
        return sol.status == 1

    lo, hi = np.deg2rad(0.2), np.deg2rad(0.45)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if hits(mid) else (mid, hi)
    want = np.arccos(A.n(0, 0.0) / A.n(0, eye[1]))
    ok = check("임계 내려다봄 각 [deg]", np.rad2deg(0.5 * (lo + hi)), np.rad2deg(want), 2e-4)
    dn = A.n(0, eye[1]) - A.n(0, 0.0)
    print(f"      Δn = {dn:.3e},  임계각 {np.rad2deg(want):.4f}° = {want * 1e3:.3f} mrad,"
          f"  근사 sqrt(2Δn) = {np.sqrt(2 * dn) * 1e3:.3f} mrad")
    return ok


OPTICORE_SNIPPET = r"""
import json, sys, numpy as np
from opticore.geometry import Sphere
from opticore.interactions import Refractive
from opticore.materials import ConstantMaterial, SellmeierMaterial
from opticore.rays import RayBundle
from opticore.surface import Surface
from opticore.system import OpticalSystem
g = SellmeierMaterial(B=(1.03961212, 0.231792344, 1.01046945),
                      C=(0.00600069867, 0.0200179144, 103.560653))
a = ConstantMaterial(1.0002778)
s = OpticalSystem(ambient=a)
s.add(Surface(Sphere(1000.0), g, Refractive(), position=(0, 0, 0)))
s.add(Surface(Sphere(-1000.0), a, Refractive(), position=(0, 0, 100.0)))
h = np.array(json.loads(sys.argv[1]))
r = RayBundle(position=np.stack([0 * h, h, -300 + 0 * h], 1),
              direction=np.tile([0, 0, 1.0], (len(h), 1)), wavelength_um=0.75)
f = s.trace(r)[-1]
p, d = np.asarray(f.position), np.asarray(f.direction)
print(json.dumps((p[:, 2] - p[:, 1] * d[:, 2] / d[:, 1]).tolist()))
"""

OPTICORE_PYTHON = os.environ.get("OPTICORE_PYTHON", os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
    "repo_opticore", "opticore", ".venv", "bin", "python"))


def verify_singlet():
    print("H. 싱글렛 (N-BK7, R=±1000 mm, t=100 mm, 750 nm)")
    lens = Singlet()
    ok = True
    # 두꺼운 렌즈 근축식: 면 굴절력 φ1, φ2 와 두께로 후초점거리를 구한다
    n, na, t = lens.n_glass, lens.n_air, lens.thickness
    p1, p2 = (n - na) / lens.r1, (na - n) / lens.r2
    phi = p1 + p2 - p1 * p2 * t / n
    z_par = t + na * (1 - p1 * t / n) / phi
    ok &= check("광축 근처 광선(h=1e-4 mm) ↔ 두꺼운 렌즈 근축 초점 [mm]",
                lens.axis_crossing([1e-4])[0], z_par, 1e-6)

    heights = [1.0, 30.0, 60.0, 90.0, 100.0]
    mine = lens.axis_crossing(heights)
    # Optiland 오라클 기록값 (repo_opticore/design_docs/PHASE0_OPTILAND_STUDY.md 절 B, Py=h/100)
    optiland = {30.0: 1059.444, 60.0: 1055.200, 90.0: 1048.055, 100.0: 1045.011}
    ok &= check("↔ Optiland 기록값 [mm] (기록이 소수 셋째 자리)",
                [mine[heights.index(h)] for h in optiland], list(optiland.values()), 1e-2)

    py = os.path.abspath(OPTICORE_PYTHON)
    if os.path.exists(py):
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1")
        out = subprocess.run([py, "-c", OPTICORE_SNIPPET, json.dumps(heights)],
                             capture_output=True, text=True, env=env, check=True)
        ok &= check("↔ opticore 3D 실광선 추적 [mm]", mine, json.loads(out.stdout), 1e-7)
    else:
        print(f"  [SKIP] opticore 가 없어 대조를 건너뛴다 ({py})")
    print(f"      근축 초점 {z_par:.4f} mm, 가장자리 광선 {mine[-1]:.4f} mm,"
          f" 종방향 구면수차 {mine[-1] - z_par:.3f} mm")
    return ok


def verify_nephroid():
    print("I. 원형 거울 화선 ↔ 네프로이드")
    h = np.linspace(-0.95, 0.95, 20001)
    hit, d = ring_reflection(1.0, h)
    ex, ey = envelope_of_lines(hit[:, 0], hit[:, 1], d[:, 1] / d[:, 0], h)
    t = np.linspace(-np.pi / 2, np.pi / 2, 200001)
    nx, ny = 0.25 * (3 * np.cos(t) - np.cos(3 * t)), 0.25 * (3 * np.sin(t) - np.sin(3 * t))
    sel = slice(10, -10, 97)
    dist = [np.min(np.hypot(nx - x, ny - y)) for x, y in zip(ex[sel], ey[sel])]
    ok = check("포락선 점에서 네프로이드까지 최대 거리 (반경 1)", max(dist), 0.0, 1e-4)
    ok &= check("첨점 위치 = R/2 (구면거울의 근축 초점)", ex[len(h) // 2], 0.5, 1e-6)
    return ok


def verify_phasor():
    print("J. 위상자 합 ↔ 정상위상 근사")
    x = np.linspace(-0.6, 0.6, 600001)
    ph, csum = phasor_sum(x)
    k = 2 * np.pi / PH_LAMBDA
    f2 = 2 * PH_H**2 / (PH_A**2 + PH_H**2) ** 1.5         # 광로길이의 x=0 2차 미분
    want = np.sqrt(2 * np.pi / (k * f2))
    ok = check("|합| ↔ sqrt(2π / k f'') (끝점 진동 몫만큼 차이)", abs(csum[-1]) / want, 1.0, 5e-2)
    near = np.abs(x) < 0.08
    part = abs(np.trapezoid(ph[near], x[near]))
    print(f"      |x| < 0.08 mm 구간(전체 폭의 {0.16 / 1.2:.0%})이 만드는 크기 / 전체 = "
          f"{part / abs(csum[-1]):.3f},  1차 프레넬 영역 반폭 = {np.sqrt(PH_LAMBDA / f2):.4f} mm")
    return ok


def verify_sa_scaling():
    print("K. 종방향 구면수차 ∝ h^2")
    lens = Singlet()
    z_par = lens.axis_crossing([1e-4])[0]
    h = np.array([10.0, 30.0, 60.0, 100.0])
    ratio = (lens.axis_crossing(h) - z_par) / h**2
    for hi, r in zip(h, ratio):
        print(f"      h = {hi:5.1f} mm   Δz / h^2 = {r:.4e} /mm")
    # 5차 이상의 항이 붙으므로 가장자리에서 조금 벗어난다
    return check("h=10 대비 Δz/h^2 의 상대 변화 (h ≤ 100 mm)", ratio / ratio[0], 1.0, 2e-2)


def verify_least_confusion():
    print("L. 최소 착란원 ↔ 가장자리 광선과 화선의 교점")
    lens = Singlet()
    z_lc, d_lc = least_confusion(lens)
    hd = np.linspace(1e-3, 100, 40001)
    p, d = lens.trace(hd)
    lz, ly = envelope_of_lines(p[:, 0], p[:, 1], d[:, 1] / d[:, 0], hd)
    pm, dm = lens.trace([100.0])
    o = np.argsort(lz)
    # 위에서 내려오는 가장자리 광선이 반대쪽(아래) 화선 y = -ly 와 만나는 z
    f = lambda z: pm[0, 1] + (z - pm[0, 0]) * dm[0, 1] / dm[0, 0] + np.interp(z, lz[o], ly[o])
    z_x = brentq(f, 1045.5, 1060.8)
    print(f"      번짐 지름 최소화: z = {z_lc:.4f} mm, 지름 {d_lc:.3f} mm")
    return check("번짐 최소 위치 ↔ 교점 [mm]", z_lc, z_x, 1e-4)


def main():
    results = [
        verify_closed_forms(), verify_invariants(), verify_snell_from_fermat(),
        verify_fermat_path(), verify_morse(), verify_wavefront_orthogonality(),
        verify_mirage(), verify_singlet(), verify_nephroid(), verify_phasor(),
        verify_sa_scaling(), verify_least_confusion(),
    ]
    print("\n전체:", "통과" if all(results) else "실패 있음")


if __name__ == "__main__":
    main()

"""retarder.py 를 서로 독립적인 두 경로로 교차검증한다.

실행: python verify_retarder.py

타원계측기 1편에서 참고문헌 식을 옮겨 적다 생긴 오독을 이 방식으로 잡았다.
여기서도 같은 원칙을 지킨다 — 한쪽 결과를 다른 쪽에서 유도하지 않는다.

  겹친 위상지연자의 유효 지연량
    경로 A: 존스 행렬 곱 -> 고유값 -> delta_eff = 2*arccos(Re(tr)/2)
    경로 B: 뮬러 행렬 곱 -> 우하단 3x3 회전 -> delta_eff = arccos((tr-1)/2)
  전반사 위상차
    경로 A: 프레넬 반사계수(복소수)의 편각 차
    경로 B: 닫힌 식 tan(d/2) = cos(t) sqrt(sin^2 t - n^2) / sin^2 t
  단일 파장판
    설계 두께 -> 설계 파장에서 지연량이 정확히 2*pi*order 인가
"""
import numpy as np

from retarder import (DN_QUARTZ, critical_angle, effective_fast_axis_jones,
                      effective_retardance_jones, effective_retardance_mueller,
                      jones_retarder, jones_stack, jones_to_stokes,
                      mueller_stack, plate_thickness, retardance,
                      stokes_to_jones, tir_phase_difference,
                      tir_phase_difference_closed)

TOL = 1e-10
rng = np.random.default_rng(0)


def check(name, got, want, tol=TOL):
    err = np.max(np.abs(np.asarray(got) - np.asarray(want)))
    ok = err < tol
    print(f"  [{'OK  ' if ok else 'FAIL'}] {name:<52s} 최대오차 {err:.3e}")
    return ok


def verify_single_plate():
    """설계 파장에서 지연량이 의도한 차수와 맞는가."""
    print("단일 파장판 — 설계 두께와 지연량")
    ok = True
    for order in (0.25, 0.5, 1.25, 2.25):
        for lam in (266e-9, 550e-9, 633e-9):
            d = plate_thickness(order, lam)
            ok &= check(f"order={order:5.2f} lambda={lam*1e9:5.1f}nm -> delta=2pi*order",
                        retardance(lam, d), 2 * np.pi * order)
    # Handbook 이 든 수치: 석영 적색광 영차 4분의1파장판은 약 18 um 다.
    d = plate_thickness(0.25, 633e-9) * 1e6
    print(f"  [참고] 석영 영차 lambda/4 @633nm 두께 = {d:.1f} um  (Handbook: 약 18 um)")
    ok &= check("  Handbook 의 18 um 와 일치", d, 18.0, tol=0.6)
    return ok


def random_polarized_stokes():
    """푸앵카레 구면 위의 무작위 완전편광 상태."""
    v = rng.normal(size=3)
    return np.array([1.0, *(v / np.linalg.norm(v))])


def verify_stack_two_paths():
    """겹친 지연자를 존스 경로와 뮬러 경로가 같게 다루는가.

    유효 지연량만 비교하면 arccos 치역(0~pi) 때문에 delta > pi 에서 2pi-delta 로
    접힌다. 그래서 두 가지로 나눠 검사한다.
      (a) 작용 — 같은 편광 상태를 넣었을 때 나오는 스토크스 벡터가 같은가
      (b) cos(delta_eff) — 접힘에 무관한 양이므로 그대로 비교할 수 있다
    """
    print("겹친 위상지연자 — 존스 경로 vs 뮬러 경로")
    ok = True
    for n_plate in (1, 2, 3):
        for _ in range(4):
            deltas = rng.uniform(0.1, 2 * np.pi - 0.1, n_plate)
            thetas = rng.uniform(-np.pi / 2, np.pi / 2, n_plate)
            J, M = jones_stack(deltas, thetas), mueller_stack(deltas, thetas)
            S = random_polarized_stokes()
            ok &= check(f"판 {n_plate}장  편광 상태 변환 일치",
                        jones_to_stokes(J @ stokes_to_jones(S)), M @ S)
            ok &= check(f"판 {n_plate}장  cos(delta_eff) 일치",
                        np.cos(effective_retardance_jones(J)),
                        np.cos(effective_retardance_mueller(M)))
    return ok


def verify_stack_sanity():
    """겹침 계산이 자명한 경우를 제대로 재현하는가."""
    print("겹친 위상지연자 — 자명한 경우")
    ok = True
    # 같은 방위각이면 지연량이 그냥 더해진다.
    d1, d2, th = 0.7, 1.1, 0.3
    ok &= check("같은 축 두 장 -> 지연량 합",
                effective_retardance_jones(jones_stack([d1, d2], [th, th])), d1 + d2)
    # 한 장짜리는 원래 값과 축이 그대로 나와야 한다.
    for _ in range(3):
        d = rng.uniform(0.2, np.pi)
        th = rng.uniform(-np.pi / 4, np.pi / 4)
        J = jones_retarder(d, th)
        ok &= check("한 장  delta_eff = delta", effective_retardance_jones(J), d)
        ok &= check("한 장  빠른 축 = theta", effective_fast_axis_jones(J), th)
    # 축이 직교한 같은 지연량 두 장은 서로 상쇄된다.
    d = 1.234
    ok &= check("직교한 같은 두 장 -> 지연량 0",
                effective_retardance_jones(jones_stack([d, d], [0.0, np.pi / 2])), 0.0)
    return ok


def verify_tir():
    """전반사 위상차 — 프레넬 계수 경로 vs 닫힌 식."""
    print("전반사 위상차 — 프레넬 계수 vs 닫힌 식")
    ok = True
    for n_in in (1.50, 1.51, 1.52, 1.75):
        tc = critical_angle(n_in)
        th = np.linspace(tc + 1e-4, np.pi / 2 - 1e-4, 400)
        ok &= check(f"n={n_in:.2f}  두 경로 일치",
                    tir_phase_difference(th, n_in),
                    tir_phase_difference_closed(th, n_in), tol=1e-9)
        # 전반사이므로 반사율의 크기는 1 이어야 한다.
        st = n_in * np.sin(th)
        ct = np.sqrt(1 - st ** 2 + 0j)
        rs = (n_in * np.cos(th) - ct) / (n_in * np.cos(th) + ct)
        ok &= check(f"n={n_in:.2f}  |r_s| = 1", np.abs(rs), np.ones_like(th), tol=1e-12)
    return ok


def report_rhomb():
    """프레넬 롬이 실제로 가능한 설계인지 확인한다(Handbook 그림 4.13)."""
    print("프레넬 롬 — 반사 1회당 22.5도가 나오는 입사각")
    n = 1.51
    tc = critical_angle(n)
    th = np.linspace(tc + 1e-6, np.pi / 2 - 1e-6, 200000)
    d = np.abs(np.rad2deg(tir_phase_difference(th, n)))   # 규약상 음수라 크기를 본다
    print(f"  임계각 {np.rad2deg(tc):.2f}°,  최대 위상차 크기 {d.max():.2f}° "
          f"(@ {np.rad2deg(th[d.argmax()]):.2f}°)")
    hits = np.rad2deg(th[np.where(np.diff(np.sign(d - 22.5)))[0]])
    print(f"  22.5° 를 주는 입사각: {np.round(hits, 2)}")
    ok = d.max() > 22.5 and len(hits) == 2
    print(f"  [{'OK  ' if ok else 'FAIL'}] 반사 4회로 90도를 만들 수 있다")
    return ok


if __name__ == "__main__":
    results = [verify_single_plate(), verify_stack_sanity(),
               verify_stack_two_paths(), verify_tir(), report_rhomb()]
    print("\n전체 결과:", "모두 통과" if all(results) else "실패 항목 있음")
    raise SystemExit(0 if all(results) else 1)

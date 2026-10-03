"""3편 계산의 교차검증. 문헌의 주장을, 그 주장에서 유도하지 않은 경로로 확인한다.

  1. PIE = ER         조명이 마스크이고 위치가 하나, beta = 1, alpha -> 0 이면 PIE 는 Fienup 의 ER 과 같다
                      (Rodenburg & Faulkner 2004). 출사파 P*O 를 2편 방식의 ER 반복과 비교한다.
  2. ePIE = 경사하강  위치별 오차 L = || |F(P O)| - sqrt(I) ||^2 의 해석적 기울기 2 conj(P)(psi - psi') 를
                      유한차분과 대조하고, ePIE 의 시편 갱신이 그 기울기 방향 1/(2 max|P|^2) 걸음임을 본다.
  3. 모호성           (cO, P/c) 와 (O e^{ig.r}, P e^{-ig.r}) 가 모든 회절 패턴을 바꾸지 않는다.
  4. 고정점           정답 (O, P) 는 PIE 와 ePIE 가 고치지 않는다.
  5. PIE 가중        보정 중 |P|/|P|max 만큼만 시편에 반영한다. 최대점에서는 나눗셈 전체, |P| -> 0 이면 0.
"""
import numpy as np

from ptycho_core import (fft2c, circular_aperture, pinhole_probe, make_object, raster, object_shape,
                         intensities, fourier_project, pie, epie, pie_weight)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
n = 64

# 1. PIE = ER
mask = circular_aperture(n, 16).real
obj1 = make_object((n, n))
I1 = np.abs(fft2c(mask * obj1)) ** 2
o = np.exp(2j * np.pi * rng.random((n, n)))
g = mask * o
worst = 0.0
for _ in range(10):
    o = pie(I1[None], np.zeros((1, 2), int), mask.astype(complex), (n, n), 1, alpha=1e-12, obj0=o)
    g = mask * fourier_project(g, np.sqrt(I1))
    worst = max(worst, np.abs(mask * o - g).max())
check("마스크 조명·위치 하나의 PIE = ER (10회)", worst < 1e-8, f"출사파 최대 차이 {worst:.1e}")

# 2. ePIE = 경사하강
P = pinhole_probe()
O = make_object((n, n))
sI = np.sqrt(np.abs(fft2c(P * make_object((n, n), phase_max=2.0))) ** 2)   # 다른 시편의 측정
psi = P * O
grad = 2 * np.conj(P) * (psi - fourier_project(psi, sI))
L = lambda oo: np.sum((np.abs(fft2c(P * oo)) - sI) ** 2)
h = 1e-6
rel = 0.0
for _ in range(3):
    d = rng.normal(size=O.shape) + 1j * rng.normal(size=O.shape)
    fd = (L(O + h * d) - L(O - h * d)) / (2 * h)
    rel = max(rel, abs(fd - np.real(np.vdot(grad, d))) / abs(fd))
check("위치별 오차의 해석적 기울기 = 유한차분", rel < 1e-6, f"상대 차이 최대 {rel:.1e}")
o_epie, _ = epie(sI[None] ** 2, np.zeros((1, 2), int), P, (n, n), 1, obj0=O)
o_grad = O - grad / (2 * np.max(np.abs(P)) ** 2)
d = np.abs(o_epie - o_grad).max()
check("ePIE 시편 갱신 = 기울기 방향 1/(2 max|P|^2) 걸음", d < 1e-12, f"최대 차이 {d:.1e}")

# 3. 모호성
pos = raster(10, span=40)
shape = object_shape(pos)
obj = make_object(shape)
I = intensities(obj, P, pos)
c = 1.7 * np.exp(0.6j)
d_scale = np.abs(intensities(obj * c, P / c, pos) - I).max() / I.max()
Y, X = np.indices(shape)
gx, gy = 0.11, -0.07
ramp_o = obj * np.exp(1j * (gx * X + gy * Y))
yy, xx = np.indices((n, n))
ramp_p = P * np.exp(-1j * (gx * xx + gy * yy))
d_ramp = np.abs(intensities(ramp_o, ramp_p, pos) - I).max() / I.max()
check("(cO, P/c) 는 세기를 바꾸지 않는다", d_scale < 1e-12, f"최대 상대 차이 {d_scale:.1e}")
check("(O e^{ig.r}, P e^{-ig.r}) 는 세기를 바꾸지 않는다", d_ramp < 1e-12, f"최대 상대 차이 {d_ramp:.1e}")

# 4. 고정점
o_fix = pie(I, pos, P, shape, 1, obj0=obj)
o2, p2 = epie(I, pos, P, shape, 3, obj0=obj)
d4 = max(np.abs(o_fix - obj).max(), np.abs(o2 - obj).max(), np.abs(p2 - P).max())
check("정답은 PIE·ePIE 의 고정점", d4 < 1e-10, f"최대 변화 {d4:.1e}")

# 5. PIE 가중: 나눗셈 1/P 중 |P|/|P|max 만큼만 적용한다 (alpha 가 무시될 만큼 |P| 가 클 때)
amps = np.array([1.0, 0.5, 0.2, 0.05, 1e-3, 1e-6])
Ptest = amps * np.exp(1j * np.linspace(0, 3, amps.size))
applied = np.abs(pie_weight(Ptest) * Ptest)          # psi' - psi 중 시편에 실제로 반영되는 비율
expect = amps * amps**2 / (amps**2 + 1e-4)
check("PIE 가중이 반영하는 비율 = |P|/|P|max x |P|^2/(|P|^2 + alpha)", np.allclose(applied, expect, rtol=1e-12),
      "비율 " + ", ".join(f"{a:g}->{v:.2g}" for a, v in zip(amps, applied)))

print("\n전부 통과" if ok else "\n실패 있음")

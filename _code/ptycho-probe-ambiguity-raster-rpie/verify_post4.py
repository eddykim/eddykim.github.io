"""4편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 격자 주기 모호성  규칙 격자(간격 a)에서 주기 a 인 복소 함수 f 로 (O f, P / f) 를 만들면 모든 세기가 같다.
                       흔든 격자와 페르마 나선에서는 깨진다. (raster grid pathology 의 원인, Thibault 등 2009)
  2. 척도·위상 기울기  (cO, P/c) 와 (O e^{ig.r}, P e^{-ig.r}) 는 세기를 바꾸지 않는다.
  3. rPIE 의 두 극한   alpha = beta = 1 이면 3편 식으로 따로 짠 ePIE 와 같다.
                       alpha -> 0 이면 한 단계 뒤 P o' = psi' (출사파를 프로브로 나눈 것)가 된다.
  4. 중첩 사영         DM 의 중첩 사영 시편은 sum_j || psi_j - P O_j ||^2 의 정류점이다 (기울기 = 0).
  5. 모멘텀 eta = 0    mPIE 가 모멘텀 없는 rPIE 와 같다.
  6. 세기 제어         곱 P O 를 바꾸지 않고, 투과율 1 인 곳에 놓인 패턴의 총세기는 sum|P|^2 와 같다 (파세발).
  7. 잔차와 흔들림     규칙 격자의 (O f, P/f) 는 잔차가 0 이라 데이터로 가려낼 수 없다.
                       ±1 화소 흔든 격자의 위치 차이는 가로세로 1 화소 이동을 만든다 (정수 조합의 최대공약수 = 1).
"""
import numpy as np

from ptycho4_core import (N, fft2c, disc, pinhole_probe, make_object, scan, object_shape, intensities,
                          fourier_project, rpie, overlap_projection, residual)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
P = pinhole_probe()
yy, xx = np.indices((N, N))

# 1. 격자 주기 모호성
a = 8
f_p = (1 + 0.3 * np.cos(2 * np.pi * xx / a)) * np.exp(0.5j * np.sin(2 * np.pi * yy / a))
for kind in ("regular", "jitter", "fermat"):
    pos = scan(kind)
    shp = object_shape(pos)
    obj = make_object(shp)
    Y, X = np.indices(shp)
    y0, x0 = pos[0]
    f_o = (1 + 0.3 * np.cos(2 * np.pi * (X - x0) / a)) * np.exp(0.5j * np.sin(2 * np.pi * (Y - y0) / a))
    I = intensities(obj, P, pos)
    d = np.abs(intensities(obj * f_o, P / f_p, pos) - I).max() / I.max()
    if kind == "regular":
        check("규칙 격자: (O f, P/f) 가 세기를 바꾸지 않는다", d < 1e-12, f"최대 상대 차이 {d:.1e}")
    else:
        check(f"{kind}: 같은 변환이 세기를 바꾼다", d > 1e-3, f"최대 상대 차이 {d:.1e}")

# 2. 척도와 위상 기울기
pos = scan("fermat"); shp = object_shape(pos); obj = make_object(shp); I = intensities(obj, P, pos)
c = 0.7 * np.exp(1.1j)
d1 = np.abs(intensities(obj * c, P / c, pos) - I).max() / I.max()
Y, X = np.indices(shp)
g = (0.09, -0.13)
d2 = np.abs(intensities(obj * np.exp(1j * (g[0] * X + g[1] * Y)), P * np.exp(-1j * (g[0] * xx + g[1] * yy)), pos) - I).max() / I.max()
check("척도 교환과 위상 기울기는 세기를 바꾸지 않는다", max(d1, d2) < 1e-12, f"최대 상대 차이 {d1:.1e}, {d2:.1e}")

# 3. rPIE 의 두 극한
P0 = disc(N, 14)


def epie_reference(I, pos, probe0, shape, n_iter, rng):
    """3편 식 그대로의 ePIE (rpie 와 독립적으로 작성)."""
    rng = np.random.default_rng(rng)
    sI = np.sqrt(I); o = np.ones(shape, complex); p = probe0.astype(complex).copy()
    for k in range(1, n_iter + 1):
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            v = o[y:y + N, x:x + N]; old = v.copy(); psi = p * old
            dpsi = fourier_project(psi, sI[j]) - psi
            v += np.conj(p) / np.max(np.abs(p)) ** 2 * dpsi
            if k >= 2:
                p += np.conj(old) / np.max(np.abs(old)) ** 2 * dpsi
    return o, p


# 같은 입력에서 한 단계(위치 하나)씩 비교한다. 수백 번 이어 돌리면 연산 순서가 다른 두 식의 반올림 차이가
# 불어나 1e-7 수준까지 벌어진다 (2편의 HIO·RAAR 와 같은 현상).
worst = 0.0
o_s, p_s = make_object(shp, phase_max=1.0), P0 * np.exp(0.3j)
for j in range(20):
    pj = pos[j:j + 1]
    oa, pa = rpie(I[j:j + 1], pj, p_s, shp, 1, alpha=1.0, beta=1.0, probe_start=1, obj0=o_s)
    y, x = pj[0]
    v = o_s.copy(); old = v[y:y + N, x:x + N].copy(); psi = p_s * old
    dpsi = fourier_project(psi, np.sqrt(I[j])) - psi
    v[y:y + N, x:x + N] += np.conj(p_s) / np.max(np.abs(p_s)) ** 2 * dpsi
    pb = p_s + np.conj(old) / np.max(np.abs(old)) ** 2 * dpsi
    worst = max(worst, np.abs(oa - v).max(), np.abs(pa - pb).max())
    o_s, p_s = oa, pa
check("alpha = beta = 1 의 rPIE 한 단계 = 3편 ePIE 한 단계 (20 위치)", worst < 1e-13, f"최대 차이 {worst:.1e}")
y, x = pos[0]
o_start = make_object(shp, phase_max=1.0)
o_new, _ = rpie(I[:1], pos[:1], P, shp, 1, alpha=1e-12, obj0=o_start, probe_start=10**9)
psi_c = fourier_project(P * o_start[y:y + N, x:x + N], np.sqrt(I[0]))
strong = np.abs(P) > 1e-2 * np.abs(P).max()
d3 = np.abs((P * o_new[y:y + N, x:x + N] - psi_c)[strong]).max() / np.abs(psi_c).max()
check("alpha -> 0 의 rPIE 한 단계는 나눗셈 (P o' = psi')", d3 < 1e-8, f"최대 상대 차이 {d3:.1e}")

# 4. 중첩 사영 = 최소자승
psi = np.stack([P * obj[y:y + N, x:x + N] for y, x in pos])
psi = psi + 0.05 * (rng.normal(size=psi.shape) + 1j * rng.normal(size=psi.shape))
o_ls, _ = overlap_projection(psi, np.ones(shp, complex), P, pos, update_probe=False, n_inner=1, eps=0.0 + 1e-30)
grad = np.zeros(shp, complex); cover = np.zeros(shp, bool)
for (y, x), ps in zip(pos, psi):
    grad[y:y + N, x:x + N] += np.conj(P) * (P * o_ls[y:y + N, x:x + N] - ps)
    cover[y:y + N, x:x + N] |= np.abs(P) > 1e-3 * np.abs(P).max()
d4 = np.abs(grad[cover]).max() / np.abs(np.conj(P) * psi[0]).max()
check("중첩 사영의 시편에서 최소자승 기울기 = 0", d4 < 1e-10, f"최대 상대 크기 {d4:.1e}")

# 5. 모멘텀 eta = 0
o3, p3 = rpie(I, pos, P0, shp, 3, alpha=0.05, rng=2)
o4, p4 = rpie(I, pos, P0, shp, 3, alpha=0.05, rng=2, momentum=dict(T=30, eta_obj=0.0, eta_prb=0.0))
check("eta = 0 인 mPIE = rPIE", np.allclose(o3, o4) and np.allclose(p3, p4), f"최대 차이 {np.abs(o3 - o4).max():.1e}")

# 6. 세기 제어
o5, p5 = rpie(I, pos, P0, shp, 2, alpha=0.05, rng=1)
o6, p6 = rpie(I, pos, P0, shp, 2, alpha=0.05, rng=1, probe_power="auto")
y, x = pos[7]
d6 = np.abs(p5 * o5[y:y + N, x:x + N] - p6 * o6[y:y + N, x:x + N]).max() / np.abs(p5 * o5[y:y + N, x:x + N]).max()
# 갱신식이 척도 변환에 등변이라 두 경로는 반올림을 빼면 같다. 반올림은 반복마다 불어난다.
check("세기 제어는 출사파 P O 를 바꾸지 않는다 (2회 뒤)", d6 < 1e-6, f"최대 상대 차이 {d6:.1e}")
flat = np.ones((N + 8, N + 8), complex)
I_flat = intensities(flat, P, np.array([[4, 4]]))
check("투과율 1 인 곳의 패턴 총세기 = sum|P|^2 (파세발)", np.isclose(I_flat.sum(), np.sum(np.abs(P) ** 2)),
      f"{I_flat.sum():.6f} vs {np.sum(np.abs(P) ** 2):.6f}")

# 7. 잔차와 흔들림
from math import gcd
from functools import reduce
pos = scan("regular"); shp = object_shape(pos); obj = make_object(shp); I = intensities(obj, P, pos)
Y, X = np.indices(shp); y0, x0 = pos[0]
f_o = (1 + 0.3 * np.cos(2 * np.pi * (X - x0) / a)) * np.exp(0.5j * np.sin(2 * np.pi * (Y - y0) / a))
r_true, r_lat = residual(obj, P, pos, I), residual(obj * f_o, P / f_p, pos, I)
r_shift = residual(obj, np.roll(P, 1, axis=1), pos, I)
check("규칙 격자의 격자 주기 해는 잔차가 0, 1화소 비켜난 프로브는 아니다",
      r_true < 1e-12 and r_lat < 1e-12 and r_shift > 1e-2, f"{r_true:.1e}, {r_lat:.1e}, {r_shift:.1e}")
pj = scan("jitter"); d = pj - pj[0]
g = [reduce(gcd, np.abs(d[:, k]).tolist()) for k in (0, 1)]
check("±1화소 흔든 격자의 위치 차이: 세로·가로 최대공약수 = 1", g == [1, 1], f"{g}")

print("\n전부 통과" if ok else "\n실패 있음")

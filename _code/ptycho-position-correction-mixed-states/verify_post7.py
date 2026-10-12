"""7편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 소수 화소 이동   shift 를 정수로 주면 np.roll 과 같고, 두 번 옮기면 합만큼 옮긴 것과 같다.
  2. 이동의 미분      dshift 가 shift 의 중앙 유한차분과 같다.
  3. 위치 기울기      grad 보정이 쓰는 dE/d위치 = 2 Re<chi, D> 가 푸리에 오차의 위치 유한차분과 같다.
  4. 모드 사영        모드가 하나면 혼합 상태의 푸리에 사영이 기존 사영과 같고, 사영 뒤 sum_k |Psi_k|^2 = n 이다.
  5. 모드 섞기        모드들을 유니터리 행렬로 섞어도 세기가 그대로다 (모드는 직교화로만 정해진다).
  6. 결맞음 모드      SVD 모드의 세기 합이 광원 점 세기의 가중 합과 같고, 그 모드로 만든 회절 세기가
                     광원 점마다 따로 만든 회절 세기의 합과 같다. 점유율의 합은 1 이다.
"""
import numpy as np

from ptycho7_core import (N, fft2c, ifft2c, pinhole_probe, make_object, scan, object_shape, fourier_project, shift,
                          dshift, split, fourier_err, source_probes, coherent_modes, modal_intensities, modal_project)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
P = pinhole_probe()

# 1. 소수 화소 이동
d1 = np.abs(shift(P, (3, -5)) - np.roll(P, (3, -5), axis=(0, 1))).max() / np.abs(P).max()
a, b = np.array([0.3, -0.7]), np.array([1.45, 0.2])
d2 = np.abs(shift(shift(P, a), b) - shift(P, a + b)).max() / np.abs(P).max()
check("정수 이동 = np.roll, 두 번 이동 = 합의 이동", max(d1, d2) < 1e-12, f"{d1:.1e}, {d2:.1e}")

# 2. 이동의 미분
d0 = np.array([0.37, 0.61]); h = 1e-5
fd = np.stack([(shift(P, d0 + h * e) - shift(P, d0 - h * e)) / (2 * h) for e in np.eye(2)])
d3 = np.abs(fd - dshift(P, d0)).max() / np.abs(fd).max()
check("dshift = 이동의 유한차분", d3 < 1e-8, f"최대 상대 차이 {d3:.1e}")

# 3. 위치 기울기
pos = scan("fermat")[:30].astype(float) + rng.uniform(0, 0.5, (30, 2))
shp = object_shape(np.ceil(pos).astype(int)) ; shp = (shp[0] + 2, shp[1] + 2)
obj = make_object(shp)
from ptycho7_core import intensities_sub
I = intensities_sub(obj, P, pos)
O1 = obj * (1 + 0.05 * rng.normal(size=shp))                           # 참 해에서 조금 벗어난 시편
p = pos[4] + np.array([0.12, -0.08])
(y, x), f = split(p)
o = O1[y:y + N, x:x + N]
Psi = fft2c(shift(P, f) * o)
chi = (1 - np.sqrt(I[4]) / np.abs(Psi)) * Psi
an = np.array([2 * np.real(np.vdot(fft2c(dP * o), chi)) for dP in dshift(P, f)])
fdp = np.array([(fourier_err(O1, P, p + h * e, np.sqrt(I[4])) - fourier_err(O1, P, p - h * e, np.sqrt(I[4]))) / (2 * h)
                for e in np.eye(2)])
d4 = np.abs(an - fdp).max() / np.abs(fdp).max()
check("위치 기울기 2Re<chi, D> = 푸리에 오차의 위치 유한차분", d4 < 1e-5, f"해석 {an.round(3)}, 유한차분 {fdp.round(3)}")

# 4. 모드 사영
psi = P * obj[4:4 + N, 4:4 + N] * 1.3
sI = np.sqrt(np.abs(fft2c(P * obj[4:4 + N, 4:4 + N])) ** 2)
one = ifft2c(modal_project(fft2c(psi)[None], sI))[0]
d5 = np.abs(one - fourier_project(psi, sI)).max() / np.abs(psi).max()
modes = np.stack([P, 0.4 * np.roll(P, 3, 1), 0.2 * np.roll(P, -2, 0)]) * obj[4:4 + N, 4:4 + N]
proj = modal_project(fft2c(modes), sI)
d6 = np.abs(np.sum(np.abs(proj) ** 2, 0) - sI ** 2).max() / (sI ** 2).max()
check("모드 하나 = 기존 사영, 사영 뒤 sum|Psi_k|^2 = n", max(d5, d6) < 1e-12, f"{d5:.1e}, {d6:.1e}")

# 5. 모드 섞기
Ps, ws = source_probes(0.012)
M, occ = coherent_modes(Ps, ws)
M5 = M[:5]
Q, _ = np.linalg.qr(rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5)))
mixed = np.einsum("kl,lyx->kyx", Q, M5)
pos_i = scan("fermat")[:20]
Ia = modal_intensities(obj, M5, pos_i)
Ib = modal_intensities(obj, mixed, pos_i)
d7 = np.abs(Ia - Ib).max() / Ia.max()
check("모드를 유니터리로 섞어도 세기는 그대로", d7 < 1e-12, f"최대 상대 차이 {d7:.1e}")

# 6. 결맞음 모드
lhs = np.sum(np.abs(M) ** 2, 0)
rhs = np.sum(ws[:, None, None] * np.abs(Ps) ** 2, 0)
d8 = np.abs(lhs - rhs).max() / rhs.max()
I_src = np.array([np.sum([w * np.abs(fft2c(Pp * obj[y:y + N, x:x + N])) ** 2 for Pp, w in zip(Ps, ws)], 0)
                  for y, x in pos_i[:6]])
I_mod = modal_intensities(obj, M, pos_i[:6])
d9 = np.abs(I_src - I_mod).max() / I_src.max()
check("SVD 모드 = 광원 점들의 세기 합 (프로브 면, 회절 면)", max(d8, d9) < 1e-12 and abs(occ.sum() - 1) < 1e-12,
      f"{d8:.1e}, {d9:.1e}, 점유율 합 {occ.sum():.12f}")

print("\n전부 통과" if ok else "\n실패 있음")

"""6편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 기울기       세 비용(amp, poisson, int)의 해석적 기울기가 중앙 유한차분의 방향 미분과 같다.
  2. 푸리에 사영  3편의 푸리에 사영 psi' - psi 가 진폭 비용의 기울기 chi 를 뺀 것과 같다 (걸음 1의 경사하강).
  3. 분산 안정화  푸아송 n 에서 sqrt(n) 의 분산은 평균이 커지면 1/4 로 간다. n 자체의 분산은 평균과 같다.
  4. 정류점       노이즈 없는 데이터에서는 참 시편·프로브에서 세 비용의 기울기가 모두 0 이다.
  5. 척도 불변    세 비용 모두 (cO, P/c) 에서 값이 같다.
  6. 단조 감소    ml() 의 비용은 반복마다 줄거나 그대로다.
"""
import numpy as np

from ptycho6_core import (N, fft2c, ifft2c, pinhole_probe, make_object, scan, object_shape, intensities,
                          fourier_project, noisy, total_cost, grad_obj, grad_probe, chi, ml)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
P = pinhole_probe()
pos = scan("fermat")[:40]
shp = object_shape(pos)
obj = make_object(shp)
I = intensities(obj, P, pos)
n, s = noisy(I, 1e5, rng)
Ps = P * np.sqrt(s)
O1 = obj * (1 + 0.1 * rng.normal(size=shp)) * np.exp(0.2j * rng.normal(size=shp))   # 참 해에서 조금 벗어난 점

# 1. 기울기
for model in ("amp", "poisson", "int"):
    g, _ = grad_obj(O1, Ps, pos, n, model)
    gp, _ = grad_probe(O1, Ps, pos, n, model)
    dO = rng.normal(size=shp) + 1j * rng.normal(size=shp)
    dP = rng.normal(size=Ps.shape) + 1j * rng.normal(size=Ps.shape)
    h = 1e-6
    fd_o = (total_cost(O1 + h * dO, Ps, pos, n, model) - total_cost(O1 - h * dO, Ps, pos, n, model)) / (2 * h)
    fd_p = (total_cost(O1, Ps + h * dP, pos, n, model) - total_cost(O1, Ps - h * dP, pos, n, model)) / (2 * h)
    an_o = 2 * np.real(np.vdot(g, dO))             # 비켤레 미분 dL/dO* 와 실수 방향 미분의 관계
    an_p = 2 * np.real(np.vdot(gp, dP))
    e = max(abs(fd_o - an_o) / abs(an_o), abs(fd_p - an_p) / abs(an_p))
    check(f"{model}: 해석적 기울기 = 유한차분 (시편·프로브)", e < 1e-5, f"최대 상대 차이 {e:.1e}")

# 2. 푸리에 사영 = 진폭 비용의 경사 한 걸음
psi = Ps * O1[pos[3][0]:pos[3][0] + N, pos[3][1]:pos[3][1] + N]
step = fourier_project(psi, np.sqrt(n[3])) - psi
grad_step = -ifft2c(chi(fft2c(psi), n[3], "amp"))
d2 = np.abs(step - grad_step).max() / np.abs(step).max()
check("푸리에 사영의 보정 = - (진폭 비용의 기울기)", d2 < 1e-12, f"최대 상대 차이 {d2:.1e}")

# 3. 분산 안정화
lam = np.array([1.0, 5.0, 30.0, 300.0, 3000.0])
x = rng.poisson(lam, size=(400000, lam.size))
v_sqrt = np.var(np.sqrt(x), 0)
v_n = np.var(x, 0)
check("sqrt(n) 의 분산 -> 1/4 (평균 30 이상에서 5 % 이내)", np.all(np.abs(v_sqrt[2:] - 0.25) < 0.0125),
      "평균 " + ", ".join(f"{l:g}: {v:.3f}" for l, v in zip(lam, v_sqrt)))
check("n 의 분산 = 평균 (1 % 이내)", np.all(np.abs(v_n / lam - 1) < 0.01), ", ".join(f"{v / l:.3f}" for l, v in zip(lam, v_n)))

# 4. 노이즈 없는 데이터에서 참 해는 정류점
I0 = I * s
for model in ("amp", "poisson", "int"):
    g, _ = grad_obj(obj, Ps, pos, I0, model)
    g1, _ = grad_obj(O1, Ps, pos, I0, model)
    r = np.abs(g).max() / np.abs(g1).max()
    check(f"{model}: 노이즈 없는 데이터, 참 해의 기울기 = 0", r < 1e-8, f"벗어난 점 대비 {r:.1e}")

# 5. 척도 불변
c = 0.7 * np.exp(0.4j)
for model in ("amp", "poisson", "int"):
    a = total_cost(O1, Ps, pos, n, model)
    b = total_cost(O1 * c, Ps / c, pos, n, model)
    check(f"{model}: 비용(cO, P/c) = 비용(O, P)", abs(a - b) / abs(a) < 1e-12, f"상대 차이 {abs(a - b) / abs(a):.1e}")

# 6. 단조 감소
for model in ("amp", "poisson"):
    cs = []
    ml(n, pos, Ps, shp, 15, model=model, update_probe=True,
       callback=lambda k, o, p: cs.append(total_cost(o, p, pos, n, model)))
    check(f"{model}: ml() 의 비용은 줄기만 한다 (15회)", np.all(np.diff(cs) <= 1e-9 * abs(cs[0])), f"{cs[0]:.6g} -> {cs[-1]:.6g}")

print("\n전부 통과" if ok else "\n실패 있음")

"""2편 계산의 교차검증. 문헌의 주장을, 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 자기상관 정리     |F x|^2 의 역변환  <->  순환 자기상관을 직접 합으로 계산
  2. 자명한 모호성     평행이동·공액 반전·전역 위상이 푸리에 크기를 바꾸지 않는다
  3. 1D 반례           Shechtman 등(2015)의 u, v 가 같은 자기상관을 갖지만 서로 자명하게 같지 않다
  4. 모듈러스 사영     멱등성, 그리고 같은 크기를 가진 무작위 점들보다 가깝다 (가장 가까운 점)
  5. ER = 최급강하     오차 E = || |F g| - m ||^2 의 해석적 경사 2(g - Pm g) 를 유한차분과 대조하고,
                       경사 반 걸음 뒤 지지 사영이 ER 한 단계와 같은지 본다 (Fienup 1982 의 주장)
  6. ER 오차 비증가    Fienup 1982 의 수렴 정리
  7. HIO = RAAR        양수 제약 없이 beta = 1 이면 한 단계가 같다 (Marchesini 2007 의 주장)
"""
import numpy as np

from phase_retrieval import (F, iF, P_m, P_s, make_object, random_start, step, twin, run)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)

# 1. 자기상관 정리 (작은 배열에서 직접 합)
n = 16
x = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
direct = np.zeros((n, n), complex)
for a in range(n):
    for b in range(n):
        direct[a, b] = np.sum(np.roll(x, (-a, -b), axis=(0, 1)) * np.conj(x))   # r[k] = sum_n x[n+k] x*[n]
via_fft = iF(np.abs(F(x)) ** 2) * n            # ortho 정규화의 척도 n 을 되돌린다
check("자기상관 = |F x|^2 의 역변환", np.allclose(direct, via_fft), f"최대 차이 {np.abs(direct - via_fft).max():.1e}")

# 2. 자명한 모호성
obj, sup = make_object(64, 24, "complex")
m = np.abs(F(obj))
for name, y in (("평행이동", np.roll(obj, (5, -3), axis=(0, 1))), ("공액 반전", twin(obj)),
                ("전역 위상", obj * np.exp(1j * 1.234))):
    d = np.abs(np.abs(F(y)) - m).max()
    check(f"{name}는 푸리에 크기를 바꾸지 않는다", d < 1e-12, f"최대 차이 {d:.1e}")

# 3. 1D 반례
u = np.array([1, 0, -2, 0, -2.])
v = np.array([1 - np.sqrt(3), 0, 1, 0, 1 + np.sqrt(3)])
au, av = np.correlate(u, u, "full"), np.correlate(v, v, "full")
trivial = any(np.allclose(s * w, v) for s in (1, -1) for w in (u, u[::-1]))
check("1D 반례: 자기상관이 같고 자명한 관계가 아니다", np.allclose(au, av) and not trivial,
      f"자기상관 {np.round(au, 6).tolist()}")

# 4. 모듈러스 사영
g = rng.normal(size=(64, 64)) + 1j * rng.normal(size=(64, 64))
pm = P_m(g, m)
check("Pm 은 멱등", np.allclose(P_m(pm, m), pm), f"최대 차이 {np.abs(P_m(pm, m) - pm).max():.1e}")
d0 = np.linalg.norm(pm - g)
others = [np.linalg.norm(iF(m * np.exp(1j * (np.angle(F(g)) + 0.3 * rng.normal(size=m.shape)))) - g)
          for _ in range(50)]
check("Pm 은 같은 크기를 가진 다른 점들보다 가깝다", d0 < min(others),
      f"Pm 까지 {d0:.4f}, 다른 점들 중 최소 {min(others):.4f}")

# 5. ER = 최급강하 (복소 물체, 양수 제약 없음)
def E(gg):
    return np.sum((np.abs(F(gg)) - m) ** 2)
grad = 2 * (g - P_m(g, m))                       # 비트링거 미분: dE = Re<grad, dg>
dirs = [rng.normal(size=g.shape) + 1j * rng.normal(size=g.shape) for _ in range(3)]
h = 1e-6
fd = [(E(g + h * d) - E(g - h * d)) / (2 * h) for d in dirs]
an = [np.real(np.vdot(grad, d)) for d in dirs]
rel = max(abs(a - b) / abs(b) for a, b in zip(fd, an))
check("해석적 경사 2(g - Pm g) = 유한차분", rel < 1e-6, f"상대 차이 최대 {rel:.1e}")
er_from_gradient = P_s(g - 0.5 * grad, sup)
d = np.abs(er_from_gradient - step("ER", g, m, sup)).max()
check("경사 반 걸음 + 지지 사영 = ER 한 단계", d < 1e-12, f"최대 차이 {d:.1e}")

# 6. ER 오차 비증가 (실수 양수 물체)
obj_r, sup_r = make_object(64, 24, "real")
m_r = np.abs(F(obj_r))
xk = random_start(m_r, np.random.default_rng(1))
errs = []
for _ in range(200):
    xk = step("ER", xk, m_r, sup_r, real=True)
    errs.append(np.linalg.norm(P_m(xk, m_r) - xk))
inc = max(np.diff(errs))
check("ER 의 오차는 늘지 않는다 (200회)", inc <= 1e-12, f"최대 증가량 {inc:.1e}, {errs[0]:.4f} -> {errs[-1]:.4f}")

# 7. beta = 1 에서 HIO = RAAR
# 같은 입력에서 한 단계씩 비교한다. 궤적 전체를 비교하면 반올림 차이가 반복마다 불어나
# (5회 1e-15 -> 50회 1e-6) 어긋난다. HIO 궤적이 작은 섭동에 민감하다는 뜻이기도 하다.
c = random_start(m, np.random.default_rng(3))
worst = 0.0
for _ in range(50):
    worst = max(worst, np.abs(step("HIO", c, m, sup, beta=1.0) - step("RAAR", c, m, sup, beta=1.0)).max())
    c = step("HIO", c, m, sup, beta=1.0)
check("beta = 1 에서 HIO 한 단계 = RAAR 한 단계 (50개 상태)", worst < 1e-14, f"최대 차이 {worst:.1e}")
a = b = random_start(m, np.random.default_rng(2))
for _ in range(50):
    a, b = step("HIO", a, m, sup, beta=1.0), step("RAAR", b, m, sup, beta=1.0)
print(f"       참고: 같은 시작점에서 50회 돌린 두 궤적의 차이 {np.abs(a - b).max():.1e} (반올림 증폭)")

print("\n전부 통과" if ok else "\n실패 있음")

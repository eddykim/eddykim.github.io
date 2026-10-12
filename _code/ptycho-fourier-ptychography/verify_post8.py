"""8편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 전방 모델     스펙트럼을 잘라 동공을 곱하는 모형이, 실공간에서 시편에 기울어진 평면파를 곱하고 동공으로
                   저역통과한 고해상 영상을 4 화소마다 찍은 것과 같다 (배율 상수 하나 차이).
  2. 겹친 넓이     두 원의 겹친 넓이 공식이 몬테카를로 추정과 같다.
  3. 상호성        푸리에 타이코그래피 rPIE 의 한 번 갱신이, 켤레를 취해 역할을 바꾼 타이코그래피 rPIE(3편 식)의
                   한 번 갱신과 같다.  ifft(x) = conj(fft(conj x)) 이므로 |F^-1[P S]| = |F[P* S*]| 이다.
  4. 동공의 모호성 동공에 위상 기울기를 실어 영상을 1 화소 옮기고 시편을 반대로 4 고해상 화소 옮기면 영상이 같다.
  5. ML 기울기     fp_ml 이 쓰는 스펙트럼 기울기가 유한차분과 같다.
"""
import numpy as np

from ptycho8_core import (NH, NL, fft2c, ifft2c, make_object, fp_pupil, leds, band, fp_spectrum, fp_images,
                          overlap_frac, fp_rpie, rpie, chi, fp_cost)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
obj = make_object((NH, NH), amp_range=(0.4, 1.0), phase_max=2.0)
R = 12
P = fp_pupil(R, (0.6, 0.4))

# 1. 전방 모델
k = np.array([7, -13])
I_model = fp_images(obj, P, [k])[0]
y, x = np.indices((NH, NH))
u = obj * np.exp(-2j * np.pi * (k[0] * y + k[1] * x) / NH)              # S(q + k) 가 가운데로 온다
Pbig = np.zeros((NH, NH), complex)
Pbig[band((0, 0))] = P
I_hr = np.abs(ifft2c(fft2c(u) * Pbig)) ** 2
I_ds = I_hr[::NH // NL, ::NH // NL]
c = np.sum(I_model * I_ds) / np.sum(I_ds ** 2)
d1 = np.abs(I_model - c * I_ds).max() / I_model.max()
check("잘라 낸 스펙트럼 모형 = 실공간 기울어진 조명 + 저역통과 + 성긴 표본", d1 < 1e-10, f"최대 상대 차이 {d1:.1e} (배율 {c:.4f})")

# 2. 겹친 넓이
pts = rng.uniform(-R, R, (2_000_000, 2))
inside = np.sum(pts ** 2, 1) <= R ** 2
dev = 0.0
for dd in (4, 10, 16, 22):
    mc = np.mean(np.sum((pts[inside] - [dd, 0]) ** 2, 1) <= R ** 2)
    dev = max(dev, abs(mc - overlap_frac(R, dd)))
check("겹친 넓이 공식 = 몬테카를로 (네 거리)", dev < 2e-3, f"최대 차이 {dev:.1e}")

# 3. 상호성
ks = np.array([[10, -20]])
I = fp_images(obj, P, ks)
o_fp, _ = fp_rpie(I, ks, P, 1, alpha=0.1)
S_fp = fft2c(o_fp) * (NL / NH)
S0 = np.zeros((NH, NH), complex)
S0[band((0, 0))] = fft2c(np.sqrt(I[0]))                                  # fp_rpie 와 같은 시작점
corner = np.array([[NH // 2 + 10 - NL // 2, NH // 2 - 20 - NL // 2]])
o_pt, _ = rpie(I, corner, np.conj(P), (NH, NH), 1, alpha=0.1, probe_start=10 ** 9, obj0=np.conj(S0))
d3 = np.abs(np.conj(o_pt) - S_fp).max() / np.abs(S_fp).max()
check("푸리에 타이코그래피 한 번 갱신 = 켤레로 옮긴 타이코그래피 rPIE 한 번 갱신", d3 < 1e-12, f"최대 상대 차이 {d3:.1e}")

# 4. 동공의 모호성
ks4 = leds(3, 10)
qy, qx = np.indices((NL, NL))
ramp = np.exp(-2j * np.pi * qx / NL)                                     # 영상을 x 로 1 화소
Ia = fp_images(obj, P, ks4)
Ib = fp_images(np.roll(obj, -NH // NL, axis=1), P * ramp, ks4)
d4 = np.abs(Ia - Ib).max() / Ia.max()
check("동공 위상 기울기 + 시편 반대 이동 = 같은 영상", d4 < 1e-10, f"최대 상대 차이 {d4:.1e}")

# 5. ML 기울기
ks5 = leds(3, 10)
n5 = rng.poisson(fp_images(obj, P, ks5) * 50).astype(float)
S1 = fp_spectrum(obj) * (1 + 0.05 * rng.normal(size=(NH, NH)))
for model in ("amp", "poisson", "int"):
    g = np.zeros_like(S1)
    for kv, nj in zip(ks5, n5):
        g[band(kv)] += np.conj(P * 1.0) * fft2c(chi(ifft2c(P * S1[band(kv)]), nj, model))
    dS = rng.normal(size=S1.shape) + 1j * rng.normal(size=S1.shape)
    h = 1e-6
    fd = (fp_cost(S1 + h * dS, P, ks5, n5, model) - fp_cost(S1 - h * dS, P, ks5, n5, model)) / (2 * h)
    an = 2 * np.real(np.vdot(g, dS))
    check(f"{model}: 스펙트럼 기울기 = 유한차분", abs(fd - an) < 1e-5 * abs(fd), f"상대 차이 {abs(fd - an) / abs(fd):.1e}")

print("\n전부 통과" if ok else "\n실패 있음")

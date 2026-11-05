"""11편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 유리 굴절      na_through_glass 의 광선이 스넬 법칙을 만족하고 옆 이동의 합이 l 이다. 유리가 없으면 l / sqrt(l^2 + H^2).
  2. 나선 순서      UConn 설명의 그림(1 가운데, 2 오른쪽, 3 그 위, 4·5 왼쪽, 6·7 아래 ...)과 같고, 15x15 를 한 번씩 덮는다.
  3. 노이즈 바닥    대역 제한 합성 영상에 알려진 진폭 노이즈를 넣으면 noise_floor 가 그 크기를 되찾는다.
  4. 밝기 배율      LED 마다 밝기가 다른 합성 데이터를 배율을 켜고 6회 복원하면 그 비를 되찾는다.
  5. 잔차           모형과 같은 데이터의 잔차는 0 이고, 영상마다 척도를 곱해도 0 이다.
  6. 켤레 모호성    동공이 실수면 모든 조명 파수의 부호를 뒤집고 시편을 켤레로 바꿔도 영상이 같다.
  7. FRC            같은 영상끼리는 1 이고, 독립 노이즈만 있는 영상끼리는 0 근처다.
"""
import numpy as np

import fpm_core as C

ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


# 1. 유리 굴절
H, t, n = 90.88, 1.0, 1.52
worst = 0.0
for l in (4.0, 20.0, 40.0, 56.0):
    na = C.na_through_glass(l, H, t, n)
    th_a = np.arcsin(na)
    th_g = np.arcsin(na / n)
    a = H * np.tan(th_a)
    worst = max(worst, abs(a + t * np.tan(th_g) - l))
na0 = C.na_through_glass(40.0, H, 0.0, n)
check("유리 굴절: 옆 이동의 합 = l, 유리 없으면 기하 각", worst < 1e-9 and abs(na0 - 40 / np.hypot(40, H)) < 1e-12,
      f"최대 차이 {worst:.1e} mm, 유리 1 mm 가 NA 를 {C.na_through_glass(56, H, t, n) - 56 / np.hypot(56, H):+.4f} 바꾼다(56 mm)")

# 2. 나선 순서
sp = C.spiral(15).astype(int)
first = [tuple(p) for p in sp[:10]]
want = [(0, 0), (1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (2, -1)]
check("나선: 설명 그림과 같고 15x15 를 한 번씩", first == want and len({tuple(p) for p in sp}) == 225 and np.abs(sp).max() == 7)

# 공통 합성 문제
rng = np.random.default_rng(0)
nl, nh, r = 64, 192, 14
y, x = np.indices((nh, nh)) - nh // 2
obj = np.exp(1j * 0.8 * np.sin(x / 5.0) * np.cos(y / 7.0)) * (0.7 + 0.3 * np.cos(x / 9.0 + y / 4.0))
S_true = C.fft2c(obj) * (nl / nh)
P = C.defocus_pupil(nl, r, 0, 1, 1)
ks = np.array([(i, j) for i in range(-20, 21, 10) for j in range(-20, 21, 10)])
amps = C.model_amps(S_true, P, ks)

# 3. 노이즈 바닥
sig = 0.02
noisy_amp = amps + sig * amps[:, :, :].mean((1, 2))[:, None, None] * rng.normal(size=amps.shape)
I_n = noisy_amp ** 2
nf = C.noise_floor(I_n, r)
bf = np.sum(ks ** 2, 1) <= r ** 2
rel_true = sig * amps.mean((1, 2)) / np.sqrt((amps ** 2).mean((1, 2)))
check("노이즈 바닥: 알려진 진폭 노이즈를 되찾는다 (명시야)", np.all(np.abs(nf[bf] / rel_true[bf] - 1) < 0.2),
      f"되찾은/참 {np.round(nf[bf] / rel_true[bf], 2)}")

# 4. 밝기 배율
gain = np.exp(0.2 * rng.normal(size=len(ks)))
I_g = (amps * np.sqrt(gain)[:, None, None]) ** 2
S, P2, c = C.fp_cal(I_g, ks, P, nh, 6, bright_from=1)
est = 1 / c ** 2                       # 측정 진폭 x c = 모형 진폭  ->  c^2 = 1 / 밝기
ratio = (est / gain) / np.median(est / gain)
check("밝기 배율: LED 밝기 비를 되찾는다 (중앙값 기준)", np.abs(ratio - 1).max() < 0.05, f"최대 상대 차이 {np.abs(ratio - 1).max():.3f}")

# 5. 잔차
I0 = amps ** 2
r0 = C.residual(S_true, P, ks, I0)
r1 = C.residual(S_true, P, ks, I0 * 3.7)
check("잔차: 모형과 같은 데이터는 0, 척도에 무관", max(r0, r1) < 1e-12, f"{r0:.1e}, {r1:.1e}")

# 6. 켤레 모호성 (실수 동공)
S_conj = C.fft2c(np.conj(obj)) * (nl / nh)
a2 = C.model_amps(S_conj, P, -ks)
d6 = np.abs(a2 - amps).max() / amps.max()
Pc = C.defocus_pupil(nl, r, 30e-6, 0.5e-6, 100e-6)            # 복소 동공이면 깨진다
d6b = np.abs(C.model_amps(S_conj, Pc, -ks) - C.model_amps(S_true, Pc, ks)).max() / amps.max()
check("켤레 모호성: 실수 동공이면 k -> -k, 시편 -> 켤레가 같은 영상, 초점 어긋남이 있으면 깨진다", d6 < 1e-12 and d6b > 1e-3,
      f"{d6:.1e}, 초점 어긋남 동공 {d6b:.2f}")

# 7. FRC
img = np.abs(obj[:128, :128])
f, v = C.frc(img, img)
f2, v2 = C.frc(rng.normal(size=(128, 128)), rng.normal(size=(128, 128)))
check("FRC: 같은 영상 1, 독립 노이즈 0 근처", np.allclose(v[1:], 1) and np.abs(v2[5:]).mean() < 0.1, f"노이즈 평균 {np.abs(v2[5:]).mean():.3f}")

print("\n전부 통과" if ok else "\n실패 있음")

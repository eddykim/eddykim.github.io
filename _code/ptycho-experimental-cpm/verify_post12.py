"""12편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다. 데이터 없이 돈다(1분 안팎).

  1. 화소 크기      주기 p 인 격자의 회절 봉우리가 검출기에서 λz/p 에 놓인다(격자 방정식) — dxo = λz/(N dxd).
  2. 이차 위상      프레넬 회절 적분을 직접 더한 세기 = 출사파에 exp(iπr²/λz) 를 곱해 FFT 한 세기.
                    그래서 원시야가 아니어도(프레넬 수 > 1) 그 위상을 프로브가 떠안으면 원시야 모형이 맞다.
  3. 척도 모호성    (sP, O/s) 는 같은 데이터를 낸다 — 데이터 오차가 똑같다.
  4. 척도 비대칭    rPIE 시편 갱신은 P → sP 에 대해 ΔO(sP) = β P*(ψ'/s − ψ)/D(P). 어두운 출발(s ≪ 1)은
                    측정 진폭 ψ' 를 1/s 배 키워 시편에 넣는다. 같은 데이터·같은 해라도 갱신이 지나는 길은 s 에 달렸다.
  5. 무게중심 맞춤  프로브와 시편을 함께 굴리면 데이터 오차가 그대로이고, 프로브만 굴리면 커진다.
  6. 직접 구현      합성 데이터(무작위 위치 169 개)에서 곡률을 대강 넣은 원 프로브로 어둡게 시작한 mPIE 가
                    참 시편을 되찾는다(상관 0.95 이상).
  7. 상관 판정      이동·180° 회전한 같은 영상은 1, 독립 노이즈끼리는 0 근처.
"""
import numpy as np

import cpm_core as C

ok = True


def check(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f": {detail}" if detail else ""))


lam, z, dxd, N = 708.8e-9, 38.6e-3, 18.16e-6, 364
dxo = lam * z / (N * dxd)

# 1. 화소 크기 — 격자 방정식
p = N * dxo / 16                                 # 창 안에 16 주기
x = (np.arange(N) - N // 2) * dxo
grating = np.ones((N, N)) * (1 + np.cos(2 * np.pi * x / p))[None, :]
I = np.abs(C.fft2c(grating)) ** 2
row = I[N // 2]
side = np.argsort(row)[-3:]                      # 0차와 ±1차
k = np.max(np.abs(side - N // 2))
check("화소 크기: ±1 차 봉우리가 검출기의 λz/p 에", abs(k * dxd - lam * z / p) < 1e-12,
      f"dxo = {dxo * 1e6:.3f} µm, 봉우리 {k} 화소 = {k * dxd * 1e3:.4f} mm, λz/p = {lam * z / p * 1e3:.4f} mm")

# 2. 이차 위상 — 프레넬 적분을 직접 더해 비교 (작은 창)
n = 48
dx_s = lam * z / (n * dxd)
rng = np.random.default_rng(0)
xs = (np.arange(n) - n // 2) * dx_s
Xs, Ys = np.meshgrid(xs, xs)
psi = (np.hypot(Xs, Ys) < 15 * dx_s) * np.exp(1j * rng.uniform(0, 2 * np.pi, (n, n)))
I_fft = C.fresnel_intensity(psi, dx_s, lam, z)
xd = (np.arange(n) - n // 2) * dxd
worst = 0.0
for (u, v) in [(0, 0), (5, -3), (-17, 11), (20, 20), (-24, 7)]:
    X0, Y0 = xd[v + n // 2], xd[u + n // 2]
    U = np.sum(psi * np.exp(1j * np.pi * ((X0 - Xs) ** 2 + (Y0 - Ys) ** 2) / (lam * z)))
    worst = max(worst, abs(abs(U) ** 2 / n ** 2 - I_fft[u + n // 2, v + n // 2]) / I_fft.max())
NF = (0.3e-3) ** 2 / (lam * z)
check("이차 위상: 프레넬 적분의 직접 합 = 이차 위상 곱 + FFT", worst < 1e-9,
      f"최대 상대 차 {worst:.1e}; 이 데이터의 프레넬 수 a²/(λz) = {NF:.2f} (원시야가 아님)")

# 합성 실공간 타이코 문제 (6·5·4·3 이 함께 쓴다)
Np = 64
pos = rng.integers(2, 103, (169, 2))           # 무작위 위치(격자 주사의 주기 모호성을 피한다, 4편)
No = int(pos.max() + Np + 4)
yy, xx = np.indices((No, No))
O_true = (0.75 + 0.25 * np.sin(xx / 4.0) * np.cos(yy / 6.0)) * np.exp(1j * 0.7 * np.cos(xx / 9.0 + yy / 13.0))
yp, xp = np.indices((Np, Np)) - Np // 2
P_true = (np.hypot(yp, xp) < 18) * np.exp(1j * 0.012 * (xp ** 2 + yp ** 2))      # 초점이 어긋난 원 프로브
I_syn = np.array([np.abs(C.fft2c(O_true[r:r + Np, c:c + Np] * P_true)) ** 2 for r, c in pos])

# 3. 척도 모호성
e1 = C.data_error(O_true, P_true, pos, I_syn)
e2 = C.data_error(O_true / 7.3, P_true * 7.3, pos, I_syn)
check("척도 모호성: (sP, O/s) 의 데이터 오차가 같다", abs(e1 - e2) < 1e-9 and e1 < 1e-9, f"{e1:.1e}, {e2:.1e}")

# 4. 척도 비대칭 — 실제 갱신 함수로
o = dict(C.DEFAULT)
patch = np.ones((Np, Np), complex)
P0 = (np.hypot(yp, xp) < 18).astype(complex)
Ij = I_syn[len(pos) // 2]
psi_new = None
worst = 0.0
sizes = {}
for s in (0.03, 1.0, 3.0):
    new_patch, _, _ = C.rpie_step(patch, s * P0, Ij, o)
    psi = patch * P0
    Psi = C.fft2c(psi)
    psi_p = C.ifft2c(Psi * np.sqrt(Ij / (np.abs(s * Psi) ** 2 + 1e-10)) * s)   # 측정 진폭, 위상은 ψ 의 것
    aP = np.abs(P0) ** 2
    want = patch + o["beta_o"] * np.conj(P0) * (psi_p / s - psi) / (o["alpha_o"] * aP.max() + (1 - o["alpha_o"]) * aP)
    worst = max(worst, np.max(np.abs(new_patch - want)))
    sizes[s] = np.sqrt(np.mean(np.abs(new_patch - patch)[np.abs(P0) > 0] ** 2))
check("척도 비대칭: ΔO(sP) = β P*(ψ'/s − ψ)/D(P)", worst < 1e-9,
      "시편 갱신의 rms 크기 " + ", ".join(f"s={s}: {v:.3g}" for s, v in sizes.items()))

# 5. 무게중심 맞춤
Pshift = np.roll(P_true, (5, -4), (0, 1))
Oshift = np.roll(O_true, (5, -4), (0, 1))
z0 = np.zeros_like(P_true)
Pa, _, _, Oa, _, _ = C.com_center(Pshift, z0, z0, Oshift, np.zeros_like(O_true), np.zeros_like(O_true))
ea = C.data_error(Oa, Pa, pos, I_syn)
eb = C.data_error(Oshift, C.com_center(Pshift)[0], pos, I_syn)
check("무게중심: 함께 굴리면 오차 그대로, 프로브만 굴리면 커진다", ea < 1e-9 and eb > 1.0,
      f"함께 {ea:.1e}, 프로브만 {eb:.1f} (패턴 {len(pos)} 장)")

# 6. 직접 구현이 합성 시편을 되찾는가
d = {"I": I_syn, "N": Np, "wavelength": lam, "zo": z, "dxd": dxd, "dxo": dxo,
     "entrancePupilDiameter": 36 * dxo}
P_start = (np.hypot(yp, xp) < 18) * np.exp(1j * 0.008 * (xp ** 2 + yp ** 2))   # 곡률은 참값의 2/3 만 안다고 둔다
O_rec, P_rec, trace, _ = C.reconstruct(d, seed=0, P0=P_start, pos=pos, **dict(C.SETTINGS["dim"], rounds=1, iters=200))
score = C.ncc(np.abs(C.scan_crop(O_true, pos, Np)), np.abs(C.scan_crop(O_rec, pos, Np)))
check("직접 구현: 어둡게 시작해 합성 시편을 되찾는다", score > 0.95 and trace[-1, 1] < 0.1 * trace[0, 1],
      f"상관 {score:.3f}, 오차 {trace[0, 1]:.1f} → {trace[-1, 1]:.2f}")

# 7. 상관 판정
a = np.abs(C.scan_crop(O_true, pos, Np))
big = np.abs(O_true)
b = np.roll(big, (3, -2), (0, 1))[::-1, ::-1]
s1 = C.ncc(a, np.abs(C.scan_crop(np.roll(O_true, (3, -2), (0, 1)), pos, Np)))
s2 = C.ncc(a, a[::-1, ::-1])
s3 = C.ncc(rng.standard_normal((80, 80)), rng.standard_normal((80, 80)))
check("상관 판정: 이동 1, 180° 회전 1, 독립 노이즈 0 근처", s1 > 0.9 and s2 > 0.999 and s3 < 0.15,
      f"{s1:.3f}, {s2:.3f}, {s3:.3f}")

print("모두 통과" if ok else "실패 있음")

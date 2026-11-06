"""acse.py 를 독립 경로와 대조한다.

acse.py 는 신호를 부품 뮬러 행렬의 곱으로만 만든다. 여기서는 그 결과를
  - a = cos 2Psi, b = sin 2Psi sin Delta 의 닫힌 식
  - 저널 식 5 (지연자 방위각 Omega_C 가 일반값일 때)
  - 5편의 식 (Omega = 0 에서 부호가 뒤집힌 1 - a cos phi - b sin phi)
  - 손으로 전개한 결과 (수직 입사 줄에 빔 스플리터 이색성이 남기는 무늬)
  - 5편 코드(channeled.py)의 두께 창 함수
와 맞춰 본다. 실행: python verify_acse.py
"""
import os
import sys

import numpy as np
from scipy.signal import hilbert

import acse as A

D = np.deg2rad
FAIL = []


def check(name, ok, detail=""):
    print(f"[{'OK ' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAIL.append(name)


def fourier_ab(signal_fn):
    """phi 를 0, pi/2, pi, 3pi/2 로 넣은 네 세기에서 정규화 계수."""
    I = [signal_fn(p) for p in (0.0, np.pi / 2, np.pi, 3 * np.pi / 2)]
    dc = 0.5 * (I[0] + I[2])
    return (I[0] - I[2]) / (2 * dc), (I[1] - I[3]) / (2 * dc)


def single_point(psi, delta, omega, pol_c, ana=D(-45.0)):
    """한 점의 세기를 뮬러 곱으로 (verify 쪽에서 따로 조립한다)."""
    def f(phi):
        M = A._chain(A.mueller_polarizer(ana), A.mueller_rotation(-omega),
                     A.mueller_sample(psi, delta), A.mueller_rotation(omega),
                     A.mueller_retarder(pol_c, phi), A.mueller_polarizer(0.0))
        return (M @ A.S_UNPOLARIZED)[0]
    return f


# 1. 닫힌 식 a = cos2Psi, b = sin2Psi sinDelta  (Omega = 90, 270)
rng = np.random.default_rng(7)
psi = rng.uniform(D(5), D(85), 40)
dl = rng.uniform(-np.pi, np.pi, 40)
worst = 0.0
for p, d in zip(psi, dl):
    for om in (D(90.0), D(270.0)):
        a, b = fourier_ab(single_point(p, d, om, D(45.0)))
        worst = max(worst, abs(a - np.cos(2 * p)), abs(b - np.sin(2 * p) * np.sin(d)))
check("Omega=90/270 에서 a=cos2Psi, b=sin2Psi sinDelta", worst < 1e-12, f"최대 차 {worst:.1e}")

# 2. 5편의 부호: Omega = 0 이면 두 항의 부호가 함께 뒤집힌다
worst = 0.0
for p, d in zip(psi, dl):
    a, b = fourier_ab(single_point(p, d, 0.0, D(45.0)))
    worst = max(worst, abs(a + np.cos(2 * p)), abs(b + np.sin(2 * p) * np.sin(d)))
check("Omega=0 에서 1 - a cos phi - b sin phi (5편 1절)", worst < 1e-12, f"최대 차 {worst:.1e}")

# 3. 저널 식 5 — 지연자 방위각 C 가 일반값일 때의 정규화 계수
def journal_eq5(p, d, C):
    s2c, c2c = np.sin(2 * C), np.cos(2 * C)
    den = 1 - s2c * np.sin(2 * p) * c2c * np.cos(d) + c2c ** 2 * np.cos(2 * p)
    al = (s2c * np.cos(2 * p) + np.sin(2 * p) * c2c * np.cos(d)) * s2c / den
    be = s2c * np.sin(d) * np.sin(2 * p) / den
    return al, be


res = {}
for ana_deg in (-45.0, 45.0):
    worst = 0.0
    for p, d in zip(psi[:15], dl[:15]):
        for C in D(np.array([10.0, 25.0, 45.0, 60.0, 80.0])):
            a, b = fourier_ab(single_point(p, d, D(90.0), C, ana=D(ana_deg)))
            al, be = journal_eq5(p, d, C)
            worst = max(worst, abs(a - al), abs(b - be))
    res[ana_deg] = worst
print(f"       저널 식 5 와의 최대 차: 분석기 -45° {res[-45.0]:.1e}, +45° {res[45.0]:.1e}")
check("저널 식 5 는 Omega=90 에서 Omega_C 일반값까지 한 분석기 각으로 맞는다",
      min(res.values()) < 1e-12,
      f"맞는 분석기 각 {'-45' if res[-45.0] < res[45.0] else '+45'}°")
check("Omega_C = 45 도에서 기저 신호에 Psi, Delta 가 남지 않는다(식 5 분모 = 1)",
      all(abs(1 - np.sin(2 * D(45)) * np.sin(2 * p) * np.cos(2 * D(45)) * np.cos(d)
              + np.cos(2 * D(45)) ** 2 * np.cos(2 * p) - 1) < 1e-12 for p, d in zip(psi, dl)))

# 4. acse.ideal_ab 가 닫힌 식과 같다 (격자 전체)
sig = A.sigma_grid(300)
th = A.theta_grid(61)
for dn in (0.0, 100.0, 1000.0):
    a, b = A.ideal_ab(sig, th, dn)
    SIG, TH = np.meshgrid(sig, th)
    ps, de, _ = A.film_psi_delta(SIG, TH, dn)
    err = max(np.max(np.abs(a - np.cos(2 * ps))), np.max(np.abs(b - np.sin(2 * ps) * np.sin(de))))
    check(f"ideal_ab 격자 전체 (막 {dn:.0f} nm)", err < 1e-12, f"{err:.1e}")

# 5. 수직 입사 줄에는 무늬가 없다, +-theta 는 같은 신호다
i0 = int(np.argmin(np.abs(th)))
for dn in (0.0, 100.0, 1000.0):
    a, b = A.ideal_ab(sig, th, dn)
    amp0 = np.max(np.hypot(a[i0], b[i0]))
    check(f"theta = 0 에서 A' = 0 (막 {dn:.0f} nm)", amp0 < 1e-12, f"{amp0:.1e}")
    img = A.detect(sig, th, dn, A.IDEAL)
    sym = np.max(np.abs(img - img[::-1])) / np.max(img)
    check(f"+-theta 대칭 (막 {dn:.0f} nm)", sym < 1e-12, f"{sym:.1e}")

# 6. 수직 입사 줄의 무늬는 장비 몫이다: 빔 스플리터 이색성 D 를 왕복하면 A' = 2D/(1+D^2)
for Dd in (0.02, 0.04, 0.10):
    e = dict(A.IDEAL, bs_diatt=Dd)
    f = lambda phi: A.detect(sig[:1], np.array([0.0]), 300.0,
                             dict(e, phi0=phi - A.retardance(sig[:1])[0]), source=False)[0, 0]
    a, b = fourier_ab(f)
    check(f"theta=0 무늬 진폭 = 2D/(1+D^2), D={Dd}", abs(np.hypot(a, b) - 2 * Dd / (1 + Dd ** 2)) < 1e-12,
          f"{np.hypot(a, b):.6f} 대 {2 * Dd / (1 + Dd ** 2):.6f}")

# 7. 힐베르트 부호: Re(Z e^{-j phi}) = a, -Im = b
n = 1024
x = np.arange(n)
phi = 2 * np.pi * 40 * x / n          # 정수 개 주기라 이산 HT 가 정확하다
for a0, b0 in ((0.5, 0.8138), (-0.3, -0.6), (0.9, -0.2)):
    z = hilbert(a0 * np.cos(phi) + b0 * np.sin(phi)) * np.exp(-1j * phi)
    check(f"HT 부호 a={a0}, b={b0}", np.allclose(z.real, a0, atol=1e-10)
          and np.allclose(-z.imag, b0, atol=1e-10),
          f"Re {z.real.mean():+.4f}, Im {z.imag.mean():+.4f}")

# 8. 위상·가시도 오차만 있으면 오차 지도가 그 오차를 그대로 담는다
#    (a) 모형 자체: 지연량을 pi/2 씩 네 번 밀어 계수를 정확히 뽑는다. 해석 경로의 수치 오차가 없다
e = dict(A.IDEAL, phi0=D(25.0), kappa=4.0, depol=0.97)
sig2 = A.sigma_grid(600)
th2 = np.deg2rad(np.array([20.0, 35.0, 50.0]))
sbar = 0.5 * (1 / A.LAM_MIN + 1 / A.LAM_MAX)
dphi = e["phi0"] + e["kappa"] * (sig2 - sbar) ** 2


def stepped_ab(d_nm, err):
    """a - j b 를 공칭 지연량 기준으로 (네 계단 위상 이동)."""
    I = [A.detect(sig2, th2, d_nm, dict(err, phi0=err["phi0"] + k * np.pi / 2))
         for k in range(4)]
    dc = 0.25 * sum(I)
    c = sum(Ik * np.exp(-1j * k * np.pi / 2) for k, Ik in enumerate(I)) / (2 * dc)
    z = c * np.exp(-1j * A.retardance(sig2))     # A' exp(-j chi_M)
    return z.real, -z.imag


am, bm = stepped_ab(0.0, e)
ar, br = A.ideal_ab(sig2, th2, 0.0)
ea, ec = A.error_maps(am, bm, ar, br)
dev_c = np.max(np.abs(np.angle(np.exp(1j * (ec + dphi)))))
dev_a = np.max(np.abs(ea - 0.97))
check("[모형] eps_chi = -(지연량 오차), 5편 delta-phi 의 지도판", dev_c < 1e-9,
      f"최대 {np.rad2deg(dev_c):.1e}°")
check("[모형] eps_A = 가시도 0.97", dev_a < 1e-9, f"최대 차 {dev_a:.1e}")
a, b = stepped_ab(300.0, e)
a0, b0 = A.ideal_ab(sig2, th2, 300.0)
ac, bc = A.apply_correction(a, b, ea, ec)
r = np.max(np.hypot(ac - a0, bc - b0))
check("[모형] 그 지도로 다른 시편(300 nm)이 정확히 돌아온다", r < 1e-9, f"최대 잔차 {r:.1e}")

#    (a') 시편과 섞이는 오차(대물렌즈 지연, 실험실 축 빔 스플리터)는 정확한 계수로도 남는다
for name, eo in (("대물렌즈 지연 6°", dict(A.IDEAL, obj_ret=D(6.0))),
                 ("빔 스플리터 이색성 0.04", dict(A.IDEAL, bs_diatt=0.04))):
    am_, bm_ = stepped_ab(0.0, eo)
    ea_, ec_ = A.error_maps(am_, bm_, ar, br)
    a_, b_ = stepped_ab(300.0, eo)
    ac_, bc_ = A.apply_correction(a_, b_, ea_, ec_)
    r_ = np.max(np.hypot(ac_ - a0, bc_ - b0))
    check(f"[모형] {name}: 기준 시편과 다른 시편에 잔차가 남는다", r_ > 1e-3, f"최대 잔차 {r_:.1e}")

#    (b) 해석 경로(EEM + HT): 이산 HT 의 경계 오차가 끝에서 1/거리 꼴로 안쪽까지 들어오고
#        (끝 10 % 지점에서 1~3 %, 가운데에서 0.02 % — 본문 5절), 무늬가 약한 줄에서는 EEM 의
#        기저 신호 오차가 더해진다. 대역 가운데 40 % 에서 그 크기 안에 드는지만 본다
k = slice(int(0.3 * len(sig2)), int(0.7 * len(sig2)))
am, bm = A.analyze_image(sig2, A.detect(sig2, th2, 0.0, e))
ea, ec = A.error_maps(am, bm, ar, br)
dev_c = np.max(np.abs(np.angle(np.exp(1j * (ec[:, k] + dphi[k])))))
dev_a = np.max(np.abs(ea[:, k] - 0.97))
check("[해석 경로] eps_chi 가 -(지연량 오차)와 1.5° 안에서 맞는다", dev_c < D(1.5),
      f"최대 {np.rad2deg(dev_c):.2f}°")
check("[해석 경로] eps_A 가 0.97 과 3 % 안에서 맞는다", dev_a < 3e-2, f"최대 차 {dev_a:.1e}")
a, b = A.analyze_image(sig2, A.detect(sig2, th2, 300.0, e))
ac, bc = A.apply_correction(a, b, ea, ec)
r = np.max(np.hypot(ac - a0, bc - b0)[:, k])
check("[해석 경로] 300 nm 시편 잔차가 1e-2 아래", r < 1e-2, f"대역 가운데 최대 잔차 {r:.1e}")

# 9. 두께 창 — 5편 코드와 대조
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "ellipsometer-channeled-spectroscopic"))
import channeled as C5  # noqa: E402

t5 = C5.thickness_window(A.LAM_MIN, A.LAM_MAX, A.DLAM_LSS, ratio=(1,))
t7 = A.thickness_window()
check("두께 창이 5편 thickness_window(지연자 하나)와 같다",
      np.allclose(t5, t7, rtol=1e-12), f"{t7[0]:.1f} ~ {t7[1]:.1f} um")
check("저널 지연자 4092.3 um 가 창 안에 있다", t7[0] < A.T_RETARDER < t7[1])

# 10. 사각 창 가시도 = sinc — 수치 적분과 대조
for lam in (0.475, 0.575, 0.675):
    offs = (np.arange(4001) + 0.5) / 4001 - 0.5
    s = 1.0 / (lam + offs * A.DLAM_LSS)
    num = abs(np.mean(np.exp(2j * np.pi * A.opd() * s)))
    check(f"가시도 sinc 근사 {lam*1e3:.0f} nm", abs(num - A.visibility_box(lam, A.T_RETARDER)) < 2e-3,
          f"수치 {num:.4f}, sinc {A.visibility_box(lam, A.T_RETARDER):.4f}")

# 11. 무늬 수 N 과 극대 개수
y = A.detect(sig2, np.array([D(50.0)]), 0.0, A.IDEAL, source=False)[0]
(xu, _), _ = A.find_extrema(sig2, y)
check("대역 안 극대 개수가 N 과 1 개 이내로 맞는다", abs(len(xu) - A.fringe_count()) <= 1.0,
      f"극대 {len(xu)}개, N = {A.fringe_count():.1f}")

print()
print("모두 통과" if not FAIL else f"실패 {len(FAIL)}건: {FAIL}")
sys.exit(1 if FAIL else 0)

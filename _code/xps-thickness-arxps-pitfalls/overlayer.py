"""XPS 기초 5편 — 오버레이어 두께의 정·역계산, 감도 분석, ARXPS 순방향 계산.

실행: python overlayer.py
      2절 두께 역산 예제, 그림 2의 감도 수치, 오염층 효과, 그림 3의 ARXPS 구분 가능성을 출력한다.
      generate_figures.py가 이 모듈을 불러 그림을 그린다.

기준 계는 SiO2/Si의 Si 2p다.
- 감쇠길이 L = 3.0 nm: 2편에서 TPP-2M으로 구한 SiO2 속 Si 2p의 IMFP(3.88 nm)에
  Powell(2020)이 정리한 Si의 EAL/IMFP 비 하한 0.77을 곱한 수준. 산화물과 기판에 같은 값을 쓴다.
- R0 = I(두꺼운 SiO2) / I(Si) = 0.88 (실측), 0.53 (계산): Seah & Spencer, SIA 33, 640 (2002).
"""
import numpy as np
from scipy.special import erfc

L = 3.0            # nm, 감쇠길이
R0 = 0.88          # 실측 (Seah & Spencer 2002)
R0_CALC = 0.53     # 계산값 (같은 논문)


# ---------------------------------------------------------------------------
# 균일 오버레이어 모델 (θ는 법선 기준)
# ---------------------------------------------------------------------------
def ratio(d, theta=0.0, lam=L, r0=R0):
    """두께 d(nm) 오버레이어의 강도비 R = I_o / I_m."""
    c = np.cos(np.radians(theta))
    return r0 * (np.exp(d / (lam * c)) - 1)


def thickness(R, theta=0.0, lam=L, r0=R0):
    """강도비 R에서 두께를 역산: d = λ cosθ ln(1 + R/R0)."""
    c = np.cos(np.radians(theta))
    return lam * c * np.log(1 + R / r0)


def signals(d, n_sub0, theta=0.0, lam=L, r0=R0):
    """기판만 있을 때 기판 피크 계수가 n_sub0이 되는 측정 조건에서, 두께 d의 두 피크 계수."""
    c = np.cos(np.radians(theta))
    i_m = n_sub0 * np.exp(-d / (lam * c))
    i_o = n_sub0 * r0 * (1 - np.exp(-d / (lam * c)))
    return i_o, i_m


def counting_sigma(d, n_sub0, theta=0.0):
    """포아송 계수 잡음에서 오는 두께의 표준편차 (오차 전파)."""
    i_o, i_m = signals(d, n_sub0, theta)
    R = i_o / i_m
    rel = np.sqrt(1 / i_o + 1 / i_m)
    c = np.cos(np.radians(theta))
    return L * c * (R * rel) / (R0 + R)


def ratio_error(d, frac=0.05, theta=0.0):
    """강도비에 frac만큼 오차가 있을 때의 두께 오차."""
    return thickness(ratio(d, theta) * (1 + frac), theta) - d


# ---------------------------------------------------------------------------
# 오염층 효과: 탄화수소층이 두 피크를 함께 감쇠
# ---------------------------------------------------------------------------
def tpp2m(E, Nv, rho, M, Eg=0.0):
    """TPP-2M IMFP (nm). 2편 코드와 같은 식."""
    Ep = 28.816 * np.sqrt(Nv * rho / M)
    U = (Ep / 28.816) ** 2
    beta = -1.0 + 9.44 / np.sqrt(Ep**2 + Eg**2) + 0.69 * rho**0.1
    gamma = 0.191 * rho**-0.5
    C, D = 19.7 - 9.1 * U, 534 - 208 * U
    return E / (Ep**2 * (beta * np.log(gamma * E) - C / E + D / E**2))


def contamination_effect(t_c=1.0, d=2.0):
    """두께 t_c의 탄화수소층(폴리에틸렌 근사)이 있을 때 역산 두께가 얼마나 변하는가.
    Si 2p의 산화물·기판 성분은 운동에너지가 약 4 eV 다르다(Al Kα, 103.3 / 99.4 eV)."""
    hv = 1486.6
    lam_ox = tpp2m(hv - 103.3, 6, 0.92, 14.027, 8.0)
    lam_m = tpp2m(hv - 99.4, 6, 0.92, 14.027, 8.0)
    i_o, i_m = signals(d, 1.0)
    R_c = (i_o * np.exp(-t_c / lam_ox)) / (i_m * np.exp(-t_c / lam_m))
    return lam_ox, lam_m, thickness(R_c) - d


# ---------------------------------------------------------------------------
# ARXPS 순방향 계산
# ---------------------------------------------------------------------------
Z = np.linspace(0, 40, 8001)
DZ = Z[1] - Z[0]


def profile_step(d):
    """산화물 비율 c(z): 깊이 d까지 1, 그 아래 0 (급격한 계면)."""
    return (Z < d).astype(float)


def profile_graded(center, width):
    """오차함수 모양의 확산 계면. width는 Gaussian 표준편차(nm)."""
    return 0.5 * erfc((Z - center) / (np.sqrt(2) * width))


def arxps_ratio(c, angles):
    """각도마다 산화물/기판 강도비. I(θ) ∝ ∫ c(z) exp(-z/(L cosθ)) dz (R0 상수는 생략)."""
    out = []
    for t in np.atleast_1d(angles):
        k = np.exp(-Z / (L * np.cos(np.radians(t))))
        out.append(np.sum(c * k) / np.sum((1 - c) * k))
    return np.array(out)


R0_REL_UNC = 0.03 / 0.88   # R0 실측 불확도(±0.03)의 상대값 — 강도비 전체에 곱해지는 상수


def angle_error_band(d_step, angles, dtheta=0.5):
    """검출 각도가 ±dtheta 틀렸을 때 강도비의 상대 변화 (Seah 2005가 꼽은 0.5°)."""
    base = arxps_ratio(profile_step(d_step), angles)
    hi = arxps_ratio(profile_step(d_step), np.asarray(angles) + dtheta)
    lo = arxps_ratio(profile_step(d_step), np.asarray(angles) - dtheta)
    return np.maximum(np.abs(hi / base - 1), np.abs(lo / base - 1))


def counting_band(d_step, angles, n_sub0=1e5):
    """계수 잡음에서 오는 강도비의 상대 표준편차. 각도가 커지면 기판 신호가 줄어 커진다."""
    out = []
    for t in np.atleast_1d(angles):
        i_o, i_m = signals(d_step, n_sub0, t)
        out.append(np.sqrt(1 / i_o + 1 / i_m))
    return np.array(out)


def total_band(d_step, angles, n_sub0=1e5):
    """각도 오차와 계수 잡음을 합친 강도비의 상대 불확도."""
    return np.hypot(angle_error_band(d_step, angles), counting_band(d_step, angles, n_sub0))


def best_graded_match(d_step, width, angles, n_sub0=1e5):
    """급격한 계면(d_step)과 각도 응답이 가장 가까운 확산 계면을 찾는다.
    확산 계면의 중심과, R0 불확도 안에서 강도비 전체에 곱해지는 상수를 함께 조정해
    '불확도 띠 대비 최대 차이'를 최소화한다. 반환: (중심, 상수, 각도별 상대차, 띠)"""
    from scipy.optimize import minimize
    target = arxps_ratio(profile_step(d_step), angles)
    band = total_band(d_step, angles, n_sub0)

    def rel(p):
        c0, s = p
        return s * arxps_ratio(profile_graded(c0, width), angles) / target - 1

    def cost(p):
        if abs(p[1] - 1) > R0_REL_UNC:
            return 1e3 + abs(p[1] - 1)
        return np.max(np.abs(rel(p)) / band)

    best = min((minimize(cost, [c, 1.0], method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6})
                for c in np.linspace(d_step - 0.3, d_step + 0.8, 8)), key=lambda r: r.fun)
    return best.x[0], best.x[1], rel(best.x), band


# ---------------------------------------------------------------------------
# 스퍼터 깊이 프로파일 개념 모사 (파라미터는 설명용 임의값)
# ---------------------------------------------------------------------------
def sputter_profile(depth, t_ox=5.0, mix_sigma=1.0, lam_info=2.0, pref_loss=0.3):
    """MO2 산화막(두께 t_ox) / 금속 M의 산소 원자 분율.
    true: 참 분율. mixed: 원자 혼합(Gaussian) 뒤. observed: 정보깊이 지수 가중까지.
    reduced: 산소 우선 제거로 산화막 영역의 산소가 pref_loss만큼 빠진 관측값."""
    x_true = np.where(depth < t_ox, 2 / 3, 0.0)
    dx = depth[1] - depth[0]
    g = np.exp(-0.5 * (np.arange(-5 * mix_sigma, 5 * mix_sigma + dx, dx) / mix_sigma) ** 2)
    n_pad = g.size // 2                    # 표면 위를 진공(0)으로 채우면 표면 농도가 인위적으로 깎이므로
    padded = np.pad(x_true, n_pad, mode="edge")   # 가장자리 값을 연장해 채운다
    mixed = np.convolve(padded, g / g.sum(), mode="same")[n_pad:n_pad + x_true.size]
    k = np.exp(-np.arange(0, 6 * lam_info, dx) / lam_info)
    k /= k.sum()
    observed = np.array([np.sum(mixed[i:i + k.size] * k[:mixed[i:].size]) / k[:mixed[i:].size].sum()
                         for i in range(depth.size)])
    reduced = observed * (1 - pref_loss)
    return x_true, mixed, observed, reduced


def main():
    print(f"=== 기준: L = {L} nm, R0 = {R0} (실측), 수직 방출")
    print(f"  3λ = {3 * L:.1f} nm에서 기판 신호는 {100 * np.exp(-3):.1f}%로 준다")

    print("\n=== 2절: 두께 역산 예제 (측정 강도비 R = 1.00, 수직 방출)")
    for r0, lab in ((R0, "실측 R0 0.88"), (R0_CALC, "계산 R0 0.53")):
        print(f"  {lab}: d = {thickness(1.0, r0=r0):.2f} nm")
    print(f"  각도 혼동: 분석기가 법선에서 30°일 때 cos30 → {thickness(1.0, 30):.2f} nm, "
          f"잘못 넣은 cos60(=sin30) → {thickness(1.0, 60):.2f} nm")

    print("\n=== 그림 2: 두께 오차 (λ = 3 nm)")
    print("    d   | 강도비 5%  | 계수 1e4   | 계수 1e5   | λ+10%   | R0 0.53 사용")
    for d in (0.5, 1, 2, 4, 6, 8, 10):
        r = ratio(d)
        print(f"  {d:5.1f} | {ratio_error(d):+.3f} nm | {counting_sigma(d, 1e4):.3f} nm | "
              f"{counting_sigma(d, 1e5):.3f} nm | {thickness(r, lam=1.1 * L) - d:+.3f} | "
              f"{thickness(r, r0=R0_CALC):.2f} nm ({100 * (thickness(r, r0=R0_CALC) / d - 1):+.0f}%)")

    lo, lm, dd = contamination_effect()
    print(f"\n=== 오염층 1 nm (Si 2p 두 성분의 탄화수소층 속 λ: {lo:.3f}, {lm:.3f} nm)")
    print(f"  두께 2 nm 역산값 변화: {dd * 1000:+.2f} pm")

    print("\n=== 그림 3: ARXPS (급격한 계면 2 nm vs 확산 계면, 0–70°)")
    angles = np.arange(0, 71, 5)
    band = total_band(2.0, angles)
    print(f"  불확도 띠(각도 ±0.5° + 계수 1e5): 0° {100 * band[0]:.2f}%, 40° {100 * band[8]:.2f}%, "
          f"70° {100 * band[-1]:.2f}%   / R0 상대 불확도 ±{100 * R0_REL_UNC:.1f}%")
    for w in (0.3, 0.5, 0.8, 1.2):
        c0, sc, rel, bnd = best_graded_match(2.0, w, angles)
        over = angles[np.abs(rel) > bnd]
        print(f"  계면 폭 {w:.1f} nm: 중심 {c0:.2f} nm, 상수 {100 * (sc - 1):+.2f}%, "
              f"최대 |차| {100 * np.max(np.abs(rel)):.2f}%, 띠를 넘는 각도 "
              + (", ".join(f"{a:.0f}" for a in over) if over.size else "없음"))


if __name__ == "__main__":
    main()

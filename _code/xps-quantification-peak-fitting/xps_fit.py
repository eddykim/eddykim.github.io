"""XPS 기초 4편 — 배경 모델, 제약 피팅, 과적합 실험.

실행: python xps_fit.py
      7절의 수치(참값 vs 추정값, 제약 유무 비교, 성분 수에 따른 변화)를 출력한다.
      generate_figures.py가 이 모듈을 불러 그림을 그린다.

합성 스펙트럼은 SiO2/Si의 Si 2p를 흉내 낸다. 스핀-궤도 분리 0.61 eV, 강도비 1/2,
중간 산화물의 화학이동(Si0 기준 +0.98, +1.82, +2.65 eV)은
Z. H. Lu et al., J. Appl. Phys. 77, 4110 (1995)를 따랐다.
참 배경은 Tougaard 범용 단면적(B·T/(C+T²)², C = 1643 eV²)으로 만든 손실 배경이다.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.special import voigt_profile

FWHM_TO_SIGMA = 1 / (2 * np.sqrt(2 * np.log(2)))
E = np.arange(95.0, 108.0, 0.05)          # 결합에너지 축 (오름차순)
SO_SPLIT, SO_RATIO = 0.61, 0.5             # Si 2p 스핀-궤도 분리(eV), 2p1/2 : 2p3/2
LORENTZ_FWHM = 0.10                        # 수명 폭은 고정 (Si 2p는 작다)
TOUGAARD_C = 1643.0                        # eV², 범용 단면적

# 참값: 성분별 (2p3/2 위치 eV, 이중선 전체 면적, Gaussian FWHM eV)
TRUE = {"Si0": (99.40, 1000.0, 0.55), "Si4+": (103.30, 600.0, 1.30)}
TRUE_OX_FRACTION = TRUE["Si4+"][1] / (TRUE["Si0"][1] + TRUE["Si4+"][1])
COUNTS_PER_UNIT = 1.0                      # 면적 1 = 계수 1 eV
BG_OFFSET = 60.0

# 중간 산화물까지 포함한 문헌 위치 (Si0 기준 화학이동, Lu 1995). 성분 판정에 쓴다.
SI0 = TRUE["Si0"][0]
REFERENCE = {"Si0": SI0, "Si1+": SI0 + 0.98, "Si2+": SI0 + 1.82,
             "Si3+": SI0 + 2.65, "Si4+": SI0 + 3.84}


# ---------------------------------------------------------------------------
# 선형
# ---------------------------------------------------------------------------
def voigt(x, center, area, g_fwhm, l_fwhm=LORENTZ_FWHM):
    """면적이 area인 Voigt (Gaussian FWHM, Lorentzian FWHM)."""
    return area * voigt_profile(x - center, g_fwhm * FWHM_TO_SIGMA, l_fwhm / 2)


def doublet(x, c32, area, g_fwhm):
    """스핀-궤도 이중선. 2p1/2는 독립 파라미터가 아니라 2p3/2에서 유도한다.
    area는 두 성분을 합친 면적, c32는 2p3/2 위치."""
    a32 = area / (1 + SO_RATIO)
    return voigt(x, c32, a32, g_fwhm) + voigt(x, c32 + SO_SPLIT, a32 * SO_RATIO, g_fwhm)


# ---------------------------------------------------------------------------
# 배경
# ---------------------------------------------------------------------------
def _cumarea(x, y):
    """낮은 결합에너지 끝에서부터 x까지의 누적 면적 (사다리꼴)."""
    return np.concatenate([[0.0], np.cumsum((y[1:] + y[:-1]) / 2 * np.diff(x))])


def shirley(x, y, n_end=5, tol=1e-7, max_iter=100):
    """Shirley 배경. x는 결합에너지 오름차순.
    어떤 결합에너지에서의 배경 높이 = 낮은 끝 값 + (높은 끝 - 낮은 끝) ×
    (그 에너지보다 낮은 결합에너지 쪽 피크 면적 / 전체 피크 면적). 자기참조라 반복한다."""
    lo, hi = y[:n_end].mean(), y[-n_end:].mean()
    bg = np.full_like(y, lo)
    for it in range(max_iter):
        cum = _cumarea(x, np.clip(y - bg, 0, None))
        new = lo + (hi - lo) * cum / cum[-1]
        if np.max(np.abs(new - bg)) < tol * hi:
            return new, it + 1
        bg = new
    return bg, max_iter


def linear(x, y, n_end=5):
    lo, hi = y[:n_end].mean(), y[-n_end:].mean()
    return lo + (hi - lo) * (x - x[0]) / (x[-1] - x[0])


def tougaard_loss(x, signal, B, C=TOUGAARD_C):
    """Tougaard 범용 단면적 F(T) = B·T/(C+T²)²로 signal이 만드는 손실 배경.
    결합에너지 b에서의 배경은 그보다 낮은 결합에너지 b'의 신호가 T = b - b'만큼 잃은 몫의 합."""
    dx = x[1] - x[0]
    T = np.arange(x.size) * dx
    kernel = B * T / (C + T**2) ** 2
    return np.convolve(signal, kernel)[: x.size] * dx


def tougaard(x, y, n_end=5, C=TOUGAARD_C):
    """데이터에서 Tougaard 배경을 추정. 낮은 끝 값을 상수로 빼고, B는 높은 끝에서 데이터와
    만나도록 정한다(CasaXPS와 같은 방식)."""
    lo, hi = y[:n_end].mean(), y[-n_end:].mean()
    unit = tougaard_loss(x, y - lo, 1.0, C)
    B = (hi - lo) / unit[-n_end:].mean()
    return lo + B * unit, B


# ---------------------------------------------------------------------------
# 합성 스펙트럼
# ---------------------------------------------------------------------------
def true_signal(x=E):
    return sum(doublet(x, *p) for p in TRUE.values())


B_UNIVERSAL = 2866.0                       # eV², 금속에 맞춘 범용값
B_MULT = 6.0                               # 합성용 배율. 1이면 이 창 안의 계단이 피크의 1% 남짓이라
                                           # 배경 모델 차이가 거의 안 보인다. 민감도는 main()에서 따로 본다.


def true_background(x=E, mult=None):
    m = B_MULT if mult is None else mult
    return BG_OFFSET + tougaard_loss(x, true_signal(x), B_UNIVERSAL * m)


def synthetic(seed, mult=None):
    """포아송 잡음을 얹은 합성 스펙트럼 (계수)."""
    lam = (true_signal() + true_background(mult=mult)) * COUNTS_PER_UNIT
    return np.random.default_rng(seed).poisson(lam).astype(float)


def area(x, y):
    return float(_cumarea(x, y)[-1])


# ---------------------------------------------------------------------------
# 피팅
# ---------------------------------------------------------------------------
def _weights(y_raw):
    return 1.0 / np.sqrt(np.clip(y_raw, 1, None))          # 포아송 가중


def fit_constrained(x, y_sub, y_raw):
    """제약 모델: Si0 이중선 + Si4+ 이중선. 파라미터 6개(성분마다 위치·면적·폭)."""
    w = _weights(y_raw)
    p0 = [99.5, 800.0, 0.8, 103.2, 500.0, 1.5]
    lo = [98.8, 0, 0.2, 102.3, 0, 0.3]
    hi = [100.0, 1e5, 2.0, 104.3, 1e5, 3.0]

    def resid(p):
        return (doublet(x, *p[:3]) + doublet(x, *p[3:]) - y_sub) * w

    r = least_squares(resid, p0, bounds=(lo, hi), method="trf")
    chi2 = float(np.sum(r.fun**2))
    dof = x.size - len(p0)
    cov = np.linalg.inv(r.jac.T @ r.jac) * chi2 / dof
    return r.x, np.sqrt(np.diag(cov)), chi2 / dof


def fit_free(x, y_sub, y_raw, n, rng, starts=12):
    """제약 없는 모델: 위치·면적·폭이 모두 자유인 단일선 n개. 초기값을 여러 번 바꿔
    가장 작은 χ²를 고른다."""
    w = _weights(y_raw)
    lo = np.tile([96.5, 0, 0.2], n)
    hi = np.tile([106.5, 1e5, 3.0], n)
    best = None
    for _ in range(starts):
        cs = np.sort(rng.uniform(98.5, 104.5, n))
        p0 = np.ravel([[c, 1600.0 / n, 0.9] for c in cs])

        def resid(p):
            return (sum(voigt(x, *p[3 * i:3 * i + 3]) for i in range(n)) - y_sub) * w

        r = least_squares(resid, p0, bounds=(lo, hi), method="trf")
        chi2 = float(np.sum(r.fun**2))
        if best is None or chi2 < best[1]:
            best = (r.x.reshape(n, 3), chi2)
    return best[0], best[1] / (x.size - 3 * n)


def assign(components):
    """성분마다 가장 가까운 문헌 위치의 화학상태를 붙이고, 상태별 면적 몫을 돌려준다.
    분석자가 제약 없는 피팅 결과를 문헌표와 대조해 읽는 방식을 흉내 낸다."""
    total = components[:, 1].sum()
    share = {k: 0.0 for k in REFERENCE}
    for c, a, _ in components:
        k = min(REFERENCE, key=lambda s: abs(REFERENCE[s] - c))
        share[k] += a / total
    return share


# ---------------------------------------------------------------------------
# 실험
# ---------------------------------------------------------------------------
def background_comparison(seed=0, mult=None):
    """그림 1: 같은 스펙트럼에 배경 모델 셋을 적용해 피크 면적을 비교."""
    y = synthetic(seed, mult)
    bgs = {"linear": linear(E, y), "Shirley": shirley(E, y)[0], "Tougaard": tougaard(E, y)[0]}
    areas = {k: area(E, y - b) for k, b in bgs.items()}
    return y, bgs, areas, area(E, true_signal())


def constrained_demo(seed=0, mult=None):
    """그림 3·7.3절: Shirley 배경을 빼고 제약 모델로 피팅.
    대조군으로 참 배경을 정확히 알고 뺐을 때의 피팅도 함께 돌려준다."""
    y = synthetic(seed, mult)
    bg, n_iter = shirley(E, y)
    p, err, chi2r = fit_constrained(E, y - bg, y)
    p_ex, err_ex, chi2r_ex = fit_constrained(E, y - true_background(mult=mult), y)
    return y, bg, n_iter, p, err, chi2r, (p_ex, err_ex, chi2r_ex)


def sweep(n_values=range(2, 8), realizations=20, starts=12, seed0=100):
    """그림 4·7.4절: 성분 수 n을 바꿔 가며 제약 없는 피팅을 잡음 실현마다 반복."""
    rng = np.random.default_rng(7)
    out = {n: {"chi2": [], "chi2r": [], "fake": [], "ox": [], "species": []} for n in n_values}
    con = {"chi2": [], "chi2r": [], "ox": []}
    for k in range(realizations):
        y = synthetic(seed0 + k)
        ys = y - shirley(E, y)[0]
        p, _, c2 = fit_constrained(E, ys, y)
        con["chi2r"].append(c2)
        con["chi2"].append(c2 * (E.size - 6))
        con["ox"].append(p[4] / (p[1] + p[4]))
        for n in n_values:
            comps, c2n = fit_free(E, ys, y, n, rng, starts)
            sh = assign(comps)
            out[n]["chi2r"].append(c2n)
            out[n]["chi2"].append(c2n * (E.size - 3 * n))
            out[n]["fake"].append(sh["Si1+"] + sh["Si2+"] + sh["Si3+"])
            out[n]["ox"].append(sh["Si4+"])
            out[n]["species"].append(sum(v > 0.02 for v in sh.values()))
    return out, con


def main():
    y, bgs, areas, a_true = background_comparison()
    print("=== 4절: 배경 모델에 따른 피크 면적 (참 면적 대비)")
    for k, a in areas.items():
        print(f"  {k:9s} {a:8.1f}  ({100 * (a / a_true - 1):+.1f}%)")
    print(f"  참 면적   {a_true:8.1f}")
    print("  배경 세기에 따른 민감도 (배율: 계단 높이/피크 높이 → 면적 오차 linear/Shirley/Tougaard)")
    for m in (1.0, 3.0, 6.0):
        bgt = true_background(mult=m)
        step = 100 * (bgt[-1] - bgt[0]) / true_signal().max()
        _, _, ar, at = background_comparison(mult=m)
        errs = "/".join(f"{100 * (ar[k] / at - 1):+.1f}%" for k in ("linear", "Shirley", "Tougaard"))
        print(f"    ×{m:.0f}: 계단 {step:.1f}% → {errs}")

    y, bg, n_iter, p, err, chi2r, (p_ex, err_ex, c2_ex) = constrained_demo()
    print(f"\n=== 7.3절: 제약 피팅 (Shirley {n_iter}회 반복 수렴, χ²_red = {chi2r:.2f})")
    names = ["Si0 위치", "Si0 면적", "Si0 FWHM", "Si4+ 위치", "Si4+ 면적", "Si4+ FWHM"]
    truth = [*TRUE["Si0"], *TRUE["Si4+"]]
    for nm, t, v, e in zip(names, truth, p, err):
        print(f"  {nm:10s} 참 {t:8.2f}  추정 {v:8.2f} ± {e:.2f}")
    print(f"  Si4+ 몫: 참 {100 * TRUE_OX_FRACTION:.1f}%  추정 {100 * p[4] / (p[1] + p[4]):.1f}%")
    print(f"  대조군(참 배경을 정확히 뺌): χ²_red {c2_ex:.2f}, Si0 면적 {p_ex[1]:.1f} ± {err_ex[1]:.1f}, "
          f"Si4+ 면적 {p_ex[4]:.1f} ± {err_ex[4]:.1f}, Si4+ FWHM {p_ex[5]:.2f}, "
          f"Si4+ 몫 {100 * p_ex[4] / (p_ex[1] + p_ex[4]):.1f}%")

    rng = np.random.default_rng(3)
    ys = y - bg
    print("\n=== 7.4절: 같은 데이터, 스핀-궤도 제약 없이 단일선 n개")
    for n in (2, 3, 5):
        comps, c2 = fit_free(E, ys, y, n, rng)
        sh = assign(comps)
        cs = ", ".join(f"{c:.2f}" for c in np.sort(comps[:, 0]))
        print(f"  n={n}: χ²_red {c2:.2f} | 위치 [{cs}]")
        print("        " + ", ".join(f"{k} {100 * v:.1f}%" for k, v in sh.items() if v > 0.001))

    print("\n=== 그림 4: 성분 수 스윕 (잡음 실현 20회, 평균 ± 표준편차)")
    out, con = sweep()
    print(f"  제약 모델: χ² {np.mean(con['chi2']):.0f}, χ²_red {np.mean(con['chi2r']):.2f} ± {np.std(con['chi2r']):.2f}, "
          f"Si4+ {100 * np.mean(con['ox']):.1f} ± {100 * np.std(con['ox']):.1f}%")
    for n, d in out.items():
        print(f"  n={n}: χ² {np.mean(d['chi2']):.0f}, χ²_red {np.mean(d['chi2r']):.2f} ± {np.std(d['chi2r']):.2f}, "
              f"가짜 중간산화물 {100 * np.mean(d['fake']):.1f} ± {100 * np.std(d['fake']):.1f}%, "
              f"Si4+ {100 * np.mean(d['ox']):.1f} ± {100 * np.std(d['ox']):.1f}%, "
              f"보고되는 화학종 {np.mean(d['species']):.1f}개")


if __name__ == "__main__":
    main()

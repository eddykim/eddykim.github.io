"""이미지처리 5편 그림 6개 생성 (한국어판·영문판)."""
import os

import matplotlib.pyplot as plt
import numpy as np
from skimage.data import coins
from skimage.feature import canny
from skimage.measure import label, regionprops
from scipy.optimize import minimize

from circle_fit import FITTERS, arc_points, geometric, kasa
from ellipse_fit import (conic_type, ellipse_params, ellipse_points,
                         fit_conic_unconstrained, fit_ellipse)
from ransac import contaminate, ransac_circle

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
BASE = os.path.join(os.path.dirname(__file__), "..", "..",
                    "assets", "img", "posts", "imgproc-geometric-fitting")

NK = {"kasa": "Kasa (대수적)", "pratt": "Pratt", "taubin": "Taubin", "geometric": "기하학적 (LM)"}
NE = {"kasa": "Kasa (algebraic)", "pratt": "Pratt", "taubin": "Taubin", "geometric": "Geometric (LM)"}
COL = {"kasa": "#c0392b", "pratt": "#e67e22", "taubin": "#1e8449", "geometric": "#2b5c9b"}
MK = {"kasa": "o", "pratt": "v", "taubin": "s", "geometric": "^"}

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT, "names": NK,
        "radius": "반경 $R$", "cost": "정규화한 목적함수", "alg": "대수적 거리", "geo": "기하학적 거리",
        "truth": "참 반경", "pts": "엣지 점", "fitted": "맞춘 원",
        "fig1": "그림1. 대수적 목적함수의 최솟값은 참 반경보다 작은 쪽에 있다",
        "span": "원호 각도 (도)", "bias": "반경 편향", "scatter": "반경 산포",
        "fig2": "그림2. 짧은 원호에서 Kasa 법은 반경을 체계적으로 작게 잡는다",
        "rmse": "반경 총 오차 RMSE", "unusable": "어느 방법도\n쓸 수 없는 구간",
        "fig3": "그림3. 20도 아래에서는 편향이 작은 방법이 오히려 더 크게 틀린다",
        "noise_ax": "점의 노이즈", "ellipse_rate": "타원이 나온 비율 (%)",
        "constrained": "Fitzgibbon (제약 있음)", "unconstrained": "제약 없음",
        "example": "제약 없는 피팅이 쌍곡선을 낸 예",
        "fig4": "그림4. 제약이 없으면 이차곡선 피팅은 타원을 보장하지 않는다",
        "outlier": "이상치 비율 (%)", "rerr": "반경 오차",
        "ls": "최소자승 (전체 점)", "rn": "RANSAC + 기하학적",
        "inlier": "정상치", "outl": "이상치", "scene": "이상치가 섞인 점 집합",
        "fig5": "그림5. 이상치가 섞이면 최소자승은 끌려가고 RANSAC 은 버틴다",
        "coin_all": "전체 윤곽 226점", "coin_arc": "45도 호 32점만 남겼을 때",
        "fig6": "그림6. 실제 영상에서도 같은 편향이 나타난다",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {}, "names": NE,
        "radius": "Radius $R$", "cost": "Normalized objective", "alg": "Algebraic distance",
        "geo": "Geometric distance", "truth": "True radius", "pts": "Edge points",
        "fitted": "Fitted circle",
        "fig1": "Fig 1. The algebraic objective has its minimum below the true radius",
        "span": "Arc span (degrees)", "bias": "Radius bias", "scatter": "Radius scatter",
        "fig2": "Fig 2. On short arcs the Kasa method systematically underestimates the radius",
        "rmse": "Total radius error RMSE", "unusable": "no method\nis usable here",
        "fig3": "Fig 3. Below 20 degrees the unbiased methods are the ones that fail",
        "noise_ax": "Point noise", "ellipse_rate": "Fraction returning an ellipse (%)",
        "constrained": "Fitzgibbon (constrained)", "unconstrained": "Unconstrained",
        "example": "An unconstrained fit returning a hyperbola",
        "fig4": "Fig 4. Without the constraint a conic fit does not guarantee an ellipse",
        "outlier": "Outlier fraction (%)", "rerr": "Radius error",
        "ls": "Least squares (all points)", "rn": "RANSAC + geometric",
        "inlier": "Inliers", "outl": "Outliers", "scene": "Point set with outliers",
        "fig6": "Fig 6. The same bias appears on a real image",
        "coin_all": "Full contour, 226 points", "coin_arc": "Only a 45-degree arc, 32 points",
        "fig5": "Fig 5. Outliers drag the least-squares fit while RANSAC holds",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
R_TRUE, CX, CY = 40.0, 50.0, 30.0
N_PTS, NOISE, TRIALS = 120, 0.4, 400

# 그림1: 두 목적함수를 반경의 함수로. 중심은 각 R 마다 다시 최적화한다.
rng1 = np.random.default_rng(11)
x1, y1 = arc_points(60, CX, CY, R_TRUE, 60.0, noise=NOISE * 2, rng=rng1)
R_scan = np.linspace(20.0, 70.0, 260)


def _cost(R, mode):
    """반경을 R 로 고정하고 중심만 최적화했을 때의 목적함수 값."""
    def obj(c):
        d = np.hypot(x1 - c[0], y1 - c[1])
        return np.sum((d ** 2 - R ** 2) ** 2) if mode == "alg" else np.sum((d - R) ** 2)
    best = np.inf
    for guess in ([CX, CY], [x1.mean(), y1.mean()]):
        r = minimize(obj, guess, method="Nelder-Mead",
                     options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 4000})
        best = min(best, float(r.fun))
    return best


cost_alg = np.array([_cost(R, "alg") for R in R_scan])
cost_geo = np.array([_cost(R, "geo") for R in R_scan])
cost_alg /= cost_alg.min()
cost_geo /= cost_geo.min()
R_alg = R_scan[int(np.argmin(cost_alg))]
R_geo = R_scan[int(np.argmin(cost_geo))]

# 그림2~3: 원호 각도 스윕
SPANS = np.array([360, 240, 180, 120, 90, 60, 45, 30, 20, 15, 10])
bias = {k: [] for k in FITTERS}
sctr = {k: [] for k in FITTERS}
rmse = {k: [] for k in FITTERS}
for sp in SPANS:
    acc = {k: [] for k in FITTERS}
    rng = np.random.default_rng(3)
    for _ in range(TRIALS):
        x, y = arc_points(N_PTS, CX, CY, R_TRUE, float(sp), noise=NOISE, rng=rng)
        for k, f in FITTERS.items():
            try:
                acc[k].append(f(x, y)[2])
            except Exception:
                acc[k].append(np.nan)
    for k in FITTERS:
        a = np.array(acc[k], float)
        bias[k].append(np.nanmean(a) - R_TRUE)
        sctr[k].append(np.nanstd(a))
        rmse[k].append(np.sqrt(np.nanmean((a - R_TRUE) ** 2)))
bias = {k: np.array(v) for k, v in bias.items()}
sctr = {k: np.array(v) for k, v in sctr.items()}
rmse = {k: np.array(v) for k, v in rmse.items()}

# 그림4: 제약의 효과
E_NOISE = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0])
rng4 = np.random.default_rng(0)
rate_c, rate_u = [], []
for nz in E_NOISE:
    c = u = 0
    for _ in range(300):
        ex, ey = ellipse_points(80, 10.0, 20.0, 60.0, 35.0, np.deg2rad(30),
                                span_deg=120, noise=nz, rng=rng4)
        c += conic_type(fit_ellipse(ex, ey)) == "타원"
        u += conic_type(fit_conic_unconstrained(ex, ey)) == "타원"
    rate_c.append(100.0 * c / 300)
    rate_u.append(100.0 * u / 300)
# 쌍곡선이 나온 실례 하나를 찾는다
rng4b = np.random.default_rng(4)
demo_x = demo_y = None
for _ in range(400):
    ex, ey = ellipse_points(80, 10.0, 20.0, 60.0, 35.0, np.deg2rad(30),
                            span_deg=120, noise=3.0, rng=rng4b)
    if conic_type(fit_conic_unconstrained(ex, ey)) != "타원":
        demo_x, demo_y = ex, ey
        break
demo_fit = ellipse_params(fit_ellipse(demo_x, demo_y)) if demo_x is not None else None
demo_unc = fit_conic_unconstrained(demo_x, demo_y) if demo_x is not None else None

# 그림5: RANSAC
FRACS = np.array([0.0, 0.05, 0.10, 0.20, 0.30, 0.40])
rng5 = np.random.default_rng(1)
err_ls, err_rn = [], []
for fr in FRACS:
    a, b = [], []
    for _ in range(120):
        x, y = arc_points(N_PTS, CX, CY, R_TRUE, 360.0, noise=NOISE, rng=rng5)
        if fr > 0:
            x, y, _ = contaminate(x, y, fr, 70.0, rng5, CX, CY)
        a.append(geometric(x, y)[2] - R_TRUE)
        b.append(ransac_circle(x, y, threshold=1.5, n_iter=200, seed=0)[2] - R_TRUE)
    err_ls.append(np.mean(a))
    err_rn.append(np.mean(b))
rng5b = np.random.default_rng(9)
sx, sy = arc_points(N_PTS, CX, CY, R_TRUE, 360.0, noise=NOISE, rng=rng5b)
sx, sy, out_idx = contaminate(sx, sy, 0.25, 70.0, rng5b, CX, CY)
s_ls = geometric(sx, sy)
s_cx, s_cy, s_R, s_inl = ransac_circle(sx, sy, threshold=1.5, n_iter=300, seed=0)

# 그림6: 실제 coins 영상
img = coins().astype(float)
edges = canny(img / 255.0, sigma=2.0, low_threshold=0.08, high_threshold=0.2)
prop = sorted(regionprops(label(edges, connectivity=2)), key=lambda q: -q.area)[0]
cy_pts = prop.coords[:, 0].astype(float)
cx_pts = prop.coords[:, 1].astype(float)
full = geometric(cx_pts, cy_pts)
ang = np.rad2deg(np.arctan2(cy_pts - full[1], cx_pts - full[0])) % 360
arc_m = (ang >= 20) & (ang < 65)
coin_fits = {k: FITTERS[k](cx_pts[arc_m], cy_pts[arc_m]) for k in FITTERS}


def circ(a, b, R, n=200):
    t = np.linspace(0, 2 * np.pi, n)
    return a + R * np.cos(t), b + R * np.sin(t)


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    F, LG, NM = L["font"], L["legend"], L["names"]

    def cap(fig, key, rect=(0, 0.06, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **F)
        fig.tight_layout(rect=rect)

    # 그림1
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 4.0))
    ax[0].plot(x1, y1, ".", ms=4, color="#555", label=L["pts"])
    for k, R in (("kasa", R_alg), ("geometric", R_geo)):
        pass
    ka, kb, kR = kasa(x1, y1)
    ga, gb, gR = geometric(x1, y1)
    cxk, cyk = circ(ka, kb, kR); ax[0].plot(cxk, cyk, "-", lw=1.4, color=COL["kasa"],
                                            label=f"{NM['kasa']} R={kR:.1f}")
    cxg, cyg = circ(ga, gb, gR); ax[0].plot(cxg, cyg, "-", lw=1.4, color=COL["geometric"],
                                            label=f"{NM['geometric']} R={gR:.1f}")
    cxt, cyt = circ(CX, CY, R_TRUE); ax[0].plot(cxt, cyt, ":", lw=1.2, color="k",
                                                label=f"{L['truth']} {R_TRUE:.0f}")
    ax[0].set_aspect("equal"); ax[0].legend(prop=LG, fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(R_scan, cost_alg, "-", lw=1.7, color=COL["kasa"], label=L["alg"])
    ax[1].plot(R_scan, cost_geo, "-", lw=1.7, color=COL["geometric"], label=L["geo"])
    ax[1].axvline(R_TRUE, ls=":", lw=1.2, color="k")
    ax[1].text(R_TRUE + 0.8, cost_alg.max() * 0.8, L["truth"], fontsize=8.5, **F)
    ax[1].set_xlabel(L["radius"], **F); ax[1].set_ylabel(L["cost"], **F)
    ax[1].set_yscale("log"); ax[1].legend(prop=LG, fontsize=9); ax[1].grid(alpha=0.3)
    cap(fig, "fig1", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig1-algebraic-vs-geometric.png"), dpi=150)
    plt.close(fig)

    # 그림2
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.0))
    for k in FITTERS:
        ax[0].plot(SPANS, bias[k], MK[k] + "-", ms=4, color=COL[k], label=NM[k])
        ax[1].loglog(SPANS, np.maximum(sctr[k], 1e-4), MK[k] + "-", ms=4,
                     color=COL[k], label=NM[k])
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_xlabel(L["span"], **F); ax[0].set_ylabel(L["bias"], **F)
    ax[0].set_xscale("log"); ax[0].legend(prop=LG, fontsize=8.5); ax[0].grid(alpha=0.3, which="both")
    ax[1].set_xlabel(L["span"], **F); ax[1].set_ylabel(L["scatter"], **F)
    ax[1].legend(prop=LG, fontsize=8.5); ax[1].grid(alpha=0.3, which="both")
    cap(fig, "fig2", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig2-arc-span-bias.png"), dpi=150)
    plt.close(fig)

    # 그림3
    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    for k in FITTERS:
        ax.loglog(SPANS, rmse[k], MK[k] + "-", ms=4, color=COL[k], label=NM[k])
    ax.axvspan(SPANS.min(), 20, color="#c0392b", alpha=0.08)
    ax.text(12, rmse["geometric"].max() * 0.05, L["unusable"], fontsize=8.5,
            color="#c0392b", ha="center", **F)
    ax.set_xlabel(L["span"], **F); ax.set_ylabel(L["rmse"], **F)
    ax.legend(prop=LG, fontsize=9); ax.grid(alpha=0.3, which="both")
    cap(fig, "fig3")
    fig.savefig(os.path.join(out, "fig3-total-error.png"), dpi=150)
    plt.close(fig)

    # 그림4
    fig, ax = plt.subplots(1, 2, figsize=(10.0, 4.0))
    ax[0].plot(E_NOISE, rate_c, "s-", ms=5, color="#1e8449", label=L["constrained"])
    ax[0].plot(E_NOISE, rate_u, "o-", ms=5, color="#c0392b", label=L["unconstrained"])
    ax[0].set_xlabel(L["noise_ax"], **F); ax[0].set_ylabel(L["ellipse_rate"], **F)
    ax[0].set_ylim(-5, 105); ax[0].legend(prop=LG, fontsize=9); ax[0].grid(alpha=0.3)
    if demo_x is not None:
        ax[1].plot(demo_x, demo_y, ".", ms=4, color="#555", label=L["pts"])
        if demo_fit is not None:
            t = np.linspace(0, 2 * np.pi, 300)
            cx0, cy0, aa, bb, th = demo_fit
            ax[1].plot(cx0 + aa * np.cos(t) * np.cos(th) - bb * np.sin(t) * np.sin(th),
                       cy0 + aa * np.cos(t) * np.sin(th) + bb * np.sin(t) * np.cos(th),
                       "-", lw=1.5, color="#1e8449", label=L["constrained"])
        gx, gy = np.meshgrid(np.linspace(demo_x.min() - 60, demo_x.max() + 60, 400),
                             np.linspace(demo_y.min() - 60, demo_y.max() + 60, 400))
        A, B, C, D, E, Fc = demo_unc
        ax[1].contour(gx, gy, A * gx ** 2 + B * gx * gy + C * gy ** 2 + D * gx + E * gy + Fc,
                      levels=[0.0], colors=["#c0392b"], linewidths=1.5)
        ax[1].plot([], [], "-", color="#c0392b", label=L["unconstrained"])
        ax[1].set_title(L["example"], fontsize=9.5, **F)
        ax[1].set_aspect("equal"); ax[1].legend(prop=LG, fontsize=8); ax[1].grid(alpha=0.3)
    cap(fig, "fig4", rect=(0, 0.07, 1, 0.95))
    fig.savefig(os.path.join(out, "fig4-ellipse-constraint.png"), dpi=150)
    plt.close(fig)

    # 그림5
    fig, ax = plt.subplots(1, 2, figsize=(10.2, 4.0))
    ax[0].plot(sx[s_inl], sy[s_inl], ".", ms=5, color="#2b5c9b", label=L["inlier"])
    ax[0].plot(sx[~s_inl], sy[~s_inl], "x", ms=5, color="#c0392b", label=L["outl"])
    a, b = circ(*s_ls); ax[0].plot(a, b, "-", lw=1.5, color="#c0392b", label=L["ls"])
    a, b = circ(s_cx, s_cy, s_R); ax[0].plot(a, b, "-", lw=1.5, color="#1e8449", label=L["rn"])
    ax[0].set_aspect("equal"); ax[0].set_title(L["scene"], fontsize=9.5, **F)
    ax[0].legend(prop=LG, fontsize=8); ax[0].grid(alpha=0.3)
    ax[1].plot(FRACS * 100, err_ls, "o-", ms=5, color="#c0392b", label=L["ls"])
    ax[1].plot(FRACS * 100, err_rn, "s-", ms=5, color="#1e8449", label=L["rn"])
    ax[1].axhline(0, color="k", lw=0.6)
    ax[1].set_xlabel(L["outlier"], **F); ax[1].set_ylabel(L["rerr"], **F)
    ax[1].legend(prop=LG, fontsize=9); ax[1].grid(alpha=0.3)
    cap(fig, "fig5", rect=(0, 0.07, 1, 1))
    fig.savefig(os.path.join(out, "fig5-ransac.png"), dpi=150)
    plt.close(fig)

    # 그림6
    fig, ax = plt.subplots(1, 2, figsize=(9.4, 4.3))
    for a in ax:
        a.imshow(img, cmap="gray")
        a.set_xlim(306, 384); a.set_ylim(224, 150)
        a.set_xticks([]); a.set_yticks([])
    ax[0].plot(cx_pts, cy_pts, ".", ms=2.5, color="#2b5c9b")
    a, b = circ(*full); ax[0].plot(a, b, "-", lw=1.4, color="#1e8449")
    ax[0].set_title(f"{L['coin_all']}  R={full[2]:.2f}", fontsize=9.5, **F)
    ax[1].plot(cx_pts[arc_m], cy_pts[arc_m], ".", ms=3.5, color="#2b5c9b")
    for k in ("kasa", "taubin", "geometric"):
        a, b = circ(*coin_fits[k])
        ax[1].plot(a, b, "-", lw=1.3, color=COL[k],
                   label=f"{NM[k]} R={coin_fits[k][2]:.2f}")
    ax[1].set_title(L["coin_arc"], fontsize=9.5, **F)
    ax[1].legend(prop=LG, fontsize=7.5, loc="upper center",
                 bbox_to_anchor=(0.5, -0.02), ncol=3, frameon=False)
    cap(fig, "fig6", rect=(0, 0.13, 1, 1))
    fig.savefig(os.path.join(out, "fig6-real-image.png"), dpi=150)
    plt.close(fig)


for lang, L in LABELS.items():
    draw(L)
    print(f"[{lang}] 그림 6개 저장: {L['dir']}")

print()
print("=== 그림1 ===")
print(f"  호 60도, 노이즈 {NOISE*2}. 대수적 목적함수 최소 R={R_alg:.2f}, "
      f"기하학적 최소 R={R_geo:.2f}, 참값 {R_TRUE}")
print(f"  Kasa R={kasa(x1,y1)[2]:.2f}, 기하학적 R={geometric(x1,y1)[2]:.2f}")
print()
print("=== 그림2~3: 원호 각도 스윕 ===")
print(f"{'호(도)':>7s} " + " ".join(f"{k:>22s}" for k in FITTERS))
for i, sp in enumerate(SPANS):
    print(f"{sp:7d} " + " ".join(f"{bias[k][i]:+9.3f}/{sctr[k][i]:7.3f}/{rmse[k][i]:7.3f}"
                                 for k in FITTERS))
print("  (편향/산포/RMSE)")
cross = [sp for i, sp in enumerate(SPANS) if rmse["kasa"][i] < rmse["taubin"][i]]
print(f"  Kasa 의 총 오차가 더 작아지는 구간: {cross} 도")
print()
print("=== 그림4: 타원 제약 ===")
for nz, c, u in zip(E_NOISE, rate_c, rate_u):
    print(f"  노이즈 {nz:4.1f}: Fitzgibbon {c:5.1f}%, 제약 없음 {u:5.1f}%")
print()
print("=== 그림5: RANSAC ===")
for fr, a, b in zip(FRACS, err_ls, err_rn):
    print(f"  이상치 {fr*100:4.0f}%: 최소자승 {a:+7.3f}, RANSAC {b:+7.3f}")
print(f"  예시 장면: 이상치 {int(round(0.25*N_PTS))}개, RANSAC 정상치 판정 {s_inl.sum()}개")
print()
print("=== 그림6: 실제 coins ===")
print(f"  전체 {len(cx_pts)}점: R={full[2]:.3f}, 중심=({full[0]:.2f},{full[1]:.2f})")
print(f"  45도 호 {arc_m.sum()}점: " + ", ".join(f"{k} R={coin_fits[k][2]:.3f}" for k in FITTERS))

"""기하광학 배경이론 2편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/geometric-optics-paraxial-abcd/      (한국어)
      ../../assets/img/posts/geometric-optics-paraxial-abcd/en/   (영문)

그림 1·5는 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from abcd import (cardinal_points, grin, iterate, periodic_cell, prop, stability, surface,
                  system, thin_lens)
from rays import Singlet, hit_sphere, n2_parabolic, refract, sample, trace_ray

KFONT = {"fontfamily": "AppleGothic"}
EFONT: dict = {}
plt.rcParams["axes.unicode_minus"] = False

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                        "assets", "img", "posts", "geometric-optics-paraxial-abcd")

C_RAY, C_PAR, C_KEY, C_GREY = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f"

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT,
        "f2t": "근축 근사의 상대오차",
        "f2x": r"각도 ($^\circ$)", "f2y": "상대오차 (%)",
        "f2sin": r"$\sin\theta \approx \theta$", "f2snell": "굴절각 (공기→유리, n = 1.5)",
        "f2mar": "1편 렌즈의\n가장자리 광선 5.7°",
        "p_a": "(a) 50 mm 진행: 가로로 밀린다", "p_b": "(b) 초점거리 50 mm 렌즈: 세로로 밀린다",
        "p_c": "(c) 앞 초점면 → 뒤 초점면: 90° 돈다",
        "p_x": "높이 y (mm)", "p_y": "각 u (mrad)", "p_area": "넓이",
        "p_before": "들어가기 전", "p_after": "나온 뒤",
        "f3a": "(a) 1편 렌즈의 주요점 (가로 축척 압축)",
        "f3b": "(b) 실제 광선이 꺾이는 것처럼 보이는 위치",
        "f3x": "광축 위치 z (mm)", "f3y": "높이 y (mm)",
        "f3par": "근축 주평면 H′", "f3meet": "연장선이 만나는 점",
        "f3bx": "근축 주평면 H′ 로부터의 거리 (mm)", "f3by": "입사 높이 h (mm)",
        "f5a": "(a) 한 점에서 나온 광선 (실선: 광선 방정식, 점선: 근축 행렬)",
        "f5b": "(b) 광축으로 돌아오는 거리",
        "f5x": "진행 거리 z", "f5y": "높이 y", "f5by": "첫 교차 거리",
        "f5bx": r"출발각 ($^\circ$)", "f5exact": "광선 방정식 적분", "f5cos": r"$(\pi/\alpha)\cos\theta_0$",
        "f5par": r"근축 행렬 $\pi/\alpha$",
        "f6a": "(a) 안정 조건", "f6ax": "렌즈 간격 d / 초점거리 f", "f6ay": "b = (A + D)/2",
        "f6stable": "안정 |b| ≤ 1",
        "f6b": "(b) 단마다 본 광선 높이", "f6bx": "단 번호 m", "f6by": r"높이 $y_m$",
        "f6d1": "d = f (6단 주기)", "f6d2": "d = 2f (4단 주기)",
        "f6d35": "d = 3.5f (안정, 비주기)", "f6d42": "d = 4.2f (불안정)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT,
        "f2t": "Relative error of the paraxial approximation",
        "f2x": r"Angle ($^\circ$)", "f2y": "Relative error (%)",
        "f2sin": r"$\sin\theta \approx \theta$", "f2snell": "Refraction angle (air→glass, n = 1.5)",
        "f2mar": "Marginal ray of\nthe Part 1 lens, 5.7°",
        "p_a": "(a) 50 mm of travel: a horizontal shear", "p_b": "(b) A 50 mm lens: a vertical shear",
        "p_c": "(c) Front to back focal plane: a quarter turn",
        "p_x": "Height y (mm)", "p_y": "Angle u (mrad)", "p_area": "area",
        "p_before": "before", "p_after": "after",
        "f3a": "(a) Cardinal points of the Part 1 lens (axial scale compressed)",
        "f3b": "(b) Where real rays appear to bend",
        "f3x": "Axial position z (mm)", "f3y": "Height y (mm)",
        "f3par": "Paraxial principal plane H′", "f3meet": "Where the extensions meet",
        "f3bx": "Distance from paraxial H′ (mm)", "f3by": "Ray height h (mm)",
        "f5a": "(a) Rays from one point (solid: ray equation, dashed: paraxial matrix)",
        "f5b": "(b) Distance to return to the axis",
        "f5x": "Propagation distance z", "f5y": "Height y", "f5by": "First crossing",
        "f5bx": r"Launch angle ($^\circ$)", "f5exact": "Ray equation", "f5cos": r"$(\pi/\alpha)\cos\theta_0$",
        "f5par": r"Paraxial matrix $\pi/\alpha$",
        "f6a": "(a) Stability condition", "f6ax": "Lens spacing d / focal length f",
        "f6ay": "b = (A + D)/2", "f6stable": "Stable |b| ≤ 1",
        "f6b": "(b) Ray height stage by stage", "f6bx": "Stage m", "f6by": r"Height $y_m$",
        "f6d1": "d = f (period 6)", "f6d2": "d = 2f (period 4)",
        "f6d35": "d = 3.5f (stable, aperiodic)", "f6d42": "d = 4.2f (unstable)",
    },
}


def legend(ax, L, **kw):
    ax.legend(prop={"family": L["font"].get("fontfamily", "DejaVu Sans"),
                    "size": kw.pop("size", 8.5)}, **kw)


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    fig.savefig(os.path.join(L["dir"], name), dpi=150, facecolor="white")
    plt.close(fig)


# ── 그림2: 근축 근사의 오차 ─────────────────────────────────
def figure2(L):
    t = np.deg2rad(np.linspace(0.5, 60, 400))
    err_sin = (t - np.sin(t)) / np.sin(t) * 100
    exact = np.arcsin(np.sin(t) / 1.5)
    err_snell = (t / 1.5 - exact) / exact * 100
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    deg = np.rad2deg(t)
    ax.plot(deg, err_sin, color=C_RAY, label=L["f2sin"])
    ax.plot(deg, np.abs(err_snell), color=C_PAR, label=L["f2snell"])
    for lv in (1, 10):
        ax.axhline(lv, color=C_GREY, lw=0.6, ls=":")
    ax.axvline(5.74, color=C_KEY, lw=1.0, ls="--")
    ax.text(6.5, 12, L["f2mar"], color=C_KEY, fontsize=8.5, **L["font"])
    ax.set_yscale("log"); ax.set_ylim(0.01, 30); ax.set_xlim(0, 60)
    ax.set_xlabel(L["f2x"], **L["font"]); ax.set_ylabel(L["f2y"], **L["font"])
    ax.set_title(L["f2t"], **L["font"])
    legend(ax, L, loc="lower right")
    fig.tight_layout()
    save(fig, L, "fig2-paraxial-error.png")


# ── 그림3: 위상공간에서 본 광선 다발 ─────────────────────────
def polygon_area(pts):
    x, y = pts[:, 0], pts[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def figure3(L):
    f = 50.0
    cases = (("p_a", prop(50.0)), ("p_b", thin_lens(1.0 / f)),
             ("p_c", system(prop(f), thin_lens(1.0 / f), prop(f))))
    # 높이 ±1 mm, 각 ±20 mrad 의 광선 다발 (공기 중이라 감소각 = 각)
    side = np.linspace(-1, 1, 41)
    edge = np.r_[np.c_[side, -np.ones_like(side)], np.c_[np.ones_like(side), side],
                 np.c_[side[::-1], np.ones_like(side)], np.c_[-np.ones_like(side), side[::-1]]]
    box = edge * np.array([1.0, 0.02])
    marks = np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]) * np.array([1.0, 0.02])
    cols = ["#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.0))
    for ax, (key, m) in zip(axes, cases):
        out = box @ m.T
        mk = marks @ m.T
        ax.fill(box[:, 0], box[:, 1] * 1e3, color="#dddddd", label=L["p_before"])
        ax.fill(out[:, 0], out[:, 1] * 1e3, color=C_RAY, alpha=0.35, label=L["p_after"])
        for (a, b), (c, d), col in zip(marks, mk, cols):
            ax.plot(a, b * 1e3, "o", color=col, ms=5, mfc="white")
            ax.plot(c, d * 1e3, "o", color=col, ms=5)
        a0, a1 = polygon_area(box) * 1e3, polygon_area(out) * 1e3
        ax.text(0.03, 0.04, f"{L['p_area']} {a0:.1f} → {a1:.1f} mm·mrad", transform=ax.transAxes,
                fontsize=8.5, **L["font"])
        ax.set_xlim(-2.3, 2.3); ax.set_ylim(-45, 45)
        ax.axhline(0, color=C_GREY, lw=0.4); ax.axvline(0, color=C_GREY, lw=0.4)
        ax.set_xlabel(L["p_x"], **L["font"])
        ax.set_title(L[key], fontsize=10, **L["font"])
    axes[0].set_ylabel(L["p_y"], **L["font"])
    legend(axes[0], L, loc="upper left", size=8)
    fig.tight_layout()
    save(fig, L, "fig3-phase-space.png")


# ── 그림4: 1편 렌즈의 주요점 ────────────────────────────────
def lens_matrix(lens):
    n, na = lens.n_glass, lens.n_air
    return system(surface(na, n, lens.r1), prop(lens.thickness, n), surface(n, na, lens.r2))


def real_ray_path(lens, h, z_start=-30.0):
    """높이 h 의 평행 광선이 렌즈 앞면·뒷면에서 꺾이는 점과 나가는 방향."""
    p = np.array([[z_start, h]])
    d = np.array([[1.0, 0.0]])
    p1, n1 = hit_sphere(p, d, 0.0, lens.r1)
    d1 = refract(d, n1, lens.n_air, lens.n_glass)
    p2, n2 = hit_sphere(p1, d1, lens.thickness, lens.r2)
    d2 = refract(d1, n2, lens.n_glass, lens.n_air)
    return p[0], p1[0], p2[0], d2[0]


def lens_outline(lens, h=100.0):
    y = np.linspace(-h, h, 200)
    z1 = lens.r1 - np.sqrt(lens.r1**2 - y**2)
    z2 = lens.thickness + lens.r2 + np.sqrt(lens.r2**2 - y**2)
    return np.r_[z1, z2[::-1], z1[:1]], np.r_[y, y[::-1], y[:1]]


def figure4(L):
    lens = Singlet()
    c = cardinal_points(lens_matrix(lens), lens.n_air, lens.n_air)
    zF, zH = c.z_front_focus, c.z_front_principal
    zHp, zFp = lens.thickness + c.z_back_principal, lens.thickness + c.z_back_focus
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.5, 4.3), gridspec_kw={"width_ratios": [1.5, 1]})

    zo, yo = lens_outline(lens)
    a1.fill(zo, yo, color="#cfe3f3", ec="k", lw=0.8)
    # (a) 전체: 평행광 → F', F 에서 나온 광선 → 평행
    for h in (60.0, -60.0):
        s, p1, p2, d2 = real_ray_path(lens, h, z_start=-1050)
        t = (zFp - p2[0]) / d2[0]
        a1.plot([s[0], p1[0], p2[0], p2[0] + t * d2[0]], [s[1], p1[1], p2[1], p2[1] + t * d2[1]],
                color=C_RAY, lw=0.9)
        # 대칭 렌즈이므로 앞 초점에서 나온 광선은 위 광선을 거울상으로 뒤집은 것이다
        a1.plot([lens.thickness - (p2[0] + t * d2[0]), lens.thickness - p2[0],
                 lens.thickness - p1[0], 1150],
                [p2[1] + t * d2[1], p2[1], p1[1], s[1]], color=C_PAR, lw=0.9)
    for z, lab, dx, ha in ((zF, "F", 0, "center"), (zH, "H", -8, "right"),
                           (zHp, "H′", 8, "left"), (zFp, "F′", 0, "center")):
        a1.axvline(z, color=C_GREY, lw=0.6, ls="--")
        a1.text(z + dx, 102, lab, ha=ha, fontsize=10)
    a1.set_xlim(-1050, 1150); a1.set_ylim(-110, 110)
    a1.set_xlabel(L["f3x"], **L["font"]); a1.set_ylabel(L["f3y"], **L["font"])
    a1.set_title(L["f3a"], fontsize=10, **L["font"])
    # (b) 실제 광선의 연장선이 만나는 위치가 근축 주평면에서 얼마나 벗어나는가
    hs = np.linspace(1, 100, 80)
    dev = []
    for h in hs:
        _, _, p2, d2 = real_ray_path(lens, h)
        dev.append(p2[0] + (h - p2[1]) * d2[0] / d2[1] - zHp)
    a2.plot(dev, hs, color=C_KEY, lw=2.0, label=L["f3meet"])
    a2.axvline(0, color=C_PAR, lw=1.2, ls="--", label=L["f3par"])
    a2.set_xlim(-1.4, 0.2); a2.set_ylim(0, 105)
    a2.set_xlabel(L["f3bx"], **L["font"]); a2.set_ylabel(L["f3by"], **L["font"])
    a2.set_title(L["f3b"], fontsize=10, **L["font"])
    legend(a2, L, loc="lower left", size=8)
    fig.tight_layout()
    save(fig, L, "fig4-cardinal-points.png")


# ── 그림6: GRIN 막대 렌즈 ──────────────────────────────────
N0, ALPHA = 1.5, 0.5


def figure6(L):
    medium = n2_parabolic(N0, ALPHA)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.0, 4.0), gridspec_kw={"width_ratios": [1.6, 1]})
    z_end = 2 * np.pi / ALPHA
    zz = np.linspace(0, z_end, 400)
    colors = plt.cm.viridis(np.linspace(0, 0.85, 4))
    for deg, col in zip((5, 15, 30, 45), colors):
        th = np.deg2rad(deg)
        sol = trace_ray(medium, (0, 0), th, 3 * z_end, stop=lambda x, y: x - z_end, max_step=0.02)
        _, x, y, *_ = sample(sol, 1500)
        a1.plot(x, y, color=col, lw=1.4, label=f"{deg}°")
        par = np.array([grin(N0, ALPHA, z) @ np.array([0.0, N0 * th]) for z in zz])
        a1.plot(zz, par[:, 0], color=col, lw=1.0, ls="--")
    a1.axvline(np.pi / ALPHA, color=C_GREY, lw=0.6, ls=":")
    a1.set_xlabel(L["f5x"], **L["font"]); a1.set_ylabel(L["f5y"], **L["font"])
    a1.set_title(L["f5a"], fontsize=10, **L["font"])
    legend(a1, L, loc="upper right", size=8)
    degs = np.linspace(1, 60, 24)
    cross = []
    for deg in degs:
        th = np.deg2rad(deg)
        sol = trace_ray(medium, (0, 0), th, 40, stop=lambda x, y: y if x > 0.5 else 1.0)
        cross.append(sol.y[0][-1])
    a2.plot(degs, cross, "o", color=C_RAY, ms=4, label=L["f5exact"])
    dd = np.linspace(0, 60, 200)
    a2.plot(dd, np.pi / ALPHA * np.cos(np.deg2rad(dd)), color=C_KEY, lw=1.2, label=L["f5cos"])
    a2.axhline(np.pi / ALPHA, color=C_PAR, ls="--", lw=1.0, label=L["f5par"])
    a2.set_xlabel(L["f5bx"], **L["font"]); a2.set_ylabel(L["f5by"], **L["font"])
    a2.set_title(L["f5b"], fontsize=10, **L["font"])
    legend(a2, L, loc="lower left", size=8)
    fig.tight_layout()
    save(fig, L, "fig6-grin-rod.png")


# ── 그림7: 주기 렌즈 도파로 ─────────────────────────────────
def figure7(L):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.0, 4.0), gridspec_kw={"width_ratios": [1, 1.5]})
    r = np.linspace(0, 5, 300)
    a1.plot(r, 1 - r / 2, color=C_RAY)
    a1.axhspan(-1, 1, color="#d8eed8", label=L["f6stable"])
    a1.axvline(4, color=C_GREY, ls=":", lw=0.8)
    a1.set_xlabel(L["f6ax"], **L["font"]); a1.set_ylabel(L["f6ay"], **L["font"])
    a1.set_title(L["f6a"], fontsize=10, **L["font"])
    a1.set_ylim(-1.6, 1.2)
    legend(a1, L, loc="lower left", size=8)
    cases = ((1.0, "f6d1", C_RAY, "o"), (2.0, "f6d2", "#2ca02c", "s"),
             (3.5, "f6d35", C_KEY, "^"), (4.2, "f6d42", C_PAR, "x"))
    for d, key, col, mk in cases:
        tr = iterate(periodic_cell(1.0, d), [1.0, 0.0], 24)[:, 0]
        m = np.arange(25)
        keep = np.cumsum(np.abs(tr) > 6) == 0      # 화면을 벗어나기 전까지만
        last = min(keep.sum() + 1, len(tr))
        a2.plot(m[:last], tr[:last], color=col, marker=mk, ms=4, lw=0.8, label=L[key])
    a2.set_ylim(-6.5, 6.5)
    a2.axhline(0, color=C_GREY, lw=0.5)
    a2.set_xlabel(L["f6bx"], **L["font"]); a2.set_ylabel(L["f6by"], **L["font"])
    a2.set_title(L["f6b"], fontsize=10, **L["font"])
    legend(a2, L, loc="upper left", size=8, ncol=2)
    fig.tight_layout()
    save(fig, L, "fig7-periodic-guide.png")


def main():
    for L in LABELS.values():
        figure2(L)
        figure3(L)
        figure4(L)
        figure6(L)
        figure7(L)
    lens = Singlet()
    c = cardinal_points(lens_matrix(lens), lens.n_air, lens.n_air)
    print(f"EFL {c.f_obj:.4f} mm, F {c.z_front_focus:.4f}, H {c.z_front_principal:.4f}, "
          f"H' {c.z_back_principal:.4f}, F' {c.z_back_focus:.4f} (각 꼭짓점 기준)")
    for d in (1.0, 2.0, 3.5, 4.2):
        b, phi = stability(periodic_cell(1.0, d))
        print(f"d = {d}f: b = {b:+.3f}, φ/2π = {phi / (2 * np.pi) if phi == phi else float('nan'):.4f}")


if __name__ == "__main__":
    main()

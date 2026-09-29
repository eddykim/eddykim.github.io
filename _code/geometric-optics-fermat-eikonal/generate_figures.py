"""기하광학 배경이론 1편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/geometric-optics-fermat-eikonal/      (한국어)
      ../../assets/img/posts/geometric-optics-fermat-eikonal/en/   (영문)

그림 1·4·7은 위키미디어 커먼즈 사진(ext-*.jpg)이라 여기서 만들지 않는다.
두 언어의 그림은 같은 계산 결과를 쓰고 표기만 다르다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from rays import (Singlet, envelope_of_lines, least_confusion, n2_linear, opl_polyline, ring_reflection,
                  road_air, air_temperature, sample, trace_ray)

KFONT = {"fontfamily": "AppleGothic"}
EFONT: dict = {}
plt.rcParams["axes.unicode_minus"] = False

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "geometric-optics-fermat-eikonal",
)

C_RAY = "#1f77b4"
C_WAVE = "#d62728"
C_HIT = "#7f7f7f"
C_KEY = "#ff7f0e"

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT,
        # 그림2
        "f2t": "굴절률이 위로 갈수록 커지는 매질 속의 파면과 광선",
        "f2ray": "광선", "f2wave": "파면 (광로길이가 같은 점)",
        "f2n": "굴절률 n",
        # 그림3
        "f3a": "(a) 도로 위 공기의 기온과 굴절률",
        "f3ax": r"기온 ($^\circ$C)", "f3ay": "지면 위 높이 (m)", "f3n": r"$n-1$ ($\times 10^{-4}$)",
        "f3b": "(b) 눈에서 거꾸로 추적한 광선 (세로 축척 과장)",
        "f3bx": "눈으로부터의 거리 (m)", "f3by": "높이 (m)",
        "f3up": "되올라가는 광선: 하늘이 보인다", "f3down": "지면에 닿는 광선: 도로가 보인다",
        "f3eye": "눈", "f3app": "겉보기 방향", "f3edge": "웅덩이가 시작되는 거리",
        # 그림5
        "f5a": "(a) 굴절: 광로길이가 최소인 경로",
        "f5ax": "경계면 위 통과점 x", "f5ay": "광로길이",
        "f5snell": "스넬 법칙이 주는 점",
        "f5b": "(b) 두 초점을 잇는 반사 경로",
        "f5bx": "반사점의 가로 위치", "f5by": "광로길이 S→C→P",
        "f5plane": "평면거울 (최소)", "f5ell": "타원거울 (일정)", "f5circ": "더 굽은 거울 (최대)",
        "f5geo": "(c) 세 거울의 배치",
        # 그림6
        "f6a": "(a) 거울 위 각 점을 거쳐 오는 파의 위상 (실수부)",
        "f6ax": "거울 위 위치 x (mm)", "f6ay": "cos(k·광로길이)",
        "f6b": "(b) 위상자를 차례로 더한 궤적",
        "f6bx": "실수부", "f6by": "허수부",
        "f6c": "정상점 근처",
        "f6start": "시작", "f6end": "끝",
        # 그림8
        "f8a": "(a) 원형 거울 안쪽의 반사 화선",
        "f8b": "(b) 양볼록 렌즈의 굴절 화선",
        "f8bx": "광축 위치 z (mm)", "f8by": "높이 y (mm)",
        "f8par": "근축 초점", "f8mar": "가장자리 광선 초점",
        "f8env": "포락선 (계산)", "f8lc": "최소 착란원", "f8neph": "네프로이드 (닫힌 식)",
        "f8c": "(c) z = 1053 mm 단면의 광선 밀도",
        "f8cx": "높이 y (mm)", "f8cy": "광선 수 (상대값)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT,
        "f2t": "Wavefronts and rays in a medium whose index increases upward",
        "f2ray": "Rays", "f2wave": "Wavefronts (equal optical path)",
        "f2n": "Refractive index n",
        "f3a": "(a) Air temperature and index above a road",
        "f3ax": r"Temperature ($^\circ$C)", "f3ay": "Height above road (m)", "f3n": r"$n-1$ ($\times 10^{-4}$)",
        "f3b": "(b) Rays traced back from the eye (vertical scale exaggerated)",
        "f3bx": "Distance from the eye (m)", "f3by": "Height (m)",
        "f3up": "Rays that turn back up: you see the sky", "f3down": "Rays that hit the road: you see the road",
        "f3eye": "Eye", "f3app": "Apparent direction", "f3edge": "Where the puddle begins",
        "f5a": "(a) Refraction: the path of least optical length",
        "f5ax": "Crossing point x on the interface", "f5ay": "Optical path length",
        "f5snell": "Point given by Snell's law",
        "f5b": "(b) Reflected paths joining two foci",
        "f5bx": "Horizontal position of the reflection point", "f5by": "Optical path S→C→P",
        "f5plane": "Plane mirror (minimum)", "f5ell": "Elliptical mirror (constant)",
        "f5circ": "More curved mirror (maximum)",
        "f5geo": "(c) The three mirrors",
        "f6a": "(a) Phase of the wave arriving via each mirror point (real part)",
        "f6ax": "Position on the mirror x (mm)", "f6ay": "cos(k · optical path)",
        "f6b": "(b) Running sum of the phasors",
        "f6bx": "Real part", "f6by": "Imaginary part",
        "f6c": "Near the stationary point",
        "f6start": "start", "f6end": "end",
        "f8a": "(a) Reflection caustic inside a circular mirror",
        "f8b": "(b) Refraction caustic of a biconvex lens",
        "f8bx": "Axial position z (mm)", "f8by": "Height y (mm)",
        "f8par": "Paraxial focus", "f8mar": "Marginal focus",
        "f8env": "Envelope (computed)", "f8lc": "Circle of least confusion", "f8neph": "Nephroid (closed form)",
        "f8c": "(c) Ray density across z = 1053 mm",
        "f8cx": "Height y (mm)", "f8cy": "Ray count (relative)",
    },
}


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    fig.savefig(os.path.join(L["dir"], name), dpi=150, facecolor="white")
    plt.close(fig)


# ── 그림2: 굴절률 분포 매질 속 파면과 광선 ─────────────────
def compute_wavefronts():
    medium = n2_linear(1.0, 1.5)           # n^2 = 1 + 1.5 y, y=0 에서 1, y=1 에서 1.58
    src = (0.0, 0.25)
    angles = np.deg2rad(np.linspace(-40, 85, 26))
    inside = lambda x, y: min(y, 1.0 - y, 1.6 - x, x + 0.1)
    rays = []
    for a in angles:
        sol = trace_ray(medium, src, a, 4.0, stop=inside, max_step=0.01)
        rays.append(sample(sol, 1500))
    levels = np.arange(0.1, 2.0, 0.1)
    fronts = []
    for lv in levels:
        pts = []
        for s, x, y, px, py, L in rays:
            if L[-1] >= lv:
                pts.append((np.interp(lv, L, x), np.interp(lv, L, y)))
            else:
                pts.append((np.nan, np.nan))
        fronts.append(np.array(pts))
    return medium, src, rays, fronts


def figure2(L, data):
    medium, src, rays, fronts = data
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    yy, xx = np.mgrid[0:1:200j, -0.1:1.6:300j]
    im = ax.imshow(medium.n(xx, yy), origin="lower", extent=(-0.1, 1.6, 0, 1),
                   cmap="Greys", alpha=0.35, aspect="equal", vmin=0.9, vmax=1.7)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label(L["f2n"], **L["font"])
    for i, (s, x, y, *_rest) in enumerate(rays):
        ax.plot(x, y, color=C_RAY, lw=0.9, label=L["f2ray"] if i == 0 else None)
    for i, f in enumerate(fronts):
        ax.plot(f[:, 0], f[:, 1], color=C_WAVE, lw=1.1, label=L["f2wave"] if i == 0 else None)
    ax.plot(*src, "o", color="k", ms=5)
    ax.set_xlim(-0.1, 1.6)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(L["f2t"], **L["font"])
    ax.legend(loc="lower right", prop={"family": L["font"].get("fontfamily", "DejaVu Sans")},
              framealpha=0.9)
    fig.tight_layout()
    save(fig, L, "fig2-wavefronts-rays.png")


# ── 그림3: 뜨거운 도로 위 신기루 ──────────────────────────
EYE = (0.0, 1.5)


def compute_mirage():
    medium = road_air()
    deps = np.deg2rad(np.array([0.05, 0.10, 0.15, 0.20, 0.25, 0.28, 0.30, 0.31, 0.32,
                                0.34, 0.37, 0.42, 0.50, 0.65]))
    rays = []
    for d in deps:
        sol = trace_ray(medium, EYE, -d, 700.0, stop=lambda x, y: min(y, 2.0 - y, 650 - x),
                        max_step=0.5)
        s, x, y, *_ = sample(sol, 3000)
        rays.append((d, x, y, y[-1] <= 1e-6))
    return medium, rays


def figure3(L, data):
    medium, rays = data
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10.5, 4.0),
                                 gridspec_kw={"width_ratios": [1, 3.2]})
    y = np.linspace(0, 1.6, 400)
    a1.plot(air_temperature(y) - 273.15, y, color=C_WAVE)
    a1.set_xlabel(L["f3ax"], color=C_WAVE, **L["font"])
    a1.set_ylabel(L["f3ay"], **L["font"])
    a1.set_title(L["f3a"], fontsize=10, **L["font"])
    a1b = a1.twiny()
    a1b.plot((medium.n(0, y) - 1) * 1e4, y, color=C_RAY, ls="--")
    a1b.set_xlabel(L["f3n"], color=C_RAY, **L["font"])
    a1.set_ylim(0, 1.6)

    up_done = down_done = False
    edge = None
    for d, x, yy, hit in rays:
        if hit:
            lab = None if down_done else L["f3down"]
            down_done = True
            a2.plot(x, yy, color=C_HIT, lw=0.9, label=lab)
        else:
            lab = None if up_done else L["f3up"]
            up_done = True
            a2.plot(x, yy, color=C_RAY, lw=1.1, label=lab)
    # 임계각 바로 아래 광선이 지면을 스치는 거리 = 웅덩이가 시작되는 거리
    theta_c = np.arccos(medium.n(0, 0.0) / medium.n(0, EYE[1]))
    sol = trace_ray(medium, EYE, -(theta_c - 1e-5), 700.0, stop=lambda x, y: y, max_step=0.2)
    _, xc, yc, *_ = sample(sol, 20000)
    edge = xc[np.argmin(yc)]
    ups = [(d, x, yy) for d, x, yy, hit in rays if not hit]
    d_last = max(r[0] for r in ups)
    # 그 광선의 처음 방향을 그대로 연장하면 지면 아래로 들어간다 = 겉보기 방향
    xs = np.linspace(0, 520, 50)
    a2.plot(xs, EYE[1] - np.tan(d_last) * xs, color=C_KEY, ls=":", lw=1.3, label=L["f3app"])
    a2.axvline(edge, color=C_KEY, lw=0.8, alpha=0.6)
    a2.text(edge + 5, 1.75, f"{L['f3edge']} ≈ {edge:.0f} m", color=C_KEY, fontsize=9, **L["font"])
    a2.fill_between([0, 650], -0.3, 0, color="#555555", alpha=0.35)
    a2.plot(*EYE, "o", color="k")
    a2.text(EYE[0] + 8, EYE[1] + 0.07, L["f3eye"], rotation=0, **L["font"])
    a2.set_xlim(0, 650)
    a2.set_ylim(-0.3, 2.0)
    a2.set_xlabel(L["f3bx"], **L["font"])
    a2.set_ylabel(L["f3by"], **L["font"])
    a2.set_title(L["f3b"], fontsize=10, **L["font"])
    a2.legend(loc="lower left", fontsize=8.5,
              prop={"family": L["font"].get("fontfamily", "DejaVu Sans"), "size": 8.5})
    fig.tight_layout()
    save(fig, L, "fig3-road-mirage-rays.png")
    return edge


# ── 그림5: 페르마 원리 — 최소, 일정, 최대 ──────────────────
ELL_A, ELL_B = 2.0, 1.0
FOCUS = np.sqrt(ELL_A**2 - ELL_B**2)
CIRC_R = 2.0                      # 꼭짓점 (0, 1) 에서 타원(곡률반경 a^2/b = 4)보다 더 굽은 원


def mirror_paths(u):
    """반사점의 가로 위치 u 에 대해 세 거울의 반사점과 광로길이 S→C→P."""
    S, P = np.array([-FOCUS, 0.0]), np.array([FOCUS, 0.0])
    pts = {
        "plane": np.stack([u, np.full_like(u, ELL_B)], -1),
        "ell": np.stack([u, ELL_B * np.sqrt(1 - (u / ELL_A) ** 2)], -1),
        "circ": np.stack([u, (ELL_B - CIRC_R) + np.sqrt(CIRC_R**2 - u**2)], -1),
    }
    opl = {k: np.linalg.norm(v - S, axis=-1) + np.linalg.norm(v - P, axis=-1)
           for k, v in pts.items()}
    return S, P, pts, opl


def figure5(L):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8),
                             gridspec_kw={"width_ratios": [1, 1, 1.6]})
    # (a) 굴절: S=(0,1) 공기, P=(2,-1) 유리
    n1, n2 = 1.0, 1.5
    x = np.linspace(-0.5, 2.5, 400)
    opl = n1 * np.hypot(x, 1.0) + n2 * np.hypot(2.0 - x, 1.0)
    i = np.argmin(opl)
    ax = axes[0]
    ax.plot(x, opl, color=C_RAY)
    ax.plot(x[i], opl[i], "o", color=C_KEY, label=L["f5snell"])
    ax.set_xlabel(L["f5ax"], **L["font"]); ax.set_ylabel(L["f5ay"], **L["font"])
    ax.set_title(L["f5a"], fontsize=10, **L["font"])
    ax.legend(prop={"family": L["font"].get("fontfamily", "DejaVu Sans"), "size": 8.5})

    # (b) 세 거울
    u = np.linspace(-0.9, 0.9, 400)
    S, P, pts, opl = mirror_paths(u)
    ax = axes[1]
    ax.plot(u, opl["plane"], color=C_RAY, label=L["f5plane"])
    ax.plot(u, opl["ell"], color="k", label=L["f5ell"])
    ax.plot(u, opl["circ"], color=C_WAVE, label=L["f5circ"])
    ax.set_xlabel(L["f5bx"], **L["font"]); ax.set_ylabel(L["f5by"], **L["font"])
    ax.set_title(L["f5b"], fontsize=10, **L["font"])
    ax.legend(prop={"family": L["font"].get("fontfamily", "DejaVu Sans"), "size": 8})

    # (c) 배치
    ax = axes[2]
    t = np.linspace(0, np.pi, 300)
    ax.plot(ELL_A * np.cos(t), ELL_B * np.sin(t), color="k", lw=1.2)
    ax.plot([-1.5, 1.5], [ELL_B, ELL_B], color=C_RAY, lw=1.2)
    uu = np.linspace(-1.3, 1.3, 200)
    ax.plot(uu, (ELL_B - CIRC_R) + np.sqrt(CIRC_R**2 - uu**2), color=C_WAVE, lw=1.2)
    for k, col in (("plane", C_RAY), ("ell", "k"), ("circ", C_WAVE)):
        c = pts[k][np.argmin(abs(u - 0.6))]
        ax.plot([S[0], c[0], P[0]], [S[1], c[1], P[1]], color=col, lw=0.7, ls="--")
    ax.plot([S[0], 0, P[0]], [S[1], ELL_B, P[1]], color=C_KEY, lw=1.4)
    ax.plot(*S, "ko"); ax.plot(*P, "ko")
    ax.text(S[0], S[1] - 0.18, "S", ha="center"); ax.text(P[0], P[1] - 0.18, "P", ha="center")
    ax.text(0.03, ELL_B + 0.06, "Q")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlim(-2.15, 2.15); ax.set_ylim(-0.3, 1.2)
    ax.set_title(L["f5geo"], fontsize=10, **L["font"])
    fig.tight_layout()
    save(fig, L, "fig5-fermat-opl.png")


# ── 그림6: 정상 경로 근처의 위상자가 보강된다 ───────────────
PH_A, PH_H, PH_LAMBDA = 10.0, 10.0, 0.5e-3     # mm. 초점 사이 반거리, 높이, 파장


def phasor_sum(x):
    k = 2 * np.pi / PH_LAMBDA
    opl = np.hypot(x + PH_A, PH_H) + np.hypot(x - PH_A, PH_H)
    ph = np.exp(1j * k * (opl - opl.min()))
    csum = np.concatenate([[0], np.cumsum(0.5 * (ph[1:] + ph[:-1]) * np.diff(x))])
    return ph, csum


def figure6(L):
    x = np.linspace(-0.6, 0.6, 60001)
    ph, csum = phasor_sum(x)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.0, 4.0),
                                 gridspec_kw={"width_ratios": [1.5, 1]})
    a1.plot(x, ph.real, color=C_RAY, lw=0.5)
    a1.axvspan(-0.08, 0.08, color=C_KEY, alpha=0.2, label=L["f6c"])
    a1.set_xlabel(L["f6ax"], **L["font"]); a1.set_ylabel(L["f6ay"], **L["font"])
    a1.set_title(L["f6a"], fontsize=10, **L["font"])
    a1.set_xlim(x[0], x[-1])
    a1.legend(loc="lower right",
              prop={"family": L["font"].get("fontfamily", "DejaVu Sans"), "size": 8.5})

    pts = np.stack([csum.real, csum.imag], -1)
    segs = np.stack([pts[:-1], pts[1:]], 1)
    lc = LineCollection(segs, cmap="coolwarm", linewidths=1.2)
    lc.set_array(0.5 * (x[1:] + x[:-1]))
    a2.add_collection(lc)
    near = np.abs(x) < 0.08
    a2.plot(csum.real[near], csum.imag[near], color=C_KEY, lw=3.5, alpha=0.5)
    a2.plot(*pts[0], "ko", ms=4); a2.plot(*pts[-1], "ks", ms=4)
    a2.annotate(L["f6start"], pts[0], xytext=(8, -4), textcoords="offset points", **L["font"])
    a2.annotate(L["f6end"], pts[-1], xytext=(8, -4), textcoords="offset points", **L["font"])
    a2.autoscale(); a2.set_aspect("equal")
    a2.set_xlabel(L["f6bx"], **L["font"]); a2.set_ylabel(L["f6by"], **L["font"])
    a2.set_title(L["f6b"], fontsize=10, **L["font"])
    fig.tight_layout()
    save(fig, L, "fig6-phasor-sum.png")


# ── 그림8: 화선 ─────────────────────────────────────────
Z_CUT = 1053.0


def compute_caustics():
    # (a) 원형 거울, 반경 1
    h = np.linspace(-0.98, 0.98, 41)
    hit, d = ring_reflection(1.0, h)
    hd = np.linspace(-0.999, 0.999, 4001)
    hit_d, d_d = ring_reflection(1.0, hd)
    # x-y 평면이므로 envelope_of_lines 의 (z, y) 자리에 (x, y) 를 넣는다
    ex, ey = envelope_of_lines(hit_d[:, 0], hit_d[:, 1], d_d[:, 1] / d_d[:, 0], hd)
    # (b) 렌즈
    lens = Singlet()
    hl = np.linspace(-100, 100, 81)
    p, dl = lens.trace(hl)
    hd2 = np.linspace(1e-3, 100, 4001)
    p2, d2 = lens.trace(hd2)
    lz, ly = envelope_of_lines(p2[:, 0], p2[:, 1], d2[:, 1] / d2[:, 0], hd2)
    z_par = lens.axis_crossing([1e-4])[0]
    z_mar = lens.axis_crossing([100.0])[0]
    z_lc, d_lc = least_confusion(lens)
    # (c) 단면의 광선 밀도: 입사 동공을 고르게 채운 광선 (2차원 단면이라 높이에 균등)
    hh = np.linspace(-100, 100, 400001)
    p3, d3 = lens.trace(hh)
    y_cut = p3[:, 1] + (Z_CUT - p3[:, 0]) * d3[:, 1] / d3[:, 0]
    return dict(hit=hit, d=d, ex=ex, ey=ey, p=p, dl=dl, lz=lz, ly=ly,
                z_par=z_par, z_mar=z_mar, z_lc=z_lc, d_lc=d_lc, y_cut=y_cut)


def figure8(L, c):
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.2),
                             gridspec_kw={"width_ratios": [1, 1.35, 1]})
    ax = axes[0]
    t = np.linspace(0, 2 * np.pi, 400)
    ax.plot(np.cos(t), np.sin(t), color="k", lw=1.5)
    for (hx, hy), (dx, dy) in zip(c["hit"], c["d"]):
        ax.plot([-1.0, hx], [hy, hy], color=C_RAY, lw=0.35, alpha=0.5)
        ax.plot([hx, hx + 2.2 * dx], [hy, hy + 2.2 * dy], color=C_RAY, lw=0.35, alpha=0.5)
    tt = np.linspace(-np.pi / 2, np.pi / 2, 400)
    neph_x = 0.25 * (3 * np.cos(tt) - np.cos(3 * tt))
    neph_y = 0.25 * (3 * np.sin(tt) - np.sin(3 * tt))
    ok = np.abs(c["ey"]) < 1
    ax.plot(c["ex"][ok], c["ey"][ok], color=C_WAVE, lw=2.5, alpha=0.6, label=L["f8env"])
    ax.plot(neph_x, neph_y, color="k", lw=0.9, ls="--", label=L["f8neph"])
    ax.set_xlim(-1.05, 1.05); ax.set_ylim(-1.05, 1.05); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(L["f8a"], fontsize=10, **L["font"])
    ax.legend(loc="lower left", prop={"family": L["font"].get("fontfamily", "DejaVu Sans"),
                                      "size": 7.5})

    ax = axes[1]
    z_end = 1075.0
    for (pz, py), (dz, dy) in zip(c["p"], c["dl"]):
        ax.plot([pz, z_end], [py, py + (z_end - pz) * dy / dz], color=C_RAY, lw=0.4, alpha=0.6)
    ok = (c["lz"] > 1035) & (c["lz"] < z_end)
    ax.plot(c["lz"][ok], c["ly"][ok], color=C_WAVE, lw=2.2, alpha=0.7, label=L["f8env"])
    ax.plot(c["lz"][ok], -c["ly"][ok], color=C_WAVE, lw=2.2, alpha=0.7)
    ax.axvline(c["z_par"], color="k", ls="--", lw=0.8)
    ax.axvline(c["z_mar"], color=C_KEY, ls="--", lw=0.8)
    ax.axvline(Z_CUT, color=C_HIT, ls=":", lw=0.8)
    ax.plot([c["z_lc"]] * 2, [-c["d_lc"] / 2, c["d_lc"] / 2], color="k", lw=3.0,
            solid_capstyle="butt", label=L["f8lc"])
    ax.text(c["z_par"] + 0.4, 3.2, f"{L['f8par']}\n{c['z_par']:.1f}", fontsize=8, **L["font"])
    ax.text(c["z_mar"] - 0.4, 3.2, f"{L['f8mar']}\n{c['z_mar']:.1f}", fontsize=8,
            ha="right", color=C_KEY, **L["font"])
    ax.set_xlim(1035, z_end); ax.set_ylim(-4, 4)
    ax.set_xlabel(L["f8bx"], **L["font"]); ax.set_ylabel(L["f8by"], **L["font"])
    ax.set_title(L["f8b"], fontsize=10, **L["font"])
    ax.legend(loc="lower left", prop={"family": L["font"].get("fontfamily", "DejaVu Sans"),
                                      "size": 8})

    ax = axes[2]
    cnt, edges = np.histogram(c["y_cut"], bins=241, range=(-2.4, 2.4))
    ax.step(0.5 * (edges[1:] + edges[:-1]), cnt / cnt.max(), color=C_RAY, where="mid")
    ax.set_xlabel(L["f8cx"], **L["font"]); ax.set_ylabel(L["f8cy"], **L["font"])
    ax.set_title(L["f8c"], fontsize=10, **L["font"])
    fig.tight_layout()
    save(fig, L, "fig8-caustics.png")


def main():
    wf = compute_wavefronts()
    mg = compute_mirage()
    cs = compute_caustics()
    for lang, L in LABELS.items():
        figure2(L, wf)
        edge = figure3(L, mg)
        figure5(L)
        figure6(L)
        figure8(L, cs)
    print(f"웅덩이가 시작되는 거리 {edge:.1f} m")
    print(f"근축 초점 {cs['z_par']:.3f} mm, 가장자리 광선 초점 {cs['z_mar']:.3f} mm, "
          f"종방향 구면수차 {cs['z_mar'] - cs['z_par']:.3f} mm")


if __name__ == "__main__":
    main()

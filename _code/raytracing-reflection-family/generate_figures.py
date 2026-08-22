"""3편 그림 2개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-reflection-family/ 에 fig1~2 저장
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from reflection_optics import (
    _vertex_normal,
    make_arbitrary_mirror,
    make_beam_splitter,
    make_flat_mirror,
    make_spherical_mirror,
    reflect_arbitrary_mirror,
    reflect_beam_splitter,
    reflect_flat_mirror,
    reflect_spherical_mirror,
    reflect_vector,
)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "raytracing-reflection-family",
)
os.makedirs(OUT_DIR, exist_ok=True)


def extend(path, length=60.0):
    x, y, vx, vy = path[-1]
    return np.vstack([path, [x + length * vx, y + length * vy, vx, vy]])


def draw_ray(ax, path, **kwargs):
    p = extend(np.atleast_2d(path))
    ax.plot(p[:, 0], p[:, 1], **kwargs)


# ── 그림1: 네 부품의 반사 경로 ──────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(9, 8))

ax = axes[0, 0]
mirror = make_flat_mirror(200, 0, 0, np.deg2rad(20))
ray = np.array([-100.0, 30.0, 1.0, 0.0])
out = reflect_flat_mirror(ray, mirror)
b = mirror["BOUNDARY"][0]
ax.plot(b[:, 0], b[:, 1], "k-", linewidth=3)
draw_ray(ax, np.vstack([ray, out[0]]), color="C3", linewidth=1.2)
ax.set_title("FlatMirror", **KFONT)

ax = axes[0, 1]
mirror = make_spherical_mirror(-300, 200, 0, 0, 0)
ray = np.array([-100.0, 40.0, 1.0, 0.0])
out = reflect_spherical_mirror(ray, mirror)
b = mirror["BOUNDARY"][0]
ax.plot(b[:, 0], b[:, 1], "k-", linewidth=3)
draw_ray(ax, np.vstack([ray, out[0]]), color="C3", linewidth=1.2)
ax.set_title("SphericalMirror", **KFONT)

ax = axes[1, 0]
mirror = make_beam_splitter(60, 0, 0, np.deg2rad(45))
ray = np.array([-100.0, 0.0, 1.0, 0.0])
r, t = reflect_beam_splitter(ray, mirror)
b = mirror["BOUNDARY"][0]
ax.plot(b[:, 0], b[:, 1], "k-", linewidth=3)
draw_ray(ax, np.vstack([ray, r[0]]), color="C3", linewidth=1.2, label="반사")
draw_ray(ax, np.vstack([ray, t[0]]), color="C0", linewidth=1.2, label="투과")
ax.plot(ray[0], ray[1], "o", color="k", markersize=3)
ax.legend(prop=LFONT, fontsize=8, loc="upper left")
ax.set_title("BeamSplitter", **KFONT)

ax = axes[1, 1]
mirror = make_arbitrary_mirror([-100, 0, 100, 200], [0, 50, 0, 50])
ray = np.array([-50.0, 100.0, 0.0, -1.0])
out = reflect_arbitrary_mirror(ray, mirror)
for seg in mirror["BOUNDARY"]:
    ax.plot(seg[:, 0], seg[:, 1], "k-", linewidth=3)
draw_ray(ax, np.vstack([ray, out[1:]]), color="C3", linewidth=1.2)
ax.set_title("ArbitraryMirror", **KFONT)

for ax in axes.flat:
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("x (mm)", **KFONT)
    ax.set_ylabel("y (mm)", **KFONT)

fig.suptitle("그림1. 네 반사 부품의 광선 경로 -- 반사 공식은 하나, 법선 찾는 법만 다르다", **KFONT)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1-reflection-family.png"), dpi=150)
plt.close(fig)

# ── 그림2: 다각형 꼭짓점 법선 -- 점을 재사용한 경우 vs 세 점을 모두 쓴 경우 ──
x0, y0, x1, y1, x2, y2 = -100, 0, 0, 50, 100, 0
n_dup = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, dup_point=True)
n_uniq = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, dup_point=False)
ang = np.rad2deg(np.arccos(np.clip(np.dot(n_dup, n_uniq), -1, 1)))
print(f"대칭 꼭짓점: n_dup={n_dup}, n_uniq={n_uniq}, 차이={ang:.3f}deg")

vx_in, vy_in = 0.0, -1.0
v_dup = reflect_vector([vx_in, vy_in], n_dup)
v_uniq = reflect_vector([vx_in, vy_in], n_uniq)

fig, ax = plt.subplots(figsize=(6, 5))
poly = np.array([[x0, y0], [x1, y1], [x2, y2]])
ax.plot(poly[:, 0], poly[:, 1], "k-", linewidth=3)
ax.plot([0, 0], [150, 50], "-", color="0.5", linewidth=1.2)
ax.annotate("입사광", (5, 120), **KFONT, fontsize=9)
L = 80
ax.plot([x1, x1 + L * v_dup[0]], [y1, y1 + L * v_dup[1]], "-", color="C3", linewidth=1.5, label="세 번째 점을 재사용한 법선")
ax.plot([x1, x1 + L * v_uniq[0]], [y1, y1 + L * v_uniq[1]], "--", color="C0", linewidth=1.5, label="세 점을 모두 쓴 법선")
ax.plot(x1, y1, "o", color="k", markersize=4)
ax.set_xlim(-120, 120)
ax.set_ylim(-20, 160)
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.legend(prop=LFONT, fontsize=9, loc="lower right")
ax.set_title(f"그림2. 세 번째 점을 재사용하면 반사 방향이 {ang:.1f}° 달라진다", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig2-vertex-normal-typo.png"), dpi=150)
plt.close(fig)

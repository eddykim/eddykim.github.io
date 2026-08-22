"""2편 그림 3개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-total-internal-reflection/ 에 fig1~3 저장
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from lens_tir import (
    exit_surface_is_tir,
    find_entry_surface,
)
from lens_tir import load_nk_file as load_nk_lens
from lens_tir import make_parallel_beam, make_spherical_lens
from lens_tir import refract_through_lens as refract_through_lens_fixed
from lens_tir import refractive_index as refractive_index_lens
from prism_optics import load_nk_file as load_nk_prism
from prism_optics import make_prism, refractive_index as refractive_index_prism
from prism_optics import trace_ray_through_prism

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "raytracing-total-internal-reflection",
)
os.makedirs(OUT_DIR, exist_ok=True)

WAVELENGTH_NM = 750.0


def extend_path(path, length=150.0):
    """마지막 점에서 마지막 방향으로 length만큼 더 이어붙인다 (그림에서
    출사 후 광선이 화면 밖으로 빠져나가는 걸 보여주기 위한 시각화 전용 연장)."""
    x, y, vx, vy = path[-1]
    return np.vstack([path, [x + length * vx, y + length * vy, vx, vy]])


def draw_prism_outline(ax, prism):
    verts = prism["verts"]
    outline = np.vstack([verts, verts[0:1]])
    ax.fill(outline[:, 0], outline[:, 1], color="C0", alpha=0.2)


# ============================================================
# 그림1: 정삼각형 프리즘, 수직 입사 광선이 한 면에서 전반사 후
#        다른 면으로 수직 출사
# ============================================================
nk_prism = load_nk_file_prism = load_nk_prism("N-BK7.nk")
n_glass = refractive_index_prism(nk_prism, WAVELENGTH_NM)
print(f"n(N-BK7, {WAVELENGTH_NM:.0f}nm) = {n_glass:.6f}")

prism = make_prism([0, 300 * np.sqrt(3), 0], [300, 0, -300])
ray0 = np.array([-150.0, 100.0, 1.0, 0.0])
path0 = trace_ray_through_prism(ray0, prism, n_glass, 1.0, verbose=True)
path0 = extend_path(path0)

fig, ax = plt.subplots(figsize=(7, 5))
draw_prism_outline(ax, prism)
ax.plot(path0[:, 0], path0[:, 1], "-o", color="C3", linewidth=1.5, markersize=4)
ax.annotate("입사\n(수직, 0°)", (-100, 105), fontsize=9, **KFONT)
ax.annotate("전반사(TIR)\n내부입사각 60°", (346, 108), fontsize=9, **KFONT)
ax.annotate("무굴절 출사\n(내부입사각 0°)", (433, -42), fontsize=9, **KFONT)
ax.set_xlim(-200, 600)
ax.set_ylim(-350, 350)
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.set_title("그림1. 수직 입사 광선이 정삼각형 프리즘 안에서 한 번 전반사한다", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1-prism-tir-path.png"), dpi=150)
plt.close(fig)

# ============================================================
# 그림2: 입사각 스윕 -- 임계각 아래에서는 전반사, 위에서는 투과
# ============================================================
theta_c = np.arcsin(1.0 / n_glass)
A_apex = np.deg2rad(60.0)
theta1_thresh = np.rad2deg(np.arcsin(n_glass * np.sin(A_apex - theta_c)))
print(f"theta_c = {np.rad2deg(theta_c):.4f}deg, theta1 임계값(이론) = {theta1_thresh:.6f}deg")


def fire(theta1_deg, y_hit=100.0, x0=-150.0):
    theta1 = np.deg2rad(theta1_deg)
    vx, vy = np.cos(theta1), np.sin(theta1)
    y0 = y_hit - (0 - x0) * np.tan(theta1)
    ray = np.array([x0, y0, vx, vy])
    return trace_ray_through_prism(ray, prism, n_glass, 1.0)


angles = [0.0, 10.0, 20.0, 28.813, 35.0, 45.0]
fig, ax = plt.subplots(figsize=(7, 5))
draw_prism_outline(ax, prism)
cmap = plt.get_cmap("coolwarm")
for i, ang in enumerate(angles):
    p = extend_path(fire(ang))
    is_tir = len(fire(ang)) >= 4
    color = cmap(i / (len(angles) - 1))
    label = f"입사각={ang:.1f}°  ({'TIR' if is_tir else '투과'})"
    ax.plot(p[:, 0], p[:, 1], "-", color=color, linewidth=1.3, label=label)
ax.set_xlim(-200, 600)
ax.set_ylim(-350, 350)
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.legend(prop=LFONT, loc="lower right", fontsize=8)
ax.set_title(f"그림2. 입사각이 임계값({theta1_thresh:.1f}°)을 넘으면 전반사 대신 투과한다", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig2-prism-angle-sweep.png"), dpi=150)
plt.close(fig)

# ============================================================
# 그림3: 평철렌즈 -- 일부 마진광선은 투과, 일부는 출사면에서 반사
# ============================================================
nk_lens = load_nk_lens("N-BK7.nk")
n_lens = refractive_index_lens(nk_lens, WAVELENGTH_NM)
lens = make_spherical_lens(R1=1e6, R2=-50.0, t=20.0, D=78.0)
source = make_parallel_beam(x=-1000, y=0, vx=1, vy=0, D=76, num=39)

fig, ax = plt.subplots(figsize=(8, 5))
b1, b2 = lens["BOUNDARY"]
outline = np.vstack([b1, b2[::-1], b1[0:1]])
ax.fill(outline[:, 0], outline[:, 1], color="C0", alpha=0.2)

n_tir = 0
for ray in source:
    entry = find_entry_surface(ray, lens)
    if entry is None:
        continue
    surf_idx, _ = entry
    inside = refract_through_lens_fixed(ray, surf_idx, lens, n_lens, 1.0)
    tir = exit_surface_is_tir(inside[0], surf_idx, lens, n_lens, 1.0)
    n_tir += tir
    full = extend_path(np.vstack([ray, inside]))
    ax.plot(full[:, 0], full[:, 1], "-", color=("C3" if tir else "C0"), linewidth=0.8, alpha=0.7)

print(f"광선 {len(source)}개 중 출사면에서 전반사한 광선 수: {n_tir}")

ax.plot([], [], "-", color="C0", label="투과")
ax.plot([], [], "-", color="C3", label="출사면에서 전반사(TIR)")
ax.set_xlim(-70, 70)
ax.set_ylim(-45, 45)
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.legend(prop=LFONT, loc="upper left", fontsize=9)
ax.set_title("그림3. 평철렌즈(R2=-50mm) 마진광선 중 가장자리는 출사면에서 전반사한다", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig3-lens-marginal-tir.png"), dpi=150)
plt.close(fig)

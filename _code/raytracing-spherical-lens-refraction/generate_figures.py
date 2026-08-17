"""1편 그림 3개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-spherical-lens-refraction/ 에 fig1~3 저장
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from sphere_optics import (
    find_entry_surface,
    load_nk_file,
    make_parallel_beam,
    make_spherical_lens,
    refract_through_lens,
    refractive_index,
    trace_parallel_beam_through_lens,
)

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "raytracing-spherical-lens-refraction",
)
os.makedirs(OUT_DIR, exist_ok=True)

# 공통 렌즈/광원/환경 설정
R1, R2, T, D = 1000.0, -1000.0, 100.0, 500.0
WAVELENGTH_NM = 750.0
N_AIR = 1.0002778
ROI = (-1500, 1500, -500, 500)

lens = make_spherical_lens(R1=R1, R2=R2, t=T, D=D)
nk = load_nk_file(os.path.join(os.path.dirname(__file__), "N-BK7.nk"))
n_lens = refractive_index(nk, WAVELENGTH_NM)
source = make_parallel_beam(x=1000, y=0, vx=-1, vy=0, D=200, num=11)
paths = trace_parallel_beam_through_lens(lens, source, WAVELENGTH_NM, N_AIR, nk, ROI)

print(f"n(N-BK7, {WAVELENGTH_NM:.0f}nm) = {n_lens:.6f}")


def draw_lens_outline(ax, lens):
    b1, b2 = lens["BOUNDARY"]
    outline = np.vstack([b1, b2[::-1], b1[0:1]])
    ax.fill(outline[:, 0], outline[:, 1], color="C0", alpha=0.25)


def draw_paths(ax, paths, **kwargs):
    for p in paths:
        ax.plot(p[:, 0], p[:, 1], **kwargs)


# ── 그림1: 평행광이 양볼록 렌즈에 입사해 굴절되는 모습 ───────────────────
fig, ax = plt.subplots(figsize=(7, 4))
draw_lens_outline(ax, lens)
draw_paths(ax, paths, color="C3", linewidth=0.8)
ax.set_xlim(ROI[0], ROI[1])
ax.set_ylim(ROI[2], ROI[3])
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.set_title("그림1. 평행광선 11개가 양볼록 렌즈를 지나 한 점으로 모인다", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1-lens-focusing.png"), dpi=150)
plt.close(fig)

# ── 그림2: 광학 가역성(reversibility) 검증 ───────────────────────────────
# R1 쪽으로 들어온 광선을 굴절시킨 뒤(정방향, entry_surface=1), 그 출사
# 광선을 반대로 뒤집어 다시 쏘면(역방향, entry_surface=2 분기를 타게 됨)
# 원래 입사 경로를 정확히 되짚어 나와야 한다 -- Snell's law 자체가 시간
# 역전에 대해 대칭이기 때문이다. 이건 refract_through_lens의 두 분기
# (entry_surface==1 / else)가 서로 모순 없이 같은 물리를 구현했는지 보는
# 좋은 교차검증이 된다.
ray_fwd = np.array([1000.0, 100.0, -1.0, 0.0])
entry_fwd = find_entry_surface(ray_fwd, lens)[0]
path_fwd = refract_through_lens(ray_fwd, entry_fwd, lens, n_lens, N_AIR)
full_fwd = np.vstack([ray_fwd, path_fwd])

# 출사 광선을 렌즈 밖으로 500mm 더 진행시킨 지점에서, 방향만 반대로 뒤집어
# 새 광선을 만든다.
x2, y2, vx2, vy2 = path_fwd[-1]
far_point = np.array([x2 + 500 * vx2, y2 + 500 * vy2])
ray_retrace = np.array([far_point[0], far_point[1], -vx2, -vy2])
entry_retrace = find_entry_surface(ray_retrace, lens)[0]
path_retrace = refract_through_lens(ray_retrace, entry_retrace, lens, n_lens, N_AIR)
# 원래 입사 광선과 같은 x=1000까지 뻗어서, 되짚은 경로가 원래 경로와
# 겹치는지 눈으로 바로 비교할 수 있게 한다.
xr, yr, vxr, vyr = path_retrace[-1]
t_tail = (1000 - xr) / vxr
retrace_tail_xy = np.array([1000.0, yr + t_tail * vyr])
full_retrace_xy = np.vstack([ray_retrace[:2], path_retrace[:, :2], retrace_tail_xy])

print(f"정방향 진입면={entry_fwd}, 역추적 진입면={entry_retrace} (반대쪽 면이어야 함)")
print(f"원래 입사 광선 시작점: ({ray_fwd[0]:.3f}, {ray_fwd[1]:.3f}), 방향 ({ray_fwd[2]:.3f}, {ray_fwd[3]:.3f})")
print(f"역추적 도착점(렌즈 첫 굴절점): ({path_retrace[-1,0]:.3f}, {path_retrace[-1,1]:.3f}), 방향 ({path_retrace[-1,2]:.3f}, {path_retrace[-1,3]:.3f})")
pos_err = np.hypot(path_retrace[-1, 0] - path_fwd[0, 0], path_retrace[-1, 1] - path_fwd[0, 1])
print(f"위치 오차: {pos_err:.2e} mm")

fig, ax = plt.subplots(figsize=(7, 4))
draw_lens_outline(ax, lens)
ax.plot(full_fwd[:, 0], full_fwd[:, 1], "-", color="C0", linewidth=3, alpha=0.6, label=f"정방향 ({entry_fwd}번 면 먼저 진입)")
ax.plot(full_retrace_xy[:, 0], full_retrace_xy[:, 1], "--", color="C3", linewidth=1.5, label="역추적 (출사 광선을 반대로)")
ax.set_xlim(-1200, 1600)
ax.set_ylim(-50, 150)
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.legend(prop=LFONT, loc="lower right")
ax.set_title("그림2. 출사 광선을 반대로 쏘면 입사 경로를 그대로 되짚는다", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig2-reversibility.png"), dpi=150)
plt.close(fig)

# ── 그림3: 시뮬레이션 초점 vs thick-lens 이론값 ──────────────────────────
n = n_lens
inv_f = (n - 1) * (1 / R1 - 1 / R2 + (n - 1) * T / (n * R1 * R2))
f = 1 / inv_f
bfd = f * (1 - (n - 1) * T / (n * R1))
rear_vertex_x = -T / 2
expected_focus_x = rear_vertex_x - bfd

# 렌즈를 나온 직후 위치(x2,y2)와 방향(vx2,vy2)을 알고 있으니, 그 직선이
# 광축(y=0)과 만나는 x좌표를 각 광선마다 구한다 (렌즈 출사점 자체는 초점이
# 아니라 렌즈 두께 절반(x=-T/2) 근방일 뿐이라 이 계산이 따로 필요하다).
# 근축(paraxial) 근사와 비교하려면 광축에 가장 가까운 광선들을 봐야 한다 -
# 가장자리 광선은 구면수차 때문에 초점 위치 자체가 다르다(그림3 자체가 이걸
# 보여준다).
exit_rows = [p[-2] for p in paths]  # [x2, y2, vx2, vy2], 광원 y좌표 오름차순
start_ys = np.array([s[1] for s in source])
def axis_crossing(y0_limit_lo, y0_limit_hi):
    xs = [
        (x2 - y2 * vx2 / vy2)
        for (x2, y2, vx2, vy2), y0 in zip(exit_rows, start_ys)
        if y0_limit_lo <= abs(y0) <= y0_limit_hi and abs(vy2) > 1e-9
    ]
    return float(np.median(xs))


sim_focus_x = axis_crossing(0, 20)  # 근축: |y0| <= 20mm (렌즈 반지름 250mm의 8%)
marginal_focus_x = axis_crossing(90, 100)  # 최외곽: |y0| = 100mm (광원 조리개 경계)

print(f"thick-lens 이론 초점(근축): f={f:.3f}mm, 후방초점거리={bfd:.3f}mm, 초점 x={expected_focus_x:.3f}")
print(f"시뮬레이션 근축 초점(|y0|<=20mm) x={sim_focus_x:.3f}, 이론 대비 {sim_focus_x - expected_focus_x:+.3f}mm ({100*(sim_focus_x-expected_focus_x)/abs(expected_focus_x):+.3f}%)")
print(f"시뮬레이션 최외곽 초점(|y0|=100mm) x={marginal_focus_x:.3f}, 근축 대비 {marginal_focus_x - sim_focus_x:+.3f}mm (구면수차)")

fig, ax = plt.subplots(figsize=(7, 4))
draw_lens_outline(ax, lens)
draw_paths(ax, paths, color="C3", linewidth=0.8, alpha=0.6)
ax.axvline(expected_focus_x, color="k", linestyle="--", linewidth=1, label=f"thick-lens 이론값 (x={expected_focus_x:.1f}mm)")
ax.axvline(sim_focus_x, color="C0", linestyle=":", linewidth=1.5, label=f"근축 광선 초점 (x={sim_focus_x:.1f}mm)")
ax.axvline(marginal_focus_x, color="C3", linestyle=":", linewidth=1.5, label=f"최외곽 광선 초점 (x={marginal_focus_x:.1f}mm)")
ax.set_xlim(-1300, -700)
ax.set_ylim(-150, 150)
ax.set_aspect("equal")
ax.set_xlabel("x (mm)", **KFONT)
ax.set_ylabel("y (mm)", **KFONT)
ax.legend(prop=LFONT, loc="upper right")
ax.set_title("그림3. 근축 이론값과 실제(유한 조리개) 광선추적 초점의 차이", **KFONT)
ax.grid(True, alpha=0.3)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig3-focus-vs-theory.png"), dpi=150)
plt.close(fig)

"""기하광학 1편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-spherical-lens-refraction/      (한국어)
      ../../assets/img/posts/raytracing-spherical-lens-refraction/en/   (영문)

광선추적 계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의
그림은 데이터가 완전히 동일하고 표기만 다르다.
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

# 한글 텍스트(제목/축이름/범례)에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면
# 눈금의 마이너스 기호가 한글 폰트에 없어 깨지기 때문이다.
# macOS 전용 설정이다. 다른 OS는 나눔고딕 등 설치된 한글 폰트 이름으로 교체할 것.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}  # legend(prop=...)는 FontProperties 인자 이름(family)을 쓴다

# 영문판은 한글 폰트가 필요 없으므로 matplotlib 기본 폰트를 그대로 쓴다.
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "raytracing-spherical-lens-refraction",
)

# 언어별 라벨. 계산 코드는 공유하고 문자열만 갈아 끼운다.
LABELS = {
    "ko": {
        "dir": BASE_DIR,
        "font": KFONT, "legend": LFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "fig1": "그림1. 평행광선 11개가 양볼록 렌즈를 지나 한 점으로 모인다",
        "forward": "정방향 ({}번 면 먼저 진입)",
        "retrace": "역추적 (출사 광선을 반대로)",
        "fig2": "그림2. 출사 광선을 반대로 쏘면 입사 경로를 그대로 되짚는다",
        "theory": "thick-lens 이론값 (x={:.1f}mm)",
        "paraxial": "근축 광선 초점 (x={:.1f}mm)",
        "marginal": "최외곽 광선 초점 (x={:.1f}mm)",
        "fig3": "그림3. 근축 이론값과 실제(유한 조리개) 광선추적 초점의 차이",
        "p_index": "n(N-BK7, {:.0f}nm) = {:.6f}",
        "p_entry": "정방향 진입면={}, 역추적 진입면={} (반대쪽 면이어야 함)",
        "p_start": "원래 입사 광선 시작점: ({:.3f}, {:.3f}), 방향 ({:.3f}, {:.3f})",
        "p_end": "역추적 도착점(렌즈 첫 굴절점): ({:.3f}, {:.3f}), 방향 ({:.3f}, {:.3f})",
        "p_err": "위치 오차: {:.2e} mm",
        "p_theory": "thick-lens 이론 초점(근축): f={:.3f}mm, 후방초점거리={:.3f}mm, 초점 x={:.3f}",
        "p_paraxial": "시뮬레이션 근축 초점(|y0|<=20mm) x={:.3f}, 이론 대비 {:+.3f}mm ({:+.3f}%)",
        "p_marginal": "시뮬레이션 최외곽 초점(|y0|=100mm) x={:.3f}, 근축 대비 {:+.3f}mm (구면수차)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"),
        "font": EFONT, "legend": ELFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "fig1": "Fig 1. Eleven parallel rays converge after a biconvex lens",
        "forward": "Forward (enters surface {} first)",
        "retrace": "Retrace (exit ray reversed)",
        "fig2": "Fig 2. Reversing the exit ray retraces the incident path",
        "theory": "Thick-lens theory (x={:.1f}mm)",
        "paraxial": "Paraxial focus (x={:.1f}mm)",
        "marginal": "Marginal focus (x={:.1f}mm)",
        "fig3": "Fig 3. Paraxial theory versus finite-aperture ray-traced focus",
        "p_index": "n(N-BK7, {:.0f}nm) = {:.6f}",
        "p_entry": "forward entry surface={}, retrace entry surface={} (must be the opposite face)",
        "p_start": "original incident ray start: ({:.3f}, {:.3f}), direction ({:.3f}, {:.3f})",
        "p_end": "retrace arrival (first refraction point): ({:.3f}, {:.3f}), direction ({:.3f}, {:.3f})",
        "p_err": "position error: {:.2e} mm",
        "p_theory": "thick-lens theory (paraxial): f={:.3f}mm, BFD={:.3f}mm, focus x={:.3f}",
        "p_paraxial": "simulated paraxial focus (|y0|<=20mm) x={:.3f}, vs theory {:+.3f}mm ({:+.3f}%)",
        "p_marginal": "simulated marginal focus (|y0|=100mm) x={:.3f}, vs paraxial {:+.3f}mm (spherical aberration)",
    },
}

# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

R1, R2, T, D = 1000.0, -1000.0, 100.0, 500.0
WAVELENGTH_NM = 750.0
N_AIR = 1.0002778
ROI = (-1500, 1500, -500, 500)

lens = make_spherical_lens(R1=R1, R2=R2, t=T, D=D)
nk = load_nk_file(os.path.join(os.path.dirname(__file__), "N-BK7.nk"))
n_lens = refractive_index(nk, WAVELENGTH_NM)
source = make_parallel_beam(x=1000, y=0, vx=-1, vy=0, D=200, num=11)
paths = trace_parallel_beam_through_lens(lens, source, WAVELENGTH_NM, N_AIR, nk, ROI)

# 광학 가역성(reversibility) 검증
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
pos_err = np.hypot(path_retrace[-1, 0] - path_fwd[0, 0], path_retrace[-1, 1] - path_fwd[0, 1])

# thick-lens 이론값과 시뮬레이션 초점 위치
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


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def draw_lens_outline(ax, lens):
    b1, b2 = lens["BOUNDARY"]
    outline = np.vstack([b1, b2[::-1], b1[0:1]])
    ax.fill(outline[:, 0], outline[:, 1], color="C0", alpha=0.25)


def draw_paths(ax, paths, **kwargs):
    for p in paths:
        ax.plot(p[:, 0], p[:, 1], **kwargs)


def render(L):
    """주어진 라벨 묶음으로 그림 3개를 그리고 검증 수치를 출력한다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    print(L["p_index"].format(WAVELENGTH_NM, n_lens))

    # 그림1: 평행광이 양볼록 렌즈에 입사해 굴절되는 모습
    fig, ax = plt.subplots(figsize=(7, 4))
    draw_lens_outline(ax, lens)
    draw_paths(ax, paths, color="C3", linewidth=0.8)
    ax.set_xlim(ROI[0], ROI[1])
    ax.set_ylim(ROI[2], ROI[3])
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.set_title(L["fig1"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-lens-focusing.png"), dpi=150)
    plt.close(fig)

    # 그림2: 광학 가역성 검증
    print(L["p_entry"].format(entry_fwd, entry_retrace))
    print(L["p_start"].format(ray_fwd[0], ray_fwd[1], ray_fwd[2], ray_fwd[3]))
    print(L["p_end"].format(*path_retrace[-1]))
    print(L["p_err"].format(pos_err))

    fig, ax = plt.subplots(figsize=(7, 4))
    draw_lens_outline(ax, lens)
    ax.plot(full_fwd[:, 0], full_fwd[:, 1], "-", color="C0", linewidth=3, alpha=0.6,
            label=L["forward"].format(entry_fwd))
    ax.plot(full_retrace_xy[:, 0], full_retrace_xy[:, 1], "--", color="C3", linewidth=1.5,
            label=L["retrace"])
    ax.set_xlim(-1200, 1600)
    ax.set_ylim(-50, 150)
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.legend(prop=LF, loc="lower right")
    ax.set_title(L["fig2"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2-reversibility.png"), dpi=150)
    plt.close(fig)

    # 그림3: 시뮬레이션 초점 vs thick-lens 이론값
    print(L["p_theory"].format(f, bfd, expected_focus_x))
    print(L["p_paraxial"].format(
        sim_focus_x, sim_focus_x - expected_focus_x,
        100 * (sim_focus_x - expected_focus_x) / abs(expected_focus_x),
    ))
    print(L["p_marginal"].format(marginal_focus_x, marginal_focus_x - sim_focus_x))

    fig, ax = plt.subplots(figsize=(7, 4))
    draw_lens_outline(ax, lens)
    draw_paths(ax, paths, color="C3", linewidth=0.8, alpha=0.6)
    ax.axvline(expected_focus_x, color="k", linestyle="--", linewidth=1,
               label=L["theory"].format(expected_focus_x))
    ax.axvline(sim_focus_x, color="C0", linestyle=":", linewidth=1.5,
               label=L["paraxial"].format(sim_focus_x))
    ax.axvline(marginal_focus_x, color="C3", linestyle=":", linewidth=1.5,
               label=L["marginal"].format(marginal_focus_x))
    ax.set_xlim(-1300, -700)
    ax.set_ylim(-150, 150)
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.legend(prop=LF, loc="upper right")
    ax.set_title(L["fig3"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3-focus-vs-theory.png"), dpi=150)
    plt.close(fig)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

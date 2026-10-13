"""광선추적 1편 생성 그림 4개 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-spherical-lens-refraction/      (한국어)
      ../../assets/img/posts/raytracing-spherical-lens-refraction/en/   (영문)

- fig4-candidates-sag.png : 직선과 원의 두 교점(후보해)과 처짐량(sag)
- fig1-lens-focusing.png  : 평행광 11가닥이 양볼록 렌즈를 지나 모인다
- fig2-reversibility.png  : 광학 가역성 검증(경로와 두 경로의 높이 차이)
- fig3-focus-vs-theory.png: 입사 높이에 따른 초점 이동(종구면수차 곡선)

파일명의 숫자는 처음 만든 순서일 뿐이고, 글의 그림 번호는 캡션이 정한다.
그림 안에는 번호를 넣지 않는다(웹 그림이 끼면 번호가 밀리기 때문이다).

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
        # 후보해·처짐량 그림
        "cand_title_a": "(a) 광선을 늘인 직선은 원과 두 번 만난다",
        "cand_title_b": "(b) 꼭짓점 접평면과 처짐량",
        "cand_circle": "R2 면을 연장한 원",
        "cand_ray": "광선 (y = 100 mm)",
        "cand_box": "렌즈 면이 놓인 범위",
        "cand_keep": "arcsin 후보: 렌즈 면 위 (채택)",
        "cand_drop": "arccos 후보: 원의 반대편 (버림)",
        "cand_tangent": "꼭짓점 접평면",
        "cand_sag_rim": "가장자리 처짐량\n{:.2f} mm",
        "cand_sag_ray": "{:.2f} mm",
        # 그림: 렌즈 집속
        "focus_title": "평행광선 11가닥이 양볼록 렌즈를 지나 모인다",
        # 그림: 가역성
        "forward": "정방향 ({}번 면 먼저 진입)",
        "retrace": "역추적 (출사 광선을 반대로)",
        "rev_title": "출사 광선을 반대로 쏘면 입사 경로를 되짚는다",
        "rev_ylabel_b": "두 경로의 높이 차이 (mm)",
        "rev_lens": "렌즈",
        # 그림: 종구면수차
        "lsa_title": "입사 높이에 따른 초점 이동 (종구면수차)",
        "lsa_xlabel": "입사 높이 $y_0$ (mm)",
        "lsa_ylabel": "근축 초점에서 렌즈 쪽으로 당겨진 거리 (mm)",
        "lsa_points": "광선추적",
        "lsa_fit": "$k\\,y_0^2$  ($k$ = {:.3f}×10$^{{-3}}$ mm$^{{-1}}$, $y_0$ = 1 mm 광선에서 구함)",
        # 출력 문구
        "p_index": "n(N-BK7, {:.0f}nm) = {:.6f}",
        "p_entry": "정방향 진입면={}, 역추적 진입면={} (반대쪽 면이어야 함)",
        "p_start": "원래 입사 광선 시작점: ({:.3f}, {:.3f}), 방향 ({:.3f}, {:.3f})",
        "p_end": "역추적 도착점(렌즈 첫 굴절점): ({:.3f}, {:.3f}), 방향 ({:.3f}, {:.3f})",
        "p_err": "위치 오차: {:.2e} mm",
        "p_theory": "두꺼운 렌즈 이론 초점(근축): f={:.3f}mm, 후초점거리={:.3f}mm, 초점 x={:.3f}",
        "p_paraxial": "시뮬레이션 근축 초점(|y0|<=20mm) x={:.3f}, 이론 대비 {:+.3f}mm ({:+.3f}%)",
        "p_marginal": "시뮬레이션 최외곽 초점(|y0|=100mm) x={:.3f}, 근축 대비 {:+.3f}mm (구면수차)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"),
        "font": EFONT, "legend": ELFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "cand_title_a": "(a) The extended ray meets the circle twice",
        "cand_title_b": "(b) Vertex tangent plane and sag",
        "cand_circle": "Full circle of surface R2",
        "cand_ray": "Ray (y = 100 mm)",
        "cand_box": "Extent of the lens surface",
        "cand_keep": "arcsin candidate: on the lens (kept)",
        "cand_drop": "arccos candidate: far side (rejected)",
        "cand_tangent": "Vertex tangent plane",
        "cand_sag_rim": "Sag at rim\n{:.2f} mm",
        "cand_sag_ray": "{:.2f} mm",
        "focus_title": "Eleven parallel rays converge after a biconvex lens",
        "forward": "Forward (enters surface {} first)",
        "retrace": "Retrace (exit ray reversed)",
        "rev_title": "Reversing the exit ray retraces the incident path",
        "rev_ylabel_b": "Height difference (mm)",
        "rev_lens": "Lens",
        "lsa_title": "Focus shift versus ray height (longitudinal spherical aberration)",
        "lsa_xlabel": "Ray height $y_0$ (mm)",
        "lsa_ylabel": "Shift of focus toward the lens (mm)",
        "lsa_points": "Ray trace",
        "lsa_fit": "$k\\,y_0^2$  ($k$ = {:.3f}×10$^{{-3}}$ mm$^{{-1}}$, from the $y_0$ = 1 mm ray)",
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
N_AIR = 1.0  # SCHOTT 카탈로그 굴절률은 공기 기준 상대값이라 공기는 1로 둔다
ROI = (-1500, 1500, -500, 500)

lens = make_spherical_lens(R1=R1, R2=R2, t=T, D=D)
nk = load_nk_file(os.path.join(os.path.dirname(__file__), "N-BK7.nk"))
n_lens = refractive_index(nk, WAVELENGTH_NM)
source = make_parallel_beam(x=1000, y=0, vx=-1, vy=0, D=200, num=11)
paths = trace_parallel_beam_through_lens(lens, source, WAVELENGTH_NM, N_AIR, nk, ROI)

# 후보해 두 개: y=100 광선을 늘인 직선이 R2 면의 원과 만나는 두 점.
# refract_through_lens와 같은 식(arccos/arcsin)으로 구한다.
Y_DEMO = 100.0
_a, _b = 0.0, -1.0  # 방향 (-1, 0)인 광선의 a=-vy, b=vx
_c = _a * 1000.0 + _b * Y_DEMO
_cx, _cy = lens["cx2"], lens["cy2"]
_q = (_c - _a * _cx - _b * _cy) / (R2 * np.hypot(_a, _b))
w_c = np.arccos(_q) + np.arctan2(_b, _a)
w_s = np.arcsin(_q) - np.arctan2(_a, _b)
cand_c = (R2 * np.cos(w_c) + _cx, R2 * np.sin(w_c) + _cy)  # 원의 반대편
cand_s = (R2 * np.cos(w_s) + _cx, R2 * np.sin(w_s) + _cy)  # 렌즈 면 위
sag_rim = abs(R2) - np.sqrt(R2**2 - (D / 2) ** 2)
sag_demo = abs(R2) - np.sqrt(R2**2 - Y_DEMO**2)

# 광학 가역성(reversibility) 검증
# R2 쪽으로 들어온 광선을 굴절시킨 뒤(정방향, entry_surface=2), 그 출사
# 광선을 반대로 뒤집어 다시 쏘면(역방향, entry_surface=1 분기를 타게 됨)
# 원래 입사 경로를 정확히 되짚어 나와야 한다 -- 스넬의 법칙은 시간을
# 거꾸로 돌려도 꼴이 바뀌지 않기 때문이다. refract_through_lens의 두 분기
# (entry_surface==1 / else)가 서로 모순 없이 같은 물리를 구현했는지 보는
# 교차검증이 된다.
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

# 두 경로의 높이 차이: 정방향 경로(출사 후 500 mm까지 포함)와 역추적 경로를
# 같은 x 위치에서 비교한다. 두 경로 모두 x에 대해 단조이므로 선분 보간으로 충분하다.
full_fwd_xy = np.vstack([full_fwd[:, :2], far_point])
fx, fy = full_fwd_xy[::-1, 0], full_fwd_xy[::-1, 1]  # x 오름차순
rx, ry = full_retrace_xy[:, 0], full_retrace_xy[:, 1]  # 이미 x 오름차순
x_cmp = np.linspace(far_point[0], 1000.0, 2001)
dy = np.abs(np.interp(x_cmp, fx, fy) - np.interp(x_cmp, rx, ry))
dy_max = float(dy.max())

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
# 가장자리 광선은 구면수차 때문에 초점 위치 자체가 다르다.
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
marginal_focus_x = axis_crossing(90, 100)  # 최외곽: |y0| = 100mm (빔 가장자리)


def focus_shift(y0):
    """높이 y0로 들어온 광선이 광축을 지나는 점이 근축 이론 초점보다 렌즈 쪽으로
    얼마나 당겨졌는지(mm, 렌즈 쪽이 +)."""
    ray = np.array([1000.0, y0, -1.0, 0.0])
    p = refract_through_lens(ray, find_entry_surface(ray, lens)[0], lens, n_lens, N_AIR)
    xe, ye, vxe, vye = p[-1]
    return (xe - ye * vxe / vye) - expected_focus_x


lsa_k = focus_shift(1.0) / 1.0**2  # 광축에 아주 가까운 광선에서 구한 y0^2 계수
lsa_y = np.arange(5.0, 100.0 + 1e-9, 5.0)
lsa_shift = np.array([focus_shift(y) for y in lsa_y])


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def draw_lens_outline(ax, lens, **kwargs):
    # BOUNDARY[0]은 위→아래, BOUNDARY[1]은 아래→위로 저장되어 있으므로
    # 그대로 이으면 닫힌 단면이 된다. 한쪽을 뒤집으면 닫는 선이 가운데서
    # X자로 엇갈려 렌즈가 잎사귀 두 장처럼 그려진다.
    b1, b2 = lens["BOUNDARY"]
    outline = np.vstack([b1, b2, b1[0:1]])
    style = dict(color="C0", alpha=0.25)
    style.update(kwargs)
    ax.fill(outline[:, 0], outline[:, 1], **style)


def draw_paths(ax, paths, **kwargs):
    for p in paths:
        ax.plot(p[:, 0], p[:, 1], **kwargs)


def render(L):
    """주어진 라벨 묶음으로 그림 4개를 그리고 검증 수치를 출력한다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    print(L["p_index"].format(WAVELENGTH_NM, n_lens))

    # 후보해 두 개와 처짐량
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(8.4, 4), gridspec_kw={"width_ratios": [2.2, 1.15]})
    th = np.linspace(0, 2 * np.pi, 721)
    axa.plot(abs(R2) * np.cos(th) + _cx, abs(R2) * np.sin(th) + _cy, "--", color="0.6",
             linewidth=1, label=L["cand_circle"])
    draw_lens_outline(axa, lens)
    b2 = lens["BOUNDARY"][1]
    off2 = lens["BD_OFFSET"][1]
    bx0, bx1 = b2[:, 0].min() - off2, b2[:, 0].max() + off2
    by0, by1 = b2[:, 1].min() - off2, b2[:, 1].max() + off2
    axa.add_patch(plt.Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False,
                                edgecolor="C2", linewidth=1.2, label=L["cand_box"]))
    axa.annotate("", xy=(-2100, Y_DEMO), xytext=(1000, Y_DEMO),
                 arrowprops=dict(arrowstyle="-|>", color="C3", linewidth=1))
    axa.plot([], [], color="C3", linewidth=1, label=L["cand_ray"])
    axa.plot(*cand_s, "o", color="C3", markersize=7, label=L["cand_keep"])
    axa.plot(*cand_c, "o", markerfacecolor="white", markeredgecolor="C3", markersize=7,
             label=L["cand_drop"])
    axa.set_xlim(-2150, 1050)
    axa.set_ylim(-1100, 1100)
    axa.set_aspect("equal")
    axa.set_xlabel(L["xlabel"], **F)
    axa.set_ylabel(L["ylabel"], **F)
    axa.set_title(L["cand_title_a"], fontsize=10, **F)
    axa.grid(True, alpha=0.3)
    axa.legend(prop=dict(LF, size=7), loc="lower left")

    # (b) R2 면 확대: 꼭짓점(x=+t/2)의 접평면과 처짐량
    yy = np.linspace(0, D / 2, 200)
    xx = (lens["cx2"] + np.sqrt(R2**2 - yy**2))  # R2 면: 중심 cx2=-950, 오른쪽 호
    axb.fill_betweenx(yy, -10, xx, color="C0", alpha=0.25)
    axb.plot(xx, yy, color="C0", linewidth=1.2)
    vx_ = T / 2
    axb.axvline(vx_, color="k", linestyle="--", linewidth=1)
    axb.text(vx_ + 2, 112, L["cand_tangent"], rotation=90, fontsize=7.5, va="bottom", **F)
    axb.annotate("", xy=(vx_ - sag_rim, D / 2), xytext=(vx_, D / 2),
                 arrowprops=dict(arrowstyle="<->", color="C2"))
    axb.text(vx_ - sag_rim / 2, D / 2 - 30, L["cand_sag_rim"].format(sag_rim),
             ha="center", va="top", fontsize=8, **F)
    axb.annotate("", xy=(vx_ - sag_demo, Y_DEMO), xytext=(vx_, Y_DEMO),
                 arrowprops=dict(arrowstyle="<->", color="C3"))
    axb.text(vx_ - sag_demo - 2, Y_DEMO + 8, L["cand_sag_ray"].format(sag_demo), fontsize=8, color="C3", ha="right", **F)
    axb.plot(*cand_s, "o", color="C3", markersize=5)
    axb.set_xlim(0, 75)
    axb.set_ylim(0, 270)
    axb.set_xlabel(L["xlabel"], **F)
    axb.set_title(L["cand_title_b"], fontsize=9, **F)
    axb.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig4-candidates-sag.png"), dpi=150)
    plt.close(fig)

    # 평행광이 양볼록 렌즈에 입사해 굴절되는 모습
    fig, ax = plt.subplots(figsize=(7, 4))
    draw_lens_outline(ax, lens)
    draw_paths(ax, paths, color="C3", linewidth=0.8)
    ax.set_xlim(ROI[0], ROI[1])
    ax.set_ylim(ROI[2], ROI[3])
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.set_title(L["focus_title"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-lens-focusing.png"), dpi=150)
    plt.close(fig)

    # 광학 가역성 검증: 위는 경로(세로 확대), 아래는 두 경로의 높이 차이
    print(L["p_entry"].format(entry_fwd, entry_retrace))
    print(L["p_start"].format(ray_fwd[0], ray_fwd[1], ray_fwd[2], ray_fwd[3]))
    print(L["p_end"].format(*path_retrace[-1]))
    print(L["p_err"].format(pos_err))
    print(f"max |dy| along path: {dy_max:.2e} mm")

    fig, (ax, axd) = plt.subplots(2, 1, figsize=(7, 5), sharex=True,
                                  gridspec_kw={"height_ratios": [2, 1]})
    ax.axvspan(-T / 2, T / 2, color="C0", alpha=0.15, label=L["rev_lens"])
    ax.plot(full_fwd_xy[:, 0], full_fwd_xy[:, 1], "-", color="C0", linewidth=4, alpha=0.5,
            label=L["forward"].format(entry_fwd))
    ax.plot(full_retrace_xy[:, 0], full_retrace_xy[:, 1], "--", color="C3", linewidth=1.5,
            label=L["retrace"])
    ax.set_xlim(-600, 1050)
    ax.set_ylim(30, 115)
    ax.set_ylabel(L["ylabel"], **F)
    ax.legend(prop=dict(LF, size=8), loc="lower right")
    ax.set_title(L["rev_title"], **F)
    ax.grid(True, alpha=0.3)
    axd.semilogy(x_cmp, np.maximum(dy, 1e-17), color="k", linewidth=1)
    axd.axvspan(-T / 2, T / 2, color="C0", alpha=0.15)
    axd.set_ylim(1e-17, 1e-11)
    axd.set_xlabel(L["xlabel"], **F)
    axd.set_ylabel(L["rev_ylabel_b"], fontsize=8, **F)
    axd.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2-reversibility.png"), dpi=150)
    plt.close(fig)

    # 시뮬레이션 초점 vs thick-lens 이론값 (출력 블록은 본문에 그대로 실린다)
    print(L["p_theory"].format(f, bfd, expected_focus_x))
    print(L["p_paraxial"].format(
        sim_focus_x, sim_focus_x - expected_focus_x,
        100 * (sim_focus_x - expected_focus_x) / abs(expected_focus_x),
    ))
    print(L["p_marginal"].format(marginal_focus_x, marginal_focus_x - sim_focus_x))
    print(f"LSA k = {lsa_k:.6e} /mm; shift(20) = {focus_shift(20.0):.3f}; shift(100) = {focus_shift(100.0):.3f}; "
          f"k-spread = {100 * np.ptp(lsa_shift / lsa_y**2) / lsa_k:.2f}%")

    # 입사 높이에 따른 초점 이동(종구면수차 곡선)
    fig, ax = plt.subplots(figsize=(7, 4))
    yfine = np.linspace(0, 100, 201)
    ax.plot(yfine, lsa_k * yfine**2, "-", color="C0", linewidth=1.2,
            label=L["lsa_fit"].format(lsa_k * 1e3))
    ax.plot(lsa_y, lsa_shift, "o", color="C3", markersize=5, label=L["lsa_points"])
    for y0, tx, ty, ha in ((20.0, 3.0, 3.2, "left"), (100.0, 96.0, 17.0, "right")):
        s = focus_shift(y0)
        ax.annotate(f"$y_0$ = {y0:.0f} mm: {s:.3f} mm", xy=(y0, s), xytext=(tx, ty),
                    ha=ha, fontsize=8,
                    arrowprops=dict(arrowstyle="-", color="0.4", linewidth=0.8))
    ax.set_xlim(0, 105)
    ax.set_ylim(0, 18)
    ax.set_xlabel(L["lsa_xlabel"], **F)
    ax.set_ylabel(L["lsa_ylabel"], fontsize=9, **F)
    ax.legend(prop=dict(LF, size=8), loc="upper left")
    ax.set_title(L["lsa_title"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3-focus-vs-theory.png"), dpi=150)
    plt.close(fig)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

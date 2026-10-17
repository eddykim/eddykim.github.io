"""광선추적 2편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-total-internal-reflection/      (한국어)
      ../../assets/img/posts/raytracing-total-internal-reflection/en/   (영문)

광선추적 계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의
그림은 데이터가 완전히 동일하고 표기만 다르다.

그림 번호는 포스트의 캡션만 가진다. 그림 안 제목에는 번호를 넣지 않는다
(웹 그림이 끼어들면 번호가 밀리기 때문이다).
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

# 한글 텍스트에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면 눈금의 마이너스
# 기호가 한글 폰트에 없어 깨진다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "raytracing-total-internal-reflection",
)

WAVELENGTH_NM = 750.0

LABELS = {
    "ko": {
        "dir": BASE_DIR,
        "font": KFONT, "legend": LFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "verbose": {"transmit": "투과", "tir": "전반사(TIR)"},
        "ann_in": "입사\n(수직, 0$^\\circ$)",
        "ann_tir": "전반사\n내부입사각 60$^\\circ$",
        "ann_out": "무굴절 출사\n(내부입사각 0$^\\circ$)",
        "fig1": "수직 입사 광선이 정삼각형 프리즘 안에서 한 번 전반사한다",
        "tir_short": "전반사", "transmit_short": "투과",
        "sweep_legend": "입사각={:g}$^\\circ$  ({})",
        "fig2": "입사각이 경계값({:.4f}$^\\circ$)을 넘으면 전반사 대신 투과한다",
        "lens_pass": "투과", "lens_tir": "출사면에서 전반사(반사 뒤 경로까지 추적)",
        "fig3": "평볼록 렌즈(R2=-50mm)의 최외곽 광선은 출사면에서 전반사한다",
        "p_index": "n(N-BK7, {:.0f}nm) = {:.6f}",
        "p_thresh": "theta_c = {:.4f}deg, theta1 임계값(이론) = {:.6f}deg",
        "p_count": "광선 {}개 중 출사면에서 전반사한 광선 수: {}",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"),
        "font": EFONT, "legend": ELFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "verbose": {"transmit": "transmitted", "tir": "total internal reflection (TIR)"},
        "ann_in": "incident\n(normal, 0°)",
        "ann_tir": "TIR\ninternal angle 60°",
        "ann_out": "exits undeviated\n(internal angle 0°)",
        "fig1": "A normally incident ray undergoes one TIR inside an equilateral prism",
        "tir_short": "TIR", "transmit_short": "transmitted",
        "sweep_legend": "incidence={:g}°  ({})",
        "fig2": "Above the threshold ({:.4f}°) the ray transmits instead of reflecting",
        "lens_pass": "transmitted", "lens_tir": "TIR at the exit surface (traced onward)",
        "fig3": "In a plano-convex lens (R2=-50mm) the outermost rays undergo TIR on exit",
        "p_index": "n(N-BK7, {:.0f}nm) = {:.6f}",
        "p_thresh": "theta_c = {:.4f}deg, theta1 threshold (theory) = {:.6f}deg",
        "p_count": "of {} rays, {} undergo TIR at the exit surface",
    },
}


def extend_path(path, length=150.0):
    """마지막 점에서 마지막 방향으로 length만큼 더 이어붙인다 (그림에서
    출사 후 광선이 화면 밖으로 빠져나가는 걸 보여주기 위한 시각화 전용 연장)."""
    x, y, vx, vy = path[-1]
    return np.vstack([path, [x + length * vx, y + length * vy, vx, vy]])


def draw_prism_outline(ax, prism):
    verts = prism["verts"]
    outline = np.vstack([verts, verts[0:1]])
    ax.fill(outline[:, 0], outline[:, 1], color="C0", alpha=0.2)


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

nk_prism = load_nk_prism("N-BK7.nk")
n_glass = refractive_index_prism(nk_prism, WAVELENGTH_NM)

prism = make_prism([0, 300 * np.sqrt(3), 0], [300, 0, -300])
ray0 = np.array([-150.0, 100.0, 1.0, 0.0])

theta_c = np.arcsin(1.0 / n_glass)
A_apex = np.deg2rad(60.0)
theta1_thresh = np.rad2deg(np.arcsin(n_glass * np.sin(A_apex - theta_c)))


def fire(theta1_deg, y_hit=100.0, x0=-150.0):
    theta1 = np.deg2rad(theta1_deg)
    vx, vy = np.cos(theta1), np.sin(theta1)
    y0 = y_hit - (0 - x0) * np.tan(theta1)
    ray = np.array([x0, y0, vx, vy])
    return trace_ray_through_prism(ray, prism, n_glass, 1.0)


SWEEP_ANGLES = [0.0, 10.0, 20.0, 28.813, 35.0, 45.0]
sweep = [(ang, extend_path(fire(ang)), len(fire(ang)) >= 4) for ang in SWEEP_ANGLES]

nk_lens = load_nk_lens("N-BK7.nk")
n_lens = refractive_index_lens(nk_lens, WAVELENGTH_NM)
lens = make_spherical_lens(R1=1e6, R2=-50.0, t=20.0, D=78.0)
source = make_parallel_beam(x=-1000, y=0, vx=1, vy=0, D=76, num=39)

def trace_inside_lens(p, v, lens, n_glass, n_air=1.0, max_hits=10, eps=1e-7):
    """fig3(포스트의 그림5) 전용: 렌즈 안의 광선을 경계와 더 이상 만나지 않을 때까지 추적한다.

    refract_through_lens는 진입 1회 + 출사 1회만 계산하므로, 출사면에서 전반사한
    광선이 그 뒤 어디로 가는지는 알려주지 않는다. 그림에서 그 경로를 보여주려고
    프리즘 함수와 같은 방식(가장 가까운 경계 -> 굴절 또는 v'=v-2(v.n)n 반사)으로
    두 구면과 옆면(y=+-D/2)을 반복해서 만난다. 옆면은 매끄러운 굴절면으로 둔다.
    첫 교차점과 방향은 refract_through_lens의 출사 행과 같아야 한다(아래에서 확인).
    """
    h = lens["D"] / 2
    arcs = [(lens["cx1"], abs(lens["R1"]), -np.sign(lens["R1"])),
            (lens["cx2"], abs(lens["R2"]), -np.sign(lens["R2"]))]
    x_rim = [c + s * np.sqrt(R * R - h * h) for c, R, s in arcs]
    p, v = np.asarray(p, float), np.asarray(v, float) / np.linalg.norm(v)
    rows = []
    for _ in range(max_hits):
        best = None
        for c, R, s in arcs:
            d = p - np.array([c, 0.0])
            b, cc = d @ v, d @ d - R * R
            disc = b * b - cc
            if disc < 0:
                continue
            for t in (-b - np.sqrt(disc), -b + np.sqrt(disc)):
                q = p + t * v
                if t > eps and abs(q[1]) <= h and np.sign(q[0] - c) == s:
                    if best is None or t < best[0]:
                        best = (t, q, (q - np.array([c, 0.0])) / R)
        if abs(v[1]) > 1e-15:
            for yy in (h, -h):
                t = (yy - p[1]) / v[1]
                q = p + t * v
                if t > eps and min(x_rim) <= q[0] <= max(x_rim):
                    if best is None or t < best[0]:
                        best = (t, q, np.array([0.0, np.sign(yy)]))
        if best is None:
            break
        _, q, n_out = best  # n_out: 유리 밖을 향하는 단위 법선
        cosi = v @ n_out
        eta = n_glass / n_air
        k = 1 - eta * eta * (1 - cosi * cosi)
        if k < 0:  # 전반사
            v = v - 2 * cosi * n_out
            rows.append([q[0], q[1], v[0], v[1]])
            p = q
            continue
        v = eta * v + (np.sqrt(k) - eta * cosi) * n_out
        v = v / np.linalg.norm(v)
        rows.append([q[0], q[1], v[0], v[1]])
        break
    return np.array(rows)


lens_rays = []
n_tir = 0
max_dev = 0.0
for ray in source:
    entry = find_entry_surface(ray, lens)
    if entry is None:
        continue
    surf_idx, _ = entry
    inside = refract_through_lens_fixed(ray, surf_idx, lens, n_lens, 1.0)
    tir = exit_surface_is_tir(inside[0], surf_idx, lens, n_lens, 1.0)
    n_tir += tir
    onward = trace_inside_lens(inside[0, :2], inside[0, 2:], lens, n_lens)
    max_dev = max(max_dev, np.max(np.abs(onward[0] - inside[1])))
    path = np.vstack([ray, inside[:1], onward]) if tir else np.vstack([ray, inside])
    lens_rays.append((extend_path(path), tir))
print(f"trace_inside_lens vs refract_through_lens (first exit row) max |diff| = {max_dev:.1e}")


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L):
    """주어진 라벨 묶음으로 그림 3개를 그리고 검증 수치를 출력한다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    print(L["p_index"].format(WAVELENGTH_NM, n_glass))

    # fig1(포스트의 그림1): 정삼각형 프리즘, 수직 입사 광선이 한 면에서 전반사 후 다른 면으로 수직 출사
    path0 = extend_path(
        trace_ray_through_prism(ray0, prism, n_glass, 1.0, verbose=True, words=L["verbose"])
    )
    fig, ax = plt.subplots(figsize=(7, 5))
    draw_prism_outline(ax, prism)
    ax.plot(path0[:, 0], path0[:, 1], "-o", color="C3", linewidth=1.5, markersize=4)
    ax.annotate(L["ann_in"], (-190, 125), fontsize=9, va="bottom", **F)
    ax.annotate(L["ann_tir"], (330, 125), fontsize=9, va="bottom", **F)
    ax.annotate(L["ann_out"], xy=(433, -50), xytext=(300, -250), fontsize=9, va="top",
                arrowprops=dict(arrowstyle="-", color="0.4", linewidth=0.8), **F)
    ax.set_xlim(-200, 600)
    ax.set_ylim(-350, 350)
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.set_title(L["fig1"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-prism-tir-path.png"), dpi=150)
    plt.close(fig)

    # fig2(포스트의 그림3): 입사각 훑기 -- 임계각 아래에서는 전반사, 위에서는 투과
    print(L["p_thresh"].format(np.rad2deg(theta_c), theta1_thresh))
    fig, ax = plt.subplots(figsize=(9.2, 5))
    draw_prism_outline(ax, prism)
    cmap = plt.get_cmap("viridis")  # 연한 유리색과 구별되도록 어두운 구간만 쓴다
    for i, (ang, p, is_tir) in enumerate(sweep):
        color = cmap(0.85 * i / (len(sweep) - 1))
        word = L["tir_short"] if is_tir else L["transmit_short"]
        ax.plot(p[:, 0], p[:, 1], "-", color=color, linewidth=1.3,
                label=L["sweep_legend"].format(ang, word))
    ax.set_xlim(-200, 600)
    ax.set_ylim(-350, 350)
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.legend(prop=dict(LF, size=8), loc="center left", bbox_to_anchor=(1.02, 0.5))
    ax.set_title(L["fig2"].format(theta1_thresh), **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2-prism-angle-sweep.png"), dpi=150)
    plt.close(fig)

    # fig3(포스트의 그림5): 평볼록 렌즈 -- 일부 광선은 투과, 가장자리는 출사면에서 전반사한 뒤
    # 반사 경로까지 그린다(trace_inside_lens)
    print(L["p_count"].format(len(source), n_tir))
    fig, ax = plt.subplots(figsize=(8, 5))
    b1, b2 = lens["BOUNDARY"]
    # b1은 위->아래, b2는 아래->위 순서라 그대로 이어야 닫힌 윤곽이 된다
    outline = np.vstack([b1, b2, b1[0:1]])
    ax.fill(outline[:, 0], outline[:, 1], color="C0", alpha=0.2)
    for full, tir in lens_rays:
        ax.plot(full[:, 0], full[:, 1], "-", color=("C3" if tir else "C0"),
                linewidth=0.8, alpha=0.7)
    ax.plot([], [], "-", color="C0", label=L["lens_pass"])
    ax.plot([], [], "-", color="C3", label=L["lens_tir"])
    ax.set_xlim(-70, 70)
    ax.set_ylim(-45, 45)
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.legend(prop=LF, loc="upper right", fontsize=9)
    ax.set_title(L["fig3"], **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3-lens-marginal-tir.png"), dpi=150)
    plt.close(fig)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

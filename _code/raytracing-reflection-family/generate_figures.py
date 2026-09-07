"""기하광학 3편 그림 2개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/raytracing-reflection-family/      (한국어)
      ../../assets/img/posts/raytracing-reflection-family/en/   (영문)

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은
데이터가 완전히 동일하고 표기만 다르다.
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

# 한글 텍스트에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면 눈금의 마이너스
# 기호가 한글 폰트에 없어 깨진다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "raytracing-reflection-family",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR,
        "font": KFONT, "legend": LFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "reflected": "반사", "transmitted": "투과",
        "fig1": "그림1. 네 반사 부품의 광선 경로 -- 반사 공식은 하나, 법선 찾는 법만 다르다",
        "incident": "입사광",
        "n_dup": "세 번째 점을 재사용한 법선",
        "n_uniq": "세 점을 모두 쓴 법선",
        "fig2": "그림2. 세 번째 점을 재사용하면 반사 방향이 {:.1f}° 달라진다",
        "p_vertex": "대칭 꼭짓점: n_dup={}, n_uniq={}, 차이={:.3f}deg",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"),
        "font": EFONT, "legend": ELFONT,
        "xlabel": "x (mm)", "ylabel": "y (mm)",
        "reflected": "reflected", "transmitted": "transmitted",
        "fig1": "Fig 1. Ray paths for the four reflective components -- one formula, four ways to find the normal",
        "incident": "incident ray",
        "n_dup": "normal from the reused third point",
        "n_uniq": "normal from all three points",
        "fig2": "Fig 2. Reusing the third point turns the reflected ray by {:.1f}°",
        "p_vertex": "symmetric vertex: n_dup={}, n_uniq={}, diff={:.3f}deg",
    },
}


def extend(path, length=60.0):
    x, y, vx, vy = path[-1]
    return np.vstack([path, [x + length * vx, y + length * vy, vx, vy]])


def draw_ray(ax, path, **kwargs):
    p = extend(np.atleast_2d(path))
    ax.plot(p[:, 0], p[:, 1], **kwargs)


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

flat_mirror = make_flat_mirror(200, 0, 0, np.deg2rad(20))
flat_ray = np.array([-100.0, 30.0, 1.0, 0.0])
flat_out = reflect_flat_mirror(flat_ray, flat_mirror)

sph_mirror = make_spherical_mirror(-300, 200, 0, 0, 0)
sph_ray = np.array([-100.0, 40.0, 1.0, 0.0])
sph_out = reflect_spherical_mirror(sph_ray, sph_mirror)

bs_mirror = make_beam_splitter(60, 0, 0, np.deg2rad(45))
bs_ray = np.array([-100.0, 0.0, 1.0, 0.0])
bs_r, bs_t = reflect_beam_splitter(bs_ray, bs_mirror)

arb_mirror = make_arbitrary_mirror([-100, 0, 100, 200], [0, 50, 0, 50])
arb_ray = np.array([-50.0, 100.0, 0.0, -1.0])
arb_out = reflect_arbitrary_mirror(arb_ray, arb_mirror)

# 꼭짓점 법선 -- 대칭 꼭짓점에 수직 입사
vx0, vy0, vx1, vy1_, vx2, vy2 = -100, 0, 0, 50, 100, 0
n_dup = _vertex_normal(vx0, vy0, vx1, vy1_, vx2, vy2, vx1, vy1_, dup_point=True)
n_uniq = _vertex_normal(vx0, vy0, vx1, vy1_, vx2, vy2, vx1, vy1_, dup_point=False)
vertex_ang = np.rad2deg(np.arccos(np.clip(np.dot(n_dup, n_uniq), -1, 1)))
v_dup = reflect_vector([0.0, -1.0], n_dup)
v_uniq = reflect_vector([0.0, -1.0], n_uniq)


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L):
    """주어진 라벨 묶음으로 그림 2개를 그리고 검증 수치를 출력한다."""
    out_dir, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out_dir, exist_ok=True)

    # 그림1: 네 부품의 반사 경로
    fig, axes = plt.subplots(2, 2, figsize=(9, 8))

    ax = axes[0, 0]
    b = flat_mirror["BOUNDARY"][0]
    ax.plot(b[:, 0], b[:, 1], "k-", linewidth=3)
    draw_ray(ax, np.vstack([flat_ray, flat_out[0]]), color="C3", linewidth=1.2)
    ax.set_title("FlatMirror", **F)

    ax = axes[0, 1]
    b = sph_mirror["BOUNDARY"][0]
    ax.plot(b[:, 0], b[:, 1], "k-", linewidth=3)
    draw_ray(ax, np.vstack([sph_ray, sph_out[0]]), color="C3", linewidth=1.2)
    ax.set_title("SphericalMirror", **F)

    ax = axes[1, 0]
    b = bs_mirror["BOUNDARY"][0]
    ax.plot(b[:, 0], b[:, 1], "k-", linewidth=3)
    draw_ray(ax, np.vstack([bs_ray, bs_r[0]]), color="C3", linewidth=1.2, label=L["reflected"])
    draw_ray(ax, np.vstack([bs_ray, bs_t[0]]), color="C0", linewidth=1.2, label=L["transmitted"])
    ax.plot(bs_ray[0], bs_ray[1], "o", color="k", markersize=3)
    ax.legend(prop=LF, fontsize=8, loc="upper left")
    ax.set_title("BeamSplitter", **F)

    ax = axes[1, 1]
    for seg in arb_mirror["BOUNDARY"]:
        ax.plot(seg[:, 0], seg[:, 1], "k-", linewidth=3)
    draw_ray(ax, np.vstack([arb_ray, arb_out[1:]]), color="C3", linewidth=1.2)
    ax.set_title("ArbitraryMirror", **F)

    for ax in axes.flat:
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.3)
        ax.set_xlabel(L["xlabel"], **F)
        ax.set_ylabel(L["ylabel"], **F)

    fig.suptitle(L["fig1"], **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig1-reflection-family.png"), dpi=150)
    plt.close(fig)

    # 그림2: 다각형 꼭짓점 법선 -- 점을 재사용한 경우 vs 세 점을 모두 쓴 경우
    print(L["p_vertex"].format(np.round(n_dup, 3), np.round(n_uniq, 3), vertex_ang))

    fig, ax = plt.subplots(figsize=(6, 5))
    poly = np.array([[vx0, vy0], [vx1, vy1_], [vx2, vy2]])
    ax.plot(poly[:, 0], poly[:, 1], "k-", linewidth=3)
    ax.plot([0, 0], [150, 50], "-", color="0.5", linewidth=1.2)
    ax.annotate(L["incident"], (5, 120), fontsize=9, **F)
    seg_len = 80
    ax.plot([vx1, vx1 + seg_len * v_dup[0]], [vy1_, vy1_ + seg_len * v_dup[1]],
            "-", color="C3", linewidth=1.5, label=L["n_dup"])
    ax.plot([vx1, vx1 + seg_len * v_uniq[0]], [vy1_, vy1_ + seg_len * v_uniq[1]],
            "--", color="C0", linewidth=1.5, label=L["n_uniq"])
    ax.plot(vx1, vy1_, "o", color="k", markersize=4)
    ax.set_xlim(-120, 120)
    ax.set_ylim(-20, 160)
    ax.set_aspect("equal")
    ax.set_xlabel(L["xlabel"], **F)
    ax.set_ylabel(L["ylabel"], **F)
    ax.legend(prop=LF, fontsize=9, loc="lower right")
    ax.set_title(L["fig2"].format(vertex_ang), **F)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig2-vertex-normal-typo.png"), dpi=150)
    plt.close(fig)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

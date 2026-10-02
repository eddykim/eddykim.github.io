"""기하광학 배경이론 3편 그림 5개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/geometric-optics-stops-pupils-etendue/      (한국어)
      ../../assets/img/posts/geometric-optics-stops-pupils-etendue/en/   (영문)

그림 1·3·8 은 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rays import N_BK7_B, N_BK7_C, sellmeier
from stops import (Element, find_aperture_stop, image_of_aperture, marginal_and_chief,
                   pupil_irradiance, trace, vignetting_fraction)

KFONT = {"fontfamily": "AppleGothic"}
EFONT: dict = {}
plt.rcParams["axes.unicode_minus"] = False

BASE_DIR = os.path.join(os.path.dirname(__file__), "..", "..",
                        "assets", "img", "posts", "geometric-optics-stops-pupils-etendue")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c"

N_GLASS = sellmeier(0.5876, N_BK7_B, N_BK7_C)
Z_OBJ, Y_OBJ = -150.0, 10.0


def test_system():
    """두 N-BK7 양볼록 렌즈와 그 사이의 조리개 (단위 mm, 공기 n = 1)."""
    n = N_GLASS
    return [Element(0.0, "surface", (1.0, n, 60.0), 15.0),
            Element(8.0, "surface", (n, 1.0, -60.0), 15.0),
            Element(33.0, "stop", (), 5.0),
            Element(55.0, "surface", (1.0, n, 80.0), 15.0),
            Element(61.0, "surface", (n, 1.0, -80.0), 15.0)]


def image_plane(elements):
    _, _, m, _ = marginal_and_chief(elements, Z_OBJ, Y_OBJ)
    return m[-1, 0] - m[-1, 1] / m[-1, 2]


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT,
        "s_a": "(a) 조리개 조리개를 조이면", "s_b": "(b) 시야 조리개를 좁히면",
        "s_x": "상면 위 높이 y (mm)", "s_y": "상면에 도착하는 각 u (mrad)",
        "s_open": "열었을 때", "s_close": "조였을 때",
        "h_t": "Hecht 예제 5.6: 어느 구멍이 조리개인가",
        "h_x": "광축 위치 z (mm)", "h_y": "높이 (mm)",
        "h_obj": "물점 S", "h_lens": "렌즈 (지름 140 mm)", "h_hole": "구멍 (지름 40 mm)",
        "h_img": "구멍의 상 = 입사동\n(지름 200 mm, 허상)",
        "h_ang_l": "렌즈 테두리 19.3°", "h_ang_h": "입사동 9.5°",
        "r_t": "두 렌즈와 조리개: 주변광선과 주광선",
        "r_x": "광축 위치 z (mm)", "r_y": "높이 (mm)",
        "r_marg": "주변광선", "r_chief": "주광선", "r_stop": "조리개",
        "r_ep": "입사동", "r_xp": "출사동", "r_obj": "물체", "r_img": "상",
        "v_a": r"(a) 원판 광원이 만드는 조도와 $\cos^4$ 법칙",
        "v_ax": r"주광선 각 $\varphi$ ($^\circ$)", "v_ay": "상대 조도",
        "v_cos4": r"$\cos^4\varphi$", "v_small": "작은 동공 (a/L = 0.05)", "v_big": "큰 동공 (a/L = 0.5)",
        "v_b": "(b) 시험계의 비네팅과 상대 조도",
        "v_bx": "물체 높이 (mm)", "v_by": "비율",
        "v_vig": "동공을 통과하는 광선 비율", "v_ri": r"상대 조도 (비네팅 × $\cos^4$)",
        "t_t": "초점이 어긋날 때 잰 크기의 오차",
        "t_x": r"물체 위치 어긋남 $\delta$ (mm)", "t_y": "잰 크기의 상대 오차 (%)",
        "t_conv": "조리개가 렌즈에 있을 때", "t_tele": "조리개가 뒤 초점면에 있을 때 (물체 쪽 텔레센트릭)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT,
        "s_a": "(a) Closing the aperture stop", "s_b": "(b) Narrowing the field stop",
        "s_x": "Height y at the image plane (mm)", "s_y": "Arrival angle u (mrad)",
        "s_open": "open", "s_close": "closed",
        "h_t": "Hecht Example 5.6: which opening is the aperture stop?",
        "h_x": "Axial position z (mm)", "h_y": "Height (mm)",
        "h_obj": "Object point S", "h_lens": "Lens (140 mm dia.)", "h_hole": "Hole (40 mm dia.)",
        "h_img": "Image of the hole = entrance pupil\n(200 mm dia., virtual)",
        "h_ang_l": "Lens rim 19.3°", "h_ang_h": "Entrance pupil 9.5°",
        "r_t": "Two lenses and a stop: marginal and chief rays",
        "r_x": "Axial position z (mm)", "r_y": "Height (mm)",
        "r_marg": "Marginal ray", "r_chief": "Chief ray", "r_stop": "Stop",
        "r_ep": "Entrance pupil", "r_xp": "Exit pupil", "r_obj": "Object", "r_img": "Image",
        "v_a": r"(a) Irradiance from a disc source and the $\cos^4$ law",
        "v_ax": r"Chief-ray angle $\varphi$ ($^\circ$)", "v_ay": "Relative irradiance",
        "v_cos4": r"$\cos^4\varphi$", "v_small": "Small pupil (a/L = 0.05)", "v_big": "Large pupil (a/L = 0.5)",
        "v_b": "(b) Vignetting and relative illumination of the test system",
        "v_bx": "Object height (mm)", "v_by": "Fraction",
        "v_vig": "Fraction of rays through the pupil", "v_ri": r"Relative illumination (vignetting × $\cos^4$)",
        "t_t": "Measured-size error when the focus is off",
        "t_x": r"Object displacement $\delta$ (mm)", "t_y": "Relative error of measured size (%)",
        "t_conv": "Stop at the lens", "t_tele": "Stop at the back focal plane (object-side telecentric)",
    },
}


def legend(ax, L, size=8.5, **kw):
    ax.legend(prop={"family": L["font"].get("fontfamily", "DejaVu Sans"), "size": size}, **kw)


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    fig.savefig(os.path.join(L["dir"], name), dpi=150, facecolor="white")
    plt.close(fig)


# ── 그림2: 위상공간에서 본 두 조리개 ─────────────────────────
def band(f, a, b):
    """초점거리 f 렌즈(조리개 반지름 a)가 무한 물체를 상면에 맺을 때, 높이 |y| ≤ b 인 상점에 도착하는 광선의 (y, u) 영역."""
    y = np.array([-b, b, b, -b])
    u = np.array([(-b - a) / f, (b - a) / f, (b + a) / f, (-b + a) / f])
    return y, -u      # 렌즈 위 높이 h 에서 상점 y 로 가는 각은 (y - h)/f


def figure2(L):
    f = 50.0
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.0), sharey=True)
    cases = ((axes[0], "s_a", ((10.0, 12.0), (4.0, 12.0))),
             (axes[1], "s_b", ((10.0, 12.0), (10.0, 6.0))))
    for ax, key, ((a1, b1), (a2, b2)) in cases:
        y, u = band(f, a1, b1)
        ax.fill(y, -u * 1e3, color="#dddddd", label=L["s_open"])
        y, u = band(f, a2, b2)
        ax.fill(y, -u * 1e3, color=C_RAY, alpha=0.45, label=L["s_close"])
        ax.set_xlim(-15, 15); ax.set_ylim(-500, 500)
        ax.axhline(0, color=C_GREY, lw=0.4); ax.axvline(0, color=C_GREY, lw=0.4)
        ax.set_xlabel(L["s_x"], **L["font"])
        ax.set_title(L[key], fontsize=10, **L["font"])
        legend(ax, L, loc="upper left", size=8)
    axes[0].set_ylabel(L["s_y"], **L["font"])
    fig.tight_layout()
    save(fig, L, "fig2-phase-space-stops.png")


# ── 그림4: Hecht 예제 5.6 ──────────────────────────────────
def figure4(L):
    els = [Element(0.0, "thin", (100.0,), 70.0), Element(80.0, "stop", (), 20.0)]
    z_s = -200.0
    _, table = find_aperture_stop(els, z_s)
    ep = table[1][1]
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    ax.plot([0, 0], [-70, 70], color="k", lw=2.0)
    ax.plot([80, 80], [20, 110], color="k", lw=3); ax.plot([80, 80], [-110, -20], color="k", lw=3)
    ax.plot([ep.z, ep.z], [-ep.radius, ep.radius], color=C_KEY, lw=2.5, ls="--")
    ax.plot(z_s, 0, "ko")
    # 물점에서 렌즈 테두리로 가는 원뿔, 입사동 가장자리로 가는 원뿔
    for sgn in (1, -1):
        ax.plot([z_s, 0], [0, sgn * 70], color=C_GREY, lw=0.9, ls=":")
        ax.plot([z_s, ep.z], [0, sgn * ep.radius], color=C_KEY, lw=1.0)
    # 실제 광선: 물점 → 렌즈 → 구멍 가장자리 → 상점
    for sgn in (1, -1):
        u = sgn * ep.radius / (ep.z - z_s)
        pts = trace(els, z_s, 0.0, u, z_end=320.0)
        ax.plot(pts[:, 0], pts[:, 1], color=C_RAY, lw=1.2)
    ax.text(z_s, -14, L["h_obj"], ha="center", **L["font"])
    ax.text(-6, -95, L["h_lens"], fontsize=9, ha="right", **L["font"])
    ax.text(85, 100, L["h_hole"], fontsize=9, **L["font"])
    ax.text(ep.z - 5, ep.radius + 8, L["h_img"], ha="right", fontsize=9, color=C_KEY, **L["font"])
    ax.text(-150, 42, L["h_ang_l"], fontsize=9, color=C_GREY, **L["font"])
    ax.text(-150, 12, L["h_ang_h"], fontsize=9, color=C_KEY, **L["font"])
    ax.set_xlim(-220, 430); ax.set_ylim(-125, 135)
    ax.set_xlabel(L["h_x"], **L["font"]); ax.set_ylabel(L["h_y"], **L["font"])
    ax.set_title(L["h_t"], fontsize=10, **L["font"])
    fig.tight_layout()
    save(fig, L, "fig4-hecht-entrance-pupil.png")


# ── 그림5: 주변광선과 주광선 ─────────────────────────────────
def lens_shape(z0, t, r1, r2, h):
    y = np.linspace(-h, h, 100)
    za = z0 + r1 - np.sqrt(r1**2 - y**2)
    zb = z0 + t + r2 + np.sqrt(r2**2 - y**2)
    return np.r_[za, zb[::-1], za[:1]], np.r_[y, y[::-1], y[:1]]


def figure5(L):
    els = test_system()
    z_img = image_plane(els)
    stop, ep, marg, chief = marginal_and_chief(els, Z_OBJ, Y_OBJ, z_end=z_img)
    xp = image_of_aperture(els, stop, "image")
    fig, ax = plt.subplots(figsize=(11.5, 4.6))
    for z0, t, r1, r2 in ((0.0, 8.0, 60.0, -60.0), (55.0, 6.0, 80.0, -80.0)):
        zz, yy = lens_shape(z0, t, r1, r2, 15.0)
        ax.fill(zz, yy, color="#cfe3f3", ec="k", lw=0.8)
    ax.plot([33, 33], [5, 17], color="k", lw=3); ax.plot([33, 33], [-17, -5], color="k", lw=3)
    ax.text(33, 18, L["r_stop"], ha="center", fontsize=9, **L["font"])
    for pts, col, lab in ((marg, C_RAY, L["r_marg"]), (chief, C_RED, L["r_chief"])):
        ax.plot(pts[:, 0], pts[:, 1], color=col, lw=1.4, label=lab)
    # 입사동: 물체 쪽 광선의 연장선이 향하는 곳, 출사동: 상 쪽 광선의 연장선이 오는 곳
    for pts, col in ((marg, C_RAY), (chief, C_RED)):
        z0, y0, u0 = pts[0, :3]
        ax.plot([z0, ep.z], [y0, y0 + (ep.z - z0) * u0], color=col, lw=0.8, ls=":")
        z1, y1, u1 = pts[-1, :3]
        ax.plot([xp.z, z1], [y1 - (z1 - xp.z) * u1, y1], color=col, lw=0.8, ls=":")
    ax.plot([ep.z, ep.z], [-ep.radius, ep.radius], color=C_KEY, lw=2.5, alpha=0.8)
    ax.text(ep.z + 2, -ep.radius - 4, L["r_ep"], color=C_KEY, fontsize=9, **L["font"])
    ax.plot([xp.z, xp.z], [-xp.radius, xp.radius], color=C_GRN, lw=2.5, alpha=0.8)
    ax.text(xp.z - 2, xp.radius + 2, L["r_xp"], color=C_GRN, fontsize=9, ha="right", **L["font"])
    ax.annotate("", xy=(Z_OBJ, Y_OBJ), xytext=(Z_OBJ, 0), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.text(Z_OBJ, -3, L["r_obj"], ha="center", fontsize=9, **L["font"])
    ax.annotate("", xy=(z_img, chief[-1, 1]), xytext=(z_img, 0), arrowprops=dict(arrowstyle="->", lw=1.5))
    ax.text(z_img, 2, L["r_img"], ha="center", fontsize=9, **L["font"])
    ax.axhline(0, color=C_GREY, lw=0.5)
    ax.set_xlim(-160, 100); ax.set_ylim(-20, 22)
    ax.set_xlabel(L["r_x"], **L["font"]); ax.set_ylabel(L["r_y"], **L["font"])
    ax.set_title(L["r_t"], fontsize=10, **L["font"])
    legend(ax, L, loc="lower left", size=8.5)
    fig.tight_layout()
    save(fig, L, "fig5-marginal-chief.png")


# ── 그림6: cos⁴ 법칙과 비네팅 ────────────────────────────────
def figure6(L):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12.0, 4.2))
    phi = np.linspace(0, 50, 26)
    Ldist = 1.0
    a1.plot(phi, np.cos(np.deg2rad(phi)) ** 4, color="k", lw=1.2, label=L["v_cos4"])
    for a, col, key in ((0.05, C_RAY, "v_small"), (0.5, C_RED, "v_big")):
        e0 = pupil_irradiance(a, Ldist, 0.0)
        e = [pupil_irradiance(a, Ldist, Ldist * np.tan(np.deg2rad(p))) / e0 for p in phi]
        a1.plot(phi, e, "o", ms=4, color=col, label=L[key])
    a1.set_xlabel(L["v_ax"], **L["font"]); a1.set_ylabel(L["v_ay"], **L["font"])
    a1.set_title(L["v_a"], fontsize=10, **L["font"])
    legend(a1, L, loc="lower left", size=8)
    els = test_system()
    z_img = image_plane(els)
    ep = image_of_aperture(els, 2, "object")
    ys = np.linspace(0, 70, 36)
    vig, ri = [], []
    for y in ys:
        v = vignetting_fraction(els, Z_OBJ, y, ep, n_grid=301)
        _, _, _, chief = marginal_and_chief(els, Z_OBJ, max(y, 1e-6), z_end=z_img)
        vig.append(v)
        ri.append(v * np.cos(np.arctan(chief[-1, 2])) ** 4)
    a2.plot(ys, vig, color=C_RAY, lw=1.5, label=L["v_vig"])
    a2.plot(ys, ri, color=C_KEY, lw=1.5, ls="--", label=L["v_ri"])
    a2.set_xlabel(L["v_bx"], **L["font"]); a2.set_ylabel(L["v_by"], **L["font"])
    a2.set_title(L["v_b"], fontsize=10, **L["font"])
    legend(a2, L, loc="lower left", size=8)
    fig.tight_layout()
    save(fig, L, "fig6-cos4-vignetting.png")


# ── 그림7: 텔레센트릭 ───────────────────────────────────────
def measured_size(stop_z, delta, f=50.0, y=5.0):
    """얇은 렌즈(z=0, 초점거리 f) + 조리개(stop_z). 물체 2f 앞, 센서 2f 뒤 고정.

    물체가 δ 만큼 어긋나면 상이 흐려지지만, 흐린 점의 중심(주광선)이 센서에 닿는 높이로 크기를 잰다.
    """
    z_obj = -2 * f + delta
    els = [Element(0.0, "thin", (f,), None)]
    # 주광선은 조리개 중심을 지난다: 물체 공간에서 조리개 중심의 상(입사동 중심)을 향한다
    if stop_z == 0.0:
        z_ep = 0.0
    else:
        els_s = els + [Element(stop_z, "stop", (), 1.0)]
        z_ep = image_of_aperture(els_s, 1, "object").z
    if np.isinf(z_ep) or abs(z_ep) > 1e12:
        u = 0.0
    else:
        u = -y / (z_ep - z_obj)
    pts = trace(els, z_obj, y, u, z_end=2 * f)
    return -pts[-1, 1]       # 상은 뒤집힌다


def figure7(L):
    f = 50.0
    d = np.linspace(-5, 5, 41)
    ref = measured_size(0.0, 0.0)
    conv = [(measured_size(0.0, x) / ref - 1) * 100 for x in d]
    tele = [(measured_size(f * (1 - 1e-12), x) / ref - 1) * 100 for x in d]
    fig, ax = plt.subplots(figsize=(7.5, 4.0))
    ax.plot(d, conv, color=C_RED, lw=1.5, label=L["t_conv"])
    ax.plot(d, tele, color=C_RAY, lw=1.5, label=L["t_tele"])
    ax.axhline(0, color=C_GREY, lw=0.5)
    ax.set_xlabel(L["t_x"], **L["font"]); ax.set_ylabel(L["t_y"], **L["font"])
    ax.set_title(L["t_t"], fontsize=10, **L["font"])
    legend(ax, L, loc="upper left", size=8)
    fig.tight_layout()
    save(fig, L, "fig7-telecentric.png")


def main():
    for L in LABELS.values():
        figure2(L)
        figure4(L)
        figure5(L)
        figure6(L)
        figure7(L)
    els = test_system()
    z_img = image_plane(els)
    stop, ep, marg, chief = marginal_and_chief(els, Z_OBJ, Y_OBJ, z_end=z_img)
    xp = image_of_aperture(els, stop, "image")
    print(f"조리개 = 요소 {stop}, 입사동 z {ep.z:.3f} r {ep.radius:.3f}, 출사동 z {xp.z:.3f} r {xp.radius:.3f}, "
          f"상면 z {z_img:.3f}, 상 높이 {chief[-1, 1]:.3f}")


if __name__ == "__main__":
    main()

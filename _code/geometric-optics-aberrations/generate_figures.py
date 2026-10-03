"""기하광학 배경이론 4편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/geometric-optics-aberrations/      (한국어)
      ../../assets/img/posts/geometric-optics-aberrations/en/   (영문)

그림 1·7·10·12 는 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from aberr import (LINE_C, LINE_D, LINE_F, axis_crossing, fit_seidel_terms, opd, paraxial_image_height,
                   paraxial_trace, psf, psf_pixel, pupil_grid, seidel_sums, seidel_surface_terms, setup,
                   spot, strehl_center, trace_real, to_plane, wave_coefficients, zernike_radial)
from systems import (F, N_D, WL, achromat, center_stop_z, bent_lens, cooke_triplet, equiconvex, plano_convex,
                     singlet_for_color)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "geometric-optics-aberrations")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"
LAM_MM = WL * 1e-3


def t(L, ko, en):
    """언어에 맞는 문구. AppleGothic 의 ° 는 간격이 어긋나므로 수식 기호로 바꾼다."""
    s = ko if L["lang"] == "ko" else en
    return s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$")


def font(L):
    return {"fontfamily": ["AppleGothic", "DejaVu Sans"]} if L["lang"] == "ko" else {}


def setfonts(L):
    plt.rcParams["font.family"] = ["AppleGothic", "DejaVu Sans"] if L["lang"] == "ko" else ["DejaVu Sans"]


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    fig.savefig(os.path.join(L["dir"], name), dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("  저장:", os.path.join(L["dir"], name))


# ------------------------------------------------------------------ 그림 2: 파면수차와 광선수차

def figure2(L):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw=dict(width_ratios=[1.15, 1]))
    # (a) 과장한 개념도: 출사동(z=0)에서 상점(z=Rr)으로 모이는 기준 구면과, ρ⁴ 항이 더해진 실제 파면
    Rr, a_p = 10.0, 3.0
    yy = np.linspace(-a_p, a_p, 301)
    ref = Rr - np.sqrt(Rr**2 - yy**2)                       # 기준 구면 (꼭짓점 z = 0)
    k4 = 0.0022
    act = ref + k4 * yy**4                                  # 가장자리가 상 쪽으로 앞선 파면 (W > 0)
    a.plot(ref, yy, color=C_GREY, lw=2, ls="--", label=t(L, "기준 구면 (상점 P*로 모이는 구면파)", "reference sphere (converges on P*)"))
    a.plot(act, yy, color=C_RAY, lw=2.2, label=t(L, "실제 파면", "actual wavefront"))
    for y0 in [-2.8, -2.0, -1.0, 0.0, 1.0, 2.0, 2.8]:
        z0 = Rr - np.sqrt(Rr**2 - y0**2) + k4 * y0**4
        dz = y0 / np.sqrt(Rr**2 - y0**2) + 4 * k4 * y0**3   # 파면의 기울기 dz/dy
        nvec = np.array([1.0, -dz])
        nvec /= np.linalg.norm(nvec)
        s = (Rr - z0) / nvec[0]
        a.plot([z0, z0 + s * nvec[0]], [y0, y0 + s * nvec[1]], color=C_RED if abs(y0) > 2.5 else C_KEY,
               lw=1.0, alpha=0.9)
    a.axvline(Rr, color="k", lw=0.8)
    a.plot([Rr], [0], "ko", ms=4)
    a.annotate("P*", (Rr, 0), (Rr + 0.25, 0.25), fontsize=10)
    i_w = 8
    a.annotate("", (act[i_w], yy[i_w]), (ref[i_w], yy[i_w]), arrowprops=dict(arrowstyle="<->", color=C_PUR, lw=1.2))
    a.text(act[i_w] + 0.15, yy[i_w] + 0.05, "W", fontsize=12, color=C_PUR, va="bottom")
    a.text(Rr + 0.15, -2.9, t(L, "가우스 상면", "Gaussian\nimage plane"), fontsize=9, **font(L))
    a.text(Rr - 1.4, -1.25, t(L, "광선 = 파면의 법선", "rays = normals\nof the wavefront"), fontsize=9,
           color=C_KEY, ha="right", **font(L))
    a.set_xlim(-0.3, Rr + 1.6)
    a.set_ylim(-3.3, 3.5)
    a.set_xlabel(t(L, "광축 방향 (임의 단위, 과장)", "along the axis (arbitrary, exaggerated)"), **font(L))
    a.set_yticks([])
    a.legend(loc="upper right", fontsize=8.5, prop=dict(family=plt.rcParams["font.family"], size=8.5))
    a.set_title(t(L, "(a) 파면이 기준 구면에서 벗어나면 광선이 상점을 비켜 간다",
                  "(a) A wavefront off the reference sphere sends rays past P*"), fontsize=10.5, **font(L))
    # (b) f/5 평볼록 렌즈: 실광선 W(ρ) 와 횡광선수차, 그리고 −(1/n′NA) dW/dρ
    lens = plano_convex(D=20.0)
    st = setup(lens, WL, 0.0)
    rho = np.linspace(-1, 1, 201)
    W, _ = opd(lens, WL, st, np.zeros_like(rho), rho)
    x, y, _ = spot(lens, WL, st, np.zeros_like(rho), rho, st.z_img)
    dW = np.gradient(W, rho)
    eps_pred = -dW / (st.n_img * st.sin_u_img)
    a2 = b.twinx()
    b.plot(rho, W / LAM_MM, color=C_RAY, lw=2, label=t(L, "파면수차 W (실광선)", "wave aberration W (real rays)"))
    a2.plot(rho, y * 1e3, color=C_RED, lw=2, label=t(L, "횡광선수차 (실광선)", "transverse ray error (real rays)"))
    a2.plot(rho[::12], eps_pred[::12] * 1e3, "o", mfc="none", color="k", ms=5,
            label=t(L, r"$-\frac{1}{n'\,\mathrm{NA}}\,\frac{dW}{d\rho}$ (파면에서 계산)", r"$-\frac{1}{n'\,\mathrm{NA}}\,\frac{dW}{d\rho}$ (from W)"))
    b.set_xlabel(t(L, r"동공 위 정규화 높이 $\rho$", r"normalized pupil height $\rho$"), **font(L))
    b.set_ylabel(t(L, r"W (파장 $\lambda$)", r"W (waves)"), color=C_RAY, **font(L))
    a2.set_ylabel(t(L, "가우스 상면에서 빗나간 높이 (µm)", "miss distance at Gaussian plane (µm)"), color=C_RED, **font(L))
    h1, l1 = b.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, loc="upper center", fontsize=8.5, prop=dict(family=plt.rcParams["font.family"], size=8.5))
    b.set_ylim(-0.5, 6.5)
    a2.set_ylim(-75, 75)
    b.set_title(t(L, "(b) f = 100 mm, F/5 평볼록 렌즈의 축상 물점", "(b) On-axis point, f = 100 mm F/5 plano-convex"),
                fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig2-wave-and-ray-aberration.png")


# ------------------------------------------------------------------ 그림 3: 다섯 가지 자이델 수차

SEIDEL_TERMS = [  # (이름 ko, 이름 en, W(x, y), ∂W/∂x, ∂W/∂y) — 시야는 +y, h = 1, 계수 = 1 파장
    ("구면수차", "spherical", lambda x, y: (x**2 + y**2) ** 2,
     lambda x, y: 4 * x * (x**2 + y**2), lambda x, y: 4 * y * (x**2 + y**2)),
    ("코마", "coma", lambda x, y: y * (x**2 + y**2),
     lambda x, y: 2 * x * y, lambda x, y: x**2 + 3 * y**2),
    ("비점수차", "astigmatism", lambda x, y: y**2,
     lambda x, y: 0 * x, lambda x, y: 2 * y),
    ("상면만곡", "field curvature", lambda x, y: x**2 + y**2,
     lambda x, y: 2 * x, lambda x, y: 2 * y),
    ("왜곡", "distortion", lambda x, y: y,
     lambda x, y: 0 * x, lambda x, y: 1 + 0 * y),
]
SEIDEL_FORMS = [r"$W_{040}\,\rho^4$", r"$W_{131}\,h\,\rho^3\cos\theta$", r"$W_{222}\,h^2\rho^2\cos^2\theta$",
                r"$W_{220}\,h^2\rho^2$", r"$W_{311}\,h^3\rho\cos\theta$"]


def figure3(L):
    n, pad = 127, 4
    X, Y, M = pupil_grid(n)
    pix = psf_pixel(n, pad)
    half = 4.5
    k = int(round(half / pix))
    fig, axes = plt.subplots(3, 5, figsize=(13, 8.1))
    rng = np.random.default_rng(1)
    # 동공을 고르게 채운 광선 (점 퍼짐용)
    rr = np.sqrt(rng.random(6000))
    th = 2 * np.pi * rng.random(6000)
    sx, sy = rr * np.cos(th), rr * np.sin(th)
    for j, (ko, en, Wf, gx, gy) in enumerate(SEIDEL_TERMS):
        W = np.where(M, Wf(X, Y), np.nan)
        ax = axes[0, j]
        lim = np.nanmax(np.abs(W))
        ax.imshow(W, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), extent=[-1, 1, -1, 1], origin="lower")
        ax.contour(X, Y, W, levels=np.linspace(-1, 1, 9), colors="k", linewidths=0.5)
        ax.set_title(f"{t(L, ko, en)}\n{SEIDEL_FORMS[j]}", fontsize=10.5, **font(L))
        ax.set_xticks([]); ax.set_yticks([])
        # 광선: ε = −∇W / NA (λ/NA 단위)
        ex, ey = -gx(sx, sy), -gy(sx, sy)
        ax = axes[1, j]
        ax.plot(ex, ey, ".", ms=0.8 if np.ptp(ey) > 0.1 else 6, color=C_RAY, alpha=0.6)
        ax.plot([0], [0], "+", color=C_RED, ms=9, mew=1.2)
        ax.set_xlim(-half, half); ax.set_ylim(-half, half); ax.set_aspect("equal")
        ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4]); ax.tick_params(labelsize=8)
        # 회절상
        I = psf(np.where(M, Wf(X, Y), 0.0), M, pad)
        c = I.shape[0] // 2
        ax = axes[2, j]
        ax.imshow(np.sqrt(I[c - k:c + k + 1, c - k:c + k + 1]), cmap="magma", origin="lower",
                  extent=[-k * pix, k * pix, -k * pix, k * pix])
        ax.plot([0], [0], "+", color="w", ms=9, mew=1.0)
        ax.set_xticks([-4, 0, 4]); ax.set_yticks([-4, 0, 4]); ax.tick_params(labelsize=8)
        ax.text(0.03, 0.04, t(L, f"정점 {I.max():.2f}", f"peak {I.max():.2f}"), color="w", fontsize=8.5,
                transform=ax.transAxes, **font(L))
    axes[0, 0].set_ylabel(t(L, "출사동의 파면수차", "wavefront in exit pupil"), fontsize=10, **font(L))
    axes[1, 0].set_ylabel(t(L, r"광선 도착점 ($\lambda$/NA)", r"ray hits ($\lambda$/NA)"), fontsize=10, **font(L))
    axes[2, 0].set_ylabel(t(L, r"회절상 ($\lambda$/NA)", r"diffraction image ($\lambda$/NA)"), fontsize=10, **font(L))
    fig.suptitle(t(L, "계수를 모두 1파장으로 둔 다섯 가지 자이델 수차 (가우스 상면, 시야는 위쪽, + 는 근축 상점)",
                   "The five Seidel aberrations, each 1 wave (Gaussian plane, field point upward, + = paraxial image)"),
                 fontsize=11.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig3-seidel-gallery.png")


# ------------------------------------------------------------------ 그림 4: 구면수차의 초점 부근

def figure4(L):
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6), gridspec_kw=dict(width_ratios=[1.25, 1]))
    lens = plano_convex(D=40.0)
    st = setup(lens, WL, 0.0)
    hs = np.linspace(-20, 20, 41)
    z_cross = axis_crossing(lens, WL, np.abs(hs[hs != 0]))
    zp = st.z_img
    zm = axis_crossing(lens, WL, [20.0])[0]
    lsa = zp - zm
    P = np.stack([np.zeros_like(hs), hs, np.full_like(hs, -10.0)], 1)
    D = np.tile([0.0, 0.0, 1.0], (len(hs), 1))
    Pf, Df, _, _, _ = trace_real(lens, P, D, WL)
    z0, z1 = zm - 1.2, zp + 0.8
    for p, d, h in zip(Pf, Df, hs):
        q0, q1 = to_plane(p[None], d[None], z0)[0], to_plane(p[None], d[None], z1)[0]
        col = plt.cm.viridis(abs(h) / 20)
        a.plot([q0[2] - zp, q1[2] - zp], [q0[1], q1[1]], color=col, lw=0.9)
    zs = np.linspace(z0, z1, 400)
    ys = []
    for z in zs:
        ys.append(np.max(np.abs(to_plane(Pf, Df, z)[:, 1])))
    ys = np.array(ys)
    z_lc = zs[np.argmin(ys)]
    marks = [(0.0, t(L, "근축 초점", "paraxial focus"), C_GREY),
             (zm - zp, t(L, "주변광선 초점", "marginal focus"), C_RED),
             (z_lc - zp, t(L, "최소 착란원", "least confusion"), C_KEY)]
    for zz, lab, col in marks:
        a.axvline(zz, color=col, lw=1.2, ls="--")
        a.text(zz - 0.04, -0.43, lab, rotation=90, ha="right", va="bottom", fontsize=9, color=col,
               bbox=dict(fc="white", ec="none", pad=1.0), **font(L))
    a.set_xlim(z0 - zp, z1 - zp)
    a.set_ylim(-0.45, 0.45)
    a.set_xlabel(t(L, "근축 초점에서 잰 위치 (mm)", "position from paraxial focus (mm)"), **font(L))
    a.set_ylabel(t(L, "높이 (mm)", "height (mm)"), **font(L))
    a.set_title(t(L, f"(a) F/2.5 평볼록 렌즈: 종구면수차 {lsa:.2f} mm", f"(a) F/2.5 plano-convex: LSA {lsa:.2f} mm"),
                fontsize=10.5, **font(L))
    print(f"    그림4: LSA {lsa:.4f} mm, 최소 착란원 {(z_lc - zm) / lsa:.3f} (주변 초점에서 근축 쪽으로, LSA 비)")
    # (b) 순수 3차 구면수차 W040 = 0.5 λ: 초점 위치에 따른 기하 번짐 반지름과 스트렐 비
    n = 201
    X, Y, M = pupil_grid(n)
    r2 = X**2 + Y**2
    w040 = 0.5
    frac = np.linspace(-0.25, 1.25, 121)      # 0: 근축 초점, 1: 주변광선 초점 (W020 = −2 W040 · frac)
    strehl, blur = [], []
    rho = np.linspace(0, 1, 2001)
    for f_ in frac:
        w020 = -2 * w040 * f_
        strehl.append(strehl_center(w040 * r2**2 + w020 * r2, M))
        blur.append(np.max(np.abs(4 * w040 * rho**3 + 2 * w020 * rho)))
    strehl, blur = np.array(strehl), np.array(blur)
    b2 = b.twinx()
    b.plot(frac, blur, color=C_KEY, lw=2, label=t(L, "기하광학: 번짐 반지름", "geometric blur radius"))
    b2.plot(frac, strehl, color=C_RAY, lw=2, label=t(L, "회절: 중심 세기 (스트렐 비)", "diffraction: central intensity (Strehl)"))
    b.axvline(0.75, color=C_KEY, ls=":", lw=1.2)
    b2.axvline(frac[np.argmax(strehl)], color=C_RAY, ls=":", lw=1.2)
    b.set_xlabel(t(L, "초점 위치 (0 = 근축 초점, 1 = 주변광선 초점)", "focus position (0 = paraxial, 1 = marginal)"), **font(L))
    b.set_ylabel(t(L, r"번짐 반지름 ($\lambda$/NA)", r"blur radius ($\lambda$/NA)"), color=C_KEY, **font(L))
    b2.set_ylabel(t(L, "스트렐 비", "Strehl ratio"), color=C_RAY, **font(L))
    b.set_ylim(0, 2.6)
    b2.set_ylim(0, 1.05)
    h1, l1 = b.get_legend_handles_labels()
    h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8.5, prop=dict(family=plt.rcParams["font.family"], size=8.5))
    b.text(0.77, 0.12, "3/4", color=C_KEY, fontsize=10)
    b.text(frac[np.argmax(strehl)] - 0.13, 0.12, "1/2", color=C_RAY, fontsize=10)
    b.set_title(t(L, r"(b) 순수 3차 구면수차 $W_{040} = 0.5\,\lambda$", r"(b) Pure third-order spherical, $W_{040} = 0.5\,\lambda$"),
                fontsize=10.5, **font(L))
    print(f"    그림4: 스트렐 최대 위치 {frac[np.argmax(strehl)]:.3f}, 최대 {strehl.max():.3f}, 번짐 최소 위치 {frac[np.argmin(blur)]:.3f}")
    fig.tight_layout()
    save(fig, L, "fig4-spherical-focus.png")


# ------------------------------------------------------------------ 그림 5: 렌즈 휘기

def lens_profile(ax, q, sx=0.011, h=22.0):
    """모양 인자 q 인 렌즈 단면(초점거리 100 mm, 반지름 h)을 x = q 자리에 그린다. 가로 1 mm = sx (q 단위)."""
    P = 1.0 / ((N_D - 1.0) * F)
    c1, c2 = P * (q + 1) / 2, P * (q - 1) / 2
    yy = np.linspace(-h, h, 80)
    s1 = lambda y: (1 / c1 - np.sign(c1) * np.sqrt(1 / c1**2 - y**2)) if abs(c1) > 1e-12 else 0 * y
    s2 = lambda y: (1 / c2 - np.sign(c2) * np.sqrt(1 / c2**2 - y**2)) if abs(c2) > 1e-12 else 0 * y
    t_c = max(2.0 - s2(h) + s1(h), 2.0)          # 가장자리 두께가 2 mm 이상이 되게
    z1, z2 = s1(yy), t_c + s2(yy)
    zc = 0.5 * (z1.mean() + z2.mean())
    zz = np.concatenate([z1, z2[::-1]]) - zc
    yv = np.concatenate([yy, yy[::-1]])
    ax.fill(q + zz * sx, yv, color="#cfe3f5", ec=C_RAY, lw=0.9)


def figure5(L):
    qs = np.linspace(-1.3, 2.3, 73)
    w040, w131, lsa = [], [], []
    for q in qs:
        lens = bent_lens(q)
        st = setup(lens, WL, np.radians(5.0))
        W = wave_coefficients(seidel_sums(lens, WL, st))
        w040.append(W["W040"] / LAM_MM)
        w131.append(W["W131"] / LAM_MM)
        z = axis_crossing(lens, WL, [1e-4, 10.0])
        lsa.append(z[0] - z[1])
    w040, w131, lsa = map(np.array, (w040, w131, lsa))
    fig = plt.figure(figsize=(9.5, 7.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.55, 1.2, 1.0], hspace=0.08)
    top = fig.add_subplot(gs[0])
    a = fig.add_subplot(gs[1], sharex=top)
    c = fig.add_subplot(gs[2], sharex=top)
    for q in [-1, 0, 1, 2]:
        lens_profile(top, q)
    top.text(-1.3, 0, t(L, "빛 →", "light →"), va="center", fontsize=9, color=C_GREY, **font(L))
    top.text(-1, -30, t(L, "평면이 앞", "flat side first"), ha="center", va="top", fontsize=8.5, **font(L))
    top.text(1, -30, t(L, "볼록면이 앞", "curved side first"), ha="center", va="top", fontsize=8.5, **font(L))
    top.text(0, -30, t(L, "대칭", "symmetric"), ha="center", va="top", fontsize=8.5, **font(L))
    top.set_ylim(-42, 26)
    top.axis("off")
    a.plot(qs, w040, color=C_RAY, lw=2.2, label=t(L, r"구면수차 $W_{040}$ (3차 이론, 왼쪽 축)", r"spherical $W_{040}$ (third order, left)"))
    a2 = a.twinx()
    a2.plot(qs[::4], lsa[::4], "s", mfc="none", color="k", ms=5,
            label=t(L, "종구면수차 (실광선, 오른쪽 축)", "longitudinal SA (real rays, right)"))
    q_min = qs[np.argmin(w040)]
    a.axvline(q_min, color=C_RAY, ls=":", lw=1)
    a.text(q_min + 0.04, w040.min() * 0.55, t(L, f"최소 q = {q_min:.2f}", f"minimum q = {q_min:.2f}"), color=C_RAY, fontsize=9, **font(L))
    a.set_ylim(0, w040.max() * 1.05)
    a2.set_ylim(0, w040.max() * 1.05 * lsa.min() / w040.min())
    a.set_ylabel(t(L, r"$W_{040}$ (파장)", r"$W_{040}$ (waves)"), **font(L))
    a2.set_ylabel(t(L, "종구면수차 (mm)", "longitudinal SA (mm)"), **font(L))
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, loc="upper center", fontsize=8.5, prop=dict(family=plt.rcParams["font.family"], size=8.5))
    plt.setp(a.get_xticklabels(), visible=False)
    plt.setp(top.get_xticklabels(), visible=False)
    c.plot(qs, w131, color=C_RED, lw=2.2, label=t(L, r"코마 $W_{131}$, 시야 5° (3차 이론)", r"coma $W_{131}$ at 5° (third order)"))
    c.axhline(0, color="k", lw=0.6)
    q_c0 = np.interp(0.0, w131, qs) if w131[-1] > w131[0] else np.interp(0.0, w131[::-1], qs[::-1])
    c.axvline(q_c0, color=C_RED, ls=":", lw=1)
    a.axvline(q_c0, color=C_RED, ls=":", lw=1)
    c.text(q_c0 + 0.04, w131.min() * 0.5, t(L, f"코마 0, q = {q_c0:.2f}", f"zero coma, q = {q_c0:.2f}"), color=C_RED, fontsize=9, **font(L))
    c.set_ylabel(t(L, r"$W_{131}$ (파장)", r"$W_{131}$ (waves)"), **font(L))
    c.set_xlabel(t(L, r"모양 인자 $q = (R_2 + R_1)/(R_2 - R_1)$", r"shape factor $q = (R_2 + R_1)/(R_2 - R_1)$"), **font(L))
    c.legend(loc="upper left", fontsize=8.5, prop=dict(family=plt.rcParams["font.family"], size=8.5))
    c.set_xlim(-1.35, 2.35)
    top.set_title(t(L, "초점거리 100 mm, F/5 인 N-BK7 렌즈를 여러 모양으로 휘었을 때 (조리개는 렌즈에)",
                    "An f = 100 mm, F/5 N-BK7 lens bent into different shapes (stop at the lens)"), fontsize=11, **font(L))
    i1, im1 = np.argmin(abs(qs - 1)), np.argmin(abs(qs + 1))
    print(f"    그림5: q_min {q_min:.3f}, q_coma0 {q_c0:.3f}, LSA(q=+1) {lsa[i1]:.4f}, LSA(q=-1) {lsa[im1]:.4f}, 비 {lsa[im1] / lsa[i1]:.2f}, "
          f"W040 비 {w040[im1] / w040[i1]:.2f}")
    save(fig, L, "fig5-lens-bending.png")


# ------------------------------------------------------------------ 그림 6: 커버글라스 두께 오차

def coverslip_strehl(NA, dt_mm, n=1.523, lam_mm=0.55e-3, m=4001):
    """두께 오차 dt 인 평판이 만드는 파면수차의 최선 초점 스트렐 비와 RMS (사인 조건을 만족하는 대물렌즈,
    동공 좌표 ρ = sinθ/NA 위에서 고른 세기). 초점 조절 항은 고NA 형태 δz·cosθ."""
    rho = np.linspace(0, 1, m)
    w = rho / np.sum(rho)
    s = NA * rho
    c, cp = np.sqrt(1 - s**2), np.sqrt(1 - (s / n) ** 2)
    W = dt_mm * (n * cp - c)
    A = np.vstack([np.ones_like(rho), c]).T * np.sqrt(w)[:, None]
    coef, *_ = np.linalg.lstsq(A, W * np.sqrt(w), rcond=None)
    res = W - coef[0] - coef[1] * c
    sig = np.sqrt(np.sum(res**2 * w))
    # 정확한 스트렐 비: 초점 조절량을 RMS 최소 근처에서 미세 탐색
    best = 0.0
    for dz in coef[1] + np.linspace(-3, 3, 1201) * lam_mm:
        r = W - dz * c
        best = max(best, abs(np.sum(w * np.exp(2j * np.pi * r / lam_mm))) ** 2)
    return sig / lam_mm, best


def figure6(L):
    NAs = np.linspace(0.1, 0.95, 69)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.3))
    for dt, col in [(0.005, C_GRN), (0.010, C_RAY), (0.020, C_RED)]:
        res = np.array([coverslip_strehl(na, dt) for na in NAs])
        a.semilogy(NAs, res[:, 0], color=col, lw=2, label=t(L, f"두께 오차 {dt * 1e3:.0f} µm", f"thickness error {dt * 1e3:.0f} µm"))
        b.plot(NAs, res[:, 1], color=col, lw=2, label=t(L, f"두께 오차 {dt * 1e3:.0f} µm", f"thickness error {dt * 1e3:.0f} µm"))
    a.axhline(1 / 14, color="k", ls="--", lw=1)
    a.text(0.12, 1 / 14 * 1.2, t(L, r"마레샬 기준 $\lambda/14$", r"Maréchal $\lambda/14$"), fontsize=9, **font(L))
    b.axhline(0.8, color="k", ls="--", lw=1)
    for ax in (a, b):
        ax.axvline(0.4, color=C_GREY, ls=":", lw=1)
        ax.set_xlabel(t(L, "건식 대물렌즈의 개구수 NA", "numerical aperture of a dry objective"), **font(L))
        ax.legend(fontsize=8.5, loc="lower left" if ax is b else "upper left",
                  prop=dict(family=plt.rcParams["font.family"], size=8.5))
    a.set_ylabel(t(L, r"최선 초점의 RMS 파면수차 ($\lambda$)", r"RMS wave aberration at best focus (waves)"), **font(L))
    b.set_ylabel(t(L, "스트렐 비", "Strehl ratio"), **font(L))
    a.set_ylim(1e-4, 2)
    a.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.set_ylim(0, 1.03)
    a.set_title(t(L, "(a) 0.17 mm 커버글라스가 조금 두껍거나 얇으면", "(a) A 0.17 mm coverslip slightly off in thickness"), fontsize=10.5, **font(L))
    b.set_title(t(L, "(b) 상의 중심 세기", "(b) Peak intensity of the image"), fontsize=10.5, **font(L))
    for na in [0.4, 0.65, 0.75, 0.85, 0.95]:
        s_, S_ = coverslip_strehl(na, 0.010)
        print(f"    그림6: NA {na}: 10 µm → RMS {s_:.4f} λ, 스트렐 {S_:.3f}")
    fig.tight_layout()
    save(fig, L, "fig6-coverslip.png")


# ------------------------------------------------------------------ 그림 8: 코마

def paraboloid_spot(f, D, theta, rings=(0.2, 0.4, 0.6, 0.8, 1.0), m=240):
    """초점거리 f 의 포물면 거울에 각 θ 로 들어온 평면파를 실광선 추적해 초점면 z = f 위 도착점을 돌려준다.
    거울: z = r²/(4f), 빛은 +z 쪽에서 −z 방향으로 들어온다. 조리개는 거울 위."""
    out = []
    d = np.array([0.0, -np.sin(theta), -np.cos(theta)])
    for rho in rings:
        th = np.linspace(0, 2 * np.pi, m, endpoint=False)
        x, y = rho * D / 2 * np.cos(th), rho * D / 2 * np.sin(th)
        P = np.stack([x, y, (x**2 + y**2) / (4 * f)], 1)        # 거울 위의 동공점에서 위로 물러난 곳에서 출발
        P = P - 3 * f * d[None]
        a = d[0] ** 2 + d[1] ** 2
        b = 2 * (P[:, 0] * d[0] + P[:, 1] * d[1]) - 4 * f * d[2]
        c = P[:, 0] ** 2 + P[:, 1] ** 2 - 4 * f * P[:, 2]
        qq = -0.5 * (b + np.sign(b) * np.sqrt(b * b - 4 * a * c))     # 수치적으로 안정한 근의 공식
        roots = np.stack([qq / a if a > 1e-15 else np.full_like(qq, np.inf), c / qq], 1)
        roots[roots <= 0] = np.inf
        tt = roots.min(1)                                           # 가장 가까운 교점
        Q = P + tt[:, None] * d
        N = np.stack([2 * Q[:, 0], 2 * Q[:, 1], -4 * f * np.ones(len(Q))], 1)
        N /= np.linalg.norm(N, axis=1)[:, None]
        R = d - 2 * np.sum(d * N, 1)[:, None] * N
        s = (f - Q[:, 2]) / R[:, 2]
        Z = Q + s[:, None] * R
        out.append((rho, th, Z[:, 0], Z[:, 1]))
    return out


def figure8(L):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 5.0))
    # (a) 순수 3차 코마의 광선 도착점 (W131 = 1 λ, λ/NA 단위): 동공의 같은 고리는 상면에서 원 하나
    for rho in [0.2, 0.4, 0.6, 0.8, 1.0]:
        th = np.linspace(0, 2 * np.pi, 400)
        x, y = rho * np.sin(th), rho * np.cos(th)     # θ 는 +y 에서 잰다
        ex, ey = -(2 * x * y), -(x**2 + 3 * y**2)
        col = plt.cm.viridis(rho)
        a.plot(ex, ey, color=col, lw=1.6)
        a.plot([0], [0], "+", color=C_RED, ms=10, mew=1.5)
        idx = [0, 50, 100]
        a.plot(ex[idx], ey[idx], "o", color=col, ms=3.5)
    yy = np.linspace(0, -3.3, 10)
    a.plot(yy * np.tan(np.radians(30)), yy, "k--", lw=0.8)
    a.plot(-yy * np.tan(np.radians(30)), yy, "k--", lw=0.8)
    a.text(0, -3.55, r"60$^\circ$", ha="center", fontsize=10)
    a.set_aspect("equal")
    a.set_xlim(-2.2, 2.2); a.set_ylim(-3.8, 0.6)
    a.set_xlabel(r"x ($\lambda$/NA)"); a.set_ylabel(r"y ($\lambda$/NA)")
    a.set_title(t(L, r"(a) 순수 3차 코마 $W_{131} = 1\,\lambda$" "\n동공의 고리 하나 → 상면의 원 하나 (두 바퀴)",
                  r"(a) Pure third-order coma, $W_{131} = 1\,\lambda$" "\none pupil ring → one circle, traced twice"),
                fontsize=10.5, **font(L))
    # (b) 포물면 거울 (f = 400 mm, F/4) 의 시야 0.5° 실광선 점 퍼짐
    f_m, D_m, th_deg = 300.0, 100.0, 0.2
    spots = paraboloid_spot(f_m, D_m, np.radians(th_deg))
    yc = spots[0][3].mean()
    chief_y = -f_m * np.tan(np.radians(th_deg))
    for rho, th, x, y in spots:
        b.plot(x * 1e3, (y - chief_y) * 1e3, ".", color=plt.cm.viridis(rho), ms=1.6)
    b.plot([0], [0], "+", color=C_RED, ms=10, mew=1.5)
    b.set_aspect("equal")
    b.set_xlabel(t(L, "x (µm)", "x (µm)")); b.set_ylabel(t(L, "y, 주광선 기준 (µm)", "y from the chief ray (µm)"), **font(L))
    b.set_title(t(L, f"(b) 포물면 거울 f = {f_m:.0f} mm, F/{f_m / D_m:.0f}, 시야 {th_deg}°\n실광선 추적 (색 = 동공 고리 반지름)",
                  f"(b) Paraboloid f = {f_m:.0f} mm, F/{f_m / D_m:.0f}, field {th_deg}°\nreal rays (color = pupil ring radius)"),
                fontsize=10.5, **font(L))
    ys = np.concatenate([s[3] for s in spots]) - chief_y
    print(f"    그림8: 포물면 코마 꼬리 길이 {np.ptp(ys) * 1e3:.2f} µm (3차식 3θf/(16N²) = {3 * np.radians(th_deg) * f_m / (16 * (f_m / D_m) ** 2) * 1e3:.2f} µm)")
    fig.tight_layout()
    save(fig, L, "fig8-coma.png")


# ------------------------------------------------------------------ 그림 9: 비점수차와 상면만곡

def focal_tangential_sagittal(lens, field, dp=0.02):
    """주광선 근처의 가는 자오·구결 광선 다발이 주광선과 만나는 z (상 공간)."""
    st = setup(lens, WL, field)
    out = []
    for px, py in [(0.0, dp), (dp, 0.0)]:
        P, D, o = launch_rays(st, [0.0, px, -px], [0.0, py, -py])
        P, D, n, opl, ok = trace_real(lens, P, D, WL)
        # 주광선과 이웃 광선의 교점: 주광선에 수직한 방향 성분이 0 이 되는 z
        c0, c1 = P[0], P[1]
        d0, d1 = D[0], D[1]
        # 자오: y–z 평면, 구결: x–z 평면에서의 교점
        if py:
            s = ((c1[1] - c0[1]) - (c1[2] - c0[2]) * d1[1] / d1[2]) / (d0[1] - d0[2] * d1[1] / d1[2])
        else:
            # 구결 광선은 x 방향으로 벌어져 있다: x 가 0 이 되는 곳(주광선의 x = 0)
            s_ = -c1[0] / d1[0]
            q = c1 + s_ * d1
            out.append(q[2])
            continue
        q = c0 + s * d0
        out.append(q[2])
    return out, st


def launch_rays(st, px, py):
    from aberr import launch
    return launch(st, np.asarray(px, float), np.asarray(py, float))


def figure9(L):
    lens = plano_convex(D=10.0)
    field = np.radians(10.0)
    st = setup(lens, WL, field)
    (zt, zs), _ = focal_tangential_sagittal(lens, field)
    fig = plt.figure(figsize=(12.5, 7.6))
    outer = fig.add_gridspec(2, 1, height_ratios=[1, 1.25], hspace=0.32)
    gs = outer[0].subgridspec(1, 5, wspace=0.38)
    gs_b = outer[1].subgridspec(1, 2, wspace=0.25)
    rng = np.random.default_rng(3)
    rr = np.sqrt(rng.random(3000)); th = 2 * np.pi * rng.random(3000)
    px, py = rr * np.cos(th), rr * np.sin(th)
    planes = [zt, zt + (zs - zt) * 0.25, (zt + zs) / 2, zt + (zs - zt) * 0.75, zs]
    names = [t(L, "자오 초점선", "tangential focal line"), "", t(L, "최소 착란원", "least confusion"), "",
             t(L, "구결 초점선", "sagittal focal line")]
    for k, (zp_, nm) in enumerate(zip(planes, names)):
        ax = fig.add_subplot(gs[k])
        x, y, ok = spot(lens, WL, st, px, py, zp_)
        xc, yc, _ = spot(lens, WL, st, np.array([0.0]), np.array([0.0]), zp_)
        ax.plot((x - xc) * 1e3, (y - yc) * 1e3, ".", ms=0.8, color=C_RAY, alpha=0.6)
        ax.set_aspect("equal"); ax.set_xlim(-160, 160); ax.set_ylim(-160, 160)
        ax.tick_params(labelsize=8)
        ax.set_title(f"{nm}\nz = {zp_ - zt:+.2f} mm" if nm else f"z = {zp_ - zt:+.2f} mm", fontsize=9.5, **font(L))
        if k == 0:
            ax.set_ylabel(t(L, "y (µm), 위 = 광축에서 먼 쪽", "y (µm), up = away from axis"), fontsize=9, **font(L))
    print(f"    그림9: 10° 에서 자오 초점 {zt:.3f}, 구결 초점 {zs:.3f}, 차 {zs - zt:.3f} mm")
    # (b) 시야에 따른 T·S·페츠발 상면 (렌즈 앞 조리개) / (c) 조리개를 옮겨 비점수차를 없앤 경우
    for col, (lens_, title) in enumerate([
            (plano_convex(D=10.0), t(L, "(b) 조리개가 렌즈에 붙어 있을 때", "(b) Stop at the lens")),
            (plano_convex(D=10.0, stop_z=center_stop_z(), flip=True),
             t(L, f"(c) 렌즈를 뒤집고 조리개를 {-center_stop_z():.1f} mm 앞(곡률 중심의 겉보기 자리)에",
               f"(c) Lens reversed, stop {-center_stop_z():.1f} mm in front (apparent center of curvature)"))]):
        ax = fig.add_subplot(gs_b[col])
        z_par = setup(lens_, WL, 0.0).z_img
        degs = np.linspace(0.5, 15, 30)
        T, Sg, Pz, Tr, Sr = [], [], [], [], []
        for d in degs:
            st_ = setup(lens_, WL, np.radians(d))
            S = seidel_sums(lens_, WL, st_)
            # 3차 이론: 상면이 근축 상면에서 벗어난 거리 Δz = −(2 W / (n′ u′²)) (W: 그 시야의 초점 항)
            u2 = st_.n_img * st_.sin_u_img**2
            Wc = wave_coefficients(S)
            T.append(-2 * (Wc["W220"] + Wc["W222"]) / u2)
            Sg.append(-2 * Wc["W220"] / u2)
            Pz.append(-2 * Wc["W220P"] / u2)
            (zt_, zs_), _ = focal_tangential_sagittal(lens_, np.radians(d))
            Tr.append(zt_ - z_par); Sr.append(zs_ - z_par)
        ax.plot(T, degs, color=C_RED, lw=2, label=t(L, "자오(T) 3차", "tangential (T), 3rd order"))
        ax.plot(Sg, degs, color=C_RAY, lw=2, label=t(L, "구결(S) 3차", "sagittal (S), 3rd order"))
        ax.plot(Pz, degs, color="k", lw=1.4, ls="--", label=t(L, "페츠발 면", "Petzval surface"))
        ax.plot(Tr[::3], degs[::3], "o", mfc="none", color=C_RED, ms=5, label=t(L, "T 실광선", "T real rays"))
        ax.plot(Sr[::3], degs[::3], "s", mfc="none", color=C_RAY, ms=5, label=t(L, "S 실광선", "S real rays"))
        ax.axvline(0, color=C_GREY, lw=0.8)
        ax.set_xlabel(t(L, "근축 상면에서 잰 초점 위치 (mm)", "focus from paraxial image plane (mm)"), **font(L))
        ax.set_ylabel(t(L, "시야각 (°)", "field angle (°)"), **font(L))
        ax.set_title(title, fontsize=10.5, **font(L))
        ax.set_xlim(-7.5, 1.0)
        ax.legend(fontsize=8, loc="lower left", prop=dict(family=plt.rcParams["font.family"], size=8))
        i = np.argmin(abs(degs - 10))
        print(f"    그림9{'bc'[col]}: 10° (T−P)/(S−P) = {(T[i] - Pz[i]) / (Sg[i] - Pz[i]) if abs(Sg[i] - Pz[i]) > 1e-9 else float('nan'):.3f}, "
              f"P {Pz[i]:.3f}, 실광선 T {Tr[i]:.3f} S {Sr[i]:.3f}")
    fig.suptitle(t(L, "(a) f = 100 mm F/10 평볼록 렌즈, 시야 10°: 초점을 지나며 바뀌는 점의 모양",
                   "(a) f = 100 mm F/10 plano-convex at 10°: the spot through focus"), fontsize=11, y=0.995, **font(L))
    fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.07)
    save(fig, L, "fig9-astigmatism-field-curvature.png")


# ------------------------------------------------------------------ 그림 11: 왜곡

def figure11(L):
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.5))
    g = np.linspace(-1, 1, 9)
    tmax = np.tan(np.radians(20.0))
    cases = [(-25.0, t(L, "(a) 조리개가 렌즈 25 mm 앞", "(a) Stop 25 mm in front")),
             (0.0, t(L, "(b) 조리개가 렌즈 한가운데", "(b) Stop at the lens")),
             (25.0, t(L, "(c) 조리개가 렌즈 25 mm 뒤", "(c) Stop 25 mm behind"))]
    for ax, (zs_, title) in zip(axes, cases):
        lens = equiconvex(zs_)
        z_img = setup(lens, WL, 0.0).z_img
        fpar = None
        for orient in range(2):
            for gv in g:
                pts_real, pts_par = [], []
                for gu in np.linspace(-1, 1, 41):
                    tx, ty = (gu, gv) if orient == 0 else (gv, gu)
                    tx, ty = tx * tmax, ty * tmax
                    th = np.arctan(np.hypot(tx, ty))
                    if th < 1e-12:
                        pts_real.append((0, 0)); pts_par.append((0, 0)); continue
                    st = setup(lens, WL, th)
                    x, y, ok = spot(lens, WL, st, np.array([0.0]), np.array([0.0]), z_img)
                    yp = paraxial_image_height(lens, WL, st)
                    phi = np.arctan2(tx, ty)
                    r_real, r_par = y[0], yp
                    pts_real.append((r_real * np.sin(phi), r_real * np.cos(phi)))
                    pts_par.append((r_par * np.sin(phi), r_par * np.cos(phi)))
                pr, pp = np.array(pts_real), np.array(pts_par)
                ax.plot(pp[:, 0], pp[:, 1], color=C_GREY, lw=0.7, ls="--")
                ax.plot(pr[:, 0], pr[:, 1], color=C_RAY, lw=1.3)
        st = setup(lens, WL, np.radians(20.0))
        x, y, ok = spot(lens, WL, st, np.array([0.0]), np.array([0.0]), z_img)
        yp = paraxial_image_height(lens, WL, st)
        ax.set_aspect("equal")
        ax.set_xlim(-60, 60); ax.set_ylim(-60, 60)
        dist = 100 * (y[0] - yp) / yp
        dist_s = "0.0" if abs(dist) < 0.05 else f"{dist:+.1f}"
        ax.set_title(f"{title}\n" + t(L, f"20° 에서 {dist_s}%", f"{dist_s}% at 20°"),
                     fontsize=10.5, **font(L))
        ax.set_xlabel("x (mm)")
        print(f"    그림11: 조리개 {zs_:+.0f} mm → 20° 왜곡 {100 * (y[0] - yp) / yp:+.2f}%")
    axes[0].set_ylabel("y (mm)")
    fig.suptitle(t(L, "f = 100 mm 양볼록 렌즈로 본 바둑판 (파란 선: 실광선, 회색 점선: 근축 상)",
                   "A grid seen through an f = 100 mm biconvex lens (blue: real rays, grey dashed: paraxial)"), fontsize=11, **font(L))
    fig.tight_layout()
    save(fig, L, "fig11-distortion.png")


# ------------------------------------------------------------------ 그림 13: 색수차

def back_focus(lens, wl):
    m = paraxial_trace(lens, 1.0, 0.0, lens[0].z, wl)
    return lens[-1].z - m[-1, 0] / m[-1, 2]


def figure13(L):
    wls = np.linspace(0.40, 0.75, 141)
    single = singlet_for_color()
    doub, info = achromat()
    fs = np.array([back_focus(single, w) for w in wls]) - back_focus(single, LINE_D)
    fd = np.array([back_focus(doub, w) for w in wls]) - back_focus(doub, LINE_D)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.4))
    a.plot(wls * 1e3, fs, color=C_RED, lw=2.2, label=t(L, "N-BK7 단렌즈", "N-BK7 singlet"))
    a.plot(wls * 1e3, fd, color=C_RAY, lw=2.2, label=t(L, "N-BK7 + N-F2 색지움 렌즈", "N-BK7 + N-F2 achromat"))
    b.plot(wls * 1e3, fd * 1e3, color=C_RAY, lw=2.2)
    for ax, sc in ((a, 1.0), (b, 1e3)):
        for lam, nm, col in [(LINE_F, "F", "#3060ff"), (LINE_D, "d", "#c8a000"), (LINE_C, "C", "#e03030")]:
            ax.axvline(lam * 1e3, color=col, lw=0.9, ls=":")
            ax.text(lam * 1e3 + 3, 0, nm, color=col, fontsize=10, va="bottom")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_xlabel(t(L, "파장 (nm)", "wavelength (nm)"), **font(L))
    a.set_ylabel(t(L, "d 선(587.6 nm) 초점에서 잰 초점 위치 (mm)", "focus shift from the d-line focus (mm)"), **font(L))
    b.set_ylabel(t(L, "초점 위치 (µm)", "focus shift (µm)"), **font(L))
    a.legend(fontsize=9, prop=dict(family=plt.rcParams["font.family"], size=9))
    a.set_title(t(L, "(a) f = 100 mm 두 렌즈의 색에 따른 초점", "(a) Focus vs color, two f = 100 mm lenses"), fontsize=10.5, **font(L))
    b.set_title(t(L, "(b) 색지움 렌즈만 확대: 남는 2차 스펙트럼", "(b) Achromat only: the secondary spectrum"), fontsize=10.5, **font(L))
    print(f"    그림13: 단렌즈 F−C {back_focus(single, LINE_F) - back_focus(single, LINE_C):.4f} mm, "
          f"색지움 F−C {1e3 * (back_focus(doub, LINE_F) - back_focus(doub, LINE_C)):.2f} µm, "
          f"F−d {1e3 * (back_focus(doub, LINE_F) - back_focus(doub, LINE_D)):.2f} µm, 최소 파장 {wls[np.argmin(fd)] * 1e3:.0f} nm, "
          f"f1 {info[0]:.3f} f2 {info[1]:.3f} R {info[2]:.3f} {info[3]:.3f} R3 얇은 {info[4]:.2f} → {info[5]:.2f}")
    fig.tight_layout()
    save(fig, L, "fig13-chromatic-focus.png")


# ------------------------------------------------------------------ 그림 14: 쿡 삼중렌즈의 자이델 장부

def draw_lens_system(ax, surfs, wl, field, hmax=None):
    st0 = setup(surfs, wl, 0.0)
    for i in range(len(surfs) - 1):
        s1, s2 = surfs[i], surfs[i + 1]
        if s1.medium(wl) == 1.0:
            continue
        h = min(s.semi for s in (s1, s2) if not s.stop) if hmax is None else hmax
        yy = np.linspace(-h, h, 80)
        z1 = s1.z + (0 if np.isinf(s1.R) else s1.R - np.sign(s1.R) * np.sqrt(s1.R**2 - yy**2))
        z2 = s2.z + (0 if np.isinf(s2.R) else s2.R - np.sign(s2.R) * np.sqrt(s2.R**2 - yy**2))
        ax.fill(np.concatenate([z1, z2[::-1]]), np.concatenate([yy, yy[::-1]]), color="#cfe3f5", ec=C_RAY, lw=0.8)
    stp = next(s for s in surfs if s.stop)
    for sgn in (1, -1):
        ax.plot([stp.z, stp.z], [sgn * stp.semi, sgn * (stp.semi + 3)], color="k", lw=2)


def figure14(L):
    tri = cooke_triplet()
    fld = np.radians(20.0)
    st = setup(tri, WL, fld)
    S = seidel_surface_terms(tri, WL, st)
    fig = plt.figure(figsize=(12.5, 4.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.5])
    a = fig.add_subplot(gs[0])
    draw_lens_system(a, tri, WL, fld)
    for fieldv, col in [(0.0, C_KEY), (fld, C_GRN)]:
        st_ = setup(tri, WL, fieldv)
        for py in ([-1, 0, 1] if fieldv else [-1, 1]):
            from aberr import launch
            P, D, o = launch(st_, np.array([0.0]), np.array([float(py)]), back=8.0)
            pts = [P[0]]
            for k in range(1, len(tri) + 1):          # 면마다의 교점: 앞쪽 k 개 면까지만 추적
                Pk, Dk, _, _, _ = trace_real(tri[:k], P, D, WL)
                pts.append(Pk[0])
            q = to_plane(Pk, Dk, st_.z_img)[0]
            pts.append(q)
            pts = np.array(pts)
            a.plot(pts[:, 2], pts[:, 1], color=col, lw=1.0)
    a.axhline(0, color=C_GREY, lw=0.6)
    a.set_aspect("equal")
    a.set_xlim(-9, st.z_img + 2)
    a.set_ylim(-14, 26)
    a.set_xlabel("z (mm)")
    a.set_title(t(L, "(a) 쿡 삼중렌즈 f = 50 mm, F/5\n주황: 축상 물점, 초록: 시야 20°", "(a) Cooke triplet f = 50 mm, F/5\norange: on axis, green: 20° field"),
                fontsize=10.5, **font(L))
    for k, s in enumerate(tri):
        a.text(s.z + (-0.8 if k == 2 else 0.8 if k == 3 else 0), -12.5, str(k + 1), ha="center", fontsize=9)
    b = fig.add_subplot(gs[1])
    names = ["$S_I$", "$S_{II}$", "$S_{III}$", "$S_{IV}$", "$S_V$"]
    labels_ko = ["구면", "코마", "비점", "페츠발", "왜곡"]
    labels_en = ["spherical", "coma", "astig.", "Petzval", "distortion"]
    nsurf = S.shape[1]
    width = 0.8 / (nsurf + 1)
    xs = np.arange(5)
    cmap = plt.cm.tab10
    for k in range(nsurf):
        b.bar(xs - 0.4 + width * (k + 0.5), S[:, k] * 1e3, width,
              color=cmap(k), label=t(L, f"{k + 1}면", f"surface {k + 1}"))
    b.bar(xs - 0.4 + width * (nsurf + 0.5), S.sum(1) * 1e3, width, color="k", label=t(L, "합", "sum"))
    b.axhline(0, color="k", lw=0.6)
    b.set_xticks(xs)
    b.set_xticklabels([f"{n_}\n{t(L, k_, e_)}" for n_, k_, e_ in zip(names, labels_ko, labels_en)], **font(L))
    b.set_ylabel(t(L, r"면마다의 기여 ($\times 10^{-3}$ mm)", r"contribution per surface ($\times 10^{-3}$ mm)"), **font(L))
    b.legend(ncol=4, fontsize=8, prop=dict(family=plt.rcParams["font.family"], size=8))
    b.set_title(t(L, "(b) 근축 두 광선으로 계산한 면별 자이델 기여와 그 합", "(b) Seidel contributions per surface from two paraxial rays, and the sums"),
                fontsize=10.5, **font(L))
    np.set_printoptions(precision=5, suppress=True)
    print("    그림14: 면별 S (×1e-3 mm)\n", S * 1e3, "\n    합", S.sum(1) * 1e3,
          "\n    합/면별 최대 절댓값", np.abs(S.sum(1)) / np.abs(S).max(1))
    fig.tight_layout()
    save(fig, L, "fig14-seidel-budget.png")


# ------------------------------------------------------------------ 그림 15: 얼마나 작아야 하나

def best_strehl(Wf, X, Y, M, amp):
    """주어진 모양의 수차(진폭 amp 파장)에서, 기울기와 초점을 조절해 얻는 최대 중심 세기와 그때의 RMS."""
    from aberr import rms
    W = amp * Wf(X, Y)
    # RMS 최소가 되는 기울기·초점을 빼면 작은 수차에서 최대 세기 위치와 같다
    w = W[M]
    A = np.stack([np.ones_like(w), X[M], Y[M], X[M] ** 2 + Y[M] ** 2], 1)
    coef, *_ = np.linalg.lstsq(A, w, rcond=None)
    Wb = W.copy()
    Wb[M] = w - (A * coef).sum(1)
    best = strehl_center(Wb, M)
    # 큰 수차에서는 RMS 최소 위치와 세기 최대 위치가 다를 수 있어 초점을 미세 탐색한다
    for d in np.linspace(-0.3, 0.3, 61) * max(amp, 0.1):
        best = max(best, strehl_center(Wb + d * (X**2 + Y**2), M))
    return best, float(np.sqrt(np.mean(Wb[M] ** 2)))


def figure15(L):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.8, 4.5))
    rho = np.linspace(0, 1, 400)
    pure = rho**4
    bal = rho**4 - rho**2
    zern = (6 * rho**4 - 6 * rho**2 + 1) / 6
    rms_ = lambda f: np.sqrt(np.trapezoid((f - np.trapezoid(f * 2 * rho, rho)) ** 2 * 2 * rho, rho))
    a.plot(rho, pure - np.trapezoid(pure * 2 * rho, rho), color=C_RED, lw=2,
           label=t(L, rf"$\rho^4$ 그대로 (RMS {rms_(pure):.3f})", rf"$\rho^4$ alone (RMS {rms_(pure):.3f})"))
    a.plot(rho, bal - np.trapezoid(bal * 2 * rho, rho), color=C_RAY, lw=2,
           label=t(L, rf"$\rho^4 - \rho^2$, 초점 조절 (RMS {rms_(bal):.3f})", rf"$\rho^4 - \rho^2$, refocused (RMS {rms_(bal):.3f})"))
    a.plot(rho, zern, "k:", lw=1.5, label=t(L, r"제르니케 $(6\rho^4 - 6\rho^2 + 1)/6$", r"Zernike $(6\rho^4 - 6\rho^2 + 1)/6$"))
    a.axhline(0, color=C_GREY, lw=0.6)
    a.set_xlabel(t(L, r"동공 반지름 $\rho$", r"pupil radius $\rho$"), **font(L))
    a.set_ylabel(t(L, r"파면 (평균을 뺌, $W_{040} = 1$)", r"wavefront (mean removed, $W_{040} = 1$)"), **font(L))
    a.legend(fontsize=8.5, loc="upper left", prop=dict(family=plt.rcParams["font.family"], size=8.5))
    a.set_title(t(L, "(a) 초점을 조금 옮기면 같은 구면수차의 RMS 가 1/4 로", "(a) Refocusing cuts the RMS of the same aberration to 1/4"),
                fontsize=10.5, **font(L))
    print(f"    그림15a: RMS 순수 {rms_(pure):.4f} (1/√{1 / rms_(pure) ** 2:.2f}), 균형 {rms_(bal):.4f} (=1/(6√5) {1 / (6 * np.sqrt(5)):.4f})")
    n = 201
    X, Y, M = pupil_grid(n)
    shapes = [(t(L, "구면수차", "spherical"), lambda x, y: (x**2 + y**2) ** 2, C_RAY, 0.94),
              (t(L, "코마", "coma"), lambda x, y: y * (x**2 + y**2), C_RED, 0.60),
              (t(L, "비점수차", "astigmatism"), lambda x, y: y**2, C_GRN, 0.35)]
    sig = np.linspace(0, 0.25, 200)
    b.plot(sig, np.exp(-(2 * np.pi * sig) ** 2), color="k", lw=1.2, ls="--",
           label=t(L, r"마레샬 근사 $e^{-(2\pi\sigma)^2}$", r"Maréchal $e^{-(2\pi\sigma)^2}$"))
    for name, Wf, col, tol in shapes:
        amps = np.linspace(0.02, 3.0, 40)
        res = np.array([best_strehl(Wf, X, Y, M, A) for A in amps])
        res = res[res[:, 1] <= 0.17]
        b.plot(res[:, 1], res[:, 0], color=col, lw=2, label=t(L, f"{name} (정확)", f"{name} (exact)"))
        s_tol, r_tol = best_strehl(Wf, X, Y, M, tol)
        b.plot([r_tol], [s_tol], "o", color=col, ms=6)
        print(f"    그림15b: {name} 계수 {tol} λ → RMS {r_tol:.4f} λ, 스트렐 {s_tol:.3f}")
    b.axhline(0.8, color=C_GREY, lw=0.8)
    b.axvline(1 / 14, color=C_GREY, lw=0.8)
    b.text(1 / 14 + 0.004, 0.95, r"$\lambda/14$", fontsize=10)
    b.text(0.15, 0.82, t(L, "스트렐 0.8", "Strehl 0.8"), fontsize=9, **font(L))
    b.annotate(t(L, "보른·울프의 허용치\n구면 0.94λ, 코마 0.60λ, 비점 0.35λ", "Born & Wolf tolerances\nspherical 0.94λ, coma 0.60λ, astig. 0.35λ"),
               (0.0715, 0.82), (0.095, 0.93), fontsize=8.5, arrowprops=dict(arrowstyle="->", lw=0.8), **font(L))
    b.set_xlim(0, 0.18); b.set_ylim(0, 1.02)
    b.set_xlabel(t(L, r"RMS 파면수차 $\sigma$ ($\lambda$, 기울기·초점을 뺌)", r"RMS wave aberration $\sigma$ (waves, tilt & focus removed)"), **font(L))
    b.set_ylabel(t(L, "최선 초점의 스트렐 비", "Strehl ratio at best focus"), **font(L))
    b.legend(fontsize=8.5, prop=dict(family=plt.rcParams["font.family"], size=8.5))
    b.set_title(t(L, "(b) 수차의 종류와 상관없이 RMS 가 정한다 (작을 때)", "(b) For small aberrations, RMS alone decides"), fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig15-strehl-rms.png")


FIGS = {2: figure2, 3: figure3, 4: figure4, 5: figure5, 6: figure6, 8: figure8, 9: figure9,
        11: figure11, 13: figure13, 14: figure14, 15: figure15}


def main():
    which = [int(a) for a in sys.argv[1:]] or sorted(FIGS)
    for lang, d in [("ko", BASE_DIR), ("en", os.path.join(BASE_DIR, "en"))]:
        L = {"lang": lang, "dir": d}
        setfonts(L)
        for k in which:
            print(f"[{lang}] 그림 {k}")
            FIGS[k](L)


if __name__ == "__main__":
    main()

"""파동광학 배경이론 5편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/wave-optics-gaussian-beams/      (한국어)
      ../../assets/img/posts/wave-optics-gaussian-beams/en/   (영문)

그림 1 은 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
그림 2 는 공용 사진(CC0) ext-green-laser-tem00.jpg 의 세기 단면에 가우스 곡선을 맞춘 것이다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import gaussian_filter
from scipy.optimize import curve_fit

from gaussian import (M_free, M_lens, W_R_from_q, abcd, angular_spectrum, fox_li, gaussian_field, grid, hermite_gaussian,
                      laguerre_gaussian, q_at_waist, rayleigh_range, resonator_mode_q, second_moment_width, thin_lens_phase,
                      waist_from_q)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "wave-optics-gaussian-beams")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"


def t(L, ko, en):
    s = ko if L["lang"] == "ko" else en
    s = s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$").replace("θ", r"$\theta$").replace("ζ", r"$\zeta$").replace("π", r"$\pi$")
    if "$" in s:
        s = s.replace("µ", r"$\mathrm{\mu}$")
    return s


def font(L):
    return {"fontfamily": ["AppleGothic", "DejaVu Sans"]} if L["lang"] == "ko" else {}


def setfonts(L):
    plt.rcParams["font.family"] = ["AppleGothic", "DejaVu Sans"] if L["lang"] == "ko" else ["DejaVu Sans"]


def legend(ax, L, **kw):
    ax.legend(prop=dict(family=plt.rcParams["font.family"], size=kw.pop("size", 8.5)), **kw)


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    path = os.path.join(L["dir"], name)
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("  저장:", path)


# ------------------------------------------------------------------ 그림 2: 레이저 포인터의 세기 단면

def figure2(L):
    img = Image.open(os.path.join(BASE_DIR, "ext-green-laser-tem00.jpg")).convert("RGB")
    g = np.asarray(img).astype(float)[:, :, 1]
    cy, cx = np.unravel_index(np.argmax(gaussian_filter(g, 15)), g.shape)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), gridspec_kw=dict(width_ratios=[1, 1.35]))
    axes[0].imshow(np.asarray(img))
    axes[0].axhline(cy, color="w", lw=0.8, ls="--")
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_title(t(L, "(a) 초록 레이저 포인터(5 mW)의 단면", "(a) Cross-section of a 5 mW green laser pointer"), fontsize=10.5, **font(L))
    row = g[cy - 5:cy + 6].mean(0)
    x = np.arange(row.size)
    f = lambda x, a, x0, W, b: a * np.exp(-2 * (x - x0) ** 2 / W**2) + b      # noqa: E731
    p, _ = curve_fit(f, x, row, p0=[row.max(), cx, 300, 0])
    b = axes[1]
    b.plot(x - p[1], row, color=C_GRN, lw=1, label=t(L, "사진의 초록 채널 (점선 가로줄)", "green channel along the dashed line"))
    b.plot(x - p[1], f(x, *p), color="k", lw=1.2, ls="--", label=t(L, "가우스 곡선 맞춤", "Gaussian fit"))
    b.set_xlabel(t(L, "가로 위치 (화소)", "horizontal position (pixels)"), **font(L))
    b.set_ylabel(t(L, "화소 값", "pixel value"), **font(L))
    resid = (row - f(x, *p)).std() / p[0]
    b.set_title(t(L, f"(b) 세기 단면과 가우스 곡선 (잔차 rms = 정점의 {resid * 100:.1f} %)", f"(b) Profile and Gaussian fit (rms residual = {resid * 100:.1f} % of peak)"), fontsize=10.5, **font(L))
    legend(b, L, loc="upper right")
    fig.tight_layout()
    save(fig, L, "fig2-tem00-profile.png")
    print(f"    그림2: 잔차 {resid:.4f}")


# ------------------------------------------------------------------ 그림 3: 가우시안 빔의 모양

def figure3(L):
    lam, W0 = 0.5, 5.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 512, 0.25
    x, (X, Y) = grid(n, dx)
    U0 = gaussian_field(X, Y, W0, lam)
    zs = np.linspace(-4, 4, 161) * z0
    sel = np.abs(x) <= 6 * W0 * 2.5
    side = np.array([np.abs(angular_spectrum(U0, dx, lam, z)[n // 2, sel]) ** 2 for z in zs])
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.3), gridspec_kw=dict(width_ratios=[1.6, 1]))
    a = axes[0]
    a.imshow(side.T, extent=[zs[0] / z0, zs[-1] / z0, x[sel][0] / W0, x[sel][-1] / W0], aspect="auto", origin="lower", cmap="inferno", vmax=1)
    zz = np.linspace(-4, 4, 400)
    Wz = np.sqrt(1 + zz**2)
    a.plot(zz, Wz, color="c", lw=1.2)
    a.plot(zz, -Wz, color="c", lw=1.2, label=r"$W(z) = W_0\sqrt{1 + (z/z_0)^2}$")
    a.plot(zz, np.abs(zz), color="w", lw=0.8, ls=":")
    a.plot(zz, -np.abs(zz), color="w", lw=0.8, ls=":", label=t(L, "먼 곳의 원뿔 (반각 θ0 = λ/πW0)", "far-field cone (half-angle θ0 = λ/πW0)"))
    for zv in (-1, 1):
        a.axvline(zv, color="c", lw=0.6, ls="--")
    a.set_xlabel("z / z0")
    a.set_ylabel("x / W0")
    a.set_ylim(-7, 7)
    a.set_title(t(L, "(a) 허리 반지름 W0 = 10λ 인 빔을 옆에서 본 세기 (각스펙트럼 전파). 세로 점선: ±z0",
                  "(a) Side view of a beam with waist W0 = 10λ (angular-spectrum propagation). Dashed: ±z0"), fontsize=10.5, **font(L))
    legend(a, L, loc="upper center", size=8)
    b = axes[1]
    zn = np.linspace(-4, 4, 17) * z0
    Wn = np.array([second_moment_width(np.abs(angular_spectrum(U0, dx, lam, z)) ** 2, x) for z in zn])
    b.plot(zz, Wz, color=C_RAY, lw=1.5, label=t(L, "폭 W/W0 (공식)", "width W/W0 (formula)"))
    b.plot(zn / z0, Wn / W0, "o", color=C_RAY, ms=4, mfc="none", label=t(L, "폭 (파동 전파에서 잰 값)", "width (measured from propagation)"))
    for side_, lab_ in [(zz[zz < -0.15], None), (zz[zz > 0.15], t(L, "곡률 반지름 R/(4z0)", "radius of curvature R/(4z0)"))]:
        b.plot(side_, side_ * (1 + 1 / side_**2) / 4, color=C_RED, lw=1.5, label=lab_)      # z = 0 에서 R → ∞ 이라 둘로 나눠 그린다
    b.axhline(0, color="k", lw=0.5)
    b.set_ylim(-1.6, 4.4)
    b.set_xlabel("z / z0")
    b.set_title(t(L, "(b) 폭과 파면의 곡률 반지름\n(R 은 ±z0 에서 가장 작은 2z0)", "(b) Width and wavefront radius of curvature\n(R is smallest, 2z0, at ±z0)"), fontsize=10.5, **font(L))
    legend(b, L, loc="upper left", size=8)
    fig.tight_layout()
    save(fig, L, "fig3-gaussian-beam-shape.png")


# ------------------------------------------------------------------ 그림 4: 구이 위상

def figure4(L):
    lam, W0 = 0.5, 5.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 512, 0.4
    x, (X, Y) = grid(n, dx)
    zz = np.linspace(-5, 5, 400)
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.2), gridspec_kw=dict(width_ratios=[1.15, 1]))
    a = axes[0]
    zn = np.linspace(-5, 5, 21) * z0
    for (l, m, col, lab) in [(0, 0, C_RAY, "HG00"), (1, 0, C_GRN, "HG10"), (1, 1, C_RED, "HG11")]:
        U0 = hermite_gaussian(X, Y, W0, l, m).astype(complex)
        ph = []
        for z in zn:
            U = angular_spectrum(U0, dx, lam, z)
            Wz = W0 * np.sqrt(1 + (z / z0) ** 2)
            curv = np.exp(1j * np.pi * (X**2 + Y**2) * z / (lam * (z**2 + z0**2)))
            T = hermite_gaussian(X, Y, Wz, l, m) * curv
            ph.append(-np.angle(np.sum(U * np.conj(T))))
        ph = np.unwrap(np.array(ph))
        ph -= ph[len(ph) // 2]
        a.plot(zz, (l + m + 1) * np.arctan(zz), color=col, lw=1.4, label=lab + r": $" + str(l + m + 1) + r"\,\tan^{-1}(z/z_0)$")
        a.plot(zn / z0, ph, "o", color=col, mfc="none", ms=4)
    a.axhline(0, color="k", lw=0.5)
    a.set_yticks(np.arange(-3, 4) * np.pi / 2)
    a.set_yticklabels([r"$-3\pi/2$", r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$"])
    a.set_xlabel("z / z0")
    a.set_ylabel(t(L, "평면파보다 늦은 위상", "phase lag behind a plane wave"), **font(L))
    a.set_title(t(L, "(a) 구이 위상: 선은 공식, 원은 파동 전파에서 잰 값", "(a) Gouy phase: lines are the formula, circles measured from propagation"), fontsize=10.5, **font(L))
    legend(a, L, loc="upper left", size=8)
    # (b) 공진기의 공명 주파수: 가로 모드가 구이 위상만큼 밀린다
    lam_r, Lc, R = 0.6328, 300e3, 600e3
    q = resonator_mode_q(Lc, R, R, lam_r)
    _, z0r = None, None
    zw, W0r = waist_from_q(abcd(M_free(Lc / 2), q), lam_r)
    z0r = rayleigh_range(W0r, lam_r)
    dzeta = 2 * np.arctan(Lc / 2 / z0r)
    nuF = 3e8 / (2 * Lc * 1e-6) / 1e6                       # MHz
    b = axes[1]
    for qn in range(0, 3):
        for N, col in zip(range(0, 4), [C_RAY, C_GRN, C_KEY, C_RED]):
            nu = qn * nuF + N * dzeta / np.pi * nuF
            if nu <= 2.3 * nuF:
                b.vlines(nu, 0, 1 / (N + 1), color=col, lw=2.2, zorder=3 + N, label=f"l + m = {N}" if qn == 0 else None)
    b.set_xlabel(t(L, "공명 주파수 − 기준 (MHz)", "resonance frequency − reference (MHz)"), **font(L))
    b.set_yticks([])
    b.set_title(t(L, f"(b) 공진기(L = 30 cm, 거울 곡률 60 cm)의 공명선\n세로 모드 간격 {nuF:.0f} MHz, 가로 모드는 {dzeta / np.pi * nuF:.0f} MHz 씩 밀린다",
                  f"(b) Resonances of a cavity (L = 30 cm, mirror radius 60 cm)\nlongitudinal spacing {nuF:.0f} MHz, transverse modes shifted by {dzeta / np.pi * nuF:.0f} MHz"), fontsize=10, **font(L))
    legend(b, L, loc="upper right", size=8)
    fig.tight_layout()
    save(fig, L, "fig4-gouy-phase.png")
    print(f"    그림4: Δζ/π = {dzeta / np.pi:.4f}, 가로 모드 간격 {dzeta / np.pi * nuF:.1f} MHz")


# ------------------------------------------------------------------ 그림 5: 렌즈로 모으기

def figure5(L):
    lam, W0, f = 0.5, 20.0, 1000.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 1024, 0.5
    x, (X, Y) = grid(n, dx)
    d1 = z0
    U = angular_spectrum(gaussian_field(X, Y, W0, lam), dx, lam, d1) * thin_lens_phase(X, Y, f, lam)
    zs = np.linspace(0, 2.2 * f, 111)
    sel = np.abs(x) <= 50
    side = np.array([np.abs(angular_spectrum(U, dx, lam, z)[n // 2, sel]) ** 2 for z in zs])
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.3), gridspec_kw=dict(width_ratios=[1.45, 1, 1]))
    a = axes[0]
    a.imshow(side.T, extent=[0, zs[-1] / 1000, x[sel][0], x[sel][-1]], aspect="auto", origin="lower", cmap="inferno", vmax=side.max() * 0.8)
    q1 = d1 + q_at_waist(W0, lam)
    zq = np.linspace(0, 2.2 * f, 300)
    Wq = np.array([W_R_from_q(abcd(M_free(z) @ M_lens(f), q1), lam)[0] for z in zq])
    a.plot(zq / 1000, Wq, color="c", lw=1.2, label=t(L, "ABCD 법칙의 W(z)", "W(z) from the ABCD law"))
    a.plot(zq / 1000, -Wq, color="c", lw=1.2)
    zw, Ww = waist_from_q(abcd(M_lens(f), q1), lam)
    a.axvline(-zw / 1000, color="w", lw=0.8, ls="--")
    a.axvline(f / 1000, color=C_GRN, lw=0.8, ls=":")
    a.text(-zw / 1000 + 0.02, 40, t(L, "허리", "waist"), color="w", fontsize=9, **font(L))
    a.text(f / 1000 - 0.02, -44, t(L, "초점 f", "focus f"), color=C_GRN, fontsize=9, ha="right", **font(L))
    a.set_xlabel(t(L, "렌즈에서 떨어진 거리 (mm)", "distance from the lens (mm)"), **font(L))
    a.set_ylabel("x (µm)")
    a.set_title(t(L, "(a) 렌즈(f = 1 mm) 뒤의 빔. 허리(W0 = 20 µm)가 렌즈 앞 z0 = 2.5 mm 에 있다", "(a) Beam after a lens (f = 1 mm); input waist W0 = 20 µm at z0 = 2.5 mm before it"), fontsize=10, **font(L))
    legend(a, L, loc="lower left", size=8)
    # (b) 출력 허리의 위치: 가우시안 빔 대 기하광학
    b = axes[1]
    s = np.linspace(-6, 6, 600)                   # 입력 허리의 렌즈 앞 거리 / f
    for zr_f, col in [(0.0, "k"), (0.25, C_RAY), (0.5, C_GRN), (1.0, C_RED)]:
        out = []
        for sv in s:
            if zr_f == 0:
                out.append(1 / (1 / f - 1 / (sv * f)) / f if abs(sv - 1) > 1e-6 else np.nan)
            else:
                W0i = np.sqrt(lam * zr_f * f / np.pi)
                zo, _ = waist_from_q(abcd(M_lens(f), sv * f + q_at_waist(W0i, lam)), lam)
                out.append(-zo / f)
        out = np.array(out)
        if zr_f == 0:
            out[np.abs(out) > 8] = np.nan
        b.plot(s, out, color=col, lw=1.4, ls="--" if zr_f == 0 else "-",
               label=t(L, "기하광학 (z0 = 0)", "geometrical optics (z0 = 0)") if zr_f == 0 else f"z0 = {zr_f:g} f")
    b.set_xlim(-1, 6)
    b.set_ylim(-4, 8)
    b.axhline(1, color=C_GREY, lw=0.5)
    b.set_xlabel(t(L, "입력 허리의 렌즈 앞 거리 / f", "input waist distance before the lens / f"), **font(L))
    b.set_ylabel(t(L, "출력 허리의 렌즈 뒤 거리 / f", "output waist distance after the lens / f"), **font(L))
    b.set_title(t(L, "(b) 출력 허리의 위치\n(빔이 가늘면 상의 자리로 가지 않는다)", "(b) Position of the output waist\n(a narrow beam does not go to the image point)"), fontsize=10, **font(L))
    legend(b, L, loc="upper right", size=8)
    # (c) F수에 따른 초점의 크기와 깊이
    c = axes[2]
    F = np.linspace(0.5, 20, 200)
    lam_c = 0.532
    spot = 4 / np.pi * lam_c * F
    dof = 2 * np.pi * (spot / 2) ** 2 / lam_c
    c.plot(F, spot, color=C_RAY, lw=1.6, label=t(L, "초점 지름 2W0' = (4/π) λ F#", "spot diameter 2W0' = (4/π) λ F#"))
    c.plot(F, dof, color=C_RED, lw=1.6, label=t(L, "초점 깊이 2z0'", "depth of focus 2z0'"))
    c.set_yscale("log")
    c.set_xscale("log")
    c.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    c.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    c.set_xlabel(t(L, "F수 F# = f/D (D = 2W0)", "f-number F# = f/D (D = 2W0)"), **font(L))
    c.set_ylabel("µm")
    c.set_title(t(L, "(c) 평행 빔을 모을 때 (532 nm)\n지름은 F#, 깊이는 $F_\\#^2$ 에 비례", "(c) Focusing a collimated beam (532 nm)\ndiameter ∝ F#, depth ∝ $F_\\#^2$"), fontsize=10, **font(L))
    legend(c, L, loc="upper left", size=8)
    fig.tight_layout()
    save(fig, L, "fig5-lens-focusing.png")
    print(f"    그림5: 허리 {-zw:.1f} µm 에서 W {Ww:.2f} µm")


# ------------------------------------------------------------------ 그림 6: 공진기와 폭스-리 반복

def figure6(L):
    lam, Lc, R = 0.6328, 300e3, 600e3
    q = resonator_mode_q(Lc, R, R, lam)
    Wm, _ = W_R_from_q(q, lam)
    n, dx = 512, 8.0
    x, _ = grid(n, dx)
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.4), gridspec_kw=dict(width_ratios=[1.15, 1.1, 0.9]))
    a = axes[0]
    rec = (1, 2, 5, 20, 1500)
    U, losses, snaps = fox_li(n, dx, lam, Lc, R, a=600.0, n_iter=1500, record=rec)
    sel = np.abs(x) <= 800
    cols = [C_GREY, C_KEY, C_GRN, C_RAY, C_RED]
    for it, col in zip(rec, cols):
        prof = np.abs(snaps[it][n // 2]) ** 2
        a.plot(x[sel], prof[sel] / prof.max(), color=col, lw=1.3, label=t(L, f"{it}번 왕복", f"after {it} pass{'es' if it > 1 else ''}"))
    a.plot(x[sel], np.exp(-2 * x[sel] ** 2 / Wm**2), "k--", lw=1, label=t(L, "자기 일치 q 의 가우시안", "self-consistent Gaussian"))
    a.axvline(600, color=C_GREY, lw=0.6, ls=":")
    a.axvline(-600, color=C_GREY, lw=0.6, ls=":")
    a.set_xlabel(t(L, "거울 위 위치 x (µm)", "position on the mirror x (µm)"), **font(L))
    a.set_ylabel(t(L, "세기 (정점 = 1)", "intensity (peak = 1)"), **font(L))
    a.set_title(t(L, "(a) 폭스-리 반복: 균일한 빛에서 출발\n(L = 30 cm, 거울 곡률 60 cm, 구멍 반지름 0.6 mm)",
                  "(a) Fox–Li iteration from uniform light\n(L = 30 cm, mirror radius 60 cm, aperture radius 0.6 mm)"), fontsize=10, **font(L))
    legend(a, L, loc="upper right", size=7.5)
    b = axes[1]
    for ap, col in [(450.0, C_RED), (600.0, C_RAY), (1000.0, C_GRN)]:
        _, ls, _ = fox_li(n, dx, lam, Lc, R, a=ap, n_iter=600)
        NF = ap**2 / (lam * Lc)
        b.semilogy(np.arange(1, 601), np.maximum(ls, 1e-9), color=col, lw=1.2, label=t(L, f"구멍 {ap / 1000:g} mm (N_F = {NF:.1f})", f"aperture {ap / 1000:g} mm (N_F = {NF:.1f})"))
    b.set_xlabel(t(L, "왕복 횟수", "number of passes"), **font(L))
    b.set_ylabel(t(L, "한 번 왕복에서 잃는 비율", "loss per pass"), **font(L))
    b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.set_title(t(L, "(b) 왕복할 때마다 잃는 빛\n구멍이 크면 고차 모드가 오래 남아 출렁인다", "(b) Light lost per pass\nwith a large aperture, higher modes linger and the loss oscillates"), fontsize=10, **font(L))
    legend(b, L, loc="upper right", size=8)
    # (c) 안정성 그림
    c = axes[2]
    g = np.linspace(-2, 2, 401)
    G1, G2 = np.meshgrid(g, g)
    stable = (G1 * G2 >= 0) & (G1 * G2 <= 1)
    c.contourf(G1, G2, stable, levels=[0.5, 1.5], colors=["#cfe3f5"])
    c.axhline(0, color="k", lw=0.5)
    c.axvline(0, color="k", lw=0.5)
    pts = [((1, 1), t(L, "평면", "planar")), ((0, 0), t(L, "공초점", "confocal")), ((-1, -1), t(L, "동심", "concentric")), ((0.5, 0.5), t(L, "(a)의 공진기", "cavity of (a)"))]
    for (g1, g2), lab in pts:
        c.plot(g1, g2, "o", color=C_RED, ms=5)
        c.text(g1 + 0.08, g2 + 0.08, lab, fontsize=9, **font(L))
    c.set_xlim(-2, 2)
    c.set_ylim(-2, 2)
    c.set_aspect("equal")
    c.set_xlabel("g1 = 1 − L/R1")
    c.set_ylabel("g2 = 1 − L/R2")
    c.set_title(t(L, "(c) 안정 영역 0 ≤ g1 g2 ≤ 1 (파랑)", "(c) Stable region 0 ≤ g1 g2 ≤ 1 (blue)"), fontsize=10, **font(L))
    fig.tight_layout()
    save(fig, L, "fig6-resonator-fox-li.png")
    print(f"    그림6: 거울 위 폭 이론 {Wm:.1f} µm, 반복 {second_moment_width(np.abs(U) ** 2, x):.1f} µm, 손실 {losses[-1]:.2e}")


# ------------------------------------------------------------------ 그림 7: 고차 모드와 M²

def figure7(L):
    lam, W0 = 0.5, 5.0
    z0 = rayleigh_range(W0, lam)
    n, dx = 256, 0.15
    x, (X, Y) = grid(n, dx)
    fig = plt.figure(figsize=(16, 4.6))
    gs = fig.add_gridspec(2, 6, width_ratios=[1, 1, 1, 1, 0.45, 2.6], hspace=0.3, wspace=0.12)
    modes = [(hermite_gaussian(X, Y, W0, 0, 0), "HG00"), (hermite_gaussian(X, Y, W0, 1, 0), "HG10"),
             (hermite_gaussian(X, Y, W0, 1, 1), "HG11"), (hermite_gaussian(X, Y, W0, 2, 1), "HG21"),
             (laguerre_gaussian(X, Y, W0, 0, 1), "LG$_0^1$"), (laguerre_gaussian(X, Y, W0, 1, 0), "LG$_1^0$"),
             (laguerre_gaussian(X, Y, W0, 0, 3), "LG$_0^3$"), (None, None)]
    sel = slice(n // 2 - 80, n // 2 + 80)
    for i, (U, lab) in enumerate(modes):
        ax = fig.add_subplot(gs[i // 4, i % 4])
        if U is None:
            Uc = laguerre_gaussian(X, Y, W0, 0, 1)
            ax.imshow(np.angle(Uc)[sel, sel], cmap="twilight", origin="lower")
            ax.set_title(t(L, "LG$_0^1$ 의 위상 (나선)", "phase of LG$_0^1$ (helix)"), fontsize=9.5, **font(L))
        else:
            ax.imshow(np.abs(U[sel, sel]) ** 2, cmap="inferno", origin="lower")
            ax.set_title(lab, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])
    # M²
    b = fig.add_subplot(gs[:, 5])
    n2, dx2 = 512, 0.5
    x2, (X2, Y2) = grid(n2, dx2)
    zs = np.linspace(-3, 3, 25) * z0
    m0 = hermite_gaussian(X2, Y2, W0, 0, 0).astype(complex)
    m1 = hermite_gaussian(X2, Y2, W0, 1, 0).astype(complex)
    m0 /= np.sqrt(np.sum(np.abs(m0) ** 2))
    m1 /= np.sqrt(np.sum(np.abs(m1) ** 2))
    Wg, Wmix = [], []
    for z in zs:
        I0 = np.abs(angular_spectrum(m0, dx2, lam, z)) ** 2
        I1 = np.abs(angular_spectrum(m1, dx2, lam, z)) ** 2
        Wg.append(second_moment_width(I0, x2))
        Wmix.append(second_moment_width(I0 + I1, x2))
    zz = np.linspace(-3, 3, 200)
    b.plot(zz, np.sqrt(1 + zz**2), color=C_RAY, lw=1.4, label=t(L, "가우시안 ($M^2$ = 1)", "Gaussian ($M^2$ = 1)"))
    b.plot(zs / z0, np.array(Wg) / W0, "o", color=C_RAY, mfc="none", ms=4)
    b.plot(zz, np.sqrt(2) * np.sqrt(1 + zz**2), color=C_RED, lw=1.4, label=t(L, "HG00 + HG10 결맞지 않게 섞음 ($M^2$ = 2)", "HG00 + HG10 mixed incoherently ($M^2$ = 2)"))
    b.plot(zs / z0, np.array(Wmix) / W0, "o", color=C_RED, mfc="none", ms=4)
    b.set_xlabel("z / z0")
    b.set_ylabel(t(L, "x 방향 2차 모멘트 폭 / W0", "second-moment width in x / W0"), **font(L))
    b.set_title(t(L, "(b) 빔 품질: 폭과 발산이 함께 $\\sqrt{M^2}$ 배\n(원: 파동 전파에서 잰 값)", "(b) Beam quality: width and divergence both × $\\sqrt{M^2}$\n(circles: measured from propagation)"), fontsize=10, **font(L))
    legend(b, L, loc="upper center", size=8.5)
    fig.text(0.02, 0.97, t(L, "(a) 고차 모드의 세기 (허리에서)", "(a) Intensity of higher-order modes (at the waist)"), fontsize=10.5, va="top", **font(L))
    save(fig, L, "fig7-higher-modes-m2.png")


FIGS = {2: figure2, 3: figure3, 4: figure4, 5: figure5, 6: figure6, 7: figure7}


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

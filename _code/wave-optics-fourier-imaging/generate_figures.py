"""파동광학 배경이론 4편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/wave-optics-fourier-imaging/      (한국어)
      ../../assets/img/posts/wave-optics-fourier-imaging/en/   (영문)

그림 1·8 은 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
모든 계산은 λ = 0.5 µm, NA = 0.5 (λ/NA = 1 µm) 로 둔다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from imaging import (coherent_image, dark_field, filter_4f, grid, incoherent_image, modulation, otf_by_overlap,
                     otf_from_psf, partially_coherent_image, phase_contrast, spoke_target)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "wave-optics-fourier-imaging")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"
LAM, NA = 0.5, 0.5


def t(L, ko, en):
    s = ko if L["lang"] == "ko" else en
    s = s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$").replace("θ", r"$\theta$").replace("σ", r"$\sigma$").replace("φ", r"$\varphi$")
    if "$" in s:
        s = s.replace("µ", r"$\mathrm{\mu}$")
    return s


def font(L):
    return {"fontfamily": ["AppleGothic", "DejaVu Sans"]} if L["lang"] == "ko" else {}


def setfonts(L):
    plt.rcParams["font.family"] = ["AppleGothic", "DejaVu Sans"] if L["lang"] == "ko" else ["DejaVu Sans"]


def legend(ax, L, **kw):
    ax.legend(prop=dict(family=plt.rcParams["font.family"], size=kw.pop("size", 8.5)), **kw)


def save(fig, L, name, palette=False):
    os.makedirs(L["dir"], exist_ok=True)
    path = os.path.join(L["dir"], name)
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    if palette:
        from PIL import Image
        Image.open(path).convert("RGB").quantize(256, dither=Image.Dither.NONE).save(path, optimize=True)
    print("  저장:", path)


def noticks(ax):
    ax.set_xticks([])
    ax.set_yticks([])


# ------------------------------------------------------------------ 그림 2: 렌즈의 푸리에 변환과 4f 계

def figure2(L):
    fig = plt.figure(figsize=(15, 4.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.55, 1, 1, 1], wspace=0.25)
    a = fig.add_subplot(gs[0])
    # 4f 계 그림
    zs = [0, 1, 2, 3, 4]
    for z, lab_ko, lab_en in [(0, "물체", "object"), (2, "푸리에 면\n(후초점면)", "Fourier plane\n(back focal plane)"), (4, "상", "image")]:
        a.plot([z, z], [-1.1, 1.1], color="k", lw=1.2)
        a.text(z, -1.45, t(L, lab_ko, lab_en), ha="center", va="top", fontsize=9, **font(L))
    for z in (1, 3):
        a.add_patch(mpatches.Ellipse((z, 0), 0.16, 2.3, fc="#cfe3f5", ec=C_RAY))
        a.text(z, 1.3, t(L, "렌즈 f", "lens f"), ha="center", fontsize=9, **font(L))
    tn = 0.4                                   # 한 회절 차수: 물체 면의 여러 높이에서 같은 각으로 나가는 평행광
    for y0 in (-0.6, 0.0, 0.6):
        a.plot([0, 1, 2], [y0, y0 + tn, tn], color=C_RED, lw=0.9)
        a.plot([0, 1, 2], [y0, y0, 0.0], color=C_GRN, lw=0.9)
    a.plot(2, tn, "o", color=C_RED, ms=4)
    a.plot(2, 0, "o", color=C_GRN, ms=4)
    for z, txt in [(0.5, "f"), (1.5, "f"), (2.5, "f"), (3.5, "f")]:
        a.annotate("", xy=(z - 0.5, -1.25), xytext=(z + 0.5, -1.25), arrowprops=dict(arrowstyle="<->", lw=0.7))
        a.text(z, -1.22, txt, ha="center", va="bottom", fontsize=9)
    a.text(2.08, -0.6, t(L, "같은 방향의\n빛은\n한 점에", "one\ndirection\n→ one\npoint"), fontsize=8.5, color=C_RED, va="center", **font(L))
    a.set_xlim(-0.3, 4.4)
    a.set_ylim(-1.9, 1.6)
    a.axis("off")
    a.set_title(t(L, "(a) 4f 계: 앞 렌즈가 변환하고 뒤 렌즈가 되돌린다", "(a) 4f system: the first lens transforms, the second undoes it"), fontsize=10.5, **font(L))
    # 물체: 세로 줄무늬(주기 4 µm), 가로 줄무늬(주기 0.8 µm), 둥근 구멍
    n, dx = 512, 0.1
    x, f = grid(n, dx)
    X, Y = np.meshgrid(x, x)
    obj = np.zeros((n, n))
    obj[(X < -2) & (Y > 2) & (X > -20) & (Y < 20)] = (np.cos(2 * np.pi * X / 4) > 1e-9)[(X < -2) & (Y > 2) & (X > -20) & (Y < 20)]
    obj[(X > 2) & (Y > 2) & (X < 20) & (Y < 20)] = (np.cos(2 * np.pi * Y / 0.8) > 1e-9)[(X > 2) & (Y > 2) & (X < 20) & (Y < 20)]
    obj[np.hypot(X, Y + 10) < 5] = 1.0
    b = fig.add_subplot(gs[1])
    ext = [x[0], x[-1], x[0], x[-1]]
    b.imshow(obj, cmap="gray", extent=ext, origin="lower")
    b.set_title(t(L, "(b) 물체 (투과율)", "(b) Object (transmittance)"), fontsize=10.5, **font(L))
    b.set_xlabel("x (µm)")
    F = np.fft.fftshift(np.abs(np.fft.fft2(obj)) ** 2)
    fs = np.fft.fftshift(f)
    c = fig.add_subplot(gs[2])
    c.imshow(np.log10(F / F.max() + 1e-7), cmap="inferno", vmin=-6, vmax=0, extent=[fs[0], fs[-1], fs[0], fs[-1]], origin="lower")
    c.set_xlim(-1.4, 1.4)
    c.set_ylim(-1.4, 1.4)
    c.add_patch(mpatches.Circle((0, 0), NA / LAM, fill=False, ec="c", ls="--", lw=1))
    c.set_xlabel(t(L, "공간 주파수 (주기/µm)", "spatial frequency (cycles/µm)"), **font(L))
    c.set_title(t(L, "(c) 후초점면의 세기 (로그)\n점선: NA 0.5 렌즈가 받는 범위", "(c) Back-focal-plane\nintensity (log); dashes:\nwhat an NA 0.5 lens accepts"), fontsize=10, **font(L))
    d = fig.add_subplot(gs[3])
    Ic = np.abs(coherent_image(obj.astype(complex), dx, NA, LAM)) ** 2
    d.imshow(Ic, cmap="gray", extent=ext, origin="lower", vmax=1.2)
    d.set_title(t(L, "(d) 점선 안만 지난 상 (결맞음)\n가로 줄무늬(주기 0.8 µm)가 사라졌다", "(d) Image from inside\nthe circle (coherent):\n0.8 µm stripes are gone"), fontsize=10, **font(L))
    d.set_xlabel("x (µm)")
    save(fig, L, "fig2-lens-fourier-4f.png")


# ------------------------------------------------------------------ 그림 3: 아베-포터 실험

def figure3(L):
    n, dx = 256, 0.25
    x, f = grid(n, dx)
    X, Y = np.meshgrid(x, x)
    per = 4.0
    wire = (np.mod(X, per) < 0.75) | (np.mod(Y, per) < 0.75)
    t_obj = (~wire).astype(float)
    rng = np.random.default_rng(4)
    for _ in range(14):                                   # 먼지 몇 알
        cx, cy = rng.uniform(-26, 26, 2)
        t_obj[np.hypot(X - cx, Y - cy) < rng.uniform(0.8, 1.6)] = 0.0
    U = t_obj.astype(complex)
    FX, FY = np.meshgrid(f, f)
    f1 = 1 / per
    on_order = (np.abs(FX / f1 - np.round(FX / f1)) < 0.02) & (np.abs(FY / f1 - np.round(FY / f1)) < 0.02)
    masks = [
        (np.ones((n, n)), t(L, "(a) 모두 통과", "(a) All pass")),
        ((np.hypot(FX, FY) < 0.5 * f1).astype(float), t(L, "(b) 가운데만\n(0차와 그 주변)", "(b) Centre only\n(0th order and its surroundings)")),
        ((np.abs(FY) < 0.5 * f1).astype(float), t(L, "(c) 가로 슬릿", "(c) Horizontal slit")),
        ((np.hypot(FX, FY) < 1.6 * f1).astype(float), t(L, "(d) 작은 조리개\n(0차와 ±1차)", "(d) Small iris\n(0th and ±1st)")),
        (1.0 - (np.hypot(FX, FY) < 0.02 * f1), t(L, "(e) 0차만 막음", "(e) 0th order blocked")),
        (1.0 - on_order, t(L, "(f) 그물의 차수를 모두 막음", "(f) All mesh orders blocked")),
    ]
    S = np.fft.fftshift(np.abs(np.fft.fft2(U)) ** 2)
    S = np.log10(S / S.max() + 1e-6)
    fig, axes = plt.subplots(2, 6, figsize=(15.5, 5.6))
    sl = slice(n // 2 - 40, n // 2 + 40)
    for i, (m, title) in enumerate(masks):
        a = axes[0, i]
        a.imshow(S[sl, sl], cmap="inferno", vmin=-5, vmax=0, origin="lower")
        a.imshow(np.ma.masked_where(np.fft.fftshift(m)[sl, sl] > 0.5, np.ones((80, 80))), cmap="Blues", vmin=0, vmax=1.5, alpha=0.75, origin="lower")
        if i == 4:
            a.plot(40, 40, "o", mfc="none", mec="c", ms=9, mew=1.5)      # 막은 0차 자리
        a.set_title(title, fontsize=10, **font(L))
        noticks(a)
        I = np.abs(filter_4f(U, m)) ** 2
        b = axes[1, i]
        b.imshow(I, cmap="gray", origin="lower")
        noticks(b)
    axes[0, 0].set_ylabel(t(L, "푸리에 면\n(파랑: 가림)", "Fourier plane\n(blue: blocked)"), fontsize=10, **font(L))
    axes[1, 0].set_ylabel(t(L, "상의 세기", "image intensity"), fontsize=10, **font(L))
    fig.tight_layout()
    save(fig, L, "fig3-abbe-porter.png", palette=True)


# ------------------------------------------------------------------ 그림 4: 결맞는 결상과 결맞지 않는 결상 — 가장자리와 두 점

def figure4(L):
    n, dx = 1024, 0.02
    x, f = grid(n, dx)
    xe = x + dx / 2
    step = ((xe > 0) & (xe < n * dx / 2 - 2)).astype(float)
    U = np.tile(step, (n, 1)).astype(complex)
    Ic = (np.abs(coherent_image(U, dx, NA, LAM)) ** 2)[n // 2]
    Ii = incoherent_image(np.abs(U) ** 2, dx, NA, LAM)[n // 2]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.0), gridspec_kw=dict(width_ratios=[1.15, 1.15, 0.85, 0.85]))
    a = axes[0]
    s = (xe > -2.5) & (xe < 3)
    a.plot(xe[s], step[s], color=C_GREY, lw=1, ls=":", label=t(L, "물체", "object"))
    a.plot(xe[s], Ic[s], color=C_RED, lw=1.6, label=t(L, "결맞는 조명", "coherent"))
    a.plot(xe[s], Ii[s], color=C_RAY, lw=1.6, label=t(L, "결맞지 않는 조명", "incoherent"))
    a.axvline(0, color="k", lw=0.6)
    for yv in (0.25, 0.5):
        a.axhline(yv, color=C_GREY, lw=0.6, ls="--")
    a.text(-2.4, 0.27, "1/4", fontsize=9)
    a.text(-2.4, 0.52, "1/2", fontsize=9)
    a.set_xlabel(t(L, "x / (λ/NA), 0 이 물체의 경계", "x / (λ/NA); 0 is the object edge"), **font(L))
    a.set_ylabel(t(L, "세기", "intensity"), **font(L))
    a.set_title(t(L, "(a) 가장자리의 상", "(a) Image of an edge"), fontsize=10.5, **font(L))
    legend(a, L, loc="lower right")
    # 두 점 (레일리 간격 0.61 λ/NA)
    b = axes[1]
    sep = 0.61
    xx = np.linspace(-1.6, 1.6, 641)
    c = NA / LAM
    g = np.linspace(-c, c, 401)
    FX, FY = np.meshgrid(g, g)
    fx = FX[np.hypot(FX, FY) <= c]
    def amp(x0):
        return np.array([np.sum(np.exp(2j * np.pi * fx * (xv - x0))) for xv in xx]) / fx.size
    A1, A2 = amp(-sep / 2), amp(sep / 2)
    Iinc = np.abs(A1) ** 2 + np.abs(A2) ** 2
    b.plot(xx, Iinc / Iinc.max(), color=C_RAY, lw=2.2, label=t(L, "결맞지 않음", "incoherent"))
    for ph, col, ls, lab_ko, lab_en in [(0, C_RED, "-", "결맞음, 위상차 0", "coherent, phase 0"),
                                        (np.pi / 2, C_KEY, "--", "결맞음, 위상차 π/2", "coherent, phase π/2"),
                                        (np.pi, C_PUR, "-", "결맞음, 위상차 π", "coherent, phase π")]:
        I = np.abs(A1 + np.exp(1j * ph) * A2) ** 2
        b.plot(xx, I / I.max(), color=col, lw=1.3, ls=ls, label=t(L, lab_ko, lab_en))
    for xv in (-sep / 2, sep / 2):
        b.axvline(xv, color=C_GREY, lw=0.6, ls=":")
    b.set_xlabel("x / (λ/NA)")
    b.set_title(t(L, "(b) 레일리 간격(0.61 λ/NA)의 두 점", "(b) Two points at the Rayleigh spacing (0.61 λ/NA)"), fontsize=10.5, **font(L))
    legend(b, L, loc="upper right", size=7.5)
    b.set_ylim(0, 1.45)
    # 막대 셋 (간격 1.2 λ/NA)
    n2, dx2 = 256, 0.04
    x2, _ = grid(n2, dx2)
    X2, Y2 = np.meshgrid(x2, x2)
    bars = np.zeros((n2, n2))
    for k in (-1, 0, 1):
        bars[(np.abs(X2 - 0.7 * k) < 0.2) & (np.abs(Y2) < 2.5)] = 1.0
    Uc = np.abs(coherent_image(bars.astype(complex), dx2, NA, LAM)) ** 2
    Ui = incoherent_image(bars, dx2, NA, LAM)
    for ax, img, ttl in [(axes[2], Uc, t(L, "(c) 막대 셋(간격 0.7 λ/NA), 결맞음", "(c) Three bars (spacing 0.7 λ/NA), coherent")),
                         (axes[3], Ui, t(L, "(d) 같은 막대, 결맞지 않음", "(d) Same bars, incoherent"))]:
        ax.imshow(img[64:192, 64:192], cmap="gray", origin="lower", vmax=img.max(), extent=[-2.56, 2.56, -2.56, 2.56])
        ax.set_title(ttl, fontsize=10.5, **font(L))
        ax.set_xlabel("x / (λ/NA)")
    fig.tight_layout()
    save(fig, L, "fig4-coherent-incoherent-edge-points.png")
    print(f"    그림4: 결맞음 넘침 {Ic.max():.3f}")


# ------------------------------------------------------------------ 그림 5: 전달함수

def figure5(L):
    fig, axes = plt.subplots(1, 4, figsize=(16, 3.9), gridspec_kw=dict(width_ratios=[1.2, 0.9, 1, 1]))
    nu = np.linspace(0, 2.2, 441)
    f = nu * NA / LAM
    a = axes[0]
    a.plot(nu, (nu <= 1).astype(float), color=C_RED, lw=1.8, label=t(L, "결맞음: 진폭 전달함수 (= 동공)", "coherent: amplitude transfer (= pupil)"))
    otf = np.real(otf_by_overlap(f, NA, LAM, m=601))
    a.plot(nu, otf, color=C_RAY, lw=1.8, label=t(L, "결맞지 않음: OTF (= MTF)", "incoherent: OTF (= MTF)"))
    a.set_xlabel(t(L, "공간 주파수 / (NA/λ)", "spatial frequency / (NA/λ)"), **font(L))
    a.set_ylabel(t(L, "전달률", "transfer"), **font(L))
    a.set_ylim(-0.02, 1.15)
    a.set_title(t(L, "(a) 수차 없는 원형 동공의 전달함수", "(a) Transfer functions of an aberration-free circular pupil"), fontsize=10.5, **font(L))
    legend(a, L, loc="upper right", size=8)
    # 동공 겹침 그림
    b = axes[1]
    s = 0.7
    b.add_patch(mpatches.Circle((-s / 2, 0), 1, fc=C_RAY, alpha=0.25, ec=C_RAY))
    b.add_patch(mpatches.Circle((s / 2, 0), 1, fc=C_RED, alpha=0.25, ec=C_RED))
    b.annotate("", xy=(s / 2, 0), xytext=(-s / 2, 0), arrowprops=dict(arrowstyle="<->", lw=1))
    b.text(0, 0.1, r"$\nu$", ha="center", fontsize=10)
    b.set_xlim(-1.6, 1.6)
    b.set_ylim(-1.3, 1.3)
    b.set_aspect("equal")
    b.axis("off")
    b.set_title(t(L, "(b) OTF = 겹친 넓이 / 동공 넓이\n(두 동공을 주파수 ν 만큼 어긋나게)", "(b) OTF = overlap area / pupil area\n(pupils offset by frequency ν)"), fontsize=10, **font(L))
    # Goodman 그림 6.17: 진폭 cos 과 |cos| (세기는 둘 다 cos²)
    n, dx = 256, 0.05
    x, _ = grid(n, dx)
    ft = 0.75 * NA / LAM
    ft = round(ft * n * dx) / (n * dx)                     # 창에 맞춘다
    A = np.tile(np.cos(2 * np.pi * ft * x), (n, 1)).astype(complex)
    B = np.abs(A)
    sel = slice(n // 2, n // 2 + int(4 / (ft * dx)))
    for ax, kind in [(axes[2], "coh"), (axes[3], "inc")]:
        ax.plot(x[sel], (np.abs(A[0]) ** 2)[sel], color=C_GREY, lw=1, ls=":", label=t(L, "물체 세기 (둘 다 같다)", "object intensity (same for both)"))
        if kind == "coh":
            IA = (np.abs(coherent_image(A, dx, NA, LAM)) ** 2)[n // 2]
            IB = (np.abs(coherent_image(B, dx, NA, LAM)) ** 2)[n // 2]
            ttl = t(L, "(c) 결맞는 조명", "(c) Coherent illumination")
        else:
            IA = incoherent_image(np.abs(A) ** 2, dx, NA, LAM)[n // 2]
            IB = incoherent_image(np.abs(B) ** 2, dx, NA, LAM)[n // 2]
            ttl = t(L, "(d) 결맞지 않는 조명", "(d) Incoherent illumination")
        ax.plot(x[sel], IA[sel], color=C_RED, lw=1.6, label=t(L, "물체 A: 진폭 cos", "object A: amplitude cos"))
        ax.plot(x[sel], IB[sel], color=C_RAY, lw=1.6, ls="--", label=t(L, "물체 B: 진폭 |cos|", "object B: amplitude |cos|"))
        ax.set_ylim(-0.05, 1.45)
        ax.set_xlabel("x / (λ/NA)")
        ax.set_title(ttl, fontsize=10.5, **font(L))
        legend(ax, L, loc="upper right", size=7.5)
    axes[2].set_ylabel(t(L, "상의 세기", "image intensity"), **font(L))
    fig.tight_layout()
    save(fig, L, "fig5-transfer-functions.png")


# ------------------------------------------------------------------ 그림 6: 부분 결맞음

def figure6(L):
    sigmas = [0.0, 0.3, 0.6, 1.0]
    cols = [C_RED, C_KEY, C_GRN, C_RAY]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.0))
    # (a) 막대 격자 (주기 0.8 λ/NA)
    p = 0.8
    n, dx = 256, p / 32
    x, _ = grid(n, dx)
    U = np.tile((np.cos(2 * np.pi * x / p) > 1e-9).astype(float), (n, 1)).astype(complex)     # cos = 0 인 표본의 반올림 잡음을 피한다
    a = axes[0]
    sel = slice(n // 2, n // 2 + 96)
    a.plot(x[sel], np.abs(U[0, sel]) ** 2, color=C_GREY, lw=1, ls=":", label=t(L, "물체", "object"))
    for sg, col in zip(sigmas, cols):
        I = partially_coherent_image(U, dx, NA, LAM, sg)[n // 2]
        a.plot(x[sel], I[sel], color=col, lw=1.5, label=f"σ = {sg:g}")
    a.set_xlabel("x / (λ/NA)")
    a.set_ylabel(t(L, "상의 세기", "image intensity"), **font(L))
    a.set_title(t(L, "(a) 주기 0.8 λ/NA 막대 격자의 상", "(a) Image of a bar grating, period 0.8 λ/NA"), fontsize=10.5, **font(L))
    legend(a, L, loc="upper right", size=8)
    a.set_ylim(-0.05, 1.5)
    # (b) 가장자리
    n2, dx2 = 512, 0.04
    x2, _ = grid(n2, dx2)
    xe = x2 + dx2 / 2
    step = ((xe > 0) & (xe < n2 * dx2 / 2 - 2)).astype(float)
    U2 = np.tile(step, (n2, 1)).astype(complex)
    b = axes[1]
    s = (xe > -2) & (xe < 2.5)
    for sg, col in zip(sigmas, cols):
        I = partially_coherent_image(U2, dx2, NA, LAM, sg)[n2 // 2]
        b.plot(xe[s], I[s], color=col, lw=1.5, label=f"σ = {sg:g}")
    b.axvline(0, color="k", lw=0.6)
    for yv in (0.25, 0.5):
        b.axhline(yv, color=C_GREY, lw=0.6, ls="--")
    b.set_xlabel("x / (λ/NA)")
    b.set_title(t(L, "(b) 가장자리의 상", "(b) Image of an edge"), fontsize=10.5, **font(L))
    legend(b, L, loc="lower right", size=8)
    # (c) 약한 사인 격자의 대비 전달
    c = axes[2]
    nu = np.linspace(0.05, 2.0, 40)
    n3 = 256
    for sg, col in zip(sigmas, cols):
        mt = []
        for v in nu:
            per = 1 / (v * NA / LAM)
            k = max(1, round(n3 * 0.05 / per))
            dx3 = k * per / n3
            x3, _ = grid(n3, dx3)
            U3 = np.tile(1 + 0.1 * np.cos(2 * np.pi * x3 / per), (n3, 1)).astype(complex)
            I = partially_coherent_image(U3, dx3, NA, LAM, sg, n_ring=8)[0]
            mt.append(modulation(I) / modulation(np.abs(U3[0]) ** 2))
        c.plot(nu, mt, color=col, lw=1.5, label=f"σ = {sg:g}")
    otf = np.real(otf_by_overlap(nu * NA / LAM, NA, LAM, m=401))
    c.plot(nu, otf, color="k", lw=1, ls="--", label=t(L, "결맞지 않음 (OTF)", "incoherent (OTF)"))
    c.set_xlabel(t(L, "공간 주파수 / (NA/λ)", "spatial frequency / (NA/λ)"), **font(L))
    c.set_ylabel(t(L, "상의 대비 / 물체의 대비", "image contrast / object contrast"), **font(L))
    c.set_title(t(L, "(c) 약한 사인 무늬의 대비 전달", "(c) Contrast transfer of a weak sinusoid"), fontsize=10.5, **font(L))
    legend(c, L, loc="upper right", size=8)
    fig.tight_layout()
    save(fig, L, "fig6-partial-coherence.png")


# ------------------------------------------------------------------ 그림 7: 초점 어긋남과 MTF

def figure7(L):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw=dict(width_ratios=[1.25, 1, 1]))
    nu = np.linspace(0, 2.0, 201)
    f = nu * NA / LAM
    a = axes[0]
    for Wd, col in zip([0, 0.25, 0.5, 0.75, 1.0], [C_GREY, C_RAY, C_GRN, C_KEY, C_RED]):
        W = (lambda rho, th, Wd=Wd: Wd * rho**2)
        O = np.real(otf_by_overlap(f, NA, LAM, W=W, m=501))
        a.plot(nu / 2, O, color=col, lw=1.6, label=f"$W_{{020}}$ = {Wd:g} λ")
    a.axhline(0, color="k", lw=0.6)
    a.set_xlabel(t(L, "공간 주파수 / 차단 주파수 (2NA/λ)", "spatial frequency / cutoff (2NA/λ)"), **font(L))
    a.set_ylabel("OTF")
    a.set_title(t(L, "(a) 초점 어긋남에 따른 OTF (원형 동공)\n음수는 명암이 뒤집힌다는 뜻", "(a) OTF with defocus (circular pupil)\nnegative means contrast reversal"), fontsize=10, **font(L))
    legend(a, L, loc="upper right", size=8)
    n, dx = 512, 0.05
    obj = spoke_target(n, 36)
    for ax, Wd, lab in [(axes[1], 0.0, "b"), (axes[2], 1.0, "c")]:
        W = (lambda rho, th, Wd=Wd: Wd * rho**2) if Wd else None
        I = incoherent_image(obj, dx, NA, LAM, W)
        ax.imshow(I, cmap="gray", origin="lower", extent=[-n * dx / 2, n * dx / 2] * 2)
        ttl = t(L, f"({lab}) 바큇살 표적, 초점", f"({lab}) Spoke target, in focus") if Wd == 0 else \
            t(L, f"({lab}) 초점 어긋남 $W_{{020}}$ = 1 λ\n주황 점선 두 고리 사이에서 명암이 뒤집힌다", f"({lab}) Defocus $W_{{020}}$ = 1 λ\ncontrast flips between the two orange dotted rings")
        ax.set_title(ttl, fontsize=10, **font(L))
        ax.set_xlabel("x (µm)")
        r_cut = 36 / (2 * np.pi * 2 * NA / LAM)
        ax.add_patch(mpatches.Circle((0, 0), r_cut, fill=False, ec="c", ls="--", lw=0.9))
        if Wd:                       # OTF 가 음수인 주파수 구간(첫 영점 ~ 마지막 음수 끝)에 해당하는 고리
            W1 = (lambda rho, th: rho**2)
            nu_f = np.linspace(0.02, 1.0, 400) * 2 * NA / LAM
            O = np.real(otf_by_overlap(nu_f * 0 + nu_f, NA, LAM, W=W1, m=501))
            neg = nu_f[O < 0]
            r_hi, r_lo = (36 / (2 * np.pi * neg.min()), 36 / (2 * np.pi * neg.max()))
            for rr in (r_lo, r_hi):
                ax.add_patch(mpatches.Circle((0, 0), rr, fill=False, ec=C_KEY, ls=":", lw=1.4))
            # 확대 삽입: 이미 계산한 영상 I 의 오른쪽 위 쪽 일부를 크게 보여 준다 (반전 고리 두 개가 걸리는 창)
            lim = 8.0                                  # 영상 안쪽(±12.8 µm)에서 첫 고리 둘레를 크게
            ins = ax.inset_axes([0.56, 0.02, 0.42, 0.42])
            ins.imshow(I, cmap="gray", origin="lower", extent=[-n * dx / 2, n * dx / 2] * 2)
            ins.set_xlim(0.0, lim)
            ins.set_ylim(0.0, lim)
            for rr in (r_lo, r_hi):
                ins.add_patch(mpatches.Circle((0, 0), rr, fill=False, ec=C_KEY, ls=":", lw=1.2))
            ins.set_xticks([]); ins.set_yticks([])
            for sp in ins.spines.values():
                sp.set_edgecolor(C_KEY)
                sp.set_linewidth(1.2)
            ins.text(0.04, 0.95, t(L, "확대", "zoom"), transform=ins.transAxes, color=C_KEY, fontsize=8, va="top", **font(L))
            print(f"    그림7(c): OTF<0 구간 반지름 {r_lo:.2f}–{r_hi:.2f} µm (주파수 {neg.min()/(2*NA/LAM):.2f}–{neg.max()/(2*NA/LAM):.2f} 차단)")
    fig.tight_layout()
    save(fig, L, "fig7-defocus-mtf-spoke.png", palette=True)


# ------------------------------------------------------------------ 그림 9: 위상 물체를 보는 다섯 가지 방법

def figure9(L):
    n, dx = 256, 0.1
    x, f = grid(n, dx)
    X, Y = np.meshgrid(x, x)
    def soft_disk(cx, cy, r, w=0.4):
        return 0.5 * (1 - np.tanh((np.hypot(X - cx, Y - cy) - r) / w))
    phi = 0.25 * soft_disk(0, 0, 8) + 0.15 * soft_disk(1.5, 1, 2.5)
    rng = np.random.default_rng(2)
    for _ in range(14):
        cx, cy = rng.uniform(-6, 6, 2)
        if np.hypot(cx - 1.5, cy - 1) > 3.2 and np.hypot(cx, cy) < 7:
            phi += 0.08 * soft_disk(cx, cy, 0.35, 0.15)
    U = np.exp(1j * phi)
    dot = 0.06 * NA / LAM                       # 위상판의 반지름 (공간 주파수): 0차와 아주 낮은 주파수를 함께 덮는다
    imgs = [
        (phi, t(L, "(a) 시편의 위상 (최대 0.4 rad)", "(a) Specimen phase (max 0.4 rad)"), "viridis"),
        (np.abs(coherent_image(U, dx, NA, LAM)) ** 2, t(L, "(b) 명시야, 초점", "(b) Brightfield, in focus"), "gray"),
        (np.abs(coherent_image(U, dx, NA, LAM, W=lambda r, th: 1.5 * r**2)) ** 2, t(L, "(c) 명시야, 초점 어긋남 1.5 λ", "(c) Brightfield, defocus 1.5 λ"), "gray"),
        (np.abs(phase_contrast(U, dx, NA, LAM, dot, a=0.5, alpha=-np.pi / 2)) ** 2, t(L, "(d) 양성 위상차\n(0차를 λ/4 앞당기고 세기 1/4)", "(d) Positive phase contrast\n(0th advanced λ/4, intensity 1/4)"), "gray"),
        (np.abs(phase_contrast(U, dx, NA, LAM, dot, a=0.5, alpha=+np.pi / 2)) ** 2, t(L, "(e) 음성 위상차\n(0차를 λ/4 늦추고 세기 1/4)", "(e) Negative phase contrast\n(0th retarded λ/4, intensity 1/4)"), "gray"),
        (np.abs(dark_field(U, dx, NA, LAM, dot)) ** 2, t(L, "(f) 암시야 (0차 제거)", "(f) Dark field (0th removed)"), "gray"),
    ]
    fig, axes = plt.subplots(1, 6, figsize=(17, 3.6))
    sl = slice(40, 216)
    for i, (ax, (img, ttl, cm)) in enumerate(zip(axes, imgs)):
        kw = dict(vmin=0.5, vmax=1.5) if i in (1, 2) else {}     # 명시야 둘은 같은 눈금 (1 ± 0.5)
        ax.imshow(img[sl, sl], cmap=cm, origin="lower", **kw)
        ax.set_title(ttl, fontsize=9.5, **font(L))
        noticks(ax)
    fig.tight_layout()
    save(fig, L, "fig9-phase-contrast.png")
    bf = imgs[1][0][sl, sl]
    pc = imgs[3][0][sl, sl]
    print(f"    그림9: 명시야 대비 {modulation(bf):.3f}, 양성 위상차 대비 {modulation(pc):.3f}")


FIGS = {2: figure2, 3: figure3, 4: figure4, 5: figure5, 6: figure6, 7: figure7, 9: figure9}


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

"""파동광학 배경이론 6편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/wave-optics-high-na-vector-focusing/      (한국어)
      ../../assets/img/posts/wave-optics-high-na-vector-focusing/en/   (영문)

그림 1 은 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다. 길이 단위는 매질 속 파장 λ 다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from vector_focus import (crossed_polarizer_leakage, fwhm_1d, half_max_area, pupil_field_after_lens, rw_intensity_map,
                          rw_linear, two_beam_fringe_visibility)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "wave-optics-high-na-vector-focusing")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"


def t(L, ko, en):
    s = ko if L["lang"] == "ko" else en
    s = s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$").replace("θ", r"$\theta$").replace("φ", r"$\varphi$")
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


def noticks(ax):
    ax.set_xticks([])
    ax.set_yticks([])


# ------------------------------------------------------------------ 그림 2: 두 평면파의 간섭과 편광

def figure2(L):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3), gridspec_kw=dict(width_ratios=[1, 1.25]))
    a = axes[0]
    th = np.radians(40)
    for sgn, col in [(1, C_RED), (-1, C_RAY)]:
        kx, kz = sgn * np.sin(th), np.cos(th)
        a.annotate("", xy=(0, 0), xytext=(-1.6 * kx, -1.6 * kz), arrowprops=dict(arrowstyle="->", color=col, lw=2))
        # p 편광: 진행 방향에 수직, 그림 평면 안
        mx, mz = -0.9 * kx, -0.9 * kz
        px, pz = np.cos(th), -sgn * np.sin(th)
        a.annotate("", xy=(mx + 0.32 * px, mz + 0.32 * pz), xytext=(mx - 0.32 * px, mz - 0.32 * pz),
                   arrowprops=dict(arrowstyle="<->", color=col, lw=1.2))
        a.plot(-1.3 * kx, -1.3 * kz, "o", mfc="none", mec=col, ms=9)
        a.plot(-1.3 * kx, -1.3 * kz, ".", color=col, ms=4)
    a.axhline(0, color=C_GREY, lw=0.6, ls=":")
    a.axvline(0, color=C_GREY, lw=0.6, ls=":")
    a.text(0.05, -1.55, "z", fontsize=10)
    a.text(1.15, 0.05, "x", fontsize=10)
    a.text(-1.55, -1.75, t(L, "⊙: s 편광 (그림 밖으로)\n↔: p 편광 (그림 평면 안)", "⊙: s polarization (out of page)\n↔: p polarization (in the page)"), fontsize=9, **font(L))
    for sgn in (1, -1):                         # 각 빔이 광축(z)과 이루는 각 θ
        th0, th1 = (270 - 40, 270) if sgn == 1 else (270, 270 + 40)
        a.add_patch(matplotlib.patches.Arc((0, 0), 0.9, 0.9, theta1=th0, theta2=th1, color="k", lw=0.9))
        mid = np.radians((th0 + th1) / 2)
        a.text(0.66 * np.cos(mid) + sgn * 0.02, 0.66 * np.sin(mid), "θ", fontsize=10, ha="center", va="center")
    a.text(0, -1.08, "2θ", fontsize=10, ha="center", va="center")
    a.set_xlim(-1.7, 1.4)
    a.set_ylim(-1.9, 0.4)
    a.set_aspect("equal")
    a.axis("off")
    a.set_title(t(L, "(a) 초점을 향해 ±θ 로 엇갈리는 두 평면파", "(a) Two plane waves crossing at ±θ toward the focus"), fontsize=10.5, **font(L))
    b = axes[1]
    deg = np.linspace(0, 89, 300)
    b.plot(deg, np.ones_like(deg), color=C_RAY, lw=1.6, label=t(L, "s 편광: 가시도 1", "s polarization: visibility 1"))
    b.plot(deg, np.abs(np.cos(2 * np.radians(deg))), color=C_RED, lw=1.6, label=t(L, "p 편광: |cos 2θ|", "p polarization: |cos 2θ|"))
    dn = np.arange(5, 90, 10)
    b.plot(dn, [two_beam_fringe_visibility(np.radians(d), "s") for d in dn], "o", color=C_RAY, mfc="none", ms=5)
    b.plot(dn, [two_beam_fringe_visibility(np.radians(d), "p") for d in dn], "o", color=C_RED, mfc="none", ms=5,
           label=t(L, "원: 벡터 합으로 직접 계산", "circles: direct vector sum"))
    for na in (0.3, 0.71, 0.9):
        d = np.degrees(np.arcsin(na))
        b.axvline(d, color=C_GREY, lw=0.6, ls=":")
        b.text(d + 0.8, {0.3: 0.4, 0.71: 0.62, 0.9: 0.5}[na], f"NA {na:g}", fontsize=8.5, color=C_GREY)
    b.set_xlabel(t(L, "반각 θ (°)", "half-angle θ (°)"), **font(L))
    b.set_ylabel(t(L, "간섭무늬의 가시도", "fringe visibility"), **font(L))
    b.set_ylim(-0.03, 1.12)
    b.set_title(t(L, "(b) 두 빛의 편광이 이루는 각이 간섭을 줄인다\n(θ = 45° 에서 p 편광은 직교해 무늬가 사라진다)", "(b) The angle between polarizations weakens interference\n(at θ = 45° p fields are orthogonal and fringes vanish)"), fontsize=10, **font(L))
    legend(b, L, loc="lower left", size=8.5)
    fig.tight_layout()
    save(fig, L, "fig2-two-beam-polarization.png")


# ------------------------------------------------------------------ 그림 3: 기준 구면 위의 편광

def figure3(L):
    NA = 0.9
    fig = plt.figure(figsize=(16, 4.4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.35, 1, 1, 1], wspace=0.15)
    a = fig.add_subplot(gs[0])
    f = 1.0
    am = np.arcsin(NA)
    a.plot([-0.9, -0.9], [-1, 1], color=C_RAY, lw=3, alpha=0.5)
    a.text(-0.9, 1.07, t(L, "렌즈", "lens"), ha="center", fontsize=9, **font(L))
    arc = np.linspace(-am, am, 100)
    cx = -0.9 + f                                                     # 초점 (기준 구면의 중심)
    a.plot(cx - f * np.cos(arc), f * np.sin(arc), color=C_GREY, lw=1, ls="--")
    a.text(cx - f * np.cos(am) - 0.05, f * np.sin(am) + 0.06, t(L, "기준 구면", "reference sphere"), fontsize=9, color=C_GREY, **font(L))
    for h in (0.0, 0.45, 0.85):
        thh = np.arcsin(h)
        sx, sy = cx - f * np.cos(thh), f * np.sin(thh)
        a.plot([-1.6, -0.9], [h, h], color=C_RED, lw=1)
        a.plot([-0.9, sx, cx], [h, sy, 0], color=C_RED, lw=1)
        # 입력 편광(세로 = x, 그림 평면 안) 과 구면 위 편광(광선에 수직)
        a.annotate("", xy=(-1.35, h + 0.1), xytext=(-1.35, h - 0.1), arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
        ux, uy = np.sin(thh), np.cos(thh)                             # 광선 방향 (cx - sx, -sy) 에 수직
        a.annotate("", xy=(sx + 0.11 * ux, sy + 0.11 * uy), xytext=(sx - 0.11 * ux, sy - 0.11 * uy),
                   arrowprops=dict(arrowstyle="<->", color="k", lw=1.2))
    a.plot(cx, 0, "ko", ms=4)
    a.text(cx + 0.04, -0.12, t(L, "초점", "focus"), fontsize=9, **font(L))
    a.text(-1.62, -0.45, t(L, "↕: 전기장 (x 편광)\n구면 위에서는 광선에\n수직으로 기운다", "↕: electric field (x-polarized)\non the sphere it tilts to stay\nperpendicular to the ray"), fontsize=8.5, **font(L))
    a.set_xlim(-1.7, 0.35)
    a.set_ylim(-0.75, 1.2)
    a.set_aspect("equal")
    a.axis("off")
    a.set_title(t(L, "(a) 무수차 렌즈: 평행광이 기준 구면을 거쳐 초점으로", "(a) Aplanatic lens: rays reach the focus via the reference sphere"), fontsize=10, **font(L))
    s = np.linspace(-1, 1, 241)
    SX, SY = np.meshgrid(s * NA, s * NA)
    rr = np.hypot(SX, SY)
    inside = rr <= NA
    TH = np.arcsin(np.clip(rr, 0, 0.999))
    PH = np.arctan2(SY, SX)
    E = pupil_field_after_lens(TH, PH, "x")
    titles = [t(L, r"(b) 구면 위 장의 x 성분 $|E_{\infty x}|$" + "\n" + r"($\sqrt{\cos\theta}$ 만큼 가장자리가 어둡다)", r"(b) x component $|E_{\infty x}|$ on the sphere" + "\n" + r"(edges dimmed by $\sqrt{\cos\theta}$)"),
              t(L, r"(c) y 성분 $|E_{\infty y}|$: 대각선 방향에 생긴다", r"(c) y component $|E_{\infty y}|$: appears along diagonals"),
              t(L, r"(d) z 성분 $|E_{\infty z}|$: 편광 방향 양쪽", r"(d) z component $|E_{\infty z}|$: on either side along x")]
    for i in range(3):
        ax = fig.add_subplot(gs[i + 1])
        comp = np.where(inside, np.abs(E[i]), np.nan)
        im = ax.imshow(comp, extent=[-NA, NA, -NA, NA], origin="lower", cmap="magma", vmin=0, vmax=1)
        ax.set_title(titles[i], fontsize=9.5, **font(L))
        ax.set_xlabel(r"$\sin\theta\cos\varphi$")
        if i == 0:
            ax.set_ylabel(r"$\sin\theta\sin\varphi$")
        else:
            ax.set_yticks([])
        ax.text(0.62, -0.85, f"max {np.nanmax(comp):.2f}", color="w", fontsize=8.5)
    fig.colorbar(im, ax=fig.axes[1:], fraction=0.02, pad=0.01)
    save(fig, L, "fig3-reference-sphere-polarization.png")


# ------------------------------------------------------------------ 그림 4: 선편광의 초점 (NA 0.9)

def figure4(L):
    NA = 0.9
    x = np.linspace(-1.0, 1.0, 301)
    Ix, Iy, Iz = rw_intensity_map("x", x, NA)
    I = Ix + Iy + Iz
    m = I.max()
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.1))
    labs = [(Ix, r"$|E_x|^2$"), (Iy, r"$|E_y|^2$"), (Iz, r"$|E_z|^2$"), (I, t(L, "전체 세기", "total intensity"))]
    for i, (ax, (img, lab)) in enumerate(zip(axes, labs)):
        ax.imshow(img / m, extent=[x[0], x[-1], x[0], x[-1]], origin="lower", cmap="inferno", vmin=0, vmax=max(img.max() / m, 1e-3))
        ax.set_title(f"({'abcd'[i]}) {lab}  " + t(L, f"(최대 {img.max() / m:.3f})", f"(max {img.max() / m:.3f})"), fontsize=10, **font(L))
        ax.set_xlabel("x / λ")
        if i == 0:
            ax.set_ylabel("y / λ")
    axes[3].contour(x, x, I / m, levels=[0.5], colors="c", linewidths=1.2)
    r_half = 0.514 / NA / 2                                           # 스칼라 에어리 무늬의 반치 반지름
    axes[3].add_patch(mpatches.Circle((0, 0), r_half, fill=False, ec="w", ls="--", lw=1))
    axes[3].text(-0.95, -0.92, t(L, "청록: 반치 등고선\n흰 점선: 스칼라 에어리", "cyan: half maximum\nwhite dashes: scalar Airy"), color="w", fontsize=8.5, **font(L))
    axes[0].annotate("", xy=(0.85, 0.85), xytext=(0.55, 0.85), arrowprops=dict(arrowstyle="<->", color="w", lw=1.4))
    axes[0].text(0.5, 0.68, t(L, "입력 편광", "input pol."), color="w", fontsize=8.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig4-linear-focus-components.png")


# ------------------------------------------------------------------ 그림 5: NA 에 따른 변화

def figure5(L):
    NAs = np.linspace(0.1, 0.95, 18)
    fx, fy, ez = [], [], []
    r = np.linspace(0, 6, 1501)
    for NA in NAs:
        rr = r / NA * 0.25
        Ex, Ey, Ez = rw_linear(rr, np.zeros_like(rr), 0.0, NA)
        Ix = np.abs(Ex) ** 2 + np.abs(Ey) ** 2 + np.abs(Ez) ** 2
        Ex2, Ey2, Ez2 = rw_linear(rr, np.full_like(rr, np.pi / 2), 0.0, NA)
        Iy = np.abs(Ex2) ** 2 + np.abs(Ey2) ** 2 + np.abs(Ez2) ** 2
        fx.append(fwhm_1d(rr, Ix))
        fy.append(fwhm_1d(rr, Iy))
        ez.append((np.abs(Ez) ** 2).max() / Ix[0])
    fx, fy, ez = map(np.array, (fx, fy, ez))
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.2))
    a = axes[0]
    a.plot(NAs, fx * NAs, "o-", color=C_RED, ms=4, label=t(L, "편광 방향 (x)", "along polarization (x)"))
    a.plot(NAs, fy * NAs, "o-", color=C_RAY, ms=4, label=t(L, "편광에 수직 (y)", "across polarization (y)"))
    a.axhline(0.514, color="k", ls="--", lw=1, label=t(L, "스칼라 에어리 0.514", "scalar Airy 0.514"))
    a.set_xlabel("NA")
    a.set_ylabel(t(L, "반치 전폭 × NA / λ", "FWHM × NA / λ"), **font(L))
    a.set_title(t(L, "(a) 초점의 반치 전폭 (선편광, 공기 중)\nNA 를 곱해 스칼라 이론이 상수가 되게 그렸다", "(a) FWHM of the focus (linear polarization, in air)\nmultiplied by NA so the scalar theory is a constant"), fontsize=10, **font(L))
    legend(a, L, loc="upper left")
    b = axes[1]
    b.plot(NAs, ez * 100, "o-", color=C_GRN, ms=4)
    b.set_xlabel("NA")
    b.set_ylabel(t(L, r"$|E_z|^2$ 최대 / 중심 세기 (%)", r"max $|E_z|^2$ / central intensity (%)"), **font(L))
    b.set_title(t(L, "(b) 세로 성분의 세기\n" + r"낮은 NA 에서는 $\mathrm{NA}^2$ 에 비례", "(b) Strength of the longitudinal component\n" + r"$\propto \mathrm{NA}^2$ at low NA"), fontsize=10, **font(L))
    fig.tight_layout()
    save(fig, L, "fig5-na-dependence.png")
    i9 = np.argmin(np.abs(NAs - 0.9))
    print(f"    그림5: NA {NAs[i9]:.2f}: FWHM x {fx[i9]:.3f}, y {fy[i9]:.3f}, Ez {ez[i9]:.3f}")


# ------------------------------------------------------------------ 그림 6: 원·방사·방위 편광

def figure6(L):
    NA = 0.9
    x = np.linspace(-1.0, 1.0, 241)
    fig = plt.figure(figsize=(16.5, 4.2))
    gs = fig.add_gridspec(1, 6, width_ratios=[1, 1, 1, 1, 0.22, 1.45], wspace=0.18)
    items = [("circular", None, t(L, "(a) 원편광: 전체", "(a) Circular: total")),
             ("radial", None, t(L, "(b) 방사 편광: 전체", "(b) Radial: total")),
             ("radial", 2, t(L, r"(c) 방사 편광: $|E_z|^2$ 만", r"(c) Radial: $|E_z|^2$ only")),
             ("azimuthal", None, t(L, "(d) 방위 편광: 전체 (도넛)", "(d) Azimuthal: total (doughnut)"))]
    for i, (kind, comp, title) in enumerate(items):
        ax = fig.add_subplot(gs[i])
        parts = rw_intensity_map(kind, x, NA)
        img = parts[comp] if comp is not None else sum(parts)
        ax.imshow(img / img.max(), extent=[x[0], x[-1], x[0], x[-1]], origin="lower", cmap="inferno")
        ax.contour(x, x, img / img.max(), levels=[0.5], colors="c", linewidths=1)
        ax.set_title(title, fontsize=10, **font(L))
        ax.set_xlabel("x / λ")
        if i == 0:
            ax.set_ylabel("y / λ")
        else:
            ax.set_yticks([])
    b = fig.add_subplot(gs[5])
    inner = np.array([0.0, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.88])
    xs = np.linspace(-1.2, 1.2, 241)
    dx = xs[1] - xs[0]
    ar_r, ar_l = [], []
    for ni in inner:
        ar_r.append(half_max_area(sum(rw_intensity_map("radial", xs, NA, NA_inner=ni)), dx))
        ar_l.append(half_max_area(sum(rw_intensity_map("x", xs, NA, NA_inner=ni)), dx))
    b.plot(inner, ar_r, "o-", color=C_RED, ms=4, label=t(L, "방사 편광", "radial"))
    b.plot(inner, ar_l, "o-", color=C_RAY, ms=4, label=t(L, "선편광", "linear"))
    b.plot([0.85, 0.85], [0.16, 0.26], "k*", ms=11, label=t(L, r"Dorn 등 (0.16, 0.26 $\lambda^2$)", r"Dorn et al. (0.16, 0.26 $\lambda^2$)"))
    b.set_xlabel(t(L, "고리 동공의 안쪽 NA (바깥 NA 0.9)", "inner NA of annular pupil (outer NA 0.9)"), **font(L))
    b.set_ylabel(t(L, r"반치 넓이 ($\lambda^2$)", r"half-maximum area ($\lambda^2$)"), **font(L))
    b.set_title(t(L, "(e) 동공 가운데를 막을수록 방사 편광이 유리해진다", "(e) Blocking the pupil centre favours radial polarization"), fontsize=10, **font(L))
    legend(b, L, loc="upper right")
    save(fig, L, "fig6-cylindrical-polarizations.png")
    print(f"    그림6: 전체 동공 방사 {ar_r[0]:.3f}, 선편광 {ar_l[0]:.3f}; 0.85: 방사 {ar_r[7]:.3f}, 선편광 {ar_l[7]:.3f}")


# ------------------------------------------------------------------ 그림 7: 공기-유리 경계가 만드는 몰타 십자

def figure7(L):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw=dict(width_ratios=[1, 1.3]))
    NA = 0.9
    s = np.linspace(-1, 1, 401)
    SX, SY = np.meshgrid(s * NA, s * NA)
    rr = np.hypot(SX, SY)
    TH = np.arcsin(np.clip(rr, 0, 0.999))
    PH = np.arctan2(SY, SX)
    leak, tot = crossed_polarizer_leakage(TH, PH, 1.5, n_surfaces=2)
    rel = np.where(rr <= NA, leak / tot, np.nan)
    a = axes[0]
    im = a.imshow(np.log10(rel + 1e-9), extent=[-NA, NA, -NA, NA], origin="lower", cmap="viridis", vmin=-6, vmax=-1)
    fig.colorbar(im, ax=a, fraction=0.046, pad=0.03, label=t(L, "log₁₀ (새는 비율)", "log₁₀ (leakage fraction)").replace("₁₀", "10"))
    a.set_xlabel(r"$\sin\theta\cos\varphi$")
    a.set_ylabel(r"$\sin\theta\sin\varphi$")
    a.set_title(t(L, "(a) 직교 편광판 사이의 새는 빛 (NA 0.9)\n공기-유리 경계(n = 1.5) 2면을 지난 광선", "(a) Leakage between crossed polarizers (NA 0.9)\nrays through 2 air–glass surfaces (n = 1.5)"), fontsize=10, **font(L))
    b = axes[1]
    NAs = np.linspace(0.05, 0.95, 60)
    for ns, col in [(2, C_RED), (8, C_RAY)]:
        ext = []
        for na in NAs:
            sm = np.linspace(-1, 1, 201) * na
            X, Y = np.meshgrid(sm, sm)
            R = np.hypot(X, Y)
            T = np.arcsin(np.clip(R, 0, 0.999))
            P = np.arctan2(Y, X)
            lk, tt = crossed_polarizer_leakage(T, P, 1.5, n_surfaces=ns)
            m = R <= na
            ext.append(tt[m].sum() / lk[m].sum())
        b.semilogy(NAs, ext, color=col, lw=1.6, label=t(L, f"공기-유리 경계 {ns}면", f"{ns} air–glass surfaces"))
    b.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    b.set_xlabel("NA")
    b.set_ylabel(t(L, "소광비 (지난 빛 / 새는 빛)", "extinction ratio (transmitted / leaked)"), **font(L))
    b.set_title(t(L, "(b) 소광비는 NA 가 커질수록 급격히 떨어진다\n(코팅 없는 면만 고려)", "(b) Extinction falls steeply with NA\n(uncoated surfaces only)"), fontsize=10, **font(L))
    legend(b, L, loc="upper right")
    fig.tight_layout()
    save(fig, L, "fig7-maltese-cross.png")


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

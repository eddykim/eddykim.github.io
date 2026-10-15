"""파동광학 배경이론 3편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/wave-optics-diffraction/      (한국어)
      ../../assets/img/posts/wave-optics-diffraction/en/   (영문)

그림 1·6·7·9 는 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from diffraction import (angular_spectrum, axial_circular_aperture, axial_ring_sum, axial_zone_plate, circle_pupil,
                         disk_shadow, encircled_energy, far_field, focal_field, propagate_stack, segmented_pupil,
                         slit_array_far_field, transfer_function, zone_edges)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "wave-optics-diffraction")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"
NA_D1, NA_D2 = 0.5889950, 0.5895924          # 나트륨 D 선 (µm)


def t(L, ko, en):
    s = ko if L["lang"] == "ko" else en
    s = s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$").replace("θ", r"$\theta$")
    if "$" in s:                     # 수식이 섞인 문자열은 글꼴 대체가 안 되므로 µ 도 수식으로
        s = s.replace("µ", r"$\mathrm{\mu}$")
    return s


def font(L):
    return {"fontfamily": ["AppleGothic", "DejaVu Sans"]} if L["lang"] == "ko" else {}


def setfonts(L):
    plt.rcParams["font.family"] = ["AppleGothic", "DejaVu Sans"] if L["lang"] == "ko" else ["DejaVu Sans"]


def legend(ax, L, **kw):
    ax.legend(prop=dict(family=plt.rcParams["font.family"], size=kw.pop("size", 8.5)), **kw)


def plain_log_ticks(ax, axis="y"):
    fmt = matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}")
    (ax.yaxis if axis == "y" else ax.xaxis).set_major_formatter(fmt)


def save(fig, L, name, palette=False):
    os.makedirs(L["dir"], exist_ok=True)
    path = os.path.join(L["dir"], name)
    fig.savefig(path, dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    if palette:                      # 색이 연속인 지도 그림은 256색 팔레트로 줄여 파일 크기를 낮춘다
        from PIL import Image
        Image.open(path).convert("RGB").quantize(256, dither=Image.Dither.NONE).save(path, optimize=True)
    print("  저장:", os.path.join(L["dir"], name))


# ------------------------------------------------------------------ 그림 2: 틈 뒤에서 빛이 퍼지는 모습

def figure2(L):
    lam, a = 0.5, 50.0                        # 폭 100 µm 의 틈 (반폭 a), 파장 0.5 µm
    n, dx = 2**17, 0.5                        # 전달함수 표본화 조건 z ≤ N dx²/λ = 65 mm 안
    x = (np.arange(n) - n // 2) * dx
    U0 = (np.abs(x) <= a).astype(complex)
    zs = np.linspace(0.05, 60, 400) * 1000    # µm
    sel = np.abs(x) <= 700
    I = np.abs(propagate_stack(U0, dx, lam, zs)[:, sel]) ** 2
    fig = plt.figure(figsize=(13, 6.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.25, 1], hspace=0.42, wspace=0.28)
    ax = fig.add_subplot(gs[0, :])
    ax.imshow(np.sqrt(I.T), aspect="auto", origin="lower", cmap="inferno",
              extent=[0, zs[-1] / 1000, x[sel][0], x[sel][-1]], vmax=1.3)
    ax.plot([0, 60], [a, a], color="w", ls="--", lw=0.8, alpha=0.7)
    ax.plot([0, 60], [-a, -a], color="w", ls="--", lw=0.8, alpha=0.7)
    NFs = [10, 1, 0.5, 0.1]
    for nf in NFs:
        z = a**2 / (lam * nf)
        ax.axvline(z / 1000, color="c", lw=0.9, ls=":")
        ax.text(z / 1000 + 0.4, 600, f"$N_F$ = {nf:g}", color="c", fontsize=9)
    ax.set_xlabel(t(L, "틈에서 떨어진 거리 z (mm)", "distance from the slit z (mm)"), **font(L))
    ax.set_ylabel("x (µm)")
    ax.set_title(t(L, "(a) 폭 100 µm 틈을 지난 빛 (λ = 0.5 µm, 세기의 제곱근). 흰 점선은 기하학적 그림자의 경계",
                   "(a) Light after a 100 µm slit (λ = 0.5 µm, square root of intensity). White dashes: geometric shadow"),
                 fontsize=10.5, **font(L))
    for i, nf in enumerate(NFs):
        z = a**2 / (lam * nf)
        U = angular_spectrum(U0, dx, lam, z)
        b = fig.add_subplot(gs[1, i])
        half = max(3 * a, 2.2 * lam * z / (2 * a))
        s2 = np.abs(x) <= half
        b.fill_between(x[s2], 0, (np.abs(x[s2]) <= a) * 1.0, color=C_GREY, alpha=0.2, lw=0,
                       label=t(L, "기하광학", "geometric optics"))
        b.plot(x[s2], np.abs(U[s2]) ** 2, color=C_RED, lw=1.2, label=t(L, "회절 계산", "diffraction"))
        if nf <= 0.5:
            th = x[s2] / z
            beta = np.pi * 2 * a * th / lam
            far = (2 * a) ** 2 / (lam * z) * np.sinc(beta / np.pi) ** 2
            b.plot(x[s2], far, color="k", ls="--", lw=0.9, label=t(L, "프라운호퍼 근사", "Fraunhofer approx."))
        b.set_title(f"({'bcde'[i]}) $N_F$ = {nf:g},  z = {z / 1000:g} mm", fontsize=10)
        b.set_xlabel("x (µm)")
        if i == 0:
            b.set_ylabel(t(L, "세기 (입사 = 1)", "intensity (incident = 1)"), **font(L))
        if i in (0, 3):
            legend(b, L, loc="upper right", size=7.5)
    save(fig, L, "fig2-slit-near-to-far.png")


# ------------------------------------------------------------------ 그림 3: 프레넬 존과 아라고 점

def figure3(L):
    fig, axes = plt.subplots(1, 4, figsize=(15.5, 4.1), gridspec_kw=dict(width_ratios=[1, 1.15, 1.25, 1.15]))
    # (a) 위상자 나선: 축 위 한 점(z = 10 µm)에서 본 띠들의 기여를 이어 붙인다
    lam, z = 0.5, 10.0
    k = 2 * np.pi / lam
    edges = zone_edges(z, lam, 60)
    fine = np.concatenate([np.linspace(edges[m], edges[m + 1], 25, endpoint=False) for m in range(60)] + [[edges[-1]]])
    c = np.concatenate([[0], np.cumsum(axial_ring_sum(fine, z, lam))]) * np.exp(-1j * k * z)
    ax = axes[0]
    for m in range(60):
        seg = c[m * 25:(m + 1) * 25 + 1]
        ax.plot(seg.real, seg.imag, color=C_RED if m % 2 == 0 else C_RAY, lw=1.0)
    full = np.exp(0j)   # 구멍 없이 다 열면 평면파 그대로 = 1
    ax.annotate("", xy=(c[25].real, c[25].imag), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color="k", lw=1.6))
    ax.annotate("", xy=(full.real, full.imag), xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=C_GRN, lw=1.6))
    box = dict(fc="white", ec="none", pad=1)
    ax.text(1.5, -0.16, t(L, "첫 존만", "zone 1 only"), fontsize=9, ha="center", bbox=box, **font(L))
    ax.text(0.5, -0.16, t(L, "다 열면", "fully open"), fontsize=9, ha="center", color=C_GRN, bbox=box, **font(L))
    ax.set_aspect("equal")
    ax.set_title(t(L, "(a) 존의 기여를 이어 붙인 위상자\n(빨강·파랑이 번갈아 한 존)", "(a) Phasor sum over the zones\n(red and blue alternate per zone)"), fontsize=10, **font(L))
    ax.set_xlabel(t(L, "실수부", "real part"), **font(L))
    ax.set_ylabel(t(L, "허수부", "imaginary part"), **font(L))
    # (b) 원형 구멍 뒤 축 위 세기와 N_F
    lam, a = 0.5, 1000.0
    NF = np.linspace(0.15, 8.5, 700)
    I = np.array([axial_circular_aperture(a, a**2 / (lam * nf), lam) for nf in NF])
    b = axes[1]
    b.plot(NF, I, color=C_RED, lw=1.5)
    b.axhline(1, color=C_GREY, ls="--", lw=0.8)
    b.text(8.4, 1.25, t(L, "구멍이 없을 때", "no aperture"), color=C_GREY, fontsize=9, ha="right",
           bbox=dict(fc="white", ec="none", pad=1), **font(L))
    for m in range(1, 9):
        b.text(m, 4.15 if m % 2 else 0.2, str(m), ha="center", fontsize=8.5, color="k")
    b.set_xlabel(t(L, "프레넬 수 $N_F = a^2/\\lambda z$ (= 열린 존의 수)", "Fresnel number $N_F = a^2/\\lambda z$ (= number of open zones)"), **font(L))
    b.set_ylabel(t(L, "축 위 세기 (입사 = 1)", "on-axis intensity (incident = 1)"), **font(L))
    b.set_ylim(0, 4.6)
    b.set_title(t(L, "(b) 원형 구멍 뒤 축 위의 세기\n(반지름 1 mm, λ = 0.5 µm, z = 0.24–13 m)", "(b) On-axis intensity behind a circular hole\n(radius 1 mm, λ = 0.5 µm, z = 0.24–13 m)"), fontsize=10, **font(L))
    # (c) 원판 그림자: 그림1 사진의 조건 (지름 5.8 mm, 광원까지 1.53 m, 화면까지 1.83 m)
    lam, a, z1, z2 = 0.55, 2900.0, 1.53e6, 1.83e6
    zeff, M = z1 * z2 / (z1 + z2), (z1 + z2) / z1          # 점광원의 구면파 = 평면파 + 거리 zeff, 화면 좌표는 M 배
    r = np.concatenate([np.linspace(0, 400, 401), np.linspace(402, 4200, 1900)])
    I = np.abs(disk_shadow(a, r, zeff, lam, n=40001)) ** 2
    cc = axes[2]
    R = r * M / 1000
    cc.plot(np.concatenate([-R[::-1], R]), np.concatenate([I[::-1], I]), color=C_RED, lw=0.9)
    cc.axvspan(-a * M / 1000, a * M / 1000, color=C_GREY, alpha=0.15, lw=0)
    cc.set_xlim(-8.5, 8.5)
    cc.set_ylim(0, 1.45)
    cc.set_xlabel(t(L, "화면 위 위치 (mm)", "position on the screen (mm)"), **font(L))
    cc.set_ylabel(t(L, "세기 (원판 없을 때 = 1)", "intensity (no disk = 1)"), **font(L))
    cc.set_title(t(L, f"(c) 원판 그림자 (사진과 같은 배치, 550 nm)\n$N_F$ ≈ {a**2 / (lam * zeff):.0f}, 회색은 기하학적 그림자", f"(c) Shadow of the disk (same geometry as the photograph, 550 nm)\n$N_F$ ≈ {a**2 / (lam * zeff):.0f}, grey: geometric shadow"), fontsize=10, **font(L))
    ins = cc.inset_axes([0.33, 0.52, 0.34, 0.4])
    s = r <= 400
    ins.plot(np.concatenate([-R[s][::-1], R[s]]), np.concatenate([I[s][::-1], I[s]]), color=C_RED, lw=1)
    ins.set_xlim(-0.6, 0.6)
    ins.tick_params(labelsize=7)
    ins.set_title(t(L, "가운데 확대", "centre, zoomed"), fontsize=8, **font(L))
    # (d) 존 플레이트의 축 위 세기
    lam, f, n_open = 0.5, 1e5, 20
    z = np.linspace(0.12, 1.35, 1400) * f
    Izp = axial_zone_plate(f, lam, n_open, z, n=60001)
    d = axes[3]
    d.plot(z / f, Izp, color=C_RAY, lw=1.2)
    d.set_yscale("log")
    d.set_ylim(1, 8000)
    for zz, lab, hy in [(1, "f", 1.3), (1 / 3, "f/3", 1.3), (1 / 5, "f/5", 1.3), (1 / 7, "f/7", 2.4)]:
        i = np.argmin(np.abs(z / f - zz))
        j = i - 20 + np.argmax(Izp[i - 20:i + 20])
        d.text(z[j] / f, Izp[j] * hy, lab, ha="center", fontsize=8.5)
    plain_log_ticks(d)
    d.set_xlabel(t(L, "거리 z / f", "distance z / f"), **font(L))
    d.set_ylabel(t(L, "축 위 세기 (입사 = 1)", "on-axis intensity (incident = 1)"), **font(L))
    d.set_title(t(L, "(d) 존 플레이트 (홀수 존 20개만 열림)\n$f = r_1^2/\\lambda$ = 10 cm", "(d) Zone plate (20 odd zones open)\n$f = r_1^2/\\lambda$ = 10 cm"), fontsize=10, **font(L))
    fig.tight_layout()
    save(fig, L, "fig3-fresnel-zones-arago.png")
    print(f"    그림3: 아라고 점 중심 세기 {I[0]:.3f}, 주 초점 {Izp.max():.0f}")


# ------------------------------------------------------------------ 그림 4: 직선 가장자리

def figure4(L):
    lam = 0.5
    n, dx = 2**18, 0.1
    x = (np.arange(n) - n // 2) * dx
    Lopen, T = 10000.0, 4000.0
    tr = np.where((x >= 0) & (x < Lopen), np.where(x < Lopen - T, 1, 0.5 * (1 + np.cos(np.pi * (x - (Lopen - T)) / T))), 0)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.1), gridspec_kw=dict(width_ratios=[1.2, 1.2, 0.9]))
    cols = [C_RAY, C_RED, C_GRN]
    zs = [200.0, 1000.0, 5000.0]
    for z, col in zip(zs, cols):
        I = np.abs(angular_spectrum(tr.astype(complex), dx, lam, z)) ** 2
        s = (x > -120) & (x < 300)
        axes[0].plot(x[s], I[s], color=col, lw=1.1, label=f"z = {z / 1000:g} mm")
        w = x * np.sqrt(2 / (lam * z))
        s2 = (w > -3) & (w < 6)
        axes[1].plot(w[s2], I[s2], color=col, lw=2.2 if z == zs[0] else 1.1, alpha=0.6 if z == zs[0] else 1,
                     label=f"z = {z / 1000:g} mm")
    for ax in axes[:2]:
        ax.axvline(0, color="k", lw=0.7)
        ax.axhline(0.25, color=C_GREY, ls=":", lw=0.8)
        ax.axhline(1, color=C_GREY, ls="--", lw=0.8)
        ax.set_ylabel(t(L, "세기 (입사 = 1)", "intensity (incident = 1)"), **font(L))
        legend(ax, L, loc="lower right", size=8)
    axes[0].set_xlabel(t(L, "x (µm), 0 이 기하학적 그림자의 경계, 음수 쪽이 그림자", "x (µm); 0 is the geometric shadow edge, negative side is shadow"), **font(L))
    axes[0].set_title(t(L, "(a) 칼날 가장자리 뒤의 세기 (λ = 0.5 µm)", "(a) Intensity behind a straight edge (λ = 0.5 µm)"), fontsize=10.5, **font(L))
    axes[1].set_xlabel(r"$w = x\sqrt{2/\lambda z}$")
    axes[1].set_title(t(L, r"(b) 가로축을 $\sqrt{\lambda z}$ 로 나누면 한 곡선이 된다", r"(b) Scaled by $\sqrt{\lambda z}$, all distances collapse"), fontsize=10.5, **font(L))
    axes[1].annotate("1/4", xy=(0, 0.25), xytext=(-2.3, 0.45), arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
    axes[1].annotate("1.37", xy=(1.22, 1.37), xytext=(2.4, 1.3), arrowprops=dict(arrowstyle="->", lw=0.8), fontsize=9)
    # (c) 코르뉴 나선: e^{iπt²/2} 의 누적 적분을 직접 쌓는다
    tt = np.linspace(-8, 8, 40001)
    ph = np.exp(1j * np.pi * tt**2 / 2)
    cum = np.concatenate([[0], np.cumsum((ph[1:] + ph[:-1]) / 2 * np.diff(tt))])
    cum -= cum[np.argmin(np.abs(tt))]
    c = axes[2]
    c.plot(cum.real, cum.imag, color="k", lw=0.8)
    eye = cum[-1]
    for w0, col, lw, al in [(1.22, C_RED, 3.2, 0.45), (0.0, C_RAY, 1.6, 1), (-1.0, C_GRN, 1.6, 1)]:
        i0 = np.argmin(np.abs(tt + w0))                 # 관측점 w 의 진폭은 −w 에서 눈까지의 화살표
        c.annotate("", xy=(eye.real, eye.imag), xytext=(cum[i0].real, cum[i0].imag),
                   arrowprops=dict(arrowstyle="->", color=col, lw=lw, alpha=al))
        c.plot(cum[i0].real, cum[i0].imag, "o", color=col, ms=4)
    c.set_aspect("equal")
    c.set_title(t(L, "(c) 코르뉴 나선. 화살표 길이가 진폭\n파랑: 경계, 빨강: 첫 봉우리, 초록: 그림자 속", "(c) Cornu spiral; arrow length = amplitude\nblue: edge, red: first maximum, green: in the shadow"), fontsize=9.5, **font(L))
    c.set_xlabel(t(L, "실수부", "real part"), **font(L))
    c.set_ylabel(t(L, "허수부", "imaginary part"), **font(L))
    fig.tight_layout()
    save(fig, L, "fig4-straight-edge.png")


# ------------------------------------------------------------------ 그림 5: 구멍 모양과 원거리 무늬

def webb_pupil(n, strut_px=2.5):
    """그림 6 사진과 같은 배치: 육각 조각의 윗변·아랫변이 수직(조각 배열을 90° 돌림), 받침대는 위로 하나, 아래 양쪽으로 둘."""
    p = segmented_pupil(n, 17.5, gap_px=2.0, strut_px=0).T.astype(bool)
    y, x = np.mgrid[:n, :n] - (n - 1) / 2
    r = np.hypot(x, y)
    for ang in (90, -120, -60):
        a = np.radians(ang)
        along = x * np.cos(a) + y * np.sin(a)
        across = -x * np.sin(a) + y * np.cos(a)
        p &= ~((along > 0) & (np.abs(across) < strut_px / 2) & (r < n / 2))
    return p.astype(float)


def figure5(L):
    n, pad = 256, 4
    y, x = np.mgrid[:n, :n] - (n - 1) / 2
    shapes = [
        (((np.abs(x) <= 12) & (np.abs(y) <= 70)).astype(float), t(L, "(a) 가늘고 긴 직사각형", "(a) Narrow rectangle")),
        (((np.abs(x) <= 50) & (np.abs(y) <= 50)).astype(float), t(L, "(b) 정사각형", "(b) Square")),
        (circle_pupil(n, 60), t(L, "(c) 원", "(c) Circle")),
        (webb_pupil(n), t(L, "(d) 육각 조각 18장 + 받침대 3개", "(d) 18 hexagonal segments + 3 struts")),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(14, 7.2))
    for i, (p, title) in enumerate(shapes):
        axes[0, i].imshow(p, cmap="gray", origin="lower")
        axes[0, i].set_title(title, fontsize=10.5, **font(L))
        I = far_field(p, pad)
        I /= I.max()
        c = I.shape[0] // 2
        h = 150
        axes[1, i].imshow(np.log10(I[c - h:c + h, c - h:c + h] + 1e-7), cmap="inferno", vmin=-6, vmax=0, origin="lower")
        for a in axes[:, i]:
            a.set_xticks([])
            a.set_yticks([])
    axes[0, 0].set_ylabel(t(L, "구멍", "aperture"), fontsize=11, **font(L))
    axes[1, 0].set_ylabel(t(L, "원거리 세기 (로그, 6자리)", "far-field intensity (log, 6 decades)"), fontsize=11, **font(L))
    fig.tight_layout()
    save(fig, L, "fig5-aperture-far-field.png", palette=True)


# ------------------------------------------------------------------ 그림 8: 에어리 무늬, 분해능, 초점 근처

def figure8(L):
    # 에어리 단면은 초점면의 반지름 적분(focal_field, u = 0)으로 촘촘하게 구한다. v = π D sinθ/λ
    v = np.linspace(1e-4, 14 * np.pi, 7001)
    A = focal_field(v, [0.0])[0]
    I_air = np.abs(A) ** 2
    s_ = np.linspace(0, 4, 4001)
    I_slit = np.abs(slit_array_far_field(s_, 1, 0, 1.0, 1.0, pts_per_slit=400)) ** 2      # 폭 D = 1, λ = 1
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.1), gridspec_kw=dict(width_ratios=[1.15, 1, 1.15, 1.25]))
    a = axes[0]
    a.semilogy(v / np.pi, I_air + 1e-9, color=C_RED, lw=1.4, label=t(L, "원 (지름 D)", "circle (diameter D)"))
    a.semilogy(s_, I_slit + 1e-9, color=C_RAY, lw=1.0, ls="--", label=t(L, "틈 (폭 D)", "slit (width D)"))
    a.set_xlim(0, 4)
    a.set_ylim(1e-4, 1.5)
    plain_log_ticks(a)
    for zz in (1.22, 2.23, 3.24):
        a.axvline(zz, color=C_RED, lw=0.5, ls=":")
    a.set_xlabel(r"$\sin\theta\,/\,(\lambda/D)$")
    a.set_ylabel(t(L, "세기 (중심 = 1)", "intensity (centre = 1)"), **font(L))
    a.set_title(t(L, "(a) 에어리 무늬와 틈의 무늬\n첫 영점: 원 1.22, 틈 1", "(a) Airy pattern vs slit pattern\nfirst zero: circle 1.22, slit 1"), fontsize=10, **font(L))
    legend(a, L, loc="upper right")
    # (b) 원 안 에너지: 세기를 넓이로 적분한다. 전체는 동공을 지난 에너지와 같다 (파세발): 정규화하면 ∫ I v dv 의 끝값이 2
    b = axes[1]
    ee = np.concatenate([[0], np.cumsum((I_air[1:] * v[1:] + I_air[:-1] * v[:-1]) / 2 * np.diff(v))]) / 2
    b.plot(v / np.pi, ee * 100, color="k", lw=1.4)
    for zz in (3.8317, 7.0156, 10.1735):
        e = np.interp(zz, v, ee)
        b.plot(zz / np.pi, e * 100, "o", color=C_RED, ms=5)
        b.text(zz / np.pi + 0.08, e * 100 - 8, f"{e * 100:.1f} %", fontsize=9)
    b.set_xlim(0, 4)
    b.set_xlabel(t(L, "반지름 r / (λf/D)", "radius r / (λf/D)"), **font(L))
    b.set_ylabel(t(L, "원 안의 에너지 (%)", "encircled energy (%)"), **font(L))
    b.set_ylim(0, 102)
    b.set_title(t(L, "(b) 반지름 안에 든 에너지\n(점은 1, 2, 3번째 어두운 고리)", "(b) Encircled energy\n(dots: 1st, 2nd, 3rd dark rings)"), fontsize=10, **font(L))
    # (c) 레일리 간격의 두 점 (서로 결맞지 않음: 세기로 더한다)
    sep = 1.2197
    g = np.linspace(-2.4, 2.4, 481)
    X, Y = np.meshgrid(g, g)
    airy = lambda r: np.interp(np.pi * r, v, I_air, right=0)          # noqa: E731
    Itot = airy(np.hypot(X - sep / 2, Y)) + airy(np.hypot(X + sep / 2, Y))
    cc = axes[2]
    cc.imshow(Itot, cmap="gray", extent=[g[0], g[-1], g[0], g[-1]], origin="lower", vmax=1.0)
    line = Itot[240]
    cc.plot(g, -2.4 + 2.0 * line, color=C_KEY, lw=1.3)
    dip = line[240] / line.max()
    cc.text(0, -2.4 + 2.0 * line[240] - 0.35, f"{dip:.3f}", color=C_KEY, ha="center", fontsize=9)
    cc.set_xlabel(r"$x\,/\,(\lambda f/D)$")
    cc.set_title(t(L, "(c) 레일리 간격(1.22 λf/D)의 두 점\n주황: 가로 단면, 숫자는 골의 높이", "(c) Two points at the Rayleigh spacing (1.22 λf/D)\norange: horizontal cut, number: dip height"), fontsize=10, **font(L))
    # (d) 초점 근처의 세기 (광학 좌표)
    u = np.linspace(-32, 32, 321)
    v = np.linspace(0, 14, 141)
    If = np.abs(focal_field(v, u)) ** 2
    full = np.concatenate([If[:, ::-1], If[:, 1:]], axis=1)
    d = axes[3]
    d.imshow(np.log10(full.T + 1e-6), extent=[u[0], u[-1], -v[-1], v[-1]], aspect="auto", cmap="inferno", vmin=-3.5, vmax=0, origin="lower")
    d.plot(u, np.abs(u), color="c", lw=0.8, ls="--")
    d.plot(u, -np.abs(u), color="c", lw=0.8, ls="--")
    d.set_ylim(-14, 14)
    for uz in (-4 * np.pi, 4 * np.pi):
        d.plot(uz, 0, "o", mfc="none", mec="w", ms=6)
    d.set_xlabel(t(L, r"축 방향 $u = k\,\mathrm{NA}^2 z$", r"axial $u = k\,\mathrm{NA}^2 z$"), **font(L))
    d.set_ylabel(t(L, r"가로 $v = k\,\mathrm{NA}\,r$", r"lateral $v = k\,\mathrm{NA}\,r$"), **font(L))
    d.set_title(t(L, "(d) 초점 근처의 세기 (로그). 점선: 기하광학의 원뿔\n흰 원: 축 위 첫 영점 u = ±4π", "(d) Intensity near focus (log); dashes: geometric cone\nwhite circles: first axial zeros u = ±4π"), fontsize=10, **font(L))
    fig.tight_layout()
    save(fig, L, "fig8-airy-resolution-focus.png")
    print(f"    그림8: 레일리 골 {dip:.3f}")


# ------------------------------------------------------------------ 그림 10: 격자

def figure10(L):
    lam, d = 0.5, 2.0
    b = 0.3 * d
    s = np.linspace(0, 0.9, 9001)
    fig = plt.figure(figsize=(14, 5.2))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.35, 1], hspace=0.15, wspace=0.22)
    env = np.abs(slit_array_far_field(s, 1, d, b, lam)) ** 2
    for i, N in enumerate([2, 5, 20]):
        a = fig.add_subplot(gs[i, 0])
        I = np.abs(slit_array_far_field(s, N, d, b, lam)) ** 2
        a.plot(s * d / lam, I, color=C_RED, lw=1.0)
        a.plot(s * d / lam, env, color="k", ls="--", lw=0.8)
        a.text(0.985, 0.8, t(L, f"틈 {N}개", f"{N} slits"), transform=a.transAxes, ha="right", fontsize=10, **font(L))
        a.set_ylim(0, 1.08)
        a.set_xlim(0, 3.6)
        if i < 2:
            a.set_xticklabels([])
        if i == 0:
            a.set_title(t(L, "(a) 틈 N개의 원거리 세기 (주봉우리를 1로). 점선은 틈 하나의 포락선",
                          "(a) Far field of N slits (principal peaks = 1); dashes: single-slit envelope"), fontsize=10.5, **font(L))
        if i == 1:
            a.set_ylabel(t(L, "세기", "intensity"), **font(L))
    a.set_xlabel(r"$d\,\sin\theta\,/\,\lambda$ " + t(L, "(정수가 회절 차수)", "(integers are diffraction orders)"), **font(L))
    # (b) 나트륨 D 선, 600 줄/mm 격자의 2차
    d2 = 1 / 0.6
    ax = fig.add_subplot(gs[:, 1])
    lam_axis = np.linspace(588.4, 590.2, 3001)                      # 그 각도에 2차로 오는 파장 (nm)
    ss = 2 * lam_axis / 1000 / d2
    for N, col, off in [(400, C_RAY, 0.0), (600, C_RED, 1.25)]:
        I1 = np.abs(slit_array_far_field(ss, N, d2, 0.3 * d2, NA_D1, 4)) ** 2
        I2 = np.abs(slit_array_far_field(ss, N, d2, 0.3 * d2, NA_D2, 4)) ** 2
        top = max(I1.max(), I2.max())
        ax.plot(lam_axis, off + I1 / top, color=col, lw=0.7, ls=":")
        ax.plot(lam_axis, off + I2 / top, color=col, lw=0.7, ls=":")
        ax.plot(lam_axis, off + (I1 + I2) / top, color=col, lw=1.6, label=t(L, f"{N}줄 (mN = {2 * N})", f"{N} lines (mN = {2 * N})"))
    ax.set_xlabel(t(L, "2차로 그 각도에 오는 파장 (nm)", "wavelength arriving at that angle in 2nd order (nm)"), **font(L))
    ax.set_yticks([])
    ax.set_title(t(L, "(b) 나트륨 D 선(0.6 nm 간격), 600 줄/mm 격자의 2차\n가르려면 mN ≥ λ/Δλ ≈ 987", "(b) Sodium D lines (0.6 nm apart), 600 lines/mm, 2nd order\nresolving needs mN ≥ λ/Δλ ≈ 987"), fontsize=10.5, **font(L))
    ax.set_xticks([588.5, 589.0, 589.5, 590.0])
    ax.set_ylim(-0.05, 2.7)
    legend(ax, L, loc="upper center", ncol=2, size=8.5)
    save(fig, L, "fig10-grating.png")


# ------------------------------------------------------------------ 그림 11: 평면파로 쪼개기

def figure11(L):
    lam = 0.5
    n, dx = 720, lam / 60                       # 창 12λ: 1.5λ 와 0.8λ 주기가 정수 번 들어간다
    x = (np.arange(n) - n // 2) * dx
    zs = np.linspace(0, 3 * lam, 181)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.1]))
    for ax, ratio, lab in [(axes[0], 1.5, "a"), (axes[1], 0.8, "b")]:
        Lam = ratio * lam
        # 격자 주기를 창에 맞춘다 (주기 경계 조건)
        nper = round(n * dx / Lam)
        Lam_ = n * dx / nper
        U0 = np.exp(1j * 2 * np.pi * x / Lam_)
        U = propagate_stack(U0, dx, lam, zs)
        sel = np.abs(x) <= 2.5 * lam
        ax.imshow(U[:, sel].real, extent=[x[sel][0] / lam, x[sel][-1] / lam, 0, zs[-1] / lam], origin="lower",
                  cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal")
        if ratio > 1:
            th = np.degrees(np.arcsin(lam / Lam_))
            ttl = t(L, f"({lab}) 주기 Λ = {ratio:g}λ: 각도 {th:.1f}° 로 진행", f"({lab}) Period Λ = {ratio:g}λ: travels at {th:.1f}°")
        else:
            dec = 1 / (2 * np.pi * np.sqrt(1 / Lam_**2 - 1 / lam**2))
            ttl = t(L, f"({lab}) 주기 Λ = {ratio:g}λ: {dec / lam:.2f}λ 마다 1/e 로 감쇠", f"({lab}) Period Λ = {ratio:g}λ: decays by 1/e every {dec / lam:.2f}λ")
        ax.set_title(ttl, fontsize=10.5, **font(L))
        ax.set_xlabel("x / λ")
        ax.set_ylabel("z / λ")
    f = np.linspace(-2.2, 2.2, 2001) / lam
    c = axes[2]
    for z, col in [(0.25 * lam, C_GRN), (0.5 * lam, C_RAY), (2 * lam, C_RED)]:
        c.plot(f * lam, np.abs(transfer_function(f, lam, z)), color=col, lw=1.4, label=f"z = {z / lam:g}λ")
    c.axvspan(-1, 1, color=C_GREY, alpha=0.12, lw=0)
    c.text(0, 0.5, t(L, "진행하는 평면파\n(|H| = 1)", "propagating\n(|H| = 1)"), ha="center", fontsize=9, **font(L))
    c.text(1.6, 0.85, t(L, "소멸파", "evanescent"), ha="center", fontsize=9, **font(L))
    c.set_xlabel(t(L, "공간 주파수 × λ (= λ/Λ)", "spatial frequency × λ (= λ/Λ)"), **font(L))
    c.set_ylabel(t(L, "전달함수의 크기 |H|", "transfer function magnitude |H|"), **font(L))
    c.set_title(t(L, "(c) 자유 공간의 전달함수", "(c) Transfer function of free space"), fontsize=10.5, **font(L))
    legend(c, L, loc="upper left", size=8)
    fig.tight_layout()
    save(fig, L, "fig11-angular-spectrum.png")


FIGS = {2: figure2, 3: figure3, 4: figure4, 5: figure5, 8: figure8, 10: figure10, 11: figure11}


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

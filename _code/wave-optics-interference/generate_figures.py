"""파동광학 배경이론 1편 그림 생성 (한국어판·영문판).

실행: python generate_figures.py [그림번호 ...]
출력: ../../assets/img/posts/wave-optics-interference/      (한국어)
      ../../assets/img/posts/wave-optics-interference/en/   (영문)

그림 1·4·6·7·10 은 위키미디어 공용 사진(ext-*)이라 여기서 만들지 않는다.
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from interference import (etalon_transmission, interferogram, michelson_outputs, multibeam_sum, newton_gap,
                          point_sources_intensity, superpose, thin_layer_reflectance, young_screen)

plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["mathtext.fontset"] = "dejavusans"
BASE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                        "assets", "img", "posts", "wave-optics-interference")
C_RAY, C_RED, C_KEY, C_GREY, C_GRN, C_PUR = "#1f77b4", "#d62728", "#ff7f0e", "#7f7f7f", "#2ca02c", "#9467bd"


def t(L, ko, en):
    """언어에 맞는 문구. AppleGothic 의 °·λ 는 모양이 어긋나므로 수식 기호로 바꾼다."""
    s = ko if L["lang"] == "ko" else en
    return s.replace("°", r"$^\circ$").replace("λ", r"$\lambda$")


def font(L):
    return {"fontfamily": ["AppleGothic", "DejaVu Sans"]} if L["lang"] == "ko" else {}


def setfonts(L):
    plt.rcParams["font.family"] = ["AppleGothic", "DejaVu Sans"] if L["lang"] == "ko" else ["DejaVu Sans"]


def legend(ax, L, **kw):
    ax.legend(prop=dict(family=plt.rcParams["font.family"], size=kw.pop("size", 8.5)), **kw)


def save(fig, L, name):
    os.makedirs(L["dir"], exist_ok=True)
    fig.savefig(os.path.join(L["dir"], name), dpi=150, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print("  저장:", os.path.join(L["dir"], name))


# ------------------------------------------------------------------ 그림 2: 진폭을 더하고 세기를 잰다

def figure2(L):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.0), gridspec_kw=dict(width_ratios=[0.8, 1.2, 1]))
    a, b, c = axes
    # (a) 위상자 합: 세 경우를 위아래로 나누어 그린다
    A1, A2 = 1.0, 0.7
    for row, (dl, col, lab) in enumerate([(0.0, C_RED, r"$\delta = 0$"), (np.pi / 2, C_GRN, r"$\delta = \pi/2$"),
                                          (np.pi, C_RAY, r"$\delta = \pi$")]):
        y0 = -1.1 * row
        tip = A1 + A2 * np.exp(1j * dl)
        a.annotate("", (A1, y0), (0, y0), arrowprops=dict(arrowstyle="->", color="k", lw=1.6))
        a.annotate("", (tip.real, y0 + tip.imag), (A1, y0), arrowprops=dict(arrowstyle="->", color=col, lw=1.6))
        a.annotate("", (tip.real, y0 + tip.imag), (0, y0), arrowprops=dict(arrowstyle="->", color=col, lw=1.2, ls="--"))
        a.text(-0.12, y0 + 0.2, lab, color=col, fontsize=10.5, ha="left")
        a.text(2.75, y0 + 0.25, t(L, f"세기 {abs(tip) ** 2:.2f}", f"intensity {abs(tip) ** 2:.2f}"), color=col, fontsize=9.5,
               ha="right", va="center", **font(L))
    a.set_xlim(-0.15, 2.8); a.set_ylim(-2.45, 0.95); a.set_aspect("equal")
    a.axis("off")
    a.set_title(t(L, "(a) 위상자로 본 두 진폭의 합\n(검정: 진폭 1, 색: 진폭 0.7, 점선: 합)",
                  "(a) Adding two amplitudes as phasors\n(black: amplitude 1, colour: 0.7, dashed: sum)"),
                fontsize=10.5, **font(L))
    # (b) 위상차에 따른 세기
    d = np.linspace(-2 * np.pi, 2 * np.pi, 801)
    for I2, col in [(1.0, C_RAY), (0.25, C_KEY), (0.01, C_GRN)]:
        b.plot(d / np.pi, superpose(1.0, I2, d), color=col, lw=2,
               label=t(L, f"I₂/I₁ = {I2:g}", f"I₂/I₁ = {I2:g}").replace("I₂/I₁", r"$I_2/I_1$"))
    b.axhline(2.0, color=C_RAY, ls=":", lw=1, label=t(L, "평균 = I₁ + I₂ (같은 세기)", "mean = I₁ + I₂ (equal beams)").replace("I₁ + I₂", r"$I_1 + I_2$"))
    b.set_xlabel(t(L, "위상차 δ (π 단위)", "phase difference δ (units of π)").replace("δ", r"$\delta$").replace("π", r"$\pi$"), **font(L))
    b.set_ylabel(t(L, "세기 (I₁ = 1)", "intensity (I₁ = 1)").replace("I₁", r"$I_1$"), **font(L))
    legend(b, L, loc="upper right")
    b.set_ylim(0, 5.3)
    b.set_title(t(L, "(b) 위상차가 정하는 세기", "(b) Intensity set by the phase difference"), fontsize=10.5, **font(L))
    # (c) 가시도
    ratio = np.logspace(-3, 0, 200)
    V = 2 * np.sqrt(ratio) / (1 + ratio)
    c.semilogx(ratio, V, color="k", lw=2)
    for r_, col in [(1.0, C_RAY), (0.25, C_KEY), (0.01, C_GRN)]:
        c.plot([r_], [2 * np.sqrt(r_) / (1 + r_)], "o", color=col, ms=7)
    c.set_xlabel(t(L, "세기 비 I₂/I₁", "intensity ratio I₂/I₁").replace("I₂/I₁", r"$I_2/I_1$"), **font(L))
    c.set_ylabel(t(L, "가시도 V", "visibility V"), **font(L))
    c.set_ylim(0, 1.05)
    c.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    c.text(0.0013, 0.9, r"$V = \frac{2\sqrt{I_1 I_2}}{I_1 + I_2}$", fontsize=13)
    c.set_title(t(L, "(c) 약한 빔도 무늬를 잘 만든다", "(c) Even a weak beam makes clear fringes"), fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig2-two-wave-superposition.png")


# ------------------------------------------------------------------ 그림 3: 두 점광원

def figure3(L):
    lam = 1.0
    d = 6.0
    fig, (a, b) = plt.subplots(1, 2, figsize=(12.5, 4.8), gridspec_kw=dict(width_ratios=[1.35, 1]))
    x = np.linspace(-25, 25, 900)
    z = np.linspace(0.3, 30, 700)
    X, Z = np.meshgrid(x, z)
    I = point_sources_intensity(X, Z, [(-d / 2, 0.0), (d / 2, 0.0)], lam)
    I_n = I * (np.hypot(X, Z) ** 2)           # 1/r² 감쇠를 빼고 무늬만 보이게
    a.imshow(I_n, extent=[x[0], x[-1], z[0], z[-1]], origin="lower", cmap="inferno", vmax=4, aspect="equal")
    for m in range(-5, 6):
        # r1 − r2 = mλ 인 쌍곡선 (점광원이 초점)
        if abs(m) * lam >= d:
            continue
        aa = abs(m) * lam / 2
        bb = np.sqrt((d / 2) ** 2 - aa**2)
        zz = np.linspace(0, 30, 300)
        xx = aa * np.sqrt(1 + (zz / bb) ** 2) * np.sign(m) if m else 0 * zz
        a.plot(xx, zz, color="w", lw=0.9, ls="--", alpha=0.9)
    a.plot([-d / 2, d / 2], [0.3, 0.3], "o", color=C_GRN, ms=6)
    a.set_xlim(-25, 25)
    a.set_xlabel(t(L, "x (파장 단위)", "x (wavelengths)"), **font(L))
    a.set_ylabel(t(L, "z (파장 단위)", "z (wavelengths)"), **font(L))
    a.set_title(t(L, "(a) 간격 6λ 인 두 점광원 (거리에 따른 감쇠를 뺀 세기)\n흰 점선: 두 광원까지의 거리 차가 mλ 인 쌍곡선", "(a) Two point sources 6λ apart (intensity with 1/r² removed)\nwhite dashed: hyperbolas of path difference mλ"),
                fontsize=10.5, **font(L))
    # (b) 먼 스크린: 이중 슬릿 크기 (d = 0.5 mm, L = 1 m, λ = 0.55 µm)
    dd, LL, lam2 = 500.0, 1.0e6, 0.55
    xs = np.linspace(-4000, 4000, 4001)
    Ie = young_screen(xs, dd, LL, lam2)
    Ie /= Ie.max()
    b.plot(xs / 1000, Ie, color=C_KEY, lw=1.8)
    p = lam2 * LL / dd
    b.annotate("", (p / 1000, 1.06), (0, 1.06), arrowprops=dict(arrowstyle="<->", lw=1))
    b.text(p / 2000, 1.09, f"{p / 1000:.2f} mm", ha="center", fontsize=9.5)
    b.set_ylim(0, 1.2)
    b.set_xlabel(t(L, "스크린 위 위치 (mm)", "position on screen (mm)"), **font(L))
    b.set_ylabel(t(L, "세기 (상대)", "intensity (relative)"), **font(L))
    b.set_title(t(L, "(b) 간격 0.5 mm, 1 m 거리, 550 nm", "(b) 0.5 mm spacing, 1 m away, 550 nm"), fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig3-two-point-sources.png")


# ------------------------------------------------------------------ 그림 5: 박막 쐐기와 뉴턴 링

def figure5(L):
    lam = 0.5893                                   # 나트륨 D 선
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(13, 4.2), gridspec_kw=dict(width_ratios=[1.2, 1, 1]))
    tt = np.linspace(0, 1.2, 1200)
    Rw2 = thin_layer_reflectance(1.0, 1.33, 1.0, tt, lam, multiple=False)
    Rwm = thin_layer_reflectance(1.0, 1.33, 1.0, tt, lam, multiple=True)
    a.plot(tt * 1000, Rw2 * 100, color=C_KEY, lw=2, label=t(L, "첫 두 반사만", "first two reflections"))
    a.plot(tt * 1000, Rwm * 100, color="k", lw=1, ls="--", label=t(L, "다중반사 모두", "all reflections"))
    a.axvline(lam / (4 * 1.33) * 1000, color=C_GREY, ls=":", lw=1)
    a.text(lam / (4 * 1.33) * 1000 + 8, 15.5, r"$\lambda/4n$", fontsize=10)
    a.set_xlabel(t(L, "막 두께 (nm)", "film thickness (nm)"), **font(L))
    a.set_ylabel(t(L, "반사율 (%)", "reflectance (%)"), **font(L))
    legend(a, L, loc="upper right")
    a.set_ylim(0, 18)
    a.set_title(t(L, "(a) 공기 중 물막(n = 1.33)의 반사율\n두께 0 에서 어둡다", "(a) Reflectance of a water film (n = 1.33) in air\ndark at zero thickness"),
                fontsize=10.5, **font(L))
    # (b) 쐐기: 두께가 x 에 비례
    xs = np.linspace(0, 10, 1000)                  # mm
    wedge_angle = 1e-4
    th = xs * 1000 * wedge_angle                   # µm
    Rx = thin_layer_reflectance(1.0, 1.33, 1.0, th, lam)
    b.imshow(np.tile(Rx, (120, 1)), extent=[0, 10, 0, 1.2], cmap="gray", aspect="auto", vmin=0, vmax=Rx.max())
    b.set_yticks([])
    b.set_xlabel(t(L, "막의 한쪽 끝에서 잰 거리 (mm)", "distance from the thin edge (mm)"), **font(L))
    b.set_title(t(L, f"(b) 기울기 {wedge_angle * 1e3:.1f} mrad 쐐기의 같은 두께 무늬\n" + r"간격 $\lambda/(2n\alpha)$" + f" = {lam / (2 * 1.33 * wedge_angle) / 1000:.2f} mm",
                  f"(b) Equal-thickness fringes of a {wedge_angle * 1e3:.1f} mrad wedge\n" + r"spacing $\lambda/(2n\alpha)$" + f" = {lam / (2 * 1.33 * wedge_angle) / 1000:.2f} mm"),
                fontsize=10.5, **font(L))
    # (c) 뉴턴 링: R = 1 m 볼록면 + 평판, 공기층, 유리 n = 1.52
    Rc = 1.0e6
    s = np.linspace(-3000, 3000, 701)
    XX, YY = np.meshgrid(s, s)
    gap = newton_gap(np.hypot(XX, YY), Rc)
    # 위 유리 → 공기층 → 아래 유리
    Rn = thin_layer_reflectance(1.52, 1.0, 1.52, gap, lam)
    c.imshow(Rn, extent=[-3, 3, -3, 3], cmap="gray", origin="lower")
    for m in [1, 2, 3, 4, 8, 12]:
        rm = np.sqrt(m * lam * Rc) / 1000
        c.add_patch(plt.Circle((0, 0), rm, fill=False, ec=C_KEY, lw=0.6, ls="--"))
    c.set_xlabel("x (mm)")
    c.set_title(t(L, "(c) 곡률반지름 1 m 렌즈의 뉴턴 링 (589 nm)\n" + r"주황 점선: 반지름 $\sqrt{m\lambda R}$", "(c) Newton's rings of an R = 1 m lens (589 nm)\n" + r"orange dashed: radius $\sqrt{m\lambda R}$"),
                fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig5-thin-film-newton.png")


# ------------------------------------------------------------------ 그림 8: 마이컬슨의 두 출력

def figure8(L):
    lam = 0.6328
    dz = np.linspace(0, 1.5 * lam, 600)            # 거울 이동 (µm), 경로차는 두 배
    delta = 2 * np.pi * 2 * dz / lam
    a_out, b_out = michelson_outputs(delta)
    fig, (a, b) = plt.subplots(1, 2, figsize=(11.5, 4.0), gridspec_kw=dict(width_ratios=[1.4, 1]))
    a.plot(dz * 1000, a_out, color=C_RED, lw=2, label=t(L, "검출기 쪽 출력", "output to the detector"))
    a.plot(dz * 1000, b_out, color=C_RAY, lw=2, label=t(L, "광원 쪽으로 돌아가는 출력", "output back to the source"))
    a.plot(dz * 1000, a_out + b_out, color="k", lw=1, ls="--", label=t(L, "합", "sum"))
    for k in range(1, 4):
        a.axvline(k * lam / 2 * 1000, color=C_GREY, lw=0.7, ls=":")
    a.set_xlabel(t(L, "거울을 민 거리 (nm), 633 nm 레이저", "mirror displacement (nm), 633 nm laser"), **font(L))
    a.set_ylabel(t(L, "세기 (입력 = 1)", "intensity (input = 1)"), **font(L))
    a.set_ylim(0, 1.15)
    legend(a, L, loc="upper center", ncol=3)
    a.set_title(t(L, "(a) 거울을 λ/2 밀 때마다 무늬 하나가 지나간다", "(a) One fringe passes for every λ/2 of mirror travel"),
                fontsize=10.5, **font(L))
    # (b) 원형 무늬: 두 거울의 상이 평행하게 거리 e 떨어져 있을 때 각 θ 의 경로차 2e cosθ
    e = 1000.0
    th = np.linspace(-0.06, 0.06, 601)
    TX, TY = np.meshgrid(th, th)
    ang = np.hypot(TX, TY)
    I = 0.5 * (1 + np.cos(2 * np.pi * 2 * e * np.cos(ang) / lam))
    b.imshow(I, extent=[th[0] * 1e3, th[-1] * 1e3, th[0] * 1e3, th[-1] * 1e3], cmap="gray", origin="lower")
    b.set_xlabel(t(L, "보는 각 (mrad)", "viewing angle (mrad)"), **font(L))
    b.set_title(t(L, "(b) 두 팔 길이가 1 mm 다를 때\n같은 기울기 무늬(원형)", "(b) Arms differing by 1 mm:\nequal-inclination (circular) fringes"),
                fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig8-michelson.png")


# ------------------------------------------------------------------ 그림 9: 수차의 간섭무늬

def figure9(L):
    n = 401
    s = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(s, s)
    M = X**2 + Y**2 <= 1
    r2 = X**2 + Y**2
    cases = [(t(L, "수차 없음", "no aberration"), 0 * X),
             (t(L, "초점 어긋남 2λ", "defocus 2λ"), 2 * r2),
             (t(L, "구면수차 3λ", "spherical 3λ"), 3 * r2**2),
             (t(L, "코마 3λ", "coma 3λ"), 3 * Y * r2),
             (t(L, "비점수차 3λ", "astigmatism 3λ"), 3 * Y**2)]
    fig, axes = plt.subplots(2, 5, figsize=(13, 5.6))
    for j, (name, W) in enumerate(cases):
        for row, tilt in enumerate([0.0, 4.0]):
            I = interferogram(W, tilt, X)
            ax = axes[row, j]
            ax.imshow(np.where(M, I, np.nan), cmap="gray", origin="lower", extent=[-1, 1, -1, 1], vmin=0, vmax=1)
            ax.set_xticks([]); ax.set_yticks([])
            if row == 0:
                ax.set_title(name, fontsize=10.5, **font(L))
    axes[0, 0].set_ylabel(t(L, "기준파와 나란히", "reference parallel"), fontsize=10, **font(L))
    axes[1, 0].set_ylabel(t(L, "기준파를 4λ 기울임", "reference tilted by 4λ"), fontsize=10, **font(L))
    fig.suptitle(t(L, "4편의 수차를 평면 기준파와 겹친 간섭무늬 (한 번 지나는 파면, 무늬 하나 = 파면 차 1λ)",
                   "Interferograms of Part 4's aberrations against a plane reference (single pass, one fringe = 1λ of wavefront)"),
                 fontsize=11, **font(L))
    fig.tight_layout()
    save(fig, L, "fig9-aberration-interferograms.png")


# ------------------------------------------------------------------ 그림 11: 여러 빔의 간섭과 파브리-페로

def figure11(L):
    phi = np.linspace(-np.pi, 3 * np.pi, 4001)
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.3))
    for M, col in [(2, C_GREY), (5, C_KEY), (20, C_RAY)]:
        I = multibeam_sum(phi, 1.0, M) / M**2
        a.plot(phi / np.pi, I, color=col, lw=1.8, label=t(L, f"같은 세기의 빔 {M}개", f"{M} equal beams"))
    a.set_xlabel(t(L, "이웃 빔 사이 위상차 (π 단위)", "phase step between beams (units of π)").replace("π", r"$\pi$"), **font(L))
    a.set_ylabel(t(L, "세기 (정점 = 1)", "intensity (peak = 1)"), **font(L))
    legend(a, L, loc="upper right")
    a.set_ylim(0, 1.4)
    a.set_title(t(L, "(a) 빔이 많을수록 봉우리가 좁아진다", "(a) More beams, narrower peaks"), fontsize=10.5, **font(L))
    for R, col in [(0.04, C_GRN), (0.5, C_KEY), (0.9, C_RAY), (0.98, C_PUR)]:
        T = etalon_transmission(phi, R)
        F = np.pi * np.sqrt(R) / (1 - R)
        b.plot(phi / np.pi, T, color=col, lw=1.8, label=t(L, f"거울 반사율 {R:g} (피네스 {F:.1f})", f"mirror reflectance {R:g} (finesse {F:.1f})"))
    b.set_xlabel(t(L, "왕복 위상 (π 단위) = 4πnd/λ", "round-trip phase (units of π) = 4πnd/λ").replace("π", r"$\pi$").replace("λ", r"$\lambda$"), **font(L))
    b.set_ylabel(t(L, "투과율", "transmittance"), **font(L))
    legend(b, L, loc="upper right", size=8)
    b.set_ylim(0, 1.45)
    b.set_title(t(L, "(b) 파브리-페로 간섭계의 투과율", "(b) Transmittance of a Fabry–Pérot interferometer"), fontsize=10.5, **font(L))
    fig.tight_layout()
    save(fig, L, "fig11-fabry-perot.png")


FIGS = {2: figure2, 3: figure3, 5: figure5, 8: figure8, 9: figure9, 11: figure11}


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

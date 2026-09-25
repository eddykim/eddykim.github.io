"""전자현미경 배경이론 2편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/electron-microscopy-electron-optics/      (한국어)
      ../../assets/img/posts/electron-microscopy-electron-optics/en/   (영문)

자기렌즈 궤적은 Glaser 종형 자기장에 대한 근축 광선 방정식을 RK4 로 적분해서 얻는다.
scipy 를 쓰지 않는 이유는 _code/requirements.txt 가 numpy 와 matplotlib 만 요구하기
때문이다. 두 언어의 그림은 같은 계산 결과를 쓰고 라벨만 다르다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "electron-microscopy-electron-optics",
)

# ── 물리 상수 ─────────────────────────────────────────────
QE = 1.602176634e-19
M0 = 9.1093837015e-31
C = 2.99792458e8
H = 6.62607015e-34

# ── 렌즈 조건 ─────────────────────────────────────────────
V_ACC = 200e3        # 가속전압 [V]
B0 = 0.30            # 종형 자기장의 최대값 [T]
A_BELL = 5.0e-3      # 종형 자기장의 반치 반폭 [m]


def relativistic_potential(v):
    """상대론 보정 전위 V* = V (1 + eV / 2 m0 c^2) [V]. 1편과 같은 보정항이다."""
    return v * (1.0 + QE * v / (2.0 * M0 * C**2))


def wavelength(v):
    """가속전압 v 에서의 전자 파장 [m]. 1편과 같은 식이다."""
    return H / np.sqrt(2.0 * M0 * QE * relativistic_potential(v))


def bz(z):
    """Glaser 종형 모델 B_z(z) = B0 / (1 + (z/a)^2) [T]."""
    return B0 / (1.0 + (z / A_BELL) ** 2)


def focusing_k(z, v=V_ACC):
    """근축 광선 방정식 r'' = -k(z) r 의 k(z) [1/m^2].

    k(z) = e B_z(z)^2 / (8 m0 V*)  — 회전 좌표계(Larmor frame)에서 성립한다.
    B_z 가 제곱으로 들어가므로 자기장의 방향을 뒤집어도 k 는 양수다.
    자기렌즈가 언제나 수렴 렌즈인 이유가 여기에 있다.
    """
    return QE * bz(z) ** 2 / (8.0 * M0 * relativistic_potential(v))


def larmor_rate(z, v=V_ACC):
    """상이 광축 둘레로 돌아가는 비율 dtheta/dz = sqrt(k(z)) [rad/m]."""
    return np.sqrt(focusing_k(z, v))


def trace_ray(r0, z_grid):
    """평행 입사 광선(r=r0, r'=0)을 RK4 로 적분한다.

    상태 y = [r, r', theta] 이고, r'' = -k(z) r, theta' = sqrt(k(z)) 이다.
    """
    def deriv(z, y):
        return np.array([y[1], -focusing_k(z) * y[0], larmor_rate(z)])

    y = np.array([r0, 0.0, 0.0])
    out = np.empty((len(z_grid), 3))
    out[0] = y
    for i in range(len(z_grid) - 1):
        z, h = z_grid[i], z_grid[i + 1] - z_grid[i]
        k1 = deriv(z, y)
        k2 = deriv(z + h / 2, y + h / 2 * k1)
        k3 = deriv(z + h / 2, y + h / 2 * k2)
        k4 = deriv(z + h, y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        out[i + 1] = y
    return out


def focal_length(r0, z_grid):
    """적분 결과에서 광선이 광축을 만나는 위치로 초점거리를 읽는다."""
    out = trace_ray(r0, z_grid)
    r, rp = out[:, 0], out[:, 1]
    # 렌즈를 빠져나온 뒤(자기장이 거의 0인 구간)는 직선이다.
    r_end, rp_end, z_end = r[-1], rp[-1], z_grid[-1]
    z_cross = z_end - r_end / rp_end
    f = -r0 / rp_end          # 주평면 기준 초점거리
    return z_cross, f


# ── 수차와 분해능 ─────────────────────────────────────────
def d_diffraction(alpha, lam):
    """회절에 의한 흐림 [m]. 1편에서 쓴 Airy 원판 기준."""
    return 0.61 * lam / alpha


def d_spherical(alpha, cs):
    """구면수차에 의한 흐림 [m]. 교재에 따라 계수가 1, 1/2, 1/4 로 갈린다."""
    return cs * alpha**3


def d_chromatic(alpha, cc, de_over_e):
    """색수차에 의한 흐림 [m]."""
    return cc * alpha * de_over_e


def d_total(alpha, lam, cs, cc, de_over_e):
    """세 항을 제곱합으로 합친 총 흐림 [m]."""
    return np.sqrt(d_diffraction(alpha, lam) ** 2
                   + d_spherical(alpha, cs) ** 2
                   + d_chromatic(alpha, cc, de_over_e) ** 2)


def optimum_alpha(lam, cs):
    """회절과 구면수차만 고려했을 때의 최적 조리개각 [rad].

    d^2 = (0.61 lam / alpha)^2 + (cs alpha^3)^2 를 alpha 로 미분해 0 으로 두면
    alpha_opt = (0.61 lam / (sqrt(3) cs))^(1/4) 가 나온다.
    """
    return (0.61 * lam / (np.sqrt(3.0) * cs)) ** 0.25


def scherzer_resolution(lam, cs):
    """Scherzer 점분해능 d = 0.66 (Cs lam^3)^(1/4) [m]. 교과서 표준 표기."""
    return 0.66 * (cs * lam**3) ** 0.25


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "fig1": "자기렌즈 안에서 전자가 지나가는 길",
        "f1a": "축상 자기장 $B_z(z)$", "f1ax": "광축 위치 z (mm)", "f1ay": "$B_z$ (T)",
        "f1b": "근축 광선의 궤적 (회전 좌표계)",
        "f1bx": "광축 위치 z (mm)", "f1by": "광축에서의 거리 r (mm)",
        "f1c": "광축에 수직인 면으로 본 같은 궤적",
        "f1cx": "x (mm)", "f1cy": "y (mm)",
        "f1focus": "초점 f = {:.1f} mm", "f1lens": "렌즈 자기장 영역",
        "f1rot": "상 회전 {:.0f}°", "f1in": "입사", "f1out": "초점",
        "fig2": "조리개각이 만드는 세 가지 흐림과 그 합",
        "f2x": "조리개 반각 α (mrad)", "f2y": "흐림의 크기 (nm)",
        "f2diff": "회절 $0.61λ/α$", "f2sph": "구면수차 $C_sα^3$",
        "f2chr": "색수차 $C_cα\\,ΔE/E$", "f2tot": "합 (제곱합)",
        "f2opt": "최적 α = {:.1f} mrad\n분해능 {:.2f} nm",
        "f2a": "수차보정 전 ($C_s$ = 0.5 mm)", "f2b": "수차보정 후 ($C_s$ = 1 $\mu$m)",
        "fig3": "구면수차 계수가 분해능을 어떻게 끌어내리는가",
        "f3x": "구면수차 계수 $C_s$ (mm)", "f3y": "Scherzer 점분해능 (nm)",
        "f3_200": "200 kV (λ = 2.51 pm)", "f3_300": "300 kV (λ = 1.97 pm)",
        "f3unc": "보정 전\n$C_s$ ≈ 0.5 mm", "f3cor": "보정 후\n$C_s$ ≈ 1 $\mu$m",
        "f3slope": "기울기 1/4\n$C_s$ 를 1000배 줄여야\n분해능이 5.6배 좋아진다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "The path an electron takes through a magnetic lens",
        "f1a": "Axial field $B_z(z)$", "f1ax": "Axial position z (mm)", "f1ay": "$B_z$ (T)",
        "f1b": "Paraxial ray trajectory (rotating frame)",
        "f1bx": "Axial position z (mm)", "f1by": "Distance from axis r (mm)",
        "f1c": "The same trajectory seen along the axis",
        "f1cx": "x (mm)", "f1cy": "y (mm)",
        "f1focus": "Focus f = {:.1f} mm", "f1lens": "Lens field region",
        "f1rot": "Image rotation {:.0f}°", "f1in": "Entry", "f1out": "Focus",
        "fig2": "Three sources of blur set by aperture angle, and their sum",
        "f2x": "Aperture semi-angle α (mrad)", "f2y": "Blur diameter (nm)",
        "f2diff": "Diffraction $0.61λ/α$", "f2sph": "Spherical $C_sα^3$",
        "f2chr": "Chromatic $C_cα\\,ΔE/E$", "f2tot": "Sum (in quadrature)",
        "f2opt": "Optimum α = {:.1f} mrad\nresolution {:.2f} nm",
        "f2a": "Before correction ($C_s$ = 0.5 mm)", "f2b": "After correction ($C_s$ = 1 $\mu$m)",
        "fig3": "How the spherical aberration coefficient pulls resolution down",
        "f3x": "Spherical aberration coefficient $C_s$ (mm)",
        "f3y": "Scherzer point resolution (nm)",
        "f3_200": "200 kV (λ = 2.51 pm)", "f3_300": "300 kV (λ = 1.97 pm)",
        "f3unc": "Uncorrected\n$C_s$ ≈ 0.5 mm", "f3cor": "Corrected\n$C_s$ ≈ 1 $\mu$m",
        "f3slope": "Slope 1/4\na 1000x smaller $C_s$\nbuys only 5.6x resolution",
    },
}

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREEN = "#27795b"
GRAY = "#7f8c8d"


def fig1_lens_trajectory(L):
    z = np.linspace(-30e-3, 60e-3, 6000)
    r0s = [0.3e-3, 0.6e-3, 0.9e-3]
    traces = [trace_ray(r0, z) for r0 in r0s]
    z_cross, f = focal_length(r0s[-1], z)

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    zmm = z * 1e3

    ax = axes[0]
    ax.plot(zmm, bz(z), color=BLUE, lw=2.0)
    ax.fill_between(zmm, 0, bz(z), color=BLUE, alpha=0.12)
    ax.set_xlabel(L["f1ax"], **L["font"])
    ax.set_ylabel(L["f1ay"], **L["font"])
    ax.set_title(L["f1a"], **L["font"])
    ax.grid(alpha=0.3)
    ax.set_xlim(zmm[0], zmm[-1])

    ax = axes[1]
    for r0, tr in zip(r0s, traces):
        ax.plot(zmm, tr[:, 0] * 1e3, color=BLUE, lw=1.6)
        ax.plot(zmm, -tr[:, 0] * 1e3, color=BLUE, lw=1.6)
    ax.axhline(0, color="#2c3e50", lw=0.8, ls="-")
    ax.axvspan(-2 * A_BELL * 1e3, 2 * A_BELL * 1e3, color=BLUE, alpha=0.10)
    ax.set_ylim(-2.6, 2.6)
    ax.text(0, 2.25, L["f1lens"], ha="center", fontsize=8.5, color=BLUE, **L["font"])
    ax.plot([z_cross * 1e3], [0], "o", color=ACCENT, ms=7, zorder=5)
    # 주석은 광선이 없는 축 아래쪽 빈 자리에 둔다.
    ax.annotate(L["f1focus"].format(f * 1e3), xy=(z_cross * 1e3, 0),
                xytext=(z_cross * 1e3 + 6, -1.9), color=ACCENT, fontsize=9,
                ha="left", arrowprops=dict(arrowstyle="->", color=ACCENT), **L["font"])
    ax.set_xlabel(L["f1bx"], **L["font"])
    ax.set_ylabel(L["f1by"], **L["font"])
    ax.set_title(L["f1b"], **L["font"])
    ax.grid(alpha=0.3)
    ax.set_xlim(zmm[0], zmm[-1])

    # 회전 좌표계를 실험실 좌표계로 되돌리면 궤적이 나선이 된다.
    # 초점을 지나면 r 의 부호가 뒤집혀 그림이 헷갈리므로 초점까지만 그린다.
    ax = axes[2]
    tr = traces[-1]
    m = zmm <= z_cross * 1e3
    r, theta = tr[m, 0], tr[m, 2]
    x, y = r * np.cos(theta) * 1e3, r * np.sin(theta) * 1e3
    pts = ax.scatter(x, y, c=zmm[m], cmap="viridis", s=4, zorder=3)
    ax.plot([x[0]], [y[0]], "o", color=BLUE, ms=7, zorder=5)
    ax.annotate(L["f1in"], xy=(x[0], y[0]), xytext=(-8, 8),
                textcoords="offset points", color=BLUE, fontsize=9,
                ha="right", **L["font"])
    ax.plot([0], [0], "o", color=ACCENT, ms=7, zorder=5)
    ax.annotate(L["f1out"], xy=(0, 0), xytext=(10, -12),
                textcoords="offset points", color=ACCENT, fontsize=9, **L["font"])
    ax.axhline(0, color=GRAY, lw=0.6)
    ax.axvline(0, color=GRAY, lw=0.6)
    ax.set_aspect("equal")
    lim = abs(r[0]) * 1e3 * 1.35
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(L["f1cx"], **L["font"])
    ax.set_ylabel(L["f1cy"], **L["font"])
    ax.set_title(L["f1c"], **L["font"])
    ax.grid(alpha=0.3)
    total_rot = np.degrees(theta[-1])
    ax.text(0.04, 0.05, L["f1rot"].format(total_rot), transform=ax.transAxes,
            fontsize=9.5, color=ACCENT, **L["font"])
    fig.colorbar(pts, ax=ax, label="z (mm)", fraction=0.046)

    fig.suptitle(L["fig1"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig, f, total_rot


def _tradeoff_panel(ax, L, lam, cs, cc, de_over_e, title):
    alpha = np.linspace(0.3e-3, 40e-3, 2000)
    amr = alpha * 1e3
    dd = d_diffraction(alpha, lam) * 1e9
    ds = d_spherical(alpha, cs) * 1e9
    dc = d_chromatic(alpha, cc, de_over_e) * 1e9
    dt = d_total(alpha, lam, cs, cc, de_over_e) * 1e9

    ax.plot(amr, dd, color=BLUE, lw=1.7, ls="--", label=L["f2diff"])
    ax.plot(amr, ds, color=ACCENT, lw=1.7, ls="--", label=L["f2sph"])
    ax.plot(amr, dc, color=GREEN, lw=1.7, ls=":", label=L["f2chr"])
    ax.plot(amr, dt, color="#2c3e50", lw=2.4, label=L["f2tot"])

    i = int(np.argmin(dt))
    ax.plot([amr[i]], [dt[i]], "o", color="#2c3e50", ms=7, zorder=6)
    ax.annotate(L["f2opt"].format(amr[i], dt[i]), xy=(amr[i], dt[i]),
                xytext=(amr[i] * 1.9, dt[i] * 2.6), fontsize=9, color="#2c3e50",
                arrowprops=dict(arrowstyle="->", color="#2c3e50"), **L["font"])

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(L["f2x"], **L["font"])
    ax.set_ylabel(L["f2y"], **L["font"])
    ax.set_title(title, **L["font"])
    ax.set_ylim(5e-3, 20)
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop=L["legend"] or None, fontsize=8.5, loc="upper left")
    return amr[i], dt[i]


def fig2_aperture_tradeoff(L):
    lam = wavelength(V_ACC)
    de_over_e = 0.7 / 200e3          # 쇼트키 전계방출총 기준 에너지폭
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    a1, d1 = _tradeoff_panel(axes[0], L, lam, 0.5e-3, 1.2e-3, de_over_e, L["f2a"])
    a2, d2 = _tradeoff_panel(axes[1], L, lam, 1.0e-6, 1.2e-3, de_over_e, L["f2b"])
    fig.suptitle(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig, (a1, d1, a2, d2)


def fig3_resolution_vs_cs(L):
    cs = np.logspace(-6, -2.3, 400)          # 1 µm ~ 5 mm
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for v, color, key in [(200e3, BLUE, "f3_200"), (300e3, ACCENT, "f3_300")]:
        lam = wavelength(v)
        ax.loglog(cs * 1e3, scherzer_resolution(lam, cs) * 1e9,
                  color=color, lw=2.1, label=L[key])

    lam200 = wavelength(200e3)
    # 주석은 두 곡선을 피해 점 기준 오프셋으로 놓는다.
    for cs_mark, key, offset, ha in [(0.5e-3, "f3unc", (-12, 16), "right"),
                                     (1.0e-6, "f3cor", (14, -26), "left")]:
        d = scherzer_resolution(lam200, cs_mark) * 1e9
        ax.plot([cs_mark * 1e3], [d], "o", color="#2c3e50", ms=7, zorder=5)
        ax.annotate(f"{L[key]}\n{d:.3f} nm", xy=(cs_mark * 1e3, d),
                    xytext=offset, textcoords="offset points", fontsize=9,
                    color="#2c3e50", ha=ha,
                    va="top" if offset[1] < 0 else "bottom", **L["font"])

    ax.set_ylim(bottom=0.014)
    ax.set_xlabel(L["f3x"], **L["font"])
    ax.set_ylabel(L["f3y"], **L["font"])
    ax.set_title(L["fig3"], fontsize=13, **L["font"])
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop=L["legend"] or None, loc="lower right")
    ax.text(0.035, 0.93, L["f3slope"], transform=ax.transAxes, fontsize=9,
            color=GRAY, va="top", **L["font"])
    fig.tight_layout()
    return fig


def main():
    stats = {}
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        fig, f, rot = fig1_lens_trajectory(L)
        fig.savefig(os.path.join(L["dir"], "fig1-magnetic-lens-trajectory.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
        stats["f"], stats["rot"] = f, rot

        fig, tr = fig2_aperture_tradeoff(L)
        fig.savefig(os.path.join(L["dir"], "fig2-aperture-tradeoff.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
        stats["tradeoff"] = tr

        fig = fig3_resolution_vs_cs(L)
        fig.savefig(os.path.join(L["dir"], "fig3-resolution-vs-cs.png"),
                    dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"[{lang}] 3 figures written to {L['dir']}")

    # 본문에 인용한 수치를 검산할 수 있게 출력한다.
    lam = wavelength(V_ACC)
    print()
    print(f"200 kV : lambda = {lam*1e12:.3f} pm, V* = {relativistic_potential(V_ACC)/1e3:.1f} kV")
    print(f"렌즈    : B0 = {B0} T, a = {A_BELL*1e3:.0f} mm "
          f"-> f = {stats['f']*1e3:.2f} mm, 상 회전 = {stats['rot']:.1f} deg")
    a1, d1, a2, d2 = stats["tradeoff"]
    print(f"Cs=0.5mm: alpha_opt = {a1:.2f} mrad, d_total = {d1:.3f} nm")
    print(f"Cs=1um  : alpha_opt = {a2:.2f} mrad, d_total = {d2:.3f} nm")
    for cs in (0.5e-3, 1e-6):
        print(f"Scherzer 점분해능 (Cs={cs*1e3:g} mm) = "
              f"{scherzer_resolution(lam, cs)*1e9:.3f} nm")
    print(f"해석식 alpha_opt (Cs=0.5mm) = {optimum_alpha(lam, 0.5e-3)*1e3:.2f} mrad")


if __name__ == "__main__":
    main()

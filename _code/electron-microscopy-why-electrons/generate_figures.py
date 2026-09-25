"""전자현미경 배경이론 1편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/electron-microscopy-why-electrons/      (한국어)
      ../../assets/img/posts/electron-microscopy-why-electrons/en/   (영문)

모든 물리량은 SI 단위로 계산하고 표시할 때만 pm/nm 로 바꾼다.
두 언어의 그림은 같은 계산 결과를 쓰고 표기만 다르다.
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 한글 텍스트에만 한글 폰트를 지정한다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "electron-microscopy-why-electrons",
)

# ── 물리 상수 (CODATA 2018) ────────────────────────────────
H = 6.62607015e-34      # 플랑크 상수 [J·s]
M0 = 9.1093837015e-31   # 전자 정지질량 [kg]
QE = 1.602176634e-19    # 기본 전하 [C]
C = 2.99792458e8        # 광속 [m/s]
KB = 1.380649e-23       # 볼츠만 상수 [J/K]


def wavelength(volts, relativistic=True):
    """가속전압 V 로 가속된 전자의 드브로이 파장 [m].

    비상대론: lambda = h / sqrt(2 m0 e V)
    상대론  : 분모에 (1 + eV / 2 m0 c^2) 보정항이 곱해진다.
    """
    v = np.asarray(volts, dtype=float)
    p2 = 2.0 * M0 * QE * v
    if relativistic:
        p2 = p2 * (1.0 + QE * v / (2.0 * M0 * C**2))
    return H / np.sqrt(p2)


def mean_free_path(pressure_pa, temperature=293.0, diameter=3.7e-10):
    """분자 지름 d 인 기체의 평균자유행로 [m]. 공기 기준 d ~ 0.37 nm."""
    p = np.asarray(pressure_pa, dtype=float)
    return KB * temperature / (np.sqrt(2.0) * np.pi * diameter**2 * p)


# ── 전자총 4종의 대표값 ────────────────────────────────────
# 제조사 사양과 교과서 값의 범위를 대표하는 자릿수 수준의 값이다.
GUNS = [
    # key,        밝기[A/cm^2·sr], 에너지폭[eV], 소스크기[m], 동작온도[K], 요구진공[Pa]
    ("w",         1e5,  2.3, 50e-6, 2800, 1e-3),
    ("lab6",      1e6,  1.5, 10e-6, 1900, 1e-4),
    ("schottky",  1e8,  0.7, 20e-9, 1800, 1e-7),
    ("cfe",       1e9,  0.3,  5e-9,  300, 1e-8),
]

# ── 분해능의 역사 ──────────────────────────────────────────
# key, 분해능[nm], 막대 색 구분용 그룹
WAVELENGTH_200KV = 0.0025   # nm, 그림2 의 기준선
MILESTONES = [
    ("optical",  200.0,   "other"),
    ("ruska",     50.0,   "tem"),
    ("tem1970",    0.3,   "tem"),
    ("hrtem",      0.2,   "tem"),
    ("corrected",  0.05,  "tem"),
    ("ptycho",     0.02,  "tem"),
]

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        # 그림1
        "fig1": "가속전압에 따른 전자의 드브로이 파장",
        "f1x": "가속전압 (kV)", "f1y": "파장 (pm)",
        "f1rel": "상대론 보정", "f1nr": "비상대론 근사",
        "f1err": "비상대론 근사의 오차",
        "f1erry": "파장 과대평가 (%)",
        "f1ann": "200 kV\n2.51 pm",
        # 그림2
        "fig2": "전자의 파장과 실제로 달성된 분해능",
        "f2y": "분해능 (nm, 로그 눈금)",
        "wavelength": "200 kV 전자의 드브로이 파장 (2.5 pm)",
        "optical": "광학현미경\n(가시광 한계)",
        "ruska": "최초의 TEM\n(1933)",
        "tem1970": "1970년대 TEM",
        "hrtem": "수차보정 이전\nHRTEM",
        "corrected": "수차보정 STEM",
        "ptycho": "전자 타이코그래피",
        "f2gap": "약 80배",
        # 그림3
        "fig3": "전자총 4종의 밝기·에너지폭·소스 크기",
        "f3a": "밝기", "f3ay": "밝기 (A/cm²·sr)",
        "f3b": "에너지폭", "f3by": "에너지폭 ΔE (eV)",
        "f3c": "소스 크기", "f3cy": "소스 지름 (m)",
        "w": "텅스텐\n열전자", "lab6": "LaB$_6$\n열전자",
        "schottky": "쇼트키\n전계방출", "cfe": "냉전계방출",
        # 그림4
        "fig4": "진공도에 따른 전자의 평균자유행로",
        "f4x": "압력 (Pa)", "f4y": "평균자유행로 (m)",
        "f4col": "경통 길이 ~1 m",
        "f4atm": "대기압", "f4hv": "고진공\n(열전자총)", "f4uhv": "초고진공\n(전계방출총)",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "De Broglie wavelength of an accelerated electron",
        "f1x": "Accelerating voltage (kV)", "f1y": "Wavelength (pm)",
        "f1rel": "Relativistic", "f1nr": "Non-relativistic",
        "f1err": "Error of the non-relativistic approximation",
        "f1erry": "Wavelength overestimate (%)",
        "f1ann": "200 kV\n2.51 pm",
        "fig2": "Electron wavelength vs. resolution actually achieved",
        "f2y": "Resolution (nm, log scale)",
        "wavelength": "De Broglie wavelength, 200 kV electron (2.5 pm)",
        "optical": "Optical microscope\n(visible light)",
        "ruska": "First TEM\n(1933)",
        "tem1970": "TEM, 1970s",
        "hrtem": "HRTEM,\nuncorrected",
        "corrected": "Aberration-\ncorrected STEM",
        "ptycho": "Electron\nptychography",
        "f2gap": "~80x",
        "fig3": "Brightness, energy spread and source size of four electron guns",
        "f3a": "Brightness", "f3ay": "Brightness (A/cm²·sr)",
        "f3b": "Energy spread", "f3by": "Energy spread ΔE (eV)",
        "f3c": "Source size", "f3cy": "Source diameter (m)",
        "w": "Tungsten\nthermionic", "lab6": "LaB$_6$\nthermionic",
        "schottky": "Schottky\nfield emission", "cfe": "Cold field\nemission",
        "fig4": "Electron mean free path vs. vacuum level",
        "f4x": "Pressure (Pa)", "f4y": "Mean free path (m)",
        "f4col": "Column length ~1 m",
        "f4atm": "Atmosphere", "f4hv": "High vacuum\n(thermionic gun)",
        "f4uhv": "Ultra-high vacuum\n(field-emission gun)",
    },
}

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GRAY = "#7f8c8d"


def fig1_wavelength(L):
    kv = np.linspace(1.0, 300.0, 600)
    volts = kv * 1e3
    lam_rel = wavelength(volts, True) * 1e12      # pm
    lam_nr = wavelength(volts, False) * 1e12
    overestimate = (lam_nr / lam_rel - 1.0) * 100.0

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax.plot(kv, lam_rel, color=BLUE, lw=2.0, label=L["f1rel"])
    ax.plot(kv, lam_nr, color=ACCENT, lw=1.6, ls="--", label=L["f1nr"])
    ax.set_yscale("log")
    ax.set_xlabel(L["f1x"], **L["font"])
    ax.set_ylabel(L["f1y"], **L["font"])
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop=L["legend"] or None)

    lam200 = wavelength(200e3) * 1e12
    ax.plot([200], [lam200], "o", color=BLUE, ms=7, zorder=5)
    ax.annotate(L["f1ann"], xy=(200, lam200), xytext=(150, lam200 * 3.2),
                arrowprops=dict(arrowstyle="->", color=GRAY), color=BLUE,
                ha="center", **L["font"])

    ax2.plot(kv, overestimate, color=ACCENT, lw=2.0)
    ax2.set_xlabel(L["f1x"], **L["font"])
    ax2.set_ylabel(L["f1erry"], **L["font"])
    ax2.set_title(L["f1err"], **L["font"])
    ax2.grid(alpha=0.3)
    ax2.axvline(200, color=GRAY, ls=":", lw=1.2)
    err200 = (wavelength(200e3, False) / wavelength(200e3) - 1) * 100
    ax2.annotate(f"200 kV: {err200:.1f} %", xy=(200, err200),
                 xytext=(60, err200 * 0.72), color=ACCENT, **L["font"],
                 arrowprops=dict(arrowstyle="->", color=GRAY))

    fig.suptitle(L["fig1"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def fig2_resolution_gap(L):
    keys = [k for k, _, _ in MILESTONES]
    vals = np.array([v for _, v, _ in MILESTONES])
    groups = [g for _, _, g in MILESTONES]
    colors = {"other": GRAY, "tem": ACCENT}

    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    x = np.arange(len(keys))
    ax.bar(x, vals, color=[colors[g] for g in groups], width=0.55, zorder=3)
    ax.set_yscale("log")
    ax.set_ylim(1e-3, 2e3)
    ax.set_ylabel(L["f2y"], **L["font"])
    ax.set_xticks(x)
    ax.set_xticklabels([L[k] for k in keys], fontsize=9, **L["font"])
    ax.grid(axis="y", alpha=0.3, which="major", zorder=0)

    for xi, v in zip(x, vals):
        ax.text(xi, v * 1.35, f"{v:g} nm", ha="center", fontsize=9, color="#2c3e50")

    # 전자의 파장은 막대가 아니라 기준선으로 둔다. 어떤 막대도 여기까지 내려오지 못한다.
    ax.axhline(WAVELENGTH_200KV, color=BLUE, ls="--", lw=1.8, zorder=4,
               label=L["wavelength"])
    ax.legend(loc="upper right", prop=L["legend"] or None, framealpha=0.95)

    # 오른쪽에 빈 칸을 만들어, 막대와 겹치지 않는 자리에 간극을 표시한다.
    xa = len(keys) + 0.35
    ax.set_xlim(-0.6, len(keys) + 1.5)
    i_hr = keys.index("hrtem")
    ax.plot([i_hr + 0.28, xa], [vals[i_hr], vals[i_hr]],
            color="#2c3e50", ls=":", lw=1.0, zorder=4)
    ax.annotate("", xy=(xa, WAVELENGTH_200KV), xytext=(xa, vals[i_hr]),
                arrowprops=dict(arrowstyle="<->", color="#2c3e50", lw=1.5))
    ax.text(xa + 0.14, np.sqrt(WAVELENGTH_200KV * vals[i_hr]), L["f2gap"],
            fontsize=10, color="#2c3e50", va="center", ha="left", **L["font"])

    ax.set_title(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout()
    return fig


def fig3_guns(L):
    keys = [g[0] for g in GUNS]
    brightness = np.array([g[1] for g in GUNS])
    spread = np.array([g[2] for g in GUNS])
    size = np.array([g[3] for g in GUNS])

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.0))
    x = np.arange(len(keys))
    bar_colors = [GRAY, GRAY, BLUE, ACCENT]
    panels = [
        (brightness, L["f3a"], L["f3ay"], True, "{:.0e}"),
        (spread, L["f3b"], L["f3by"], False, "{:.1f} eV"),
        (size, L["f3c"], L["f3cy"], True, None),
    ]
    for ax, (vals, title, ylab, logy, fmt) in zip(axes, panels):
        ax.bar(x, vals, color=bar_colors, width=0.6, zorder=3)
        if logy:
            ax.set_yscale("log")
        ax.set_title(title, **L["font"])
        ax.set_ylabel(ylab, **L["font"])
        ax.set_xticks(x)
        ax.set_xticklabels([L[k] for k in keys], fontsize=8, **L["font"])
        ax.grid(axis="y", alpha=0.3, which="both", zorder=0)
        for xi, v in zip(x, vals):
            if fmt is None:
                txt = f"{v * 1e6:.0f} µm" if v > 1e-6 else f"{v * 1e9:.0f} nm"
            else:
                txt = fmt.format(v)
            off = v * 1.5 if logy else v + 0.12
            ax.text(xi, off, txt, ha="center", fontsize=8.5, color="#2c3e50")
        ax.set_ylim(top=(vals.max() * (8 if logy else 1.35)))

    fig.suptitle(L["fig3"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def fig4_mean_free_path(L):
    p = np.logspace(5, -8, 400)
    mfp = mean_free_path(p)

    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    ax.loglog(p, mfp, color=BLUE, lw=2.2, zorder=3)
    ax.set_xlabel(L["f4x"], **L["font"])
    ax.set_ylabel(L["f4y"], **L["font"])
    ax.invert_xaxis()
    ax.grid(alpha=0.3, which="both")

    ax.axhline(1.0, color=ACCENT, ls="--", lw=1.4)
    ax.text(3e4, 1.5, L["f4col"], color=ACCENT, fontsize=9.5, **L["font"])

    # 점마다 라벨이 축 밖으로 나가지 않도록 위치를 따로 준다.
    markers = [
        (1e5,  "f4atm", (12, 16),   "left"),
        (1e-4, "f4hv",  (-14, 12),  "right"),
        (1e-8, "f4uhv", (-16, 14),  "right"),
    ]
    for pressure, key, offset, ha in markers:
        m = mean_free_path(pressure)
        ax.plot([pressure], [m], "o", color="#2c3e50", ms=6, zorder=5)
        ax.annotate(L[key], xy=(pressure, m), xytext=offset,
                    textcoords="offset points", ha=ha, fontsize=9,
                    color="#2c3e50", **L["font"])
    ax.set_ylim(top=mean_free_path(1e-8) * 60)

    ax.set_title(L["fig4"], fontsize=13, **L["font"])
    fig.tight_layout()
    return fig


def main():
    builders = [
        ("fig1-wavelength-vs-voltage", fig1_wavelength),
        ("fig2-resolution-gap", fig2_resolution_gap),
        ("fig3-electron-guns", fig3_guns),
        ("fig4-mean-free-path", fig4_mean_free_path),
    ]
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for name, builder in builders:
            fig = builder(L)
            path = os.path.join(L["dir"], name + ".png")
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"[{lang}] {path}")

    # 본문에 인용한 수치를 콘솔에 남겨 검산할 수 있게 한다.
    print()
    for kv in (1, 30, 100, 200, 300):
        lam = wavelength(kv * 1e3) * 1e12
        lam_nr = wavelength(kv * 1e3, False) * 1e12
        print(f"{kv:3d} kV : lambda = {lam:6.3f} pm "
              f"(비상대론 {lam_nr:6.3f} pm, +{(lam_nr / lam - 1) * 100:4.1f} %)")
    for pressure in (1e5, 1e-3, 1e-4, 1e-7, 1e-8):
        print(f"{pressure:8.0e} Pa : mfp = {mean_free_path(pressure):.3e} m")


if __name__ == "__main__":
    main()

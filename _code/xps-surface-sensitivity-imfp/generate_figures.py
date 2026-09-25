"""XPS 기초 2편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/xps-surface-sensitivity-imfp/      (한국어)
      ../../assets/img/posts/xps-surface-sensitivity-imfp/en/   (영문)

IMFP는 TPP-2M 예측식으로 계산한다. 계수는 C. J. Powell, J. Vac. Sci. Technol. A 38,
023209 (2020) 식 (1)–(2)를 따랐고, 유효 범위는 50 eV 이상이다. 50 eV 아래는
Seah–Dench(1979) 유니버설 커브로 모양만 보인다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

# 한글 텍스트에만 한글 폰트를 지정한다. 전역 폰트를 바꾸면 눈금의 마이너스
# 기호가 한글 폰트에 없어 깨진다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "xps-surface-sensitivity-imfp",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "ke": "운동에너지 (eV)", "imfp": r"IMFP $\lambda$ (nm)",
        "sd": "Seah–Dench 유니버설 커브 (원소, a = 0.25 nm)",
        "alka": "Al K$\\alpha$ XPS 영역",
        "fig1": "운동에너지에 따른 비탄성 평균자유행로 (TPP-2M)",
        "depth": r"깊이 $z$ (nm)", "weight": "상대 기여 (단위 깊이당)",
        "cum": "누적 기여 (%)", "d": r"표면에서 깊이 $d$까지 (nm)",
        "fig2a": r"(a) 깊이별 기여 $e^{-z/\lambda}$", "fig2b": r"(b) 누적 기여 $1-e^{-d/\lambda}$",
        "fig2": r"깊이별 신호 기여 (Si 2p, $\lambda$ = {:.2f} nm)",
        "bin": "1 nm 구간별 몫",
        "surface": "표면", "normal": "법선", "analyzer": "분석기 방향",
        "sample": "시료", "path": r"탈출 경로 $z/\cos\theta$",
        "fig3a": "(a) 검출 각도와 탈출 경로", "fig3b": r"(b) 각도별 누적 기여 ($\lambda$ = {:.2f} nm)",
        "info95": r"95% 깊이 $3\lambda\cos\theta$",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "ke": "Kinetic energy (eV)", "imfp": r"IMFP $\lambda$ (nm)",
        "sd": "Seah–Dench universal curve (elements, a = 0.25 nm)",
        "alka": "Al K$\\alpha$ XPS range",
        "fig1": "Inelastic mean free path versus kinetic energy (TPP-2M)",
        "depth": r"Depth $z$ (nm)", "weight": "Relative contribution (per unit depth)",
        "cum": "Cumulative contribution (%)", "d": r"Depth $d$ from the surface (nm)",
        "fig2a": r"(a) Contribution by depth $e^{-z/\lambda}$", "fig2b": r"(b) Cumulative $1-e^{-d/\lambda}$",
        "fig2": r"Signal contribution by depth (Si 2p, $\lambda$ = {:.2f} nm)",
        "bin": "share per 1 nm slice",
        "surface": "surface", "normal": "normal", "analyzer": "to analyzer",
        "sample": "sample", "path": r"escape path $z/\cos\theta$",
        "fig3a": "(a) Emission angle and escape path", "fig3b": r"(b) Cumulative by angle ($\lambda$ = {:.2f} nm)",
        "info95": r"95% depth $3\lambda\cos\theta$",
    },
}

HV_AL = 1486.6  # Al Kα 광자 에너지 (eV)


def tpp2m(E, Nv, rho, M, Eg=0.0):
    """TPP-2M IMFP (nm). E: 운동에너지(eV), Nv: 원자·분자당 가전자 수,
    rho: 밀도(g/cm^3), M: 원자량·분자량, Eg: 밴드갭(eV)."""
    Ep = 28.816 * np.sqrt(Nv * rho / M)            # 자유전자 플라즈몬 에너지
    U = (Ep / 28.816) ** 2
    beta = -1.0 + 9.44 / np.sqrt(Ep**2 + Eg**2) + 0.69 * rho**0.1
    gamma = 0.191 * rho**-0.5
    C, D = 19.7 - 9.1 * U, 534 - 208 * U
    return E / (Ep**2 * (beta * np.log(gamma * E) - C / E + D / E**2))


def seah_dench(E, a=0.25):
    """Seah–Dench(1979) 원소용 유니버설 커브 (nm). a: 단원자층 두께(nm)."""
    return a * (538 / E**2 + 0.41 * np.sqrt(a * E))


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------
MATERIALS = {          # (Nv, rho, M, Eg, 색)
    "SiO$_2$": (16, 2.19, 60.08, 9.0, "tab:orange"),
    "Si": (4, 2.33, 28.086, 1.1, "tab:blue"),
    "Au": (11, 19.3, 196.97, 0.0, "tab:red"),
}
E_tpp = np.logspace(np.log10(50), np.log10(5000), 400)
E_sd = np.logspace(np.log10(5), np.log10(5000), 400)
imfp_curves = {k: tpp2m(E_tpp, *v[:4]) for k, v in MATERIALS.items()}
sd_curve = seah_dench(E_sd)
E_sd_min = E_sd[np.argmin(sd_curve)]

# Al Kα 여기 대표 피크: (라벨, 결합에너지, 물질 키)
PEAKS = [("Si 2p", 99.4, "Si"), ("O 1s", 533.0, "SiO$_2$"), ("Au 4f", 84.0, "Au")]
peak_pts = []
for label, be, mat in PEAKS:
    ke = HV_AL - be
    peak_pts.append((label, ke, float(tpp2m(ke, *MATERIALS[mat][:4])), MATERIALS[mat][4]))

# 그림2·3에 쓰는 λ: Si 속 Si 2p 광전자
LAM = float(tpp2m(HV_AL - 99.4, *MATERIALS["Si"][:4]))
z = np.linspace(0, 16, 800)
weight = np.exp(-z / LAM)
slices = np.arange(0, 12)                                   # 0–1, 1–2, ... nm
slice_share = np.exp(-slices / LAM) - np.exp(-(slices + 1) / LAM)
cum = 1 - np.exp(-z / LAM)
ANGLES = [0, 45, 60, 75]
cum_by_angle = {t: 1 - np.exp(-z / (LAM * np.cos(np.radians(t)))) for t in ANGLES}


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L_):
    """주어진 라벨 묶음으로 그림 3개를 그린다."""
    out, F, LF = L_["dir"], L_["font"], L_["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. IMFP vs 운동에너지 (로그-로그)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.axvspan(HV_AL - 1350, HV_AL, color="0.92", zorder=0)
    ax.text(np.sqrt((HV_AL - 1350) * HV_AL), 0.22, L_["alka"], ha="center",
            fontsize=9, color="0.35", **F)
    ax.plot(E_sd, sd_curve, color="0.4", lw=1.2, ls="--", label=L_["sd"])
    for name, (*_, col) in MATERIALS.items():
        ax.plot(E_tpp, imfp_curves[name], color=col, lw=2, label=f"{name} (TPP-2M)")
    offsets = {"Si 2p": ((10, -2), "left"), "O 1s": ((-10, 6), "right"),
               "Au 4f": ((10, -8), "left")}
    for label, ke, lam, col in peak_pts:
        ax.plot(ke, lam, "o", ms=7, mfc="white", mec=col, mew=2, zorder=5)
        off, ha = offsets[label]
        ax.annotate(f"{label}\n{lam:.1f} nm", xy=(ke, lam), xytext=off,
                    textcoords="offset points", ha=ha, va="center", fontsize=8.5, color=col)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(5, 5000)
    ax.set_ylim(0.2, 20)
    ax.set_xlabel(L_["ke"], fontsize=11, **F)
    ax.set_ylabel(L_["imfp"], fontsize=11, **F)
    ax.grid(alpha=0.3, which="both")
    ax.legend(prop={**LF, "size": 8.5}, loc="upper left")
    ax.set_title(L_["fig1"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-imfp-universal-curve.png"), dpi=150)
    plt.close(fig)

    # 그림 2. 깊이별 기여와 누적 기여
    fig2, (a1, a2) = plt.subplots(1, 2, figsize=(10.0, 4.2))
    a1.bar(slices + 0.5, slice_share / slice_share.max(), width=0.9, color="tab:blue",
           alpha=0.25, label=L_["bin"])
    for i in range(4):
        a1.text(i + 0.5, slice_share[i] / slice_share.max() + 0.03,
                f"{100 * slice_share[i]:.0f}%", ha="center", fontsize=8.5, color="tab:blue")
    a1.plot(z, weight, color="tab:blue", lw=2)
    a1.set_xlim(0, 12)
    a1.set_ylim(0, 1.15)
    a1.set_xlabel(L_["depth"], fontsize=11, **F)
    a1.set_ylabel(L_["weight"], fontsize=11, **F)
    a1.set_title(L_["fig2a"], fontsize=11, **F)
    a1.legend(prop={**LF, "size": 9}, loc="upper right")
    a1.grid(alpha=0.3)

    a2.plot(z, 100 * cum, color="tab:blue", lw=2)
    for n, col in zip((1, 2, 3), ("tab:green", "tab:orange", "tab:red")):
        d, c = n * LAM, 100 * (1 - np.exp(-n))
        a2.plot([d, d], [0, c], color=col, lw=1, ls="--")
        a2.plot([0, d], [c, c], color=col, lw=1, ls="--")
        a2.plot(d, c, "o", color=col)
        a2.annotate(f"{n}λ = {d:.1f} nm\n{c:.1f}%", xy=(d, c), xytext=(8, -26),
                    textcoords="offset points", fontsize=9, color=col)
    a2.set_xlim(0, 16)
    a2.set_ylim(0, 105)
    a2.set_xlabel(L_["d"], fontsize=11, **F)
    a2.set_ylabel(L_["cum"], fontsize=11, **F)
    a2.set_title(L_["fig2b"], fontsize=11, **F)
    a2.grid(alpha=0.3)
    fig2.suptitle(L_["fig2"].format(LAM), fontsize=12, **F)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "fig2-depth-contribution.png"), dpi=150)
    plt.close(fig2)

    # 그림 3. 검출 각도와 표면민감성
    fig3, (b1, b2) = plt.subplots(1, 2, figsize=(10.4, 4.4),
                                  gridspec_kw={"width_ratios": [1, 1.25]})
    b1.fill_between([-1.4, 2.6], -1.3, 0, color="0.88")
    b1.axhline(0, color="black", lw=1.5)
    b1.text(-1.3, -1.2, L_["sample"], fontsize=10, **F)
    b1.text(2.55, 0.06, L_["surface"], fontsize=9, ha="right", **F)
    b1.plot([0, 0], [-1.3, 1.6], color="gray", lw=1, ls=(0, (5, 4)))
    b1.text(0.06, 1.45, L_["normal"], fontsize=9, color="gray", **F)
    depth = 0.8                                      # 방출 지점의 깊이 z (도식 단위)
    src = np.array([0.0, -depth])
    b1.plot(*src, "o", color="black", zorder=5)
    b1.text(-0.12, -depth - 0.05, "$z$", fontsize=11, ha="right", va="center")
    for t, col in ((0, "tab:blue"), (60, "tab:red")):
        u = np.array([np.sin(np.radians(t)), np.cos(np.radians(t))])
        exit_pt = src + u * (depth / np.cos(np.radians(t)))
        b1.plot(*zip(src, exit_pt), color=col, lw=2.5)
        b1.annotate("", xy=exit_pt + 0.7 * u, xytext=exit_pt,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=1.8))
        b1.text(*(exit_pt + 0.8 * u), rf"$\theta$ = {t}°", color=col, fontsize=10,
                ha="center", va="bottom")
    ex60 = src + np.array([np.sin(np.radians(60)), np.cos(np.radians(60))]) * depth / 0.5
    b1.plot([ex60[0]] * 2, [0, 0.9], color="gray", lw=0.8, ls=":")
    arc = np.linspace(np.pi / 2 - np.radians(60), np.pi / 2, 30)
    b1.plot(ex60[0] + 0.4 * np.cos(arc), 0.4 * np.sin(arc), color="tab:red", lw=1)
    b1.text(ex60[0] + 0.12, 0.45, r"$\theta$", color="tab:red", fontsize=11)
    b1.text(0.55, -0.72, L_["path"], color="tab:red", fontsize=9, **F)
    b1.set_xlim(-1.4, 2.6)
    b1.set_ylim(-1.3, 1.7)
    b1.set_aspect("equal")
    b1.axis("off")
    b1.set_title(L_["fig3a"], fontsize=11, **F)

    colors = ["tab:blue", "tab:green", "tab:orange", "tab:red"]
    for t, col in zip(ANGLES, colors):
        d95 = 3 * LAM * np.cos(np.radians(t))
        b2.plot(z, 100 * cum_by_angle[t], color=col, lw=2,
                label=rf"$\theta$ = {t}°  (3$\lambda\cos\theta$ = {d95:.1f} nm)")
        b2.plot(d95, 95, "o", color=col)
    b2.axhline(95, color="gray", lw=1, ls=":")
    b2.text(15.8, 91, L_["info95"], ha="right", fontsize=9, color="gray", **F)
    b2.set_xlim(0, 16)
    b2.set_ylim(0, 105)
    b2.set_xlabel(L_["d"], fontsize=11, **F)
    b2.set_ylabel(L_["cum"], fontsize=11, **F)
    b2.set_title(L_["fig3b"].format(LAM), fontsize=11, **F)
    b2.legend(prop={**LF, "size": 9}, loc="lower right", bbox_to_anchor=(1.0, 0.08))
    b2.grid(alpha=0.3)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, "fig3-takeoff-angle.png"), dpi=150)
    plt.close(fig3)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

print()
print(f"Seah–Dench 최소: E = {E_sd_min:.0f} eV, λ = {sd_curve.min():.2f} nm")
for label, ke, lam, _ in peak_pts:
    print(f"{label}: KE {ke:.1f} eV, λ = {lam:.2f} nm, 3λ = {3 * lam:.1f} nm, "
          f"최표면 1 nm 몫 {100 * (1 - np.exp(-1 / lam)):.0f}%")
print(f"그림2 λ(Si 2p in Si) = {LAM:.2f} nm, 1 nm 구간 몫:",
      " ".join(f"{100 * s:.1f}%" for s in slice_share[:4]))
for t in ANGLES:
    print(f"θ = {t:2d}°: 3λcosθ = {3 * LAM * np.cos(np.radians(t)):.2f} nm, "
          f"최표면 1 nm 몫 {100 * cum_by_angle[t][np.searchsorted(z, 1.0)]:.0f}%")

# ---------------------------------------------------------------------------
# 본문 수치 (그림에는 없음)
# ---------------------------------------------------------------------------
# X선 감쇠길이: CXRO(Henke) atten2, 1486.6 eV, 수직 입사. 단위 µm.
XRAY_ATT_UM = {"Si": 7.89, "SiO$_2$": 4.11, "Au": 0.230}
for label, ke, lam, _ in peak_pts:
    mat = next(m for lb, _, m in PEAKS if lb == label)
    print(f"{label}: X선 감쇠길이/IMFP = {XRAY_ATT_UM[mat] * 1000 / lam:.0f}")

# 흡착 탄소층: 폴리에틸렌 (CH2)n 으로 근사 (Nv=6, rho=0.92, M=14.027, Eg≈8 eV 가정)
lam_ch2 = float(tpp2m(HV_AL - 99.4, 6, 0.92, 14.027, 8.0))
print(f"Si 2p 광전자의 탄화수소층 속 λ = {lam_ch2:.2f} nm; 남는 하부 신호:",
      " ".join(f"{t} nm→{100 * np.exp(-t / lam_ch2):.0f}%" for t in (1.0, 1.5, 2.0)))

"""엘립소메트리 배경이론 2편 그림 4개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometry-polarization-mueller-matrix/      (한국어)
      ../../assets/img/posts/ellipsometry-polarization-mueller-matrix/en/   (영문)

무작위 앙상블(그림2)은 같은 seed 를 쓰고 계산도 한 번만 하므로, 두 언어의
그림은 데이터가 완전히 동일하고 표기만 다르다.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

# 한글 텍스트에만 한글 폰트를 지정한다. macOS 전용 설정이다.
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
EFONT: dict = {}
ELFONT: dict = {}

BASE_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometry-polarization-mueller-matrix",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "linear": "선형편광", "circular": "원형편광", "elliptical": "타원편광",
        "fig1": "전기장 성분의 진폭비·위상차에 따른 편광 상태",
        "full": "완전편광", "partial": "부분편광", "unpol": "무편광",
        "dop": "{}\n(측정된 P={:.2f})",
        "fig2": "존스 벡터로는 왼쪽 한 장만 표현 가능 — 스토크스 벡터는 셋 다 표현",
        "source": "광원", "detector": "검출기",
        "polarizer": "편광자\n$M_P(\\Omega_P)$", "retarder": "위상지연자\n$M_C(\\Omega_C,\\phi)$",
        "sample": "시편\n$M_{sample}$", "retarder2": "위상지연자\n$M_C'$",
        "analyzer": "편광자\n$M_A$",
        "chain": "PSG → 시편 → PSA",
        "psg": "편광 생성단 (PSG)", "psa": "편광 분석단 (PSA)",
        "fig3": "PSG-시편-PSA 구조와 뮬러 행렬의 순차 곱 "
                "$S_{out}=M_A M_C' M_{sample} M_C M_P S_{in}$",
        "phi_label": r"위상 지연량 $\phi$ (rad)",
        "stokes_norm": "정규화된 스토크스 성분",
        "fig4": "편광자(0°)+위상지연자(45°) 통과 후 출력 편광: 선형 → 타원 → 원형 → 타원 → 선형",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "linear": "linear", "circular": "circular", "elliptical": "elliptical",
        "fig1": "Polarization state by amplitude ratio and phase difference of the field components",
        "full": "fully polarized", "partial": "partially polarized", "unpol": "unpolarized",
        "dop": "{}\n(measured P={:.2f})",
        "fig2": "A Jones vector describes only the left panel — a Stokes vector describes all three",
        "source": "source", "detector": "detector",
        "polarizer": "polarizer\n$M_P(\\Omega_P)$", "retarder": "retarder\n$M_C(\\Omega_C,\\phi)$",
        "sample": "sample\n$M_{sample}$", "retarder2": "retarder\n$M_C'$",
        "analyzer": "analyzer\n$M_A$",
        "chain": "PSG → sample → PSA",
        "psg": "Polarization State Generator", "psa": "Polarization State Analyzer",
        "fig3": "PSG-sample-PSA chain as a product of Mueller matrices "
                "$S_{out}=M_A M_C' M_{sample} M_C M_P S_{in}$",
        "phi_label": r"Retardance $\phi$ (rad)",
        "stokes_norm": "Normalized Stokes components",
        "fig4": "Output polarization after a polarizer (0°) and a retarder (45°): "
                "linear → elliptical → circular → elliptical → linear",
    },
}


# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소
# ---------------------------------------------------------------------------

def mueller_rotation(omega):
    """회전 행렬 M_R(omega), omega: rad."""
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([
        [1, 0, 0, 0],
        [0, c, s, 0],
        [0, -s, c, 0],
        [0, 0, 0, 1],
    ])


def mueller_polarizer(omega_p):
    """이상적 편광자 M_P(omega_p)."""
    base = 0.5 * np.array([
        [1, 1, 0, 0],
        [1, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ])
    return mueller_rotation(-omega_p) @ base @ mueller_rotation(omega_p)


def mueller_retarder(omega_c, phi):
    """위상 지연자 M_C(omega_c, phi)."""
    c, s = np.cos(phi), np.sin(phi)
    base = np.array([
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, c, s],
        [0, 0, -s, c],
    ])
    return mueller_rotation(-omega_c) @ base @ mueller_rotation(omega_c)


def stokes_to_ellipse(S, n_pts=200):
    """정규화된 스토크스 벡터로부터 편광 타원 (Ex(t), Ey(t)) 좌표를 복원."""
    S0, S1, S2, S3 = S
    s1, s2, s3 = S1 / S0, S2 / S0, S3 / S0
    psi = 0.5 * np.arctan2(s2, s1)
    s3c = np.clip(s3, -1, 1)
    chi = 0.5 * np.arcsin(s3c)
    a, b = np.cos(chi), np.sin(chi)
    t_ = np.linspace(0, 2 * np.pi, n_pts)
    X, Y = a * np.cos(t_), b * np.sin(t_)
    c, s = np.cos(psi), np.sin(psi)
    Ex = c * X - s * Y
    Ey = s * X + c * Y
    return np.sqrt(S0) * Ex, np.sqrt(S0) * Ey


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

# 그림1: 선형/원형/타원 편광의 전기장 궤적
t = np.linspace(0, 2 * np.pi, 400)
CASES = [("linear", 1.0, 1.0, 0.0), ("circular", 1.0, 1.0, np.pi / 2),
         ("elliptical", 1.0, 0.55, np.pi / 4)]

# 그림2: 완전편광 -> 부분편광 -> 무편광
# 매 순간의 편광 상태는 (완전편광이므로) 포앵카레 구면 위의 한 점으로 표현된다.
# 부분편광은 "고정된 한 점(신호)"과 "구면 위에 고르게 흩어진 무작위 점(잡음)"을
# p:(1-p) 비율로 섞은 앙상블의 시간평균으로 모델링한다.
rng = np.random.default_rng(0)
N = 2000

Ax0, Ay0, dphi0 = 1.0, 0.65, np.pi / 4
S0_ref = Ax0 ** 2 + Ay0 ** 2
s_signal = np.array([
    (Ax0 ** 2 - Ay0 ** 2) / S0_ref,
    (2 * Ax0 * Ay0 * np.cos(dphi0)) / S0_ref,
    (2 * Ax0 * Ay0 * np.sin(dphi0)) / S0_ref,
])


def random_sphere_points(n):
    """포앵카레 구면 위에 고르게 분포하는 무작위 (s1,s2,s3) n개."""
    z = rng.uniform(-1, 1, n)
    az = rng.uniform(0, 2 * np.pi, n)
    r = np.sqrt(1 - z ** 2)
    return np.stack([r * np.cos(az), r * np.sin(az), z], axis=1)


PANELS = []
for key, p in [("full", 1.0), ("partial", 0.5), ("unpol", 0.0)]:
    is_signal = rng.uniform(size=N) < p
    s_rand = random_sphere_points(N)
    s_each = np.where(is_signal[:, None], s_signal[None, :], s_rand)
    dop = float(np.linalg.norm(s_each.mean(axis=0)))
    PANELS.append((key, p, s_each, dop))

# 그림4: 무편광 광원 -> 편광자(0deg) -> 위상지연자(45deg)
S_in = np.array([1.0, 0.0, 0.0, 0.0])
omega_p, omega_c = 0.0, np.pi / 4
phis = np.linspace(0, np.pi, 300)
S_out_all = np.array([
    mueller_retarder(omega_c, phi) @ mueller_polarizer(omega_p) @ S_in
    for phi in phis
])
s_norm = S_out_all[:, 1:] / S_out_all[:, :1]
SNAP_PHIS = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
SNAP_LABELS = [r"$\phi=0$", r"$\phi=\pi/4$", r"$\phi=\pi/2$", r"$\phi=3\pi/4$", r"$\phi=\pi$"]

BOX_KEYS = ["source", "polarizer", "retarder", "sample", "retarder2", "analyzer", "detector"]
BOX_X = [0.5, 2.0, 3.7, 5.6, 7.5, 9.2, 10.7]
BOX_COLOR = [None, "tab:blue", "tab:blue", "tab:orange", "tab:green", "tab:green", None]


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L):
    """주어진 라벨 묶음으로 그림 4개를 그린다."""
    out, F, LF = L["dir"], L["font"], L["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. 선형/원형/타원 편광의 전기장 궤적
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    for ax, (key, Ax, Ay, dphi) in zip(axes, CASES):
        Ex = Ax * np.cos(t)
        Ey = Ay * np.cos(t - dphi)
        ax.plot(Ex, Ey, color="tab:blue", lw=2)
        idx = len(t) // 5
        ax.annotate("", xy=(Ex[idx + 1], Ey[idx + 1]), xytext=(Ex[idx], Ey[idx]),
                    arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=1.8))
        ax.scatter([Ex[0]], [Ey[0]], color="black", zorder=5, s=20)
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_aspect("equal")
        ax.axhline(0, color="gray", lw=0.6)
        ax.axvline(0, color="gray", lw=0.6)
        ax.set_xlabel("$E_x$", fontsize=10)
        ax.set_ylabel("$E_y$", fontsize=10)
        ax.set_title(L[key], fontsize=12, **F)
    fig.suptitle(L["fig1"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-polarization-trajectories.png"), dpi=150)
    plt.close(fig)

    # 그림 2. 완전편광 -> 부분편광 -> 무편광
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    for ax, (key, p, s_each, dop) in zip(axes, PANELS):
        if p < 1.0:
            for s1, s2, s3 in s_each[:400]:
                Ex, Ey = stokes_to_ellipse((1.0, s1, s2, s3), n_pts=120)
                ax.plot(Ex, Ey, color="tab:blue", alpha=0.06, lw=1)
        else:
            Ex, Ey = stokes_to_ellipse((1.0, *s_signal), n_pts=120)
            ax.plot(Ex, Ey, color="tab:blue", alpha=1.0, lw=2)
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.set_aspect("equal")
        ax.axhline(0, color="gray", lw=0.6)
        ax.axvline(0, color="gray", lw=0.6)
        ax.set_xlabel("$E_x$", fontsize=10)
        ax.set_ylabel("$E_y$", fontsize=10)
        ax.set_title(L["dop"].format(L[key], dop), fontsize=11, **F)
    fig.suptitle(L["fig2"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig2-partial-polarization.png"), dpi=150)
    plt.close(fig)

    # 그림 3. PSG - 시편 - PSA 구조와 뮬러 행렬 체인
    fig, ax = plt.subplots(figsize=(11, 3.6))
    y = 0.5
    for x, key, color in zip(BOX_X, BOX_KEYS, BOX_COLOR):
        fc = "white" if color is None else color
        ec = "black" if color is None else color
        ax.add_patch(plt.Rectangle((x - 0.55, y - 0.32), 1.1, 0.64, fill=True,
                                   facecolor=fc, alpha=1.0 if color is None else 0.15,
                                   edgecolor=ec, lw=1.6))
        ax.text(x, y, L[key], ha="center", va="center", fontsize=9.5, **F)
    for i in range(len(BOX_X) - 1):
        ax.annotate("", xy=(BOX_X[i + 1] - 0.55, y), xytext=(BOX_X[i] + 0.55, y),
                    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.3))
    ax.text(1.25, y + 0.55, "$S_{in}$", fontsize=11, ha="center")
    ax.text(6.55, y + 0.55, L["chain"], fontsize=10, ha="center", color="gray", **F)
    ax.text(9.95, y + 0.55, "$S_{out}$", fontsize=11, ha="center")
    ax.annotate(L["psg"], xy=(2.85, y - 0.55), fontsize=9.5, ha="center",
                color="tab:blue", **F)
    ax.annotate(L["psa"], xy=(8.35, y - 0.55), fontsize=9.5, ha="center",
                color="tab:green", **F)
    ax.set_xlim(-0.3, 11.3)
    ax.set_ylim(-0.5, 1.4)
    ax.axis("off")
    ax.set_title(L["fig3"], fontsize=11, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig3-psg-psa-diagram.png"), dpi=150)
    plt.close(fig)

    # 그림 4. 위상지연량에 따른 출력 편광 변화
    fig = plt.figure(figsize=(11, 6.5))
    for i, (phi, lab) in enumerate(zip(SNAP_PHIS, SNAP_LABELS)):
        ax = fig.add_subplot(2, 5, i + 1)
        S_out = mueller_retarder(omega_c, phi) @ mueller_polarizer(omega_p) @ S_in
        Ex, Ey = stokes_to_ellipse(S_out)
        ax.plot(Ex, Ey, color="tab:purple", lw=2)
        ax.set_xlim(-0.8, 0.8)
        ax.set_ylim(-0.8, 0.8)
        ax.set_aspect("equal")
        ax.axhline(0, color="gray", lw=0.5)
        ax.axvline(0, color="gray", lw=0.5)
        ax.set_title(lab, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

    ax_bottom = fig.add_subplot(2, 1, 2)
    ax_bottom.plot(phis, s_norm[:, 0], label="$S_1/S_0$", lw=2)
    ax_bottom.plot(phis, s_norm[:, 1], label="$S_2/S_0$", lw=2)
    ax_bottom.plot(phis, s_norm[:, 2], label="$S_3/S_0$", lw=2)
    for phi in SNAP_PHIS:
        ax_bottom.axvline(phi, color="gray", lw=0.6, ls="--")
    ax_bottom.set_xlabel(L["phi_label"], fontsize=11, **F)
    ax_bottom.set_ylabel(L["stokes_norm"], fontsize=11, **F)
    ax_bottom.set_xlim(0, np.pi)
    ax_bottom.legend(fontsize=10)
    ax_bottom.grid(alpha=0.3)
    fig.suptitle(L["fig4"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig4-retarder-phase-sweep.png"), dpi=150)
    plt.close(fig)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

# ---------------------------------------------------------------------------
# 본문에 인용할 수치 출력
# ---------------------------------------------------------------------------
print()
print("=== 그림2: 편광도 P ===")
for key, p, _, dop in PANELS:
    print(f"  {key:8s} (신호 비율 p={p:.2f}) -> 측정된 P={dop:.4f}")
print()
print("=== 그림4: S_out = M_C(45deg, phi) M_P(0deg) S_in ===")
for phi, lab in zip(SNAP_PHIS, SNAP_LABELS):
    S = mueller_retarder(omega_c, phi) @ mueller_polarizer(omega_p) @ S_in
    print(f"  phi={phi:.4f} rad -> S=[{S[0]:.3f}, {S[1]:+.3f}, {S[2]:+.3f}, {S[3]:+.3f}]")
print(f"  S_2 는 모든 phi 에서 0 인가: {np.allclose(S_out_all[:, 2], 0)}")

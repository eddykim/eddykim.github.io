"""2편 그림 4개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometry-polarization-mueller-matrix/ 에 fig1~4 저장
"""
import os

import matplotlib.pyplot as plt
import numpy as np

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometry-polarization-mueller-matrix",
)
os.makedirs(OUT_DIR, exist_ok=True)


# ────────────────────────────────────────────────────────────
# 뮬러 계산법 기본 요소 (식 2.21~2.24)
# ────────────────────────────────────────────────────────────
def mueller_rotation(omega):
    """회전 행렬 M_R(omega), omega: rad (식 2.23)."""
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([
        [1, 0, 0, 0],
        [0, c, s, 0],
        [0, -s, c, 0],
        [0, 0, 0, 1],
    ])


def mueller_polarizer(omega_p):
    """이상적 편광자 M_P(omega_p) (식 2.21)."""
    base = 0.5 * np.array([
        [1, 1, 0, 0],
        [1, 1, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0],
    ])
    return mueller_rotation(-omega_p) @ base @ mueller_rotation(omega_p)


def mueller_retarder(omega_c, phi):
    """위상 지연자 M_C(omega_c, phi) (식 2.22)."""
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
    t = np.linspace(0, 2 * np.pi, n_pts)
    X, Y = a * np.cos(t), b * np.sin(t)
    c, s = np.cos(psi), np.sin(psi)
    Ex = c * X - s * Y
    Ey = s * X + c * Y
    scale = np.sqrt(S0)
    return scale * Ex, scale * Ey


# ────────────────────────────────────────────────────────────
# 그림 1. 선형/원형/타원 편광의 전기장 궤적 (식 2.18)
# ────────────────────────────────────────────────────────────
t = np.linspace(0, 2 * np.pi, 400)
cases = [
    ("선형편광", 1.0, 1.0, 0.0),
    ("원형편광", 1.0, 1.0, np.pi / 2),
    ("타원편광", 1.0, 0.55, np.pi / 4),
]

fig, axes = plt.subplots(1, 3, figsize=(11, 4))
for ax, (name, Ax, Ay, dphi) in zip(axes, cases):
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
    ax.set_title(name, fontsize=12, **KFONT)

fig.suptitle("전기장 성분의 진폭비·위상차에 따른 편광 상태 (식 2.18)", fontsize=12, **KFONT)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1-polarization-trajectories.png"), dpi=150)
plt.close(fig)


# ────────────────────────────────────────────────────────────
# 그림 2. 완전편광 -> 부분편광 -> 무편광: 존스 벡터가 표현 못 하는 영역
#
# 매 순간의 편광 상태는 (완전편광이므로) 포앵카레 구면 위의 한 점으로
# 표현된다. 부분편광은 "고정된 한 점(신호)"과 "구면 위에 고르게 흩어진
# 무작위 점(잡음)"을 p:(1-p) 비율로 섞은 앙상블의 시간평균으로 모델링한다.
# ────────────────────────────────────────────────────────────
rng = np.random.default_rng(0)
N = 2000

# 기준이 되는 순수 타원편광 (Ax=1, Ay=0.65, dphi=pi/4)을 정규화 스토크스 성분으로 변환
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


panels = [("완전편광", 1.0), ("부분편광", 0.5), ("무편광", 0.0)]

fig, axes = plt.subplots(1, 3, figsize=(11, 4))
for ax, (label, p) in zip(axes, panels):
    is_signal = rng.uniform(size=N) < p
    s_rand = random_sphere_points(N)
    s_each = np.where(is_signal[:, None], s_signal[None, :], s_rand)

    if p < 1.0:
        for s1, s2, s3 in s_each[:400]:
            Ex, Ey = stokes_to_ellipse((1.0, s1, s2, s3), n_pts=120)
            ax.plot(Ex, Ey, color="tab:blue", alpha=0.06, lw=1)
    else:
        Ex, Ey = stokes_to_ellipse((1.0, *s_signal), n_pts=120)
        ax.plot(Ex, Ey, color="tab:blue", alpha=1.0, lw=2)

    s_avg = s_each.mean(axis=0)
    dop = np.linalg.norm(s_avg)

    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axhline(0, color="gray", lw=0.6)
    ax.axvline(0, color="gray", lw=0.6)
    ax.set_xlabel("$E_x$", fontsize=10)
    ax.set_ylabel("$E_y$", fontsize=10)
    ax.set_title(f"{label}\n(측정된 P={dop:.2f})", fontsize=11, **KFONT)

fig.suptitle("존스 벡터로는 왼쪽 한 장만 표현 가능 — 스토크스 벡터는 셋 다 표현", fontsize=12, **KFONT)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig2-partial-polarization.png"), dpi=150)
plt.close(fig)


# ────────────────────────────────────────────────────────────
# 그림 3. PSG - 시편 - PSA 구조와 뮬러 행렬 체인 (그림 1.1, 2.3 재구성)
# ────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 3.6))
boxes = [
    (0.5, "광원", None),
    (2.0, "편광자\n$M_P(\\Omega_P)$", "tab:blue"),
    (3.7, "위상지연자\n$M_C(\\Omega_C,\\phi)$", "tab:blue"),
    (5.6, "시편\n$M_{sample}$", "tab:orange"),
    (7.5, "위상지연자\n$M_C'$", "tab:green"),
    (9.2, "편광자\n$M_A$", "tab:green"),
    (10.7, "검출기", None),
]
y = 0.5
for x, label, color in boxes:
    fc = "white" if color is None else color
    alpha = 1.0 if color is None else 0.15
    ec = "black" if color is None else color
    ax.add_patch(plt.Rectangle((x - 0.55, y - 0.32), 1.1, 0.64, fill=True,
                                facecolor=fc, alpha=alpha if color else 1.0,
                                edgecolor=ec, lw=1.6))
    ax.text(x, y, label, ha="center", va="center", fontsize=9.5, **KFONT)

for i in range(len(boxes) - 1):
    x0, x1 = boxes[i][0] + 0.55, boxes[i + 1][0] - 0.55
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.3))

ax.text(1.25, y + 0.55, "$S_{in}$", fontsize=11, ha="center")
ax.text(6.55, y + 0.55, "PSG → 시편 → PSA", fontsize=10, ha="center", color="gray", **KFONT)
ax.text(9.95, y + 0.55, "$S_{out}$", fontsize=11, ha="center")

ax.annotate("편광 생성단 (PSG)", xy=(2.85, y - 0.55), fontsize=9.5, ha="center",
            color="tab:blue", **KFONT)
ax.annotate("편광 분석단 (PSA)", xy=(8.35, y - 0.55), fontsize=9.5, ha="center",
            color="tab:green", **KFONT)

ax.set_xlim(-0.3, 11.3)
ax.set_ylim(-0.5, 1.4)
ax.axis("off")
ax.set_title("PSG-시편-PSA 구조와 뮬러 행렬의 순차 곱 $S_{out}=M_A M_C' M_{sample} M_C M_P S_{in}$",
             fontsize=11, **KFONT)

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig3-psg-psa-diagram.png"), dpi=150)
plt.close(fig)


# ────────────────────────────────────────────────────────────
# 그림 4. 편광자(0°) + 위상지연자(45°) 통과 후, 위상지연량에 따른 출력 편광 변화
# ────────────────────────────────────────────────────────────
S_in = np.array([1.0, 0.0, 0.0, 0.0])  # 무편광 광원
omega_p, omega_c = 0.0, np.pi / 4  # 편광자 0deg, 위상지연자 45deg

phis = np.linspace(0, np.pi, 300)
S_out_all = np.array([
    mueller_retarder(omega_c, phi) @ mueller_polarizer(omega_p) @ S_in
    for phi in phis
])
s_norm = S_out_all[:, 1:] / S_out_all[:, :1]

fig = plt.figure(figsize=(11, 6.5))
snap_phis = [0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]
snap_labels = [r"$\phi=0$", r"$\phi=\pi/4$", r"$\phi=\pi/2$", r"$\phi=3\pi/4$", r"$\phi=\pi$"]

for i, (phi, lab) in enumerate(zip(snap_phis, snap_labels)):
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
    ax.set_xticks([]); ax.set_yticks([])

ax_bottom = fig.add_subplot(2, 1, 2)
ax_bottom.plot(phis, s_norm[:, 0], label="$S_1/S_0$", lw=2)
ax_bottom.plot(phis, s_norm[:, 1], label="$S_2/S_0$", lw=2)
ax_bottom.plot(phis, s_norm[:, 2], label="$S_3/S_0$", lw=2)
for phi in snap_phis:
    ax_bottom.axvline(phi, color="gray", lw=0.6, ls="--")
ax_bottom.set_xlabel(r"위상 지연량 $\phi$ (rad)", fontsize=11, **KFONT)
ax_bottom.set_ylabel("정규화된 스토크스 성분", fontsize=11, **KFONT)
ax_bottom.set_xlim(0, np.pi)
ax_bottom.legend(fontsize=10)
ax_bottom.grid(alpha=0.3)

fig.suptitle("편광자(0°)+위상지연자(45°) 통과 후 출력 편광: 선형 → 타원 → 원형 → 타원 → 선형",
             fontsize=12, **KFONT)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig4-retarder-phase-sweep.png"), dpi=150)
plt.close(fig)

print("저장 완료:", OUT_DIR)

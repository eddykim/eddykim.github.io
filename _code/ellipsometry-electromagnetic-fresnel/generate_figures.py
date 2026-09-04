"""1편 그림 2개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometry-electromagnetic-fresnel/ 에 fig1~2 저장
"""
import os

import matplotlib.pyplot as plt
import numpy as np

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometry-electromagnetic-fresnel",
)
os.makedirs(OUT_DIR, exist_ok=True)


def fresnel_coefficients(n1, n2, theta_i):
    """s/p 편광 프레넬 반사·투과계수. theta_i: 입사각(rad), n1/n2: 실수/복소 굴절률."""
    sin_t = (n1 / n2) * np.sin(theta_i)
    cos_t = np.sqrt(1 - sin_t ** 2 + 0j)
    cos_i = np.cos(theta_i)

    r_p = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    return r_p, r_s


# ────────────────────────────────────────────────────────────
# 그림 1. 입사/반사/굴절 벡터와 s/p 편광 정의 개념도
# ────────────────────────────────────────────────────────────
theta_i_deg = 35.0
theta_t_deg = 22.0  # n2 > n1 이므로 법선 쪽으로 굽음 (도식용 예시각)
ti, tt = np.radians(theta_i_deg), np.radians(theta_t_deg)

fig, ax = plt.subplots(figsize=(6.4, 5.4))

# 매질 경계면과 법선
ax.axhline(0, color="black", lw=1.5)
ax.plot([0, 0], [-1.3, 1.3], color="gray", lw=1, ls=(0, (5, 4)))
ax.text(1.35, -0.08, "경계면", ha="left", va="top", fontsize=10, **KFONT)
ax.text(0.05, 1.28, "법선", ha="left", va="top", fontsize=10, color="gray", **KFONT)

origin = np.array([0, 0])
L = 1.15

# 입사광: 왼쪽 위(-x,+y)에서 출발해 원점(경계면)에 도달
p_inc = L * np.array([-np.sin(ti), np.cos(ti)])
# 반사광: 원점에서 오른쪽 위(+x,+y)로 진행
p_ref = L * np.array([np.sin(ti), np.cos(ti)])
# 굴절광: 원점에서 오른쪽 아래(+x,-y)로 진행 (n2>n1이므로 법선 쪽으로 굽음)
p_tra = L * np.array([np.sin(tt), -np.cos(tt)])

ax.annotate("", xy=origin, xytext=p_inc,
            arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=2))
ax.annotate("", xy=p_ref, xytext=origin,
            arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=2))
ax.annotate("", xy=p_tra, xytext=origin,
            arrowprops=dict(arrowstyle="-|>", color="tab:green", lw=2))

ax.text(*(p_inc * 1.12), r"입사광 ($\theta_i$)", color="tab:blue", fontsize=10,
        ha="center", va="bottom", **KFONT)
ax.text(*(p_ref * 1.12), r"반사광 ($\theta_r$)", color="tab:red", fontsize=10,
        ha="center", va="bottom", **KFONT)
ax.text(*(p_tra * 1.16), r"굴절광 ($\theta_t$)", color="tab:green", fontsize=10,
        ha="center", va="top", **KFONT)

# 각도 호 (법선 기준, 각 광선이 있는 사분면에 맞춰 그린다)
arc_ti = np.linspace(np.pi / 2, np.pi / 2 + ti, 30)         # 입사광: 좌상단
ax.plot(0.35 * np.cos(arc_ti), 0.35 * np.sin(arc_ti), color="tab:blue", lw=1.2)

arc_tr = np.linspace(np.pi / 2 - ti, np.pi / 2, 30)         # 반사광: 우상단
ax.plot(0.45 * np.cos(arc_tr), 0.45 * np.sin(arc_tr), color="tab:red", lw=1.2)

arc_tt = np.linspace(-np.pi / 2, -np.pi / 2 + tt, 30)       # 굴절광: 우하단
ax.plot(0.35 * np.cos(arc_tt), 0.35 * np.sin(arc_tt), color="tab:green", lw=1.2)

ax.text(-0.13, 0.30, r"$\theta_i$", color="tab:blue", fontsize=11)
ax.text(0.16, 0.42, r"$\theta_r$", color="tab:red", fontsize=11)
ax.text(0.10, -0.30, r"$\theta_t$", color="tab:green", fontsize=11)

# 매질 표시
ax.text(-1.45, 0.9, r"매질 1 ($N_1$)", fontsize=11, **KFONT)
ax.text(-1.45, -0.9, r"매질 2 ($N_2$)", fontsize=11, **KFONT)

# s/p 편광 기호 (입사광 중간 지점, 각도 표시와 겹치지 않는 위치에 표시)
mid_inc = p_inc * 0.72
dir_inc = p_inc / np.linalg.norm(p_inc)
perp_inc = np.array([-dir_inc[1], dir_inc[0]])

# s편광: 지면에서 튀어나오는 방향 (동그라미+점), 광선 왼쪽에 배치
s_pt = mid_inc - 0.16 * perp_inc
ax.scatter(*s_pt, s=90, facecolors="none", edgecolors="black", zorder=5)
ax.scatter(*s_pt, s=8, color="black", zorder=5)
ax.text(*(s_pt + np.array([-0.10, 0.0])), "s-편광\n(지면에 수직)", fontsize=8.5,
        ha="right", va="center", **KFONT)

# p편광: 광선에 수직하고 지면 내에 있는 화살표(양방향), 광선 오른쪽에 배치
p_base = mid_inc + 0.16 * perp_inc
ax.annotate("", xy=p_base + 0.14 * perp_inc, xytext=p_base - 0.14 * perp_inc,
            arrowprops=dict(arrowstyle="<->", color="black", lw=1.4))
ax.text(*(p_base + 0.24 * perp_inc), "p-편광\n(입사면 내)", fontsize=8.5,
        ha="left", va="center", **KFONT)

ax.set_xlim(-1.7, 1.7)
ax.set_ylim(-1.5, 1.5)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("입사면 위의 s/p 편광과 반사·굴절 벡터 정의", fontsize=12, **KFONT)

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1-concept-diagram.png"), dpi=150)
plt.close(fig)


# ────────────────────────────────────────────────────────────
# 그림 2. 프레넬 반사율 R_p, R_s vs 입사각, 브루스터각 표시
# ────────────────────────────────────────────────────────────
n1, n2 = 1.0, 1.5  # 공기 -> 유리(예시)
angles_deg = np.linspace(0, 89.9, 500)
angles_rad = np.radians(angles_deg)

r_p, r_s = fresnel_coefficients(n1, n2, angles_rad)
R_p, R_s = np.abs(r_p) ** 2, np.abs(r_s) ** 2

theta_brewster_deg = np.degrees(np.arctan(n2 / n1))

fig2, ax2 = plt.subplots(figsize=(6.6, 4.6))
ax2.plot(angles_deg, R_p, color="tab:green", lw=2, label=r"$R_p$")
ax2.plot(angles_deg, R_s, color="tab:purple", lw=2, label=r"$R_s$")
ax2.axvline(theta_brewster_deg, color="gray", lw=1, ls="--")
ax2.annotate(
    f"브루스터각\n{theta_brewster_deg:.1f}°",
    xy=(theta_brewster_deg, 0.0), xytext=(theta_brewster_deg + 6, 0.18),
    fontsize=9, color="gray", **KFONT,
    arrowprops=dict(arrowstyle="->", color="gray", lw=1),
)

ax2.set_xlabel("입사각 " + r"$\theta_i$ (deg)", fontsize=11, **KFONT)
ax2.set_ylabel("반사율 R", fontsize=11, **KFONT)
ax2.set_title(f"프레넬 반사율 ($N_1$={n1}, $N_2$={n2})", fontsize=12, **KFONT)
ax2.set_xlim(0, 90)
ax2.set_ylim(0, 1)
ax2.legend(prop=LFONT, fontsize=11)
ax2.grid(alpha=0.3)

fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "fig2-fresnel-reflectance.png"), dpi=150)
plt.close(fig2)

print(f"브루스터각: {theta_brewster_deg:.3f} deg (R_p={np.interp(theta_brewster_deg, angles_deg, R_p):.2e})")
print("저장 완료:", OUT_DIR)

"""엘립소메트리 배경이론 1편 그림 2개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometry-electromagnetic-fresnel/      (한국어)
      ../../assets/img/posts/ellipsometry-electromagnetic-fresnel/en/   (영문)

계산은 한 번만 수행하고 라벨 문자열만 갈아 끼우므로, 두 언어의 그림은 데이터가
완전히 동일하고 표기만 다르다.
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
    "assets", "img", "posts", "ellipsometry-electromagnetic-fresnel",
)

LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "interface": "경계면", "normal": "법선",
        "incident": r"입사광 ($\theta_i$)", "reflected": r"반사광 ($\theta_r$)",
        "transmitted": r"굴절광 ($\theta_t$)",
        "medium1": r"매질 1 ($N_1$)", "medium2": r"매질 2 ($N_2$)",
        "s_pol": "s-편광\n(지면에 수직)", "p_pol": "p-편광\n(입사면 내)",
        "fig1": "입사면 위의 s/p 편광과 반사·굴절 벡터 정의",
        "brewster": "브루스터각\n{:.1f}°",
        "xlabel": "입사각 " + r"$\theta_i$ (deg)", "ylabel": "반사율 R",
        "fig2": "프레넬 반사율 ($N_1$={}, $N_2$={})",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "interface": "interface", "normal": "normal",
        "incident": r"incident ($\theta_i$)", "reflected": r"reflected ($\theta_r$)",
        "transmitted": r"refracted ($\theta_t$)",
        "medium1": r"medium 1 ($N_1$)", "medium2": r"medium 2 ($N_2$)",
        "s_pol": "s-polarized\n(normal to the page)", "p_pol": "p-polarized\n(in the plane of incidence)",
        "fig1": "s/p polarization and the reflected and refracted vectors on the plane of incidence",
        "brewster": "Brewster angle\n{:.1f}°",
        "xlabel": r"Angle of incidence $\theta_i$ (deg)", "ylabel": "Reflectance R",
        "fig2": "Fresnel reflectance ($N_1$={}, $N_2$={})",
    },
}


def fresnel_coefficients(n1, n2, theta_i):
    """s/p 편광 프레넬 반사·투과계수. theta_i: 입사각(rad), n1/n2: 실수/복소 굴절률."""
    sin_t = (n1 / n2) * np.sin(theta_i)
    cos_t = np.sqrt(1 - sin_t ** 2 + 0j)
    cos_i = np.cos(theta_i)

    r_p = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t)
    r_s = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t)
    return r_p, r_s


# ---------------------------------------------------------------------------
# 계산 (언어와 무관하게 한 번만 수행한다)
# ---------------------------------------------------------------------------

# 그림1: 개념도용 각도 (도식용 예시각, n2 > n1 이므로 법선 쪽으로 굽는다)
theta_i_deg, theta_t_deg = 35.0, 22.0
ti, tt = np.radians(theta_i_deg), np.radians(theta_t_deg)
L = 1.15
origin = np.array([0, 0])
p_inc = L * np.array([-np.sin(ti), np.cos(ti)])   # 왼쪽 위에서 원점으로
p_ref = L * np.array([np.sin(ti), np.cos(ti)])    # 원점에서 오른쪽 위로
p_tra = L * np.array([np.sin(tt), -np.cos(tt)])   # 원점에서 오른쪽 아래로

# 그림2: 공기 -> 유리(예시)에서의 프레넬 반사율
n1, n2 = 1.0, 1.5
angles_deg = np.linspace(0, 89.9, 500)
angles_rad = np.radians(angles_deg)
r_p, r_s = fresnel_coefficients(n1, n2, angles_rad)
R_p, R_s = np.abs(r_p) ** 2, np.abs(r_s) ** 2
theta_brewster_deg = np.degrees(np.arctan(n2 / n1))


# ---------------------------------------------------------------------------
# 그리기
# ---------------------------------------------------------------------------

def render(L_):
    """주어진 라벨 묶음으로 그림 2개를 그린다."""
    out, F, LF = L_["dir"], L_["font"], L_["legend"]
    os.makedirs(out, exist_ok=True)

    # 그림 1. 입사/반사/굴절 벡터와 s/p 편광 정의 개념도
    fig, ax = plt.subplots(figsize=(6.4, 5.4))

    ax.axhline(0, color="black", lw=1.5)
    ax.plot([0, 0], [-1.3, 1.3], color="gray", lw=1, ls=(0, (5, 4)))
    ax.text(1.35, -0.08, L_["interface"], ha="left", va="top", fontsize=10, **F)
    ax.text(0.05, 1.28, L_["normal"], ha="left", va="top", fontsize=10, color="gray", **F)

    ax.annotate("", xy=origin, xytext=p_inc,
                arrowprops=dict(arrowstyle="-|>", color="tab:blue", lw=2))
    ax.annotate("", xy=p_ref, xytext=origin,
                arrowprops=dict(arrowstyle="-|>", color="tab:red", lw=2))
    ax.annotate("", xy=p_tra, xytext=origin,
                arrowprops=dict(arrowstyle="-|>", color="tab:green", lw=2))

    ax.text(*(p_inc * 1.12), L_["incident"], color="tab:blue", fontsize=10,
            ha="center", va="bottom", **F)
    ax.text(*(p_ref * 1.12), L_["reflected"], color="tab:red", fontsize=10,
            ha="center", va="bottom", **F)
    ax.text(*(p_tra * 1.16), L_["transmitted"], color="tab:green", fontsize=10,
            ha="center", va="top", **F)

    # 각도 호 (법선 기준, 각 광선이 있는 사분면에 맞춰 그린다)
    arc_ti = np.linspace(np.pi / 2, np.pi / 2 + ti, 30)
    ax.plot(0.35 * np.cos(arc_ti), 0.35 * np.sin(arc_ti), color="tab:blue", lw=1.2)
    arc_tr = np.linspace(np.pi / 2 - ti, np.pi / 2, 30)
    ax.plot(0.45 * np.cos(arc_tr), 0.45 * np.sin(arc_tr), color="tab:red", lw=1.2)
    arc_tt = np.linspace(-np.pi / 2, -np.pi / 2 + tt, 30)
    ax.plot(0.35 * np.cos(arc_tt), 0.35 * np.sin(arc_tt), color="tab:green", lw=1.2)

    ax.text(-0.13, 0.30, r"$\theta_i$", color="tab:blue", fontsize=11)
    ax.text(0.16, 0.42, r"$\theta_r$", color="tab:red", fontsize=11)
    ax.text(0.10, -0.30, r"$\theta_t$", color="tab:green", fontsize=11)

    ax.text(-1.45, 0.9, L_["medium1"], fontsize=11, **F)
    ax.text(-1.45, -0.9, L_["medium2"], fontsize=11, **F)

    # s/p 편광 기호 (입사광 중간 지점, 각도 표시와 겹치지 않는 위치)
    mid_inc = p_inc * 0.72
    dir_inc = p_inc / np.linalg.norm(p_inc)
    perp_inc = np.array([-dir_inc[1], dir_inc[0]])

    s_pt = mid_inc - 0.16 * perp_inc  # 지면에서 튀어나오는 방향 (동그라미+점)
    ax.scatter(*s_pt, s=90, facecolors="none", edgecolors="black", zorder=5)
    ax.scatter(*s_pt, s=8, color="black", zorder=5)
    ax.text(*(s_pt + np.array([-0.10, 0.0])), L_["s_pol"], fontsize=8.5,
            ha="right", va="center", **F)

    p_base = mid_inc + 0.16 * perp_inc  # 광선에 수직하고 지면 내에 있는 양방향 화살표
    ax.annotate("", xy=p_base + 0.14 * perp_inc, xytext=p_base - 0.14 * perp_inc,
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.4))
    ax.text(*(p_base + 0.24 * perp_inc), L_["p_pol"], fontsize=8.5,
            ha="left", va="center", **F)

    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(L_["fig1"], fontsize=12, **F)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "fig1-concept-diagram.png"), dpi=150)
    plt.close(fig)

    # 그림 2. 프레넬 반사율 R_p, R_s vs 입사각, 브루스터각 표시
    fig2, ax2 = plt.subplots(figsize=(6.6, 4.6))
    ax2.plot(angles_deg, R_p, color="tab:green", lw=2, label=r"$R_p$")
    ax2.plot(angles_deg, R_s, color="tab:purple", lw=2, label=r"$R_s$")
    ax2.axvline(theta_brewster_deg, color="gray", lw=1, ls="--")
    ax2.annotate(
        L_["brewster"].format(theta_brewster_deg),
        xy=(theta_brewster_deg, 0.0), xytext=(theta_brewster_deg + 6, 0.18),
        fontsize=9, color="gray", **F,
        arrowprops=dict(arrowstyle="->", color="gray", lw=1),
    )
    ax2.set_xlabel(L_["xlabel"], fontsize=11, **F)
    ax2.set_ylabel(L_["ylabel"], fontsize=11, **F)
    ax2.set_title(L_["fig2"].format(n1, n2), fontsize=12, **F)
    ax2.set_xlim(0, 90)
    ax2.set_ylim(0, 1)
    ax2.legend(prop=LF, fontsize=11)
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, "fig2-fresnel-reflectance.png"), dpi=150)
    plt.close(fig2)


for lang, labels in LABELS.items():
    print(f"--- [{lang}] {os.path.normpath(labels['dir'])} ---")
    render(labels)

print()
print(f"브루스터각: {theta_brewster_deg:.3f} deg "
      f"(R_p={np.interp(theta_brewster_deg, angles_deg, R_p):.2e})")
print(f"수직 입사 반사율 R(0도) = {R_p[0]:.5f} "
      f"(= ((n2-n1)/(n2+n1))^2 = {((n2 - n1) / (n2 + n1)) ** 2:.5f})")

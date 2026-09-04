"""3편 그림 3개 생성.

실행: python generate_figures.py
출력: ../../assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/ 에 fig1~3 저장
"""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from tmm import reflectance_spectrum
from smm_tensor import smm_reflectance_tensor

KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}

OUT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "..",
    "assets", "img", "posts", "ellipsometry-thin-film-multilayer-reflectance",
)
os.makedirs(OUT_DIR, exist_ok=True)


# ────────────────────────────────────────────────────────────
# 그림 1. 단층 박막 다중반사 광선 경로 개념도 (논문 그림 2.5 재구성)
# ────────────────────────────────────────────────────────────
n0, n1, n2 = 1.0, 1.5, 3.5
theta1_deg = 22.0
theta1 = np.radians(theta1_deg)
theta0 = np.arcsin((n1 / n0) * np.sin(theta1))
theta2 = theta1 - np.radians(6.0)  # 기판에서 조금 더 굽는 것처럼 보이도록(도식용)

d = 1.0  # 박막 두께 (도식 단위)
dx = d * np.tan(theta1)  # 내부 반사 한 번당 수평 이동량

# 상부(공기/박막) 계면 위의 점들: O, B, D  |  하부(박막/기판) 계면 위의 점들: A, C
O = np.array([0.0, 0.0])
B = np.array([2 * dx, 0.0])
D = np.array([4 * dx, 0.0])
A = np.array([dx, -d])
C = np.array([3 * dx, -d])

fig, ax = plt.subplots(figsize=(7.2, 5.6))

# 매질 영역 표시
ax.axhspan(0, 1.4, color="#eaf4ff", zorder=0)
ax.axhspan(-d, 0, color="#fff3d6", zorder=0)
ax.axhspan(-d - 0.9, -d, color="#e9e9e9", zorder=0)
ax.axhline(0, color="black", lw=1.6)
ax.axhline(-d, color="black", lw=1.6)

ax.text(4.6 * dx, 1.15, r"Air : $N_0$", fontsize=11, **KFONT)
ax.text(4.6 * dx, -d / 2, r"Layer : $N_1$", fontsize=11, va="center", **KFONT)
ax.text(4.6 * dx, -d - 0.55, r"Substrate : $N_2$", fontsize=11, **KFONT)

L_in = 1.0 * d / np.cos(theta0)  # 입사광 표시 길이


def ray(p_from, p_to, color, lw=1.8, style="-"):
    ax.annotate("", xy=p_to, xytext=p_from,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=style))


# 입사광 -> O, O에서의 1차 반사 Er1, O에서 A로의 굴절
p_inc_start = O + L_in * np.array([-np.sin(theta0), np.cos(theta0)])
ray(p_inc_start, O, "tab:blue")
ax.text(*(p_inc_start + [0.06, 0.03]), r"$E_i$", color="tab:blue", fontsize=12, **KFONT)

er1_end = O + L_in * np.array([np.sin(theta0), np.cos(theta0)])
ray(O, er1_end, "tab:red")
ax.text(*(er1_end + [0.03, 0.02]), r"$E_{r1}$", color="tab:red", fontsize=11, **KFONT)

ray(O, A, "tab:green")  # O -> A (박막 내부, 최초 굴절)

# A: 박막/기판 계면 -> 기판으로 투과(Et1) + 박막 내부로 반사(-> B)
et1_end = A + 0.85 * np.array([np.sin(theta2), -np.cos(theta2)])
ray(A, et1_end, "tab:orange")
ax.text(*(et1_end + [0.03, -0.03]), r"$E_{t1}$", color="tab:orange", fontsize=10, **KFONT)
ray(A, B, "tab:green")

# B: 박막/공기 계면 -> 공기로 투과(Er2, Er1과 평행) + 박막 내부로 반사(-> C)
er2_end = B + L_in * np.array([np.sin(theta0), np.cos(theta0)])
ray(B, er2_end, "tab:red")
ax.text(*(er2_end + [0.03, 0.02]), r"$E_{r2}$", color="tab:red", fontsize=11, **KFONT)
ray(B, C, "tab:green")

# C: 박막/기판 계면 -> 기판으로 투과(Et2) + 박막 내부로 반사(-> D)
et2_end = C + 0.85 * np.array([np.sin(theta2), -np.cos(theta2)])
ray(C, et2_end, "tab:orange")
ax.text(*(et2_end + [0.03, -0.03]), r"$E_{t2}$", color="tab:orange", fontsize=10, **KFONT)
ray(C, D, "tab:green")

# D: 박막/공기 계면 -> 공기로 투과(Er3, Er1과 평행)
er3_end = D + L_in * np.array([np.sin(theta0), np.cos(theta0)])
ray(D, er3_end, "tab:red")
ax.text(*(er3_end + [0.03, 0.02]), r"$E_{r3}$", color="tab:red", fontsize=11, **KFONT)

# 광경로차 주석: O->A->B(박막 내부, 실제 경로) vs O->B(상부계면을 따르는 기준 경로)
ref_y = 0.22
ax.plot([O[0], B[0]], [ref_y, ref_y], color="gray", lw=1.2, ls=(0, (4, 3)))
ax.text((O[0] + B[0]) / 2, ref_y + 0.06, "기준 경로 (OB, 공기 중 환산)", color="gray",
        fontsize=8.5, ha="center", **KFONT)
ax.text((O[0] + A[0] + B[0]) / 3, -d * 0.62,
        "실제 경로 O→A→B\n(박막 내부, 길이 $\\propto N_1$)", color="tab:green",
        fontsize=8.5, ha="center", **KFONT)
ax.annotate(r"위상차 $2\beta_{phase}$", xy=(B[0] + 0.15, 0.42),
            fontsize=10.5, color="black", **KFONT)

label_offsets = {"O": (-0.10, 0.07), "A": (-0.10, -0.10), "B": (-0.10, 0.07), "C": (-0.10, -0.10)}
for pt, name in [(O, "O"), (A, "A"), (B, "B"), (C, "C")]:
    ax.scatter(*pt, s=22, color="black", zorder=5)
    ox, oy = label_offsets[name]
    ax.text(pt[0] + ox, pt[1] + oy, name, fontsize=11, ha="center", zorder=6)

ax.text(-0.55, 0.25, r"$\theta_0$", fontsize=11, color="tab:blue")
ax.text(dx * 0.35, -0.28, r"$\theta_1$", fontsize=11, color="tab:green")

ax.set_xlim(-1.2, 4.6 * dx + 0.6)
ax.set_ylim(-d - 1.0, 1.4)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("단층 박막에서의 다중반사와 광경로차", fontsize=12, **KFONT)

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "fig1-concept-diagram.png"), dpi=150)
plt.close(fig)


# ────────────────────────────────────────────────────────────
# 그림 2. 로아드/TMM/SMM 세 방법의 계산 구조 비교
# ────────────────────────────────────────────────────────────
def box(ax, xy, w, h, text, fc, ec="black", fontsize=10):
    b = FancyBboxPatch(xy, w, h, boxstyle="round,pad=0.02,rounding_size=0.04",
                        linewidth=1.3, edgecolor=ec, facecolor=fc, zorder=3)
    ax.add_patch(b)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center", va="center",
            fontsize=fontsize, zorder=4, **KFONT)


def arrow(ax, p0, p1, color="black", style="-|>"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, color=color,
                                  lw=1.6, mutation_scale=14, zorder=2))


fig2, axes = plt.subplots(1, 3, figsize=(12.6, 5.2))

# (a) 로아드 방법: 아래(기판)에서 위(공기)로 등가계면을 순차 병합
ax = axes[0]
labels_bottom_up = [r"$r_{m-1,m}$" + "\n(기판 경계)", r"$\rho_m$", r"$\rho_{m-1}$", r"$r_{total}$" + "\n(공기 경계)"]
ys = [0.05, 0.35, 0.65, 0.95]
for y, lab in zip(ys, labels_bottom_up):
    box(ax, (0.32, y), 0.36, 0.16, lab, fc="#fff3d6", fontsize=9.5)
for y0, y1 in zip(ys[:-1], ys[1:]):
    arrow(ax, (0.5, y0 + 0.16), (0.5, y1), color="tab:orange")
ax.text(0.5, 1.16, "로아드 방법", ha="center", fontsize=13, weight="bold", **KFONT)
ax.text(0.5, -0.12, "하부층 → 상부층\n등가계면으로 순차 치환", ha="center",
        fontsize=9, color="#555", **KFONT)

# (b) TMM: D0^-1, Q1, Q2, ..., Dsub 를 한 줄로 곱해나가는 체인
ax = axes[1]
chain = [r"$D_0^{-1}$", r"$Q_1$", r"$Q_2$", r"$\cdots$", r"$D_{sub}$"]
xs = np.linspace(0.05, 0.75, len(chain))
for x, lab in zip(xs, chain):
    box(ax, (x, 0.45), 0.16, 0.18, lab, fc="#eaf4ff", fontsize=10.5)
for x0, x1 in zip(xs[:-1], xs[1:]):
    arrow(ax, (x0 + 0.16, 0.54), (x1, 0.54), color="tab:blue")
ax.text(0.5, 1.16, "전달행렬법 (TMM)", ha="center", fontsize=13, weight="bold", **KFONT)
ax.text(0.5, 0.20, "행렬을 한 줄로 순서대로 곱함\n(위상+계면이 한 행렬 $Q_m$에 결합)",
        ha="center", fontsize=9, color="#555", **KFONT)

# (c) SMM: 계면행렬(I)과 층행렬(L)이 분리되어 번갈아 곱해짐
ax = axes[2]
chain = [(r"$I_{01}$", "#ffe0e0"), (r"$L_1$", "#e6ffe6"), (r"$I_{12}$", "#ffe0e0"),
         (r"$L_2$", "#e6ffe6"), (r"$I_{2,sub}$", "#ffe0e0")]
xs = np.linspace(0.02, 0.78, len(chain))
for x, (lab, fc) in zip(xs, chain):
    box(ax, (x, 0.45), 0.18, 0.18, lab, fc=fc, fontsize=10)
for x0, x1 in zip(xs[:-1], xs[1:]):
    arrow(ax, (x0 + 0.18, 0.54), (x1, 0.54), color="tab:purple")
ax.text(0.5, 1.16, "산란행렬법 (SMM)", ha="center", fontsize=13, weight="bold", **KFONT)
ax.text(0.5, 0.20, "계면(빨강)과 층(초록)이\n독립된 행렬로 분리됨",
        ha="center", fontsize=9, color="#555", **KFONT)

for ax in axes:
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.2, 1.25)
    ax.axis("off")

fig2.tight_layout()
fig2.savefig(os.path.join(OUT_DIR, "fig2-methods-comparison.png"), dpi=150)
plt.close(fig2)


# ────────────────────────────────────────────────────────────
# 그림 3. TMM으로 계산한 두께별 반사율 스펙트럼 (간섭 프린지 비교)
# ────────────────────────────────────────────────────────────
n_air, n_sio2, n_si = 1.0, 1.46, 3.88 - 0.02j  # 상수 근사(비분산), 400-1000nm 대역 대략값
wavelengths = np.linspace(400, 1000, 500)
thicknesses_nm = [100, 300, 600, 1000]

fig3, ax3 = plt.subplots(figsize=(7.2, 4.8))
colors = ["tab:blue", "tab:green", "tab:orange", "tab:red"]
for d_nm, c in zip(thicknesses_nm, colors):
    R = reflectance_spectrum([n_air, n_sio2, n_si], [d_nm], wavelengths, theta0=0.0, pol="s")
    ax3.plot(wavelengths, R, color=c, lw=1.8, label=f"d = {d_nm} nm")

ax3.set_xlabel("파장 (nm)", fontsize=11, **KFONT)
ax3.set_ylabel("반사율 R", fontsize=11, **KFONT)
ax3.set_title("SiO$_2$/Si 박막의 TMM 반사율 스펙트럼 (수직입사)", fontsize=12, **KFONT)
ax3.set_xlim(400, 1000)
ax3.set_ylim(0, 1)
ax3.legend(prop=LFONT, fontsize=10)
ax3.grid(alpha=0.3)

fig3.tight_layout()
fig3.savefig(os.path.join(OUT_DIR, "fig3-thickness-fringes.png"), dpi=150)
plt.close(fig3)


# ────────────────────────────────────────────────────────────
# 그림 4. SMM 텐서로 계산한 입사각×파장 Psi, Delta 맵
# ────────────────────────────────────────────────────────────
d_fixed = 300.0  # nm, 그림3의 d=300nm 곡선과 동일 구조
angles_deg = np.linspace(0.1, 85, 300)
wavelengths_2d = np.linspace(400, 1000, 300)

r_s = smm_reflectance_tensor([n_air, n_sio2, n_si], [d_fixed],
                              np.radians(angles_deg), wavelengths_2d, pol="s")
r_p = smm_reflectance_tensor([n_air, n_sio2, n_si], [d_fixed],
                              np.radians(angles_deg), wavelengths_2d, pol="p")
rho = r_p / r_s
psi_deg = np.degrees(np.arctan(np.abs(rho)))
delta_deg = np.degrees(np.angle(rho)) % 360

fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(11.5, 4.8), sharey=True)

im_a = ax4a.pcolormesh(wavelengths_2d, angles_deg, psi_deg, cmap="viridis", shading="auto")
ax4a.set_title(r"$\Psi$ (deg)", fontsize=12, **KFONT)
fig4.colorbar(im_a, ax=ax4a)

im_b = ax4b.pcolormesh(wavelengths_2d, angles_deg, delta_deg, cmap="twilight", shading="auto")
ax4b.set_title(r"$\Delta$ (deg)", fontsize=12, **KFONT)
fig4.colorbar(im_b, ax=ax4b)

for ax in (ax4a, ax4b):
    ax.set_xlabel("파장 (nm)", fontsize=11, **KFONT)
ax4a.set_ylabel("입사각 (deg)", fontsize=11, **KFONT)
fig4.suptitle(f"SiO$_2$({d_fixed:.0f}nm)/Si — 입사각×파장에 대한 " + r"$\Psi, \Delta$" + " (SMM)",
              fontsize=12, **KFONT)

fig4.tight_layout()
fig4.savefig(os.path.join(OUT_DIR, "fig4-angle-wavelength-map.png"), dpi=150)
plt.close(fig4)


# ────────────────────────────────────────────────────────────
# 그림 5. 단일 입사각으로는 구별되지 않는 두 모델이 각도축에서 갈리는 모습
#   좁은 파장 대역(630-670nm)에서 입사각 65도 기준으로 축퇴된 (d, n) 두 조합을
#   (a) 단일 입사각 스펙트럼과 (b) 입사각 스캔으로 각각 비교한다.
# ────────────────────────────────────────────────────────────
MODEL_A = (300.00, 1.460)   # 기준
MODEL_B = (311.12, 1.428)   # 65도 좁은 대역에서 A와 거의 구별되지 않는 조합


def psi_delta(d, n_film, angles_deg, wavelengths_nm):
    """Psi, Delta (deg). angles_deg, wavelengths_nm 은 1D 배열."""
    th = np.radians(np.atleast_1d(angles_deg))
    wl_arr = np.atleast_1d(wavelengths_nm)
    r_s = smm_reflectance_tensor([n_air, n_film, n_si], [d], th, wl_arr, pol="s")
    r_p = smm_reflectance_tensor([n_air, n_film, n_si], [d], th, wl_arr, pol="p")
    rho = r_p / r_s
    return np.degrees(np.arctan(np.abs(rho))), np.degrees(np.angle(rho))


wl_narrow = np.linspace(630, 670, 200)   # 좁은 대역
ang_fixed = 65.0                          # 통상적인 단일 입사각
wl_fixed = 650.0                          # 각도 스캔에 쓸 단일 파장
ang_scan = np.linspace(40, 80, 300)

pA_wl, dA_wl = psi_delta(*MODEL_A, ang_fixed, wl_narrow)
pB_wl, dB_wl = psi_delta(*MODEL_B, ang_fixed, wl_narrow)
pA_ang, dA_ang = psi_delta(*MODEL_A, ang_scan, wl_fixed)
pB_ang, dB_ang = psi_delta(*MODEL_B, ang_scan, wl_fixed)

fig5, axes5 = plt.subplots(2, 2, figsize=(11.5, 6.4))
cA, cB = "tab:blue", "tab:red"
labA = f"모델 A: d={MODEL_A[0]:.0f}nm, n={MODEL_A[1]:.3f}"
labB = f"모델 B: d={MODEL_B[0]:.0f}nm, n={MODEL_B[1]:.3f}"

# (좌) 단일 입사각 65도에서 파장 스캔 — 두 모델이 겹친다
for ax, (yA, yB), name in [(axes5[0, 0], (pA_wl.ravel(), pB_wl.ravel()), r"$\Psi$"),
                            (axes5[1, 0], (dA_wl.ravel(), dB_wl.ravel()), r"$\Delta$")]:
    ax.plot(wl_narrow, yA, color=cA, lw=2.2, label=labA)
    ax.plot(wl_narrow, yB, color=cB, lw=1.4, ls="--", label=labB)
    ax.set_ylabel(name + " (deg)", fontsize=11)
    gap = np.abs(yA - yB).max()
    ax.text(0.03, 0.10, f"최대 차이 {gap:.2f}°", transform=ax.transAxes,
            fontsize=10, color="#333", **KFONT)
axes5[1, 0].set_xlabel("파장 (nm)", fontsize=11, **KFONT)
axes5[0, 0].set_title(f"(a) 단일 입사각 {ang_fixed:.0f}° · 파장 스캔", fontsize=12, **KFONT)

# (우) 단일 파장에서 입사각 스캔 — 두 모델이 갈린다
for ax, (yA, yB), name in [(axes5[0, 1], (pA_ang.ravel(), pB_ang.ravel()), r"$\Psi$"),
                            (axes5[1, 1], (dA_ang.ravel(), dB_ang.ravel()), r"$\Delta$")]:
    ax.plot(ang_scan, yA, color=cA, lw=2.2, label=labA)
    ax.plot(ang_scan, yB, color=cB, lw=1.4, ls="--", label=labB)
    ax.set_ylabel(name + " (deg)", fontsize=11)
    diff = np.abs(yA - yB)
    i = np.argmax(diff)
    ax.axvline(ang_scan[i], color="gray", lw=1, ls=":")
    ax.text(0.03, 0.10, f"최대 차이 {diff[i]:.2f}° @ {ang_scan[i]:.1f}°",
            transform=ax.transAxes, fontsize=10, color="#333", **KFONT)
axes5[1, 1].set_xlabel("입사각 (deg)", fontsize=11, **KFONT)
axes5[0, 1].set_title(f"(b) 단일 파장 {wl_fixed:.0f}nm · 입사각 스캔", fontsize=12, **KFONT)

for ax in axes5.ravel():
    ax.grid(alpha=0.3)
axes5[0, 0].legend(prop=LFONT, fontsize=9, loc="upper right")

fig5.suptitle("좁은 파장 대역에서 서로 구별되지 않는 두 (d, n) 조합", fontsize=12, **KFONT)
fig5.tight_layout()
fig5.savefig(os.path.join(OUT_DIR, "fig5-single-vs-angle-resolved.png"), dpi=150)
plt.close(fig5)

print("저장 완료:", OUT_DIR)

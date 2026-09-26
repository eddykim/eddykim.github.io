"""전자현미경 배경이론 6편 그림 3개 생성 (한국어판·영문판).

실행: python generate_figures.py
출력: ../../assets/img/posts/electron-microscopy-tem-phase-contrast-ctf/      (한국어)
      ../../assets/img/posts/electron-microscopy-tem-phase-contrast-ctf/en/   (영문)

길이는 nm, 공간주파수는 nm^-1 로 계산한다. 그림3 은 위상물체 근사로
격자상을 실제로 합성해, 디포커스만 바꿔 대비가 뒤집히는 것을 보인다.
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
    "assets", "img", "posts", "electron-microscopy-tem-phase-contrast-ctf",
)

ACCENT = "#c0392b"
BLUE = "#2c6fbb"
GREEN = "#27795b"
ORANGE = "#d68910"
GRAY = "#7f8c8d"
DARK = "#2c3e50"

# ── 물리 상수와 장비 조건 ─────────────────────────────────
H, M0, QE, C = 6.62607015e-34, 9.1093837015e-31, 1.602176634e-19, 2.99792458e8
CS = 5.0e5           # nm, 구면수차 계수 0.5 mm — 2편에서 쓴 값과 같다
FOCUS_SPREAD = 5.0   # nm, 초점 퍼짐 (색수차·전원 안정도에서 온다)
CONVERGENCE = 5.0e-4  # rad, 빔 수렴 반각 0.5 mrad
LATTICE = 0.25       # nm, 그림3 의 격자 간격


def wavelength(volts=200e3):
    """상대론 보정을 넣은 전자 파장 [nm]. 1편과 같은 식이다."""
    v = float(volts)
    p2 = 2.0 * M0 * QE * v * (1.0 + QE * v / (2.0 * M0 * C**2))
    return H / np.sqrt(p2) * 1e9


LAM = wavelength()


def chi(k, df, cs=CS, lam=LAM):
    """수차 함수. 디포커스가 k^2, 구면수차가 k^4 로 들어간다.

        chi(k) = pi df lam k^2 + (1/2) pi Cs lam^3 k^4

    df 는 관례대로 언더포커스에서 음수다.
    """
    k = np.asarray(k, dtype=float)
    return np.pi * df * lam * k**2 + 0.5 * np.pi * cs * lam**3 * k**4


def scherzer_defocus(cs=CS, lam=LAM):
    """확장 Scherzer 디포커스 -1.2 sqrt(Cs lambda) [nm]."""
    return -1.2 * np.sqrt(cs * lam)


def scherzer_resolution(cs=CS, lam=LAM):
    """Scherzer 점분해능 0.66 (Cs lambda^3)^(1/4) [nm]. 2편과 같은 식이다."""
    return 0.66 * (cs * lam**3) ** 0.25


def envelope_temporal(k, spread=FOCUS_SPREAD, lam=LAM):
    """시간 코히런스 감쇠 포락선. 초점이 퍼진 만큼 고주파가 죽는다."""
    k = np.asarray(k, dtype=float)
    return np.exp(-0.5 * (np.pi * lam * spread) ** 2 * k**4)


def envelope_spatial(k, df, alpha=CONVERGENCE, cs=CS, lam=LAM):
    """공간 코히런스 감쇠 포락선. 빔이 완전한 평면파가 아니어서 생긴다."""
    k = np.asarray(k, dtype=float)
    return np.exp(-((np.pi * alpha) ** 2) * (cs * lam**2 * k**3 + df * k) ** 2)


def ctf(k, df, damped=True):
    """대비 전달 함수 sin(chi). 약위상물체의 대비는 -2 sin(chi) 에 비례한다."""
    t = np.sin(chi(k, df))
    if damped:
        t = t * envelope_temporal(k) * envelope_spatial(k, df)
    return t


def first_zero(df, kmax=20.0):
    """Scherzer 통과대역이 끝나는 첫 영점 [nm^-1]. 원점 근처의 영점은 건너뛴다."""
    k = np.linspace(0.05, kmax, 200000)
    v = np.sin(chi(k, df))
    sign_change = np.where(np.diff(np.sign(v)))[0]
    return float(k[sign_change[0]]) if len(sign_change) else np.nan


def information_limit(threshold=np.exp(-2.0), kmax=20.0):
    """포락선이 1/e^2 로 떨어지는 공간주파수 [nm^-1]."""
    k = np.linspace(0.05, kmax, 200000)
    e = envelope_temporal(k)
    below = np.where(e < threshold)[0]
    return float(k[below[0]]) if len(below) else np.nan


# ── 그림3: 위상물체 격자상 합성 ───────────────────────────
def simulate_lattice_image(df, n=320, field=4.0, sigma_v=0.30, atom_sigma=0.045):
    """정사각 격자의 위상물체를 만들고 CTF 를 통과시켜 상을 합성한다.

    약위상물체 근사를 쓰지 않고 exit wave 를 그대로 전파시킨다.
      psi_exit = exp(i sigma V),  Psi' = FT(psi) exp(-i chi) E,  I = |IFT(Psi')|^2
    """
    x = np.linspace(-field / 2, field / 2, n, endpoint=False)
    xx, yy = np.meshgrid(x, x)

    # 격자점마다 가우시안 투영 퍼텐셜을 올린다.
    v = np.zeros_like(xx)
    m = int(np.ceil(field / LATTICE)) + 1
    for i in range(-m, m + 1):
        for j in range(-m, m + 1):
            cx, cy = i * LATTICE, j * LATTICE
            if abs(cx) > field or abs(cy) > field:
                continue
            v += np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * atom_sigma**2))
    v = v / v.max() * sigma_v

    psi = np.exp(1j * v)
    kx = np.fft.fftfreq(n, d=field / n)
    kxx, kyy = np.meshgrid(kx, kx)
    kr = np.sqrt(kxx**2 + kyy**2)

    # numpy 의 정방향 FFT 는 exp(-2 pi i k x) 를 쓰는데, TEM 교재의 광학 관례는
    # 부호가 반대다. 그래서 전달 함수에 exp(+i chi) 를 곱해야 약위상물체의
    # 예측(Scherzer 디포커스에서 원자가 어둡게 나온다)과 맞는다.
    transfer = np.exp(1j * chi(kr, df)) * envelope_temporal(kr) \
        * envelope_spatial(kr, df)
    img = np.abs(np.fft.ifft2(np.fft.fft2(psi) * transfer)) ** 2
    return x, v, img


LABELS = {
    "ko": {
        "dir": BASE_DIR, "font": KFONT, "legend": LFONT,
        "fig1": "디포커스가 통과시킬 공간주파수를 고른다",
        "f1x": "공간주파수 k (nm$^{-1}$)", "f1y": "$\\sin χ(k)$",
        "f1sch": "Scherzer 디포커스 ({:.1f} nm)", "f1u": "{:.0f} nm (언더포커스)",
        "f1o": "{:.0f} nm (오버포커스)",
        "f1band": "Scherzer 통과대역",
        "f1zero": "첫 영점 {:.2f} nm$^{{-1}}$\n= 점분해능 {:.3f} nm",
        "f1note": "$\\sin χ$ 가 음수인 구간과 양수인 구간에서\n같은 구조가 반대 명암으로 나온다",
        "fig2": "감쇠 포락선이 정보한계를 정한다",
        "f2x": "공간주파수 k (nm$^{-1}$)", "f2y": "전달 함수",
        "f2ctf": "$\\sin χ$ (감쇠 없음)", "f2et": "시간 코히런스 포락선",
        "f2es": "공간 코히런스 포락선", "f2prod": "실제 전달 함수",
        "f2point": "점분해능\n{:.3f} nm", "f2info": "정보한계\n{:.3f} nm",
        "f2gap": "이 사이의 정보는 전달되지만\n부호가 뒤집혀 그대로 읽을 수 없다",
        "fig3": "같은 시편, 디포커스만 바꾼 상",
        "f3obj": "실제 원자 위치\n(투영 퍼텐셜)",
        "f3panel": "Δf = {:.1f} nm\n$\\sin χ$({:.0f} nm$^{{-1}}$) = {:+.2f}",
        "f3dark": "원자가 어둡다", "f3none": "격자가 거의 사라진다\n(진폭이 30분의 1)",
        "f3bright": "원자가 밝다 — 대비가 뒤집혔다",
    },
    "en": {
        "dir": os.path.join(BASE_DIR, "en"), "font": EFONT, "legend": ELFONT,
        "fig1": "Defocus selects which spatial frequencies get through",
        "f1x": "Spatial frequency k (nm$^{-1}$)", "f1y": "$\\sin χ(k)$",
        "f1sch": "Scherzer defocus ({:.1f} nm)", "f1u": "{:.0f} nm (underfocus)",
        "f1o": "{:.0f} nm (overfocus)",
        "f1band": "Scherzer passband",
        "f1zero": "First zero {:.2f} nm$^{{-1}}$\n= point resolution {:.3f} nm",
        "f1note": "Where $\\sin χ$ is negative and where it is positive,\nthe same structure appears with opposite contrast",
        "fig2": "Damping envelopes set the information limit",
        "f2x": "Spatial frequency k (nm$^{-1}$)", "f2y": "Transfer function",
        "f2ctf": "$\\sin χ$ (undamped)", "f2et": "Temporal coherence envelope",
        "f2es": "Spatial coherence envelope", "f2prod": "Effective transfer function",
        "f2point": "Point resolution\n{:.3f} nm", "f2info": "Information limit\n{:.3f} nm",
        "f2gap": "Information here is transmitted but\nsign-reversed, so it cannot be read directly",
        "fig3": "Same specimen, defocus alone changed",
        "f3obj": "True atom positions\n(projected potential)",
        "f3panel": "Δf = {:.1f} nm\n$\\sin χ$({:.0f} nm$^{{-1}}$) = {:+.2f}",
        "f3dark": "Atoms appear dark", "f3none": "The lattice nearly vanishes\n(amplitude down 30x)",
        "f3bright": "Atoms appear bright — contrast inverted",
    },
}


def fig1_ctf_curves(L):
    # k > 8 nm^-1 에서는 진동이 너무 빨라 선이 뭉개지기만 한다.
    k = np.linspace(0, 8, 4000)
    df_s = scherzer_defocus()
    fig, ax = plt.subplots(figsize=(10.5, 5.4))

    for df, color, key, ls in [(df_s, ACCENT, "f1sch", "-"),
                               (-20.0, BLUE, "f1u", "--"),
                               (20.0, GREEN, "f1o", ":")]:
        ax.plot(k, np.sin(chi(k, df)), color=color, lw=2.0, ls=ls,
                label=L[key].format(df))

    kz = first_zero(df_s)
    ax.axvspan(0, kz, color=ACCENT, alpha=0.10, zorder=0)
    ax.text(kz / 2, 1.72, L["f1band"], ha="center", fontsize=9.5, color=ACCENT,
            **L["font"])
    ax.plot([kz], [0], "o", color=ACCENT, ms=8, zorder=5)
    # 주석은 곡선이 없는 위쪽 여백에 둔다.
    ax.annotate(L["f1zero"].format(kz, 1.0 / kz), xy=(kz, 0),
                xytext=(kz + 0.25, 1.38), ha="left", fontsize=9, color=ACCENT,
                **L["font"], arrowprops=dict(arrowstyle="->", color=ACCENT))

    ax.axhline(0, color=DARK, lw=0.9)
    ax.set_xlabel(L["f1x"], **L["font"])
    ax.set_ylabel(L["f1y"], **L["font"])
    ax.set_ylim(-1.15, 1.95)
    ax.set_xlim(0, 8)
    ax.grid(alpha=0.3)
    ax.set_title(L["fig1"], fontsize=13, **L["font"])
    # 범례는 축 아래로 빼서 곡선을 가리지 않게 한다.
    ax.legend(prop=L["legend"] or None, fontsize=9, ncol=3,
              loc="upper center", bbox_to_anchor=(0.5, -0.13), frameon=False)
    fig.tight_layout()
    return fig


def fig2_envelopes(L):
    k = np.linspace(0, 12, 4000)
    df_s = scherzer_defocus()
    fig, ax = plt.subplots(figsize=(10.5, 5.0))

    ax.plot(k, np.sin(chi(k, df_s)), color=GRAY, lw=1.2, alpha=0.75, label=L["f2ctf"])
    ax.plot(k, envelope_temporal(k), color=BLUE, lw=1.8, ls="--", label=L["f2et"])
    ax.plot(k, envelope_spatial(k, df_s), color=GREEN, lw=1.8, ls=":", label=L["f2es"])
    ax.plot(k, ctf(k, df_s), color=ACCENT, lw=2.3, label=L["f2prod"])

    kz, ki = first_zero(df_s), information_limit()
    for kv, key, col in [(kz, "f2point", DARK), (ki, "f2info", ORANGE)]:
        ax.axvline(kv, color=col, ls="-.", lw=1.3)
        # sin(chi) 진동 위에 얹히므로 흰 배경을 깔아 읽히게 한다.
        ax.text(kv + 0.12, 0.88, L[key].format(1.0 / kv), fontsize=9, color=col,
                **L["font"],
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5))
    ax.axvspan(kz, ki, color=ORANGE, alpha=0.10, zorder=0)
    ax.text((kz + ki) / 2, 1.28, L["f2gap"], ha="center", fontsize=9,
            color="#9c6500", **L["font"])

    ax.axhline(0, color=DARK, lw=0.9)
    ax.set_xlabel(L["f2x"], **L["font"])
    ax.set_ylabel(L["f2y"], **L["font"])
    ax.set_ylim(-1.05, 1.62)
    ax.set_xlim(0, 12)
    ax.grid(alpha=0.3)
    ax.legend(prop=L["legend"] or None, loc="lower right", fontsize=8.5)
    ax.set_title(L["fig2"], fontsize=13, **L["font"])
    fig.tight_layout()
    return fig


def fig3_lattice_images(L):
    k_lat = 1.0 / LATTICE
    cases = [(scherzer_defocus(), "f3dark"), (-50.1, "f3none"), (-62.5, "f3bright")]
    fig, axes = plt.subplots(1, 4, figsize=(14, 4.2))

    x, v, _ = simulate_lattice_image(cases[0][0])
    ext = [x[0], x[-1], x[0], x[-1]]
    axes[0].imshow(v, cmap="magma", extent=ext, origin="lower")
    axes[0].set_title(L["f3obj"], fontsize=9.5, **L["font"])

    # 세 상을 같은 명암 범위로 그린다. 패널마다 따로 정규화하면 대비가 거의
    # 사라진 상까지 끝까지 늘어나, 실제로는 30 분의 1 인 진폭이 똑같아 보인다.
    imgs = [simulate_lattice_image(df)[2] for df, _ in cases]
    lo = min(im.min() for im in imgs)
    hi = max(im.max() for im in imgs)

    for ax, (df, key), img in zip(axes[1:], cases, imgs):
        ax.imshow(img, cmap="gray", extent=ext, origin="lower", vmin=lo, vmax=hi)
        ax.set_title(L["f3panel"].format(df, k_lat, np.sin(chi(k_lat, df))),
                     fontsize=9.5, **L["font"])
        ax.set_xlabel(L[key], fontsize=9, **L["font"])

    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([])

    fig.suptitle(L["fig3"], fontsize=13, **L["font"])
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def main():
    for lang, L in LABELS.items():
        os.makedirs(L["dir"], exist_ok=True)
        for name, builder in [("fig1-ctf-curves", fig1_ctf_curves),
                              ("fig2-envelopes", fig2_envelopes),
                              ("fig3-lattice-images", fig3_lattice_images)]:
            fig = builder(L)
            fig.savefig(os.path.join(L["dir"], name + ".png"), dpi=150,
                        bbox_inches="tight")
            plt.close(fig)
        print(f"[{lang}] 3 figures written to {L['dir']}")

    df_s = scherzer_defocus()
    kz, ki = first_zero(df_s), information_limit()
    print()
    print(f"lambda = {LAM*1000:.3f} pm, Cs = {CS/1e6:.1f} mm")
    print(f"Scherzer 디포커스 = {df_s:.2f} nm")
    print(f"첫 영점 k = {kz:.3f} nm^-1 -> 점분해능 {1/kz:.4f} nm")
    print(f"Scherzer 점분해능 공식값       = {scherzer_resolution():.4f} nm")
    print(f"정보한계 k = {ki:.3f} nm^-1 -> {1/ki:.4f} nm")
    print(f"\n격자 {LATTICE} nm (k = {1/LATTICE:.0f} nm^-1) 에서 sin(chi):")
    for df in (df_s, -50.1, -62.5):
        print(f"  df = {df:7.2f} nm : {np.sin(chi(1/LATTICE, df)):+.3f}")


if __name__ == "__main__":
    main()

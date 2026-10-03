"""계산 이미징과 타이코그래피 2편 그림 8개 생성 (한국어판·영문판)."""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from skimage import data, transform

from phase_retrieval import (F, iF, P_m, P_s, make_object, random_start, step, estimate,
                             align_error, twin, _gray, _resize)
import studies

SLUG = "ptycho-phase-problem-iterative-retrieval"
BASE = os.path.join(os.path.dirname(__file__), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
C_TRUE, C_A, C_B, C_C = "#555555", "#c8553d", "#2f6690", "#6a994e"

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "camA": "영상 A", "camB": "영상 B", "swap1": "A의 크기 + B의 위상", "swap2": "B의 크기 + A의 위상",
        "fig1": "그림1. 푸리에 크기와 위상을 맞바꾸면, 상은 위상을 준 쪽을 닮는다",
        "orig": "원래 물체", "shift": "평행이동", "twin": "쌍둥이 상 (공액 반전)",
        "inten": "세 경우의 회절 세기 (로그)",
        "fig2": "그림2. 서로 다른 세 물체가 똑같은 회절 세기를 만든다",
        "idx": "위치 $n$", "lag": "어긋남 $k$", "sig_t": "지지 영역이 같은 두 신호", "ac_t": "자기상관 (= 세기의 역변환)",
        "fig3": "그림3. 1차원에서는 자명하지 않은 두 신호가 같은 세기를 낸다",
        "box_t": "물체와 그 자기상관", "obj_l": "물체 지지", "ac_l": "자기상관 지지",
        "sig_ax": "오버샘플링 비 $\\sigma$ (상자 면적 / 지지 면적)", "succ": "성공 횟수 (8회 중)",
        "real600": "실수 물체, 600회", "real3000": "실수 물체, 3000회",
        "cplx600": "복소 물체, 600회", "cplx3000": "복소 물체, 3000회",
        "fig4": "그림4. 자기상관은 물체의 두 배 크기다. σ가 2 근처면 해가 있어도 찾기 어렵다",
        "est": "실공간 추정 $g$", "fourier": "푸리에 공간 $G$", "mod": "크기만 측정값으로\n(위상은 유지)",
        "supp": "지지 밖을 0으로\n(실수면 음수도 0)", "fft": "FFT", "ifft": "역 FFT",
        "fig5": "그림5. Gerchberg–Saxton과 ER: 두 공간의 제약을 번갈아 강제한다",
        "iters": "반복 횟수", "err": "상대 오차", "truth": "참 물체",
        "er600": "ER 600회", "hio600": "HIO 600회",
        "fig6": "그림6. ER은 금방 0.19 근처에서 멈추고, HIO는 한동안 오르내리다 정답(여기서는 쌍둥이 상)에 닿는다",
        "line": "지지 집합 S (볼록)", "curve": "모듈러스 집합 M (비볼록)", "start": "시작",
        "er": "ER", "hio": "HIO", "sol": "교점 (정답)", "trap": "ER이 멈춘 곳",
        "fig7": "그림7. ER은 두 집합이 가장 가까운 곳에 갇히고, HIO는 틈의 방향으로 밀려나 교점을 찾는다",
        "ph_ax": "총 광자 수", "med": "오차 중앙값 (6회)", "bestt": "가장 좋은 시도",
        "amp": "진폭", "phs": "위상", "phot": "광자",
        "fig8": "그림8. 패턴 한 장에는 여유분이 없어서, 광자가 줄면 복원이 빠르게 무너진다",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "camA": "Image A", "camB": "Image B", "swap1": "|A| with phase of B", "swap2": "|B| with phase of A",
        "fig1": "Fig 1. Swap Fourier magnitude and phase, and the image follows the phase",
        "orig": "Original object", "shift": "Translated", "twin": "Twin (conjugate inversion)",
        "inten": "Diffraction intensity of all three (log)",
        "fig2": "Fig 2. Three different objects produce exactly the same diffraction intensity",
        "idx": "Position $n$", "lag": "Lag $k$", "sig_t": "Two signals with the same support",
        "ac_t": "Autocorrelation (= inverse transform of intensity)",
        "fig3": "Fig 3. In 1D, two non-trivially different signals give the same intensity",
        "box_t": "Object and its autocorrelation", "obj_l": "Object support", "ac_l": "Autocorrelation support",
        "sig_ax": "Oversampling ratio $\\sigma$ (box area / support area)", "succ": "Successes (out of 8)",
        "real600": "Real object, 600 it.", "real3000": "Real object, 3000 it.",
        "cplx600": "Complex object, 600 it.", "cplx3000": "Complex object, 3000 it.",
        "fig4": "Fig 4. The autocorrelation is twice the object's size; near σ = 2 a solution exists but is hard to find",
        "est": "Real-space estimate $g$", "fourier": "Fourier space $G$", "mod": "Replace magnitude\n(keep phase)",
        "supp": "Zero outside support\n(and negatives, if real)", "fft": "FFT", "ifft": "Inverse FFT",
        "fig5": "Fig 5. Gerchberg–Saxton and ER: enforce the constraints of each space in turn",
        "iters": "Iterations", "err": "Relative error", "truth": "True object",
        "er600": "ER, 600 it.", "hio600": "HIO, 600 it.",
        "fig6": "Fig 6. ER stalls near 0.19 almost at once; HIO wanders before reaching the answer (here, the twin)",
        "line": "Support set S (convex)", "curve": "Modulus set M (non-convex)", "start": "Start",
        "er": "ER", "hio": "HIO", "sol": "Intersection (answer)", "trap": "Where ER stops",
        "fig7": "Fig 7. ER is trapped where the sets come closest; HIO is pushed along the gap and finds the intersection",
        "ph_ax": "Total photons", "med": "Median error (6 trials)", "bestt": "Best trial",
        "amp": "Amplitude", "phs": "Phase", "phot": "photons",
        "fig8": "Fig 8. A single pattern has no redundancy, so the reconstruction collapses quickly as photons drop",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
N = 128

# 그림1: 위상 맞바꾸기
IA = _resize(data.camera().astype(float), 256)
IB = _resize(_gray(data.astronaut()), 256)
FA, FB = np.fft.fft2(IA), np.fft.fft2(IB)
SW1 = np.real(np.fft.ifft2(np.abs(FA) * np.exp(1j * np.angle(FB))))
SW2 = np.real(np.fft.ifft2(np.abs(FB) * np.exp(1j * np.angle(FA))))

# 그림2: 자명한 모호성 (실수 물체라 공액 반전 = 뒤집기)
OBJ_R, SUP_R = make_object(N, 48, "real")
SHIFTED = np.roll(OBJ_R, (14, -20), axis=(0, 1))
TWIN = twin(OBJ_R)
INT_LOG = np.log10(np.fft.fftshift(np.abs(F(OBJ_R)) ** 2) + 1e-6)

# 그림3: 1D 반례
U = np.array([1, 0, -2, 0, -2.])
V = np.array([1 - np.sqrt(3), 0, 1, 0, 1 + np.sqrt(3)])

# 그림4: 자기상관 지지와 오버샘플링
AC = np.abs(np.fft.fftshift(iF(np.abs(F(OBJ_R)) ** 2)))
OVS = studies.oversampling()

# 그림6: ER 과 HIO 의 오차 궤적 (같은 시작점)
MAG_R = np.abs(F(OBJ_R))
X0 = random_start(MAG_R, np.random.default_rng(3))
TRAJ = {}
RECS = {}
for alg in ("ER", "HIO"):
    x = X0.copy(); errs = []
    for k in range(600):
        x = step(alg, x, MAG_R, SUP_R, real=True)
        errs.append(align_error(estimate(alg, x, MAG_R, SUP_R, real=True), OBJ_R))
    TRAJ[alg] = np.array(errs); RECS[alg] = estimate(alg, x, MAG_R, SUP_R, real=True)

# 그림7: 집합 장난감
XS = np.linspace(-4, 12, 200001)
CURVE = np.stack([XS, 0.6 - 0.5 * np.cos(XS) - 0.15 * XS], 1)
X_CROSS = XS[np.where(np.diff(np.sign(CURVE[:, 1])))[0][0]]


def toy(alg, n=600, beta=0.9):
    def pm(p):
        return CURVE[np.argmin(((CURVE - p) ** 2).sum(1))]
    x = np.array([0.8, 0.0]); path = [x.copy()]
    for _ in range(n):
        q = pm(x)
        x = np.array([q[0], 0.0]) if alg == "ER" else np.array([q[0], x[1] - beta * q[1]])
        path.append(x.copy())
    return np.array(path)


TOY = {a: toy(a) for a in ("ER", "HIO")}

# 그림8: 노이즈
NOISE = studies.noise()
OBJ_C, _ = make_object(N, 48, "complex")

# 본문 표용
ALG = studies.algorithm_table()
RBETA = studies.raar_beta()

SUMMARY = {"traj_ER": {k: TRAJ["ER"][k - 1] for k in (1, 10, 50, 100, 200, 600)},
           "traj_HIO": {k: TRAJ["HIO"][k - 1] for k in (1, 10, 50, 100, 200, 600)},
           "toy_cross": X_CROSS, "toy_ER_end": TOY["ER"][-1], "toy_HIO_end": TOY["HIO"][-1],
           "oversampling": OVS, "noise": {k: (v[0], v[2]) for k, v in NOISE.items()},
           "alg": ALG, "raar_beta": RBETA}


def crop(a, c=N // 2, h=40):
    return a[c - h:c + h, c - h:c + h]


# ── 그리기 ───────────────────────────────────────────────────
def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.06, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **Fk)
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    def img(ax, a, title, cmap="gray", **kw):
        ax.imshow(a, cmap=cmap, **kw); ax.set_title(title, fontsize=9.5, **Fk); ax.set_axis_off()

    # 그림1
    fig, ax = plt.subplots(1, 4, figsize=(12, 3.5))
    for a, im, k in zip(ax, (IA, IB, SW1, SW2), ("camA", "camB", "swap1", "swap2")):
        img(a, im, L[k])
    cap(fig, "fig1", rect=(0, 0.07, 1, 1)); save(fig, "fig1-phase-swap.png")

    # 그림2
    fig, ax = plt.subplots(1, 4, figsize=(12, 3.5))
    for a, im, k in zip(ax[:3], (OBJ_R, SHIFTED, TWIN), ("orig", "shift", "twin")):
        img(a, np.abs(im), L[k], vmin=0, vmax=1)
    img(ax[3], INT_LOG, L["inten"], cmap="magma")
    cap(fig, "fig2", rect=(0, 0.07, 1, 1)); save(fig, "fig2-trivial-ambiguities.png")

    # 그림3
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.4))
    n = np.arange(5)
    ax[0].stem(n - 0.08, U, linefmt=C_A, markerfmt="o", basefmt=" ", label="u")
    ax[0].stem(n + 0.08, V, linefmt=C_B, markerfmt="s", basefmt=" ", label="v")
    ax[0].axhline(0, color="#999", lw=0.8); ax[0].set_xticks(n); ax[0].set_xlabel(L["idx"], **Fk)
    ax[0].set_title(L["sig_t"], fontsize=10, **Fk); ax[0].legend()
    k = np.arange(-4, 5)
    ax[1].plot(k, np.correlate(U, U, "full"), "o-", color=C_A, label="u ⋆ u")
    ax[1].plot(k, np.correlate(V, V, "full"), "s--", color=C_B, mfc="none", ms=9, label="v ⋆ v")
    ax[1].set_xlabel(L["lag"], **Fk); ax[1].set_title(L["ac_t"], fontsize=10, **Fk); ax[1].legend()
    for a in ax: a.grid(alpha=0.3)
    cap(fig, "fig3", rect=(0, 0.07, 1, 1)); save(fig, "fig3-1d-counterexample.png")

    # 그림4
    fig, ax = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1, 1.5]})
    ax[0].imshow(np.log10(AC + 1e-6), cmap="magma", vmin=-2)
    c, s = N // 2, 48
    ax[0].add_patch(plt.Rectangle((c - s / 2, c - s / 2), s, s, fill=False, ec="#7fd3ff", ls="--", lw=1.4,
                                  label=L["obj_l"]))
    ax[0].add_patch(plt.Rectangle((c - s, c - s), 2 * s, 2 * s, fill=False, ec=C_C, lw=1.5, label=L["ac_l"]))
    ax[0].set_axis_off(); ax[0].set_title(L["box_t"], fontsize=10, **Fk)
    ax[0].legend(prop=LG, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False)
    sizes = sorted({k[1] for k in OVS}, reverse=True)
    sig = [N * N / s_ ** 2 for s_ in sizes]
    for kind, it, key, c_, ls, mk, ms in (("real", 600, "real600", C_A, "--", "o", 9), ("real", 3000, "real3000", C_A, "-", "o", 9),
                                          ("complex", 600, "cplx600", C_B, "--", "s", 5), ("complex", 3000, "cplx3000", C_B, "-", "s", 5)):
        ax[1].plot(sig, [OVS[(kind, s_, it)] for s_ in sizes], ls, marker=mk, ms=ms, color=c_,
                   mfc="none" if it == 600 else c_, label=L[key])
    ax[1].axvline(2, color="#999", ls=":", lw=1); ax[1].set_ylim(-0.5, 8.5)
    ax[1].set_xlabel(L["sig_ax"], **Fk); ax[1].set_ylabel(L["succ"], **Fk)
    ax[1].legend(prop=LG, fontsize=8, loc="lower right"); ax[1].grid(alpha=0.3)
    cap(fig, "fig4", rect=(0, 0.07, 1, 1)); save(fig, "fig4-oversampling.png")

    # 그림5: 순환 도식
    fig, ax = plt.subplots(figsize=(9, 4.2)); ax.set_axis_off(); ax.set_xlim(0, 10); ax.set_ylim(0, 5)

    def box(x, y, w, h, t, c):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.06", fc=c, ec="none", alpha=0.9))
        ax.text(x + w / 2, y + h / 2, t, ha="center", va="center", fontsize=9.5, color="white", **Fk)

    box(0.6, 3.2, 3.0, 1.0, L["est"], C_TRUE); box(6.4, 3.2, 3.0, 1.0, L["fourier"], C_B)
    box(6.4, 0.6, 3.0, 1.1, L["mod"], C_A); box(0.6, 0.6, 3.0, 1.1, L["supp"], C_C)
    kw = dict(arrowprops=dict(arrowstyle="->", color="#333", lw=1.5))
    ax.annotate("", (6.3, 3.7), (3.7, 3.7), **kw); ax.text(5.0, 3.9, L["fft"], ha="center", fontsize=9, **Fk)
    ax.annotate("", (7.9, 1.8), (7.9, 3.1), **kw)
    ax.annotate("", (3.7, 1.15), (6.3, 1.15), **kw); ax.text(5.0, 1.35, L["ifft"], ha="center", fontsize=9, **Fk)
    ax.annotate("", (2.1, 3.1), (2.1, 1.8), **kw)
    fig.text(0.5, 0.02, L["fig5"], ha="center", fontsize=9, **Fk)
    save(fig, "fig5-er-loop.png")

    # 그림6
    fig = plt.figure(figsize=(12, 4))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.6, 1, 1, 1])
    a0 = fig.add_subplot(gs[0])
    it = np.arange(1, 601)
    a0.semilogy(it, TRAJ["ER"], color=C_A, lw=1.4, label="ER")
    a0.semilogy(it, TRAJ["HIO"], color=C_B, lw=1.4, label="HIO")
    a0.set_xlabel(L["iters"], **Fk); a0.set_ylabel(L["err"], **Fk); a0.legend(); a0.grid(alpha=0.3, which="both")
    for j, (im, key, e) in enumerate(((OBJ_R, "truth", None), (RECS["ER"], "er600", TRAJ["ER"][-1]),
                                      (RECS["HIO"], "hio600", TRAJ["HIO"][-1]))):
        a = fig.add_subplot(gs[j + 1])
        img(a, crop(np.abs(im)), L[key] + (f" ({e:.3f})" if e is not None else ""), vmin=0, vmax=1)
    cap(fig, "fig6", rect=(0, 0.06, 1, 1)); save(fig, "fig6-er-vs-hio.png")

    # 그림7
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={"width_ratios": [1.4, 1]})
    for a in ax:
        a.axhline(0, color=C_C, lw=2.2, label=L["line"])
        a.plot(CURVE[::200, 0], CURVE[::200, 1], color=C_TRUE, lw=1.6, label=L["curve"])
        a.plot(*TOY["ER"].T, ".-", color=C_A, ms=3, lw=0.8, label=L["er"])
        a.plot(*TOY["HIO"].T, ".-", color=C_B, ms=2, lw=0.6, alpha=0.8, label=L["hio"])
        a.plot(X_CROSS, 0, "*", color="gold", ms=15, mec="k", label=L["sol"])
        a.plot(0.8, 0, "ko", ms=5)
        a.grid(alpha=0.3)
    ax[0].set_xlim(-2, 8); ax[0].set_ylim(-16, 2.5); ax[0].legend(prop=LG, fontsize=8, loc="lower left")
    ax[1].set_xlim(-0.5, 5.5); ax[1].set_ylim(-0.5, 1.2)
    ax[1].annotate(L["trap"], TOY["ER"][-1], (1.2, 0.8), color=C_A, fontsize=9,
                   arrowprops=dict(arrowstyle="->", color=C_A), **Fk)
    ax[1].annotate(L["start"], (0.8, 0), (0.2, -0.35), fontsize=9, arrowprops=dict(arrowstyle="->"), **Fk)
    cap(fig, "fig7", rect=(0, 0.06, 1, 1)); save(fig, "fig7-sets-toy.png")

    # 그림8
    fig = plt.figure(figsize=(12, 4.6))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.6, 1, 1, 1])
    a0 = fig.add_subplot(gs[:, 0])
    ph = np.array(sorted(NOISE))
    a0.semilogx(ph, [NOISE[p][0] for p in ph], "o-", color=C_B, label=L["med"])
    a0.semilogx(ph, [NOISE[p][2] for p in ph], "s--", color=C_A, mfc="none", label=L["bestt"])
    a0.invert_xaxis(); a0.set_xlabel(L["ph_ax"], **Fk); a0.set_ylabel(L["err"], **Fk)
    a0.legend(prop=LG, fontsize=8); a0.grid(alpha=0.3, which="both")
    for j, p in enumerate((1e9, 1e7, 1e5)):
        rec = NOISE[p][1]
        # 표시용으로 모호성을 맞춘 상을 다시 만든다
        best = None
        for cand in (rec, twin(rec)):
            xc = iF(F(OBJ_C) * np.conj(F(cand)))
            sh = np.unravel_index(np.argmax(np.abs(xc)), xc.shape)
            c2 = np.roll(cand, sh, axis=(0, 1)); a_ = np.vdot(c2, OBJ_C) / np.vdot(c2, c2)
            e = np.linalg.norm(a_ * c2 - OBJ_C)
            if best is None or e < best[0]: best = (e, a_ * c2)
        r = crop(best[1], h=30)
        img(fig.add_subplot(gs[0, j + 1]), np.abs(r), f"{p:.0e} {L['phot']}: {L['amp']}", vmin=0, vmax=1)
        img(fig.add_subplot(gs[1, j + 1]), np.angle(r) * (np.abs(r) > 0.2), L["phs"], cmap="twilight",
            vmin=-np.pi, vmax=np.pi)
    cap(fig, "fig8", rect=(0, 0.06, 1, 1)); save(fig, "fig8-noise.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 8개 저장: {L['dir']}")
    print("\n본문 수치")
    for k, v in SUMMARY.items():
        print(f"  {k}: {v}")

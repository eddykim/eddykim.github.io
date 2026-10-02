"""계산 이미징과 타이코그래피 1편 그림 10개 생성 (한국어판·영문판)."""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch
from skimage import data, transform

from forward_models import (K, RASKAR_CODE, box_kernel, coded_kernel, smear_matrix,
                            transfer_function, test_signal, measure)
from regularize import (least_squares, tikhonov, landweber_svd, tv_admm, difference_matrix,
                        rel_error, best_lambda, l_curve, discrepancy_lambda)
import noise_study

SLUG = "ptycho-computational-imaging-inverse-problem"
BASE = os.path.join(os.path.dirname(__file__), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
C_TRUE, C_BOX, C_CODE, C_ACC = "#555555", "#c8553d", "#2f6690", "#6a994e"

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "scene": "장면", "lens": "렌즈", "image": "상", "done": "끝",
        "optics": "부호화 광학", "meas": "측정\n(사람이 못 알아봄)", "algo": "복원 계산",
        "trad": "전통적 이미징: 광학이 상을 완성한다",
        "comp": "계산 이미징: 광학은 부호화하고, 계산이 상을 만든다",
        "fig1": "그림1. 계산 이미징에서 측정은 상이 아니라 부호다",
        "pos": "위치 (화소)", "val": "밝기", "orig": "원 신호 $x$",
        "chip": "노출 칩", "kbox": "상자 셔터 (52칩 모두 열림)", "kcode": "부호 셔터 (26칩 열림)",
        "ybox": "상자 번짐", "ycode": "부호 번짐", "meas_t": "측정 $y = Ax$ (노이즈 없음)",
        "kern_t": "번짐 커널", "sig_t": "원 신호",
        "fig2": "그림2. 같은 길이만큼 움직여도 셔터를 어떻게 여닫느냐에 따라 측정이 달라진다",
        "freq": "공간주파수 (주기/화소)", "mag": "전달함수 크기 $|H|$",
        "zeros": "영점: 이 주파수의 정보는 사라진다",
        "fig3": "그림3. 상자 셔터의 전달함수는 영점을 갖고, 부호 셔터는 빛이 절반인 대신 영점이 없다",
        "s0": "노이즈 0", "s1": "노이즈 $\\sigma = 10^{-3}$", "s2": "노이즈 $\\sigma = 10^{-2}$",
        "err": "상대 오차",
        "fig4": "그림4. 노이즈가 없으면 완벽히 되돌아오지만, 노이즈가 1 %면 되돌린 신호가 원 신호보다 커진다",
        "idx": "특잇값 번호 $i$", "sv": "특잇값 $s_i$", "coef": "측정 성분 $|u_i \\cdot y|$",
        "ratio": "복원한 해의 성분 $|u_i \\cdot y| / s_i$", "floor": "노이즈 바닥",
        "coef0": "노이즈가 없을 때의 측정 성분", "xtrue": "참 신호의 성분 $|v_i \\cdot x|$",
        "p_data": "측정: 노이즈 바닥 아래로는 못 내려간다", "p_sol": "해: 나누면 노이즈가 불어난다",
        "fig5": "그림5. 측정 성분이 노이즈 바닥에 닿는 곳부터 복원한 해가 참 신호에서 떨어져 나간다 (피카르 그림)",
        "small": "작은 $\\lambda$", "bestl": "최적 $\\lambda$", "large": "큰 $\\lambda$",
        "res": "잔차 $\\|Ax - y\\|$", "soln": "해의 크기 $\\|x\\|$",
        "disc": "불일치 원리", "noiselvl": "노이즈의 크기", "bestpt": "오차 최소", "ff": "필터 인자 $f_i$",
        "recon_t": "$\\lambda$ 에 따른 복원", "lc_t": "잔차가 노이즈만큼 남는 $\\lambda$", "ff_t": "특잇값별로 남기는 비율",
        "fig6": "그림6. 티호노프 정규화는 작은 특잇값 쪽 성분만 눌러 노이즈를 막는다",
        "iters": "반복 횟수", "stop": "여기서 멈추면 최선", "times": "회",
        "fig7": "그림7. 최소자승 경사하강은 처음엔 좋아지다가, 오래 돌리면 노이즈를 복원하기 시작한다",
        "tikd": "매끈함 가정 (미분 티호노프)", "tv": "계단 가정 (TV)",
        "fig8": "그림8. 신호가 계단형이면 기울기가 드물다는 가정이 모서리를 살린다",
        "sig_ax": "가우시안 노이즈 $\\sigma$", "ph_ax": "칩당 광자 수 (밝기 1 기준)",
        "g_t": "읽기 노이즈 (신호와 무관)", "p_t": "샷 노이즈 (광자 수의 제곱근)",
        "bt": "상자, 티호노프", "ct": "부호, 티호노프", "bv": "상자, TV", "cv": "부호, TV",
        "fig9": "그림9. 부호 셔터의 이득은 노이즈의 성격에 달렸다. 샷 노이즈에서는 빛 손해가 덜 든다",
        "i_orig": "원 영상", "i_box": "상자 번짐 측정", "i_code": "부호 번짐 측정",
        "i_bls": "상자: 최소자승", "i_cls": "부호: 최소자승",
        "i_btk": "상자: 티호노프", "i_ctk": "부호: 티호노프",
        "fig10": "그림10. 칩당 광자 1000개의 수평 모션 블러. 괄호는 원 영상 대비 상대 오차",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "scene": "Scene", "lens": "Lens", "image": "Image", "done": "Done",
        "optics": "Coding\noptics", "meas": "Measurement\n(unreadable)", "algo": "Recon-\nstruction",
        "trad": "Conventional imaging: optics finish the image",
        "comp": "Computational imaging: optics encode, computation forms the image",
        "fig1": "Fig 1. In computational imaging the measurement is a code, not an image",
        "pos": "Position (pixel)", "val": "Brightness", "orig": "Original $x$",
        "chip": "Exposure chip", "kbox": "Box shutter (all 52 chips open)", "kcode": "Coded shutter (26 chips open)",
        "ybox": "Box blur", "ycode": "Coded blur", "meas_t": "Measurement $y = Ax$ (noiseless)",
        "kern_t": "Blur kernel", "sig_t": "Original signal",
        "fig2": "Fig 2. Same motion, different measurements, depending on how the shutter flutters",
        "freq": "Spatial frequency (cycles/pixel)", "mag": "Transfer function $|H|$",
        "zeros": "Zeros: information at these frequencies is gone",
        "fig3": "Fig 3. The box shutter has spectral zeros; the coded shutter has none, at the cost of half the light",
        "s0": "No noise", "s1": "Noise $\\sigma = 10^{-3}$", "s2": "Noise $\\sigma = 10^{-2}$",
        "err": "relative error",
        "fig4": "Fig 4. Without noise the inversion is perfect; with 1 % noise the error exceeds the signal",
        "idx": "Singular value index $i$", "sv": "Singular value $s_i$", "coef": "Data $|u_i \\cdot y|$",
        "ratio": "Recovered $|u_i \\cdot y| / s_i$", "floor": "Noise floor",
        "coef0": "Noiseless data", "xtrue": "True signal $|v_i \\cdot x|$",
        "p_data": "Data cannot go below the noise floor", "p_sol": "Dividing by $s_i$ amplifies the noise",
        "fig5": "Fig 5. Where the data hit the noise floor, the recovered coefficients depart from the true signal (Picard plot)",
        "small": "Small $\\lambda$", "bestl": "Best $\\lambda$", "large": "Large $\\lambda$",
        "res": "Residual $\\|Ax - y\\|$", "soln": "Solution norm $\\|x\\|$",
        "disc": "Discrepancy principle", "noiselvl": "Noise level", "bestpt": "Minimum error", "ff": "Filter factor $f_i$",
        "recon_t": "Reconstruction vs $\\lambda$", "lc_t": "$\\lambda$ that leaves a noise-sized residual", "ff_t": "Fraction kept per singular value",
        "fig6": "Fig 6. Tikhonov regularization damps only the small-singular-value components",
        "iters": "Iterations", "stop": "Best place to stop", "times": " it.",
        "fig7": "Fig 7. Least-squares gradient descent first improves, then starts reconstructing the noise",
        "tikd": "Smoothness prior (derivative Tikhonov)", "tv": "Piecewise-constant prior (TV)",
        "fig8": "Fig 8. For a piecewise-constant signal, a sparse-gradient prior keeps the edges",
        "sig_ax": "Gaussian noise $\\sigma$", "ph_ax": "Photons per chip (at brightness 1)",
        "g_t": "Read noise (signal-independent)", "p_t": "Shot noise (square root of counts)",
        "bt": "Box, Tikhonov", "ct": "Coded, Tikhonov", "bv": "Box, TV", "cv": "Coded, TV",
        "fig9": "Fig 9. The coded shutter's gain depends on the noise; under shot noise the light loss costs less",
        "i_orig": "Original", "i_box": "Box-blurred", "i_code": "Coded-blurred",
        "i_bls": "Box: least squares", "i_cls": "Coded: least squares",
        "i_btk": "Box: Tikhonov", "i_ctk": "Coded: Tikhonov",
        "fig10": "Fig 10. Horizontal motion blur at 1000 photons per chip; relative error in parentheses",
    },
}

# ── 계산 ─────────────────────────────────────────────────────
rng = np.random.default_rng(0)
X = test_signal()
N = X.size
POS = np.arange(N)
A_BOX = smear_matrix(box_kernel(), N)
A_CODE = smear_matrix(coded_kernel(), N)
D = difference_matrix(N)

# 그림4: 최소자승 역변환
Y_CLEAN = A_BOX @ X
LS = {}
for key, s in (("s0", 0.0), ("s1", 1e-3), ("s2", 1e-2)):
    xs = least_squares(A_BOX, measure(A_BOX, X, s, rng) if s else Y_CLEAN)
    LS[key] = (xs, rel_error(xs, X))

# 그림5: 피카르
SIG = 1e-2
Y_N = measure(A_BOX, X, SIG, np.random.default_rng(5))
U, S, Vt = np.linalg.svd(A_BOX, full_matrices=False)
COEF = np.abs(U.T @ Y_N)
COEF0 = np.abs(U.T @ Y_CLEAN)          # 잡음 없는 측정의 성분
XCOEF = np.abs(Vt @ X)                 # 참 신호의 성분 (해가 가져야 할 크기)


def smooth(v, w=9):
    """로그 척도에서 이웃 w 개의 평균. 흩어진 점들의 추세만 보려는 것이다."""
    lv = np.log(np.maximum(v, 1e-12))
    return np.exp(np.convolve(lv, np.ones(w) / w, mode="same"))

# 그림6: 티호노프
LAMS = np.logspace(-4, 1, 60)
LBEST, EBEST, ERRS = best_lambda(lambda A, y, l: tikhonov(A, y, l), A_BOX, Y_N, X, LAMS)
RES, SOL = l_curve(A_BOX, Y_N, LAMS)
LDISC = discrepancy_lambda(A_BOX, Y_N, SIG)
EDISC = rel_error(tikhonov(A_BOX, Y_N, LDISC), X)
LSMALL, LLARGE = 1e-3, 0.5
TIK = {k: tikhonov(A_BOX, Y_N, l) for k, l in (("small", LSMALL), ("bestl", LBEST), ("large", LLARGE))}

# 그림7: 랜드웨버 (닫힌 해로 반복 횟수 스캔; verify 에서 반복 계산과 일치 확인)
ITERS = np.unique(np.logspace(0, np.log10(2e5), 80).astype(int))
LW_ERR = np.array([rel_error(landweber_svd(A_BOX, Y_N, k), X) for k in ITERS])
K_BEST = int(ITERS[np.argmin(LW_ERR)])
LW_SNAP = {k: landweber_svd(A_BOX, Y_N, k) for k in (10, K_BEST, 10000, 200000)}

# 그림8: TV vs 미분 티호노프 (같은 측정)
LD, ED, _ = best_lambda(lambda A, y, l: tikhonov(A, y, l, D), A_BOX, Y_N, X, np.logspace(-4, 1, 40))
LT, ET, _ = best_lambda(lambda A, y, l: tv_admm(A, y, l), A_BOX, Y_N, X, np.logspace(-5, -1, 20))
X_TIKD, X_TV = tikhonov(A_BOX, Y_N, LD, D), tv_admm(A_BOX, Y_N, LT)

# 그림9: 잡음 조건별 비교
G_LEVELS = np.array([1e-3, 3e-3, 1e-2, 3e-2])
P_LEVELS = np.array([10, 100, 1000, 10000])
GRES = noise_study.gaussian(tuple(G_LEVELS), methods=("tikd", "tv"))
PRES = noise_study.poisson(tuple(P_LEVELS), methods=("tikd", "tv"))

# 그림10: 2D 모션 블러, 샷 잡음 칩당 1000 광자
IMG = transform.resize(data.camera().astype(float) / 255, (256, 256), anti_aliasing=True)
PH = 1000
IM = {}
r2 = np.random.default_rng(3)
for name, kern in (("box", box_kernel()), ("code", coded_kernel())):
    A = smear_matrix(kern, 256)
    Yn = r2.poisson(np.clip(IMG @ A.T, 0, None) * PH * K) / (PH * K)
    ls = least_squares(A, Yn.T).T
    lams = np.logspace(-4, 0, 25)
    errs = [rel_error(tikhonov(A, Yn.T, l).T, IMG) for l in lams]
    tk = tikhonov(A, Yn.T, lams[int(np.argmin(errs))]).T
    IM[name] = dict(meas=Yn, ls=ls, tk=tk, e_ls=rel_error(ls, IMG), e_tk=min(errs))

SUMMARY = {
    "ls_err": {k: v[1] for k, v in LS.items()}, "lam_best": LBEST, "err_best": EBEST,
    "lam_disc": LDISC, "err_disc": EDISC, "lcurve_sol_range": (SOL.max(), SOL[int(np.argmin(ERRS))]),
    "lw_best_iter": K_BEST, "lw_best_err": LW_ERR.min(), "lw_final_err": LW_ERR[-1],
    "lw_snap_err": {k: rel_error(v, X) for k, v in LW_SNAP.items()},
    "tikd_err": ED, "tv_err": ET, "gauss": GRES, "poisson": PRES,
    "img": {k: (v["e_ls"], v["e_tk"]) for k, v in IM.items()},
    "cond": (S[0] / S[-1]),
}


# ── 그리기 ───────────────────────────────────────────────────
def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    F, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.06, 1, 1)):
        fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **F)
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림1: 개념도
    fig, ax = plt.subplots(figsize=(10, 3.6)); ax.set_axis_off()
    ax.set_xlim(0, 10); ax.set_ylim(0, 4)

    def box(x, y, w, t, c):
        ax.add_patch(FancyBboxPatch((x, y), w, 0.8, boxstyle="round,pad=0.05", fc=c, ec="none", alpha=0.9))
        ax.text(x + w / 2, y + 0.4, t, ha="center", va="center", fontsize=9, color="white", **F)

    def arrow(x0, x1, y):
        ax.annotate("", (x1, y), (x0, y), arrowprops=dict(arrowstyle="->", color="#333", lw=1.3))
    y1, y2 = 2.7, 0.6
    ax.text(0.1, 3.75, L["trad"], fontsize=10, **F)
    for x, t, c in ((0.3, L["scene"], C_TRUE), (2.6, L["lens"], C_BOX), (4.9, L["image"], C_ACC)):
        box(x, y1, 1.6, t, c)
    arrow(1.95, 2.55, y1 + 0.4); arrow(4.25, 4.85, y1 + 0.4)
    ax.text(0.1, 1.65, L["comp"], fontsize=10, **F)
    for x, t, c in ((0.3, L["scene"], C_TRUE), (2.3, L["optics"], C_BOX), (4.3, L["meas"], "#888888"),
                    (6.3, L["algo"], C_CODE), (8.3, L["image"], C_ACC)):
        box(x, y2, 1.5, t, c)
    for a in (1.85, 3.85, 5.85, 7.85):
        arrow(a, a + 0.4, y2 + 0.4)
    fig.text(0.5, 0.015, L["fig1"], ha="center", fontsize=9, **F)
    save(fig, "fig1-concept.png")

    # 그림2: 전방 모델
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.4), gridspec_kw={"width_ratios": [1.2, 1, 1.3]})
    ax[0].plot(POS, X, color=C_TRUE, lw=1.5); ax[0].set_title(L["sig_t"], fontsize=10, **F)
    ax[0].set_xlabel(L["pos"], **F); ax[0].set_ylabel(L["val"], **F)
    chips = np.arange(K)
    ax[1].bar(chips, np.ones(K), color=C_BOX, alpha=0.35, width=1.0, label=L["kbox"])
    ax[1].bar(chips, RASKAR_CODE, color=C_CODE, width=0.8, label=L["kcode"])
    ax[1].set_ylim(0, 1.6); ax[1].set_yticks([0, 1]); ax[1].set_xlabel(L["chip"], **F)
    ax[1].set_title(L["kern_t"], fontsize=10, **F); ax[1].legend(prop=LG, fontsize=7.5, loc="upper left")
    pe = np.arange(A_BOX.shape[0])
    ax[2].plot(pe, A_BOX @ X, color=C_BOX, lw=1.5, label=L["ybox"])
    ax[2].plot(pe, A_CODE @ X, color=C_CODE, lw=1.5, label=L["ycode"])
    ax[2].set_title(L["meas_t"], fontsize=10, **F); ax[2].set_xlabel(L["pos"], **F)
    ax[2].legend(prop=LG, fontsize=8)
    for a in ax: a.grid(alpha=0.3)
    cap(fig, "fig2", rect=(0, 0.07, 1, 1)); save(fig, "fig2-forward-model.png")

    # 그림3: 전달함수
    f, Hb = transfer_function(box_kernel(), 4096)
    _, Hc = transfer_function(coded_kernel(), 4096)
    fig, ax = plt.subplots(figsize=(8.5, 3.6))
    ax.semilogy(f, np.maximum(Hb, 1e-5), color=C_BOX, lw=1.2, label=L["kbox"])
    ax.semilogy(f, Hc, color=C_CODE, lw=1.2, label=L["kcode"])
    zf = np.arange(1, 6) / K
    ax.plot(zf, np.full_like(zf, 1.3e-5), "v", color=C_BOX, ms=6)
    ax.text(zf[-1] + 0.01, 1.3e-5, L["zeros"], fontsize=8.5, va="center", color=C_BOX, **F)
    ax.set_ylim(1e-5, 1.5); ax.set_xlim(0, 0.5)
    ax.set_xlabel(L["freq"], **F); ax.set_ylabel(L["mag"], **F)
    ax.legend(prop=LG, fontsize=8.5, loc="lower right"); ax.grid(alpha=0.3, which="both")
    cap(fig, "fig3", rect=(0, 0.07, 1, 1)); save(fig, "fig3-transfer-function.png")

    # 그림4: 역변환 폭주
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.3), sharey=False)
    for a, key in zip(ax, ("s0", "s1", "s2")):
        xs, e = LS[key]
        a.plot(POS, X, color=C_TRUE, lw=1.0, ls="--")
        a.plot(POS, xs, color=C_BOX, lw=1.0)
        a.set_title(f"{L[key]}  ({L['err']} {e:.2g})", fontsize=9.5, **F)
        a.set_xlabel(L["pos"], **F); a.grid(alpha=0.3)
    cap(fig, "fig4", rect=(0, 0.07, 1, 1)); save(fig, "fig4-inverse-blowup.png")

    # 그림5: 피카르
    i = np.arange(1, S.size + 1)
    sl = slice(4, -4)                   # 이동 평균의 가장자리는 버린다
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    ax[0].semilogy(i, S, color=C_TRUE, lw=1.8, label=L["sv"])
    ax[0].semilogy(i, COEF, ".", color=C_CODE, ms=2.5, alpha=0.35)
    ax[0].semilogy(i[sl], smooth(COEF)[sl], color=C_CODE, lw=1.6, label=L["coef"])
    ax[0].semilogy(i[sl], smooth(COEF0)[sl], color=C_CODE, lw=1.2, ls="--", label=L["coef0"])
    ax[0].axhline(SIG * np.sqrt(2 / np.pi), color=C_CODE, ls=":", lw=1.2)
    ax[0].text(215, SIG * np.sqrt(2 / np.pi) * 0.5, L["floor"], color=C_CODE, fontsize=8.5, **F)
    ax[0].set_title(L["p_data"], fontsize=10, **F)
    ax[1].semilogy(i, COEF / S, ".", color=C_BOX, ms=2.5, alpha=0.35)
    ax[1].semilogy(i[sl], smooth(COEF / S)[sl], color=C_BOX, lw=1.6, label=L["ratio"])
    ax[1].semilogy(i[sl], smooth(XCOEF)[sl], color=C_TRUE, lw=1.2, ls="--", label=L["xtrue"])
    ax[1].set_title(L["p_sol"], fontsize=10, **F)
    for a in ax:
        a.set_xlabel(L["idx"], **F); a.legend(prop=LG, fontsize=8.5, loc="lower left")
        a.grid(alpha=0.3, which="both")
    cap(fig, "fig5", rect=(0, 0.06, 1, 1)); save(fig, "fig5-picard.png")

    # 그림6: 티호노프
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.7), gridspec_kw={"width_ratios": [1.4, 1, 1]})
    ax[0].plot(POS, X, color=C_TRUE, lw=1.0, ls="--")
    for key, c in (("small", "#bbbbbb"), ("bestl", C_CODE), ("large", C_BOX)):
        lam = {"small": LSMALL, "bestl": LBEST, "large": LLARGE}[key]
        ax[0].plot(POS, TIK[key], color=c, lw=1.1, label=f"{L[key]} = {lam:.2g}")
    ax[0].set_ylim(-0.6, 1.6); ax[0].set_title(L["recon_t"], fontsize=10, **F)
    ax[0].set_xlabel(L["pos"], **F); ax[0].legend(prop=LG, fontsize=8)
    ax[1].loglog(LAMS, RES, color=C_TRUE, lw=1.3, label=L["res"])
    ax[1].axhline(SIG * np.sqrt(A_BOX.shape[0]), color=C_TRUE, ls=":", lw=1)
    ax[1].text(LAMS[1], SIG * np.sqrt(A_BOX.shape[0]) * 1.15, L["noiselvl"], fontsize=8, color=C_TRUE, **F)
    ax[1].axvline(LDISC, color=C_BOX, ls="--", lw=1, label=f"{L['disc']} = {LDISC:.2g}")
    ax2 = ax[1].twinx()
    ax2.semilogx(LAMS, ERRS, color=C_CODE, lw=1.3)
    ax2.plot(LBEST, EBEST, "s", mfc="none", color=C_CODE, ms=8)
    ax2.set_ylabel(L["err"], color=C_CODE, **F); ax2.tick_params(axis="y", colors=C_CODE)
    ax[1].set_xlabel("$\\lambda$"); ax[1].set_title(L["lc_t"], fontsize=10, **F)
    ax[1].legend(prop=LG, fontsize=7.5, loc="center left")
    for key, c in (("small", "#bbbbbb"), ("bestl", C_CODE), ("large", C_BOX)):
        lam = {"small": LSMALL, "bestl": LBEST, "large": LLARGE}[key]
        ax[2].semilogx(S, S**2 / (S**2 + lam**2), ".", color=c, ms=3)
    ax[2].set_xlabel(L["sv"], **F); ax[2].set_ylabel(L["ff"], **F)
    ax[2].set_title(L["ff_t"], fontsize=10, **F)
    for a in ax: a.grid(alpha=0.3, which="both")
    cap(fig, "fig6", rect=(0, 0.06, 1, 1)); save(fig, "fig6-tikhonov.png")

    # 그림7: 랜드웨버
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.7), gridspec_kw={"width_ratios": [1, 1.3]})
    ax[0].semilogx(ITERS, LW_ERR, color=C_TRUE, lw=1.4)
    ax[0].plot(K_BEST, LW_ERR.min(), "o", color=C_CODE)
    ax[0].annotate(L["stop"], (K_BEST, LW_ERR.min()), (1.2, 0.2),
                   arrowprops=dict(arrowstyle="->", color=C_CODE), color=C_CODE, fontsize=9, **F)
    ax[0].set_xlabel(L["iters"], **F); ax[0].set_ylabel(L["err"], **F)
    ax[1].plot(POS, X, color=C_TRUE, lw=1.0, ls="--")
    for k, c, lw in ((200000, C_BOX, 0.7), (10, "#999999", 1.4), (K_BEST, C_CODE, 1.4)):
        xs = LW_SNAP[k]
        ax[1].plot(POS, xs, color=c, lw=lw, label=f"{k:,}{L['times']} ({L['err']} {rel_error(xs, X):.2f})")
    ax[1].set_ylim(-0.8, 1.8); ax[1].set_xlabel(L["pos"], **F)
    ax[1].legend(prop=LG, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)
    for a in ax: a.grid(alpha=0.3, which="both")
    cap(fig, "fig7", rect=(0, 0.06, 1, 1)); save(fig, "fig7-landweber.png")

    # 그림8: TV
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.plot(POS, X, color=C_TRUE, lw=1.0, ls="--")
    ax.plot(POS, X_TIKD, color=C_BOX, lw=1.1, label=f"{L['tikd']} ({ED:.3f})")
    ax.plot(POS, X_TV, color=C_CODE, lw=1.3, label=f"{L['tv']} ({ET:.3f})")
    ax.set_xlabel(L["pos"], **F); ax.legend(prop=LG, fontsize=8.5); ax.grid(alpha=0.3)
    cap(fig, "fig8", rect=(0, 0.07, 1, 1)); save(fig, "fig8-tv.png")

    # 그림9: 잡음 조건별
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.8))
    for a, res, lv, xl, tt in ((ax[0], GRES, G_LEVELS, "sig_ax", "g_t"), (ax[1], PRES, P_LEVELS, "ph_ax", "p_t")):
        for m, mk, (lb, lc) in (("tikd", "o", ("bt", "ct")), ("tv", "s", ("bv", "cv"))):
            a.loglog(lv, [res[v][m][0] for v in lv], mk + "--", color=C_BOX, mfc="none", label=L[lb])
            a.loglog(lv, [res[v][m][1] for v in lv], mk + "-", color=C_CODE, label=L[lc])
        a.set_xlabel(L[xl], **F); a.set_ylabel(L["err"], **F); a.set_title(L[tt], fontsize=10, **F)
        a.grid(alpha=0.3, which="both")
    ax[1].invert_xaxis()
    ax[0].legend(prop=LG, fontsize=8)
    cap(fig, "fig9", rect=(0, 0.06, 1, 1)); save(fig, "fig9-coded-noise.png")

    # 그림10: 2D
    fig = plt.figure(figsize=(13, 6.6))
    gs = fig.add_gridspec(2, 4)
    a0 = fig.add_subplot(gs[:, 0]); a0.imshow(IMG, cmap="gray", vmin=0, vmax=1)
    a0.set_title(L["i_orig"], fontsize=10, **F); a0.set_axis_off()
    panels = [(0, 1, IM["box"]["meas"], "i_box", None), (0, 2, IM["box"]["ls"], "i_bls", IM["box"]["e_ls"]),
              (0, 3, IM["box"]["tk"], "i_btk", IM["box"]["e_tk"]), (1, 1, IM["code"]["meas"], "i_code", None),
              (1, 2, IM["code"]["ls"], "i_cls", IM["code"]["e_ls"]), (1, 3, IM["code"]["tk"], "i_ctk", IM["code"]["e_tk"])]
    for r, c, im, key, e in panels:
        a = fig.add_subplot(gs[r, c]); a.set_axis_off()
        a.imshow(im, cmap="gray", vmin=0, vmax=1 if e is not None else im.max())
        a.set_title(L[key] + (f" ({e:.3f})" if e is not None else ""), fontsize=9.5, **F)
    cap(fig, "fig10", rect=(0, 0.04, 1, 1)); save(fig, "fig10-2d-motion.png", dpi=110)


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 10개 저장: {L['dir']}")
    print("\n본문 수치")
    for k, v in SUMMARY.items():
        print(f"  {k}: {v}")

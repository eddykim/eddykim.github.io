"""계산 이미징과 타이코그래피 9편 그림 5개 생성 (한국어판·영문판). 수치는 studies9.py 의 캐시를 읽는다."""
import os
import textwrap

import matplotlib.pyplot as plt
import numpy as np

import studies9 as S
from ptycho9_core import align

SLUG = "ptycho-autodiff-multislice"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
COL = {"hand": "#555555", "adam_0.003": "#9fc5e8", "adam_0.01": "#2f6690", "adam_0.03": "#0b3954", "lbfgs": "#c8553d",
       1: "#c8553d", 2: "#2f6690"}

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "amp": "진폭 비용", "poisson": "푸아송 비용", "int": "세기 비용", "hand_g": "손으로 유도한 기울기 (Re)", "ad_g": "자동미분 기울기 / 2 (Re)",
        "hand": "손 구현 ML (6편)", "adam_0.003": "Adam, 걸음 0.003", "adam_0.01": "Adam, 걸음 0.01", "adam_0.03": "Adam, 걸음 0.03",
        "lbfgs": "L-BFGS", "iters": "반복 횟수", "err": "시편 오차", "sec": "시간 (초)",
        "rms0": "처음 위치 오차 (rms, 화소)", "rms": "남은 위치 오차 (rms, 화소)", "same": "보정 없음 (y = x)",
        "dz": "두 장 사이 간격 ($\\mu$m)", "resid": "데이터 잔차", "one": "한 장 모형 (곱 근사)", "two": "두 장 모형",
        "e_two1": "두 장 모형: 첫 장 오차", "e_two2": "두 장 모형: 둘째 장 오차", "e_one": "한 장 모형: 곱 O1 O2 와의 오차",
        "t1": "참 첫 장 (위상)", "t2": "참 둘째 장 (위상)", "r1_30": "두 장, 30 $\\mu$m: 첫 장", "r2_30": "두 장, 30 $\\mu$m: 둘째 장",
        "one_30": "한 장, 30 $\\mu$m", "r1_1": "두 장, 1 $\\mu$m: 첫 장", "r2_1": "두 장, 1 $\\mu$m: 둘째 장", "prod": "참 곱 O1 O2 (위상)",
        "fig1": "그림1. 세 비용에서 자동미분 기울기(PyTorch 관례의 1/2)와 손으로 유도한 기울기. 화소 2000개가 모두 대각선 위에 있다",
        "fig2": "그림2. 프로브를 모를 때 rPIE 50회 결과에서 이어 간 수렴(노이즈 표본 3개의 중앙값, 오른쪽은 표본 하나). 손 구현 ML이 가장 빠르고, Adam은 걸음에 따라 늦게 같은 답에 닿으며, L-BFGS는 멈췄다",
        "fig3": "그림3. 위치까지 미지수로 둔 자동미분의 위치 보정(오차 표본 3개). 큰 오차에서는 다 고치지 못한다",
        "fig4": "그림4. 두 위상 장 사이 간격에 따른 데이터 잔차(왼쪽)와 시편 오차(오른쪽). 간격이 좁으면 두 장이 섞인다",
        "fig5": "그림5. 두 장 모형은 간격이 넉넉하면 두 장을 따로 되찾고(위), 가까우면 데이터는 맞추지만 서로 섞인다(아래 오른쪽)",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "amp": "Amplitude cost", "poisson": "Poisson cost", "int": "Intensity cost", "hand_g": "Hand-derived gradient (Re)",
        "ad_g": "Autodiff gradient / 2 (Re)",
        "hand": "Hand-coded ML (Post 6)", "adam_0.003": "Adam, step 0.003", "adam_0.01": "Adam, step 0.01", "adam_0.03": "Adam, step 0.03",
        "lbfgs": "L-BFGS", "iters": "Iterations", "err": "Specimen error", "sec": "Time (s)",
        "rms0": "Initial position error (rms, px)", "rms": "Remaining position error (rms, px)", "same": "No correction (y = x)",
        "dz": "Separation between the two slices (µm)", "resid": "Data residual", "one": "One-slice model (multiplicative)",
        "two": "Two-slice model",
        "e_two1": "Two slices: error of slice 1", "e_two2": "Two slices: error of slice 2", "e_one": "One slice: error against O1 O2",
        "t1": "True slice 1 (phase)", "t2": "True slice 2 (phase)", "r1_30": "Two slices, 30 µm: slice 1", "r2_30": "Two slices, 30 µm: slice 2",
        "one_30": "One slice, 30 µm", "r1_1": "Two slices, 1 µm: slice 1", "r2_1": "Two slices, 1 µm: slice 2", "prod": "True product O1 O2 (phase)",
        "fig1": "Fig 1. Autodiff gradient (half, per the PyTorch convention) against the hand-derived gradient for three costs. All 2000 pixels lie on the diagonal",
        "fig2": "Fig 2. Convergence from a 50-iteration rPIE result with the probe unknown (median of 3 noise draws; right: one draw). Hand-coded ML is fastest, Adam reaches the same answer later depending on its step, L-BFGS stalled",
        "fig3": "Fig 3. Position correction by autodiff with the positions as unknowns (3 error draws). Large errors are not fully corrected",
        "fig4": "Fig 4. Data residual (left) and specimen error (right) versus the separation between two phase slices. Close slices mix",
        "fig5": "Fig 5. With enough separation the two-slice model recovers each slice (top); with close slices it fits the data but the slices mix (bottom right)",
    },
}

R = S.load()


def crop(a, m):
    ys, xs = np.where(m)
    return a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0.07, 1, 1)):
        t = fig.text(0.5, 0.015, L[key], ha="center", fontsize=9, **Fk)
        fig.canvas.draw()
        if t.get_window_extent().width > 0.98 * fig.bbox.width:  # 영문 캡션이 폭을 넘으면 줄바꿈
            text = textwrap.fill(L[key], width=int(len(L[key]) * 0.95 * fig.bbox.width / t.get_window_extent().width))
            t.set_text(text); t.set_va("bottom"); t.set_y(0.012)
            rect = (rect[0], rect[1] + 0.03 * text.count("\n"), rect[2], rect[3])
        fig.tight_layout(rect=rect)

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림1: 기울기
    G = R["grads"]
    fig, ax = plt.subplots(1, 3, figsize=(12.5, 4.2))
    for a, kind in zip(ax, ("amp", "poisson", "int")):
        h, d = np.real(G[kind]["hand"]), np.real(G[kind]["ad"]) / 2
        a.plot(h, d, ".", ms=2.5, color="#2f6690")
        lim = np.abs(h).max() * 1.05
        a.plot([-lim, lim], [-lim, lim], "k-", lw=0.6)
        a.set_title(f"{L[kind]}: {G[kind]['rel']:.0e}", fontsize=9.5, **Fk)
        a.set_xlabel(L["hand_g"], **Fk); a.set_ylabel(L["ad_g"], **Fk)
        a.ticklabel_format(style="sci", scilimits=(-2, 2))
    cap(fig, "fig1"); save(fig, "fig1-gradient-check.png")

    # 그림2: 옵티마이저
    O = R["optim"]
    CK = np.array(S.CK)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    for name in ("hand", "adam_0.003", "adam_0.01", "adam_0.03", "lbfgs"):
        e = np.array([[O[sd][name][k][0] for k in S.CK] for sd in S.SEEDS])
        ax[0].loglog(CK, np.median(e, 0), "o-", ms=3.5, color=COL[name], label=L[name])
        tt = np.array([O[0][name][k][1] for k in S.CK])
        ax[1].loglog(tt, e[0], "o-", ms=3.5, color=COL[name], label=L[name])
    for a, key in zip(ax, ("iters", "sec")):
        a.set_xlabel(L[key], **Fk); a.set_ylabel(L["err"], **Fk); a.grid(alpha=0.3, which="both")
    ax[0].legend(prop=LG, fontsize=8.5)
    cap(fig, "fig2"); save(fig, "fig2-optimizers.png")

    # 그림3: 위치
    P = R["positions"]
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    r0 = [P[(s, sd, True)]["rms0"] for s in (0.5, 1.0, 2.0) for sd in S.SEEDS]
    r1 = [P[(s, sd, True)]["rms"] for s in (0.5, 1.0, 2.0) for sd in S.SEEDS]
    ax.loglog(r0, r1, "o", color=COL[2])
    xx = np.array([0.5, 4])
    ax.loglog(xx, xx, "k--", lw=0.8, label=L["same"])
    ax.set_xlabel(L["rms0"], **Fk); ax.set_ylabel(L["rms"], **Fk); ax.grid(alpha=0.3, which="both"); ax.legend(prop=LG, fontsize=9)
    cap(fig, "fig3", rect=(0, 0.08, 1, 1)); save(fig, "fig3-autodiff-positions.png")

    # 그림4: 다층 간격
    M = R["multislice"]["runs"]
    dzs = sorted({dz for dz, ns in M})
    x = np.array(dzs) * 1e6
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    for ns, key in ((1, "one"), (2, "two")):
        ax[0].loglog(x, [M[(dz, ns)]["resid"][500] for dz in dzs], "o-", color=COL[ns], label=L[key])
    ax[1].semilogx(x, [M[(dz, 1)]["e_prod"] for dz in dzs], "o-", color=COL[1], label=L["e_one"])
    ax[1].semilogx(x, [M[(dz, 2)]["e1"] for dz in dzs], "s-", color=COL[2], label=L["e_two1"])
    ax[1].semilogx(x, [M[(dz, 2)]["e2"] for dz in dzs], "^--", color=COL[2], label=L["e_two2"])
    for a, key in zip(ax, ("resid", "err")):
        a.set_xlabel(L["dz"], **Fk); a.set_ylabel(L[key], **Fk); a.grid(alpha=0.3, which="both"); a.legend(prop=LG, fontsize=8.5)
    cap(fig, "fig4"); save(fig, "fig4-multislice-separation.png")

    # 그림5: 다층 영상
    MS = R["multislice"]
    m, O1, O2 = MS["mask"], MS["O1"], MS["O2"]
    fig, ax = plt.subplots(2, 4, figsize=(13, 7))
    panels = [(O1, O1, "t1"), (O2, O2, "t2"), (M[(30e-6, 2)]["slices"][0], O1, "r1_30"), (M[(30e-6, 2)]["slices"][1], O2, "r2_30"),
              (O1 * O2, O1 * O2, "prod"), (M[(30e-6, 1)]["slices"][0], O1 * O2, "one_30"),
              (M[(1e-6, 2)]["slices"][0], O1, "r1_1"), (M[(1e-6, 2)]["slices"][1], O2, "r2_1")]
    for a, (img, ref, key) in zip(ax.ravel(), panels):
        ph = np.angle(align(img, ref, m))
        a.imshow(crop(np.where(m, ph, np.nan), m), cmap="gray", vmin=0, vmax=1.6 if key in ("prod", "one_30") else 0.8)
        a.set_title(L[key], fontsize=9.5, **Fk); a.set_axis_off()
    cap(fig, "fig5", rect=(0, 0.04, 1, 1)); save(fig, "fig5-multislice-images.png", dpi=120)


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 5개 저장: {L['dir']}")

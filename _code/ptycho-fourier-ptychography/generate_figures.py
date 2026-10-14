"""계산 이미징과 타이코그래피 8편 그림 6개 생성 (한국어판·영문판). 수치는 studies8.py 의 캐시를 읽는다."""
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

import studies8 as S
from ptycho8_core import NH, NL, align

SLUG = "ptycho-fourier-ptychography"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
COL = {"rpie": "#c8553d", "amp": "#2f6690", "poisson": "#6a994e", "int": "#9a6fb0", "center": "#2f6690", "c": "#e09f3e", 1: "#999999"}

LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
        "spec": "시편 스펙트럼 (로그)", "center": "수직 조명 (명시야)", "bright": "비스듬한 조명 (명시야)", "dark": "더 비스듬한 조명 (암시야)",
        "t_amp": "참 시편: 진폭", "t_ph": "참 시편: 위상", "raw": "카메라 영상 한 장", "lowph": "렌즈 개구수만큼의 위상 (참고)",
        "fp_amp": "푸리에 타이코그래피: 진폭", "fp_ph": "푸리에 타이코그래피: 위상",
        "ovl": "이웃한 동공 원이 겹친 넓이", "err": "영상 오차 (덮인 대역 기준)", "c_order": "가운데 LED 부터", "r_order": "무작위 순서",
        "pt": "참 동공 위상 (수차)", "pr": "복원한 동공 위상", "no_p": "동공을 모른 채", "with_p": "동공도 함께",
        "no_pupil": "동공을 모른 채 복원", "epry": "동공도 함께 복원", "small": "수차 0.63 rad rms", "large": "수차 1.27 rad rms",
        "photons": "가장 밝은 영상의 광자 수", "rpie": "rPIE", "amp": "ML (진폭)", "poisson": "ML (푸아송)", "int": "ML (세기)",
        "fig1": "그림1. 조명 각도마다 렌즈의 동공 원(흰 원)이 스펙트럼의 다른 부분을 통과시킨다. 원이 중심을 벗어나면 암시야 영상이 된다",
        "fig2": "그림2. LED 289개의 저해상 영상으로 복원한 고해상 시편. 담을 수 있는 가장 높은 공간 주파수가 축 방향으로 렌즈의 약 7.7배다",
        "fig3": "그림3. 이웃한 LED의 동공 원이 겹친 넓이에 따른 영상 오차. 40 % 근처 아래로 내려가면 갑자기 실패한다",
        "fig4": "그림4. 수차가 큰 렌즈(1.27 rad rms)에서 동공을 함께 찾으면 수차와 시편이 함께 복원된다(겹침 79 %)",
        "fig5": "그림5. 겹침에 따라 동공을 모른 채 / 함께 찾을 때의 영상 오차. 동공까지 찾으려면 겹침이 훨씬 커야 한다",
        "fig6": "그림6. 노이즈 아래 비용 함수별 영상 오차(노이즈 표본 3개). 진폭·푸아송 ML이 가장 낮고, 세기 비용과 rPIE는 1.2–1.6배 높다",
    },
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
        "spec": "Specimen spectrum (log)", "center": "Normal illumination (brightfield)", "bright": "Oblique illumination (brightfield)",
        "dark": "More oblique illumination (darkfield)",
        "t_amp": "True specimen: amplitude", "t_ph": "True specimen: phase", "raw": "One camera image", "lowph": "Phase within the lens NA (reference)",
        "fp_amp": "Fourier ptychography: amplitude", "fp_ph": "Fourier ptychography: phase",
        "ovl": "Overlap area of neighbouring pupil circles", "err": "Image error (within the covered band)", "c_order": "Centre LED first",
        "r_order": "Random order",
        "pt": "True pupil phase (aberration)", "pr": "Recovered pupil phase", "no_p": "Pupil unknown", "with_p": "Pupil recovered too",
        "no_pupil": "Pupil unknown", "epry": "Pupil recovered too", "small": "Aberration 0.63 rad rms", "large": "Aberration 1.27 rad rms",
        "photons": "Photons in the brightest image", "rpie": "rPIE", "amp": "ML (amplitude)", "poisson": "ML (Poisson)", "int": "ML (intensity)",
        "fig1": "Fig 1. For each illumination angle the lens pupil (white circle) passes a different part of the spectrum. Off-centre, the image becomes darkfield",
        "fig2": "Fig 2. High-resolution specimen recovered from 289 low-resolution images. The highest spatial frequency along an axis is about 7.7× the lens's",
        "fig3": "Fig 3. Image error versus the overlap area of neighbouring LED pupils. Below about 40 % it fails abruptly",
        "fig4": "Fig 4. For a strongly aberrated lens (1.27 rad rms), recovering the pupil too restores both aberration and specimen (79 % overlap)",
        "fig5": "Fig 5. Image error with the pupil unknown or recovered, versus overlap. Recovering the pupil needs much more overlap",
        "fig6": "Fig 6. Image error by cost function under noise (3 noise draws). Amplitude and Poisson ML are lowest; the intensity cost and rPIE are 1.2–1.6× higher",
    },
}

R = S.load()


def crop(a, m):
    ys, xs = np.where(m)
    return a[ys.min():ys.max() + 1, xs.min():xs.max() + 1]


def draw(L):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def cap(fig, key, rect=(0, 0, 1, 1)):
        # 캡션은 본문에 둔다. 이미지에는 넣지 않는다.
        fig.tight_layout(rect=(0, 0, 1, rect[3]))

    def save(fig, name, dpi=150):
        fig.savefig(os.path.join(out, name), dpi=dpi); plt.close(fig)

    # 그림1: 스펙트럼과 원본 영상
    St = R["setup"]
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.9))
    ax[0].imshow(np.log10(St["spectrum"] + 1e-6 * St["spectrum"].max()), cmap="magma")
    for key, (ky, kx) in St["picks"].items():
        ax[0].add_patch(Circle((NH // 2 + kx, NH // 2 + ky), St["R"], fill=False, ec="w", lw=1.4))
    ax[0].set_title(L["spec"], fontsize=9.5, **Fk); ax[0].set_axis_off()
    for a, key in zip(ax[1:], ("center", "bright", "dark")):
        a.imshow(St["images"][key], cmap="gray"); a.set_title(L[key], fontsize=9.5, **Fk); a.set_axis_off()
    cap(fig, "fig1", rect=(0, 0.06, 1, 0.94)); save(fig, "fig1-fp-principle.png")

    # 그림2: 합성 개구수
    Sy = R["synth"]
    m = Sy["mask"]; obj = Sy["obj"]; rec = align(Sy["rec"], obj, m)
    fig, ax = plt.subplots(2, 3, figsize=(11, 7.6))
    ph0 = np.angle(obj)[m]; pv = (ph0.min(), ph0.max())
    av = (0, np.abs(obj).max())
    ax[0, 0].imshow(crop(np.abs(obj), m), cmap="gray", vmin=av[0], vmax=av[1]); ax[0, 0].set_title(L["t_amp"], fontsize=9.5, **Fk)
    raw = np.sqrt(Sy["I0"]); s = NH // NL
    ax[0, 1].imshow(crop(np.kron(raw, np.ones((s, s))), m), cmap="gray"); ax[0, 1].set_title(L["raw"], fontsize=9.5, **Fk)
    ax[0, 2].imshow(crop(np.abs(rec), m), cmap="gray", vmin=av[0], vmax=av[1]); ax[0, 2].set_title(L["fp_amp"], fontsize=9.5, **Fk)
    ax[1, 0].imshow(crop(np.angle(obj), m), cmap="gray", vmin=pv[0], vmax=pv[1]); ax[1, 0].set_title(L["t_ph"], fontsize=9.5, **Fk)
    ax[1, 1].imshow(crop(np.angle(align(Sy["low"], obj, m)), m), cmap="gray", vmin=pv[0], vmax=pv[1]); ax[1, 1].set_title(L["lowph"], fontsize=9.5, **Fk)
    ax[1, 2].imshow(crop(np.angle(rec), m), cmap="gray", vmin=pv[0], vmax=pv[1]); ax[1, 2].set_title(L["fp_ph"], fontsize=9.5, **Fk)
    for a in ax.ravel():
        a.set_axis_off()
    cap(fig, "fig2", rect=(0, 0.04, 1, 1)); save(fig, "fig2-synthetic-aperture.png", dpi=120)

    # 그림3: 겹침
    O = R["overlap"]
    steps = sorted({st for st, od in O})
    x = np.array([O[(st, None)]["ovl"] for st in steps]) * 100
    fig, ax = plt.subplots(figsize=(7.8, 4.3))
    ax.semilogy(x, [O[(st, None)]["err"] for st in steps], "o-", color=COL["center"], label=L["c_order"])
    for od in (1, 2):
        ax.semilogy(x, [O[(st, od)]["err"] for st in steps], "x", color="#888888", label=L["r_order"] if od == 1 else None)
    ax.set_xlabel(L["ovl"] + " (%)", **Fk); ax.set_ylabel(L["err"], **Fk); ax.grid(alpha=0.3, which="both"); ax.legend(prop=LG, fontsize=9)
    cap(fig, "fig3", rect=(0, 0.08, 1, 1)); save(fig, "fig3-overlap.png")

    # 그림4: 수차 복원 (겹침 79 %, 큰 수차)
    Pp = R["pupil"]
    v = Pp[(4, S.ABS[1])]
    sup = np.abs(v["Pt"]) > 0
    Prel = v["P"] * np.exp(-1j * np.angle(np.mean(v["P"][sup] * np.conj(v["Pt"][sup]))))
    fig, ax = plt.subplots(1, 4, figsize=(14, 3.9))
    for a, P_, t in ((ax[0], v["Pt"], "pt"), (ax[1], Prel, "pr")):
        a.imshow(np.where(sup, np.angle(P_), np.nan)[20:44, 20:44], cmap="twilight", vmin=-np.pi, vmax=np.pi)
        a.set_title(L[t], fontsize=9.5, **Fk); a.set_axis_off()
    for a, o, t, e in ((ax[2], v["rec_no"], "no_p", v["no_pupil"]), (ax[3], v["rec"], "with_p", v["epry"])):
        a.imshow(crop(np.angle(align(o, R["synth"]["obj"], R["synth"]["mask"])), R["synth"]["mask"]), cmap="gray", vmin=pv[0], vmax=pv[1])
        a.set_title(f"{L[t]}: {e:.2g}", fontsize=9.5, **Fk); a.set_axis_off()
    cap(fig, "fig4", rect=(0, 0.06, 1, 0.94)); save(fig, "fig4-pupil-recovery.png")

    # 그림5: 겹침에 따른 동공 복원
    fig, ax = plt.subplots(figsize=(8, 4.4))
    sts = (10, 6, 4)
    for ab, mk, key in ((S.ABS[0], "o", "small"), (S.ABS[1], "s", "large")):
        xo = [Pp[(st, ab)]["ovl"] * 100 for st in sts]
        ax.semilogy(xo, [Pp[(st, ab)]["no_pupil"] for st in sts], mk + "--", color=COL["rpie"], label=f"{L['no_pupil']}, {L[key]}")
        ax.semilogy(xo, [Pp[(st, ab)]["epry"] for st in sts], mk + "-", color=COL["amp"], label=f"{L['epry']}, {L[key]}")
    ax.set_xlabel(L["ovl"] + " (%)", **Fk); ax.set_ylabel(L["err"], **Fk); ax.grid(alpha=0.3, which="both"); ax.legend(prop=LG, fontsize=8.5)
    cap(fig, "fig5", rect=(0, 0.08, 1, 1)); save(fig, "fig5-pupil-vs-overlap.png")

    # 그림6: 노이즈와 비용 함수
    N_ = R["noise"]
    doses = sorted({ph for ph, mo, sd in N_})
    fig, ax = plt.subplots(figsize=(7.8, 4.4))
    for mo in ("rpie", "amp", "poisson", "int"):
        e = np.array([[N_[(ph, mo, sd)]["err"] for sd in S.SEEDS] for ph in doses])
        ax.loglog(doses, np.median(e, 1), "o-", color=COL[mo], label=L[mo])
        ax.fill_between(doses, e.min(1), e.max(1), color=COL[mo], alpha=0.2)
    ax.set_xlabel(L["photons"], **Fk); ax.set_ylabel(L["err"], **Fk); ax.grid(alpha=0.3, which="both"); ax.legend(prop=LG, fontsize=9)
    cap(fig, "fig6", rect=(0, 0.08, 1, 1)); save(fig, "fig6-noise-costs.png")


if __name__ == "__main__":
    for lang, L in LABELS.items():
        draw(L)
        print(f"[{lang}] 그림 6개 저장: {L['dir']}")

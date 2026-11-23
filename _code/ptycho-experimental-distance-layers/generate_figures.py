"""계산 이미징과 타이코그래피 13편 그림 (한국어판·영문판). 캡션은 본문에 둔다.

입력: DZ_DATA 의 ptyLab_zPIE.h5·synth_real.npz, DZ_RUNS(기본 DZ_DATA/runs)의 results13*.csv 와 *.npz
      (studies13.py run/two/score 의 출력). 입력이 없는 그림은 건너뛰고 알린다.
복원한 USAF 는 180° 돌아간 채 나온다(시편·프로브를 함께 뒤집고 켤레를 취한 쌍둥이 해와 같은 데이터를 준다).
숫자를 읽을 수 있게 그림에서만 180° 돌려 보인다.
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

import dz_core as C  # noqa: E402

SLUG = "ptycho-experimental-distance-layers"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
RUNS = os.environ.get("DZ_RUNS", os.path.join(C.DATA, "runs"))
_have = {f.name for f in font_manager.fontManager.ttflist}
KNAME = next((n for n in ("AppleGothic", "Malgun Gothic", "NanumGothic") if n in _have), None)
KFONT = {"fontfamily": KNAME} if KNAME else {}
LFONT = {"family": KNAME} if KNAME else {}
BLUE, RED, GREY, GREEN, ORANGE, PURPLE = "#2f6690", "#c8553d", "#555555", "#6a994e", "#dd8a2c", "#7b5ea7"
PAPER = (34.92, 34.99)          # zPIE 논문이 200 반복 안에 수렴한 z 의 범위 [mm]
NOISE_FLOOR = 3.48              # 합성 데이터에서 정답 시편·프로브(세기 척도 맞춤)의 데이터 오차
plt.rcParams["axes.unicode_minus"] = False

L10N = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
           "pattern": "회절 패턴 하나 (√세기)", "det_mm": "검출기 위치 (mm)",
           "scan": "주사 위치 201개와 빔 크기(FWHM)", "x_um": "x (µm)", "y_um": "y (µm)",
           "beam": "빔 FWHM 572 µm",
           "R70_O": "곡률 반경 70 mm 로 시작: |O|", "R70_P": "|P| — 한 점으로 무너짐",
           "R7_O": "곡률 반경 7 mm 로 시작: |O|", "R7_P": "|P|",
           "grid": "주사 격자: 참(●)과 z 를 3 % 크게 잡은 것(○)", "arrow10": "어긋남 ×10",
           "eps": "가정한 z 의 오차 δz/z (%)", "err": "데이터 오차", "err50": "데이터 오차 (50반복)", "ncc": "정답 시편과의 상관",
           "truth": "정답", "synth_img": "δz/z = {e:+.0f} % 로 복원",
           "pred_tv": "예측 εR", "meas_tv": "가장 선명한 거리", "pred_R": "예측 R(1-ε)", "meas_R": "복원 프로브 곡률",
           "probeR": "프로브 곡률 반경 (mm)", "err_oracle": "데이터 오차 (정답 프로브로 출발)", "floor": "잡음 바닥",
           "dz": "시편을 전파한 거리 (µm)", "tv": "TV / 최솟값", "tvcurves": "복원 시편의 선명도 (가정한 z 별)",
           "zassumed": "가정한 z (mm)", "tvmin": "가장 선명한 전파 거리 (µm)",
           "fit": "맞춤: 기울기 {k:.3f}, 0 지점 {z0:.2f} mm", "pred": "예측 R/z = {k:.3f} (R = {R:.1f} mm)",
           "zpie_pts": "zPIE 끝", "paper": "논문 34.92–34.99 mm", "errmin": "데이터 오차 최소 35.00 mm",
           "zpie_line": "zPIE 35.30 mm", "step": "보폭 {s}", "seed0": "시드 0", "seed1": "시드 1",
           "iter": "반복", "z_mm": "z (mm)", "zpie_real": "실데이터: zPIE 의 z (시작점 넷, 보폭 30·100)",
           "zpie_syn": "합성(실데이터 프로브, 정답 35 mm)", "true": "정답 35 mm",
           "barz": "막대 주기로 거꾸로 구한 z (mm)", "identity": "y = x",
           "bars_at": "가정 z = {z:.0f} mm · 4군 1요소 주기 {p:.1f} px", "usaf35": "가정 z = 35.0 mm 복원 |O| (180° 돌려 보임)",
           "gap": "두 장 사이 간격 Δz (µm)", "one": "Δz = 0 (한 장과 같음)", "front": "앞 장 |O₁|", "back": "뒤 장 |O₂|",
           "gap_img": "Δz = {g:.0f} µm"},
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
           "pattern": "One diffraction pattern (√intensity)", "det_mm": "Detector position (mm)",
           "scan": "201 scan positions and beam size (FWHM)", "x_um": "x (µm)", "y_um": "y (µm)",
           "beam": "Beam FWHM 572 µm",
           "R70_O": "Start with curvature radius 70 mm: |O|", "R70_P": "|P| — collapsed to a point",
           "R7_O": "Start with curvature radius 7 mm: |O|", "R7_P": "|P|",
           "grid": "Scan grid: true (●) and with z taken 3 % too large (○)", "arrow10": "offset ×10",
           "eps": "Error in assumed z, δz/z (%)", "err": "Data error", "err50": "Data error (50 iterations)", "ncc": "Correlation with true object",
           "truth": "Truth", "synth_img": "Reconstructed with δz/z = {e:+.0f} %",
           "pred_tv": "Prediction εR", "meas_tv": "Sharpest distance", "pred_R": "Prediction R(1-ε)", "meas_R": "Recovered probe curvature",
           "probeR": "Probe curvature radius (mm)", "err_oracle": "Data error (start from true probe)", "floor": "Noise floor",
           "dz": "Propagation distance of the object (µm)", "tv": "TV / minimum", "tvcurves": "Sharpness of the reconstructed object (per assumed z)",
           "zassumed": "Assumed z (mm)", "tvmin": "Sharpest propagation distance (µm)",
           "fit": "Fit: slope {k:.3f}, zero at {z0:.2f} mm", "pred": "Prediction R/z = {k:.3f} (R = {R:.1f} mm)",
           "zpie_pts": "zPIE end points", "paper": "Paper 34.92–34.99 mm", "errmin": "Data-error minimum 35.00 mm",
           "zpie_line": "zPIE 35.30 mm", "step": "Step {s}", "seed0": "Seed 0", "seed1": "Seed 1",
           "iter": "Iteration", "z_mm": "z (mm)", "zpie_real": "Measured data: z under zPIE (four starts, step 30/100)",
           "zpie_syn": "Synthetic (measured probe, truth 35 mm)", "true": "Truth 35 mm",
           "barz": "z inferred from bar period (mm)", "identity": "y = x",
           "bars_at": "Assumed z = {z:.0f} mm · group 4/1 period {p:.1f} px", "usaf35": "Reconstruction at z = 35.0 mm, |O| (rotated 180° for display)",
           "gap": "Gap between the two slices Δz (µm)", "one": "Δz = 0 (same as single slice)", "front": "Front |O₁|", "back": "Back |O₂|",
           "gap_img": "Δz = {g:.0f} µm"},
}


def rows(name="results13.csv"):
    path = os.path.join(RUNS, name)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        out = list(csv.DictReader(fh))
    for r in out:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return out


def run(name):
    path = os.path.join(RUNS, name + ".npz")
    return np.load(path) if os.path.exists(path) else None


def show_obj(ax, O, pos=None, N=364, frac=0.25, flip=True, vmax_pct=99.7):
    A = np.abs(O)
    if pos is not None:
        A = C.scan_region(A, pos, N, frac)
    if flip:
        A = A[::-1, ::-1]
    ax.imshow(A, cmap="gray", vmax=np.percentile(A, vmax_pct))
    ax.set_xticks([]); ax.set_yticks([])


def save(fig, L, name, dpi=110):
    os.makedirs(L["dir"], exist_ok=True)
    path = os.path.join(L["dir"], name)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.1)
    plt.close(fig)
    print("  ", path)


# ── 그림 1: 데이터와 무너짐 ─────────────────────────────────────────────────

def fig1(L, d):
    bad, good = run("start_R70"), run("start_R7")
    if bad is None or good is None:
        print("  그림1: start_R70 / start_R7 없음 — 건너뜀"); return
    Fk = L["font"]
    fig, ax = plt.subplots(2, 3, figsize=(13.5, 8.6))
    half = d["N"] / 2 * d["dxd"] * 1e3
    ax[0, 0].imshow(np.sqrt(d["I"][150]), cmap="gray", extent=[-half, half, half, -half])
    ax[0, 0].set_title(L["pattern"], fontsize=10, **Fk); ax[0, 0].set_xlabel(L["det_mm"], **Fk)
    e = d["encoder"] * 1e6
    ax[1, 0].plot(e[:, 1], e[:, 0], ".", ms=3, color=BLUE)
    ax[1, 0].add_patch(Circle((0, 0), 286, fill=False, color=RED, lw=1.2, label=L["beam"]))
    ax[1, 0].set_aspect("equal"); ax[1, 0].invert_yaxis()
    ax[1, 0].set_xlabel(L["x_um"], **Fk); ax[1, 0].set_ylabel(L["y_um"], **Fk)
    ax[1, 0].set_title(L["scan"], fontsize=10, **Fk); ax[1, 0].legend(prop=L["legend"], fontsize=8, loc="upper right")
    for col, (r, kO, kP) in enumerate(((bad, "R70_O", "R70_P"), (good, "R7_O", "R7_P")), 1):
        O = r["O"][200:670, 200:670]
        show_obj(ax[0, col], O)
        ax[0, col].set_title(L[kO], fontsize=10, **Fk)
        ax[1, col].imshow(np.abs(r["P"]), cmap="magma"); ax[1, col].set_xticks([]); ax[1, col].set_yticks([])
        ax[1, col].set_title(L[kP], fontsize=10, **Fk)
    fig.tight_layout()
    save(fig, L, "fig1-data-collapse.png")


# ── 그림 2: 합성 — 거리 오차가 하는 일 ─────────────────────────────────────

def fig2(L, R):
    syn = [r for r in R if r["name"].startswith("synr1_z")] if R else []          # l2 없이 100 반복
    orc = sorted([r for r in R if r["name"].startswith("synr1_oracle_z")], key=lambda r: r["z0"]) if R else []
    truth_path = os.path.join(C.DATA, "synth_real.npz")
    if len(syn) < 5 or not os.path.exists(truth_path):
        print("  그림2: synthr 결과 없음 — 건너뜀"); return
    from studies13 import probe_R, get_data
    Fk, LG = L["font"], L["legend"]
    T = np.load(truth_path)
    d = C.load("zpie")
    R0 = probe_R(T["P"].astype(complex), get_data("synthr")) * 1e3      # 정답 프로브의 곡률 반경 [mm]
    fig = plt.figure(figsize=(15, 8.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.05])
    a = fig.add_subplot(gs[0, 0])
    No = C.grid_size(d, 0.028)
    p_bad = d["encoder"] / C.dxo(d, 0.035 * 1.03) + No // 2          # 실수 화소(반올림 전)로 어긋남을 보인다
    p_true = d["encoder"] / C.dxo(d, 0.035) + No // 2
    a.plot(p_true[:, 1], p_true[:, 0], "o", ms=3, color=BLUE)
    a.plot(p_bad[:, 1], p_bad[:, 0], "o", ms=3, mfc="none", color=RED)
    for t, b in zip(p_true, p_bad):
        a.plot([t[1], t[1] + 10 * (b[1] - t[1])], [t[0], t[0] + 10 * (b[0] - t[0])], color=GREY, lw=0.6)
    a.set_aspect("equal"); a.invert_yaxis(); a.set_xticks([]); a.set_yticks([])
    a.set_title(L["grid"] + "\n(" + L["arrow10"] + ")", fontsize=9, **Fk)
    syn.sort(key=lambda r: r["z0"])
    eps = np.array([(r["z0"] / 35 - 1) * 100 for r in syn])
    b = fig.add_subplot(gs[0, 1:3])
    b.plot(eps, [r["err_final"] for r in syn], "o-", color=BLUE, label=L["err"])
    if orc:
        b.plot([(r["z0"] / 35 - 1) * 100 for r in orc], [r["err_final"] for r in orc], "^:", color=PURPLE, label=L["err_oracle"])
    b.axhline(NOISE_FLOOR, color=GREY, lw=0.8, ls="--", label=L["floor"])
    b.set_ylim(0, max(r["err_final"] for r in syn + orc) * 1.3)
    b.set_xlabel(L["eps"], **Fk); b.set_ylabel(L["err"], color=BLUE, **Fk); b.grid(alpha=0.3)
    b2 = b.twinx()
    b2.plot(eps, [r["ncc"] for r in syn], "s--", color=RED, label=L["ncc"])
    b2.set_ylabel(L["ncc"], color=RED, **Fk); b2.set_ylim(0, 1)
    h1, l1 = b.get_legend_handles_labels(); h2, l2 = b2.get_legend_handles_labels()
    b.legend(h1 + h2, l1 + l2, prop=LG, fontsize=8, loc="lower center")
    c = fig.add_subplot(gs[0, 3])
    ee = np.linspace(eps.min(), eps.max(), 50)
    c.plot(ee, ee / 100 * R0 * 1e3, "--", color=GREY, label=L["pred_tv"])
    c.plot(eps, [r["tv_min"] for r in syn], "o", color=BLUE, label=L["meas_tv"])
    c.set_xlabel(L["eps"], **Fk); c.set_ylabel(L["tvmin"], **Fk); c.grid(alpha=0.3)
    c2 = c.twinx()
    c2.plot(ee, R0 * (1 - ee / 100), ":", color=RED, label=L["pred_R"])
    c2.plot(eps, [r["probe_R"] for r in syn], "s", color=RED, mfc="none", label=L["meas_R"])
    c2.set_ylabel(L["probeR"], color=RED, **Fk)
    h1, l1 = c.get_legend_handles_labels(); h2, l2 = c2.get_legend_handles_labels()
    c.legend(h1 + h2, l1 + l2, prop=LG, fontsize=7, loc="upper left")
    def zoom(O, h=90):                                    # 배열 가운데 ±h 화소: 배율과 흐림이 보이게
        c = O.shape[-1] // 2
        return O[c - h:c + h, c - h:c + h]
    ax_t = fig.add_subplot(gs[1, 0])
    show_obj(ax_t, zoom(T["O"]), flip=False, vmax_pct=100)
    ax_t.set_title(L["truth"], fontsize=10, **Fk)
    for k, (e, name) in enumerate(((-3, "synr1_z33.95"), (0, "synr1_z35.0"), (3, "synr1_z36.05"))):
        r = run(name)
        axk = fig.add_subplot(gs[1, k + 1])
        if r is not None:
            show_obj(axk, zoom(r["O"]), flip=False, vmax_pct=100)
        axk.set_title(L["synth_img"].format(e=e), fontsize=10, **Fk)
    fig.tight_layout()
    save(fig, L, "fig2-synthetic-error.png")


# ── 그림 3: 실데이터 — 선명도와 데이터 오차 ────────────────────────────────

def fig3(L, R):
    if not R:
        print("  그림3: results13.csv 없음 — 건너뜀"); return
    Fk, LG = L["font"], L["legend"]
    fix0 = sorted([r for r in R if r["data"] == "zpie" and r["zpie"] == 0 and r["seed"] == 0], key=lambda r: r["z0"])
    fix1 = sorted([r for r in R if r["data"] == "zpie" and r["zpie"] == 0 and r["seed"] == 1], key=lambda r: r["z0"])
    zp = [r for r in R if r["data"] == "zpie" and r["zpie"] == 1]
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.8))
    cmap = plt.get_cmap("coolwarm")
    for zz in (33.0, 34.0, 35.0, 36.0, 37.0):
        r = run(f"zfix_z{zz}")
        if r is None:
            continue
        s, dz = r["sharp"], r["sharp_dz"] * 1e6
        ax[0].plot(dz, s / s.min(), color=cmap((zz - 33) / 4), label=f"{zz:.0f} mm")
        ax[0].plot(dz[np.argmin(s)], 1.0, "v", color=cmap((zz - 33) / 4), ms=6)
    ax[0].set_xlabel(L["dz"], **Fk); ax[0].set_ylabel(L["tv"], **Fk); ax[0].grid(alpha=0.3)
    ax[0].set_title(L["tvcurves"], fontsize=10, **Fk); ax[0].legend(fontsize=8, title=L["zassumed"], title_fontproperties=LG)
    z = np.array([r["z0"] for r in fix0]); t = np.array([r["tv_min"] for r in fix0])
    k, b = np.polyfit(z, t, 1)
    z0 = -b / k
    zz = np.linspace(32.8, 37.2, 10)
    ax[1].plot(z, t, "o", color=BLUE, label=L["seed0"])
    ax[1].plot(zz, k * zz + b, "-", color=BLUE, lw=1, label=L["fit"].format(k=k / 1e3, z0=z0))
    Rm = np.mean([r["probe_R"] for r in fix0 if abs(r["z0"] - 35) < 0.01])
    kp = Rm / 35.0
    ax[1].plot(zz, kp * 1e3 * (zz - z0), "--", color=RED, lw=1, label=L["pred"].format(k=kp, R=Rm))
    ax[1].plot([r["z_final"] for r in zp], [r["tv_min"] for r in zp], "*", color=ORANGE, ms=10, label=L["zpie_pts"])
    ax[1].axhline(0, color=GREY, lw=0.6)
    ax[1].set_xlabel(L["zassumed"], **Fk); ax[1].set_ylabel(L["tvmin"], **Fk); ax[1].grid(alpha=0.3)
    ax[1].legend(prop=LG, fontsize=8, loc="upper left")
    ax[2].axvspan(*PAPER, color=GREEN, alpha=0.2, label=L["paper"])
    ax[2].plot([r["z0"] for r in fix0], [r["err_round1"] for r in fix0], "o-", color=BLUE, label=L["seed0"])
    ax[2].plot([r["z0"] for r in fix1], [r["err_round1"] for r in fix1], "s", color=PURPLE, label=L["seed1"])
    zpm = np.mean([r["z_final"] for r in zp])
    ax[2].axvline(zpm, color=ORANGE, ls="--", label=L["zpie_line"])
    ax[2].set_xlabel(L["zassumed"], **Fk); ax[2].set_ylabel(L["err50"], **Fk); ax[2].grid(alpha=0.3)
    ax[2].legend(prop=LG, fontsize=8)
    fig.tight_layout()
    save(fig, L, "fig3-sharpness-error.png")
    return dict(slope=k / 1e3, zero=z0, R=Rm, pred=kp, zpie=zpm)


# ── 그림 4: zPIE 의 궤적 ──────────────────────────────────────────────────

def fig4(L, R):
    names = [r["name"] for r in R if r["zpie"] == 1 and r["data"] == "zpie"] if R else []
    if not names:
        print("  그림4: zPIE 결과 없음 — 건너뜀"); return
    Fk, LG = L["font"], L["legend"]
    syn = [r["name"] for r in R if r["zpie"] == 1 and r["data"] == "synthr"]
    fig, ax = plt.subplots(1, 2 if syn else 1, figsize=(13 if syn else 6.5, 4.4), squeeze=False)
    a = ax[0, 0]
    a.axhspan(*PAPER, color=GREEN, alpha=0.25, label=L["paper"])
    a.axhline(35.0, color=BLUE, ls=":", label=L["errmin"])
    seen = set()
    for n in sorted(names):
        tr = run(n)["trace"]
        big = "step100" in n
        lab = None if big in seen else L["step"].format(s=100 if big else 30)
        seen.add(big)
        a.plot(np.arange(len(tr)), tr[:, 2] * 1e3, color=ORANGE if big else RED, lw=1.2,
               ls="--" if big else "-", label=lab)
    a.set_xlabel(L["iter"], **Fk); a.set_ylabel(L["z_mm"], **Fk); a.grid(alpha=0.3)
    a.set_title(L["zpie_real"], fontsize=10, **Fk); a.legend(prop=LG, fontsize=8, loc="upper right")
    if syn:
        b = ax[0, 1]
        b.axhline(35.0, color="k", lw=0.8, label=L["true"])
        for n in sorted(syn):
            tr = run(n)["trace"]
            b.plot(np.arange(len(tr)), tr[:, 2] * 1e3, color=RED, lw=1.2)
        b.set_xlabel(L["iter"], **Fk); b.set_ylabel(L["z_mm"], **Fk); b.grid(alpha=0.3)
        b.set_title(L["zpie_syn"], fontsize=10, **Fk); b.legend(prop=LG, fontsize=8)
    fig.tight_layout()
    save(fig, L, "fig4-zpie.png")


# ── 그림 5: 막대 선폭은 왜 검증이 되지 못하나 ─────────────────────────────

def fig5(L, R):
    fix0 = sorted([r for r in R if r["data"] == "zpie" and r["zpie"] == 0 and r["seed"] == 0 and r["bar_z"] != ""],
                  key=lambda r: r["z0"]) if R else []
    r33, r37, r35 = run("zfix_z33.0"), run("zfix_z37.0"), run("zfix_z35.0")
    if not fix0 or r33 is None or r37 is None or r35 is None:
        print("  그림5: 입력 없음 — 건너뜀"); return
    Fk, LG = L["font"], L["legend"]
    fig = plt.figure(figsize=(15, 5.2))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.25, 0.8, 0.8, 1.25])
    a = fig.add_subplot(gs[0])
    z = [r["z0"] for r in fix0]
    a.plot([32.8, 37.2], [32.8, 37.2], color=GREY, lw=0.8, label=L["identity"])
    a.errorbar(z, [r["bar_z"] for r in fix0], yerr=[r["bar_sem"] for r in fix0], fmt="o", color=BLUE, capsize=2)
    a.set_xlabel(L["zassumed"], **Fk); a.set_ylabel(L["barz"], **Fk); a.grid(alpha=0.3); a.legend(prop=LG, fontsize=8)
    a.set_aspect("equal")
    win = (slice(300, 400), slice(260, 420))           # 4군 1요소 (복원 화소, 뒤집기 전)
    for k, (r, zz) in enumerate(((r33, 33), (r37, 37))):
        b = fig.add_subplot(gs[1 + k])
        A = np.abs(r["O"][win])[::-1, ::-1]
        b.imshow(A, cmap="gray", vmax=np.percentile(A, 99.5))
        b.set_xticks([]); b.set_yticks([])
        import usaf
        ms = [m for m in usaf.measure(r["O"], (270, 610, 230, 600), float(r["dxo_final"]), groups=(4,)) if m["e"] == 1]
        p = float(np.mean([m["p_px"] for m in ms])) if ms else float("nan")
        b.set_title(L["bars_at"].format(z=zz, p=p), fontsize=9, **Fk)
    c = fig.add_subplot(gs[3])
    show_obj(c, r35["O"][280:600, 250:590], flip=True)
    c.set_title(L["usaf35"], fontsize=10, **Fk)
    fig.tight_layout()
    save(fig, L, "fig5-bars-circular.png")


# ── 그림 6: 두 장 모형 ─────────────────────────────────────────────────────

def fig6(L):
    R2 = rows("results13_two.csv")
    if not R2:
        print("  그림6: results13_two.csv 없음 — 건너뜀"); return
    Fk, LG = L["font"], L["legend"]
    R2.sort(key=lambda r: r["gap_um"])
    R2s1 = [r for r in R2 if r.get("seed", 0) == 1]
    R2 = [r for r in R2 if r.get("seed", 0) == 0]
    fig = plt.figure(figsize=(14, 7.6))
    gs = fig.add_gridspec(2, 4)
    a = fig.add_subplot(gs[:, 0:2])
    g = [max(r["gap_um"], 10) for r in R2]
    a.plot(g, [r["err_final"] for r in R2], "o-", color=BLUE, label=L["seed0"])
    if R2s1:
        a.plot([max(r["gap_um"], 10) for r in R2s1], [r["err_final"] for r in R2s1], "s", color=PURPLE, label=L["seed1"])
    a.annotate(L["one"], (10, R2[0]["err_final"]), textcoords="offset points", xytext=(4, 10), fontsize=8, **Fk)
    a.set_xscale("log"); a.set_xlabel(L["gap"], **Fk); a.set_ylabel(L["err"], **Fk); a.grid(alpha=0.3, which="both")
    a.legend(prop=LG, fontsize=8)
    best = min(R2[1:] or R2, key=lambda r: r["err_final"])
    for k, gap in enumerate((best["gap_um"], R2[-1]["gap_um"])):
        r = run(f"tl2_gap{gap:.0f}")
        if r is None:
            continue
        pos = C.positions(C.load("twolayer"), 0.0386, r["O1"].shape[-1])
        for j, key in enumerate(("O1", "O2")):
            ax = fig.add_subplot(gs[k, 2 + j])
            show_obj(ax, r[key], pos, flip=False)
            ax.set_title(L["gap_img"].format(g=gap) + " · " + L["front" if key == "O1" else "back"], fontsize=9, **Fk)
    fig.tight_layout()
    save(fig, L, "fig6-two-slices.png")


if __name__ == "__main__":
    R = rows()
    d = C.load("zpie")
    for lang, L in L10N.items():
        print(f"[{lang}] {L['dir']}")
        fig1(L, d)
        fig2(L, R)
        info = fig3(L, R)
        fig4(L, R)
        fig5(L, R)
        fig6(L)
    if info:
        print("그림3 수치:", {k: round(v, 4) for k, v in info.items()})

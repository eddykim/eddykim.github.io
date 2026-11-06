"""계산 이미징과 타이코그래피 12편 그림 6개 생성 (한국어판·영문판). 캡션은 본문에 둔다. 그림 안에는 넣지 않는다.
입력: data/TwoLayer_bin4.hdf5 (그림1), results.csv(data/ 에 없으면 이 폴더의 사본)와 data/runs/*.npz
(그림2–6, studies12.py 와 ptylab/run_ptylab.py 의 결과).
입력이 없는 그림은 건너뛰고 알린다.
"""
import csv
import os

import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Circle

import cpm_core as C

SLUG = "ptycho-experimental-cpm"
BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "assets", "img", "posts", SLUG)
RUNS = os.path.join(C.DATA, "runs")
KFONT = {"fontfamily": "AppleGothic"}
LFONT = {"family": "AppleGothic"}
BLUE, RED, GREY, GREEN, ORANGE = "#2f6690", "#c8553d", "#555555", "#6a994e", "#dd8a2c"

VARIANTS = ["base", "nosmooth", "nocom", "nol2", "flatprobe", "sequential", "epie",
            "nopower", "nopower_flatprobe", "nopower_sequential", "nopower_epie"]
LABELS = {
    "ko": {"dir": BASE, "font": KFONT, "legend": LFONT,
           "pattern": "회절 패턴 하나 (log 세기)", "det_mm": "검출기 위치 (mm)", "scan": "주사 위치 202개와 프로브 지름",
           "x_um": "x ($\\mu$m)", "y_um": "y ($\\mu$m)", "probe": "프로브 지름 0.6 mm",
           "sums": "패턴별 총세기", "sum_ax": "총세기 / 가장 밝은 패턴", "count": "패턴 수",
           "median": "중앙값", "fixed": "고정값 (가장 밝은 패턴)",
           "seed": "시드", "err": "데이터 오차 (202장 합)", "bright": "기본 설정 (밝게 시작)", "free": "세기 맞춤 끔",
           "ok_O": "성공 (시드 {s}) |O|", "ok_P": "성공 |P|", "bad_O": "실패 (시드 {s}) |O|", "bad_P": "실패 |P|",
           "var": {"base": "기본 설정", "nosmooth": "프로브 매끄럽게 끔", "nocom": "무게중심 맞춤 끔",
                   "nocoupling": "결합 끔", "nol2": "l2 끔", "flatprobe": "곡률 없는 초기 프로브",
                   "sequential": "순차 순서", "epie": "ePIE (모멘텀 없음)", "nopower": "세기 맞춤 끔",
                   "nopower_flatprobe": "세기 맞춤 끔 + 곡률 없음", "nopower_sequential": "세기 맞춤 끔 + 순차",
                   "nopower_epie": "세기 맞춤 끔 + ePIE"},
           "corr": "기준 복원과의 상관", "start": "시작 세기 / 고정값", "rate": "성공 비율",
           "ptylab": "PtyLab (고정 없음)", "own": "직접 구현", "iters": "반복", "trace_title": "시드 0, 첫 라운드",
           "tr": {"base": "시작 맞춤 + 매 반복 고정", "initonly": "시작 맞춤만", "iteronly": "매 반복 고정만",
                  "nopower": "둘 다 끔"},
           "dim": "어둡게 시작", "brt": "밝게 시작", "it": "반복 {i}", "curve": "데이터 오차",
           "absO": "|O| (시편 진폭)", "argO": "arg O (시편 위상)", "absP": "|P| (프로브 진폭)", "argP": "arg P (프로브 위상)"},
    "en": {"dir": os.path.join(BASE, "en"), "font": {}, "legend": {},
           "pattern": "One diffraction pattern (log intensity)", "det_mm": "Detector position (mm)",
           "scan": "202 scan positions and probe diameter",
           "x_um": "x (µm)", "y_um": "y (µm)", "probe": "Probe diameter 0.6 mm",
           "sums": "Total intensity per pattern", "sum_ax": "Total intensity / brightest pattern", "count": "Patterns",
           "median": "Median", "fixed": "Fixed value (brightest pattern)",
           "seed": "Seed", "err": "Data error (sum over 202)", "bright": "Default (bright start)", "free": "Power matching off",
           "ok_O": "Success (seed {s}) |O|", "ok_P": "Success |P|", "bad_O": "Failure (seed {s}) |O|", "bad_P": "Failure |P|",
           "var": {"base": "Default", "nosmooth": "Probe smoothing off", "nocom": "Centre-of-mass off",
                   "nocoupling": "Coupling off", "nol2": "l2 off", "flatprobe": "Initial probe without curvature",
                   "sequential": "Sequential order", "epie": "ePIE (no momentum)", "nopower": "Power matching off",
                   "nopower_flatprobe": "Power off + no curvature", "nopower_sequential": "Power off + sequential",
                   "nopower_epie": "Power off + ePIE"},
           "corr": "Correlation with reference", "start": "Starting power / fixed value", "rate": "Success rate",
           "ptylab": "PtyLab (no fixing)", "own": "Own implementation", "iters": "Iteration", "trace_title": "Seed 0, first round",
           "tr": {"base": "Start matching + per-iteration fix", "initonly": "Start matching only",
                  "iteronly": "Per-iteration fix only", "nopower": "Both off"},
           "dim": "Dim start", "brt": "Bright start", "it": "iter {i}", "curve": "Data error",
           "absO": "|O| (object amplitude)", "argO": "arg O (object phase)", "absP": "|P| (probe amplitude)",
           "argP": "arg P (probe phase)"},
}


def results():
    path = os.path.join(C.DATA, "results.csv")
    if not os.path.exists(path):
        path = os.path.join(C.HERE, "results.csv")     # 글에 쓴 채점표의 사본 (studies12.py score 의 출력)
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["seed"] = int(r["seed"]); r["error"] = float(r["error"]); r["ncc"] = float(r["ncc"])
        r["success"] = int(r["success"]); r["power_ratio"] = float(r["power_ratio"])
    return rows


def pick(rows, impl, setting):
    return sorted((r for r in rows if r["impl"] == impl and r["setting"] == setting), key=lambda r: r["seed"])


def run(name):
    p = os.path.join(RUNS, name + ".npz")
    return np.load(p, allow_pickle=True) if os.path.exists(p) else None


def crop_O(z, d):
    O = np.squeeze(z["object"])
    if "positions" in z.files:
        pos = z["positions"]
    else:
        pos = np.round(d["encoder"] / d["dxo"]).astype(int) + O.shape[-1] // 2 - d["N"] // 2
    r0, c0 = pos.min(0) + d["N"] // 4
    r1, c1 = pos.max(0) + 3 * d["N"] // 4
    return O[r0:r1, c0:c1]


def draw(L, d, rows):
    out = L["dir"]; os.makedirs(out, exist_ok=True)
    Fk, LG = L["font"], L["legend"]

    def save(fig, name, dpi=150, colors=None):
        """colors 를 주면 저장 뒤 그만큼의 색으로 줄인다(영상이 든 그림의 용량을 줄이려고)."""
        fig.tight_layout()
        path = os.path.join(out, name)
        fig.savefig(path, dpi=dpi, bbox_inches="tight", pad_inches=0.1); plt.close(fig)
        if colors:
            Image.open(path).convert("RGB").quantize(colors=colors, method=Image.Quantize.MEDIANCUT).save(path, optimize=True)

    # 그림1: 데이터 — 패턴, 주사 위치, 패턴별 총세기
    I = d["I"]; sums = I.sum((-1, -2)); mp2 = sums.max()
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.4))
    j = int(np.argsort(sums)[len(sums) // 2])
    ext = d["N"] * d["dxd"] * 1e3 / 2
    ax[0].imshow(np.log10(I[j] + 1), cmap="magma", extent=[-ext, ext, ext, -ext])
    ax[0].set_title(L["pattern"], fontsize=10, **Fk); ax[0].set_xlabel(L["det_mm"], **Fk)
    e = d["encoder"] * 1e6
    ax[1].plot(e[:, 1], e[:, 0], ".", color=BLUE, ms=4)
    k = int(np.argmin(np.hypot(*(e - e.mean(0)).T)))         # 주사 영역 가운데 위치
    ax[1].add_patch(Circle((e[k, 1], e[k, 0]), d["entrancePupilDiameter"] * 1e6 / 2, fill=False, color=RED, lw=1.5,
                           label=L["probe"]))
    ax[1].set_aspect("equal"); ax[1].set_xlabel(L["x_um"], **Fk); ax[1].set_ylabel(L["y_um"], **Fk)
    ax[1].set_title(L["scan"], fontsize=10, **Fk); ax[1].legend(prop=LG, fontsize=9, loc="upper right")
    ax[2].hist(sums / mp2, bins=30, color=GREY, alpha=0.8)
    ax[2].axvline(np.median(sums) / mp2, color=BLUE, ls="--", label=L["median"])
    ax[2].axvline(1.0, color=RED, ls="-", label=L["fixed"])
    ax[2].set_xlabel(L["sum_ax"], **Fk); ax[2].set_ylabel(L["count"], **Fk)
    ax[2].set_title(L["sums"], fontsize=10, **Fk); ax[2].legend(prop=LG, fontsize=9)
    save(fig, "fig1-data-scan.png", dpi=110, colors=128)
    if rows is None:
        print("  results.csv 가 없어 그림2–6 을 건너뛴다")
        return

    # 그림2: 시드별 성패 + 성공·실패 한 쌍 (PtyLab 기본 설정)
    br, fr = pick(rows, "ptylab", "base"), pick(rows, "ptylab", "nopower")
    ok = next(r["seed"] for r in br if r["success"]); bad = next(r["seed"] for r in br if not r["success"])
    fig = plt.figure(figsize=(15, 4.2))
    gs = fig.add_gridspec(1, 5, width_ratios=[1.8, 1, 1, 1, 1])
    a = fig.add_subplot(gs[0])
    a.plot([r["seed"] for r in br], [r["error"] for r in br], "o", color=RED, label=L["bright"])
    a.plot([r["seed"] for r in fr], [r["error"] for r in fr], "s", mfc="none", color=BLUE, label=L["free"])
    a.set_xlabel(L["seed"], **Fk); a.set_ylabel(L["err"], **Fk); a.grid(alpha=0.3); a.legend(prop=LG, fontsize=9)
    a.set_xticks(range(0, 20, 2))
    for i, (s, key) in enumerate([(ok, "ok"), (bad, "bad")]):
        z = run(f"ptylab_base_s{s}")
        O = np.abs(crop_O(z, d)); P = np.abs(np.squeeze(z["probe"]))
        b = fig.add_subplot(gs[1 + 2 * i]); b.imshow(O, cmap="gray", vmax=np.percentile(O, 99.5))
        b.set_title(L[key + "_O"].format(s=s), fontsize=10, **Fk); b.axis("off")
        c = fig.add_subplot(gs[2 + 2 * i]); c.imshow(P, cmap="magma"); c.set_title(L[key + "_P"], fontsize=10, **Fk)
        c.axis("off")
    save(fig, "fig2-seeds.png", dpi=110, colors=256)

    # 그림3: 소거 행렬 (PtyLab, 시드 0–2)
    M = np.full((len(VARIANTS), 3), np.nan)
    for i, v in enumerate(VARIANTS):
        for r in pick(rows, "ptylab", v):
            if r["seed"] < 3:
                M[i, r["seed"]] = r["ncc"]
    fig, ax = plt.subplots(figsize=(6.5, 6))
    im = ax.imshow(M, cmap="RdYlGn", norm=TwoSlopeNorm(vmin=0, vcenter=0.85, vmax=1), aspect="auto")   # 0.7 도 실패로 보이게
    for i in range(M.shape[0]):
        for k in range(3):
            if np.isfinite(M[i, k]):
                ax.text(k, i, f"{M[i, k]:.2f}", ha="center", va="center", fontsize=9)
    ax.set_yticks(range(len(VARIANTS))); ax.set_yticklabels([L["var"][v] for v in VARIANTS], **Fk)
    ax.set_xticks(range(3)); ax.set_xticklabels(["0", "1", "2"]); ax.set_xlabel(L["seed"], **Fk)
    fig.colorbar(im, ax=ax, fraction=0.04).set_label(L["corr"], **Fk)
    save(fig, "fig3-ablation.png")

    # 그림4: 시작 세기와 성공 비율 + 시드 0 의 반복별 오차
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
    p0 = 1.64e4 / 1.773e7
    pts = [(p0, "nopower")] + [(f, f"scale{f:g}") for f in (0.01, 0.1, 0.3)] + [(1.0, "initonly"), (3.0, "scale3")]
    xs, ys, ns = [], [], []
    for f, v in pts:
        rr = [r for r in pick(rows, "ptylab", v) if r["seed"] < 6]
        xs.append(f); ys.append(np.mean([r["success"] for r in rr])); ns.append(len(rr))
    ax[0].semilogx(xs, ys, "o-", color=BLUE, label=L["ptylab"])
    for x, y, n in zip(xs, ys, ns):
        ax[0].annotate(f"{round(y * n)}/{n}", (x, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
    own = [(p0, "free"), (3.0, "bright3")]
    for i, (x, v) in enumerate(own):          # 직접 구현은 두 점만 쟀으므로 선으로 잇지 않는다
        rr = pick(rows, "own", v)
        y = np.mean([r["success"] for r in rr])
        ax[0].semilogx([x], [y], "s", color=RED, mfc="none", ms=9, label=L["own"] if i == 0 else None)
        ax[0].annotate(f"{sum(r['success'] for r in rr)}/{len(rr)}", (x, y), textcoords="offset points",
                       xytext=(0, -14), ha="center", fontsize=8, color=RED)
    ax[0].set_xlabel(L["start"], **Fk); ax[0].set_ylabel(L["rate"], **Fk); ax[0].set_ylim(-0.05, 1.12)
    ax[0].grid(alpha=0.3); ax[0].legend(prop=LG, fontsize=9, loc="lower left")
    for v, col in zip(["base", "initonly", "iteronly", "nopower"], [RED, ORANGE, GREEN, BLUE]):
        z = run(f"ptylab_{v}_s0")
        if z is not None:
            t = np.asarray(z["trace"])[:50, 3]
            ax[1].plot(t, color=col, label=L["tr"][v], ls="--" if v == "initonly" else "-")
    ax[1].set_yscale("log"); ax[1].set_xlabel(L["iters"], **Fk); ax[1].set_ylabel(L["err"], **Fk)
    ax[1].set_title(L["trace_title"], fontsize=10, **Fk); ax[1].grid(alpha=0.3, which="both"); ax[1].legend(prop=LG, fontsize=9)
    save(fig, "fig4-start-power.png")

    # 그림5: 같은 시드에서 어둡게/밝게 시작한 앞쪽 반복의 프로브
    show = [0, 1, 2, 4, 8, 13, 20, 29]
    zs = {m: run(f"diag_{m}_s4") for m in ("dim", "bright")}
    if all(v is not None for v in zs.values()):
        fig = plt.figure(figsize=(16, 7))
        gs = fig.add_gridspec(3, len(show), height_ratios=[1, 1, 1.25])
        for row, m in enumerate(("dim", "bright")):
            z = zs[m]; its = list(z["its"])
            for c, i in enumerate(show):
                a = fig.add_subplot(gs[row, c]); a.imshow(np.abs(z["P"][its.index(i)]), cmap="magma"); a.axis("off")
                a.set_title((L["dim" if m == "dim" else "brt"] + ", " if c == 0 else "") + L["it"].format(i=i),
                            fontsize=9, **Fk)
        a = fig.add_subplot(gs[2, :])
        for s, ls in ((4, "-"), (0, "--"), (1, ":")):
            for m, col in (("dim", BLUE), ("bright", RED)):
                z = run(f"diag_{m}_s{s}")
                if z is not None:
                    a.plot(np.asarray(z["trace"])[:, 1], color=col, ls=ls,
                           label=f"{L['dim' if m == 'dim' else 'brt']}, {L['seed']} {s}")
        a.set_yscale("log"); a.set_xlabel(L["iters"], **Fk); a.set_ylabel(L["err"], **Fk); a.grid(alpha=0.3, which="both")
        a.legend(prop=LG, fontsize=8, ncol=3)
        save(fig, "fig5-early-probe.png", dpi=100, colors=128)

    # 그림6: 직접 구현(어둡게 시작 + 매 반복 고정)의 최종 복원
    z = run("own_dim_s0")
    O = crop_O(z, d); P = np.squeeze(z["probe"])
    fig, ax = plt.subplots(1, 4, figsize=(16, 4.3))
    ax[0].imshow(np.abs(O), cmap="gray", vmax=np.percentile(np.abs(O), 99.5)); ax[0].set_title(L["absO"], fontsize=10, **Fk)
    ax[1].imshow(np.angle(O), cmap="twilight"); ax[1].set_title(L["argO"], fontsize=10, **Fk)
    ax[2].imshow(np.abs(P), cmap="magma"); ax[2].set_title(L["absP"], fontsize=10, **Fk)
    ax[3].imshow(np.angle(P), cmap="twilight"); ax[3].set_title(L["argP"], fontsize=10, **Fk)
    for a in ax:
        a.axis("off")
    save(fig, "fig6-own-final.png", dpi=100, colors=256)


if __name__ == "__main__":
    d = C.load()
    rows = results()
    for lang, L in LABELS.items():
        draw(L, d, rows)
        print(f"[{lang}] 저장: {L['dir']}")

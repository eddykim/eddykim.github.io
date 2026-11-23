"""13편의 복원 실행. 한 번 부를 때 복원 하나를 돌리고 결과를 npz 로 남긴다.

  python studies13.py run --data zpie --z0 35.0 [--zpie] [--seed 0] [--rounds 2] --out runs/zpie_z35.0.npz
      --data   zpie | synth | synthr (synth.py 가 만든 합성 데이터: 완만한 빔 | 실데이터 프로브)
      --z0     시작 거리 [mm]. --zpie 가 없으면 이 거리로 끝까지 고정한다
      --zpie   Algorithm 1 로 z 를 고친다 (--step, --sign)
      --curvR  시작 프로브의 곡률 반경 [mm] (zPIE 데이터 기본 7, 합성 완만한 빔 70)
  설정은 12편이 권한 '어둡게 시작 + 매 반복 Σ|P|² 고정'(dim)이 기본이다.

결과 npz: O(가운데 crop 없이 전체, complex64), P, trace[반복, (Σ|P|², 오차, z)], pos, merits(zPIE 의 11 면 TV),
          sharp_dz·sharp(끝난 시편을 ±5 DoF 로 전파한 TV 곡선, 41 점), z_final, dxo_final.

  python studies13.py two --gap 200 --init runs/tl_single.npz --out runs/tl2_gap200.npz
      TwoLayer 에 두 장 모형(층 간격 gap µm), 한 장 복원에서 이어 시작

  python studies13.py score --runs <runs 폴더>
      모든 결과를 채점해 results13.csv(한 장)와 results13_two.csv(두 장)를 그 폴더에 쓴다. 글의 수치는 이 표에서 온다.
      tv_min   끝난 시편이 가장 선명해지는(TV 최소) 전파 거리 [µm], 41 점 곡선의 최소 근처 7 점에 포물선
      probe_R  복원 프로브의 곡률 반경 [mm]. 위상 기울기에 1차식을 맞추고, 화소는 z 35 mm 의 것으로 환산한다
               (프레넬 축척 정리의 예측과 견주려면 가정한 z 가 아니라 고정된 참 화소로 재야 한다)
      bar_z    USAF 4군 막대 12 묶음의 주기로 거꾸로 구한 z 의 평균과 표준오차 [mm] (zpie 데이터만)
      ncc      정답 막대 표적과의 상관 (합성 데이터만)
      err_round1 첫 라운드 끝(50 반복)의 오차. 둘째 라운드는 l2 정규화(매 반복 O·P × 0.99)가 오차를 고르게 끌어올려
               거리마다의 차이를 묻으므로, 거리 비교에는 이 값(또는 --rounds 1 로 돌린 결과)을 쓴다
"""
import argparse
import csv
import glob
import os

import numpy as np

from dz_core import DATA, dof, dxo, load, reconstruct, reconstruct_two, sharpness_curve
from synth import load_synth


def get_data(name):
    if name == "synth":                                  # 완만하게 퍼지는 빔(R ≈ 2z)
        return load_synth(os.path.join(DATA, "synth_zpie.npz"))
    if name == "synthr":                                 # 실데이터에서 복원한 프로브(R ≈ 5 mm)
        return load_synth(os.path.join(DATA, "synth_real.npz"))
    return load(name)


def cmd_run(a):
    d = get_data(a.data)
    if a.curvR is not None:                              # 시작 프로브의 곡률 반경을 바꿔 본다(그림 1: 70 mm vs 7 mm)
        d["curv_R"] = a.curvR * 1e-3
    zpie = dict(step=a.step, sign=a.sign, start=a.start) if a.zpie else None
    opt = dict(power=True, power_init=False, rounds=a.rounds, iters=a.iters)
    P0 = np.load(a.probe0)["P"] if a.probe0 else None     # 시작 프로브를 파일에서(합성의 정답 프로브로 출발하는 대조군)
    O, P, trace, pos, merits = reconstruct(d, a.z0 * 1e-3, seed=a.seed, zpie=zpie, zmin=a.zmin * 1e-3, P0=P0, **opt)
    z = float(trace[-1, 2])
    dzs = np.linspace(-5, 5, 41) * dof(d, z)
    sharp = sharpness_curve(O, d, z, dzs)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, O=O.astype(np.complex64), P=P.astype(np.complex64), trace=trace, pos=pos,
                        merits=merits, sharp_dz=dzs, sharp=sharp, z_final=z, dxo_final=dxo(d, z),
                        z0=a.z0 * 1e-3, seed=a.seed, zpie=bool(a.zpie))
    print("저장:", a.out, f"z {z * 1e3:.4f} mm, 오차 {trace[-1, 1]:.4g}")


def cmd_two(a):
    """TwoLayer 에 두 장 모형을 간격 gap 으로. gap 0 은 한 장과 같은 대조군."""
    d = load(a.data)
    opt = dict(power=True, power_init=False, rounds=a.rounds, iters=a.iters)
    init = np.load(a.init)
    O1, O2, P, trace = reconstruct_two(d, a.gap * 1e-6, seed=a.seed, O0=init["O"], P0=init["P"], **opt)[:4]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, O1=O1.astype(np.complex64), O2=O2.astype(np.complex64), P=P.astype(np.complex64),
                        trace=trace, gap=a.gap * 1e-6, seed=a.seed)
    print("저장:", a.out, f"gap {a.gap:.0f} um, 오차 {trace[-1, 1]:.4g}")


Z_REF = 0.035   # 곡률 환산에 쓰는 참 화소의 기준 거리 (합성 데이터의 정답, 실데이터 기록값)


def tv_min(r):
    s, dz = r["sharp"], r["sharp_dz"] * 1e6
    i = int(np.argmin(s))
    sl = slice(max(i - 3, 0), min(i + 4, len(s)))
    a = np.polyfit(dz[sl], s[sl], 2)
    return float(-a[1] / (2 * a[0])) if a[0] > 0 else float(dz[i])


def probe_R(P, d):
    """프로브 위상의 x·y 기울기를 진폭 가중 1차식으로 맞춰 곡률 반경을 낸다(두 방향 평균의 역수)."""
    A, n = np.abs(P), P.shape[0]
    dx = d["wavelength"] * Z_REF / d["D"]
    x = np.arange(n - 1) + 0.5 - n / 2
    ks = []
    for g, w, X in ((np.angle(P[:, 1:] * np.conj(P[:, :-1])), A[:, 1:] * A[:, :-1], np.broadcast_to(x, (n, n - 1))),
                    (np.angle(P[1:, :] * np.conj(P[:-1, :])), A[1:, :] * A[:-1, :], np.broadcast_to(x[:, None], (n - 1, n)))):
        m = w > 0.2 * w.max()
        ks.append(np.polyfit(X[m], g[m], 1, w=w[m])[0] / 2 / dx ** 2)
    return float(np.pi / (d["wavelength"] * np.mean(ks)))


def cmd_score(a):
    import usaf
    from dz_core import ncc, scan_region
    cache = {}
    rows = []
    for f in sorted(glob.glob(os.path.join(a.runs, "*.npz"))):
        name = os.path.splitext(os.path.basename(f))[0]
        kind = next((k for k in ("zfix_", "zpie_", "synr1_", "synr_", "synrzpie_", "synfix_", "synzpie_") if name.startswith(k)), None)
        if kind is None:
            continue
        data = {"zfix_": "zpie", "zpie_": "zpie", "synr1_": "synthr", "synr_": "synthr", "synrzpie_": "synthr",
                "synfix_": "synth", "synzpie_": "synth"}[kind]
        if data not in cache:
            cache[data] = get_data(data)
        d = cache[data]
        r = np.load(f)
        tr = r["trace"]
        row = dict(name=name, data=data, zpie=int(bool(r["zpie"])), seed=int(r["seed"]), z0=float(r["z0"]) * 1e3,
                   z_final=float(r["z_final"]) * 1e3, err_final=float(tr[-1, 1]), err_round1=float(tr[min(49, len(tr) - 1), 1]),
                   tv_min=tv_min(r), probe_R=probe_R(r["P"].astype(complex), d) * 1e3, bar_z="", bar_sem="", ncc="")
        if data == "zpie":
            ms = [m for m in usaf.measure(r["O"], (270, 610, 230, 600), float(r["dxo_final"]), groups=(4,)) if m["g"] == 4]
            zb = np.array([m["dx"] * d["D"] / d["wavelength"] * 1e3 for m in ms])
            if len(zb):
                row.update(bar_z=float(zb.mean()), bar_sem=float(zb.std(ddof=1) / np.sqrt(len(zb))), bar_n=len(zb))
        else:
            T = d["truth"]
            row["ncc"] = ncc(np.abs(scan_region(r["O"], r["pos"], d["N"])), np.abs(scan_region(T["O"], T["pos"], d["N"])))
        rows.append(row)
        print(name, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items() if k != "name"})
    keys = ["name", "data", "zpie", "seed", "z0", "z_final", "err_final", "err_round1", "tv_min", "probe_R",
            "bar_z", "bar_sem", "bar_n", "ncc"]
    with open(os.path.join(a.runs, "results13.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    two = []
    for f in sorted(glob.glob(os.path.join(a.runs, "tl2_gap*.npz"))):
        r = np.load(f)
        tr = r["trace"]
        two.append(dict(name=os.path.splitext(os.path.basename(f))[0], gap_um=float(r["gap"]) * 1e6, seed=int(r["seed"]),
                        err_final=float(tr[-1, 1]), err_round1=float(tr[min(49, len(tr) - 1), 1]), err_min=float(tr[:, 1].min())))
    if two:
        with open(os.path.join(a.runs, "results13_two.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, list(two[0]))
            w.writeheader()
            w.writerows(sorted(two, key=lambda t: (t["seed"], t["gap_um"])))
    print("저장:", os.path.join(a.runs, "results13.csv"), len(rows), "줄;", "두 장", len(two), "줄")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--data", default="zpie")
    r.add_argument("--z0", type=float, default=35.0)
    r.add_argument("--zmin", type=float, default=28.0)
    r.add_argument("--curvR", type=float, default=None, help="시작 프로브 곡률 반경 [mm]")
    r.add_argument("--probe0", default=None, help="시작 프로브를 담은 npz (키 P)")
    r.add_argument("--zpie", action="store_true")
    r.add_argument("--step", type=float, default=1.0)
    r.add_argument("--sign", type=float, default=1.0)
    r.add_argument("--start", type=int, default=1)
    r.add_argument("--seed", type=int, default=0)
    r.add_argument("--rounds", type=int, default=2)
    r.add_argument("--iters", type=int, default=50)
    r.add_argument("--out", required=True)
    t = sub.add_parser("two")
    t.add_argument("--data", default="twolayer")
    t.add_argument("--gap", type=float, required=True, help="층 간격 [µm]")
    t.add_argument("--init", required=True, help="한 장 복원 결과(run --data twolayer 의 npz)")
    t.add_argument("--seed", type=int, default=0)
    t.add_argument("--rounds", type=int, default=2)
    t.add_argument("--iters", type=int, default=50)
    t.add_argument("--out", required=True)
    sc = sub.add_parser("score")
    sc.add_argument("--runs", default=os.path.join(DATA, "runs"))
    a = ap.parse_args()
    {"run": cmd_run, "two": cmd_two, "score": cmd_score}[a.cmd](a)

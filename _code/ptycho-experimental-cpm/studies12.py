"""12편의 표와 그림에 들어가는 실행과 집계. python3 get_data.py 로 데이터를 받은 뒤 쓴다.
직접 구현 한 번(2 라운드 × 50 반복, 202 장)은 데스크톱 한 코어에서 8–15 분 걸린다. 여러 개를 동시에 돌린다.

  python3 studies12.py run <설정> <시드>       data/runs/own_<설정>_s<시드>.npz
        설정: bright(밝게 시작 + 매 반복 고정, PtyLab 기본), dim(원 프로브 그대로 + 매 반복 고정),
              free(원 프로브 그대로, 고정 없음), bright3(고정값의 3 배로 시작, 고정 없음)
  python3 studies12.py diag <시드> <dim|bright>  data/runs/diag_<dim|bright>_s<시드>.npz (1 라운드 30 반복, 앞쪽 프로브 스냅숏)
  python3 studies12.py score [파일 ...]           data/results.csv — 데이터 오차, 기준 복원(own_dim_s0)과의 |O| 상관
        파일을 주지 않으면 data/runs/*.npz 전부. PtyLab 참고 실행(ptylab/)의 결과도 같은 표에 넣는다.

글의 판정: 성공 = 데이터 오차 70 이하 이면서 상관 0.9 이상. 실제 값은 두 봉우리로 갈렸다(64–66 대 103 이상).
"""
import csv
import glob
import os
import re
import sys

import numpy as np

import cpm_core as C

RUNS = os.path.join(C.DATA, "runs")
RESULTS = os.path.join(C.DATA, "results.csv")
REF = os.path.join(RUNS, "own_dim_s0.npz")
KEEP = (0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 29)


def run(setting, seed):
    d = C.load()
    O, P, trace, pos = C.reconstruct(d, seed=seed, **C.SETTINGS[setting])
    os.makedirs(RUNS, exist_ok=True)
    np.savez_compressed(os.path.join(RUNS, f"own_{setting}_s{seed}.npz"), object=O.astype(np.complex64),
                        probe=P.astype(np.complex64), trace=trace, positions=pos, seed=seed, setting=setting)


def diag(seed, mode):
    """어둡게/밝게 시작의 차이는 첫 훑기뿐이다(둘 다 매 반복 고정). 앞쪽 반복의 프로브와 시편 가운데를 남긴다."""
    d = C.load()
    snaps = {"its": [], "P": [], "O": []}

    def cb(rnd, it, O, P, trace):
        if it in KEEP:
            c = O.shape[0] // 2
            snaps["its"].append(it)
            snaps["P"].append(P.astype(np.complex64).copy())
            snaps["O"].append(O[c - 300:c + 300, c - 300:c + 300].astype(np.complex64).copy())

    O, P, trace, pos = C.reconstruct(d, seed=seed, callback=cb, rounds=1, iters=30,
                                     **C.SETTINGS["dim" if mode == "dim" else "bright"])
    os.makedirs(RUNS, exist_ok=True)
    np.savez_compressed(os.path.join(RUNS, f"diag_{mode}_s{seed}.npz"), its=np.array(snaps["its"]),
                        P=np.array(snaps["P"]), O=np.array(snaps["O"]), trace=trace, object=O.astype(np.complex64),
                        probe=P.astype(np.complex64), positions=pos, seed=seed, mode=mode)


def _final(z):
    """마지막 오차와 Σ|P|². 직접 구현의 기록은 (Σ|P|², 오차), PtyLab 기록은 (Σ|P|², 무게중심, 가장자리, 오차, |O|²)."""
    P = np.squeeze(z["probe"])
    if "trace" in z.files:
        t = np.asarray(z["trace"])
        return float(t[-1, 3] if t.shape[1] == 5 else t[-1, 1]), float(np.sum(np.abs(P) ** 2))
    return float(np.asarray(z["error"])[-1]), float(np.sum(np.abs(P) ** 2))


def _name(path):
    """파일 이름 → (구현, 설정, 시드). own_dim_s3, ptylab_base_s12, ptylab_scale0.1_s4 꼴."""
    m = re.match(r"(own|ptylab|diag)_(.+)_s(\d+)\.npz$", os.path.basename(path))
    return (m.group(1), m.group(2), int(m.group(3))) if m else (os.path.basename(path), "", -1)


def score(files):
    d = C.load()
    N = d["N"]
    mp2 = float(np.max(d["I"].sum((-1, -2))))
    pos_default, _ = C.geometry(d)
    zr = np.load(REF)
    ref = np.abs(C.scan_crop(np.squeeze(zr["object"]), zr["positions"], N))
    rows = []
    for f in sorted(files):
        z = np.load(f, allow_pickle=True)
        pos = z["positions"] if "positions" in z.files else None
        O = np.squeeze(z["object"])
        if pos is None:      # PtyLab 실행: 위치 = round(encoder/dxo) + No//2 − N//2, No 는 시편 배열의 크기
            pix = np.round(d["encoder"] / d["dxo"]).astype(int)
            pos = pix + O.shape[-1] // 2 - N // 2
        err, psum = _final(z)
        impl, setting, seed = _name(f)
        s = C.ncc(ref, np.abs(C.scan_crop(O, pos, N)))
        rows.append(dict(impl=impl, setting=setting, seed=seed, error=round(err, 3), ncc=round(s, 4),
                         power_ratio=round(psum / mp2, 4), success=int(err < 70 and s > 0.9)))
        print(f"{os.path.basename(f):32s} 오차 {err:8.2f}  상관 {s:.3f}  Σ|P|²/고정값 {psum / mp2:8.3g}", flush=True)
    with open(RESULTS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("저장:", RESULTS)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "run":
        run(sys.argv[2], int(sys.argv[3]))
    elif cmd == "diag":
        diag(int(sys.argv[2]), sys.argv[3])
    elif cmd == "score":
        score(sys.argv[2:] or glob.glob(os.path.join(RUNS, "*.npz")))

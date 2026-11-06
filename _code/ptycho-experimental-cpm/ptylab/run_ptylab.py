"""12편의 참고 실행: PtyLab(Loetgering 등 2023, pip install ptylab, 0.3.7 에서 확인)의 예제 설정으로 TwoLayer 를 복원한다.
본문 2–4절의 PtyLab 결과(시드 20개, 소거 행렬, 출발 세기 훑기)가 이 스크립트에서 나왔다. 한 번에 10분 안팎(CPU).

  python run_ptylab.py <변형> <시드> [라운드 = 2]        → data/runs/ptylab_<변형>_s<시드>.npz
변형은 + 로 잇는다(파일 이름에서는 _ 로 바뀐다).
  base        예제 설정 그대로 (mPIE 50 반복 × 라운드, 짝수 라운드에 l2, 세기 맞춤 켬)
  nosmooth    프로브 진폭 매끄럽게 하기 끔
  nocom       무게중심 맞춤 끔
  nol2        l2 정규화 끔
  flatprobe   초기 프로브에 곡률 위상을 넣지 않음
  sequential  패턴 순서를 늘 차례대로
  epie        엔진을 ePIE 로 (모멘텀 없음)
  nopower     세기 맞춤 끔 (출발도 원 프로브 그대로, 고정값의 0.00093 배)
  initonly    첫 라운드 시작 때만 고정값에 맞춤, 매 반복 고정 없음
  iteronly    매 반복 고정만, 출발은 원 프로브 그대로
  scale<f>    첫 라운드 시작 때 고정값의 f 배로 맞춤, 고정 없음 (예: scale0.1)
반복마다 (Σ|P|², 프로브 무게중심 반지름, 조리개 가장자리 고리의 세기 비율, 데이터 오차, 주사 영역 평균 |O|²)를 남긴다.
"""
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import numpy as np
import PtyLab
from PtyLab import Engines

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from cpm_core import DATA, HDF5  # noqa: E402

var, seed = sys.argv[1], int(sys.argv[2])
rounds = int(sys.argv[3]) if len(sys.argv) > 3 else 2
parts = set(var.split("+"))
scale = next((float(p[5:]) for p in parts if p.startswith("scale")), None)
np.random.seed(seed)        # PtyLab 은 전역 난수로 순서를 섞고 모멘텀을 건다

data, rec, params, monitor, _ = PtyLab.easyInitialize(HDF5, operationMode="CPM", dummyMonitor=True)
params.TV_autofocus = False
params.l2reg = False
params.l2reg_probe_aleph = 1e-2
params.l2reg_object_aleph = 1e-2
rec.No = int(rec.No * 0.9)
rec.npsm = rec.nosm = rec.nlambda = rec.nslice = 1
rec.initialProbe = "circ"
rec.initialObject = "ones"
rec.initializeObjectProbe()
if "flatprobe" not in parts:
    rec.probe *= np.exp(1.0j * 2 * np.pi / (rec.wavelength * rec.zo * 2) * (rec.Xp ** 2 + rec.Yp ** 2) / 2)
probe0 = np.squeeze(np.asarray(rec.probe)).copy()
params.positionOrder = "sequential" if "sequential" in parts else "random"
params.propagatorType = "Fraunhofer"
params.positionCorrectionSwitch = False
params.modulusEnforcedProbeSwitch = False
params.probeSmoothenessSwitch = "nosmooth" not in parts
params.probeSmoothnessAleph = 1e-2
params.probeSmoothenessWidth = 10
params.comStabilizationSwitch = 10 ** 9 if "nocom" in parts else 10     # 0 을 주면 PtyLab 이 0 으로 나눈다
params.orthogonalizationSwitch = True
params.orthogonalizationFrequency = 10
params.absorbingProbeBoundary = False
params.objectContrastSwitch = False
params.absObjectSwitch = False
params.backgroundModeSwitch = False
params.couplingSwitch = True
params.couplingAleph = 1
params.probePowerCorrectionSwitch = not (parts & {"nopower", "initonly"} or scale is not None)
params.gpuSwitch = False
params.fftshiftSwitch = False
params.intensityConstraint = "standard"

eng = (Engines.ePIE if "epie" in parts else Engines.mPIE)(rec, data, params, monitor)
eng.numIterations = 50
if "iteronly" in parts:
    eng._initialProbePowerCorrection = lambda: None      # 라운드 시작 때의 맞춤만 막는다

mp2 = float(data.maxProbePower) ** 2                       # 고정값 = 가장 밝은 패턴의 총세기
Np = rec.Np
yy, xx = np.indices((Np, Np)) - Np // 2
rr = np.hypot(yy, xx)
r0 = np.sqrt(np.count_nonzero(np.abs(probe0) > 0) / np.pi)
rim = (rr > 0.8 * r0) & (rr < 1.2 * r0)
scan = np.zeros(np.squeeze(np.asarray(rec.object)).shape[-2:], bool)
for r, c in np.asarray(rec.positions):
    scan[r + Np // 4:r + 3 * Np // 4, c + Np // 4:c + 3 * Np // 4] = True
trace = []


def record(loop, reconstruction):
    P = np.squeeze(np.asarray(reconstruction.probe))
    O = np.squeeze(np.asarray(reconstruction.object))
    I = np.abs(P) ** 2
    pw = float(I.sum())
    err = np.asarray(reconstruction.error)
    trace.append((pw, float(np.hypot(np.sum(I * yy), np.sum(I * xx)) / pw), float(I[rim].sum() / pw),
                  float(err[-1]) if err.size else np.nan, float(np.mean(np.abs(O[scan]) ** 2))))


t0 = time.time()
for i in range(1, rounds + 1):
    params.l2reg = ("nol2" not in parts) and i % 2 == 0
    if i == 1 and ("initonly" in parts or scale is not None):
        f = 1.0 if scale is None else scale
        rec.probe = rec.probe / np.sqrt(np.sum(np.abs(rec.probe) ** 2)) * np.sqrt(f * mp2)
    eng.reconstruct(vis_after_each_iteration=record)
    pw, com, rimf, err, om = trace[-1]
    print(f"round {i}: error {err:.4g}, sum|P|^2/mp2 {pw / mp2:.3g}, rim {rimf:.2f}, {time.time() - t0:.0f}s", flush=True)
os.makedirs(os.path.join(DATA, "runs"), exist_ok=True)
np.savez_compressed(os.path.join(DATA, "runs", f"ptylab_{var.replace('+', '_')}_s{seed}.npz"),
                    object=np.squeeze(np.asarray(rec.object)), probe=np.squeeze(np.asarray(rec.probe)), probe0=probe0,
                    positions=np.asarray(rec.positions), trace=np.array(trace), mp2=mp2, var=var, seed=seed)
print("done", flush=True)

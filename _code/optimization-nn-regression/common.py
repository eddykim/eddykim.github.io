"""실험 공통 설정 — 최적화 5편.

측정 데이터 규약은 1~4편과 동일하다(generate_figures.py 에서 그대로 옮김).

numpy 2.2.6 + macOS Accelerate BLAS 조합에서 정상적인 행렬곱에도
'divide by zero / overflow / invalid value encountered in matmul' 경고가 뜬다.
완전히 유한한 난수 행렬곱에서도 재현되고 결과·가중치는 모두 유한함을 확인했으므로
허위 경고로 판단해 억제한다(수치 문제를 감추는 것이 아니다).
"""
import time
import warnings

import numpy as np

warnings.filterwarnings("ignore", message=".*encountered in matmul.*")

from mlp import MLP, Scaler
from train_data import make_dataset, make_spectra, augment, measured_reference
from levenberg_marquardt import levenberg_marquardt, objective, residual
from reflectance_model import reflectance, N_SIO2

# --- 1~4편과 동일한 규약 -------------------------------------------------
WL = np.linspace(450, 750, 300)
TRUE_D = 1490.0
NOISE_STD = 0.004
D_RANGE = (1200.0, 1800.0)
HARD_D0 = 1300.0          # 3편에서 LM 이 1293nm 에 갇히던 시작점
SUCCESS_TOL = 5.0         # 4편 관례: |d - 1490| < 5


def success(d):
    return abs(float(d) - TRUE_D) < SUCCESS_TOL


def lm_fit(d0, wl, R, n_iter=30):
    """LM 한 번 실행 → (최종 d, 반복 횟수, 최종 J)"""
    d_hist, J_hist, _ = levenberg_marquardt(d0, wl, R, n_iter=n_iter)
    return float(d_hist[-1]), len(d_hist) - 1, float(J_hist[-1])


def train_net(X, d, scaler, layers, rng, epochs, lr=1e-3, batch_size=64,
              X_val=None, d_val=None, verbose=False):
    """스케일링 → MLP 학습. 반환: (model, history, 학습 시간[s])"""
    Xs = scaler.x(X)
    ys = scaler.y(d)
    Xv = scaler.x(X_val) if X_val is not None else None
    yv = scaler.y(d_val) if d_val is not None else None
    net = MLP(layers, rng=rng)
    t0 = time.perf_counter()
    hist = net.fit(Xs, ys, Xv, yv, epochs=epochs, batch_size=batch_size,
                   lr=lr, rng=rng, verbose=verbose)
    return net, hist, time.perf_counter() - t0


def predict_d(net, scaler, X):
    """스펙트럼(들) → 두께(nm) 배열"""
    X = np.atleast_2d(X)
    return scaler.y_inv(net.predict(scaler.x(X)))


def rmse_nm(net, scaler, X, d_true):
    return float(np.sqrt(np.mean((predict_d(net, scaler, X) - d_true) ** 2)))

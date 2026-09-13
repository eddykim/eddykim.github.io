"""numpy 만으로 구현한 소형 MLP — 최적화 5편 실험용.

이 환경에는 torch/tensorflow/sklearn 이 없다. 1~4편이 LM 을 Madsen·Nielsen·Tingleff
2004 의 Algorithm 3.16 으로 직접 옮긴 것과 같은 방식으로, 여기서도 순전파·역전파·Adam 을
원 논문 수식 그대로 numpy 로 옮긴다.

- 역전파:     Rumelhart, Hinton & Williams 1986, Nature 323, 533-536
- 가중치 초기화: Glorot & Bengio 2010, AISTATS (PMLR 9:249-256) 의 normalized initialization
                W ~ U[-sqrt(6/(fan_in+fan_out)), +sqrt(6/(fan_in+fan_out))]
- 옵티마이저:  Kingma & Ba 2015, arXiv:1412.6980 의 Algorithm 1

활성함수는 tanh 를 쓴다. Glorot 초기화가 상정한 대칭 포화 함수이고, 회귀 출력의
부호를 살리기에도 편하다. 출력층은 선형이다.
"""
import numpy as np


class MLP:
    """전결합 MLP. 은닉층은 tanh, 출력층은 선형. 손실은 MSE.

    layer_sizes: 예) [300, 256, 1] 이면 입력 300 → 은닉 256 → 출력 1
    """

    def __init__(self, layer_sizes, rng=None):
        self.layer_sizes = list(layer_sizes)
        rng = np.random.default_rng(0) if rng is None else rng

        self.W, self.b = [], []
        for fan_in, fan_out in zip(self.layer_sizes[:-1], self.layer_sizes[1:]):
            # Glorot & Bengio 2010 식 (16)
            limit = np.sqrt(6.0 / (fan_in + fan_out))
            self.W.append(rng.uniform(-limit, limit, size=(fan_in, fan_out)))
            self.b.append(np.zeros(fan_out))

        # Adam 의 1차·2차 모멘트 (Kingma & Ba Algorithm 1 의 m, v)
        self._mW = [np.zeros_like(w) for w in self.W]
        self._vW = [np.zeros_like(w) for w in self.W]
        self._mb = [np.zeros_like(b) for b in self.b]
        self._vb = [np.zeros_like(b) for b in self.b]
        self._t = 0

    # ------------------------------------------------------------------
    def forward(self, X, cache=False):
        """X: (N, n_in) → (N, n_out). cache=True 면 역전파용 중간값을 함께 반환."""
        a = X
        acts = [a]
        for i, (w, b) in enumerate(zip(self.W, self.b)):
            z = a @ w + b
            a = z if i == len(self.W) - 1 else np.tanh(z)   # 마지막 층만 선형
            acts.append(a)
        return (a, acts) if cache else a

    def _backward(self, acts, y_true):
        """MSE 손실의 기울기. 반환: dW, db (리스트)."""
        n = y_true.shape[0]
        dW = [None] * len(self.W)
        db = [None] * len(self.b)

        # dL/dz_out : L = mean((y - t)^2) 이므로 2/n 배가 붙는다
        delta = 2.0 * (acts[-1] - y_true) / n

        for i in reversed(range(len(self.W))):
            dW[i] = acts[i].T @ delta
            db[i] = delta.sum(axis=0)
            if i > 0:
                # tanh' = 1 - tanh^2, acts[i] 가 이미 tanh 출력이다
                delta = (delta @ self.W[i].T) * (1.0 - acts[i] ** 2)
        return dW, db

    def _adam_step(self, dW, db, lr, b1=0.9, b2=0.999, eps=1e-8):
        """Kingma & Ba 2015 Algorithm 1 그대로."""
        self._t += 1
        bc1 = 1.0 - b1 ** self._t     # bias correction 분모
        bc2 = 1.0 - b2 ** self._t
        for i in range(len(self.W)):
            self._mW[i] = b1 * self._mW[i] + (1 - b1) * dW[i]
            self._vW[i] = b2 * self._vW[i] + (1 - b2) * dW[i] ** 2
            self.W[i] -= lr * (self._mW[i] / bc1) / (np.sqrt(self._vW[i] / bc2) + eps)

            self._mb[i] = b1 * self._mb[i] + (1 - b1) * db[i]
            self._vb[i] = b2 * self._vb[i] + (1 - b2) * db[i] ** 2
            self.b[i] -= lr * (self._mb[i] / bc1) / (np.sqrt(self._vb[i] / bc2) + eps)

    # ------------------------------------------------------------------
    def _snapshot(self):
        return ([w.copy() for w in self.W], [b.copy() for b in self.b])

    def _restore(self, snap):
        self.W = [w.copy() for w in snap[0]]
        self.b = [b.copy() for b in snap[1]]

    def fit(self, X, y, X_val=None, y_val=None, epochs=300, batch_size=64,
            lr=1e-3, rng=None, verbose=False, restore_best=True):
        """미니배치 Adam 학습. 반환: history dict (epoch 별 train/val MSE).

        restore_best=True 면 검증 오차가 최소였던 epoch 의 가중치를 되돌린다.
        Fried & Masa 1994 가 권고한 절차다 — 훈련 오차는 계속 줄어도 시험 오차는
        어느 시점부터 증가하므로, 학습 중 '낯선' 데이터로 주기적으로 시험해
        최적 지점을 찾아야 한다.
        """
        rng = np.random.default_rng(0) if rng is None else rng
        n = X.shape[0]
        hist = {"epoch": [], "train": [], "val": []}
        best_val, best_snap, best_ep = np.inf, None, -1

        for ep in range(epochs):
            idx = rng.permutation(n)
            for s in range(0, n, batch_size):
                bi = idx[s:s + batch_size]
                _, acts = self.forward(X[bi], cache=True)
                dW, db = self._backward(acts, y[bi])
                self._adam_step(dW, db, lr)

            tr = float(np.mean((self.forward(X) - y) ** 2))
            va = float(np.mean((self.forward(X_val) - y_val) ** 2)) if X_val is not None else np.nan
            hist["epoch"].append(ep)
            hist["train"].append(tr)
            hist["val"].append(va)
            if X_val is not None and va < best_val:
                best_val, best_snap, best_ep = va, self._snapshot(), ep
            if verbose and (ep % max(1, epochs // 10) == 0 or ep == epochs - 1):
                print(f"    epoch {ep:4d}  train {tr:.3e}  val {va:.3e}")

        if restore_best and best_snap is not None:
            self._restore(best_snap)
        hist["best_epoch"] = best_ep
        hist["best_val"] = float(best_val)
        return hist

    def predict(self, X):
        return self.forward(X)


class Scaler:
    """입력은 파장점별 표준화, 출력 두께는 [-1,1] 선형 정규화.

    정규화를 하지 않으면 학습이 되지 않는다(입력 R 은 0.05~0.4, 출력 d 는 1200~1800).
    """

    def __init__(self, d_lo, d_hi):
        self.d_lo, self.d_hi = float(d_lo), float(d_hi)
        self.mu = None
        self.sigma = None

    def fit_inputs(self, X):
        self.mu = X.mean(axis=0)
        self.sigma = X.std(axis=0)
        self.sigma[self.sigma < 1e-12] = 1.0
        return self

    def x(self, X):
        return (X - self.mu) / self.sigma

    def y(self, d):
        """두께 → [-1, 1]"""
        mid = 0.5 * (self.d_lo + self.d_hi)
        half = 0.5 * (self.d_hi - self.d_lo)
        return ((np.asarray(d, dtype=float) - mid) / half).reshape(-1, 1)

    def y_inv(self, yn):
        mid = 0.5 * (self.d_lo + self.d_hi)
        half = 0.5 * (self.d_hi - self.d_lo)
        return np.asarray(yn, dtype=float).ravel() * half + mid

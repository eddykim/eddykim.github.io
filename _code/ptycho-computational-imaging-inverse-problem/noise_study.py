"""상자 번짐과 부호 번짐을 잡음 조건별로 비교한다 — 1편 7절.

부호 셔터는 빛을 절반만 받는다. 그 손해와 전달함수 영점이 사라지는 이득 중 어느 쪽이
큰지는 잡음이 얼마나 크냐에 달렸다. 방법마다 정답과의 오차가 최소가 되는 lam 을 골라,
각 방법이 낼 수 있는 가장 좋은 결과끼리 비교한다.

  gaussian(): 신호와 무관한 가우시안 잡음 (읽기 잡음이 지배할 때)
  poisson():  광자 수에 비례하는 샷 잡음 (빛이 약할 때)
"""
import numpy as np

from forward_models import (box_kernel, coded_kernel, smear_matrix, test_signal,
                            measure, measure_poisson)
from regularize import (least_squares, tikhonov, tv_admm, difference_matrix, rel_error, best_lambda)

X = test_signal()
N = X.size
A_BOX = smear_matrix(box_kernel(), N)
A_CODE = smear_matrix(coded_kernel(), N)
D = difference_matrix(N)

METHODS = {
    "ls": (lambda A, y, l: least_squares(A, y), [0.0]),
    "tik": (lambda A, y, l: tikhonov(A, y, l), np.logspace(-5, 0, 30)),
    "tikd": (lambda A, y, l: tikhonov(A, y, l, D), np.logspace(-4, 1, 30)),
    "tv": (lambda A, y, l: tv_admm(A, y, l, n_iter=400), np.logspace(-5, -1, 14)),
}


def compare(make_y, levels, trials=5, seed=1, methods=("ls", "tik", "tikd", "tv")):
    """levels 마다 {방법: (상자 오차 중앙값, 부호 오차 중앙값)}."""
    rng = np.random.default_rng(seed)
    out = {}
    for lv in levels:
        row = {}
        for m in methods:
            solver, lams = METHODS[m]
            med = []
            for A in (A_BOX, A_CODE):
                errs = [best_lambda(solver, A, make_y(A, lv, rng), X, lams)[1] for _ in range(trials)]
                med.append(float(np.median(errs)))
            row[m] = tuple(med)
        out[lv] = row
    return out


def gaussian(levels=(1e-3, 3e-3, 1e-2, 3e-2), **kw):
    return compare(lambda A, s, rng: measure(A, X, s, rng), levels, **kw)


def poisson(levels=(10, 100, 1000, 10000), **kw):
    return compare(lambda A, p, rng: measure_poisson(A, X, p, rng), levels, **kw)


if __name__ == "__main__":
    for name, res in (("가우시안 sigma", gaussian()), ("샷 잡음, 칩당 광자", poisson())):
        print(f"\n{name} | " + " | ".join(f"{m} 상자/부호" for m in METHODS))
        for lv, row in res.items():
            print(f"{lv:g} | " + " | ".join(f"{b:.3f}/{c:.3f}" for b, c in row.values()))

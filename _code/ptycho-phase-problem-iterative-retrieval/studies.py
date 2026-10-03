"""2편의 표와 그림에 들어가는 수치 실험.

모든 실험은 무작위 위상에서 출발한 여러 번의 시도를 세고, 자명한 모호성을 맞춘 오차가
0.02 아래면 성공으로 본다 (phase_retrieval.align_error).
"""
import numpy as np

from phase_retrieval import make_object, run, align_error, F, poisson_magnitude

N = 128
SUCCESS = 0.02


def trials(obj, sup, alg, n_iter, n_trials=8, real=False, beta=0.9, er_tail=0, seed0=0, mag=None):
    mag = np.abs(F(obj)) if mag is None else mag
    errs = [align_error(run(alg, mag, sup, n_iter, seed0 + r, real, beta, er_tail=er_tail)[0], obj)
            for r in range(n_trials)]
    return np.array(errs)


def algorithm_table(n_iter=500):
    """물체 종류별 ER, HIO, RAAR 성공 횟수와 오차 중앙값 (지지 48, sigma = 7.1)."""
    out = {}
    for kind in ("real", "complex"):
        obj, sup = make_object(N, 48, kind)
        for alg, beta in (("ER", None), ("HIO", 0.9), ("RAAR", 0.9), ("RAAR", 0.99)):
            e = trials(obj, sup, alg, n_iter, real=kind == "real", beta=beta or 0.9)
            out[(kind, alg, beta)] = (int((e < SUCCESS).sum()), float(np.median(e)))
    return out


def raar_beta(betas=(0.75, 0.85, 0.9, 0.95, 0.99), n_iter=1000):
    obj, sup = make_object(N, 48, "complex")
    return {b: int((trials(obj, sup, "RAAR", n_iter, beta=b) < SUCCESS).sum()) for b in betas}


def oversampling(sizes=(88, 80, 72, 64), iters=(600, 3000)):
    """지지 크기 s 마다 sigma = N^2 / s^2 와 HIO(+마지막 50회 ER) 성공 횟수."""
    out = {}
    for kind in ("real", "complex"):
        for s in sizes:
            obj, sup = make_object(N, s, kind)
            for it in iters:
                e = trials(obj, sup, "HIO", it, real=kind == "real", er_tail=50, seed0=100)
                out[(kind, s, it)] = int((e < SUCCESS).sum())
    return out


def noise(photons=(1e9, 1e8, 1e7, 1e6, 1e5), n_trials=6):
    """총 광자 수별 오차 중앙값과 가장 좋은 시도의 복원 (복소 물체)."""
    obj, sup = make_object(N, 48, "complex")
    out = {}
    for ph in photons:
        mag = poisson_magnitude(obj, ph, rng=7)
        recs = [run("HIO", mag, sup, 600, r, er_tail=50)[0] for r in range(n_trials)]
        errs = np.array([align_error(r, obj) for r in recs])
        out[ph] = (float(np.median(errs)), recs[int(np.argmin(errs))], float(errs.min()))
    return out


if __name__ == "__main__":
    import time
    t = time.time()
    for k, v in algorithm_table().items():
        print("알고리즘", k, v)
    print("RAAR beta", raar_beta())
    for k, v in oversampling().items():
        print("오버샘플링", k, f"sigma={N * N / k[1] ** 2:.2f}", v)
    for k, v in noise().items():
        print("노이즈", f"{k:.0e}", round(v[0], 3), round(v[2], 3))
    print(f"{time.time() - t:.0f} s")

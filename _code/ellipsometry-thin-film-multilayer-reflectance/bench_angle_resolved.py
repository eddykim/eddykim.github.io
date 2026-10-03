"""입사각×파장 격자에서 세 방법의 계산 속도 비교 (본문 6절의 "로아드가 SMM보다 4~5배 빠르다").

로아드 재귀도 (입사각, 파장) 격자의 각 점에서 독립적인 스칼라 연산이라,
rouard.py 의 재귀를 numpy 배열에 그대로 적용하면(반복문은 층 수만큼) 벡터화된다.
TMM 도 같은 방식으로 배치 행렬곱으로 벡터화해 함께 잰다. SMM 은 smm_tensor.py 그대로다.

실행: python bench_angle_resolved.py
세 방법의 결과가 1e-10 안에서 일치하는지 확인한 뒤, 7회 측정 중 최솟값을 출력한다.
절대 시간은 기계마다 다르고, 비율을 보는 것이 목적이다.
"""
import time

import numpy as np

from rouard import _fresnel_r
from smm_tensor import _cos_theta_grid, smm_reflectance_tensor


def rouard_tensor(n_list, d_list, theta0, wavelength, pol="s"):
    """rouard.py 의 재귀를 (n_angle, n_wav) 배열에 그대로 적용한 버전."""
    c = _cos_theta_grid(n_list, theta0)              # 매질별 (n_angle, 1)
    M = len(d_list)
    rho = _fresnel_r(n_list[M], c[M], n_list[M + 1], c[M + 1], pol)
    for m in range(M, 0, -1):
        beta = 2 * np.pi / wavelength[None, :] * n_list[m] * d_list[m - 1] * c[m]
        r_top = _fresnel_r(n_list[m - 1], c[m - 1], n_list[m], c[m], pol)
        phase = np.exp(-2j * beta)
        rho = (r_top + rho * phase) / (1 + r_top * rho * phase)
    return rho


def tmm_tensor(n_list, d_list, theta0, wavelength, pol="s"):
    """tmm.py 의 특성행렬 곱을 (n_angle, n_wav) 배치 행렬곱으로 옮긴 버전."""
    c = _cos_theta_grid(n_list, theta0)
    eta = [n * ci if pol == "s" else n / ci for n, ci in zip(n_list, c)]
    shape = (len(theta0), len(wavelength))

    def boundary(e):
        e = np.broadcast_to(e, shape)
        D = np.empty(shape + (2, 2), dtype=complex)
        D[..., 0, 0], D[..., 0, 1] = 1, 1
        D[..., 1, 0], D[..., 1, 1] = e, -e
        return D

    Mx = np.linalg.inv(boundary(eta[0]))
    for j in range(1, len(d_list) + 1):
        beta = 2 * np.pi / wavelength[None, :] * n_list[j] * d_list[j - 1] * c[j]
        e = np.broadcast_to(eta[j], shape)
        Q = np.empty(shape + (2, 2), dtype=complex)
        Q[..., 0, 0], Q[..., 0, 1] = np.cos(beta), 1j / e * np.sin(beta)
        Q[..., 1, 0], Q[..., 1, 1] = 1j * e * np.sin(beta), np.cos(beta)
        Mx = Mx @ Q
    Mx = Mx @ boundary(eta[-1])
    r = Mx[..., 1, 0] / Mx[..., 0, 0]
    return -r if pol == "p" else r   # 어드미턴스 규약 → 1편 프레넬 규약


def best_of(f, n=7):
    times = []
    for _ in range(n):
        t0 = time.perf_counter()
        f()
        times.append(time.perf_counter() - t0)
    return min(times)


if __name__ == "__main__":
    n_si = 3.88 - 0.02j
    stacks = {
        "1층": ([1.0, 1.46, n_si], [300.0]),
        "3층": ([1.0, 1.46, 2.35, 2.02, n_si], [120.0, 80.0, 45.0]),
        "20층": ([1.0] + [1.46, 2.35] * 10 + [n_si], [100.0, 60.0] * 10),
    }
    theta = np.radians(np.linspace(0.1, 85, 300))   # 그림4와 같은 격자
    wl = np.linspace(400, 1000, 300)

    print("== 입사각 300 x 파장 300 격자, p-편광 ==")
    for name, (n, d) in stacks.items():
        a = rouard_tensor(n, d, theta, wl, "p")
        b = smm_reflectance_tensor(n, d, theta, wl, "p")
        c = tmm_tensor(n, d, theta, wl, "p")
        assert np.abs(a - b).max() < 1e-10 and np.abs(c - b).max() < 1e-10
        tr = best_of(lambda: rouard_tensor(n, d, theta, wl, "p"))
        tt = best_of(lambda: tmm_tensor(n, d, theta, wl, "p"))
        ts = best_of(lambda: smm_reflectance_tensor(n, d, theta, wl, "p"))
        print(f"  {name:>4}: 로아드 {tr * 1e3:7.2f} ms | TMM {tt * 1e3:7.2f} ms | "
              f"SMM {ts * 1e3:7.2f} ms | SMM/로아드 = {ts / tr:.1f}배")

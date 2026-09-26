"""원 피팅 — 이미지처리 5편.

엣지 점 수백 개에서 중심과 반경을 뽑는 방법은 무엇을 최소화하느냐로 갈린다.

  기하학적 거리  점에서 원까지의 실제 거리 sqrt((x-a)^2+(y-b)^2) - R 의 제곱합.
                 재고 싶은 것이 바로 이것이지만 파라미터에 비선형이라 반복법이 필요하다.
  대수적 거리    x^2+y^2+Dx+Ey+F 처럼 원의 방정식에 좌표를 넣은 값. 파라미터에
                 선형이라 연립방정식 한 번으로 풀리지만, 재는 것이 거리가 아니다.

대수적 거리는 반경이 큰 쪽에서 더 큰 값을 내므로, 점들이 원 전체에 고루 퍼져
있지 않으면 해가 한쪽으로 끌린다. 짧은 원호에서 Kasa 법이 반경을 작게 잡는 편향이
여기서 나온다. Pratt 와 Taubin 은 같은 대수적 거리를 쓰되 제약을 바꿔 이 편향을
줄인 것이고, 기하학적 피팅은 아예 거리를 재러 간다.
"""
import numpy as np
import scipy.linalg
from scipy.optimize import least_squares


def kasa(x, y):
    """가장 단순한 대수적 피팅. x^2+y^2+Dx+Ey+F = 0 을 최소자승으로 푼다."""
    A = np.column_stack([x, y, np.ones_like(x)])
    b = -(x ** 2 + y ** 2)
    D, E, F = np.linalg.lstsq(A, b, rcond=None)[0]
    a, c = -D / 2.0, -E / 2.0
    R = np.sqrt(max(a ** 2 + c ** 2 - F, 0.0))
    return a, c, R


def pratt(x, y):
    """Pratt 의 제약 B^2+C^2-4AD = 1 을 쓴 대수적 피팅.

    A(x^2+y^2)+Bx+Cy+D = 0 에서 A 를 1 로 고정하지 않으므로, 점들이 직선에
    가까워져도(A -> 0) 해가 발산하지 않는다. 제약이 붙은 최소자승이라
    일반화 고유값 문제가 되고, 음이 아닌 가장 작은 고유값의 고유벡터가 해다.
    """
    xm, ym = x.mean(), y.mean()
    u, v = x - xm, y - ym                      # 무게중심으로 옮겨 수치 안정화
    M = np.column_stack([u ** 2 + v ** 2, u, v, np.ones_like(u)])
    S = M.T @ M
    N = np.array([[0.0, 0.0, 0.0, -2.0],
                  [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0],
                  [-2.0, 0.0, 0.0, 0.0]])
    try:
        evals, evecs = scipy.linalg.eig(S, N)
    except Exception:
        return kasa(x, y)
    ok = np.isfinite(evals.real) & (np.abs(evals.imag) < 1e-8) & (evals.real > -1e-9)
    if not ok.any():
        return kasa(x, y)
    A, B, C, D = evecs[:, np.flatnonzero(ok)[np.argmin(evals.real[ok])]].real
    if abs(A) < 1e-12:
        return kasa(x, y)
    disc = B ** 2 + C ** 2 - 4.0 * A * D
    return (-B / (2 * A) + xm, -C / (2 * A) + ym,
            np.sqrt(max(disc, 0.0)) / (2 * abs(A)))


def taubin(x, y):
    """Taubin 의 피팅. 대수적 거리를 그 기울기 크기로 정규화해 편향을 줄인다.

    Chernov 가 정리한 표준 구현을 그대로 따른다. 좌표를 무게중심으로 옮긴 뒤
    특성 4차 다항식의 근을 뉴턴법으로 찾고, 그 근으로 중심과 반경을 계산한다.
    """
    xm, ym = x.mean(), y.mean()
    u, v = x - xm, y - ym
    z = u ** 2 + v ** 2

    Mxx, Myy, Mxy = (u * u).mean(), (v * v).mean(), (u * v).mean()
    Mxz, Myz, Mzz = (u * z).mean(), (v * z).mean(), (z * z).mean()
    Mz = Mxx + Myy
    cov_xy = Mxx * Myy - Mxy * Mxy
    var_z = Mzz - Mz * Mz

    a3 = 4.0 * Mz
    a2 = -3.0 * Mz * Mz - Mzz
    a1 = var_z * Mz + 4.0 * cov_xy * Mz - Mxz * Mxz - Myz * Myz
    a0 = (Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy)
          - var_z * cov_xy)
    a22, a33 = a2 + a2, a3 + a3 + a3

    root, f = 0.0, a0
    for _ in range(99):
        df = a1 + root * (a22 + a33 * root)
        if df == 0.0:
            break
        step = f / df
        new = root - step
        if new == root or not np.isfinite(new):
            break
        f_new = a0 + new * (a1 + new * (a2 + new * a3))
        if abs(f_new) >= abs(f) and root != 0.0:
            break
        root, f = new, f_new
        if abs(step) < 1e-14:
            break

    det = root * root - root * Mz + cov_xy
    if abs(det) < 1e-14:
        return kasa(x, y)
    cu = (Mxz * (Myy - root) - Myz * Mxy) / det / 2.0
    cv = (Myz * (Mxx - root) - Mxz * Mxy) / det / 2.0
    return cu + xm, cv + ym, np.sqrt(max(cu * cu + cv * cv + Mz, 0.0))


def geometric(x, y, init=None):
    """기하학적 거리의 제곱합을 Levenberg-Marquardt 로 최소화한다.

    최적화 시리즈 3편에서 다룬 그 LM 이다. 초기값은 대수적 피팅으로 잡는다.
    비선형이지만 잔차가 작고 야코비안이 잘 정의되어 있어 몇 번 만에 수렴한다.
    """
    a0, b0, R0 = init if init is not None else kasa(x, y)

    def residual(p):
        a, b, R = p
        return np.hypot(x - a, y - b) - R

    sol = least_squares(residual, [a0, b0, R0], method="lm")
    return sol.x[0], sol.x[1], sol.x[2]


FITTERS = {"kasa": kasa, "pratt": pratt, "taubin": taubin, "geometric": geometric}


def arc_points(n, cx, cy, R, span_deg, start_deg=0.0, noise=0.0, rng=None):
    """원호 위에 점을 고르게 놓고 반경 방향으로 노이즈를 준다.

    노이즈를 반경 방향으로만 주는 이유는 4편에서 얻은 엣지 점의 오차가 엣지에
    수직인 방향으로 생기기 때문이다. 접선 방향 오차는 원 피팅에 거의 영향이 없다.
    """
    th = np.deg2rad(start_deg + np.linspace(0.0, span_deg, n, endpoint=False))
    r = np.full(n, float(R))
    if noise > 0.0:
        r = r + noise * (rng or np.random.default_rng(0)).standard_normal(n)
    return cx + r * np.cos(th), cy + r * np.sin(th)

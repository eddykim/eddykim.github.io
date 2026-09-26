"""타원 피팅 — 이미지처리 5편.

원이 다섯 개가 아니라 세 개의 파라미터를 갖는 특수한 경우라면, 일반적인 이차곡선은
여섯 개(스케일을 빼면 다섯 개)를 갖는다. 문제는 이차곡선의 대수적 피팅이 타원만
내놓지 않는다는 데 있다. 점들이 조금만 흐트러져도 포물선이나 쌍곡선이 나온다.

Fitzgibbon 등이 1999년에 제안한 방법은 제약 4AC - B^2 = 1 을 걸어 이 문제를 푼다.
이 값은 이차곡선의 판별식이고, 양수라는 것이 곧 타원이라는 뜻이다. 제약을 정규화
인자로 흡수하면 일반화 고유값 문제 하나로 풀리며, 어떤 입력에도 반드시 타원이 나온다.

구현은 Halir 와 Flusser(1998)의 수치적으로 안정한 형태를 따른다. 설계 행렬을
이차항과 일차항으로 나눠 6x6 대신 3x3 고유값 문제로 줄인다.
"""
import numpy as np


def fit_ellipse(x, y):
    """Fitzgibbon 의 직접 타원 피팅. 계수 (A,B,C,D,E,F) 를 돌려준다.

    A x^2 + B xy + C y^2 + D x + E y + F = 0
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    # 좌표를 무게중심으로 옮기고 크기를 맞춰야 6차 행렬의 조건수가 견딜 만해진다.
    mx, my = x.mean(), y.mean()
    s = max(np.hypot(x - mx, y - my).max(), 1e-12)
    u, v = (x - mx) / s, (y - my) / s

    D1 = np.column_stack([u * u, u * v, v * v])      # 이차항
    D2 = np.column_stack([u, v, np.ones_like(u)])    # 일차항
    S1, S2, S3 = D1.T @ D1, D1.T @ D2, D2.T @ D2
    try:
        T = -np.linalg.solve(S3, S2.T)
    except np.linalg.LinAlgError:
        return None
    M = S1 + S2 @ T
    # 제약 4AC - B^2 = 1 을 행렬로 옮긴 뒤 좌변에 곱해 둔 형태
    M = np.array([M[2] / 2.0, -M[1], M[0] / 2.0])
    evals, evecs = np.linalg.eig(M)
    cond = 4.0 * evecs[0] * evecs[2] - evecs[1] ** 2   # 타원 조건
    idx = np.flatnonzero(cond > 0)
    if len(idx) == 0:
        return None
    a1 = evecs[:, idx[0]].real
    coef = np.concatenate([a1, T @ a1])                # (A,B,C,D,E,F), 정규화 좌표계

    # 정규화를 되돌린다. u=(x-mx)/s 를 대입해 계수를 원좌표로 옮긴다.
    A, B, C, D, E, F = coef
    return np.array([
        A / s ** 2,
        B / s ** 2,
        C / s ** 2,
        (-2 * A * mx - B * my) / s ** 2 + D / s,
        (-B * mx - 2 * C * my) / s ** 2 + E / s,
        (A * mx ** 2 + B * mx * my + C * my ** 2) / s ** 2 - (D * mx + E * my) / s + F,
    ])


def ellipse_params(coef):
    """계수에서 중심, 반축, 기울기를 뽑는다. 타원이 아니면 None.

    이차곡선의 표준 공식을 쓴다. 반축은
        a,b = -sqrt(2*num*((A+C) +- sqrt((A-C)^2+B^2))) / (B^2-4AC)
    이고 num = A E^2 + C D^2 - B D E + (B^2-4AC) F 다.
    """
    if coef is None:
        return None
    A, B, C, D, E, F = coef
    disc = B * B - 4.0 * A * C
    if disc >= 0:
        return None                      # 포물선 또는 쌍곡선
    cx = (2.0 * C * D - B * E) / disc
    cy = (2.0 * A * E - B * D) / disc
    num = A * E * E + C * D * D - B * D * E + disc * F
    root = np.sqrt(max((A - C) ** 2 + B * B, 0.0))
    t1 = 2.0 * num * ((A + C) + root)
    t2 = 2.0 * num * ((A + C) - root)
    if t1 < 0 or t2 < 0:
        return None
    ax1 = -np.sqrt(t1) / disc
    ax2 = -np.sqrt(t2) / disc
    if not (np.isfinite(ax1) and np.isfinite(ax2)) or ax1 <= 0 or ax2 <= 0:
        return None
    theta = 0.5 * np.arctan2(-B, C - A)
    if ax2 > ax1:                        # 장축이 theta 방향이 되도록 맞춘다
        ax1, ax2 = ax2, ax1
        theta += np.pi / 2.0
    return cx, cy, ax1, ax2, (theta + np.pi) % np.pi


def conic_type(coef):
    """판별식으로 이차곡선의 종류를 판정한다."""
    if coef is None:
        return "none"
    disc = coef[1] ** 2 - 4 * coef[0] * coef[2]
    return "타원" if disc < 0 else ("포물선" if abs(disc) < 1e-12 else "쌍곡선")


def fit_conic_unconstrained(x, y):
    """제약 없는 이차곡선 피팅. 계수 벡터의 크기만 1로 고정한다.

    Fitzgibbon 의 제약이 왜 필요한지 보이기 위한 비교 대상이다. 최소 고유벡터를
    그대로 쓰므로 결과가 타원이라는 보장이 없다.
    """
    x = np.asarray(x, float); y = np.asarray(y, float)
    mx, my = x.mean(), y.mean()
    s = max(np.hypot(x - mx, y - my).max(), 1e-12)
    u, v = (x - mx) / s, (y - my) / s
    M = np.column_stack([u * u, u * v, v * v, u, v, np.ones_like(u)])
    _, _, Vt = np.linalg.svd(M, full_matrices=False)
    A, B, C, D, E, F = Vt[-1]
    return np.array([
        A / s ** 2, B / s ** 2, C / s ** 2,
        (-2 * A * mx - B * my) / s ** 2 + D / s,
        (-B * mx - 2 * C * my) / s ** 2 + E / s,
        (A * mx ** 2 + B * mx * my + C * my ** 2) / s ** 2 - (D * mx + E * my) / s + F,
    ])


def ellipse_points(n, cx, cy, a, b, theta, span_deg=360.0, start_deg=0.0,
                   noise=0.0, rng=None):
    """타원 위의 점. 노이즈는 법선 방향 대신 반축 방향으로 간단히 준다."""
    t = np.deg2rad(start_deg + np.linspace(0.0, span_deg, n, endpoint=False))
    ex, ey = a * np.cos(t), b * np.sin(t)
    if noise > 0.0:
        g = (rng or np.random.default_rng(0)).standard_normal((2, n))
        ex, ey = ex + noise * g[0], ey + noise * g[1]
    ct, st = np.cos(theta), np.sin(theta)
    return cx + ct * ex - st * ey, cy + st * ex + ct * ey

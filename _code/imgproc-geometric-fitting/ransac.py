"""이상치에 견디는 피팅 — 이미지처리 5편.

최소자승은 모든 점을 똑같이 믿는다. 잔차를 제곱해 더하므로 크게 벗어난 점 하나가
나머지 수백 개보다 큰 영향을 준다. 계측 영상에서 이런 점은 드물지 않다. 다른
구조의 엣지가 섞여 들어오거나, 4편에서 본 대로 창 안에 엣지가 둘일 때 생긴다.

RANSAC 은 모든 점을 쓰는 대신 최소 개수만 뽑아 모델을 세우고, 그 모델에 가까운
점이 몇 개인지 센다. 이 과정을 여러 번 되풀이해 가장 많은 점이 동의한 모델을
고른 뒤, 그 점들만으로 다시 정밀하게 피팅한다.
"""
import numpy as np

from circle_fit import geometric, kasa


def circle_from_three(p1, p2, p3):
    """세 점을 지나는 원. 세 점이 거의 일직선이면 None."""
    (x1, y1), (x2, y2), (x3, y3) = p1, p2, p3
    d = 2.0 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-9:
        return None
    s1, s2, s3 = x1 ** 2 + y1 ** 2, x2 ** 2 + y2 ** 2, x3 ** 2 + y3 ** 2
    cx = (s1 * (y2 - y3) + s2 * (y3 - y1) + s3 * (y1 - y2)) / d
    cy = (s1 * (x3 - x2) + s2 * (x1 - x3) + s3 * (x2 - x1)) / d
    return cx, cy, np.hypot(x1 - cx, y1 - cy)


def ransac_circle(x, y, threshold, n_iter=300, seed=0):
    """원을 정하는 최소 개수인 세 점을 뽑아 반복한다.

    threshold 는 원에서 이 거리 안에 들면 동의한 것으로 치는 값이다. 엣지 점의
    산포보다 두세 배로 잡는다. 4편에서 잰 서브픽셀 산포가 그 기준이 된다.
    """
    rng = np.random.default_rng(seed)
    n = len(x)
    best_inliers = np.zeros(n, bool)
    for _ in range(n_iter):
        idx = rng.choice(n, 3, replace=False)
        model = circle_from_three(*[(x[i], y[i]) for i in idx])
        if model is None:
            continue
        cx, cy, R = model
        inliers = np.abs(np.hypot(x - cx, y - cy) - R) < threshold
        if inliers.sum() > best_inliers.sum():
            best_inliers = inliers
    if best_inliers.sum() < 3:
        return kasa(x, y) + (best_inliers,)
    # 동의한 점들만으로 기하학적 피팅을 다시 한다.
    cx, cy, R = geometric(x[best_inliers], y[best_inliers])
    inliers = np.abs(np.hypot(x - cx, y - cy) - R) < threshold
    if inliers.sum() >= 3:
        cx, cy, R = geometric(x[inliers], y[inliers])
    return cx, cy, R, inliers


def contaminate(x, y, frac, spread, rng, cx=0.0, cy=0.0):
    """점 일부를 이상치로 바꾼다. 다른 구조에서 온 엣지를 흉내 낸 것이다."""
    n = len(x)
    k = int(round(frac * n))
    idx = rng.choice(n, k, replace=False)
    xo, yo = x.copy(), y.copy()
    ang = rng.uniform(0, 2 * np.pi, k)
    rad = spread * rng.uniform(0.4, 1.0, k)
    xo[idx] = cx + rad * np.cos(ang)
    yo[idx] = cy + rad * np.sin(ang)
    return xo, yo, idx

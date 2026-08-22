"""반사 계열 부품 4종 -- 평면거울/구면거울/빔스플리터/임의형상거울.

MATLAB 광선추적기의 MakeFlatMirror.m / MakeSphericalMirror.m /
MakeBeamSplitter.m / MakeArbitraryMirror.m과 그 CalculateRayPath_*.m을
이 글에 필요한 만큼 재구성했다. 네 함수 모두 충돌점과 법선을 찾는 방식만
다르고, 반사 자체는 reflect_vector() 하나로 끝난다.

실행: python reflection_optics.py
"""
import numpy as np


def reflect_vector(v, n):
    """법선 n 기준 완전 반사. v, n은 (vx, vy) 튜플/배열, n은 단위벡터.

    부호가 뒤집힌 n을 넣어도 결과는 같다 -- (v.n)이 부호와 함께 뒤집혀서
    상쇄되기 때문이다. 그래서 아래 네 함수는 법선 부호를 굳이 맞추지 않는다.
    """
    v = np.asarray(v, dtype=float)
    n = np.asarray(n, dtype=float)
    return v - 2 * np.dot(v, n) * n


# ---------------------------------------------------------------------------
# 1. 평면거울 -- 법선은 직선의 방향벡터에 수직인 상수 하나
# ---------------------------------------------------------------------------

def make_flat_mirror(D, x, y, theta):
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    local = np.array([[0, D / 2], [0, -D / 2]])
    boundary = (rot @ local.T).T + np.array([x, y])
    return dict(BOUNDARY=[boundary])


def reflect_flat_mirror(ray, mirror):
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y
    (x1, y1), (x2, y2) = mirror["BOUNDARY"][0]
    a_e, b_e = y2 - y1, x1 - x2
    c_e = a_e * x1 + b_e * y1
    xc, yc = np.linalg.solve([[a_r, b_r], [a_e, b_e]], [c_r, c_e])

    n = np.array([a_e, b_e]) / np.hypot(a_e, b_e)
    vout = reflect_vector([vx, vy], n)
    return np.array([[xc, yc, vout[0], vout[1]], [xc + 0.5 * vout[0], yc + 0.5 * vout[1], vout[0], vout[1]]])


# ---------------------------------------------------------------------------
# 2. 구면거울 -- 법선은 곡률중심에서 충돌점으로의 반지름 방향
# ---------------------------------------------------------------------------

def make_spherical_mirror(R, D, x, y, theta):
    w = np.linspace(np.arcsin(D / (2 * R)), -np.arcsin(D / (2 * R)), 21)
    t = R - R * np.cos(np.abs(np.arcsin(D / (2 * R))))
    Xc, Yc = -R * np.cos(w) + R - t / 2, R * np.sin(w)
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    boundary = (rot @ np.vstack([Xc, Yc])).T + np.array([x, y])
    cx, cy = rot @ np.array([R - t / 2, 0]) + np.array([x, y])
    bd_offset = 0.01 * (abs(R) - np.sqrt(R**2 - (D / 2) ** 2))
    return dict(R=R, D=D, cx=cx, cy=cy, BOUNDARY=[boundary], BD_OFFSET=bd_offset)


def reflect_spherical_mirror(ray, mirror):
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y
    phi_c, phi_s = np.arctan2(b_r, a_r), np.arctan2(a_r, b_r)
    cx, cy, R = mirror["cx"], mirror["cy"], mirror["R"]

    w_c = np.arccos((c_r - a_r * cx - b_r * cy) / (R * np.hypot(a_r, b_r))) + phi_c
    w_s = np.arcsin((c_r - a_r * cx - b_r * cy) / (R * np.hypot(a_r, b_r))) - phi_s
    xc_c, yc_c = R * np.cos(w_c) + cx, R * np.sin(w_c) + cy
    xc_s, yc_s = R * np.cos(w_s) + cx, R * np.sin(w_s) + cy

    b, off = mirror["BOUNDARY"][0], mirror["BD_OFFSET"]

    def inside(px, py):
        return (
            px <= np.max(b[:, 0]) + off and px >= np.min(b[:, 0]) - off
            and py <= np.max(b[:, 1]) + off and py >= np.min(b[:, 1]) - off
        )

    in_c, in_s = inside(xc_c, yc_c), inside(xc_s, yc_s)
    if in_c and in_s:
        xc, yc = (xc_c, yc_c) if np.hypot(xc_c - x, yc_c - y) < np.hypot(xc_s - x, yc_s - y) else (xc_s, yc_s)
    elif in_c:
        xc, yc = xc_c, yc_c
    else:
        xc, yc = xc_s, yc_s

    n = np.array([xc - cx, yc - cy]) / np.hypot(xc - cx, yc - cy)
    vout = reflect_vector([vx, vy], n)
    return np.array([[xc, yc, vout[0], vout[1]], [xc + 0.5 * vout[0], yc + 0.5 * vout[1], vout[0], vout[1]]])


# ---------------------------------------------------------------------------
# 3. 빔스플리터 -- 평면거울과 같은 법선, 반사광 + 투과광 둘 다 반환
# ---------------------------------------------------------------------------

def make_beam_splitter(D, x, y, theta):
    return make_flat_mirror(2 * D, x, y, theta)


def reflect_beam_splitter(ray, splitter):
    x, y, vx, vy = ray
    reflected = reflect_flat_mirror(ray, splitter)
    xc, yc = reflected[0, 0], reflected[0, 1]
    transmitted = np.array([[xc, yc, vx, vy], [xc + 0.5 * vx, yc + 0.5 * vy, vx, vy]])
    return reflected, transmitted


# ---------------------------------------------------------------------------
# 4. 임의형상거울 -- 변은 수선, 꼭짓점은 인접 3점을 지나는 원의 국소 곡률
# ---------------------------------------------------------------------------

def make_arbitrary_mirror(ex, ey, cx=0.0, cy=0.0, tx=0.0, ty=0.0, theta=0.0):
    ex, ey = np.asarray(ex, dtype=float), np.asarray(ey, dtype=float)
    rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    XY = (rot @ np.vstack([ex - cx, ey - cy]) + np.array([[tx], [ty]])).T
    boundary = [XY[i : i + 2] for i in range(len(XY) - 1)]
    return dict(BOUNDARY=boundary)


def _vertex_normal(x0, y0, x1, y1, x2, y2, xc, yc, dup_point):
    """꼭짓점(xc,yc)에서의 국소 곡률 법선 -- 인접 3점을 지나는 원의 중심에서 구한다.

    dup_point=True는 P에서 세 번째 점(x2,y2) 대신 두 번째 점(x1,y1)을 다시 쓴다.
    세 점이 서로 달라야 원이 유일하게 정해지므로, 이 경우 그 조건이 깨진다.
    dup_point=False가 세 번째 점을 올바르게 쓰는 버전이다.
    """
    M = np.array([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]], dtype=float)
    if dup_point:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x1**2 + y1**2], dtype=float)
    else:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x2**2 + y2**2], dtype=float)
    J = np.linalg.solve(M, P)
    acx, acy = -J[0] / 2, -J[1] / 2
    n = np.array([acx - xc, acy - yc])
    return n / np.linalg.norm(n)


def reflect_arbitrary_mirror(ray, mirror, max_bounces=20, dup_point=True):
    xs, ys, vx, vy = map(float, ray)
    boundary = mirror["BOUNDARY"]
    n_seg = len(boundary)
    path = [[xs, ys, vx, vy]]

    for _ in range(max_bounces):
        a_r, b_r = -vy, vx
        c_r = a_r * xs + b_r * ys

        best = None
        for i, seg in enumerate(boundary):
            (x1, y1), (x2, y2) = seg
            a_e, b_e = y2 - y1, x1 - x2
            c_e = a_e * x1 + b_e * y1
            M = np.array([[a_r, b_r], [a_e, b_e]])
            if abs(np.linalg.det(M)) < 1e-8:
                continue
            J = np.linalg.solve(M, [c_r, c_e])
            seg_len_sq = (x1 - x2) ** 2 + (y1 - y2) ** 2
            on_segment = (
                np.sum((J - [x1, y1]) ** 2) <= seg_len_sq
                and np.sum((J - [x2, y2]) ** 2) <= seg_len_sq
            )
            forward = vx * (J[0] - xs) + vy * (J[1] - ys) > 0
            if on_segment and forward:
                dist = np.hypot(xs - J[0], ys - J[1])
                if best is None or dist < best[0]:
                    best = (dist, J, i)

        if best is None:
            break
        _, (xc, yc), i = best
        (x1, y1), (x2, y2) = boundary[i]
        x0, y0 = boundary[i - 1][0] if i > 0 else boundary[-1][1]
        x3, y3 = boundary[i + 1][1] if i < n_seg - 1 else boundary[0][0]

        if np.isclose([xc, yc], [x1, y1]).all():
            n = _vertex_normal(x0, y0, x1, y1, x2, y2, xc, yc, dup_point)
        elif np.isclose([xc, yc], [x2, y2]).all():
            n = _vertex_normal(x1, y1, x2, y2, x3, y3, xc, yc, dup_point)
        else:
            seg_n = np.array([-(y2 - y1), x2 - x1])
            n = seg_n / np.linalg.norm(seg_n)

        vout = reflect_vector([vx, vy], n)
        vout /= np.linalg.norm(vout)
        vx, vy = vout
        xs, ys = xc + 0.5 * vx, yc + 0.5 * vy
        path.append([xc, yc, vx, vy])

    return np.array(path)


if __name__ == "__main__":
    flat = make_flat_mirror(200, 0, 0, np.deg2rad(20))
    out = reflect_flat_mirror(np.array([-100.0, 30.0, 1.0, 0.0]), flat)
    print("FlatMirror:", out.tolist())

    sph = make_spherical_mirror(-300, 200, 0, 0, 0)
    out = reflect_spherical_mirror(np.array([-100.0, 40.0, 1.0, 0.0]), sph)
    print("SphericalMirror:", out.tolist())

    bs = make_beam_splitter(100, 0, 0, np.deg2rad(45))
    r, t = reflect_beam_splitter(np.array([-100.0, 0.0, 1.0, 0.0]), bs)
    print("BeamSplitter reflected:", r.tolist())
    print("BeamSplitter transmitted:", t.tolist())

    zigzag = make_arbitrary_mirror([-100, 0, 100, 200], [0, 50, 0, 50])
    out = reflect_arbitrary_mirror(np.array([-50.0, 100.0, 0.0, -1.0]), zigzag)
    print("ArbitraryMirror (segment hit):", out.tolist())

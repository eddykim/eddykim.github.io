"""기하광학 배경이론 3편의 계산 모듈: 조리개, 동공, 주광선·주변광선, 비네팅, 조도.

광학계는 광축 위에 놓인 요소들의 목록이다. 요소마다 z 위치(mm)를 갖고, 종류는
  ("surface", n1, n2, R, 반경)   구면 굴절면 (반경은 구경 반지름, 없으면 None)
  ("thin", f, 반경)              얇은 렌즈
  ("stop", 반경)                 구멍만 있는 조리개
이다. 근축 계산은 abcd.py 의 감소각 규약 (y, n·u) 을 따른다.
문헌의 닫힌 식은 여기에 쓰지 않는다. verify_stops.py 가 대조할 독립 경로여야 하기 때문이다.
"""
from dataclasses import dataclass

import numpy as np

from abcd import prop, surface, system, thin_lens


@dataclass
class Element:
    z: float
    kind: str
    params: tuple
    radius: float | None   # 구경 반지름 (None 이면 제한 없음)


def element_matrix(e):
    if e.kind == "surface":
        n1, n2, R = e.params
        return surface(n1, n2, R)
    if e.kind == "thin":
        return thin_lens(1.0 / e.params[0])
    return np.eye(2)


def index_after(elements, i, n_object=1.0):
    """i 번째 요소 바로 뒤 매질의 굴절률."""
    n = n_object
    for e in elements[: i + 1]:
        if e.kind == "surface":
            n = e.params[1]
    return n


def matrix_between(elements, z_from, z_to, n_object=1.0):
    """평면 z_from 에서 z_to 까지 (z_from < z_to) 요소들을 지나는 행렬. 요소 위의 평면은 그 요소를 지난 뒤로 본다."""
    m = np.eye(2)
    z, n = z_from, n_object
    for e in elements:
        if e.z <= z_from:
            n = e.params[1] if e.kind == "surface" else n
            continue
        if e.z > z_to:
            break
        m = system(m, prop(e.z - z, n), element_matrix(e))
        z = e.z
        if e.kind == "surface":
            n = e.params[1]
    return system(m, prop(z_to - z, n)), n


@dataclass
class Pupil:
    z: float        # 광축 위 위치 (전역 좌표)
    radius: float
    magnification: float


def image_of_aperture(elements, i, side, n_object=1.0):
    """i 번째 요소의 구경을 앞쪽(side='object') 또는 뒤쪽(side='image') 요소들로 본 상.

    앞쪽: 물체 공간의 평면 z_p 에서 출발해 구경면까지의 행렬이 B = 0 이 되는 z_p 를 찾는다.
    뒤쪽: 구경면에서 출발해 상 공간의 평면 z_p 까지의 행렬이 B = 0 이 되는 z_p 를 찾는다.
    """
    e = elements[i]
    z_a = e.z
    if side == "object":
        z0 = elements[0].z
        pre, _ = matrix_between(elements[:i], z0 - 1e-9, z_a, n_object)   # 첫 요소부터 구경면까지
        a, b, c, d = pre.ravel()
        # 전체 = pre · T(d0), d0 = z0 - z_p (물체 공간 굴절률 n_object)
        d0 = -n_object * b / a
        return Pupil(z=z0 - d0, radius=e.radius / abs(a), magnification=a)
    # 상 쪽: 구경면 바로 뒤에서 마지막 요소까지
    z_last = elements[-1].z
    post, n_img = matrix_between(elements, z_a, z_last, n_object)
    if i == len(elements) - 1:
        post, n_img = np.eye(2), index_after(elements, i, n_object)
    a, b, c, d = post.ravel()
    d1 = -n_img * b / d
    mag = a + d1 * c / n_img
    return Pupil(z=z_last + d1, radius=e.radius * abs(mag), magnification=mag)


def find_aperture_stop(elements, z_object, n_object=1.0):
    """각 구경을 물체 공간으로 옮긴 상 가운데 축상 물점에서 가장 작은 각으로 보이는 것이 입사동이다."""
    best, best_angle, table = None, np.inf, []
    for i, e in enumerate(elements):
        if e.radius is None:
            continue
        p = image_of_aperture(elements, i, "object", n_object) if i > 0 else Pupil(e.z, e.radius, 1.0)
        angle = np.arctan2(p.radius, p.z - z_object)
        table.append((i, p, angle))
        if angle < best_angle:
            best, best_angle = i, angle
    return best, table


def trace(elements, z_start, y, u, n_start=1.0, z_end=None):
    """근축 광선 하나를 요소들을 따라 추적한다. 반환: 각 요소를 지난 직후의 (z, y, 실제 각 u, 굴절률 n)."""
    pts = [(z_start, y, u, n_start)]
    v, n, z = n_start * u, n_start, z_start
    for e in elements:
        if e.z < z_start:
            continue
        y = y + (e.z - z) * v / n
        z = e.z
        y, v = element_matrix(e) @ np.array([y, v])
        if e.kind == "surface":
            n = e.params[1]
        pts.append((z, y, v / n, n))
    if z_end is not None:
        y = y + (z_end - z) * v / n
        pts.append((z_end, y, v / n, n))
    return np.array(pts)


def marginal_and_chief(elements, z_object, y_object, n_object=1.0, z_end=None):
    """주변광선: 축상 물점 → 입사동 가장자리. 주광선: 물체 끝점 → 입사동 중심."""
    stop, _ = find_aperture_stop(elements, z_object, n_object)
    ep = image_of_aperture(elements, stop, "object", n_object) if stop > 0 else \
        Pupil(elements[0].z, elements[0].radius, 1.0)
    u_m = ep.radius / (ep.z - z_object)
    u_c = (0.0 - y_object) / (ep.z - z_object)
    marg = trace(elements, z_object, 0.0, u_m, n_object, z_end)
    chief = trace(elements, z_object, y_object, u_c, n_object, z_end)
    return stop, ep, marg, chief


def lagrange_invariant(marg, chief):
    """H = n (ȳ u - y ū): 주광선 (ȳ, ū) 와 주변광선 (y, u). 굴절률은 그 지점의 매질 값."""
    return marg[:, 3] * (chief[:, 1] * marg[:, 2] - marg[:, 1] * chief[:, 2])


def vignetting_fraction(elements, z_object, y_object, ep, n_object=1.0, n_grid=401):
    """축외 물점에서 입사동을 고르게 채운 광선 가운데 모든 구경을 통과하는 비율 (근축, 2차원 동공)."""
    s = np.linspace(-1, 1, n_grid)
    px, py = np.meshgrid(s, s)
    inside = px**2 + py**2 <= 1
    px, py = px[inside] * ep.radius, py[inside] * ep.radius
    L = ep.z - z_object
    ux, uy = px / L, (py - y_object) / L      # 물점 (0, y_object) 에서 동공점으로 가는 방향
    ok = np.ones_like(px, dtype=bool)
    for i, e in enumerate(elements):
        if e.radius is None:
            continue
        # x, y 방향은 서로 독립인 근축 사상이다: 높이 = A·y0 + B·u0 (감소각 규약)
        m, _ = matrix_between(elements[: i + 1], z_object, e.z, n_object)
        a, b = m[0, 0], m[0, 1]
        x_at = a * 0.0 + b * n_object * ux
        y_at = a * y_object + b * n_object * uy
        ok &= x_at**2 + y_at**2 <= e.radius**2 * (1 + 1e-12)
    return ok.mean()


def pupil_irradiance(a, L, h, n=300):
    """반지름 a 인 균일 람베르트 원판(복사휘도 1)이 거리 L 의 평행한 평면 위 높이 h 인 점에 만드는 조도.

    E = ∫ cosθ cosθ' / r^2 dA, 두 면이 평행하므로 cosθ = cosθ' = L/r.
    """
    r = np.linspace(0, a, n)
    phi = np.linspace(0, 2 * np.pi, 2 * n, endpoint=False)
    R, P = np.meshgrid(r, phi)
    x, y = R * np.cos(P), R * np.sin(P)
    rr2 = (x - h) ** 2 + y**2 + L**2
    f = L**2 / rr2**2 * R
    # 방위각은 주기 함수라 닫힌 구간의 단순 합(사각형 공식)이 정확하다
    return np.sum(np.trapezoid(f, r, axis=1)) * (2 * np.pi / (2 * n))

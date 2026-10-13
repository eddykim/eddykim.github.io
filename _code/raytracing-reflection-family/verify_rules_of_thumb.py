"""광선추적 3편의 경험칙 수치.

1. 곡면을 N개 변으로 쪼갠 다각형 거울: 변에 맞은 광선의 반사 방향 오차
   (매끈한 SphericalMirror 대비)와 한 변의 호각 비교. 원 위에서 뽑은 점이면
   세 점 원 꼭짓점 법선은 실제 반지름 방향과 같다.
2. 서로 수직인 두 거울(2차원 코너 반사경): 두 번 반사하면 v -> -v.
   두 거울 사이 각을 100°로 벌리면 되돌아오는 방향이 2 x 10° = 20° 어긋난다.
3. 판형 빔 분할기의 투과광 옆밀림: d = t sinθ (1 - cosθ / sqrt(n^2 - sin^2θ)).

실행: python verify_rules_of_thumb.py  (한국어·영문 결과를 차례로 출력한다)
"""
import numpy as np

from reflection_optics import (
    _vertex_normal,
    make_arbitrary_mirror,
    make_spherical_mirror,
    reflect_arbitrary_mirror,
    reflect_spherical_mirror,
)

R, D = 300.0, 200.0
ARC = np.rad2deg(2 * np.arcsin(D / (2 * R)))   # 전체 호각


def polygon_like(N):
    """make_spherical_mirror 와 같은 원호 위의 N+1개 점을 이은 다각형 거울."""
    sph = make_spherical_mirror(-R, D, 0, 0, 0)
    Rr = -R
    w = np.linspace(np.arcsin(D / (2 * Rr)), -np.arcsin(D / (2 * Rr)), N + 1)
    t = Rr - Rr * np.cos(abs(np.arcsin(D / (2 * Rr))))
    X, Y = -Rr * np.cos(w) + Rr - t / 2, Rr * np.sin(w)
    return sph, make_arbitrary_mirror(X, Y)


def max_dir_err(N):
    sph, arb = polygon_like(N)
    worst = 0.0
    for h in np.linspace(-95, 95, 381):
        ray = np.array([-200.0, h, 1.0, 0.0])
        v1 = reflect_spherical_mirror(ray, sph)[0, 2:4]
        o2 = reflect_arbitrary_mirror(ray, arb, max_bounces=1, dup_point=False)
        if len(o2) < 2:
            continue
        v2 = o2[1, 2:4]
        c = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
        worst = max(worst, np.rad2deg(np.arccos(np.clip(c, -1, 1))))
    return worst


POLY = [(N, ARC / N, max_dir_err(N)) for N in [10, 20, 39, 50, 100, 390]]

sph20, arb20 = polygon_like(20)
B = arb20["BOUNDARY"]
(x0, y0), (x1, y1) = B[4]
x2, y2 = B[5][1]
n_u = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, False)
n_true = np.array([sph20["cx"] - x1, sph20["cy"] - y1])
n_true /= np.linalg.norm(n_true)
VERTEX_EXACT = np.rad2deg(np.arccos(np.clip(abs(n_u @ n_true), -1, 1)))


def corner(opening_deg, tilt_deg):
    half = np.deg2rad(opening_deg / 2)
    yy = 100 / np.tan(half)
    m = make_arbitrary_mirror([-100, 0, 100], [yy, 0, yy])
    a = np.deg2rad(tilt_deg)
    v = np.array([np.sin(a), -np.cos(a)])
    path = reflect_arbitrary_mirror(np.r_[-20.0, 150.0, v], m, max_bounces=5, dup_point=False)
    vout = path[-1, 2:4]
    return len(path) - 1, float(np.dot(vout, v)), np.rad2deg(np.arccos(np.clip(-vout @ v, -1, 1)))


CORNER = [(90, t, *corner(90, t)) for t in [0, 10, 25]] + [(100, t, *corner(100, t)) for t in [0, 10]]

n_glass, th = 1.5168, np.deg2rad(45)   # N-BK7, 587.6 nm
SHIFT = [(t, t * np.sin(th) * (1 - np.cos(th) / np.sqrt(n_glass**2 - np.sin(th) ** 2)))
         for t in [1.0, 3.0, 6.0]]

LABELS = {
    "ko": {
        "poly_head": "다각형 거울 (R={:.0f}mm, 지름 {:.0f}mm, 전체 호각 {:.3f}°)",
        "poly": "N={:>3}변: 한 변의 호각={:.3f}°  반사 방향 최대 오차={:.3f}°",
        "vertex": "원 위 점 20변 다각형의 꼭짓점: 세 점 원 법선과 실제 반지름 방향의 차이={:.1f}°",
        "corner": "{}° V자 거울, 기울기 {:>2}°: 반사 {}번  v_out·v_in={:.12f}  되돌아오는 방향과의 차이={:.3f}°",
        "shift": "판형 빔 분할기 (n={}, 45°) 두께 {:.0f}mm: 투과광 옆밀림={:.3f}mm",
    },
    "en": {
        "poly_head": "polygon mirror (R={:.0f}mm, diameter {:.0f}mm, total arc {:.3f}°)",
        "poly": "N={:>3} edges: arc per edge={:.3f}°  max reflected-direction error={:.3f}°",
        "vertex": "vertex of 20-edge polygon on the circle: three-point normal vs true radius={:.1f}°",
        "corner": "{}° V mirror, tilt {:>2}°: {} reflections  v_out·v_in={:.12f}  off retro-direction by {:.3f}°",
        "shift": "plate beam splitter (n={}, 45°) thickness {:.0f}mm: lateral shift of transmitted ray={:.3f}mm",
    },
}

for lang, L in LABELS.items():
    print(f"--- [{lang}] ---")
    print(L["poly_head"].format(R, D, ARC))
    for N, per, err in POLY:
        print(L["poly"].format(N, per, err))
    print(L["vertex"].format(VERTEX_EXACT))
    for op, t, nb, dot, dev in CORNER:
        print(L["corner"].format(op, t, nb, dot, dev))
    for t, d in SHIFT:
        print(L["shift"].format(n_glass, t, d))

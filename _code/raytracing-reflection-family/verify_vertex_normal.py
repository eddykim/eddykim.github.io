"""3편 3절 표: 세 번째 점을 재사용한 법선(dup_point=True)과 세 점을 모두 쓴
법선(dup_point=False)이 꼭짓점의 좌우 비대칭 정도에 따라 얼마나 달라지는지
세 가지 꼭짓점 형태로 비교한다.

대칭 꼭짓점은 인접한 두 변의 길이·기울기가 같고, 비대칭 꼭짓점으로 갈수록
두 변의 길이·기울기 차이를 키웠다.

dup_point=True 로 얻은 원이 세 번째 점을 얼마나 빗나가는지도 함께 출력한다.
법선이 틀어지는 정도는 결국 이 빗나감의 크기를 따라간다.

실행: python verify_vertex_normal.py  (한국어·영문 결과를 차례로 출력한다)
"""
import numpy as np

from reflection_optics import _vertex_normal

LABELS = {
    "ko": {
        "cases": ["대칭 꼭짓점      ", "비대칭 꼭짓점     ", "더 비대칭한 꼭짓점"],
        "row": "{}: n_dup=[{:.3f}, {:.3f}]  n_uniq=[{:.3f}, {:.3f}]  차이={:.3f}deg  "
               "세 번째 점 빗나감={:.1f}mm",
    },
    "en": {
        "cases": ["symmetric vertex ", "asymmetric vertex", "more asymmetric  "],
        "row": "{}: n_dup=[{:.3f}, {:.3f}]  n_uniq=[{:.3f}, {:.3f}]  diff={:.3f}deg  "
               "third point missed by {:.1f}mm",
    },
}

# 대칭: 두 변의 길이(100)와 기울기(|dy/dx|=0.5)가 좌우로 같다
# 비대칭: 왼쪽 변(60, 상승)보다 오른쪽 변(130, 완만한 하강)이 더 길고 기울기도 다르다
# 더 비대칭: 왼쪽 변을 더 짧고 가파르게, 오른쪽 변을 더 길고 완만하게
CASES = [
    (-100, 0, 0, 50, 100, 0),
    (-60, 0, 0, 50, 130, -10),
    (-30, 0, 0, 50, 180, -40),
]


def circle_center(x0, y0, x1, y1, x2, y2, dup_point):
    """_vertex_normal 이 내부에서 푸는 것과 같은 원의 중심."""
    M = np.array([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]], dtype=float)
    if dup_point:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x1**2 + y1**2], dtype=float)
    else:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x2**2 + y2**2], dtype=float)
    J = np.linalg.solve(M, P)
    return np.array([-J[0] / 2, -J[1] / 2])


def measure(x0, y0, x1, y1, x2, y2):
    n_dup = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, dup_point=True)
    n_uniq = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, dup_point=False)
    ang = np.rad2deg(np.arccos(np.clip(np.dot(n_dup, n_uniq), -1, 1)))
    # dup 원은 앞의 두 점은 지나고 세 번째 점만 빗나간다. 그 빗나간 거리.
    c = circle_center(x0, y0, x1, y1, x2, y2, dup_point=True)
    radius = np.hypot(x1 - c[0], y1 - c[1])
    miss = abs(np.hypot(x2 - c[0], y2 - c[1]) - radius)
    return n_dup, n_uniq, ang, miss


RESULTS = [measure(*case) for case in CASES]

for lang, L in LABELS.items():
    print(f"--- [{lang}] ---")
    for name, (n_dup, n_uniq, ang, miss) in zip(L["cases"], RESULTS):
        print(L["row"].format(name, n_dup[0], n_dup[1], n_uniq[0], n_uniq[1], ang, miss))

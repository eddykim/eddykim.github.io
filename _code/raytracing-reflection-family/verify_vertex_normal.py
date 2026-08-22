"""3편 2절 표: 세 번째 점을 재사용한 법선(dup_point=True)과 세 점을 모두 쓴
법선(dup_point=False)이 꼭짓점의 좌우 비대칭 정도에 따라 얼마나 달라지는지
세 가지 꼭짓점 형태로 비교한다.

대칭 꼭짓점은 인접한 두 변의 길이·기울기가 같고, 비대칭 꼭짓점으로 갈수록
두 변의 길이·기울기 차이를 키웠다.

실행: python verify_vertex_normal.py
"""
import numpy as np

from reflection_optics import _vertex_normal, reflect_vector


def compare(label, x0, y0, x1, y1, x2, y2):
    n_dup = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, dup_point=True)
    n_uniq = _vertex_normal(x0, y0, x1, y1, x2, y2, x1, y1, dup_point=False)
    ang = np.rad2deg(np.arccos(np.clip(np.dot(n_dup, n_uniq), -1, 1)))
    print(f"{label}: n_dup=[{n_dup[0]:.3f}, {n_dup[1]:.3f}]  "
          f"n_uniq=[{n_uniq[0]:.3f}, {n_uniq[1]:.3f}]  차이={ang:.3f}deg")
    return n_dup, n_uniq, ang


# 대칭 꼭짓점: 두 변의 길이(100)와 기울기(|dy/dx|=0.5)가 좌우로 같다
compare("대칭 꼭짓점      ", -100, 0, 0, 50, 100, 0)

# 비대칭 꼭짓점: 왼쪽 변(60, 상승)보다 오른쪽 변(130, 완만한 하강)이 더 길고 기울기도 다르다
compare("비대칭 꼭짓점     ", -60, 0, 0, 50, 130, -10)

# 더 비대칭한 꼭짓점: 왼쪽 변을 더 짧고 가파르게, 오른쪽 변을 더 길고 완만하게
compare("더 비대칭한 꼭짓점", -30, 0, 0, 50, 180, -40)

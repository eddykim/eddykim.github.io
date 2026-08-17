"""BD_OFFSET IndexError 재현과 수정.

make_spherical_lens()는 BD_OFFSET을 1차원 배열로 만든다:
    BD_OFFSET = 0.01 * np.array([offset_R1, offset_R2])   # shape (2,)

그런데 MATLAB 원본은 BD_OFFSET(1,1), BD_OFFSET(1,2)처럼 2차원(1xN 행벡터)으로
인덱싱한다. 이걸 그대로 옮기면서 처음 짠 Python 코드는 BD_OFFSET[0, 0]
(2차원 인덱싱)을 썼다. 실행하면 다음 에러가 난다.

실행: python bug_demo.py
"""
import sys
import traceback

import numpy as np

BD_OFFSET = 0.01 * np.array([1.234, 2.345])  # make_spherical_lens가 실제로 만드는 형태와 동일

print("=== 처음 짠 코드 (MATLAB의 BD_OFFSET(1,1)을 그대로 옮김) ===")
try:
    offset = BD_OFFSET[0, 0]
    print(f"BD_OFFSET[0, 0] = {offset}")
except IndexError:
    traceback.print_exc(file=sys.stdout)

print()
print("=== 수정 (1차원 배열이므로 인덱스 하나만 필요) ===")
offset = BD_OFFSET[0]
print(f"BD_OFFSET[0] = {offset}")

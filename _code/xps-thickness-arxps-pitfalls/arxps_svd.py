"""XPS 기초 5편 4절: 각도 데이터에 독립된 정보가 몇 개 들어 있는가.

실행: python arxps_svd.py

ARXPS의 측정식 I(θ) ∝ ∫ c(z) exp(-z / (λ cosθ)) dz 를 행렬 K c 로 적는다.
행은 그림 3과 같은 29개 각도(0–70°, 2.5° 간격), 열은 overlayer.py의 깊이 격자다.
각 행을 노름 1로 맞춰 측정마다 상대 잡음이 같다고 보고 특이값을 구한다.
첫 특이값에 대한 비가 측정의 상대 잡음보다 작은 성분은 잡음에 묻힌다.
"""
import numpy as np

import overlayer as O

angles = np.arange(0, 71, 2.5)
dz = O.Z[1] - O.Z[0]
K = np.exp(-O.Z[None, :] / (O.L * np.cos(np.radians(angles))[:, None])) * dz
K /= np.linalg.norm(K, axis=1, keepdims=True)
s = np.linalg.svd(K, compute_uv=False)
rel = s / s[0]

print(f"각도 {angles.size}개 (0–70°, 2.5° 간격), λ = {O.L} nm, 깊이 0–{O.Z[-1]:.0f} nm")
print("상대 특이값: " + ", ".join(f"{v:.2g}" for v in rel[:6]))
band = 100 * O.total_band(2.0, angles)
print(f"그림 3 불확도 띠: {band.min():.2f}% ~ {band.max():.2f}%")
for noise in (0.05, 0.01, 0.005):
    print(f"  상대 잡음 {100 * noise:.1f}%보다 큰 성분: {(rel > noise).sum()}개")

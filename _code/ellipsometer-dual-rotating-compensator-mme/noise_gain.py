"""조건수의 비와 실제 오차의 비가 어떻게 다른지 계산한다 (본문 6절).

실행: python noise_gain.py

조건수는 측정 잡음이 뮬러 요소로 번질 때 상대오차가 커질 수 있는 최악의 한도다.
같은 잡음에서 뮬러 요소 오차가 몇 배 줄어드는지는 이와 다른 양이다.
여기서는 5:3 배치에서 광학주기당 36개 표본을 등간격으로 잡고, 표본마다 같은 크기의
독립 잡음이 실린다고 둔다. 최소제곱 해의 공분산은 sigma^2 (W^T W)^-1 이므로

  열여섯 요소 오차의 합   sqrt(trace (W^T W)^-1)
  요소별 표준편차         sqrt(diag (W^T W)^-1)
  최악 방향의 증폭        1 / sigma_min(W)

를 90°와 128°에서 비교한다. W 는 표본 영역 행렬(36 x 16)이다. 세기 계산은
mme.intensity(부품 뮬러 행렬의 곱)를 쓴다. 팩트체크 때 존스 행렬 경로로 따로 계산한
값(오차 합 14.00 -> 6.98, 2.0배)과 비가 같은지 확인한다. 절대값은 세기의 상수 배율에
따라 달라지지만 비는 배율과 무관하다.
"""
import numpy as np

from mme import RATIO, N_SCAN, intensity, optical_cycle

LABELS = [f"M{j}{k}" for j in range(1, 5) for k in range(1, 5)]


def sample_matrix(delta, ratio=RATIO, n_pts=N_SCAN):
    """표본 세기 = W @ (뮬러 요소 16개) 의 W. 단위행렬 E_jk 를 넣어 열을 만든다."""
    C = optical_cycle(n_pts)
    cols = []
    for j in range(4):
        for k in range(4):
            E = np.zeros((4, 4))
            E[j, k] = 1.0
            cols.append(intensity(E, C, delta, delta, ratio))
    return np.array(cols).T


def noise_stats(delta_deg, ratio=RATIO, n_pts=N_SCAN):
    W = sample_matrix(np.deg2rad(delta_deg), ratio, n_pts)
    s = np.linalg.svd(W, compute_uv=False)
    cov = np.linalg.inv(W.T @ W)
    return {
        "cond": s[0] / s[-1],
        "total": np.sqrt(np.trace(cov)),
        "per": np.sqrt(np.diag(cov)),
        "worst": 1.0 / s[-1],
    }


def optimum(key, lo=100.0, hi=150.0, step=0.1):
    deg = np.arange(lo, hi + 1e-9, step)
    val = np.array([noise_stats(d)[key] for d in deg])
    i = int(np.argmin(val))
    return deg[i], val[i]


if __name__ == "__main__":
    a, b = noise_stats(90.0), noise_stats(128.0)
    print(f"회전비 {RATIO[0]}:{RATIO[1]}, 주기당 표본 {N_SCAN}개, 표본마다 같은 잡음")
    print(f"{'':24s}{'90°':>10s}{'128°':>10s}{'비(90/128)':>14s}")
    for key, name in (("cond", "조건수(표본 영역)"),
                      ("total", "16요소 오차 합"),
                      ("worst", "최악 방향 증폭")):
        print(f"  {name:22s}{a[key]:10.3f}{b[key]:10.3f}{a[key]/b[key]:12.2f}배")
    ratio = a["per"] / b["per"]
    print("  요소별 표준편차 비(90/128):")
    for row in range(4):
        print("    " + "  ".join(f"{LABELS[4*row+c]} {ratio[4*row+c]:.2f}"
                                 for c in range(4)))
    print(f"  요소별 비의 범위: {ratio[:-1].min():.2f}~{ratio.max():.2f}배, "
          f"M44 만 {ratio[-1]:.2f}배(128°에서 {1/ratio[-1]:.1f}배 커진다)")
    dk, vk = optimum("cond")
    de, ve = optimum("total")
    print(f"  최적 지연량: 조건수 기준 {dk:.1f}°, 오차 합 기준 {de:.1f}°")

    # 존스 행렬 경로(팩트체크 독립 계산)와 비가 같은지 확인한다
    ok = (abs(a["cond"] / b["cond"] - 3.84) < 0.01
          and abs(a["total"] / b["total"] - 14.000 / 6.976) < 0.01)
    print(f"\n[{'OK  ' if ok else 'FAIL'}] 조건수 비 3.84배, 오차 합 비 2.01배"
          f"(존스 행렬 경로와 일치)")
    raise SystemExit(0 if ok else 1)

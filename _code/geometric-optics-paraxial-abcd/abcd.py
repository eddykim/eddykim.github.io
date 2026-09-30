"""기하광학 배경이론 2편의 근축 행렬 모듈.

광선 벡터는 (y, v) 이고 v = n·u 는 감소각(reduced angle)이다 (Goodman 부록 B 규약).
이 규약에서는 모든 기본 행렬의 행렬식이 1 이다. 보통 각 (y, u) 규약(Saleh & Teich, Fowles)으로
바꾸는 함수는 to_angle_convention() 이다.

행렬은 빛이 지나는 순서대로 system(M1, M2, ...) 에 넘긴다. 곱은 M_N ··· M_2 M_1 이다.
문헌의 두꺼운 렌즈 닫힌 식은 여기에 쓰지 않는다. verify_abcd.py 가 대조할 독립 경로여야 하기 때문이다.
"""
from dataclasses import dataclass

import numpy as np


# ── 기본 행렬 ─────────────────────────────────────────────
def prop(d, n=1.0):
    """굴절률 n 인 균질 매질을 거리 d 만큼 진행."""
    return np.array([[1.0, d / n], [0.0, 1.0]])


def surface(n1, n2, radius):
    """반경 radius 인 구면 경계 (n1 → n2). 곡률 중심이 뒤쪽이면 radius > 0."""
    power = 0.0 if np.isinf(radius) else (n2 - n1) / radius
    return np.array([[1.0, 0.0], [-power, 1.0]])


def thin_lens(power):
    """굴절력 power 인 얇은 렌즈."""
    return np.array([[1.0, 0.0], [-power, 1.0]])


def system(*mats):
    """빛이 지나는 순서대로 받은 행렬들의 곱."""
    m = np.eye(2)
    for x in mats:
        m = x @ m
    return m


def to_angle_convention(m, n_in, n_out):
    """(y, n·u) 규약의 행렬을 (y, u) 규약으로 바꾼다. 행렬식은 n_in / n_out 이 된다."""
    return np.diag([1.0, 1.0 / n_out]) @ m @ np.diag([1.0, n_in])


# ── 주요점 ───────────────────────────────────────────────
@dataclass
class Cardinal:
    """입력면·출력면 기준 주요점 위치 (빛의 진행 방향이 +)."""

    power: float        # 굴절력 Φ = -C
    f_obj: float        # 물체 쪽 초점거리 n_in / Φ
    f_img: float        # 상 쪽 초점거리 n_out / Φ
    z_front_focus: float   # 입력면 기준 앞 초점 F 위치
    z_back_focus: float    # 출력면 기준 뒤 초점 F' 위치
    z_front_principal: float  # 입력면 기준 앞 주평면 H 위치
    z_back_principal: float   # 출력면 기준 뒤 주평면 H' 위치
    z_front_nodal: float      # 입력면 기준 앞 절점 N
    z_back_nodal: float       # 출력면 기준 뒤 절점 N'


def cardinal_points(m, n_in=1.0, n_out=1.0):
    """계 행렬에서 주요점을 읽는다.

    광축에 평행한 광선 (1, 0) 은 (A, C) 로 나온다 → 뒤 초점은 출력면에서 -n_out·A/C.
    출력 쪽에서 평행해지는 광선 (C y + D v = 0) 을 거꾸로 따라가면 앞 초점은 입력면에서 n_in·D/C.
    """
    a, b, c, d = m.ravel()
    power = -c
    f_obj, f_img = n_in / power, n_out / power
    z_f = n_in * d / c
    z_fp = -n_out * a / c
    z_h = z_f + f_obj
    z_hp = z_fp - f_img
    # 절점: 주점에서 (n_out - n_in)/Φ 만큼 이동 (앞뒤 매질이 같으면 주점과 겹친다)
    shift = (n_out - n_in) / power
    return Cardinal(power, f_obj, f_img, z_f, z_fp, z_h, z_hp, z_h + shift, z_hp + shift)


def image_distance(m, s_obj, n_in=1.0, n_out=1.0):
    """입력면 앞 s_obj 에 놓인 물체의 상이 출력면 뒤 어디에 맺히는지 (전체 B = 0 조건).

    반환: (출력면 기준 상 거리, 전체 행렬)
    """
    a, b, c, d = m.ravel()
    x = s_obj / n_in
    s_img = -n_out * (a * x + b) / (c * x + d)
    total = system(prop(s_obj, n_in), m, prop(s_img, n_out))
    return s_img, total


# ── 굴절률 분포 매질 ─────────────────────────────────────
def grin(n0, alpha, length):
    """n^2 = n0^2(1 - α^2 y^2) 매질을 length 만큼 지나는 근축 행렬 (감소각 규약).

    근축 광선 방정식 y'' = -α^2 y 의 해 y = y0 cos αz + (u0/α) sin αz 에서 나온다.
    """
    k = alpha * length
    return np.array([[np.cos(k), np.sin(k) / (n0 * alpha)],
                     [-n0 * alpha * np.sin(k), np.cos(k)]])


def grin_sliced(n0, alpha, length, n_slices):
    """같은 매질을 얇은 균질 판 + 얇은 렌즈의 반복으로 근사한 행렬 (독립 경로)."""
    dz = length / n_slices
    # 두께 dz 의 판 하나가 갖는 굴절력: n(y) ≈ n0 (1 - α²y²/2) → 광로차 -n0 α² y² dz / 2
    lens = thin_lens(n0 * alpha**2 * dz)
    half = prop(dz / 2, n0)
    cell = system(half, lens, half)
    return np.linalg.matrix_power(cell, n_slices)


# ── 주기 광학계 ─────────────────────────────────────────
def periodic_cell(f, d):
    """거리 d 진행 후 초점거리 f 인 얇은 렌즈 (Saleh & Teich 예제 1.4-1 의 단위계)."""
    return system(prop(d), thin_lens(1.0 / f))


def stability(m):
    """b = (A + D)/2, |b| ≤ 1 이면 안정. 안정이면 한 단마다 위상 φ = arccos b 만큼 돈다."""
    b = 0.5 * np.trace(m)
    phi = np.arccos(b) if abs(b) <= 1 else np.nan
    return b, phi


def iterate(m, ray, n):
    """행렬 m 을 n 번 적용한 광선 궤적 (n+1 개 점)."""
    out = [np.asarray(ray, dtype=float)]
    for _ in range(n):
        out.append(m @ out[-1])
    return np.array(out)

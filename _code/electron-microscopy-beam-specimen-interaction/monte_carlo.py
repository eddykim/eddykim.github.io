"""전자-시료 상호작용의 단일산란 몬테카를로 모형.

전자현미경 배경이론 3편에서 쓰는 계산을 모아 둔 모듈이다. 교과서적인
단일산란(single-scattering) 모형이며, 구성은 다음 네 조각뿐이다.

1. 차폐 러더퍼드 단면적으로 평균자유행로를 구한다
2. 지수분포에서 다음 산란까지의 거리를 뽑는다
3. 차폐 러더퍼드 분포에서 산란각을 뽑고 방위각은 균일분포에서 뽑는다
4. 이동 구간에서 Bethe 식으로 에너지를 잃는다

전자가 표면(z < 0)을 넘어 되돌아 나오면 후방산란으로 집계하고, 에너지가
E_MIN 아래로 떨어지면 그 자리에서 멈춘 것으로 본다.

단위는 관례를 따라 에너지 keV, 길이 cm 로 계산하고 표시할 때만 µm/nm 로 바꾼다.
"""
from dataclasses import dataclass

import numpy as np

N_A = 6.02214076e23
E_MIN = 0.5          # keV, 이 아래로 떨어지면 추적을 멈춘다
MAX_STEPS = 20000


@dataclass(frozen=True)
class Material:
    name: str
    Z: float          # 원자번호
    A: float          # 원자량 [g/mol]
    rho: float        # 밀도 [g/cm^3]


C = Material("C", 6, 12.011, 2.26)
SI = Material("Si", 14, 28.086, 2.329)
CU = Material("Cu", 29, 63.546, 8.96)
AG = Material("Ag", 47, 107.868, 10.49)
W = Material("W", 74, 183.84, 19.25)
AU = Material("Au", 79, 196.967, 19.30)


def screening(mat, e_kev):
    """차폐 상수 alpha. 원자핵의 전하를 궤도전자가 가리는 정도를 나타낸다."""
    return 3.4e-3 * mat.Z**0.67 / e_kev


def cross_section(mat, e_kev):
    """차폐 러더퍼드 탄성산란 단면적 [cm^2/atom].

    Z^2 에 비례하고 E^2 에 반비례한다 — 무거운 원소일수록, 느린 전자일수록
    더 자주 크게 휜다는 3편의 논지가 이 한 줄에서 나온다.
    마지막 괄호는 상대론 보정이다.
    """
    a = screening(mat, e_kev)
    rel = ((e_kev + 511.0) / (e_kev + 1022.0)) ** 2
    return 5.21e-21 * (mat.Z**2 / e_kev**2) * (4.0 * np.pi / (a * (1.0 + a))) * rel


def mean_free_path(mat, e_kev):
    """탄성산란 평균자유행로 [cm]."""
    return mat.A / (N_A * mat.rho * cross_section(mat, e_kev))


def mean_ionization_potential(mat):
    """평균 이온화 퍼텐셜 J [keV] (Berger-Seltzer)."""
    return (9.76 * mat.Z + 58.5 * mat.Z**-0.19) * 1e-3


def stopping_power(mat, e_kev):
    """Bethe 저지능 dE/ds [keV/cm]. 음수로 돌려준다.

    저에너지에서 발산하지 않도록 Joy-Luo 의 수정항(0.85 J)을 넣었다.
    """
    j = mean_ionization_potential(mat)
    return -78500.0 * mat.rho * (mat.Z / (mat.A * e_kev)) * np.log(
        1.166 * (e_kev + 0.85 * j) / j)


def kanaya_okayama_range(mat, e0_kev):
    """Kanaya-Okayama 침투 깊이 [cm]. 상호작용 부피의 크기를 가늠하는 경험식."""
    return 0.0276 * mat.A * e0_kev**1.67 / (mat.Z**0.89 * mat.rho) * 1e-4


def reuter_eta(z):
    """Reuter 의 경험식으로 얻는 후방산란계수. 몬테카를로 결과와 비교할 기준선이다."""
    z = np.asarray(z, dtype=float)
    return -0.0254 + 0.016 * z - 1.86e-4 * z**2 + 8.3e-7 * z**3


def _scatter_direction(c, cos_t, phi):
    """방향코사인 c 를 극각 theta, 방위각 phi 만큼 회전시킨다."""
    cx, cy, cz = c
    sin_t = np.sqrt(max(0.0, 1.0 - cos_t * cos_t))
    cos_p, sin_p = np.cos(phi), np.sin(phi)
    perp = np.sqrt(max(0.0, 1.0 - cz * cz))
    if perp < 1e-8:                      # 진행 방향이 z 축과 거의 나란할 때
        return np.array([sin_t * cos_p, sin_t * sin_p, np.sign(cz) * cos_t])
    return np.array([
        cx * cos_t + sin_t * (cx * cz * cos_p - cy * sin_p) / perp,
        cy * cos_t + sin_t * (cy * cz * cos_p + cx * sin_p) / perp,
        cz * cos_t - sin_t * perp * cos_p,
    ])


def trace_electron(mat, e0_kev, rng, keep_path=True):
    """전자 한 개를 추적한다.

    반환: dict
      path         (n,3) 궤적 [cm] — keep_path 가 False 면 None
      backscattered 표면으로 되돌아 나왔는가
      depths/losses 구간 중점의 깊이 [cm] 와 그 구간에서 잃은 에너지 [keV]
      energies     각 구간에 들어갈 때의 전자 에너지 [keV]
      exit_depth   후방산란 직전 마지막 산란의 깊이 [cm], 아니면 None
      max_depth    도달한 최대 깊이 [cm]
    """
    pos = np.zeros(3)
    direction = np.array([0.0, 0.0, 1.0])     # 표면에 수직으로 입사
    e = float(e0_kev)
    path = [pos.copy()]
    depths, losses, energies = [], [], []
    backscattered = False
    last_pos = pos.copy()

    for _ in range(MAX_STEPS):
        lam = mean_free_path(mat, e)
        step = -lam * np.log(max(rng.random(), 1e-12))
        new_pos = pos + direction * step

        if new_pos[2] < 0.0:                  # 표면을 넘어섰다 = 후방산란
            backscattered = True
            last_pos = pos.copy()
            if keep_path:
                t = -pos[2] / direction[2]    # 표면과의 교점까지만 그린다
                path.append(pos + direction * t)
            break

        # 이동 구간에서 잃은 에너지를 구간 중점의 깊이에 쌓는다.
        loss = min(-stopping_power(mat, e) * step, e - E_MIN * 0.5)
        depths.append(0.5 * (pos[2] + new_pos[2]))
        losses.append(max(loss, 0.0))
        energies.append(e)
        e -= loss
        pos = new_pos
        if keep_path:
            path.append(pos.copy())
        if e <= E_MIN:
            break

        # 차폐 러더퍼드 분포에서 산란각을 뽑는다.
        a = screening(mat, e)
        r = rng.random()
        cos_t = 1.0 - 2.0 * a * r / (1.0 + a - r)
        direction = _scatter_direction(direction, cos_t, 2.0 * np.pi * rng.random())

    arr = np.array(path) if keep_path else None
    all_z = arr[:, 2] if keep_path else np.array(depths or [0.0])
    return {
        "path": arr,
        "backscattered": backscattered,
        "depths": np.array(depths),
        "losses": np.array(losses),
        "energies": np.array(energies),
        "exit_depth": float(last_pos[2]) if backscattered else None,
        "max_depth": float(np.max(all_z)) if len(all_z) else 0.0,
    }


def simulate(mat, e0_kev, n, seed=0, keep_path=True):
    """전자 n 개를 추적해 결과 목록을 돌려준다."""
    rng = np.random.default_rng(seed)
    return [trace_electron(mat, e0_kev, rng, keep_path) for _ in range(n)]


def backscatter_coefficient(mat, e0_kev, n, seed=0):
    """몬테카를로로 얻은 후방산란계수."""
    res = simulate(mat, e0_kev, n, seed=seed, keep_path=False)
    return sum(r["backscattered"] for r in res) / len(res)


if __name__ == "__main__":
    print(f"{'물질':>4} {'E0':>6} {'K-O 침투깊이':>14} {'MC eta':>8} {'Reuter eta':>11}")
    for mat in (C, SI, CU, W):
        for e0 in (5.0, 30.0):
            r_ko = kanaya_okayama_range(mat, e0) * 1e4      # µm
            eta = backscatter_coefficient(mat, e0, 2000, seed=1)
            print(f"{mat.name:>4} {e0:5.0f}kV {r_ko:11.3f} µm "
                  f"{eta:8.3f} {reuter_eta(mat.Z):11.3f}")

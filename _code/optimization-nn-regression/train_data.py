"""학습 데이터 합성과 augmentation — 최적화 5편 실험용.

전방 모델(Airy 공식)이 사실상 공짜이므로 학습 데이터를 원하는 만큼 합성할 수 있다.
1~4편과 같은 Air/SiO2/Si 단층 구조를 쓴다.

augmentation 축은 Kwak et al. 2021 (Light: Adv. Manuf. 2(1), 9-19) 의 설정을 따랐다.
그 논문은 3D 다층 반도체 계측에서 시료당 40개씩 125→5,000 개로 늘렸고, 교란 축이
"세로 오프셋"과 "가로(파장축) 오프셋"이었다. 임의의 잡음이 아니라 실제 분광기의
오차 모드라는 점이 핵심이다.
"""
import numpy as np

from reflectance_model import reflectance, N_SIO2


def make_spectra(d_values, wavelength_nm, n1=N_SIO2):
    """두께 배열 → 반사율 스펙트럼 행렬 (N, n_wl)."""
    d_values = np.atleast_1d(np.asarray(d_values, dtype=float))
    return np.stack([reflectance(d, wavelength_nm, n1=n1) for d in d_values])


def make_dataset(n_samples, d_range, wavelength_nm, rng, noise_std=0.0, n1=N_SIO2):
    """[d_lo, d_hi] 에서 두께를 균일 추출해 스펙트럼을 만든다.

    noise_std > 0 이면 고정 표준편차의 가우시안 노이즈를 더한다(단순 노이즈 훈련).
    """
    d_lo, d_hi = d_range
    d = rng.uniform(d_lo, d_hi, size=n_samples)
    X = make_spectra(d, wavelength_nm, n1=n1)
    if noise_std > 0:
        X = X + rng.normal(0.0, noise_std, size=X.shape)
    return X, d


def augment(d_values, wavelength_nm, rng, n_per_sample=1,
            noise_std_range=(0.0, 0.02),
            vertical_offset=0.04, vertical_scale=0.05,
            lateral_shift_nm=6.0, n1=N_SIO2,
            use_noise=True, use_vertical=True, use_lateral=True):
    """실측 분광기의 오차 모드를 흉내 낸 증강.

    - noise    : 표본마다 노이즈 표준편차를 noise_std_range 에서 뽑아 주입
    - vertical : R -> (1+s)*R + c  (광량 스케일 s, 베이스라인 c)
    - lateral  : 파장축이 delta 만큼 어긋난 분광기 — reflectance 를 shifted 파장에서 평가

    가로 오프셋은 보간 대신 어긋난 파장에서 전방 모델을 직접 평가한다(보간 오차 없음).
    반환: X (N*n_per_sample, n_wl), d (N*n_per_sample,)
    """
    d_values = np.atleast_1d(np.asarray(d_values, dtype=float))
    Xs, ds = [], []

    for d in d_values:
        for _ in range(n_per_sample):
            wl = wavelength_nm
            if use_lateral and lateral_shift_nm > 0:
                delta = rng.uniform(-lateral_shift_nm, lateral_shift_nm)
                wl = wavelength_nm + delta
            r = reflectance(d, wl, n1=n1)

            if use_vertical:
                s = rng.uniform(-vertical_scale, vertical_scale)
                c = rng.uniform(-vertical_offset, vertical_offset)
                r = (1.0 + s) * r + c
            if use_noise:
                sd = rng.uniform(*noise_std_range)
                if sd > 0:
                    r = r + rng.normal(0.0, sd, size=r.shape)

            Xs.append(r)
            ds.append(d)

    return np.stack(Xs), np.asarray(ds)


def measured_reference(wavelength_nm, true_thickness=1490.0, noise_std=0.004, seed=0,
                       n1=N_SIO2):
    """1~4편과 동일한 '측정 데이터' 생성 규약.

    _code/optimization-global-heuristics/generate_figures.py 의 블록을 그대로 옮겼다.
    np.random.seed(0) 기반 legacy RandomState 를 써야 1~4편과 같은 배열이 나온다.
    """
    np.random.seed(seed)
    R_true = reflectance(true_thickness, wavelength_nm, n1=n1)
    return R_true + np.random.normal(0, noise_std, size=wavelength_nm.shape)

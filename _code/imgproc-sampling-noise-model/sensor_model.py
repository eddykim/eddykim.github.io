"""합성 이미지 센서 — 광자에서 DN 까지의 사슬을 그대로 따라간다.

이미지처리 1편에서 쓰는 모델이다. 실제 카메라에서 픽셀 값이 만들어지는 순서를
그대로 코드로 옮겨두면, 뒤에서 광전달곡선(photon transfer curve)이나 프레임
평균의 한계를 실험할 때 "정답"을 알고 있는 상태로 검증할 수 있다.

사슬:
    광자 -> (양자효율) -> 신호 전자 -> (+암전류 전자, +고정패턴)
         -> (+읽기 노이즈) -> (시스템 이득 K) -> DN -> (양자화, 포화)

노이즈원 세 가지를 구분해서 넣는다.
  - 샷 노이즈(shot noise): 광자 도착이 포아송이라 신호 자체에 붙는다. 분산 = 평균.
  - 읽기 노이즈(read noise): 신호와 무관한 가우시안. 판독 회로에서 온다.
  - 고정패턴 노이즈(fixed pattern noise, FPN): 픽셀마다 고정된 편차라
    시간 평균으로 줄지 않는다. 곱셈형(PRNU)과 덧셈형(DSNU)으로 나눈다.
"""
import numpy as np

# EMVA 1288 에 맞춘 기호. 전자 단위는 e-, 디지털 단위는 DN 으로 적는다.
DEFAULTS = dict(
    qe=0.60,            # 양자효율 [e-/photon]
    gain_K=0.25,        # 시스템 이득 [DN/e-]
    read_noise_e=3.0,   # 읽기 노이즈 [e- rms]
    dark_rate_e=5.0,    # 암전류 [e-/s]
    full_well_e=15000,  # 최대 전하 용량 [e-]
    bit_depth=12,
    offset_dn=100.0,    # 검은 영역이 0 에 붙지 않도록 띄워두는 값
    prnu=0.010,         # 곱셈형 고정패턴 (이득 편차) [비율]
    dsnu_e=2.0,         # 덧셈형 고정패턴 (암전류 편차) [e- rms]
)


class Sensor:
    """픽셀 하나하나가 독립인 단순 센서. 크로스토크와 블루밍은 다루지 않는다."""

    def __init__(self, shape=(256, 256), seed=0, **kwargs):
        p = {**DEFAULTS, **kwargs}
        for k, v in p.items():
            setattr(self, k, v)
        self.shape = shape
        self.max_dn = 2 ** self.bit_depth - 1

        # 고정패턴은 센서를 만들 때 한 번만 뽑는다. 이후 몇 장을 찍든 그대로다.
        # 이 한 줄이 "평균으로는 못 지우는 노이즈"의 정체다.
        fp_rng = np.random.default_rng(seed)
        self.gain_map = 1.0 + self.prnu * fp_rng.standard_normal(shape)
        self.dark_map = self.dsnu_e * fp_rng.standard_normal(shape)

    def capture(self, photons, exposure_s, rng, fpn=True, quantize=True):
        """광자 수 photons(픽셀당 평균)를 받아 DN 배열 한 장을 내놓는다."""
        mean_e = self.qe * photons * (self.gain_map if fpn else 1.0)
        mean_e = np.broadcast_to(np.asarray(mean_e, float), self.shape)
        signal_e = rng.poisson(np.maximum(mean_e, 0.0)).astype(float)

        mean_dark = self.dark_rate_e * exposure_s
        dark_e = rng.poisson(np.full(self.shape, max(mean_dark, 0.0))).astype(float)
        if fpn:
            dark_e = dark_e + self.dark_map

        read_e = self.read_noise_e * rng.standard_normal(self.shape)

        total_e = np.clip(signal_e + dark_e + read_e, 0.0, self.full_well_e)
        dn = self.gain_K * total_e + self.offset_dn
        if quantize:
            dn = np.floor(dn)
        return np.clip(dn, 0.0, self.max_dn)

    def expected_variance_dn(self, mean_dn):
        """이론 분산 곡선. 광전달곡선 직선의 정답에 해당한다.

        sigma_y^2 = K^2 * sigma_dark^2 + K * (mu_y - y_0)
        """
        sigma_dark_sq = self.read_noise_e ** 2
        return (self.gain_K ** 2) * sigma_dark_sq + self.gain_K * (mean_dn - self.offset_dn)


def photon_transfer_curve(sensor, photon_levels, exposure_s, rng, n_pairs=1):
    """노출을 바꿔가며 (평균 DN, 분산 DN^2) 쌍을 모은다.

    분산은 한 장이 아니라 **두 장의 차이**에서 구한다. 같은 고정패턴이 두 장에
    똑같이 실리므로 빼면 사라지고, 시간적으로 변하는 노이즈만 남는다. 차이의
    분산은 원래의 두 배이므로 2 로 나눈다 (EMVA 1288 의 표준 절차다).
    """
    means, variances = [], []
    for p in photon_levels:
        acc_mean, acc_var = [], []
        for _ in range(n_pairs):
            a = sensor.capture(p, exposure_s, rng)
            b = sensor.capture(p, exposure_s, rng)
            acc_mean.append(0.5 * (a.mean() + b.mean()))
            acc_var.append(0.5 * (a - b).var())
        means.append(np.mean(acc_mean))
        variances.append(np.mean(acc_var))
    return np.array(means), np.array(variances)


def fit_ptc(signal_dn, var_dn, fit_range, sigma_q_sq=0.0):
    """광전달곡선의 선형 구간에 직선을 맞춰 K 와 어두운 영역 노이즈를 되찾는다.

    signal_dn 은 오프셋을 뺀 값(mu_y - y_0)이어야 한다. 그래야 기울기가 K [DN/e-],
    절편이 K^2 * sigma_dark^2 로 곧장 읽힌다. 최소자승 직선 하나로 데이터시트 없이
    센서의 두 상수를 얻는 셈이다.

    절편은 기울기 항에 비해 몇 자릿수 작으므로, 맞춤 구간을 신호가 작은 쪽으로
    잡아야 절편이 분해된다. 밝은 쪽까지 통째로 맞추면 절편은 잡음에 묻힌다.

    엄밀히는 절편에 양자화 노이즈 sigma_q^2 도 함께 들어 있다. 어두운 영역 노이즈만
    떼어내려면 이 값을 빼고 환산해야 한다. 균등 분포를 가정하면 sigma_q^2 = 1/12 DN^2
    이고, 절편 자체가 1 DN^2 이하인 저노이즈 센서에서는 무시할 수 없는 크기다.
    """
    lo, hi = fit_range
    m = (signal_dn >= lo) & (signal_dn <= hi)
    slope, intercept = np.polyfit(signal_dn[m], var_dn[m], 1)
    K = slope
    sigma_dark_e = np.sqrt(max(intercept - sigma_q_sq, 0.0)) / K
    return K, sigma_dark_e, slope, intercept

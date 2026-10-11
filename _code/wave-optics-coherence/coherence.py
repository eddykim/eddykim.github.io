"""파동광학 배경이론 2편의 계산 모듈: 무작위 빛, 결맞음도, 간섭무늬와 푸리에 변환 분광, 박막·OCT, 공간 결맞음, 스페클.

길이는 µm, 파수 σ = 1/λ 는 µm⁻¹ 로 쓴다. 광로차 Δ 에서 위상은 2πσΔ 다.
문헌의 닫힌 식(λ²/Δλ, 2J1(x)/x, 1.22λ/θ, 1/√N, OCT 분해능 식)은 여기서 쓰지 않는다.
verify_coherence.py 가 그것들과 대조할 독립 경로여야 하기 때문이다.
"""
import numpy as np


# ---------------------------------------------------------------- 스펙트럼

def gaussian_line(sigma, sigma0, fwhm):
    """파수 sigma0, 반치폭 fwhm(파수) 인 가우스 선 모양 (넓이 1)."""
    s = fwhm / (2 * np.sqrt(2 * np.log(2)))
    return np.exp(-0.5 * ((sigma - sigma0) / s) ** 2) / (s * np.sqrt(2 * np.pi))


def lorentz_line(sigma, sigma0, fwhm):
    g = fwhm / 2
    return (g / np.pi) / ((sigma - sigma0) ** 2 + g**2)


def wavelength_band(lam0, dlam_fwhm):
    """중심 파장 lam0, 파장 반치폭 dlam 을 파수 (sigma0, 반치폭) 으로 (좁은 띠 근사 없이 양 끝으로 정한다)."""
    s_hi, s_lo = 1 / (lam0 - dlam_fwhm / 2), 1 / (lam0 + dlam_fwhm / 2)
    return 1 / lam0, s_hi - s_lo


# ---------------------------------------------------------------- 무작위 빛

def random_wave(S_of_nu, nu, n_t, dt, rng):
    """파워 스펙트럼 S(ν) 를 갖는 정상 가우스 무작위 복소 파동 U(t) 한 표본.
    주파수 성분마다 독립인 복소 가우스 진폭을 주고 역 FFT 한다 (각 성분의 위상이 무작위)."""
    freqs = np.fft.fftfreq(n_t, dt)
    S = np.interp(freqs, nu, S_of_nu, left=0.0, right=0.0)
    amp = np.sqrt(S / 2) * (rng.standard_normal(n_t) + 1j * rng.standard_normal(n_t))
    return np.fft.ifft(amp) * n_t


def autocorrelation(U, max_lag):
    """시간 평균으로 잰 정규화 자기상관 g(τ) = <U*(t) U(t+τ)> / <|U|²> (원형 FFT, 표본이 충분히 길다고 본다)."""
    F = np.fft.fft(U)
    G = np.fft.ifft(np.abs(F) ** 2) / len(U)
    return G[:max_lag] / G[0]


def degree_from_spectrum(sigma, S, delays):
    """스펙트럼에서 복소 결맞음도 g(Δ) = ∫ S(σ) e^{i2πσΔ} dσ / ∫ S dσ (광로차 Δ 에서)."""
    w = S / np.trapezoid(S, sigma)
    d = np.atleast_1d(np.asarray(delays, float))
    out = np.empty(len(d), complex)
    for i in range(0, len(d), 256):                       # 지연 256 개씩 행렬로 계산
        blk = d[i:i + 256]
        out[i:i + 256] = np.trapezoid(w[None, :] * np.exp(2j * np.pi * blk[:, None] * sigma[None, :]), sigma, axis=1)
    return out


# ---------------------------------------------------------------- 간섭무늬와 분광

def michelson_interferogram(sigma, S, path_diff):
    """스펙트럼 S(σ) 인 빛의 마이컬슨 출력 (입력 세기 1로 정규화): I(Δ) = ½ ∫ s(σ) [1 + cos 2πσΔ] dσ."""
    w = S / np.trapezoid(S, sigma)
    return np.array([0.5 * np.trapezoid(w * (1 + np.cos(2 * np.pi * sigma * d)), sigma) for d in np.atleast_1d(path_diff)])


def envelope_visibility(path_diff, I, period_samples):
    """간섭무늬에서 국소 가시도: 한 무늬 주기 남짓의 창 안 (최대 − 최소)/(최대 + 최소)."""
    h = max(int(period_samples), 1)
    V = np.empty_like(I)
    for i in range(len(I)):
        seg = I[max(0, i - h):i + h + 1]
        V[i] = (seg.max() - seg.min()) / (seg.max() + seg.min())
    return V


def spectrum_from_interferogram(path_diff, I, sigma):
    """푸리에 변환 분광: 대칭으로 잰 간섭무늬의 교류 성분을 코사인 변환해 스펙트럼을 되찾는다."""
    ac = I - I.mean()
    return np.array([np.trapezoid(ac * np.cos(2 * np.pi * s * path_diff), path_diff) for s in sigma])


# ---------------------------------------------------------------- 박막과 OCT

def film_reflectance_two_beam(sigma, n, t, r1=None, r2=None):
    """공기 중 막(굴절률 n, 두께 t)의 두 빔 반사율. 아랫면 반사에 π 가 붙는다."""
    r = (1 - n) / (1 + n) if r1 is None else r1
    rb = -r if r2 is None else r2
    T = 1 - r**2
    U = r + T * rb * np.exp(1j * 4 * np.pi * n * t * sigma)
    return np.abs(U) ** 2


def observed_film_reflectance(n, t, sigma_grid, S):
    """스펙트럼 S 로 비춘 막을 파장을 가리지 않는 검출기로 볼 때의 반사율 (파장별 반사율의 가중 평균)."""
    w = S / np.trapezoid(S, sigma_grid)
    return np.array([np.trapezoid(w * film_reflectance_two_beam(sigma_grid, n, tt), sigma_grid) for tt in np.atleast_1d(t)])


def oct_ascan(sigma, S, mirror_z, layers):
    """시간 영역 OCT: 기준 거울 위치 z 를 움직이며 잰 세기. layers = [(깊이(공기 환산), 진폭 반사계수)].
    광로차는 왕복이라 2(z − z_i)."""
    w = S / np.trapezoid(S, sigma)
    out = []
    for z in np.atleast_1d(mirror_z):
        U = np.exp(1j * 4 * np.pi * sigma * z)
        for zi, ri in layers:
            U = U + ri * np.exp(1j * 4 * np.pi * sigma * zi)
        out.append(np.trapezoid(w * np.abs(U) ** 2, sigma))
    return np.array(out)


# ---------------------------------------------------------------- 공간 결맞음

def young_extended_source(x, d, L, lam, theta_s, n_pts=4001, disk=True):
    """각 지름 theta_s 인 고른 광원(원판 또는 띠)의 서로 독립인 점들이 각각 만든 영의 무늬를 세기로 더한다.
    광원 점이 각 θ 에서 오면 두 구멍 사이 위상차에 2π d θ / λ 가 더해진다."""
    th = np.linspace(-theta_s / 2, theta_s / 2, n_pts)
    w = np.sqrt(np.maximum(1 - (2 * th / theta_s) ** 2, 0)) if disk else np.ones_like(th)
    w /= w.sum()
    X = np.atleast_1d(x)[:, None]
    phase = 2 * np.pi * (d * X / (lam * L) + d * th[None, :] / lam)
    return (w[None, :] * (1 + np.cos(phase))).sum(1)


def visibility_extended_source(d, lam, theta_s, n_pts=4001, disk=True):
    """위 무늬의 가시도를 한 주기 안에서 직접 잰다."""
    x = np.linspace(0, lam * 1e6 / d, 400, endpoint=False)
    I = young_extended_source(x, d, 1e6, lam, theta_s, n_pts, disk)
    return (I.max() - I.min()) / (I.max() + I.min())


# ---------------------------------------------------------------- 스페클

def speckle_pattern(n, aperture_frac, rng, roughness_waves=3.0):
    """거친 면(무작위 높이)에서 흩어진 결맞는 빛을 원형 동공의 렌즈로 본 상의 세기.
    거친 면 → 무작위 위상, 렌즈의 동공 → 공간 주파수 저역 통과."""
    phase = 2 * np.pi * roughness_waves * rng.standard_normal((n, n))
    field = np.exp(1j * phase)
    F = np.fft.fft2(field)
    fx = np.fft.fftfreq(n)
    FX, FY = np.meshgrid(fx, fx)
    F[FX**2 + FY**2 > (aperture_frac / 2) ** 2] = 0
    I = np.abs(np.fft.ifft2(F)) ** 2
    return I / I.mean()


def contrast(I):
    return I.std() / I.mean()

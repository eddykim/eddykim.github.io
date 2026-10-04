"""파동광학 배경이론 1편의 계산 모듈: 두 파동의 간섭, 점광원 간섭장, 박막·뉴턴 링, 빔스플리터, 파브리-페로.

길이 단위는 µm 다. 복소 진폭은 exp(+i(kr − ωt)) 규약으로 시간 인자를 생략해 쓴다.
문헌의 닫힌 식(무늬 간격 λL/d, 뉴턴 링 반지름 √(mλR), 에어리 함수, 피네스 π√R/(1−R))은
여기서 쓰지 않는다. verify_interference.py 가 그것들과 대조할 독립 경로여야 하기 때문이다.
"""
import numpy as np


# ---------------------------------------------------------------- 두 파동

def superpose(I1, I2, delta):
    """세기 I1, I2, 위상차 delta 인 두 파동을 진폭으로 더한 뒤 세기를 잰다."""
    U = np.sqrt(I1) + np.sqrt(I2) * np.exp(1j * np.asarray(delta))
    return np.abs(U) ** 2


def visibility_of(I):
    """무늬 세기 배열에서 가시도 (Imax − Imin)/(Imax + Imin)."""
    return (I.max() - I.min()) / (I.max() + I.min())


# ---------------------------------------------------------------- 점광원

def point_sources_intensity(X, Z, sources, wavelength, amplitudes=None, phases=None):
    """같은 진동수의 점광원들이 평면 (x, z) 위에 만드는 세기. 구면파 e^{ikr}/r 를 그대로 더한다."""
    k = 2 * np.pi / wavelength
    U = np.zeros_like(X, dtype=complex)
    for j, (xs, zs) in enumerate(sources):
        r = np.hypot(X - xs, Z - zs)
        a = 1.0 if amplitudes is None else amplitudes[j]
        p = 0.0 if phases is None else phases[j]
        U += a * np.exp(1j * (k * r + p)) / r
    return np.abs(U) ** 2


def young_screen(x, d, L, wavelength):
    """간격 d 의 두 바늘구멍(광축에 대칭)에서 거리 L 의 스크린 위 세기. 근사 없이 두 구면파를 더한다."""
    X, Z = np.asarray(x, float), np.full_like(np.asarray(x, float), L)
    return point_sources_intensity(X, Z, [(-d / 2, 0.0), (d / 2, 0.0)], wavelength)


# ---------------------------------------------------------------- 박막

def fresnel_r(n1, n2):
    """수직 입사 진폭 반사계수 (n1 → n2)."""
    return (n1 - n2) / (n1 + n2)


def thin_layer_reflectance(n0, n1, n2, t, wavelength, multiple=True):
    """두께 t, 굴절률 n1 인 층(위 n0, 아래 n2)의 수직 입사 반사율.
    multiple=False 면 첫 두 반사만 더한 두 빔 근사, True 면 층 안의 다중반사를 모두 더한 값 (급수를 직접 합산)."""
    r01, r12 = fresnel_r(n0, n1), fresnel_r(n1, n2)
    t01, t10 = 1 + r01, 1 - r01
    r10 = -r01
    phase = np.exp(1j * 4 * np.pi * n1 * np.asarray(t, float) / wavelength)
    if not multiple:
        return np.abs(r01 + t01 * t10 * r12 * phase) ** 2
    total = r01 + 0j
    term = t01 * t10 * r12 * phase
    for _ in range(400):                     # 층 안을 한 번 오갈 때마다 r10·r12·phase 가 곱해진다
        total = total + term
        term = term * r10 * r12 * phase
    return np.abs(total) ** 2


def newton_gap(r, R):
    """곡률반지름 R 인 볼록면이 평판에 닿아 있을 때 반지름 r 에서의 공기층 두께 (정확한 구면)."""
    return R - np.sqrt(R**2 - r**2)


# ---------------------------------------------------------------- 빔스플리터와 마이컬슨

def beamsplitter(t=1 / np.sqrt(2)):
    """손실 없는 대칭 빔스플리터의 산란 행렬. 에너지 보존(유니타리)이 반사 위상 i 를 강제한다."""
    r = 1j * np.sqrt(1 - t**2)
    return np.array([[t, r], [r, t]])


def michelson_outputs(delta, t=1 / np.sqrt(2)):
    """마이컬슨: 빔스플리터로 나눈 두 빔이 위상차 delta 를 얻고 같은 빔스플리터로 돌아와 합쳐진다.
    반환: (검출기 쪽 출력 세기, 광원 쪽으로 되돌아가는 출력 세기), 입력 세기 1."""
    B = beamsplitter(t)
    delta = np.atleast_1d(np.asarray(delta, float))
    out_a, out_b = [], []
    for d in delta:
        arms = B @ np.array([1.0, 0.0])              # 두 팔로 나뉜다
        arms = arms * np.array([1.0, np.exp(1j * d)])  # 둘째 팔이 위상 d 만큼 더 간다
        out = B @ arms                                # 다시 합쳐진다
        out_a.append(abs(out[1]) ** 2)
        out_b.append(abs(out[0]) ** 2)
    return np.array(out_a), np.array(out_b)


# ---------------------------------------------------------------- 여러 빔

def multibeam_sum(phi, h, M=None):
    """진폭이 |h| 배씩 줄고 위상이 phi 씩 늘어나는 파동들의 합의 세기 (첫 파동 세기 1).
    M 개만 더하거나(M 지정) 수렴할 때까지 더한다."""
    phi = np.asarray(phi, float)
    n = M if M is not None else int(np.ceil(np.log(1e-12) / np.log(abs(h)))) + 1
    U = np.zeros_like(phi, dtype=complex)
    term = np.ones_like(phi, dtype=complex)
    for _ in range(n):
        U += term
        term = term * h * np.exp(1j * phi)
    return np.abs(U) ** 2


def etalon_transmission(phi, R, M=None):
    """반사율 R 인 두 거울 사이(손실 없음)를 오가며 새어 나온 빔들의 합: 투과율.
    한 번 왕복할 때 진폭 R, 위상 phi. 첫 빔의 투과 진폭은 (1 − R)."""
    return (1 - R) ** 2 * multibeam_sum(phi, R, M)


def fwhm_of_peak(phi, I):
    """phi = 0 근처 봉우리의 반치폭 (수치)."""
    half = I.max() / 2
    i0 = np.argmax(I)
    right = i0 + np.argmax(I[i0:] < half)
    left = i0 - np.argmax(I[:i0 + 1][::-1] < half)
    def interp(i_in, i_out):
        return np.interp(half, [I[i_out], I[i_in]], [phi[i_out], phi[i_in]])
    return interp(right - 1, right) - interp(left + 1, left)


# ---------------------------------------------------------------- 전달행렬 (검증용 독립 경로)

def tmm_normal(ns, ds, wavelength):
    """수직 입사 전달행렬법: 굴절률 ns = [n0, n1, ..., nN], 두께 ds = [d1, ..., d_{N-1}] → (r, t)."""
    M = np.eye(2, dtype=complex)
    for j in range(len(ns) - 1):
        n1, n2 = ns[j], ns[j + 1]
        D = 0.5 / n2 * np.array([[n2 + n1, n2 - n1], [n2 - n1, n2 + n1]])   # 경계면
        M = D @ M
        if j + 1 < len(ns) - 1:
            b = 2 * np.pi * n2 * ds[j] / wavelength
            M = np.diag([np.exp(1j * b), np.exp(-1j * b)]) @ M
    r = -M[1, 0] / M[1, 1]
    t = M[0, 0] + M[0, 1] * r
    return r, t


# ---------------------------------------------------------------- 수차의 간섭무늬

def interferogram(W_waves, tilt_waves=0.0, X=None):
    """파면수차 W(파장 단위)를 평면 기준파와 겹친 간섭무늬. tilt 는 동공 가로 방향 기울기(파장/반지름)."""
    tilt = 0.0 if X is None else tilt_waves * X
    return 0.5 * (1 + np.cos(2 * np.pi * (W_waves + tilt)))

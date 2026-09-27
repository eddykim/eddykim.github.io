"""후초점면 마이크로 타원계측기의 신호 합성과 고리 복조.

고개구수 대물렌즈로 수직 입사 광학계를 꾸미면 후초점면에서 두 가지가 동시에 일어난다.

  반경  r     -> 입사각   theta = arcsin(r/r_max * NA)      (아베 사인조건)
  방위각 phi  -> 편광축 회전. 입사면이 phi 만큼 돌아간 것과 같다

그래서 같은 반경의 고리를 따라 방위각별 세기를 읽으면 부품을 돌리지 않고도 변조가
일어난다. 다만 편광자와 분석기가 둘 다 실험실 좌표계에 고정돼 있으므로, 시편 좌표계
에서 보면 **둘이 함께** 도는 배치다 — 3편이 분류한 "편광자·분석기 동시 회전"이다.
따라서 신호에 2phi 와 4phi 성분이 함께 실린다(1편의 RAE 는 2omega 뿐이다).
반경을 바꾸면 다른 입사각이 나온다.

신호 합성은 부품 뮬러 행렬의 곱으로만 한다. Ye 2007 의 닫힌 식과 문헌 수치는
verify_bfp.py 가 대조할 독립 경로로 남겨둔다.

부호 규약은 배경이론 2·3편, 타원계측기 1~5편과 일치시킨다.
"""
import numpy as np

D = np.deg2rad

# 본문에 쓰는 기준 시편과 대물렌즈
N_SIO2, N_SI = 1.460, 3.875 - 0.0156j      # Ye 2007 의 값
N_FUSED_SILICA = 1.458                      # 김영준 2025 4.3 절의 기준 시편
NA_DEFAULT = 0.9


# ---------------------------------------------------------------------------
# 뮬러 계산법 기본 요소 (배경이론 2편, 타원계측기 3~5편과 동일한 정의)
# ---------------------------------------------------------------------------


def mueller_rotation(omega):
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])


def mueller_polarizer(theta):
    base = 0.5 * np.array([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
    return mueller_rotation(-theta) @ base @ mueller_rotation(theta)


def mueller_sample(psi, delta):
    """등방성 시편. 3~5편의 mueller_sample 과 같다."""
    c2, s2 = np.cos(2 * psi), np.sin(2 * psi)
    cd, sd = np.cos(delta), np.sin(delta)
    return np.array([[1, -c2, 0, 0], [-c2, 1, 0, 0],
                     [0, 0, s2 * cd, s2 * sd], [0, 0, -s2 * sd, s2 * cd]])


S_UNPOLARIZED = np.array([1.0, 0.0, 0.0, 0.0])


# ---------------------------------------------------------------------------
# 시편 — 프레넬 반사와 단층막
# ---------------------------------------------------------------------------


def fresnel_rp_rs(theta, n1, n2):
    """매질 n1 에서 n2 로 입사할 때의 진폭 반사계수. theta 는 입사각(rad)."""
    theta = np.asarray(theta, dtype=complex)
    c1 = np.cos(theta)
    c2 = np.sqrt(1.0 - (n1 * np.sin(theta) / n2) ** 2)
    rp = (n2 * c1 - n1 * c2) / (n2 * c1 + n1 * c2)
    rs = (n1 * c1 - n2 * c2) / (n1 * c1 + n2 * c2)
    return rp, rs


def psi_delta_bare(theta, n_sub, n_amb=1.0):
    """맨 기판의 Psi, Delta."""
    rp, rs = fresnel_rp_rs(theta, n_amb, n_sub)
    rho = rp / rs
    return np.arctan(np.abs(rho)), np.angle(rho)


def psi_delta_film(theta, d_nm, lam_nm=632.8, n_film=N_SIO2, n_sub=N_SI, n_amb=1.0):
    """기판 위 단층막의 Psi, Delta. 5편의 film_psi_delta 와 같은 모형이다."""
    theta = np.asarray(theta, dtype=complex)
    c1 = np.sqrt(1.0 - (n_amb * np.sin(theta) / n_film) ** 2)
    r01p, r01s = fresnel_rp_rs(theta, n_amb, n_film)
    th1 = np.arcsin(n_amb * np.sin(theta) / n_film)
    r12p, r12s = fresnel_rp_rs(th1, n_film, n_sub)
    beta = 2 * np.pi * (d_nm / lam_nm) * n_film * c1
    z = np.exp(-2j * beta)
    rp = (r01p + r12p * z) / (1.0 + r01p * r12p * z)
    rs = (r01s + r12s * z) / (1.0 + r01s * r12s * z)
    rho = rp / rs
    return np.arctan(np.abs(rho)), np.angle(rho)


def brewster_angle(n_sub, n_amb=1.0):
    """유전체의 브루스터각. 흡수가 있으면 |rp| 최소점을 수치로 찾는다."""
    if np.imag(n_sub) == 0:
        return np.arctan(np.real(n_sub) / n_amb)
    th = np.linspace(1e-4, np.pi / 2 - 1e-4, 200001)
    rp, _ = fresnel_rp_rs(th, n_amb, n_sub)
    return th[int(np.argmin(np.abs(rp)))]


# ---------------------------------------------------------------------------
# 후초점면 — 반경이 입사각, 방위각이 편광축 회전
# ---------------------------------------------------------------------------


def radius_to_angle(r_norm, na=NA_DEFAULT, n_amb=1.0):
    """정규화 반경(0~1)을 입사각으로. 아베 사인조건이다."""
    return np.arcsin(np.clip(np.asarray(r_norm) * na / n_amb, -1.0, 1.0))


def angle_to_radius(theta, na=NA_DEFAULT, n_amb=1.0):
    """입사각을 정규화 반경으로 (위의 역)."""
    return n_amb * np.sin(theta) / na


def max_angle(na, n_amb=1.0):
    """대물렌즈가 닿을 수 있는 최대 입사각."""
    return np.arcsin(np.clip(na / n_amb, -1.0, 1.0))


def bfp_intensity(psi, delta, phi, pol=0.0, ana=np.pi / 2):
    """후초점면 한 점의 검출 세기.

    방위각 phi 인 점에서는 입사면이 phi 만큼 돌아가 있다. 곧 편광자와 분석기가
    시편 좌표계에서 -phi 만큼 회전한 것과 같다. 부품을 돌리는 대신 좌표가 돈다.

    기본값은 편광자 0도, 분석기 90도(교차 편광)다. Ye 2007 은 둘 다 x 축에
    맞췄으므로 ana=0 으로 두면 그 배치가 된다.
    """
    phi = np.asarray(phi, dtype=float)
    ms = mueller_sample(psi, delta)
    out = np.empty_like(phi)
    for i, p in enumerate(np.atleast_1d(phi).ravel()):
        mp = mueller_polarizer(pol - p)
        ma = mueller_polarizer(ana - p)
        out.ravel()[i] = (ma @ (ms @ (mp @ S_UNPOLARIZED)))[0]
    return out


def bfp_map(d_nm, na=NA_DEFAULT, n_pix=241, pol=0.0, ana=0.0, lam_nm=632.8,
            n_film=N_SIO2, n_sub=N_SI):
    """후초점면 전체의 세기 지도. 반환: (x, y, I). 동공 밖은 NaN."""
    ax = np.linspace(-1.0, 1.0, n_pix)
    X, Y = np.meshgrid(ax, ax)
    R = np.hypot(X, Y)
    PHI = np.arctan2(Y, X)
    I = np.full_like(R, np.nan)
    inside = R <= 1.0
    th = radius_to_angle(R[inside], na)
    psi, dl = psi_delta_film(th, d_nm, lam_nm, n_film, n_sub)
    for k, (t_psi, t_dl, t_phi) in enumerate(zip(psi, dl, PHI[inside])):
        I[np.nonzero(inside)[0][k], np.nonzero(inside)[1][k]] = bfp_intensity(
            t_psi, t_dl, np.array([t_phi]), pol, ana)[0]
    return X, Y, I


def annulus(d_nm, r_norm, na=NA_DEFAULT, n_phi=360, pol=0.0, ana=0.0,
            lam_nm=632.8, n_film=N_SIO2, n_sub=N_SI):
    """반경 r_norm 인 고리를 따라 방위각별 세기를 읽는다.

    반환: (phi 배열, 세기 배열, 그 반경의 (Psi, Delta))
    """
    phi = np.linspace(0.0, 2 * np.pi, n_phi, endpoint=False)
    th = radius_to_angle(r_norm, na)
    psi, dl = psi_delta_film(np.array([th]), d_nm, lam_nm, n_film, n_sub)
    return phi, bfp_intensity(psi[0], dl[0], phi, pol, ana), (psi[0], dl[0])


# ---------------------------------------------------------------------------
# 고리 복조 — 회전편광자 타원계측기와 등가다
# ---------------------------------------------------------------------------


def fourier_alpha(signal, phi):
    """I = dc[1 + a2 cos 2phi + a4 cos 4phi] 의 정규화 계수 (a2, a4).

    Ye 2007 이 {alpha_2, alpha_4} 라 부르는 값이다. 편광자와 분석기가 함께 도는
    배치라 사인 성분이 항등적으로 0 이고 코사인 두 개만 남는다.
    """
    dc = np.mean(signal)
    return (2 * np.mean(signal * np.cos(2 * phi)) / dc,
            2 * np.mean(signal * np.cos(4 * phi)) / dc)


def invert_annulus(a2, a4):
    """고리의 푸리에 계수에서 Psi, Delta 를 읽는다.

    뮬러 곱을 전개하면 (verify_bfp.py 가 대조한다)

      I ∝ (3 + sin2Psi cosDelta)/2 - 2 cos2Psi cos2phi
          + (1 - sin2Psi cosDelta)/2 cos4phi

    이고, 정규화 계수로 풀면 아래가 된다.

      cos 2Psi        = -a2 / (1 + a4)
      sin 2Psi cosDelta = (1 - 3 a4) / (1 + a4)

    RAE 와 마찬가지로 cos(Delta) 만 얻으므로 Delta 의 부호는 잃는다.
    """
    c2 = np.clip(-np.asarray(a2) / (1.0 + np.asarray(a4)), -1.0, 1.0)
    psi = 0.5 * np.arccos(c2)
    s2 = np.sin(2 * psi)
    u = (1.0 - 3.0 * np.asarray(a4)) / (1.0 + np.asarray(a4))
    cd = np.clip(u / np.where(np.abs(s2) < 1e-12, np.nan, s2), -1.0, 1.0)
    return psi, np.arccos(cd)


# ---------------------------------------------------------------------------
# 집속이 만드는 평균 — 평면파가 아니다
# ---------------------------------------------------------------------------


def jones_sample(psi, delta):
    """등방성 시편의 존스 행렬(정규화). diag(rho, 1) 꼴이다."""
    return np.array([[np.tan(psi) * np.exp(1j * delta), 0.0], [0.0, 1.0]])


def jones_to_mueller(J):
    """존스 행렬을 뮬러 행렬로. 편광 소멸이 없는 행렬만 나온다."""
    A = np.array([[1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1j, -1j, 0]])
    return np.real(A @ np.kron(J, np.conj(J)) @ np.linalg.inv(A))


def averaged_mueller(theta_c, half_width, d_nm, n_ray=201, coherent=False,
                     lam_nm=632.8, n_film=N_SIO2, n_sub=N_SI):
    """입사각 폭에 걸쳐 평균한 시편 행렬.

    coherent=True  : 존스 행렬을 평균한 뒤 뮬러로 (공초점 검출)
    coherent=False : 뮬러 행렬을 평균 (통상 검출)

    Munro & Torok 2008 이 가른 두 경우다. 뮬러 평균은 편광 소멸을 만들고
    존스 평균은 만들지 않는다 — verify 가 확인한다.
    """
    th = np.linspace(theta_c - half_width, theta_c + half_width, n_ray)
    th = th[(th > 1e-6) & (th < np.pi / 2 - 1e-6)]
    psi, dl = psi_delta_film(th, d_nm, lam_nm, n_film, n_sub)
    if coherent:
        J = np.mean([jones_sample(p, d) for p, d in zip(psi, dl)], axis=0)
        return jones_to_mueller(J) / jones_to_mueller(J)[0, 0]
    M = np.mean([jones_to_mueller(jones_sample(p, d)) for p, d in zip(psi, dl)],
                axis=0)
    return M / M[0, 0]


def depolarization_index(M):
    """편광소멸 지수. 1 이면 편광을 유지하고 0 이면 완전히 소멸시킨다."""
    M = np.asarray(M, dtype=float)
    num = np.sum(M ** 2) - M[0, 0] ** 2
    return float(np.sqrt(num / (3.0 * M[0, 0] ** 2)))


# ---------------------------------------------------------------------------
# 입사각 눈금 보정 — 브루스터각을 기준점으로
# ---------------------------------------------------------------------------


def tan2_psi_profile(r_norm, na, n_sub, n_amb=1.0):
    """반경을 따라 본 tan^2(Psi) = Rp/Rs. 브루스터각에서 최소가 된다."""
    th = radius_to_angle(r_norm, na)
    psi, _ = psi_delta_bare(th, n_sub, n_amb)
    return np.tan(psi) ** 2


def calibrate_radius(r_norm, na_true, n_sub, n_amb=1.0, n_scan=4001):
    """tan^2(Psi) 의 최소점에서 브루스터 반경을 찾아 반경-각도 눈금을 세운다.

    반환: (r_b, theta_b, 그 눈금으로 환산한 입사각)
    식은 김영준 2025 식 4.2 다:  theta = arcsin(r sin(theta_b) / r_b)
    """
    rr = np.linspace(1e-4, 1.0, n_scan)
    prof = tan2_psi_profile(rr, na_true, n_sub, n_amb)
    r_b = rr[int(np.argmin(prof))]
    th_b = brewster_angle(n_sub, n_amb)
    theta = np.arcsin(np.clip(np.asarray(r_norm) * np.sin(th_b) / r_b, -1.0, 1.0))
    return r_b, th_b, theta

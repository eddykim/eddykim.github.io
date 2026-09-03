"""Levenberg-Marquardt법: (J^T J + mu I) h = -J^T r, gain ratio로 mu를 매 스텝 조정.

Madsen, Nielsen, Tingleff (2004) Algorithm 3.16 및 식 2.21을 그대로 구현했다.
gauss_newton.py와 마찬가지로 파라미터가 두께 d 하나뿐이라 J^T J, J^T r 모두 스칼라다.
mu=0으로 고정하면 gauss_newton.py의 스텝과 정확히 같아진다.
"""
import numpy as np
from reflectance_model import reflectance


def residual(thickness_nm, wavelength_nm, measured_R):
    return reflectance(thickness_nm, wavelength_nm) - measured_R


def objective(thickness_nm, wavelength_nm, measured_R):
    r = residual(thickness_nm, wavelength_nm, measured_R)
    return 0.5 * np.sum(r ** 2)


def numerical_jacobian(thickness_nm, wavelength_nm, measured_R, h=1e-3):
    r_plus = residual(thickness_nm + h, wavelength_nm, measured_R)
    r_minus = residual(thickness_nm - h, wavelength_nm, measured_R)
    return (r_plus - r_minus) / (2 * h)


def levenberg_marquardt(d0, wavelength_nm, measured_R, n_iter=30, h=1e-3,
                         tau=1e-3, eps1=1e-12, eps2=1e-12):
    """반환: d_hist, J_hist, mu_hist. mu_hist[k]는 k번째 스텝에 사용한 mu.

    tau는 초기 댐핑 mu0 = tau * (J^T J) 를 정하는 스케일 인자다(Madsen 3.2절).
    매 스텝 gain ratio rho = (F(d)-F(d+h)) / (L(0)-L(h)) 를 계산해
    rho>0 이면 스텝을 받아들이고 mu를 줄이며(식 2.21), rho<=0 이면 스텝을
    기각하고 mu를 nu배 늘린 뒤 nu를 2배로 키운다(Marquardt 1963의 flutter를
    피하기 위한 Nielsen 1999 권장 규칙).
    """
    d = d0
    r = residual(d, wavelength_nm, measured_R)
    Jr = numerical_jacobian(d, wavelength_nm, measured_R, h)
    A = np.sum(Jr ** 2)   # J^T J (파라미터 1개라 1x1 행렬 = 스칼라)
    g = np.sum(Jr * r)    # J^T r

    mu = tau * A
    nu = 2.0

    d_hist = [d]
    J_hist = [objective(d, wavelength_nm, measured_R)]
    mu_hist = []

    found = abs(g) <= eps1

    for _ in range(n_iter):
        if found:
            break
        mu_hist.append(mu)
        h_lm = -g / (A + mu)

        d_new = d + h_lm
        F_old = objective(d, wavelength_nm, measured_R)
        F_new = objective(d_new, wavelength_nm, measured_R)
        L0_minus_Lh = 0.5 * h_lm * (mu * h_lm - g)
        rho = (F_old - F_new) / L0_minus_Lh if L0_minus_Lh > 0 else -1.0

        if rho > 0:
            d = d_new
            r = residual(d, wavelength_nm, measured_R)
            Jr = numerical_jacobian(d, wavelength_nm, measured_R, h)
            A = np.sum(Jr ** 2)
            g = np.sum(Jr * r)
            found = (abs(g) <= eps1) or (abs(h_lm) <= eps2 * (abs(d) + eps2))
            mu = mu * max(1.0 / 3.0, 1.0 - (2 * rho - 1) ** 3)
            nu = 2.0
        else:
            mu = mu * nu
            nu = 2.0 * nu

        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))

    return np.array(d_hist), np.array(J_hist), np.array(mu_hist)

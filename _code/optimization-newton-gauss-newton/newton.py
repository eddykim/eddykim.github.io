"""목적함수 J(d)의 1차·2차 미분(수치미분)을 이용한 Newton법."""
import numpy as np
from reflectance_model import reflectance


def objective(thickness_nm, wavelength_nm, measured_R):
    """J(d) = 0.5 * sum((model(d) - measured)^2)"""
    model_R = reflectance(thickness_nm, wavelength_nm)
    residual = model_R - measured_R
    return 0.5 * np.sum(residual ** 2)


def numerical_grad_hess(thickness_nm, wavelength_nm, measured_R, h=1e-3):
    """중심차분으로 dJ/dd, d^2J/dd^2 를 근사."""
    J0 = objective(thickness_nm, wavelength_nm, measured_R)
    J_plus = objective(thickness_nm + h, wavelength_nm, measured_R)
    J_minus = objective(thickness_nm - h, wavelength_nm, measured_R)
    grad = (J_plus - J_minus) / (2 * h)
    hess = (J_plus - 2 * J0 + J_minus) / (h ** 2)
    return grad, hess


def newton(d0, wavelength_nm, measured_R, n_iter=100, h=1e-3):
    """반환: iteration별 두께 추정값, 목적함수 값."""
    d = d0
    d_hist = [d]
    J_hist = [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        grad, hess = numerical_grad_hess(d, wavelength_nm, measured_R, h)
        d = d - grad / hess
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)

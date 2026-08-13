"""최소자승 목적함수 + 수치미분 기반 gradient descent (두께 1개 파라미터)."""
import numpy as np
from reflectance_model import reflectance


def objective(thickness_nm, wavelength_nm, measured_R):
    """J(d) = 0.5 * sum((model(d) - measured)^2)"""
    model_R = reflectance(thickness_nm, wavelength_nm)
    residual = model_R - measured_R
    return 0.5 * np.sum(residual ** 2)


def numerical_gradient(thickness_nm, wavelength_nm, measured_R, h=1e-3):
    """중심차분으로 dJ/dd 근사."""
    J_plus = objective(thickness_nm + h, wavelength_nm, measured_R)
    J_minus = objective(thickness_nm - h, wavelength_nm, measured_R)
    return (J_plus - J_minus) / (2 * h)


def gradient_descent(d0, wavelength_nm, measured_R, alpha, n_iter=100):
    """반환: iteration별 두께 추정값, 목적함수 값."""
    d = d0
    d_hist = [d]
    J_hist = [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        grad = numerical_gradient(d, wavelength_nm, measured_R)
        d = d - alpha * grad
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)

"""잔차 벡터 r(d)의 Jacobian(수치미분)을 이용한 Gauss-Newton법.

Newton법(newton.py)은 J(d)의 2차미분 d^2J/dd^2 를 직접 수치미분으로 구하지만,
Gauss-Newton은 잔차 r_i(d) = model_i(d) - measured_i 의 1차미분(Jacobian)만으로
d^2J/dd^2 를 근사한다:

    d^2J/dd^2 = sum_i [ (dr_i/dd)^2 + r_i * d^2r_i/dd^2 ]  (Newton, 정확한 Hessian)
    d^2J/dd^2 ~= sum_i (dr_i/dd)^2                          (Gauss-Newton, 둘째 항을 버림)

파라미터가 두께 d 하나뿐이라 정규방정식 (J^T J) h = -J^T r 이 스칼라 나눗셈으로 풀린다.
"""
import numpy as np
from reflectance_model import reflectance


def residual(thickness_nm, wavelength_nm, measured_R):
    return reflectance(thickness_nm, wavelength_nm) - measured_R


def objective(thickness_nm, wavelength_nm, measured_R):
    r = residual(thickness_nm, wavelength_nm, measured_R)
    return 0.5 * np.sum(r ** 2)


def numerical_jacobian(thickness_nm, wavelength_nm, measured_R, h=1e-3):
    """중심차분으로 파장별 잔차 각각에 대한 dr_i/dd 를 근사."""
    r_plus = residual(thickness_nm + h, wavelength_nm, measured_R)
    r_minus = residual(thickness_nm - h, wavelength_nm, measured_R)
    return (r_plus - r_minus) / (2 * h)


def gauss_newton(d0, wavelength_nm, measured_R, n_iter=100, h=1e-3):
    """반환: iteration별 두께 추정값, 목적함수 값."""
    d = d0
    d_hist = [d]
    J_hist = [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        r = residual(d, wavelength_nm, measured_R)
        Jr = numerical_jacobian(d, wavelength_nm, measured_R, h)
        step = -np.sum(Jr * r) / np.sum(Jr ** 2)  # (J^T J) h = -J^T r
        d = d + step
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)

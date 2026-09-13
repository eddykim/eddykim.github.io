"""2파라미터 Levenberg-Marquardt — (두께 d, SiO2 굴절률 n1) 동시 피팅.

복사해 온 levenberg_marquardt.py 는 파라미터 1개 전용이다(A = sum(Jr**2) 가 스칼라).
실험 5-2 에서 파라미터 상호상관을 보려면 미지수가 최소 2개여야 하므로,
같은 Nielsen 댐핑 규칙(Madsen·Nielsen·Tingleff 2004)을 2x2 정규방정식으로 확장했다.
원본 파일은 손대지 않았다.
"""
import numpy as np

from reflectance_model import reflectance


def residual2(p, wl, R_meas):
    d, n1 = p
    return reflectance(d, wl, n1=n1) - R_meas


def objective2(p, wl, R_meas):
    r = residual2(p, wl, R_meas)
    return 0.5 * float(np.sum(r ** 2))


def jacobian2(p, wl, R_meas, h=(1e-3, 1e-6)):
    """중심차분 Jacobian (n_wl, 2). 열: [dr/dd, dr/dn1]"""
    cols = []
    for i, hi in enumerate(h):
        pp = np.array(p, dtype=float); pp[i] += hi
        pm = np.array(p, dtype=float); pm[i] -= hi
        cols.append((residual2(pp, wl, R_meas) - residual2(pm, wl, R_meas)) / (2 * hi))
    return np.column_stack(cols)


def lm2(p0, wl, R_meas, n_iter=60, tau=1e-3, eps1=1e-14, eps2=1e-14):
    """반환: (p_final, n_iter_used, J_final)"""
    p = np.array(p0, dtype=float)
    r = residual2(p, wl, R_meas)
    Jm = jacobian2(p, wl, R_meas)
    A = Jm.T @ Jm
    g = Jm.T @ r

    mu = tau * float(np.max(np.diag(A)))
    nu = 2.0
    used = 0

    for _ in range(n_iter):
        if np.linalg.norm(g, np.inf) <= eps1:
            break
        used += 1
        try:
            h_lm = np.linalg.solve(A + mu * np.eye(2), -g)
        except np.linalg.LinAlgError:
            mu *= nu; nu *= 2.0
            continue

        p_new = p + h_lm
        F_old = objective2(p, wl, R_meas)
        F_new = objective2(p_new, wl, R_meas)
        denom = 0.5 * float(h_lm @ (mu * h_lm - g))
        rho = (F_old - F_new) / denom if denom > 0 else -1.0

        if rho > 0:
            p = p_new
            r = residual2(p, wl, R_meas)
            Jm = jacobian2(p, wl, R_meas)
            A = Jm.T @ Jm
            g = Jm.T @ r
            mu *= max(1.0 / 3.0, 1.0 - (2 * rho - 1) ** 3)
            nu = 2.0
            if np.linalg.norm(h_lm) <= eps2 * (np.linalg.norm(p) + eps2):
                break
        else:
            mu *= nu
            nu *= 2.0

    return p, used, objective2(p, wl, R_meas)


def correlation_and_condition(p, wl, R_meas):
    """참 파라미터 근방의 파라미터 상호상관계수와 J^T J 조건수.

    공분산 C = (J^T J)^{-1} 에 대해 corr = C01 / sqrt(C00 * C11).
    Johnson 2000 (Methods Enzymol. 321, 424) 이 -0.999933 을 보인 그 값과 같은 정의다.
    """
    Jm = jacobian2(p, wl, R_meas)
    A = Jm.T @ Jm
    cond = float(np.linalg.cond(A))
    C = np.linalg.inv(A)
    corr = float(C[0, 1] / np.sqrt(C[0, 0] * C[1, 1]))
    return corr, cond

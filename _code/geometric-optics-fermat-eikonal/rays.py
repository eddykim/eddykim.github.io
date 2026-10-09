"""기하광학 배경이론 1편의 계산 모듈.

광선 방정식  d/ds (n dr/ds) = grad n  을 2차원(x 수평, y 수직)에서 적분하고,
광로길이(OPL), 신기루 대기 모델, 반사·굴절 코스틱을 계산한다.

그림(generate_figures.py)과 검증(verify_rays.py)이 모두 이 파일을 쓴다.
문헌의 닫힌 식은 여기에 쓰지 않는다. 검증 파일이 대조할 독립 경로여야 하기 때문이다.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp


# ── 굴절률 분포 ────────────────────────────────────────────
@dataclass
class Medium:
    """굴절률 n(x, y) 와 그 기울기를 돌려주는 매질."""

    n: callable        # (x, y) -> n
    grad: callable     # (x, y) -> (dn/dx, dn/dy)


def stratified(n_of_y, dn_dy):
    """굴절률이 높이 y 에만 의존하는 층상 매질."""
    return Medium(
        n=lambda x, y: n_of_y(y),
        grad=lambda x, y: (0.0 * np.asarray(x), dn_dy(y)),
    )


def n2_linear(a, b):
    """n^2 = a + b y 인 매질."""
    return stratified(lambda y: np.sqrt(a + b * y),
                      lambda y: 0.5 * b / np.sqrt(a + b * y))


def n2_parabolic(n0, alpha):
    """n^2 = n0^2 (1 - alpha^2 y^2) 인 매질 (GRIN 렌즈·광섬유의 모형)."""
    return stratified(lambda y: n0 * np.sqrt(1.0 - (alpha * y) ** 2),
                      lambda y: -n0 * alpha**2 * y / np.sqrt(1.0 - (alpha * y) ** 2))


# ── 광선 방정식 적분 ───────────────────────────────────────
def trace_ray(medium, r0, angle, s_max, stop=None, max_step=np.inf, rtol=1e-11, atol=1e-13):
    """시작점 r0 에서 수평과 angle [rad] 을 이루는 광선을 호 길이 s 로 적분한다.

    상태는 (x, y, px, py, L). p = n dr/ds 는 광선 벡터, L 은 누적 광로길이다.
      dr/ds = p / n,   dp/ds = grad n,   dL/ds = n
    stop(x, y) 가 0 을 지나면 멈춘다 (예: 지면에 닿음).
    """
    x0, y0 = r0
    n0 = medium.n(x0, y0)
    state0 = [x0, y0, n0 * np.cos(angle), n0 * np.sin(angle), 0.0]

    def rhs(s, u):
        x, y, px, py, _ = u
        n = medium.n(x, y)
        gx, gy = medium.grad(x, y)
        return [px / n, py / n, gx, gy, n]

    events = None
    if stop is not None:
        def ev(s, u):
            return stop(u[0], u[1])
        ev.terminal = True
        events = [ev]

    sol = solve_ivp(rhs, (0.0, s_max), state0, method="DOP853", rtol=rtol, atol=atol,
                    events=events, dense_output=True, max_step=max_step)
    return sol


def sample(sol, n=2000):
    """적분 결과를 호 길이에 대해 고르게 뽑는다. 반환: s, x, y, px, py, L."""
    s = np.linspace(sol.t[0], sol.t[-1], n)
    x, y, px, py, L = sol.sol(s)
    return s, x, y, px, py, L


# ── 페르마 원리: 경로를 잘게 나눠 광로길이를 직접 최소화 ─────
def opl_polyline(medium, xs, ys):
    """꺾은선 경로의 광로길이. 각 선분의 굴절률은 중점 값(중점 적분)."""
    dx, dy = np.diff(xs), np.diff(ys)
    xm, ym = 0.5 * (xs[1:] + xs[:-1]), 0.5 * (ys[1:] + ys[:-1])
    return np.sum(medium.n(xm, ym) * np.hypot(dx, dy))


def opl_hessian(medium, xs, ys, h=1e-6):
    """내부 점들의 y 좌표에 대한 광로길이의 헤시안 (중심차분)."""
    m = len(ys) - 2
    H = np.zeros((m, m))
    base = ys.copy()

    def f(v):
        y = base.copy()
        y[1:-1] = v
        return opl_polyline(medium, xs, y)

    v0 = base[1:-1].copy()
    f0 = f(v0)
    for i in range(m):
        for j in range(i, m):
            if i == j:
                e = np.zeros(m); e[i] = h
                H[i, i] = (f(v0 + e) - 2 * f0 + f(v0 - e)) / h**2
            else:
                ei = np.zeros(m); ei[i] = h
                ej = np.zeros(m); ej[j] = h
                H[i, j] = H[j, i] = (f(v0 + ei + ej) - f(v0 + ei - ej)
                                     - f(v0 - ei + ej) + f(v0 - ei - ej)) / (4 * h**2)
    return H


# ── 신기루: 뜨거운 도로 위 공기 ────────────────────────────
# 기준값: 15 °C, 1 기압 건조공기의 가시광 굴절률 n - 1 ≈ 2.77e-4 (자릿수 수준)
AIR_N1_REF = 2.77e-4
AIR_T_REF = 288.15


def air_temperature(y, t_far=303.15, t_road=323.15, height=0.10):
    """지면 위 높이 y [m] 의 기온 [K]. 지면 근처만 뜨거운 지수형 경계층."""
    return t_far + (t_road - t_far) * np.exp(-np.asarray(y) / height)


def road_air(t_far=303.15, t_road=323.15, height=0.10):
    """글래드스톤-데일 관계 (n-1) ∝ 밀도 ∝ 1/T (등압) 로 만든 층상 매질."""
    def n_of_y(y):
        return 1.0 + AIR_N1_REF * AIR_T_REF / air_temperature(y, t_far, t_road, height)

    def dn_dy(y):
        T = air_temperature(y, t_far, t_road, height)
        dT = -(t_road - t_far) / height * np.exp(-np.asarray(y) / height)
        return -AIR_N1_REF * AIR_T_REF * dT / T**2

    return stratified(n_of_y, dn_dy)


# ── 2차원 벡터 굴절·반사 ──────────────────────────────────
def reflect(d, nrm):
    """단위 방향 d 를 단위 법선 nrm 인 면에서 반사."""
    return d - 2.0 * np.sum(d * nrm, axis=-1, keepdims=True) * nrm


def refract(d, nrm, n1, n2):
    """벡터 스넬 법칙. nrm 은 입사 쪽을 향하는 단위 법선."""
    mu = n1 / n2
    cos_i = -np.sum(d * nrm, axis=-1, keepdims=True)
    sin2_t = mu**2 * (1.0 - cos_i**2)
    cos_t = np.sqrt(1.0 - sin2_t)
    return mu * d + (mu * cos_i - cos_t) * nrm


def hit_sphere(p, d, vertex_z, radius):
    """(z, y) 평면에서 광선 p + t d 와 꼭짓점 vertex_z, 반경 radius 인 원의 교점.

    광축(y=0) 근처의 교점(꼭짓점 쪽)을 고른다. 성분 순서는 (z, y).
    """
    c = np.array([vertex_z + radius, 0.0])
    oc = p - c
    b = np.sum(oc * d, axis=-1)
    q = np.sum(oc * oc, axis=-1) - radius**2
    disc = np.sqrt(b**2 - q)
    t1, t2 = -b - disc, -b + disc
    # 꼭짓점에 가까운 쪽: 반경이 양수면 앞쪽 교점, 음수면 뒤쪽 교점
    t = np.where(radius > 0, t1, t2)
    hit = p + t[..., None] * d
    nrm = (hit - c) / abs(radius)                    # 중심에서 바깥으로
    nrm = np.where(nrm[..., :1] > 0, -nrm, nrm)      # 입사 쪽(-z)을 향하게
    return hit, nrm


def sellmeier(wl_um, B, C):
    l2 = wl_um**2
    return np.sqrt(1.0 + sum(b * l2 / (l2 - c) for b, c in zip(B, C)))


N_BK7_B = (1.03961212, 0.231792344, 1.01046945)
N_BK7_C = (0.00600069867, 0.0200179144, 103.560653)


@dataclass
class Singlet:
    """양볼록 싱글렛. opticore 가 Optiland 와 대조해 둔 처방 그대로다."""

    r1: float = 1000.0
    r2: float = -1000.0
    thickness: float = 100.0
    wavelength_um: float = 0.750
    n_air: float = 1.0  # 유리는 카탈로그(공기 기준 상대) 굴절률이므로 공기는 1

    @property
    def n_glass(self):
        return sellmeier(self.wavelength_um, N_BK7_B, N_BK7_C)

    def trace(self, heights, z_start=-300.0):
        """광축에 평행한 광선을 높이 heights [mm] 에서 쏜다. 출사 광선 (점, 방향) 을 돌려준다."""
        h = np.asarray(heights, dtype=float)
        p = np.stack([np.full_like(h, z_start), h], axis=-1)
        d = np.tile([1.0, 0.0], (len(h), 1))
        p, nrm = hit_sphere(p, d, 0.0, self.r1)
        d = refract(d, nrm, self.n_air, self.n_glass)
        p, nrm = hit_sphere(p, d, self.thickness, self.r2)
        d = refract(d, nrm, self.n_glass, self.n_air)
        return p, d

    def axis_crossing(self, heights):
        """출사 광선이 광축(y=0)을 지나는 z."""
        p, d = self.trace(heights)
        return p[:, 0] - p[:, 1] * d[:, 0] / d[:, 1]


def envelope_of_lines(z0, y0, slope, u):
    """직선족 y = y0(u) + slope(u) (z - z0(u)) 의 포락선.

    y = a(u) + m(u) z 꼴로 바꾼 뒤 a'(u) + m'(u) z = 0 을 푼다 (미분은 수치 기울기).
    """
    a = y0 - slope * z0
    da, dm = np.gradient(a, u), np.gradient(slope, u)
    zc = -da / dm
    return zc, a + slope * zc


def ring_reflection(radius, heights):
    """+x 방향 평행광이 반경 radius 인 원형 거울 안쪽(오른쪽 벽)에서 한 번 반사된다.

    반사점과 반사 방향을 돌려준다. 성분 순서는 (x, y).
    """
    h = np.asarray(heights, dtype=float)
    x = np.sqrt(radius**2 - h**2)
    hit = np.stack([x, h], axis=-1)
    nrm = -hit / radius                     # 안쪽을 향하는 법선
    d = np.tile([1.0, 0.0], (len(h), 1))
    return hit, reflect(d, nrm)


def least_confusion(lens, h_max=100.0, n=20001, bounds=(1040.0, 1062.0)):
    """광선 다발의 번짐 지름이 가장 작은 z (최소 착란원) 와 그 지름.

    동공을 고르게 채운 광선을 추적하고, 단면에서 광선 높이의 최대-최소를 직접 최소화한다.
    """
    from scipy.optimize import minimize_scalar

    h = np.linspace(-h_max, h_max, n)
    p, d = lens.trace(h)
    slope = d[:, 1] / d[:, 0]

    def blur(z):
        y = p[:, 1] + (z - p[:, 0]) * slope
        return y.max() - y.min()

    r = minimize_scalar(blur, bounds=bounds, method="bounded", options={"xatol": 1e-9})
    return r.x, r.fun

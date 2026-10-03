"""글에 나오는 광학계들 (단위 mm, 파장 µm). generate_figures.py 와 verify_aberr.py 가 함께 쓴다."""
import numpy as np

from aberr import AIR, N_BK7, N_F2, N_SK16, Surface, glass, setup

WL = 0.5876                       # 헬륨 d 선 근처 (µm)
BK7, F2, SK16 = glass(N_BK7), glass(N_F2), glass(N_SK16)
N_D = BK7(WL)
F = 100.0                         # 단렌즈·색지움 렌즈의 초점거리 목표
T_LENS = 4.0


def plano_convex(D=10.0, stop_z=0.0, flip=False, semi=30.0):
    """초점거리 약 100 mm 의 N-BK7 평볼록 렌즈. flip=False 면 볼록면이 물체 쪽.
    조리개(반지름 D/2)는 stop_z 에 둔다 (0 이면 렌즈 첫 면 바로 앞)."""
    R = (N_D - 1.0) * F
    lens = ([Surface(0.0, R, BK7, semi), Surface(T_LENS, np.inf, AIR, semi)] if not flip else
            [Surface(0.0, np.inf, BK7, semi), Surface(T_LENS, -R, AIR, semi)])
    stop = Surface(stop_z, np.inf, AIR, D / 2, True)
    return [stop] + lens if stop_z <= 0.0 else lens + [stop]


def center_stop_z():
    """평면이 앞인 평볼록 렌즈에서, 뒤 곡면의 곡률 중심을 평면 너머로 본 겉보기 자리 (근축).
    여기에 조리개를 두면 주광선이 곡면에 수직으로 들어간다: 곡률 중심은 꼭짓점 앞 R − t 이고,
    유리 속에서 그 점을 향하던 광선은 공기 중에서 (R − t)/n 을 향한다."""
    R = (N_D - 1.0) * F
    return -(R - T_LENS) / N_D


def bent_lens(q, D=20.0, t=T_LENS):
    """모양 인자 q = (R2 + R1)/(R2 − R1) 로 휜, 얇은 렌즈 식으로 초점거리 100 mm 인 N-BK7 렌즈.
    q = +1: 볼록면이 물체 쪽인 평볼록, q = 0: 양볼록 대칭, q = −1: 평면이 물체 쪽."""
    P = 1.0 / ((N_D - 1.0) * F)
    c1, c2 = P * (q + 1) / 2, P * (q - 1) / 2
    R1 = np.inf if abs(c1) < 1e-15 else 1.0 / c1
    R2 = np.inf if abs(c2) < 1e-15 else 1.0 / c2
    return [Surface(0.0, np.inf, AIR, D / 2, True), Surface(0.0, R1, BK7, 40), Surface(t, R2, AIR, 40)]


def equiconvex(stop_z, stop_r=3.0, t=5.0):
    """초점거리 약 100 mm 의 N-BK7 양볼록 렌즈와 그 앞(음수) 또는 뒤(양수)의 조리개. stop_z=0 이면 렌즈 한가운데."""
    R = 2 * (N_D - 1.0) * F
    lens = [Surface(-t / 2, R, BK7, 60), Surface(t / 2, -R, AIR, 60)]
    if stop_z == 0.0:
        # 렌즈 한가운데의 조리개: 유리 속 평면
        return [lens[0], Surface(0.0, np.inf, BK7, stop_r, True), lens[1]]
    stop = Surface(stop_z, np.inf, AIR, stop_r, True)
    return [stop] + lens if stop_z < 0 else lens + [stop]


def achromat(D=10.0, exact=True):
    """N-BK7 양볼록 + N-F2 붙임 색지움 렌즈.
    1) 얇은 렌즈 조건 f1·V1 + f2·V2 = 0 으로 굴절력을 나누고 BK7 을 대칭으로 한다.
    2) 두께(5 mm, 2 mm)를 주면 F·C 선 초점이 조금 어긋나므로, exact=True 면 마지막 면 R3 만 고쳐
       두 선의 근축 뒤초점을 다시 맞춘다. 반환: 면 목록, (f1, f2, R1, R2, R3 얇은 렌즈값, R3 최종)."""
    from aberr import LINE_C, LINE_F, abbe_number, LINE_D, paraxial_trace
    V1, V2 = abbe_number(BK7), abbe_number(F2)
    P1, P2 = V1 / (F * (V1 - V2)), -V2 / (F * (V1 - V2))
    R1 = 2 * (BK7(LINE_D) - 1.0) / P1
    R2 = -R1
    R3 = 1.0 / (1.0 / R2 - P2 / (F2(LINE_D) - 1.0))

    def build(r3):
        return [Surface(0.0, np.inf, AIR, D / 2, True), Surface(0.0, R1, BK7, 15),
                Surface(5.0, R2, F2, 15), Surface(7.0, r3, AIR, 15)]

    def bf(lens, wl):
        m = paraxial_trace(lens, 1.0, 0.0, 0.0, wl)
        return lens[-1].z - m[-1, 0] / m[-1, 2]

    r3 = R3
    if exact:
        lo, hi = 1.0 / R3 - 2e-3, 1.0 / R3 + 2e-3          # 곡률로 이분법
        g = lambda c: bf(build(1.0 / c), LINE_F) - bf(build(1.0 / c), LINE_C)
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if np.sign(g(mid)) == np.sign(g(lo)):
                lo = mid
            else:
                hi = mid
        r3 = 1.0 / (0.5 * (lo + hi))
    return build(r3), (1 / P1, 1 / P2, R1, R2, R3, r3)


def singlet_for_color(D=10.0):
    """색지움 렌즈와 비교할 N-BK7 양볼록 단렌즈 (d 선 초점거리 약 100 mm)."""
    from aberr import LINE_D
    R = 2 * (BK7(LINE_D) - 1.0) * F
    return [Surface(0.0, np.inf, AIR, D / 2, True), Surface(0.0, R, BK7, 15), Surface(5.0, -R, AIR, 15)]


TRIPLET_RX = [  # (반지름, 뒤 두께, 유리) — 널리 쓰이는 쿡 삼중렌즈 예제 (f ≈ 50 mm, F/5, 반화각 20°)
    (22.01359, 3.25896, "SK16"), (-435.76044, 6.00755, None), (-22.21328, 0.99997, "F2"),
    (20.29192, 4.75041, None), (79.68360, 2.95208, "SK16"), (-18.39533, None, None)]
TRIPLET_STOP = 3   # 넷째 면(F2 렌즈 뒷면)에 조리개


def cooke_triplet(epd=10.0, wl=WL):
    media = {"SK16": SK16, "F2": F2, None: AIR}
    surfs, z = [], 0.0
    for i, (R, t, g) in enumerate(TRIPLET_RX):
        surfs.append(Surface(z, R, media[g], 9.0, i == TRIPLET_STOP))
        z += t or 0.0
    surfs[TRIPLET_STOP].semi = 1.0
    st = setup(surfs, wl, 0.0)
    surfs[TRIPLET_STOP].semi = (epd / 2) / st.r_ep
    return surfs

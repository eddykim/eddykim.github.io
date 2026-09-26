"""간섭 무늬 신호 — 이미지처리 6편.

간섭계든 구조광 투사든, 검출기가 담아내는 것은 언제나 같은 꼴이다.

    I(x) = a(x) + b(x) cos(phi(x))

a 는 배경(background), b 는 변조 진폭(modulation), phi 가 재고 싶은 위상이다.
문제는 세 가지가 한 장의 영상에 겹쳐 있고, 검출기는 그 합만 내놓는다는 것이다.
미지수 셋에 방정식 하나이므로 그대로는 풀리지 않는다.

이 시리즈의 앞 다섯 편이 전부 밝기를 읽는 이야기였다면, 여기서는 밝기에서
위상을 꺼내는 방법을 다룬다. 방법은 크게 둘로 갈린다. 한 장으로 푸는 쪽은
공간적인 가정(위상이 반송파를 타고 천천히 변한다)을 얹고, 여러 장으로 푸는
쪽은 위상을 아는 만큼 옮겨가며 찍어 미지수를 줄인다.
"""
import numpy as np


def background(x):
    """배경. 조명 불균일이라 위상보다 훨씬 천천히 변한다."""
    s = (x - x.min()) / (x.max() - x.min())
    return 0.55 + 0.30 * np.exp(-0.5 * ((s - 0.45) / 0.35) ** 2)


def modulation(x):
    """변조 진폭. 가장자리에서 대비가 떨어지는 실제 상황을 흉내 낸다."""
    s = (x - x.min()) / (x.max() - x.min())
    return 0.18 + 0.22 * np.sin(np.pi * s) ** 0.7


def true_phase(x, carrier, shape=1.0):
    """참 위상. 반송파에 시편 형상이 더해진 꼴이다.

    carrier 는 화소당 위상 증가량이고, 이 값이 언래핑 가능 여부를 정한다.
    shape 가 만드는 완만한 언덕이 재고 싶은 형상에 해당한다.
    """
    s = (x - x.min()) / (x.max() - x.min())
    return carrier * np.arange(len(x)) + shape * (
        6.0 * np.exp(-0.5 * ((s - 0.35) / 0.12) ** 2)
        - 4.0 * np.exp(-0.5 * ((s - 0.72) / 0.16) ** 2))


def fringe(x, carrier=0.45, shape=1.0, noise=0.0, rng=None, phase_shift=0.0):
    """간섭 무늬 한 줄. phase_shift 를 주면 위상천이 한 장이 된다."""
    a, b = background(x), modulation(x)
    phi = true_phase(x, carrier, shape)
    img = a + b * np.cos(phi + phase_shift)
    if noise > 0.0:
        img = img + noise * (rng or np.random.default_rng(0)).standard_normal(len(x))
    return img, a, b, phi


def fringe_2d(shape_xy=(200, 200), carrier=0.35, amp=14.0, noise=0.0, seed=0,
              phase_shift=0.0):
    """2차원 간섭 무늬. 언래핑의 경로 의존성을 보이려면 2차원이 필요하다.

    amp 는 형상 봉우리의 높이다. 이 값을 키우면 봉우리 옆면의 위상 기울기가
    화소당 pi 를 넘어서고, 그때부터 잔여점이 생긴다.
    """
    ny, nx = shape_xy
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    sx, sy = xx / nx, yy / ny
    # amp 가 커지면 봉우리 기울기가 화소당 pi 를 넘어 언샘플링이 일어난다.
    phi = carrier * xx + amp * (
        np.exp(-0.5 * (((sx - 0.35) / 0.10) ** 2 + ((sy - 0.40) / 0.10) ** 2))
        - 0.7 * np.exp(-0.5 * (((sx - 0.70) / 0.09) ** 2 + ((sy - 0.68) / 0.09) ** 2)))
    a = 0.55 + 0.25 * np.exp(-0.5 * (((sx - 0.5) / 0.6) ** 2 + ((sy - 0.5) / 0.6) ** 2))
    b = 0.20 + 0.18 * np.sin(np.pi * sx) * np.sin(np.pi * sy)
    img = a + b * np.cos(phi + phase_shift)
    if noise > 0.0:
        img = img + noise * np.random.default_rng(seed).standard_normal(img.shape)
    return img, a, b, phi

"""합성 분광 간섭 신호 — 이미지처리 3편.

학위논문 3.2~3.3절이 다루는 채널 분광 엘립소메트리의 신호를 재현한다. 다중차
위상 지연자를 지난 빛은 파수에 선형인 위상 지연을 얻으므로, 검출된 스펙트럼은
느리게 변하는 기저 신호(baseline signal, BLS) 위에 하나의 조화함수가 얹힌 꼴이 된다.

    I(sigma) = I_BLS(sigma) + A(sigma) * cos(2*pi*dnd*sigma - B(sigma))

기저 신호는 광원 스펙트럼과 시편 반사도를, 변조 신호(modulated signal, MDS)는
편광 정보를 담는다. 두 성분을 정확히 갈라내는 것이 해석의 출발점이고, 3편은
그 분리를 주파수 영역 필터로 할 때와 극값 포락선으로 할 때를 비교한다.

파수 단위는 1/um 다. 가시광 400~770 nm 가 1.30~2.50 1/um 에 해당한다.
"""
import numpy as np

# 다중차 위상 지연자의 위상 지연량 dn*d [um]. 대역 안에 줄무늬가 약 25개 생긴다.
RETARDER_DND = 20.8


def baseline(sigma):
    """기저 신호. 광원 스펙트럼과 박막 반사도가 섞인 완만한 곡선이다.

    양 끝 값이 서로 다르다는 점이 중요하다. 실제 스펙트럼도 400 nm 와 770 nm 에서
    반사도가 같을 이유가 없고, 이 불일치가 푸리에 기반 필터에서 문제를 일으킨다.
    """
    s = (sigma - sigma.min()) / (sigma.max() - sigma.min())     # 0~1 로 정규화
    return (0.62 + 0.30 * s                                      # 완만한 기울기
            - 0.22 * np.exp(-0.5 * ((s - 0.34) / 0.13) ** 2)     # 흡수 골
            + 0.09 * np.exp(-0.5 * ((s - 0.78) / 0.10) ** 2))    # 간섭 봉우리


def amplitude(sigma):
    """변조 진폭. 편광 소멸과 대역 가장자리의 감도 저하로 양 끝에서 줄어든다."""
    s = (sigma - sigma.min()) / (sigma.max() - sigma.min())
    return 0.26 * (0.55 + 0.45 * np.sin(np.pi * s) ** 0.6)


def phase_offset(sigma):
    """시편이 실어 보내는 위상 B. 파수에 따라 천천히 변한다."""
    s = (sigma - sigma.min()) / (sigma.max() - sigma.min())
    return 0.9 + 1.7 * s - 1.2 * s ** 2


def make_signal(n=1400, sigma_min=1.30, sigma_max=2.50, noise=0.0, seed=0):
    """분광 간섭 신호 한 벌을 만든다. 참 기저 신호도 함께 돌려준다."""
    sigma = np.linspace(sigma_min, sigma_max, n)
    bls = baseline(sigma)
    amp = amplitude(sigma)
    mds = amp * np.cos(2.0 * np.pi * RETARDER_DND * sigma - phase_offset(sigma))
    total = bls + mds
    if noise > 0.0:
        total = total + noise * np.random.default_rng(seed).standard_normal(n)
    return sigma, total, bls, amp


def fringe_count(sigma):
    """대역 안의 줄무늬 개수. 극값 사이 간격이 몇 화소인지 가늠할 때 쓴다."""
    return RETARDER_DND * (sigma.max() - sigma.min())

---
title: "XPS 기초 2편 — 왜 표면 10 nm만 보이는가"
categories: [표면분석, 깊이분석]
date: 2026-10-12 20:00:00 +0900
page_id: xps-surface-sensitivity-imfp
tags: [xps, imfp, surface-analysis, information-depth, arxps]
description: "X선은 시료 속으로 마이크로미터를 파고드는데 XPS는 왜 표면 10 nm 기법인가. 비탄성 평균자유행로와 정보깊이를 정리한다."
math: true
---

[1편](/posts/xps-photoemission-binding-energy/)은 에너지 보존식 $E_K = h\nu - E_B - \phi_{sp}$에서 결합에너지가 나오고, 코어 준위의 결합에너지가 원소 지문이 된다는 데까지 왔다. 그런데 XPS를 왜 "표면" 분석 기법이라고 부르는지는 아직 설명하지 않았다. 1편의 식에는 깊이라는 변수가 어디에도 없다.

이상한 점은 X선 쪽에 있다. Al K$\alpha$ X선(1486.6 eV)은 실리콘 속으로 수 마이크로미터를 들어간다. 광전효과는 X선이 닿는 곳이면 어디서든 일어나므로, 광전자도 그 깊이 전체에서 만들어진다. 그런데도 XPS 스펙트럼의 피크는 표면 수 nm에서 온 신호로 채워진다. 표면민감성을 만드는 것은 X선이 아니다. 그렇다면 무엇인가.

이번 편은 이 질문에 답하면서 "XPS는 표면 10 nm를 본다"는 익숙한 문장을 분해한다. 결론부터 말하면 이 10 nm는 물리 상수가 아니다. 지수 감쇠 모델, 95%라는 관례적 기준, 예측식으로 계산한 평균자유행로, 수직 검출이라는 기하 조건을 곱해 얻은 숫자다. 네 가지 중 어느 하나가 바뀌면 숫자도 바뀐다.

## 1. 표면민감성을 만드는 것은 X선이 아니라 전자다

X선은 깊이 들어가는데 왜 신호는 표면에서만 오는가. 먼저 두 입자가 고체 속에서 얼마나 멀리 가는지 비교해 보자. X선 데이터베이스(CXRO)에서 Al K$\alpha$ 광자의 감쇠길이를 찾으면 실리콘에서 7.9 µm, SiO₂에서 4.1 µm, 금에서 0.23 µm다. 같은 에너지 영역의 광전자가 에너지를 잃지 않고 가는 거리는 3절에서 볼 비탄성 평균자유행로로 재는데, 실리콘의 Si 2p 광전자는 약 3.1 nm, 금의 Au 4f 광전자는 약 1.6 nm다. 비를 내면 실리콘에서 약 2,500배, 금에서 약 150배다. 두세 자릿수 차이다.

이 차이가 뜻하는 바는 이렇다. X선의 감쇠길이에 비하면 표면 수 nm는 무시할 만큼 얇다. 그래서 신호가 나올 수 있는 깊이 안에서는 광전자가 깊이에 상관없이 거의 균일하게 만들어진다. 문제는 밖으로 나오는 길이다. 깊은 곳에서 출발한 광전자는 고체 안의 전자들과 부딪치며 비탄성 산란(inelastic scattering)을 겪는다. 플라즈몬을 들뜨게 하거나 가전자를 들뜨게 하면서 운동에너지를 조금씩 잃는다. 에너지를 잃은 전자도 표면을 빠져나와 분석기에 도달할 수는 있다. 그러나 운동에너지가 달라졌으므로 원래 피크 자리에 나타나지 못한다. 1편 그림3에서 피크마다 높은 결합에너지 쪽에 붙어 있던 계단형 배경이 바로 이 전자들이다.

따라서 피크 면적에 기여하는 전자는 **한 번도 비탄성 산란을 겪지 않고 탈출한 전자**뿐이다. 표면민감성의 정체는 이것이다. 깊은 곳의 정보가 사라지는 것이 아니라, 피크에서 배경으로 자리를 옮긴다. 이 관점은 4편에서 배경을 어떤 모양으로 빼야 하는지 따질 때 다시 쓰인다.

## 2. IMFP의 정의와 감쇠식

"무손실로 탈출할 확률"을 어떻게 정량화하는가. 비탄성 산란이 서로 독립인 사건이고 전자가 지나는 거리마다 같은 확률로 일어난다고 가정하면, 거리 $s$를 가는 동안 한 번도 산란하지 않을 확률은 지수함수가 된다.

$$ P(s) = e^{-s/\lambda} $$

여기서 $\lambda$가 비탄성 평균자유행로(inelastic mean free path, IMFP)다. 국제표준 ISO 18115는 IMFP를 "주어진 에너지의 전자가 연속한 두 비탄성 충돌 사이에 이동하는 평균 거리"로 정의한다. 식으로 보면 무손실 확률이 $1/e$로 떨어지는 거리다.

깊이 $z$에서 만들어진 광전자가 표면 법선에서 $\theta$만큼 기운 방향으로 분석기를 향한다면, 고체 안에서 지나는 거리는 $z/\cos\theta$다. 따라서 이 전자가 무손실로 탈출할 확률은 다음과 같다.

$$ P(z,\theta) = \exp\!\left(-\frac{z}{\lambda\cos\theta}\right) $$

이 글에서 $\theta$는 항상 **표면 법선에서 잰 각**이다. 5절에서 이 약속이 왜 중요한지 다룬다.

식의 모양은 [엘립소메트리 배경이론 1편](/posts/ellipsometry-electromagnetic-fresnel/)에서 본 Beer–Lambert 법칙과 같다. 그러나 감쇠하는 주체가 다르다. Beer–Lambert는 광자가 흡수되어 사라지는 과정이고, 여기서는 전자가 사라지지 않고 에너지만 잃는 과정이다. 광자는 흡수되면 신호에서 빠지지만, 전자는 배경으로 옮겨 가 스펙트럼에 계속 남는다.

실무에서는 IMFP와 비슷한 이름의 양이 셋 더 쓰이고, 이 넷을 섞어 쓰는 것이 흔한 오차의 원인이다. Powell(2020)의 정리를 따라 구분하면 다음과 같다.

| 양 | 정의 | 무엇에 의존하는가 |
|---|---|---|
| IMFP $\lambda$ | 비탄성 충돌 사이의 평균 거리 | 물질, 전자 에너지 |
| 실효 감쇠길이(effective attenuation length, EAL) | 탄성 산란을 무시한 식에 IMFP 대신 넣으면 탄성 산란 효과까지 보정되도록 정한 값 | 위 둘 + 측정 기하, 용도 |
| 평균 탈출 깊이(mean escape depth, MED) | 검출된 전자가 출발한 깊이의 평균 | 위와 같음 |
| 정보깊이(information depth, ID) | 신호의 지정 비율(95% 또는 99%)이 나오는 깊이 | 위와 같음 |

IMFP만 물질과 에너지로 정해지는 물성값이고, 나머지 셋은 측정 조건까지 끌어들인다. 차이를 만드는 것은 탄성 산란(elastic scattering)이다. 전자가 원자핵 근처에서 에너지를 잃지 않고 방향만 꺾이면 실제 경로가 직선보다 길어지고, 같은 깊이에서 출발해도 비탄성 산란을 겪을 기회가 늘어난다. 그래서 EAL은 IMFP보다 짧다. Powell이 정리한 계산에 따르면 XPS에서 주로 쓰는 200 eV~1.5 keV 구간에서 EAL/IMFP 비는 알루미늄과 실리콘에서 0.77~0.92, 구리·은·인듐·금에서 0.68~0.85다. 탄성 산란은 EAL, MED, ID를 최대 40% 가까이 바꿀 수 있다. 이 글의 계산은 탄성 산란을 무시하고 IMFP로 하므로, 실제 깊이는 이보다 얕다고 읽어야 한다.

## 3. 유니버설 커브 — λ는 운동에너지에 어떻게 의존하는가

어느 피크가 더 표면을 보는가. $\lambda$는 물질에 따라 다르지만, 그보다 더 크게 움직이는 것은 전자의 운동에너지다. 1970년대 여러 물질에서 측정한 값을 모아 운동에너지를 가로축으로 그리면 대부분 비슷한 V자 모양으로 모인다. 이것을 유니버설 커브(universal curve)라고 부른다. Seah와 Dench(1979)는 원소의 경우 이 곡선을 단원자층 단위의 경험식 하나로 정리했다.

$$ \lambda_m = \frac{538}{E^2} + 0.41\,(aE)^{1/2} $$

$E$는 eV 단위의 운동에너지, $a$는 단원자층 두께(nm), $\lambda_m$은 단원자층 수다. 첫 항은 저에너지에서 가전자를 들뜨게 할 에너지가 모자라 산란이 줄어드는 효과, 둘째 항은 고에너지에서 전자가 빨라 상호작용 시간이 짧아지는 효과를 나타낸다. $a = 0.25$ nm로 두고 계산하면 최소값은 약 41 eV에서 0.41 nm다. 수십 eV의 전자가 가장 표면에 민감하다.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/fig1-imfp-universal-curve.png" alt="Si, SiO2, Au의 TPP-2M IMFP와 Seah–Dench 유니버설 커브, Al Kα 여기 대표 피크 위치" width="640">
_그림1. 운동에너지에 따른 IMFP. 실선은 TPP-2M 예측식(50 eV 이상에서 유효), 점선은 Seah–Dench 경험식이다. 원은 Al K$\alpha$로 여기했을 때 Si 2p(Si), O 1s(SiO₂), Au 4f(Au) 광전자의 위치다._

오늘날 실무에서 쓰는 값은 경험식이 아니라 계산값이다. Tanuma, Powell, Penn은 물질의 광학 상수에서 에너지 손실 함수를 구해 IMFP를 계산하고, 그 결과를 물질 파라미터 네 개(가전자 수, 밀도, 분자량, 밴드갭)만으로 재현하는 예측식을 만들었다. 이것이 TPP-2M이다. 그림1의 실선은 이 식으로 직접 계산한 것이다.

```python
def tpp2m(E, Nv, rho, M, Eg=0.0):
    """TPP-2M IMFP (nm). E: 운동에너지(eV), Nv: 가전자 수, rho: 밀도(g/cm^3)."""
    Ep = 28.816 * np.sqrt(Nv * rho / M)            # 자유전자 플라즈몬 에너지
    U = (Ep / 28.816) ** 2
    beta = -1.0 + 9.44 / np.sqrt(Ep**2 + Eg**2) + 0.69 * rho**0.1
    gamma = 0.191 * rho**-0.5
    C, D = 19.7 - 9.1 * U, 534 - 208 * U
    return E / (Ep**2 * (beta * np.log(gamma * E) - C / E + D / E**2))
```

계수는 Powell(2020)에 실린 식을 그대로 옮겼다. 전체 스크립트는 [`_code/xps-surface-sensitivity-imfp/generate_figures.py`](https://github.com/eddykim/eddykim.github.io/blob/main/_code/xps-surface-sensitivity-imfp/generate_figures.py)에 있다.

그림1에서 읽을 수 있는 수치 감각은 다음과 같다. 50~100 eV 근처에서 $\lambda$는 0.4~0.8 nm로, 곡선에서 가장 짧은 영역이다. Al K$\alpha$로 여기한 광전자가 주로 놓이는 수백 eV~1.4 keV 구간에서는 에너지가 올라갈수록 커져 1.4 keV에서 1.6 nm(금)~4 nm(SiO₂)에 이른다. 수십 nm에 닿으려면 운동에너지가 10 keV를 넘어야 하는데, 이는 고에너지 X선을 쓰는 HAXPES(hard X-ray photoelectron spectroscopy)의 영역이다. 고에너지 쪽은 대략 $E^{\,p}$에 비례해 커진다. Seah–Dench 식의 고에너지 극한은 $p = 0.5$이고, 여러 물질의 실측 감쇠길이에 맞춘 $p$는 0.54~0.81로 보고되어 있다.

여기서 이 시리즈의 축이 처음으로 구체적인 모습을 드러낸다. 두께나 조성을 계산할 때 넣는 $\lambda$는 대개 측정값이 아니라 TPP-2M이나 NIST 데이터베이스(SRD 71)에서 가져온 계산값이다. Powell에 따르면 광학 데이터로 계산한 IMFP의 불확도는 최대 10% 정도이고, TPP-2M은 그 계산값과 제곱평균으로 9% 정도 어긋나며, 탄성 피크 측정으로 얻은 실험값과 계산값은 제곱평균 12~15% 차이 난다. 두께 결과는 $\lambda$에 정비례하므로, $\lambda$의 불확도가 그대로 두께의 불확도가 된다.

한 가지 더 눈여겨볼 점이 있다. 같은 시료 안에서도 피크마다 운동에너지가 다르므로 **피크마다 보는 깊이가 다르다**. SiO₂에서 O 1s 광전자(운동에너지 약 950 eV)의 $\lambda$는 2.9 nm인데, Si 2p 광전자(약 1,380 eV)는 이보다 길다. 두 피크의 세기 비로 조성을 계산하면 서로 다른 깊이 범위를 평균한 두 값을 나누는 셈이 된다. 이 문제는 4편의 정량과 5편의 두께 역산에 그대로 걸린다.

## 4. 정보깊이 3λ가 95%인 이유

"정보깊이 10 nm"라는 숫자는 어디서 나왔는가. 수직 검출($\theta = 0$)을 가정하고 깊이별 기여를 적분해 보자. 광전자가 깊이와 상관없이 균일하게 만들어진다면, 깊이 $z$에서 나오는 신호는 $e^{-z/\lambda}$에 비례한다. 표면부터 깊이 $d$까지가 전체 신호에서 차지하는 몫은 다음과 같다.

$$ \frac{\int_0^{d} e^{-z/\lambda}\,dz}{\int_0^{\infty} e^{-z/\lambda}\,dz} = 1 - e^{-d/\lambda} $$

$d$에 $\lambda$의 정수배를 넣으면 익숙한 숫자가 나온다. $d = \lambda$에서 63.2%, $2\lambda$에서 86.5%, $3\lambda$에서 95.0%다.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/fig2-depth-contribution.png" alt="Si 2p 광전자의 깊이별 기여와 누적 기여 곡선, λ·2λ·3λ 지점 표시" width="720">
_그림2. 실리콘의 Si 2p 광전자($\lambda$ = 3.09 nm, TPP-2M)의 깊이별 신호 기여. (a) 1 nm 구간마다의 몫, (b) 누적 기여._

실리콘 속 Si 2p 광전자의 $\lambda$ = 3.09 nm를 넣으면 $3\lambda$ = 9.3 nm다. "XPS는 표면 10 nm를 본다"는 말은 여기서 나온다. 그러나 이 숫자에 들어간 가정을 하나씩 풀어 보면 사정이 달라진다.

첫째, 95%는 물리가 정한 값이 아니라 관례다. ISO 정의도 정보깊이를 "지정한 비율"의 신호가 나오는 깊이로 두고, 그 비율로 95%와 99%를 함께 예로 든다. 99%를 기준으로 잡으면 $4.6\lambda$, 약 14 nm가 된다. 둘째, $\lambda$는 피크와 물질마다 다르다. 같은 방식으로 계산하면 SiO₂의 O 1s는 $3\lambda$ = 8.8 nm, 금의 Au 4f는 4.7 nm다. 금 위의 박막이라면 "10 nm"의 절반밖에 보지 못한다. 셋째, 2절에서 본 대로 탄성 산란까지 넣으면 실제 정보깊이는 이보다 얕다. 넷째, 수직 검출을 가정했다. 이것은 5절에서 다룬다.

더 중요한 따름결과는 신호가 깊이를 지수 함수로 가중한 **가중 평균**이라는 점이다. 그림2(a)에서 1 nm 구간마다의 몫을 보면 최표면 1 nm가 28%, 그다음 1 nm가 20%, 그다음이 14%, 10%로 줄어든다. 첫 1 nm가 셋째 nm보다 두 배 가까이 무겁다. 금이라면 최표면 1 nm의 몫이 47%에 이른다. XPS가 보고하는 조성은 "표면 10 nm의 평균 조성"이 아니라 최표면에 무게가 쏠린 가중 평균이다. 깊이에 따라 조성이 변하는 시료라면 이 차이가 결과를 바꾼다.

## 5. 각도를 바꾸면 보는 깊이가 바뀐다 — 그리고 문헌의 각도 혼동

더 얕게 보고 싶다면 무엇을 바꿔야 하는가. X선 에너지를 바꾸기는 어렵지만, 시료를 기울이는 것은 쉽다. 2절의 식에서 $\lambda$ 자리에 $\lambda\cos\theta$가 들어 있었다. 검출 방향을 법선에서 멀리 기울이면 같은 깊이에서 출발한 전자가 더 긴 경로를 지나야 하므로, 무손실로 나오는 전자의 출발 깊이가 얕아진다. 평균 탈출 깊이도 정보깊이도 $\cos\theta$에 비례해 줄어든다.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/fig3-takeoff-angle.png" alt="검출 각도에 따른 탈출 경로 개념도와 각도별 누적 기여 곡선" width="720">
_그림3. (a) 같은 깊이 $z$에서 출발해도 $\theta$ = 60°로 검출하면 탈출 경로가 두 배가 된다. (b) 각도별 누적 기여. 점은 95% 깊이 $3\lambda\cos\theta$다._

그림3(b)는 실리콘의 Si 2p를 기준으로 각도를 바꿔 가며 누적 기여를 그린 것이다. 95% 깊이는 $\theta$ = 0°에서 9.3 nm, 45°에서 6.6 nm, 60°에서 4.6 nm, 75°에서 2.4 nm다. 최표면 1 nm의 몫도 28%에서 60°에서는 48%, 75°에서는 71%로 커진다. 수직 검출과 경사 검출의 스펙트럼을 비교하면 어떤 성분이 표면 쪽에 있는지 판별할 수 있다. 이 원리를 여러 각도로 확장한 것이 각도 분해 XPS(angle-resolved XPS, ARXPS)이며, 5편에서 본격적으로 다룬다. 다만 이 $\cos\theta$ 관계에는 범위가 있다. Powell은 방출각 60° 정도까지는 EAL 하나로 충분하지만 그보다 기울이면 EAL 자체가 각도에 따라 달라진다고 정리한다. 75°처럼 표면에 가깝게 누우면 탄성 산란과 표면 거칠기의 영향이 커지기 때문이다.

여기서 실무에서 자주 나는 실수 하나를 짚어야 한다. 검출 각도를 문헌마다 다른 기준으로 잰다. 이 글과 Powell(2020)은 표면 법선에서 잰 방출각(emission angle)을 쓰므로 유효 깊이가 $\lambda\cos\theta$로 적힌다. 반면 표면 평면에서 잰 출사각(take-off angle)을 쓰는 문헌에서는 같은 물리가 $\lambda\sin\theta$로 적힌다. 수직 검출이 한쪽에서는 0°이고 다른 쪽에서는 90°다.

대표적인 예가 산화막 두께 계산에 널리 쓰이는 Strohmeier 식이다.

$$ d = \lambda_{ox}\sin\theta\,\ln\!\left(\frac{N_m\lambda_m I_{ox}}{N_{ox}\lambda_{ox} I_m} + 1\right) $$

$N$은 금속 원자의 부피 밀도, $I$는 산화물과 금속 성분의 피크 면적이다. 이 식이 법선 기준의 $\cos\theta$형 식과 같은 물리를 나타내려면 $\theta$는 표면에서 잰 각이어야 한다. 수직 검출에서 $\sin 90° = 1$이 되어 $\lambda_{ox}$가 그대로 남기 때문이다. 그런데 이 식을 소개하는 인기 있는 실무 자료들을 찾아보면, $\theta$를 "take-off angle"이라고만 부르고 어느 면에서 재는지는 적지 않은 경우가 많다. 법선 기준으로 45°가 아닌 각도를 쓰는 장비에서 이 식에 그 각도를 그대로 넣으면 두께가 틀린다. 예를 들어 법선에서 30° 기운 분석기라면 넣어야 할 값은 $\sin 60°$ = 0.87인데, $\sin 30°$ = 0.5를 넣으면 두께를 42% 작게 얻는다. 식을 가져다 쓸 때는 각도의 기준면부터 확인해야 한다.

## 6. 표면민감성의 대가

표면만 보는 것은 언제 손해인가. 표면민감성은 장점이지만, 정보깊이 전체가 알고 싶은 대상으로 채워져 있을 때만 그렇다.

대기에 노출된 시료의 표면에는 거의 예외 없이 흡착 탄소(adventitious carbon)가 덮여 있다. 공기 중의 탄화수소와 산소 함유 유기물이 달라붙은 층으로, 두께는 보통 1~2 nm다. 이 두께는 IMFP와 같은 자릿수다. 탄화수소층을 폴리에틸렌 조성으로 근사해 TPP-2M으로 계산하면, 그 안에서 Si 2p 광전자의 $\lambda$는 약 4.7 nm다. 흡착 탄소층이 1 nm면 하부 신호의 81%, 2 nm면 65%만 남는다. 4절에서 본 대로 신호는 최표면에 무게가 쏠려 있으므로, 가장 무거운 첫 1~2 nm를 알고 싶지 않은 층이 차지하는 셈이다. 그렇다고 이 층을 이온 식각으로 걷어내면 그 과정에서 하부의 화학상태가 바뀔 수 있다. 이 문제는 5편의 주제다.

벌크 정보는 원리적으로 얻을 수 없다. 합금이나 도핑된 재료에서는 특정 원소가 표면으로 몰리는 표면 편석(surface segregation)이 흔하다. 이때 XPS 조성이 벌크 조성과 다르게 나오는 것은 측정 오차가 아니다. XPS는 정의상 표면 조성을 재는 기법이고, 그 결과를 벌크 조성으로 읽는 쪽이 해석을 잘못한 것이다.

XPS 장비가 초고진공(ultra-high vacuum, UHV)을 요구하는 이유도 표면민감성에서 나온다. 이유는 두 가지다. 하나는 재오염이다. 기체 분자가 표면에 부딪치는 빈도는 압력에 비례하며, $10^{-6}$ Torr에서는 흡착 확률이 1일 때 약 1초면 단분자층 하나를 덮을 만큼 부딪친다(이 노출량을 1 랭뮤어(langmuir)라 한다). 분석실 압력이 $10^{-10}$ mbar 수준이면 같은 노출에 수 시간이 걸려, 측정하는 동안 표면이 유지된다. 다른 하나는 광전자의 기체 산란이다. 광전자가 분석기까지 가는 긴 경로에서 기체 분자와 부딪쳐 에너지를 잃으면, 고체 안에서 잃은 것과 마찬가지로 피크에서 빠진다.

## 7. 광학 계측과 비교하면 무엇이 다른가

같은 박막을 광학으로 재면 무엇이 보이는가. [엘립소메트리 배경이론 1편](/posts/ellipsometry-electromagnetic-fresnel/)에서 빛의 침투 깊이를 $\delta = \lambda/4\pi k$로 적었다. 흡수가 약한 파장대에서 $k$가 작으면 $\delta$는 수백 nm에서 수 µm에 이른다. 엘립소메트리는 다층 박막 전체를 한 번에 통과하며 각 층의 두께와 광학 상수를 모델로 역산한다. XPS는 그 가운데 최상단 수 nm만 본다.

두 기법이 잘 주는 것도 다르다. XPS는 원소와 화학상태를 직접 구분하지만, 두께는 가장 간접적으로 얻는 양이다. 엘립소메트리는 두께와 광학 상수를 서브나노미터 감도로 주지만, 그 층이 어떤 화학 결합으로 이루어졌는지는 광학 상수 모델에 미리 넣어 줘야 한다. 그래서 두 기법은 경쟁 관계라기보다 서로를 보완한다. 얇은 산화막 하나를 두고 XPS로 화학상태와 대략의 두께를, 엘립소메트리로 정밀한 두께를 얻어 교차 확인하는 방식이 흔하다. 5편에서 XPS 두께 계측의 한계를 다룰 때 이 비교로 돌아온다.

## 정리 및 다음 편 예고

이번 편에서는 XPS의 표면민감성이 어디서 오는지를 따라갔다. X선은 고체 속으로 마이크로미터를 들어가지만, 피크에 기여하는 것은 비탄성 산란 없이 탈출한 광전자뿐이고, 그 확률이 $\exp(-z/\lambda\cos\theta)$로 떨어진다. IMFP $\lambda$는 운동에너지에 따라 유니버설 커브를 그리며, Al K$\alpha$ XPS 영역에서 1.6~4 nm다. 실무에서 쓰는 $\lambda$는 대개 TPP-2M 같은 예측식이 준 계산값이고, 10% 안팎의 불확도를 안고 있다. 정보깊이 "10 nm"는 95%라는 관례, 실리콘 계열의 $\lambda$, 수직 검출, 탄성 산란 무시라는 가정을 모두 곱한 숫자다. 신호는 깊이를 지수 함수로 가중한 평균이어서 최표면에 무게가 쏠리고, 검출 각도를 기울이면 이 쏠림이 더 강해진다. 각도를 인용할 때는 기준면부터 확인해야 한다.

다음 편은 스펙트럼 자체로 돌아간다. 1편의 서베이 스펙트럼에는 원소마다 피크가 하나씩 있는 것처럼 보였지만, 고분해능으로 보면 사정이 다르다. 스핀-궤도 분리로 갈라진 짝 피크, 셰이크업(shake-up) 위성, 다중항 분리, 플라즈몬 손실, 오제 계열이 함께 있다. 각각이 초기상태 효과인지 최종상태 효과인지 가려야 화학상태를 올바르게 읽을 수 있다.

## 참고자료

- C. J. Powell, "Practical guide for inelastic mean free paths, effective attenuation lengths, mean escape depths, and information depths in x-ray photoelectron spectroscopy," *J. Vac. Sci. Technol. A* 38, 023209 (2020). <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=929257>
- S. Tanuma, C. J. Powell, and D. R. Penn, "Calculations of electron inelastic mean free paths. V. Data for 14 organic compounds over the 50–2000 eV range," *Surf. Interface Anal.* 21, 165–176 (1994). <https://doi.org/10.1002/sia.740210302>
- M. P. Seah and W. A. Dench, "Quantitative electron spectroscopy of surfaces: A standard data base for electron inelastic mean free paths in solids," *Surf. Interface Anal.* 1, 2–11 (1979). <https://doi.org/10.1002/sia.740010103>
- B. R. Strohmeier, "An ESCA method for determining the oxide thickness on aluminum alloys," *Surf. Interface Anal.* 15, 51–56 (1990). <https://doi.org/10.1002/sia.740150109>
- NIST Standard Reference Database 71: Electron Inelastic-Mean-Free-Path Database. <https://www.nist.gov/srd/nist-standard-reference-database-71>
- NIST Standard Reference Database 82: Electron Effective-Attenuation-Length Database. <https://www.nist.gov/srd/nist-standard-reference-database-82>
- Center for X-Ray Optics (CXRO), X-Ray Attenuation Length. <https://henke.lbl.gov/optical_constants/atten2.html>
- Thermo Fisher Scientific, XPS Periodic Table — Carbon. <https://www.thermofisher.com/us/en/home/materials-science/learning-center/periodic-table/non-metal/carbon.html>
- D. N. G. Krishna and J. Philip, "Review on surface-characterization applications of X-ray photoelectron spectroscopy (XPS): Recent developments and challenges," *Appl. Surf. Sci. Adv.* 12, 100332 (2022). <https://doi.org/10.1016/j.apsadv.2022.100332>
- D. J. Morgan, "X-Ray Photoelectron Spectroscopy (XPS): An Introduction," Cardiff Catalysis Institute. <https://sites.cardiff.ac.uk/xpsaccess/files/2014/07/AccessXPS_Primer_Paper.pdf>

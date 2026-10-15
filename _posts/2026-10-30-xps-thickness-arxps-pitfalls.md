---
title: "XPS 기초 5편 — XPS로 잰 산화막 두께는 얼마나 믿을 수 있는가"
categories: [표면분석, 깊이분석]
date: 2026-10-30 20:00:00 +0900
page_id: xps-thickness-arxps-pitfalls
tags: [xps, arxps, film-thickness, depth-profiling, beam-damage, metrology]
description: "기판 신호가 줄어든 만큼으로 산화막 두께를 잴 수 있다. 그 값이 감쇠길이, R₀, 각도 규약, 스퍼터링과 X선 손상에 얼마나 흔들리는지 본다."
math: true
---

[4편](/posts/xps-quantification-peak-fitting/)까지는 시편이 정보깊이 안에서 깊이 방향으로 균일하다고 가정했다. 조성비 식에 들어간 $\lambda$가 그 가정을 품고 있었다. 그런데 실제 시편은 거의 항상 층 구조다. 실리콘 웨이퍼 위의 자연산화막, 모든 표면을 덮은 흡착 탄소층, 일부러 올린 증착막이 그렇다. 4편의 합성 스펙트럼도 얇은 SiO₂가 덮인 실리콘이었다.

층 구조는 방해 요인이기만 한 것이 아니다. [2편](/posts/xps-surface-sensitivity-imfp/)에서 본 지수 감쇠는 깊이 정보를 담고 있다. 산화막이 두꺼울수록 그 아래 기판에서 나오는 신호는 더 많이 줄어든다. 그러니 신호가 얼마나 줄었는지를 거꾸로 읽으면 산화막의 두께를 알 수 있다. 이번 편은 이 원리로 두께를 재는 방법과, 그 방법이 어디서 틀어지는지를 다룬다.

결론부터 말하면 두께는 XPS가 내놓는 양 가운데 가장 간접적이다. 강도비 하나에서 모델을 거쳐 나오며, 층이 균일하다는 가정과 $\lambda$ 추정값에 전적으로 기댄다. 그리고 이 편의 후반부는 한 걸음 더 나간다. 깊이를 보려고 시편을 깎거나 X선을 오래 쬐면, 측정이 측정 대상 자체를 바꾼다.

## 1. 기판 신호가 줄어든 정도로 두께를 알 수 있는가

산화막 아래 기판의 신호가 아직 보인다면, 얼마나 남았는지로 두께를 알 수 있는가. 안개 너머의 등불을 떠올리면 된다. 안개가 두꺼울수록 불빛이 흐려지므로, 불빛이 얼마나 줄었는지에서 안개의 두께를 거꾸로 짐작할 수 있다. 가장 단순한 기하에서 시작하자. 두께 $d$의 균일한 오버레이어(overlayer), 예를 들어 SiO₂가 무한히 두꺼운 기판(Si) 위에 있고, 분석기는 표면 법선에서 $\theta$만큼 기운 방향에 있다. 2편과 같이 이 글에서도 $\theta$는 **법선 기준**으로 잰 방출각이다.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/fig1-overlayer-geometry.png" alt="균일 오버레이어 모델의 기하와 두께에 따른 강도비 곡선" width="720">
_그림1. (a) 균일 오버레이어 모델. 기판 신호는 산화막을 지나며 감쇠한다. (b) 두께에 따른 강도비 (합성 계산, $\lambda$ = 3 nm). 실선은 수직 방출($\theta$ = 0°), 검은 점선은 $\theta$ = 60°이고, 빨간 점선은 수직 방출에서 기판 신호가 5%로 떨어지는 $3\lambda$다. 60°에서는 같은 두께의 강도비가 커지는 대신 잴 수 있는 두께가 절반으로 준다._

오버레이어 안의 깊이 $z$에서 나온 광전자가 무손실로 빠져나올 확률은 $\exp(-z/\lambda\cos\theta)$다. 이것을 표면부터 $d$까지 적분하면 오버레이어의 신호가 된다.

$$ I_o = I_o^{\infty} \left[ 1 - \exp\!\left( -\frac{d}{\lambda_o \cos\theta} \right) \right] $$

$I_o^{\infty}$는 오버레이어가 무한히 두꺼울 때의 세기다. 기판의 신호는 두께 $d$의 오버레이어를 통째로 지나야 하므로 한 번 더 감쇠한다.

$$ I_m = I_m^{\infty} \exp\!\left( -\frac{d}{\lambda_m \cos\theta} \right) $$

즉, 기판 신호는 산화막을 $\lambda_m\cos\theta$만큼 지날 때마다 $1/e$, 약 37%로 준다. 수직 방출에서 감쇠길이가 3 nm라면 3 nm 막 아래에서는 37%, 9 nm 막 아래에서는 5%가 남는다.

두 식을 나누면 강도비 $R = I_o/I_m$이 나온다. SiO₂/Si의 Si 2p처럼 두 피크의 운동에너지가 거의 같으면 $\lambda_o \approx \lambda_m = \lambda$로 둘 수 있고, 이때 식이 깔끔하게 풀린다.

$$ R = R_0 \left[ \exp\!\left( \frac{d}{\lambda\cos\theta} \right) - 1 \right], \qquad R_0 = \frac{I_o^{\infty}}{I_m^{\infty}} $$

이를 $d$로 정리하면 강도비 하나에서 두께가 나온다.

$$ d = \lambda \cos\theta \, \ln\!\left( 1 + \frac{R}{R_0} \right) $$

즉, 강도비가 $R_0$의 몇 배인지가 두께를 정한다. 막이 아주 얇으면($d \ll \lambda$) 로그를 펼쳐 $R \approx R_0\, d/(\lambda\cos\theta)$가 되므로 강도비가 두께에 비례한다. $\lambda$ = 3 nm와 2절에서 볼 SiO₂/Si의 실측값 $R_0$ = 0.88을 쓰면, 0.5 nm 막의 강도비는 정확한 식으로 0.160, 이 직선 근사로 0.147이다. 반대로 두 피크의 면적이 같으면($R$ = 1) 수직 방출에서 두께는 약 2.3 nm다. 스펙트럼을 보고 바로 대 볼 수 있는 기준이다.

그림1(b)는 $\lambda$ = 3 nm에서 이 관계를 그린 것이다. 이 값은 2편의 TPP-2M 식으로 계산한 SiO₂ 속 Si 2p의 비탄성 평균자유행로(inelastic mean free path, IMFP) 3.90 nm에, 2편에서 인용한 Powell의 실효 감쇠길이(effective attenuation length, EAL)와 IMFP의 비 범위(Al·Si, 0.77~0.92)의 하한을 곱한 어림값이다. Si 2p의 운동에너지(약 1.4 keV)에서는 이 비가 0.9를 넘으므로 실제 감쇠길이는 이보다 15%쯤 길다. 이 편의 결론은 $\lambda$를 3.5 nm로 바꿔도 달라지지 않는다.

그림1(b)를 읽어 보자. 두께가 늘면 강도비는 지수적으로 커지고, $3\lambda$ = 9 nm에서 기판 신호는 원래의 5%로 줄어든다. 이 너머에서는 기판 피크가 배경에 묻혀 두께를 잴 수 없다. 검출 각도를 60°로 기울이면(검은 점선) 같은 두께에서 강도비가 훨씬 커져 얇은 막에 민감해진다. 대신 잴 수 있는 범위는 절반인 4.5 nm로 줄어든다.

이 식이 성립하려면 다섯 가지 가정이 필요하다.

1. 오버레이어가 균일하고 평탄하다.
2. 오버레이어와 기판의 계면이 급격하다.
3. $\lambda$와 $R_0$ 값을 안다.
4. 탄성 산란을 무시한다(또는 $\lambda$ 자리에 실효 감쇠길이를 쓴다).
5. 표면 오염층이 없거나, 있어도 두 신호에 같은 영향을 준다.

가정 4는 앞 문단처럼 $\lambda$ 자리에 IMFP 대신 EAL을 넣어 처리한다. 이후의 절은 나머지 네 가정이 하나씩 어디서 깨지는지를 따라간다.

## 2. 실무 공식으로 두께를 계산할 때 어디서 틀리는가

실무에서 실제로 쓰는 식은 어떤 모양이고, 가정 3($\lambda$와 $R_0$를 안다)은 거기서 얼마나 무거운가. 산화막 두께 계산에 가장 널리 쓰이는 것은 Strohmeier(1990)의 형태다.

$$ d = \lambda_o \sin\theta' \, \ln\!\left( \frac{N_m \lambda_m I_o}{N_o \lambda_o I_m} + 1 \right) $$

$N_m, N_o$는 금속과 산화물 속 금속 원자의 부피 밀도이고, $\lambda_m$은 1절과 달리 금속 기판 안의 감쇠길이다. 여기서 $$\theta'$$는 **표면 평면에서 잰** 출사각이다. 2편 5절에서 경고한 각도 혼동이 바로 여기서 나온다. 1절의 식과 비교하면 $$\sin\theta' = \cos\theta$$, 곧 $$\theta' = 90° - \theta$$이고, 로그 안의 $N_m\lambda_m/(N_o\lambda_o)$는 $1/R_0$의 이론적 추정이다. 두 식은 같은 물리다. 그러나 분석기가 법선에서 30° 기운 장비에서 이 식의 $$\theta'$$에 30°를 그대로 넣으면, 강도비 1.00에서 두께가 1.97 nm가 아니라 1.14 nm로 나온다. 식을 가져다 쓸 때 각도의 기준면부터 확인해야 하는 이유다.

더 깊은 함정은 $R_0$에 있다. Strohmeier 식은 $R_0$를 원자 밀도와 감쇠길이의 비로 계산한다. Seah와 Spencer(2002)는 SiO₂/Si에서 이 값을 두 가지 경로로 직접 측정했는데, 결과는 0.88 ± 0.03이었다. 계산값 0.53 ± 0.05와 크게 다르다. 원인으로는 IMFP와 SiO₂ 밀도의 계통 불확도, 그리고 3편에서 본 셰이크업 위성이 주 피크 밖으로 가져가는 몫을 꼽는다. 이 차이가 두께에 주는 영향은 크다. 같은 강도비 1.00에서 실측 $R_0$로 계산하면 2.28 nm, 계산 $R_0$로 계산하면 3.18 nm다.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/fig2-thickness-sensitivity.png" alt="두께 역산의 통계 오차와 계통 오차 비교: 강도비 오차와 계수 잡음, λ와 R0의 영향" width="720">
_그림2. 두께 역산의 정밀도와 정확도 (합성 계산, $\lambda$ = 3 nm, $$R_0$$ = 0.88). (a) 통계 오차: 강도비에 5% 오차를 넣은 경우(보라)와 기판 계수 10⁴(파란 점선)·10⁵(파란 실선)의 계수 잡음에서 오는 두께의 표준편차. (b) 계통 오차: $\lambda$를 ±10% 잘못 알 때(녹색 띠)와 실측 $$R_0$$ 대신 계산값 0.53을 쓸 때(빨간 선). 통계 오차는 12 nm까지 0.25 nm 아래에 머물지만, 계산 $$R_0$$에서 오는 계통 오차는 1 nm 막에서 이미 약 0.5 nm다._

그림2는 두 종류의 오차를 나란히 놓는다. 왼쪽의 통계 오차는 작다. 기판 피크에서 10⁵ 계수를 모으는 측정이라면, 계수 잡음에서 오는 두께의 표준편차는 0.5 nm 막에서 0.004 nm, 10 nm 막에서도 0.049 nm다. 기판 계수가 10⁴로 열 배 적으면(파란 점선) 잡음은 약 3배로 커져 10 nm 막에서 0.156 nm가 된다. 강도비에 5%의 오차를 일부러 넣어도 두께 오차는 0.023 nm에서 시작해 0.15 nm 근처에서 포화한다. 두께가 커지면 강도비가 지수적으로 커져, 같은 비율의 오차가 두께에서는 점점 작은 몫이 되기 때문이다.

오른쪽의 계통 오차는 이야기가 다르다. $\lambda$를 10% 크게 잡으면 두께도 정확히 10% 커진다. 식에서 $d$가 $\lambda$에 정비례하기 때문이다. 계산 $R_0$를 쓰면 1 nm 막이 1.51 nm(+51%), 4 nm 막이 5.19 nm(+30%), 8 nm 막이 9.44 nm(+18%)로 나온다. 얇은 막일수록 상대 오차가 크다.

이 절의 수치와 4절의 불확도 띠는 [`_code/xps-thickness-arxps-pitfalls/`](https://github.com/eddykim/eddykim.github.io/tree/main/_code/xps-thickness-arxps-pitfalls)의 `overlayer.py`가 출력하고, 띠를 넘는 각도의 수와 특이값은 각각 `generate_figures.py`와 `arxps_svd.py`가 출력한다.

Seah(2005)는 SiO₂/Si에서 이 대비를 실측으로 정리했다. XPS의 두께 정밀도는 1표준편차 0.025 nm까지 좋아질 수 있지만, 과거에는 감쇠길이 불확도에서 생기는 두께의 약 20% 계통 오차가 그 정밀도를 가렸다. 감쇠길이를 새로 보정한 뒤에는 95% 신뢰도에서 약 2%까지 줄일 수 있으나, 일반 실험실의 절차 편차 때문에 실제로는 0.4 nm까지 나빠지는 일이 흔하다. 앞의 Seah와 Spencer(2002)는 기판의 결정 구조 때문에 강도가 방향에 따라 달라져, 수직 방출로 얻은 데이터가 두께를 18% 작게 준다는 점도 보고했다. 가정 1과 3이 실제로 얼마나 무거운지 보여 주는 숫자다.

반대로 가정 5, 오염층은 이 방법에서 거의 문제가 되지 않는다. 흡착 탄소층이 두 신호를 함께 감쇠시키지만, Si 2p의 산화물 성분과 기판 성분은 운동에너지가 4 eV밖에 차이 나지 않는다. 탄화수소층 속 감쇠길이가 4.449 nm와 4.460 nm로 사실상 같아서, 1 nm 오염층이 2 nm 두께의 역산값에 주는 영향은 1 pm도 안 된다. 같은 원소의 두 화학상태를 비교하는 방법이 강한 이유다. 이 성질은 7절에서 다시 중요해진다.

## 3. 막과 기판이 다른 원소이거나 입자 모양이면 어떻게 재는가

오버레이어와 기판이 아예 다른 원소면 어떻게 하는가. 예를 들어 실리콘 위의 금속막이나 금속 위의 고분자 코팅이다. 이 경우 두 피크의 운동에너지가 달라 $\lambda_o \ne \lambda_m$이고, 1절의 깔끔한 역산식을 쓸 수 없다. 두께를 구하려면 식을 수치적으로 풀어야 한다. 판단부터 말하면, 계산 자체는 어렵지 않게 해결되지만 형상이 평면이 아니면 같은 강도비가 전혀 다른 두께를 뜻한다.

Cumpson(2000)의 두께 도표(Thickogram)는 이 계산을 그래프 한 장으로 바꾼 방법이다. 두 피크의 강도비와 감도인자 비, 운동에너지 비를 그래프 위에 놓고 두께를 읽는다. 강도비를 쓰므로 장비 인자가 상쇄되고, 두 피크의 운동에너지가 가까우면 오염층의 영향도 줄어든다. 다만 정확도는 결국 입력한 감쇠길이 값의 정확도가 정한다. 2절에서 본 대로 그 값에 두께가 정비례하기 때문이다.

평면이 아닌 형상으로 가면 문제가 더 어려워진다. 가정 1(균일하고 평탄한 층)이 깨지는 대표적인 경우다. 코어-셸(core–shell) 나노입자처럼 구형 껍질이 코어를 감싼 구조(그림3)에서는 광전자가 지나는 경로 길이가 위치마다 다르다. Shard(2012)는 이런 코어-셸 입자의 XPS 데이터에서 껍질 두께를 구하는 간단한 방법을 제시했다. 그림3의 껍질은 XPS 정보깊이(약 10 nm)보다 훨씬 두꺼워, 이 영상은 구조를 보여 주는 개념 예시일 뿐 Shard의 방법을 적용하는 사례가 아니다.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/ext-core-shell-nanoparticle-tem.jpg" alt="실리카 껍질이 산화철 나노입자 덩어리를 감싼 코어-셸 입자의 투과전자현미경 영상" width="600">
_그림3. 실리카 껍질이 마그헤마이트(γ-Fe₂O₃) 나노입자 덩어리를 감싼 코어-셸 입자의 투과전자현미경(transmission electron microscope, TEM) 영상. 가운데의 어두운 알갱이 덩어리가 코어이고, 그것을 둘러싼 회색 고리가 껍질이다. 구형 껍질에서는 코어 한가운데에서 분석기 쪽으로 나오는 광전자는 껍질을 그 두께만큼만 지나지만, 코어 가장자리에서 나오는 광전자는 껍질을 비스듬히 길게 지난다. (출처: Marko Petek, [Maghemite silica nanoparticle cluster](https://commons.wikimedia.org/wiki/File:Maghemite_silica_nanoparticle_cluster.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

## 4. ARXPS(Angle-Resolved X-ray Photoelectron Spectroscopy)로 깊이 분포를 되찾을 수 있는가

각도를 여러 개 재면 층 구조를 모두 되찾을 수 있는가. 1절의 식은 층이 하나이고 계면이 급격하다는 가정 2 위에 서 있었다. 층 구조를 모를 때는 2편 5절의 원리를 넓혀, 검출 각도를 바꿔 가며 여러 번 잰다. 이것이 각도 분해 XPS(angle-resolved XPS, ARXPS)다. 깊이 $z$에서 산화물 비율을 $c(z)$라 하면, 각도 $\theta$에서 측정한 강도는 다음과 같다.

$$ I(\theta) \propto \int_0^{\infty} c(z) \exp\!\left( -\frac{z}{\lambda\cos\theta} \right) dz $$

이 식은 $c(z)$의 라플라스 변환 꼴이다. 각도를 바꾸는 것은 변환 변수 $1/(\lambda\cos\theta)$를 바꾸는 것이다. 그러니 여러 각도의 $I(\theta)$에서 $c(z)$를 되찾는 것은 역라플라스 변환을 푸는 문제다.

그런데 지수 커널을 가진 이런 적분방정식은 전형적인 불량조건(ill-posed) 역문제다. 지수 함수는 $c(z)$의 세밀한 구조를 뭉개서 적분하므로, 구조가 꽤 다른 두 분포가 거의 같은 $I(\theta)$를 낸다. 거꾸로 풀면 데이터의 작은 잡음이 해에서 크게 증폭된다. Cumpson은 ARXPS의 깊이 분해능 한계를 오래 다뤄 왔다. 최근 프리프린트(2026)에서는 자유 형태의 깊이 프로파일 복원을 "심각한 불량조건의 역라플라스 변환 문제"로 규정한다. 1995년 논문에서 Cumpson이 정리한 바로는, 실질적인 깊이 분해능은 주로 신호 대 잡음비가 정하며 깊이의 상당 부분에 이른다.

얼마나 나쁜 조건인지는 숫자 하나로 볼 수 있다. 그림4와 같은 29개 각도(0~70°, 2.5° 간격)에서 위의 적분을 행렬로 적고, 측정마다 상대 잡음이 같도록 행을 맞춘 뒤 특이값 분해를 했다. 첫 특이값으로 나눈 값은 1, 0.15, 0.021, 0.0029, 0.0004로 단계마다 약 7분의 1씩 떨어진다. 측정 하나하나의 상대 잡음이 1%라면 잡음 위로 올라오는 성분은 셋뿐이다. 즉, 각도를 29개 재도 독립된 정보는 세 개 남짓이다. 두께가 조금씩 다른 흐림 필터 여러 장으로 같은 사진을 보는 것과 같다. 필터끼리 서로 닮아서 장수를 늘려도 새로 보이는 것이 거의 없다. 계산은 `arxps_svd.py`에 있다.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/fig3-arxps.png" alt="급격한 계면과 폭이 다른 확산 계면의 깊이 분포, 각도별 강도비 차이와 불확도 띠" width="720">
_그림4. ARXPS로 계면 구조를 가를 수 있는가 (합성 계산). (a) 2 nm의 급격한 계면과 폭이 다른 확산 계면. 폭은 오차함수를 이루는 가우시안(Gaussian)의 표준편차 $\sigma$다. (b) 각도별 강도비의 차이. 회색 띠는 검출 각도 ±0.5° 오차와 계수 잡음(기판 10⁵ 계수)을 합친 불확도다. $\sigma$ = 0.5 nm 이하의 확산 계면은 모든 각도에서 띠 안에 든다._

그림4가 이 불량조건성을 숫자로 보여 준다. 2 nm의 급격한 계면과, 계면이 오차함수 모양으로 퍼진 확산 계면을 비교했다. 계면 폭은 오차함수를 이루는 가우시안의 표준편차 $\sigma$로 나타냈고, 농도가 10%에서 90%로 바뀌는 구간은 그 약 2.6배다. 확산 계면마다 중심 위치를 조정했다. 또 $R_0$의 실측 불확도(±0.03/0.88, 약 ±3.4%) 안에서 강도비 전체에 곱해지는 상수도 조정해, 급격한 계면과 가장 비슷해지도록 맞췄다. 불확도 띠는 Seah(2005)가 한계 요인으로 꼽은 각도 설정 오차 0.5°와 계수 잡음을 합친 것으로, 0°에서 ±0.65%, 40°에서 ±1.29%, 70°에서 ±5.80%다.

결과는 이렇다. $\sigma$가 0.3 nm나 0.5 nm인 확산 계면은 0~70°의 모든 각도에서 띠 안에 들어온다. 급격한 계면과 구분할 수 없다. $\sigma$가 0.8 nm로 커져야 29개 각도 가운데 21개에서, 1.2 nm에서는 25개에서 띠를 넘는다. 2 nm 깊이의 계면에서 0.5 nm 이하의 구조는 각도 데이터에 사실상 흔적을 남기지 않는 셈이다. 이 결론은 $R_0$ 불확도를 끝까지 활용했을 때의 것이다. $R_0$를 정확히 안다면 구분 한계는 조금 내려가지만, 2절에서 본 대로 그 값 자체가 문헌마다 크게 달랐다.

그래서 실무에서는 임의의 $c(z)$를 자유롭게 되찾으려 하지 않는다. 층의 개수와 순서를 미리 가정한 층 모델을 피팅하거나, 해가 매끄럽도록 정규화(regularization)를 건다. [최적화 방법론 3편](/posts/optimization-levenberg-marquardt/)의 Levenberg-Marquardt법이 $J^T J$에 대각 성분을 더해 나쁜 조건수를 누그러뜨린 것과 같은 발상이다. 층 모델 피팅은 [최적화 방법론 4편](/posts/optimization-global-heuristics/)에서 본 것처럼 초기값에 따라 다른 국소 최솟값으로 수렴하기 쉽다. 각도 데이터만으로 복잡한 깊이 분포를 자유롭게 되찾았다는 결과를 보면, 어떤 가정과 정규화가 들어갔는지부터 확인해야 한다.

비파괴로 깊이를 보는 방법이 하나 더 있다. [3편](/posts/xps-spectrum-structure-satellites/) 7절과 4편 4절에서 예고한 Tougaard의 배경 분석이다. 광전자가 에너지를 잃을 확률은 고체 안에서 지나온 거리에 달려 있으므로, 피크 아래 배경의 모양에는 그 원소가 어느 깊이에 어떻게 분포하는지가 담겨 있다. 표면에만 있는 원소는 에너지 손실이 적어 피크 뒤 배경이 낮고, 깊이 묻힌 원소는 손실 꼬리가 크다. Tougaard는 비탄성 산란 단면적으로 여러 깊이 분포의 배경 모양을 계산해 측정 스펙트럼과 맞춰 보는 방법을 발전시켰다. 이렇게 하면 각도를 바꾸지 않고 한 번의 측정에서 나노 규모의 깊이 구조를 추정할 수 있다. 이 방법도 결국은 배경 모델과 단면적이라는 모델 위에 서 있으며, 깊이 분포를 거꾸로 푸는 역문제라는 점은 각도 방법과 같다.

ARXPS와 배경 분석의 공통 장점은 시편을 깎지 않는다는 점이다. 한계는 보는 깊이다. ARXPS는 에너지를 잃지 않은 광전자의 정보깊이, 곧 약 $3\lambda$(Si 2p에서 약 10 nm) 안만 본다. 배경 분석은 에너지를 잃은 전자까지 쓰므로 조금 더 깊이 닿는다. Tougaard(1998)에 따르면 주된 정보는 약 $5\lambda$ 안에서, 구조 파라미터 수준의 정보는 $5\lambda$~$10\lambda$까지에서 온다. 그보다 깊은 구조는 어느 방법으로도 닿지 않는다.

## 5. 이온빔으로 깎아 내려가면 더 깊이 볼 수 있는가

더 깊이 보려면 깎으면 되지 않는가. 가장 직접적인 방법은 아르곤 이온빔으로 표면을 조금씩 깎아 내며 XPS를 반복 측정하는 것이다. 깎은 시간을 깊이로 환산하면 조성 대 깊이 곡선이 나온다. 정보깊이의 한계를 넘어 수백 nm까지 볼 수 있고, 2편에서 걷어내고 싶다고 했던 흡착 탄소층도 걷어낼 수 있다.

문제는 이온빔이 조성을 바꾼다는 점이다. 구성 원소가 같은 비율로 떨어져 나가지 않는 현상을 우선 스퍼터링(preferential sputtering)이라 한다. 산화물에서는 산소가 먼저 빠지는 경향이 있다. Counsell 등(2014)은 비정질 TiO₂를 5 kV 단원자 Ar⁺ 이온으로 깎으면 산소가 우선 제거되며, Ti가 +4 상태에서 +3과 +2 상태로 환원된다고 보고했다. 원래 시편에 없던 Ti³⁺와 Ti²⁺가 스펙트럼에 나타나는 것이다. 이것을 "계면에 환원된 산화물이 있다"고 읽으면 측정이 만든 인공물(artifact)을 시편의 성질로 보고하는 셈이다.

그 밖에도 인공물이 여럿 있다. 입사 이온이 원자를 밀어 넣어 층이 섞이는 원자 혼합(atomic mixing), 깎을수록 거칠어지는 표면, 층마다 다른 스퍼터 속도 때문에 깎은 시간을 깊이로 바꾸는 눈금이 층마다 달라지는 문제가 있다. 여기에 XPS 자체의 정보깊이가 더해진다. 각 단계에서 측정하는 것은 새로 드러난 표면 한 점이 아니라, 그 아래 수 nm를 지수 가중 평균한 값이다.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/fig4-sputter-artifact.png" alt="참 산소 분포와 원자 혼합, 정보깊이 가중, 산소 우선 제거를 거친 관측 프로파일" width="640">
_그림5. 스퍼터 깊이 프로파일이 계면을 뭉개는 과정 (개념 모사). 두께 5 nm의 MO₂ 산화막 아래 금속 M. 원자 혼합 폭 1 nm, 정보깊이 가중의 감쇠길이 2 nm(2편의 정의로 정보깊이는 약 6 nm), 산소 우선 제거 30%는 설명을 위한 임의의 값이다. 빨간 선이 세 효과를 모두 거친 최종 관측 프로파일이다._

그림5는 이 효과들을 차례로 얹은 개념 모사다. 참 분포는 5 nm에서 급격히 끝나지만, 원자 혼합이 계면을 퍼뜨리고, 정보깊이의 지수 가중이 그 위를 한 번 더 뭉개며, 산소 우선 제거가 산화막 전체의 산소 분율을 끌어내린다. 최종 관측 프로파일(빨간 선)만 보면 산화막이 실제보다 얇고 계면이 넓으며, 산화막이 원래 조성보다 산소가 모자란 것처럼 보인다.

완화책도 있다. 이온 에너지를 낮추고, 시편을 돌려 가며 깎아 거칠어짐을 줄인다. 유기물이나 고분자, 그리고 산화물의 화학상태를 보존하고 싶을 때는 수백에서 수천 개의 아르곤 원자가 뭉친 가스 클러스터 이온(gas cluster ion, Ar$_n^+$)을 쓴다. 원자 하나당 에너지가 작아 깊이 침투하지 않기 때문이다. Counsell 등도 클러스터 이온이 TiO₂의 손상과 이온 혼입을 크게 줄인다고 보고했다. 그래도 핵심 메시지는 바뀌지 않는다. 스퍼터 깊이 프로파일은 조성의 추세를 보는 데는 쓸 만하지만, 깎은 뒤에 읽은 화학상태를 원래 시편의 화학상태로 보고해서는 안 된다.

## 6. X선은 시편을 정말 건드리지 않는가

X선은 비파괴 아닌가. XPS를 소개하는 글에는 흔히 비파괴 분석이라는 말이 붙는다. Morgan(2023)은 이 통념을 반박한다. X선 광자와 그 뒤에 이어지는 전자 폭포(electron cascade)가 분석 영역을 상당히 바꿀 수 있다는 것이다. Morgan이 보인 PVC 손상 데이터는 단색화 X선 장비에서 얻은 것이다.

특히 환원되기 쉬운 화학종이 취약하다. Cardiff 대학의 XPS 참고 자료는 Au(III) 화합물이 광방출 과정에서 광환원되어, 10분 조사 뒤에는 Au(III)와 Au(0)가 함께 보인다고 적는다. Biesinger 등(2010)은 V₂O₅가 X선 아래에서 서서히 V(IV)로 바뀌어, 210 W 광원에서 24시간 동안 15% 넘게 전환되었다고 보고했다. 3편에서 위성 구조로 판별했던 Cu(II)도 X선에 의한 환원을 고려해야 하는 화학종이다. 고분자에서는 결합 절단이 흔하다.

실무 절차는 단순하다. 민감한 영역은 측정의 처음과 끝에 두 번 찍는다. Cardiff 자료가 Au 4f 영역에 권하는 방식도 같다. 그 영역을 먼저 따로 찍고, 나머지 영역을 다 찍은 뒤 다시 찍어 환원 정도를 확인한다. 두 스펙트럼이 다르면 그 데이터로 화학상태를 결론짓지 않는다. 조사 시간을 줄이거나 위치를 옮겨 가며 찍는 것도 방법이다.

1편부터 4편까지는 모델이 어디서 깨지는지를 봤다. 여기서는 모델 이전에 측정 대상이 측정 도중에 바뀐다. 가장 근본적인 형태의 모델 붕괴다. 어떤 모델도 측정하는 동안 바뀐 시편에서 원래 시편을 되찾아 주지 않는다.

## 7. 엘립소메트리와 XPS는 서로 무엇을 보완하는가

같은 SiO₂/Si 두께를 엘립소메트리로도 XPS로도 잴 수 있는데, 무엇이 다른가. [2편](/posts/xps-surface-sensitivity-imfp/) 7절에서 두 기법이 보는 깊이와 화학 정보가 어떻게 다른지 비교했다. 여기서는 두께를 내는 데 각각 무엇이 필요한지를 놓고, 같은 산화막의 두께 값이 왜 서로 어긋나는지를 본다.

| 항목 | 엘립소메트리 | XPS |
|---|---|---|
| 측정하는 것 | 반사광의 편광 상태 변화 ($\Psi, \Delta$) | 광전자의 운동에너지 분포 |
| 측정 환경 | 대기 중 | 초고진공 |
| 속도 | 초 단위, 공정 중 측정 가능 | 분 단위 이상 |
| 두께에 필요한 사전 지식 | 층 구조와 각 층의 광학 상수 | 층 구조, 감쇠길이, $R_0$ |

마지막 행이 핵심이다. 엘립소메트리는 광학 상수를 알아야 두께가 나오고, XPS는 $\lambda$와 층 구조를 알아야 두께가 나온다. 둘 다 모델 역산이고, 필요한 사전 지식의 종류가 다를 뿐이다.

이 차이는 실측으로도 드러났다. 국제도량형위원회 산하 물질량 자문위원회(CCQM)의 예비 비교 연구에서, Seah 등(2004)은 1.5~8 nm의 열산화막을 10가지 방법으로 잰 45개 결과를 비교했다. 모든 방법이 기준값과 잘 맞는 직선 관계를 보였지만, 방법마다 고유한 영점 오프셋이 0~1 nm 있었다. 엘립소메트리는 표면의 물과 탄소 오염을 약 1 nm의 산화막 두께로 보았고, 중에너지 이온 산란(medium-energy ion scattering, MEIS), 핵반응 분석(nuclear reaction analysis, NRA), 러더퍼드 후방산란 분광(Rutherford backscattering spectrometry, RBS) 같은 이온빔 방법은 흡착된 산소를 약 0.5 nm로 보았다. Si 2p를 쓰는 XPS만 오프셋이 0이었다. 2절에서 본 대로 오염층이 두 Si 2p 성분을 똑같이 감쇠시키기 때문이다. 대신 XPS는 척도 상수, 곧 감쇠길이의 불확도가 컸다. 그래서 이 연구는 XPS의 영점과 다른 방법들의 척도를 조합해, 기준 감쇠길이에 0.986 ± 0.009를 곱하라고 권고했다. 두 계측이 서로의 약점을 메운 셈이다.

## 정리

이번 편에서는 감쇠를 거꾸로 이용해 두께를 재는 방법을 따라갔다. 균일 오버레이어 모델은 강도비 하나에서 두께를 준다. 통계 정밀도는 기판 10⁵ 계수에서 0.004~0.05 nm로 좋지만, 정확도는 $\lambda$와 $R_0$라는 모델 상수가 정하며, $R_0$를 계산값으로 쓰면 얇은 막의 두께가 50% 넘게 틀린다. 각도 규약을 헷갈리면 같은 데이터에서 전혀 다른 두께가 나온다. ARXPS는 지수 커널의 역문제라서, 이 글의 합성 조건(감쇠길이 3 nm, 오차함수 계면, 계수 10⁵)에서는 2 nm 깊이 계면의 폭($\sigma$) 0.5 nm 이하, 곧 10~90% 전이 구간 약 1.3 nm 이하의 계면 구조는 각도 데이터에 사실상 흔적을 남기지 않는다. 스퍼터 깊이 프로파일은 더 깊이 보여 주지만 우선 스퍼터링으로 원래 없던 화학상태를 만들고, X선 자체도 환원되기 쉬운 화학종을 바꾼다.

## 시리즈를 돌아보며

다섯 편을 돌아보면 같은 질문이 매번 모양을 바꿔 나왔다. [1편](/posts/xps-photoemission-binding-energy/)에서는 운동에너지를 결합에너지로 바꾸는 데 에너지 보존식과 페르미 준위 정렬이라는 가정이 필요했고, 측정된 결합에너지는 이완 때문에 궤도에너지와도 달랐다. [2편](/posts/xps-surface-sensitivity-imfp/)에서는 "표면 10 nm"가 95%라는 관례와 예측식으로 계산한 $\lambda$의 곱이었다. [3편](/posts/xps-spectrum-structure-satellites/)에서는 스펙트럼의 어떤 구조를 성분으로 셀지가 이미 해석이었고, 위성과 다중항을 잘못 세면 없는 화학종이 생겼다. [4편](/posts/xps-quantification-peak-fitting/)에서는 조성 숫자 하나가 대전 보정, 배경, 선형, 감도인자라는 네 모델을 통과했고, 배경 모델이 틀리면 오차 막대가 실제 오차를 전혀 담지 못했다. 그리고 이번 편에서 두께는 그 모든 모델 위에 층 구조와 $\lambda$라는 가정을 하나 더 얹은, 가장 간접적인 양이었다.

XPS가 실제로 재는 것은 처음부터 끝까지 광전자의 운동에너지 분포 하나였다. 원소도, 화학상태도, 조성도, 두께도 그 분포에 모델을 씌워 거꾸로 계산한 값이다. 이 구조는 XPS만의 것이 아니다. [엘립소메트리 배경이론 1편](/posts/ellipsometry-electromagnetic-fresnel/)의 출발점은 검출기가 세기만 재기 때문에 전기장의 위상이 사라진다는 사실이었다. 엘립소메트리는 편광 상태를 비교해 그 정보를 되찾았다. 측정 원리는 다르지만, 측정량과 알고 싶은 양 사이에 모델이 끼어 있고 그 모델의 가정이 결과를 좌우한다는 점은 같다.

각 모델의 가정을 결과와 함께 적어 두면, XPS의 숫자를 어디까지 믿을지 판단할 수 있다. 모델이 메우는 그 간극을 잊는 순간, 그럴듯한 숫자만 남는다.

## 참고자료

- M. P. Seah and S. J. Spencer, "Ultrathin SiO₂ on Si II. Issues in quantification of the oxide thickness," *Surf. Interface Anal.* 33, 640–652 (2002). <https://doi.org/10.1002/sia.1433>
- M. P. Seah et al., "Critical review of the current status of thickness measurements for ultrathin SiO₂ on Si Part V: Results of a CCQM pilot study," *Surf. Interface Anal.* 36, 1269–1303 (2004). <https://doi.org/10.1002/sia.1909>
- M. P. Seah, "Ultrathin SiO₂ on Si. VI. Evaluation of uncertainties in thickness measurement using XPS," *Surf. Interface Anal.* 37, 300–309 (2005). <https://doi.org/10.1002/sia.2020>
- B. R. Strohmeier, "An ESCA method for determining the oxide thickness on aluminum alloys," *Surf. Interface Anal.* 15, 51–56 (1990). <https://doi.org/10.1002/sia.740150109>
- HarwellXPS, "Technical Note #1: Oxide Thickness Determination Using the Strohmeier Equation." <https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf>
- P. J. Cumpson, "The Thickogram: a method for easy film thickness measurement in XPS," *Surf. Interface Anal.* 29, 403–406 (2000). <https://doi.org/10.1002/1096-9918(200006)29:6%3C403::AID-SIA884%3E3.0.CO;2-8>
- A. G. Shard, "A straightforward method for interpreting XPS data from core–shell nanoparticles," *J. Phys. Chem. C* 116, 16806–16813 (2012). <https://doi.org/10.1021/jp305267d>
- P. J. Cumpson, "Angle-resolved XPS and AES: Depth-resolution limits and a general comparison of properties of depth-profile reconstruction methods," *J. Electron Spectrosc. Relat. Phenom.* 73, 25–52 (1995). <https://doi.org/10.1016/0368-2048(94)02270-4>
- P. J. Cumpson, "Resolution Limits in Constrained Angle-Resolved XPS Depth-Profile Reconstruction," ChemRxiv (2026). <https://doi.org/10.26434/chemrxiv.15007123/v1>
- S. Tougaard, "Accuracy of the non-destructive surface nanostructure quantification technique based on analysis of the XPS or AES peak shape," *Surf. Interface Anal.* 26, 249–269 (1998). <https://doi.org/10.1002/(SICI)1096-9918(199804)26:4%3C249::AID-SIA368%3E3.0.CO;2-A>
- J. D. P. Counsell, A. J. Roberts, W. Boxford, C. Moffitt, and K. Takahashi, "Reduced Preferential Sputtering of TiO₂ using Massive Argon Clusters," *J. Surf. Anal.* 20, 211–215 (2014). <https://www.jstage.jst.go.jp/article/jsa/20/3/20_211/_article>
- D. J. Morgan, "XPS insights: Sample degradation in X-ray photoelectron spectroscopy," *Surf. Interface Anal.* 55, 331–335 (2023). <https://doi.org/10.1002/sia.7205>
- Cardiff University XPS Access, "Gold." <https://sites.cardiff.ac.uk/xpsaccess/reference/gold/>
- M. C. Biesinger, L. W. M. Lau, A. R. Gerson, and R. St. C. Smart, "Resolving surface chemical states in XPS analysis of first row transition metals, oxides and hydroxides: Sc, Ti, V, Cu and Zn," *Appl. Surf. Sci.* 257, 887–898 (2010). <https://doi.org/10.1016/j.apsusc.2010.07.086>
- C. J. Powell, "Practical guide for inelastic mean free paths, effective attenuation lengths, mean escape depths, and information depths in x-ray photoelectron spectroscopy," *J. Vac. Sci. Technol. A* 38, 023209 (2020). <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=929257>

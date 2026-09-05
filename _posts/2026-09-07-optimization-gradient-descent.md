---
title: 최적화 방법론 1편 — 최소자승과 Gradient Descent
date: 2026-09-07 20:00:00 +0900
categories: [계산과 알고리즘, 최적화방법]
page_id: optimization-gradient-descent
tags: [optimization, gradient-descent, least-squares, thin-film, python]
description: 수치해석과 최적화 이론이 무엇을 다루는 분야인지부터 시작해, 박막 반사율 피팅을 예로 gradient descent가 어떻게 동작하고 어디서 실패하는지를 다룬다.
math: true
---

계측 데이터에서 원하는 물리량을 뽑아내려면 결국 최적화 문제를 풀어야 하는 경우가 많다. 이 시리즈는 그러한 반복법을 gradient descent부터 순서대로 정리한다. 본론에 앞서 수치해석과 최적화 이론이 각각 무엇을 다루는 분야인지 정의한다.

## 1. 수치해석과 최적화 이론

수치해석(numerical analysis)은 해석적으로 정확한 해를 구하기 어렵거나 아예 존재하지 않는 수학 문제를, 유한한 계산으로 근사해를 구하는 방법론을 다루는 분야다. 방정식의 근을 찾는 문제, 함수를 적분하거나 미분방정식을 푸는 문제, 데이터 사이를 보간하는 문제가 모두 여기 속한다. 반도체 소자의 특성을 미리 계산하는 것부터 유체의 흐름을 예측하는 것까지, 대수적으로 정확한 해가 없는 문제를 다루는 공학 전반이 수치해석에 기대고 있다.

최적화(optimization) 이론은 이 중에서도 "어떤 목적함수(objective function) $J(x)$를 최소화(또는 최대화)하는 파라미터 $x$를 찾는" 문제를 다루는 갈래다. 일반형으로 쓰면 다음과 같다.

$$ x^* = \arg\min_x J(x) $$

$J$가 선형이거나 2차식이면 $x^*$를 대수적으로 바로 구할 수 있는 경우도 있다. 하지만 $J$가 비선형이거나 $x$의 차원이 커지면, 대수적으로 풀리는 경우가 오히려 드물다. 이때 사용하는 것이 반복법(iterative method)이다. 초기 추정값 $x_0$에서 출발해서, 매 단계 $J(x)$가 줄어드는 방향으로 $x$를 갱신해가며 최솟값에 접근하는 방식으로, 한 번의 계산으로 정확한 해를 산출하는 직접법(direct method)과 대비된다.

이 시리즈에서는 gradient descent, Newton법, Gauss-Newton법, 그리고 계측 분야에서 사실상 표준으로 쓰이는 Levenberg-Marquardt(LM)법을 차례로 다룬다. 이번 편의 대상은 그중 가장 단순한 gradient descent이다.

## 2. 간접측정과 반복법이 필요한 경우

수치적 최적화가 필요한 이유는 하나로 수렴한다. 모델이 비선형이면 측정값으로부터 그 원인이 되는 파라미터를 구하는 역함수가 일반적으로 존재하지 않기 때문이다. 이는 광학에 국한된 문제가 아니다. 물리학이든 통계학이든 현상을 기술하는 모델은 비선형성을 갖는 경우가 흔하며, 그때마다 동일한 문제에 직면한다. 선형 모델은 행렬의 역을 구해 바로 풀리지만 비선형 모델은 그렇지 않다.

대표적인 예가 박막(thin film) 반사율 스펙트럼으로부터 두께를 구하는 문제다. 모델식이 주어져 있으므로 역함수를 구하면 된다고 생각하기 쉽고, 실제로 전혀 근거가 없는 기대는 아니다. 소광계수(extinction coefficient)가 $k=0$인 흡수 없는 투명 박막에 한해서는, 반사율 스펙트럼의 간섭 무늬(interference fringe) 극값으로부터 두께를 닫힌 형태로 구하는 방법(envelope method 계열)이 존재한다. 그러나 이는 특수한 경우에 그친다. 소광계수가 0이 아니면 굴절률이 복소수 $N = n - ik$가 되며, 무엇보다 구하려는 두께 $d$가 측정량인 반사율 $R$에 대해 비선형으로 결합되어 있어 $R$로부터 $d$를 직접 얻는 닫힌 형태의 역함수는 일반적으로 존재하지 않는다.

이 문제는 간접측정(indirect measurement)의 한 예다. 간접측정이란 알고 싶은 양을 직접 잴 수 없을 때, 측정 가능한 다른 물리량과 그 파라미터 사이의 관계(모델)를 세우고, 모델을 거꾸로 풀어 파라미터를 추정하는 방식이다. 두께를 자로 재는 대신 반사율 스펙트럼 $R(\lambda)$를 측정하고, 이론 모델 $R_{model}(d,\lambda)$을 거꾸로 풀어 두께 $d$를 구하는 식이다. 모델이 선형이고 닫힌 형태의 역함수가 존재하면 대수적으로 해를 구할 수 있으나(직접법), 그렇지 않은 경우가 훨씬 많다. 이때 1절에서 정의한 반복법이 필요하다.

구체적인 구조로, 공기(air, ambient layer — 빛이 입사하는 매질) 위에 SiO2(관심층, layer of interest — 두께를 구하고 싶은 층) 하나가 있고, 그 아래를 Si 기판(substrate layer)이 받치는 구조를 쓴다. 층이 하나뿐인 이런 구조를 단층 샘플(single-layer sample)이라고 부른다. 구조는 단순하지만 반사율은 두께에 대해 지수함수(간섭항)를 포함하는 비선형식으로 주어진다. 수직입사 조건에서는 다음과 같은 닫힌 형태(Airy 공식)로 쓸 수 있다.

$$ r = \frac{r_{01} + r_{12} e^{-2i\beta}}{1 + r_{01} r_{12} e^{-2i\beta}}, \qquad \beta = \frac{2\pi n_1 d}{\lambda} $$

$r_{01}, r_{12}$는 각 계면(ambient/관심층, 관심층/substrate)의 Fresnel 반사계수, $R = \lvert r \rvert^2$이 실제로 측정하는 반사율이다(구현은 `_code/optimization-gradient-descent/reflectance_model.py`). 이 식은 $d$에 대해 닫힌 형태로 역을 취할 수 없다. 또한 측정값에는 노이즈가 포함되어 있으므로 모델과 정확히 일치하는 $d$는 존재하지 않으며, 결국 가장 잘 맞는 $d$를 찾는 문제가 된다.

여기서 최소자승(least squares) 목적함수가 필요해진다.

$$ J(d) = \frac{1}{2}\sum_i \left( R_{model}(d, \lambda_i) - R_{meas}(\lambda_i) \right)^2 $$

잔차(residual)를 제곱하여 합하는 이유는 무엇인가. 절댓값 합을 비롯한 다른 손실함수도 사용할 수 있으나, 측정 노이즈가 가우시안 분포를 따른다고 가정하면 최소자승해가 최대우도추정(maximum likelihood estimation)과 일치한다는 통계적 근거가 있다. 계측 노이즈는 여러 독립적인 잡음원의 합이므로 중심극한정리에 의해 가우시안에 근사하는 경우가 많으며, 이러한 이유로 최소자승이 이 시리즈 전체의 출발점이 된다.

목표는 모델과 측정값의 차이를 제곱하여 합한 값을 최소화하는 $d$를 찾는 것이다. 초기 모델과 측정값의 불일치를 그림으로 나타내면 다음과 같다.

<img src="/assets/img/posts/optimization-gradient-descent/fig1-model-vs-measurement.png" alt="모델 초기값 vs 측정치" width="600">
_그림1. 두께 1540nm로 가정한 초기 모델과 실제 1490nm 샘플의 측정치_

측정치는 실제 스펙트로미터 raw 데이터가 아니라, 과거 실험실에서 실측했던 SiO2/Si 두께 계열(10~190nm, 1490nm) 중 1490nm 샘플을 가정하고 같은 물리 모델에 가우시안 노이즈를 더해 합성한 것이다. 1490nm는 가시광 대역 안에 간섭 무늬(fringe)가 여러 개 들어갈 만큼 광학적으로 두꺼워, 두께 변화가 반사율 스펙트럼에 뚜렷하게 반영된다. 그림 1에서 초기 모델(파란 선, $d=1540$nm)과 측정값(회색 점)의 fringe 위치가 어긋나 있는 것이 이를 보여준다.

이 불일치를 두께 $d$의 함수로 나타내면 목적함수 $J(d)$의 형태를 얻는다.

<img src="/assets/img/posts/optimization-gradient-descent/fig2-objective-landscape.png" alt="두께에 따른 목적함수" width="600">
_그림2. 두께에 따른 목적함수 J(d) — 회색 점선이 전역 최솟값(global minimum)_

다만 그림 2와 같이 전 구간을 그려 눈으로 최솟값을 찾는 것은 실제 절차가 아니다(그림 2는 원리를 보이기 위해 넓은 구간을 미리 탐색한 결과일 뿐이다). $d$ 하나를 평가할 때마다 전 파장 대역의 반사율을 계산해야 하므로 이러한 격자 탐색(grid search)은 계산 비용이 크다. 실제 목표는 $J(d)$를 최소로 만드는 $d$를 최소한의 연산으로 찾아내는 것이다.

## 3. 목적함수와 기울기 조건

2절에서 반복법이 필요한 이유는 확인했다. 반복법이 성립하려면 매 단계마다 현재 추정값 $d_k$에서 다음 추정값 $d_{k+1}$을 어느 방향으로 얼마나 옮길지가 정해져야 한다. 임의의 방향으로 이동하면 $J$가 오히려 증가할 수 있으므로 이 결정에 근거가 필요하다. 그 근거를 제공하는 것이 기울기(gradient)다.

$J(d)$를 $d$로 미분한 $dJ/dd$는 부호로 "어느 방향"을, 크기로 "지금 얼마나 급격하게 안 맞는가"를 알려준다. 그림 1의 초기값($d=1540$nm)에서 $J$가 크다는 사실은 그 지점이 최적점이 아님을 의미하며, 해당 지점의 기울기는 $J$를 감소시키는 이동 방향을 지시한다. 이번 편에서는 파라미터가 두께 하나뿐이므로 기울기가 스칼라 미분값이지만, 다음 편들에서 파라미터가 여러 개(다층 두께, 굴절률 등)로 늘어나면 이 스칼라는 기울기 벡터(gradient vector)로, 2차 미분은 헤시안 행렬(Hessian)로 확장된다.

기울기를 따라 이동을 반복하면 결국 더 이상 이동할 방향이 없는 지점, 즉 기울기가 0인 지점에 도달한다. 따라서 $J(d)$가 최소가 되는 지점에서는 기울기가 0이어야 한다는 필요조건이 성립한다.

$$ \frac{dJ}{dd} = 0 $$

이는 필요조건일 뿐 충분조건이 아니다. $dJ/dd=0$인 지점은 전역 최솟값(global minimum)일 수도, 국소 최솟값(local minimum)일 수도 있다. 이 문제는 5절에서 다시 언급하고 다음 편(Newton법, Gauss-Newton법)에서 상세히 다룬다.

## 4. Gradient Descent

가장 단순한 접근은 현재 위치에서 기울기의 반대 방향으로 조금씩 이동하는 것이다. 반대 방향이어야 하는 근거는 1차 테일러 전개에서 직접 얻어진다.

$$ J(d + \Delta d) \approx J(d) + \frac{dJ}{dd}\Delta d $$

$\Delta d$를 기울기와 반대 부호로 잡으면($\Delta d = -\alpha \, dJ/dd$, $\alpha>0$) 우변의 둘째 항이 항상 음수가 되어 $J$가 적어도 국소적으로는 감소함이 보장된다. 이 과정을 반복하는 것이 gradient descent다.

$$ d_{k+1} = d_k - \alpha \, \frac{dJ}{dd}(d_k) $$

$\alpha$는 스텝 사이즈(step size, learning rate)다. 해석적 미분을 유도하는 대신 중심차분(central difference)으로 수치미분하였다.

```python
# gradient_descent.py 핵심부 (전체 코드: _code/optimization-gradient-descent/)
import numpy as np
from reflectance_model import reflectance

def gradient_descent(d0, wavelength_nm, measured_R, alpha, n_iter=100, h=1e-3):
    def J(d):
        return 0.5 * np.sum((reflectance(d, wavelength_nm) - measured_R) ** 2)

    d = d0
    d_hist, J_hist = [d], [J(d)]
    for _ in range(n_iter):
        grad = (J(d + h) - J(d - h)) / (2 * h)  # 중심차분으로 dJ/dd 근사
        d = d - alpha * grad
        d_hist.append(d)
        J_hist.append(J(d))
    return np.array(d_hist), np.array(J_hist)
```

먼저 $\alpha$를 충분히 크게 잡아 `alpha=1500`, 초기값 `d0=1540`nm으로 20회 반복하였다. 이 선택의 결과를 적절한 스텝 사이즈와 비교하기 위해, 동일한 초기값과 반복 횟수로 `alpha=300`을 함께 계산하여 겹쳐 그렸다.

<img src="/assets/img/posts/optimization-gradient-descent/fig3-thickness-vs-iteration.png" alt="두께 추정값 vs iteration - step size 비교" width="600">
_그림3. 두께 추정값 vs iteration - step size 비교_

<img src="/assets/img/posts/optimization-gradient-descent/fig4-objective-vs-iteration.png" alt="목적함수 J vs iteration - step size 비교" width="600">
_그림4. 목적함수 J vs iteration - step size 비교_

파란 선(`alpha=300`)은 5회 반복 만에 실제 값(1490nm, 오차 0.3nm 이내)에 도달한 뒤 그 값을 유지한다. 목적함수 역시 노이즈 바닥(noise floor)에 해당하는 약 0.0024까지 단조 감소한 뒤 평탄해진다(그림 4).

반면 빨간 선(`alpha=1500`)은 1540 → 1423 → 1527 → 1416 → 1503 순으로 매 스텝 수십에서 100nm 이상 진동하며, 20회를 모두 반복할 때까지 수렴하지 않는다(그림 3). 최솟값 부근을 지나가는 경우도 있으나(7회째 반복에서 $J\approx0.03$까지 감소, 그림 4) 다음 스텝에서 곧바로 벗어난다. 계곡을 내려갈 때 보폭이 지나치게 크면 바닥을 지나쳐 반대편 사면으로 올라가는 것과 같은 상황이다. 반대편의 기울기가 더 가파르므로 다음 스텝의 이탈 폭은 더욱 커진다. 스텝 사이즈가 과도할 때 나타나는 gradient descent의 대표적인 실패 양상이며, darkpgmr의 최적화 관련 글([darkpgmr.tistory.com/133](https://darkpgmr.tistory.com/133))에서도 같은 문제를 지적하고 있다.

스텝 사이즈 하나의 차이가 수렴과 발산을 가른 것이다.

## 5. 한계

이번 편에서 사용한 gradient descent는 두 가지 조건에 의존한다.

첫째, 초기값이 전역 최솟값으로 이어지는 영역(basin) 안에 있어야 한다. 그림 2에서 확인할 수 있듯 $d=1490$ 부근이 유일한 최솟값은 아니다. 간섭 무늬 하나만큼 떨어진 $d \approx 1290$과 $d \approx 1690$ 부근에도 각각 국소 최솟값이 존재하며, basin 경계(대략 $d=1390$과 $d=1590$, 그림 2의 두 봉우리) 바깥에서 출발하면 알고리즘은 다른 국소 최솟값으로 수렴한다. 이 국소 최솟값 문제는 이번 편에서 다루지 않으며, 다음 편(Newton법, Gauss-Newton법)에서 논의한다.

둘째, 수렴 속도가 $\alpha$에 크게 의존하며, 시행착오로 매번 $\alpha$를 결정하는 것은 비효율적이다. 최솟값에 근접할수록 $\lvert dJ/dd \rvert$가 0에 가까워져 스텝이 점차 작아지고, 반대로 초기부터 과도하게 크게 설정하면 그림 3과 같이 발산한다. 다음 편에서 다룰 Newton법은 이 스텝 사이즈를 헤시안(Hessian)으로부터 자동으로 결정하지만, 그 대가로 헤시안 계산과 역행렬 연산이라는 비용이 발생한다.

---
title: 최적화 방법론 1편 — 최소자승과 Gradient Descent
date: 2026-08-13 16:00:00 +0900
categories: [최적화, 수치해법]
tags: [optimization, gradient-descent, least-squares, thin-film, python]
description: 수치해석과 최적화 이론이 무엇을 다루는 분야인지부터 시작해, 박막 반사율 피팅을 예로 gradient descent가 어떻게 동작하고 어디서 실패하는지를 다룬다.
math: true
---

계측 데이터에서 원하는 물리량을 뽑아내려면 결국 최적화 문제를 풀어야 하는 경우가 많다. 이 시리즈에서는 그 반복법들을 gradient descent부터 하나씩 정리한다. 이번 편에 들어가기 전에, 수치해석과 최적화 이론이 정확히 뭘 다루는 분야인지부터 짚고 간다.

## 1. 수치해석과 최적화 이론

수치해석(numerical analysis)은 해석적으로 정확한 해를 구하기 어렵거나 아예 존재하지 않는 수학 문제를, 유한한 계산으로 근사해를 구하는 방법론을 다루는 분야다. 방정식의 근을 찾는 문제, 함수를 적분하거나 미분방정식을 푸는 문제, 데이터 사이를 보간하는 문제가 모두 여기 속한다. 반도체 소자의 특성을 미리 계산하는 것부터 유체의 흐름을 예측하는 것까지, 대수적으로 정확한 해가 없는 문제를 다루는 공학 전반이 수치해석에 기대고 있다.

최적화(optimization) 이론은 이 중에서도 "어떤 목적함수(objective function) $J(x)$를 최소화(또는 최대화)하는 파라미터 $x$를 찾는" 문제를 다루는 갈래다. 일반형으로 쓰면 다음과 같다.

$$ x^* = \arg\min_x J(x) $$

$J$가 선형이거나 2차식이면 $x^*$를 대수적으로 바로 구할 수 있는 경우도 있다. 하지만 $J$가 비선형이거나 $x$의 차원이 커지면, 대수적으로 풀리는 경우가 오히려 드물다. 이때 쓰는 게 반복법(iterative method)이다. 초기 추정값 $x_0$에서 출발해서, 매 단계 $J(x)$가 줄어드는 방향으로 $x$를 갱신해가며 최솟값에 접근하는 방식이다. 한 번의 계산으로 정확한 해를 내놓는 직접법(direct method)과 대비되는 개념이다.

이 시리즈는 이 반복법들 — gradient descent, Newton, Gauss-Newton, 그리고 계측 분야에서 사실상 표준으로 쓰는 Levenberg-Marquardt(LM) — 을 순서대로 다룬다. 이번 편은 가장 단순한 gradient descent부터 시작한다.

## 2. 간접측정과 반복법이 필요한 경우

수치적으로 최적화를 풀어야 하는 이유는 결국 하나로 모인다. 모델이 비선형이면, 측정한 결과값으로부터 그 원인이 되는 파라미터를 직접 구하는 역함수가 대개 존재하지 않는다는 것이다. 이건 광학에 국한된 이야기가 아니다. 물리학이든 통계학이든, 어떤 현상을 설명하려는 모델을 세우다 보면 비선형성을 가진 경우가 흔하고, 그럴 때마다 같은 문제에 부딪힌다. 모델이 선형이면 행렬을 한 번 뒤집어 바로 풀리지만, 비선형 모델은 그렇게 안 된다.

내가 오래 다뤄온 예시로는 박막(thin film) 반사율 스펙트럼에서 두께를 뽑아내는 문제가 있다. 처음 이 문제를 접했을 때는 "모델식이 있으니까 역함수를 구하면 되는 거 아닌가" 하고 생각했다. 아주 틀린 생각은 아니다 — 박막의 소광계수(extinction coefficient) $k=0$인 흡수 없는 투명한 박막이라는 특수한 경우에는, 반사율 스펙트럼의 간섭 무늬(interference fringe) 극값들로부터 두께를 닫힌 형태로 바로 구하는 방법(envelope method 계열)이 실제로 있다. 하지만 이건 특수한 경우다. 소광계수가 0이 아니면 굴절률부터 복소수($n=n-ik$)가 되고, 무엇보다 구하고 싶은 두께 $d$가 측정치인 반사율 $R$에 대해 비선형 함수의 변수로 들어가 있어서, $R$을 재고 그로부터 $d$를 직접 얻는 닫힌 형태의 역함수는 일반적으로 존재하지 않는다.

이 문제는 간접측정(indirect measurement)의 한 예다. 간접측정이란 알고 싶은 양을 직접 잴 수 없을 때, 측정 가능한 다른 물리량과 그 파라미터 사이의 관계(모델)를 세우고, 모델을 거꾸로 풀어 파라미터를 추정하는 방식이다. 두께를 자로 재는 대신 반사율 스펙트럼 $R(\lambda)$를 측정하고, 이론 모델 $R_{model}(d,\lambda)$을 거꾸로 풀어 두께 $d$를 구하는 식이다. 모델이 선형이고 닫힌 형태의 역함수가 있으면 대수적으로 바로 풀리지만(직접법), 그렇지 않은 경우가 훨씬 많다. 그럴 때 1절에서 정의한 반복법이 필요해진다.

구체적인 예로, 공기(air, ambient layer — 빛이 입사하는 매질) 위에 SiO2(관심층, layer of interest — 두께를 구하고 싶은 층) 하나가 있고, 그 아래를 Si 기판(substrate layer)이 받치는 구조를 쓴다. 층이 하나뿐인 이런 구조를 단층 샘플(single-layer sample)이라고 부른다. 단순한 구조인데도 반사율은 두께에 대해 지수함수(간섭항)가 들어간 비선형식이다. 수직입사 조건에서는 다음과 같은 닫힌 형태(Airy 공식)로 쓸 수 있다.

$$ r = \frac{r_{01} + r_{12} e^{-2i\beta}}{1 + r_{01} r_{12} e^{-2i\beta}}, \qquad \beta = \frac{2\pi n_1 d}{\lambda} $$

$r_{01}, r_{12}$는 각 계면(ambient/관심층, 관심층/substrate)의 Fresnel 반사계수, $R = \lvert r \rvert^2$이 실제로 측정하는 반사율이다(구현은 `_code/optimization-gradient-descent/reflectance_model.py`). 이 식을 $d$에 대해 닫힌 형태로 뒤집을 수가 없다. 게다가 측정치엔 노이즈가 껴 있어서 "정확히 맞는" $d$란 애초에 존재하지 않고, "가장 잘 맞는" $d$를 찾는 문제가 된다.

그래서 필요한 게 최소자승(least squares) 목적함수다.

$$ J(d) = \frac{1}{2}\sum_i \left( R_{model}(d, \lambda_i) - R_{meas}(\lambda_i) \right)^2 $$

왜 하필 잔차(residual)를 제곱해서 더하는가. 절댓값 합이나 다른 손실도 쓸 수 있지만, 측정 노이즈가 가우시안 분포를 따른다고 가정하면 최소자승해가 최대우도추정(maximum likelihood estimation)과 정확히 일치한다는 통계적 근거가 있다. 계측 노이즈는 실제로 여러 독립적 잡음원의 합이라 중심극한정리에 의해 가우시안에 가까운 경우가 많고, 그래서 최소자승이 이 시리즈 전체의 출발점이 된다.

모델과 측정치의 차이를 제곱해서 더한 값을 최소화하는 $d$를 찾는 게 목표다. 실제로 얼마나 안 맞는지 그림으로 보면 이렇다.

![모델 초기값 vs 측정치](/assets/img/posts/optimization-gradient-descent/fig1-model-vs-measurement.png){: width="600" }
_그림1. 두께 1540nm로 가정한 초기 모델과 실제 1490nm 샘플의 측정치_

측정치는 실제 스펙트로미터 raw 데이터가 아니라, 과거 실험실에서 실측했던 SiO2/Si 두께 계열(10~190nm, 1490nm) 중 1490nm 샘플을 가정하고 같은 물리 모델에 가우시안 노이즈를 더해 합성한 것이다. 1490nm는 가시광 대역 안에 간섭 무늬(fringe)가 여러 개 들어갈 만큼 광학적으로 두꺼워서, 두께 변화가 반사율 스펙트럼에 뚜렷하게 반영된다(그림 1에서 초기 모델과 측정치의 fringe 위치가 어긋나 있는 게 바로 보인다). 초기 모델(파란 선, $d=1540$nm)과 측정치(회색 점)가 얼마나 안 맞는지가 그림 1에서 바로 보인다.

이 안 맞는 정도를 두께 $d$의 함수로 그리면 목적함수 $J(d)$의 모양이 나온다.

![두께에 따른 목적함수](/assets/img/posts/optimization-gradient-descent/fig2-objective-landscape.png){: width="600" }
_그림2. 두께에 따른 목적함수 J(d) — 회색 점선이 전역 최솟값(global minimum)_

목표는 이 $J(d)$ 그래프를 전 구간에 대해 다 그려보고 눈으로 최솟값을 찾는 게 아니다(그림 2는 지금 원리를 보여주려고 미리 넓은 구간을 훑어본 것일 뿐이다). 실제로는 $d$ 하나를 계산할 때마다 전체 파장 대역에 대해 반사율을 계산해야 하므로, 이렇게 촘촘하게 격자 탐색(grid search)을 하는 건 계산량이 크다. 목표는 이 local minimum의 위치, 즉 $J(d)$가 최소가 되는 $d$를 최소한의 연산으로 찾아내는 것이다.

## 3. 목적함수와 기울기 조건

2절에서 반복법이 왜 필요한지는 확인했다. 그런데 반복법이라는 게 매 단계 "지금 추정값 $d_k$에서 다음 추정값 $d_{k+1}$을 어느 방향으로, 얼마나 옮길 것인가"를 정해야 굴러간다. 아무 방향으로나 옮기면 $J$가 오히려 커질 수도 있으니, 이 질문에 답을 줄 무언가가 필요하다. 그 힌트가 기울기(gradient)다.

$J(d)$를 $d$로 미분한 $dJ/dd$는 부호로 "어느 방향"을, 크기로 "지금 얼마나 급격하게 안 맞는가"를 알려준다. 그림 1의 초기값($d=1540$nm)에서 $J$가 크다는 사실 자체가 "여기는 최적점이 아니다"라는 신호이고, 그 지점에서의 기울기가 어느 쪽으로 움직여야 $J$가 줄어드는지를 알려준다. 이번 편은 파라미터가 두께 하나뿐이라 기울기가 스칼라 미분값이지만, 뒤 편에서 파라미터가 여러 개(두께 여러 층, 굴절률 등)로 늘어나면 이 스칼라는 벡터(gradient vector)로, 다음 단계인 2차 미분은 행렬(Hessian)로 확장된다.

이 기울기를 따라 계속 움직이다 보면 결국 어딘가에 도달한다. 더 이상 움직일 방향이 없는, 즉 기울기가 0인 지점이다. 그래서 $J(d)$가 최소가 되는 지점에서는 기울기가 0이어야 한다는 필요조건이 성립한다.

$$ \frac{dJ}{dd} = 0 $$

이건 필요조건이지 충분조건이 아니다. $dJ/dd=0$인 지점이 전역 최솟값(global minimum)일 수도, 국소 최솟값(local minimum)일 수도 있다 — 이 문제는 5절에서 다시 짚고 다음 편(Newton, Gauss-Newton)에서 자세히 다룬다.

## 4. Gradient Descent

가장 단순한 접근은 현재 위치에서 기울기의 반대 방향으로 조금씩 이동하는 것이다. 왜 반대 방향인지는 1차 테일러 전개로 바로 나온다.

$$ J(d + \Delta d) \approx J(d) + \frac{dJ}{dd}\Delta d $$

$\Delta d$를 기울기와 반대 부호로 잡으면($\Delta d = -\alpha \, dJ/dd$, $\alpha>0$) 우변의 둘째 항이 항상 음수가 되어 $J$가 (적어도 국소적으로는) 줄어드는 게 보장된다. 이걸 반복하는 게 gradient descent다.

$$ d_{k+1} = d_k - \alpha \, \frac{dJ}{dd}(d_k) $$

$\alpha$는 스텝 사이즈(step size, learning rate)다. 미분을 해석적으로 구하기 번거로워서 중심차분(central difference)으로 수치미분했다.

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

처음엔 $\alpha$를 별생각 없이 크게 잡았다. `alpha=1500`, 초기값 `d0=1540`nm으로 20번 반복했다. 이게 얼마나 나쁜 선택이었는지 보려면 잘 골랐을 때와 나란히 비교하는 게 낫겠다 싶어서, 똑같은 초기값·반복 횟수로 `alpha=300`도 같이 돌려서 겹쳐 그렸다.

![두께 추정값 vs iteration - step size 비교](/assets/img/posts/optimization-gradient-descent/fig3-thickness-vs-iteration.png){: width="600" }
_그림3. 두께 추정값 vs iteration - step size 비교_

![목적함수 J vs iteration - step size 비교](/assets/img/posts/optimization-gradient-descent/fig4-objective-vs-iteration.png){: width="600" }
_그림4. 목적함수 J vs iteration - step size 비교_

파란 선(`alpha=300`)은 5번째 반복 만에 실제 값(1490nm, 오차 0.3nm 이내)에 도달해서 그대로 눌러앉는다. 목적함수 값도 노이즈 바닥(noise floor) 수준인 약 0.0024까지 매끄럽게 떨어진 뒤 평평해진다(그림 4).

빨간 선(`alpha=1500`)은 딴판이다. 1540 → 1423 → 1527 → 1416 → 1503 → ... 로 매 스텝 수십~100nm 이상씩 요동치다가 20번을 다 돌 때까지 한 번도 안정되지 않는다(그림 3). 최솟값 근처를 몇 번 스치듯 지나가기도 하지만(예: 7번째 반복에서 $J\approx0.03$까지 내려감, 그림 4) 다음 스텝에서 바로 다시 튕겨 나간다. 계곡을 미끄러져 내려가는데 보폭이 너무 크면 바닥을 지나쳐 반대편 비탈을 타고 올라가 버리는 것과 같다 — 한 번 넘어가면 반대편 기울기는 더 가파르고, 그다음 스텝은 더 크게 튕겨 나간다. step size가 너무 크면 나타나는 gradient descent의 대표적인 실패 유형이다. 이 함정은 darkpgmr의 최적화 관련 글([darkpgmr.tistory.com/133](https://darkpgmr.tistory.com/133))에서도 지적하는 내용이다.

$\alpha$ 하나 차이로 완전히 다른 결과가 나온 것이다.

## 5. 한계

이번 편에서 쓴 gradient descent는 두 가지 조건에 기대고 있다.

첫째, 초기값이 정답의 basin(전역 최솟값으로 이어지는 영역) 안에 있어야 한다. 그림 2를 다시 보면 $d=1490$ 근처가 유일한 최솟값이 아니다. 간섭 무늬 하나만큼 떨어진 $d \approx 1290$, $d \approx 1690$ 부근에도 각각 국소 최솟값이 있고, basin 경계(대략 $d=1390$과 $d=1590$, 그림 2의 두 봉우리) 바깥으로 초기값을 두면 알고리즘이 엉뚱한 쪽으로 수렴해버린다. 이 local minimum 문제는 이번 편에서는 건드리지 않았다. 다음 편(Newton, Gauss-Newton)에서 다룬다.

둘째, 수렴 속도가 $\alpha$에 크게 의존하는데 매번 시행착오로 $\alpha$를 찾는 건 비효율적이다. 최솟값에 가까워질수록 $\lvert dJ/dd \rvert$가 0에 가까워지므로 스텝이 점점 작아지고, 반대로 처음부터 너무 크게 잡으면 그림 3처럼 발산한다. 다음 편에서 다룰 Newton 방법은 이 스텝 사이즈를 헤시안(Hessian)으로부터 자동으로 정하는데, 그 대신 헤시안 계산과 역행렬이라는 다른 비용이 붙는다.

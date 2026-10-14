---
title: "최적화 방법론 2편 — Newton법은 왜 빠르게 수렴하다가도 엉뚱한 곳으로 튀는가"
date: 2026-09-10 20:00:00 +0900
categories: [계산과 알고리즘, 최적화방법]
page_id: optimization-newton-gauss-newton
tags: [optimization, newton-method, gauss-newton, least-squares, thin-film, python]
description: Newton법과 Gauss-Newton법으로 박막 두께 피팅을 풀며, 왜 Gauss-Newton이 계측 실무의 표준인지 확인한다.
math: true
---

[1편](/posts/optimization-gradient-descent/)에서는 SiO2/Si 단층 박막의 두께를 반사율 스펙트럼으로부터 구하는 문제를 경사 하강법(gradient descent)으로 풀었다. 스텝 크기(step size) `alpha`를 잘 고르면(예: 300) 문제없이 수렴했지만, 다섯 배로 키우면(1500) 계곡을 넘나들며 수렴하지 못하고 진동했다. 이번 편은 그 스텝 크기 딜레마를 스스로 해결한다는 Newton법부터 시작해서, 계측 분야에서 실제로 훨씬 더 많이 쓰이는 Gauss-Newton법까지 다룬다.

결론부터 말하면, "2차 미분까지 쓰는 Newton법이 이론적으로 더 빠르니 무조건 낫다"는 기대는 이번 문제에서 깨진다. 그 이유를 들여다보면 왜 계측 소프트웨어들이 대부분 Newton법이 아니라 Gauss-Newton(또는 그 변형인 LM(Levenberg-Marquardt)법)을 표준으로 쓰는지 알 수 있다.

## 1. Newton법은 스텝 크기를 어떻게 스스로 정하는가

1편 4절에서 경사 하강법은 1차 테일러 전개 $$J(d+\Delta d) \approx J(d) + J'(d)\Delta d$$ 에서 $\Delta d$를 기울기 반대 방향으로 잡는 방법이었다. 문제는 "얼마나" 움직일지, 즉 스텝 크기 $\alpha$를 사람이 정해야 한다는 점이었다. 곡선의 기울기뿐 아니라 휘어진 정도까지 보면 이 값을 스스로 정할 수 있지 않을까? Newton법의 발상이 이것이다.

Newton법은 한 차수 더 나아가 2차 테일러 전개를 쓴다. 그 자리에서 곡선을 포물선으로 근사하고, 그 포물선으로 스텝을 자동으로 정한다.

$$ J(d+\Delta d) \approx J(d) + J'(d)\Delta d + \frac{1}{2}J''(d)\Delta d^2 $$

여기서 $$J'$$, $$J''$$는 $J$를 두께로 한 번, 두 번 미분한 값이고, $\Delta d$는 이번 스텝에서 움직일 거리다. 우변을 $\Delta d$의 함수로 보고 최솟값을 구하면(우변을 $\Delta d$로 미분해 0으로 놓으면) $\Delta d$가 바로 나온다.

$$ J'(d) + J''(d)\Delta d = 0 \quad\Longrightarrow\quad \Delta d = -\frac{J'(d)}{J''(d)} $$

$$ d_{k+1} = d_k - \frac{J'(d_k)}{J''(d_k)} $$

즉, Newton 스텝은 그 자리에 맞춘 포물선의 꼭짓점으로 곧장 뛰는 것이다. 경사 하강법과 비교하면 $\alpha$ 자리에 $$1/J''(d_k)$$가 들어간 셈이다. 국소적인 곡률(curvature)이 클수록(계곡이 좁고 가파를수록) 작은 스텝을, 곡률이 작을수록(완만할수록) 큰 스텝을 자동으로 잡는다. 그림1은 두께 $d_k=1500$nm 지점에서 실제 $J(d)$를 국소적으로 근사하는 포물선을 보여준다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig1-quadratic-approximation.png" alt="목적함수의 국소 2차 근사와 Newton 스텝" width="600">
_그림1. d=1500nm에서 맞춘 목적함수의 2차 근사(빨간 점선)와 Newton 스텝. 포물선이 실제 J(d)(파란선)와 거의 겹치고, 한 스텝 만에 최솟값 근방인 1489.8nm에 닿는다._

이 지점에서는 포물선이 실제 곡선을 거의 겹칠 정도로 잘 근사하고, 한 번의 스텝($d=1500 \to d\approx1489.8$)만으로 최솟값 근방까지 도달한다. $$J'(d)$$, $$J''(d)$$ 모두 1편처럼 중심차분(central difference)으로 수치미분했다.

```python
# newton.py 핵심부 (전체 코드: _code/optimization-newton-gauss-newton/)
def numerical_grad_hess(thickness_nm, wavelength_nm, measured_R, h=1e-3):
    J0 = objective(thickness_nm, wavelength_nm, measured_R)
    J_plus = objective(thickness_nm + h, wavelength_nm, measured_R)
    J_minus = objective(thickness_nm - h, wavelength_nm, measured_R)
    grad = (J_plus - J_minus) / (2 * h)
    hess = (J_plus - 2 * J0 + J_minus) / (h ** 2)
    return grad, hess


def newton(d0, wavelength_nm, measured_R, n_iter=100, h=1e-3):
    d = d0
    d_hist, J_hist = [d], [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        grad, hess = numerical_grad_hess(d, wavelength_nm, measured_R, h)
        d = d - grad / hess
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)
```

이론적으로 Newton법은 해 근방에서 2차 수렴(quadratic convergence)한다. 오차가 매 스텝 직전 오차의 제곱에 비례해 줄어든다는 뜻이다. 이 문제에서도 그렇다. $d_0=1500$nm에서 시작하면 최솟값(1490.12nm)까지의 오차가 9.88 → 0.288 → 0.0001nm로 줄고, 그다음에는 수치미분의 정밀도 한계에 닿는다. 즉, 해에 가까워질수록 한 스텝에 줄어드는 비율이 커진다. 반면 경사 하강법은 선형 수렴(linear convergence)에 그친다. 그렇다면 1편의 초기값 $d_0=1540$nm에 그대로 Newton법을 적용하면 경사 하강법보다 훨씬 빨리 수렴해야 맞다.

## 2. 같은 초기값에서 Newton법을 돌려 보면

1편과 완전히 같은 조건(seed=0, 실제 두께 1490nm, 노이즈, $d_0=1540$nm)에서 경사 하강법($\alpha=300$), Newton법, 그리고 5절에서 다룰 Gauss-Newton법을 함께 돌렸다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig2-thickness-vs-iteration.png" alt="반복에 따른 두께 추정값" width="600">
_그림2. 초기값 1540nm에서 세 방법의 반복(iteration)에 따른 두께 추정값. Newton법(빨간선)만 최솟값에서 수백 nm 떨어진 곳을 맴돌며 수렴하지 않는다._

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig3-objective-vs-iteration.png" alt="반복에 따른 목적함수 J" width="600">
_그림3. 같은 실행의 목적함수 J(로그 눈금). 경사 하강법은 5스텝, Gauss-Newton법은 3스텝 만에 노이즈 바닥에 닿고, Newton법은 출발점의 값 근처에 머문다._

예상과 다른 결과가 나왔다. 경사 하강법(파란선)은 1편에서 본 대로 매끄럽게 수렴한다. 그런데 Newton법(빨간선)은 첫 스텝부터 1876nm까지 튀어 오르더니, 1~9스텝 동안 1824~1921nm 사이를 오르내리고 10스텝째에는 그 범위를 벗어나 1637nm로 내려온다. 어느 쪽이든 최솟값(1490nm)에서 100nm 넘게 떨어져 있다. 다만 이 개별 값들은 반올림 오차에 민감해 인용할 수치가 아니다. 그 이유는 3절 끝에서 본다.

목적함수로 보면(그림3) 차이가 더 분명하다. 경사 하강법은 5스텝 만에 노이즈 바닥($J\approx0.0024$)에 닿는다. 측정 노이즈 때문에 어떤 두께로도 그 아래로는 내려갈 수 없는 값이다. 반면 Newton법의 $J$는 10스텝 내내 2.4~3.1 사이에 머물러, 출발점의 값(2.57) 근처를 벗어나지 못한다. 초록선(Gauss-Newton법)은 6절에서 본다. "2차 수렴이 더 빠르다"는 이론이 무색하게, 이 초기값에서는 Newton법이 아예 작동하지 않는다.

## 3. Newton법은 왜 실패했는가

1편 2절에서 그린 목적함수 $J(d)$ 그래프(1편 그림2)를 떠올려 보면 이유가 보인다. 2차 근사가 무너지는 지점이 있기 때문이다. 박막의 반사율은 두께가 바뀌면 간섭 무늬(fringe)를 따라 오르내리는 진동함수이고, $J(d)$도 하나의 전역 최솟값 주변에 여러 개의 국소 최솟값이 늘어선 비볼록(non-convex) 함수였다. 2차 근사 자체는 어느 지점에서나 만들 수 있다. 다만 그 포물선이 바닥을 가지려면 그 지점의 곡률, 즉 $$J''(d)$$가 양(+)이어야 한다(포물선이 아래로 볼록해야 스텝이 최솟값 쪽으로 향한다). 실제로 두 지점에서 $$J''(d)$$를 계산해보면 다음과 같다.

| 지점 | $$J''(d)$$ | 의미 |
|---|---|---|
| $d=1500$ | $+0.002487$ | 아래로 볼록, 정상적인 Newton 스텝 |
| $d=1540$ | $-0.000233$ | 위로 볼록(오목), 스텝이 반대 방향으로 튐 |

$d=1540$에서는 곡률의 부호 자체가 뒤집혀 있다. 곡률의 부호는 약 1537nm에서 바뀌고, 간섭 무늬가 만드는 국소적인 "언덕"의 꼭대기는 1589nm에 있다. 1540nm는 그 언덕을 막 오르기 시작한 사면인 셈이다. Newton 스텝 공식 $$-J'/J''$$은 $$J''$$가 음수일 때 오히려 목적함수가 커지는 방향으로 움직인다. 그림4는 같은 방법(Newton법)이 시작점만 40nm 다를 뿐인데도 완전히 다른 결과를 맞는 모습을 보여준다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig4-newton-success-vs-failure.png" alt="Newton법의 두 얼굴" width="600">
_그림4. Newton법의 두 얼굴. 초기값 1500nm(초록선)에서는 두 스텝 만에 노이즈 바닥에 닿고, 1540nm(빨간선)에서는 10스텝이 지나도 출발점의 목적함수 근처에 머문다._

$d_0=1500$(초록선)은 곡률이 정상이라 2회 만에 노이즈 바닥까지 떨어지는 교과서적인 2차 수렴을 보여준다. $d_0=1540$(빨간선)은 10번을 돌려도 초기값보다 먼 곳에 머물러 있다.

실패한 첫 스텝을 그림1과 같은 방식으로 그려 보면 무슨 일이 일어났는지 보인다(그림5). $d=1540$에서 맞춘 포물선은 위로 볼록하다. 바닥이 없으니, 기울기가 0인 점을 찾는 Newton 스텝은 포물선의 꼭대기로 간다. 그 꼭대기는 336nm 떨어진 1876nm에 있고, 첫 스텝이 거기로 튄다. 포물선은 그 자리의 $J$를 15.7로 예측하지만 실제 $J$는 2.83이다. 출발점(2.57)보다 오히려 크다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig5-concave-quadratic-approximation.png" alt="곡률이 음수인 지점의 2차 근사" width="600">
_그림5. 곡률이 음수인 지점(d=1540nm)에서 맞춘 2차 근사. 포물선(빨간 점선)이 위로 볼록해서 Newton 스텝은 그 꼭대기인 1876nm로 튀고, 그 근처에서 포물선은 실제 J(d)(파란선)와 전혀 맞지 않는다._

즉, Newton법은 포물선의 꼭짓점이 바닥인지 꼭대기인지 가리지 않는다. 이는 구현 버그가 아니라 Newton법 자체의 구조적인 약점이다. Newton법의 스텝이 목적함수를 줄이는 방향(descent direction)이라는 보장은 그 지점의 헤시안(Hessian)이 양의 정부호(positive definite)일 때만 성립한다. 이 문제처럼 목적함수가 여러 굴곡을 가진 비볼록 함수라면, 초기값이 어느 "골짜기"에 있느냐에 따라 Newton법의 운명이 갈린다.

한 가지 덧붙일 것이 있다. 그림2·그림3의 발산 궤적에서 개별 반복값은 인용할 만한 수치가 아니다. 첫 스텝 뒤 궤적이 $$J''$$가 0에 가까운 구간을 지나면서, 중심차분으로 구한 $$J''$$의 반올림 오차(약 $10^{-9}$)를 크게 증폭시키기 때문이다. 초기값을 $10^{-12}$nm만 흔들어도 이 반올림 오차의 양상이 바뀌어 10스텝째 두께가 200nm 넘게 달라진다.

```text
발산 궤적의 초기값 민감도 (10스텝째 두께):
  d0=1540.000000000000 -> 1636.55 nm
  d0=1540.000000000001 -> 1855.29 nm
  d0=1540.000000000100 -> 1854.98 nm
```

재현되는 것은 "수렴하지 않고 최솟값에서 수백 nm 떨어진 곳을 헤맨다"는 사실이지, 특정 반복값이 아니다.

## 4. Newton법의 한계

이번 편에서는 파라미터가 두께 하나뿐이라 $$J''(d)$$가 스칼라였다. 파라미터가 여러 개가 되면 무엇이 더 어려워지는가? 실제 계측에서는 두께·굴절률·소광계수 등 여러 파라미터를 동시에 피팅하는 경우가 많다. 파라미터가 $n$개면 $$J''$$는 $n \times n$ 헤시안 행렬이 되고, 매 스텝 이 행렬을 계산하고 역행렬을 구해야 한다. 대칭 행렬이라도 서로 다른 2차 도함수가 $n(n+1)/2$개다. 두께 둘과 굴절률 분산 계수 셋, 곧 파라미터가 5개면 매 스텝 15개를 수치미분해야 한다. 5절의 Gauss-Newton법은 잔차의 1차 도함수만 있으면 된다.

비용만 문제가 아니다. 3절에서 본 음의 곡률 위험도 그대로 남는다. 그래서 Newton법은 실무에서 그대로 쓰지 않고 안전장치를 반드시 곁들인다. 선 탐색(line search)은 방향은 Newton 스텝대로 두고 스텝 크기를 줄여 가며 $J$가 실제로 줄었는지 확인한다. 신뢰영역(trust region)은 2차 모델을 믿을 수 있는 반경을 정해 그 안에서만 움직인다. 다음 편의 LM법도 이 신뢰영역 계열이다.

## 5. Gauss-Newton법은 헤시안을 어떻게 근사하는가

헤시안을 통째로 계산하지 않으면서 음의 곡률도 피할 길은 없는가? 목적함수가 제곱합이라는 사실을 쓰면 된다. 비선형 최소자승(nonlinear least squares) 문제는 구조가 조금 특별하다. 목적함수가 잔차(residual) $r_i(d) = R_{model}(d,\lambda_i) - R_{meas}(\lambda_i)$ 들의 제곱합이라는 형태로 고정되어 있다.

$$ J(d) = \frac{1}{2}\sum_i r_i(d)^2 $$

이 형태를 그대로 미분하면 1차·2차 도함수가 다음과 같이 나온다.

$$ J'(d) = \sum_i r_i(d)\, r_i'(d), \qquad J''(d) = \sum_i \Big[ r_i'(d)^2 + r_i(d)\, r_i''(d) \Big] $$

여기서 $$r_i'(d)$$는 파장 $\lambda_i$의 반사율이 두께에 얼마나 민감한가, 곧 감도이고, $$r_i''(d)$$는 그 감도가 두께에 따라 변하는 정도다. Newton법은 $$J''(d)$$ 전체를 쓴다. Gauss-Newton법은 여기서 둘째 항 $$\sum_i r_i(d) r_i''(d)$$ 을 통째로 버리고 첫째 항만 근사 헤시안으로 쓴다.

$$ J''(d) \;\approx\; \sum_i r_i'(d)^2 $$

(파라미터가 여러 개인 일반적인 표기로는 이 근사 헤시안을 $J^TJ$라고 쓴다. $J$는 잔차 벡터의 야코비안(Jacobian)이다. 목적함수 $J(d)$와 글자가 같지만 관례를 따랐다. 이 글에서 $J^TJ$, $J^T\mathbf{r}$의 $J$는 야코비안이고, $J(d)$, $$J'$$, $$J''$$의 $J$는 목적함수다.) 버려지는 둘째 항은 "잔차 $r_i$"와 "잔차의 곡률 $$r_i''$$"의 곱이다. 모델이 데이터에 잘 맞아서 잔차가 작다면(또는 모델이 국소적으로 거의 직선이라 $$r_i''$$가 작다면) 이 항은 원래도 작아서 버려도 큰 차이가 없다.

같은 근사를 다른 쪽에서 볼 수도 있다. 매 스텝 잔차를 그 자리에서 직선으로 펴는 것이다.

$$ r_i(d+\Delta d) \approx r_i(d) + r_i'(d)\,\Delta d $$

이 직선들의 제곱합 $$\frac{1}{2}\sum_i \big[ r_i(d) + r_i'(d)\,\Delta d \big]^2$$을 가장 작게 하는 $\Delta d$를 고르면 아래의 스텝 공식이 그대로 나온다. 즉, Gauss-Newton법은 매 스텝 직선 맞춤(선형 최소자승)을 한 번 푸는 방법이다. 굽은 산길을 그 자리의 평평한 지도로 바꿔 갈 곳을 정하고, 한 걸음 간 뒤 지도를 다시 그리는 것과 같다. 직선 모델의 제곱합은 $\Delta d$의 2차식이고 $\Delta d^2$의 계수가 $$\frac{1}{2}\sum_i r_i'^2 \ge 0$$이므로, 언제나 아래로 볼록한 그릇이다. 뒤에서 볼 근사 헤시안의 좋은 성질이 여기서 나온다. 버린 항 $$r_i r_i''$$도 같은 그림으로 읽힌다. 잔차가 직선에서 벗어난 정도($$r_i''$$)에 잔차의 크기($r_i$)를 곱한 것이다.

이 근사를 쓰면 스텝 공식은 다음과 같다.

$$ d_{k+1} = d_k - \frac{\sum_i r_i'(d_k)\, r_i(d_k)}{\sum_i r_i'(d_k)^2} $$

파라미터가 여러 개인 일반적인 경우엔 이 식이 정규방정식(normal equation) $(J^TJ)\,\Delta\mathbf{d} = -J^T\mathbf{r}$ 형태가 된다. 코드로는 다음과 같다.

```python
# gauss_newton.py 핵심부 (전체 코드: _code/optimization-newton-gauss-newton/)
def gauss_newton(d0, wavelength_nm, measured_R, n_iter=100, h=1e-3):
    d = d0
    d_hist, J_hist = [d], [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        r = residual(d, wavelength_nm, measured_R)
        Jr = numerical_jacobian(d, wavelength_nm, measured_R, h)
        step = -np.sum(Jr * r) / np.sum(Jr ** 2)  # (J^T J) Δd = -J^T r
        d = d + step
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)
```

이 근사 헤시안 $J^TJ$(이번 문제에서는 $$\sum_i r_i'^2$$)에는 Newton법의 실제 헤시안에 없는 중요한 성질이 있다. 제곱의 합이므로 $J^TJ$는 언제나 양의 준정부호(positive semidefinite)이고, 야코비안이 최대 계수(full rank)이면 양의 정부호다. 이번 문제로 말하면 $$\sum_i r_i'^2$$는 감도의 제곱합이므로, 두께에 반응하는 파장이 하나라도 있으면 양수다. 즉 Gauss-Newton법의 스텝에는 곡률의 부호가 뒤집히는 지점이 없어서, 3절에서 Newton법이 발산했던 바로 그 이유(음의 곡률)로부터 원천적으로 자유롭다.

## 6. 같은 실패 지점에서 Gauss-Newton법을 돌려 보면

그림2·그림3의 초록선이 Gauss-Newton법이다. $d_0=1540$nm, Newton법이 발산했던 바로 그 시작점에서 Gauss-Newton법은 3번째 스텝 만에 노이즈 바닥까지 떨어졌다. 경사 하강법($\alpha=300$)이 이 조건에서 5스텝 걸린 것보다도 빠르다.

| 방법 | $d_0=1540$nm 결과 |
|---|---|
| 경사 하강법 ($\alpha=300$) | 1490.12nm, 5스텝 수렴 |
| Newton법 | 10스텝까지 미수렴, 최솟값에서 100nm 이상 떨어진 곳 |
| Gauss-Newton법 | 1490.12nm, 3스텝 수렴 |

수렴의 모양도 볼 만하다. 최솟값까지의 오차는 49.9 → 20.3 → 1.77 → 0.0036nm로 줄고, 그 뒤로는 매 스텝 약 $7.3\times10^{-4}$배씩 준다. 비율이 일정하니 엄밀히는 1차 수렴이다. 이 비율은 해에서 버린 항과 남긴 항의 비 $$\lvert \sum_i r_i r_i'' \rvert / \sum_i r_i'^2$$와 같은 값이다. 잔차가 노이즈뿐이라 버린 항이 작고, 그래서 1차 수렴인데도 2차 수렴과 구별하기 어려울 만큼 빠르다.

헤시안의 실제 곡률 부호가 무엇이든 상관없이 $$\sum_i r_i'^2$$은 항상 0 이상이라는 사실 하나가 이런 차이를 만든다. 계측 소프트웨어들이 일반 Newton법 대신 Gauss-Newton법(혹은 다음 편에서 다룰 LM법)을 표준으로 쓰는 이유가 바로 여기 있다. 2차 미분 계산 비용을 아끼는 것도 있지만, 그보다 이 구조적인 안정성이 더 크다.

## 7. Gauss-Newton도 만능은 아니다

그렇다고 Gauss-Newton법이 모든 상황에서 안전한 것은 아니다. 버려진 둘째 항 $$\sum_i r_i r_i''$$이 항상 무시할 만큼 작은 것은 아니기 때문이다. 이론적으로 알려진 실패 조건은 크게 세 가지다.

첫째, 잔차 자체가 크면(모델이 데이터에 잘 안 맞으면) 버려진 항이 커져서 2차 수렴은커녕 선형 수렴에 그치거나 발산할 수 있다. 해에서 잔차가 0인 문제(영잔차 문제, zero-residual problem)에서는 Gauss-Newton법도 Newton법과 같은 2차 수렴을 회복한다. 잔차가 큰지는 피팅 뒤 잔차의 RMS를 측정 노이즈와 견주어 보면 알 수 있다. 이번 문제의 최솟값에서 잔차 RMS는 0.0040으로 노이즈 표준편차(0.004)와 같다. 이렇게 잔차가 노이즈 수준이면 영잔차 문제에 가깝다고 보아도 무난하다.

둘째, 야코비안의 계수(rank)가 모자라거나 특이에 가까우면(이번 문제에서는 $$\sum_i r_i'^2$$이 0에 가까우면) 근사 헤시안 $J^TJ$가 거의 0이 된다. 그러면 버려진 항을 상대적으로 무시할 수 없게 되고, 정규방정식의 해, 곧 스텝이 지나치게 커진다. [3편](/posts/optimization-levenberg-marquardt/) 6절에서 굴절률을 1.02로 잘못 둔 모델로 이 상황을 재현한다.

셋째, 시작점이 다른 국소 최솟값의 분지(basin)에 있으면 Gauss-Newton법도 당연히 그 지점으로 수렴한다. 실제로 $d_0$를 1300nm 근방, 1690nm 근방으로 바꿔서 돌려보면 각각 1293nm, 1689nm로, 1편 그림2에서 본 이웃 국소 최솟값으로 정확히 수렴한다. Gauss-Newton법의 안정성은 "전역 최솟값을 찾아준다"는 뜻이 아니라 "한 번 방향을 잡으면 헤매지 않고 그 방향으로 곧장 간다"는 뜻에 가깝다.

## 정리

Newton법은 이론상 2차 수렴을 약속한다. 하지만 그 약속은 목적함수가 국소적으로 볼록(convex)할 때만 유효하다. 간섭 무늬로 굴곡진 이 문제에서는 초기값에 따라 아예 수렴하지 않는다. Gauss-Newton법은 비선형 최소자승이라는 문제의 구조를 이용해 헤시안을 $J^TJ$로 근사함으로써 이 실패 조건 자체를 없앤다. 대신 잔차가 큰 경우의 수렴 속도 저하라는 다른 대가를 치른다. 즉, 최소자승 문제라면 Newton법을 그대로 쓸 이유가 거의 없다. 목적함수가 제곱합이고 잔차가 작으면 Gauss-Newton법을, 스텝이 튀거나 감도가 약하면 다음 편의 LM법을 쓴다. 어느 쪽도 어느 분지로 갈지는 고르지 못한다. 그 문제는 [4편](/posts/optimization-global-heuristics/)에서 다룬다.

다음 편에서는 Gauss-Newton법의 스텝이 지나치게 커지거나 야코비안이 특이에 가까워 불안정해지는 경우를 댐핑 계수(damping parameter)로 보완하는 LM법을 다룬다. 계측 분야에서 사실상 표준으로 쓰이는 방법으로, 경사 하강법과 Gauss-Newton법 사이를 매끄럽게 오가는 구조를 가지고 있다. 7절의 첫째 조건에 비추어 말하면, LM법이 막아 주는 것은 잔차가 클 때 생기는 발산이고, 남는 것은 느린 수렴이다. 해 근처에서는 Gauss-Newton법과 같은 근사 헤시안을 쓰기 때문이다.

## 참고자료

- D. Bindel, "Gauss-Newton and nonlinear least squares" (lecture notes, 2023-04-10), CS 4220 Numerical Analysis, Cornell University, Spring 2023. [cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf](https://www.cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf)
- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk 원문](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

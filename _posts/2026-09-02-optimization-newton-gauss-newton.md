---
layout: post
title: "최적화 방법론 2편 — Newton법과 Gauss-Newton법"
date: 2026-09-02 09:00:00 +0900
categories: [최적화 방법]
tags: [optimization, newton-method, gauss-newton, least-squares, thin-film, python]
description: "Newton법과 Gauss-Newton법으로 박막 두께 피팅을 풀며, 왜 Gauss-Newton이 계측 실무의 표준인지 확인한다."
math: true
---

[1편](/posts/optimization-gradient-descent/)에서는 SiO2/Si 단층 박막의 두께를 반사율 스펙트럼으로부터 구하는 문제를 gradient descent로 풀었다. 스텝 사이즈 `alpha`를 잘 고르면(예: 300) 문제없이 수렴했지만, 조금만 크게 잡아도(예: 1500) 계곡을 넘나들며 발산했다. 이번 편은 그 스텝 사이즈 딜레마를 스스로 해결한다는 Newton법부터 시작해서, 계측 분야에서 실제로 훨씬 더 많이 쓰이는 Gauss-Newton법까지 다룬다.

결론부터 말하면, "2차 미분까지 쓰는 Newton법이 이론적으로 더 빠르니 무조건 낫다"는 기대는 이번 문제에서 깨진다. 그 이유를 들여다보면 왜 계측 소프트웨어들이 대부분 Newton법이 아니라 Gauss-Newton(또는 그 변형인 Levenberg-Marquardt)을 표준으로 쓰는지 알 수 있다.

## 1. Newton법 — 2차 근사로 스텝을 자동으로 정하기

1편 4절에서 gradient descent는 1차 테일러 전개 $J(d+\Delta d) \approx J(d) + J'(d)\Delta d$ 에서 $\Delta d$를 기울기 반대 방향으로 잡는 방법이었다. 문제는 "얼마나" 움직일지, 즉 스텝 크기 $\alpha$를 사람이 정해야 한다는 점이었다.

Newton법은 한 차수 더 나아가 2차 테일러 전개를 쓴다.

$$ J(d+\Delta d) \approx J(d) + J'(d)\Delta d + \frac{1}{2}J''(d)\Delta d^2 $$

우변을 $\Delta d$에 대한 함수로 보고 최솟값을 구하면(우변을 $\Delta d$로 미분해 0으로 놓으면) $\Delta d$가 바로 나온다.

$$ J'(d) + J''(d)\Delta d = 0 \quad\Longrightarrow\quad \Delta d = -\frac{J'(d)}{J''(d)} $$

$$ d_{k+1} = d_k - \frac{J'(d_k)}{J''(d_k)} $$

gradient descent와 비교하면 $\alpha$ 자리에 $1/J''(d_k)$가 들어간 셈이다. 국소적인 곡률(curvature)이 클수록(계곡이 좁고 가파를수록) 작은 스텝을, 곡률이 작을수록(완만할수록) 큰 스텝을 자동으로 잡는다. 그림1은 두께 $d_k=1500$nm 지점에서 실제 $J(d)$를 국소적으로 근사하는 포물선을 보여준다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig1-quadratic-approximation.png" alt="목적함수의 국소 2차 근사와 Newton 스텝" width="600">
_그림1. 목적함수의 국소 2차 근사와 Newton 스텝_

이 지점에서는 포물선이 실제 곡선을 거의 겹칠 정도로 잘 근사하고, 한 번의 스텝($d=1500 \to d\approx1489.8$)만으로 최솟값 근방까지 도달한다. $J'(d)$, $J''(d)$ 모두 1편처럼 중심차분(central difference)으로 수치미분했다.

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

이론적으로 Newton법은 해 근방에서 2차 수렴(quadratic convergence)한다 — 오차가 매 스텝 제곱으로 줄어든다는 뜻이다. 반면 gradient descent는 선형 수렴(linear convergence)에 그친다. 그렇다면 1편의 초기값 $d_0=1540$nm에 그대로 Newton법을 적용하면 gradient descent보다 훨씬 빨리 수렴해야 맞다.

## 2. 실험 — 같은 초기값에서 돌려보면

1편과 완전히 같은 조건(seed=0, 실제 두께 1490nm, 노이즈, $d_0=1540$nm)에서 gradient descent(alpha=300), Newton, 그리고 5절에서 다룰 Gauss-Newton을 함께 돌렸다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig2-thickness-vs-iteration.png" alt="두께 추정값 vs iteration" width="600">
_그림2. 두께 추정값 vs iteration (d0=1540nm)_

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig3-objective-vs-iteration.png" alt="목적함수 J vs iteration" width="600">
_그림3. 목적함수 J vs iteration (d0=1540nm, log scale)_

예상과 다른 결과가 나왔다. Gradient descent(파란선)는 1편에서 본 대로 매끄럽게 수렴한다. 그런데 Newton(빨간선)은 첫 스텝부터 1875nm까지 튀어 오르더니, 10번을 돌려도 1660~1920nm 사이를 오르내리며 전혀 정착하지 못한다. "2차 수렴이 더 빠르다"는 이론이 무색하게, 이 초기값에서는 Newton법이 아예 작동하지 않는다.

## 3. 왜 실패했는가 — 2차 근사가 무너지는 지점

1편 5절에서 그린 목적함수 $J(d)$ 그래프(1편 그림2)를 떠올려 보면 이유가 보인다. 박막의 반사율은 두께에 대해 간섭 무늬(fringe)를 그리는 진동함수이고, $J(d)$도 하나의 전역 최솟값 주변에 여러 개의 국소 최솟값이 늘어선 비볼록(non-convex) 함수였다. Newton법이 기대는 2차 근사가 성립하려면 그 지점의 곡률, 즉 $J''(d)$가 양(+)이어야 한다(포물선이 아래로 볼록해야 스텝이 최솟값 쪽으로 향한다). 실제로 두 지점에서 $J''(d)$를 계산해보면 다음과 같다.

| 지점 | $J''(d)$ | 의미 |
|---|---|---|
| $d=1500$ | $+0.00249$ | 아래로 볼록 — 정상적인 Newton 스텝 |
| $d=1540$ | $-0.00023$ | 위로 볼록(오목) — 스텝이 반대 방향으로 튐 |

$d=1540$에서는 곡률의 부호 자체가 뒤집혀 있다. 이 지점은 간섭 무늬가 만드는 국소적인 "언덕" 근처였던 셈이고, Newton 스텝 공식 $-J'/J''$은 $J''$가 음수일 때 오히려 목적함수가 커지는 방향으로 움직인다. 그림4는 같은 방법(Newton)이 시작점만 40nm 다를 뿐인데도 완전히 다른 운명을 맞는 모습을 보여준다.

<img src="/assets/img/posts/optimization-newton-gauss-newton/fig4-newton-success-vs-failure.png" alt="Newton법의 두 얼굴" width="600">
_그림4. Newton법의 두 얼굴 - 초기값에 따른 성공/실패_

$d_0=1500$(초록선)은 곡률이 정상이라 2회 만에 노이즈 바닥까지 떨어지는 교과서적인 2차 수렴을 보여준다. $d_0=1540$(빨간선)은 10번을 돌려도 초기값과 비슷한 수준에 머물러 있다. 이건 구현 버그가 아니라 Newton법 자체의 구조적인 약점이다. Newton법의 스텝이 매번 목적함수를 줄이는 방향(descent direction)이라는 보장은 그 지점의 Hessian이 양의 정부호(positive definite)일 때만 성립한다고 알려져 있다. 우리 문제처럼 목적함수가 여러 굴곡을 가진 비볼록 함수라면, 초기값이 어느 "골짜기"에 있느냐에 따라 Newton법의 운명이 갈린다.

## 4. Newton법의 한계

이번 편에서는 파라미터가 두께 하나뿐이라 $J''(d)$가 스칼라였지만, 실제 계측에서는 두께·굴절률·흡수계수 등 여러 파라미터를 동시에 피팅하는 경우가 많다. 파라미터가 $n$개면 $J''$는 $n \times n$ Hessian 행렬이 되고, 매 스텝 이 행렬을 계산하고 역행렬을 구해야 한다. 계산량이 $n$이 커질수록 빠르게 늘어나는 것은 물론이고, 이번 실험처럼 Hessian이 양의 정부호가 아닌 지점을 지나가면 스텝이 최솟값과 무관한 방향으로 튈 위험도 여전히 남는다. Newton법을 실무에서 그대로 쓰지 않고 line search나 신뢰영역(trust region) 같은 안전장치를 반드시 곁들이는 이유다.

## 5. Gauss-Newton법 — Hessian을 근사해서 안전하게

비선형최소자승 문제는 구조가 조금 특별하다. 목적함수가 잔차(residual) $r_i(d) = R_{model}(d,\lambda_i) - R_{meas}(\lambda_i)$ 들의 제곱합이라는 형태로 고정되어 있다.

$$ J(d) = \frac{1}{2}\sum_i r_i(d)^2 $$

이 형태를 그대로 미분하면 1차·2차 도함수가 다음과 같이 나온다.

$$ J'(d) = \sum_i r_i(d)\, r_i'(d), \qquad J''(d) = \sum_i \Big[ r_i'(d)^2 + r_i(d)\, r_i''(d) \Big] $$

Newton법은 $J''(d)$ 전체를 쓴다. Gauss-Newton법은 여기서 둘째 항 $\sum_i r_i(d) r_i''(d)$ 을 통째로 버리고 첫째 항만 근사 Hessian으로 쓴다.

$$ J''(d) \;\approx\; \sum_i r_i'(d)^2 $$

(파라미터가 여러 개인 일반적인 표기로는 이 근사 Hessian을 $J^TJ$라고 쓴다. $J$는 잔차 벡터의 Jacobian이다.) 버려지는 둘째 항은 "잔차 $r_i$"와 "잔차의 곡률 $r_i''$"의 곱이다. 모델이 데이터에 잘 맞아서 잔차가 작다면(또는 모델이 국소적으로 거의 직선이라 $r_i''$가 작다면) 이 항은 원래도 작아서 버려도 큰 차이가 없다. 이 근사를 쓰면 스텝 공식은 다음과 같다.

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
        step = -np.sum(Jr * r) / np.sum(Jr ** 2)  # (J^T J) h = -J^T r
        d = d + step
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)
```

이 근사 Hessian $J^TJ$(우리 문제에서는 $\sum_i r_i'^2$)에는 Newton의 실제 Hessian에 없는 중요한 성질이 있다. 제곱의 합이라 $J^TJ$는 (Jacobian이 rank를 유지하는 한) 항상 양의 정부호이거나 최소한 준정부호(positive semidefinite)다. 즉 Gauss-Newton의 스텝에는 곡률의 부호가 뒤집히는 지점이 없어서, 3절에서 Newton이 발산했던 바로 그 이유(음의 곡률)로부터 원천적으로 자유롭다.

## 6. 실험 — 같은 실패 지점에서 Gauss-Newton은?

그림2·그림3의 초록선이 Gauss-Newton이다. $d_0=1540$nm, Newton이 발산했던 바로 그 시작점에서 Gauss-Newton은 3번째 스텝 만에 노이즈 바닥까지 떨어졌다 — gradient descent(alpha=300)가 이 조건에서 4~5스텝 걸린 것보다도 빠르다.

| 방법 | $d_0=1540$nm 결과 |
|---|---|
| Gradient Descent (alpha=300) | 1490.12nm, 4~5스텝 근방 수렴 |
| Newton | 1666.31nm, 10스텝까지 미수렴(발산) |
| Gauss-Newton | 1490.12nm, 3스텝 수렴 |

Hessian의 실제 곡률 부호가 무엇이든 상관없이 $\sum_i r_i'^2$은 항상 0 이상이라는 사실 하나가 이런 차이를 만든다. 계측 소프트웨어들이 일반 Newton법 대신 Gauss-Newton(혹은 다음 편에서 다룰 Levenberg-Marquardt)을 표준으로 쓰는 이유가 바로 여기 있다. 2차 미분 계산 비용을 아끼는 것도 있지만, 그보다 이 구조적인 안정성이 더 크다.

## 7. Gauss-Newton도 만능은 아니다

그렇다고 Gauss-Newton이 모든 상황에서 안전한 것은 아니다. 버려진 둘째 항 $\sum_i r_i r_i''$이 항상 무시할 만큼 작은 것은 아니기 때문이다. 이론적으로 알려진 실패 조건은 크게 두 가지다.

첫째, 잔차 자체가 크면(모델이 데이터에 잘 안 맞으면) 버려진 항이 커져서 2차 수렴은커녕 선형 수렴에 그치거나 발산할 수 있다. 실제로 해에서 잔차가 정확히 0인 경우(consistent problem)에만 Gauss-Newton도 Newton과 같은 2차 수렴을 회복한다.

둘째, 시작점이 다른 국소 최솟값의 basin에 있으면 Gauss-Newton도 당연히 그 지점으로 수렴한다. 실제로 $d_0$를 1300nm 근방, 1690nm 근방으로 바꿔서 돌려보면 각각 1293nm, 1688nm — 1편 그림2에서 본 이웃 국소 최솟값으로 정확히 수렴한다. Gauss-Newton의 안정성은 "전역 최솟값을 찾아준다"는 뜻이 아니라 "한 번 방향을 잡으면 헤매지 않고 그 방향으로 곧장 간다"는 뜻에 가깝다.

## 정리

이번 편에서 확인한 것을 정리하면, Newton법은 이론상 가장 빠른 2차 수렴을 약속하지만 그 약속은 목적함수가 국소적으로 볼록(convex)할 때만 유효하고, 우리 문제처럼 간섭 무늬로 굴곡진 비볼록 목적함수에서는 초기값에 따라 완전히 실패할 수 있다. Gauss-Newton법은 비선형최소자승이라는 문제의 구조를 이용해 Hessian을 $J^TJ$로 근사함으로써 이 실패 조건 자체를 없앤다 — 대신 잔차가 큰 경우의 수렴 속도 저하라는 다른 대가를 치른다.

다음 편에서는 Gauss-Newton법의 이 약점(잔차가 클 때의 불안정성)을 damping parameter로 보완하는 Levenberg-Marquardt법을 다룬다. 계측 분야에서 사실상 표준으로 쓰이는 방법으로, gradient descent와 Gauss-Newton 사이를 매끄럽게 오가며 두 방법의 장점만 취하는 구조를 가지고 있다.

## 참고자료

- D. Bindel, "Nonlinear Least Squares (Newton and Gauss-Newton)," Numerical Analysis lecture notes, Cornell University, 2023. [cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf](https://www.cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf)
- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk 원문](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

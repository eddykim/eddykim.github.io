---
layout: post
title: "최적화 방법론 3편 — Levenberg-Marquardt법"
date: 2026-09-02 20:00:00 +0900
categories: [계산과 알고리즘, 최적화방법]
tags: [optimization, levenberg-marquardt, gauss-newton, least-squares, thin-film, python]
description: "damping parameter μ로 Gauss-Newton과 steepest descent를 매끄럽게 오가는 Levenberg-Marquardt법을 실제 실패 사례로 검증한다."
math: true
---

[2편](/posts/optimization-newton-gauss-newton/)에서는 같은 두께 피팅 문제에 Newton법과 Gauss-Newton법을 적용했다. Newton법은 초기값 $d_0=1540$nm에서 곡률의 부호가 뒤집혀($J''(d)<0$) 완전히 발산했고, Gauss-Newton법은 Hessian을 $J^TJ$로 근사한 덕분에 그 문제로부터 자유로워 3스텝 만에 수렴했다. 다만 2편 말미에서 Gauss-Newton에도 이론적인 한계가 있다고 짚었다 — 버려진 항 $\sum_i r_i r_i''$이 무시할 수 없을 만큼 크면(잔차가 크거나 Jacobian이 특이에 가까우면) 수렴이 느려지거나 아예 실패할 수 있다는 것이었다.

이번 편은 그 한계를 damping parameter $\mu$로 보완하는 Levenberg-Marquardt(LM)법을 다룬다. 계측 소프트웨어 대부분이 실제로 표준으로 쓰는 방법이며, 실험을 통해 Gauss-Newton이 정말로 실패하는 상황을 직접 만들어 LM이 어떻게 구제하는지 확인한다.

## 1. Levenberg-Marquardt법 — GN과 steepest descent를 잇는 다리

2편에서 Gauss-Newton의 스텝은 정규방정식 $(J^TJ)\,\Delta d = -J^Tr$을 풀어서 구했다(우리 문제에서는 $J^TJ = \sum_i r_i'^2$, $J^Tr=\sum_i r_i' r_i$인 스칼라). Levenberg-Marquardt는 이 정규방정식의 대각에 $\mu$를 더한다.

$$ (J^TJ + \mu I)\, h_{lm} = -J^Tr $$

파라미터가 두께 하나뿐이라 $J^TJ$가 스칼라 $A=\sum_i r_i'^2$이므로, 스텝은 다음과 같이 간단해진다.

$$ h_{lm} = -\frac{g}{A+\mu}, \qquad g = \sum_i r_i' r_i $$

$\mu$의 역할은 극단을 보면 뚜렷하다. $\mu \to 0$이면 $h_{lm}$은 그대로 Gauss-Newton 스텝이 된다. 반대로 $\mu \to \infty$이면 $A$는 무시할 만큼 작아져 $h_{lm} \approx -g/\mu$, 즉 기울기 반대 방향으로 $1/\mu$만큼 움직이는 steepest descent가 된다. $\mu$ 하나가 "안전하지만 느린 방법"과 "빠르지만 곡률에 취약한 방법" 사이를 매끄럽게 오가는 손잡이인 셈이다.

초기 damping은 $\mu_0 = \tau \cdot A$로 정한다. $\tau$는 사용자가 고르는 스케일 인자로, 초기값이 해에 가깝다고 믿으면 작게(예: $10^{-6}$), 자신이 없으면 크게(예: $10^{-3}\sim1$ 자리) 잡으라고 알려져 있다.

```python
# levenberg_marquardt.py 핵심부 (전체 코드: _code/optimization-levenberg-marquardt/)
mu = tau * A   # A = sum(Jr**2), g = sum(Jr*r)
nu = 2.0
for _ in range(n_iter):
    h_lm = -g / (A + mu)
    d_new = d + h_lm
    F_old, F_new = objective(d, ...), objective(d_new, ...)
    L0_minus_Lh = 0.5 * h_lm * (mu * h_lm - g)
    rho = (F_old - F_new) / L0_minus_Lh if L0_minus_Lh > 0 else -1.0

    if rho > 0:                    # 스텝 채택
        d = d_new
        Jr, r = numerical_jacobian(d, ...), residual(d, ...)
        A, g = np.sum(Jr ** 2), np.sum(Jr * r)
        mu = mu * max(1/3, 1 - (2 * rho - 1) ** 3)
        nu = 2.0
    else:                           # 스텝 기각, 더 보수적으로
        mu = mu * nu
        nu = 2.0 * nu
```

## 2. 실험 — $\mu_0$을 크게 잡을 때 vs 작게 잡을 때

1·2편과 같은 조건($d_0=1540$nm, 실제 두께 1490nm)에서 $\tau=10^{-6}$(Gauss-Newton에 가까움)과 $\tau=10^{6}$(steepest descent에 가까움) 두 경우를 돌려봤다. 이 지점에서 $A \approx 0.00264$이므로 각각 $\mu_0 \approx 2.64\times10^{-9}$, $\mu_0 \approx 2644$가 된다.

<img src="/assets/img/posts/optimization-levenberg-marquardt/fig1-mu0-comparison.png" alt="mu0 크기에 따른 LM 수렴 궤적" width="600">
_그림1. μ0 크기에 따른 LM 수렴 궤적 (d0=1540nm)_

$\tau=10^{-6}$(파란선)은 사실상 Gauss-Newton과 같은 궤적으로 3스텝 만에 노이즈 바닥까지 떨어진다. $\tau=10^{6}$(빨간선)은 처음 10스텝 가까이 1540nm 근방에서 거의 움직이지 않다가 — $\mu$가 워낙 커서 스텝이 $-g/\mu$로 짓눌려 있기 때문이다 — 11번째 스텝 무렵부터 갑자기 가속해 20스텝 안에 같은 지점에 도달한다. "느리지만 절대 튀지 않는 시작"에서 "빠른 Gauss-Newton"으로의 전환이 한 그래프 안에 다 보인다.

## 3. gain ratio가 $\mu$를 조정하는 방식

이 전환은 매 스텝 gain ratio $\rho$를 계산해 자동으로 일어난다.

$$ \rho = \frac{F(d) - F(d+h_{lm})}{L(0) - L(h_{lm})}, \qquad L(0)-L(h_{lm}) = \frac{1}{2}h_{lm}(\mu h_{lm} - g) $$

분모는 국소 2차 모델이 예측한 감소량, 분자는 실제 감소량이다. $\rho>0$이면(모델이 잘 맞으면) 스텝을 받아들이고 $\mu \leftarrow \mu \cdot \max\{1/3,\ 1-(2\rho-1)^3\}$로 줄인다. $\rho \le 0$이면 스텝을 기각하고 $\mu \leftarrow \nu\mu$로 키운 뒤 $\nu$를 2배로 늘린다(Marquardt의 원래 규칙은 단순히 곱하기/나누기만 반복해 flutter가 있었는데, Nielsen이 제안한 이 규칙은 그 진동을 줄인다).

<img src="/assets/img/posts/optimization-levenberg-marquardt/fig2-mu-adaptation.png" alt="gain ratio 기반 mu 자동조정 추이" width="600">
_그림2. gain ratio 기반 μ 자동조정 추이_

두 경우 모두 스텝이 성공하는 동안은 $\mu$가 정확히 매번 $1/3$씩 줄어든다(로그 스케일에서 직선). 국소 2차 근사가 실제로 잘 맞아 $\rho$가 1에 가깝기 때문이다. 흥미로운 건 $\tau=10^{-6}$ 쪽(파란선)이 노이즈 바닥에 도달한 뒤의 움직임이다. $d$가 더는 유의미하게 줄지 않는 지점(노이즈가 만든 바닥)에 도달하면 이후 스텝들은 $\rho \le 0$이 되어 기각되고, $\mu$가 $\times2, \times4, \times8, \dots$로 다시 커진다 — 알고리즘이 "더 줄일 게 없다"는 걸 감지하고 스스로 보수적으로 바뀌는 모습이다.

## 4. 네 가지 방법 종합 비교

같은 $d_0=1540$nm, 10스텝 조건에서 1~2편의 세 방법에 Levenberg-Marquardt($\tau=10^{-6}$)를 더해 나란히 돌렸다.

<img src="/assets/img/posts/optimization-levenberg-marquardt/fig3-four-methods-comparison.png" alt="네 방법 비교" width="600">
_그림3. 네 방법 비교 (d0=1540nm)_

| 방법 | $d_0=1540$nm 결과 (10스텝) |
|---|---|
| Gradient Descent (alpha=300) | 1490.12nm, 4~5스텝 근방 수렴 |
| Newton | 1666.31nm, 발산 |
| Gauss-Newton | 1490.12nm, 3스텝 수렴 |
| Levenberg-Marquardt ($\tau=10^{-6}$) | 1490.12nm, 3스텝 수렴 (GN과 사실상 동일) |

$\tau$를 충분히 작게 잡으면 LM은 GN과 구분할 수 없을 정도로 같은 궤적을 그린다. 당연한 결과다 — $\mu_0$이 $A$보다 몇 자릿수 작으면 정규방정식이 사실상 그대로이기 때문이다. LM이 진가를 발휘하는 지점은 지금부터다.

## 5. LM도 만능은 아니다 — basin 문제는 그대로

$\mu$는 정규방정식의 조건을 바꿀 뿐, 목적함수 자체의 모양을 바꾸지는 않는다. 그래서 2편에서 확인한 국소 최솟값 basin 문제는 LM에도 그대로 남는다. $d_0$를 1300nm, 1690nm 근방으로 바꿔서 LM($\tau=1$)을 돌려보면 GN과 똑같은 이웃 극소점으로 수렴한다.

| 초기값 $d_0$ | Gauss-Newton | Levenberg-Marquardt ($\tau=1$) |
|---|---|---|
| 1300nm | 1293.12nm | 1293.12nm |
| 1690nm | 1688.55nm | 1688.55nm |

$\tau$를 바꿔봐도(1e-3, 1, 100) 결과는 같았다. Damping은 "한 basin 안에서 안전하게 그 최솟값까지 가는 것"을 보장할 뿐, "어느 basin으로 갈지"는 여전히 초기값이 결정한다.

## 6. 실전 실패 사례 — 굴절률을 잘못 가정하면

2편이 예고한 GN의 이론적 한계(잔차가 크거나 Jacobian이 특이에 가까우면 실패)를 실제로 재현해봤다. 피팅 모델에서 SiO2 굴절률을 실제 값(1.46)이 아니라 1.02로 잘못 가정했다고 하자 — 공기(1.0)에 가까운 값이라 SiO2/Si 경계의 간섭 콘트라스트가 거의 사라진다. $d=1490$nm에서 $\sum_i r_i'^2$를 계산해보면 올바른 모델의 약 $2.58\times10^{-3}$에서 $4.34\times10^{-6}$로, 약 600배 작아진다. Jacobian이 거의 0에 가까워진 것이다.

이 상황에서 $d_0=1540$nm부터 Gauss-Newton을 돌리면 1666nm → 1543nm → 1417nm 세 값 사이를 20스텝을 넘도록 영원히 순환하며 전혀 수렴하지 않는다. 정규방정식의 분모($\sum_i r_i'^2$)가 지나치게 작아서 스텝이 매번 크게 오버슈트하기 때문이다. 같은 조건에서 Levenberg-Marquardt는 $\tau$ 값과 무관하게 10~20스텝 안에 이 (잘못된) 모델의 실제 최솟값 $d\approx1486.6$nm으로 안정적으로 수렴한다.

<img src="/assets/img/posts/optimization-levenberg-marquardt/fig4-gn-failure-lm-rescue.png" alt="GN 실패 vs LM 성공" width="600">
_그림4. 굴절률을 잘못 가정한 모델에서: GN 실패 vs LM 성공_

$\mu$가 분모에 더해지면서 정규방정식이 더는 특이에 가깝지 않게 되고, 스텝 크기가 자동으로 억제된다. Madsen이 이론으로 경고했던 실패 조건을, 딱 그 조건을 만족하는 예제를 만들어서 직접 확인한 셈이다.

## 정리

1편(gradient descent)부터 이번 3편까지 같은 두께 피팅 문제로 네 가지 방법을 살펴봤다.

| 방법 | 스텝 결정 방식 | 강점 | 약점 |
|---|---|---|---|
| Gradient Descent | $-\alpha J'(d)$ | 구현이 가장 간단 | 스텝 크기를 사람이 골라야 함 |
| Newton | $-J'(d)/J''(d)$ | 해 근방에서 2차 수렴 | 비볼록 구간에서 곡률 부호가 뒤집히면 발산 |
| Gauss-Newton | $-g/A$ ($J^TJ$ 근사) | $J^TJ$가 항상 준정부호라 곡률 문제에서 자유로움 | 잔차가 크거나 Jacobian이 특이에 가까우면 실패 |
| Levenberg-Marquardt | $-g/(A+\mu)$ | $\mu$가 GN↔steepest descent를 보간, GN의 실패 조건을 구제 | basin 문제는 여전히 남음, $\mu_0$ 튜닝 필요 |

Levenberg-Marquardt가 계측 소프트웨어의 사실상 표준인 이유는 이 표 한 줄로 요약된다 — Gauss-Newton의 속도를 대부분 그대로 가져가면서, Gauss-Newton이 무너지는 바로 그 지점에서 damping이 안전장치로 작동한다. 다만 어떤 방법도 "어느 국소 최솟값으로 갈지"까지는 정해주지 않는다는 점은 네 방법 모두에 공통으로 남는 숙제다.

## 참고자료

- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk 원문](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

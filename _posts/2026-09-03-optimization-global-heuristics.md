---
layout: post
title: "최적화 방법론 4편 — Global Optimization 휴리스틱: Basin-Hopping과 Simulated Annealing"
date: 2026-09-03 21:00:00 +0900
categories: [최적화 방법]
tags: [optimization, simulated-annealing, basin-hopping, global-optimization, metropolis, thin-film]
description: "국소 최적화가 못 넘는 basin 경계를 Simulated Annealing과 Basin-Hopping이 온도 하나로 어떻게 넘는지 실험으로 확인한다."
math: true
---

[3편](/posts/optimization-levenberg-marquardt/)에서 Levenberg-Marquardt(LM)법은 damping parameter $\mu$로 Gauss-Newton의 실패 조건(잔차가 크거나 Jacobian이 거의 특이할 때)을 구제한다는 걸 확인했다. 하지만 3편 말미에서 실험으로 확인한 대로, LM에는 $\mu$로도 못 고치는 한계가 있다 — 초기값 $d_0=1300$nm에서 시작하면 항상 $d=1293$nm 근처로, $d_0=1690$nm에서 시작하면 항상 $d=1688$nm 근처로 수렴했다. $\tau$(초기 damping 스케일)를 아무리 바꿔도 결과는 같았다. 국소 최적화 방법들은 전부 "지금 서 있는 골짜기"만 내려갈 뿐, 옆 골짜기에 더 낮은 바닥이 있는지는 알지 못한다.

이번 편은 이 "basin 문제"를 정면으로 다룬다. 국소 최적화기 자체를 고치는 게 아니라, 그 위에 무작위성과 확률적 수락 기준을 얹어서 나쁜 골짜기를 스스로 벗어나게 만드는 두 가지 고전적인 휴리스틱 — Simulated Annealing(SA)과 Basin-Hopping(BH) — 을 구현하고, 실제로 얼마나 잘 작동하는지 실험으로 확인한다. 1~3편이 국소 방법의 시리즈였다면, 4편은 그 시리즈를 감싸는 전역 탐색 레이어인 셈이다.

## 1. 지형을 다시 본다 — basin 세 개

지금까지 쓴 예제(SiO2/Si 단층, 실제 두께 1490nm, seed=0 노이즈)의 목적함수 $J(d)$를 1편 그림2보다 넓은 구간(1200~1800nm)에서 다시 그려보면, basin이 정확히 몇 개이고 경계가 어디인지 눈으로 확인할 수 있다.

![그림1](/assets/img/posts/optimization-global-heuristics/fig1-objective-landscape.png)
_그림1. 목적함수 지형 — 세 개의 basin (점선: basin 경계)_

이 구간 안에는 국소 최솟값이 정확히 세 개 있다: $d\approx1293$nm($J\approx0.867$), $d\approx1490$nm($J\approx0.0024$, 전역 최솟값), $d\approx1689$nm($J\approx1.059$). basin 경계(지형의 봉우리)는 $d\approx1391$nm와 $d\approx1589$nm에 있다 — 1편에서 어림잡았던 "대략 1390과 1590"과 정확히 일치한다. 전역 최솟값의 $J$가 다른 두 국소 최솟값보다 두 자릿수 이상 작다는 점도 눈여겨볼 만하다. 초기 위치가 이 200nm 폭짜리 경계 밖에 있으면, LM은 아무리 정교하게 스텝을 조정해도 자기가 서 있는 골짜기를 벗어날 방법이 없다.

## 2. Simulated Annealing — 랜덤워크 자체를 담금질한다

### 2.1 원리

Kirkpatrick, Gelatt, Vecchi(1983)가 제안한 SA는 금속을 천천히 식히면 원자가 전역적으로 가장 안정된(에너지가 가장 낮은) 결정 구조를 찾아간다는 물리적 담금질(annealing) 과정에서 이름을 땄다. 알고리즘은 단순하다 — 현재 위치 $d$에서 무작위로 제안한 $d'$을 Metropolis 기준으로 수락하거나 기각한다.

$$ P(\text{accept}) = \min\!\left(1,\ \exp(-\Delta J / T)\right), \qquad \Delta J = J(d') - J(d) $$

$J$가 줄어드는 제안은 항상 받아들이고, 늘어나는 제안도 확률적으로 받아들인다 — 이 "가끔 나쁜 방향으로도 움직이는" 성질 덕분에 얕은 골짜기에 갇히지 않고 빠져나올 수 있다. $T$(온도)는 이 확률을 조절하는 손잡이다. $T\to0$이면 개선되는 제안만 받아들이는 순수한 언덕 내려가기(greedy descent)가 되고, $T\to\infty$이면 사실상 무작위 워크(random walk)가 된다. 매 스텝 $T$를 기하급수적으로 냉각시켜($T_{k+1}=\text{cooling}\cdot T_k$), 처음에는 넓게 탐색하다가 점점 한 골짜기로 정착하도록 유도한다.

### 2.2 우리 문제에 적용

파라미터가 두께 $d$ 하나뿐이므로, 제안 분포도 단순한 가우시안 랜덤워크로 충분하다.

```python
# simulated_annealing.py 핵심부 (전체 코드: _code/optimization-global-heuristics/)
T = T0
for _ in range(n_iter):
    d_prop = d + rng.normal(0, step_sigma)
    J_prop = objective(d_prop, wavelength_nm, measured_R)
    dJ = J_prop - J
    if dJ < 0 or rng.random() < np.exp(-dJ / T):
        d, J = d_prop, J_prop
    if J < best_J:
        best_d, best_J = d, J
    T *= cooling
```

3편의 LM 코드는 여기서는 쓰지 않는다 — SA는 국소 최적화기 없이 담금질 자체로 basin을 넘나든다. 대신 현재 위치가 아니라 지금까지 방문한 것 중 최솟값(`best_d`, `best_J`)을 따로 추적한다. Metropolis 기준이 확률적으로 나쁜 스텝도 받아들이는 이상, "지금 서 있는 자리가 최선"이라는 보장이 없기 때문이다.

### 2.3 실험 — 온도 스케줄이 정말 문제인가

3편에서 LM이 실패했던 바로 그 시작점 $d_0=1300$nm에서, 냉각 속도(`cooling=0.97`)는 고정하고 초기 온도 $T_0$만 세 가지로 바꿔 300 iteration을 돌렸다.

![그림2](/assets/img/posts/optimization-global-heuristics/fig2-sa-temperature-schedule.png)
_그림2. SA 온도 스케줄에 따른 탐색 궤적 (d0=1300nm)_

- $T_0$가 매우 낮으면(0.0005): 초반 240여 스텝 동안 $d\approx1293$에 거의 그대로 머문다. 개선되지 않는 제안은 확률적으로도 거의 받아들이지 않기 때문이다. 하지만 낮은 확률로도 언젠가는 우연히 basin 경계를 넘는 제안이 나오고, 실제로 245번째 스텝 근처에서 전역 최솟값으로 건너뛴 뒤 바로 정착한다. 운이 나빴다면 300 스텝 안에 못 넘었을 수도 있다.
- $T_0$가 적당하면(0.02): 초반부터 옆 basin으로의 이동을 훨씬 자주 시도하고, 60번째 스텝 근처에서 전역 최솟값을 찾아 안정적으로 정착한다.
- $T_0$가 매우 높으면(5.0): $d$가 1240~1722nm 사이를 계속 넓게 방황한다. 결국은 전역 최솟값 근처에 도달하지만(최종 $J\approx0.0025$), 궤적이 지저분하고 여러 basin을 반복해서 들락거린다 — 냉각이 끝날 때까지 한 골짜기에 제대로 정착하지 못할 위험이 크다.

세 경우 모두 이번 300-스텝 예산 안에서는 전역 최솟값을 찾아냈지만, $T_0$가 너무 낮으면 "우연히 탈출할 때까지 기다리는" 운에 의존하게 되고, 너무 높으면 "이미 찾은 곳에 머무르지 못하고 계속 떠도는" 비효율이 생긴다. 사실상 다음 편으로 미룰 것 없이 이번 실험만으로도 확인된다 — 적당한 $T_0$는 "탐색"과 "수렴"의 균형점을 찾는 문제이지, 아무 값이나 크게 잡으면 되는 게 아니다.

## 3. Basin-Hopping — 국소 최적화를 담금질한다

### 3.1 원리

Wales와 Doye(1997)가 제안한 basin-hopping은 SA와 목적이 같지만 구조가 다르다. 매번 좌표 $d$ 자체를 담금질하는 대신, "지금 어느 basin의 바닥에 있는가"를 담금질한다.

1. 현재 해 $d$에 큰 무작위 섭동을 준다: $d_0' = d + \mathcal{N}(0,\sigma^2)$
2. $d_0'$에서 국소 최적화를 수렴할 때까지 돌려 그 basin의 바닥 $d'$을 찾는다
3. 두 바닥 $J(d)$와 $J(d')$ 사이를 Metropolis 기준으로 비교해 수락/기각한다

SA가 매 스텝 임의의(반드시 바닥이 아닌) 위치를 비교하는 반면, basin-hopping은 항상 "바닥과 바닥"을 비교한다. 그래서 지형의 미세한 굴곡(노이즈로 인한 잔물결 등)에 흔들리지 않고, basin 단위로 훨씬 깨끗하게 움직인다.

### 3.2 구현 — 3편 LM을 로컬 스텝으로 재사용

여기서 2단계의 "국소 최적화"가 정확히 3편의 `levenberg_marquardt()`다. 시리즈 전체가 여기서 만난다 — 1~2편이 실패했던 지점을 3편이 damping으로 구제했고, 이번 편은 그 3편 함수를 통째로 가져다 "한 basin을 빠르고 정확하게 바닥까지 내려가는 도구"로 재사용한다.

```python
# basin_hopping.py 핵심부 (전체 코드: _code/optimization-global-heuristics/)
d_hist_lm, J_hist_lm, _ = levenberg_marquardt(d0, wavelength_nm, measured_R, **lm_kwargs)
d, J = d_hist_lm[-1], J_hist_lm[-1]
T = T0
for _ in range(n_hops):
    d_trial0 = d + rng.normal(0, perturb_sigma)
    d_hist_trial, J_hist_trial, _ = levenberg_marquardt(d_trial0, wavelength_nm, measured_R, **lm_kwargs)
    d_trial, J_trial = d_hist_trial[-1], J_hist_trial[-1]
    dJ = J_trial - J
    if dJ < 0 or rng.random() < np.exp(-dJ / T):
        d, J = d_trial, J_trial
    if J < best_J:
        best_d, best_J = d, J
```

섭동 크기(`perturb_sigma`)는 basin 간격(~200nm, 그림1)보다 커야 다른 basin으로 건너뛸 여지가 생긴다. 여기서는 150nm를 썼다.

### 3.3 실험 — 온도가 너무 높으면 오히려 실패한다

같은 $d_0=1300$nm에서 25번의 hop을 시도했다.

![그림3](/assets/img/posts/optimization-global-heuristics/fig3-basinhopping-temperature-schedule.png)
_그림3. Basin-hopping 온도 스케줄에 따른 궤적 (d0=1300nm)_

- $T_0$가 매우 낮거나(0.0005) 적당하면(0.05): 두 경우 궤적이 완전히 같다(파란 선이 초록 선 뒤에 가려져 있다). 25번 중 17번의 hop을 수락했고, 14번째 hop에서 $d=1293$에서 $d=1490$으로 건너뛴 뒤 그대로 정착한다. 이 문제는 basin이 세 개뿐이고 간격이 넓어서, 낮은 온도의 "개선되는 hop만 받아들이는" 전략만으로도 충분히 잘 작동한다.
- $T_0$가 매우 높으면(5.0): 24번 중 24번, 즉 사실상 모든 hop을 받아들인다. 그 결과 $d=1293 \to 1022 \to 696 \to 530 \to 441$로 계속 더 나쁜 basin으로 튕겨 나간다. 전역 최솟값 근처를 지나간 적조차 없어서, `best_J` 역시 출발 basin의 값(0.867)에 머물러 있다.

이 결과는 애초에 예상했던 것과 다르다 — "온도가 너무 낮으면 multi-start처럼 될 것"이라는 예상은 이 문제에서는 틀렸다(오히려 낮은 온도 쪽이 더 잘 작동했다). 실제로 걸림돌이 된 쪽은 반대였다: basin-hopping에서 온도가 너무 높으면 나쁜 basin으로의 이동도 거의 다 받아들이게 되어, 이미 찾은 좋은 해를 스스로 걷어차고 떠나 버린다. SA와 basin-hopping은 같은 "온도"라는 이름을 쓰지만, SA에서는 고온이 "탐색 범위 확대"로 끝나는 반면 basin-hopping에서는 고온이 "이미 찾은 답을 잃어버리는" 훨씬 직접적인 손해로 이어진다. Metropolis 비교 대상이 basin의 바닥이라 스텝 하나하나의 낙폭이 훨씬 크기 때문이다.

## 4. 정량 비교 — 몇 번 만에, 얼마나 자주 찾는가

시작점 하나만 보는 대신, $d_0$를 $[1200, 1800]$ 구간에서 50번 무작위로 뽑아 세 가지 방법의 성공률을 비교했다. "성공"은 전역 최솟값(1490nm)으로부터 5nm 이내에 도달한 경우로 정의했다.

![그림4](/assets/img/posts/optimization-global-heuristics/fig4-success-rate-comparison.png)
_그림4. 무작위 시작 50회 성공률 비교_

- 단일 LM(재시도 없이 1회): 18/50 = 36%. 대략 basin 폭(200nm)이 전체 탐색 구간(600nm)에서 차지하는 비율과 맞아떨어진다 — 초기값이 정답의 basin 안에 있을 때만 성공하는 국소 방법의 한계를 그대로 보여준다.
- SA(300 iteration 예산): 50/50 = 100%. 처음 5nm 이내에 도달하기까지 평균 61.8 iteration이 걸렸다.
- Basin-hopping(15 hop 예산): 50/50 = 100%. 처음 도달하기까지 평균 2.2 hop이 걸렸다.

두 전역 탐색법 모두 이번 예산 안에서는 실패가 없었다. 다만 "iteration 수"와 "hop 수"를 그대로 비교하면 안 된다 — SA의 한 스텝은 목적함수 한 번 계산이지만, basin-hopping의 한 hop은 그 안에서 LM이 여러 번 반복(Jacobian 계산 포함)해 바닥까지 내려가는 과정을 통째로 포함한다. 그래서 이번 문제처럼 국소 구조(매끄러운 골짜기, 정확한 로컬 옵티마이저)가 뚜렷할 때는 basin-hopping이 훨씬 적은 "바깥쪽 스텝"만으로 끝나고, SA는 골짜기 바닥 근처에서도 계속 잔걸음으로 미세조정해야 해서 스텝 수 자체는 더 많이 든다.

## 5. 지형 위에서 보는 basin-hopping

그림3의 $T_0=0.05$ 궤적을 그림1의 지형 위에 그대로 겹쳐보면, 이 알고리즘이 하는 일이 "국소 최솟값들 사이를 직접 건너뛰는 것"이라는 게 한눈에 보인다.

![그림5](/assets/img/posts/optimization-global-heuristics/fig5-trajectory-on-landscape.png)
_그림5. 지형 위에 겹친 basin-hopping 궤적_

파란 사각형이 $d_0=1300$에서 LM으로 수렴한 첫 번째 바닥($d=1293$)이고, 초록 별이 14번의 hop 뒤 도달해 계속 유지한 전역 최솟값이다. LM 혼자였다면 이 그림에서 파란 사각형 지점에 영원히 머물렀을 것이다. Basin-hopping은 그 위에 무작위 섭동으로 "다른 봉우리를 넘어가 보는" 시도를 계속 반복하고, 넘어간 곳이 더 낫다고 판단되면(Metropolis) 그쪽으로 옮겨간다 — 지형을 훑는 게 아니라, 봉우리 몇 개를 건너뛰며 바닥끼리 비교하는 것이다.

## 6. 정리 — 언제 국소 방법으로 충분한가

1~3편에서 gradient descent, Newton, Gauss-Newton, LM 순서로 다듬어 온 것은 전부 "이미 맞는 basin 안에 있다"는 전제 위에서의 이야기였다. 이번 편에서 확인했듯, 그 전제가 깨지는 순간 국소 방법은 damping을 아무리 정교하게 조정해도 답이 없다 — 애초에 다른 골짜기를 볼 방법이 없기 때문이다.

실무적으로는 다음 기준이 쓸모 있다.

- **초기값에 대한 사전 정보가 있는가.** 이전 측정값, 공정 스펙, 대략적인 하드웨어 설계값처럼 초기값을 basin 폭 안으로 좁힐 수 있는 정보가 있다면 LM 한 번으로 충분하다. 이번 예제에서 basin 폭은 200nm 남짓이었다 — 이 정도 정밀도의 사전 정보조차 없는 경우는 생각보다 흔치 않다.
- **목적함수 계산이 얼마나 비싼가.** SA와 basin-hopping 모두 국소 방법보다 목적함수(및 그 이상)를 훨씬 더 많이 계산한다. 여기서는 반사율 계산이 사실상 공짜라 문제가 안 됐지만, 다층 박막의 전달행렬 계산이나 FDTD처럼 한 번 평가에 초 단위가 걸리는 모델이라면 전역 탐색의 비용은 무시할 수 없다.
- **basin이 몇 개이고 얼마나 벌어져 있는지 아는가.** 이번 실험은 basin이 세 개뿐이고 폭도 넓어서 낮은 온도로도 잘 풀렸다. basin이 훨씬 많고 촘촘하다면(파라미터가 여러 개로 늘어나는 경우 특히) 온도와 섭동 크기를 더 신중하게 잡아야 하고, 그림3에서 본 "고온이 답을 잃어버리는" 실패 모드도 더 쉽게 나타날 것이다.
- **SA vs basin-hopping.** 국소 구조가 뚜렷하고 신뢰할 만한 로컬 옵티마이저(LM 같은)가 있다면 basin-hopping 쪽이 적은 바깥쪽 스텝으로 끝난다. 반대로 지형이 매끄럽지 않거나 좋은 로컬 옵티마이저를 설계하기 어렵다면, 로컬 최적화 단계 없이도 동작하는 SA가 더 단순하고 안전한 선택이다.

1편의 gradient descent에서 시작해 4편의 basin-hopping까지, 이 시리즈는 결국 "국소 정보(기울기, 곡률)를 얼마나 똑똑하게 쓸 것인가"와 "전역 구조를 얼마나 넓게 볼 것인가"라는 두 축 사이의 트레이드오프를 순서대로 밟아온 셈이다. 실제 계측 소프트웨어에서는 흔히 이 둘을 함께 쓴다 — 대략적인 위치를 전역 탐색으로 찾고, 마지막 정밀도는 LM 같은 국소 방법에 맡기는 식이다.

## 참고자료

- S. Kirkpatrick, C. D. Gelatt, M. P. Vecchi, "Optimization by Simulated Annealing," _Science_, 220(4598), 671-680, 1983.
- D. J. Wales, J. P. K. Doye, "Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms," _J. Phys. Chem. A_, 101(28), 5111-5116, 1997.
- C. P. Chang, Y. H. Lee, S. Y. Wu, "Optimization of a thin-film multilayer design by use of the generalized simulated-annealing method," _Opt. Lett._, 15(11), 595-597, 1990.
- [scipy.optimize.basinhopping 문서](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html)

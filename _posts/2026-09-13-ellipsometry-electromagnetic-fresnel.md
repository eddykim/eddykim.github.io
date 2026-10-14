---
title: "엘립소메트리 배경이론 1편 — 검출기는 세기만 재는데 어떻게 빛의 위상을 알아내는가"
date: 2026-09-13 20:00:00 +0900
categories: [광학, 전자기학]
page_id: ellipsometry-electromagnetic-fresnel
tags: [ellipsometry, electromagnetics, maxwell-equations, fresnel-equations, polarization, thin-film]
description: 맥스웰 방정식에서 프레넬 반사계수까지의 수학 구조와, 검출기가 세기만 측정한다는 사실이 왜 엘립소메트리라는 계측법을 낳는지를 정리한다.
math: true
---

엘립소메트리(ellipsometry)는 빛의 편광 상태가 시편 표면에서 어떻게 바뀌는지를 측정해 박막의 두께와 굴절률을 알아내는 계측법이다. 이 계측법을 이해하려면 결국 "빛이 매질 경계면에서 어떻게 반사·굴절되는가"라는 전자기학의 가장 기본적인 질문으로 돌아가야 한다. 이번 시리즈는 그 배경이론을 3편에 걸쳐 정리한다. 1편(이번 글)은 맥스웰 방정식(Maxwell's equations)에서 출발해 프레넬 방정식까지의 수학적 구조를 따라간다. [2편](/posts/ellipsometry-polarization-mueller-matrix/)은 편광 상태를 기술하는 스토크스 벡터(Stokes vector)와 뮬러 행렬(Mueller matrix)을, [3편](/posts/ellipsometry-thin-film-multilayer-reflectance/)은 박막 다층 구조의 반사와 엘립소메트리 파라미터 $\Psi, \Delta$를 다룬다.

이번 글에서 특히 눈여겨볼 지점은 하나다. 프레넬 방정식 자체는 광학 교재 어디에나 나오는 표준적인 결과다. 하지만 "그래서 왜 굳이 편광 상태를 비교하는 엘립소메트리라는 계측법이 필요한가"라는 질문에 답하려면, 그 방정식이 유도되는 과정 중간의 포인팅 벡터(Poynting vector) 대목을 짚고 넘어가야 한다. 검출기가 실제로 측정하는 물리량은 전기장 자체가 아니다. 그 때문에 잃어버리는 정보가 무엇인지가 이 시리즈 전체를 관통하는 동기다.

## 1. 맥스웰 방정식에서 어떻게 파동이 나오는가

빛이 전자기파라는 말은 어디에서 나오는가. 출발점은 맥스웰 방정식이다. 자유전하와 자유전류가 없는 진공에서 맥스웰 방정식은 다음과 같은 형태를 갖는다.

$$ \nabla \cdot \mathbf{E} = 0, \qquad \nabla \cdot \mathbf{B} = 0 $$

$$ \nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t}, \qquad \nabla \times \mathbf{B} = \mu_0 \epsilon_0 \frac{\partial \mathbf{E}}{\partial t} $$

여기서 $\mathbf{E}$는 전기장, $\mathbf{B}$는 자기장, $\mu_0$와 $\epsilon_0$는 각각 진공의 투자율(permeability)과 유전율(permittivity)이다. 이 네 식을 결합하면 전기장과 자기장이 각각 만족하는 파동 방정식(wave equation)이 나온다. 두 회전 방정식 각각에 다시 회전 연산자를 취하고, 항등식 $\nabla\times(\nabla\times\mathbf{A}) = \nabla(\nabla\cdot\mathbf{A}) - \nabla^2\mathbf{A}$에 발산 조건($\nabla\cdot\mathbf{E}=\nabla\cdot\mathbf{B}=0$)을 대입하면 된다.

$$ \nabla^2 \mathbf{E} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2} = 0, \qquad \nabla^2 \mathbf{B} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{B}}{\partial t^2} = 0 $$

즉, 전기장과 자기장은 진공의 두 상수만으로 정해지는 속도 $c = 1/\sqrt{\mu_0\epsilon_0} \approx 3.00\times10^{8}\,\mathrm{m/s}$로 퍼지는 파동이다. 이 값이 다름 아닌 빛의 속도다.

이 파동 방정식의 해는 평면파(plane wave)의 형태로 주어진다.

$$ \mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0\, e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}, \qquad \mathbf{B}(\mathbf{r}, t) = \mathbf{B}_0\, e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})} $$

여기서 $\mathbf{k}$는 파수 벡터(wave vector), $\omega$는 각진동수(angular frequency)이며, 둘 사이에는 $\lvert \mathbf{k} \rvert = \omega/c$ 관계가 성립한다. 시간 인자를 $e^{j\omega t}$로 두는 것은 1968년 네브래스카 학회 이후 엘립소메트리 문헌에서 널리 쓰는 관례(Azzam & Bashara, Fujiwara)를 따른 것이다. 다만 문헌 전체가 이 관례로 통일되어 있지는 않다. 이 부호 규약은 3절에서 복소굴절률을 정의할 때 다시 문제가 되므로 미리 고정해 둔다.

지수함수로 쓰는 것은 계산이 편하기 때문이다. 실제 물리량은 이 복소수 표현의 실수부다. 하지만 미분과 위상 계산은 지수함수 꼴에서 훨씬 간단하므로, 뒤에 나올 프레넬 방정식까지 이 표기를 그대로 끌고 간다. 이 대목에서 이미 중요한 사실 하나가 드러난다. 전기장은 진폭 $\mathbf{E}_0$뿐 아니라 위상 $\omega t - \mathbf{k}\cdot\mathbf{r}$ 정보를 함께 담고 있는 복소수 값이라는 점이다.

## 2. 검출기는 전기장의 무엇을 보는가

전기장이 진폭과 위상을 함께 담는다면, 검출기는 그중 무엇을 읽는가. 답을 얻으려면 전자기파가 나르는 에너지부터 봐야 한다. 그 에너지 흐름은 포인팅 벡터 $\mathbf{S}$로 나타낸다.

$$ \mathbf{S} = \frac{1}{\mu_0}(\mathbf{E} \times \mathbf{B}) $$

가시광의 각진동수 $\omega$는 $10^{15}\,\mathrm{rad/s}$ 자릿수다. 633 nm 빛이라면 한 주기가 약 2.1 fs다. 어떤 검출기(포토다이오드, CCD, CMOS 무엇이든)도 이 속도로 진동하는 전기장 자체를 실시간으로 따라갈 수 없다. 검출기가 실제로 반응하는 물리량은 노출 시간 동안 적분된, 즉 시간 평균된 에너지 흐름이다.

$$ \langle \mathbf{S} \rangle = \frac{1}{2}\,\mathrm{Re}(\mathbf{E} \times \mathbf{H}^*) $$

여기서 $\mathbf{H} = \mathbf{B}/\mu_0$이고 $\mathbf{E}, \mathbf{H}$는 복소 진폭이다. 진공의 평면파로 이 식을 전개하면, 검출기가 측정하는 세기(intensity) $I$는 전기장 진폭의 절댓값 제곱에 비례하는 스칼라량으로 정리된다.

$$ I = \lvert \langle \mathbf{S} \rangle \rvert = \frac{1}{2}\sqrt{\frac{\epsilon_0}{\mu_0}}\, \lvert \mathbf{E}_0 \rvert^2 $$

굴절률 $n$인 매질 안에서는 우변에 $n$이 곱해지지만, 세기가 $\lvert \mathbf{E}_0 \rvert^2$에 비례한다는 점은 같다.

이 식이 이 시리즈 전체의 출발점이다. 세기 $I$는 진폭 $\lvert \mathbf{E}_0 \rvert$만 남기고, 위상은 제곱과 시간 평균 과정에서 완전히 지워 버린다. 즉, 진폭이 같으면 위상이 다른 두 빛도 검출기에는 똑같이 보인다. 반사율이나 투과율처럼 세기만 재는 계측을 흔히 포토메트리(photometry)라 부른다. 이런 계측은 한 파장에서 잰 세기만으로는 위상 정보에 접근할 수 없다. 넓은 대역의 반사율 스펙트럼에 크라머스-크로니히 관계(Kramers–Kronig relations)를 적용해 위상을 되살리는 방법이 있기는 하다. 하지만 측정 대역 밖의 외삽과 벌크 시편이라는 조건이 붙는다. 이 지점으로는 7절에서 다시 돌아온다.

## 3. 빛은 경계면에서 왜 꺾이는가

이제 이 파동이 매질 경계면을 만나면 어떻게 되는지 보자. 평면파가 굴절률이 다른 두 매질의 경계면에 입사하면 일부는 반사되고 일부는 굴절되어 두 번째 매질로 투과한다. 이 현상은 새로운 물리 법칙이 아니다. 경계면에서 $\mathbf{E}$와 $\mathbf{H}$의 접선 성분이 연속이어야 한다는 맥스웰 방정식의 경계 조건(boundary condition)에서 그대로 유도된다. 이 글처럼 비자성 매질을 다룰 때는 $\mathbf{H}$ 대신 자기장 $\mathbf{B}$의 접선 성분이 연속이라고 해도 같다.

입사광의 진행 방향과 경계면 법선이 이루는 평면을 입사면(plane of incidence)이라 한다. 모든 편광 성분은 이 평면에 수직인 s편광(s-polarization)과 평행한 p편광(p-polarization)으로 나누어 다룬다. 입사각 $\theta_i$, 반사각 $\theta_r$, 굴절각 $\theta_t$를 각각 법선과 이루는 각으로 정의한다. 그림1은 이 정의를 한 장에 모은 것이다. ⊙ 기호는 지면에서 튀어나오는 방향의 진동(s편광)을, 양방향 화살표는 지면 안에서 광선에 수직인 진동(p편광)을 뜻한다.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/fig1-concept-diagram.png" alt="입사면 위의 s/p 편광과 입사·반사·굴절 광선" width="600">
_그림1. 입사면(지면) 위의 입사·반사·굴절 광선. s편광은 지면에 수직으로(⊙), p편광은 입사면 안에서 광선에 수직으로(↔) 진동한다. 매질 2의 굴절률이 더 커서 굴절광이 법선 쪽으로 꺾인다(도식용 각도 $$\theta_i = 35°$$, $$\theta_t = 22°$$)._

반사각은 언제나 입사각과 같고, 굴절각은 스넬의 법칙(Snell's law)으로 결정된다.

$$ \theta_i = \theta_r, \qquad N_1 \sin\theta_i = N_2 \sin\theta_t $$

여기서 $N_1, N_2$는 각 매질의 굴절률이다. 흡수가 있는 실제 물질은 복소굴절률(complex refractive index) $N = n - jk$로 표현한다. 실수부 $n$은 통상적인 굴절률이고, 허수부 $k$는 물질이 빛을 얼마나 흡수하는지를 나타내는 소광계수(extinction coefficient)다. 굵은 글자 $\mathbf{k}$(파수 벡터)와 글자가 같으니 혼동하지 않기 바란다. 부호가 마이너스인 것은 1절에서 고정한 $e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}$ 규약과 짝을 이룬다. 이 규약에서 파수의 크기 $\lvert \mathbf{k} \rvert = N\omega/c$를 대입하면 진행 방향으로 $e^{-k\omega z/c}$의 감쇠가 나오기 때문이다. 시간 인자를 $e^{-j\omega t}$로 두는 물리학 쪽 관례를 따른다면 같은 감쇠를 얻기 위해 $N = n + jk$로 써야 한다. 어느 쪽이든 물리는 같다. 하지만 두 규약을 섞으면 흡수가 증폭으로 뒤집히므로, 글 전체에서 하나로 통일한다.

이 스넬의 법칙은 [광선추적 1편](/posts/raytracing-spherical-lens-refraction/)에서 렌즈 굴절 계산에 그대로 가져다 썼던 바로 그 식이다. 그 글에서는 이 식을 주어진 결과로만 쓰고 왜 성립하는지는 다루지 않았다. 답은 여기 있다. 경계면 위의 한 점에서 위쪽 파동과 아래쪽 파동이 서로 다른 위상을 가질 수는 없다. 그러니 경계면을 따라 잰 파면의 마루 간격은 양쪽이 같아야 한다. 그런데 굴절률이 큰 매질 2에서는 파장이 짧아진다. 파면이 기울어야만 같은 간격을 맞출 수 있고, 그 기울기가 곧 굴절각이다.

식으로 쓰면 이렇다. 평면파 해 $e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}$가 경계면 위의 모든 점, 모든 시각에서 입사파·반사파·굴절파 세 파동에 걸쳐 같은 위상으로 이어져야 한다. 그러려면 각 파수 벡터를 경계면 방향으로 투영한 성분이 서로 같아야 한다. 이 조건을 위상 정합(phase matching)이라 한다. 각도로 풀어 쓰면 $\lvert \mathbf{k}_1 \rvert \sin\theta_i = \lvert \mathbf{k}_2 \rvert \sin\theta_t$이고, $\lvert \mathbf{k} \rvert = N\omega/c$를 넣으면 정확히 $N_1\sin\theta_i = N_2\sin\theta_t$가 된다. 즉, 위상 정합 자체가 스넬의 법칙이다. 그림2가 이 과정을 파면으로 보여 준다.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/ext-wavefront-refraction.png" alt="경계면을 지나며 꺾이는 평면파의 파면" width="450">
_그림2. 굴절률 $$n_1$$인 매질에서 $$n_2$$인 매질로 들어가는 평면파. 굵은 회색 선은 빛줄기의 양쪽 가장자리이고, 가는 선(위는 파랑, 아래는 초록)은 파면이다. 파면은 경계면의 노란 점에서 끊기지 않고 이어지며, 아래 매질에서 파면 간격이 좁아진 만큼 진행 방향(점선 화살표)이 법선 쪽으로 꺾인다. (출처: 老陳, [Wavefront Refraction](https://commons.wikimedia.org/wiki/File:Wavefront_Refraction.svg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

## 4. 경계면에서 빛은 얼마나 반사되는가

스넬의 법칙은 굴절이 "어느 각도로" 일어나는지만 알려 준다. 경계면에서 빛이 "얼마나" 반사되고 투과되는지는 알려 주지 않는다. 그 비율을 정하는 것이 프레넬 방정식(Fresnel equations)이다. 이 방정식이 주는 계수는 진폭비와 위상차를 함께 담는다는 점에서 이 시리즈의 중심에 놓인다. 3절의 접선 성분 연속 조건을 s편광과 p편광에서 각각 풀면, 반사파와 입사파의 전기장 진폭비로 정의되는 반사계수(reflection coefficient) $r$이 다음과 같이 얻어진다.

$$ r_p = \frac{N_2\cos\theta_i - N_1\cos\theta_t}{N_2\cos\theta_i + N_1\cos\theta_t}, \qquad r_s = \frac{N_1\cos\theta_i - N_2\cos\theta_t}{N_1\cos\theta_i + N_2\cos\theta_t} $$

투과계수 $t_p, t_s$도 같은 방식으로 얻어지며, 3편에서 다룰 박막 내부의 다중 반사를 계산할 때 다시 등장한다. 여기서 주목할 부분은 $N_1, N_2$가 일반적으로 복소수라는 점이다. 흡수가 있는 매질이라면 $N=n-jk$이므로 $r_p, r_s$ 역시 복소수다. 복소수인 반사계수는 크기(진폭비)와 편각(위상차) 두 정보를 동시에 담는다. 이것이 2절에서 짚은 "세기 측정은 위상을 잃는다"는 문제를 풀 첫 실마리다. 즉, 반사계수 $r$에는 위상이 처음부터 살아 있다.

극한부터 보자. 수직 입사($\theta_i = 0$)에서는 $\cos\theta_i = \cos\theta_t = 1$이므로 $r_p = (N_2-N_1)/(N_2+N_1)$, $r_s = -(N_2-N_1)/(N_2+N_1)$이다. 공기($N_1=1.0$)에서 유리($N_2=1.5$)로 들어가면 $r_p = +0.2$, $r_s = -0.2$이고, 반사율은 둘 다 $0.2^2 = 0.04$다. 수직 입사에서는 s와 p를 가를 방향이 없으니 두 편광은 물리적으로 똑같이 반사된다. $r_p = -r_s$라는 부호 차이는 p편광의 기준 방향을 잡는 규약에서 온다. 이 글은 Born & Wolf, Fujiwara와 같은 규약을 쓰고, 2·3편도 이를 따른다. 이 규약에서는 수직 입사일 때 두 반사파의 위상차 180°가 기준값이 된다.

## 5. p편광 반사율은 왜 한 각도에서 0이 되는가

반사계수가 위상을 담는다는 것은 알았다. 그렇다면 세기로 보이는 반사율은 입사각에 따라 어떻게 바뀌는가. 실제로 측정 가능한 반사율(reflectance) $R$은 반사계수의 절댓값 제곱이다.

$$ R_p = \lvert r_p \rvert^2, \qquad R_s = \lvert r_s \rvert^2 $$

공기($N_1=1.0$)에서 유리($N_2=1.5$)로 들어가는 빛의 $R_p, R_s$를 입사각에 따라 계산하면 그림3과 같다.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/fig2-fresnel-reflectance.png" alt="프레넬 반사율 R_p, R_s를 입사각에 따라 그린 그래프, 브루스터각 표시" width="600">
_그림3. 공기($$N_1=1.0$$)에서 유리($$N_2=1.5$$)로 입사할 때의 프레넬 반사율. $$R_s$$(보라)는 입사각에 따라 단조증가하고, $$R_p$$(초록)는 브루스터각 56.3°(점선)에서 0이 된다. 수직 입사에서는 둘 다 0.04다._

두 곡선은 수직 입사의 0.04에서 함께 출발한다. $R_s$는 입사각이 커질수록 단조증가하지만, $R_p$는 줄어들다가 특정 각도에서 0을 찍고 다시 올라간다. 이 각도를 브루스터각(Brewster's angle)이라 한다. $r_p$의 분자가 0이 되는 조건 $N_2\cos\theta_i = N_1\cos\theta_t$를 스넬의 법칙과 함께 풀면 $\theta_B = \arctan(N_2/N_1)$를 얻는다. 공기-유리 조합에서는 56.3°다. 브루스터각에서는 p편광 성분이 전혀 반사되지 않고 반사광에는 s편광만 남는다. 이때 반사광과 굴절광은 서로 직각을 이룬다(유리에서 굴절각은 33.7°라 56.3° + 33.7° = 90°). 이 성질은 실제 계측 현장에서 편광자 방위각(azimuth)의 기준, 곧 입사면 방향을 잡는 보정에 활용된다. 그림4가 이 상황을 보여 준다.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/ext-brewster-angle.png" alt="브루스터각으로 입사한 빛의 반사와 굴절, 반사광은 s편광만 남는다" width="500">
_그림4. 브루스터각 $$\theta_B$$로 입사한 무편광 빛(incident ray)의 반사(reflected ray)와 굴절(refracted ray). 점은 지면에 수직인 진동(s편광), 양방향 화살표는 입사면 안의 진동(p편광)이다. 반사광에는 점만 남아 s편광이 되고, 굴절광은 두 성분을 모두 지녀 약간만 편광되며, 반사광과 굴절광 사이는 직각이다. (출처: Pajs, [Brewsters-angle](https://commons.wikimedia.org/wiki/File:Brewsters-angle.svg), 퍼블릭 도메인)_

그림3이 보여 주지 못하는 것이 하나 있다. $r_p$의 부호다. 유리에서 $r_s$는 입사각이 커질수록 음수 쪽으로 커진다(45°에서 −0.303, 70°에서 −0.547). 반면 $r_p$는 줄어들다가 브루스터각에서 0을 지나 음수가 된다(45°에서 +0.092, 70°에서 −0.206). 즉, 브루스터각을 지나는 순간 두 반사파의 위상차가 180°에서 0°로 뒤집힌다. 그런데 반사율로는 이 변화가 보이지 않는다. 45°와 63.6°의 $R_p$는 0.0085로 같지만, 그 두 각에서 $r_p$는 +0.092와 −0.092로 부호가 반대다. 세기는 위상을 지운다는 2절의 결론이 여기서 숫자로 드러난다.

흡수가 있는 시편에서는 사정이 조금 다르다. $N$이 복소수이면 $r_p$의 분자가 정확히 0이 되지 않으므로, $R_p$는 0이 아니라 최솟값만 갖는다. 이 각을 유사 브루스터각(pseudo-Brewster angle)이라 한다. Si(633 nm, $N \approx 3.88 - 0.016j$)에서는 약 75.6°에서 $R_p \approx 4\times10^{-6}$이 최솟값이다. 그 근처에서 위상차가 가장 빠르게 바뀐다. 70°에서 179.4°, 75°에서 173.9°이던 위상차가 75.6° 근처에서 90°를 지나고, 80°에서는 0.6°가 된다. 즉, 유사 브루스터각 근처는 시편의 작은 변화가 위상차에 가장 크게 실리는 자리다. 엘립소미터가 Si 기판 위의 시편을 잴 때 입사각을 70°대에 두는 이유가 이것이다. 경사 입사가 왜 필요한지, 장비가 이 각도를 어떻게 쓰는지는 [타원계측기 1편](/posts/ellipsometer-oblique-rotating-element/) 1절에서 다룬다.

## 6. 흡수하는 기판은 왜 반사로 재는가

지금까지는 경계면에서 갈라지는 빛을 보았다. 이번에는 매질 안으로 들어간 빛이 어떻게 되는지 보자. 흡수가 있는 매질 안을 진행하는 빛의 세기는 지수적으로 줄어들고, 그 빠르기는 소광계수 $k$가 정한다(Beer–Lambert 법칙).

$$ I(z) = I_0\, e^{-\alpha_{abs} z}, \qquad \alpha_{abs} = \frac{4\pi k}{\lambda} $$

여기서 $\lambda$는 진공 파장이다. 입사 세기가 $1/e$로 줄어드는 깊이를 침투 깊이(penetration depth) $\delta$라 부른다.

$$ \delta = \frac{1}{\alpha_{abs}} = \frac{\lambda}{4\pi k} $$

이 식은 박막 계측이 왜 대부분 반사형으로 이루어지는지를 설명해 준다. 흡수가 있는 파장 영역에서는 $\delta$가 기판 두께보다 훨씬 짧아진다. 실리콘은 633 nm에서 약 3 µm인데 웨이퍼는 수백 µm다. 빛이 시편을 투과하지 못하므로 투과형 측정 자체가 불가능해진다. 실리콘 위 박막처럼 흡수가 있는 기판을 다루는 엘립소메트리가 반사형 구성을 기본으로 하는 이유 중 하나가 여기에 있다.

## 7. 세기 대신 편광을 재면 무엇을 얻는가

이제 2절에서 미뤄 둔 이야기를 마저 할 수 있다. 검출기가 재는 세기 $I \propto \lvert E_0 \rvert^2$은 반사계수 $r$의 크기 $\lvert r \rvert$만 담는다. $r$의 위상은 세기 측정 과정에서 사라진다. 반사율 $R=\lvert r \rvert^2$을 한 파장에서 아무리 정밀하게 측정해도 s편광과 p편광 반사파 사이의 위상차에는 직접 접근할 수 없다. 그런데 5절에서 보았듯이 이 위상차가 바로 박막의 두께와 광학 상수에 민감하게 반응하는 정보다.

엘립소메트리가 세기 대신 편광 상태 자체의 변화를 재는 이유가 여기서 나온다. 입사광을 이미 편광 상태를 아는 빛(보통 선편광)으로 준비하고, 반사 후 편광 상태가 어떻게 바뀌었는지를 측정하면 s편광과 p편광 반사계수의 비율

$$ \rho = \frac{r_p}{r_s} = \tan\Psi \, e^{j\Delta} $$

을 얻을 수 있다. 여기서 $\tan\Psi$는 두 반사계수의 진폭비, $\Delta$는 위상차다. 4·5절에서 따라간 두 반사파의 위상차가 바로 이 $\Delta$다. 반사율 측정으로는 $\lvert r_p \rvert$와 $\lvert r_s \rvert$를 따로 얻을 수 있을 뿐이다. 반면 편광 상태의 변화는 두 편광 성분을 서로 견주는 상대적인 양이므로, 크기비와 위상차를 한 번에 준다. 즉, 위상에 직접 접근할 수 없다는 2절의 한계를 두 편광 성분의 비교로 우회하는 것이다. 엘립소메트리는 절대 세기를 버리고 두 편광의 비율을 택한 계측이다.

이 $\Psi, \Delta$ 한 쌍이 엘립소메트리가 실제로 내놓는 결과다. 여기서 박막의 두께와 굴절률을 거꾸로 구하려면 두 단계가 필요하다. 먼저 두께와 굴절률이 주어졌을 때 $\Psi, \Delta$를 계산하는 박막 반사 모델을 세운다. 이것은 [3편](/posts/ellipsometry-thin-film-multilayer-reflectance/)에서 한다. 그다음 이 모델을 측정값에 맞추는 비선형 최소자승 피팅으로 역산한다. 이 과정은 [최적화 방법론 3편](/posts/optimization-levenberg-marquardt/)의 Levenberg-Marquardt법이 다룬다.

## 정리 및 다음 편 예고

이번 편에서는 진공 맥스웰 방정식에서 출발해 파동 방정식, 평면파 해, 포인팅 벡터, 경계면의 스넬의 법칙과 프레넬 방정식까지 하나로 이어지는 수학적 구조를 따라갔다. 계측 관점에서 가장 중요한 대목은 포인팅 벡터의 시간 평균, 곧 세기가 전기장의 위상 정보를 지워 버린다는 사실이었다. 프레넬 반사계수 $r_p, r_s$는 복소수로서 위상을 담고 있다. 공기-유리 경계에서 두 반사파의 위상차는 브루스터각 56.3°를 지나며 180°에서 0°로 뒤집히고, Si에서는 유사 브루스터각 75.6° 근처에서 가장 빠르게 바뀐다. 엘립소메트리는 세기 대신 두 편광 성분의 비율을 재서 이 위상차를 되찾는다.

다음 편에서는 편광 상태를 정량적으로 기술하는 도구인 스토크스 벡터와 뮬러 행렬을 다룬다. 이번 편에서 $\rho = \tan\Psi\, e^{j\Delta}$로만 짚고 넘어간 편광 상태의 변화가 실제로 어떤 수학적 대상으로 표현되는지가 그 다음 단계다. 세기 하나만 주는 검출기로 $\Psi, \Delta$를 어떻게 끌어내는가 하는 장비 쪽 질문은 [타원계측기 1편](/posts/ellipsometer-oblique-rotating-element/)이 이어받는다. 이 배경이론 시리즈가 계산의 틀을 세우고, 타원계측기 시리즈가 그 계산을 실제 측정으로 옮기는 장비를 다룬다.

## 참고자료

- 김영준, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발," 서울대학교 대학원 박사학위논문, 2025, 2.1–2.2절.
- E. Hecht, *Optics*, 5th ed., Pearson, 2017, ch. 3–4 (전자기파와 경계 조건).
- M. Born and E. Wolf, *Principles of Optics*, 7th ed., Cambridge University Press, 1999, ch. 1 (프레넬 방정식의 표준 유도).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, §1.4 (시간 인자 규약), §5.3.2 (브루스터각을 이용한 편광자 방위각 보정).
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 2, 4 (프레넬 방정식은 2.3절, 브루스터각과 유사 브루스터각은 2.3.4절, 엘립소메트리 파라미터와의 연결은 4.1절).
- M. A. Green, "Self-consistent optical parameters of intrinsic silicon at 300 K including temperature coefficients," *Solar Energy Materials and Solar Cells* 92, 1305–1310 (2008) (Si의 633 nm 광학 상수).

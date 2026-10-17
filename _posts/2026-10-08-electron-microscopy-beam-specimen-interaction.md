---
layout: post
title: "전자현미경 배경이론 3편 — 전자빔을 1 nm로 좁히면 1 nm까지 보이는가"
categories: [전자현미경, 전자산란]
date: 2026-10-08 20:00:00 +0900
page_id: electron-microscopy-beam-specimen-interaction
tags: [electron-microscopy, monte-carlo, scattering, interaction-volume, backscattered-electron]
description: "빔을 1 nm로 좁혀도 신호는 1 µm에서 나온다. 산란 단면적부터 몬테카를로까지 그 간극을 따라간다."
math: true
---

[2편](/posts/electron-microscopy-electron-optics/)까지가 전자를 시편 위 한 점에 모으는 이야기였다. 셰르처(Scherzer) 정리 때문에 구면수차는 지울 수 없다. 그래서 조리개각은 회절과의 균형점(2편의 조건에서 6.5 mrad)에 묶이고, 분해능은 $C_s^{1/4}\lambda^{3/4}$에 갇힌다는 것까지 확인했다. 그렇게 애써서 빔을 1 nm로 좁혔다고 하자. 그러면 1 nm 분해능의 상을 얻는가?

얻지 못한다. 전자는 시편 표면에서 멈추지 않는다. 안으로 들어가 원자와 수백 번 부딪히며 옆으로 퍼지고, 그렇게 퍼진 영역 전체에서 신호가 나온다. 30 kV로 가속한 전자를 실리콘에 쏘면 신호가 나오는 영역은 깊이와 폭이 모두 수 µm에 이른다. 이 영역을 상호작용 부피(interaction volume)라 한다. 빔을 아무리 좁혀도 이 부피가 상의 분해능을 제한한다.

그런데 주사전자현미경(scanning electron microscope, SEM) 상이 실제로 1 nm 수준의 구조를 보여주는 것도 사실이다. 두 이야기가 어떻게 양립하는가? 답은 "신호"를 하나로 뭉뚱그린 데 있다. 이번 편에서는 전자가 시편 안에서 겪는 일을 산란 단면적에서 출발해 세우고, 몬테카를로(Monte Carlo) 계산으로 상호작용 부피를 직접 구한다. 그다음 그 부피 안에서 이차전자(secondary electron, SE)와 후방산란전자(backscattered electron, BSE)와 특성 X선(characteristic X-ray)이 각각 어느 깊이에서 나오는지를 갈라 본다. 4편에서 다룰 SEM의 대비는 전부 이 갈림에서 생긴다.

## 1. 전자는 무엇과 부딪히는가

시편 안에는 원자핵과 궤도전자가 있다. 입사 전자는 둘 모두와 부딪히지만 결과는 다르다.

원자핵과의 충돌은 탄성 산란(elastic scattering)이다. 원자핵이 전자보다 수만~수십만 배 무거우므로(실리콘 핵은 약 5만 배) 전자는 에너지를 거의 잃지 않고 방향만 크게 바뀐다. 30 keV 전자가 실리콘 핵에 한 번에 넘길 수 있는 에너지는 많아야 2.4 eV로, 가진 에너지의 1만분의 1에도 못 미친다. 당구공이 벽에 맞고 튀는 것에 가깝다. 이 산란이 전자를 옆으로 퍼뜨리고, 충분히 여러 번 겹치면 전자를 아예 되돌려 보낸다.

궤도전자와의 충돌은 비탄성 산란(inelastic scattering)이다. 질량이 같으므로 에너지를 실질적으로 주고받지만 방향은 거의 바뀌지 않는다. 이때 넘어간 에너지가 시편 안에서 이차전자를 만들고, 내각 전자를 때려내면 그 자리를 채우는 과정에서 특성 X선이 나온다. 그림1은 이 세 가지 신호가 원자 하나에서 생기는 모습을 나란히 그린 것이다.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/ext-electron-emission-mechanisms.png" alt="입사 전자가 원자와 부딪혀 이차전자, 후방산란전자, 특성 X선을 내는 세 과정" width="380">
_그림1. 입사 전자(PE)가 원자 하나와 부딪혀 신호를 내는 세 방식. 위는 바깥 궤도전자를 떼어 내는 이차전자(SE), 가운데는 원자핵 곁에서 크게 휘어 되돌아 나오는 후방산란전자(BSE), 아래는 내각의 빈자리를 바깥 전자가 채우며 내는 특성 X선(Characteristic X-ray)이다. (출처: Rob Hurt, [Electron emission mechanisms](https://commons.wikimedia.org/wiki/File:Electron_emission_mechanisms.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

정리하면 방향을 바꾸는 것은 탄성 산란이고 에너지를 잃게 하는 것은 비탄성 산란이다. 전자의 궤적은 이 둘이 번갈아 작용한 결과이고, 이번 편의 계산도 이 두 조각으로 이루어진다.

## 2. 얼마나 자주, 얼마나 크게 휘는가

탄성 산란을 정량화하려면 확률이 필요하다. 원자 하나가 지나가는 전자를 휘게 할 확률을 면적으로 환산한 것이 산란 단면적(scattering cross section)이다. 차폐 러더퍼드(screened Rutherford) 모형에서는 이 면적이 다음과 같이 쓰인다.

$$ \sigma \;\propto\; \frac{Z^2}{E^2} \cdot \frac{1}{\alpha_s(1+\alpha_s)}, \qquad \alpha_s = \frac{3.4\times10^{-3}\,Z^{0.67}}{E} $$

$Z$는 원자번호, $E$는 전자의 에너지(keV)다. 차폐 상수(screening parameter) $\alpha_s$는 궤도전자가 원자핵의 전하를 가리는 정도를 나타낸다. 식 앞쪽의 $Z^2/E^2$ 한 덩어리가 이번 편에서 나올 정성적인 결론 대부분을 담고 있다.

식의 두 덩어리, 곧 $Z^2/E^2$와 $\alpha_s$를 담은 뒷부분이 소제목의 두 질문에 따로 답한다. 먼저 얼마나 크게 휘는가다. 큰 각으로 휘는 산란의 단면적은 앞 인자 그대로 $Z^2/E^2$를 따른다. 그래서 큰 각으로 휠 확률은 텅스텐($Z=74$)이 실리콘($Z=14$)보다 $Z^2$ 기준으로 28배 크다(90°를 넘게 휘는 산란만 세면 27.7배).

다음은 얼마나 자주 휘는가다. 이쪽은 작은 각까지 모두 센 총단면적이 답한다. $\alpha_s$에도 $Z^{0.67}/E$가 들어 있어, $\alpha_s \ll 1$인 이 에너지 영역에서는 앞 인자의 일부가 상쇄된다. 총단면적은 대략 $Z^{4/3}/E$를 따르고, 산란 빈도로 따지면 텅스텐 원자 하나가 실리콘의 약 10배다. 즉, 무거운 원소에서 전자는 더 자주 휘기도 하지만 무엇보다 한 번에 크게 휜다. 분모의 $E^2$은 전자가 느릴수록 잘 휜다는 뜻이다. 가속전압을 낮추면 전자가 깊이 들어가기 전에 방향을 잃는다.

단면적을 알면 평균자유행로(mean free path)가 나온다. 단위 부피당 원자 수가 $N$이고, $A$는 원자량, $N_A$는 아보가드로 수, $\rho$는 밀도일 때

$$ \ell_{el} = \frac{1}{N\sigma} = \frac{A}{N_A\,\rho\,\sigma} $$

이고, 다음 산란까지의 거리는 이 값을 평균으로 하는 지수분포에서 뽑으면 된다. 예를 들어 15 keV 전자는 실리콘에서 평균 18 nm마다 한 번 꺾이고, 텅스텐에서는 1.6 nm마다 꺾인다.

에너지 손실 쪽은 베테(Bethe)의 저지능(stopping power) 식을 쓴다. 저지능은 전자가 단위 길이를 지날 때 잃는 에너지다.

$$ \frac{dE}{ds} = -78500\,\frac{\rho\,Z}{A\,E}\,\ln\!\left(\frac{1.166\,E}{J}\right) \quad [\mathrm{keV/cm}] $$

$J$는 평균 이온화 퍼텐셜(mean ionization potential)이고, 실리콘에서 172 eV다. 이 식은 $E$가 작아지면 로그 안이 1 아래로 내려가 부호가 뒤집힌다. 그래서 저에너지에서는 로그 안을 $1.166(E+0.85J)/J$로 바꾼 조이-루오(Joy-Luo) 수정을 쓴다. 아래 계산도 그렇게 했다. 이 식으로 셈하면 15 keV 전자는 실리콘 1 µm를 지날 때 약 2.8 keV를 잃는다. 즉, 15 keV를 다 쓰기까지 µm 단위의 길을 가고, 그동안 평균 320번쯤 꺾인다.

## 3. 난수로 전자의 궤적을 어떻게 만드는가

두 식이 있으면 전자 한 개의 궤적을 만들 수 있다. 해석적으로 풀 필요가 없다. 난수를 써서 하나씩 따라가면 된다. 규칙은 네 가지뿐이다.

```python
# monte_carlo.py 핵심부 (전체 코드: _code/electron-microscopy-beam-specimen-interaction/)
lam = mean_free_path(mat, e)                     # 1. 평균자유행로
step = -lam * np.log(rng.random())               # 2. 지수분포에서 이동 거리
loss = -stopping_power(mat, e) * step            # 3. Bethe 식으로 에너지 손실
e -= loss
a = screening(mat, e)                            # 4. 차폐 러더퍼드 분포에서 산란각
cos_t = 1.0 - 2.0 * a * r / (1.0 + a - r)        #    (r 은 균일 난수)
direction = _scatter_direction(direction, cos_t, 2 * np.pi * rng.random())
```

네 단계가 전부다. 이 과정을 전자가 표면 밖으로 나가거나 에너지가 0.5 keV 아래로 떨어질 때까지 되풀이한다. 표면을 넘어 되돌아 나오면 후방산란전자로 집계한다.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/fig1-interaction-volume.png" alt="실리콘과 텅스텐에서 5 kV와 30 kV 전자의 몬테카를로 궤적과 상호작용 부피" width="820">
_그림2. 물질과 가속전압에 따른 상호작용 부피. 패널마다 전자 220개의 궤적을 그렸고, 파란 궤적은 시편 안에서 멈춘 전자, 붉은 궤적은 되돌아 나온 후방산란전자다. 점선 반원은 카나야-오카야마 침투 깊이이고, η는 전자 3,000개로 따로 센 후방산란계수다. 네 패널의 축척이 서로 다르며, 가속전압을 올리면 부피가 커지고 텅스텐에서는 표면에 붙은 반구가 된다._

네 패널이 2절의 두 인자를 그대로 보여준다.

가속전압을 5 kV에서 30 kV로 올리면(왼쪽에서 오른쪽) 부피가 통째로 커진다. 실리콘에서 침투 깊이(penetration depth)가 0.47 µm에서 9.31 µm로 20배 늘어난다. 에너지가 크면 단면적이 작아 덜 휘고, 잃을 에너지도 많으니 멀리 간다.

원자번호를 올리면(위에서 아래) 부피가 작아지고 모양이 바뀐다. 30 kV에서 실리콘은 9.31 µm인데 텅스텐은 1.68 µm로 5분의 1이다. 다만 이 차이의 대부분은 밀도(2.33 대 19.25 g/cm³, 8.3배)에서 온다. 이 모형에서 밀도는 평균자유행로와 저지능에 똑같이 들어가 모든 길이를 $1/\rho$로 줄인다. 그래서 텅스텐의 밀도를 실리콘과 같게 두면 최대 깊이는 실리콘보다 16 % 얕을 뿐이다. 원자번호가 직접 바꾸는 것은 크기보다 모양이다. 실리콘의 부피는 표면 아래에서 옆으로 부푸는 물방울 모양인 반면, 텅스텐은 표면에 바짝 붙은 반구에 가깝다. 실리콘에서는 전자가 한동안 직진하다가 뒤늦게 방향을 잃지만, 텅스텐에서는 들어가자마자 큰 각도로 휘기 때문이다. $Z^2$의 결과다.

그림에 그린 점선 반원은 카나야-오카야마(Kanaya-Okayama) 침투 깊이다.

$$ R_{KO} = \frac{0.0276\,A\,E_0^{1.67}}{Z^{0.89}\,\rho} \quad [\mu\mathrm{m}] $$

여기서 $E_0$는 입사 에너지(keV)다. 궤적은 대부분 이 반원의 안쪽 절반에 몰려 있다. $R_{KO}$는 가장 멀리 가는 소수의 전자까지 포함한 최대 도달 거리의 경험적 추정이지 전형적인 값이 아니다. 실리콘 15 kV에서 에너지의 90 %가 쌓이는 깊이는 1.55 µm로 $R_{KO}$(2.93 µm)의 53 %이고, 30 kV에서도 5.45 µm로 59 %다. 즉, $R_{KO}$의 절반 남짓한 깊이 안에 에너지의 90 %가 쌓이고, 신호의 대부분도 거기서 나온다.

## 4. 신호마다 어느 깊이에서 나오는가

상호작용 부피가 µm 규모라는 것은 알았다. 그런데 SEM은 어떻게 nm 구조를 보여주는가?

부피 전체에서 나온 신호를 다 쓰는 것이 아니기 때문이다. 신호마다 시편 밖으로 빠져나올 수 있는 깊이가 다르고, 그 차이가 세 자릿수에 이른다.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/fig2-depth-distribution.png" alt="실리콘 15 kV에서 에너지 침착, 특성 X선 발생, 후방산란전자의 도달 깊이 분포와 이차전자 탈출 깊이" width="820">
_그림3. 실리콘 15 kV에서 전자 1,200개로 센 신호별 깊이 분포. 왼쪽의 파란 실선은 에너지 침착, 초록 파선은 특성 X선 발생, 붉은 점선은 후방산란전자가 되돌아서기 전에 도달한 최대 깊이이고, 주황 띠는 이차전자 탈출 깊이(5 nm)다. 오른쪽은 표면 20 nm를 확대한 것으로, 왼쪽 축척에서는 이차전자 탈출층이 선 하나로 보일 만큼 얇다._

왼쪽 그림의 파란 실선은 깊이별 에너지 침착(energy deposition)이다. 이차전자는 비탄성 산란으로 넘어간 에너지가 만드는 것이므로 생성량이 이 곡선을 따른다. 0.5 µm 근처에서 최대가 되고 2.5 µm까지 이어진다.

초록 파선은 특성 X선을 만들 수 있는 조건, 즉 전자가 아직 실리콘 K 껍질(K shell)의 이온화 에너지 1.84 keV보다 큰 에너지를 지닌 구간만 추린 것이다. 파란 실선과 거의 겹친다. 15 keV로 들어온 전자는 궤적의 대부분에서 1.84 keV보다 훨씬 큰 에너지를 지니므로, X선은 상호작용 부피 거의 전체에서 나온다는 뜻이다. 에너지 분산 X선 분광(energy dispersive X-ray spectroscopy, EDS)으로 조성을 잴 때 공간 분해능이 µm 수준에 머무는 이유가 여기에 있다. 거꾸로 말하면 EDS의 공간 분해능을 높이려면 가속전압을 낮추면 된다. 실리콘에서 1.84 keV 넘는 전자가 쌓는 에너지의 90 %가 들어가는 깊이는 15 kV에서 1.51 µm, 10 kV에서 0.74 µm, 5 kV에서 0.20 µm다. 물론 반대급부가 있다. 가속전압이 이온화 에너지에 가까워질수록 X선을 만들 수 있는 구간이 짧아져 신호가 약해진다.

붉은 점선은 되돌아 나온 전자들이 되돌아서기 전에 도달했던 최대 깊이의 분포다. 0.3~0.5 µm에서 최대가 되고 1 µm를 넘는 것은 거의 없다. 후방산란전자는 상호작용 부피의 위쪽 3분의 1 정도만 훑고 나온다.

그리고 주황색 띠가 이차전자의 탈출 깊이(escape depth)다. 이차전자의 에너지는 50 eV 아래로 매우 낮아서 시편 안에서 금방 다시 산란되어 멈춘다. 금속과 반도체에서 실제로 표면 밖으로 나올 수 있는 것은 표면에서 수 nm 안쪽에서 생긴 것뿐이다(절연체는 수십 nm). 대표값 5 nm는 이 조건에서 침투 깊이 2.93 µm의 약 0.2 %, 곧 600분의 1에 지나지 않는다. 왼쪽 그림에서 이 띠가 세로선 하나로 보이는 이유이고, 오른쪽에 20 nm 구간을 따로 확대한 이유이기도 하다.

그림3에 그리지 않은 신호가 하나 더 있다. 오제 전자(Auger electron)다. 입사 전자가 내각 전자를 때려내면 빈자리를 바깥 궤도의 전자가 채운다. 이때 남는 에너지가 X선으로 나가는 대신 다른 궤도전자를 밖으로 떼어 내기도 한다(그림4). 그렇게 튀어나온 것이 오제 전자이고, 에너지가 원소마다 정해져 있어 조성 정보를 싣는다. 탈출 깊이는 수 nm 이하이고, 에너지가 낮은 오제선은 1 nm 안팎이다(수십~수백 eV 전자의 평균자유행로는 [XPS 기초 2편](/posts/xps-surface-sensitivity-imfp/)의 그림1에서 읽을 수 있다). 전자빔으로 최상층 원자층의 조성만 골라 보는 오제 전자 분광(Auger electron spectroscopy, AES)이 성립하는 근거가 이것이고, 같은 이유로 표면이 조금만 오염되어도 신호가 묻힌다. 전자 대신 X선으로 때려서 같은 종류의 문제를 푸는 쪽이 [X선 광전자 분광(X-ray photoelectron spectroscopy, XPS)](/posts/xps-photoemission-binding-energy/)이다.
<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/ext-auger-process.png" alt="오제 과정을 원자 궤도와 에너지 준위로 그린 도식" width="640">
_그림4. 오제 과정. (a) 입사 전자(빨간 화살표)가 1s 전자를 때려낸 뒤(Electron collision), 2s 전자가 그 빈자리로 내려가고 남는 에너지를 받은 2p 전자가 오제 전자로 튀어나온다(Auger electron emission). (b) 같은 과정을 에너지 준위로 그리면 K 껍질의 빈자리를 L₁ 전자가 채우고, L₂,₃ 전자가 진공 준위(Vac) 위로 나간다. (출처: A. Carlson, [Auger Process](https://commons.wikimedia.org/wiki/File:Auger_Process.svg), 퍼블릭 도메인)_

처음의 모순은 이렇게 풀린다. 상호작용 부피는 µm지만 이차전자가 실어 나르는 정보는 표면 수 nm에서 온다. 입사점 바로 아래, 빔이 옆으로 퍼지기 전의 좁은 기둥에서 생긴 이차전자(SE1)가 빔 지름에 가까운 분해능을 준다. 후방산란전자가 표면을 빠져나가며 만든 이차전자(SE2)도 같은 얕은 층에서 나오지만, 후방산란전자가 나오는 넓은 영역에 퍼져 있어 분해능 없는 배경으로 깔린다. 그림5와 같은 몬테카를로로 실리콘 15 kV에서 전자 3,000개를 추적하면, 후방산란전자가 표면으로 나오는 자리는 입사점에서 중앙값 약 0.8 µm, 90 %가 1.4 µm 안쪽이다. 즉, SE1은 빔 지름 정도의 점에서, SE2는 지름 약 3 µm의 원에서 나온다. 세기로는 SE2가 SE1의 $\beta\eta$배쯤 되는데, Seiler에 따르면 $\beta$는 2~3이고 실리콘의 $\eta$는 0.17이니 SE2는 SE1의 3~5할이다. SEM이 지형을 nm 수준으로 보여주는 것은 상호작용 부피를 이겨서가 아니라 그 부피 중 아주 얇은 껍질만 보기 때문이다.

## 5. 무거운 원소는 왜 후방산란전자 영상에서 밝은가

이차전자가 표면 형상을 보여준다면 후방산란전자는 무엇을 보여주는가? 2절의 $Z^2$이 여기서 다시 나온다.

입사한 전자 중 되돌아 나오는 비율을 후방산란계수(backscattering coefficient) $\eta$라 한다. 큰 각 산란이 $Z^2$을 따르니 $\eta$도 원자번호를 따라 올라갈 것이라 예상할 수 있다. 20 kV에서 원소 여섯 개를 몬테카를로로 계산해 확인해 보자.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/fig3-backscatter-vs-z.png" alt="원자번호에 따른 후방산란계수를 몬테카를로 결과와 Reuter 경험식으로 비교" width="750">
_그림5. 원자번호에 따른 후방산란계수. 점은 20 kV에서 원소마다 전자 3,000개를 추적한 몬테카를로 결과이고, 회색 파선은 로이터 경험식이다. C와 Si에서는 둘이 잘 맞고, Cu부터 몬테카를로가 경험식 위로 올라가기 시작해 W에서 12 %가량 높다._

탄소($Z=6$)의 0.06에서 금($Z=79$)의 0.54까지 단조증가한다. 가벼운 원소 영역에서 가파르고 무거운 원소로 갈수록 완만해진다. 회색 파선은 실측을 맞춘 로이터(Reuter)의 경험식이다.

이 곡선이 조성 대비(compositional contrast)의 근거다. 시편 위에 원자번호가 다른 두 상이 나란히 있으면 되돌아 나오는 전자 수가 다르고, BSE 검출기로 보면 무거운 원소 쪽이 밝게 보인다(그림6은 색을 입힌 영상이라 밝기 대신 색으로 갈린다). 가벼운 원소와 중간 원소에서는 기울기가 가팔라 원자번호 차이가 작아도 대비가 비교적 크다. 곡선이 눕는 무거운 원소 쪽에서는 같은 차이가 훨씬 작은 대비로 줄어든다.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/ext-gold-arsenopyrite-bse.jpg" alt="황비철석 결정에 붙어 자란 자연금의 후방산란전자 영상" width="600">
_그림6. 황비철석(FeAsS) 결정에 붙어 자란 자연금을 후방산란전자로 찍은 영상으로, 노란 부분이 금, 회청색 부분이 황비철석이다. 금($Z=79$)과 평균 원자번호가 약 27인 황비철석은 후방산란계수가 크게 달라 두 광물이 또렷이 갈린다. 색은 원작자가 입힌 것이다. (출처: Boris Lobastov, [Gold on arsenopyrite SEM image](https://commons.wikimedia.org/wiki/File:Gold_on_arsenopyrite_SEM_image.png), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0))_

다만 그림5에는 짚고 넘어갈 어긋남이 있다. 원자번호가 큰 쪽에서 몬테카를로 점들이 경험식 위로 계통적으로 올라간다. 실리콘에서는 0.168 대 0.164로 잘 맞지만, 텅스텐에서는 0.532 대 0.476으로 12 %가량 높다. 통계 요동은 아니다. 차이(0.056)가 표준편차(약 0.009)의 6배다. 먼저 의심할 곳은 단면적이다. 차폐 러더퍼드 단면적은 첫째 보른 근사(핵전하를 약한 섭동으로 다룬다)에 궤도전자의 차폐를 경험 상수 하나로 덧붙인 것이다. 핵전하가 큰 무거운 원소일수록 이 근사가 나빠진다. 하지만 이 계산만으로는 원인이 단면적인지, 저지능 식이나 0.5 keV에서 추적을 멈춘 것 같은 다른 조각인지 가를 수 없었다. 단면적을 확인하려면 상대론적 파동 계산에서 얻는 모트(Mott) 단면적을 표로 넣어 다시 돌려 보면 된다.

이 어긋남을 남겨두는 것이 더 정직하다고 본다. 몬테카를로가 보여주는 것은 경향이지 측정값이 아니다. 상호작용 부피의 모양이 왜 그렇게 생겼는지, 신호가 어느 깊이에서 나오는지를 이해하는 데는 이 정도 모형으로 충분하지만, 정량 분석의 보정 계수로 쓰려면 단면적부터 모형의 조각들을 하나씩 검증해야 한다.

## 정리 및 다음 편 예고

전자가 시편 안에서 겪는 일은 두 갈래였다. 원자핵과의 탄성 산란이 방향을 바꾸고 궤도전자와의 비탄성 산란이 에너지를 빼앗는다. 탄성 산란, 특히 큰 각 산란의 단면적이 $Z^2/E^2$에 비례한다는 사실에서 상호작용 부피의 모양이 대부분 따라 나온다. 가속전압을 올리면 부피가 커지고, 원자번호가 큰 원소는 대개 밀도도 커서 부피가 작아지며, 원자번호 자체는 물방울 모양을 반구로 바꾼다.

그 부피가 µm 규모인데도 SEM이 nm 구조를 보여주는 이유는 신호마다 탈출 깊이가 다르기 때문이었다. 특성 X선은 부피 거의 전체에서, 후방산란전자는 위쪽 3분의 1에서, 이차전자는 표면 5 nm 안쪽(절연체는 수십 nm)에서, 오제 전자는 수 nm 이하에서만 나온다. 이차전자의 5 nm는 침투 깊이의 약 0.2 %에 불과하다. 같은 빔으로 같은 지점을 때려도 어느 검출기로 받느냐에 따라 보는 깊이가 세 자릿수 달라진다.

다음 편은 그 검출기 이야기다. 주사 방식이 어떻게 상을 만드는지부터 세우고, 이차전자 검출기와 후방산란전자 검출기가 경통 안 어디에 어떤 각도로 놓이는지가 어떻게 지형 대비(topographic contrast)와 조성 대비를 각각 만들어내는지를 본다. 표면이 기울어진 곳에서 이차전자 수율이 올라가 윤곽이 밝게 보이는 엣지 효과(edge effect), 그리고 절연체를 볼 때 부딪히는 대전(charging) 문제와 저전압 SEM이 그 해법이 되는 과정도 함께 다룬다.

## 참고자료

- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 1–5 (전자-시편 상호작용, 상호작용 부피, 후방산란계수).
- D. C. Joy, *Monte Carlo Modeling for Electron Microscopy and Microanalysis*, Oxford University Press, 1995 (단일산란 몬테카를로의 표준 구현).
- D. C. Joy and S. Luo, "An empirical stopping power relationship for low-energy electrons," *Scanning*, vol. 11, pp. 176–180, 1989 (저에너지 Bethe 수정).
- K. Kanaya and S. Okayama, "Penetration and energy-loss theory of electrons in solid targets," *Journal of Physics D*, vol. 5, pp. 43–58, 1972 (침투 깊이 경험식).
- W. Reuter, in *Proc. 6th Int. Conf. on X-Ray Optics and Microanalysis*, G. Shinoda, K. Kohra, T. Ichinokawa (eds.), University of Tokyo Press, 1972, p. 121 (후방산란계수 경험식, Goldstein 4판 식 2.2).
- L. Reimer, *Scanning Electron Microscopy: Physics of Image Formation and Microanalysis*, 2nd ed., Springer, 1998, ch. 3–4 (산란 이론과 신호 발생 깊이).
- H. Seiler, "Secondary electron emission," in *Electron Beam Interactions with Solids*, SEM Inc., 1984, pp. 33–42 (이차전자 에너지 경계 50 eV, 탈출 깊이 T ≈ 5λ).
- [Scanning Electron Microscopy, Carleton SERC](https://serc.carleton.edu/research_education/geochemsheets/techniques/SEM.html) (SEM 신호 종류의 정성적 개요).

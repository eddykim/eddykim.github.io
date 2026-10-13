---
layout: post
title: "전자현미경 배경이론 7편 — 전자현미경으로 원자의 종류까지 가려낼 수 있는가"
categories: [전자현미경, 분석전자현미경]
date: 2026-11-05 20:00:00 +0900
page_id: electron-microscopy-stem-analytical
tags: [electron-microscopy, stem, haadf, eels, 4d-stem]
description: "빔을 다시 한 점에 모아 훑는다. 검출기 내각 하나가 원자번호 의존성을 정하고, 6편의 부호 반전 문제가 사라진다."
math: true
---

[6편](/posts/electron-microscopy-tem-phase-contrast-ctf/)은 좋은 소식과 나쁜 소식을 함께 남겼다. 원자 배열이 상에 직접 나타난다는 것이 좋은 소식이고, 그 상의 밝은 점이 원자라는 보장이 없다는 것이 나쁜 소식이었다. 디포커스만 바꿔도 명암이 뒤집히므로 상을 해석하려면 계산이 필요했다.

이번 편에서는 빔을 다시 한 점에 모은다. [4편](/posts/electron-microscopy-sem-imaging/)의 주사전자현미경(scanning electron microscope, SEM)처럼 시편 위를 훑되, 되돌아 나오는 신호가 아니라 **뚫고 나간** 신호를 받는다. 주사투과전자현미경(scanning transmission electron microscope, STEM)이다.

이렇게 하면 6편의 문제가 상당 부분 사라진다. 상을 만드는 방식이 간섭이 아니라 산란된 전자를 세는 것으로 바뀌기 때문이다. 세는 것에는 부호가 없다. 그 대신 어느 각도의 전자를 셀 것인가라는 새로운 선택이 생기고, 그 선택 하나가 무엇을 보게 될지를 정한다.

이번 편은 시리즈의 마지막이다. STEM의 원리에서 출발해 분석 기법과 4D-STEM까지 훑고, 시편 준비와 반도체 계측 응용으로 마무리한다.

## 1. 빔을 모아 훑는 STEM(Scanning Transmission Electron Microscope)은 왜 투과전자현미경과 같은 상을 주는가

빔을 한 점에 모아 훑는 STEM과 넓게 비추는 투과전자현미경(transmission electron microscope, TEM)은 완전히 다른 장비처럼 보인다. 그런데 왜 비슷한 상이 나오는가?

상반성 원리(principle of reciprocity) 때문이다. 에너지를 잃지 않는 탄성 산란이라면 전자 파동의 전파에는 시간 역전 대칭이 있어, 광선의 진행 방향을 모두 뒤집어도 같은 경로가 성립한다. 광원이 있던 자리에 검출기를 두고 검출기가 있던 자리에 광원을 두면 결과가 같다는 뜻이다.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/fig1-reciprocity.png" alt="TEM과 STEM의 광선도를 상반성 원리에 따라 나란히 비교" width="800">
_그림1. TEM(왼쪽)과 STEM(오른쪽)의 광선도. 그려진 선은 같고 화살표 방향과 광원·검출기의 자리만 바뀐다. 그래서 TEM의 조명 각은 STEM의 수집 각에, TEM의 수집 각은 STEM의 프로브 수렴 각에 대응한다._

그림1의 두 광선도는 그려진 선이 완전히 같다. 다른 것은 어느 끝을 광원이라 부르고 어느 끝을 검출기라 부르느냐뿐이다. STEM은 관례대로 뒤집어 그렸다. 광선은 시편의 한 점에서 만난다. STEM으로 읽으면 아래 렌즈가 넓은 빔을 시편 위 한 점으로 모으는 집속 프로브이고, TEM으로 읽으면 시편의 한 점에서 나온 산란파를 대물렌즈가 받아 내는 그림이다.

이 대응에는 구체적인 내용이 있다. TEM에서 시편을 비추는 조명 각은 STEM에서 검출기가 받는 수집 각에 대응하고, TEM의 대물렌즈 수집 각은 STEM의 프로브 수렴 각에 대응한다. 그래서 STEM에서 작은 명시야 검출기로 얻은 상은 TEM 상과 같아진다.

대응이 성립한다면 굳이 STEM을 쓸 이유는 무엇인가? 두 가지다. 첫째, 검출기를 여러 개 놓고 각도별로 나눠 받을 수 있다. 둘째, 빔이 한 점에 모여 있으므로 그 점에서 나오는 X선이나 에너지를 잃은 전자를 동시에 받아 조성을 잴 수 있다. 상과 분석을 같은 주사에서 동시에 얻는 것이 STEM의 실질적인 이점이다.

## 2. HAADF(High-Angle Annular Dark Field) 상의 밝기는 원자번호를 얼마나 따르는가

각도별로 나눠 받는다면, 어느 각도를 받아야 무엇이 보이는가?

결론부터 말하면, 높은 각도로 산란된 전자만 고리 모양 검출기로 모아 받으면 상의 밝기가 그 자리 원자의 원자번호를 거의 그대로 따른다. 6편과 달리 밝은 점이 원자다. 그 이유와 "거의"의 정도는 3편에서 쓴 차폐 러더퍼드 모형으로 따져 볼 수 있다. 그 미분 단면적은 다음과 같고, 3편의 총단면적은 이것을 모든 방향으로 적분한 값이다($E$는 keV 단위).

$$ \frac{d\sigma}{d\Omega} \;\propto\; \frac{Z^2}{\left(\sin^2(\theta/2) + \alpha_s\right)^2}, \qquad \alpha_s = \frac{3.4\times10^{-3}\,Z^{0.67}}{E} $$

여기서 $\theta$는 산란각, $Z$는 원자번호, $\alpha_s$는 궤도전자가 원자핵의 전하를 가리는 정도를 나타내는 차폐 상수다. 내각 $\theta_1$, 외각 $\theta_2$인 환형 검출기가 받는 양은 이것을 입체각으로 적분하면 되고, 닫힌 형태로 나온다.

$$ \sigma \;\propto\; Z^2\left[\frac{1}{\sin^2(\theta_1/2)+\alpha_s} - \frac{1}{\sin^2(\theta_2/2)+\alpha_s}\right] $$

즉, 검출기가 받는 양은 내각에서 센 값과 바깥각에서 센 값의 차이다. $Z^2$이 앞에 붙어 있으니 원자번호 제곱에 비례할 것 같다. 그런데 차폐 상수 $\alpha_s$ 자체도 $Z^{0.67}$에 비례한다. 무거운 원소일수록 차폐각이 커서 대괄호 안의 값이 작아지고, 그만큼 $Z^2$이 깎인다.

깎이는 정도는 내각이 정한다. 내각이 차폐각보다 훨씬 크면 $\alpha_s$가 무시되어 $Z^2$이 살아남고, 내각이 차폐각 근처면 $\alpha_s$가 크게 작용해 지수가 내려간다.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/fig2-z-exponent.png" alt="검출기 내각에 따른 원자번호 유효 지수와 원소별 단면적" width="820">
_그림2. 3편의 차폐 러더퍼드 단면적을 환형 검출기 범위로 적분한 결과(200 keV, 바깥각 200 mrad). 왼쪽은 내각에 따른 유효 지수 n으로, C·Si·Cu·Ag·W·Au 여섯 원소에 직선을 맞춰 얻었다. 배경 띠는 명시야(BF), 환형 명시야(ABF), 환형 암시야(ADF), HAADF 검출기 구역을 나타내는 도식이고, 점선은 Si(20 mrad)와 Au(36 mrad)의 차폐각이다. 오른쪽은 내각 20 mrad(n = 1.60)와 100 mrad(n = 1.95)에서 계산한 원소별 단면적이다. 내각이 차폐각을 넘으면 지수가 2에 다가간다._

그림2의 왼쪽이 그 계산이다. 내각 20 mrad에서 유효 지수가 1.60이고, 80 mrad에서 1.93, 150 mrad에서 1.97까지 올라간다. 교과서와 실측이 HAADF 대비의 지수를 2보다 작게(흔히 1.6~1.9) 적는 이유 가운데 하나가 이 차폐 효과이고, 3편의 차폐 상수만으로 그 경향이 나온다. 다만 고립 원자의 탄성 산란만 본 모형이라, 실제 지수는 열확산산란(thermal diffuse scattering, TDS)과 채널링(channelling), 비교하는 원소 조합에 따라 더 갈린다.

점선으로 표시한 차폐각 $\theta_0 = 2\sqrt{\alpha_s}$는 차폐가 산란을 누그러뜨리기 시작하는 각으로, 교과서의 TEM용 근사식과 계수까지 같은 비상대론 근사다. 200 keV에서 실리콘은 20 mrad, 금은 36 mrad다. 내각이 이 값들을 넘어가야 지수가 2에 가까워진다는 것이 곡선에서 보인다. 고각 환형 암시야(high-angle annular dark field, HAADF) 검출기가 보통 50 mrad 이상을 쓰는 데에는 이 점도 작용한다. 내각을 그만큼 키우는 일차 목적은 브래그 회절과 간섭성 산란을 걸러 상을 비간섭성으로 만드는 것이다.

오른쪽 그림은 같은 계산을 원소별로 본 것이다. 로그-로그에서 직선의 기울기가 곧 지수이고, 내각 20 mrad와 100 mrad의 두 직선이 눈에 띄게 다른 기울기를 갖는다.

지수만으로는 상이 어떻게 보일지 감이 오지 않는다. 같은 계산을 원자 하나당 밝기의 비로 바꿔 보자. 내각 80 mrad에서 하프늄 원자는 실리콘 원자보다 약 23배($Z^2$이면 26배), 금은 약 28배($Z^2$이면 32배), 저마늄은 약 5.0배 밝고, 산소는 실리콘의 0.33배에 그친다. 내각을 20 mrad로 낮추면 금과 실리콘의 비는 약 15배로 줄어든다. 즉, 실리콘 위의 산화하프늄 게이트 절연막은 환하게 빛나고, 실리콘과 실리콘저마늄이 번갈아 쌓인 층은 갈려 보이며, 산소는 거의 보이지 않는다. 가벼운 원소를 보려면 다른 신호가 필요하다는 뜻이고, 그 신호를 3절에서 다룬다. 이 값들도 고립 원자 모형의 것이라 실제 비는 시편 두께와 결정 방향에 따라 달라진다.

여기서 6편과의 차이가 드러난다. HAADF 상의 밝기는 그 자리에 있는 원자의 원자번호와 개수를 거의 단조롭게 반영한다. $\sin\chi$처럼 부호가 뒤집히지 않으므로 밝은 점이 원자다. 초점을 바꿔도 대비가 반전되지 않는다. Z-대비(Z-contrast)라 부르는 이유이고, 조성이 다른 층이 겹친 구조를 볼 때 HAADF를 먼저 찍는 이유다. 다만 시편이 두꺼우면 밝기가 두께에 비례하지 않고, 원자 기둥을 따라 전자가 모여 흐르는 채널링의 영향도 받으므로, 밝기를 원자 수로 바로 옮길 때는 시편 두께를 함께 따져야 한다. 그림3은 실제 원자 분해능 상에서 이 차이를 보여 준다.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/ext-srtio3-adf-abf.jpg" alt="SrTiO3의 원자 분해능 환형 암시야 상과 환형 명시야 상 비교" width="600">
_그림3. 페로브스카이트 산화물 SrTiO₃를 [110] 방향에서 찍은 원자 분해능 STEM 상. 왼쪽(Annular Dark Field)은 환형 암시야 상으로 무거운 스트론튬 기둥이 가장 밝고 산소는 보이지 않으며, 오른쪽(Annular Bright Field)은 환형 명시야 상으로 산소 기둥까지 어두운 점으로 드러난다. 겹쳐 그린 원은 스트론튬(초록), 티타늄(회색), 산소(빨강)의 자리이고, 원본에 검출기 각도는 적혀 있지 않다. (출처: Magnunor, [Scanning transmission electron microscopy srtio3 compare adf abf](https://commons.wikimedia.org/wiki/File:Scanning_transmission_electron_microscopy_srtio3_compare_adf_abf.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

고각 산란의 주성분은 순수한 탄성 러더퍼드 산란이 아니라 열로 흔들리는 원자와의 산란, 즉 열확산산란이다. 그 덕분에 신호가 결정 방향에 덜 민감해져 회절 대비가 섞이지 않고, 결정 방향이 조금만 바뀌어도 회절파 세기가 크게 바뀌는 문제에서도 비교적 자유롭다. 5편에서 본 회절 대비는 이 점에서 정반대였다. 편위 파라미터 $s$가 회절파 세기를 정하는데, 두께 90 nm에서는 $s$가 0.009 nm⁻¹만 벗어나도 세기가 절반으로 줄었다. 물론 HAADF 신호도 채널링과 격자 변형에는 여전히 민감하다.

## 3. EDS(Energy Dispersive X-ray Spectroscopy)와 EELS(Electron Energy Loss Spectroscopy)는 무엇을 읽는가

빔이 한 점에 머무는 동안, 그 점이 무엇으로 이루어졌는지도 알 수 있는가? 3편에서 본 특성 X선과 비탄성 산란이 여기서 쓰인다.

에너지 분산 X선 분광(energy dispersive X-ray spectroscopy, EDS)은 시편 위쪽에 놓은 검출기로 특성 X선의 에너지를 잰다. 원소마다 정해진 에너지가 나오므로 조성을 알 수 있다. 3편에서 X선은 상호작용 부피 거의 전체에서 나온다고 했다. 그런데 TEM용 박편은 두께가 100 nm도 안 되므로 그 부피 자체가 작다. SEM에서 µm 수준이던 공간 분해능이 nm 수준으로 좋아지는 것이 박편화의 이득이다. 그림4는 HAADF 상과 EDS를 같은 자리에서 함께 얻은 예다.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/ext-stem-haadf-eds-cacarbonate.jpg" alt="남세균 세포 속 알갱이의 HAADF 상, EDS 원소 지도, EDS 스펙트럼" width="760">
_그림4. 남세균(cyanobacteria) 세포 속 무기물 알갱이를 STEM으로 분석한 결과. (a)(b)는 HAADF 상, (c)는 같은 영역의 EDS 원소 지도(파랑 탄소, 빨강 칼슘, 초록 인), (d)는 (b)의 빨간 원과 초록 원 자리에서 뽑은 EDS 스펙트럼(Counts는 계수, Energy는 X선 에너지)이다. 탄소가 주성분인 세포보다 칼슘과 인을 품은 알갱이가 HAADF에서 훨씬 밝고, 스펙트럼의 Ca Kα(3.7 keV)와 P Kα(2.0 keV) 봉우리가 두 알갱이를 탄산칼슘과 폴리인산염으로 가른다. (출처: Jinhua Li 외, [Minerals-06-00010-g002](https://commons.wikimedia.org/wiki/File:Minerals-06-00010-g002.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), 일부 잘라냄)_

전자 에너지 손실 분광(electron energy loss spectroscopy, EELS)은 시편을 지나면서 에너지를 잃은 전자를 분광기로 분산시켜 잰다. 잃은 에너지가 어느 껍질(shell)을 때렸는지를 알려주므로 조성이 나오고, 흡수단 바로 뒤의 미세한 구조는 결합 상태와 배위까지 알려준다. 가벼운 원소에 특히 강해서 산소나 질소의 분포, 산화 상태를 보는 데 쓰인다.

둘은 경쟁 관계라기보다 상보적이다. 차이를 표로 모으면 다음과 같다.

| 항목 | EDS | EELS | 이유 |
|---|---|---|---|
| 재는 것 | 시편에서 나온 특성 X선의 에너지 | 시편을 지난 전자가 잃은 에너지 | 같은 이온화 사건을 X선 쪽에서 보느냐 전자 쪽에서 보느냐의 차이다 |
| 강한 원소 | 무거운 원소 | 가벼운 원소(산소, 질소 등) | 가벼운 원소는 들뜬 뒤 X선보다 오제 전자를 주로 내 X선 수율이 낮다. EELS는 X선이 나오든 말든 전자가 잃은 에너지를 직접 잰다 |
| 에너지 분해능 | 넓어서 이웃 봉우리가 겹치기 쉽다 | 0.3~1 eV 수준 | EDS는 검출기 안에서 생기는 전하 수의 통계적 요동이 봉우리를 넓히고, EELS는 자기 분광기로 전자를 직접 분산시킨다 |
| 얻는 정보 | 원소 조성 | 조성에 더해 결합 상태, 배위, 산화 상태 | 흡수단 근처의 미세구조가 비어 있는 전자 상태를 반영한다 |
| 시편 두께 | 얇을수록 공간 분해능이 좋다 | 얇아야 해석된다 | 두꺼우면 EELS에서 다중 산란으로 손실이 겹친다 |

경계가 칼같지는 않다. EDS도 창이 없는 검출기를 쓰면 붕소에서 플루오린까지 정량할 수 있다. 그래도 산소의 분포나 산화 상태를 물으면 EELS를, 무거운 원소의 조성을 빠르게 훑으려면 EDS를 먼저 고르는 것이 보통이다.

분석에서 주의할 점은 분해능이 빔 지름으로 정해지지 않는다는 것이다. 빔이 시편을 지나는 동안 조금씩 퍼지고(빔 퍼짐, beam broadening), 에너지 손실이 작은 비탄성 산란은 원자 하나보다 넓은 범위에서 일어난다. 그래서 원자 하나의 조성을 읽으려면 얇은 시편과 충분한 신호가 함께 필요하다. 신호를 모으려면 시간이 걸리고, 그동안 빔이 시편을 망가뜨린다. [2편](/posts/electron-microscopy-electron-optics/)에서 본 가속전압과 손상의 저울질이 여기서도 되풀이된다.

## 4. 회절 패턴을 통째로 저장하면 무엇을 더 얻는가

지금까지의 검출기는 모두 신호를 하나의 숫자로 뭉갰다. 명시야 검출기는 가운데 원 안의 전자를 모두 더하고, HAADF는 고리 안의 전자를 모두 더한다. 각도 정보가 그 합 속에서 사라진다.

버리지 않으면 어떻게 되는가?

<img src="/assets/img/posts/electron-microscopy-stem-analytical/fig3-4dstem.png" alt="4D-STEM 데이터 구조와 거기서 뽑아내는 상들" width="850">
_그림5. 4D-STEM 데이터의 구조. 주사 위치 (x, y)마다 회절 패턴 (k_x, k_y) 한 장을 저장한다. 측정이 끝난 뒤 가운데를 더하면 명시야, 고리를 더하면 환형 암시야 상이 되고, 무게중심 이동으로 미분 위상 대비를, 패턴 전체로 타이코그래피를 얻는다._

주사 위치 $(x, y)$마다 2차원 회절 패턴 $(k_x, k_y)$를 통째로 저장하면 네 축짜리 데이터가 된다. 4차원 STEM(four-dimensional STEM, 4D-STEM)이다. 화소마다 상이 아니라 패턴 한 장이 들어 있으므로, 측정이 끝난 뒤에 어떤 상을 만들지 정할 수 있다. 크기를 어림해 보면, 256×256 위치를 훑으며 128×128 화소 검출기로 패턴을 받으면 값이 약 10억 개이고 32비트로 저장하면 4 GiB다. 네 축짜리 데이터가 측정 문제이자 계산 문제인 이유다.

그림5에서 보듯 가운데만 더하면 명시야 상이 나오고, 고리만 더하면 환형 암시야 상이 나온다. 기존 검출기가 하드웨어로 하던 일을 계산으로 대신하는 셈이다.

여기서 더 나아갈 수 있다. 회절 원판의 무게중심이 얼마나 치우쳤는지를 재면, 빔이 시편을 지나며 받은 전기장을 두께 방향으로 더한 양(투영 전기장)에 비례하는 값이 나온다. 이렇게 시편 내부의 전기장이나 자기장을 지도로 그리는 것이 미분 위상 대비(differential phase contrast, DPC)다.

가장 멀리 나간 것이 타이코그래피(ptychography)다. 이웃한 주사 위치에서 프로브가 비추는 영역이 서로 겹치도록 촘촘히 찍으면, 그 중복으로부터 시편의 위상을 역산할 수 있다. 측정한 것은 세기뿐인데 위상을 되찾는 것이므로 6편에서 위상이 세기에서 떨어져 나간다고 한 문제를 정면으로 푸는 방법이다.

역산은 최적화 문제로 푼다. 시편의 위상과 프로브의 모양을 미지수로 두고, 그것들로 계산한 회절 패턴이 측정값과 맞도록 반복해서 고쳐 나간다. [최적화 시리즈](/posts/optimization-gradient-descent/)에서 다룬 반복 최적화가 여기에 그대로 들어간다. 계산량이 큰 대신 분해능이 6편에서 본 정보한계 너머까지 오른다. 렌즈의 정보한계를 넘는 분해능은 1995년 Nellist, McCallum, Rodenburg가 처음 보였고, 1편에서 인용한 0.02 nm대의 분해능도 이 방법에서 나온 값이다. 역산 알고리즘 자체는 이 블로그의 계산 이미징과 타이코그래피 시리즈에서 따로 다룬다.

## 5. FIB(Focused Ion Beam)로 깎은 박편에는 무엇이 남는가

TEM 계열이 전부 전제하는 조건이 하나 있다. 시편이 충분히 얇아야 한다는 것이다. [5편](/posts/electron-microscopy-tem-diffraction-contrast/)의 회절 대비도, 6편의 위상물체 근사도, 이번 편의 분석도 모두 100 nm 아래를 가정한다.

그런데 실제 시편은 웨이퍼이고 소자다. 원하는 지점을 골라 그 두께로 깎아야 한다.

집속 이온빔(focused ion beam, FIB)이 그 일을 한다. 갈륨 이온을 가속해 원하는 자리를 파내고, 얇은 판 모양으로 남긴 뒤 그것을 떼어내 격자에 옮겨 붙인다. 리프트아웃(lift-out)이라 부르는 이 과정으로 관심 영역에서 두께 100 nm 아래, 고분해능 분석용이면 수십 nm의 박편을 만든다. 반도체에서 트랜지스터 하나를 골라 그 단면을 볼 수 있게 된 것은 이 기술 덕분이다. 그림6이 이렇게 만든 박편이다.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/ext-fib-lamella.jpg" alt="구리 격자에 붙인 FIB 박편과 그 박편의 원자 분해능 STEM 상" width="760">
_그림6. FIB로 만든 TEM 박편을 배율을 높여 가며 본 모습. 왼쪽 두 장은 FIB 장비의 이차전자 상으로, 구리 격자(Cu-grid)에 붙인 박편(Sample)과 그 위쪽의 박막(Thin film), 보호층(Protective layer), 기판(Substrate)을 보여 주며, 보호층은 깎기 전에 표면을 덮어 이온빔에서 지키는 층이다. 오른쪽은 같은 박편의 박막과 기판 경계를 원자 분해능 STEM으로 찍은 상이다. (출처: Magnunor, [Fib sample preparation tem](https://commons.wikimedia.org/wiki/File:Fib_sample_preparation_tem.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

대가가 있다. 갈륨 이온이 박편 표면을 때려 결정을 망가뜨리고 비정질층을 남긴다. 갈륨 자체가 시편에 박히기도 한다. 30 kV로만 깎으면 표면마다 20 nm 남짓한 비정질층이 생기는데, 박편 전체가 50 nm라면 양면을 합쳐 대부분이 손상층인 셈이다. 그래서 마무리는 2~5 kV의 낮은 에너지로 다듬어 손상층을 몇 nm로 줄인다.

이것이 전자현미경 관찰에서 가장 흔한 함정이다. 상에 보이는 비정질층이 시편에 원래 있던 것인지 FIB가 만든 것인지, 보이는 결함이 소자의 결함인지 준비 과정에서 생긴 것인지를 구분해야 한다. 상을 읽기 전에 그것이 어떻게 만들어진 시편인지를 먼저 알아야 한다는 뜻이다.

## 6. 반도체 공정에서 전자현미경은 무엇을 재는가

시리즈를 계측 쪽에서 마무리하자. 전자현미경은 반도체 공정에서 어디에 쓰이는가?

선폭 측정에는 선폭 측정 SEM(critical-dimension SEM, CD-SEM)이 쓰인다. 4편에서 본 이차전자 상의 엣지 효과를 이용해 패턴의 폭과 선폭 거칠기를 재는 전용 SEM이다. 비파괴이고 웨이퍼를 통째로 넣을 수 있어 양산 라인에 들어간다. 다만 4편에서 짚었듯 엣지에서 신호가 솟는 탓에 경계를 어디로 정할지가 알고리즘의 몫이고, 그 정의가 장비마다 다르면 값도 달라진다.

단면 구조를 봐야 하면 TEM으로 간다. 게이트 스택의 층 두께, 계면의 거칠기 같은 것은 박편을 만들어 보는 수밖에 없다(도핑 분포는 이차이온질량분석(secondary ion mass spectrometry, SIMS)이나 주사 탐침 기법이 맡는다). 파괴 검사이고 한 장을 얻는 데 시간이 오래 걸리므로, 주로 공정 개발과 불량 분석, 그리고 인라인 계측을 교정하는 참조 계측에 쓰인다.

최근에는 정량 쪽으로 더 밀고 들어간다. 게이트 올 어라운드(gate-all-around, GAA) 트랜지스터의 5 nm 실리콘 채널에서 원자층 사이 간격이 계면 근처에서 어떻게 풀리는지, 계면이 얼마나 거친지를 전자 타이코그래피로 3차원으로 잰 연구가 나와 있다. 상을 보는 것을 넘어 변형을 숫자로 원자 수준에서 읽는 단계다.

이 시리즈가 광학 계측을 다루는 블로그에서 나온 만큼, 마지막으로 역할 분담을 짚고 마친다. 엘립소메트리 같은 광학 계측은 비파괴이고 빠르며 웨이퍼 전체를 훑을 수 있지만, 측정값이 모델을 거쳐야 두께나 조성이 되고([엘립소메트리 배경이론 3편](/posts/ellipsometry-thin-film-multilayer-reflectance/)) 그 모델이 맞는지를 스스로 증명하지 못한다. 전자현미경은 느리고 파괴적이지만 구조를 직접 보여준다.

그래서 둘은 경쟁하지 않는다. 전자현미경으로 몇 개의 시편을 잘라 모델이 맞는지 확인하고, 그 모델로 광학 계측이 나머지 전부를 빠르게 잰다. 1편에서 파장과 분해능의 간극으로 시작한 이야기가 결국 여기에 닿는다.

## 정리

시리즈를 관통한 질문은 [1편](/posts/electron-microscopy-why-electrons/)의 것이었다. 200 kV 전자의 파장은 2.51 pm인데 분해능은 왜 0.2 nm인가.

[2편](/posts/electron-microscopy-electron-optics/)이 답했다. 셰르처(Scherzer) 정리 탓에 회전 대칭 전자렌즈는 구면수차를 피할 수 없고, 그래서 조리개각이 6.5 mrad 근처에 묶이며 분해능이 $C_s^{1/4}\lambda^{3/4}$에 갇힌다. 수차보정기가 그 금지를 우회했지만 네제곱근이라는 지수 탓에 $C_s$를 500배 줄여도 분해능은 4.7배 좋아지는 데 그쳤다.

[3편](/posts/electron-microscopy-beam-specimen-interaction/)부터는 전자가 시편에 닿은 뒤의 이야기였다. 큰 각으로 휘는 산란의 단면적은 $Z^2/E^2$에 비례하지만, 차폐 때문에 총단면적은 대략 $Z^{4/3}/E$에 비례한다. 이 사실에서 상호작용 부피의 성질이 따라 나왔고, 신호마다 탈출 깊이가 세 자릿수 다르다는 점이 [4편](/posts/electron-microscopy-sem-imaging/)의 SEM 대비를 갈랐다. [5편](/posts/electron-microscopy-tem-diffraction-contrast/)과 [6편](/posts/electron-microscopy-tem-phase-contrast-ctf/)에서는 시편을 얇게 깎아 투과시켰고, 회절파를 고르면 회절 대비가, 간섭시키면 위상 대비가 되었다. 그리고 위상 대비에서는 밝은 점이 원자라는 보장이 없었다.

이번 편에서 그 문제가 풀렸다. 간섭시키는 대신 산란된 전자를 각도별로 세면 부호가 사라진다. 그 각도를 정하는 검출기 내각 하나가 원자번호 의존성을 $Z^{1.6}$에서 $Z^{2}$ 사이 어디로 둘지 결정했고, 그 계산은 3편의 차폐 상수에서 그대로 나왔다. 무거운 원소는 HAADF로, 가벼운 원소와 화학 상태는 EELS로, 조성 지도는 EDS로 읽는다. 원자의 종류까지 가려낼 수 있는가라는 제목의 물음에는 이렇게 답할 수 있다. 그렇다, 다만 무엇으로 세느냐에 따라 보이는 원소가 다르다.

일곱 편에 걸쳐 되풀이된 구조가 하나 있다. 무엇을 보게 되는지는 장비가 정하는 것이 아니라 **무엇을 고르느냐**가 정한다는 것이다. 어느 검출기로 받을지(4편), 조리개를 어디에 둘지(5편), 초점을 어디에 맞출지(6편), 어느 각도의 전자를 셀지(7편). 전자현미경을 다루는 일은 그 선택의 목록을 알고, 각각이 무엇을 보이게 하고 무엇을 가리는지를 아는 일에 가깝다.

## 참고자료

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, §6.2.A (상반성), §9.4 (STEM 검출기), ch. 22 (진폭 대비), ch. 32–36 (XEDS), ch. 37–40 (EELS).
- S. J. Pennycook and P. D. Nellist (eds.), *Scanning Transmission Electron Microscopy: Imaging and Analysis*, Springer, 2011 (HAADF Z-대비와 열확산산란).
- R. F. Egerton, *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed., Springer, 2011 (EELS 표준 정리).
- P. D. Nellist, B. C. McCallum, J. M. Rodenburg, "Resolution beyond the 'information limit' in transmission electron microscopy," *Nature* 374, 630–632 (1995) (렌즈의 정보한계를 넘은 첫 타이코그래피).
- [C. Ophus, "Four-Dimensional Scanning Transmission Electron Microscopy (4D-STEM): From Scanning Nanodiffraction to Ptychography and Beyond," *Microscopy and Microanalysis* 25(3), 563–582 (2019)](https://academic.oup.com/mam/article/25/3/563/6887544) (4D-STEM 개관).
- [3D Atomic-Scale Metrology of Strain Relaxation and Roughness in Gate-All-Around Transistors via Electron Ptychography, *Nature Communications* 17, 3561 (2026)](https://arxiv.org/abs/2507.07265) (GAA 트랜지스터 채널의 원자층 간격 완화와 계면 거칠기 계측, arXiv:2507.07265).
- [J. Li et al., "Biomineralization Patterns of Intracellular Carbonatogenesis in Cyanobacteria: Molecular Hypotheses," *Minerals* 6(1), 10 (2016)](https://www.mdpi.com/2075-163X/6/1/10) (그림4의 HAADF·EDS 원 논문, CC BY 4.0).
- [TEM Metrology, Thermo Fisher Scientific](https://www.thermofisher.com/us/en/home/semiconductors/advanced-logic-devices/applications/tem-metrology.html) (반도체 단면 계측 응용).
- [ITRS Metrology 2015](https://www.semiconductors.org/wp-content/uploads/2018/06/2_2015-ITRS-2.0-Metrology.pdf) (CD-SEM 과 계측 로드맵).

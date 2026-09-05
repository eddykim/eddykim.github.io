# XPS 시리즈 5편 작성 프롬프트

> 사용법: 이 레포(`eddykim.github.io`) 루트에서 Claude Code를 열고, 아래 `---` 아래 전체를 붙여넣으세요.
> **선행 조건: 1~4편이 이미 작성되어 있어야 합니다.**

---

`write-post` 스킬 워크플로우로 아래 포스트를 작성해줘.

**단, 0~4단계는 이미 끝났다.** 아래 "확정 사항"이 그 결과다. 5단계(계획 제시)는 아래 목차를 그대로 확인만 받고 바로 6단계 초안으로 들어가라. 질의응답을 다시 하지 마라.

**작성 시작 전에 반드시 1~4편을 읽어라.** 이 편은 시리즈의 마지막 편이므로 앞 네 편의 실을 모두 회수해야 한다. 또한 엘립소메트리 1편(`_posts/2026-09-04-ellipsometry-electromagnetic-fresnel.md`)도 읽어라 — 7절의 비교가 이 블로그 전체를 잇는 대목이다.

## 확정 사항

- **시리즈**: "XPS 기초" 전 5편 중 **5편(마지막)**
- **slug**: `xps-thickness-arxps-pitfalls`
- **포스트 경로**: `_posts/2026-09-09-xps-thickness-arxps-pitfalls.md`
- **제목**: `XPS 기초 5편 — 두께 계측과 함정`
- **카테고리**: `[표면분석, 깊이분석]`
- **대상 독자**: 대학원생 + 계측/공정 실무 엔지니어
- **분량**: 9,000자 이상
- **수식**: 있음 → `math: true`
- **그림**: 4장
- **코드**: `_code/xps-thickness-arxps-pitfalls/`에 `generate_figures.py` + `overlayer.py`

### front matter

```yaml
---
layout: post
title: "XPS 기초 5편 — 두께 계측과 함정"
date: 2026-09-09 20:00:00 +0900
categories: [표면분석, 깊이분석]
tags: [xps, arxps, film-thickness, depth-profiling, beam-damage, metrology]
description: "감쇠 자체를 이용하면 박막 두께를 잴 수 있다. 오버레이어 식과 ARXPS의 역문제, 그리고 우선 스퍼터링과 빔 손상이라는 함정까지 정리한다."
math: true
---
```

## 시리즈를 관통하는 메시지

**XPS가 실제로 측정하는 것은 광전자의 운동에너지 분포 하나뿐이다. 원소·화학상태·조성비·두께는 전부 모델을 거쳐 역산한 값이며, 그 모델의 가정이 어디서 깨지는지가 이 시리즈의 축이다.**

5편은 이 축을 닫는 편이다. **두께는 XPS가 가장 간접적으로 주는 양이다 — 강도비 하나에서 모델을 통해 나오며, 균일층 가정과 $\lambda$ 추정값에 통째로 의존한다.** 마지막 절에서 시리즈 전체를 이 문장으로 마무리하라.

## 스타일 규칙 (반드시 준수)

1. **어조**: 평서문 "~이다".
2. **각 절은 답을 주기 전에 질문을 먼저 던진다.**
3. **수식 안에서 리터럴 `|` 금지.** `\lvert ... \rvert`를 써라.
4. 소제목 `## 1.` … 번호 매기기. 말미 `## 정리` (마지막 편이므로 "다음 편 예고" 대신 시리즈 전체 정리) + `## 참고자료`.
5. 기술 용어 영문 병기(첫 등장 시).
6. 30줄 넘는 코드는 발췌만, 전체는 `_code/<slug>/`.
7. 이미지 경로: `/assets/img/posts/xps-thickness-arxps-pitfalls/...`
8. 도입부는 4편 말미를 이어받아 시작.

## 목차와 각 절에서 다룰 내용

### 도입 (제목 없음)
4편까지는 시료가 깊이 방향으로 균일하다고 가정했다. 실제 시료는 거의 항상 층 구조다 — 자연산화막, 오염층, 증착막. 그런데 2편에서 본 지수 감쇠는 이 층 구조를 방해 요인으로만 만들지 않는다. 감쇠 자체가 두께에 대한 정보를 담고 있으므로, 거꾸로 두께를 재는 데 쓸 수 있다.

### ## 1. 균일 오버레이어 모델
- 질문: 산화막 위에서 금속 신호가 살아 있다면, 그 살아 있는 정도가 두께를 말해주는가?
- 기하 설정: 두께 $d$의 균일한 오버레이어(산화막) 아래에 무한히 두꺼운 기판(금속).
- 오버레이어 신호(두께 $0$에서 $d$까지 적분):
  $$ I_o = I_o^{\infty} \left[ 1 - \exp\!\left( -\frac{d}{\lambda_o \cos\theta} \right) \right] $$
- 기판 신호(오버레이어를 통과하며 감쇠):
  $$ I_m = I_m^{\infty} \exp\!\left( -\frac{d}{\lambda_m \cos\theta} \right) $$
- 두 식을 나눠 $d$에 대해 풀면 강도비 하나에서 두께가 나온다. 유도 과정을 직접 보여라.
- **가정을 명시적으로 나열하라**: (i) 오버레이어가 균일하고 평탄하다, (ii) 계면이 급격하다, (iii) $\lambda$ 값을 안다, (iv) 탄성산란을 무시한다, (v) 표면 오염층이 없거나 무시할 수 있다.

### ## 2. Strohmeier 식과 자연산화막
- 질문: 실무에서 실제로 쓰는 식은 어떤 형태인가?
- Carlson·Strohmeier 형태:
  $$ d = \lambda_o \sin\theta \, \ln\!\left( \frac{N_m \lambda_m I_o}{N_o \lambda_o I_m} + 1 \right) $$
  ($N_m, N_o$: 금속·산화물에서의 금속 원자 부피 밀도, $\theta$: **표면 평면에서 잰** take-off angle)
- **2편 5절의 각도 혼동을 여기서 회수하라**: 1절에서 $\cos\theta$(법선 기준)로 쓴 것과 이 식의 $\sin\theta$(표면 기준)는 같은 물리다. 문헌 식을 가져다 쓸 때 각도 정의를 확인하지 않으면 두께가 통째로 틀린다.
- 적용 조건: 금속과 그 산화물처럼 두 피크의 운동에너지가 비슷할 때. 대략 0~9 nm 범위에서 유효.
- 대표 사례: Si 웨이퍼 위 자연 SiO2 (1 nm 이하 수준까지 측정 가능), Al 위 Al2O3.
- 오염층 보정: 흡착 탄소층이 두 신호를 함께 감쇠시키지만 운동에너지가 다르면 감쇠 정도가 달라 오차가 된다.

### ## 3. Thickogram과 서로 다른 원소로 이루어진 층
- 질문: 오버레이어와 기판이 아예 다른 원소면 어떻게 하는가?
- Cumpson의 Thickogram: 두 피크의 강도비와 감도인자 비에서 두께를 그래프로 읽는 방법. 흡착 탄소 오염에 둔감하고 장비 인자가 상쇄되는 것이 장점. 정확도는 감쇠길이 값의 정확도에 좌우되며 대략 ±10% 수준으로 보고된다.
- Shard의 topofactor: 평면이 아닌 형상(구, 원기둥 — 코어-셸 나노입자 등)으로 확장.
- 이 절은 짧게. "균일 평면"이라는 가정이 얼마나 강한 가정인지를 보여주는 데 의의가 있다.

### ## 4. ARXPS — 각도 시리즈에서 깊이 프로파일로, 그리고 그 역문제
- 질문: 각도를 여러 개 재면 층 구조를 통째로 복원할 수 있는가?
- 2편 5절의 원리를 확장: 각도 $\theta$마다 측정한 강도는 깊이 방향 농도 분포 $c(z)$의 라플라스 변환 형태로 얻어진다:
  $$ I(\theta) \propto \int_0^{\infty} c(z) \exp\!\left( -\frac{z}{\lambda\cos\theta} \right) dz $$
- 이 적분방정식을 $c(z)$에 대해 푸는 것이 ARXPS 깊이 프로파일이다. **지수 커널을 갖는 제1종 적분방정식이므로 전형적인 불량조건(ill-posed) 역문제다** — 데이터의 작은 노이즈가 해에서 크게 증폭된다.
- 실무적 귀결: 층 개수와 순서를 미리 가정한 모델 피팅(층 모델)을 쓰거나, 정규화를 건다. 각도 데이터만으로 임의의 $c(z)$를 자유롭게 복원했다는 결과는 의심해야 한다.
- 최적화 시리즈 링크: [LM법](/posts/optimization-levenberg-marquardt/)의 damping이 사실상 정규화 역할을 한다는 점, [전역 최적화](/posts/optimization-global-heuristics/)에서 본 초기값 의존성이 층 모델 피팅에서도 그대로 나타난다는 점을 연결.
- 장점: **비파괴**다. 시료를 깎지 않고 깊이 정보를 얻는다.
- 한계: 정보깊이(약 10 nm) 밖은 못 본다.

### ## 5. 스퍼터 깊이 프로파일과 우선 스퍼터링
- 질문: 더 깊이 보려면 깎으면 되지 않는가?
- Ar$^+$ 이온빔으로 표면을 깎고 XPS를 반복 측정 → 조성 vs 깊이.
- **우선 스퍼터링(preferential sputtering)**: 구성 원소가 같은 비율로 제거되지 않는다. 가벼운 원소가 먼저 떨어져 나가는 경향이 있다. 대표 사례로 TiO2를 스퍼터하면 산소가 우선 제거되어 환원된 Ti 상태들이 나타난다 — 원래 시료에 없던 화학상태가 측정 결과에 생긴다.
- 그 밖의 아티팩트: 계면 broadening, 스퍼터 유도 혼합(atomic mixing), 표면 거칠기 증가, 깊이 축 보정 문제(스퍼터 속도는 물질마다 다르다).
- 완화: 저에너지 이온, 시료 회전, 유기물/고분자에는 클러스터 이온빔(Ar$_n^+$).
- **핵심 메시지**: 스퍼터 깊이 프로파일은 조성 추세에는 쓸 만하지만, 깎은 뒤에 읽은 화학상태를 원래 시료의 화학상태로 보고하면 안 된다.

### ## 6. 측정이 시료를 바꾼다 — 빔 손상
- 질문: X선은 비파괴 아닌가?
- 단색화된 X선이라도 손상은 일어난다. 특히 환원되기 쉬운 화학종이 취약하다.
- 대표 사례: Au(III)가 분석 중 10분 정도의 노출로 금속 Au(0)로 환원된다는 보고. Cu(II)의 환원도 잘 알려져 있다. 고분자의 결합 절단도 흔하다.
- **실무 절차**: 민감한 영역은 측정 시작과 끝에 두 번 스캔해서 변화 여부를 확인한다. 변하면 그 스펙트럼으로 화학상태를 결론짓지 않는다.
- 시리즈 축과 연결: 지금까지 "모델이 깨지는 곳"을 봤는데, 여기서는 **측정 대상 자체가 측정 도중에 바뀐다.** 가장 근본적인 형태의 모델 붕괴다.

### ## 7. 광학 계측과 XPS는 서로 무엇을 보완하는가
- 질문: 같은 SiO2/Si 두께를 엘립소메트리로도 XPS로도 잴 수 있는데, 무엇이 다른가?
- 비교 표를 하나 넣어라: 측정 원리 / 정보 깊이 / 두께 범위 / 화학상태 정보 유무 / 측정 환경(대기 vs UHV) / 속도 / 필요한 사전 정보(광학 모델 vs $\lambda$ 값과 층 모델).
- 엘립소메트리는 광학 상수를 알아야 두께가 나오고, XPS는 $\lambda$와 층 구조를 알아야 두께가 나온다 — **둘 다 모델 역산이고, 필요로 하는 사전 지식의 종류가 다를 뿐이다.**
- 엘립소메트리 1편(`/posts/ellipsometry-electromagnetic-fresnel/`)의 "세기만 재면 위상이 사라진다"와, 이 시리즈의 "운동에너지만 재면 나머지는 전부 모델이다"가 같은 문제의식이라는 것으로 마무리 문단을 써라.

### ## 정리 — 시리즈를 닫으며
- 5편 요약 + 시리즈 전체 5편을 한 문단씩 회수.
- 마지막 문장은 시리즈 축(측정량과 알고 싶은 양의 간극)으로 닫아라.

### ## 참고자료

## 그림 명세 (4장)

`_code/xps-thickness-arxps-pitfalls/generate_figures.py`로 생성. 합성 데이터임을 캡션에 명시.

- **fig1-overlayer-geometry.png** — (a) 오버레이어 기하 개념도(오버레이어 두께 $d$, 기판, 검출 각도, 탈출 경로). (b) 두께 $d$에 대한 강도비 $I_o/I_m$ 곡선. 두께가 커질수록 곡선이 포화되어 감도를 잃는 구간이 보이도록.
- **fig2-thickness-sensitivity.png** — 두께 역산의 감도. 강도비에 ±5% 오차가 있을 때 역산 두께의 불확실성이 두께에 따라 어떻게 커지는지, 그리고 $\lambda$를 10% 잘못 알았을 때 두께가 얼마나 밀리는지를 함께 보일 것. **이 편의 핵심 그림이다.**
- **fig3-arxps.png** — ARXPS. (a) 두 개의 서로 다른 깊이 분포 $c(z)$(예: 급격한 계면 vs 완만한 확산 계면)를 놓고, (b) 각도별 예측 강도 곡선이 노이즈 수준 안에서 거의 구분되지 않는 것을 보여라. 역문제의 불량조건성을 시각적으로 드러내는 그림이다.
- **fig4-sputter-artifact.png** — 스퍼터 깊이 프로파일 개념도. 참 조성 프로파일(급격한 계면)과, 우선 스퍼터링 + 원자 혼합 + 정보깊이 컨볼루션을 거쳐 관측되는 프로파일을 겹쳐 그려 계면이 어떻게 뭉개지는지 보일 것.

## 코드 자산

```
_code/xps-thickness-arxps-pitfalls/
  generate_figures.py   # 그림 4장
  overlayer.py          # 오버레이어 두께 정·역계산 + ARXPS 순방향 시뮬 + 감도 분석
```

- `overlayer.py`는 단독 실행 시 2절 두께 역산 예제와 fig2의 감도 수치를 출력해야 한다.
- `_code/requirements.txt`에 이미 `scipy`가 추가되어 있어야 한다(4편에서 추가). 없으면 추가하라.

## 내부 링크

- XPS 2편 `/posts/xps-surface-sensitivity-imfp/` — 지수 감쇠(1절), 각도 정의(2절), ARXPS 원리(4절)
- XPS 3편 `/posts/xps-spectrum-structure-satellites/` — Tougaard 배경(필요 시)
- XPS 4편 `/posts/xps-quantification-peak-fitting/` — 균일 가정(도입), 금속/산화물 성분 분해(2절)
- 최적화 3편 `/posts/optimization-levenberg-marquardt/`, 4편 `/posts/optimization-global-heuristics/` — 역문제와 정규화(4절)
- 엘립소메트리 1편 `/posts/ellipsometry-electromagnetic-fresnel/` — 7절 비교

## 참고자료 (말미 섹션)

- D. J. Morgan, "X-Ray Photoelectron Spectroscopy (XPS): An Introduction," Cardiff Catalysis Institute. <https://sites.cardiff.ac.uk/xpsaccess/files/2014/07/AccessXPS_Primer_Paper.pdf>
- B. R. Strohmeier, "An ESCA method for determining the oxide thickness on aluminum alloys," *Surf. Interface Anal.* 15, 51 (1990).
- P. J. Cumpson, "The Thickogram: a method for easy film thickness measurement in XPS," *Surf. Interface Anal.* 29, 403 (2000). <https://www.researchgate.net/publication/230247607_The_Thickogram_a_method_for_easy_film_thickness_measurement_in_XPS>
- "XPS and angle resolved XPS, in the semiconductor industry: Characterization and metrology control of ultra-thin films," *Vacuum* (2010). <https://www.sciencedirect.com/science/article/abs/pii/S0368204810000538>
- "Assessing advanced methods in XPS and HAXPES for determining the thicknesses of high-k oxide materials," *Appl. Surf. Sci.* (2023). <https://www.sciencedirect.com/science/article/abs/pii/S0169433222028458>
- "XPS on corrosion products of ZnCr coated steel: on the reliability of Ar+ ion depth profiling for multi component material analysis." <https://arxiv.org/pdf/1310.7725>
- S. Tougaard, *Surf. Interface Anal.* (비탄성 배경에서 깊이 분포 추출).

## 사실 확인 주의사항

- **각도 규약을 본문 안에서 일관되게 유지하라.** 1절에서 법선 기준 $\cos\theta$를 쓰기로 했으면, 2절에서 Strohmeier 식의 $\sin\theta$를 인용할 때 각도 기준이 다르다는 것을 반드시 명시하라. 이건 이 편에서 가장 실수하기 쉬운 지점이고, 동시에 좋은 교육 포인트다.
- Strohmeier 식의 유효 두께 범위(대략 0~9 nm)와 Thickogram의 정확도(대략 ±10%)는 문헌에 보고된 값이다. 더 구체적인 수치를 쓰려면 웹 검색으로 확인하고 출처를 남겨라.
- "X선은 비파괴"라는 통념을 그대로 쓰지 마라. 6절이 그 통념을 교정하는 절이다.
- 코드는 실제로 실행해 그림과 수치를 확인하고, 본문 수치를 실행 결과와 일치시켜라.
- 시리즈 마지막 편이므로 1~4편 링크가 전부 올바른 slug를 가리키는지 확인하라.
- 작성 후 `grep -n '\$[^$]*|' _posts/2026-09-09-*.md`로 수식 내 리터럴 `|` 확인.

## 커밋

사용자의 명시적 지시 전에는 `git commit` / `git push` 금지.
시리즈 전체를 마친 뒤 커밋을 원하면 편별로 나눠 커밋하는 것을 권한다 (`post: XPS 기초 N편 — ...`).

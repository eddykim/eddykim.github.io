# XPS 시리즈 4편 작성 프롬프트

> 사용법: 이 레포(`eddykim.github.io`) 루트에서 Claude Code를 열고, 아래 `---` 아래 전체를 붙여넣으세요.
> **선행 조건: 1·2·3편이 이미 작성되어 있어야 합니다.**

---

`write-post` 스킬 워크플로우로 아래 포스트를 작성해줘.

**단, 0~4단계는 이미 끝났다.** 아래 "확정 사항"이 그 결과다. 5단계(계획 제시)는 아래 목차를 그대로 확인만 받고 바로 6단계 초안으로 들어가라. 질의응답을 다시 하지 마라.

**작성 시작 전에 반드시 1~3편을 읽어라.** 또한 이 편은 기존 최적화 시리즈와 직접 연결되므로 `_posts/2026-09-02-optimization-levenberg-marquardt.md`도 읽어라 — LM법 설명을 반복하지 말고 링크로 넘겨야 한다.

## 확정 사항

- **시리즈**: "XPS 기초" 전 5편 중 **4편**
- **slug**: `xps-quantification-peak-fitting`
- **포스트 경로**: `_posts/2026-09-08-xps-quantification-peak-fitting.md`
- **제목**: `XPS 기초 4편 — 정량과 피크 피팅`
- **카테고리**: `[표면분석, 스펙트럼피팅]`
- **대상 독자**: 대학원생 + 계측/공정 실무 엔지니어
- **분량**: 9,000자 이상 (시리즈에서 가장 긴 편)
- **수식**: 있음 → `math: true`
- **그림**: 4장
- **코드**: `_code/xps-quantification-peak-fitting/`에 `generate_figures.py` + `xps_fit.py` (실행 가능한 피팅 구현)

### front matter

```yaml
---
layout: post
title: "XPS 기초 4편 — 정량과 피크 피팅"
date: 2026-09-08 20:00:00 +0900
categories: [표면분석, 스펙트럼피팅]
tags: [xps, quantification, peak-fitting, shirley-background, voigt, charge-referencing]
description: "피크 면적에서 조성으로 가는 길에는 감도인자, 배경 모델, 선형, 대전 보정이 줄줄이 끼어 있다. 각각이 어디서 깨지는지 정리하고 직접 피팅해본다."
math: true
---
```

## 시리즈를 관통하는 메시지

**XPS가 실제로 측정하는 것은 광전자의 운동에너지 분포 하나뿐이다. 원소·화학상태·조성비·두께는 전부 모델을 거쳐 역산한 값이며, 그 모델의 가정이 어디서 깨지는지가 이 시리즈의 축이다.**

4편은 이 축이 가장 노골적으로 드러나는 편이다. **"조성 30%"라는 숫자 하나가 나오기까지 최소 네 개의 모델(배경, 선형, 감도인자, 대전 보정)을 통과한다.**

## 스타일 규칙 (반드시 준수)

1. **어조**: 평서문 "~이다".
2. **각 절은 답을 주기 전에 질문을 먼저 던진다.**
3. **수식 안에서 리터럴 `|` 금지.** `\lvert ... \rvert`를 써라. 특히 잔차 $\lvert r_i \rvert$, 적분 상한 표기에서 주의.
4. 소제목 `## 1.` … 번호 매기기. 하위 절이 필요하면 `### 7.1` 형식(최적화 4편 패턴). 말미 `## 정리 및 다음 편 예고` + `## 참고자료`.
5. 기술 용어 영문 병기(첫 등장 시).
6. **30줄 넘는 코드는 본문에 통째로 넣지 마라.** 핵심 함수만 발췌하고 전체는 `_code/xps-quantification-peak-fitting/`에 저장 후 링크.
7. 이미지 경로: `/assets/img/posts/xps-quantification-peak-fitting/...`
8. 도입부는 3편 말미를 이어받아 시작.

## 목차와 각 절에서 다룰 내용

### 도입 (제목 없음)
3편에서 스펙트럼의 각 구조가 무엇인지는 알았다. 이제 "O가 30%, Si가 25%" 같은 숫자를 뽑아야 한다. 그런데 그 숫자는 피크 면적을 재는 순간 정해지는 것이 아니라, 배경을 어디에 긋고, 어떤 함수로 나누고, 어떤 감도인자를 쓰고, 어디를 에너지 기준으로 삼았는지에 전부 의존한다. 이 편은 그 의존성을 하나씩 뜯는다.

### ## 1. 피크 면적에서 조성으로
- 질문: 피크가 크면 그 원소가 많다 — 이 직관은 어디까지 맞는가?
- 첫원리 신호식:
  $$ I_x = A \, n_x \, f \, \sigma_x \, \theta \, \lambda(E_K) \, T(E_K) $$
  ($A$: 분석 면적, $n_x$: 원자 밀도, $f$: X선 플럭스, $\sigma_x$: 광이온화 단면적(Scofield), $\theta$: 각도 인자, $\lambda$: IMFP, $T$: 분광기 투과함수)
- 같은 측정 조건에서 $A, f$는 상수 → 조성비:
  $$ \frac{n_1}{n_2} = \frac{I_1 / (\sigma_1 \lambda_1 T_1)}{I_2 / (\sigma_2 \lambda_2 T_2)} $$
- **핵심 짚기**: $\lambda$가 들어 있다는 것은 2편에서 본 대로 **피크마다 보는 깊이가 다르다**는 뜻이다. 시료가 깊이 방향으로 균일하다는 가정이 이미 여기 들어간다. 균일하지 않으면(오염층, 산화막) 이 조성비는 물리적 조성이 아니라 "지수 가중 평균 조성"이다.

### ## 2. 감도인자(RSF)와 그 불일치
- 질문: 실무에서는 왜 저 식을 직접 쓰지 않는가?
- RSF(relative sensitivity factor)는 $\sigma \lambda T$와 각도 인자를 하나로 묶은 실험적/반경험적 계수. 소프트웨어가 원소별 표를 내장한다.
- **문제**: 출처(장비 제조사, Wagner 계열, Scofield 기반)마다 RSF 값이 다르다. 같은 스펙트럼을 다른 소프트웨어로 정량하면 조성이 달라진다.
- 투과함수 $T(E_K)$는 장비마다 다르고 pass energy 설정에 따라서도 달라진다 → RSF 표를 장비 간에 그대로 옮기면 안 된다.
- 실무 결론: XPS 정량은 **상대 비교에는 강하고 절대값에는 약하다**. 같은 장비·같은 조건으로 잰 시료들 사이의 추세는 믿을 만하지만, 문헌의 절대 조성과 소수점 이하를 비교하는 것은 의미가 없다.

### ## 3. 대전 보정 — C 1s 284.8 eV라는 관행과 그 논쟁
- 질문: 절연체 시료의 에너지 축은 무엇을 기준으로 잡는가?
- 1편 1절로 돌아가라: 도체에서는 시료-분광기 페르미 준위 정렬 덕분에 축이 정해졌다. 절연체에서는 광방출로 표면이 양전위가 되어 피크가 통째로 고결합에너지 쪽으로 밀린다. 전하 중화기(charge neutraliser)를 쓰면 이번엔 반대로 밀릴 수 있다.
- 관행: 모든 대기 노출 표면에 있는 흡착 탄화수소의 C 1s를 284.8 eV로 놓고 축 전체를 평행이동한다(고분자 등에서는 285.0 eV를 쓰기도 한다).
- **논쟁을 정직하게 서술하라**:
  - 이 값의 불확실성은 대개 ±0.2~0.3 eV 수준이다. 화학이동이 0.5~4 eV 자릿수임을 생각하면 무시할 수 없다.
  - 흡착 탄소의 실제 C 1s 위치는 기판에 따라 달라진다. 예를 들어 알루미늄 자연산화막 위에서는 286 eV 부근으로 나타난다는 보고가 있다.
  - 반대로, 다수 시료를 통계적으로 본 최근 연구들은 전기적으로 고립된 시료에서 284.8 eV가 대체로 적절하다고 보고한다.
  - 대안: 도체는 페르미 에지, 내부 기준(예: SiO2의 Si 2p), 그리고 3편에서 본 **Auger 파라미터**(대전에 무관).
- 시리즈 축과 연결: 결합에너지라는 "측정값"조차 모델 선택의 산물이다.

### ## 4. 배경 모델 — 어디에 선을 긋는가
- 질문: 피크 면적을 재려면 먼저 배경을 정해야 하는데, 배경은 측정되지 않는다.
- **linear**: 구간 양 끝을 직선으로 잇는다. 간단하지만 계단 구조를 무시한다.
- **Shirley**: 3편 7절의 계단 배경을 물리적으로 반영한다. 어떤 에너지에서의 배경이 그보다 높은 운동에너지(= 낮은 결합에너지) 쪽 피크 면적에 비례한다고 놓고 반복 계산한다:
  $$ B(E) = k \int_{E}^{E_{\max}} \left[ I(E') - B(E') \right] dE' $$
  자기참조 형태이므로 수렴할 때까지 반복한다. 구현이 간단해 가장 널리 쓰인다.
- **Tougaard**: 비탄성산란 단면적을 직접 모델링한다. 물리적으로 가장 엄밀하고 깊이 분포 정보까지 담지만, 넓은 에너지 구간이 필요하고 사용이 까다롭다. 5편에서 다시 다룸.
- **실무 포인트**: 배경 구간의 양 끝점을 어디에 두느냐만으로 피크 면적이 수 % 움직인다. 시료 간 비교를 할 때는 반드시 같은 구간·같은 배경 모델을 써야 한다.

### ## 5. 선형(lineshape) — Voigt는 왜 Voigt인가
- 질문: 피크는 왜 하필 그 모양인가?
- **Lorentzian**: 코어홀의 유한한 수명에서 오는 자연선폭. 불확정성 원리로 $\Gamma \sim \hbar/\tau$. 원소·준위 고유값.
- **Gaussian**: 장비 기여(X선 선폭, 분석기 분해능)와 시료 기여(포논 broadening, 불균일 대전, 화학적 불균일).
- 두 과정은 독립이므로 관측 선형은 둘의 합성곱 = **Voigt**. 소프트웨어에서는 계산이 싼 pseudo-Voigt(GL 혼합)를 쓰는 경우가 많다.
- **금속은 비대칭**(3편 4절): LA 계열이나 Doniach–Šunjić 형태를 써야 한다. 금속에 대칭 함수를 쓰면 꼬리를 메우려 가짜 산화물 성분이 생긴다.
- 어느 폭이 어디서 오는지 알면, FWHM을 제약할 때 무엇을 묶고 무엇을 풀어야 하는지가 정해진다 → 다음 절.

### ## 6. 제약조건이 필수인 이유 — 피팅 해는 유일하지 않다
- 질문: 잔차가 작으면 옳은 분해인가?
- 성분 수를 늘리면 잔차는 **항상** 줄어든다. 자유 파라미터가 늘어나기 때문이다. 따라서 잔차나 $\chi^2$만으로는 성분 수를 정할 수 없다.
- 물리에서 오는 제약을 반드시 걸어야 한다:
  - 스핀-궤도 이중선의 **면적비 고정**(p 1:2, d 2:3, f 3:4) 및 **분리 폭 고정**
  - 같은 화학종의 이중선 두 성분은 **동일 선형·동일 FWHM**
  - 성분 위치는 문헌값 부근으로 범위 제한
  - 다중항 분리를 갖는 원소는 표준 다중항 모델 사용
- 자유 파라미터가 줄면 잔차는 커지지만 해는 안정된다. 이 맞교환을 명시적으로 설명하라.
- 최적화 시리즈와 연결: 이건 [3편 LM법](/posts/optimization-levenberg-marquardt/)에서 본 비선형 최소자승 문제이고, 성분 위치 초기값에 따라 다른 해로 수렴하는 것은 [4편 전역 최적화](/posts/optimization-global-heuristics/)에서 본 local minimum 문제와 같은 구조다. LM 자체를 다시 설명하지 말고 링크로 넘겨라.

### ## 7. 직접 해보기 — 합성 스펙트럼을 되찾을 수 있는가
- 질문: 정답을 아는 스펙트럼을 만들어 피팅하면 파라미터가 제대로 복원되는가?
- `### 7.1 합성 스펙트럼 만들기`: 알려진 성분(예: 금속 이중선 + 산화물 이중선)과 계단형 배경, 그리고 포아송 노이즈로 스펙트럼을 생성. 참값을 명시.
- `### 7.2 Shirley 배경 구현`: 위 반복식을 numpy로 구현. 수렴 판정 포함.
- `### 7.3 제약조건을 걸고 피팅`: `scipy.optimize.least_squares`(method='lm' 또는 'trf')로 피팅. 제약조건을 파라미터화로 구현(이중선의 두 번째 성분을 독립 파라미터로 두지 않고 첫 성분에서 유도).
- `### 7.4 제약을 풀면 어떻게 되는가`: 같은 데이터에 제약 없이 성분을 하나 더 얹으면 잔차는 줄지만 참값과 멀어진다는 것을 수치로 보여라. 이 절이 이 편의 핵심 실험이다.
- 본문에는 핵심 함수(Shirley 반복, 제약 파라미터화)만 발췌하고 전체 코드는 링크.

### ## 8. 정량값을 보고할 때 지켜야 할 것
- 짧은 실무 체크리스트 절: 배경 모델과 구간, 선형과 제약조건, 감도인자 출처, 대전 보정 기준을 함께 보고해야 결과가 재현 가능하다. 이것이 최근 XPS 커뮤니티가 피팅 관행 문제를 제기하는 이유다.

### ## 정리 및 다음 편 예고
4편 요약 + 5편 예고: 지금까지는 시료가 깊이 방향으로 균일하다고 가정했다. 5편에서는 그 가정을 놓고, 감쇠 자체를 이용해 박막 두께를 재는 방법과 그 함정(ARXPS의 역문제, 우선 스퍼터링, 빔 손상)을 다룬다.

### ## 참고자료

## 그림 명세 (4장)

`_code/xps-quantification-peak-fitting/generate_figures.py`로 생성. 모든 스펙트럼은 **합성 데이터**임을 캡션에 명시.

- **fig1-backgrounds.png** — 같은 합성 스펙트럼에 linear / Shirley / (단순화한) Tougaard형 배경을 각각 적용하고, 그 결과 피크 면적이 얼마나 달라지는지 수치로 함께 표시.
- **fig2-voigt-composition.png** — Gaussian, Lorentzian, 그리고 둘의 합성곱인 Voigt를 같은 축에 겹쳐 그려 꼬리 차이를 보일 것. FWHM 표시.
- **fig3-fit-result.png** — 7절 피팅 결과. 상단: 데이터 점 + 전체 피팅 곡선 + 각 성분 + 배경. 하단: 잔차. 참값과 추정값을 표로 또는 주석으로 병기.
- **fig4-overfitting.png** — 성분 수를 1개씩 늘려가며 잔차 제곱합이 단조 감소하는 것과, 동시에 추정 조성비가 참값에서 멀어지는 것을 두 축으로 보일 것. 이 편의 메시지를 요약하는 그림이다.

## 코드 자산

```
_code/xps-quantification-peak-fitting/
  generate_figures.py     # 그림 4장 생성
  xps_fit.py              # Shirley 배경 + Voigt/비대칭 선형 + 제약조건 피팅
```

- `xps_fit.py`는 단독 실행 가능해야 하고, 실행하면 7절의 수치 결과(참값 vs 추정값, 제약 유무 비교)를 출력해야 한다.
- **`_code/requirements.txt`에 `scipy>=1.11`을 추가하라.** 현재 numpy, matplotlib만 등록되어 있다.

## 내부 링크

- XPS 1편 `/posts/xps-photoemission-binding-energy/` — 대전 보정(3절)
- XPS 2편 `/posts/xps-surface-sensitivity-imfp/` — IMFP가 정량식에 들어오는 지점(1절)
- XPS 3편 `/posts/xps-spectrum-structure-satellites/` — 계단 배경(4절), 비대칭 선형(5절), 스핀-궤도 제약(6절)
- 최적화 3편 `/posts/optimization-levenberg-marquardt/` — 비선형 최소자승(6·7절)
- 최적화 4편 `/posts/optimization-global-heuristics/` — 초기값 의존성(6절)

## 참고자료 (말미 섹션)

- G. H. Major et al., "Practical guide for curve fitting in x-ray photoelectron spectroscopy," *J. Vac. Sci. Technol. A* 38, 061203 (2020). <https://pubs.aip.org/avs/jva/article/38/6/061203/1023652/Practical-guide-for-curve-fitting-in-x-ray>
- "Defining the nature of adventitious carbon and improving its merit as a charge correction reference for XPS," *Appl. Surf. Sci.* (2024). <https://www.sciencedirect.com/science/article/pii/S0169433224000333>
- D. J. Morgan, "The Utility of Adventitious Carbon for Charge Correction: A Perspective From a Second Multiuser Facility," *Surf. Interface Anal.* (2025). <https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/sia.7360>
- "Referencing to adventitious carbon in XPS: Can differential charging explain C 1s peak shifts?," *Appl. Surf. Sci.* (2023). <https://www.sciencedirect.com/science/article/pii/S0169433222023832>
- D. A. Shirley, *Phys. Rev. B* 5, 4709 (1972).
- S. Tougaard, *Surf. Interface Anal.* (비탄성 배경 해석).
- J. H. Scofield, *J. Electron Spectrosc. Relat. Phenom.* 8, 129 (1976) (광이온화 단면적).

## 사실 확인 주의사항

- Shirley 배경 반복식의 적분 방향(운동에너지가 높은 쪽 = 결합에너지가 낮은 쪽)을 헷갈리지 마라. 구현과 본문 서술이 일치해야 한다.
- 284.8 eV 논쟁은 **양쪽 입장을 다 적어라.** 한쪽만 쓰면 사실을 왜곡하게 된다.
- Voigt는 Gaussian과 Lorentzian의 **합성곱**이다. 곱이나 선형결합이 아니다(pseudo-Voigt는 선형결합 근사라는 점을 구분해 서술).
- 코드는 반드시 실제로 실행해서 그림과 수치가 나오는지 확인하고, 본문의 수치는 실행 결과와 일치시켜라. 결과를 지어내지 마라.
- 작성 후 `grep -n '\$[^$]*|' _posts/2026-09-08-*.md`로 수식 내 리터럴 `|` 확인.

## 커밋

사용자의 명시적 지시 전에는 `git commit` / `git push` 금지.

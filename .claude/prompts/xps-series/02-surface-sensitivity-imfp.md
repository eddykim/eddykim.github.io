# XPS 시리즈 2편 작성 프롬프트

> 사용법: 이 레포(`eddykim.github.io`) 루트에서 Claude Code를 열고, 아래 `---` 아래 전체를 붙여넣으세요.
> **선행 조건: 1편(`_posts/2026-09-05-xps-photoemission-binding-energy.md`)이 이미 작성되어 있어야 합니다.**

---

`write-post` 스킬 워크플로우로 아래 포스트를 작성해줘.

**단, 0~4단계(입력 수집 / 분류 판단 / 질의응답 / 이전 포스트 참고 확인)는 이미 끝났다.** 아래 "확정 사항"이 그 결과다. 5단계(계획 제시)는 아래 목차를 그대로 확인만 받고 바로 6단계 초안으로 들어가라. 질의응답을 다시 하지 마라.

**작성 시작 전에 반드시 1편(`_posts/2026-09-05-xps-photoemission-binding-energy.md`)을 읽어라.** 용어·기호·문체를 그대로 이어야 하고, 1편에서 예고한 내용을 받아서 시작해야 한다.

## 확정 사항

- **시리즈**: "XPS 기초" 전 5편 중 **2편**
- **slug**: `xps-surface-sensitivity-imfp`
- **포스트 경로**: `_posts/2026-09-06-xps-surface-sensitivity-imfp.md`
- **제목**: `XPS 기초 2편 — 왜 표면 10 nm만 보이는가`
- **카테고리**: `[표면분석, 깊이분석]`
- **대상 독자**: 대학원생 + 계측/공정 실무 엔지니어
- **분량**: 8,000자 이상
- **수식**: 있음 → `math: true`
- **그림**: 개념도 + 시뮬레이션, 3장
- **코드**: `_code/xps-surface-sensitivity-imfp/generate_figures.py`

### front matter

```yaml
---
layout: post
title: "XPS 기초 2편 — 왜 표면 10 nm만 보이는가"
date: 2026-09-06 20:00:00 +0900
categories: [표면분석, 깊이분석]
tags: [xps, imfp, surface-analysis, information-depth, arxps]
description: "X선은 시료 속으로 마이크로미터를 파고드는데 XPS는 왜 표면 10 nm 기법인가. 비탄성 평균자유행로와 정보깊이를 정리한다."
math: true
---
```

## 시리즈를 관통하는 메시지

**XPS가 실제로 측정하는 것은 광전자의 운동에너지 분포 하나뿐이다. 원소·화학상태·조성비·두께는 전부 모델을 거쳐 역산한 값이며, 그 모델의 가정이 어디서 깨지는지가 이 시리즈의 축이다.**

2편에서 이 축이 걸리는 지점: **"표면 10 nm"라는 숫자 자체가 감쇠 모델에서 나온 값이며, IMFP 값은 대개 직접 측정한 것이 아니라 예측식(TPP-2M)으로 계산한 값이다.**

## 스타일 규칙 (반드시 준수)

1. **어조**: 평서문 "~이다". 1편 및 최근 포스트와 동일하게.
2. **각 절은 답을 주기 전에 질문을 먼저 던진다.** 모든 `##` 절이 "무엇이 문제인가"로 시작한 뒤 식이나 정의를 내놓는다.
3. **수식 안에서 리터럴 `|` 금지.** 절댓값은 `\lvert ... \rvert`. kramdown 렌더링이 깨진다.
4. 소제목 `## 1.` `## 2.` … 번호 매기기. 말미는 `## 정리 및 다음 편 예고` + `## 참고자료`.
5. 기술 용어 영문 병기(첫 등장 시).
6. 30줄 넘는 코드는 발췌만, 전체는 `_code/<slug>/`.
7. 이미지 경로: `/assets/img/posts/xps-surface-sensitivity-imfp/...`
8. 도입부는 1편 말미를 이어받는 방식으로 시작하라 — 기존 시리즈 포스트들의 확립된 패턴이다. 예: "[1편](/posts/xps-photoemission-binding-energy/)에서는 …로 마무리했다. 그런데 …"

## 목차와 각 절에서 다룰 내용

### 도입 (제목 없음)
1편에서 코어 준위 결합에너지가 원소 지문이 된다는 데까지 왔다. 그런데 XPS를 "표면 분석 기법"이라고 부르는 이유는 아직 설명되지 않았다. Al Kα X선은 고체 속으로 마이크로미터 단위까지 침투한다 — 그러니 표면민감성은 X선이 만드는 것이 아니다. 그럼 무엇이 만드는가.

### ## 1. 표면민감성을 만드는 것은 X선이 아니라 전자다
- 질문: X선은 깊이 들어가는데 왜 신호는 표면에서만 오는가?
- 광자의 감쇠길이와 전자의 감쇠길이가 서로 세 자릿수 다르다.
- 깊은 곳에서 생성된 광전자도 존재하지만, 밖으로 나오는 길에 비탄성산란(inelastic scattering)을 겪으면 운동에너지를 잃는다. 잃은 전자는 원래 피크 위치에 나타나지 못하고 배경(계단형 배경)에 섞인다.
- 즉 **"피크에 기여하는 전자"는 무손실로 탈출한 전자뿐**이라는 것이 표면민감성의 정체다. 이 관점은 4편 배경 모델에서 다시 쓰인다.

### ## 2. IMFP의 정의와 감쇠식
- 질문: "무손실로 탈출할 확률"을 어떻게 정량화하는가?
- IMFP(inelastic mean free path) $\lambda$의 정의: 전자가 비탄성산란을 겪지 않을 확률이 $1/e$가 되는 이동 거리.
- 깊이 $z$에서 생성된 전자가 무손실로 탈출할 확률: $\exp\!\left(-\dfrac{z}{\lambda\cos\theta}\right)$ — 여기서 $\theta$는 표면 법선에서 잰 검출 각도.
- Beer–Lambert와 형태는 같지만 물리적 주체가 다르다(광자의 흡수 vs 전자의 비탄성산란)는 점을 짚어라.
- IMFP와 EAL(effective attenuation length)의 구분: 탄성산란까지 포함하면 실효 감쇠길이는 IMFP보다 짧다. 실무에서 두 값을 섞어 쓰다 오차가 생긴다는 점을 언급.

### ## 3. 유니버설 커브 — $\lambda$는 운동에너지에 어떻게 의존하는가
- 질문: 어느 피크가 더 표면을 보는가?
- 유니버설 커브(universal curve): 대부분의 물질에서 $\lambda$가 운동에너지에 대해 비슷한 모양을 그린다. 최소값이 대략 50~100 eV 부근.
- 수치 감각: $E_K \approx 50$ eV에서 약 0.5 nm, 100~1400 eV 구간에서 대략 1~3 nm, 1 keV 이상에서 수 nm~수십 nm.
- 고에너지 쪽 근사: $\lambda \propto E_K^{\,n}$, $n \approx 0.5\text{–}0.75$.
- 실무: IMFP는 보통 직접 측정하지 않고 **TPP-2M 예측식**(Tanuma–Powell–Penn)이나 NIST SRD 71 데이터베이스에서 가져온다. 즉 두께·조성 계산에 들어가는 $\lambda$는 이미 모델값이다 — 시리즈 축을 여기서 한 번 짚어라.
- 같은 시료에서도 피크마다 $E_K$가 다르므로 **피크마다 보는 깊이가 다르다**. 이 사실이 5편의 두께 역산과 3편의 강도비 해석에 그대로 걸린다.

### ## 4. 정보깊이 $3\lambda$가 95%인 이유
- 질문: "정보깊이 10 nm"라는 숫자는 어디서 나왔는가?
- 깊이별 기여도를 적분해서 유도하라. 표면부터 깊이 $d$까지의 누적 기여:
  $$ \frac{\int_0^{d} e^{-z/\lambda}\,dz}{\int_0^{\infty} e^{-z/\lambda}\,dz} = 1 - e^{-d/\lambda} $$
- $d = \lambda$ → 63.2%, $d = 2\lambda$ → 86.5%, $d = 3\lambda$ → 95.0%.
- 따라서 $\lambda \approx 3$ nm이면 정보깊이 $3\lambda \approx 10$ nm. "10 nm"는 물리 상수가 아니라 95%라는 임의 기준과 $\lambda$ 추정값의 곱이다.
- 여기서 중요한 따름결과: 신호는 깊이에 대해 **지수 가중 평균**이지 균일 평균이 아니다. 최표면 1 nm가 전체 신호의 30% 가까이를 차지한다.

### ## 5. 각도를 바꾸면 보는 깊이가 바뀐다 — 그리고 문헌상의 각도 혼동
- 질문: 더 얕게 보고 싶으면 무엇을 바꿔야 하는가?
- 유효 깊이 $d = \lambda\cos\theta$ ($\theta$ = 표면 법선 기준 검출 각도). $\theta$를 키우면(시료를 기울이면) 탈출 경로가 길어져 표면민감성이 올라간다.
- **문헌 혼동을 명시적으로 정리하라**: "take-off angle"이 표면 평면에서 잰 각으로 정의된 문헌에서는 같은 물리가 $\lambda\sin\theta$로 적힌다(예: Strohmeier 식). $\cos$이냐 $\sin$이냐는 각도 기준의 차이일 뿐 물리는 같다. 남의 식을 가져다 쓸 때 반드시 각도 정의를 확인해야 하며, 이건 실제로 두께 계산에서 자주 나는 실수다.
- 이 원리를 각도 시리즈로 확장한 것이 ARXPS이며, 5편에서 본격적으로 다룬다고 예고.

### ## 6. 표면민감성의 대가
- 질문: 표면만 보는 것은 언제 손해인가?
- 대기 노출된 모든 표면에는 흡착 탄화수소층(adventitious carbon)이 있다. 두께가 IMFP와 같은 자릿수(대략 1~2 nm)이므로 하부 신호를 유의하게 감쇠시킨다. 정보깊이의 상당 부분이 "알고 싶지 않은 것"에 쓰인다.
- 벌크 정보는 원리적으로 얻을 수 없다. 표면 편석(surface segregation)이 있으면 XPS 조성은 벌크 조성과 다르다 — 이건 오차가 아니라 기법의 정의다.
- 시료 취급·이송·UHV(약 $10^{-10}$ mbar) 요구가 여기서 나온다. 왜 초고진공이 필요한가: 잔류 가스에 의한 재오염 시간 + 광전자의 기체 중 산란.

### ## 7. 광학 계측과 비교하면 무엇이 다른가
- 엘립소메트리는 흡수가 없는 파장대에서 수백 nm~수 µm를 투과하며 다층 구조 전체를 본다. XPS는 최상단 10 nm만 본다.
- 그래서 둘은 경쟁 기법이 아니라 상보적이다: XPS는 화학상태를, 엘립소메트리는 두께·광학상수를 잘 준다. 5편에서 두께 계측을 다룰 때 이 비교로 돌아온다고 예고.
- 엘립소메트리 1편(`/posts/ellipsometry-electromagnetic-fresnel/`)의 침투 깊이 $\delta$ 논의로 내부 링크.

### ## 정리 및 다음 편 예고
2편 요약 + 3편 예고: 스펙트럼에는 원소당 피크 하나만 있는 것이 아니다. 스핀-궤도 분리, shake-up 위성, 다중항 분리, 플라즈몬 손실, Auger 계열이 함께 있고, 각각이 초기상태 효과인지 최종상태 효과인지 구분해야 화학상태를 읽을 수 있다.

### ## 참고자료

## 그림 명세 (3장)

`_code/xps-surface-sensitivity-imfp/generate_figures.py` 하나로 생성. 1편 스크립트의 헤더·폰트(`AppleGothic`)·출력경로 패턴을 따를 것.

- **fig1-imfp-universal-curve.png** — IMFP vs 운동에너지 로그-로그 곡선. TPP-2M 형태의 근사식이나 유니버설 커브 경험식으로 곡선을 그리고, Al Kα로 여기했을 때 대표 피크들의 $E_K$ 위치를 점으로 표시. 최소값 영역과 $\lambda\propto E^{0.5\text{–}0.75}$ 영역을 구분해 표시.
- **fig2-depth-contribution.png** — 깊이별 신호 기여. (a) 지수 감쇠 $e^{-z/\lambda}$ 곡선과 깊이 구간별 기여 막대, (b) 누적 기여 $1-e^{-d/\lambda}$ 곡선에 $\lambda, 2\lambda, 3\lambda$ 지점의 63.2/86.5/95.0% 표시.
- **fig3-takeoff-angle.png** — (a) ARXPS 기하 개념도: 시료 표면, 법선, 검출 방향, 탈출 경로 길이 비교. (b) 각도별 유효 깊이 $\lambda\cos\theta$와 각도별 누적 기여 곡선을 겹쳐 그려, 각도를 키우면 표면 가중이 커지는 것을 보일 것.

## 코드 자산

- `_code/xps-surface-sensitivity-imfp/generate_figures.py`
- 추가 의존성 없음(numpy, matplotlib).
- 본문 발췌: 누적 기여 계산과 각도 의존 부분 정도만.

## 내부 링크

- XPS 1편: `/posts/xps-photoemission-binding-energy/` — 도입부
- 엘립소메트리 1편: `/posts/ellipsometry-electromagnetic-fresnel/` — 7절 침투 깊이 비교

## 참고자료 (말미 섹션)

- D. J. Morgan, "X-Ray Photoelectron Spectroscopy (XPS): An Introduction," Cardiff Catalysis Institute. <https://sites.cardiff.ac.uk/xpsaccess/files/2014/07/AccessXPS_Primer_Paper.pdf>
- S. Tanuma, C. J. Powell, D. R. Penn, "Calculations of electron inelastic mean free paths," *Surf. Interface Anal.* (TPP-2M 계열 논문)
- NIST Standard Reference Database 71: Electron Inelastic-Mean-Free-Path Database
- C. J. Powell, A. Jablonski, NIST Electron Effective-Attenuation-Length Database (SRD 82)
- "Review on surface-characterization applications of X-ray photoelectron spectroscopy (XPS)," *Appl. Surf. Sci. Adv.* (2022). <https://www.sciencedirect.com/science/article/pii/S2666523922001222>

## 사실 확인 주의사항

- IMFP와 EAL(effective attenuation length)을 혼용하지 마라. 본문에서 구분해 설명한 뒤 일관되게 써라.
- TPP-2M의 정확한 계수식을 본문에 옮길 필요는 없다. 필요하면 웹 검색으로 확인하고, 확인 못한 계수를 지어내지 마라. 그림은 유니버설 커브의 정성적 형태를 재현하는 수준이면 충분하며, 그 사실을 그림 캡션이나 본문에 명시하라.
- 63.2 / 86.5 / 95.0% 숫자는 $1-e^{-n}$에서 직접 나오는 값이므로 그대로 써도 된다.
- 작성 후 자체 점검: front matter 유효성, 이미지 경로 실재, 수식 내 리터럴 `|` 없음, 30줄 초과 코드 블록 없음.

## 커밋

사용자의 명시적 지시 전에는 `git commit` / `git push` 금지.

# XPS 기초 시리즈 — 채팅 핸드오프 프롬프트

각 편을 **별도 채팅에서 순차적으로** 작성하기 위한 프롬프트 모음이다.
파일 하나가 채팅 하나에 붙여넣을 프롬프트 전문이며, 각각 독립적으로 동작하도록 스타일 규칙·조사 결과·그림 명세·참고자료를 모두 담고 있다.

## 진행 순서

| # | 프롬프트 파일 | 포스트 | 상태 |
|---|---|---|---|
| 1 | [01-photoemission-binding-energy.md](01-photoemission-binding-energy.md) | `2026-09-05-xps-photoemission-binding-energy.md` | ☐ |
| 2 | [02-surface-sensitivity-imfp.md](02-surface-sensitivity-imfp.md) | `2026-09-06-xps-surface-sensitivity-imfp.md` | ☐ |
| 3 | [03-spectrum-structure-satellites.md](03-spectrum-structure-satellites.md) | `2026-09-07-xps-spectrum-structure-satellites.md` | ☐ |
| 4 | [04-quantification-peak-fitting.md](04-quantification-peak-fitting.md) | `2026-09-08-xps-quantification-peak-fitting.md` | ☐ |
| 5 | [05-thickness-arxps-pitfalls.md](05-thickness-arxps-pitfalls.md) | `2026-09-09-xps-thickness-arxps-pitfalls.md` | ☐ |

**앞 편이 커밋(또는 최소한 작성 완료)된 뒤에 다음 편을 시작하라.** 2편부터는 앞 편 본문을 읽어 용어·기호·문체를 잇도록 프롬프트에 지시되어 있다.

## 사전 준비 (완료됨)

- `categories.yaml`에 대분류 `표면분석` 추가 완료.
  소분류: 광전자분광 / 장비구성 / 정량분석 / 스펙트럼피팅 / 깊이분석
- 4편에서 `_code/requirements.txt`에 `scipy>=1.11`을 추가해야 한다(현재 numpy, matplotlib만 등록).

## 시리즈를 관통하는 메시지

> XPS가 실제로 측정하는 것은 광전자의 운동에너지 분포 하나뿐이다.
> 원소·화학상태·조성비·두께는 전부 모델을 거쳐 역산한 값이며,
> 그 모델의 가정이 어디서 깨지는지가 이 시리즈의 축이다.

편별로 이 축이 걸리는 지점:

1. **1편** — 결합에너지 스케일 자체가 에너지 보존식 + 일함수 보정의 산물
2. **2편** — "표면 10 nm"는 감쇠 모델과 예측식으로 계산한 $\lambda$의 곱
3. **3편** — 무엇을 "성분"으로 볼 것인가가 이미 해석
4. **4편** — 조성 숫자 하나가 배경·선형·감도인자·대전 보정 네 모델을 통과
5. **5편** — 두께는 가장 간접적인 양이고, 측정이 시료를 바꾸기까지 한다

## 전 편 공통 스타일 규칙 (각 프롬프트에도 포함되어 있음)

- 평서문 "~이다" (CLAUDE.md에는 존댓말로 적혀 있으나 실제 최근 포스트는 전부 평서문 — 실제 포스트를 따른다)
- 모든 `##` 절은 답을 주기 전에 질문을 먼저 던진다
- 수식 안 리터럴 `|` 금지 → `\lvert` / `\rvert` (kramdown 표 렌더링 깨짐)
- `## 1.` 번호 매기기, 말미 `## 정리 및 다음 편 예고` + `## 참고자료`
- 30줄 넘는 코드는 발췌만, 전체는 `_code/<slug>/`
- 이미지는 루트 기준 경로 `/assets/img/posts/<slug>/...`
- 사용자의 명시적 지시 없이 `git commit` / `git push` 금지

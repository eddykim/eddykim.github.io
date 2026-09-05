# 기존 포스트 재발행 — 채팅 핸드오프 프롬프트

> 사용법: 이 레포(`eddykim.github.io`) 루트에서 Claude Code를 열고, 아래 `---` 아래 전체를 붙여넣으세요.
> 글 하나를 끝낼 때마다 아래 진행표의 체크박스를 갱신하고 커밋하면, 다음 채팅이 상태를 이어받습니다.

---

`revise-post` 스킬 워크플로우로 재발행 작업을 이어서 해줘. 대상은 **진행표에서 아직 체크되지 않은 가장 위의 글 한 편**이다.

## 지금까지의 경위

블로그를 처음부터 다시 시작하기로 하고, 2026-09-05에 **전체 포스트 11편을 `_drafts/`로 내렸다.** 한국어 문장을 논문체로 다듬고 영문판을 함께 만들어 하나씩 재발행하는 중이다.

일부만 내리면 남은 글이 내려간 글을 링크해 htmlproofer가 실패하고 배포 전체가 막히기 때문에 전부 내렸다. 되돌릴 때도 같은 이유로 **순서를 지켜야 한다.**

## 진행표 (이 순서를 지킬 것)

```
[x]  1  start-blog                                    2026-09-05 발행
[x]  2  optimization-gradient-descent                 2026-09-07 예약
[ ]  3  raytracing-spherical-lens-refraction          ← 다음
[ ]  4  raytracing-total-internal-reflection
[ ]  5  raytracing-reflection-family
[ ]  6  optimization-newton-gauss-newton              ※ 링크 복원 필요
[ ]  7  optimization-levenberg-marquardt
[ ]  8  optimization-global-heuristics                ※ 링크 복원 필요
[ ]  9  ellipsometry-electromagnetic-fresnel
[ ] 10  ellipsometry-polarization-mueller-matrix
[ ] 11  ellipsometry-thin-film-multilayer-reflectance
```

**순서를 어기면 배포가 막힌다.** 글 사이의 링크가 전부 과거 글을 향하므로, 이 순서대로 올리면 매 시점에 링크가 유효하다. 예를 들어 6번을 2번보다 먼저 올리면 1편으로 가는 링크가 404가 된다.

## 확정 사항 (다시 묻지 말 것)

- **어투**: 한국어 논문체. 평서문 `~이다`. 구어체 축약(`~는 게`)과 1인칭은 제거한다.
  단, 수사적 질문("~하는 이유는 무엇인가")은 이 블로그 글의 강점이므로 남긴다.
  소개글만 예외적으로 존댓말이며 이미 발행됐다.
- **영문판**: 모든 글에 함께 만든다. 번역이 아니라 영어로 다시 쓴다.
- **발행 주기**: 주 3회 **월·수·금 20:00 KST**. 다음 빈 슬롯에 배정한다.
  이미 9/7(월)에 최적화 1편이 잡혀 있으므로 3번 글은 9/9(수)부터다.
- **카테고리**: 한국어는 `categories.yaml`의 분류를, 영문은 같은 파일의 `en_names`
  대응표를 쓴다. 이미 draft 파일들의 카테고리는 새 체계로 갱신해 두었다.
- **커밋·푸시**: 사용자의 명시적 지시가 있을 때만. CLAUDE.md 절대 규칙이다.

## 글마다 할 일

`revise-post` 스킬의 단계를 그대로 따르되, 아래는 이 작업에 특히 해당한다.

### 1) 사실 오류부터 찾는다

어투보다 먼저다. 최적화 1편에서 `$n = n - ik$`처럼 좌우변에 같은 기호가 들어가
성립하지 않는 식을 찾아 `$N = n - ik$`로 고쳤다. **다른 편에도 있을 수 있으니
수식을 그냥 넘기지 마라.** 시리즈 간 기호 일관성, 중복 서술, `http://` 링크도 함께 본다.

### 2) 그림 영문화

그림 라벨이 한국어이면 영문판에 쓸 수 없다. `_code/<slug>/generate_figures.py`를
언어별 라벨 사전 구조로 고쳐 `assets/img/posts/<slug>/en/`에 다시 생성한다.

**`_code/optimization-gradient-descent/generate_figures.py`가 이미 그 구조로 되어 있으니
그대로 본떠라.** 핵심은 데이터 계산을 한 번만 하고 `render(labels)`로 라벨만 갈아
끼우는 것이다. 실행 후 `git status`로 **한국어 그림이 변하지 않았는지 확인한다.**
같은 seed를 쓰므로 바이트 단위로 동일해야 한다.

### 3) front matter

한국어판에 `page_id: <slug>` 한 줄을 추가하고, 영문판에는 네 줄을 넣는다.

```yaml
lang: en
lang-exclusive: ["en"]
permalink: /posts/<slug>/
page_id: <slug>
```

`page_id`가 양쪽에 있어야 hreflang이 붙는다. 없으면 검색엔진이 두 글을 중복
콘텐츠로 보고 서로 순위를 깎는다.

### 4) 링크 복원 (6번과 8번에만 해당)

전체를 내릴 때 최적화 2편·4편이 1편을 링크하고 있어 htmlproofer가 실패했다.
임시로 일반 텍스트 `1편`로 바꿔둔 상태다. 해당 글을 재발행할 때 되돌려라.

```markdown
1편  →  [1편](/posts/optimization-gradient-descent/)
```

1편은 9/7에 발행되므로 그 뒤라면 링크가 유효하다.

### 5) 로컬 검증 (반드시)

푸시 전에 돌린다. 이 검증이 없어서 배포가 사흘간 멈춘 적이 있다.

```bash
export PATH="$HOME/.rbenv/shims:$PATH"
JEKYLL_ENV=production bundle exec jekyll b -d _site --future
bundle exec htmlproofer _site --disable-external \
  --ignore-urls "/^http:\/\/127.0.0.1/,/^http:\/\/0.0.0.0/,/^http:\/\/localhost/"
rm -rf _site
```

`--future`가 있어야 예약 날짜 글까지 검증된다. 확인할 것:

- `_site/posts/<slug>/`와 `_site/en/posts/<slug>/`가 **둘 다** 생겼는가
- hreflang이 서로를 가리키는가
- 영문 카테고리 경로가 영문인가 (`/en/categories/optimization/`)
- htmlproofer 통과

### 6) 진행표 갱신

작업을 마치면 이 파일의 체크박스와 발행일을 갱신해서 함께 커밋한다.
다음 채팅이 상태를 이어받는다.

## 환경

- 로컬 Ruby는 rbenv에 3.4.10이 있다. 셸이 shim을 안 잡으므로 위처럼 PATH를 직접 넣는다.
  시스템 기본 `ruby`는 2.6이라 Chirpy가 돌지 않는다.
- 의존성은 `vendor/bundle`에 설치돼 있다.
- 사이트는 이중언어다. 한국어는 루트(`/posts/...`), 영문은 `/en/posts/...`.
- 참고자료가 필요하면 `research-sources` 스킬을 쓴다. `_references/`에 이미 받아둔
  자료가 있을 수 있으니 새로 찾기 전에 먼저 본다.

## 하지 말 것

- **원문의 목소리를 지우지 마라.** 소개글을 과하게 다듬었다가 전부 되돌린 적이 있다.
  문법이 깨진 곳과 구어체를 고치되, 멀쩡한 문장을 더 좋게 만들려 하지 마라.
- 굵은 글씨로 항목 라벨을 다는 구성을 쓰지 마라. 전형적인 AI 글 문법이다.
- 한국어판이 확정되기 전에 번역하지 마라. 곧 고칠 문장을 옮기면 작업이 두 번 된다.
- 순서를 건너뛰지 마라.

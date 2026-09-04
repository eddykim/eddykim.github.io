# CLAUDE.md — 블로그 레포 운영 규칙

이 레포는 GitHub Pages(Jekyll) 기반 기술 블로그입니다.
블로그 포스트 작성 요청을 받으면 반드시 `.claude/skills/write-post/SKILL.md`의 워크플로우를 따르세요.

## 절대 규칙 (위반 금지)

1. **사용자의 명시적 확인 없이 `git commit` / `git push`를 실행하지 않는다.**
   "커밋해줘", "푸시해줘", "올려줘" 같은 명시적 지시가 있을 때만 실행한다.
2. **`git push --force`, 히스토리 재작성(rebase, reset --hard)은 절대 하지 않는다.**
3. `_config.yml`, 테마 파일(`_layouts/`, `_includes/`, `_sass/`)은 사용자가 직접 요청하지 않는 한 수정하지 않는다.
4. 기존 포스트를 수정·삭제하지 않는다. 새 포스트 작성이 기본이다.
5. 글 계획(outline) 승인을 받기 전에는 초안 작성을 시작하지 않는다.

## 디렉토리 구조와 저장 경로

```
_posts/                          # 포스트 본문
  YYYY-MM-DD-slug.md
assets/
  img/posts/<slug>/               # 해당 포스트 전용 이미지·미디어
  files/posts/<slug>/             # 다운로드용 첨부 파일 (데이터, PDF 등)
_code/
  <slug>/                         # 포스트에 사용된 전체 코드 (실행 가능한 형태)
```

- `slug`는 영문 소문자와 하이픈만 사용 (예: `snapshot-ellipsometry-intro`)
- 포스트 파일명(`_posts/`)에만 날짜를 붙이고, `assets/`와 `_code/`의 slug 폴더명에는 날짜를 붙이지 않는다.
- 포스트 하나에 딸린 자산은 반드시 같은 slug 폴더에 모아 저장한다.
- 이미지 삽입은 상대 경로가 아닌 루트 기준 경로 사용:
  `![설명](/assets/img/posts/slug/figure1.png)`

## Front matter 템플릿

```yaml
---
layout: post
title: "포스트 제목"
date: YYYY-MM-DD HH:MM:SS +0900
categories: [대분류, 소분류]  # categories.yaml 참고, 소분류 없으면 대분류만
tags: [태그1, 태그2]        # 소문자, 3~6개
description: "검색·미리보기용 한 줄 요약 (80자 이내)"
image: /assets/img/posts/slug/thumbnail.png   # 있을 경우만
---
```

카테고리 체계는 `categories.yaml`이 단일 출처(source of truth)다. 대분류·소분류·분류별 질문 세트가
여기 정의되어 있으며, 새 대분류를 추가할 때는 `categories.yaml`의 `fallback.new_category_checklist`를
거쳐 사용자 승인 후 이 파일에 반영한다. 분류 판단·질의응답 절차는
`.claude/skills/write-post/SKILL.md`를 따른다.

## 글 스타일

- 언어: 한국어 본문 + 기술 용어는 영문 병기 (예: "후초점면(back focal plane, BFP)")
- 어조: 평서문 존댓말("~입니다"), 과장 없이 담백하게
- 구조: 도입(왜 이 주제인가) → 본문(소제목 `##` 단위) → 정리/다음 글 예고
- 분량 기준: 일반 포스트 1,500~3,000자, 튜토리얼은 코드 포함 유동적
- 수식은 MathJax 문법 (`$...$`, `$$...$$`)
- 코드는 언어 지정 코드 블록 사용. **30줄이 넘는 코드는 본문에 전부 넣지 말고**
  핵심 부분만 발췌하고 전체는 `_code/slug/`에 저장 후 링크한다.
- 참고자료는 글 말미에 "참고자료" 섹션으로 출처(제목, 링크) 명시
- 외부 자료의 문장을 그대로 옮기지 않는다. 반드시 자신의 문장으로 재구성한다.

## 영문판 작성 (jekyll-polyglot)

한국어가 기본 언어이며 루트(`/posts/...`)를 쓴다. 영문판은 `/en/posts/...`에 생성된다.
기존 한국어 URL은 바뀌지 않는다.

- **번역이 없는 한국어 글은 아무것도 손대지 않는다.** `_config.yml`의 defaults가
  `lang: ko`와 `lang-exclusive: ["ko"]`를 자동으로 넣어 `/en/`에 노출되지 않게 한다.
- **영문판을 새로 쓸 때**는 파일명을 `YYYY-MM-DD-<slug>-en.md`로 하고 front matter에
  아래 네 줄을 반드시 넣는다.

```yaml
lang: en
lang-exclusive: ["en"]          # 빠뜨리면 기본값 ko가 남아 /en/ 에서 사라진다
permalink: /posts/<slug>/       # 파일명의 -en 접미사가 URL에 새지 않게 고정
page_id: <slug>                 # 한국어판과 짝을 맺는 식별자
```

- **그리고 한국어 원문에도 같은 `page_id: <slug>` 한 줄을 추가한다.** 이게 있어야
  두 글이 서로의 번역본으로 인식되어 hreflang이 붙는다. 없으면 검색엔진이 두 글을
  중복 콘텐츠로 보고 서로 순위를 깎는다.

`_includes/metadata-hook.html`이 `page.available_languages`를 보고 hreflang을 자동 생성하므로
직접 손댈 일은 없다. canonical은 jekyll-seo-tag가 담당하므로 polyglot의 `I18n_Headers`
태그는 쓰지 않는다 (쓰면 canonical이 중복된다).

## 커밋 컨벤션

- 포스트 추가: `post: <포스트 제목>`
- 포스트 수정: `fix(post): <slug> — 수정 내용 한 줄`
- 자산만 추가: `assets: <slug> 이미지/코드 추가`
- 커밋 전 반드시 `git status`와 `git diff --stat`을 사용자에게 보여주고 확인받는다.
- 커밋 대상 파일은 개별 지정한다 (`git add <파일>...`). `git add -A` 금지.

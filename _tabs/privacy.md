---
title: 개인정보·라이선스
# 기본 permalink 는 /:title/ 이라 한국어 제목이 슬러그가 될 수 있어 명시적으로 고정합니다.
permalink: /privacy/
icon: fas fa-shield-halved
order: 5
---

이 문서는 본 블로그(`https://eddykim.github.io`)가 방문자의 정보를 어떻게 다루는지, 그리고 이곳에 실린
글과 코드를 어떤 조건으로 쓸 수 있는지를 밝힙니다.

본 블로그는 개인이 운영하는 기술 기록용 사이트이며, 회원가입·로그인 절차가 없고 이름·연락처 같은
개인정보를 직접 수집하는 양식(form)을 두지 않습니다. 아래에 적은 것은 모두 외부 서비스를 통해
자동으로 처리되는 항목입니다.

## 1. 방문 통계

방문 추이를 파악하기 위해 아래 도구를 사용합니다.

| 도구 | 수집 항목 | 쿠키 |
| --- | --- | --- |
| [GoatCounter](https://www.goatcounter.com/) | 페이지 주소, 유입 경로, 브라우저·화면 크기 등 익명 집계값 | 사용하지 않음 |
| [Google Analytics 4](https://policies.google.com/privacy) | 페이지 주소, 유입 경로, 대략적 위치, 기기 정보, 방문 식별자 | 사용함 |

GoatCounter는 쿠키를 심지 않고 개인을 식별하지 않습니다. 반면 Google Analytics 4는 방문자를
구분하기 위한 쿠키를 사용하며, 수집된 정보는 Google의 개인정보처리방침에 따라 처리됩니다.
브라우저의 쿠키 차단 설정이나 [Google Analytics 차단 부가기능](https://tools.google.com/dlpage/gaoptout)으로
수집을 거부할 수 있으며, 거부해도 글을 읽는 데에는 아무런 제약이 없습니다.

이 통계는 어떤 글이 읽히는지를 보기 위한 용도로만 쓰며, 제3자에게 판매하거나 제공하지 않습니다.

## 2. 댓글

댓글은 [giscus](https://giscus.app)를 통해 GitHub Discussions에 저장됩니다. 본 블로그는 댓글을 담는
별도의 데이터베이스를 두지 않으며, 작성된 내용과 계정 정보는 GitHub에 귀속되어
[GitHub 개인정보처리방침](https://docs.github.com/en/site-policy/privacy-policies/github-privacy-statement)의
적용을 받습니다.

댓글을 쓰려면 GitHub 계정으로 giscus 앱을 인가해야 합니다. 이 인가는 방문자의 계정으로 댓글을
게시하기 위한 것이며, 운영자는 방문자의 GitHub 자격증명에 접근하지 않습니다. 본인이 작성한 댓글은
해당 저장소의 Discussions에서 직접 수정하거나 삭제할 수 있습니다.

## 3. 외부 리소스

페이지를 표시하는 과정에서 글꼴, 수식 렌더링(MathJax), 다이어그램(Mermaid) 등의 자원을
Google Fonts와 jsDelivr CDN에서 내려받습니다. 이 과정에서 해당 서비스에 방문자의 IP 주소가 전달될
수 있습니다.

## 4. 광고

현재 이 블로그에는 광고가 없습니다. 이후 광고를 도입할 경우 이 문서를 먼저 갱신하고, 광고 사업자가
쿠키를 사용하는지 여부와 거부 방법을 함께 안내하겠습니다.

## 5. 콘텐츠 이용 조건

본 블로그의 글·그림·코드는 별도 표시가 없는 한 저작자에게 권리가 있습니다. 기계 판독이 가능한 형태의
이용 조건을 [RSL(Really Simple Licensing)](https://rslstandard.org/) 표준에 따라
[`/license.xml`](/license.xml)에 게시하고 있으며, 요약하면 다음과 같습니다.

- **검색엔진 색인**: 자유롭게 허용합니다. 별도의 대가나 절차가 필요하지 않습니다.
- **사람이 읽고 인용하는 경우**: 출처와 원문 링크를 밝히면 자유롭게 인용할 수 있습니다.
- **AI 학습·추론·검색 색인에 사용하는 경우**: 가시적인 출처 표기와 원문 링크가 필요하며,
  `/license.xml`에 명시한 조건에 따른 대가 지급 대상입니다. 학습 코퍼스 포함, 그리고
  생성 결과에 기여하는 추론·그라운딩(RAG) 사용이 모두 여기에 해당합니다.

조건 협의나 별도 라이선스가 필요하시면 아래 연락처로 문의해 주십시오.

## 6. 문의

- 이메일: [eddyoptics@gmail.com](mailto:eddyoptics@gmail.com)
- GitHub: [github.com/eddykim](https://github.com/eddykim)

## 7. 변경 이력

- 2026-09-04: 최초 게시

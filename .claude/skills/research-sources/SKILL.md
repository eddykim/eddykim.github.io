---
name: research-sources
description: >
  포스트에 쓸 자료를 조사하고 모으는 워크플로우. 논문·특허·공개 학술자료를 찾고,
  유료 논문의 합법적인 무료 전문을 확보하고, 접근이 막히면 사용자에게 요청한다.
  "이 주제 자료 좀 찾아줘", "관련 논문 조사해줘", "특허 있나 봐줘", "레퍼런스 모아줘",
  "이 논문 구할 수 있어?" 같은 요청이면 이 스킬을 쓴다.
  _references/ 로컬 자료 확인 → 공개 DB 검색 → 전문 확보 → 저장·기록 순서로 진행한다.
---

# 자료 조사 워크플로우

CLAUDE.md의 규칙을 따른다. 특히 **외부 자료의 문장을 그대로 옮기지 않고 자기 문장으로
재구성**하며, 글 말미의 "참고자료" 절에 출처를 명시한다.

아래 접근 가능 여부는 이 환경에서 실제로 확인한 결과다. 추측이 아니다.

## 0단계 — 로컬 자료 먼저

새로 찾기 전에 `_references/`에 이미 있는지 본다. 같은 논문을 두 번 받는 일이 흔하다.

```bash
ls _references/                                   # slug 별 폴더
find _references -iname "*<키워드>*"
grep -ril "<키워드>" _references --include="*.md"  # notes.md 에 남긴 메모
```

PDF는 Read 도구로 읽는다. 10쪽이 넘으면 `pages` 인자가 필수다 (`pages: "1-10"`, 한 번에
최대 20쪽). 학위논문 등 상시 참고자료는 `_references/_thesis/`에 있다.

## 1단계 — 공개 DB 검색

### OpenAlex — 첫 번째 선택

키가 필요 없고 OA 여부와 PDF 링크까지 한 번에 준다. 검색은 여기서 시작한다.

```bash
curl -s "https://api.openalex.org/works?search=<검색어>&per-page=5&mailto=eddyoptics@gmail.com" \
| python3 -c "
import sys,json
for w in json.load(sys.stdin)['results']:
    oa=(w.get('best_oa_location') or {})
    print(w.get('publication_year'), '| OA=', w.get('open_access',{}).get('is_oa'), '|', w.get('doi'))
    print('   ', (w.get('title') or '')[:80])
    if oa.get('pdf_url'): print('    PDF:', oa['pdf_url'])
"
```

`mailto` 파라미터를 붙이면 속도 제한이 완화된다.

### Crossref — 서지정보 확정

DOI를 알 때 저자·저널·권호·연도를 정확히 확인한다. 참고자료 절을 쓸 때 쓴다.

```bash
curl -s "https://api.crossref.org/works/<DOI>"
```

### arXiv — 프리프린트

물리·광학 분야는 여기 있는 경우가 많다. abs 페이지, API, PDF 모두 직접 접근된다.

```bash
curl -s "https://export.arxiv.org/api/query?search_query=all:<검색어>&max_results=5"
curl -sL -o _references/<slug>/<이름>.pdf "https://arxiv.org/pdf/<arXiv-ID>"
```

### Google Patents — 특허

**WebFetch 로 명세서와 청구항 전문까지 읽힌다.** 검증된 경로다.

```
https://patents.google.com/patent/<특허번호>/en     특정 특허
https://patents.google.com/?q=<검색어>              검색
```

계측 장비 주제는 특허가 논문보다 구현 세부를 잘 담고 있는 경우가 많다. 광학 계측이라면
Therma-Wave, KLA, Nanometrics, J.A. Woollam 같은 곳의 특허를 훑을 만하다.

### 그 밖에 열리는 곳

| 대상 | 주소 |
| --- | --- |
| 생명·의학 전문 | `https://www.ncbi.nlm.nih.gov/pmc/` |
| 오픈액세스 저널 색인 | `https://doaj.org/api/search/articles/<검색어>` |
| XPS 결합에너지 표준 DB | `https://srdata.nist.gov/xps/` |

## 2단계 — 유료 논문의 무료 전문 확보

**출판사 사이트로 직접 가지 않는다.** ACS·Elsevier·Optica는 봇을 차단한다(403). 대신
Unpaywall로 합법적인 무료 사본을 찾는다. 저자가 기관 리포지토리에 올린 판본(green OA)이
있는 경우가 많다.

```bash
curl -s "https://api.unpaywall.org/v2/<DOI>?email=eddyoptics@gmail.com" \
| python3 -c "
import sys,json
d=json.load(sys.stdin)
print('is_oa:', d.get('is_oa'), '| status:', d.get('oa_status'))
for l in (d.get('oa_locations') or []):
    print(' ', l.get('host_type'), (l.get('url_for_pdf') or l.get('url')))
"
```

`is_oa: true`면 해당 URL에서 PDF를 받아 `_references/<slug>/`에 저장하고 Read로 읽는다.

## 3단계 — 접근이 막혔을 때 사용자에게 요청

아래는 이 환경에서 **확인된 제약**이다. 우회하려고 시간을 쓰지 말고 바로 요청한다.

| 상황 | 증상 |
| --- | --- |
| 출판사 직접 접근 | `doi.org` 리졸브 후 403 (ACS 확인) |
| Optica / OSA | `opg.optica.org` 접근 불가 |
| Semantic Scholar | 429 레이트 리밋, 키 없이는 불안정 |
| Zenodo API | 403 |
| KIPRIS(한국 특허) | 홈은 열리나 검색이 JS 기반이라 결과를 못 읽음 |
| Google Scholar | 응답은 오지만 봇 차단 페이지일 수 있어 신뢰하지 않는다 |
| 유료 표준 (ISO, SEMI, ASTM) | 구매 필요 |
| 교재 | Fujiwara, Born & Wolf, Hecht 등 |

요청할 때는 **무엇이 왜 필요한지와 대안을 함께** 제시한다. 막연히 "구해달라"고 하지 않는다.

> Fujiwara 2007 5.3절이 필요합니다. 회전보상자 배치에서 Δ 부호 모호성이 어떻게
> 해소되는지가 이 글의 핵심인데, 공개 자료로는 결론만 있고 유도 과정이 없습니다.
>
> - 소장하고 계시면 해당 절을 `_references/<slug>/`에 넣어 주세요.
> - 기관 접근이 되시면 DOI 10.1002/9780470060193 로 받으실 수 있습니다.
> - 없으면 대안으로 Azzam & Bashara 1977 을 쓰겠습니다. 표기가 달라 기호 대조가 필요합니다.

**추측으로 메우지 않는다.** 확인 못 한 내용은 글에 넣지 말고, 넣어야 한다면 무엇이
확인되지 않았는지 사용자에게 알린다.

## 4단계 — 저장과 기록

받은 자료는 포스트 slug 폴더에 넣는다.

```
_references/<slug>/
  fujiwara-ch5.pdf
  us6829049-therma-wave.pdf
  notes.md
```

`_references/`는 **git 에 올라가지 않는다.** 용량 문제와 재배포 권한 때문이다. 따라서
파일을 잃으면 되찾을 방법이 필요하고, 그 역할을 `notes.md`가 한다. 자료를 받을 때마다
아래를 남긴다.

```markdown
## Survey of methods to characterize thin absorbing films (2008)
- DOI: 10.1016/j.tsf.2008.04.060
- 출처: Thin Solid Films 516(22) / OA 아님, EPFL 리포지토리 사본
- URL: http://infoscience.epfl.ch/record/205389
- 쓸 곳: 3편 다층 계산법 비교 절
- 메모: TMM 과 Rouard 의 수치 안정성 비교가 표 2에 있음
```

## 5단계 — 포스트에 반영

- 본문에는 **자기 문장으로 재구성**해 쓴다. 원문 문장을 옮기지 않는다.
- 글 말미 "참고자료" 절에 제목·저자·출처·링크를 남긴다. 여기에 링크가 남으면
  `notes.md`와 함께 이중으로 보존된다.
- **링크는 반드시 `https://`로 쓴다.** `http://` 하나면 htmlproofer 가 실패해 배포 전체가
  막힌다. 실제로 사흘간 멈춘 적이 있다.
- 특허를 인용할 때는 번호·명칭·출원인·등록일을 함께 적는다.

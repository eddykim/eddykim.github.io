#!/usr/bin/env python3
"""포스트 한 편의 사실 점검 가운데 기계로 할 수 있는 부분.

  factcheck.py refs   <포스트.md> [--json] [--no-links]   참고자료·본문 인용을 Crossref·DataCite·arXiv 와 대조
  factcheck.py pair   <한국어판.md> <영문판.md>            두 판의 수식과 숫자가 같은지 대조
  factcheck.py lint   <포스트.md>                          렌더링·배포를 깨는 함정 (수식 속 |, http://, 없는 그림 파일)
  factcheck.py claims <포스트.md>                          주장 장부의 후보 줄 (연도, 단위 붙은 숫자, 별행 수식, 인명)

refs 의 판정은 "서지정보가 맞는가"까지다. 인용한 논문이 본문의 주장을 실제로 뒷받침하는지는
스크립트가 알 수 없다. 그건 SKILL.md 4단계의 독립 검증이 한다.

네트워크는 전부 curl 로 한다. 이 맥의 Python urllib 은 SSL 인증서 검증에 실패한다.
"""
import argparse
import html
import json
import os
import re
import subprocess
import sys
import time
import unicodedata
import urllib.parse
import xml.etree.ElementTree as ET

MAILTO = "eddyoptics@gmail.com"
UA = f"eddykim-blog-factcheck/1.0 (https://eddykim.github.io; mailto:{MAILTO})"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))

DOI = re.compile(r"(?:doi\.org/|doi:\s?|DOI\s)(10\.\d{4,9}/[^\s<>\"]+)", re.I)
ARXIV = re.compile(r"(?:arxiv\.org/(?:abs|pdf)/|arXiv:\s?)(\d{4}\.\d{4,5}|[a-z\-]+(?:\.[A-Z]{2})?/\d{7})(?:v\d+)?", re.I)
URL = re.compile(r"https?://[^\s<>)\"]+")
MDLINK = re.compile(r"\[([^\]]*)\]\(((?:[^()\s]|\([^()\s]*\))+)\)")
AUTHORISH = re.compile(r"[A-Z]\.\s?(?:[A-Z]\.\s?)*[A-Z][\w'\-]+|et al|[가-힣]{2,4},")
TITLE = re.compile(r"[\"“]([^\"”]{8,}?),?[\"”]")
YEAR = re.compile(r"(?<!\d)(1[6-9]\d\d|20[0-3]\d)(?!\d)")
THESIS = re.compile(r"학위논문|[Tt]hesis|[Dd]issertation")
BOOKISH = re.compile(r"학위논문|\bed\.|\bEd\.|edition|판\)|Press|Wiley|Springer|Pergamon|McGraw|Pearson|Elsevier\)|"
                     r"Academic|Holt|Dover|CRC|교재|Handbook|\bpp\.\s|§|절\b|장\b")


# ── 공통 ──────────────────────────────────────────────────────────────────

def curl_json(url, timeout=40):
    try:
        r = subprocess.run(["curl", "-sS", "-L", "-m", str(timeout), "-A", UA,
                            "-w", "\n%{http_code}", url], capture_output=True, text=True)
    except OSError as e:
        return None, f"curl 실행 실패: {e}"
    body, _, code = r.stdout.rpartition("\n")
    if r.returncode != 0:
        return None, f"네트워크 오류 ({r.stderr.strip()[:80]})"
    if code != "200":
        return None, code
    try:
        return json.loads(body), "200"
    except json.JSONDecodeError:
        return body, "200"


def http_status(url):
    r = subprocess.run(["curl", "-sS", "-L", "-m", "30", "-A", UA, "-o", "/dev/null",
                        "-w", "%{http_code}", url], capture_output=True, text=True)
    return r.stdout.strip() or "ERR"


def norm(s):
    s = html.unescape(re.sub(r"<[^>]+>", " ", s or ""))
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    return s.lower()


def tokens(s):
    return {t for t in re.findall(r"[a-z0-9]+", norm(s)) if len(t) > 1 or t.isdigit()}


def title_sim(a, b):
    A, B = tokens(a), tokens(b)
    if not A or not B:
        return 0.0
    return len(A & B) / max(len(A), len(B))


def read(path):
    with open(path, encoding="utf-8") as f:
        return f.read().split("\n")


def body_lines(lines):
    """front matter 와 코드 블록을 뺀 (줄 번호, 내용)."""
    fm = bool(lines) and lines[0].strip() == "---"
    fence = None
    for n, line in enumerate(lines, 1):
        s = line.strip()
        if fm:
            if n > 1 and s == "---":
                fm = False
            continue
        if fence:
            if s.startswith(fence):
                fence = None
            continue
        m = re.match(r"(```+|~~~+)", s)
        if m:
            fence = m.group(1)
            continue
        yield n, line


# ── refs ──────────────────────────────────────────────────────────────────

def ref_items(lines):
    """참고자료 절의 항목. (줄 번호, 한 줄로 합친 텍스트)."""
    items, inside = [], False
    for n, line in body_lines(lines):
        if re.match(r"^##\s+(참고자료|참고 문헌|References|Further reading)\s*$", line.strip()):
            inside = True
            continue
        if inside and re.match(r"^#{1,2}\s", line):
            break
        if not inside or not line.strip():
            continue
        if re.match(r"^\s*(?:[-*]|\d+\.)\s+", line):
            items.append([n, re.sub(r"^\s*(?:[-*]|\d+\.)\s+", "", line).strip()])
        elif items and line.startswith((" ", "\t")):
            items[-1][1] += " " + line.strip()
    return items


def citation_part(text):
    """설명(— 뒤)과 URL 을 뺀 서지 부분."""
    t = re.split(r"\s[—–]\s", text, maxsplit=1)[0]
    t = MDLINK.sub(lambda m: m.group(1) if m.start() and m.string[m.start() - 1] in "\"“"
                   else f'"{m.group(1)}"' if m.group(1) else "", t)
    t = re.sub(r"<https?://[^>]+>", "", t)
    t = URL.sub("", t)
    return re.sub(r"[*<>]", "", t).strip()


def crossref_record(msg):
    def year_of(key):
        try:
            return msg[key]["date-parts"][0][0]
        except (KeyError, IndexError, TypeError):
            return None
    years = {y for y in (year_of(k) for k in ("published-print", "published-online", "issued", "posted")) if y}
    authors = [a.get("family") or a.get("name") or "" for a in msg.get("author", [])]
    return {
        "source": "Crossref",
        "doi": msg.get("DOI"),
        "title": " ".join(html.unescape(re.sub(r"<[^>]+>", "", ": ".join(
            [(msg.get("title") or [""])[0]] + (msg.get("subtitle") or [])[:1]))).split()),
        "authors": authors,
        "container": html.unescape((msg.get("short-container-title") or msg.get("container-title") or [""])[0]),
        "volume": msg.get("volume"),
        "page": msg.get("page") or msg.get("article-number"),
        "years": sorted(years),
        "type": msg.get("type"),
    }


def clean_doi(doi):
    doi = urllib.parse.unquote(doi).rstrip(".,;>")
    while doi.endswith(")") and doi.count(")") > doi.count("("):
        doi = doi[:-1]
    return doi


def lookup_doi(doi):
    doi = clean_doi(doi)
    q = urllib.parse.quote(doi, safe="/")
    data, code = curl_json(f"https://api.crossref.org/works/{q}?mailto={MAILTO}")
    if isinstance(data, dict) and "message" in data:
        return crossref_record(data["message"]), None
    data2, code2 = curl_json(f"https://api.datacite.org/dois/{q}")
    if isinstance(data2, dict) and "data" in data2:
        a = data2["data"]["attributes"]
        return {
            "source": "DataCite",
            "doi": doi,
            "title": (a.get("titles") or [{}])[0].get("title", ""),
            "authors": [c.get("familyName") or c.get("name", "") for c in a.get("creators", [])],
            "container": a.get("publisher") if isinstance(a.get("publisher"), str) else "",
            "volume": None, "page": None,
            "years": [a.get("publicationYear")] if a.get("publicationYear") else [],
            "type": (a.get("types") or {}).get("resourceTypeGeneral"),
        }, None
    if code in ("404",) and code2 in ("404",):
        return None, "DEAD_DOI"
    return None, f"NET_ERROR (Crossref {code}, DataCite {code2})"


def lookup_arxiv(aid):
    raw = subprocess.run(["curl", "-sS", "-m", "40", "-A", UA,
                          f"https://export.arxiv.org/api/query?id_list={aid}"],
                         capture_output=True, text=True).stdout
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return None, "NET_ERROR (arXiv)"
    ns = {"a": "http://www.w3.org/2005/Atom"}
    e = root.find("a:entry", ns)
    if e is None or e.find("a:title", ns) is None:
        return None, "NOT_FOUND"
    title = " ".join(e.find("a:title", ns).text.split())
    if title.lower() == "error":
        return None, "NOT_FOUND"
    authors = [a.find("a:name", ns).text.split()[-1] for a in e.findall("a:author", ns)]
    year = int(e.find("a:published", ns).text[:4])
    return {"source": "arXiv", "doi": f"arXiv:{aid}", "title": title, "authors": authors,
            "container": "arXiv", "volume": None, "page": None, "years": [year], "type": "preprint"}, None


def lookup_biblio(cite):
    q = urllib.parse.quote(cite[:300])
    data, code = curl_json(f"https://api.crossref.org/works?query.bibliographic={q}&rows=3&mailto={MAILTO}")
    if not isinstance(data, dict):
        return [], f"NET_ERROR (Crossref {code})"
    return [crossref_record(m) for m in data["message"].get("items", [])], None


def compare(line, rec):
    """블로그 표기와 레코드의 필드별 대조. (판정 코드 목록, 메모 목록, 일치한 필드 수).

    어긋난 필드가 없다는 것만으로는 OK 가 아니다. 대조할 필드가 없으면 아무것도 안 어긋나기
    때문이다. 일치한 필드 수를 함께 돌려주고, 호출하는 쪽이 그것으로 OK 를 준다."""
    codes, notes, hits = [], [], 0
    cite = citation_part(line)
    n_line = norm(cite)
    m = TITLE.search(cite)
    if m and rec["title"]:
        s = title_sim(m.group(1), rec["title"])
        hits += 2 if s >= 0.75 else 0
        if s < 0.4:
            codes.append("TITLE_MISMATCH")
        elif s < 0.75:
            codes.append("TITLE_DIFF")
            notes.append(f"제목 유사도 {s:.2f}")
    if rec["authors"] and AUTHORISH.search(cite):
        first = (norm(rec["authors"][0]).split() or [""])[-1]
        if first and first in n_line:
            hits += 1
        elif first:
            codes.append("AUTHOR_MISMATCH")
            notes.append(f"제1저자 레코드={rec['authors'][0]}")
    yrs = [int(y) for y in YEAR.findall(cite)]
    if rec["years"] and yrs and set(yrs) & set(rec["years"]):
        hits += 1
    elif rec["years"] and yrs and rec["source"] == "arXiv":
        codes.append("YEAR_DIFF")
        notes.append(f"연도 표기={yrs}, arXiv 첫 게시={rec['years']} — 학회·저널 연도라면 맞을 수 있다")
    elif rec["years"] and yrs:
        codes.append("YEAR_MISMATCH")
        notes.append(f"연도 표기={yrs} 레코드={rec['years']}")
    vol = str(rec.get("volume") or "")
    if vol and re.search(rf"(?<![\d.]){re.escape(vol)}(?![\d.])", cite):
        hits += 1
    elif vol:
        codes.append("VOLUME_MISMATCH")
        notes.append(f"권 레코드={vol}")
    page = str(rec.get("page") or "").split("-")[0].split("–")[0]
    others = [x for x in re.findall(r"(?<![\w.])\d+(?![\w.])", cite) if x not in map(str, yrs) and x != vol]
    if page and re.search(rf"(?<![\d.]){re.escape(page)}(?![\d.])", cite):
        hits += 1
    elif page and others:
        codes.append("PAGE_MISMATCH")
        notes.append(f"쪽 레코드={rec['page']}")
    if not codes and hits < 2:
        codes.append("UNVERIFIED")
        notes.append("대조할 수 있는 필드가 부족하다")
    return codes or ["OK"], notes, hits


def match_score(codes):
    bad = {"TITLE_MISMATCH": 3, "AUTHOR_MISMATCH": 2, "YEAR_MISMATCH": 1, "YEAR_DIFF": 0.5,
           "VOLUME_MISMATCH": 1, "PAGE_MISMATCH": 1, "TITLE_DIFF": 0.5}
    return sum(bad.get(c, 0) for c in codes)


def check_item(n, text, check_links):
    res = {"line": n, "text": text, "codes": [], "notes": [], "record": None}
    dois = DOI.findall(text)
    arx = ARXIV.findall(text)
    if dois or arx:
        rec, err = lookup_doi(dois[0]) if dois else lookup_arxiv(arx[0])
        if err:
            res["codes"] = [err]
        else:
            res["record"] = rec
            res["codes"], res["notes"], _ = compare(text, rec)
    else:
        cite = citation_part(text)
        if THESIS.search(text) or (BOOKISH.search(text) and not TITLE.search(cite)):
            res["codes"] = ["MANUAL"]
            res["notes"].append("교재·학위논문 — 소장본이나 Research Library 로 판·연도·장절을 확인")
        elif URL.search(text) and not re.search(r"\d+\s?,\s?\d+", cite):
            res["codes"] = ["MANUAL"]
            res["notes"].append("웹 자료 — 링크를 열어 제목·발행처를 확인")
        else:
            cands, err = lookup_biblio(cite)
            if err:
                res["codes"] = [err]
            else:
                best = None
                for rec in cands:
                    codes, notes, hits = compare(text, rec)
                    sc = match_score(codes) - hits
                    if best is None or sc < best[0]:
                        best = (sc, rec, codes, notes, hits)
                # 서지 검색은 늘 무언가를 돌려준다. 제목이 맞거나 세 필드 이상 맞아야 같은 논문으로 본다.
                # 서지 검색이 고른 레코드의 제목이 다르면, 대개 블로그 제목이 틀린 게 아니라 다른 논문을
                # 집어 온 것이다. 권과 쪽까지 맞을 때만 "같은 논문인데 제목이 틀렸다"고 본다.
                if best and "TITLE_MISMATCH" in best[2] and (
                        "VOLUME_MISMATCH" in best[2] or "PAGE_MISMATCH" in best[2]
                        or not (best[1].get("volume") and best[1].get("page"))):
                    best = (best[0], best[1], best[2], best[3], 0)
                if best is None or best[4] < 3:
                    res["codes"] = ["NOT_FOUND"]
                    res["notes"].append("Crossref 서지 검색에서 맞는 논문을 못 찾음 — DB에 없는 것일 수도, "
                                        "지어낸 인용일 수도 있다. 다른 경로로 존재부터 확인")
                    if best:
                        res["record"] = best[1]
                else:
                    _, res["record"], res["codes"], res["notes"], _ = best
                    res["notes"].append("DOI 없음 — 서지 검색으로 찾은 레코드와 대조")
    if len(URL.findall(text)) > 1 and res["codes"] != ["OK"]:
        res["notes"].append("링크가 둘 이상 — 여러 문헌을 한 항목에 묶었을 수 있다. 제목과 DOI가 다른 문헌을 가리키는지 볼 것")
    if check_links:
        for u in URL.findall(text):
            if "doi.org/" in u or "arxiv.org/" in u:
                continue
            code = http_status(u.rstrip(">.,"))
            if code != "200":
                res["notes"].append(f"링크 {code}: {u}" + (" (봇 차단일 수 있음)" if code in ("403", "429") else ""))
                if code in ("404", "410", "000", "ERR"):
                    res["codes"].append("DEAD_LINK")
    time.sleep(0.2)
    return res


def cmd_refs(a):
    lines = read(a.post)
    items = ref_items(lines)
    ref_lines = {n for n, _ in items}
    results = [check_item(n, t, not a.no_links) for n, t in items]
    inline = []
    for n, line in body_lines(lines):
        if n in ref_lines:
            continue
        for d in DOI.findall(line):
            rec, err = lookup_doi(d)
            inline.append({"line": n, "id": d, "codes": [err] if err else ["EXISTS"], "record": rec})
        for x in ARXIV.findall(line):
            rec, err = lookup_arxiv(x)
            inline.append({"line": n, "id": f"arXiv:{x}", "codes": [err] if err else ["EXISTS"], "record": rec})
    if a.json:
        print(json.dumps({"refs": results, "inline": inline}, ensure_ascii=False, indent=1))
        return
    if not items:
        print("참고자료 절을 못 찾았다 (## 참고자료 / ## References).")
    for r in results:
        flag = "  " if r["codes"] == ["OK"] else "!!" if any(c.endswith("MISMATCH") or c.startswith("DEAD") for c in r["codes"]) else " ?"
        print(f"{flag} L{r['line']:<4} {' '.join(r['codes'])}")
        print(f"        표기: {citation_part(r['text'])[:150]}")
        rec = r["record"]
        if rec and r["codes"] != ["OK"]:
            au = ", ".join(rec["authors"][:3]) + (" 외" if len(rec["authors"]) > 3 else "")
            print(f"        {rec['source']}: {au} | {rec['title'][:90]} | {rec['container']} "
                  f"{rec.get('volume') or ''} {rec.get('page') or ''} {rec['years']} | {rec['doi']}")
        for note in r["notes"]:
            print(f"        - {note}")
    for i in inline:
        rec = i["record"]
        t = f"{rec['title'][:90]} ({rec['years']})" if rec else ""
        print(f"   L{i['line']:<4} 본문 링크 {i['id']} {' '.join(i['codes'])} {t}")
    bad = sum(1 for r in results if r["codes"] != ["OK"])
    print(f"\n항목 {len(results)}개, 확인 필요 {bad}개, 본문 링크 {len(inline)}개")


# ── pair ──────────────────────────────────────────────────────────────────

MATH = re.compile(r"\$\$(.+?)\$\$|\$(.+?)\$", re.S)
NUM = re.compile(r"(?<![A-Za-z0-9_.])[-−]?\d+(?:[.,]\d+)?(?:[eE][-−]?\d+)?(?![A-Za-z0-9_])")


def strip_links(s):
    s = re.sub(r"!\[[^\]]*\]\([^)]*\)|<img[^>]*>", " ", s)
    s = re.sub(r"\]\([^)]*\)", "]", s)
    return URL.sub(" ", s)


def split_post(path):
    text = "\n".join(l for _, l in body_lines(read(path)))
    text = strip_links(text)
    maths = [re.sub(r"\s+", "", (m.group(1) or m.group(2))) for m in MATH.finditer(text)]
    prose = MATH.sub(" ", text)
    nums = [n.replace("−", "-").replace(",", "") for n in NUM.findall(prose)]
    return maths, nums


def multiset_diff(a, b):
    from collections import Counter
    ca, cb = Counter(a), Counter(b)
    return ca - cb, cb - ca


def cmd_pair(a):
    mk, nk = split_post(a.ko)
    me, ne = split_post(a.en)
    only_k, only_e = multiset_diff(mk, me)
    print(f"수식 {len(mk)} / {len(me)}개")
    for m, c in only_k.items():
        print(f"  한국어판에만: {m[:120]}" + (f" ×{c}" if c > 1 else ""))
    for m, c in only_e.items():
        print(f"  영문판에만:  {m[:120]}" + (f" ×{c}" if c > 1 else ""))
    only_k, only_e = multiset_diff(nk, ne)
    print(f"\n본문 숫자 {len(nk)} / {len(ne)}개 (편 번호·서수 차이는 무시해도 된다)")
    if only_k:
        print("  한국어판에만:", " ".join(f"{k}×{v}" if v > 1 else k for k, v in only_k.items()))
    if only_e:
        print("  영문판에만: ", " ".join(f"{k}×{v}" if v > 1 else k for k, v in only_e.items()))
    if not any([only_k, only_e]) and mk == me:
        print("차이 없음")


# ── lint ──────────────────────────────────────────────────────────────────

def cmd_lint(a):
    lines = read(a.post)
    found = 0
    for n, line in body_lines(lines):
        for m in MATH.finditer(line):
            if "|" in m.group(0):
                found += 1
                print(f"L{n}: 수식 속 '|' — kramdown 이 표로 렌더한다. \\lvert \\rvert 로: {m.group(0)[:80]}")
        for m in re.finditer(r"(\]\(|src=\"|href=\"|<)http://(?!127\.0\.0\.1|0\.0\.0\.0|localhost)", line):
            found += 1
            print(f"L{n}: http:// 링크 — htmlproofer 가 배포를 막는다")
        for m in re.finditer(r"(?:src=\"|\]\()(/assets/[^\")\s]+)", line):
            if not os.path.exists(os.path.join(REPO, m.group(1).lstrip("/"))):
                found += 1
                print(f"L{n}: 없는 파일 {m.group(1)}")
    text = "\n".join(l for _, l in body_lines(lines))
    if text.count("$$") % 2:
        found += 1
        print("$$ 개수가 홀수 — 닫히지 않은 별행 수식이 있다")
    print(f"\n{found}건" if found else "문제 없음")


# ── claims ────────────────────────────────────────────────────────────────

UNIT = re.compile(r"\d(?:[\d.,]*)\s?(?:nm|µm|um|mm|cm|m|pm|Å|eV|keV|kV|V|mrad|rad|°|도|deg|%|Hz|kHz|MHz|K|s|ms|µs|ns|fs|W|mW|J|T|Pa|배)(?![A-Za-z가-힣])")
EPONYM = re.compile(r"[A-Z][a-z]+(?:[-–][A-Z][a-z]+)?(?:의| 법칙| 식| 방정식| 원리| 조건| 각| 근사| 변환|'s)|"
                    r"[가-힣]{2,}(?:가|이|는|은) \d{4}년")


CITE = re.compile(r"\[\d+(?:\s?[,–-]\s?\d+)*\]")
HISTORY = re.compile(r"처음|최초|시초|시조|창시|고안|발견|명명|이름을 따|first|introduced|coined|named after", re.I)


def cmd_claims(a):
    """주장 장부를 채울 후보 줄. 놓치는 것이 있으므로 본문 정독을 대신하지 않는다."""
    lines = read(a.post)
    for n, line in body_lines(lines):
        s = line.strip()
        if re.match(r"^##\s+(참고자료|참고 문헌|References)", s):
            break
        if not s or s.startswith(("<img", "![")):
            continue
        tags = []
        if s.startswith("$$") or s.endswith("$$"):
            tags.append("식")
        if YEAR.search(s):
            tags.append("연도")
        if HISTORY.search(s):
            tags.append("역사")
        if CITE.search(s):
            tags.append("인용")
        if UNIT.search(s):
            tags.append("수치")
        if EPONYM.search(s):
            tags.append("인명")
        if s.startswith("_") and re.match(r"_(그림|Fig)", s):
            tags.append("캡션")
        if tags:
            print(f"L{n:<4} [{','.join(tags)}] {s[:160]}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("refs"); r.add_argument("post"); r.add_argument("--json", action="store_true")
    r.add_argument("--no-links", action="store_true"); r.set_defaults(f=cmd_refs)
    q = sub.add_parser("pair"); q.add_argument("ko"); q.add_argument("en"); q.set_defaults(f=cmd_pair)
    l = sub.add_parser("lint"); l.add_argument("post"); l.set_defaults(f=cmd_lint)
    c = sub.add_parser("claims"); c.add_argument("post"); c.set_defaults(f=cmd_claims)
    a = p.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()

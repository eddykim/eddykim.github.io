#!/usr/bin/env python3
"""포스트 한 편의 사실 점검 가운데 기계로 할 수 있는 부분.

  factcheck.py refs   <포스트.md> [--json] [--no-links]   참고자료·본문 인용을 Crossref·DataCite·arXiv 와 대조 (공통 명령에 위임)
  factcheck.py pair   <한국어판.md> <영문판.md>            두 판의 수식과 숫자가 같은지 대조
  factcheck.py lint   <포스트.md>                          렌더링·배포를 깨는 함정 (수식 속 |, http://, 없는 그림 파일)
  factcheck.py claims <포스트.md>                          주장 장부의 후보 줄 (연도, 단위 붙은 숫자, 별행 수식, 인명)

refs 의 판정은 "서지정보가 맞는가"까지다. 인용한 논문이 본문의 주장을 실제로 뒷받침하는지는
스크립트가 알 수 없다. 그건 SKILL.md 4단계의 독립 검증이 한다. refs 의 구현은 공통 팩트체커
(repo_kms/tools/factcheck.py)로 옮겼고 여기서는 그 명령을 부른다.
"""
import argparse
import os
import re
import subprocess
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))

URL = re.compile(r"https?://[^\s<>)\"]+")
MDLINK = re.compile(r"\[([^\]]*)\]\(((?:[^()\s]|\([^()\s]*\))+)\)")
YEAR = re.compile(r"(?<!\d)(1[6-9]\d\d|20[0-3]\d)(?!\d)")


# ── 공통 ──────────────────────────────────────────────────────────────────

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

KMS = os.path.normpath(os.path.join(REPO, "..", "repo_kms", "tools", "factcheck.py"))


def cmd_refs(a):
    """서지 대조는 공통 팩트체커(repo_kms/tools/factcheck.py refs)가 한다. 인자를 그대로 넘긴다."""
    if not os.path.exists(KMS):
        sys.exit(f"공통 팩트체커를 찾지 못했다: {KMS}")
    args = ["python3", KMS, "refs", a.post] + (["--json"] if a.json else []) + (["--no-links"] if a.no_links else [])
    sys.exit(subprocess.run(args).returncode)


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

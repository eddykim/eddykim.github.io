#!/usr/bin/env python3
"""포스트의 그림 번호 · 캡션 · 출처 표기를 점검하고, 그림 번호를 등장 순서대로 다시 매긴다.

  figures.py check    <포스트.md> [<짝 포스트.md>]   한·영 두 파일을 주면 그림 순서도 대조한다
  figures.py renumber <포스트.md> [--apply]          --apply 없이 돌리면 바뀔 곳만 보여준다

그림은 이 블로그의 관례를 따른다고 가정한다.

  <img src="/assets/img/posts/<slug>/fig1-....png" alt="..." width="600">
  _그림1. 캡션._            (영문판은 _Fig 1. 캡션._ 또는 _Figure 1. 캡션._)

웹 이미지는 파일명이 ext- 로 시작한다. 그 캡션에는 출처와 라이선스가 있어야 한다.
"""
import argparse
import os
import re
import sys

IMG = re.compile(r'<img\s[^>]*?src="([^"]+)"|!\[[^\]]*\]\(([^)\s]+)')
CAP = re.compile(r"^_(그림|Figure|Fig\.?)\s?(\d+)\.")
REF = re.compile(r"(그림|Figure|Fig\.?)(\s?)(\d+)")
# "그림 2, 3", "Figures 3, 4 and 5", "그림 가운데 3·4·5" 처럼 번호를 묶어 부르는 줄.
# 첫 번호만 REF 에 걸리거나 아예 안 걸리므로 스크립트가 고치지 않고 알리기만 한다.
GROUP = re.compile(r"(그림|Fig)[^.$\n]{0,15}?(?<![\d.])\d{1,2}\s?(·|,|~|–|-|/|와|과|및|and|to)\s?\d{1,2}(?![\d.%°])")
FOREIGN = re.compile(r"et al|원문|논문의|of \[|in \[")
SOURCE = re.compile(r"출처:|Source:")
LICENSE = re.compile(r"CC0|CC BY|퍼블릭 도메인|[Pp]ublic domain")
PLACEHOLDER = re.compile(r"<작성자 확인 필요>|<author: check file page>")
HTTP = re.compile(r'(\]\(|src="|href=")http://(?!127\.0\.0\.1|0\.0\.0\.0|localhost)')
MAX_BYTES = 800_000


def body_lines(lines):
    """front matter 와 코드 블록을 뺀 (줄 번호, 내용)."""
    fm = lines[0].strip() == "---" if lines else False
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
        if s.startswith(("```", "~~~")):
            fence = s[:3]
            continue
        yield n, line


def is_group(line):
    return bool(GROUP.search(line)) and not CAP.match(line.strip())


def repo_root(path):
    d = os.path.dirname(os.path.abspath(path))
    while d != os.path.dirname(d):
        if os.path.exists(os.path.join(d, "_config.yml")):
            return d
        d = os.path.dirname(d)
    return os.getcwd()


def parse(path):
    lines = open(path, encoding="utf-8").read().split("\n")
    body = list(body_lines(lines))
    figs, caps, refs = [], [], []
    for k, (n, line) in enumerate(body):
        m = IMG.search(line)
        if m:
            src = m.group(1) or m.group(2)
            cap = None
            for n2, nxt in body[k + 1:k + 3]:     # 캡션은 바로 다음 줄(빈 줄 하나까지 허용)
                c = CAP.match(nxt.strip())
                if c:
                    cap = (n2, int(c.group(2)), nxt, c.group(1))
                    break
                if nxt.strip():
                    break
            figs.append(dict(line=n, src=src, cap=cap))
        c = CAP.match(line.strip())
        if c:
            caps.append((n, int(c.group(2))))
        for r in REF.finditer(line):
            if c and r.start() == 1:              # 캡션 자신의 번호
                continue
            refs.append((n, int(r.group(3)), line))
    return lines, figs, caps, refs


def check_one(path, problems):
    lines, figs, caps, refs = parse(path)
    root = repo_root(path)
    nums = [num for _, num in caps]
    name = os.path.basename(path)
    print(f"\n{name}")

    for i, f in enumerate(figs, 1):
        src = f["src"]
        web = os.path.basename(src).startswith("ext-")
        tag = "웹" if web else "생성"
        cap = f["cap"]
        info = []
        if src.startswith("http"):
            problems.append(f"{name} L{f['line']}: 외부 이미지를 직접 링크한다 — 내려받아 assets 에 둔다")
            info.append("핫링크")
        else:
            fp = os.path.join(root, src.lstrip("/"))
            if not os.path.exists(fp):
                problems.append(f"{name} L{f['line']}: 이미지 파일 없음 {src}")
                info.append("파일 없음")
            else:
                size = os.path.getsize(fp)
                info.append(f"{size // 1000} KB")
                if web and size > MAX_BYTES:
                    problems.append(f"{name} L{f['line']}: {size // 1000} KB — 800KB 이하로 줄인다")
        if not cap:
            problems.append(f"{name} L{f['line']}: 캡션 없음 ({os.path.basename(src)})")
        elif web:
            text = cap[2]
            if not SOURCE.search(text):
                problems.append(f"{name} L{cap[0]}: 웹 이미지 캡션에 '출처:'/'Source:' 가 없다")
            if not LICENSE.search(text):
                problems.append(f"{name} L{cap[0]}: 웹 이미지 캡션에 라이선스 표기가 없다")
            if PLACEHOLDER.search(text):
                problems.append(f"{name} L{cap[0]}: 작성자 자리표시가 남아 있다 — 파일 페이지에서 채운다")
        label = f"{cap[3]}{cap[1]}" if cap else "캡션 없음"
        print(f"  {i:>2}. {label:<8} [{tag}] {os.path.basename(src)}  {' '.join(info)}")

    if nums != list(range(1, len(nums) + 1)):
        problems.append(f"{name}: 캡션 번호가 1부터 차례대로가 아니다 {nums} — renumber 를 돌린다")
    for n, num, line in refs:
        if num not in nums:
            warn = " (외부 문헌 그림일 수 있음)" if FOREIGN.search(line) else ""
            problems.append(f"{name} L{n}: 없는 그림 {num}번을 부른다{warn}")
    for n, line in body_lines(lines):
        if is_group(line):
            print(f"  ⚠ L{n}: 묶음 참조 — 번호를 바꿨다면 직접 확인: {line.strip()[:70]}")
        if HTTP.search(line):
            problems.append(f"{name} L{n}: http:// 링크 — htmlproofer 가 배포를 막는다")
    return figs


def cmd_check(a):
    problems = []
    for p in a.posts:
        if not os.path.exists(p):
            print(f"{p} 없음 — 건너뛴다 (영문판을 아직 안 썼다면 정상)")
    figs = [check_one(p, problems) for p in a.posts if os.path.exists(p)]
    if len(figs) == 2:
        seq = [[os.path.basename(f["src"]) for f in fs] for fs in figs]
        if seq[0] != seq[1]:
            problems.append("한·영 그림 순서가 다르다:\n    " + "\n    ".join(" | ".join(s) for s in seq))
        web = [[f["src"] for f in fs if os.path.basename(f["src"]).startswith("ext-")] for fs in figs]
        if web[0] != web[1]:
            problems.append("한·영이 서로 다른 웹 이미지 경로를 쓴다 (웹 이미지는 en/ 사본 없이 같은 파일을 쓴다)")
    print()
    if problems:
        print("문제:")
        for p in problems:
            print("  - " + p)
        sys.exit(1)
    print("이상 없음")


def cmd_renumber(a):
    lines, figs, caps, refs = parse(a.post)
    nums = [num for _, num in caps]
    dup = sorted({x for x in nums if nums.count(x) > 1})
    if dup:
        sys.exit(f"캡션 번호가 겹친다 {dup}. 새 그림에는 기존과 겹치지 않는 임시 번호(101, 102 …)를 붙인다.")
    mapping = {old: i for i, old in enumerate(nums, 1)}
    if all(k == v for k, v in mapping.items()):
        print("이미 차례대로다. 바꿀 것 없음.")
        return
    print("번호 대응: " + ", ".join(f"{k}→{v}" for k, v in mapping.items() if k != v))

    body = {n for n, _ in body_lines(lines)}
    out = list(lines)
    for n in sorted(body):
        line = lines[n - 1]
        changes = []

        def sub(m):
            old = int(m.group(3))
            new = mapping.get(old, old)
            if new != old:
                changes.append(f"{m.group(0)} → {m.group(1)}{m.group(2)}{new}")
            return f"{m.group(1)}{m.group(2)}{new}"

        new_line = REF.sub(sub, line)
        if changes:
            ctx = line.strip()[:70]
            flag = "  ⚠ 외부 문헌 그림일 수 있음 — 확인" if FOREIGN.search(line) else ""
            print(f"  L{n}: {', '.join(changes)}{flag}\n        {ctx}")
            out[n - 1] = new_line
        for r in REF.finditer(line):
            if int(r.group(3)) not in mapping and not CAP.match(line.strip()):
                print(f"  L{n}: {r.group(0)} 는 캡션에 없는 번호라 그대로 둔다")

    groups = [(n, lines[n - 1]) for n in sorted(body) if is_group(lines[n - 1])]
    if groups:
        print("\n묶음 참조 — 스크립트가 고치지 않는다. 위 번호 대응대로 직접 고친다:")
        for n, line in groups:
            print(f"  L{n}: {line.strip()[:90]}")

    if a.apply:
        open(a.post, "w", encoding="utf-8").write("\n".join(out))
        print("적용함. check 로 다시 확인한다.")
    else:
        print("\n미리보기다. 맞으면 --apply 를 붙여 다시 돌린다.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("posts", nargs="+")
    c.set_defaults(fn=cmd_check)
    r = sub.add_parser("renumber")
    r.add_argument("post")
    r.add_argument("--apply", action="store_true")
    r.set_defaults(fn=cmd_renumber)
    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()

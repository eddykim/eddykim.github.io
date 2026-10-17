"""빌드된 사이트의 포스트에서 수식을 꺼내 렌더링 문제를 찾는다.
사용: python3 mathcheck.py <_site 경로> [출력.md]
- kramdown이 수식 안을 마크다운으로 해석한 흔적(<em>, <strong>, <a>, <td> 등 태그가 수식 안에 섞임)
- 짝이 안 맞는 $ (문단 단위)
- 표 칸에 쪼개진 수식
- MathJax 파서 오류 (정의 안 된 명령, 괄호 불일치 등)"""
import html, json, re, subprocess, sys
from pathlib import Path

SITE = Path(sys.argv[1])
OUT = Path(sys.argv[2]) if len(sys.argv) > 2 else None
PARSE = Path(__file__).parent / "mjx/parse.js"

ART = re.compile(r'<div class="content">(.*?)</div>\s*<div class="post-tail', re.S)
STRIP = re.compile(r"<(script|pre|code|style)\b.*?</\1>", re.S)
BLOCK = re.compile(r"(?<!\\)\\\[(.+?)\\\]", re.S)
INL = re.compile(r"(?<!\\)\\\((.+?)\\\)", re.S)
DD = re.compile(r"(?<![\\$])\$\$(.+?)\$\$", re.S)
SD = re.compile(r"(?<![\\$])\$(?!\$)(.+?)(?<![\\$])\$(?!\$)", re.S)
TAG = re.compile(r"</?(em|strong|a|sup|sub|del|td|tr|table|p|br|li|span)\b[^>]*>")

items, issues = [], []


def ctx(s):
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", s))[:120]


for f in sorted(SITE.glob("**/posts/*/index.html")):
    page = f.relative_to(SITE).parent.as_posix()
    m = ART.search(f.read_text(encoding="utf-8"))
    if not m:
        continue
    body = STRIP.sub("", m.group(1))
    found = []
    rest = body
    for rx, disp in ((BLOCK, True), (DD, True), (INL, False)):
        for mm in rx.finditer(rest):
            found.append((mm.group(1), disp))
        rest = rx.sub(" ", rest)
    for para in re.split(r"</p>|</li>|</td>|</th>|<br\s*/?>|</h\d>|</figcaption>", rest):
        text = re.sub(r"<[^>]+>", "", para)
        if (text.count("$") - text.count(r"\$")) % 2:
            issues.append((page, "짝 안 맞는 $", ctx(para)))
        for mm in SD.finditer(para):
            found.append((mm.group(1), False))
    if "<td" in body:
        for cell in re.findall(r"<td[^>]*>(.*?)</td>", body, re.S):
            t = re.sub(r"<[^>]+>", "", cell)
            if (t.count("$") - t.count(r"\$")) % 2:
                issues.append((page, "표 칸에서 쪼개진 수식", ctx(cell)))
    for tex, disp in found:
        if TAG.search(tex):
            issues.append((page, "수식 안에 HTML 태그 (마크다운으로 해석됨)", ctx(tex) + " ‖ " + TAG.search(tex).group(0)))
            continue
        k = len(items)
        items.append({"k": k, "tex": html.unescape(tex), "display": disp, "page": page})

r = subprocess.run(["node", str(PARSE)], input=json.dumps([{k: v for k, v in it.items() if k != "page"} for it in items]),
                   capture_output=True, text=True, cwd=PARSE.parent)
if r.returncode:
    print(r.stderr[-2000:]); sys.exit(1)
for e in json.loads(r.stdout):
    it = items[e["k"]]
    issues.append((it["page"], "MathJax 오류: " + e["error"], it["tex"][:120].replace("\n", " ")))

lines = [f"# 수식 렌더링 점검 — 수식 {len(items)}개, 문제 {len(issues)}건", ""]
cur = None
for page, kind, s in sorted(issues):
    if page != cur:
        lines += ["", f"## {page}"]; cur = page
    lines.append(f"- **{kind}** — `{s}`")
txt = "\n".join(lines) + "\n"
(OUT.write_text(txt) if OUT else print(txt))
print(f"수식 {len(items)}개, 문제 {len(issues)}건, 쪽 {len({p for p,_,_ in issues})}개")

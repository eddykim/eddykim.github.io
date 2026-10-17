r"""원문의 인라인 수식($...$)이 빌드 결과에 글자 그대로 남아 있는지 본다.
kramdown은 $...$를 수식으로 보지 않아서 \{ \} \\ \| \_ \* 같은 이스케이프를 풀거나 _ * 를 강조로 바꾼다.
사용: python3 mathdiff.py <저장소> <_site> [출력.md]"""
import html, re, sys
from pathlib import Path

REPO, SITE = Path(sys.argv[1]), Path(sys.argv[2])
OUT = Path(sys.argv[3]) if len(sys.argv) > 3 else None
FENCE = re.compile(r"^```.*?^```", re.S | re.M)
CODE = re.compile(r"`[^`\n]+`")
DD = re.compile(r"\$\$.+?\$\$", re.S)
SD = re.compile(r"(?<![\\$])\$(?!\$)([^$\n]+?)(?<![\\$])\$(?!\$)")
ART = re.compile(r'<div class="content">(.*?)</div>\s*<div class="post-tail', re.S)
TAGS = re.compile(r"<[^>]+>")

rows = []
for md in sorted(REPO.glob("_posts/*.md")):
    name = md.stem
    slug = re.sub(r"^\d{4}-\d{2}-\d{2}-", "", name)
    en = slug.endswith("-en")
    slug = slug[:-3] if en else slug
    page = SITE / ("en/posts" if en else "posts") / slug / "index.html"
    if not page.exists():
        continue
    m = ART.search(page.read_text(encoding="utf-8"))
    if not m:
        continue
    built = html.unescape(m.group(1))
    plain = html.unescape(TAGS.sub("", m.group(1)))
    src = md.read_text(encoding="utf-8").split("---", 2)[-1]
    src = DD.sub(" ", CODE.sub(" ", FENCE.sub(" ", src)))
    for ln, line in enumerate(md.read_text(encoding="utf-8").splitlines(), 1):
        pass
    seen = set()
    for mm in SD.finditer(src):
        tex = mm.group(1)
        if tex in seen:
            continue
        seen.add(tex)
        if f"${tex}$" in built:
            continue
        why = "태그가 끼어듦" if f"${tex}$" in plain or tex in plain else "글자가 바뀜"
        lno = md.read_text(encoding="utf-8")[: md.read_text(encoding="utf-8").find(f"${tex}$")].count("\n") + 1
        rows.append((md.name, lno, why, tex))

lines = [f"# 인라인 수식 원문·빌드 대조 — 어긋남 {len(rows)}건", ""]
cur = None
for f, ln, why, tex in rows:
    if f != cur:
        lines += ["", f"## {f}"]; cur = f
    lines.append(f"- L{ln} {why}: `{tex[:110]}`")
txt = "\n".join(lines) + "\n"
(OUT.write_text(txt) if OUT else print(txt))
print(f"어긋남 {len(rows)}건, 파일 {len({r[0] for r in rows})}개")

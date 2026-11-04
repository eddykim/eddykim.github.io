#!/usr/bin/env python3
"""웹 이미지 후보 검색 · 라이선스 확인 · 내려받기.

  webimg.py wiki   "<위키백과 문서 제목>" [--lang en] [--thumbs DIR]
  webimg.py search "<영문 검색어>" [--source commons|openverse] [-n 10] [--thumbs DIR]
  webimg.py info   "<Commons 파일명 또는 URL>"
  webimg.py fetch  "<Commons 파일명 또는 URL>" <저장경로> [--width 1400]
  webimg.py fetch  --url <이미지 직접 URL> <저장경로> [--width 1400]

라이선스 등급: A 사용 가능 / B 조건부(사용자 판단) / C 사용 불가 / ? 직접 확인.
등급은 라이선스 이름만 보고 매긴 것이다. 최종 확인은 파일 페이지에서 한다.

네트워크는 전부 curl 로 한다. 이 맥의 Python urllib 은 SSL 인증서 검증에 실패한다.
"""
import argparse
import html
import json
import os
import re
import subprocess
import sys
import urllib.parse

UA = "eddykim-blog-figures/1.0 (https://eddykim.github.io)"
COMMONS_API = "https://commons.wikimedia.org/w/api.php"
OPENVERSE_API = "https://api.openverse.org/v1/images/"
META = "LicenseShortName|LicenseUrl|Artist|Credit|ImageDescription|Restrictions|AttributionRequired"
# 위키백과 문서에 딸린 아이콘·틀 그림. 후보가 아니다.
ICON = re.compile(r"(?i)commons-logo|wik(i|t)(pedia|data|books|quote|source|versity|tionary)|"
                  r"question[_ ]book|ambox|edit-clear|symbol|oojs|padlock|lock-|icon|portal|"
                  r"folder|crystal|nuvola|increase|decrease|steady|red pencil|disambig")
MAX_BYTES = 800_000


# ── 네트워크 ──────────────────────────────────────────────────────────────

def curl(url, out=None):
    cmd = ["curl", "-sSL", "--fail", "-m", "90", "-A", UA]
    if out:
        cmd += ["-o", out]
    r = subprocess.run(cmd + [url], capture_output=True)
    if r.returncode:
        sys.exit(f"curl 실패 (exit {r.returncode}): {url}\n{r.stderr.decode(errors='replace').strip()}")
    return r.stdout


def api(base, **params):
    return json.loads(curl(base + "?" + urllib.parse.urlencode(params)))


# ── 라이선스와 출처 표기 ──────────────────────────────────────────────────

def text(s):
    """extmetadata 의 HTML 조각을 한 줄 텍스트로."""
    s = re.sub(r"<[^>]+>", "", s or "")
    return re.sub(r"\s+", " ", html.unescape(s)).strip()


def https(u):
    return re.sub(r"^(https?:)?//", "https://", u) if u else ""


def is_pd(lic):
    return bool(re.match(r"(?i)public domain|pd\b|pd-", lic or ""))


def tier(lic):
    s = (lic or "").lower()
    if not s:
        return "?"
    if re.search(r"\bnc\b|non-?commercial|fair use|non-free|all rights", s):
        return "C"
    if re.search(r"\bnd\b|no-?deriv", s) or (s.startswith("gfdl") and "cc" not in s):
        return "B"
    if re.search(r"cc0|public domain|^pd|cc[ -]by|^attribution|copyrighted free use|no restrictions", s):
        return "A"
    return "?"


def md(s):
    """캡션(_..._) 안에 넣어도 강조·표·링크 문법이 깨지지 않게 이스케이프."""
    return re.sub(r"([\\`*_\[\]|])", r"\\\1", s)


def tails(who, title, page, lic, lic_url):
    """캡션 끝에 붙일 출처 표기(TASL: 제목·작성자·원본·라이선스) 한·영 두 벌."""
    lic_url = https(lic_url)
    out = []
    # 작성자가 메타데이터에 없으면 파일 페이지의 Author 란을 직접 읽어 채운다.
    # "미상"으로 적어 두면 틀린 표기가 그대로 게시되므로 일부러 눈에 띄는 자리표시를 남긴다.
    for lang, unknown, pd, head in (("KO", "<작성자 확인 필요>", "퍼블릭 도메인", "출처"),
                                    ("EN", "<author: check file page>", "public domain", "Source")):
        name = pd if is_pd(lic) else lic
        lic_md = f"[{name}]({lic_url})" if lic_url else name
        out.append(f"{lang}: ({head}: {md(who) if who else unknown}, [{md(title)}]({https(page)}), {lic_md})")
    return out


def short_title(file_title):
    t = re.sub(r"^(File|파일):", "", file_title)
    return os.path.splitext(t)[0]


def commons_page(file_title):
    return "https://commons.wikimedia.org/wiki/" + urllib.parse.quote(file_title.replace(" ", "_"), safe=":,")


# ── Commons ──────────────────────────────────────────────────────────────

def commons_target(name):
    """파일명·파일 페이지 URL·upload URL·curid URL 을 API 파라미터로."""
    s = urllib.parse.unquote(name.strip())
    if s.startswith("http"):
        u = urllib.parse.urlparse(s)
        q = urllib.parse.parse_qs(u.query)
        if "curid" in q:
            return {"pageids": q["curid"][0]}
        if "/wiki/" in u.path:
            s = u.path.split("/wiki/", 1)[1]
        else:  # upload.wikimedia.org/.../thumb/a/ab/X.jpg/800px-X.jpg 또는 .../a/ab/X.jpg
            parts = u.path.split("/")
            s = parts[-2] if "/thumb/" in u.path else parts[-1]
    s = s.replace("_", " ")
    if not re.match(r"(?i)(file|파일):", s):
        s = "File:" + s
    return {"titles": s}


def commons_info(name, width=None):
    p = dict(action="query", format="json", prop="imageinfo",
             iiprop="url|size|mime|extmetadata", iiextmetadatafilter=META, **commons_target(name))
    if width:
        p["iiurlwidth"] = width
    page = next(iter(api(COMMONS_API, **p)["query"]["pages"].values()))
    if "imageinfo" not in page:
        sys.exit(f"Commons 에 없는 파일: {name}")
    return page["title"], page["imageinfo"][0]


def describe(title, ii, repo="shared"):
    m = {k: v.get("value", "") for k, v in ii.get("extmetadata", {}).items()}
    lic = text(m.get("LicenseShortName"))
    if repo == "local" and tier(lic) == "A":
        grade = "?"   # 위키백과 로컬 파일은 대개 공정이용 그림이다. 파일 페이지를 직접 본다.
    else:
        grade = tier(lic)
    page = commons_page(title) if repo == "shared" else https(ii.get("descriptionurl"))
    who = text(m.get("Artist")) or text(m.get("Credit"))
    return dict(title=title, page=page, lic=lic, grade=grade, lic_url=text(m.get("LicenseUrl")),
                who=who, credit=text(m.get("Credit")), restrictions=text(m.get("Restrictions")),
                desc=text(m.get("ImageDescription")), w=ii.get("width"), h=ii.get("height"),
                mime=ii.get("mime", ""), thumb=ii.get("thumburl"), repo=repo)


def print_row(i, d):
    local = "  [위키백과 로컬 파일]" if d["repo"] == "local" else ""
    print(f"{i:>2}. [{d['grade']}] {d['title']}{local}")
    print(f"    {d['w']}x{d['h']} {d['mime']} | {d['lic'] or '라이선스 표기 없음'} | {d['who'][:70]}")
    if d["restrictions"]:
        print(f"    제한: {d['restrictions']}")
    if d["desc"]:
        print(f"    설명: {d['desc'][:140]}")
    print(f"    {d['page']}")


def save_thumbs(rows, thumbs):
    """후보 섬네일을 받아 둔다. 제목만 보고 고르지 말고 Read 로 직접 본다."""
    os.makedirs(thumbs, exist_ok=True)
    for i, d in rows:
        if not d.get("thumb"):
            continue
        ext = os.path.splitext(urllib.parse.urlparse(d["thumb"]).path)[1].lower() or ".jpg"
        if ext not in (".jpg", ".jpeg", ".png", ".gif", ".webp"):
            ext = ".jpg"
        slug = re.sub(r"[^A-Za-z0-9]+", "-", short_title(d["title"]))[:50].strip("-")
        path = os.path.join(thumbs, f"{i:02d}-{slug}{ext}")
        curl(d["thumb"], path)
    print(f"\n섬네일 → {thumbs}")


def cmd_wiki(a):
    base = f"https://{a.lang}.wikipedia.org/w/api.php"
    data = api(base, action="query", format="json", generator="images", titles=a.article,
               gimlimit=50, redirects=1, prop="imageinfo", iiprop="url|size|mime|extmetadata",
               iiextmetadatafilter=META, iiurlwidth=320)
    pages = data.get("query", {}).get("pages", {})
    rows = []
    for p in pages.values():
        ii = (p.get("imageinfo") or [{}])[0]
        if not ii or ICON.search(p["title"]) or (ii.get("width", 0) < 200 and ii.get("height", 0) < 200):
            continue
        rows.append(describe(p["title"], ii, p.get("imagerepository", "shared")))
    if not rows:
        sys.exit(f"'{a.article}' 문서에서 쓸 만한 그림을 찾지 못했다 (문서 제목·언어 확인).")
    rows.sort(key=lambda d: "ABC?".index(d["grade"]))
    for i, d in enumerate(rows, 1):
        print_row(i, d)
    if a.thumbs:
        save_thumbs(list(enumerate(rows, 1)), a.thumbs)


def search_commons(a):
    data = api(COMMONS_API, action="query", format="json", generator="search", gsrnamespace=6,
               gsrsearch=a.query, gsrlimit=a.n, prop="imageinfo", iiprop="url|size|mime|extmetadata",
               iiextmetadatafilter=META, iiurlwidth=320)
    pages = sorted(data.get("query", {}).get("pages", {}).values(), key=lambda p: p.get("index", 0))
    rows = [(i, describe(p["title"], p["imageinfo"][0])) for i, p in enumerate(pages, 1) if p.get("imageinfo")]
    if not rows:
        sys.exit("검색 결과 없음. 검색어를 바꾸거나 --source openverse 를 써 본다.")
    for i, d in rows:
        print_row(i, d)
    if a.thumbs:
        save_thumbs(rows, a.thumbs)


def ov_license(code, ver):
    code = (code or "").lower()
    if code == "cc0":
        return "CC0 1.0"
    if code == "pdm":
        return "Public domain"
    return f"CC {code.upper()} {ver or ''}".strip()


def search_openverse(a):
    # NC·ND 는 처음부터 뺀다. 결과의 라이선스는 원본 페이지에서 다시 확인해야 한다.
    data = api(OPENVERSE_API, q=a.query, license="by,by-sa,cc0,pdm", page_size=a.n)
    rows = []
    for i, r in enumerate(data.get("results", []), 1):
        lic = ov_license(r.get("license"), r.get("license_version"))
        d = dict(title=r.get("title") or "untitled", page=r.get("foreign_landing_url", ""), lic=lic,
                 grade=tier(lic), lic_url=r.get("license_url", ""), who=r.get("creator") or "",
                 thumb=r.get("thumbnail"))
        rows.append((i, d))
        print(f"{i:>2}. [{d['grade']}] {d['title'][:80]}  ({r.get('source')})")
        print(f"    {r.get('width')}x{r.get('height')} {r.get('filetype') or ''} | {lic} | {d['who'][:70]}")
        print(f"    원본 페이지: {d['page']}")
        print(f"    이미지: {r.get('url')}")
        if r.get("source") == "wikimedia":
            print("    → Commons 파일이다. info 에 원본 페이지 URL 을 넣어 메타데이터를 다시 받는다.")
        else:
            for t in tails(d["who"], d["title"], d["page"], lic, d["lic_url"]):
                print(f"    {t}")
    if not rows:
        sys.exit("검색 결과 없음.")
    if a.thumbs:
        save_thumbs(rows, a.thumbs)


def cmd_search(a):
    (search_openverse if a.source == "openverse" else search_commons)(a)


def cmd_info(a):
    title, ii = commons_info(a.name)
    d = describe(title, ii)
    print(f"파일     {d['title']}")
    print(f"페이지   {d['page']}")
    print(f"크기     {d['w']}x{d['h']} {d['mime']}")
    print(f"라이선스 {d['lic'] or '(표기 없음)'} [{d['grade']}] {https(d['lic_url'])}")
    print(f"작성자   {d['who'] or '(표기 없음)'}")
    if d["credit"] and d["credit"] != d["who"]:
        print(f"크레딧   {d['credit']}")
    if d["restrictions"]:
        print(f"제한     {d['restrictions']}")
    if d["desc"]:
        print(f"설명     {d['desc'][:400]}")
    print()
    for t in tails(d["who"], short_title(title), d["page"], d["lic"], d["lic_url"]):
        print(t)
    if d["grade"] != "A":
        print(f"\n⚠ {d['grade']}등급이다. SKILL.md 의 라이선스 표를 보고 판단한다.")


# ── 내려받기 ─────────────────────────────────────────────────────────────

def shrink(path, width):
    """흰 바탕으로 펴고, 가로 width 를 넘거나 용량이 크면 줄인다. PIL 이 없으면 경고만 한다.

    Commons 의 SVG 렌더링은 배경이 투명하다. 그대로 두면 블로그 다크 모드에서 검은 글씨가
    배경에 묻힌다. 생성 그림(matplotlib)이 전부 흰 바탕이므로 웹 이미지도 흰 바탕으로 맞춘다.
    """
    try:
        from PIL import Image
    except ImportError:
        print("⚠ PIL 이 없어 투명 배경 처리와 크기 조정을 못 했다.")
        return
    im = Image.open(path)
    fmt = im.format
    if fmt == "GIF" and getattr(im, "is_animated", False):
        print("⚠ 움직이는 GIF 는 손대지 않았다.")
        return
    changed = False
    if im.mode in ("RGBA", "LA", "P") and (im.mode != "P" or "transparency" in im.info):
        rgba = im.convert("RGBA")
        bg = Image.new("RGB", rgba.size, "white")
        bg.paste(rgba, mask=rgba.split()[-1])
        im, changed = bg, True
    if im.width > width:
        if im.mode not in ("RGB", "L"):
            im = im.convert("RGB")
        im = im.resize((width, round(im.height * width / im.width)), Image.LANCZOS)
        changed = True
    if changed:
        im.save(path, fmt, **({"quality": 88} if fmt == "JPEG" else {"optimize": True}))
    if fmt == "JPEG" and os.path.getsize(path) > MAX_BYTES:
        im.save(path, "JPEG", quality=80, optimize=True)


def fix_ext(out, src_url):
    """받은 형식에 맞게 확장자를 고친다 (SVG 는 Commons 가 PNG 로 렌더링해 준다)."""
    got = os.path.splitext(urllib.parse.urlparse(src_url).path)[1].lower().replace(".jpeg", ".jpg")
    base, want = os.path.splitext(out)
    if got in (".jpg", ".png", ".gif", ".webp") and want.lower().replace(".jpeg", ".jpg") != got:
        print(f"확장자 {want} → {got} (받은 형식에 맞춤)")
        return base + got
    return out


def cmd_fetch(a):
    if a.url:
        src, out = a.url, a.target
    else:
        _, ii = commons_info(a.target, width=a.width)
        src, out = ii.get("thumburl") or ii["url"], a.out
    if not out:
        sys.exit("저장 경로가 필요하다.")
    out = fix_ext(out, src)
    if not os.path.basename(out).startswith("ext-"):
        print("⚠ 파일명이 ext- 로 시작하지 않는다. 생성 그림과 구분되도록 ext-<내용> 을 권한다.")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    curl(src, out)
    shrink(out, a.width)
    try:
        from PIL import Image
        with Image.open(out) as im:
            w, dims = im.width, f"{im.width}x{im.height}"
    except Exception:
        w, dims = None, "?"
    kb = os.path.getsize(out) // 1000
    print(f"저장 {out}  {dims}  {kb} KB")
    if kb * 1000 > MAX_BYTES:
        print("⚠ 800KB 초과 — --width 를 줄인다.")
    if w and w < 800:
        print(f"⚠ 가로 {w}px — 본문 폭(600~760px)을 채우면 흐리게 보인다. 더 큰 판을 찾는다.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    w = sub.add_parser("wiki", help="위키백과 문서에 실린 그림 목록")
    w.add_argument("article")
    w.add_argument("--lang", default="en")
    w.add_argument("--thumbs")
    w.set_defaults(fn=cmd_wiki)

    s = sub.add_parser("search", help="Commons 또는 Openverse 검색")
    s.add_argument("query")
    s.add_argument("--source", choices=["commons", "openverse"], default="commons")
    s.add_argument("-n", type=int, default=10)
    s.add_argument("--thumbs")
    s.set_defaults(fn=cmd_search)

    i = sub.add_parser("info", help="Commons 파일의 라이선스·작성자·출처 표기")
    i.add_argument("name")
    i.set_defaults(fn=cmd_info)

    f = sub.add_parser("fetch", help="줄여서 내려받기")
    f.add_argument("target", help="Commons 파일명/URL, 또는 --url 일 때 저장 경로")
    f.add_argument("out", nargs="?")
    f.add_argument("--url", help="Commons 가 아닌 곳의 이미지 직접 URL")
    f.add_argument("--width", type=int, default=1400)
    f.set_defaults(fn=cmd_fetch)

    a = ap.parse_args()
    a.fn(a)


if __name__ == "__main__":
    main()

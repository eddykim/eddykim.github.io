#!/usr/bin/env bash
#
# SessionStart hook — Claude Code on the web 세션에서 이 블로그를 빌드/검증할 수 있게
# 환경을 준비한다. 로컬(맥) 개발에는 영향을 주지 않는다.
#
#   1. Ruby gem 설치 (jekyll, html-proofer)
#   2. ruby 실행 파일 디렉터리를 PATH 에 추가 — 컨테이너 기본 PATH 에는 없다
#   3. UTF-8 로케일 고정 — POSIX 로케일이면 html-proofer 가 한글 페이지를 열지 못해
#      링크를 0 개 검사하고도 "성공" 으로 끝난다 (거짓 통과)
#   4. _code/ 그림 생성 스크립트용 파이썬 가상환경(.venv) 구성
#
# 원격 세션에서만 동작한다.

set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

cd "${CLAUDE_PROJECT_DIR:-$(git rev-parse --show-toplevel)}"

env_file="${CLAUDE_ENV_FILE:-/dev/null}"

# ── 1. 로케일 ─────────────────────────────────────────────────────────────
# 본문이 한글이므로 US-ASCII 로는 빌드 산출물을 파싱할 수 없다.
export LANG=C.UTF-8 LC_ALL=C.UTF-8
{
  echo 'export LANG=C.UTF-8'
  echo 'export LC_ALL=C.UTF-8'
} >>"$env_file"

# ── 2. Ruby PATH ─────────────────────────────────────────────────────────
# bundle 은 PATH 에 있지만 gem 이 설치한 실행 파일(jekyll 등) 디렉터리는 빠져 있다.
ruby_bin="$(ruby -e 'puts Gem.bindir' 2>/dev/null || true)"
if [ -d "$ruby_bin" ]; then
  export PATH="$ruby_bin:$PATH"
  echo "export PATH=\"$ruby_bin:\$PATH\"" >>"$env_file"
fi

# ── 3. gem 설치 ──────────────────────────────────────────────────────────
echo "==> bundle install"
bundle install --quiet
bundle check

# ── 4. 파이썬 환경 ────────────────────────────────────────────────────────
# _code/<slug>/generate_figures.py 를 돌려 그림을 다시 뽑을 때 필요하다.
venv=".venv"
if [ ! -x "$venv/bin/python" ]; then
  echo "==> creating $venv"
  python3 -m venv "$venv"
fi
echo "==> installing python deps"
"$venv/bin/python" -m pip install --quiet --upgrade pip
"$venv/bin/python" -m pip install --quiet -r _code/requirements.txt

venv_bin="$PWD/$venv/bin"
export PATH="$venv_bin:$PATH"
echo "export PATH=\"$venv_bin:\$PATH\"" >>"$env_file"

echo "==> ready: $(jekyll --version 2>/dev/null), $("$venv/bin/python" -c 'import numpy,matplotlib;print("numpy",numpy.__version__,"matplotlib",matplotlib.__version__)')"

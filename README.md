# Ed Kim — 광학 계측 연구노트

[eddykim.github.io](https://eddykim.github.io/)에서 운영하는 개인 기술 블로그입니다.
광학·물리학, 최적화 및 프로그래밍, 이미지 프로세싱, 머신러닝과 관련한 글을 다룹니다.

[Chirpy](https://github.com/cotes2020/jekyll-theme-chirpy) Jekyll 테마로 만들었습니다.

## 구조

```shell
.
├── _posts/                포스트 본문 (YYYY-MM-DD-slug.md)
├── _code/<slug>/           각 포스트의 시뮬레이션·분석 코드
├── assets/img/posts/<slug>/  각 포스트의 이미지
└── categories.yaml          분류 체계
```

## 개발

로컬:

```shell
bundle install
bash tools/run.sh            # 로컬 서버
bash tools/test.sh           # 빌드 + html-proofer 링크 검사
```

그림 생성 스크립트는 numpy/matplotlib 를 쓴다 (`_code/requirements.txt`).

```shell
python -m venv .venv && .venv/bin/pip install -r _code/requirements.txt
cd _code/<slug> && python generate_figures.py
```

Claude Code on the web 세션에서는 `.claude/hooks/session-start.sh` 가
gem·파이썬 의존성 설치와 `PATH`·`LANG` 설정을 자동으로 처리한다.
본문이 한글이므로 로케일이 UTF-8 이어야 한다 — POSIX 로케일에서는
html-proofer 가 페이지를 열지 못한 채 링크를 0 개만 검사하고 통과한다.

`generate_figures.py` 들은 한글 라벨 폰트로 `AppleGothic` 을 지정하고 있어
맥이 아닌 환경에서는 한글이 깨진 채로 그려진다.

## License

사이트 인프라(Chirpy 테마 스캐폴드)는 [MIT][mit] 라이선스를 따릅니다.
포스트 본문은 [CC BY 4.0][cc-by]을 따릅니다.

[mit]: https://github.com/cotes2020/chirpy-starter/blob/master/LICENSE
[cc-by]: https://creativecommons.org/licenses/by/4.0/

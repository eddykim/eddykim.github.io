"""13편이 쓰는 공개 데이터 두 벌을 data/ 에 받는다 (약 214 MB). 데이터는 이 폴더에 싣지 않는다.

  ptyLab_zPIE.h5     L. Loetgering, "PtyLab zPIE test data", figshare 10.6084/m9.figshare.21804093.v1, CC BY 4.0
                     USAF 표적, 201 장(364²), λ 708.8 nm, 기록된 z 35 mm. 실험은 zPIE 논문(Optics Letters 45(7), 2030, 2020).
  TwoLayer_bin4.hdf5 M. Du, "ptyLab example data", figshare 10.6084/m9.figshare.14915628.v1, CC BY 4.0
                     두 층 시편(반사 배치), 202 장(364²), λ 708.8 nm, z 38.6 mm. 실험은 POCT 논문(Optics Letters 46(6), 1337, 2021).

받는 곳은 환경 변수 DZ_DATA 로 바꿀 수 있다. 내려받기 주소는 figshare API 에서 그때그때 묻고, MD5 로 확인한다.
"""
import hashlib
import json
import os
import urllib.request

from dz_core import DATA

FILES = {"ptyLab_zPIE.h5": (21804093, "f4604c251ac45103133764d8ac927908"),
         "TwoLayer_bin4.hdf5": (14915628, "a3c36e2899010db92440f64013cdd194")}


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    for name, (article, digest) in FILES.items():
        path = os.path.join(DATA, name)
        if not os.path.exists(path):
            with urllib.request.urlopen(f"https://api.figshare.com/v2/articles/{article}") as r:
                files = json.load(r)["files"]
            url = next(f["download_url"] for f in files if f["name"] == name)
            print("받는 중:", name, url)
            urllib.request.urlretrieve(url, path)
        print(name, "MD5", "일치" if md5(path) == digest else "다름")

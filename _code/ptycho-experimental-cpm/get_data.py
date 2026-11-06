"""12편이 쓰는 공개 데이터를 data/ 에 받는다 (107 MB). 데이터는 이 폴더에 싣지 않는다.

  PtyLab 예제 데이터, M. Du, figshare 10.6084/m9.figshare.14915628, CC BY 4.0
      TwoLayer_bin4.hdf5 — 실공간 타이코그래피 202 장(364², λ 708.8 nm, z 38.6 mm).
      실험은 Du 등, Optics Letters 46(6), 1337 (2021).
받는 곳은 환경 변수 CPM_DATA 로 바꿀 수 있다.
"""
import hashlib
import os
import urllib.request

from cpm_core import DATA, HDF5

URL = "https://api.figshare.com/v2/file/download/28720086"
SHA256 = "0c4234f78dfd77c82e1a0c7e0b861b161918ee10fe76a38a23a0ecd2a212a5ce"

if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    if os.path.exists(HDF5):
        print("있음:", HDF5)
    else:
        print("받는 중:", URL)
        urllib.request.urlretrieve(URL, HDF5)
    h = hashlib.sha256(open(HDF5, "rb").read()).hexdigest()
    print("sha256", "일치" if h == SHA256 else f"다름 ({h})")

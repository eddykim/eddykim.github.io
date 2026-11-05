"""11편이 쓰는 공개 데이터를 data/ 에 받는다 (약 22 MB). 데이터는 이 폴더에 싣지 않는다.

  UConn Smart Imaging Lab (Song & Zheng), https://github.com/SmartImagingLabUConn/Fourier-Ptychography
      USAF_red.mat, blood_aberration.mat, bloodsmear_green.mat — 레포가 GPL-2.0 이라 받아서 쓰기만 한다.
  PtyLab FPM 데이터 (Loetgering 등 2023), figshare 19111472, CC BY 4.0
      USAFTargetFPM.hdf5
"""
import os
import urllib.request

from fpm_core import DATA

UCONN = "https://raw.githubusercontent.com/SmartImagingLabUConn/Fourier-Ptychography/master/Data/"
FILES = {f: UCONN + f for f in ("USAF_red.mat", "blood_aberration.mat", "bloodsmear_green.mat")}
FILES["USAFTargetFPM.hdf5"] = "https://api.figshare.com/v2/file/download/33956984"

if __name__ == "__main__":
    os.makedirs(DATA, exist_ok=True)
    for name, url in FILES.items():
        path = os.path.join(DATA, name)
        if os.path.exists(path):
            print("있음:", name)
            continue
        print("받는 중:", name)
        urllib.request.urlretrieve(url, path)

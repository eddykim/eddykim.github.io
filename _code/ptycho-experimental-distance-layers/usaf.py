"""복원한 USAF 표적에서 세 줄 막대의 주기를 재어 화소 크기(따라서 z)를 거꾸로 구한다.

USAF 1951: 군 g, 요소 e 의 선 쌍 수 = 2^(g + (e−1)/6) lp/mm, 선폭 = 1/(2·선 쌍 수) mm, 주기 = 선폭의 두 배.
막대 하나하나를 어두운 덩어리로 찾고(문턱 + 연결 성분), 같은 방향·같은 크기로 나란한 셋을 한 묶음으로 본다.
묶음마다 막대 길이의 가운데 60 % 를 평균한 단면에 '세 개의 사각형 막대 ⊗ 가우스 흐림' 모형을 맞춰
주기 p [화소] 를 얻고, 가장 가까운 공칭 주기와 짝지어 화소 크기 = 공칭 주기 / p 를 낸다.
"""
import numpy as np
from scipy import ndimage
from scipy.optimize import curve_fit
from scipy.special import erf


def usaf_period(g, e):
    """공칭 주기 [m]."""
    return 1e-3 / 2 ** (g + (e - 1) / 6)


TABLE = [(g, e, usaf_period(g, e)) for g in range(2, 8) for e in range(1, 7)]


def bars_model(x, base, depth, x0, p, w, s):
    """밝기 base 위에 폭 w, 주기 p 인 어두운 막대 셋(가운데 x0), 가장자리 흐림 s."""
    y = np.full_like(x, base, dtype=float)
    for k in (-1, 0, 1):
        c = x0 + k * p
        y -= depth * 0.5 * (erf((x - c + w / 2) / (np.sqrt(2) * s)) - erf((x - c - w / 2) / (np.sqrt(2) * s)))
    return y


def find_triples(A, region, thresh=0.5, min_area=12):
    """어두운 막대 덩어리를 찾아 나란한 셋끼리 묶는다. region = (r0, r1, c0, c1)."""
    r0, r1, c0, c1 = region
    sub = A[r0:r1, c0:c1]
    bright = np.percentile(sub, 90)
    lab, n = ndimage.label(sub < thresh * bright)
    bars = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        area = int((lab[sl] == i).sum())
        if area < min_area or max(h, w) < 2.5 * min(h, w) or max(h, w) > 6.5 * min(h, w):
            continue                                     # 막대 모양(길이 ≈ 5 × 폭)이 아닌 것은 버린다
        vert = h > w
        cy, cx = ndimage.center_of_mass(lab[sl] == i)
        bars.append(dict(vert=vert, L=max(h, w), W=min(h, w), cy=cy + sl[0].start + r0, cx=cx + sl[1].start + c0,
                         sl=(slice(sl[0].start + r0, sl[0].stop + r0), slice(sl[1].start + c0, sl[1].stop + c0))))
    triples, used = [], set()
    for i, b in enumerate(bars):
        if i in used:
            continue
        group = [i]
        for j, c in enumerate(bars):
            if j == i or j in used or c["vert"] != b["vert"] or abs(c["L"] - b["L"]) > 0.25 * b["L"]:
                continue
            along = abs(c["cy"] - b["cy"]) if b["vert"] else abs(c["cx"] - b["cx"])
            across = abs(c["cx"] - b["cx"]) if b["vert"] else abs(c["cy"] - b["cy"])
            if along < 0.3 * b["L"] and across < 1.0 * b["L"]:
                group.append(j)
        if len(group) == 3:
            used.update(group)
            triples.append([bars[k] for k in group])
    return triples


def measure_triple(A, tri, dx_guess):
    """세 막대 묶음의 주기를 맞춤으로 잰다. 반환: dict(주기 px, 짝지은 공칭, 화소 크기)."""
    vert = tri[0]["vert"]
    rs = min(b["sl"][0].start for b in tri), max(b["sl"][0].stop for b in tri)
    cs = min(b["sl"][1].start for b in tri), max(b["sl"][1].stop for b in tri)
    L = tri[0]["L"]
    pad = int(round(L / 5))
    if vert:                                             # 세로 막대 → 열 방향 단면, 행은 막대 길이의 가운데 60 %
        rr = slice(rs[0] + int(0.2 * (rs[1] - rs[0])), rs[1] - int(0.2 * (rs[1] - rs[0])))
        prof = A[rr, cs[0] - pad:cs[1] + pad].mean(0)
        x0 = cs[0] - pad
    else:
        cc = slice(cs[0] + int(0.2 * (cs[1] - cs[0])), cs[1] - int(0.2 * (cs[1] - cs[0])))
        prof = A[rs[0] - pad:rs[1] + pad, cc].mean(1)
        x0 = rs[0] - pad
    x = np.arange(prof.size, dtype=float)
    centers = sorted((b["cx"] if vert else b["cy"]) - x0 for b in tri)
    p0 = (centers[2] - centers[0]) / 2
    guess = [np.percentile(prof, 90), np.ptp(prof), centers[1], p0, p0 / 2, 0.8]
    popt, pcov = curve_fit(bars_model, x, prof, p0=guess, maxfev=20000)
    p_px, p_err = float(popt[3]), float(np.sqrt(pcov[3, 3]))
    g, e, p_nom = min(TABLE, key=lambda t: abs(np.log(t[2] / (p_px * dx_guess))))
    return dict(vert=bool(vert), g=g, e=e, p_px=p_px, p_err=p_err, p_nom=p_nom, dx=p_nom / p_px,
                width_px=float(popt[4]), blur_px=float(popt[5]),
                center=(float(np.mean([b["cy"] for b in tri])), float(np.mean([b["cx"] for b in tri]))))


def measure(O, region, dx_guess, groups=(4, 5)):
    """region 안의 모든 세 막대 묶음을 재고, 지정한 군의 것만 남긴다."""
    A = np.abs(O)
    out = []
    for tri in find_triples(A, region):
        try:
            m = measure_triple(A, tri, dx_guess)
        except (RuntimeError, ValueError):
            continue
        if m["g"] in groups:
            out.append(m)
    return sorted(out, key=lambda m: (m["g"], m["e"], m["vert"]))

"""10편 본문의 계산을 두 경로로 맞대어 확인한다.  python3 verify_post10.py"""
import numpy as np

import ptycho10_core as C

SI, SIO2, SIN = (C.nk(m, C.WL) for m in ("Si", "SiO2", "Si3N4"))
ok = []


def check(name, cond, detail=""):
    ok.append(cond)
    print(("PASS " if cond else "FAIL ") + name + (f"  ({detail})" if detail else ""))


# 1. Parratt 재귀와 특성 행렬(Abeles)이 같은 반사계수를 주는가
worst = 0.0
for layers in ([], [(SIO2, 1.5)], [(SIN, 10.0), (SIO2, 5.0)], [(SIO2, 3.0), (SIN, 50.0), (SIO2, 5.0)]):
    for g in (3, 10, 22, 45, 70):
        a, b = C.parratt(layers, SI, C.WL, g), C.transfer_matrix_refl(layers, SI, C.WL, g)
        worst = max(worst, abs(a - b))
check("Parratt = 특성 행렬", worst < 1e-12, f"최대 차이 {worst:.1e}")

# 2. 층이 없으면 s 편광 Fresnel 식
g = 22
th = np.deg2rad(90 - g)
n2 = SI
ct = np.sqrt(1 - (np.sin(th) / n2) ** 2 + 0j)
fres = (np.cos(th) - n2 * ct) / (np.cos(th) + n2 * ct)
check("단일 계면 = Fresnel r_s", abs(C.parratt([], SI, C.WL, g) - fres) < 1e-12)

# 3. 공통 기준면: 같은 Si 로 받침대를 h 만큼 쌓아 표면만 올리면 위상만 -2 k sin g h 바뀐다(코드 규약)
#    (맨 위에 진공층을 얹는 것은 표면을 올리는 것이 아니므로 위상이 바뀌지 않아야 한다)
h = 0.2
r0 = C.refl([(SIO2, 1.5)], SI, C.WL, 25)
r1 = C.refl([(SIO2, 1.5), (SI, h)], SI, C.WL, 25)
rv = C.refl([(1.0 + 0j, h), (SIO2, 1.5)], SI, C.WL, 25)
check("진공층을 얹어도 반사계수는 그대로", abs(rv - r0) < 1e-12)
expect = -2 * 2 * np.pi / C.WL * np.sin(np.deg2rad(25)) * h
check("높이 2 A = 위상 0.036 rad (25도)", abs(np.angle(r1 / r0) - expect) < 1e-12 and abs(abs(expect) - 0.0363) < 5e-4,
      f"{abs(expect):.4f} rad")

# 4. 비균일 DFT 를 균일 격자에서 구하면 FFT 와 같다
G = C.Geometry(22)
rng = np.random.default_rng(0)
psi = rng.normal(size=(C.N, C.N)) + 1j * rng.normal(size=(C.N, C.N))
a = C.nudft(psi, C.dft_factors(G, G.q_uniform())).reshape(C.N, C.N)
check("비균일 DFT = FFT (균일 격자)", np.abs(a - C.fft2c(psi)).max() < 1e-10)

# 5. 역사상 왕복과 중심의 선형 근사
qe = G.q_exact()
u, v = G.q_to_detector(*qe)
check("q -> 검출기 역사상 왕복", np.nanmax(np.abs(u - G.U)) < 1e-12 and np.nanmax(np.abs(v - G.V)) < 1e-12)
c = C.N // 2
dq = (qe[0][c + 1, c] - qe[0][c - 1, c]) / 2
check("중심 화소 간격 = k sin g p", abs(dq / G.dqx - 1) < 1e-3, f"비 {dq / G.dqx:.5f}")
check("화소 비 dx/dy = 1/sin 22도 = 2.67", abs(G.dx / G.dy - 1 / np.sin(np.deg2rad(22))) < 1e-12, f"{G.dx / G.dy:.3f}")
lo, hi = qe[0].min() / G.k, qe[0].max() / G.k
check("수치구경 0.2, 22도에서 q_x 범위 -0.056k ~ +0.107k", abs(lo + 0.0562) < 5e-4 and abs(hi - 0.1066) < 5e-4,
      f"{lo:.4f} ~ {hi:.4f}")

# 6. 산란 지수: 두께와 무관, 표값이 Eschen 2022 의 92 eV 값(SiO2 2.0, Si3N4 2.9, Al -0.1)과 맞는가
fq = {m: C.delta_beta(m, 13.5)[0] / C.delta_beta(m, 13.5)[1] for m in ("SiO2", "Si3N4", "Al")}
check("delta/beta 표값 = Eschen", abs(fq["SiO2"] - 2.0) < 0.05 and abs(fq["Si3N4"] - 2.9) < 0.05 and abs(fq["Al"] + 0.1) < 0.05,
      ", ".join(f"{k} {v:.2f}" for k, v in fq.items()))
k = 2 * np.pi / 13.5
d, b = C.delta_beta("SiO2", 13.5)
ratios = [(k * d * t) / (k * b * t) for t in (20, 80, 150)]
check("위상/로그 진폭 비는 두께와 무관", np.ptp(ratios) < 1e-12)

# 7. 침투 깊이 (세기 1/e): 29.3 nm 에서 SiO2·Si3N4 약 28 nm
pen = {m: C.WL / (4 * np.pi * C.nk(m, C.WL).imag) for m in ("SiO2", "Si3N4", "Si")}
check("침투 깊이 SiO2 28 nm, Si3N4 28 nm", abs(pen["SiO2"] - 28) < 0.5 and abs(pen["Si3N4"] - 28.4) < 0.5,
      ", ".join(f"{m} {v:.0f} nm" for m, v in pen.items()))

# 8. 본문이 인용하는 계산 결과 (studies10.py 의 캐시)
import os
import studies10 as S
if os.path.exists(S.CACHE):
    R = S.load()
    t = {(r["na"], r["mode"]): r["err"] for r in R["tilt"]}
    check("보정 0.059 / 보정 없음 0.69 (NA 0.2)", round(t[(0.2, "corrected")], 3) == 0.059 and round(t[(0.2, "uncorrected")], 2) == 0.69)
    check("보정 없음 NA 0.05/0.10: 0.14/0.34", round(t[(0.05, "uncorrected")], 2) == 0.14 and round(t[(0.1, "uncorrected")], 2) == 0.34)
    b = {(r["g"], r["os"]): max(abs(r["dphB"]), abs(r["dphC"])) for r in R["bias"]}
    check("보간 편향: 1배 0.04–0.11, 4배 0.002 이하", 0.04 <= min(b[(g, 1)] for g in (22, 30, 40)) and
          max(b[(g, 1)] for g in (22, 30, 40)) <= 0.115 and max(b[(g, 4)] for g in (22, 30, 40)) <= 0.002)
    dev = [(abs(abs(r["rB"]) / abs(S.model_ratio((10, 5, 1), [r["g"]])[0]) - 1),
            abs(np.angle(r["rB"] / S.model_ratio((10, 5, 1), [r["g"]])[0]))) for r in R["pipeline"]]
    check("영역 비: 진폭 최대 0.3%, 위상 최대 0.010 rad", max(d[0] for d in dev) < 0.0035 and max(d[1] for d in dev) < 0.0105)
    sd = {(k, u): R["fit"][k]["crlb"][u]["sd"][0] for k in ("narrow", "wide") for u in ("amp", "phase", "both")}
    want = {("narrow", "amp"): 1.96, ("narrow", "phase"): 11.6, ("narrow", "both"): 0.53,
            ("wide", "amp"): 0.28, ("wide", "phase"): 0.76, ("wide", "both"): 0.13}
    check("하한 표 (질화막 두께 sd)", all(abs(sd[k] - v) / v < 0.02 for k, v in want.items()))
    fq = R["transmission"]["cal"]["fq"]
    check("투과 산란 지수 2.88 / 2.03 / -0.10", round(fq[0], 2) == 2.88 and round(fq[1], 2) == 2.03 and round(fq[2], 2) == -0.10)

print(f"\n{sum(ok)}/{len(ok)} 통과")

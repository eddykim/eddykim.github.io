"""9편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. 전방 모델     PyTorch 의 세기가 6–8편 numpy 전방 모델의 세기와 같다 (정수 위치, 소수 위치).
  2. 자동미분      세 비용의 자동미분 기울기가 6편에서 손으로 유도한 기울기의 정확히 두 배다 (PyTorch 복소 관례).
  3. 전파 연산자   각스펙트럼 전파는 세기를 보존하고(유니터리), dz1 과 dz2 를 이어 전파하면 dz1 + dz2 와 같다.
  4. 다층의 극한   두 장 사이의 간격이 0 이면 다층 모델이 곱 O1 O2 한 장의 모델과 같다.
  5. 다층의 모호성 장마다 상수 c1, c2 를 곱해도 c1 c2 = 1 이면 세기가 같다.
"""
import numpy as np
import torch

import ad_torch as A
from ptycho9_core import (N, pinhole_probe, focused_probe, make_object, scan, object_shape, intensities, noisy,
                          grad_obj, propagate)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


rng = np.random.default_rng(0)
P = pinhole_probe()
pos = scan("fermat")[:40]
shp = object_shape(pos)
obj = make_object(shp)

# 1. 전방 모델
I_np = intensities(obj, P, pos)
I_t = A.intensities(A.t(obj), A.t(P), pos).numpy()
d1 = np.abs(I_np - I_t).max() / I_np.max()
dsub = rng.uniform(-0.5, 0.5, (len(pos), 2))
F = np.fft.fftfreq(N); FY, FX = np.meshgrid(F, F, indexing="ij")
I_np2 = np.array([np.abs(np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(
    np.fft.ifft2(np.fft.fft2(P) * np.exp(-2j * np.pi * (FY * dd[0] + FX * dd[1]))) * obj[y:y + N, x:x + N]), norm="ortho"))) ** 2
    for (y, x), dd in zip(pos, dsub)])
I_t2 = A.intensities(A.t(obj), A.t(P), pos, torch.as_tensor(dsub)).numpy()
d1b = np.abs(I_np2 - I_t2).max() / I_np2.max()
check("PyTorch 전방 모델 = numpy 전방 모델 (정수·소수 위치)", max(d1, d1b) < 1e-12, f"{d1:.1e}, {d1b:.1e}")

# 2. 자동미분 = 손 기울기의 두 배
n, s = noisy(I_np, 1e5, rng)
Ps = P * np.sqrt(s)
O1 = obj * (1 + 0.1 * rng.normal(size=shp))
for kind in ("amp", "poisson", "int"):
    Ot = A.t(O1).requires_grad_(True)
    A.loss(A.intensities(Ot, A.t(Ps), pos), torch.as_tensor(n), kind).backward()
    g_an, _ = grad_obj(O1, Ps, pos, n, kind)
    r = np.abs(Ot.grad.numpy() - 2 * g_an).max() / np.abs(2 * g_an).max()
    check(f"{kind}: 자동미분 기울기 = 2 x 손 기울기", r < 1e-7, f"최대 상대 차이 {r:.1e}")

# 3. 전파 연산자
dx = 0.5e-6
Pf = focused_probe(na=0.15, defocus=15e-6, dx=dx)
H1, H2, H12 = A.propagator(7e-6, dx), A.propagator(11e-6, dx), A.propagator(18e-6, dx)
u = A.t(Pf)
v = A.tifft(A.tfft(u) * H1)
d3a = abs(torch.sum(torch.abs(v) ** 2).item() - torch.sum(torch.abs(u) ** 2).item()) / torch.sum(torch.abs(u) ** 2).item()
w1 = A.tifft(A.tfft(v) * H2)
w2 = A.tifft(A.tfft(u) * H12)
d3b = (torch.abs(w1 - w2).max() / torch.abs(w2).max()).item()
w3 = A.t(propagate(Pf, dx, 18e-6))
d3c = (torch.abs(w3 - w2).max() / torch.abs(w2).max()).item()
check("전파: 세기 보존, 7 + 11 um = 18 um, numpy propagate 와 같음", max(d3a, d3b, d3c) < 1e-12, f"{d3a:.1e}, {d3b:.1e}, {d3c:.1e}")

# 4. 다층의 극한
posf = scan("fermat", step=4)[:30]
shpf = object_shape(posf)
O_a = np.exp(1j * 0.8 * rng.uniform(size=shpf))
O_b = np.exp(1j * 0.8 * rng.uniform(size=shpf))
I_ms = A.multislice([A.t(O_a), A.t(O_b)], A.t(Pf), posf, A.propagator(0.0, dx)).numpy()
I_1 = A.intensities(A.t(O_a * O_b), A.t(Pf), posf).numpy()
d4 = np.abs(I_ms - I_1).max() / I_1.max()
check("간격 0 의 두 장 = 곱 한 장", d4 < 1e-12, f"최대 상대 차이 {d4:.1e}")

# 5. 다층의 모호성
H = A.propagator(20e-6, dx)
c1 = 0.6 * np.exp(0.9j)
Ia = A.multislice([A.t(O_a), A.t(O_b)], A.t(Pf), posf, H).numpy()
Ib = A.multislice([A.t(O_a * c1), A.t(O_b / c1)], A.t(Pf), posf, H).numpy()
d5 = np.abs(Ia - Ib).max() / Ia.max()
check("장마다 상수 c1, 1/c1 을 곱해도 세기가 같다", d5 < 1e-12, f"최대 상대 차이 {d5:.1e}")

print("\n전부 통과" if ok else "\n실패 있음")

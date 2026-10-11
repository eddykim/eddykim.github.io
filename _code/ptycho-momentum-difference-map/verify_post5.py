"""5편 계산의 교차검증. 주장마다 그 주장에서 유도하지 않은 경로로 확인한다.

  1. Maiden 형식   rpie(momentum=maiden()) 가 Maiden 2017 식 19-22 를 그대로 옮긴 참조 구현과 같다.
  2. PtyLab 형식   마찰·피드백을 나눈 모멘텀이 PtyLab mPIE 의 갱신(기울기 = 버퍼 - 현재, 현재 -= feed * m)과 같다.
  3. 피드백 0      모멘텀의 피드백이 0 이면 rPIE 와 같다.
  4. DM = RAAR 꼴  프로브를 고정하면 DM 한 단계 Psi + P_F(2P_O Psi - Psi) - P_O Psi 가 (R_F R_O + I)/2 와 같다.
                   2편 RAAR(beta=1) = (R_s R_m + I)/2 에서 두 집합의 순서를 바꾼 꼴이다.
  5. 중첩 사영     프로브를 고정한 P_O 는 사영이다 (두 번 해도 한 번과 같다).
  6. DM 의 고정점  참 출사파 묶음은 DM 한 단계에 움직이지 않는다.
  7. 세기 제어     모멘텀이 있어도 곱 P O 를 바꾸지 않는다.
  8. 걸음 함수     gamma_obj 를 상수 함수로 주면 같은 상수를 준 것과 같다.
"""
import numpy as np

from ptycho5_core import (N, disc, pinhole_probe, diffuser_probe, make_object, scan, object_shape, intensities,
                          fourier_project, rpie, overlap_projection, dm, maiden, ptylab)

ok = True


def check(name, cond, detail):
    global ok
    ok &= bool(cond)
    print(f"[{'PASS' if cond else 'FAIL'}] {name}: {detail}")


P = diffuser_probe()
P0 = disc(N, 12)
pos = scan("fermat")
shp = object_shape(pos)
obj = make_object(shp)
I = intensities(obj, P, pos)
sI = np.sqrt(I)


def reference(n_iter, alpha, g_o, g_p, T, fric, feed, sign_ptylab, rng):
    """식을 그대로 옮긴 mPIE. sign_ptylab 이면 PtyLab 의 부호 관례(m = 버퍼 - 현재 + fric m, 현재 -= feed m)로 쓴다."""
    rng = np.random.default_rng(rng)
    O = np.ones(shp, complex); Pr = P0.astype(complex).copy()
    vO = np.zeros_like(O); vP = np.zeros_like(Pr); bO = O.copy(); bP = Pr.copy(); c = 0
    for k in range(1, n_iter + 1):
        for j in rng.permutation(len(pos)):
            y, x = pos[j]
            o = O[y:y + N, x:x + N].copy()
            psi = Pr * o
            d = fourier_project(psi, sI[j]) - psi
            O[y:y + N, x:x + N] = o + g_o * np.conj(Pr) * d / ((1 - alpha) * np.abs(Pr) ** 2 + alpha * np.max(np.abs(Pr) ** 2))
            if k >= 2:
                Pr = Pr + g_p * np.conj(o) * d / np.max(np.abs(o) ** 2)
            c += 1
            if c % T == 0:
                if sign_ptylab:
                    vO = (bO - O) + fric * vO; O = O - feed * vO
                else:
                    vO = fric * vO + (O - bO); O = O + feed * vO
                bO = O.copy()
                if k >= 2:
                    if sign_ptylab:
                        vP = (bP - Pr) + fric * vP; Pr = Pr - feed * vP
                    else:
                        vP = fric * vP + (Pr - bP); Pr = Pr + feed * vP
                bP = Pr.copy()
    return O, Pr


def rel(a, b):
    return float(np.abs(a - b).max() / np.abs(b).max())


# 1. Maiden 형식
o1, p1 = rpie(I, pos, P0, shp, 2, alpha=0.05, gamma_obj=0.2, momentum=maiden(), rng=4)
o2, p2 = reference(2, 0.05, 0.2, 1.0, 30, 0.9, 0.9, False, 4)
check("rpie(momentum=maiden()) = 식 19-22 참조 구현 (2회)", max(rel(o1, o2), rel(p1, p2)) < 1e-12,
      f"최대 상대 차이 {max(rel(o1, o2), rel(p1, p2)):.1e}")

# 2. PtyLab 형식
o1, p1 = rpie(I, pos, P0, shp, 2, alpha=0.05, gamma_obj=0.25, gamma_prb=0.5, momentum=ptylab(), rng=5)
o2, p2 = reference(2, 0.05, 0.25, 0.5, 20, 0.7, 0.3, True, 5)
check("마찰·피드백 분리 = PtyLab 부호 관례의 갱신 (2회)", max(rel(o1, o2), rel(p1, p2)) < 1e-12,
      f"최대 상대 차이 {max(rel(o1, o2), rel(p1, p2)):.1e}")

# 3. 피드백 0
mom0 = dict(T=30, fric_obj=0.9, feed_obj=0.0, fric_prb=0.9, feed_prb=0.0)
o1, p1 = rpie(I, pos, P0, shp, 3, gamma_obj=0.25, momentum=mom0, rng=6)
o2, p2 = rpie(I, pos, P0, shp, 3, gamma_obj=0.25, rng=6)
check("피드백 0 인 모멘텀 = rPIE", np.array_equal(o1, o2) and np.array_equal(p1, p2), f"최대 차이 {np.abs(o1 - o2).max():.1e}")

# 4. DM = RAAR(beta = 1), 프로브 고정
rng = np.random.default_rng(0)
psi = np.stack([P * obj[y:y + N, x:x + N] for y, x in pos]) * (1 + 0.3 * rng.normal(size=(len(pos), N, N)))


def P_O(ps, eps=1e-15):
    # 엔진은 분모에 eps * max 를 더해 조명이 거의 없는 화소를 누른다(기본 1e-8). 사영 성질은 eps -> 0 에서 본다.
    o, _ = overlap_projection(ps, np.ones(shp, complex), P, pos, update_probe=False, n_inner=1, eps=eps)
    return np.stack([P * o[y:y + N, x:x + N] for y, x in pos])


def P_F(ps):
    return np.stack([fourier_project(a, s) for a, s in zip(ps, sI)])


po = P_O(psi)
dm_step = psi + P_F(2 * po - psi) - po
R_O = 2 * po - psi
R_F_R_O = 2 * P_F(R_O) - R_O
raar = 0.5 * (R_F_R_O + psi)
check("DM 한 단계 = (R_F R_O + I)/2, 2편 RAAR(beta=1)의 꼴에서 두 집합 순서를 바꾼 것", rel(dm_step, raar) < 1e-12, f"최대 상대 차이 {rel(dm_step, raar):.1e}")

# 5. 중첩 사영은 사영
d5 = rel(P_O(po), po)
d5e = rel(P_O(P_O(psi, 1e-8), 1e-8), P_O(psi, 1e-8))
check("P_O(P_O(Psi)) = P_O(Psi)", d5 < 1e-9, f"최대 상대 차이 {d5:.1e} (엔진 기본 eps 1e-8 에서는 {d5e:.1e})")

# 6. DM 고정점
psi_true = np.stack([P * obj[y:y + N, x:x + N] for y, x in pos])
po_t2 = P_O(psi_true)
step_t = psi_true + P_F(2 * po_t2 - psi_true) - po_t2
check("참 출사파 묶음은 DM 한 단계에 움직이지 않는다 (프로브 고정)", rel(step_t, psi_true) < 1e-9,
      f"최대 상대 변화 {rel(step_t, psi_true):.1e}")

# 7. 세기 제어 + 모멘텀
o1, p1 = rpie(I, pos, P0, shp, 2, gamma_obj=0.2, momentum=maiden(), rng=1)
o2, p2 = rpie(I, pos, P0, shp, 2, gamma_obj=0.2, momentum=maiden(), probe_power="auto", rng=1)
y, x = pos[11]
d7 = rel(p2 * o2[y:y + N, x:x + N], p1 * o1[y:y + N, x:x + N])
check("모멘텀이 있어도 세기 제어는 출사파 P O 를 바꾸지 않는다 (2회)", d7 < 1e-6, f"최대 상대 차이 {d7:.1e}")

# 8. 걸음 함수
o1, _ = rpie(I, pos, P0, shp, 2, gamma_obj=lambda k: 0.25, rng=2)
o2, _ = rpie(I, pos, P0, shp, 2, gamma_obj=0.25, rng=2)
check("gamma_obj = (k -> 0.25) 는 0.25 와 같다", np.array_equal(o1, o2), f"최대 차이 {np.abs(o1 - o2).max():.1e}")

print("\n전부 통과" if ok else "\n실패 있음")

"""그림4(입사각×파장 Psi, Delta 맵)에서 본문 10절이 인용하는 자리를 찾는다.

- Psi 최대점(밝은 점): r_s 가 거의 0 이 되는 자리
- Psi 최소점(어두운 점): r_p 가 거의 0 이 되는 유사-브루스터각 자리
- 그 자리의 박막 왕복 위상 2*beta_phase / 2pi
- 맨 Si 기판(박막 없음)의 유사-브루스터각
- 파장을 고정하고 70~80도를 훑을 때 Delta 변화폭

실행: python map_features.py  (격자와 물성값은 generate_figures.py 의 그림4와 같다)
"""
import numpy as np

from smm_tensor import smm_reflectance_tensor

n_air, n_sio2, n_si = 1.0, 1.46, 3.88 - 0.02j
d_film = 300.0
ang = np.linspace(0.1, 85, 300)
wl = np.linspace(400, 1000, 300)

stack = [n_air, n_sio2, n_si]
r_s = smm_reflectance_tensor(stack, [d_film], np.radians(ang), wl, pol="s")
r_p = smm_reflectance_tensor(stack, [d_film], np.radians(ang), wl, pol="p")
rho = r_p / r_s
psi = np.degrees(np.arctan(np.abs(rho)))
delta = np.degrees(np.angle(rho))


def round_trip(theta_deg, wl_nm):
    """박막 왕복 위상 2*beta_phase 를 2pi 단위로."""
    cos1 = np.sqrt(1 - (n_air * np.sin(np.radians(theta_deg)) / n_sio2) ** 2)
    return 2 * n_sio2 * d_film * cos1 / wl_nm


for label, idx in [("Psi 최대(밝은 점)", psi.argmax()), ("Psi 최소(어두운 점)", psi.argmin())]:
    i, j = np.unravel_index(idx, psi.shape)
    print(f"{label}: Psi = {psi[i, j]:.2f} deg @ {ang[i]:.1f} deg, {wl[j]:.0f} nm | "
          f"|r_s| = {abs(r_s[i, j]):.4f}, |r_p| = {abs(r_p[i, j]):.4f} | "
          f"2*beta/2pi = {round_trip(ang[i], wl[j]):.4f}")

# 맨 Si 기판의 유사-브루스터각 (|r_p| 최소, 1편 프레넬 규약)
th = np.radians(np.linspace(60, 85, 25001))
cos_t = np.sqrt(1 - (np.sin(th) / n_si) ** 2 + 0j)
r_p_bare = (n_si * np.cos(th) - cos_t) / (n_si * np.cos(th) + cos_t)
print(f"맨 Si({n_si}) 유사-브루스터각: {np.degrees(th[np.argmin(abs(r_p_bare))]):.2f} deg")

# 파장 고정, 70~80도 스캔에서 Delta 변화폭 (위상 언랩 후)
sel = (ang >= 70) & (ang <= 80)
swing = np.array([np.ptp(np.degrees(np.unwrap(np.radians(delta[sel, k])))) for k in range(len(wl))])
k = swing.argmax()
print(f"고정 파장 70~80deg 스캔의 Delta 변화폭: 최대 {swing[k]:.1f} deg @ {wl[k]:.0f} nm, "
      f"650nm 부근 {swing[np.searchsorted(wl, 650)]:.1f} deg")

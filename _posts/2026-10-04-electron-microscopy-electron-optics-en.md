---
title: "Electron Microscopy Foundations 2 — Electron Optics: a World Without Concave Lenses"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-electron-optics/
date: 2026-10-04 20:00:00 +0900
page_id: electron-microscopy-electron-optics
categories: [Electron Microscopy, Electron Optics]
tags: [electron-microscopy, electron-optics, magnetic-lens, spherical-aberration, scherzer]
description: Why a magnetic lens can only ever converge, and how Scherzer's theorem fixes the optimum aperture angle and the resolution limit.
math: true
---

[Post 1](/en/posts/electron-microscopy-why-electrons/) ended by naming a suspect. The wavelength of a 200 kV electron is 2.51 pm, yet resolution stalls at 0.2 nm, and the reason is not the wavelength but the aperture angle. Electron microscopes work at 10 mrad, one hundredth of an optical objective, and what prevents that value from growing is spherical aberration.

But glass lenses suffer spherical aberration too. [Geometrical Optics 1](/en/posts/raytracing-spherical-lens-refraction/) traced rays through a spherical surface and found exactly this: the further a ray sits from the axis, the earlier it crosses. Optical designers cancel it with a concave element. Why can an electron lens not do the same? This post answers that. It builds up how a magnetic field acts as a lens, arrives at Scherzer's conclusion that the very structure of such a lens forbids the concave counterpart, and derives the optimum aperture angle and resolution limit from there.

## 1. How a magnetic field becomes a lens

A glass lens bends light at a boundary between two refractive indices. Inside an electron column there is no such boundary. What bends the electron?

The Lorentz force.

$$ \mathbf{F} = -e\,(\mathbf{v} \times \mathbf{B}) $$

The force is perpendicular to both velocity and field, and that matters. If the electron flies straight along the axis and the field points along the axis too, then $\mathbf{v} \times \mathbf{B} = 0$ and nothing happens. A lens needs a radial field component.

In an axially symmetric field, that component appears on its own. Solving $\nabla\cdot\mathbf{B}=0$ near the axis gives

$$ B_r \approx -\frac{r}{2}\frac{dB_z}{dz} $$

The field of a coil strengthens at the lens entrance and weakens at the exit, so $dB_z/dz$ is non-zero and a radial field exists at any distance $r$ off the axis.

The electron now bends in two stages. An electron entering parallel to the axis carries only $v_z$. Meeting $B_r$, it feels a force along $v_z \times B_r$ — azimuthal, around the axis — which gives it a rotational velocity $v_\theta$. That new $v_\theta$ then meets the main field $B_z$, and the force along $v_\theta \times B_z$ points inward, toward the axis. This second force is what collects electrons onto the axis.

Converging by way of two successive bends determines everything that follows. The converging force comes from the product of $B_r$ and $B_z$, and $B_r$ is itself proportional to the gradient of $B_z$. The converging force therefore scales with the *square* of the field. Keep the square in mind.

## 2. The paraxial ray equation, and why it matches light

The electron bends. Can the imaging equations of a glass lens be carried over unchanged?

They can, after one change of coordinates. An electron spirals as it advances, so its path is a helix in the laboratory frame. Viewed from a frame that rotates with it — the Larmor frame — the rotation disappears and pure convergence remains. In that frame a near-axis electron obeys

$$ \frac{d^2 r}{dz^2} + k(z)\,r = 0, \qquad k(z) = \frac{e\,B_z(z)^2}{8\,m_0\,V^{*}} $$

where $V^{*}$ is the accelerating potential with the relativistic correction from post 1, $V(1 + eV/2m_0c^2)$, equal to 239.1 kV at 200 kV. This equation has exactly the form of the paraxial ray equation for light. Focal length, magnification and principal planes therefore carry over intact. When the field region is short enough to treat the lens as thin, the focal length takes a familiar integral form.

$$ \frac{1}{f} = \frac{e}{8\,m_0\,V^{*}} \int B_z(z)^2\,dz $$

Here the paths diverge. $B_z$ enters the integral squared. A square is never negative, so $1/f$ is always positive and $f > 0$. Reversing the coil current and flipping the whole field changes nothing.

A glass lens becomes diverging when the sign of the curvature or the index ratio flips. A magnetic lens has no such freedom. **An axially symmetric magnetic lens can only converge.** The square set aside in section 1 collects its debt here.

## 3. Simulation — an electron in a bell-shaped field

Equations alone are hard to feel, so integrate them. For the field, the analytically convenient Glaser bell model $B_z(z) = B_0 / (1 + (z/a)^2)$ serves, with $B_0 = 0.3$ T and $a = 5$ mm, close to a real objective lens.

```python
# core of generate_figures.py (full code: _code/electron-microscopy-electron-optics/)
def focusing_k(z, v=V_ACC):
    """k(z) in r'' = -k(z) r. B_z is squared, so k > 0 always — hence always converging."""
    return QE * bz(z) ** 2 / (8.0 * M0 * relativistic_potential(v))

def larmor_rate(z, v=V_ACC):
    """Rate at which the image rotates about the axis, dtheta/dz."""
    return np.sqrt(focusing_k(z, v))

def deriv(z, y):                       # y = [r, r', theta]
    return np.array([y[1], -focusing_k(z) * y[0], larmor_rate(z)])
```

<img src="/assets/img/posts/electron-microscopy-electron-optics/en/fig1-magnetic-lens-trajectory.png" alt="Bell-shaped axial field, paraxial electron trajectories converging within it, and the same path seen along the axis" width="850">
_Fig 1. Axial field (left), paraxial rays converging in the Larmor frame (centre), and the same trajectory viewed along the axis (right)_

In the centre panel, three electrons entering parallel to the axis meet at one point. The focal length is 18.0 mm, the same order as a real objective. That they meet at precisely one point regardless of entry height follows from solving the paraxial equation; it is not a statement about real lenses. The paraxial approximation keeps only the first-order term in $r$, and as the height grows the discarded third-order term revives and pulls the focus forward. That discrepancy is spherical aberration, and everything from section 4 onward is the price of the discarded term.

The right-hand panel shows what the Larmor frame cost. Seen along the axis, the electron spirals into the axis rather than running straight, turning through 70° by the time it reaches focus. This is not a calculational device but something that physically happens. Changing magnification in an electron microscope changes the lens current, and a change in current changes the rotation, so the whole image rotates. Nothing like it occurs in an optical microscope, and every electron microscope user meets it whenever they adjust magnification.

## 4. Scherzer's theorem — a list of prohibitions

Back to the opening question. Glass optics cancels spherical aberration by pairing convex with concave. Section 2 showed that a magnetic lens cannot be made diverging, so that route is clearly blocked — but is this merely "no one has found the trick yet", or is it impossible in principle?

Scherzer answered in 1936. Any electron lens satisfying all four of the following conditions has spherical and chromatic aberration coefficients $C_s$ and $C_c$ that are **necessarily positive**.

1. It is rotationally symmetric.
2. Its electric and magnetic fields are static.
3. There is no space charge on the axis.
4. It actually converges.

Positive means uncorrectable. Stacking any number of lenses with positive spherical aberration leaves a positive sum. The cancellation strategy of glass optics is blocked at the level of principle.

The theorem need not read as purely pessimistic. Listing four premises also hands over a list of escape routes: abandon one of them and negative spherical aberration becomes possible. That is precisely what an aberration corrector does, and section 6 returns to it.

## 5. The optimum aperture angle — diffraction against aberration

If $C_s$ cannot be removed, the only remaining choice is where to set the aperture angle. What decides that value?

Three terms pull in opposing directions. The first is diffraction, from post 1, which worsens as the aperture closes.

$$ d_d = 0.61\,\frac{\lambda}{\alpha} $$

The second is spherical aberration. Rays further off-axis miss the focus, so this worsens as the aperture opens — and it worsens as the cube.

$$ d_s = C_s\,\alpha^3 $$

The third is chromatic aberration. Electrons of slightly different energy focus at slightly different distances, and the blur scales with the energy spread.

$$ d_c = C_c\,\alpha\,\frac{\Delta E}{E} $$

Post 1 compared electron guns and noted that a cold field emission gun holds its energy spread to 0.3 eV against 2.3 eV for tungsten. This is where that number enters.

Adding the three in quadrature and setting the derivative with respect to $\alpha$ to zero gives the optimum. Keeping only diffraction and spherical aberration yields a closed form.

$$ \alpha_{opt} = \left(\frac{0.61\,\lambda}{\sqrt{3}\,C_s}\right)^{1/4}, \qquad d_{min} \propto C_s^{1/4}\,\lambda^{3/4} $$

<img src="/assets/img/posts/electron-microscopy-electron-optics/en/fig2-aperture-tradeoff.png" alt="Diffraction, spherical and chromatic blur against aperture semi-angle, compared before and after aberration correction" width="800">
_Fig 2. The three blurs set by aperture angle and their sum. Left, before correction; right, after_

The left panel is the uncorrected case. At $C_s = 0.5$ mm the optimum is 6.5 mrad, the same order as the 10 mrad quoted loosely in post 1. The sum reaches its minimum near where the blue dashed curve (diffraction) crosses the red one (spherical), and the value there is 0.27 nm. Chromatic aberration sits far below the other two in this range and hardly matters.

A word about coefficients is needed here. Whether the spherical blur is written $C_s\alpha^3$, $\frac{1}{2}C_s\alpha^3$ or $\frac{1}{4}C_s\alpha^3$ varies between textbooks, as does whether the terms add in quadrature or directly. A figure such as 0.27 nm should therefore not be taken as absolute. What carries meaning is the proportionality $d_{min} \propto C_s^{1/4}\lambda^{3/4}$. For a conventional form, use the Scherzer point resolution.

$$ d = 0.66\,(C_s\,\lambda^3)^{1/4} $$

Substituting $C_s = 0.5$ mm and $\lambda = 2.51$ pm returns 0.197 nm, matching exactly the 0.2 nm quoted in post 1 for uncorrected HRTEM.

The right panel reduces $C_s$ to 1 µm. The red curve shifts far to the right, the optimum widens almost threefold to 19 mrad, and resolution improves to 0.11 nm. The more consequential change is in the green curve. Chromatic aberration now contributes meaningfully to the sum. Once spherical aberration is removed, chromatic aberration becomes the next wall, which is why gun energy spread and power-supply stability matter increasingly in corrected instruments.

## 6. Why buying down $C_s$ costs so much

Reducing $C_s$ improves resolution. By how much, for how much reduction?

<img src="/assets/img/posts/electron-microscopy-electron-optics/en/fig3-resolution-vs-cs.png" alt="Scherzer point resolution against spherical aberration coefficient at 200 kV and 300 kV" width="750">
_Fig 3. Scherzer point resolution against the spherical aberration coefficient_

The log-log slope is $1/4$. Cutting $C_s$ from 0.5 mm to 1 µm — a factor of 500 — improves resolution from 0.197 nm to 0.042 nm, a factor of only 4.7. A fourth root is that unforgiving. It explains why an aberration corrector, for all its cost and complexity, returns so little, and equally why that little was still judged worth building for.

Wavelength fares better, with an exponent of $3/4$. Raising the accelerating voltage from 200 kV to 300 kV shortens the wavelength from 2.51 pm to 1.97 pm and resolution follows. That road carries its own toll. The more energy an electron carries, the more readily it knocks atoms out of their sites, and specimens made of light elements degrade while being observed. What the beam does to the specimen is the subject of post 3.

So how does a corrector reduce $C_s$? It abandons the first condition on Scherzer's list, rotational symmetry. Stacking multipole elements that are not rotationally symmetric — quadrupoles, hexapoles — over several stages produces negative spherical aberration, which then cancels the positive value of the objective. Scherzer himself pointed out the possibility in 1947, but aligning dozens of elements to micrometre precision and measuring aberrations in real time to feed back took half a century to realize.

## 7. Astigmatism — the aberration users meet daily

Spherical and chromatic aberration follow from principle, and the user cannot touch them. Sitting down at a microscope, however, one meets another aberration constantly.

Astigmatism. It arises because the lens is not perfectly rotationally symmetric. Slight ellipticity from machining the pole piece, magnetic inhomogeneity in the material, local charging from contamination on an aperture — any of these suffice. Focal lengths differ between two perpendicular directions, so focusing one way stretches the image along one axis and focusing the other way stretches it along the other.

Since the cause is broken symmetry, the cure breaks symmetry too. A stigmator, a quadrupole coil, generates opposing astigmatism to cancel it. Unlike the previous two aberrations, astigmatism is not among Scherzer's prohibitions, which is why it yields so easily. It does drift back as contamination accumulates, so it has to be reset repeatedly.

## Summary and what comes next

A magnetic field becomes a lens in two stages. The radial component of an axially symmetric field gives the electron a rotational velocity, and that rotational velocity, meeting the main field, produces a force toward the axis. The converging force scales with the square of the field, which is why a magnetic lens converges regardless of current direction. Moving to the Larmor frame makes the paraxial ray equation identical in form to the one for light, so geometrical optics carries over — at the cost of the image rotating whenever magnification changes.

That no concave counterpart exists is confirmed at the level of principle by Scherzer's theorem: a rotationally symmetric, static, space-charge-free converging lens has strictly positive spherical and chromatic aberration. The aperture angle is therefore pinned near 6.5 mrad, where diffraction and aberration balance, and resolution is trapped at a value proportional to $C_s^{1/4}\lambda^{3/4}$. Evaluated as the Scherzer point resolution with $C_s = 0.5$ mm, that is 0.197 nm — which is where the 0.2 nm of post 1 comes from. A corrector escapes by giving up the rotational symmetry Scherzer assumed, but the fourth-root exponent means a factor of 500 in effort returns a factor of 5 in resolution.

That covers getting electrons onto the specimen. The next post takes up what happens after they arrive. It separates elastic from inelastic collisions, computes the interaction volume the beam spreads into using a Monte Carlo calculation, and establishes the depth from which secondary electrons, backscattered electrons and characteristic X-rays each emerge. Where the contrast in an electron micrograph comes from is decided there.

## References

- P. W. Hawkes and E. Kasper, *Principles of Electron Optics*, 2nd ed., Academic Press, 2018, vol. 1 (paraxial ray equation and magnetic lenses), vol. 2 (aberration theory).
- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 6 (lenses, apertures, resolution).
- O. Scherzer, "Über einige Fehler von Elektronenlinsen," *Zeitschrift für Physik*, vol. 101, pp. 593–603, 1936 (aberration theorem for rotationally symmetric electron lenses).
- O. Scherzer, "Sphärische und chromatische Korrektur von Elektronen-Linsen," *Optik*, vol. 2, pp. 114–132, 1947 (correction using multipole elements).
- M. Haider et al., "Electron microscopy image enhanced," *Nature*, vol. 392, pp. 768–769, 1998 (demonstration of a hexapole corrector).
- [Aberration correction for TEM, *Materials Today*](https://www.sciencedirect.com/science/article/pii/S1369702104005711) (corrector designs and how they came into use).

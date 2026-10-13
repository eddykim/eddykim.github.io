---
title: "Electron Microscopy Foundations 5 — Why Do Defects in a Crystal Show Up as Dark Lines?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-tem-diffraction-contrast/
date: 2026-10-26 20:00:00 +0900
page_id: electron-microscopy-tem-diffraction-contrast
categories: [Electron Microscopy, TEM Imaging]
tags: [electron-microscopy, tem, diffraction, ewald-sphere, bright-field]
description: The electron wavelength is short enough to flatten the Ewald sphere, and TEM diffraction contrast follows from that one fact.
math: true
---

Everything through [post 4](/en/posts/electron-microscopy-sem-imaging/) read signals coming back out of the specimen surface. Secondary electrons came from the top 5 nm, backscattered electrons from the upper third of the interaction volume, and which detector received them separated topographic from compositional contrast.

Now the electrons that pass through the specimen get used. Thin the specimen below 100 nm and a good fraction emerge on the far side without major energy loss. What they carry is not surface information but the crystal structure inside.

What does a crystal do to them? Regularly spaced atomic planes send the electron wave back strongly in particular directions only. This is Bragg diffraction. Which electrons to build an image from, the diffracted ones or the ones that passed straight through, is the subject of this post. The contrast produced by that choice is called diffraction contrast, and it is what makes defects such as dislocations and grain boundaries visible.

The argument starts from one number in post 1. A 200 kV electron has a wavelength of 2.51 pm, sixty times shorter than an X-ray, and that alone changes the character of diffraction completely.

## 1. How does electron diffraction differ from X-ray diffraction?

Electrons obey the same Bragg equation as X-rays. What changes when the wavelength is sixty times shorter? The starting point is the same. For waves reflected from planes spaced $d_{hkl}$ apart to interfere constructively, the path difference must be an integer number of wavelengths.

$$ 2\,d_{hkl}\sin\theta_B = n\lambda $$

Familiar from X-ray diffraction. Put numbers in it. For a typical metal or semiconductor spacing of $d_{hkl} = 0.2$ nm, copper $K\alpha$ X-rays ($\lambda = 0.154$ nm) give a Bragg angle of $\theta_B = 22.65°$. The same planes with 200 kV electrons ($\lambda = 2.51$ pm) give $\theta_B = 0.36°$, smaller by more than sixty.

A Bragg angle below one degree means the diffracted electrons travel almost parallel to the transmitted ones. The diffraction pattern therefore crowds close to the optical axis, and a single objective aperture suffices to separate the two.

More than the angle changes. The rule deciding which reflections switch on changes too, and seeing that requires reciprocal space.

The Ewald construction puts the Bragg condition into geometry. Draw a sphere of radius $1/\lambda$, placed so that the radius from its centre to the origin of the reciprocal lattice is the wavevector of the incident beam. Only reciprocal lattice points lying on that surface satisfy the diffraction condition. The radius varies inversely with wavelength, so a shorter wavelength makes a larger sphere, and a larger sphere grows flatter near the origin. It is the same reason the Earth looks flat to someone standing on a sports field.

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/fig1-ewald-sphere.png" alt="Ewald spheres for X-rays and electrons drawn on the same reciprocal lattice" width="820">
_Fig 1. Ewald spheres for X-rays (Cu Kα, radius 6.5 nm⁻¹) and 200 kV electrons (radius 399 nm⁻¹) drawn on the same reciprocal lattice (spacing 5 nm⁻¹). The X-ray sphere passes through no lattice point except the origin, while the nearly flat electron sphere passes close to a whole row of points. Red dots mark reflections within 0.35 nm⁻¹ of the sphere._

The left panel is the X-ray case. A radius of $1/\lambda = 6.5$ nm$^{-1}$ is only 1.3 times the reciprocal lattice spacing of 5 nm$^{-1}$, so the sphere curves away immediately past the origin. No reciprocal lattice point besides the origin lies on it. That is why X-ray diffraction requires rotating the specimen or scanning the wavelength.

The right panel is the electron case. A radius of 399 nm$^{-1}$ is eighty times the lattice spacing. Near the origin the sphere is effectively a plane: at $g_x = 5$ nm$^{-1}$ it has risen by only 0.031 nm$^{-1}$. Several reciprocal lattice points along a row therefore lie on it, or nearly so, at the same time. Many reflections switch on together without rotating the specimen. How much "nearly" is allowed depends on the specimen thickness, which Section 4 takes up.

This difference sets the character of a transmission electron microscope (TEM) diffraction pattern. Align one crystal direction with the optical axis and the points of the reciprocal lattice plane through the origin perpendicular to it (the zero-order Laue zone) meet the sphere together, within the tolerance of Section 4. A regular grid of spots appears on the screen at once, rather than one spot at a time as in an X-ray measurement.

## 2. Where below the lens does the diffraction pattern form?

Where do the diffracted electrons converge? Consider what the objective lens does.

Electrons leaving different points of the specimen in the same direction converge to one point of the back focal plane (BFP). Conversely, electrons leaving one point in many directions converge to one point of the image plane. Direction maps to position in the back focal plane; origin maps to position in the image plane.

[Ellipsometry Foundations 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/) made the same statement about light: different angles of incidence separate by position in the objective's back focal plane. Light or electrons, a lens converts direction into position, and the back focal plane is where that happens.

So the diffraction pattern forms in the back focal plane and the image forms in the image plane. Figure 2 sets the two side by side as recorded on a real instrument. Focusing the intermediate lens on the back focal plane puts a diffraction pattern on the screen; focusing it on the image plane puts an image there. That is why one control switches between diffraction mode and imaging mode.

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/ext-lens-image-diffraction.png" alt="Lens ray diagram flanked by an image from the image plane and a spot pattern from the back focal plane" width="760">
_Fig 2. The two planes formed by one lens, shown with real TEM recordings. In the central ray diagram, rays leaving in the same direction (same colour) meet at one point of the back focal plane (upper grey line), and rays from one point meet at one point of the image plane (lower grey line). On the left is the image of a polycrystalline specimen in the image plane (scale bar 500 nm), on the right a spot diffraction pattern in the back focal plane (scale bar 5 nm⁻¹); the grain-to-grain differences in brightness on the left are the diffraction contrast of Section 5. (Source: Klingm01, [ElmagLensScheme](https://commons.wikimedia.org/wiki/File:ElmagLensScheme.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

One more thing is needed. The pattern in the back focal plane is made by the entire illuminated region together, so as it stands there is no telling which grain a reflection came from. A second aperture in the image plane solves this. Masking everything but the region of interest lets only electrons from that region continue downward, yielding the diffraction pattern of that region alone. This is selected area electron diffraction (SAED), and it is what allows picking one grain to measure its orientation, or pairing an image with its diffraction pattern.

## 3. Why does moving the aperture invert the contrast?

If the diffraction pattern spreads out in the back focal plane, an aperture placed there can pass only the beam wanted. Making that choice is the job of the objective aperture.

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/fig2-bright-dark-field.png" alt="Ray diagrams with the objective aperture on the transmitted beam and on the diffracted beam" width="800">
_Fig 3. With the objective aperture on the transmitted beam 000 the image is bright field (left); on the diffracted beam g it is dark field (right). Strongly diffracting regions are dark in bright field and bright in dark field. Dark field is drawn here by moving the aperture; in practice the beam is usually tilted instead (centred dark field)._

Place the aperture on the axis, at the transmitted 000 beam, and the diffracted beams are blocked. This is bright field (BF). The more strongly a region diffracts, the fewer electrons remain in the transmitted beam, so that region appears **dark**.

Move the aperture aside onto a particular diffracted beam and the situation reverses. This is dark field (DF). A diffracted beam is named after the reciprocal lattice point that produces it. The vector from the origin to that point is the diffraction vector $g$, and the beam is simply called "the diffracted beam $g$". The axes $g_x, g_y$ of Figure 1 are the components of this vector. In dark field, strongly diffracting regions appear **bright**. Same specimen, same point, and the contrast inverts with nothing but the aperture position.

Recording both is standard practice. A region dark in bright field and bright in the dark field of that $g$ owes its contrast to diffraction into that reflection. One dark in both is either scattering by thickness or atomic number (mass–thickness contrast) or diffraction into some other reflection that was not selected. It is a simple first test for narrowing down the cause.

Dark field carries one caution. Moving the aperture alone means imaging with an off-axis beam, which raises aberration and blurs the image, because off-axis rays suffer more aberration, as [post 2](/en/posts/electron-microscopy-electron-optics/) showed. In practice the aperture stays on the axis and the incident beam itself is tilted so that the wanted diffracted beam travels along the axis. This is called centred dark field. Tilting the beam also changes which reflection sits at the Bragg condition, so $-g$ is set strongly excited before tilting; after the tilt, $g$ is the one at the Bragg condition. On recent instruments with small aberrations, moving the aperture is often good enough.

## 4. How exactly must the Bragg condition be met?

Section 1 said the flattened Ewald sphere switches many reflections on at once. But is "lying exactly on the sphere" really so exact a condition?

It is not. The excitation error $s$ measures how far a reciprocal lattice point departs from the sphere. At $s = 0$ the Bragg condition holds exactly, and the diffracted beam weakens as $s$ grows.

Thinness enters here. In a slab of thickness $t$ each reciprocal lattice point is not a point but a rod extended along the thickness direction over a length of roughly $2/t$. Thinner means longer. The sphere therefore need not pass through the centre of a point; passing anywhere through the rod suffices to produce diffraction. Electron diffraction is forgiving not only because of the wavelength but because the specimen is thin.

The formula simplifies when only the transmitted beam and one diffracted beam are strong. This is the two-beam approximation. Neglecting absorption, the diffracted intensity in this approximation is

$$ I_g = \frac{\sin^2(\pi t\, s_{\text{eff}})}{(\xi_g\, s_{\text{eff}})^2}, \qquad s_{\text{eff}} = \sqrt{s^2 + \frac{1}{\xi_g^2}} $$

where $\xi_g$ is the extinction distance, a length set by the crystal, the reflection and the accelerating voltage. For the silicon 111 reflection it is 60 nm at 100 kV, growing with voltage to about 76 nm at 200 kV; the calculation below uses the round figure of 60 nm.

The two limits show what the formula means. At $s = 0$, $s_{\text{eff}} = 1/\xi_g$ and $I_g = \sin^2(\pi t/\xi_g)$. At a thickness of $\xi_g/2$ the intensity reaches 1: all of the transmitted energy has passed into the diffracted beam. When $s$ is much larger than $1/\xi_g$, $s_{\text{eff}} \approx s$ and $I_g \approx \sin^2(\pi t s)/(\xi_g s)^2$. The period shortens to $1/s$ and the intensity falls as $1/s^2$.

```python
# core of generate_figures.py (full code: _code/electron-microscopy-tem-diffraction-contrast/)
def two_beam_intensity(t, s, xi_g=XI_G):
    """Diffracted intensity in the two-beam approximation.
    At s = 0 this oscillates with a period of exactly xi_g."""
    s_eff = np.sqrt(np.asarray(s, dtype=float) ** 2 + 1.0 / xi_g**2)
    return np.sin(np.pi * np.asarray(t, dtype=float) * s_eff) ** 2 / (xi_g * s_eff) ** 2
```

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/fig3-excitation-error.png" alt="Oscillation of diffracted intensity with thickness and its dependence on excitation error" width="820">
_Fig 4. Diffracted intensity in the two-beam approximation (extinction distance 60 nm). Left, intensity against thickness: at s = 0 (solid) the period is exactly 60 nm, and as s grows (dashed 0.01 nm⁻¹, dotted 0.02 nm⁻¹) period and amplitude both shrink. Right, intensity against s at a thickness of 90 nm: it has already halved at s = 0.009 nm⁻¹._

In the left panel the red curve at $s = 0$ oscillates periodically with thickness, and the period is exactly the extinction distance of 60 nm. Energy passes back and forth between the transmitted and diffracted beams as the electron crosses the specimen, and the length of one round trip is the extinction distance. It resembles two pendulums joined by a string trading their swing, which is why Ewald called this solution the Pendellösung (pendulum solution).

Seen in an image, that oscillation is a thickness fringe. Where the specimen tapers to a wedge, as at its edge, the contrast repeats each time the thickness grows by $\xi_g$, producing stripes. Counting the stripes reads off the thickness in units of $\xi_g$. Because $\xi_g$ depends on the accelerating voltage, the value for the instrument's own voltage must be used: for silicon 111, the 60 nm at 100 kV and the 76 nm at 200 kV differ by a little over 25%.

Departing from $s = 0$ shortens the period and lowers the amplitude, as the blue dashed ($s$ = 0.01 nm$^{-1}$) and green dotted (0.02 nm$^{-1}$) curves show. The right panel fixes the thickness at 90 nm and sweeps $s$ instead. The intensity peaks at $s = 0$ and falls away sharply on both sides, with small subsidiary maxima beyond.

How small a departure is enough to matter? At a thickness of 90 nm the intensity halves at $s$ = 0.009 nm$^{-1}$ and reaches zero at 0.015 nm$^{-1}$. Tilting the crystal by a small angle $\Delta\theta$ gives $s \approx g\,\Delta\theta$, and for planes with $d$ = 0.2 nm, $g = 1/d$ = 5 nm$^{-1}$. A tilt of only 0.1° therefore halves the diffracted intensity, and 0.17° extinguishes it. The Bragg condition need not be met exactly for diffraction to occur, but the price is this much sensitivity, and that sensitivity is what produces the contrast of the next section.

## 5. Why do dislocations and grain boundaries become visible?

In a perfect crystal the whole specimen has the same $s$, so the contrast is uniform. What changes when there is a dislocation?

Around a dislocation the atomic planes are bent. Bent planes mean a locally different $s$. Defects are therefore imaged with the crystal tilted deliberately a little off the Bragg condition, to a small positive $s$ (the reciprocal lattice point just inside the Ewald sphere). On one side of the dislocation the bent planes then rotate back towards the Bragg condition ($s \approx 0$), climbing the steep curve of Figure 4 so that the diffracted beam strengthens there and only there. As the previous section showed, a bend of about 0.1° changes the intensity by half. The transmitted beam loses that much, and the dislocation draws a dark line on a bright background in a bright field image (Figure 5). Because the planes return to the Bragg condition on one side of the dislocation, the line sits to one side of the core rather than directly over it.

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/ext-metadislocation-bright-field.jpg" alt="TEM image in which defects appear as dark curves on a bright background" width="600">
_Fig 5. TEM image of a plastically deformed complex metallic alloy (ε6-Al-Pd-Mn); the field is about 3.5 µm wide, converted from the 200 nm scale bar on the photographic plate. Each dark curve on the bright background is diffraction contrast from the bent planes around a defect, and the path of each defect can be followed wherever it bends or splits. The defects are metadislocations, a kind of dislocation found in complex metallic alloys, and they mostly occur in pairs. (Source: Dr. Michael Feuerbacher, [Metadislocation Network in a Complex Metallic Alloy](https://commons.wikimedia.org/wiki/File:Metadislocation_Network_in_a_Complex_Metallic_Alloy.tif), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), cropped)_

A dislocation core is a few atoms across, yet under strong-beam conditions the line in the image is roughly $\xi_g/3$ wide, some 20 nm here, because what the image records is not the dislocation but the distorted strain field around it. For the same reason, the same dislocation appears or vanishes depending on which $g$ is used. If the displacement is perpendicular to $g$ ($g\cdot R = 0$; for a dislocation usually $g\cdot b = 0$, with $b$ the Burgers vector), the planes behind that reflection merely slide within themselves, neither bending nor shifting along their normal. That reflection is then unaffected and the contrast disappears. Imaging with several $g$ and noting when it vanishes is the standard way of determining a Burgers vector: two vanishing $g$ fix the direction of $b$ as their cross product.

Grain boundaries are simpler. Crystal orientation differs across the boundary, so the two grains start from different $s$. Whichever sits closer to the diffraction condition goes dark in bright field, and grain-to-grain differences in brightness reveal the microstructure, as in the left image of Figure 2. Telling grains apart one by one rests on this.

Diffraction contrast, then, is a map of $s$. How far atomic planes have bent inside the specimen, and where the crystal orientation changes, become differences in diffracted intensity on the screen.

## Summary and what comes next

Most of this post followed from the single fact that electrons have a wavelength sixty times shorter than X-rays. The Bragg angle shrinks to $0.36°$, so the diffracted beams travel nearly parallel to the transmitted one and a single objective aperture separates them. The Ewald sphere radius reaches eighty times the reciprocal lattice spacing, flattening the sphere so that many reflections switch on together.

Because a lens converts direction into position at the back focal plane, that is where the diffraction pattern forms. An aperture placed there selects the transmitted beam for bright field or a diffracted beam for dark field, inverting the contrast. The excitation error, measuring departure from the diffraction condition, sets the diffracted intensity. The dependence is steep: under the conditions of this post, a 0.1° bend of the atomic planes halves the intensity, and that becomes contrast. Dislocations and grain boundaries become visible through it, while the periodic variation with thickness appears as thickness fringes with the extinction distance as their period.

So far the diffracted beams have only been switched on or off. The next post lets them interfere with the transmitted beam instead. Starting from the weak phase object approximation (WPOA), in which the specimen alters only the phase of the electron wave and not its amplitude, it follows how the aberration function built from the spherical aberration and defocus of post 2 becomes the contrast transfer function (CTF). Why a bright spot in a high-resolution image is not guaranteed to be an atom emerges there.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 11–13 (Ewald construction, two-beam equations and extinction distance), ch. 22–26 (diffraction contrast, thickness fringes, dislocation contrast).
- P. B. Hirsch et al., *Electron Microscopy of Thin Crystals*, Butterworths, 1965 (classic treatment of the two-beam approximation and defect contrast).
- J. M. Cowley, *Diffraction Physics*, 3rd ed., North-Holland, 1995, ch. 6–9 (kinematical diffraction by crystals and two-beam dynamical theory).
- L. Reimer and H. Kohl, *Transmission Electron Microscopy: Physics of Image Formation*, 5th ed., Springer, 2008, ch. 7–9 (extinction distance and excitation error).

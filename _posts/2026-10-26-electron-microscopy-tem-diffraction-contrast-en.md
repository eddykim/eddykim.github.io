---
title: "Electron Microscopy Foundations 5 — TEM Diffraction Contrast: the First Way of Seeing a Crystal"
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

What does a crystal do to them? Regularly spaced atomic planes return the electron wave strongly in particular directions. That is Bragg diffraction. Which to build an image from — the diffracted electrons or the ones that passed straight through — is the subject of this post, and that choice is what makes defects such as dislocations and grain boundaries visible.

The argument starts from one number in post 1. A 200 kV electron has a wavelength of 2.51 pm, sixty times shorter than an X-ray, and that alone changes the character of diffraction completely.

## 1. The Bragg condition and the Ewald construction

For waves reflected from planes spaced $d_{hkl}$ apart to interfere constructively, the path difference must be an integer number of wavelengths.

$$ 2\,d_{hkl}\sin\theta_B = n\lambda $$

Familiar from X-ray diffraction. Put numbers in it. For a typical metal or semiconductor spacing of $d_{hkl} = 0.2$ nm, copper $K\alpha$ X-rays ($\lambda = 0.154$ nm) give $\theta_B = 22.65°$. The same planes with 200 kV electrons ($\lambda = 2.51$ pm) give $\theta_B = 0.36°$, smaller by more than sixty.

A Bragg angle below one degree means the diffracted electrons travel almost parallel to the transmitted ones. The diffraction pattern therefore crowds close to the optical axis, and a single objective aperture suffices to separate the two.

More than the angle changes. The rule deciding which reflections switch on changes too, and seeing that requires reciprocal space.

The Ewald construction puts the Bragg condition into geometry. Draw a sphere of radius $1/\lambda$ passing through the origin of the reciprocal lattice; only reciprocal lattice points lying on that surface satisfy the diffraction condition. The radius varies inversely with wavelength, so a shorter wavelength makes a larger sphere, and a larger sphere grows flatter near the origin.

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/fig1-ewald-sphere.png" alt="Ewald spheres for X-rays and electrons drawn on the same reciprocal lattice" width="820">
_Fig 1. Ewald spheres for X-rays and electrons on the same reciprocal lattice_

The left panel is the X-ray case. A radius of $1/\lambda = 6.5$ nm$^{-1}$ is only 1.3 times the reciprocal lattice spacing of 5 nm$^{-1}$, so the sphere curves away immediately past the origin. No reciprocal lattice point besides the origin lies on it. This is why X-ray diffraction requires rotating the specimen or scanning the wavelength.

The right panel is the electron case. A radius of 399 nm$^{-1}$ is eighty times the lattice spacing. Near the origin the sphere is effectively a plane: at $g_x = 5$ nm$^{-1}$ it has risen by only 0.031 nm$^{-1}$. Several reciprocal lattice points along a row therefore lie on it, or nearly so, at the same time. Many reflections switch on together without rotating anything.

This difference sets the character of a TEM diffraction pattern. Align one crystal direction with the optical axis and the whole reciprocal lattice plane perpendicular to it meets the sphere, putting a regular grid of spots on the screen at once, rather than one spot at a time as in an X-ray measurement.

## 2. The diffraction pattern forms at the back focal plane

Where do the diffracted electrons converge? Consider what the objective lens does.

Electrons leaving different points of the specimen in the same direction converge to one point of the back focal plane (BFP). Conversely, electrons leaving one point in many directions converge to one point of the image plane. Direction maps to position in the back focal plane; origin maps to position in the image plane.

[Ellipsometry 1](/en/posts/ellipsometry-electromagnetic-fresnel/) made the same statement about light: different angles of incidence separate by position in the objective's back focal plane. Light or electrons, a lens converts direction into position, and the back focal plane is where that happens.

So the diffraction pattern forms in the back focal plane and the image forms in the image plane. Focusing the intermediate lens on the back focal plane puts a diffraction pattern on the screen; focusing it on the image plane puts an image there. That is the substance of switching between diffraction mode and imaging mode with one control.

One more thing is needed. The pattern in the back focal plane is made by the entire illuminated region together, so as it stands there is no telling which grain a reflection came from. A second aperture in the image plane solves this. Masking everything but the region of interest lets only electrons from that region continue downward, yielding the diffraction pattern of that region alone. This is selected area diffraction (SAED), and it is what allows picking one grain to measure its orientation, or pairing an image with its diffraction pattern.

## 3. Bright field and dark field — choosing with the aperture

If the diffraction pattern spreads out in the back focal plane, an aperture placed there can pass only the beam wanted. That is the objective aperture.

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/fig2-bright-dark-field.png" alt="Ray diagrams with the objective aperture on the transmitted beam and on the diffracted beam" width="800">
_Fig 2. Bright field and dark field according to objective aperture position_

Place the aperture on the axis, at the transmitted 000 beam, and the diffracted beams are blocked. This is bright field (BF). The more strongly a region diffracts, the fewer electrons remain in the transmitted beam, so that region appears **dark**.

Move the aperture aside onto a particular diffracted beam $g$ and the situation reverses. This is dark field (DF), where strongly diffracting regions appear **bright**. Same specimen, same point, and the contrast inverts with nothing but the aperture position.

Recording both is standard practice. A region dark in bright field and bright in dark field owes its contrast to diffraction; one dark in both owes it to absorption or thickness. It is the simplest test for separating the causes of contrast.

Dark field carries one caution. Moving the aperture alone means working with an off-axis beam, which raises aberration and blurs the image. In practice the aperture stays on the axis and the incident beam itself is tilted so that the wanted diffracted beam travels along the axis. Centered dark field, this is called, and the fact from post 2 that off-axis rays suffer more aberration is what drives it.

## 4. The excitation error — tolerance, and its price

Section 1 said the flattened Ewald sphere switches many reflections on at once. But is "lying exactly on the sphere" really so exact a condition?

It is not. The excitation error $s$ measures how far a reciprocal lattice point departs from the sphere. At $s = 0$ the Bragg condition holds exactly, and the diffracted beam weakens as $s$ grows.

Thinness enters here. In a slab of thickness $t$ each reciprocal lattice point is not a point but a rod extended along the thickness direction over a length of roughly $2/t$. Thinner means longer. The sphere therefore need not pass through the centre of a point; passing anywhere through the rod suffices to produce diffraction. Electron diffraction is forgiving not only because of the wavelength but because the specimen is thin.

In the two-beam approximation the diffracted intensity comes out as

$$ I_g = \frac{\sin^2(\pi t\, s_{\text{eff}})}{(\xi_g\, s_{\text{eff}})^2}, \qquad s_{\text{eff}} = \sqrt{s^2 + \frac{1}{\xi_g^2}} $$

where $\xi_g$ is the extinction distance, a length set by the crystal, the reflection and the accelerating voltage. For a typical silicon reflection at 200 kV it is around 60 nm.

```python
# core of generate_figures.py (full code: _code/electron-microscopy-tem-diffraction-contrast/)
def two_beam_intensity(t, s, xi_g=XI_G):
    """Diffracted intensity in the two-beam approximation.
    At s = 0 this oscillates with a period of exactly xi_g."""
    s_eff = np.sqrt(np.asarray(s, dtype=float) ** 2 + 1.0 / xi_g**2)
    return np.sin(np.pi * np.asarray(t, dtype=float) * s_eff) ** 2 / (xi_g * s_eff) ** 2
```

<img src="/assets/img/posts/electron-microscopy-tem-diffraction-contrast/en/fig3-excitation-error.png" alt="Oscillation of diffracted intensity with thickness and its dependence on excitation error" width="820">
_Fig 3. Diffracted intensity computed in the two-beam approximation_

In the left panel the red curve at $s = 0$ oscillates periodically with thickness, and the period is exactly the extinction distance of 60 nm. Energy passes back and forth between the transmitted and diffracted beams as the electron crosses the specimen, and the length of one round trip is the extinction distance.

Seen in an image, that oscillation is a thickness fringe. Where the specimen tapers to a wedge, as at its edge, the contrast repeats each time the thickness grows by $\xi_g$, producing stripes. Counting the stripes reads off the thickness in units of $\xi_g$.

Departing from $s = 0$ shortens the period and lowers the amplitude, as the blue dashed and green dotted curves show. The right panel fixes the thickness at 90 nm and sweeps $s$ instead: the intensity peaks at $s = 0$ and falls away sharply on both sides, with small subsidiary maxima beyond. A slight departure from the diffraction condition costs a great deal of intensity, and that sensitivity is what produces the contrast of the next section.

## 5. Why defects become visible

Assembling the pieces explains how a TEM shows dislocations and grain boundaries.

A dislocation in an otherwise perfect crystal bends the atomic planes around it. Bent planes mean a locally different $s$. If the surrounding crystal is set to $s = 0$, then near the dislocation $s \ne 0$, and the diffracted intensity swings along the steep curve of Figure 3. Hence the dark line a dislocation draws in a bright field image.

A dislocation core is a few atoms across, yet the line in the image runs several nanometres wide, because what is seen is not the dislocation but the distorted strain field around it. For the same reason, the same dislocation appears or vanishes depending on which $g$ is used. If the distortion leaves the spacing along $g$ unchanged, that reflection is unaffected and the contrast disappears. Imaging with several $g$ and noting when it vanishes is the standard way of determining a Burgers vector.

Grain boundaries are simpler. Crystal orientation differs across the boundary, so the two grains start from different $s$. Whichever sits closer to the diffraction condition goes dark in bright field, and grain-to-grain differences in brightness reveal the microstructure. Measuring grain size or examining texture rests on this.

Diffraction contrast, then, is a map of $s$. How far atomic planes have bent inside the specimen, and where the crystal orientation changes, are translated into variations in diffracted intensity and appear on the screen.

## Summary and what comes next

Most of this post followed from the single fact that electrons have a wavelength sixty times shorter than X-rays. The Bragg angle shrinks to $0.36°$, so the diffracted beams travel nearly parallel to the transmitted one and a single objective aperture separates them. The Ewald sphere radius reaches eighty times the reciprocal lattice spacing, flattening the sphere so that many reflections switch on together.

Because a lens converts direction into position at the back focal plane, that is where the diffraction pattern forms. An aperture placed there selects the transmitted beam for bright field or a diffracted beam for dark field, inverting the contrast. The excitation error, measuring departure from the diffraction condition, sets the diffracted intensity, and its steep dependence means a slight bending of atomic planes suffices to generate contrast. Dislocations and grain boundaries become visible through that, while the periodic variation with thickness appears as thickness fringes with the extinction distance as their period.

So far the diffracted beams have only been switched on or off. The next post lets them interfere with the transmitted beam instead. Starting from the weak phase object approximation, in which the specimen alters only the phase of the electron wave and not its amplitude, it follows how the aberration function built from the spherical aberration and defocus of post 2 becomes the contrast transfer function. Why a bright spot in a high-resolution image is not guaranteed to be an atom emerges there.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 11–12 (Ewald construction), ch. 22–26 (diffraction contrast, thickness fringes, dislocation contrast).
- P. B. Hirsch et al., *Electron Microscopy of Thin Crystals*, Butterworths, 1965 (classic treatment of the two-beam approximation and defect contrast).
- J. M. Cowley, *Diffraction Physics*, 3rd ed., North-Holland, 1995, ch. 5–6 (kinematical and dynamical diffraction theory).
- L. Reimer and H. Kohl, *Transmission Electron Microscopy: Physics of Image Formation*, 5th ed., Springer, 2008, ch. 7–9 (extinction distance and excitation error).

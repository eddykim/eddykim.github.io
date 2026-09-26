---
title: "Electron Microscopy Foundations 7 — STEM and Analysis: a Microscope That Counts"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-stem-analytical/
date: 2026-11-05 20:00:00 +0900
page_id: electron-microscopy-stem-analytical
categories: [Electron Microscopy, Analytical EM]
tags: [electron-microscopy, stem, haadf, eels, 4d-stem]
description: Converge the beam to a point again and scan. One detector inner angle sets the atomic number dependence, and the sign-reversal problem of post 6 disappears.
math: true
---

[Post 6](/en/posts/electron-microscopy-tem-phase-contrast-ctf/) left good news and bad news together. The good news was that atomic arrangements appear directly in the image; the bad news was that a bright spot is not guaranteed to be an atom. Contrast inverts with defocus alone, so interpreting such an image requires computation.

This post converges the beam back to a point. As in the SEM of post 4 it scans across the specimen, but it collects what passes **through** rather than what comes back. This is scanning transmission electron microscopy (STEM).

Doing so removes much of post 6's difficulty, because images now form by counting scattered electrons rather than by interference. Counting has no sign. In exchange a new choice appears — which angular range to count — and that one choice decides what becomes visible.

This post closes the series. It starts from the principle behind STEM, covers the analytical techniques and 4D-STEM, and ends with specimen preparation and semiconductor metrology.

## 1. Reciprocity — STEM is TEM turned around

Converging a beam and scanning looks like a wholly different instrument from illuminating broadly. Why then do the two produce similar images?

Because of the principle of reciprocity. Electromagnetic propagation carries time-reversal symmetry, so reversing the direction of every ray leaves the same paths valid. Put a detector where the source was and a source where the detector was, and the result is unchanged.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/en/fig1-reciprocity.png" alt="TEM and STEM ray diagrams compared side by side under the principle of reciprocity" width="800">
_Fig 1. One set of ray paths, read once as TEM and once as STEM_

The lines drawn in the two diagrams are identical. What differs is only which end is called the source and which the detector. STEM is drawn inverted, as is conventional.

The correspondence carries specific content. The illumination angle onto the specimen in TEM corresponds to the collection angle of the detector in STEM, and the collection angle of the TEM objective corresponds to the convergence angle of the STEM probe. An image obtained in STEM with a small bright field detector therefore matches a TEM image.

If they correspond, why use STEM at all? Two reasons. First, several detectors can be placed to divide the signal by angle. Second, since the beam sits at one point, X-rays and energy-lost electrons leaving that point can be collected at the same time to measure composition. Obtaining image and analysis from one scan is the practical advantage.

## 2. HAADF — the detector inner angle sets the atomic number dependence

If the signal divides by angle, which angles show what?

The screened Rutherford cross section of post 3 answers. The differential cross section was

$$ \frac{d\sigma}{d\Omega} \;\propto\; \frac{Z^2}{\left(\sin^2(\theta/2) + \alpha_s\right)^2}, \qquad \alpha_s = \frac{3.4\times10^{-3}\,Z^{0.67}}{E} $$

What an annular detector of inner angle $\theta_1$ and outer angle $\theta_2$ receives follows from integrating this over solid angle, and the integral closes.

$$ \sigma \;\propto\; Z^2\left[\frac{1}{\sin^2(\theta_1/2)+\alpha_s} - \frac{1}{\sin^2(\theta_2/2)+\alpha_s}\right] $$

The leading $Z^2$ suggests a square dependence, but the screening constant $\alpha_s$ itself scales as $Z^{0.67}$. Heavier elements screen over a larger angle, shrinking the bracket and eating into the $Z^2$.

How much it eats depends on the inner angle. An inner angle well above the screening angle makes $\alpha_s$ negligible and the $Z^2$ survives; an inner angle near the screening angle lets $\alpha_s$ bite and the exponent drops.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/en/fig2-z-exponent.png" alt="Effective atomic number exponent against detector inner angle, and cross section by element" width="820">
_Fig 2. The screened Rutherford cross section of post 3 integrated over annular detector ranges_

The left panel is that calculation. The effective exponent is 1.60 at an inner angle of 20 mrad, 1.93 at 80 mrad, and 1.97 at 150 mrad. This is the basis for textbooks quoting HAADF contrast as somewhere between $Z^{1.7}$ and $Z^2$, and it follows from post 3's screening constant with no further assumption.

The dotted lines mark the screening angles: 20 mrad for silicon, 36 mrad for gold. The curve shows that the inner angle must exceed these before the exponent approaches two, which is why high-angle annular dark field (HAADF) detectors typically start above 50 mrad.

The right panel runs the same calculation element by element. The slope on log-log axes is the exponent, and the two inner angles give visibly different slopes.

Here the difference from post 6 emerges. Brightness in a HAADF image reflects the atomic number and the number of atoms beneath the probe almost monotonically. Nothing flips sign the way $\sin\chi$ does, so a bright spot is an atom, and changing focus does not invert the contrast. Hence the name Z-contrast, and hence the habit of recording HAADF first when examining stacks of differing composition.

One qualification belongs here. High-angle scattering is not purely elastic Rutherford scattering; the dominant contribution comes from scattering off atoms displaced by thermal vibration, known as thermal diffuse scattering. That makes the signal far less sensitive to crystal orientation, so diffraction contrast does not intrude and the intensity does not swing with the excitation error of post 5.

## 3. EDS and EELS — reading composition at the same point

Since the beam sits at one point, other signals leaving that point can be collected too. The characteristic X-rays and inelastic scattering of post 3 come into use here.

Energy dispersive X-ray spectroscopy measures the energy of characteristic X-rays with a detector mounted above the specimen. Each element emits at fixed energies, which gives composition. Post 3 noted that X-rays emerge from nearly the whole interaction volume; a TEM foil, however, is under 100 nm thick, so that volume is small to begin with. Spatial resolution that stood at micrometres in the SEM improves to nanometres, and that improvement is the reward for thinning.

Electron energy loss spectroscopy (EELS) disperses the electrons that lost energy crossing the specimen. The energy lost identifies which shell was struck, giving composition, and fine structure just past an absorption edge reveals bonding and coordination as well. It is particularly strong on light elements, making it the tool for mapping oxygen or nitrogen and reading oxidation states.

The two complement rather than compete. EDS favours heavy elements and EELS light ones; EELS offers far better energy resolution but becomes hard to interpret in thick specimens where multiple scattering sets in.

One caution applies to analysis generally: resolution is not set by the probe diameter. The beam broadens as it crosses the specimen, and inelastic scattering is not localized to a single atom in the first place. Reading the composition of one atom therefore demands both a thin specimen and adequate signal, and accumulating signal takes time during which the beam damages the specimen. The trade between accelerating voltage and damage from post 2 returns here.

## 4. 4D-STEM — storing the whole diffraction pattern

Every detector so far has collapsed its signal to one number. A bright field detector sums all electrons inside a central disc; HAADF sums all electrons within a ring. Angular information disappears into those sums.

What if nothing is discarded?

<img src="/assets/img/posts/electron-microscopy-stem-analytical/en/fig3-4dstem.png" alt="The 4D-STEM data structure and the images extracted from it" width="850">
_Fig 3. What comes out of a single four-axis dataset_

Storing the full two-dimensional diffraction pattern $(k_x, k_y)$ at every scan position $(x, y)$ yields a four-axis dataset. This is 4D-STEM. Each pixel holds a pattern rather than an image, so which image to form can be decided after the measurement is over.

Summing the centre gives a bright field image; summing the ring gives an annular dark field image. What conventional detectors did in hardware is now done by computation.

More becomes possible. Measuring how far the centre of mass of the diffraction disc has shifted yields a quantity proportional to the local electric field. Mapping electric or magnetic fields inside a specimen this way is differential phase contrast (DPC).

Ptychography goes furthest. Scanning finely enough that neighbouring diffraction patterns overlap allows the specimen phase to be recovered from that redundancy. Only intensities were measured, yet phase comes back — a direct attack on the problem post 6 described, where phase detaches from intensity.

The recovery is posed as an optimization. The specimen phase and the probe shape become unknowns, iteratively adjusted until the diffraction patterns they predict match those measured. The iterative optimization of the [optimization series](/en/posts/optimization-gradient-descent/) enters here unchanged. The computation is heavy, but resolution beyond the information limit of post 6 repays it. The 0.02 nm figure quoted in post 1 comes from this method.

## 5. Specimen preparation — is that the specimen or the preparation?

Every TEM technique assumes one thing: a sufficiently thin specimen. The diffraction contrast of post 5, the phase object approximation of post 6 and the analysis here all presume something under 100 nm.

Real samples, however, are wafers and devices. A chosen location must be thinned to that dimension.

A focused ion beam (FIB) does this. Accelerated gallium ions mill away the surroundings, leaving a thin slab that is then detached and mounted on a grid. This lift-out procedure produces a 10–20 nm foil from the region of interest. Sectioning one specific transistor became possible because of it.

There is a price. Gallium ions batter the foil surfaces, destroying crystallinity and leaving an amorphous layer, and gallium itself implants into the specimen. Milling at 30 kV leaves damage layers more than 20 nm thick on each face — nearly half the material if the foil is 50 nm. Final polishing at 2–5 kV therefore reduces the damage to a few nanometres.

This is the most common trap in electron microscopy. Whether an amorphous layer in the image was in the specimen or created by the FIB, and whether a defect belongs to the device or to the preparation, has to be settled. Reading an image means first knowing how the specimen was made.

## 6. The role in semiconductor metrology

Let the series close on the metrology side. Where does electron microscopy sit in semiconductor manufacturing?

Linewidth measurement uses CD-SEM, a dedicated instrument that exploits the edge effect in secondary electron images from post 4 to measure pattern width and line width roughness. It is non-destructive and takes whole wafers, so it runs in production lines. As post 4 noted, though, the signal peaks at edges, leaving the placement of the boundary to an algorithm — and differing definitions across instruments give differing numbers.

Cross-sectional structure calls for TEM. Gate stack layer thicknesses, interface roughness and dopant distribution can only be seen by making a foil. Being destructive and slow per image, it serves process development and failure analysis rather than production metrology.

Recently the push has been toward the quantitative. Work has appeared measuring in three dimensions how far the nanosheets of a gate-all-around transistor are strained, using electron ptychography. That is past viewing structure and into reading strain as a number at atomic scale.

Since this series appeared on a blog about optical metrology, a word on the division of labour is fitting. Optical techniques such as ellipsometry are non-destructive, fast and able to cover a whole wafer, but their measurements become thickness or composition only by passing through a model, and they cannot verify that model themselves. Electron microscopy is slow and destructive but shows the structure directly.

The two therefore do not compete. A few specimens are sectioned under the electron microscope to confirm the model, and optical metrology then measures everything else quickly with it. The question that opened post 1, about the gap between wavelength and resolution, ends here.

## Summary — looking back over seven posts

One question ran through the series, posed in post 1: a 200 kV electron has a wavelength of 2.51 pm, so why is the resolution 0.2 nm?

Post 2 answered. Scherzer's theorem denies a rotationally symmetric electron lens any escape from spherical aberration, pinning the aperture angle near 6.5 mrad and trapping resolution at a value proportional to $C_s^{1/4}\lambda^{3/4}$. Aberration correctors circumvented the prohibition, but the fourth-root exponent returned a factor of five for a factor of five hundred in effort.

From post 3 the story moved inside the specimen. That the elastic cross section scales as $Z^2/E^2$ determined most properties of the interaction volume, and that escape depths differ by three orders of magnitude between signals separated the contrast mechanisms of the SEM in post 4. Posts 5 and 6 thinned the specimen to transmit electrons: selecting a diffracted beam gave diffraction contrast, letting it interfere gave phase contrast. And in phase contrast, a bright spot was not guaranteed to be an atom.

This post resolved that. Counting scattered electrons by angle instead of letting them interfere removes the sign. The detector inner angle that sets the angular range decided whether the atomic number dependence lands nearer $Z^{1.6}$ or $Z^{2}$, and that calculation followed directly from post 3's screening constant.

One structure recurred across all seven posts. What becomes visible is decided not by the instrument but by **what gets selected**: which detector receives (post 4), where the aperture sits (post 5), where the focus is set (post 6), which angular range is counted (post 7). Working with an electron microscope is largely a matter of knowing that list of choices, and knowing what each one reveals and what it hides.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 22 (reciprocity), ch. 32–37 (STEM, EDS, EELS).
- S. J. Pennycook and P. D. Nellist (eds.), *Scanning Transmission Electron Microscopy*, Springer, 2011 (HAADF Z-contrast and thermal diffuse scattering).
- R. F. Egerton, *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed., Springer, 2011 (standard treatment of EELS).
- [Four-Dimensional Scanning Transmission Electron Microscopy (4D-STEM), *Microscopy and Microanalysis*](https://academic.oup.com/mam/article/25/3/563/6887544) (overview of 4D-STEM).
- [3D Atomic-Scale Metrology of Strain Relaxation and Roughness in Gate-All-Around Transistors via Electron Ptychography, arXiv](https://arxiv.org/pdf/2507.07265) (strain metrology on GAA transistors).
- [TEM Metrology, Thermo Fisher Scientific](https://www.thermofisher.com/us/en/home/semiconductors/advanced-logic-devices/applications/tem-metrology.html) (cross-sectional metrology in semiconductors).
- [ITRS Metrology 2015](https://www.semiconductors.org/wp-content/uploads/2018/06/2_2015-ITRS-2.0-Metrology.pdf) (CD-SEM and the metrology roadmap).

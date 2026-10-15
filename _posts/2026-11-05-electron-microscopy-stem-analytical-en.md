---
title: "Electron Microscopy Foundations 7 — Can an Electron Microscope Tell Which Atom Is Which?"
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

This post converges the beam back to a point. As in the scanning electron microscope (SEM) of [post 4](/en/posts/electron-microscopy-sem-imaging/), the beam scans across the specimen, but what is collected passes **through** it rather than coming back. This is the scanning transmission electron microscope (STEM).

Doing so removes much of post 6's difficulty, because images now form by counting scattered electrons rather than by interference. Counting has no sign. In exchange a new choice appears: which angular range to count. That one choice decides what becomes visible.

This post closes the series. It starts from the principle behind STEM, covers the analytical techniques and 4D-STEM, and ends with specimen preparation and semiconductor metrology.

## 1. Why does a scanned, focused STEM (Scanning Transmission Electron Microscope) give an image that corresponds to that of a conventional TEM?

A STEM, which converges the beam and scans it, and a transmission electron microscope (TEM), which illuminates broadly, look like wholly different instruments. Why then do they produce similar images?

Because of the principle of reciprocity. For elastic scattering, electron wave propagation carries time-reversal symmetry, so reversing the direction of every ray leaves the same paths valid. Put a detector where the source was and a source where the detector was, and the result is unchanged.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/en/fig1-reciprocity.png" alt="TEM and STEM ray diagrams compared side by side under the principle of reciprocity" width="800">
_Fig 1. Ray diagrams for TEM (left) and STEM (right). The lines are identical; only the arrow and the places of source and detector change. The TEM illumination angle therefore corresponds to the STEM collection angle, and the TEM collection angle to the STEM probe convergence angle._

The lines drawn in the two diagrams of Fig 1 are identical. What differs is only which end is called the source and which the detector. STEM is drawn inverted here so that it matches the TEM layout and the reciprocity is visible. In both diagrams the rays meet at one point in the specimen. Read as STEM, the lower lens focuses a broad beam onto that point as a probe; read as TEM, the objective lens collects the waves scattered from that point.

The correspondence carries specific content. The illumination angle onto the specimen in TEM corresponds to the collection angle of the detector in STEM, and the collection angle of the TEM objective corresponds to the convergence angle of the STEM probe. An image obtained in STEM with a small bright field detector therefore matches a TEM image.

If they correspond, why use STEM at all? Two reasons. First, several detectors can be placed to divide the signal by angle. Second, since the beam sits at one point, X-rays and energy-lost electrons leaving that point can be collected at the same time to measure composition. Obtaining image and analysis from one scan is the practical advantage.

## 2. How closely does HAADF (High-Angle Annular Dark Field) brightness follow atomic number?

If the signal divides by angle, which angles show what?

The conclusion first. Collect only electrons scattered to high angles with a ring-shaped detector, and image brightness follows the atomic number of the atoms beneath the probe almost directly. Unlike post 6, a bright spot is an atom. Why, and the way and extent to which screening lowers the exponent, can be worked out with the screened Rutherford model of post 3. This model's differential cross section is as follows, and the total cross section of post 3 is this integrated over all directions ($E$ in keV):

$$ \frac{d\sigma}{d\Omega} \;\propto\; \frac{Z^2}{\left(\sin^2(\theta/2) + \alpha_s\right)^2}, \qquad \alpha_s = \frac{3.4\times10^{-3}\,Z^{0.67}}{E} $$

Here $\theta$ is the scattering angle, $Z$ the atomic number and $\alpha_s$ the screening constant, which measures how far the orbital electrons hide the nuclear charge. What an annular detector of inner angle $\theta_1$ and outer angle $\theta_2$ receives follows from integrating this over solid angle, and the integral closes.

$$ \sigma \;\propto\; Z^2\left[\frac{1}{\sin^2(\theta_1/2)+\alpha_s} - \frac{1}{\sin^2(\theta_2/2)+\alpha_s}\right] $$

In other words, the detector receives the value for everything inside the inner angle minus the value for everything inside the outer angle. The leading $Z^2$ suggests a square dependence. But the screening constant $\alpha_s$ itself scales as $Z^{0.67}$. Heavier elements screen over a larger angle, shrinking the bracket and eating into the $Z^2$.

How much it eats depends on the inner angle. An inner angle well above the screening angle makes $\alpha_s$ negligible and the $Z^2$ survives; an inner angle near the screening angle lets $\alpha_s$ bite and the exponent drops.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/en/fig2-z-exponent.png" alt="Effective atomic number exponent against detector inner angle, and cross section by element" width="820">
_Fig 2. The screened Rutherford cross section of post 3 integrated over annular detector ranges (200 keV, outer angle 200 mrad). Left: the effective exponent n against inner angle, from a straight-line fit over six elements (C, Si, Cu, Ag, W, Au). The shaded bands are a schematic of the bright field (BF), annular bright field (ABF), annular dark field (ADF) and HAADF detector ranges, and the dotted lines mark the screening angles of Si (20 mrad) and Au (36 mrad). Right: cross section by element at inner angles of 20 mrad (n = 1.60) and 100 mrad (n = 1.95), normalised to carbon = 1. Once the inner angle passes the screening angle, the exponent approaches two._

The left panel of Fig 2 is that calculation. The effective exponent is 1.60 at an inner angle of 20 mrad, 1.93 at 80 mrad, and 1.97 at 150 mrad. Screening is one reason textbooks and measurements put the HAADF exponent below two (often 1.6–1.9), and post 3's screening constant alone reproduces the trend. The model covers only elastic scattering by isolated atoms, though; real exponents shift further with thermal diffuse scattering (TDS), channelling and the pair of elements compared.

The dotted lines mark the screening angle $\theta_0 = 2\sqrt{\alpha_s}$, the angle at which screening starts to soften the scattering. It is the same non-relativistic approximation as the textbook formula for TEM, down to the coefficient. At 200 keV it is 20 mrad for silicon and 36 mrad for gold. The curve shows that the inner angle must exceed these before the exponent approaches two. This is one reason high-angle annular dark field (HAADF) detectors typically start above 50 mrad. The primary purpose of so large an inner angle is to exclude Bragg and other coherent scattering so that the image becomes incoherent.

The right panel runs the same calculation element by element. The slope on log-log axes is the exponent, and the two inner angles give visibly different slopes.

An exponent alone gives little feel for what the image looks like, so convert it to brightness ratios per atom. At an inner angle of 80 mrad a hafnium atom is about 23 times brighter than a silicon atom (26 for pure $Z^2$), gold about 28 times (32 for $Z^2$), germanium about 5.0 times, and oxygen only 0.33 times. Lowering the inner angle to 20 mrad cuts the gold-to-silicon ratio to about 15. In practice, a hafnium oxide gate dielectric on silicon glows, alternating silicon and silicon–germanium layers are told apart, and oxygen nearly vanishes. Seeing light elements takes another signal, which section 3 covers. These are isolated-atom values too, and real ratios change with specimen thickness and orientation.

Here the difference from post 6 emerges. Brightness in a HAADF image reflects the atomic number and the number of atoms beneath the probe almost monotonically. Nothing flips sign the way $\sin\chi$ does, so a bright spot is an atom, and changing focus does not invert the contrast. Hence the name Z-contrast, and hence the habit of recording HAADF first when examining stacks of differing composition. In thick specimens, however, brightness is no longer proportional to thickness, and channelling, in which electrons concentrate along atomic columns, also plays a part. Converting brightness into a number of atoms therefore requires the thickness as well. The annular bright field (ABF) detector, as the bands of Fig 2 show, collects a low-angle ring just outside the bright field. Fig 3 shows the same specimen imaged in annular dark field and annular bright field. In dark field the heavy columns are bright and oxygen is invisible; in annular bright field even oxygen appears. The angles a detector collects decide which elements show.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/ext-srtio3-adf-abf.jpg" alt="Atomic-resolution annular dark field and annular bright field images of SrTiO3 compared" width="600">
_Fig 3. Atomic-resolution STEM images of the perovskite oxide SrTiO₃ viewed along [110]. In the annular dark field image (left) the heavy strontium columns are brightest and oxygen is invisible; in the annular bright field image (right) even the oxygen columns appear as dark dots. The overlaid circles mark strontium (green), titanium (grey) and oxygen (red); the source does not give the detector angles. (Source: Magnunor, [Scanning transmission electron microscopy srtio3 compare adf abf](https://commons.wikimedia.org/wiki/File:Scanning_transmission_electron_microscopy_srtio3_compare_adf_abf.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

Most high-angle scattering is not pure elastic Rutherford scattering. Its main part is scattering off atoms displaced by thermal vibration, which is thermal diffuse scattering. That makes the signal less sensitive to crystal orientation, so diffraction contrast does not intrude, and the signal is largely free of the problem that a slight change in orientation changes the diffracted intensity sharply. Diffraction contrast in post 5 was the opposite case. There the excitation error $s$ set the diffracted intensity, and at a thickness of 90 nm a departure of only 0.009 nm⁻¹ in $s$ halved it. The HAADF signal is still sensitive to channelling and lattice strain, of course.

## 3. What do EDS (Energy Dispersive X-ray Spectroscopy) and EELS (Electron Energy Loss Spectroscopy) read?

While the beam stays at one point, can we also learn what that point is made of? The characteristic X-rays and inelastic scattering of post 3 come into use here.

Energy dispersive X-ray spectroscopy (EDS) measures the energy of characteristic X-rays with a detector mounted above the specimen. Each element emits at fixed energies, which gives composition. Post 3 noted that X-rays emerge from nearly the whole interaction volume. A TEM foil, however, is under 100 nm thick, so that volume is small to begin with. Spatial resolution that stood at micrometres in the SEM improves to nanometres, and that improvement is the reward for thinning. Fig 4 is an example not of resolution but of obtaining a HAADF image and EDS together at the same place.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/ext-stem-haadf-eds-cacarbonate.jpg" alt="HAADF images, an EDS element map and EDS spectra of inclusions in cyanobacteria" width="760">
_Fig 4. STEM analysis of mineral inclusions inside cyanobacteria cells. (a) and (b) are HAADF images, (c) an EDS element map of the same area (blue carbon, red calcium, green phosphorus), and (d) EDS spectra taken at the red and green circles in (b). Inclusions holding calcium and phosphorus are much brighter in HAADF than the carbon-rich cells, and the Ca Kα (3.7 keV) and P Kα (2.0 keV) peaks separate the two inclusions into a calcium one and a phosphorus one, which the source paper interpreted as calcium carbonate and polyphosphate respectively. (Source: Jinhua Li et al., [Minerals-06-00010-g002](https://commons.wikimedia.org/wiki/File:Minerals-06-00010-g002.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), cropped)_

Electron energy loss spectroscopy (EELS) disperses the electrons that lost energy crossing the specimen. The energy lost identifies which shell was struck, giving composition, and fine structure just past an absorption edge reveals bonding and coordination as well. It is particularly strong on light elements, making it the tool for mapping oxygen or nitrogen and reading oxidation states.

The two complement rather than compete. The table collects the differences.

| Item | EDS | EELS | Reason |
|---|---|---|---|
| Measures | Energy of characteristic X-rays leaving the specimen | Energy lost by electrons crossing the specimen | The same ionization event, seen from the X-ray side or the electron side |
| Strong on | Heavy elements | Light elements (oxygen, nitrogen and so on) | Excited light atoms mostly emit Auger electrons rather than X-rays, so X-ray yield is low. EELS measures the energy loss directly, whether or not an X-ray follows |
| Energy resolution | Broad; neighbouring peaks overlap easily | Around 0.3–1 eV | In EDS, statistical fluctuation in the charge created inside the detector broadens each peak; EELS disperses the electrons directly with a magnetic spectrometer |
| Information | Elemental composition | Composition plus bonding, coordination and oxidation state | Fine structure near the edge reflects the empty electronic states |
| Specimen thickness | Thinner gives better spatial resolution | Must be thin to interpret | In thick specimens, multiple scattering makes the losses overlap in EELS |

The boundary is not sharp. With a windowless detector, EDS can quantify boron through fluorine. Still, the usual first choice is EELS for oxygen distribution or oxidation state, and EDS for a quick survey of heavier elements.

One caution applies to analysis generally: resolution is not set by the probe diameter. The beam spreads as it crosses the specimen (beam broadening), and inelastic scattering with small energy loss takes place over a region wider than one atom. Reading the composition of one atom therefore demands both a thin specimen and adequate signal. Accumulating signal takes time, and during that time the beam damages the specimen. The trade between accelerating voltage and damage from [post 2](/en/posts/electron-microscopy-electron-optics/) returns here.

## 4. What more do we gain by storing the whole diffraction pattern?

Every detector so far has collapsed its signal to one number. A bright field detector sums all electrons inside a central disc; HAADF sums all electrons within a ring. Angular information disappears into those sums.

What if nothing is discarded?

<img src="/assets/img/posts/electron-microscopy-stem-analytical/en/fig3-4dstem.png" alt="The 4D-STEM data structure and the images extracted from it" width="850">
_Fig 5. The structure of 4D-STEM data. One diffraction pattern $$(k_x, k_y)$$ is stored at every scan position $(x, y)$. After the measurement, summing the centre gives a bright field image and summing the ring an annular dark field image, while the centre-of-mass shift gives differential phase contrast and the full patterns give ptychography._

Storing the full two-dimensional diffraction pattern $(k_x, k_y)$ at every scan position $(x, y)$ yields a four-axis dataset. This is four-dimensional STEM (4D-STEM). Each pixel holds a pattern rather than an image, so which image to form can be decided after the measurement is over. For scale: scanning 256×256 positions with a 128×128-pixel detector gives about a billion values, or 4 GiB at 32 bits each. That is why a four-axis dataset is a measurement problem and a computing problem at once.

As Fig 5 shows, summing the centre gives a bright field image and summing the ring gives an annular dark field image. What conventional detectors did in hardware is now done by computation.

More becomes possible. Measuring how far the centre of mass of the diffraction disc has shifted gives a value proportional to the electric field the beam met while crossing the specimen, summed through the thickness (the projected field). Mapping electric or magnetic fields inside a specimen this way is differential phase contrast (DPC).

Ptychography goes furthest. Scanning finely enough that the areas illuminated at neighbouring probe positions overlap allows the specimen phase to be recovered from that redundancy. Only intensities were measured, yet phase comes back. This attacks the problem post 6 described head-on, where phase detaches from intensity.

The recovery is posed as an optimization. The specimen phase and the probe shape become unknowns, iteratively adjusted until the diffraction patterns they predict match those measured. The iterative optimization of the [optimization series](/en/posts/optimization-gradient-descent/) enters here unchanged. The computation is heavy, but resolution rises beyond the information limit of post 6. Resolution beyond the lens information limit was first shown by Nellist, McCallum and Rodenburg in 1995, and the 0.02 nm figure quoted in post 1 also comes from this method. The reconstruction algorithms themselves are the subject of this blog's Computational Imaging and Ptychography series.

## 5. What does milling a foil with a FIB (Focused Ion Beam) leave behind?

Every TEM technique assumes one thing: a sufficiently thin specimen. The diffraction contrast of [post 5](/en/posts/electron-microscopy-tem-diffraction-contrast/), the phase object approximation of post 6 and the analysis here all presume something under 100 nm.

Real specimens, however, are wafers and devices. A chosen location must be thinned to that dimension.

A focused ion beam (FIB) does this. Accelerated gallium ions mill away the surroundings, leaving a thin slab that is then detached and mounted on a grid. This lift-out procedure produces a foil under 100 nm thick from the region of interest, or a few tens of nanometres for high-resolution work. Picking out one transistor and sectioning it became possible because of it. Fig 6 shows a foil made this way.

<img src="/assets/img/posts/electron-microscopy-stem-analytical/ext-fib-lamella.jpg" alt="A FIB foil mounted on a copper grid and an atomic-resolution STEM image of that foil" width="760">
_Fig 6. A TEM foil made by FIB, seen at increasing magnification. The left two images are secondary electron images from the FIB, showing the foil on its copper grid and the thin film, protective layer and substrate near its top; the protective layer is deposited before milling to shield the surface from the ion beam. The right image is an atomic-resolution STEM image of the film–substrate interface in the same foil. (Source: Magnunor, [Fib sample preparation tem](https://commons.wikimedia.org/wiki/File:Fib_sample_preparation_tem.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

There is a price. Gallium ions batter the foil surfaces, destroying crystallinity and leaving an amorphous layer, and gallium itself implants into the specimen. Milling at 30 kV alone leaves an amorphous layer of some 20 nm on each face. If the foil is 50 nm, that is most of the material, counting both faces. Final polishing at 2–5 kV therefore reduces the damage to a few nanometres.

This is one of the common traps in electron microscopy. Whether an amorphous layer in the image was in the specimen or created by the FIB, and whether a defect belongs to the device or to the preparation, has to be settled. Reading an image means first knowing how the specimen was made.

## 6. What does electron microscopy measure in semiconductor manufacturing?

Let the series close on the metrology side. Where does electron microscopy sit in semiconductor manufacturing?

Linewidth measurement uses the critical-dimension SEM (CD-SEM), a dedicated instrument that exploits the edge effect in secondary electron images from post 4 to measure pattern width and line width roughness. It is non-destructive and takes whole wafers, so it runs in production lines. As post 4 noted, though, the signal peaks at edges, leaving the placement of the boundary to an algorithm, and differing definitions across instruments give differing numbers.

Cross-sectional structure calls for TEM. Gate stack layer thicknesses and interface roughness can only be seen by making a foil (dopant profiles are left to secondary ion mass spectrometry (SIMS) and scanning-probe methods). Being destructive and slow per image, it serves mainly process development, failure analysis and the reference metrology that calibrates inline tools.

Recently the push has been toward the quantitative. Work has appeared that uses electron ptychography to measure, in three dimensions, how the spacing between atomic layers relaxes near the interfaces of 5 nm silicon channels in gate-all-around (GAA) transistors, and how rough those interfaces are. That is past viewing structure and into reading strain as a number at atomic scale.

Since this series appeared on a blog about optical metrology, it closes on the division of labour. Optical techniques such as ellipsometry are non-destructive, fast and able to cover a whole wafer. But their measurements become thickness or composition only by passing through a model ([Ellipsometry Foundations 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/)), and they cannot verify that model themselves. Electron microscopy is slow and destructive but shows the structure directly.

The two therefore do not compete. A few specimens are sectioned under the electron microscope to confirm the model, and optical metrology then measures everything else quickly with it. The question that opened post 1, about the gap between wavelength and resolution, ends here.

## Summary

One question ran through the series, posed in [post 1](/en/posts/electron-microscopy-why-electrons/): a 200 kV electron has a wavelength of 2.51 pm, so why is the resolution 0.2 nm?

[Post 2](/en/posts/electron-microscopy-electron-optics/) answered. Scherzer's theorem denies a rotationally symmetric electron lens any escape from spherical aberration, pinning the aperture angle near 6.5 mrad and trapping resolution at a value proportional to $C_s^{1/4}\lambda^{3/4}$. Aberration correctors circumvented the prohibition, but because of the fourth-root exponent, cutting $C_s$ by a factor of 500 improved resolution only 4.7 times.

From [post 3](/en/posts/electron-microscopy-beam-specimen-interaction/) the story moved inside the specimen. The cross section for large-angle scattering scales as $Z^2/E^2$, but because of screening the total cross section scales roughly as $Z^{4/3}/E$. This determined most properties of the interaction volume, and the three-orders-of-magnitude spread in escape depth between signals separated the contrast mechanisms of the SEM in [post 4](/en/posts/electron-microscopy-sem-imaging/). [Post 5](/en/posts/electron-microscopy-tem-diffraction-contrast/) and [post 6](/en/posts/electron-microscopy-tem-phase-contrast-ctf/) thinned the specimen to transmit electrons: selecting a diffracted beam gave diffraction contrast, letting it interfere gave phase contrast. And in phase contrast, a bright spot was not guaranteed to be an atom.

This post resolved that. Counting scattered electrons by angle instead of letting them interfere removes the sign. The detector inner angle that sets the angular range decided whether the atomic number dependence lands nearer $Z^{1.6}$ or $Z^{2}$, and that calculation followed directly from post 3's screening constant. Heavy elements are read with HAADF, light elements and chemical state with EELS, and composition maps with EDS. So can an electron microscope tell which atom is which? Yes, depending on what it counts.

One structure recurred across all seven posts. What becomes visible is decided not by the instrument but by **what gets selected**: which detector receives (post 4), where the aperture sits (post 5), where the focus is set (post 6), which angular range is counted (post 7). Working with an electron microscope is largely a matter of knowing that list of choices, and knowing what each one reveals and what it hides.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, §6.2.A (reciprocity), §9.4 (STEM detectors), ch. 22 (amplitude contrast), ch. 32–36 (XEDS), ch. 37–40 (EELS).
- S. J. Pennycook and P. D. Nellist (eds.), *Scanning Transmission Electron Microscopy: Imaging and Analysis*, Springer, 2011 (HAADF Z-contrast and thermal diffuse scattering).
- R. F. Egerton, *Electron Energy-Loss Spectroscopy in the Electron Microscope*, 3rd ed., Springer, 2011 (standard treatment of EELS).
- P. D. Nellist, B. C. McCallum, J. M. Rodenburg, "Resolution beyond the 'information limit' in transmission electron microscopy," *Nature* 374, 630–632 (1995) (the first ptychography beyond the lens information limit).
- [C. Ophus, "Four-Dimensional Scanning Transmission Electron Microscopy (4D-STEM): From Scanning Nanodiffraction to Ptychography and Beyond," *Microscopy and Microanalysis* 25(3), 563–582 (2019)](https://academic.oup.com/mam/article/25/3/563/6887544) (overview of 4D-STEM).
- [3D Atomic-Scale Metrology of Strain Relaxation and Roughness in Gate-All-Around Transistors via Electron Ptychography, *Nature Communications* 17, 3561 (2026)](https://arxiv.org/abs/2507.07265) (layer-spacing relaxation and interface roughness in GAA channels, arXiv:2507.07265).
- [J. Li et al., "Biomineralization Patterns of Intracellular Carbonatogenesis in Cyanobacteria: Molecular Hypotheses," *Minerals* 6(1), 10 (2016)](https://www.mdpi.com/2075-163X/6/1/10) (source paper of Fig 4, CC BY 4.0).
- [TEM Metrology, Thermo Fisher Scientific](https://www.thermofisher.com/us/en/home/semiconductors/advanced-logic-devices/applications/tem-metrology.html) (cross-sectional metrology in semiconductors).
- [ITRS Metrology 2015](https://www.semiconductors.org/wp-content/uploads/2018/06/2_2015-ITRS-2.0-Metrology.pdf) (CD-SEM and the metrology roadmap).

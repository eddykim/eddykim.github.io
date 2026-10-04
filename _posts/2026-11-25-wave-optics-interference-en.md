---
layout: post
date: 2026-11-25 20:00:00 +0900
title: "Wave Optics Foundations 1 — Why Can Light Plus Light Be Darker? Interference"
lang: en
lang-exclusive: ["en"]
permalink: /posts/wave-optics-interference/
page_id: wave-optics-interference
categories: [Optics, Wave Optics]
tags: [wave-optics, interference, interferometer, thin-film, fabry-perot]
description: "Light adds as amplitudes and is measured as intensity. From Young's experiment to thin films, the Michelson and the Fabry–Pérot, this post follows how interference works and what it is used for."
math: true
---

When a drop of oil lands on wet asphalt after rain, rainbow-coloured bands appear (Figure 1). The oil is transparent and the asphalt is black; nothing here has a colour of its own, yet colours appear. Stranger still, the film looks black at its edge, where it becomes very thin. Light reflected from the top and from the bottom of the film both reach the eye, and yet two beams added together are darker than one.

<img src="/assets/img/posts/wave-optics-interference/ext-oil-film.jpg" alt="Rainbow interference bands from an oil film on a wet road" width="600">
_Figure 1. A film of diesel spread on a wet road. The film thickness varies from place to place, and depending on the thickness different colours are strengthened or weakened. (Source: John (ex-user Guinnog), [Dieselrainbow](https://commons.wikimedia.org/wiki/File:Dieselrainbow.jpg), [CC BY-SA 2.5](https://creativecommons.org/licenses/by-sa/2.5))_

In the [geometrical optics foundations](/en/posts/geometric-optics-fermat-eikonal/), rays were lines carrying energy, and where two bundles of rays overlapped their brightnesses added. That picture cannot explain Figure 1. We need the fact that light is a wave, that a wave has crests and troughs, and that a crest meeting a trough cancels it. This is interference. At the end of [Geometrical Optics Foundations 4](/en/posts/geometric-optics-aberrations/), once aberrations shrank to the size of a wavelength, the ray picture gave way to diffraction. This series starts from there and takes in turn interference, coherence, diffraction, Fourier imaging, Gaussian beams and high-NA focusing. The first part is interference, the basis of all wave optics.

## 1. Add amplitudes, measure intensity

There is one key to interference: light adds as amplitudes, and we measure intensity.

Where two waves meet, the electric field is the sum of the two fields. This is the principle of superposition, which holds because the wave equation is linear [3]. But neither the eye nor any detector sees the electric field directly. The field of visible light oscillates more than $4\times10^{14}$ times per second, too fast for any detector to follow; a detector records only the time average of the squared field, the irradiance [1, 2]. Since what adds is the amplitude and what is measured is its square, the intensity of two combined beams is not the sum of their intensities [3]. Two waves of intensities $I_1$ and $I_2$ meeting with a phase difference $\delta$ give

$$ I = I_1 + I_2 + 2\sqrt{I_1 I_2}\,\cos\delta $$

[1, 2, 3]. The third term is the interference term. At zero phase difference the intensity exceeds the sum (constructive interference); at $\pi$ it falls below it (destructive interference). If the two intensities are equal, $I = 4I_1\cos^2(\delta/2)$: four times a single beam at the bright points and exactly zero at the dark ones [1, 2].

Phasors make this intuitive. Draw each wave as an arrow whose length is the amplitude and whose direction is the phase, and place the arrows head to tail. The squared length of the resulting arrow is the intensity (Figure 2(a)). Arrows pointing the same way add their lengths; arrows pointing opposite ways subtract.

<img src="/assets/img/posts/wave-optics-interference/en/fig2-two-wave-superposition.png" alt="Phasor sum of two waves, intensity versus phase difference, and visibility versus intensity ratio" width="780">
_Figure 2. (a) Phasor sum of two waves of amplitude 1 and 0.7. At zero phase difference the arrows line up and lengthen; at $\pi$ they cut into each other. (b) Intensity versus phase difference. Averaged evenly over phase, it is always $I_1 + I_2$. (c) Fringe visibility versus the ratio of the two intensities._

Where did the energy go? The energy missing from the dark places has not vanished; it has moved to the bright ones. Averaging over a full cycle of phase in Figure 2(b), the cosine in the interference term averages to zero, so the mean intensity is always $I_1 + I_2$ [1, 3]. Interference neither creates nor destroys light; it only moves it around.

How clear the fringes are is measured by the visibility

$$ V = \frac{I_{\max} - I_{\min}}{I_{\max} + I_{\min}} = \frac{2\sqrt{I_1 I_2}}{I_1 + I_2} $$

With equal intensities $V = 1$ and the dark fringes are completely black [1]. What stands out in Figure 2(c) is how gently the curve falls. Even when one beam is only 1/100 as intense as the other, the visibility is 0.2. Its amplitude is the square root, 1/10, and adding to and subtracting from the strong beam it swings the intensity by ±20%. Techniques that read a very weak signal by mixing it with a strong reference beam and amplifying it through interference exploit exactly this.

## 2. When can interference be seen?

Why, then, don't two lamps in a room produce interference fringes? The interference formula assumes that the phase difference $\delta$ is constant. Three things are needed to see fringes.

First, the two beams must have almost the same frequency. If the frequencies differ, the phase difference rotates rapidly in time and the interference term averages to zero over the detection time [1]. Second, the phase difference must stay constant. Even if the phases differ, a constant difference only shifts the pattern sideways. Two such sources are called coherent [1]. Light from an ordinary source oscillates as a well-behaved sine wave for less than about 10 ns before its phase changes at random [1]. The pattern from two lamps jumps to a new position every such interval and averages away [1]. So to see interference, one splits a single beam in two and recombines it. The two halves change phase together, so their difference stays constant [1, 2]. How long and over how wide a region the phase stays orderly is coherence, the subject of the next part.

Third, the polarizations of the two beams must be parallel. Strictly, the interference term comes from the dot product of the two electric field vectors [1, 2]. If the fields are perpendicular the dot product is zero and the interference term vanishes. Fresnel and Arago found experimentally that two beams polarized at right angles do not interfere, which became evidence that light is a transverse wave vibrating perpendicular to its direction of travel [2]. These are the Fresnel–Arago laws.

There are two ways to split light in two [2]. Division of wavefront takes different parts of a wavefront; division of amplitude splits the intensity at a partially reflecting surface. Division of wavefront requires a small source, while division of amplitude works with extended sources and gives brighter fringes [2]. Sections 3 and 4 take them in turn.

## 3. Division of wavefront: Young's experiment

In 1665 Grimaldi let sunlight through two pinholes to see whether the region where the two discs of light overlapped would darken [1]. A hundred and forty years later the physician Thomas Young repeated the experiment, and the key was a third pinhole placed before the two [1]. The same wavefront spreading from the first hole illuminated both of the others, so the light leaving them formed phase-locked twins. Young published the result in 1804, and it became the experiment that established the wave theory of light against Newton's widely accepted corpuscular theory [1].

Treat the two holes as two point sources. The phase difference at a point in space is set by the difference in its distances from the two sources. Points where that difference is a whole number of wavelengths are bright, and points where it is a half-integer number are dark. Surfaces of constant distance difference from two points are hyperboloids with the points as foci [1]. Figure 3(a) shows the intensity around two point sources six wavelengths apart: bright and dark bands fan out along hyperbolas.

<img src="/assets/img/posts/wave-optics-interference/en/fig3-two-point-sources.png" alt="Interference field of two point sources and the fringes on a distant screen" width="780">
_Figure 3. (a) Intensity around two point sources (green) six wavelengths apart, with the fall-off with distance removed. White dashed lines are the hyperbolas on which the path difference is a whole number of wavelengths. (b) Intensity on a screen 1 m from two holes 0.5 mm apart (wavelength 550 nm), computed by adding the two spherical waves without approximation._

On a distant screen the hyperbolas become nearly straight, and the pattern becomes evenly spaced stripes. With hole separation $d$, screen distance $L$ and wavelength $\lambda$, the spacing between adjacent bright fringes is

$$ \Delta x = \frac{\lambda L}{d} $$

[1, 3]. Two holes 0.5 mm apart lit with 550 nm light give fringes 1.1 mm apart on a screen 1 m away (Figure 3(b)). The wavelength is less than a micrometre, yet the fringes are millimetres wide; the distance $L$ magnifies the wavelength. Measuring the fringe spacing and the two distances gives the wavelength in return.

A more general form hides in this formula. Seen from afar, the light from the two holes is two plane waves meeting at an angle $\theta \approx d/L$. Two plane waves meeting at an angle $\theta$ produce sinusoidal fringes of period $\lambda/\sin\theta$ [3]. At 30° the period is twice the wavelength. This is how fine diffraction gratings are written with light, and holography, which records the fringes between a reference beam and the light from an object, starts from the same idea [3].

<img src="/assets/img/posts/wave-optics-interference/ext-double-slit-sunlight.jpg" alt="Coloured interference fringes from sunlight through two slits" width="480">
_Figure 4. Interference fringes from sunlight passing through two slits about 0.5 mm apart. The central fringe is white because every colour is bright there; further out the spacing $\lambda L/d$ differs from colour to colour and the fringes split into rainbows. (Source: Aleksandr Berdnikov, [Double slit interference](https://commons.wikimedia.org/wiki/File:Double_slit_interference.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

Figure 4 is Young's experiment done directly with sunlight. Sunlight contains every colour. Red interferes with red and blue with blue, and the patterns of all the colours overlap into a single white-light pattern [1]. At the central fringe, where the path difference is zero, every colour is bright together, so it is white. As the path difference grows, the patterns of different colours drift apart, and after a few fringes the colours mix and the fringes vanish. This is why white-light fringes extend over a much narrower range than monochromatic ones, and it returns in the next part as the coherence length.

## 4. Division of amplitude: thin films and Newton's rings

Back to the oil film of Figure 1. Light reflected from the top of the film overlaps light reflected from the bottom. The two surfaces have split the same light by intensity, so this is division of amplitude. The optical path difference between the two is one round trip inside the film, $2nt$ at normal incidence, where $n$ is the film index and $t$ its thickness.

One more thing enters. Reflection at a boundary from the low-index side and reflection from the high-index side have opposite signs. The top reflection of an oil or soap film (air to film) and the bottom reflection (film to air) are exactly these two cases, so an extra phase difference of $\pi$ appears between them [1]. That is why a film thin enough to approach zero thickness looks black: the two reflections cancel. Hecht explains that this $\pi$ is needed for the reflection to go smoothly to zero as the film thins and disappears [1]. A vertical soap film turns black at the top as it drains and thins for exactly this reason [1]. Figure 5(a) computes the reflectance of a water film in air as a function of thickness. It is zero at zero thickness and reaches its first maximum at a thickness of $\lambda/4n$ (111 nm at 589 nm).

<img src="/assets/img/posts/wave-optics-interference/en/fig5-thin-film-newton.png" alt="Reflectance of a water film versus thickness, equal-thickness fringes of a wedge, and computed Newton's rings" width="780">
_Figure 5. (a) Normal-incidence reflectance of a water film (n = 1.33) in air at 589 nm. Orange: only the reflections from the two surfaces added; black dashed: all multiple reflections inside the film added. (b) Reflected light from a wedge-shaped water film tilted by 0.1 mrad. Places of equal thickness form stripes of equal brightness. (c) Newton's rings formed by the air gap under a convex lens of 1 m radius of curvature resting on a flat. The orange dashed circles have radius $\sqrt{m\lambda R}$._

In Figure 5(a) the calculation with only two reflections (orange) and the one with all the multiply reflected light (black dashed) are almost identical. A water–air boundary reflects only 2% of the light, so the third reflection onwards is negligible. Things change when the reflectance is high, as Section 6 shows. Reflection from a stack of many layers is computed with the transfer-matrix method of [Ellipsometry Part 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/).

### Fringes of equal thickness

When the film thickness varies, places of equal thickness form stripes of equal brightness. These are fringes of equal thickness [1, 2]. The coloured bands of Figure 1 are contour lines of the oil film. In white light, different thicknesses strengthen different colours and give coloured bands; in monochromatic light they give bright and dark stripes. In a wedge whose thickness grows uniformly in one direction, as in Figure 5(b), the stripes are evenly spaced. The thickness increases by $\lambda/2n$ from one stripe to the next, so a water film tilted by 0.1 mrad gives stripes 2.2 mm apart. Press two clean microscope slides together and you can see the coloured bands of the air film between them shift as you press [1].

### Newton's rings

Place a convex lens on a flat plate and the air gap between them thickens outwards in proportion to the square of the radius. The fringes of equal thickness become concentric circles, Newton's rings (Figure 5(c), Figure 6). For a lens of radius of curvature $R$, the radius of the $m$-th dark ring is

$$ r_m = \sqrt{m\lambda R} $$

[1]. The centre, where the thickness is zero, is dark. A lens with a 1 m radius of curvature seen in sodium light (589 nm) has its first dark ring at a radius of 0.77 mm. Measuring the ring radii gives the curvature of the lens, and distorted rings mean the surface is uneven [1]. This is the check done in lens shops when a lens surface is laid against a reference.

<img src="/assets/img/posts/wave-optics-interference/ext-newton-rings.jpg" alt="Newton's rings in sodium light" width="560">
_Figure 6. Newton's rings formed in sodium light by the air gap between a 20 cm focal-length convex lens and a flat, seen through a microscope. (Source: Warrencarpani, [20cm Air 1](https://commons.wikimedia.org/wiki/File:20cm_Air_1.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

The rings bear Newton's name, but Hooke and Newton each studied the colours of such thin films in detail, from soap bubbles to the air between lenses [1].

### Films that cancel reflection

Interference can remove reflection as well as create it. Coat a glass surface with a film whose index lies between those of air and glass, with a thickness of $\lambda/4n$: the two reflections from the top and bottom of the film then differ by half a wavelength and cancel. This is an antireflection coating. A bare surface of glass with index 1.52 reflects 4.3%; a quarter-wave coating of magnesium fluoride (MgF₂, index 1.38) reduces that to 1.3% at the design wavelength. The coating is usually tuned to suppress reflection most strongly near green, where the eye is most sensitive. What remains is reflected at the two ends of the spectrum, violet and red, which is why coated camera lenses and spectacle lenses show a purplish tint when viewed at an angle.

## 5. Interferometers

An interferometer splits light with a beamsplitter, sends the two parts along different paths and recombines them [3]. The slightest change in the optical path of either arm shifts the phase difference and moves the fringes, so an interferometer is a ruler that measures length, refractive index and wavefronts in units of the wavelength.

### The Michelson interferometer

The Michelson interferometer splits light with a single half-silvered mirror, sends the two parts to mirrors at the ends of two arms, and recombines the returning light at the same half-silvered mirror (Figure 7) [1, 3]. Moving one mirror by $\lambda/2$ lengthens that arm's round trip by $\lambda$, and one fringe passes. With a 633 nm helium–neon laser, one fringe passes for every 316 nm of mirror travel (Figure 8(a)). Counting fringes measures the mirror's motion to a fraction of a wavelength.

<img src="/assets/img/posts/wave-optics-interference/ext-michelson.jpg" alt="A laboratory Michelson interferometer with the light path" width="520">
_Figure 7. A laboratory Michelson interferometer. Light from the source on the left is split at the central half-silvered mirror, reflected by the mirrors at the top and on the right, recombined, and sent out through the viewing port at the bottom. The compensating plate makes both arms pass through the same thickness of glass. (Source: Warren Leywon, [Michelson interferometer](https://commons.wikimedia.org/wiki/File:Michelson_interferometer.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

<img src="/assets/img/posts/wave-optics-interference/en/fig8-michelson.png" alt="The two outputs of a Michelson interferometer and circular fringes" width="760">
_Figure 8. (a) Intensities of the two outputs of a Michelson interferometer as a mirror moves (lossless 50:50 beamsplitter, 633 nm). When one brightens the other darkens, and their sum always equals the input. (b) Circular fringes when the arms differ in length by 1 mm, as a function of the viewing angle across an extended source._

Figure 8(a) shows two outputs. The recombined light leaves both towards the viewing port and back towards the source. When one is bright the other is dark, and their sum always equals the incoming light, as energy conservation requires. Energy conservation alone shows that this complementarity needs a $\pi/2$ phase difference between the light reflected and transmitted at the beamsplitter [3]. In one output each beam has been reflected once and transmitted once, so the two are in phase; in the other, one beam has been reflected twice and the other transmitted twice, putting them $\pi$ apart [3].

With an extended source, concentric fringes appear as in Figure 8(b). Light entering at an angle from different points of the source has a path difference between the arms that depends on its tilt, and light with the same tilt forms the same ring. These are fringes of equal inclination [1, 2]. If the fringes of equal thickness in Section 4 are a map of thickness, these are a map of angle.

Other interferometers follow from the same idea [1, 3]. The Mach–Zehnder interferometer has separate splitting and recombining mirrors, so its two arms are physically apart and a sample can be placed in one of them. A typical use is to view density variations in gas flow in a wind tunnel [1]. The Sagnac interferometer sends two beams round the same loop in opposite directions; if the whole device rotates, the two optical paths differ, and the rotation rate can be measured. Optical gyroscopes work on this principle [1, 3].

### Seeing aberrations as fringes

Interferometers are used to test lenses. [Geometrical Optics Foundations 4](/en/posts/geometric-optics-aberrations/) defined aberration as the departure $W$ of the wavefront in the exit pupil from a reference sphere. Overlap the wavefront that has passed through a test lens with a perfect reference wave, and the fringes are the contour lines of $W$: from one dark fringe to the next, the wavefront changes by one wavelength. The Twyman–Green interferometer is a Michelson adapted for optical testing. One arm contains the lens under test, followed by a spherical mirror whose centre of curvature coincides with the focus of the lens [1]. If the lens is free of aberrations, the returning light is again a plane wave; if not, the aberration shows up as fringes [1]. Since the light passes through the test lens twice, one fringe corresponds to half a wavelength of wavefront error.

<img src="/assets/img/posts/wave-optics-interference/en/fig9-aberration-interferograms.png" alt="Interferograms of the Seidel aberrations" width="780">
_Figure 9. Interferograms of the aberrations of Geometrical Optics Foundations 4 against a plane reference wave (light passing the test optics once). Top row: reference parallel. Bottom row: reference tilted. One fringe is one wavelength of wavefront difference._

Figure 9 computes how the aberrations of Part 4 look in an interferometer. Two waves of defocus give two concentric rings; three waves of spherical aberration give rings that crowd together towards the edge. Coma gives rings pushed to one side, and astigmatism a saddle-shaped pattern. Deliberately tilting the reference (bottom row) turns the pattern into straight stripes bent by the aberration, which is easier to read. Reading the phase quantitatively from fringes (phase shifting, unwrapping) was the subject of [Image Processing for Metrology Part 6](/en/posts/imgproc-phase-hilbert-unwrapping/). Interferometers and wavefront sensors report the measured wavefront as the Zernike coefficients met in Part 4 [4], and lens makers use those coefficients to polish or align their optics.

### A length smaller than a proton: LIGO

The most extreme use of the Michelson interferometer is the detection of gravitational waves. LIGO in the United States is a Michelson interferometer with 4 km arms (Figure 10). A passing gravitational wave stretches space in one direction and squeezes it in the perpendicular one, changing the difference between the arm lengths. The first gravitational wave detected, in 2015, changed the arm length by about $5\times10^{-22}$ of itself, roughly 2 am (attometres, $10^{-18}$ m) — some 400 times smaller than the radius of a proton [3].

<img src="/assets/img/posts/wave-optics-interference/ext-ligo-hanford.jpg" alt="Aerial view of the LIGO Hanford Observatory" width="680">
_Figure 10. The LIGO Hanford Observatory in Washington State. A 4 km arm runs across the desert from the building at lower right; the other arm runs at right angles to it. (Source: LIGO Laboratory, [LIGO Hanford aerial 05](https://commons.wikimedia.org/wiki/File:LIGO_Hanford_aerial_05.jpg), public domain)_

To catch such a small change, LIGO places a Fabry–Pérot interferometer of Section 6 inside each arm. The light travels up and down the arm hundreds of times, multiplying the phase difference produced by the same change in length. An arm with a finesse of about 450 amplifies phase changes by $2\mathcal{F}/\pi \approx 286$ [3]. This is multiple-beam interference, the subject of the next section.

## 6. Multiple-beam interference: the Fabry–Pérot

So far only two beams have overlapped. Two-beam fringes have a $\cos^2$ shape, with bright and dark bands of equal width. Adding $M$ beams that each differ by the same phase step changes this. Only when the step is a multiple of $2\pi$ do the phasors line up and the intensity reach $M^2$ times one beam; for a slight departure the phasors curl round and cancel. The first zero comes at a step of $2\pi/M$ [3]. The more beams, the narrower the bright peaks (Figure 11(a)). The sharp spectral lines of a diffraction grating come from this.

<img src="/assets/img/posts/wave-optics-interference/en/fig11-fabry-perot.png" alt="Interference peaks for different numbers of beams and the transmittance of a Fabry–Pérot interferometer" width="780">
_Figure 11. (a) Intensity from adding 2, 5 and 20 equal beams with the same phase step between neighbours (normalized to the peak). (b) Transmittance of light bouncing between two mirrors of reflectance $R$. Light passes only when the round-trip phase is a multiple of $2\pi$, and the higher the reflectance, the narrower the peaks._

Trapping light between two parallel mirrors produces such multiple beams naturally. On each round trip a little light leaks out; the leaked beams differ in phase by the round-trip path $2nd$ and weaken by a factor $R$ each time. Adding them all gives the transmittance

$$ T = \frac{(1-R)^2}{1 + R^2 - 2R\cos\phi} = \frac{1}{1 + \frac{4R}{(1-R)^2}\sin^2(\phi/2)}, \qquad \phi = \frac{4\pi n d}{\lambda} $$

[1, 3], known as the Airy function. The device is the Fabry–Pérot interferometer. Fabry and Pérot built it at the end of the nineteenth century; it is a spectroscope of very high resolving power and the basic structure of a laser resonator [1].

The sharpness of the peaks is expressed by the finesse, the spacing between neighbouring peaks divided by their full width at half maximum. For high reflectance

$$ \mathcal{F} \approx \frac{\pi\sqrt{R}}{1-R} $$

[1, 3]. A reflectance of 0.9 gives a finesse of 30, and 0.98 gives 156 (Figure 11(b)). According to Hecht, ordinary Fabry–Pérot instruments have a finesse of about 30, and curved mirrors with dielectric multilayer coatings reach about 1000 [1]. Conversely, with the 4% reflectance of uncoated glass the finesse is below 1, and the pattern is no different from the gentle two-beam fringes of Figure 5(a). The higher the finesse, the more times light bounces between the mirrors, and the more sensitive the phase becomes to a given change in length. That is why LIGO puts Fabry–Pérot cavities in its arms.

With its length fixed, a Fabry–Pérot becomes a filter that selects wavelengths. Peaks appear only at wavelengths for which the round-trip phase is a multiple of $2\pi$, so scanning the mirror spacing finely while measuring the transmitted light reveals the fine structure of spectral lines [1]. Inside a laser, the same principle lets only particular wavelengths survive and be amplified.

## Summary

- Light adds as amplitudes and is measured as intensity. Two beams therefore give $I_1 + I_2 + 2\sqrt{I_1I_2}\cos\delta$, and the third term is interference. Interference neither creates nor destroys energy; it only moves it.
- Fringes require equal frequencies, a constant phase difference and parallel polarizations. Ordinary light loses its phase quickly, so a single beam is split and recombined.
- Young's experiment is the classic division of wavefront. The fringe spacing $\lambda L/d$ magnifies a sub-micrometre wavelength into millimetre-scale fringes.
- Thin films are the classic division of amplitude. The $\pi$ phase on reflection makes very thin films black, and places of equal thickness form the same fringe. Newton's rings measure lens curvature, and antireflection coatings control reflection.
- An interferometer is a ruler graduated in wavelengths. A Michelson interferometer shows one fringe for every $\lambda/2$ of mirror travel, maps the aberrations of a lens as fringes, and with 4 km arms detects length changes of attometres.
- Overlapping many beams sharpens the peaks. The sharpness of a Fabry–Pérot interferometer is set by its finesse, $\pi\sqrt{R}/(1-R)$.

Throughout this post light has been treated as a perfect sine wave whose phase stays orderly forever. Real light is not like that. That the white-light fringes of Figure 4 vanished after a few fringes, and that two lamps make no fringes, both depend on how orderly the phase of light is in time and in space. The next part deals with the concept that measures this: coherence.

## References

1. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Chapter 9 — the interference term and conditions for interference, coherence, the Fresnel–Arago laws, Young's experiment and its history, thin films and fringes of equal thickness, soap films, Newton's rings, the Michelson, Mach–Zehnder and Sagnac interferometers, fringes of equal inclination, multiple-beam interference and the Fabry–Pérot interferometer, the Twyman–Green interferometer (Section 9.8.2).
2. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Sections 7.1–7.6 — interference of two monochromatic waves, the Fresnel–Arago laws and transversality, division of wavefront and of amplitude, localized fringes, multiple-beam interference.
3. B. E. A. Saleh and M. C. Teich, *Fundamentals of Photonics*, 3rd ed. (Wiley, 2019), Section 2.5 — the interference equation and phasors, interferometers and beamsplitter phases, fringes of two tilted plane waves, interference of many waves and finesse, the Fabry–Pérot interferometer, LIGO (Example 2.5-1).
4. R. K. Tyson and B. W. Frazier, *Principles of Adaptive Optics*, 5th ed. (CRC Press, 2022), Section 1.3 — interference and wavefront measurement, the Zernike series.

---
date: 2026-11-09 20:00:00 +0900
layout: post
title: "Ellipsometer Instrumentation 6 — Can One Microscope Objective Measure Many Angles of Incidence at Once?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-back-focal-plane-micro/
page_id: ellipsometer-back-focal-plane-micro
categories: [Optics, Instrumentation]
tags: [ellipsometry, back-focal-plane, angle-resolved, microscopy, instrumentation]
description: A pattern that polarization microscopy treats as a nuisance becomes the signal here. This post follows what it costs to spread the angles of incidence without a goniometer.
math: true
---

[Post 5](/en/posts/ellipsometer-channeled-spectroscopic/) moved the modulation from the time axis to the wavelength axis. A thick birefringent plate becomes a variable retarder once it meets a spectrometer, and the polarization went into the fringes riding on a single spectrum.

This post does the same thing along a **spatial** axis. And it is not only the polarization modulation that moves. **The angle of incidence is spread out as well.**

Every instrument from [Post 1](/en/posts/ellipsometer-oblique-rotating-element/) onward sent light in obliquely. Near the Brewster angle the difference between p- and s-polarized reflection is largest, and that difference is the signal. Sending light in obliquely brings its own problems, though. Changing the angle of incidence means swinging an arm on a **goniometer**, and a beam that strikes the sample at a slant stretches the spot along one direction.

Send the light in **normally** through an objective with a high numerical aperture (NA) and those problems disappear. And the back focal plane (BFP) of that objective holds something unexpected.

## 1. What goes wrong with oblique incidence on a small area?

Cohn's 1994 chapter in the edited volume *Microanalysis of Solids* divides microellipsometry in two. One approach adds focusing optics to a standard ellipsometer and measures a region 10–200 μm across. This is the **microspot** approach. The other replaces the detector with a CCD, which is **full-field imaging**.

Both pay for it. The microspot approach keeps its precision but has to move the sample mechanically from point to point. Cohn's numbers show the difference. A commercial microspot instrument needs **more than 11 hours** to acquire a $100 \times 100$ point image, while a comparatively slow full-field instrument acquires $480 \times 512$ points in **55 seconds**. In exchange, full-field imaging carries a larger uncertainty.

Cohn's table of instruments reported through the 1980s shows something else. An instrument with a $10 \times 30$ μm spot has an uncertainty of $0.1°$ in $\Delta$. A transmission instrument with a 100 μm spot and **eight-zone averaging** reaches $0.005°$. [Post 4](/en/posts/ellipsometer-calibration/) covered four-zone averaging; this instrument went to eight. **Spatial resolution and accuracy trade against each other within a single table.**

Full-field imaging carries one more chronic problem. Cohn puts it this way.

> The **defocusing** encountered when observing at a 45 to 70° angle of incidence affects most two-dimensional imaging microellipsometers.

The sample plane is tilted with respect to the optical axis, so the image plane tilts too and one end falls out of focus. The remedies on offer are wedge prisms for focus correction or scanning stepwise and keeping only the rows with acceptable resolution. Neither is clean.

Normal incidence removes both the defocus and the working-distance constraint. The image plane does not tilt, so there is no defocus. Objective and sample face each other squarely, so **a short working distance no longer risks a collision.** A high-NA objective has a short focal length, and in an oblique arrangement that was exactly the danger.

## 2. Where is the back focal plane, and what does it show?

The back focal plane is the rear focal plane of the objective. **Rays that leave the front focal plane at the same angle converge to a single point there.** It is the counterpart of the property that light from one position on the sample converges to one point in the image plane.

In other words, **position in the back focal plane means angle.** Where the image plane shows where light came from, the back focal plane shows which way it went. The two planes are not conjugate planes. They are related by a Fourier transform, so what is a point in one is a spread distribution in the other. Microscopists look into the back focal plane through a Bertrand lens to align the illumination (lamp filament, condenser aperture) for the same reason: that plane shows angles. Figure 1 shows the relation with rays.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/ext-conoscopy-relay.png" alt="Ray diagram in which two lenses form a directions image and then an image of the object" width="760">
_Figure 1. Ray diagram of a two-lens imaging system. Rays leaving different points of the object in the same direction (same colour) meet at one point in the back focal plane of Lens 1, the "Directions Image", and Lens 2 turns that directions image back into an image of the object. (Source: Panjasan, [Cono-2](https://commons.wikimedia.org/wiki/File:Cono-2.png), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

The difficulty is knowing where that plane is. The author's dissertation records the situation.

> For commercial objectives the manufacturer **often does not provide the exact position of the back focal plane.**

So it is found experimentally. Send a collimated beam into the objective and every ray converges to one small focal spot at the back focal plane, growing by defocus on either side. Mount an imaging system on a micrometer stage, scan along the optical axis, and find where the spot diameter is smallest.

The result explains the layout of these instruments. In the objective the author used, the back focal plane sat **6.1 mm inside the rear face of the objective**, inside the barrel. A camera cannot be placed there, so **a relay lens must image it outside.** That is why instruments built around the back focal plane usually carry one.

## 3. Which angle of incidence does each radius in the back focal plane correspond to?

So which point of the back focal plane belongs to which angle? Radius from the centre maps to angle of incidence, and the mapping comes from the **Abbe sine condition**.

$$\theta = \sin^{-1}\!\big(r\,\mathrm{NA}\big)$$

Here $r$ is the pupil radius normalized to 1. The centre ($r = 0$) is normal incidence and the rim ($r = 1$) is the maximum angle $\sin^{-1}(\mathrm{NA})$. The angle grows with radius, and it grows faster towards the rim: at NA 0.90 the slope at the rim is about 2.3 times that at the centre.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig1-bfp-mapping.png" alt="Radius in the back focal plane maps to angle of incidence" width="780">
_Figure 2. Normalized radius in the back focal plane versus angle of incidence. Left, the curves for several numerical apertures, which steepen towards the rim. Right, iso-angle rings drawn on the back focal plane of an NA 0.90 objective; the rim is at $64.2°$ and the Brewster angle of fused silica (green ring, $55.6°$) lies inside it._

The numerical aperture therefore fixes **the range of angles available**. NA 0.55 reaches $33.4°$, NA 0.90 reaches $64.2°$, NA 0.95 reaches $71.8°$. The numerical aperture is engraved on the barrel next to the magnification (Figure 3).

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/ext-objective-100x.jpg" alt="An oil-immersion objective engraved with its magnification and numerical aperture" width="600">
_Figure 3. An oil-immersion objective engraved with its magnification and numerical aperture. "100×/1,25 Oil" means 100× magnification, NA 1.25 and use in immersion oil; "∞/0.17" means infinity-corrected optics and a 0.17 mm cover glass. A dry objective, used without oil, cannot exceed NA 1, and its practical limit is about 0.95. (Source: QuodScripsiScripsi, [Objective zeiss 100x](https://commons.wikimedia.org/wiki/File:Objective_zeiss_100x.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

That range matters because of the Brewster angle. As [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/) showed, it is the angle at which p-polarized reflection is minimal and where ellipsometry is most sensitive. Fused silica (index 1.458) has a Brewster angle of $55.55°$, so **NA 0.825** reaches it. Glass is similar at $56.31°$.

Silicon is not. With an index of 3.875 at 632.8 nm its Brewster angle is $75.53°$, and reaching it takes **NA 0.968**. Dry objectives stop at about 0.95, so this is practically out of reach in air. The most sensitive angle is unavailable on the main substrate of semiconductor metrology, and that is one practical limit of this method.

For strongly absorbing metals the reference is the pseudo-Brewster angle. Ye and colleagues give $67.9°$ for gold ($N = 0.35 - j2.45$). The three common definitions give $65.55°$ (minimum of $\lvert r_p\rvert$), $69.94°$ (minimum of $\lvert r_p/r_s\rvert$) and $70.38°$ ($\Delta = 90°$), and none of them gives $67.9°$.

## 4. How does the polarization change as you go around in azimuth?

If radius gives the angle, what does azimuth give? This is where the method turns.

The plane of incidence **differs from ray to ray.** For a ray through the point at azimuth $\phi$ in the back focal plane, the plane of incidence is rotated by $\phi$. (In Post 5 $\phi$ was the retardance of a retarder; here it is an azimuth. Each field uses the symbol its own way, so both are left as they are.) Linearly polarized light that has passed a polarizer **fixed** in the laboratory frame therefore sits at a different angle in the sample's p-s frame at every azimuth.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig2-azimuth-rotation.png" alt="Azimuth rotates the polarization axis" width="780">
_Figure 4. How azimuth divides the polarization. Left, the local p direction (red, radial) and s direction (blue, tangential) across the pupil. Right, how light from a polarizer fixed along $x$ splits between p and s. It is entirely p at azimuth $0°$ and entirely s at $90°$._

At azimuth $0°$ the polarization axis lies in the plane of incidence, so the light is **entirely p-polarized**; at $90°$ it is perpendicular, so it is **entirely s-polarized**. In between it splits as $\cos^2\phi$ and $\sin^2\phi$.

In other words, **the polarization axis rotates as the azimuth changes, with no polarizer turned at all.**

This phenomenon was not newly discovered. Polarization microscopy had known it for a long time, and **had a name for it.** The sentence Ye and colleagues wrote in 2007 is the thesis of this post.

> In polarization microscopy this phenomenon is referred to as a Fresnel effect and is **treated as a nuisance because it causes the four-corners problem.** In a microellipsometer system, however, **the incident plane rotating effect plays the role of a rotating polarizer.** This is the key algorithm of real-time acquisition without moving parts.

Microscopy read the pattern as something that spoils image contrast; metrology read it as a signal carrying angle and polarization.

Post 5 already performed this inversion once. In [Post 2](/en/posts/ellipsometer-components/) the wavelength dependence of retardance was a defect to be avoided. In Post 5 the plate was deliberately made thick so that the same dependence carried the modulation. **Turning a defect into the operating principle** is a structure shared by the two posts.

The pattern itself had been analyzed by 1995. Shatalin and colleagues applied **conoscopy**, imaging the exit pupil between crossed polarizers, to isotropic thin films. They explained why the image is not angularly symmetric even for an isotropic sample through the difference in Fresnel coefficients between the TM and TE (that is, p and s) components, and derived the expressions. The paper also notes that ellipsometric data can be extracted from the pattern.

## 5. What can one ring tell you?

If azimuth rotates the polarization axis, then **reading the intensity around a ring of fixed radius** is one measurement taken while the polarization is modulated. Change the radius and it becomes a measurement at a different angle of incidence. That is what "annular" means in the title of Ye's paper.

Which arrangement it is equivalent to needs care, however. Polarizer and analyzer are **both** fixed in the laboratory frame, so in the sample frame **they turn together.** Among the arrangements discussed in Section 2 of [Post 3](/en/posts/ellipsometer-dual-rotating-compensator-mme/), this is the one that rotates polarizer and analyzer together, not the rotating-analyzer arrangement of Post 1.

Multiplying out the Mueller matrices gives

$$I(\phi) \propto \frac{3 + u}{2} - 2\cos 2\Psi\,\cos 2\phi + \frac{1-u}{2}\cos 4\phi, \qquad u = \sin 2\Psi\,\cos\Delta$$

**The sine terms vanish identically and only the $2\phi$ and $4\phi$ cosines survive.** The $2\phi$ term carries $\cos 2\Psi$ alone, and $\cos\Delta$ enters only the constant and the $4\phi$ term. A rotating analyzer, by contrast, produces $2\omega$ alone. This is why Ye's paper writes the Fourier coefficients as $\{\alpha_2, \alpha_4\}$ rather than $\alpha$ and $\beta$. Ye's equations take $\Delta$ $180°$ away from this series, so the sign of $\cos\Delta$ is reversed; keep that in mind when comparing.

With $\alpha_2$ and $\alpha_4$ normalized by the constant term, the inversion is straightforward.

$$\cos 2\Psi = \frac{-\alpha_2}{1+\alpha_4}, \qquad \sin 2\Psi\,\cos\Delta = \frac{1-3\alpha_4}{1+\alpha_4}$$

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig3-annular-readout.png" alt="Pupil intensity and annular demodulation" width="780">
_Figure 5. A synthesized example of annular acquisition (32.9 nm SiO$_2$ on Si, 632.8 nm, NA 0.90). Left, pupil intensity with polarizer and analyzer both along $x$. Centre, intensity around the white ring ($46.1°$ incidence) and its $2\phi$ and $4\phi$ components. Right, $\Psi$ and $\Delta$ recovered at each radius. The data are noise-free, so the recovered values (markers) sit on the true values (pale lines)._

Putting the values in Figure 5(b), $\alpha_2 = -0.705$ and $\alpha_4 = 0.904$, into the inversion gives $\Psi = 34.13°$, the true value. The two coefficients of one ring fix $\Psi$ and $\cos\Delta$.

The four-lobed pattern on the left is the one Ye's paper describes. Dark and bright axes cross. Since polarizer and analyzer are both aligned to $x$, the vertical axis is the brighter one.

The right panel shows what the method buys. **A single image yields $\Psi$ and $\Delta$ from $6°$ to $64°$ of incidence at once.** Work that meant swinging a goniometer and measuring again at each angle becomes a matter of reading a different radius.

What is obtained, though, is $\cos\Delta$. **The sign of $\Delta$ is lost here too.** This is the very arrangement that Section 2 of Post 3 said misses the fourth row and column when polarizer and analyzer turn together. With no component to produce circular polarization, it cannot read $S_3$, just as the rotating analyzer of Post 1 could not. In addition, $\cos\Delta$ has zero slope near $0°$ and $180°$, so sensitivity to $\Delta$ drops where $\Delta$ hugs $180°$, as it does at small angles in the right panel of Figure 5.

Three branches have grown from the question of what to read in the back focal plane.

| Branch | What is read | Representative |
| --- | --- | --- |
| Conoscopy | intensity pattern between crossed polarizers | Shatalin 1995 |
| Annular acquisition | intensity around a ring of fixed radius | Ye 2007 |
| Interferometry | phase distribution in the pupil | Feke 1998 |

The third has an interesting origin. Feke and colleagues measured the **phase** distribution in the back focal plane with a Michelson-type phase-shifting interferometer and obtained $\Delta$ directly. Their motivation was not thin-film metrology. A precision topography interferometer picks up an error on heterogeneous samples because the phase change on reflection differs between materials. **They wanted to measure that phase change and subtract it.** A topography-interferometer manufacturer among the co-authors is the evidence.

## 6. If one camera reads the pupil, where does wavelength go?

That is the gain. The cost appears in the same place as in Post 5.

Reading the back focal plane requires a two-dimensional imaging camera, and **radius and azimuth already occupy both of its axes.** Once angle of incidence and polarization are stored, no axis is left for wavelength. The author's dissertation states the limitation directly.

> An ordinary imaging sensor **cannot resolve wavelength directly.**

Three remedies exist and each costs something: a band-pass filter, an RGB camera with a Bayer filter, or several single-wavelength sources switched in turn. The author's dissertation lists three prices as well. The system becomes more complex, it costs more, and **measurement range and speed are constrained.** In a study that photographed the back focal plane through a Bayer filter, the limited spectral resolution left the fitted parameters markedly less distinguishable.

The data volume grows too. Cohn notes that a full-field imaging instrument needs at least three irradiance images to determine $\Psi$ and $\Delta$ completely, which at $480 \times 512$ pixels and 8-bit depth comes to at least 0.74 MByte for the set (about 0.25 MByte per image). Annular acquisition in the back focal plane gets $\Psi$ and $\cos\Delta$ from a single image. Adding wavelengths in sequence, however, multiplies the image count by the number of wavelengths, and the data needed for one point grow quickly.

Set beside Post 5, the structure is clear.

| | Axis the modulation moved to | What is given up |
| --- | --- | --- |
| Post 5, channeled spectroscopy | wavelength | the spectrometer's resolvable points are shared among channels |
| Post 6, back focal plane | space | radius and azimuth consume both camera axes |

Getting everything at once means giving up an axis. And because the two methods give up **different** axes, there is a way to combine them. The summary returns to this.

## 7. Why does high-NA focusing break the plane-wave assumption?

Normal incidence and a high numerical aperture removed the geometric problems and introduced a new one. What Cohn identified in the same chapter is fundamental.

> The theory of ellipsometry is based on the interaction of electromagnetic **plane waves** at one or more interfaces. When lenses are used to focus the incident light, **the plane wave case no longer directly applies.**

Erman and Theeten, as presented in Cohn's chapter, treated this with Fourier optics. Focusing through a finite aperture gives the angle of incidence a spread, and $r_p$ and $r_s$ are **each convolved** over that spread before their ratio sets the detected $\tan\Psi$. Under coherent illumination amplitude functions are convolved; under incoherent illumination irradiance functions are.

Munro and Török treated the same problem in Mueller form in 2008. The two lines of work do not cite each other, but both say that what gets averaged depends on detection, so they read as the same distinction. A focused field is a sum of plane waves, each meeting the sample matrix once, and **what gets averaged depends on how it is detected.** Confocal detection averages the sample's **Jones** matrices; conventional detection averages the **Mueller** matrices.

Why the distinction matters: the average of Jones matrices is still a single Jones matrix, so it produces no depolarization. The average of Mueller matrices does.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig4-focusing-average.png" alt="The averaging that focusing performs" width="780">
_Figure 6. Averages computed by widening the angular spread about $45°$ incidence (32.9 nm SiO$_2$ on Si, 632.8 nm). Left, the depolarization index: the Jones average (green, solid) stays at 1 and only the Mueller average (red, dashed) falls below it. Right, the error the Mueller average leaves in $\Psi$, which stays below $0.01°$ even at a $12°$ half-width._

The Jones average holds an index of $1.000000000$ exactly, and **only the Mueller average falls below 1**, further as the spread widens. Not one component in the system depolarizes. Yet with conventional detection **the act of focusing itself produces depolarization.**

How much should this worry you? Not much. At a half-width of $\pm 2°$ the index is 0.99978 and the error in $\Psi$ is $0.0004°$. Even at $\pm 10°$ the index is 0.99445 and the error is $0.0075°$. In a microspot instrument that collects a cone this wide onto one detector, that error is of the same order as the uncertainties in the Section 1 table ($0.005$–$0.1°$) and cannot be ignored. Read the back focal plane pixel by pixel, though, and each pixel averages only the narrow spread it covers. At $\pm 0.5°$ the error in $\Psi$ is just $0.00002°$. Focusing-induced averaging is a problem of integrating a wide angular range at once.

The author's dissertation points at a similar gap when criticizing earlier work. Modelling the error from the temperature of a birefringent material and an inaccurate retardance alone **cannot account for the total error of the system**, and the amplitude error of the spectral interference signal was left out. The dissertation traces depolarization to components such as the retarder, the beam splitter and the lenses; focusing-induced averaging is a further source added here. Post 5 covered the retardance side, the drift that rotates the $(a, b)$ plane; the amplitude and depolarization side is what remained.

The remedy already appeared in Posts 4 and 5. Do not try to remove it. **Put it in the forward model.** Post 4's regression calibration said a model cannot account for what is not in it, and in Post 5 Okabe folded the convolution of the window function into the theoretical model. In 2023 Liu and colleagues parameterized the Mueller matrix of the focusing lens pair and fitted it together with the sample.

What the instrument averaged, the model must average too. The idea recurs across Posts 4, 5 and 6.

## 8. If the nominal NA cannot be trusted, how is the angle scale set?

The thesis of [Post 4](/en/posts/ellipsometer-calibration/) was that the angle a dial reports is not the true angle. In Post 5 it was the retardance and the wavelength scale that drifted. Here it is the **angle-of-incidence scale**.

Converting radius to angle needs only the numerical aperture. Yet the number engraved on the barrel, as in Figure 3, cannot be trusted. As the author's dissertation puts it:

> Because of the internal aperture of the objective, **a discrepancy appears between the NA specified by the manufacturer and the NA actually measured.**

Mapping the centre to the maximum-NA radius through the sine condition inherits that discrepancy directly. So what serves as a reference?

**The Brewster angle.** By definition $\tan^2\Psi = R_p/R_s$, and the Brewster angle minimizes $R_p$, hence also $\tan^2\Psi$. Finding that radius $r_b$ in the pupil image fixes one absolute angular reference, and the remaining pixels follow from the sine condition.

$$\theta = \sin^{-1}\!\left(\frac{r\,\sin\theta_b}{r_b}\right)$$

Here $\theta_b$ is the Brewster angle of the reference sample. With one reference point fixed, the numerical aperture is no longer needed.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig5-brewster-calibration.png" alt="Calibrating the angle scale with the Brewster angle" width="780">
_Figure 7. Calibrating the angle scale with the Brewster angle. Left, the $\tan^2\Psi$ profile of a fused silica reference; the minimum at $r_b = 0.916$ corresponds to the Brewster angle of $55.55°$. Right, for a lens whose true NA is 0.90, the scale obtained by trusting the nominal NA of 0.95 (red, dashed) and the scale calibrated with the Brewster angle (green markers), which agrees with the true scale (grey) to about $0.005°$._

The objective the author used was rated NA 0.95 but measured 0.90. Using the rated value as it stands puts the rim **off by up to $7.6°$.** Solving for a thickness from $\Psi$ and $\Delta$ attributed to an angle that is $7.6°$ wrong cannot give the right answer. Setting the scale by the Brewster angle reduces that error to about $0.005°$.

Measuring $\tan^2\Psi$ is simple. Remove the retarder from the polarization state generator, set the polarizer alternately to $0°$ and $90°$, and hold the analyzer at $45°$. The profile is read along the $90°$ azimuth line, where a polarizer along the laboratory $x$ axis ($0°$) delivers s-polarized light, so

$$I_{0°} = I_{in}(1 + \cos 2\Psi) = 2I_{in}\cos^2\Psi, \qquad I_{90°} = I_{in}(1 - \cos 2\Psi) = 2I_{in}\sin^2\Psi$$

The ratio $I_{90°}/I_{0°}$ is $\tan^2\Psi$. **The source intensity cancels, so no absolute radiometric calibration is needed.**

It is the same idea as Post 4's residual calibration, which read three numbers off the minimum of a single curve. Take a ratio to cancel what is unknown, and use only the **position** of the minimum.

The author used a fused silica flat 12.70 mm thick as the reference. It is thick enough for internal reflections to be ignored, and its index of 1.458 gives a Brewster angle of $55.55°$. That index is the value near 588 nm; at the 633 nm measurement wavelength the Brewster angle is about $0.02°$ smaller.

## 9. What is missing from existing back-focal-plane work?

This problem is the subject of the author's doctoral work. Several groups have used the back focal plane and each gained something. The map drawn in the author's dissertation shows where the method stands. "Single-shot" in the table means completing a measurement in one exposure, the same as "snapshot" elsewhere in this series. ○ means satisfied, △ partly satisfied, ✗ not satisfied.

| | Single-shot | Spectral and angular together | Error calibration | Full $\Psi,\Delta$ |
| --- | --- | --- | --- | --- |
| Bayer-filter RGB camera | ✗ | △ | ✗ | ○ |
| Structured-illumination ring | ✗ | △ | ✗ | ✗ |
| Channeled spectroscopy combined | ○ | ✗ | △ | ○ |
| Line-scan plus polarization filter | ○ | ○ | ✗ | ✗ |
| Back focal plane with polarization modulation | ✗ | ○ | ✗ | ○ |

Illuminating a shaped ring along the radial direction gives high spatial resolution. It requires measuring at many radii, however, and as the ring narrows the angular resolution improves while the signal weakens. A line-scan spectrograph with a polarization filter achieves single-shot operation and a microspot. It **cannot measure the phase change fully** and stays reflectance-centred, which makes Mueller-matrix analysis difficult.

The author's dissertation sums up the map as follows.

> Existing studies have not presented an integrated methodology that satisfies **microspot measurement, single-shot analysis, angular resolution, and systematic error calibration** at the same time.

## Summary

In the back focal plane of an objective with a high numerical aperture, radius becomes angle of incidence and azimuth becomes a rotation of the polarization axis. A single image therefore yields $\Psi$ and $\Delta$ at every angle, with no goniometer and no rotating polarization component. A pattern that polarization microscopy treats as a nuisance is the signal here.

Three costs follow. Radius and azimuth consume both camera axes, so **there is no axis left for wavelength.** Focusing breaks the plane-wave assumption and, with conventional detection, **produces depolarization in a system that has none**, though pixel-by-pixel reading of the back focal plane keeps it negligible. And the nominal NA cannot be trusted, so **the angle scale itself has to be rebuilt.**

The three ideas that carried this series are all recovered here. A slot that should be empty watches the instrument; a defect in one context becomes the operating principle in another; and what the instrument averaged, the model must average too.

One link remains. Post 5 loaded the polarization onto the wavelength axis so as not to give that axis up; Post 6 used the spatial axes to spread the angles and lost the wavelength axis instead. **What each method gives up and what it keeps fit together.** Provide a separate wavelength axis with a line-scan spectrograph, carry the polarization modulation on that axis through channeled spectroscopy, and spread the angles in the back focal plane, and all four requirements can be met at once. The author's dissertation attempts that combination, and it is why this series placed Posts 5 and 6 side by side.

Across six posts the question was a single one. **What cannot be measured, and how do you work with that once you know it?**

## References

- S.-H. Ye, S. H. Kim, Y. K. Kwak, H. M. Cho, Y. J. Cho, W. Chegal, "Angle-resolved annular data acquisition method for microellipsometry," *Opt. Express* **15**(26), 18056 (2007) — the origin of annular acquisition.
- R. F. Cohn, "Microellipsometry," in *Microanalysis of Solids*, B. G. Yacobi, D. B. Holt, L. L. Kazmerski, eds. (Plenum, 1994), Ch. 11 — microspot versus full-field imaging, the failure of the plane-wave assumption, and an account of Erman and Theeten's Fourier-optics analysis.
- S. V. Shatalin, R. Juskaitis, J. B. Tan, T. Wilson, "Reflection conoscopy and micro-ellipsometry of isotropic thin film structures," *J. Microsc.* **179**(3), 241 (1995) — the asymmetry of the pupil pattern.
- G. D. Feke, D. P. Snow, R. D. Grober, P. J. de Groot, L. Deck, "Interferometric back focal plane microellipsometry," *Appl. Opt.* **37**, 1796 (1998) — the branch that reads phase.
- P. R. T. Munro and P. Török, "Properties of high-numerical-aperture Mueller-matrix polarimeters," *Opt. Lett.* **33**(21), 2428 (2008) — Jones averaging versus Mueller averaging.
- J. Liu, Z. Jiang, S. Zhang, T. Huang, H. Jiang, S. Liu, "Calibration of polarization effects for the focusing lens pair in a micro-spot Mueller matrix ellipsometer," *Thin Solid Films* **766**, 139656 (2023) — a parametric model for the focusing lens pair.
- Young Joon Kim, "Development of snapshot angle-resolved ellipsometry using a line-scan spectrograph and back focal plane spectral interference," Ph.D. dissertation, Seoul National University, 2025, §1.3 and §4.1–4.3.

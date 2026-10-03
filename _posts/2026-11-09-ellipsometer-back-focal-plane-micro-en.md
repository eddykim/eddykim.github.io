---
date: 2026-11-09 20:00:00 +0900
layout: post
title: "Ellipsometer Instrumentation 6 — The Back Focal Plane: Spreading Angles Without Turning Anything"
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

Send the light in **normally** through a high-numerical-aperture objective and those problems disappear. And the back focal plane (BFP) of that objective holds something unexpected.

## 1. The geometry oblique incidence leaves behind

Cohn's 1994 chapter in the edited volume *Microanalysis of Solids* divides microellipsometry in two. One approach adds focusing optics to a standard ellipsometer and measures a region 10–200 μm across — the **microspot** approach. The other replaces the detector with a CCD — **full-field imaging**.

Both pay for it. The microspot approach keeps its precision but has to move the sample mechanically from point to point. Cohn's numbers show the difference. A commercial microspot instrument needs **more than 11 hours** to acquire a $100 \times 100$ point image, while a comparatively slow full-field instrument acquires $480 \times 512$ points in **55 seconds**. In exchange, full-field imaging carries a larger uncertainty.

Cohn's table of instruments reported through the 1980s shows something else. An instrument with a $10 \times 30$ μm spot has an uncertainty of $0.1°$ in $\Delta$, while a transmission instrument with a 100 μm spot and **eight-zone averaging** reaches $0.005°$. Post 4 covered four-zone averaging; this goes to eight. **Spatial resolution and accuracy trade against each other within a single table.**

Full-field imaging carries one more chronic problem.

> The **defocusing** encountered when observing at a 45 to 70° angle of incidence affects most two-dimensional imaging microellipsometers.

The sample plane is tilted with respect to the optical axis, so the image plane tilts too and one end falls out of focus. The remedies on offer are wedge prisms for focus correction or scanning stepwise and keeping only the rows with acceptable resolution. Neither is clean.

Normal incidence removes both at once. The image plane does not tilt, so there is no defocus, and because objective and sample face each other squarely, **a short working distance no longer risks a collision.** A high-NA objective has a short focal length, and in an oblique arrangement that was exactly the danger.

## 2. The back focal plane sits inside the lens

The back focal plane is the rear focal plane of the objective. **Rays that leave the front focal plane at the same angle converge to a single point there.** It is the counterpart of the property that light from one position on the sample converges to one point in the image plane.

Because of it, **position in the back focal plane means angle.** Where the image plane shows where light came from, the back focal plane shows which way it went. The two planes are not conjugate but related by a Fourier transform, so what is a point in one is a spread distribution in the other. Microscopists look into the back focal plane through a Bertrand lens to align the illumination (lamp filament, condenser aperture) for the same reason: that plane shows angles.

The difficulty is knowing where that plane is. Kim's dissertation records the practice.

> For commercial objectives the manufacturer **often does not provide the exact position of the back focal plane.**

So it is found experimentally. Send a collimated beam into the objective and every ray converges to one small focal spot at the back focal plane, growing by defocus on either side. Mount an imaging system on a micrometer stage, scan along the optical axis, and find where the spot diameter is smallest.

The result explains the layout of these instruments. The measured back focal plane sat **6.1 mm inside the rear face of the objective** — inside the barrel. A camera cannot be placed there, so **a relay lens must image it outside.** That is why every instrument built around the back focal plane carries one.

## 3. Radius is angle of incidence

Radius from the centre of the back focal plane maps to angle of incidence. The mapping comes from the **Abbe sine condition**.

$$\theta = \sin^{-1}\!\big(r\,\mathrm{NA}\big)$$

with $r$ the pupil radius normalized to 1. The centre ($r = 0$) is normal incidence and the rim ($r = 1$) is the maximum angle $\sin^{-1}(\mathrm{NA})$.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig1-bfp-mapping.png" alt="Radius in the back focal plane maps to angle of incidence" width="780">
_Figure 1. Left, the radius-to-angle mapping for several numerical apertures. Right, iso-angle rings drawn on the back focal plane of an NA 0.90 objective._

The numerical aperture therefore fixes **the range of angles available**. NA 0.55 reaches $33.4°$, NA 0.90 reaches $64.2°$, NA 0.95 reaches $71.8°$.

That range matters because of the Brewster angle. As [the background post on Fresnel reflection](/en/posts/ellipsometry-electromagnetic-fresnel/) showed, it is the angle at which p-polarized reflection is minimal and where ellipsometry is most sensitive. Fused silica has a Brewster angle of $55.55°$, so **NA 0.825** reaches it. Glass is similar at $56.31°$.

Silicon is not. Its Brewster angle is $75.53°$ and reaching it takes **NA 0.968**, practically out of reach for objectives used in air. The most sensitive angle is unavailable on the main substrate of semiconductor metrology, and that is one practical limit of this method.

## 4. Azimuth rotates the polarization axis

If radius gives the angle, what does azimuth give? This is where the method turns.

The plane of incidence **differs from ray to ray.** For a ray through the point at azimuth $\phi$ in the back focal plane, the plane of incidence is rotated by $\phi$. (In Post 5 $\phi$ was the retardance of a retarder; here it is an azimuth. Each field uses the symbol its own way, so both are left as they are.) Linearly polarized light that has passed a polarizer **fixed** in the laboratory frame therefore sits at a different angle in the sample's p-s frame at every azimuth.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig2-azimuth-rotation.png" alt="Azimuth rotates the polarization axis" width="780">
_Figure 2. Left, the local p direction (radial) and s direction (tangential) across the pupil. Right, how a fixed polarizer splits between the two components._

At azimuth $0°$ the polarization axis lies in the plane of incidence, so the light is **entirely p-polarized**; at $90°$ it is perpendicular, so it is **entirely s-polarized**. In between it splits as $\cos^2\phi$ and $\sin^2\phi$.

In other words, **the polarization axis rotates as the azimuth changes, with no polarizer turned at all.**

This phenomenon was not newly discovered. Polarization microscopy had known it for a long time, and **had a name for it.** The sentence Ye and colleagues wrote in 2007 is the thesis of this post.

> In polarization microscopy this phenomenon is referred to as a Fresnel effect and is **treated as a nuisance because it causes the four-corners problem.** In a microellipsometer system, however, **the incident plane rotating effect plays the role of a rotating polarizer.** This is the key algorithm of real-time acquisition without moving parts.

Microscopy read the pattern as something that spoils image contrast; metrology read it as angle and polarization encoded into a picture.

Post 5 already performed this inversion once. In [Post 2](/en/posts/ellipsometer-components/) the wavelength dependence of retardance was a defect to be avoided, and in Post 5 the plate was deliberately made thick so that the same dependence carried the modulation. **Turning a defect into the operating principle** is a structure shared by the two posts.

The pattern itself had been analyzed by 1995. Shatalin and colleagues turned **conoscopy** — imaging the exit pupil between crossed polarizers — onto isotropic thin films, explained why the image is not angularly symmetric even for an isotropic sample through the difference in Fresnel coefficients between the TM and TE components, and derived the expressions. That ellipsometric data can be extracted from the pattern is written there too.

## 5. What one ring gives

If azimuth rotates the polarization axis, then **reading the intensity around a ring of fixed radius** is one measurement taken while the polarization is modulated. Change the radius and it becomes a measurement at a different angle of incidence. That is what "annular" means in the title of Ye's paper.

Which arrangement it is equivalent to needs care, however. Polarizer and analyzer are **both** fixed in the laboratory frame, so in the sample frame **they turn together.** Among the arrangements [Post 3](/en/posts/ellipsometer-dual-rotating-compensator-mme/) classified, this is the one that rotates polarizer and analyzer simultaneously — not the rotating-analyzer arrangement of Post 1.

Multiplying out the Mueller matrices gives

$$I(\phi) \propto \frac{3 + u}{2} - 2\cos 2\Psi\,\cos 2\phi + \frac{1-u}{2}\cos 4\phi, \qquad u = \sin 2\Psi\,\cos\Delta$$

**The sine terms vanish identically and only the $2\phi$ and $4\phi$ cosines survive**, in contrast to a rotating analyzer, which produces $2\omega$ alone. This is why Ye's paper writes the Fourier coefficients as $\{\alpha_2, \alpha_4\}$ rather than $\alpha$ and $\beta$.

Inverting the normalized coefficients is then straightforward.

$$\cos 2\Psi = \frac{-\alpha_2}{1+\alpha_4}, \qquad \sin 2\Psi\,\cos\Delta = \frac{1-3\alpha_4}{1+\alpha_4}$$

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig3-annular-readout.png" alt="Pupil intensity and annular demodulation" width="780">
_Figure 3. Left, synthesized pupil intensity for a 32.9 nm SiO$_2$ film. Centre, intensity around the white ring and its two components. Right, $\Psi$ and $\Delta$ recovered at each radius._

The four-lobed pattern on the left is the one Ye's paper describes. Dark and bright axes cross, and since polarizer and analyzer are both aligned to $x$, the vertical axis is the brighter one.

The right panel shows what the method buys. **A single image yields $\Psi$ and $\Delta$ from $6°$ to $64°$ of incidence at once.** Work that meant swinging a goniometer and measuring again at each angle becomes a matter of reading a different radius.

What is obtained, though, is $\cos\Delta$. **The sign of $\Delta$ is lost here too.** The limitation Post 1 traced to a rotating analyzer missing $S_3$ follows along even after the arrangement changes.

Three branches have grown from the question of what to read in the back focal plane.

| Branch | What is read | Representative |
| --- | --- | --- |
| Conoscopy | intensity pattern between crossed polarizers | Shatalin 1995 |
| Annular acquisition | intensity around a ring of fixed radius | Ye 2007 |
| Interferometry | phase distribution in the pupil | Feke 1998 |

The third has an interesting origin. Feke and colleagues measured the **phase** distribution in the back focal plane with a Michelson-type phase-shifting interferometer and obtained $\Delta$ directly, but their motivation was not thin-film metrology. A precision topography interferometer picks up an error on heterogeneous samples because the phase change on reflection differs between materials, and **they wanted to measure that phase change and subtract it.** A topography-interferometer manufacturer among the co-authors is the evidence.

## 6. The cost: the camera has no wavelength axis

That is the gain. The cost appears in the same place as in Post 5.

Reading the back focal plane requires a two-dimensional imaging camera, and **radius and azimuth already occupy both of its axes.** Once angle of incidence and polarization are stored, no axis is left for wavelength.

> An ordinary imaging sensor **cannot resolve wavelength directly.**

Three remedies exist and each costs something: a band-pass filter, an RGB camera with a Bayer filter, or several single-wavelength sources switched in turn. Kim's dissertation records the price — system complexity rises, cost rises, and **measurement range and speed are constrained.** In a study that photographed the back focal plane through a Bayer filter, the limited spectral resolution left the fitted parameters markedly less distinguishable.

The data volume grows too. Cohn notes that a full-field imaging instrument needs at least three irradiance images to determine $\Psi$ and $\Delta$ completely, which at $480 \times 512$ pixels and 8-bit depth comes to at least 0.74 MByte for the set (about 0.25 MByte per image). Annular acquisition in the back focal plane gets $\Psi$ and $\cos\Delta$ from a single image, but adding wavelengths in sequence multiplies the image count by the number of wavelengths, and what a single point costs grows quickly.

Set beside Post 5, the structure is clear.

| | Axis the modulation moved to | What is given up |
| --- | --- | --- |
| Post 5, channeled spectroscopy | wavelength | the spectrometer's resolvable points are shared among channels |
| Post 6, back focal plane | space | radius and azimuth consume both camera axes |

Getting everything at once means giving up an axis. And because the two methods give up **different** axes, there is a way to combine them. More on that at the end.

## 7. The bill for high NA: this is not a plane wave

Normal incidence and a high numerical aperture removed the geometric problems and introduced a new one. What Cohn identified in the same chapter is fundamental.

> The theory of ellipsometry is based on the interaction of electromagnetic **plane waves** at one or more interfaces. When lenses are used to focus the incident light, **the plane wave case no longer directly applies.**

Erman and Theeten treated this with Fourier optics. Focusing through a finite aperture gives the angle of incidence a spread, and $r_p$ and $r_s$ are **each convolved** over that spread before their ratio sets the detected $\tan\Psi$. Under coherent illumination amplitude functions are convolved; under incoherent illumination irradiance functions are.

Munro and Török rebuilt the same distinction in Mueller form in 2008. A focused field is a sum of plane waves, each meeting the sample matrix once, and **what gets averaged depends on how it is detected.** Confocal detection averages the sample's **Jones** matrices; conventional detection averages the **Mueller** matrices.

Why the distinction matters: the average of Jones matrices is still a single Jones matrix, so it does not depolarize. The average of Mueller matrices does.

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig4-focusing-average.png" alt="The averaging that focusing performs" width="780">
_Figure 4. Computed by widening the spread about $45°$ incidence. Left, the depolarization index. Right, the error the Mueller average leaves in $\Psi$._

The Jones average holds an index of $1.000000000$ exactly, and **only the Mueller average falls below 1**, further as the spread widens. Not one component in the system depolarizes, yet **the act of focusing manufactures depolarization.**

Kim's dissertation points at a similar gap when criticizing earlier work: modelling the error from the temperature of a birefringent material and an inaccurate retardance alone **cannot account for the total error of the system**, and the amplitude error of the spectral interference signal was left out. The dissertation traces depolarization to components such as the retarder, the beam splitter and the lenses; focusing-induced averaging is a further source added here. Post 5 covered the retardance side, the drift that rotates the $(a, b)$ plane; the amplitude and depolarization side is what remained.

The remedy already appeared in Posts 4 and 5. Do not try to remove it — **put it in the forward model.** Post 4's regression calibration said a model cannot account for what is not in it, and in Post 5 Okabe folded the convolution of the window function into the theoretical model. Here the Mueller matrix of the focusing lens pair is parameterized and fitted alongside.

**What the instrument averaged, the model must average too** — an idea that recurs across three posts.

## 8. Calibrating the angle scale with the Brewster angle

The thesis of [Post 4](/en/posts/ellipsometer-calibration/) was that the angle a dial reports is not the true angle. In Post 5 it was the retardance and the wavelength scale that drifted. Here it is the **angle-of-incidence scale**.

Converting radius to angle needs only the numerical aperture, and that number cannot be trusted.

> Because of the internal aperture of the objective, **a discrepancy appears between the NA specified by the manufacturer and the NA actually measured.**

Mapping the centre to the maximum-NA radius through the sine condition inherits that discrepancy directly. So what serves as a reference?

**The Brewster angle.** By definition $\tan^2\Psi = R_p/R_s$, and the Brewster angle minimizes $R_p$, hence also $\tan^2\Psi$. Finding that radius $r_b$ in the pupil image fixes one absolute angular reference, and the remaining pixels follow from the sine condition.

$$\theta = \sin^{-1}\!\left(\frac{r\,\sin\theta_b}{r_b}\right)$$

<img src="/assets/img/posts/ellipsometer-back-focal-plane-micro/en/fig5-brewster-calibration.png" alt="Calibrating the angle scale with the Brewster angle" width="780">
_Figure 5. Left, the $\tan^2\Psi$ profile of a fused silica reference. Right, the difference between trusting the nominal NA and calibrating with the Brewster angle._

If the true NA is 0.90 and the specified 0.95 is used as it stands, the rim is **off by up to $7.6°$.** Solving for a thickness from $\Psi$ and $\Delta$ attributed to an angle that is $7.6°$ wrong cannot give the right answer. Setting the scale by the Brewster angle removes that error.

Measuring $\tan^2\Psi$ is simple. Remove the retarder from the generator, set the polarizer alternately to $0°$ and $90°$, and hold the analyzer at $45°$. The profile is read along the $90°$ azimuth line, where a polarizer along the laboratory $x$ axis ($0°$) delivers s-polarized light, so

$$I_{0°} = I_{in}(1 + \cos 2\Psi) = 2I_{in}\cos^2\Psi, \qquad I_{90°} = I_{in}(1 - \cos 2\Psi) = 2I_{in}\sin^2\Psi$$

The ratio $I_{90°}/I_{0°}$ is $\tan^2\Psi$. **The source intensity cancels, so no absolute radiometric calibration is needed.**

It is the same idea as Post 4's residual calibration, which read three numbers off the minimum of a single curve. Take a ratio to cancel what is unknown, and use only the **position** of the minimum.

The reference is a fused silica flat 12.70 mm thick — thick enough for internal reflections to be ignored — with a refractive index of 1.458 and hence a Brewster angle of $55.55°$.

## 9. No study satisfied all four at once

Several groups have used the back focal plane and each gained something. The map Kim's dissertation draws shows where the method stands.

| | Snapshot | Spectral and angular together | Error calibration | Full $\Psi,\Delta$ |
| --- | --- | --- | --- | --- |
| Bayer-filter RGB camera | ✗ | △ | ✗ | ○ |
| Structured-illumination ring | ○ | △ | ✗ | ✗ |
| Channeled spectroscopy combined | ○ | ✗ | △ | ○ |
| Line-scan plus polarization filter | ○ | ○ | ✗ | ✗ |
| Back focal plane with polarization modulation | ✗ | ○ | ✗ | ○ |

Illuminating a shaped ring along the radial direction gives high spatial resolution but requires measuring at many radii, and as the ring narrows the angular resolution improves while the signal weakens. A line-scan spectrograph with a polarization filter achieves snapshot operation and a micro spot, but **cannot measure the phase change fully** and stays reflectance-centred, which makes Mueller-matrix analysis difficult.

> Existing studies have not presented an integrated methodology that satisfies **micro-spot measurement, snapshot analysis, angular resolution, and systematic error calibration** at the same time.

## Summary

In the back focal plane of a high-NA objective, radius becomes angle of incidence and azimuth becomes a rotation of the polarization axis. A single image therefore yields $\Psi$ and $\Delta$ at every angle, with no goniometer and no rotating polarization component. A pattern that polarization microscopy treats as a nuisance is the signal here.

Three costs follow. Radius and azimuth consume both camera axes, so **there is no axis left for wavelength.** Focusing breaks the plane-wave assumption and **manufactures depolarization in a system that has none.** And the nominal NA cannot be trusted, so **the angle scale itself has to be rebuilt.**

The three ideas that carried this series are all recovered here. A slot that should be empty watches the instrument; a defect in one context becomes the operating principle in another; and what the instrument averaged, the model must average too.

One link remains. Post 5 loaded the polarization onto the wavelength axis so as not to give that axis up; Post 6 used the spatial axes to spread the angles and lost the wavelength axis instead. **What each method gives up and what it keeps fit together.** Provide a separate wavelength axis with a line-scan spectrograph, carry the polarization modulation on that axis through channeled spectroscopy, and spread the angles in the back focal plane, and all four requirements can be met at once. Kim's dissertation is that combination, and it is why this series placed Posts 5 and 6 side by side.

Across six posts the question was a single one. **What cannot be measured, and how do you work with that once you know it?**

## Reproducing the calculations

Every number and figure in this post is computed from products of component Mueller matrices alone. The closed forms and published numbers are kept out of the synthesis and reserved as an independent path for the verification script.

- `bfp.py` — pupil intensity synthesis, annular demodulation, focusing average, Brewster calibration
- `verify_bfp.py` — checks against Ye's 2007 numbers, the round trip of the Abbe sine condition, the harmonic decomposition of the ring signal, Munro and Török's Jones and Mueller averages, and Eq. 4.2 of Kim's dissertation. All nine groups pass
- `generate_figures.py` — Figures 1–5 in both languages

The full code sits in `_code/ellipsometer-back-focal-plane-micro/`.

This cross-check caught **a misidentified arrangement.** The first attempt treated one ring as equivalent to the rotating-analyzer measurement of Post 1 and expected the form $I = I_0[1 + \alpha\cos 2\phi + \beta\sin 2\phi]$, which did not fit. The coordinate relation worked out in Section 5 had been overlooked at the planning stage. Expanding the Mueller product confirmed that only the $2\phi$ and $4\phi$ cosines survive, and the inversion was solved again. Why Ye's paper writes $\{\alpha_2, \alpha_4\}$ only became clear at that point.

One number was not reproduced. Ye's paper gives the pseudo-Brewster angle of gold ($N = 0.35 - j2.45$) as $67.9°$, whereas sweeping for the minimum of $\lvert r_p\rvert$ gives $65.6°$. The paper may use a different definition, so the value is left out of the body. The $56.3°$ for glass at $n = 1.5$ and the $64.16°$ for NA 0.9 reproduce exactly.

## References

- S.-H. Ye, S. H. Kim, Y. K. Kwak, H. M. Cho, Y. J. Cho, W. Chegal, "Angle-resolved annular data acquisition method for microellipsometry," *Opt. Express* **15**(26), 18056 (2007) — the origin of annular acquisition.
- R. F. Cohn, "Microellipsometry," in *Microanalysis of Solids*, B. G. Yacobi, D. B. Holt, L. L. Kazmerski, eds. (Plenum, 1994), Ch. 11 — microspot versus full-field imaging, and the failure of the plane-wave assumption.
- S. V. Shatalin, R. Juskaitis, J. B. Tan, T. Wilson, "Reflection conoscopy and micro-ellipsometry of isotropic thin film structures," *J. Microsc.* **179**(3), 241 (1995) — the asymmetry of the pupil pattern.
- G. D. Feke, D. P. Snow, R. D. Grober, P. J. de Groot, L. Deck, "Interferometric back focal plane microellipsometry," *Appl. Opt.* **37**, 1796 (1998) — the branch that reads phase.
- P. R. T. Munro and P. Török, "Properties of high-numerical-aperture Mueller-matrix polarimeters," *Opt. Lett.* **33**(21), 2428 (2008) — Jones averaging versus Mueller averaging.
- J. Liu, Z. Jiang, S. Zhang, T. Huang, H. Jiang, S. Liu, "Calibration of polarization effects for the focusing lens pair in a micro-spot Mueller matrix ellipsometer," *Thin Solid Films* **766**, 139656 (2023) — a parametric model for the focusing lens pair.
- Young Joon Kim, "Development of snapshot angle-resolved ellipsometry using a line-scan spectrograph and back focal plane spectral interference," Ph.D. dissertation, Seoul National University, 2025, §1.3 and §4.2–4.3.

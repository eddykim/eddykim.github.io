---
layout: post
title: "Computational Imaging and Ptychography 10 — Can Phase Measure Thickness? EUV Imaging Reflectometry"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-euv-imaging-reflectometry/
page_id: ptycho-euv-imaging-reflectometry
date: 2026-12-19 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, euv, reflectometry, thin-film, metrology, phase-retrieval]
description: "Following, by calculation, an EUV imaging reflectometer that reads height, thickness and composition of thin films from the complex reflectance recovered by reflection ptychography."
math: true
---

Up to post 9 this series looked at specimens. We measured how closely a phase image resembled the true specimen with a single error, and called an error of 0.01 a good reconstruction. Semiconductor metrology asks something else. How many nanometres thick is this nitride film? How many ångströms is the error bar? Is the oxide next to it the same thickness? The question is not whether the image looks right, but whether a number can be read from each of its pixels.

Ptychography already holds the material to answer this, because its result is a complex number, amplitude and phase, not an intensity. In 2021 Tanksalvala et al. used this in extreme ultraviolet (EUV) reflection. They tilted the sample through several angles, recovered complex reflectance images by reflection ptychography, and solved for layer thicknesses and dopant levels pattern by pattern. They called the instrument a phase-sensitive imaging reflectometer.

This post asks four questions. Why EUV, and why reflection? When the sample is tilted, how is the diffraction pattern distorted, and how is it straightened out? What can be measured with the phase that cannot be measured with the amplitude alone? And in transmission, can a single complex transmittance identify a material? The calculations follow the conditions of Tanksalvala et al. (wavelength 29.3 nm, grazing angles 21–25°), with optical constants from the Henke tables at CXRO.

The sign convention first. As in [Ellipsometry Foundations 1](/posts/ellipsometry-electromagnetic-fresnel/), the time factor is $e^{j\omega t}$. The EUV and X-ray communities write the complex refractive index as $n = 1 - \delta - j\beta$, where $\delta$ is how far the index falls below 1 and $\beta$ is the absorption. This is the same convention as $N = n - jk$ in the ellipsometry series. At these wavelengths both $\delta$ and $\beta$ are around 0.01, so writing the index as a small departure from 1 is convenient.

## 1. Why EUV, and Why Reflection?

What goes wrong when a visible-light microscope looks at a 10 nm film? Resolution stops at about half the wavelength. A shorter wavelength would help, but lenses cannot be made for EUV. Every material absorbs this light strongly, so there is no transparent lens material, and single-surface mirrors reflect so little that multilayer mirrors are needed. Even a Mo/Si multilayer at 13.5 nm reaches at most 70%. With image-forming optics this awkward, ptychography, which computes the image from diffraction patterns without a lens, suits the band well. The source is available too. Focusing a femtosecond laser into a gas produces EUV light by high harmonic generation (HHG), giving a coherent beam on a tabletop. Why coherence matters was covered in [Wave Optics Foundations 2](/posts/wave-optics-coherence/).

Why reflection rather than transmission? EUV light does not get far into matter. The penetration depth $\lambda/4\pi\beta$, over which the intensity falls by $1/e$, is 28 nm for both SiO$_2$ and Si$_3$N$_4$ at 29.3 nm, and 274 nm for Si. A wafer cannot be seen through. Reflection, on the other hand, uses this short depth to pick up only the top few tens of nanometres. That is exactly the depth thin-film metrology cares about. Figure 1 shows the instrument of Tanksalvala et al.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/ext-tanksalvala2021-setup.jpg" alt="Layout of the EUV imaging reflectometer, the sample structure, phase reconstructions before and after improvements, a wide-field amplitude reconstruction, and reflectivity and reflection phase curves for several materials" width="950">
_Fig 1. Phase-sensitive EUV imaging reflectometer. (A) Multilayer and ellipsoidal mirrors bring the beam to the sample at a glancing angle; sample and detector rotate in a θ-2θ geometry. (B, C) Phase reconstructions before and after improvements. (D) Wide-field amplitude reconstruction. (E, F) Reflectivity and reflection phase of several materials at 30 nm. Source: M. Tanksalvala et al., Science Advances 7, eabd9667 (2021), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)_

In Figure 1(A) the beam meets the surface at a few tens of degrees. In this post, as in Tanksalvala et al., the angle is **measured from the surface**: the grazing angle $\theta$. The ellipsometry series measured from the normal, so take care. Why such a glancing angle? EUV reflectivity drops quickly as the beam approaches normal incidence. Figure 2 shows the calculation.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/en/fig2-reflectivity.png" alt="s-polarized reflectivity and reflection phase of Si, Si with native oxide, SiO2 and Si3N4 at 29.3 nm as functions of grazing angle" width="950">
_Fig 2. Computed s-polarized reflectivity and reflection phase at 29.3 nm. Between 21° and 25° the values separate strongly by material_

Si reflects 77% at a grazing angle of 15°, but only 9% at 25° and 1% at 40°. Adding 1.5 nm of native oxide changes the reflectivity at 21° from 35% to 23%. A nanometre of film moves the reflectivity by more than ten percentage points, which is why EUV reflection is so sensitive to thin films. The curves agree in magnitude and shape with Figure 1(E, F). The Henke tables, however, were built for 50–30,000 eV, and 29.3 nm (42 eV) lies below that range, so their accuracy is reduced. Real measurements, like those of Tanksalvala et al., take the measured values as the reference. The numbers here show trends.

Now the phase on the right. At a grazing angle of 0° every material starts at 180°: a beam at a very glancing angle is almost totally reflected, with its sign flipped, whatever it hits. As the angle rises the materials separate. At 23°, Si is at 17°, SiO$_2$ at 90° and Si$_3$N$_4$ at 98°. Tanksalvala et al. chose 21–25° for this separation. The reflectivity there is still a few percent, enough to collect photons, while the differences between materials are large.

## 2. How Does Tilting the Sample Distort the Diffraction Pattern?

From post 3 to post 9 the specimen faced the beam squarely. A detector pixel then corresponded directly to one cell of the specimen's spatial frequencies. The diffraction pattern was the Fourier transform of the exit wave, the detector pixels formed a uniform grid, and one FFT took us between the two. What changes when the sample is tilted?

In reflection, what a pixel receives is the **component parallel to the surface** of the scattering vector $\mathbf{q} = \mathbf{k}_{\text{out}} - \mathbf{k}_{\text{in}}$, because the reflectance map $r(x, y)$ lives on the surface. Place the detector at distance $L$, normal to the specular direction, with coordinates $u$ in the plane of incidence and $v$ across it. Near the centre,

$$ q_x \approx k \sin\theta\, \frac{u}{L}, \qquad q_y \approx k\, \frac{v}{L} $$

where $k = 2\pi/\lambda$ and $x$ is the incidence direction. Along the incidence direction, then, a given distance on the detector spans a spatial-frequency range smaller by $\sin\theta$. Since pixel size is the inverse of that range, reconstructed pixels are $1/\sin\theta$ times longer along $x$. At 22°, $1/\sin 22° = 2.67$. Tanksalvala et al. reported pixels of 64 nm by 172 nm, a ratio of 2.69, which is this relation at work. For the same reason the illuminated footprint is stretched by $1/\sin\theta$ along the incidence direction. This anisotropy, called conical diffraction, widens the field of view at the cost of resolution in that direction.

The trouble is that the equation above holds only near the centre. Towards the edge of the detector, $q_x$ is no longer proportional to $u$. Figure 3 shows the exact mapping at 22° for a detector NA of 0.2.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/en/fig3-geometry.png" alt="Layout of incident beam, specular beam and detector in reflection ptychography, and the curved grid of surface spatial frequencies received by the detector pixels" width="950">
_Fig 3. Geometry of reflection ptychography (left) and the surface spatial frequencies received by the detector pixels at 22° (right). The grey square is the extent of the uniform grid used in reconstruction_

The blue lines on the right are the rows and columns of detector pixels. The square grid is bent, asymmetrically, along the incidence direction. Under the linear approximation $q_x$ would be spread evenly between $\pm 0.075k$; in fact it runs from $-0.056k$ to $+0.107k$. Pixels at the edge receive spatial frequencies 25–43% away from where the linear approximation places them.

The remedy is to compute, before reconstruction, where each point of a uniform grid lands on the detector, and to interpolate the measured intensity there. This is the tilted-plane correction. Points that fall off the detector have no data, so those cells are left free in the reconstruction. In Figure 3 they are the left strip of the grey square outside the blue grid; at 22° and NA 0.2 they make up 21% of the uniform grid. On the right, conversely, the detector receives spatial frequencies beyond the uniform grid, which the reconstruction cannot use.

How much does the correction matter? The same test pattern was reconstructed three ways: measured directly on a uniform grid (the ideal case), measured with the exact mapping and then corrected, and measured with the exact mapping but interpreted with the linear approximation. The probe was taken as known here, so that only the geometry differs.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/en/fig4-tilted-plane.png" alt="Phase of the true sample and of reconstructions with and without tilted-plane correction, the phase-step bias versus detector oversampling, and reconstruction error versus detector NA" width="950">
_Fig 4. Reconstructed phase from reflection ptychography at 22° (top); without correction it fails. Bottom left: bias that intensity interpolation in the correction leaves in region phase steps, shrinking as the detector is sampled more finely. Bottom right: reconstruction error versus NA_

At NA 0.2 the corrected reconstruction has an error of 0.059 and the uncorrected one 0.69. As in the top right of Figure 4, structures slide along the incidence direction and the phase is scrambled. The bottom-right panel varies the NA, with the test pattern kept the same in pixel units. The corrected error stays at 0.055–0.069 regardless of NA. The uncorrected error grows from 0.14 at NA 0.05 to 0.34 at 0.10 and 0.69 at 0.20.

Why does the linear approximation fail even at small NA? The relative size of the nonlinear term at the edge works out to roughly $\mathrm{NA}\cot\theta/2$, and $\cot\theta$ grows as the grazing angle falls. At 22°, $\cot\theta = 2.5$, so even at NA 0.05 the edge frequencies are off by 6%. With a 32-pixel half-width, that is a shift of about two pixels. In reflection ptychography, the tilted-plane correction is not optional. Tanksalvala et al. likewise reported that an exact correction including all three sample rotation angles was critical to image quality. In Figure 1, (B) is before and (C) after; (C) also includes the total variation (TV) regularization of [post 1](/posts/ptycho-computational-imaging-inverse-problem/).

Correction is not the end of the story. Metrology cares less about whether the image looks right than about whether region averages are accurate. The bottom-left panel of Figure 4 shows how far the region phase step (defined in section 3) departs from the truth after reconstruction from corrected, noise-free data with a known probe. With detector pixels at the same spacing as the reconstruction grid, a bias of 0.04–0.11 rad remains at 22–40°. The cause is the intensity interpolation. A sharp-edged probe pupil gives the diffraction intensity sharp edges too, and cubic interpolation fails near them. Unlike noise, this error lands on the same side however many times one measures. Sampling the detector twice as finely cuts the bias to 0.006–0.02 rad, and four times as finely to 0.002 rad or less. An accurate correction needs a detector finer than the reconstruction itself requires. All calculations from section 3 on use the fourfold detector.

## 3. What Does a Phase Step Measure?

Suppose the reconstruction has succeeded. What does its phase measure? Here the ambiguities of [post 4](/posts/ptycho-probe-ambiguity-raster-rpie/) return. There we called scale exchange, phase ramp and translation effectively harmless. They are harmless for looking at an image, but not for metrology. If the recovered reflectance carries an arbitrary complex constant $c$, its magnitude and phase on their own mean nothing.

What survives is the **ratio between regions**. Dividing two regions A and B of the same image removes $c$:

$$ \frac{r_B}{r_A} = \left\lvert \frac{r_B}{r_A} \right\rvert e^{j\Delta\varphi} $$

Here $\lvert r_B/r_A \rvert$ is the reflectance ratio of the two regions and $\Delta\varphi$ the phase step. The phase ramp is removed by fitting a plane so that a broad flat region (here the substrate under its native oxide) has flat phase. The reference for the measurement must lie inside the sample. It is the same idea as ellipsometry measuring the ratio $r_p/r_s$ to cancel the source intensity and the absolute phase.

A phase step mixes two things. One is geometric height. If a surface is higher by $h$, the round-trip path is shorter by $2h\sin\theta$, so the phase advances by

$$ \Delta\varphi_{\text{height}} = 2k h \sin\theta = \frac{4\pi h \sin\theta}{\lambda} $$

The factor 2 is there because the light goes down and comes back up. At 29.3 nm and 25°, a height of 2 Å gives 0.036 rad. Bringing the phase noise down to around 0.01 rad makes ångström heights visible, and the axial precision of 2 Å reported by Tanksalvala et al. is of this order. The other is the reflection phase of the layer stack. Where nitride and oxide are stacked, reflections from the top surface and from the buried interfaces interfere, giving a phase that depends on material and thickness as in the right of Figure 2. This is computed with Rouard's method from [Ellipsometry Foundations 3](/posts/ellipsometry-thin-film-multilayer-reflectance/), known on the X-ray side as Parratt's (1954) recursion: start at the bottom interface and fold the reflection coefficient upward one layer at a time. Here every region is referred to a common plane, the top of the substrate, so that height differences and stack phases enter one expression together.

Now follow the whole measurement. The test sample has three regions: the reference substrate (Si with 1.5 nm native oxide), a structure (5 nm SiO$_2$ on Si with 10 nm Si$_3$N$_4$ on top) and a pad (15 nm SiO$_2$ on Si). Poisson noise was added at a scale of $10^8$ photons per position; the detected photons are reduced by the reflectivity, to $1.6\times10^7$–$6\times10^7$ in the brightest pattern. After tilted-plane correction, rPIE recovered object and probe together over 200 iterations. The phase ramp was then removed using the substrate, and each region, trimmed by two pixels at its edges, was averaged to give the ratios. The angles were the five of Tanksalvala et al. (21–25°) and, for comparison, five spread more widely (15–35°). There were three noise draws per angle.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/en/fig5-phase-steps.png" alt="Phase images reconstructed at 21° and 25° relative to the substrate, and reflectance ratio and phase step of the structure region versus angle" width="950">
_Fig 5. Reflectance ratio and phase step of the structure relative to the substrate, from reconstructions at each angle. Dots are 3 noise draws; lines are the model with the true parameters_

The two images at the top are the same sample reconstructed at 21° and 25°. With the substrate (dark) set to zero, brightness is the phase step. The structure and pad look different at the two angles; this angle dependence is the fingerprint of the layer stack. In the plots below, the recovered ratios (dots) sit on the model (lines). Over 27 reconstructions the reflectance ratio deviated from the truth by at most 0.3% and the phase step by at most 0.010 rad, with medians of 0.07% and 0.003 rad. Reconstruction errors were 0.021–0.035. A phase noise of about 0.002 rad corresponds to 0.1 Å of height at 25°.

The brightness of the reference region deserves a word, and it is why the wide range stops at 35°. With the same photon budget at 40° and 45°, the phase step is off by up to 0.022 rad for the structure and 0.090 rad for the pad, against at most 0.010 and 0.032 rad over 15–35°. The substrate reflectivity falls to 0.4% at 40° and 0.2% at 45°, and the reference region becomes dark. A dark region collects few photons, and, as in post 6, rPIE follows the noise where photons are scarce. When the reference average wobbles, every ratio wobbles with it. The reference region should reflect at least a few percent at every measurement angle. This is part of why Tanksalvala et al. chose 21–25°, where the reflectivity is 3–15%.

## 4. Can Reflectance Ratios and Phase Steps Separate Thickness from Composition?

This is where a second inverse problem begins. The first turned diffraction patterns into a reflectance image. The second turns reflectance ratios and phase steps at several angles into layer thicknesses and composition. Take three unknowns for the structure: the nitride thickness, the oxide thickness and a density factor for the nitride. The density factor scales $\delta$ and $\beta$ together, standing in for a film whose composition has shifted. Each angle gives one complex number, two real values. Five angles give ten measurements for three unknowns. By count, that looks ample.

Is it? As in [post 1](/posts/ptycho-computational-imaging-inverse-problem/), having more measurements than unknowns does not make a problem easy. What matters is how differently the measurements see different directions in the space of unknowns. The tool for this is the Fisher information matrix. From the Jacobian $J$ of the measurements with respect to the unknowns and the measurement noise $\sigma$, form $F = J^{\mathsf T}J/\sigma^2$. The diagonal of $F^{-1}$ is the lower bound on the variance of any unbiased estimate, the Cramér–Rao lower bound. The noise $\sigma$ was set for each angle from the spread of the three noise draws in Figure 5. With the same noise, actual fits were also run: least squares from 20 starting points, keeping the lowest cost.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/en/fig6-fit.png" alt="Confidence ellipses and fit results for nitride thickness and density factor using amplitude ratio only, phase step only and both, for narrow and wide angle ranges" width="950">
_Fig 6. 1σ confidence ellipses (information matrix) for nitride thickness and density, with the actual fits. Over a narrow angle range the amplitude ratio alone ties the two together_

The ellipses are the 1σ regions given by the Cramér–Rao bound; the symbols are actual fits to the three noise draws (× for one kind of data, ● for both). The table collects the results as the lower bound on the standard deviation of the nitride thickness.

| Angles | Amplitude ratio only | Phase step only | Both |
|---|---|---|---|
| Narrow (21–25°, 5) | 1.96 nm | 11.6 nm | 0.53 nm |
| Wide (15–35°, 5) | 0.28 nm | 0.76 nm | 0.13 nm |

Start with the narrow range on the left. With the amplitude ratio alone, the correlation between thickness and density is −0.993; the ellipse lies along a diagonal. A thicker nitride and a lower density are nearly indistinguishable in these data, and the actual fits scatter to 8.9, 12.0 and 9.7 nm. The phase step alone is worse: correlation −1.000 and a standard deviation of 11.6 nm. The fits land at 14.0, 5.9 and 6.1 nm, and some fall outside the plot once the density factor is included.

Combined, the bound drops to 0.53 nm and the fits gather at 9.2–10.2 nm. Why are two weak measurements strong together? The two ellipses lie in nearly the same direction, but at slightly different slopes. Where the amplitude ratio says "add this much thickness and remove this much density and nothing changes", the phase step says "not quite this much". The two valleys cross at only one place. The phase is useful not because it is a strong measurement on its own, but because it looks in a different direction from the amplitude.

Over the wide range on the right, the amplitude ratio alone already gives 0.28 nm. With a wide spread of angles, the shape of the reflectance curve itself separates thickness from density. The number of measurements is still five. Making measurements see different things matters more than making more of them. Together the bound is 0.13 nm, and the actual fits, at 10.15–10.26 nm, are slightly offset. Once the bound is down to 0.13 nm, the noise estimate from only three draws and small biases left by interpolation begin to show.

For comparison, fitting the pad, which has a single unknown (its SiO$_2$ thickness), gives 15.02, 14.96 and 15.02 nm even over the narrow range, against a true 15 nm. The same data are easy with one unknown and hard with three. This is why Tanksalvala et al. limited the number of unknowns to the number of phase steps they measured and fixed the rest at design values. They expected that more angles would allow more parameters to be solved together.

## 5. In Transmission, Can One Complex Transmittance Identify a Material?

In reflection we had to change the angle several times. Transmission offers a simpler route. In 2022 Eschen et al. imaged a lamella from an integrated circuit in transmission with tabletop ptychography at 13.5 nm and identified the material pixel by pixel. How can one know the material without knowing the thickness?

The transmittance of a film of thickness $t$ is $T = e^{-jk(n-1)t}$. Substituting $n = 1 - \delta - j\beta$ separates phase and log amplitude:

$$ \varphi = k\delta t, \qquad -\ln\lvert T \rvert = k\beta t, \qquad \frac{\varphi}{-\ln\lvert T \rvert} = \frac{\delta}{\beta} $$

Here $\varphi$ is the phase advance through the film. Both quantities are proportional to thickness, so their ratio cancels it and leaves only the material property $\delta/\beta$. Eschen et al. called this ratio the scattering quotient. From the Henke tables at 13.5 nm it is 2.88 for Si$_3$N$_4$, 2.04 for SiO$_2$ and $-0.10$ for Al, matching the 2.9, 2.0 and $-0.1$ tabulated by Eschen et al. at 92 eV. Al is negative because its L absorption edge lies just below, at 73 eV, pushing $\delta$ below zero.

To check this by calculation, three material strips were placed side by side, each a wedge whose thickness increases in one direction (30–150 nm for Si$_3$N$_4$ and SiO$_2$, 10–60 nm for Al). The leftmost strip was left empty as vacuum. Noise was added at the same photon scale as in reflection, the probe was recovered too, and the phase and log amplitude of each pixel were plotted. Transmission needs no tilted-plane correction, so the data were measured directly on a uniform grid.

<img src="/assets/img/posts/ptycho-euv-imaging-reflectometry/en/fig7-scattering-quotient.png" alt="Layout of three material wedges and vacuum, and scatter plots of pixel phase against log amplitude for the vacuum-calibrated and the unreferenced reconstruction" width="950">
_Fig 7. Pixel-wise phase and log amplitude of a 13.5 nm transmission reconstruction. Calibrated to vacuum, each material lies on a line through the origin (centre); without a reference the points scatter (right)_

In the centre panel, pixels of different thickness line up on one straight line per material. The slope of each line is the scattering quotient, and the dashed lines are the tabulated values. The median pixel ratios are 2.88 for Si$_3$N$_4$, 2.03 for SiO$_2$ and $-0.10$ for Al, in agreement with the tables. The reconstruction error was 0.005.

The key is the vacuum. The lines pass through the origin only after the scale, global phase and phase ramp of the reconstruction have been set so that the vacuum region transmits exactly 1. Multiply the same reconstruction by an arbitrary complex constant (right) and the lines leave the origin; the ratios become 0.85, 0.47 and $-0.91$, which are meaningless. Just as the substrate served as the reference in reflection, vacuum is the reference in transmission. The practical advice that follows: leave an empty hole in the field of view when designing the sample. Eschen et al. reported that the materials separated cleanly only when the model also allowed the probe to drift during the scan. Their method added orthogonal probe relaxation, which lets the probe vary slightly from position to position, to the mixed states of [post 7](/posts/ptycho-position-correction-mixed-states/). The probe here was fixed, so that effect is not covered.

## 6. How Does This Compare with an Ellipsometer?

Readers of the ellipsometer series will recognize the flow. Measure a ratio, compute reflection coefficients with a multilayer model, and find the parameters that minimize the difference between measurement and model. Side by side:

| | Spectroscopic ellipsometer | EUV imaging reflectometer |
|---|---|---|
| Ratio measured | $r_p/r_s$ (same spot, two polarizations) | $r_B/r_A$ (same polarization, two regions) |
| What the ratio cancels | Source intensity, absolute phase | Scale and global phase of the reconstruction |
| Variable scanned | Wavelength (hundreds) | Angle (about five) |
| Spatial resolution | Averaged over the spot (10 µm or more even for microspot) | Per pixel (about 100 nm) |
| Depth | Up to µm for transparent films | Penetration depth, tens of nm |
| Model fitting | Hundreds of measurements | About ten measurements |

The imaging reflectometer trades the ellipsometer's number of measurements for spatial resolution. An ellipsometer averages over a large area and collects hundreds of spectral values, so it can afford many unknowns. An imaging reflectometer gives a value for every 100 nm pixel, but a handful of angles yields few measurements. That is why Tanksalvala et al. limited their unknowns to the number of measured phase steps and fixed the rest from design values or other measurements. Section 4 shows why: when measurements are few, the angle range and the presence of phase decide what can be solved.

The back-focal-plane micro-ellipsometer of [Ellipsometer Instrumentation 6](/posts/ellipsometer-back-focal-plane-micro/) sits between the two. It reads many angles of incidence at once in the back focal plane of an objective, measuring a small area while turning angles into extra measurements. The ellipsometer is moving towards finer spatial resolution, the imaging reflectometer towards more measurements, and they are approaching each other.

## Summary

EUV light is hard to focus with lenses and penetrates only tens of nanometres. These two weaknesses make lensless reflection ptychography a natural tool for thin-film metrology. Imaging in reflection, though, means the tilt of the sample must enter the reconstruction. The spatial frequencies received by the detector pixels are nonlinear along the incidence direction, and without a tilted-plane correction that maps them onto a uniform grid the reconstruction fails (error 0.69 against 0.059).

The absolute magnitude and phase of the recovered reflectance cannot be trusted because of the ambiguities of post 4. Ratios between regions can. A phase step measures height as $4\pi h\sin\theta/\lambda$ and also carries the reflection phase of the layer stack. When thickness and composition are solved together over a narrow range of angles, neither the reflectance ratio nor the phase step separates them alone (nitride thickness bounds of 1.96 nm and 11.6 nm). Together the bound drops to 0.53 nm, because the phase looks in a different direction from the amplitude. A wider range of angles gives 0.28 nm with the amplitude alone. And as the bounds shrink, errors that are not noise, such as interpolation in the correction and a dark reference region, come to decide the result. In transmission, calibrating to vacuum turns the ratio of phase to log amplitude into a material fingerprint independent of thickness.

## Closing the Series

Post 1 started from the inverse problem of undoing a blurred photograph. Post 2 moved to phase retrieval, recovering phase from intensity alone, and post 3 showed the core of ptychography: overlapping, shifted illumination solves the phase problem. Posts 4 and 5 dealt with the ambiguities that appear when the probe is also unknown, and with rPIE, momentum and the difference map. Posts 6 and 7 brought the imperfections of real experiments, photon noise, position errors and partial coherence, into the cost function and the model. Post 8 turned the same mathematics towards illumination angles in Fourier ptychography, and post 9 tied everything together in one framework of automatic differentiation and multislice models.

In this last post the image became a measurement. Looking back, the same questions recurred through all ten posts. What does the measurement contain, and what has it erased? What fills in what was erased? Is the mathematically correct answer the physically correct one? Ptychography is one elaborate answer to these questions: recover the phase from the redundancy of overlap, put the imperfections of the experiment into the model, and read numbers against a reference inside the sample.

## References

- M. Tanksalvala, C. L. Porter, Y. Esashi, B. Wang, N. W. Jenkins, Z. Zhang, G. P. Miley, J. L. Knobloch, B. McBennett, N. Horiguchi, S. Yazdi, J. Zhou, M. N. Jacobs, C. S. Bevis, R. M. Karl Jr., P. Johnsen, D. Ren, L. Waller, D. E. Adams, S. L. Cousin, C.-T. Liao, J. Miao, M. Gerrity, H. C. Kapteyn, M. M. Murnane, "[Nondestructive, high-resolution, chemically specific 3D nanostructure characterization using phase-sensitive EUV imaging reflectometry](https://doi.org/10.1126/sciadv.abd9667)," *Science Advances* 7(5), eabd9667 (2021). Phase-sensitive imaging reflectometer, tilted-plane correction, fitting phase steps across angles, 64×172 nm pixels and 2 Å axial precision.
- W. Eschen, L. Loetgering, V. Schuster, R. Klas, A. Kirsche, L. Berthold, M. Steinert, T. Pertsch, H. Gross, M. Krause, J. Limpert, J. Rothhardt, "[Material-specific high-resolution table-top extreme ultraviolet microscopy](https://doi.org/10.1038/s41377-022-00797-6)," *Light: Science & Applications* 11, 117 (2022). Transmission ptychography at 13.5 nm, material identification by the scattering quotient δ/β, mixed-state orthogonal probe relaxation.
- L. G. Parratt, "[Surface studies of solids by total reflection of X-rays](https://doi.org/10.1103/PhysRev.95.359)," *Physical Review* 95(2), 359–369 (1954). Recursive calculation of multilayer reflection.
- B. L. Henke, E. M. Gullikson, J. C. Davis, "[X-ray interactions: photoabsorption, scattering, transmission, and reflection at E = 50–30,000 eV, Z = 1–92](https://doi.org/10.1006/adnd.1993.1013)," *Atomic Data and Nuclear Data Tables* 54(2), 181–342 (1993). Optical constant tables; values obtained from [CXRO](https://henke.lbl.gov/optical_constants/).

---
title: "Computational Imaging and Ptychography 8 — Can We See More Sharply Without Changing the Lens?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-fourier-ptychography/
page_id: ptycho-fourier-ptychography
date: 2026-12-11 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, fourier-ptychography, synthetic-aperture, phase-retrieval, microscopy]
description: "Fourier ptychography, which turns many low-resolution images taken at different illumination angles into an image beyond the lens's numerical aperture, solved with the same mathematics as ptychography."
math: true
---

Choosing a microscope objective always means giving something up. A lens with a large numerical aperture (NA) resolves well but sees a small field; a lens with a small NA sees widely but blurs. The number of pixels a single image can hold, the space-bandwidth product, is set by the lens. In 2013, Zheng, Horstmeyer and Yang showed how to get around this limit by computation. They replaced the illumination of a 2×, 0.08 NA objective with an LED matrix and took many low-resolution images while changing the illumination angle. Combining those images gave a resolution equivalent to 0.5 NA over the wide field of the 2× lens, a space-bandwidth product of about 0.9 gigapixels.

This method is called Fourier ptychography. As the name suggests, it uses the same mathematics as the ptychography of the previous posts. This post asks three questions. Why does changing the illumination angle bring in information the lens missed? How can the algorithms so far be reused as they are? And how do the lessons of the earlier posts, overlap, probe recovery and noise models, reappear here? The calculations use a 256×256 high-resolution specimen recorded by a camera with a coarser 64×64 grid.

## 1. Does Tilting the Illumination Bring In What the Lens Missed?

Consider what a lens does in Fourier space. Splitting the light that has passed through the specimen by angle gives the specimen's spatial-frequency spectrum. An objective collects only light leaving at angles smaller than its numerical aperture, so it passes only a circle of the spectrum whose radius is proportional to the NA. This circle is the pupil. Higher spatial frequencies outside the circle, the fine structure, do not reach the image.

What happens when the specimen is illuminated obliquely? Multiplying by a tilted plane wave $e^{-i \mathbf{k}\cdot\mathbf{r}}$ is the same as shifting the spectrum by $\mathbf{k}$ in Fourier space. With the Fourier transform defined by $e^{-i\mathbf{q}\cdot\mathbf{r}}$, this plane wave turns $\tilde{O}(\mathbf{q})$ into $\tilde{O}(\mathbf{q}+\mathbf{k})$; the formula below and the code follow this convention. The pupil circle stays put while the spectrum moves, so parts that used to lie outside now fall inside the circle.

$$ I_{\mathbf{k}}(\mathbf{r}) = \left\lvert \mathcal{F}^{-1}\left\{ P(\mathbf{q})\, \tilde{O}(\mathbf{q} + \mathbf{k}) \right\} \right\rvert^2 $$

Here $\tilde{O}$ is the specimen spectrum, $P$ the pupil, and $\mathbf{k}$ the shift introduced by the illumination. For each illumination angle, the camera records a different circle of the spectrum as a low-resolution image. That this formula equals multiplying the specimen by a tilted plane wave in real space, filtering with the lens, and sampling coarsely with the camera was checked separately (they agree to $10^{-15}$).

<img src="/assets/img/posts/ptycho-fourier-ptychography/en/fig1-fp-principle.png" alt="Three pupil circles on the specimen spectrum, and low-resolution images under normal, oblique brightfield, and darkfield illumination" width="950">
_Fig 1. For each illumination angle the lens pupil (white circle) passes a different part of the spectrum. Off-centre, the image becomes darkfield_

Fig 1 shows images for three illuminations. Under normal and slightly oblique illumination, the centre of the spectrum, the light that passed straight through the specimen, stays inside the pupil and the image is brightfield with a bright background. When the illumination tilts more than the pupil allows, the centre leaves the circle, only light scattered by the specimen gets in, and the image is darkfield with a dark background. A darkfield image carries high frequencies the lens could not see on its own, but by itself it is not higher in resolution, because a single circle is still the same size.

## 2. How Is This the Same Problem as Ptychography?

A single circle does not raise the resolution, so the circles need to be stitched together. The camera records only intensity, though, so each circle loses its phase. This is the phase problem of Post 2 onwards. And if the illumination angles are chosen so that neighbouring circles overlap, the same thing happens as when overlapping positions solved the phase problem in Post 3. Placing the formulas side by side makes the correspondence clear.

| Ptychography | Fourier ptychography |
|---|---|
| Specimen $O(\mathbf{r})$ | Specimen spectrum $\tilde{O}(\mathbf{q})$ |
| Probe $P(\mathbf{r})$ | Pupil $P(\mathbf{q})$ |
| Scan position $\mathbf{R}_j$ | Illumination angle $\mathbf{k}$ |
| Diffraction pattern | Camera image |
| Overlap of neighbouring positions | Overlap of neighbouring pupil circles |

It is the same problem with the roles of real space and Fourier space exchanged. The authors of PtyLab stress that this reciprocity lets both methods be solved with the same code. In fact the rPIE of Post 4 can be used almost unchanged: put the spectrum in place of the specimen and the pupil in place of the probe, and replace the Fourier transform by its inverse. The inverse transform is the forward transform of the conjugate, so the two are mathematically the same update. Comparing one update computed both ways, they agree to $10^{-16}$.

## 3. How Far Does Stitching the Circles Widen the Aperture?

The pupil radius was set to 12 (Fourier pixels), and the LEDs were placed on a 17×17 square grid with spacing 10. The outermost LED shifts the spectrum by 80, so the highest frequency that can be captured along an axis is $80 + 12 = 92$, about 7.7 times what the lens alone reaches. This enlarged aperture is called a synthetic aperture. rPIE was run for 50 iterations, starting from the normally illuminated image.

<img src="/assets/img/posts/ptycho-fourier-ptychography/en/fig2-synthetic-aperture.png" alt="Amplitude and phase of the true specimen, one camera image, the phase within the lens NA, and the amplitude and phase recovered by Fourier ptychography" width="900">
_Fig 2. High-resolution specimen recovered from 289 low-resolution images. The highest spatial frequency along an axis is about 7.7× the lens's_

One camera image (top middle) has 64×64 pixels and is as blurred as the lens's aperture allows. Even if that image had recorded phase, it could hold no more than the bottom middle. The reconstruction (right) is close to the true specimen in both amplitude and phase. Compared within the part of the spectrum covered by the LED circles, the error is 0.0046. Against the whole true specimen it is 0.043; the difference is the highest frequencies the circles did not cover. Fourier ptychography recovers exactly as much as the circles cover, and within that it is nearly exact. Note also that it recovered the phase, which the camera never recorded.

## 4. How Much Must the Pupil Circles Overlap?

In Post 3, positions in real space had to overlap enough. How much must the circles in Fourier space overlap? The LED spacing was varied to reduce the overlap area of neighbouring pupil circles from 79 % to 0 %. The shift produced by the outermost LED was kept similar, between 72 and 84.

<img src="/assets/img/posts/ptycho-fourier-ptychography/en/fig3-overlap.png" alt="Image error versus overlap area of neighbouring pupil circles, for three LED processing orders" width="750">
_Fig 3. Image error versus the overlap area of neighbouring LED pupils. Below about 40 % it fails abruptly_

The result comes down to a single threshold. With an overlap of 39 % or more the reconstruction succeeds. Processing the LEDs from the centre outwards gives errors all below 0.005; with a random order the error at 39 % rises slightly to 0.01–0.02. At 79 % it is $10^{-8}$. At 35 % it starts to break down, 0.16–0.41, and at 30 % or less it fails, 0.54–0.70. Changing the order in which the LEDs were processed did not move the threshold. The LEDs should be placed so that around 40 % of the circle area overlaps. The overlapping part of neighbouring circles is the link that ties two images to the same phase; if that link is too narrow, the phases of the circles cannot be brought into agreement.

## 5. Can the Lens Aberration Be Recovered Too?

Real lenses, especially inexpensive wide-field low-magnification ones, have aberrations: a phase on the pupil. What happens when the aberration is ignored and a perfect pupil is assumed? With a mix of defocus and astigmatism at two strengths (0.63 and 1.27 rad rms), the reconstruction error was 0.49–0.96: every case failed.

Fortunately the cure is already known. The pupil sits in the place of the ptychographic probe. Just as ePIE in Post 3 recovered the probe along with the specimen, the pupil can be corrected along with the spectrum. Ou and colleagues called this embedded pupil recovery (EPRY). In this post's implementation the update is the probe update of rPIE with the roles of specimen and pupil exchanged.

<img src="/assets/img/posts/ptycho-fourier-ptychography/en/fig4-pupil-recovery.png" alt="True and recovered pupil phase, and the phase reconstructed with the pupil unknown and with the pupil recovered" width="950">
_Fig 4. For a strongly aberrated lens (1.27 rad rms), recovering the pupil too restores both aberration and specimen (79 % overlap)_

| Pupil overlap | Aberration 0.63 rad rms | Aberration 1.27 rad rms |
|---|---|---|
| 49 % | 0.49 → 0.28 (pupil 0.46 rad) | 0.80 → 0.63 (pupil 1.2 rad) |
| 69 % | 0.88 → 0.28 (pupil 0.30 rad) | 0.63 → 0.48 (pupil 0.20 rad) |
| 79 % | 0.96 → 0.052 (pupil 0.049 rad) | 0.86 → 0.083 (pupil 0.075 rad) |

In Fig 4, the phase reconstructed with the pupil unknown (third panel) has lost its structure, while recovering the pupil too (last panel) gives a sharp specimen. The recovered pupil phase (second panel) has the same shape as the true one (first panel).

The table gives the image error with the pupil unknown → recovered, and in parentheses the phase error left in the recovered pupil. A phase tilt of the pupil produces the same images as a shift of the specimen and cannot be fixed by the data, so it was removed before measuring. It is the Fourier-space counterpart of the phase-tilt ambiguity of Post 4.

<img src="/assets/img/posts/ptycho-fourier-ptychography/en/fig5-pupil-vs-overlap.png" alt="Image error with the pupil unknown and with the pupil recovered, versus overlap of the pupil circles" width="750">
_Fig 5. Image error with the pupil unknown or recovered, versus overlap. Recovering the pupil needs much more overlap_

In Fig 5, the blue lines (pupil recovered too) drop sharply only at 79 % overlap. Pupil recovery depends strongly on overlap. With the pupil known, 39 % was enough, but recovering the pupil as well needed 79 %. At 49 % the pupil was hardly recovered. This is the lesson of Post 3 again, though Post 3 measured overlap as a linear fraction of the diameter and this post as a fraction of the circle area, so the numbers cannot be compared directly. There, too, 56 % overlap was enough with the probe known, while ePIE, which also recovers the probe, needed 72 % or more. With more unknowns, neighbouring measurements have to tie each other together more tightly. The review by Konda and colleagues also points out that EPRY easily falls into local minima when the aberrations are large.

## 6. Which Cost Function Is Better When Bright and Dim Images Are Mixed?

Fourier ptychography data have one feature of their own: the brightness of brightfield and darkfield images differs enormously. In a setting where the brightest image receives $10^5$ photons, brightfield images average $9.5\times10^4$ photons while darkfield images average 179. That is more than a factor of 500. Almost all of the high-frequency information lies in these dim darkfield images.

In [Post 6](/en/posts/ptycho-noise-maximum-likelihood/), the cost function decided the outcome for such data. Yeh and colleagues compared several Fourier ptychography algorithms on experimental data and reported that amplitude-based costs are more robust to noise and model mismatch than intensity-based costs, because the large fluctuations of bright brightfield images dominate an intensity cost and bury the darkfield signal. The same comparison was run here.

<img src="/assets/img/posts/ptycho-fourier-ptychography/en/fig6-noise-costs.png" alt="Image error of rPIE and ML with amplitude, Poisson and intensity costs versus dose" width="750">
_Fig 6. Image error by cost function under noise (3 noise draws). Amplitude and Poisson ML are lowest; the intensity cost and rPIE are 1.2–1.7× higher (relative to the lower of amplitude and Poisson ML)_

| Photons in the brightest image | rPIE | ML (amplitude) | ML (Poisson) | ML (intensity) |
|---|---|---|---|---|
| $10^5$ | 0.215 | 0.165 | 0.188 | 0.20 |
| $10^6$ | 0.107 | 0.084 | 0.086 | 0.101 |
| $10^7$ | 0.046 | 0.037 | 0.034 | 0.057 |

Amplitude and Poisson ML are lowest, and the intensity cost and rPIE are 1.2–1.7 times higher. This factor is relative to the lower of amplitude and Poisson ML; relative to amplitude ML alone it is 1.2–1.6. Amplitude and Poisson ML being lowest matches Post 6, but the gap between the intensity cost and rPIE is much smaller than in Post 6 (where the intensity cost was more than four times worse), and their order changes with dose. There, intensity least squares stayed near 0.4; here, starting from rPIE, it comes down to around 0.2.

Another source of error Yeh and colleagues identified is wrong LED positions. If the LED matrix is slightly rotated or at the wrong height, the illumination angles $\mathbf{k}$ are off. This is the same problem as the scan-position errors of [Post 7](/en/posts/ptycho-position-correction-mixed-states/), and the angles are found from the data in the same way.

## Summary

Tilting the illumination shifts the specimen spectrum, and the lens's pupil circle passes a different part of it. Low-resolution images taken at different illumination angles are overlapping circles of the spectrum, and stitching them together gives a high-resolution image beyond the lens's numerical aperture. Stitching phase-less circles through their overlap is the same problem as ptychography with the roles of real and Fourier space exchanged, so the algorithms so far carry over directly.

The lessons of the earlier posts came back as well. The pupil circles must overlap by about 40 % of their area. Lens aberrations are corrected by recovering the pupil like a probe, which needs close to 80 % overlap. For data mixing bright brightfield and dim darkfield images, amplitude or Poisson costs beat the intensity cost. In this setting, the reconstruction recovered frequencies 7.7 times higher along an axis than the lens alone, and was nearly exact within the covered band.

## Next

So far, PIE, DM, ML, position correction, mixed states and Fourier ptychography have been covered one at a time. Post 9 views them all as a single optimization problem. Once the forward model is written to be differentiable, automatic differentiation supplies the gradients, and the earlier methods fit into one framework simply by changing the cost and the unknowns. The multislice model for thick specimens, where the multiplicative approximation breaks down, is treated in the same framework.

## References

- G. Zheng, R. Horstmeyer, C. Yang, "[Wide-field, high-resolution Fourier ptychographic microscopy](https://doi.org/10.1038/nphoton.2013.187)," *Nature Photonics* 7, 739-745 (2013). LED-matrix illumination, a synthetic NA of 0.5 from a 0.08 NA lens, a space-bandwidth product of about 0.9 gigapixels, and the need for overlap for convergence.
- L.-H. Yeh, J. Dong, J. Zhong, L. Tian, M. Chen, G. Tang, M. Soltanolkotabi, L. Waller, "[Experimental robustness of Fourier ptychography phase retrieval algorithms](https://doi.org/10.1364/OE.23.033214)," *Optics Express* 23(26), 33214 (2015). Noise, aberrations and LED misalignment, and the greater robustness of amplitude-based costs.
- P. C. Konda, L. Loetgering, K. C. Zhou, S. Xu, A. R. Harvey, R. Horstmeyer, "[Fourier ptychography: current applications and future promises](https://doi.org/10.1364/OE.386168)," *Optics Express* 28(7), 9603 (2020). Reciprocity with ptychography, embedded pupil recovery (EPRY) and its limits.
- L. Loetgering et al., "[PtyLab.m/py/jl: a cross-platform, open-source inverse modeling toolbox for conventional and Fourier ptychography](https://doi.org/10.1364/OE.485370)," *Optics Express* 31(9), 13763 (2023). Reciprocity that lets both methods be solved with the same code.
- X. Ou, G. Zheng, C. Yang, "[Embedded pupil function recovery for Fourier ptychographic microscopy](https://doi.org/10.1364/OE.22.004960)," *Optics Express* 22(5), 4960 (2014). The original EPRY paper, checked against the open full text (PMC4086333).

---
title: "Computational Imaging and Ptychography 3 — Why Moving the Illumination Solves the Phase Problem: PIE and ePIE"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-overlap-pie-epie/
page_id: ptycho-overlap-pie-epie
date: 2026-11-21 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, pie, epie, phase-retrieval, overlap]
description: "Why recording overlapping diffraction patterns makes the phase problem easy, checked by writing PIE and ePIE from scratch."
math: true
---

[Post 2](/en/posts/ptycho-phase-problem-iterative-retrieval/) ended with the problem of recovering phase from a single diffraction pattern. With enough measurements an answer exists, but it was hard to find and fragile against noise. The original and the twin image could also "fight" each other. The post closed with one hint: recording separate regions of an object with a little overlap makes the twin image disappear.

This post follows that hint. Recording diffraction patterns from overlapping regions while moving the illumination, and using them together to reconstruct the object, is ptychography. The name comes from the Greek ptycho, "to fold". In 1969 Hoppe noted that when a narrow beam illuminates a crystal, the diffraction pattern is the crystal's Bragg reflections convolved with the transform of the illumination, a convolution being *Faltung* (folding) in German. Neighbouring reflections overlap and interfere, so their phase difference can be read from the interference, and he argued that the remaining conjugate ambiguity is removed by moving the illumination once more. The name itself was coined by Hoppe with Hegerl in 1970.

The idea is old, but a practical algorithm came much later: the PIE (ptychographical iterative engine) published by Rodenburg and Faulkner in 2004. According to Rodenburg and Maiden, the name also teases an eminent scientist who, in the early 1990s, called ptychography "pie in the sky". This post writes PIE and its successor ePIE from scratch and quantifies why overlapping measurements make such a difference.

## 1. What Ptychography Measures

Start with what is measured. If the specimen is thin enough, the exit wave, the wave just after the specimen, is the product of the illumination and the specimen's transmission function. The illumination is called the probe. At position $j$, with the probe shifted by $R_j$,

$$ \psi_j(r) = P(r - R_j)\, O(r), \qquad I_j(k) = \lvert \mathcal{F}\{\psi_j\}(k) \rvert^2 $$

where $P$ is the probe, $O$ the specimen, and $I_j$ the diffraction intensity recorded by a far-field detector at that position. One such pattern accumulates per position.

<img src="/assets/img/posts/ptycho-overlap-pie-epie/en/fig1-setup.png" alt="Scan positions over the specimen and diffraction patterns at individual positions" width="900">
_Fig 1. The probe is stepped with overlap, and one diffraction pattern is recorded at each position_

The computations follow the 2009 simulation conditions of Maiden and Rodenburg. The probe is the field 200 µm behind a pinhole of radius 12 µm (wavelength 633 nm, pixel 1 µm). It spreads as it propagates and develops Fresnel fringes at its edge. Since 90 % of its energy lies within a diameter of 36 pixels, that is used as the probe diameter here. The specimen carries different images in its amplitude (0.2–1) and phase (0–π). The scan is a regular grid with each position jittered by ±1 pixel; Maiden and Rodenburg warned that a perfectly regular grid leaves a periodic pattern in the reconstruction, a topic for Post 4. Fig. 1 shows an 8-pixel step, 12 × 12 = 144 positions.

In the language of [Post 1](/en/posts/ptycho-computational-imaging-inverse-problem/), this is again an inverse problem. The difference is that each pattern draws on a different set of unknowns. The same specimen pixel appears in many patterns, each time multiplied by a different part of the probe. That repetition is the protagonist of this post.

## 2. Dividing the Exit Wave by the Probe

As in Post 2's ER, pick a position, form the exit wave from the current estimate, and replace its magnitude in the diffraction plane with the measurement. This gives an exit wave $\psi_j'$ consistent with the data. How is the specimen extracted from it?

The first idea is to divide: since $\psi' = P\,O$, $O = \psi' / P$. But as Rodenburg and Maiden point out, the probe is not a sharp disc. It weakens towards its edge and drops almost to zero in the troughs of its Fresnel fringes, and dividing there blows up. The red curve in Fig. 2 shows how much division amplifies noise: a factor of 50 just 30 pixels from the centre.

<img src="/assets/img/posts/ptycho-overlap-pie-epie/en/fig2-pie-weight.png" alt="Probe amplitude and phase, and the gain of plain division and the PIE weight versus distance from the probe centre" width="900">
_Fig 2. Instead of dividing by the probe, PIE corrects only as much as the probe is bright_

This is exactly the inverse-filter problem of Post 1, and the remedy is the same. Multiply numerator and denominator by $P^*$ to make the denominator real, then add a small constant $\alpha$. Rodenburg and Faulkner themselves call this effectively a Wiener filter. One more ingredient is added: strongly illuminated regions put more information into the diffraction pattern, so they are corrected more, and weakly illuminated regions less. The result is the PIE update.

$$ O(r) \leftarrow O(r) + \frac{\lvert P \rvert}{\lvert P \rvert_{\max}} \frac{P^*}{\lvert P \rvert^2 + \alpha}\, \beta \left( \psi_j' - \psi_j \right) $$

What it does is simple. Of the correction that division would give, it applies a fraction $\lvert P \rvert / \lvert P \rvert_{\max}$. Where the probe is brightest it divides fully; as the probe weakens it corrects less; where nothing was illuminated it changes nothing. That is why the applied fraction (blue dashes) in Fig. 2 nearly coincides with the $\lvert P \rvert / \lvert P \rvert_{\max}$ curve; they part only in very dark regions where $\lvert P \rvert^2$ becomes comparable to $\alpha$. Rodenburg and Faulkner used $\alpha = 10^{-4}$ and $\beta = 1$, and so does this post.

If the illumination is a sharp disc, the update amounts to pasting the data-consistent values into the disc. With a single position as well, the original paper states, PIE is mathematically identical to Post 2's ER. Running PIE and ER ten times each with one disc illumination, the exit waves agree to $10^{-12}$. PIE is Post 2's method extended to many positions.

## 3. Running PIE

Assume the probe is known exactly and run PIE. One iteration uses each of the 144 positions once, in random order.

<img src="/assets/img/posts/ptycho-overlap-pie-epie/en/fig3-pie-convergence.png" alt="Error versus iteration for PIE and ePIE, and the amplitude and phase of the true specimen and the PIE reconstruction" width="900">
_Fig 3. With a known probe PIE converges in tens of iterations; with a wrong probe it stalls_

The specimen error falls from 0.47 after one iteration to 0.039 after 10, 0.0033 after 50, and 0.00096 after 100. The two images placed in amplitude and phase come back separately, without leaking into each other. The contrast with Post 2 is large. No support had to be known in advance, no twin image appeared, and what HIO needed hundreds of iterations for took tens. Nor did it stagnate the way ER did.

The reason is the overlap. Taken alone, any single pattern is matched equally well by the twin image or a shifted image. But where it overlaps the neighbouring pattern, both must use the same specimen values. If one position drifts towards the twin, it no longer agrees with its neighbours. Overlap lets neighbours hold each other's ambiguities in place.

## 4. How Much Overlap?

If overlap is the key, how much is needed? The step was varied from 4 to 24 pixels while covering the same area. Wider steps reduce the number of positions from 484 to 25, and with it the overlap.

<img src="/assets/img/posts/ptycho-overlap-pie-epie/en/fig4-overlap.png" alt="PIE and ePIE error versus overlap, and reconstructions at low and high overlap" width="900">
_Fig 4. With a known probe half overlap suffices; recovering the probe as well needs more_

| Step (pixels) | Overlap (90 % diameter) | Overlap (pinhole diameter) | PIE, 100 it. | ePIE, 200 it. |
|---|---|---|---|---|
| 4 | 0.89 | 0.83 | 0.0001 | 0.006 |
| 8 | 0.78 | 0.67 | 0.0010 | 0.048 |
| 10 | 0.72 | 0.58 | 0.0022 | 0.048 |
| 12 | 0.67 | 0.50 | 0.0048 | 0.18 |
| 16 | 0.56 | 0.33 | 0.020 | 0.66 |
| 20 | 0.44 | 0.17 | 0.079 | 0.81 |
| 24 | 0.33 | 0.00 | 0.24 | 0.88 |

Overlap is one minus the step divided by the probe diameter. PIE with an exactly known probe has plenty of margin: the error is 0.02 even at 56 % overlap. When overlap nearly vanishes it breaks down. The first image in Fig. 4 is at 33 % overlap. It splits into disconnected islands of different brightness, because each island was solved on its own, like Post 2's problem, with its own scale and phase.

The two overlap columns measure the same step in two ways. Taking the pinhole diameter (24 pixels) as the probe size, a step of 24 means no overlap; but the real probe spreads as it propagates, and by the 90 %-energy criterion it still overlaps by 33 %. Maiden and Rodenburg wrote that good results need an overlap of around 60–70 %. Such figures should be read together with how the probe size was defined.

## 5. When the Probe Is Unknown Too

So far the probe has been known exactly. What about real experiments? Maiden and Rodenburg note that modelling the probe accurately, phase included, is at best very time-consuming and at worst impossible to measure accurately enough. What happens to PIE if the probe is slightly wrong?

Telling PIE that the probe is a featureless disc 2 pixels wider than the pinhole gives the red curve in Fig. 3: the error stays near 0.6. Whatever error the wrong probe introduces through the division goes straight into the specimen.

The ePIE (extended PIE) of Maiden and Rodenburg (2009) corrects the probe as well. The probe and the specimen enter the exit wave $\psi = P\,O$ symmetrically, so swapping their roles in the specimen update gives the probe update.

$$ O \leftarrow O + \frac{P^*}{\lvert P \rvert_{\max}^2} (\psi' - \psi), \qquad P \leftarrow P + \frac{O^*}{\lvert O \rvert_{\max}^2} (\psi' - \psi) $$

The weighting is simpler than PIE's: multiply by the conjugate and divide by the squared maximum, with no division by the probe itself. This form has a clean interpretation. The gradient of the error at one position, $\lVert\, \lvert \mathcal{F}(P\,O) \rvert - \sqrt{I_j} \,\rVert^2$, with respect to the specimen is $2P^*(\psi - \psi')$. The ePIE specimen update is exactly a step of $1/(2\lvert P \rvert_{\max}^2)$ against that gradient. Checked against a finite-difference gradient, the gradient agrees to about $10^{-9}$, and the updates agree to $10^{-16}$. Unlike [gradient descent](/en/posts/optimization-gradient-descent/), which reduces the total error at once, it picks positions one at a time in random order and reduces only that position's error. This is called stochastic gradient descent. Post 8 returns to this view.

<img src="/assets/img/posts/ptycho-overlap-pie-epie/en/fig5-epie-probe.png" alt="Amplitude and phase of the initial disc probe, the probe found by ePIE, and the true probe" width="700">
_Fig 5. Starting from a disc of roughly the right size, ePIE recovers the probe's fringes and phase_

ePIE started from the same disc follows the green curve in Fig. 3. It is not smooth: the error rises from 0.42 at 10 iterations to 0.49 at 20, then falls to 0.19 at 50, 0.076 at 100, and 0.048 at 200. The probe error reaches 0.065 by 200 iterations. As Fig. 5 shows, the probe ePIE finds matches the true probe down to its Fresnel fringes and concentric phase, although all it was given was a disc of roughly the right size. Maiden and Rodenburg also started from a disc 2 pixels too large and held the probe fixed for the first iteration, and so does this post.

The price is overlap. The ePIE column of the table and Fig. 4 show that recovering the probe as well needs much more overlap than PIE. With 72 % overlap or more (90 % diameter) the error is below 0.05, but at 67 % it jumps to 0.18, and at 56 % it is 0.66, effectively a failure. With more unknowns, the measurements need more grip on each other.

More unknowns also bring new ambiguities. Multiplying the specimen by a constant $c$ and dividing the probe by $c$ leaves the exit wave unchanged. Multiplying the specimen by a phase ramp $e^{i g \cdot r}$ and removing the same ramp from the probe adds only a constant phase per position, leaving the intensities intact. Shifting specimen and probe together does the same. Every error in this post was measured after aligning these three. When they cause real trouble in a reconstruction is the subject of Post 4.

## 6. The Redundancy Overlap Creates

Post 2 ended by saying that a single pattern has no redundancy. Does overlap create it? Counting makes the difference plain. Each of the 144 positions in Fig. 1 gives 64 × 64 = 4096 intensities, about 590,000 measurements in total. The well-illuminated specimen has about 11,000 pixels, or about 23,000 unknowns since they are complex. Measurements outnumber unknowns about 26 to 1, and still 19 to 1 counting the probe as unknown. Each specimen pixel appears in 12 different patterns on average.

Does this redundancy help when the photon budget is the same? The number of photons per illuminated specimen pixel was matched. For the single-pattern case, the central 48 × 48 pixels of the same specimen were reconstructed with Post 2's method (HIO followed by ER, median of six starts).

<img src="/assets/img/posts/ptycho-overlap-pie-epie/en/fig6-dose.png" alt="Error of single-pattern CDI, PIE, and ePIE versus photons per illuminated pixel, with reconstructions" width="900">
_Fig 6. At equal dose, overlapping measurements give a smaller error (2-4 times, for PIE)_

| Photons per pixel | Single-pattern CDI | PIE | ePIE |
|---|---|---|---|
| 125 | 0.36 | 0.16 | 0.15 |
| 1,255 | 0.20 | 0.061 | 0.072 |
| 12,549 | 0.073 | 0.019 | 0.049 |

At equal dose, the PIE error is 2–4 times smaller than the single-pattern error. At 1,255 photons per pixel the single pattern gives 0.20 with visible noise texture, while PIE gives a clean 0.061. Spreading the photons over many patterns makes each pattern dimmer, but the same specimen pixel is measured more than ten times under different conditions, and the noise averages down. Robustness to noise was also what Maiden and Rodenburg emphasized when introducing ePIE.

At high dose ePIE has a larger error than PIE, because 200 iterations are not yet enough to finish recovering the probe. At 12,549 photons per pixel PIE reaches 0.019 while ePIE stays at 0.049. Not knowing the probe also costs iterations.

## Summary

Ptychography moves the probe with overlap and records many diffraction patterns. PIE corrects the exit wave at one position at a time to agree with the data and feeds the correction back into the specimen. Plain division by the probe blows up in dark regions, so a small constant is added as in a Wiener filter and the correction is applied in proportion to the probe's brightness. With a disc illumination and a single position it is identical to Post 2's ER. With many overlapping positions there is no need for a support or for worry about the twin image, and it converges in tens of iterations.

When the probe is unknown, ePIE corrects specimen and probe symmetrically, as stochastic gradient descent on the error at each position. The cost is more overlap: in these computations PIE needed 56 % and ePIE 72 % or more. Overlap gives the measurement redundancy, and for the same photons a reconstruction 2–4 times more accurate than a single pattern.

## Next

Making the probe unknown in ePIE introduced new ambiguities. Scale exchange and phase ramps look harmless, something to align only when measuring error, but on an overly regular scan grid they combine to imprint a grid pattern on the reconstruction. Post 4 looks at when these ambiguities cause real trouble and introduces rPIE, a refinement of ePIE. Post 5 then compares mPIE, which adds momentum, with the difference map (DM), Post 2's algorithm applied to ptychography.

## References

- J. M. Rodenburg, H. M. L. Faulkner, "[A phase retrieval algorithm for shifting illumination](https://doi.org/10.1063/1.1823034)," *Applied Physics Letters* 85(20), 4795-4797 (2004).
- A. M. Maiden, J. M. Rodenburg, "[An improved ptychographical phase retrieval algorithm for diffractive imaging](https://doi.org/10.1016/j.ultramic.2009.05.012)," *Ultramicroscopy* 109, 1256-1262 (2009).
- J. M. Rodenburg, A. M. Maiden, "[Ptychography](https://eprints.whiterose.ac.uk/id/eprint/127795/)," in *Springer Handbook of Microscopy*, P. W. Hawkes, J. C. H. Spence (eds.), Springer, 2019, pp. 819-904, [doi:10.1007/978-3-030-00069-1_17](https://doi.org/10.1007/978-3-030-00069-1_17). Origin of the name and Hoppe's idea (Sec. 2), intuitive derivation of the PIE update (Sec. 3.4), and the example of overlap removing ambiguities (Sec. 3.3).
- W. Hoppe, "[Beugung im inhomogenen Primärstrahlwellenfeld. I. Prinzip einer Phasenmessung von Elektronenbeugungsinterferenzen](https://doi.org/10.1107/S0567739469001045)," *Acta Crystallographica A* 25(4), 495-501 (1969).
- R. Hegerl, W. Hoppe, "[Dynamische Theorie der Kristallstrukturanalyse durch Elektronenbeugung im inhomogenen Primärstrahlwellenfeld](https://doi.org/10.1002/bbpc.19700741112)," *Berichte der Bunsengesellschaft für physikalische Chemie* 74(11), 1148-1154 (1970).

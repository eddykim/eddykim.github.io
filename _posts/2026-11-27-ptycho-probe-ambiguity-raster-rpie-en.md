---
title: "Computational Imaging and Ptychography 4 — Why Does a Regular Scan Leave a Pattern? Probe Retrieval and Its Pitfalls"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-probe-ambiguity-raster-rpie/
page_id: ptycho-probe-ambiguity-raster-rpie
date: 2026-11-27 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, raster-grid-pathology, rpie, fermat-spiral, phase-retrieval]
description: "When the ambiguities of an unknown probe imprint a pattern on the reconstruction, and how rPIE makes ePIE faster."
math: true
---

ePIE in [Post 3](/en/posts/ptycho-overlap-pie-epie/) recovered the probe along with the specimen. The price was a new set of solutions that the measurement cannot tell apart. Multiplying the specimen by a constant and dividing the probe by the same constant left every diffraction pattern unchanged. So did giving the specimen and the probe opposite phase ramps. Post 3 simply aligned these away when measuring the error.

Is that really safe? This post looks at these ambiguities more closely. Most of them are indeed harmless. On a regular scan grid, however, the number of ambiguities becomes infinite, and they are imprinted on the reconstruction as a grid-shaped pattern. After that, the post returns to ePIE, which Post 3 found to need many iterations once the probe is unknown. rPIE changes a single weight in the update and solves the same problem more than ten times faster.

## 1. What the Measurement Cannot Distinguish

Restating the measurement from Post 3: at position $j$ the exit wave is $\psi_j(r) = P(r - R_j)\,O(r)$, and the detector records $\lvert \mathcal{F}\{\psi_j\} \rvert^2$. In what ways can $O$ and $P$ be changed without changing this measurement?

The first is a scale exchange. Multiply $O$ by a complex number $c$ and divide $P$ by $c$; the product is unchanged. The second is a phase ramp. Multiply the specimen by $e^{i g \cdot r}$ and the probe by $e^{-i g \cdot r}$, and the exit wave picks up only a constant phase $e^{i g \cdot R_j}$. The constant differs from one position to the next, but within each pattern it is a constant, so the intensity does not change. The third is to shift the specimen and the probe together in the same direction. The scan positions fix only their relative position, so the measurement is the same. When Maiden and Rodenburg introduced ePIE, they listed the scale exchange and the shift among these as ambiguities of the solution.

<img src="/assets/img/posts/ptycho-probe-ambiguity-raster-rpie/en/fig1-ambiguities.png" alt="Specimen phase and probe amplitude and phase for the original pair and for pairs transformed by a scale exchange, a phase ramp, and a lattice-periodic function" width="850">
_Fig 1. All four give identical diffraction patterns; the last only on a regular grid_

The first three columns of Fig 1 are the original pair and the two transformations above. The specimen and probe phases look clearly different, yet when all 144 diffraction patterns are computed and compared, the largest relative difference is of order $10^{-16}$, which is rounding error.

These three are effectively harmless. Each acts on the whole specimen at once, so it does not destroy any information in the image. They amount to a choice of brightness reference, a tilt of the phase baseline, and a choice of coordinate origin. In a real measurement one takes empty space without specimen as the reference, or aligns them when measuring the error. The fourth column of Fig 1 is a different matter.

## 2. The Trap of the Regular Grid

Suppose the scan positions lie on a square grid of spacing $a$, so every $R_j$ is an integer multiple of $a$. Now pick any complex function $f(r)$ that is periodic with period $a$ in both directions, multiply the specimen by it and divide the probe by it.

$$ \frac{P(r - R_j)}{f(r - R_j)}\, O(r)\, f(r) = P(r - R_j)\, O(r)\, \frac{f(r)}{f(r - R_j)} = \psi_j(r) $$

Since $f$ has period $a$ and $R_j$ is a multiple of $a$, $f(r - R_j) = f(r)$. The exit wave does not change at all. The three ambiguities above were fixed by a few numbers, a constant $c$ or a ramp $g$. Here a whole periodic function is free. It can take any shape within one grid cell, so there are as many degrees of freedom as there are pixels in a cell. The fourth column of Fig 1 is such a pair on a grid of 8-pixel spacing. The specimen has stripes and the probe is distorted, yet the diffraction patterns agree to within rounding error. On a jittered grid or a Fermat spiral, the same transformation changes the intensities by up to 20 %.

This is the raster grid pathology. Thibault and colleagues pointed it out in 2009, and tools such as PtyLab have since taken it as a basic premise that ptychography needs a non-periodic scan to avoid ambiguities. The equation also shows that the problem exists only when the probe is unknown. If the probe is known exactly, $P/f$ is not allowed.

Does a pattern really remain in the reconstruction? The same 144 positions were arranged in three ways and reconstructed without knowledge of the probe (rPIE from the next section, 300 iterations). Each case was run with three different orders of processing the patterns, and the table and figure show the run with the middle error.

<img src="/assets/img/posts/ptycho-probe-ambiguity-raster-rpie/en/fig2-raster-grid.png" alt="Scan positions, reconstructed amplitude, error image, and error spectrum for a regular grid, a jittered grid, and a Fermat spiral" width="900">
_Fig 2. Only on the regular grid does the error form a grid pattern, concentrated at grid frequencies_

| Scan | Specimen error | Probe error | Data residual | Probe known | Error share at grid frequencies |
|---|---|---|---|---|---|
| Regular grid (8-pixel step) | 0.057 | 0.077 | $8.1\times10^{-4}$ | 0.00035 | 34 % |
| Grid jittered by ±1 pixel | 0.00053 | 0.00058 | $9.0\times10^{-4}$ | 0.00038 | 1.0 % |
| Fermat spiral | 0.00062 | 0.00047 | $1.5\times10^{-3}$ | 0.00062 | 1.0 % |

The error image of the regular grid shows a clear dot pattern at 8-pixel spacing. In the error spectrum, the bins at the grid frequency and its harmonics make up only 1 % of all bins, yet they hold 34 % of the error power. For the jittered grid and the Fermat spiral that share is 1 %, the same as the bin fraction, so nothing is concentrated there. The size of the error also differs by a factor of a hundred. When the probe is given, all three scans reach the same error of 0.0004–0.0006. This shows that the pattern comes from the freedom in the probe.

The data residual column is the relative difference between the measured diffraction amplitudes and those recomputed from the reconstructed specimen and probe. The regular grid's residual is at the same level as the other two scans. The reconstruction explains the data as well as the others do, yet the specimen is a hundred times more wrong. Nothing in the measurement reveals that this reconstruction is wrong, and that is exactly what an ambiguity is. Indeed, two runs with different orders stopped at different answers, with errors of 0.045 and 0.057. Since the data cannot choose between them, the starting point and the path of the updates decide the answer.

How much, then, does the grid need to be jittered? The jitter amplitude was set to ±1, ±2 and ±4 pixels, and each amplitude was reconstructed nine times, with three jitter seeds and three pattern orders. Among the runs that succeeded in explaining the data, the median errors are 0.00053, 0.00058 and 0.00068, nearly independent of the amplitude, and the same as the 0.0006–0.0009 of the Fermat spiral over twelve orders. Jittering by just one pixel, an eighth of the step, is enough.

The equation makes this unsurprising. In the derivation above, $f$ had to stay the same when shifted by every $R_j$. On a regular grid every $R_j$ is a multiple of 8 pixels, so any function with an 8-pixel period satisfies the condition. Once each position carries an offset of a pixel or so, the differences between positions generate shifts of one pixel horizontally and vertically. Then $f$ must be unchanged by a one-pixel shift, and on a pixel grid the only such function is a constant. What remains is the scale exchange of Section 1, which is harmless.

The jittered grid did show a different kind of failure, though. In about four of every nine runs the reconstruction stopped near an error of 0.018 or 0.034. These runs have residuals of $3\times10^{-2}$ or more, so they do not explain the data, and running them to 2000 iterations did not get them out. Looking closer, the centre of mass of the probe is off by 1–2 pixels from the true probe. The probe appears to have slid along the shift direction of Section 1, while the specimen, whose region outside the illumination still holds its initial value, could not follow as a whole, leaving the two locked in a misaligned state. This kind of stagnation differs in nature from an ambiguity. The residual is large, so the data does report the failure. In these calculations, rerunning with a different order got out of it. The Fermat spiral never stagnated over twelve orders, but these calculations do not explain why stagnation is frequent on the jittered grid alone.

In practice the Fermat spiral is the common choice. Odstrčil and colleagues and PtyLab use it, citing Huang and colleagues. The points on the spiral have nearly uniform spacing to their neighbours, yet do not repeat in any direction. It removes periodicity while keeping the overlap as even as a grid. Kandel and colleagues showed an example where, with few photons, the pattern from a regular grid spoiled the reconstruction.

## 3. Why Is ePIE Slow?

In Post 3, ePIE needed far more iterations than PIE with a known probe. What makes it slow? In a 2017 paper, Maiden, Johnson and Li answered this question by reading the ePIE specimen update in three ways.

$$ o' = o + \frac{P^*}{\lvert P \rvert_{\max}^2} (\psi' - \psi) $$

First, it can be read as a weighted average. The exit wave corrected to the measurement, divided by the probe, $\psi'/P$, is a new estimate, and it is mixed with the old estimate $o$ in the ratio $w = \lvert P \rvert^2 / \lvert P \rvert_{\max}^2$. Where the probe is bright the new value is trusted; where it is dark the old value is kept. Second, it can be read as gradient descent. As shown in Post 3, one ePIE step moves along the gradient of the per-position error with a step proportional to $1/\lvert P \rvert_{\max}^2$. That step comes from the Lipschitz constant of the error gradient, the same kind of step used by the Landweber iteration of [Post 1](/en/posts/ptycho-computational-imaging-inverse-problem/). It is stable but conservative. Third, it can be read as a regularized cost. Minimising the exit-wave error plus a penalty saying "do not change the specimen much" gives the ePIE update. The strength of that penalty, however, differs from pixel to pixel. It is as if the $\lambda$ of the Tikhonov regularization in Post 1 varied with position.

All three views point to the same problem. The weight $w = \lvert P \rvert^2 / \lvert P \rvert_{\max}^2$ falls off too quickly. Where the probe is half as bright as its maximum, ePIE applies only 25 % of the correction, and where it is 0.2 only 4 %. In effect it corrects little apart from the probe centre.

<img src="/assets/img/posts/ptycho-probe-ambiguity-raster-rpie/en/fig3-update-weights.png" alt="Fraction of the correction applied to the specimen by ePIE, PIE and rPIE as a function of probe brightness" width="650">
_Fig 3. ePIE barely corrects anything but the bright region; rPIE lifts the curve through α_

So Maiden and colleagues proposed rPIE (regularized PIE), which changes only the denominator.

$$ o' = o + \frac{P^*}{(1 - \alpha)\lvert P \rvert^2 + \alpha \lvert P \rvert_{\max}^2} (\psi' - \psi) $$

With $\alpha = 1$ it is ePIE. As $\alpha$ approaches 0 the denominator becomes $\lvert P \rvert^2$, and the update simply divides the exit wave by the probe. In between, it is close to division like the PIE of Post 3, while a floor of $\alpha\lvert P \rvert_{\max}^2$ keeps it from blowing up where the probe vanishes. As Fig 3 shows, with $\alpha = 0.05$ the update applies 87 % of the correction even where the probe is half as bright. The probe update is the same formula with the roles of specimen and probe exchanged, and Maiden and colleagues recommended keeping it at $\beta = 1$, the same form as ePIE. The reason is that the specimen usually has a fairly uniform, high transmission, so the weight in the probe update does not vary much. This post also uses $\alpha = 0.05$ and $\beta = 1$.

<img src="/assets/img/posts/ptycho-probe-ambiguity-raster-rpie/en/fig4-convergence.png" alt="Convergence of ePIE and rPIE for a pinhole probe, a defocused focused probe, and a diffuser probe" width="900">
_Fig 4. rPIE converges more than ten times faster than ePIE; neither solves the diffuser probe_

Fig 4 shows reconstructions from a Fermat spiral scan, repeated with three pattern orders. The line is the median and the band is the range of the three. With the same pinhole probe as Post 3, rPIE reaches an error of 0.014 within 20 iterations, better than the 0.024 that ePIE reaches after 200. At 200 iterations rPIE is at 0.0012, a twentieth of ePIE. rPIE traces almost the same curve whatever the order, while for ePIE the iteration at which a given error is reached varies by tens depending on the order.

The second panel is the case Maiden and colleagues gave as one where ePIE struggles. When a beam focused by a lens is used as the probe, it is hard to know exactly how far the specimen is from focus. The curvature of the probe wavefront is therefore estimated wrongly. Here the true probe is 100 µm out of focus, and the initial guess was a probe exactly in focus. rPIE goes down to 0.0005 in 200 iterations, while ePIE stays between 0.008 and 0.07 depending on the order.

<img src="/assets/img/posts/ptycho-probe-ambiguity-raster-rpie/en/fig5-defocus-probe.png" alt="Amplitude and phase of the defocused true probe, the in-focus initial guess, and the probes found by ePIE and rPIE after 200 iterations" width="800">
_Fig 5. From an in-focus guess, both ePIE and rPIE recover the curved wavefront; rPIE is far faster_

In the phase of Fig 5, the initial guess is in focus and its phase is nearly flat, whereas the true probe is out of focus and its phase is curved into concentric rings. Both algorithms found that curved wavefront within 200 iterations, with probe errors of 0.0003–0.0006 for rPIE and 0.009–0.023 for ePIE.

The third panel is different. The diffuser probe, made from a pinhole with a random phase plate, has a complicated structure. Starting from a single disc as the guess, both ePIE and rPIE sit at an error of 0.99 without moving. Maiden and colleagues also reported that ePIE fails completely with a diffuser probe. In their results rPIE still converged in this case, but in these calculations rPIE failed too. The strength of the diffuser or the initial guess may have differed. Either way, there are cases that changing the weight alone cannot escape.

## Summary

With the probe unknown, solutions appear that the measurement cannot distinguish. The scale exchange, the phase ramp and the shift are questions of a reference applied to the whole image, and are harmless. On a regular scan grid, however, any function with the grid's period can be passed between specimen and probe, and that freedom remains in the reconstruction as a grid-shaped pattern. In these calculations 34 % of the error power was concentrated at grid frequencies. The cure is to remove periodicity from the scan, and jittering by one pixel was enough. What makes this ambiguity dangerous is that the data residual does not reveal it. The occasional stagnation on the jittered grid, by contrast, has a large residual, so the data reports it.

ePIE is slow because its weight barely corrects anything outside the bright part of the probe. rPIE changes the denominator so that pixels of intermediate brightness are corrected substantially too, and it solved the same problem more than ten times faster. It was also stable where ePIE struggles, such as with a defocused probe.

## Next

rPIE could not solve the diffuser probe either. Post 5 starts from there. The momentum that Maiden and colleagues added to rPIE (mPIE) solves this problem, but brings a new risk with it. The scale exchange, called harmless in this post, runs away once it meets momentum. Post 5 then compares how the method that applies the difference map of Post 2 to ptychography (DM) handles all of these cases.

## References

- A. Maiden, D. Johnson, P. Li, "[Further improvements to the ptychographical iterative engine](https://doi.org/10.1364/OPTICA.4.000736)," *Optica* 4(7), 736-745 (2017). rPIE and mPIE, the three views of the ePIE update, and the defocused and diffuser probe examples.
- A. M. Maiden, J. M. Rodenburg, "[An improved ptychographical phase retrieval algorithm for diffractive imaging](https://doi.org/10.1016/j.ultramic.2009.05.012)," *Ultramicroscopy* 109, 1256-1262 (2009).
- P. Thibault, M. Dierolf, O. Bunk, A. Menzel, F. Pfeiffer, "[Probe retrieval in ptychographic coherent diffractive imaging](https://doi.org/10.1016/j.ultramic.2008.12.011)," *Ultramicroscopy* 109(4), 338-343 (2009). The original source of the raster grid pathology, cited by the two works below. The original was not consulted for this post; the lattice-periodic ambiguity was derived independently and checked numerically.
- X. Huang, H. Yan, R. Harder, Y. Hwu, I. K. Robinson, Y. S. Chu, "[Optimization of overlap uniformness for ptychography](https://doi.org/10.1364/OE.22.012634)," *Optics Express* 22(10), 12634 (2014). Cited by Odstrčil and colleagues and by PtyLab as the source of Fermat spiral scanning.
- M. Odstrčil, A. Menzel, M. Guizar-Sicairos, "[Iterative least-squares solver for generalized maximum-likelihood ptychography](https://doi.org/10.1364/OE.26.003108)," *Optics Express* 26(3), 3108 (2018).
- L. Loetgering et al., "[PtyLab.m/py/jl: a cross-platform, open-source inverse modeling toolbox for conventional and Fourier ptychography](https://doi.org/10.1364/OE.485370)," *Optics Express* 31(9), 13763 (2023).
- S. Kandel, S. Maddali, M. Allain, S. O. Hruszkewycz, C. Jacobsen, Y. S. G. Nashed, "[Using automatic differentiation as a general framework for ptychographic reconstruction](https://doi.org/10.1364/OE.27.018653)," *Optics Express* 27(13), 18653 (2019).

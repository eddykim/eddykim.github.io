---
title: "Computational Imaging and Ptychography 6 — What Should We Trust When Photons Are Scarce?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-noise-maximum-likelihood/
page_id: ptycho-noise-maximum-likelihood
date: 2026-12-03 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, maximum-likelihood, poisson-noise, phase-retrieval, photon-counting]
description: "How the cost function and the update scheme change a ptychographic reconstruction under the Poisson noise of photon counting, compared by calculation."
math: true
---

The data in the previous five posts had no noise. A real detector counts photons one by one. When the specimen is easily damaged, as with X-rays or electrons, photons have to be rationed, and the outer pixels of a diffraction pattern often receive only one or two. [Post 1](/en/posts/ptycho-computational-imaging-inverse-problem/) said that the noise model of photon-counting measurements would return in Post 6. This is that post.

There are three questions. What happens to the algorithms so far when photons are scarce? What cost function reflects the noise properly, and how much better does it do? And what noise were the Fourier projections used so far assuming? The answer to the last question turns out to be the centre of this post. All calculations use the pinhole probe and the 144-position Fermat spiral used since Post 3, each repeated with three noise draws.

## 1. How Much Does a Photon Count Fluctuate?

When a detector pixel receives $I$ photons on average, the number actually counted, $n$, differs every time. This fluctuation follows a Poisson distribution.

$$ p(n \mid I) = \frac{I^{\,n}\, e^{-I}}{n!} $$

Here $I$ is the expected value set by the specimen and probe, the intensity written so far as $\lvert \mathcal{F}\{\psi_j\} \rvert^2$, and $n$ is the measured photon count. The variance of this distribution equals its mean. Where 100 photons are expected, the count varies by about ±10; where 1 is expected, 0 or 2 are common. The relative error is $1/\sqrt{I}$, so brighter pixels are more trustworthy.

<img src="/assets/img/posts/ptycho-noise-maximum-likelihood/en/fig1-photon-patterns.png" alt="The same diffraction pattern drawn with 10 to the 4th, 6th and 8th photons per pattern" width="900">
_Fig 1. The same diffraction pattern drawn at three doses (log colour). At low dose most outer pixels hold 0 or 1 photon_

Diffraction patterns span a very wide range of brightness. Fig 1 shows the pattern at the central position drawn at three doses. In this post the dose is defined as the number of photons in the brightest single pattern. At $10^8$ photons the centre receives close to a million, and only 0.4 % of pixels hold 0–1 photons. At $10^6$ that fraction is 25 %; at $10^4$ it is 94 %. The fine structure of the specimen lives at high spatial frequencies, in the outer part of the pattern, and that information arrives as one or two photons.

## 2. Which Cost Function Matches Photon Noise?

With noise, no specimen fits the measurement exactly. What should be minimised then? Maximum likelihood (ML) picks the specimen and probe that make the measured data most probable. Taking the logarithm of the Poisson distribution and changing its sign gives the cost to minimise.

$$ L_{\text{Poisson}} = \sum_{j,q} \left( I_{jq} - n_{jq} \log I_{jq} \right) $$

Here $j$ is the position and $q$ the detector pixel, and constants depending only on $n$ are dropped. Two other costs are compared with it. One is intensity least squares, the sum of squared intensity differences.

$$ L_{\text{intensity}} = \frac{1}{2}\sum_{j,q} \left( I_{jq} - n_{jq} \right)^2 $$

It is the obvious choice, but it ignores the fact that the variance equals the mean. It weights the large fluctuations of bright pixels and the small ones of dark pixels equally. The other is amplitude least squares, which compares square roots of the intensities.

$$ L_{\text{amplitude}} = \sum_{j,q} \left( \sqrt{I_{jq}} - \sqrt{n_{jq}} \right)^2 $$

Which of the three suits the noise? In a 2012 paper, Thibault and Guizar-Sicairos showed that expanding the Poisson cost around $\sqrt{n}$ gives the amplitude cost. The reason is that the square root evens out the variance. The variance of the count $n$ grows with its mean, but the variance of $\sqrt{n}$ settles at about 1/4 regardless of the mean. This is called variance stabilization.

<img src="/assets/img/posts/ptycho-noise-maximum-likelihood/en/fig2-variance-stabilization.png" alt="Variance of n and of the square root of n under a Poisson distribution, as a function of the mean photon count" width="650">
_Fig 2. The variance of the photon count n equals its mean, but the variance of √n settles near 1/4 once the mean exceeds a few_

Fig 2 checks this by Monte Carlo. The variance of $\sqrt{n}$ is 0.253 at a mean of 30 and 0.251 at 300, close to 1/4. At a mean of 5 it is 0.285 and at 1 it is 0.40, above 1/4, and at even lower means it drops below. The amplitude cost is a Gaussian assumption that every pixel fluctuates by the same amount, and for pixels with more than a few photons that assumption is nearly the same as Poisson. It breaks down for pixels holding 0–1 photons.

One point deserves attention here. The Fourier projection used since Post 3 replaced only the amplitude with the measured value and kept the phase. Differentiating the amplitude cost with respect to the conjugate exit wave gives $\left(1 - \sqrt{n}/\lvert\Psi\rvert\right)\Psi$, and the correction the Fourier projection adds to the exit wave is exactly this with a minus sign. Computed both ways, they agree to $10^{-15}$. Every algorithm in Posts 1–5 has been walking down the gradient of the amplitude cost, and so has been using a noise model close to Poisson from the start. Thibault and colleagues noted that this connection may partly explain why phase-retrieval algorithms are fairly robust to noise.

## 3. How Much Does the Cost Function Change the Reconstruction of the Same Data?

Does the theory hold up in calculation? To isolate the effect of noise, the probe was given. The three costs were minimised the same way, by gradient descent: the gradient is divided by the illumination each pixel received, $\sum \lvert P \rvert^2$, to set the direction, and a line search picks a step that actually lowers the cost. rPIE from Post 3 and the difference map (DM) from Post 5 are included alongside.

<img src="/assets/img/posts/ptycho-noise-maximum-likelihood/en/fig3-error-vs-dose.png" alt="Final error versus dose with the probe known, for rPIE, rPIE with a decaying step, DM, amplitude ML, Poisson ML and intensity least squares" width="750">
_Fig 3. Final error versus dose with the probe known (300 iterations, median and range of 3 noise draws). Amplitude and Poisson ML are lowest_

| Photons in the brightest pattern | rPIE | DM | ML (amplitude) | ML (Poisson) | Intensity least squares |
|---|---|---|---|---|---|
| $10^4$ | 0.22 | 0.71 | 0.143 | 0.148 | 0.46 |
| $10^5$ | 0.098 | 0.17 | 0.056 | 0.052 | 0.41 |
| $10^6$ | 0.031 | 0.085 | 0.017 | 0.0165 | 0.44 |
| $10^7$ | 0.0094 | 0.059 | 0.0053 | 0.0052 | 0.42 |

First, above $10^5$ photons the ML error falls by about $\sqrt{10} \approx 3.2$ for every tenfold increase in dose (2.6 between $10^4$ and $10^5$). The error is inversely proportional to the square root of the photon count, as expected for a measurement whose relative noise is $1/\sqrt{N}$. Reaching a given error takes photons in proportion to the square of the improvement: halving the error costs four times the dose.

Amplitude ML and Poisson ML are practically the same at every dose. The largest difference is 7 %, and even at $10^4$, where 94 % of pixels hold 0–1 photons, they give 0.143 and 0.148. The variance stabilization of Section 2 holds well enough for this problem. Intensity least squares, by contrast, stays near 0.4 however high the dose.

Why is intensity least squares so poor? The cause can be split in two. Starting from the true specimen, each method stays near the minimum of its own cost, and that error measures the quality of the estimator itself. At $10^6$ photons, amplitude ML stays at 0.0172, Poisson ML at 0.0165, and intensity least squares at 0.038. Its minimum alone is more than twice as bad, because it trusts the large fluctuations of bright pixels as much as the dark ones. Starting from a blank specimen, it is still at 0.31 after 1000 iterations. The intensity enters its gradient once more, so the bright central pixels dominate the step size and the structure carried by the dark pixels barely moves. The large wrong patches in the rightmost reconstruction of Fig 4 are the trace of this.

<img src="/assets/img/posts/ptycho-noise-maximum-likelihood/en/fig4-reconstructions.png" alt="True specimen phase and phases reconstructed by rPIE, Poisson ML and intensity least squares at 10 to the 5th photons per pattern" width="950">
_Fig 4. Reconstructed phase at $10^5$ photons per pattern with the probe known. Same data; only the cost and update differ, and only intensity least squares gets large patches wrong_

## 4. Same Cost, So Why Does rPIE Lag?

The table has one odd feature. As Section 2 showed, rPIE follows the amplitude cost, yet at every dose it is 1.5–1.8 times worse than amplitude ML. The cost is the same, so why is the answer different?

The update scheme is different. ML looks at all 144 positions at once and only moves in directions that lower the total cost. rPIE goes through the positions one at a time and corrects the specimen to fit that position's pattern almost completely. Without noise this is fine, because every pattern points to the same answer. With noise each pattern points to a slightly different answer, and rPIE is pulled toward whichever pattern it just saw. Near the minimum it keeps fluctuating by the noise of the last few patterns. [Post 3](/en/posts/ptycho-overlap-pie-epie/) read ePIE as stochastic gradient descent, and stochastic gradient descent with a fixed step wanders within a range proportional to the step size near the minimum. Starting from the true specimen still lands at 0.031, so the gap is not a matter of incomplete convergence.

Then reducing the step should help. With a step of 1 for the first 30 iterations and $1/(1 + (k - 30)/10)$ afterwards, the error at $10^6$ photons fell from 0.031 to 0.022. It did not reach ML's 0.017 within 300 iterations. DM did worse. With noise its projections bounce between two sets that no longer meet, and the result varied widely between noise draws. Thibault and colleagues also used DM with a relaxed Fourier constraint on noisy data.

Projections that fit the measurement exactly are tools for noise-free data. With noise, the right move is not to fit each pattern completely but to step gradually in the direction that lowers the total cost.

## 5. What Happens at Low Dose When the Probe Is Unknown?

So far the probe was given. What happens when it is unknown? Following the conclusion of Post 5, rPIE with a small step for the first 10 iterations was run for 50 iterations to get a starting point, then either continued as rPIE for 400 more iterations or refined with ML.

<img src="/assets/img/posts/ptycho-noise-maximum-likelihood/en/fig5-low-dose-drift.png" alt="Error and amplitude of an unilluminated pixel for rPIE, rPIE without power control, amplitude ML and Poisson ML at 10 to the 5th photons per pattern with the probe unknown" width="950">
_Fig 5. $10^5$ photons per pattern, probe unknown. Continued rPIE drifts along the scale direction and collapses; ML refinement from iteration 50 drops at once and stays down_

At $10^5$ photons, rPIE stays near 0.11 until about iteration 200 and then collapses, reaching 0.95–0.97 in all three runs by iteration 450. Why, with no momentum? The right side of Fig 5 gives the clue. The amplitude of a pixel that is never illuminated grows by more than $10^5$ times. This is the same picture as in Post 5. The factors that power control applies each time it pulls the probe back have piled up on the specimen side, which means the specimen and probe drifted along the scale ambiguity. Without power control, the probe power grows to $6\times10^9$–$2\times10^{11}$ times its true value and the reconstruction collapses the same way. This did not happen with the noise-free rPIE of Post 5.

One way to understand this is as follows. rPIE's per-position updates are not steps that lower a single cost. With noise the patterns point to different answers, and the combined effect of their updates pushes slowly but steadily along the scale direction, where the cost is flat. Because the cost does not change in that direction, nothing pushes back. At $10^6$ photons the same drift was slow enough not to cause trouble within 450 iterations.

ML has no such drift. The line search accepts only steps that actually lower the total cost, and the gradient along the scale direction is zero, so it does not move that way. Refining the rPIE result at iteration 50 with ML brings the error to 0.054–0.056 at $10^5$ photons and 0.017–0.020 at $10^6$, where it stays. These are almost the same as with the probe given in Section 3. With the right noise model, not knowing the probe costs essentially nothing in this problem. Thibault and colleagues proposed ML as a refinement step for a DM result, which is the same order.

## 6. Should the Photons Be Spread Over Many Positions or Concentrated in Few?

Given a fixed dose, is it better to place positions densely and use few photons per pattern, or to use fewer positions with more photons each? Rodenburg and Maiden note in their handbook chapter that what matters is not the noise in any single pattern but the total number of photons the reconstruction gathers.

<img src="/assets/img/posts/ptycho-noise-maximum-likelihood/en/fig6-dose-fractionation.png" alt="Specimen error as the number of positions changes from 47 to 369 with the total photon count fixed" width="650">
_Fig 6. Total photons fixed at $1.44\times10^8$ while the number of positions changes. The total, not the number of positions, sets the error_

The spacing of a Fermat spiral covering the same area was varied to use 47 to 369 positions, with the total photon count fixed at $1.44 \times 10^8$. While the number of positions changes by nearly a factor of eight, the amplitude-ML error stays at 0.016–0.020. It was slightly lower with fewer positions, but small compared with the change from varying the dose. The total photon count sets the error. The number of positions can be chosen for overlap, scan time, and the algorithmic convergence covered up to Post 5.

## Summary

The noise of photon counting follows a Poisson distribution, whose variance equals its mean. The Poisson cost, which uses this directly, and the amplitude cost, which evens out the variance with a square root, gave almost the same answer in these calculations. The Fourier projection used since Post 3 is one gradient step on the amplitude cost, so the algorithms so far were already using a sensible noise model. Intensity least squares, which fits the intensities directly, had an estimator twice as bad and failed to converge. That is the one to avoid.

Even with the same cost, the update scheme changes the result. rPIE, which fits each position's measurement almost completely, fluctuated with the noise and was 1.5–1.8 times worse than ML; at low dose with the probe unknown it drifted along the scale direction and collapsed. DM oscillated with noise. ML refinement, which accepts only steps that lower the total cost, has no such drift and reaches the same error as with the probe known. The error is inversely proportional to the square root of the total photon count, and how those photons are split among positions matters little.

In practice: start quickly with rPIE or DM as before. When photons are scarce, use that result as the starting point and finish with ML on the amplitude or Poisson cost. Do not use a cost that fits the intensities directly.

## Next

Besides noise, real experiments break two more assumptions: scan positions are not exact, and the light is not fully coherent. Post 7 covers correcting the scan positions along with the specimen, and representing partial coherence as a mixture of several probe modes.

## References

- P. Thibault, M. Guizar-Sicairos, "[Maximum-likelihood refinement for coherent diffractive imaging](https://doi.org/10.1088/1367-2630/14/6/063004)," *New Journal of Physics* 14, 063004 (2012). The Poisson, Gaussian and amplitude costs (Eqs. 4, 8–12), the amplitude cost as an approximation to Poisson and a Gaussian with variance 1/4, and ML refinement after DM.
- J. M. Rodenburg, A. M. Maiden, "[Ptychography](https://eprints.whiterose.ac.uk/id/eprint/127795/)," in *Springer Handbook of Microscopy*, P. W. Hawkes, J. C. H. Spence (eds.), Springer, 2019, pp. 819-904, [doi:10.1007/978-3-030-00069-1_17](https://doi.org/10.1007/978-3-030-00069-1_17). Shot noise and total photon count (Sec. 4.7), and the development of algorithms using noise models.
- M. Odstrčil, A. Menzel, M. Guizar-Sicairos, "[Iterative least-squares solver for generalized maximum-likelihood ptychography](https://doi.org/10.1364/OE.26.003108)," *Optics Express* 26(3), 3108 (2018). A least-squares solver using Gaussian and Poisson likelihoods, with step sizes.

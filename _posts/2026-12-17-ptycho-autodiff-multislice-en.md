---
layout: post
title: "Computational Imaging and Ptychography 9 — Is a Forward Model All We Need?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-autodiff-multislice/
page_id: ptycho-autodiff-multislice
date: 2026-12-17 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, automatic-differentiation, multislice, optimization, phase-retrieval]
description: "Posts 3 to 8 recast as one cost minimization: what automatic differentiation does and does not do for us, and a multislice model for thick specimens."
math: true
---

Posts 3 through 8 covered many algorithms: PIE, ePIE, rPIE, momentum, difference map, maximum likelihood, position correction, mixed states and Fourier ptychography. The names differ, but the skeleton is the same. A forward model describes the measurement, a cost measures the mismatch between model and data, and the unknowns are updated to reduce that cost. What changed from post to post was only which quantities were unknown, which cost was used, and how each step was taken.

Yet the most effort always went into the gradient. We had to derive the derivative of the cost with respect to each unknown by hand: one expression for position, another for the modes, each new. Kandel et al. (2019) pointed to exactly this. Analytical derivatives are hard to obtain, which makes new models hard to try. Their remedy is automatic differentiation (AD): write only the forward model as a program, and the computer produces the gradient.

This post asks two questions. Does automatic differentiation really replace the hand derivation? And what does it not replace? Finally, we use the same framework on the thick specimens we have so far avoided. The calculations were done in PyTorch.

## 1. What Do the Methods So Far Look Like in One Table?

| Post | Unknowns | Cost | One step |
|---|---|---|---|
| 3 ePIE | object, probe | amplitude | one position at a time, step size $1/\lvert P \rvert_{\max}^2$ |
| 4 rPIE | object, probe | amplitude | one position at a time, per-pixel step size |
| 5 mPIE, DM | object, probe | amplitude | momentum / combination of two projections |
| 6 ML | object, probe | amplitude, Poisson, intensity | all positions at once, line search |
| 7 | + positions, + modes | amplitude | annealing, gradient / mode projection |
| 8 | spectrum, pupil | amplitude | rPIE with the roles swapped |

Seen this way, every new capability added one unknown, and each time we rewrote the gradient and the step. Automatic differentiation takes over the gradient.

## 2. What Automatic Differentiation Does

Automatic differentiation is neither numerical differentiation (nudge a value and look at the difference) nor symbolic differentiation (simplify an expression). A computation is a chain of elementary operations: addition, multiplication, Fourier transform, absolute value. The derivative of each elementary operation is known. So we record the order of the operations during the forward pass, then start from the cost and walk the record backwards, multiplying the derivatives by the chain rule. It is like leaving a marker at every fork on the way out and reading the markers in reverse on the way back. The cost of obtaining the derivative therefore never exceeds a small multiple of one forward pass.

Is the result the same as the hand derivation? We rewrote the three costs of post 6 in PyTorch and compared the AD gradients with the gradients derived by hand there.

<img src="/assets/img/posts/ptycho-autodiff-multislice/en/fig1-gradient-check.png" alt="Pixel-by-pixel comparison of AD and hand-derived gradients for three costs" width="950">
_Fig 1. Autodiff gradient (half, per the PyTorch convention) against the hand-derived gradient for three costs. All 2000 pixels lie on the diagonal_

For all three costs the maximum relative difference is below $2\times10^{-9}$. One point differs, by exactly a factor of two. PyTorch differentiates a complex variable through its real and imaginary parts, so it returns twice the gradient $\partial L / \partial O^*$ taken with respect to the conjugate variable. The direction is identical; only the step size needs rescaling. In short, writing a few lines of forward model reproduces the derivative worked out by hand in post 6 automatically.

## 3. Does Changing the Optimizer Change the Result?

If the gradients are the same, are the results? We used the problem of post 6 ($10^6$ photons per pattern, unknown probe), started with 50 iterations of rPIE, and continued with three methods: the hand-coded ML of post 6, Adam, the optimizer common in machine learning, and L-BFGS, a quasi-Newton method. All reduce the same Poisson cost.

<img src="/assets/img/posts/ptycho-autodiff-multislice/en/fig2-optimizers.png" alt="Object error versus iteration count and time for hand-coded ML, Adam with three step sizes, and L-BFGS" width="950">
_Fig 2. Convergence from a 50-iteration rPIE result with the probe unknown (median of 3 noise draws; right: one draw). Hand-coded ML is fastest, Adam reaches the same answer later depending on its step, L-BFGS stalled (at 0.026 in the draw on the right)_

| Method | 10 iter. | 50 iter. | 100 iter. | 300 iter. | Time per iteration |
|---|---|---|---|---|---|
| Hand-coded ML (post 6) | 0.0169 | 0.0168 | 0.0168 | 0.0168 | 0.36 s |
| Adam, step 0.03 | 0.029 | 0.017 | 0.0168 | 0.0168 | 0.18 s |
| Adam, step 0.01 | 0.024 | 0.018 | 0.017 | 0.0168 | 0.18 s |
| Adam, step 0.003 | 0.024 | 0.020 | 0.019 | 0.017 | 0.18 s |
| L-BFGS | 0.028 | 0.019 | 0.019 | 0.020 | 0.37 s |

The starting error is 0.036, and the entries are medians over three noise realizations. Hand-coded ML was fastest: it reached the final value in 10 iterations, 3.7 s. Adam takes half the time per iteration, but needed 50–100 iterations (9–18 s) to reach the same answer. With step 0.003 it had not arrived by 300 iterations. L-BFGS stopped at 0.026 in one of the three realizations.

Why is the hand-coded version faster when the gradients agree? The difference lies in how a step is taken. The ML of post 6 divided the gradient by the illumination each pixel received, $\sum \lvert P \rvert^2$. A well-lit pixel has plenty of data and should change little; a poorly lit pixel should change a lot. This scaling comes from the physics. On top of it sat a line search that accepts only steps that actually reduce the cost. Adam does something similar by dividing each pixel's step by the magnitude of its past gradients, but it has to learn this from the data, which takes time, and the overall step size is still chosen by a person. This is the same situation as the momentum value in post 5, which differed from problem to problem. The LSQ-ML of Odstrčil et al. pushes the idea further and computes the best step along each update direction in closed form by least squares. In short, automatic differentiation replaces the gradient, but not a good step. Since the gradients themselves agree (Fig 1), putting the same illumination scaling and line search on the autodiff gradient would likely narrow the gap; this post did not test that combination. What separated in this comparison was not automatic differentiation but the optimizer that chooses the step.

The starting point matters in the same way. On the same data, starting directly from an empty object and a disk probe without rPIE, Adam descended only from 0.78 to 0.62 in 300 iterations, and L-BFGS stayed at 0.69–0.80. rPIE reached 0.036 in 50 iterations. In this calculation, in the early phase of a reconstruction with an unknown probe, the projection methods refined since post 3 were still strong. The order of post 6, starting with rPIE and finishing with ML, holds here too.

## 4. What Happens When Positions Become Unknowns Too?

The real advantage of AD is how easily it adds unknowns. Consider the position correction of post 7 again. Shifting the probe by a sub-pixel amount at each position is a Fourier phase ramp, so it is differentiable with respect to the shift. We put the shift into the probe in the forward model and add the shifts to the list of unknowns. We did not write a single line of the position derivative. As in post 7, the ambiguity of shifting all positions together was removed by resetting the mean of the corrections to zero at every iteration.

<img src="/assets/img/posts/ptycho-autodiff-multislice/en/fig3-autodiff-positions.png" alt="Initial position error and the error remaining after AD-based correction" width="650">
_Fig 3. Position correction by autodiff with the positions as unknowns (3 error draws). Large errors are not fully corrected_

| Initial position error (rms) | Start (rPIE 50 iter., no correction) | AD (Adam 300 iter.) | Post 7 annealing | Post 7 gradient |
|---|---|---|---|---|
| 0.71 px | 0.10–0.14 | 0.0085–0.010 (position 0.014–0.020 px) | 0.0026 | 0.0038 |
| 1.4 px | 0.22–0.23 | 0.013–0.18 (2 of 3 succeeded) | 0.011 | 0.0065 |
| 2.8 px | 0.38–0.43 | 0.19–0.24 (position 0.57–0.77 px left) | 0.010 | 0.036 |

From the rPIE 50-iteration result we updated the object, probe and positions together with Adam for 300 iterations. Removing the positions from the unknowns and running the same Adam leaves the 1.4 px case at 0.19–0.22. Small errors were corrected well: the 0.71 px error shrank to 0.014–0.020 px. The error remained larger than with the dedicated methods of post 7, however, and at 2.8 px, like the gradient method of post 7, it fell outside the capture range and was not fully corrected. Annealing is a random search rather than a derivative, so it is robust to large errors. A derivative sees only its neighborhood.

In short, automatic differentiation nearly removes the cost of adding a new unknown. How well that unknown is found still depends on the shape of the cost and on the step.

## 5. How Do We Model a Thick Specimen?

In every post so far the exit wave was the product of the probe and the object. That is an approximation for thin specimens. In a thick specimen, light spreads inside, and the front and back parts are crossed at different angles. A single product cannot express this. In 2012 Maiden, Humphry and Rodenburg brought the multislice model, used in electron microscopy, to ptychography. The specimen is divided into several thin slices; the wave is multiplied by each slice and then propagated through free space to the next.

$$ \psi_1 = P\,O_1, \qquad \psi_{k+1} = \mathcal{P}_{\Delta z}\{\psi_k\}\; O_{k+1}, \qquad I = \left\lvert \mathcal{F}\{\psi_{\text{last}}\} \right\rvert^2 $$

Here $\mathcal{P}_{\Delta z}$ is angular-spectrum propagation over a distance $\Delta z$. Deriving the gradient by hand requires unwinding this chain slice by slice, which is what the 3PIE of Maiden et al. does. With AD, it takes one added line of propagation in the forward model.

We computed two phase slices (phase 0.8 rad) separated by $\Delta z$, with a focused probe of numerical aperture 0.15. At a wavelength of 633 nm, $\lambda/\mathrm{NA}^2$ (half the distance over which the axial intensity first drops to zero) is about 28 µm. With the probe known, we reconstructed from an empty object with Adam for 500 iterations, using a one-slice model and a two-slice model.

<img src="/assets/img/posts/ptycho-autodiff-multislice/en/fig4-multislice-separation.png" alt="Data residual and object error of the one- and two-slice models versus the spacing of two phase slices" width="950">
_Fig 4. Data residual (left) and specimen error (right) versus the separation between two phase slices. Close slices mix_

| Spacing | One-slice model: residual / error against the product $O_1O_2$ | Two-slice model: residual / per-slice error |
|---|---|---|
| 1 µm | $9.6\times10^{-3}$ / 0.033 | $1.7\times10^{-3}$ / 0.078 |
| 3 µm | $2.5\times10^{-2}$ / 0.063 | $1.1\times10^{-3}$ / 0.035 |
| 10 µm | $5.2\times10^{-2}$ / 0.075 | $2.4\times10^{-4}$ / 0.002 |
| 30 µm | $9.2\times10^{-2}$ / 0.13 | $6\times10^{-6}$ / $10^{-6}$ |
| 80 µm | 0.12 / 0.25 | $1\times10^{-5}$ / $2\times10^{-6}$ |

The one-slice model explains the data less and less as the spacing grows. The residual rises from $10^{-2}$ at 1 µm to 0.12 at 80 µm. The product approximation breaks down. The two-slice model recovers the two slices separately once the spacing is 10 µm or more. At 30 µm and beyond, the error is of order $10^{-6}$.

<img src="/assets/img/posts/ptycho-autodiff-multislice/en/fig5-multislice-images.png" alt="True first and second slices, the two slices reconstructed by the two-slice model at 30 µm and 1 µm, and the one-slice result at 30 µm" width="950">
_Fig 5. With enough separation the two-slice model recovers each slice (top); with close slices it fits the data but the slices mix (bottom third)_

The situation differs when the spacing is small. At 1 µm the two-slice model fits the data well (residual $1.7\times10^{-3}$). The per-slice error, however, is 0.078, and as in the third lower panel of Fig 5 the pattern of the second slice shows through in the first. This is the ambiguity pointed out by Maiden et al. When two slices are close enough for the product approximation to hold, replacing them by $(f\,O_1,\ O_2/f)$ for any function $f$ leaves the measurement nearly unchanged. It has the same form as the periodic function that a regular grid let the object and probe trade in [post 4](/posts/ptycho-probe-ambiguity-raster-rpie/).

In short, the power of the multislice model to separate slices, its depth resolution, is set by the probe. The larger the numerical aperture, the more the spreading of light differs even over short distances, and the better two slices can be told apart. In this calculation the two slices separated clearly at 10 µm, about a third of $\lambda/\mathrm{NA}^2$; among the spacings tested, the boundary lay between 3 µm and 10 µm.

## 6. What Does This Framework Cost, and Which Tools Implement It?

Automatic differentiation has a cost. Walking backwards requires remembering every intermediate value of the forward pass, which takes a lot of memory. With real data of tens of thousands of positions and several slices, this becomes the limit. Practical tools therefore process positions in batches and use GPUs. The calculations here ran on a laptop CPU, with 0.18 s per Adam iteration.

Tools in this direction have indeed taken hold. PtyRAD by Lee et al. optimizes positions, mixed modes, specimen thickness and local tilt together on top of PyTorch, and is reported to be up to 24 times faster than earlier packages (PtychoShelves, py4DSTEM) on GPU for experimental data. The features we attached one by one by hand up to post 7 become one line each in the list of unknowns. The electron ptychography of Chen et al. (2021), which pushed resolution to the limit set by atomic lattice vibrations, was likewise a reconstruction that combined mixed states with a multislice model.

## Summary

The algorithms of posts 3 through 8 were all ways of reducing one cost. Automatic differentiation takes over the gradient. It matched the hand-derived gradient exactly, and it turned adding a new unknown or a new model, such as positions or slices, into one or two lines of forward model.

What it does not take over was equally clear. In this problem the hand-coded version, which divides the step by the illumination, was 2–5 times faster than Adam, and Adam's step still had to be chosen by a person. From an empty start with an unknown probe, both Adam and L-BFGS lagged far behind rPIE. Position correction was less precise than the dedicated methods and could not fix large errors. The practical order is therefore unchanged: start quickly with a projection method and finish with a differentiable model. Automatic differentiation nearly removes the cost of changing the model in that finishing stage.

Thick specimens were handled with the multislice model. The product approximation explained the data less as the spacing between the two slices grew, and the two-slice model recovered the two slices separately when their spacing was wide enough for the probe to tell them apart. When they were too close, it fit the data while the two slices mixed.

## Next Post

Post 10 closes this series. So far imaging has meant looking at a specimen; the same technology can also measure it. We move to imaging reflectometry, which uses extreme ultraviolet (EUV) reflection ptychography to measure the thickness and composition of thin films and patterns. It is where the reflection measurements of the ellipsometer series meet the phase retrieval of this one.

## References

- S. Kandel, S. Maddali, M. Allain, S. O. Hruszkewycz, C. Jacobsen, Y. S. G. Nashed, "[Using automatic differentiation as a general framework for ptychographic reconstruction](https://doi.org/10.1364/OE.27.018653)," *Optics Express* 27(13), 18653 (2019). The problem that hand-derived derivatives block model exploration; automatic differentiation with Adam.
- A. M. Maiden, M. J. Humphry, J. M. Rodenburg, "[Ptychographic transmission microscopy in three dimensions using a multi-slice approach](https://doi.org/10.1364/JOSAA.29.001606)," *Journal of the Optical Society of America A* 29(8), 1606 (2012). Multislice ptychography (3PIE), per-slice constants and the ambiguity between close slices.
- C.-H. Lee, S. E. Zeltmann, D. Yoon, D. Ma, D. A. Muller, "[PtyRAD: A high-performance and flexible ptychographic reconstruction framework with automatic differentiation](https://doi.org/10.1093/mam/ozaf070)," *Microscopy and Microanalysis* 31(4), ozaf070 (2025). An AD tool that jointly optimizes positions, mixed modes, thickness and tilt.
- M. Odstrčil, A. Menzel, M. Guizar-Sicairos, "[Iterative least-squares solver for generalized maximum-likelihood ptychography](https://doi.org/10.1364/OE.26.003108)," *Optics Express* 26(3), 3108 (2018). LSQ-ML, which sets the step along each update direction by least squares.
- Z. Chen, Y. Jiang, Y.-T. Shao, M. E. Holtz, M. Odstrčil, M. Guizar-Sicairos, I. Hanke, S. Ganschow, D. G. Schlom, D. A. Muller, "[Electron ptychography achieves atomic-resolution limits set by lattice vibrations](https://doi.org/10.1126/science.abg2533)," *Science* 372, 826-831 (2021). Mixed-state multislice electron ptychography.
- D. P. Kingma, J. Ba, "[Adam: A method for stochastic optimization](https://arxiv.org/abs/1412.6980)," arXiv:1412.6980 (2014).

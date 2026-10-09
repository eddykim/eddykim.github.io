---
title: "Computational Imaging and Ptychography 1 — Can a Blurred Photo Be Undone? Inverse Problems and Regularization"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-computational-imaging-inverse-problem/
page_id: ptycho-computational-imaging-inverse-problem
date: 2026-11-15 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [computational-imaging, inverse-problem, regularization, deconvolution, coded-exposure]
description: "Why inverting a blur collapses, and how regularization and measurement design share the work of preventing it."
math: true
---

If the camera moves during an exposure, the photo smears sideways. Suppose we know the direction and length of that smear exactly. Can we compute our way back to the sharp photo?

It seems it should work. The smear is a computation that followed a fixed rule, and a known rule can be run backwards. In practice it almost always fails: the "restored" photo comes out as a field of noise far rougher than the original. This post traces where that failure comes from and follows two ways of preventing it. One works on the computing side, filling in what the data cannot tell us with assumptions; this is regularization. The other works on the measuring side, designing the measurement so that it is easy to invert in the first place.

Where the two meet is the field called computational imaging. In a conventional camera the lens finishes the image and the sensor merely records it. In computational imaging the optics encode the scene into something no one could recognize, and the image is produced by computation.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig1-concept.png" alt="Flow of conventional imaging compared with computational imaging" width="800">
_Fig 1. In computational imaging the measurement is a code, not an image_

The destination of this series is ptychography. As previewed in [Electron Microscopy Foundations 7](/en/posts/electron-microscopy-stem-analytical/), ptychography measures nothing but diffraction intensities and still recovers the phase of the specimen. Every tool needed to get there comes down to one pattern: write an equation for how the measurement is formed, then solve that equation backwards. This first post starts from the simplest case, where the equation is linear.

## 1. A Measurement Is a Function of the Unknowns

What do we need before we can compute an image back? An equation stating how the measured values were produced from the scene. This is called the forward model.

Reduce the problem to a single line. One row of a photo is a 1D signal $x$ of 300 brightness values. Let the camera drift 52 pixels to the right at constant speed during the exposure. Split the exposure into 52 chips of time; during each chip the sensor sees the scene shifted by one more pixel, and what it records is the sum of all 52. A single point becomes a 52-pixel streak, and for the whole signal this is a convolution with a 52-sample box kernel.

$$ y = A x + n $$

$A$ is a 351 × 300 matrix. Column $j$ is the trace a single bright pixel $j$ leaves on the sensor, and the smear makes the measurement 51 pixels longer than the signal. $n$ is the noise added by the sensor. Convolution, in other words, is a matrix product. If the edges are assumed to wrap around periodically, the matrix becomes circulant and the Fourier transform diagonalizes it. This is the familiar statement that convolution in space is multiplication in frequency, from [Image Processing 2](/en/posts/imgproc-convolution-kernels/) and [3](/en/posts/imgproc-frequency-domain-filtering/), rewritten as linear algebra.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig2-forward-model.png" alt="Original signal, blur kernels of the box and coded shutters, and the measurements from each" width="850">
_Fig 2. Same motion, different measurements, depending on how the shutter flutters_

Fig. 2 includes a second shutter. Instead of staying open for all 52 chips, it opens and closes in a fixed sequence. This is the coded exposure of Raskar and colleagues (2006), which opens 26 of the 52 chips. What it changes is the subject of Section 7. Until then we work with the box shutter only.

## 2. Dividing Backwards Breaks Down

We know the forward model, so why not invert $A$? The frequency domain makes it look even simpler. The spectrum of the measurement is the spectrum of the signal multiplied by the blur's transfer function $H$, so dividing the measured spectrum by $H$ should undo it. This is the inverse filter.

Goodman lists three defects of this naive approach. First, at frequencies the optics never passed, $H = 0$ and the division is undefined. Second, even inside the passband $H$ may vanish at isolated frequencies; severe defocus and motion blur are the standard examples. Third, every measurement carries noise, and the inverse filter amplifies most strongly exactly those frequencies where the signal-to-noise ratio is worst.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig3-transfer-function.png" alt="Transfer-function magnitude of the box and coded shutters" width="750">
_Fig 3. The box shutter has spectral zeros; the coded shutter has none, at the cost of half the light_

The box shutter is the second case. A 52-sample box has a sinc-shaped transfer function that is exactly zero at 1/52, 2/52, … cycles per pixel. Whatever the signal held at those frequencies is gone from the measurement.

Without noise, then, can we invert? Because a linear (non-wrapping) convolution produces a measurement longer than the signal, the matrix $A$ itself remains invertible even with spectral zeros. Indeed, noiseless least squares recovers the signal to a relative error of $3.6 \times 10^{-14}$. The trouble is the noise.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig4-inverse-blowup.png" alt="Least-squares inversion at three noise levels" width="850">
_Fig 4. Without noise the inversion is perfect; with 1 % noise the error exceeds the signal_

With a noise standard deviation of 0.1 % of the peak signal the relative error is 0.13; at 1 % it jumps to 1.5, meaning the error is larger than the signal itself. Hadamard called a problem well-posed if a solution exists, is unique, and changes only slightly when the data change slightly. Violate any of the three and the problem is ill-posed. Deblurring violates the third, stability.

[Electron Microscopy Foundations 6](/en/posts/electron-microscopy-tem-phase-contrast-ctf/) showed the same structure. At frequencies where the contrast transfer function (CTF) crosses zero, specimen information never reaches the image. Whether the optics are made of glass or of magnetic fields, a frequency erased during image formation cannot be computed back.

## 3. What Amplifies the Noise

How does 1 % noise become a 150 % error? The singular value decomposition (SVD) of the matrix makes it visible.

$$ A = \sum_i s_i \, u_i v_i^{T} $$

The $v_i$ are directions in signal space and the $u_i$ directions in measurement space; the singular value $s_i$ says how strongly the signal component along $v_i$ is carried into the measurement. Least squares splits the measurement along each $u_i$ and divides by $s_i$.

$$ x_{\text{LS}} = \sum_i \frac{u_i \cdot y}{s_i} \, v_i $$

Now follow the noise. White noise favours no direction, so its component along any $u_i$ has the same average size. For Gaussian noise with standard deviation 0.01 that size is $0.01\sqrt{2/\pi} \approx 0.008$; Hansen describes this as the noise floor. The components the true signal leaves in the measurement, $s_i (v_i \cdot x)$, shrink as the index grows. Past the point where the two meet, the measured component is noise and nothing else.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig5-picard.png" alt="Singular values, data coefficients, and recovered versus true solution coefficients" width="850">
_Fig 5. Where the data hit the noise floor, the recovered coefficients depart from the true signal (Picard plot)_

On the left, the noiseless data coefficients (dashed) drop below the noise floor around index 50. The actual data coefficients (solid) are held up by the floor and go no lower. On the right, dividing them by $s_i$ gives the solution coefficients. Beyond index 50 they are noise divided by small singular values, and they exceed the true coefficients by more than a factor of ten. This is a Picard plot, named after the Picard condition: for a solution to be meaningful, the data coefficients must decay faster than the singular values.

The ratio of the largest to the smallest singular value of the box blur, its condition number, is 231. By the standards of numerical analysis that is not a terrible matrix, yet 1 % noise still ruins the solution. Along the weakest direction the noise is amplified nearly 230-fold, and there are hundreds of such directions stacked together.

## 4. Filling the Unknown with an Assumption

What should we do with components that hold nothing but noise? The measurement cannot say what the signal was along those directions, so something else has to fill the gap. The humblest assumption is that the signal is not absurdly large. Attaching that assumption to the objective as a penalty gives Tikhonov regularization.

$$ \hat{x} = \arg\min_x \; \lVert Ax - y \rVert^2 + \lambda^2 \lVert x \rVert^2 $$

The first term asks the solution to explain the measurement, the second to stay small, and $\lambda$ sets the balance. Written with the SVD, its effect is plain: each term of the least-squares solution is multiplied by a filter factor.

$$ \hat{x} = \sum_i f_i \, \frac{u_i \cdot y}{s_i} \, v_i, \qquad f_i = \frac{s_i^2}{s_i^2 + \lambda^2} $$

When $s_i \gg \lambda$, $f_i \approx 1$ and the component is kept; when $s_i \ll \lambda$, $f_i \approx 0$ and it is discarded. Only the components drowned in noise in Fig. 5 are suppressed.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig6-tikhonov.png" alt="Tikhonov reconstructions for several lambda, lambda chosen by the discrepancy principle, and filter factors per singular value" width="900">
_Fig 6. Tikhonov regularization damps only the small-singular-value components_

Too small a $\lambda$ (grey) leaves the noise in; too large (red) smears the signal as well. At the error-minimizing $\lambda = 0.052$ the relative error is 0.198, almost eight times smaller than the 1.5 of the unregularized solution.

The same answer arrives by other routes. With periodic edges, the Tikhonov solution becomes a filter that multiplies each frequency by $\overline{H} / (\lvert H \rvert^2 + \lambda^2)$. This is Goodman's Wiener filter with the noise-to-signal power spectral ratio set to a constant $\lambda^2$. The Wiener filter is derived as the linear filter of least mean-square error; it reduces to the inverse filter where the signal-to-noise ratio is high and attenuates strongly where it is low. Computing both side by side, they agree to within $4 \times 10^{-14}$.

There is also a probabilistic reading. With Gaussian noise and a Gaussian prior of zero mean on the signal, the Tikhonov solution is the maximum a posteriori (MAP) estimate. The penalty term is prior knowledge. Goodfellow and colleagues explain weight decay (L2 regularization) in deep learning the same way. The $\mu I$ added to $J^T J$ in the Levenberg-Marquardt (LM) method of [Optimization 3](/en/posts/optimization-levenberg-marquardt/) has the same structure; there it stabilized each step.

That leaves the choice of $\lambda$. The "error-minimizing $\lambda$" above can only be found when the answer is known, which is never the case with real data. The most intuitive criterion is the discrepancy principle. The measurement contains noise, so a solution that fits it more closely than the noise level is fitting the noise. Choose the $\lambda$ at which the residual $\lVert Ax - y \rVert$ equals the size of the noise, $\sigma \sqrt{m}$, with $m$ the number of measurements. As the middle panel of Fig. 6 shows, this gives $\lambda = 0.054$, practically the optimum, with a relative error of 0.200. Over 15 combinations of random seed and noise level, the worst case was only 6 % above the optimal error.

The price is that $\sigma$ must be known. Measuring a detector's noise was the topic of [Image Processing 1](/en/posts/imgproc-sampling-noise-model/). When it is unknown, a common alternative is to look for the corner of the L-curve, a log-log plot of residual against solution norm. As Hansen points out, though, that shape is pronounced only when the singular values decay gradually to zero. In this example, with a condition number of only 231, the solution norm changes merely from 15.8 to 9.7 and no clear corner forms.

## 5. Stopping Early Is Also Regularization

What if the matrix is millions by millions, too large for an SVD or an inverse? In image restoration this is the usual situation, and iterative methods take over. The simplest one runs [gradient descent](/en/posts/optimization-gradient-descent/) on the least-squares objective.

$$ x_{k+1} = x_k + \omega A^{T}(y - A x_k) $$

This is the Landweber iteration. Run to convergence, it reaches the least-squares solution, the ruined answer of Section 2. The path there is a different story.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig7-landweber.png" alt="Landweber error versus iteration count and reconstructions at three iteration counts" width="850">
_Fig 7. Least-squares gradient descent first improves, then starts reconstructing the noise_

The relative error is 0.40 after 10 iterations, reaches its minimum of 0.18 at 303, then climbs back to 0.76 at 10,000 and 1.17 at 200,000. Improving first and deteriorating later is called semi-convergence.

Writing the $k$-th iterate with the SVD explains it. It equals a solution with filter factors $f_i = 1 - (1 - \omega s_i^2)^k$. Components with large singular values settle within a few steps; those with small singular values fill in very slowly. Stopping early therefore leaves the small-singular-value components mostly empty, a state much like the Tikhonov filter. The filter factor switches from 0 to 1 around $s_i^2 \approx 1/(k\omega)$. Here $\omega \approx 1$, so 303 iterations correspond to $\lambda \approx 1/\sqrt{303} \approx 0.057$, the same order as the optimal 0.052 of Section 4. The iteration count acts as the regularization strength. Hansen classifies this as iterative regularization, and Goodfellow and colleagues show that for a quadratic objective early stopping is equivalent to L2 regularization. The discrepancy principle also tells us when to stop: once the residual has fallen to the noise level.

The Richardson–Lucy iteration, long used in astronomy and fluorescence microscopy, shares this behaviour. It replaces the Gaussian with the Poisson statistics of photon counts. How the noise model changes when photons are counted will come back in Post 6, together with ptychography. The core algorithm of Post 3 is an iteration of this kind too, so it is worth remembering that when to stop matters.

## 6. Change the Assumption, Change the Result

Tikhonov assumed the solution is small. To assume instead that it is smooth, penalize the sum of squared differences between neighbouring pixels, $\lVert Dx \rVert^2$, rather than $\lVert x \rVert^2$. But this signal has several steps. Does a smoothness assumption hold at a step edge?

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig8-tv.png" alt="Derivative Tikhonov and total variation reconstructions compared" width="800">
_Fig 8. For a piecewise-constant signal, a sparse-gradient prior keeps the edges_

Derivative Tikhonov (red) rounds the edges and leaves ripples beside them, with a relative error of 0.166. Replacing the sum of squares with a sum of absolute values, $\lVert Dx \rVert_1$, changes the picture. Squares punish a large difference heavily and spread it across many pixels; absolute values charge the same whether the change happens in one place or many. The penalty therefore favours solutions whose gradient is zero almost everywhere and large only occasionally, which is exactly what a step signal looks like. This is the total variation (TV) regularization proposed by Rudin, Osher and Fatemi in 1992, and it brings the relative error down to 0.100.

The cost is visible too. The smooth sinusoidal stretch has been cut into steps. This staircasing is the characteristic bias of TV: where the assumption fits it wins, and where it does not, the shape of the assumption is stamped onto the result.

The absolute value has no derivative at zero, so plain gradient descent struggles. The common remedy is the alternating direction method of multipliers (ADMM). The gradient $Dx$ is split off as a separate variable $z$, and the requirement that the two agree is enforced with a penalty and a multiplier. Each iteration then takes three lines: solve one linear system for $x$, apply soft thresholding (shrinking small values to zero) for $z$, and accumulate the constraint violation into the multiplier. Solving the same problem with a general-purpose constrained optimizer gives an objective value that agrees to eight digits.

Pushed to its limit, this direction can solve problems with more unknowns than measurements. If the assumption that the signal is mostly zero in some basis is strong enough, a shortage of equations can still pin down a single solution. That is the idea of compressed sensing. DiffuserCam, by Antipa and colleagues, replaces the lens with an inexpensive diffuser and reconstructs 100 million 3D voxels from a single 1.3-megapixel image. As the authors themselves note, however, the effective resolution varies strongly with scene content. The ability to recover more unknowns than measurements comes entirely from the assumption, and it works only as far as the assumption fits the scene.

## 7. Designing the Measurement

So far the measurement matrix $A$ was given, and only the computation changed. What if $A$ itself could be changed? In Section 2 the root of the trouble was the zeros of the transfer function. The fix is a blur without zeros.

That is what the coded exposure of Raskar and colleagues does. Instead of leaving the shutter open for 52 chips, it opens 26 of them according to a binary code found by computer search. The camera moves just as far, but the blur kernel becomes a sparse code rather than a box. As the blue curve in Fig. 3 shows, its transfer function never touches zero; its lowest point is 0.026, against a maximum of 0.5. The condition number drops from 231 to 17.

It is not free. With the shutter open half the time, only half the light arrives. Which wins, the zeros removed or the light lost? The answer turned out to depend on what kind of noise dominates. For each shutter and each method, the best $\lambda$ was chosen and the errors compared.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig9-coded-noise.png" alt="Reconstruction error of box and coded shutters under read noise and shot noise" width="850">
_Fig 9. The coded shutter's gain depends on the noise; under shot noise the light loss costs less_

| Noise condition | Derivative Tikhonov (box → coded) | TV (box → coded) |
|---|---|---|
| Read noise σ = 0.001 | 0.078 → 0.027 | 0.027 → 0.016 |
| Read noise σ = 0.01 | 0.169 → 0.156 | 0.131 → 0.091 |
| Read noise σ = 0.03 | 0.236 → 0.258 | 0.220 → 0.211 |
| Shot noise, 1000 photons per chip | 0.117 → 0.050 | 0.050 → 0.027 |
| Shot noise, 10 photons per chip | 0.217 → 0.212 | 0.199 → 0.156 |

Read noise is constant regardless of signal level, so losing half the light halves the signal-to-noise ratio outright. At low noise the coded shutter wins by about a factor of three; as the noise grows the gain shrinks and, with derivative Tikhonov, reverses. Shot noise, by contrast, scales with the square root of the photon count, so halving the light worsens the signal-to-noise ratio only by $\sqrt{2}$. Here the derivative-Tikhonov error falls by more than half at 1000 photons per chip, and even at an extremely dark 10 photons per chip the coded shutter ties or wins. Which noise dominates depends on the light level, and the crossover can be read from the photon transfer curve of [Image Processing 1](/en/posts/imgproc-sampling-noise-model/).

In a 2D photo the difference is visible to the eye. The same horizontal blur was applied to every row, with shot noise at 1000 photons per chip.

<img src="/assets/img/posts/ptycho-computational-imaging-inverse-problem/en/fig10-2d-motion.png" alt="A photo blurred by the box and coded shutters, restored by least squares and by Tikhonov" width="900">
_Fig 10. Horizontal motion blur at 1000 photons per chip; relative error in parentheses_

Inverting the box-shutter measurement without regularization gives a relative error of 0.364, with the whole photo covered in a noise pattern; even the best Tikhonov leaves a soft smear at 0.092. The coded-shutter measurement inverted without any regularization reaches 0.062, and 0.056 with Tikhonov. A well-designed measurement, borrowing almost nothing from assumptions, beat the box shutter borrowing as much as it could.

This is the point of computational imaging. Measurement design and prior knowledge share the same job. The less information the measurement loses, the less the computation has to lean on assumptions, which is why the optics and the algorithm are designed together rather than separately. Examples abound. Nagahara and colleagues moved the sensor along the optical axis during the exposure so that defocus blur became independent of scene depth; a single kernel then restores the whole photo without knowing the depth. The diffuser in DiffuserCam is likewise a code, chosen so that each point leaves a distinct pattern across the whole sensor.

## Summary

Undoing a blur looks solvable once the forward model is known, but it collapses the moment noise is divided by small singular values. There were two broad defences. On the computing side, the missing components are filled with assumptions: that the solution is small (Tikhonov), smooth, or has sparse gradients (TV), and stopping an iteration early has the same effect. The strength of the assumption is set by the size of the noise. On the measuring side, the optics are designed to lose less information in the first place. How the work is divided between the two is decided by the kind and size of the noise.

## Next

The forward model in this post stayed linear throughout. With $y = Ax$, superposition holds and the SVD explains everything. The story changes when the detector records intensity. A sensor stores not the complex amplitude but its squared magnitude, so the forward model becomes $y = \lvert Ax \rvert^2$. The phase vanishes from the measurement entirely, and the equation is no longer linear. This is where the phase-contrast difficulty of Electron Microscopy Foundations 6 comes from. The next post takes up phase retrieval, recovering phase from intensity alone. It starts from the iterative method of Gerchberg and Saxton (1972), follows how Fienup repaired it, and asks why it still occasionally stalls.

## References

- J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed., McGraw-Hill, 1996, Sec. 8.8 (three defects of the inverse filter; the Wiener filter).
- P. C. Hansen, "[Regularization Tools: A Matlab package for analysis and solution of discrete ill-posed problems](https://www2.imm.dtu.dk/~pcha/Regutools/RTv4manual.pdf)," manual v4.1, Ch. 2 (Hadamard's definition, filter factors, Picard condition, L-curve, discrepancy principle, iterative regularization). Package paper: [Numerical Algorithms 46, 189-194 (2007)](https://doi.org/10.1007/s11075-007-9136-9).
- G. Wetzstein, "[Image Deconvolution](https://stanford.edu/class/ee367/reading/lecture6_notes.pdf)," Stanford EE367 Computational Imaging lecture notes (inverse filter, Wiener filter, TV with ADMM).
- M. P. Deisenroth, A. A. Faisal, C. S. Ong, *[Mathematics for Machine Learning](https://mml-book.github.io/)*, Cambridge University Press, 2020, Sec. 4.5 (singular value decomposition).
- I. Goodfellow, Y. Bengio, A. Courville, *[Deep Learning](https://www.deeplearningbook.org/)*, MIT Press, 2016, Secs. 7.1 and 7.8 (L2 penalty as MAP estimation; equivalence of early stopping and L2).
- L. I. Rudin, S. Osher, E. Fatemi, "[Nonlinear total variation based noise removal algorithms](https://doi.org/10.1016/0167-2789(92)90242-F)," *Physica D* 60, 259-268 (1992).
- S. Boyd, N. Parikh, E. Chu, B. Peleato, J. Eckstein, "[Distributed optimization and statistical learning via the alternating direction method of multipliers](https://doi.org/10.1561/2200000016)," *Foundations and Trends in Machine Learning* 3(1), 1-122 (2011).
- W. H. Richardson, "[Bayesian-based iterative method of image restoration](https://doi.org/10.1364/JOSA.62.000055)," *J. Opt. Soc. Am.* 62(1), 55 (1972); L. B. Lucy, "[An iterative technique for the rectification of observed distributions](https://doi.org/10.1086/111605)," *Astron. J.* 79, 745 (1974).
- R. Raskar, A. Agrawal, J. Tumblin, "[Coded exposure photography: motion deblurring using fluttered shutter](https://doi.org/10.1145/1141911.1141957)," *ACM Trans. Graph.* 25(3), 795-804 (2006). The 52-chip code is taken directly from this paper.
- H. Nagahara, S. Kuthirummal, C. Zhou, S. K. Nayar, "Flexible depth of field photography," *ECCV* 2008.
- N. Antipa, G. Kuo, R. Heckel, B. Mildenhall, E. Bostan, R. Ng, L. Waller, "[DiffuserCam: lensless single-exposure 3D imaging](https://doi.org/10.1364/OPTICA.5.000001)," *Optica* 5(1), 1 (2018).

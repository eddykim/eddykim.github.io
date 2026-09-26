---
title: "Metrology Image Processing 6 — Reading Phase: Analytic Signals, Phase Shifting, Unwrapping"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-phase-hilbert-unwrapping/
page_id: imgproc-phase-hilbert-unwrapping
date: 2026-11-03 20:00:00 +0900
categories: [Computation, Structured Light]
tags: [image-processing, hilbert-transform, phase-shifting, phase-unwrapping, residue, interferometry]
description: "In a fringe pattern the information sits in the phase, not the intensity. Single-frame and multi-frame routes to it, and where recovering it breaks."
math: true
---

The previous five posts were all about reading intensity. Fixing what a pixel value counts ([Post 1](/en/posts/imgproc-sampling-noise-model/)), mixing neighbours ([Post 2](/en/posts/imgproc-convolution-kernels/)), filtering in the frequency domain ([Post 3](/en/posts/imgproc-frequency-domain-filtering/)), locating an edge below the pixel ([Post 4](/en/posts/imgproc-edge-subpixel/)), and extracting one number from many points ([Post 5](/en/posts/imgproc-geometric-fitting/)).

Interferometry and structured-light projection work differently. The height or thickness being measured lives in the phase rather than the intensity, and no detector reports phase directly. This post covers how phase is extracted from intensity, and where recovering it breaks down.

## 1. One frame holds three unknowns

Whether the instrument is an interferometer or a projector, the detector records the same form.

$$ I(x) = a(x) + b(x)\cos\phi(x) $$

Here $a$ is the background, $b$ the modulation amplitude, and $\phi$ the phase of interest. The difficulty is immediate: three unknowns, and the detector supplies only their combination.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig1-where-is-phase.png" alt="The three components of a fringe pattern, and intensity versus shift at one pixel" width="750">
_Fig 1. One frame holds three unknowns and offers one equation_

The left panel is a single frame. The background $a$ sets the centreline of the fringes and $a \pm b$ forms the envelope. Countless combinations of $(a, b, \phi)$ produce the same intensity.

Two routes diverge from here. One solves with a single frame by adding a spatial assumption: if the phase varies slowly on a carrier, the background and the sideband separate in the spectrum. The other captures several frames with the phase deliberately shifted by known amounts. The right panel shows that principle. Holding one pixel and sweeping the shift $\delta$ traces a sinusoid, and the phase of that sinusoid is $\phi$.

## 2. The analytic signal — inventing the missing imaginary part

Take the single-frame route first. If only $\cos\phi$ is available and $\sin\phi$ is missing, the answer is to manufacture $\sin\phi$.

The Hilbert transform does this. It is a linear operation that rotates the phase of every frequency component by $-\pi/2$, turning a cosine into a sine. Attaching the result to the original as an imaginary part gives the analytic signal.

$$ z(x) = I_{ac}(x) + j\,\mathcal{H}[I_{ac}(x)] = b(x)\,e^{j\phi(x)} $$

Its magnitude is the instantaneous amplitude and its argument the instantaneous phase. This is the transform mentioned only in passing in [Post 3](/en/posts/imgproc-frequency-domain-filtering/) as part of the dissertation's signal analysis.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig2-analytic-signal.png" alt="Background-removed signal with its Hilbert transform, instantaneous amplitude, and phase error with and without the background" width="750">
_Fig 2. The analytic signal gives amplitude as magnitude and phase as argument_

In the centre panel the magnitude of the analytic signal tracks the true modulation $b(x)$ closely, with a recovery error of 0.0009 RMS.

Two conditions attach. First, the background $a$ must be removed beforehand. The right panel shows the difference: leaving it in gives a phase error of 0.235 rad RMS, removing it gives 0.061 rad. The argument is otherwise pulled toward the background.

That remaining 0.061 comes almost entirely from the two ends. Excluding 30 pixels at each border gives 0.0034 rad, while those 30 pixels alone give 0.178 rad. The Hilbert transform computes as though the signal continued forever, which is the boundary problem from section 6 of Post 2 appearing again. So is the remedy: shrink the region of interest.

Second, the spectra of $b$ and $\cos\phi$ must not overlap. This is the Bedrosian condition, and in practice it means the modulation must vary far more slowly than the fringes. Images whose contrast changes within a few pixels break the method.

## 3. One frame or several

The second single-frame route uses a tool from Post 3. Keeping only the $+1$ sideband in the spectrum and transforming back yields a complex signal whose argument is the phase. Takeda proposed it in 1982, on the same principle as the notch filter of Post 3.

The multi-frame route is phase shifting. Since the shifts $\delta_k$ are known,

$$ I_k = a + b\cos(\phi+\delta_k) = a + (b\cos\phi)\cos\delta_k - (b\sin\phi)\sin\delta_k $$

is **linear** in the unknowns $(a,\ b\cos\phi,\ b\sin\phi)$. Three frames suffice. That the background $a$ is solved for as one of the unknowns is the crux: it never has to be estimated separately.

```python
# core of phase_extract.py (full code: _code/imgproc-phase-hilbert-unwrapping/)
def phase_shift_n(frames, shifts):
    d = np.asarray(shifts, float)
    A = np.column_stack([np.ones_like(d), np.cos(d), -np.sin(d)])
    sol, *_ = np.linalg.lstsq(A, frames, rcond=None)
    a, c, s = sol
    return np.arctan2(s, c), np.hypot(c, s), a
```

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig3-three-methods.png" alt="Phase error of three methods against noise, and the sideband width trade of the Fourier method" width="750">
_Fig 3. The Fourier method's noise immunity comes from its bandpass, whose width is another trade_

Without noise, phase shifting errs by $6\times10^{-15}$ rad — floating-point precision. Adding no spatial assumption makes it exact in principle. The Hilbert and Fourier methods carry error even without noise: 0.0022 and 0.0037 rad respectively where the shape is steep, which is exactly how far the slow-phase assumption has been violated.

The left panel reverses that order once noise grows. At a noise level of 0.10, Fourier gives 0.127 rad against 0.214 for phase shifting. The Fourier method carries a bandpass that admits only the sideband and discards noise outside it, whereas phase shifting is a per-pixel operation that averages nothing spatially.

This is a difference in what is paid, not in quality. The right panel shows the payment. Narrowing the sideband rejects more noise but smears steep shape, driving the noiseless error to 0.23 rad at a width of $0.25f_0$. Widening it preserves shape but admits noise. At a noise level of 0.05 the optimum sits near $0.4f_0$. The same trade as the smoothing width of Post 2 and the blur width of Post 4.

## 4. When the shift is wrong

Phase shifting assumes the shifts $\delta_k$ are known exactly. In practice piezo non-linearity and calibration error make them drift. Taking the actual shift as $(1+\varepsilon)$ times nominal, the resulting error was measured.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig4-detuning.png" alt="Phase error of four algorithms as the shift error grows" width="700">
_Fig 4. Only the Hariharan 5-step responds quadratically to shift error_

The 3-step and 4-step errors are proportional to $\varepsilon$: tenfold in $\varepsilon$ is tenfold in phase error. The 5-step least-squares solution has the same slope with a smaller coefficient.

Only the algorithm of Schwider and Hariharan behaves differently.

$$ \phi = \operatorname{atan2}\!\big(2(I_4-I_2),\ I_1-2I_3+I_5\big) $$

As $\varepsilon$ grows tenfold from 0.01 to 0.10, its error grows **a hundredfold**, from 0.000044 to 0.004379 rad. The coefficients were chosen so the first-order term cancels, which makes it 40 times more accurate than the 3-step at $\varepsilon = 0.05$.

The practical reading is this. Where the shifting hardware can be calibrated precisely, fewer steps are better, since a shorter capture leaves less room for the specimen to move. Where the calibration cannot be trusted, the answer is not more frames but an algorithm insensitive to detuning. Section 3.4 of the dissertation, which separately corrects amplitude and phase errors arising from retarder non-linearity and depolarization, addresses the same problem by another route.

## 5. Phase breaks at $\pi$

Every phase obtained so far is folded into $(-\pi, \pi]$, because all the formulas above use a four-quadrant arctangent. An ordinary arctangent returns only $(-\pi/2, \pi/2)$ and loses half the phase, so the four-quadrant version, which takes the signs separately, is required. The price is that the result folds every $2\pi$.

Undoing that fold is phase unwrapping. In one dimension it amounts to wrapping the difference between neighbours and accumulating. The method Itoh set out rests on a single assumption: **the true phase difference between adjacent pixels is smaller than $\pi$**.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig5-unwrap-limit.png" alt="Unwrapping error against phase increment per pixel, and the rate of 2-pi jumps against noise" width="750">
_Fig 5. Past $\pi$ the error jumps from $10^{-13}$ to hundreds of radians_

The cliff in the left panel is that assumption. Up to a phase increment of 3.108 rad per pixel the unwrapping error is $10^{-13}$ rad, which is to say exact. At 3.158 rad it leaps to 427 rad. The boundary is exactly $\pi = 3.1416$. Once broken, the result stays off by a multiple of $2\pi$ for the rest of the line.

The structure matches the Nyquist limit of section 6 of Post 1. There, patterns finer than half a period per pixel folded down in frequency; here, phase differences beyond $\pi$ per pixel cannot be undone. They are in fact the same condition. Packing fringes more densely raises phase sensitivity but eventually makes unwrapping impossible.

The right panel is the noise limit. Even with the phase difference comfortably below $\pi$, noise eats the margin. Up to a noise level of 0.08 not one of 200 trials failed; at 0.14, 52% failed; at 0.20, all of them did. The steepness of that transition matters. A system that was working moments ago fails entirely after a small increase in noise.

## 6. In two dimensions the path matters

One dimension offers a single integration path, so there is nothing to choose. Two dimensions offer many paths between two points, and those paths can disagree.

The test is local. Walking a $2\times2$ loop and summing the wrapped differences gives a non-zero result exactly when the unwrapped value inside that loop depends on the path. Such a point is called a residue.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig6-residues-2d.png" alt="Wrapped phase, residues from undersampling, and residues from noise" width="750">
_Fig 6. With residues present, unwrapping depends on the path_

The left panel is the healthy case, with no residues at all and a maximum phase gradient of 0.77 rad per pixel.

The centre panel raises the shape peak until the gradient on its flank reaches 4.59 rad per pixel, past $\pi$. Two hundred and twenty-six residues appear in pairs around the peaks. As with the short arcs of Post 5, this is not a defect of the algorithm but an absence of information in the data. Integrating between the same two points along the two drawn paths gives results differing by exactly $2\times2\pi$.

The right panel keeps the gradient at a healthy 0.77 rad per pixel and raises only the noise. Residues scatter across the whole image, 4658 of them. The two origins call for different responses. Undersampling requires sparser fringes or splitting the shape into separately captured ranges; noise can be reduced by averaging and filtering. The distribution tells them apart at a glance: clustered in pairs around peaks means undersampling, spread evenly means noise.

A phase map with residues can still be forced open with branch cuts or least-norm methods, but those do not create the missing information — they choose where the error will go. In metrology the better response is to treat the appearance of residues as a failed measurement and change the conditions.

## Summary — closing the series

A single fringe frame carries three unknowns — background, modulation, phase — against one equation. Solving from one frame requires the spatial assumption that phase varies slowly, which leaves error even without noise while the accompanying bandpass rejects noise. Solving from several frames is exact without that assumption, but requires known shifts, and when those are wrong the choice of algorithm sets the order of the resulting error. In the recovery step $\pi$ per pixel is an absolute limit, and in two dimensions residues make the result path-dependent.

One thread runs through all six posts. Every operation in image processing rests on an assumption, and the way that assumption breaks is the bias of the method. The sampling of Post 1, the smoothing width of Post 2, the periodic extension of Post 3, the grid-locked bias of Post 4, the algebraic distance of Post 5, and the spatial assumption and shift calibration of this post were all the same story. Choosing a method in metrology means choosing which assumption to live with.

## What comes next

This series covered tools. What follows is measuring real shape with them.

Structured-light projection and moiré form one branch. Attaching the phase shifting of this post to a projected grating gives shape measurement, and it explains why variants such as high-precision computer-generated moiré profilometry (HCGMP), which cancels background light with complementary gratings, came about. The other branch is defect and contaminant detection trained on normal data alone, where the question is how far a reconstruction error can be trusted as an indicator on metrology images.

## References

- M. Takeda, H. Ina, S. Kobayashi, "[Fourier-transform method of fringe-pattern analysis](https://doi.org/10.1364/JOSA.72.000156)," _J. Opt. Soc. Am._, 72(1), 156-160, 1982.
- P. Hariharan, B. F. Oreb, T. Eiju, "[Digital phase-shifting interferometry: a simple error-compensating phase calculation algorithm](https://doi.org/10.1364/AO.26.002504)," _Applied Optics_, 26(13), 2504-2506, 1987.
- K. Itoh, "[Analysis of the phase unwrapping algorithm](https://doi.org/10.1364/AO.21.002470)," _Applied Optics_, 21(14), 2470, 1982.
- R. M. Goldstein, H. A. Zebker, C. L. Werner, "[Satellite radar interferometry: Two-dimensional phase unwrapping](https://doi.org/10.1029/RS023i004p00713)," _Radio Science_, 23(4), 713-720, 1988 (residues and branch cuts).
- E. Bedrosian, "A product theorem for Hilbert transforms," _Proc. IEEE_, 51(5), 868-869, 1963.
- Y. Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," Ph.D. dissertation, Seoul National University, 2025, secs. 3.2 and 3.4 (Hilbert-transform signal analysis and correction of amplitude and phase errors).

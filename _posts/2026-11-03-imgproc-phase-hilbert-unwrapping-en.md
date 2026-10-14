---
title: "Metrology Image Processing 6 — How Do You Read Height from a Photo of Interference Fringes?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-phase-hilbert-unwrapping/
page_id: imgproc-phase-hilbert-unwrapping
date: 2026-11-03 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, hilbert-transform, phase-shifting, phase-unwrapping, residue, interferometry]
description: "In a fringe pattern the information sits in the phase, not the intensity. Single-frame and multi-frame routes to it, and where recovering it breaks."
math: true
---

The previous five posts were all about reading intensity. [Post 1](/en/posts/imgproc-sampling-noise-model/) fixed what a pixel value counts, [Post 2](/en/posts/imgproc-convolution-kernels/) mixed neighbours, [Post 3](/en/posts/imgproc-frequency-domain-filtering/) filtered in the frequency domain, [Post 4](/en/posts/imgproc-edge-subpixel/) located an edge below the pixel, and [Post 5](/en/posts/imgproc-geometric-fitting/) extracted one number from many points.

Interferometry and structured-light projection work differently. The height or thickness being measured lives in the phase rather than the intensity, and no detector reports phase directly. This post covers how phase is extracted from intensity, and where recovering it breaks down.

## 1. Why can the phase not be read straight off one frame?

A camera records only brightness. What has to be known to get height out of that brightness? Whether the instrument is an interferometer or a projector, the detector records the same form.

$$ I(x) = a(x) + b(x)\cos\phi(x) $$

Here $a$ is the background, $b$ the modulation amplitude, and $\phi$ the phase of interest. The difficulty is immediate: three unknowns, and the detector supplies only their combination.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig1-where-is-phase.png" alt="The three components of a fringe pattern, and intensity versus shift at one pixel" width="750">
_Fig 1. Pixels 220 to 380 of a 512-pixel fringe line (carrier 0.45 rad per pixel, period about 14 pixels). Left: the fringe (grey), the background $a$ (blue dashed) and the envelope $a \pm b$ (red dotted); the dash-dot line marks pixel 300, shown on the right. Right: intensity at that pixel as the shift $\delta$ runs from 0 to 2π (green), with the four 4-step samples (red dots). The intensity peaks at $\delta = -\phi$ (dotted, 5.21 rad), so the phase at this pixel is 1.07 rad._

The left panel is a single frame. The background $a$ sets the centreline of the fringes and $a \pm b$ forms the envelope. Countless combinations of $(a, b, \phi)$ produce the same intensity.

Real images are worse. Fig 2 is a photograph of fringes from a Michelson interferometer lit by a sodium lamp. The centre is bright enough to saturate, and towards the right both the background and the fringe contrast fall away. In other words, $a$ and $b$ are not constants but unknowns that change across the frame.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/ext-michelson-rings.jpg" alt="Orange concentric interference rings formed by a Michelson interferometer on a screen" width="600">
_Fig 2. Concentric fringes formed on a screen by a Michelson interferometer in sodium D light. The brightness centreline (background) and the fringe amplitude vary across the screen, and because the fringes are closed rings the phase does not increase in one direction only. (Source: Warren Leywon, [The Fringe formed by the Michelson interferometer](https://commons.wikimedia.org/wiki/File:The_Fringe_formed_by_the_Michelson_interferometer.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), cropped)_

Two routes diverge from here. One solves with a single frame by adding a spatial assumption: if the phase varies slowly on a carrier, the background and the sideband separate in the spectrum. The other captures several frames with the phase deliberately shifted by known amounts. The right panel of Fig 1 shows that principle. Holding one pixel and sweeping the shift $\delta$ traces a sinusoid whose peak falls at $\delta = -\phi$ (dotted), so reading the peak position is reading $\phi$.

## 2. How can a missing sine component be manufactured?

Take the single-frame route first. Once the background is removed the fringe holds only $\cos\phi$. Reading the phase as the argument of a complex number needs $\sin\phi$ as well. If it is missing, manufacture it.

An analogy helps. Seeing $\cos\phi$ alone is like watching a needle turn on a dial through its shadow on one wall. A single shadow cannot tell whether the needle points up or down, that is, whether the angle is $\phi$ or $-\phi$. Add its shadow on the side wall, $\sin\phi$, and the angle can be read directly.

The Hilbert transform casts that second shadow. It is a linear operation that rotates the phase of every frequency component by $-\pi/2$, turning a cosine into a sine. Attaching the result to the original as an imaginary part gives the analytic signal.

$$ z(x) = I_{ac}(x) + j\,\mathcal{H}[I_{ac}(x)] \approx b(x)\,e^{j\phi(x)} $$

Here $I_{ac}$ is the background-removed signal and $j$ the imaginary unit. Its magnitude is the instantaneous amplitude and its argument the instantaneous phase. The relation is approximate for the three reasons below. The author's dissertation (section 3.2) analysed the channelled spectroscopic signal from section 5 of [Post 3](/en/posts/imgproc-frequency-domain-filtering/) with exactly this analytic signal; Post 3 covered only the step of separating out the baseline.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig2-analytic-signal.png" alt="Background-removed signal with its Hilbert transform, instantaneous amplitude, and phase error with and without the background" width="750">
_Fig 3. The analytic signal applied to the same 512-pixel fringe as Fig 1. Left: the background-removed signal (grey) and its Hilbert transform (orange); the background was estimated with a 15-pixel Gaussian smoothing and subtracted. Centre: the magnitude of the analytic signal (green) and the true modulation $b(x)$ (black dashed). Right: the phase error computed with the background left in (red) and removed (blue); once the background is removed, error remains only at the ends._

In the centre panel the magnitude of the analytic signal tracks the true modulation $b(x)$ closely, with a recovery error of 0.0009 root mean square (RMS) away from the outer 40 pixels.

Three conditions attach. First, the background $a$ must be removed beforehand. The right panel shows the difference: leaving it in gives a phase error of 0.235 rad RMS, removing it gives 0.061 rad. The argument is otherwise pulled toward the background.

Yet that remaining 0.061 rad comes almost entirely from the two ends. Excluding 30 pixels at each border gives 0.0034 rad, while those 30 pixels alone give 0.178 rad. The Hilbert transform computes as though the signal continued forever, which is the boundary problem from section 6 of Post 2 appearing again. So is the remedy: shrink the region of interest.

Second, the spectra of $b$ and $\cos\phi$ must not overlap. This is the Bedrosian condition, and in practice it means the modulation must vary far more slowly than the fringes. Images whose modulation changes within a few pixels break the method.

Third, the phase must keep increasing in one direction, which the carrier guarantees. Where the phase decreases, the argument comes out as $-\phi$: the very ambiguity a single shadow could not resolve. Ring fringes like those of Fig 2 reverse direction as they pass through the centre, so without a carrier the method cannot be applied as it stands.

## 3. One frame or several?

The second single-frame route uses a tool from Post 3. Keeping only the $+1$ sideband in the spectrum and transforming back yields a complex signal whose argument is the phase. Takeda proposed it in 1982. Like the notch filter of Post 3 it masks the spectrum, but where the notch erases points of noise, this keeps a single sideband.

The multi-frame route is phase shifting. Since the shifts $\delta_k$ are known,

$$ I_k = a + b\cos(\phi+\delta_k) = a + (b\cos\phi)\cos\delta_k - (b\sin\phi)\sin\delta_k $$

is **linear** in the unknowns $(a,\ b\cos\phi,\ b\sin\phi)$. Three frames suffice. That the background $a$ is solved for as one of the unknowns is the crux: it never has to be estimated separately. Why, then, take four or five frames? More frames let least squares average down the noise, and above all they change how well the result survives a wrong shift. Section 4 takes up the second point.

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
_Fig 4. Phase error of the three methods on the same 512-pixel fringe (30 pixels excluded at each end). Left: RMS error against frame noise for the Hilbert method (red), the Fourier method (orange) and 4-step phase shifting (blue), averaged over 60 runs per noise level; as noise grows the Fourier method becomes the most accurate. Right: as the Fourier sideband width is varied (passband $$f_0 \pm$$ width), the error under noise 0.05 (orange) and the noiseless error where the shape is steep (green)._

Without noise, phase shifting errs by $6\times10^{-15}$ rad, which is floating-point precision. Adding no spatial assumption makes it exact in principle. The Hilbert and Fourier methods carry error even without noise: 0.0022 and 0.0037 rad respectively where the shape is steep, which is exactly how far the slow-phase assumption has been violated.

The left panel reverses that order once noise grows. Noise here is a standard deviation on the same scale as the intensity; with a modulation amplitude of 0.18 to 0.40, a noise of 0.10 is 25 to 56% of the modulation. At that level Fourier gives 0.127 rad against 0.214 for phase shifting. The Fourier method carries a bandpass that admits only the sideband and discards noise outside it, whereas phase shifting is a per-pixel operation that averages nothing spatially.

The difference, then, is in what each method pays, not in quality. The right panel shows the payment. Narrowing the sideband rejects more noise but smears steep shape, driving the noiseless error to 0.23 rad at a width of $0.25f_0$. Widening it preserves shape but admits noise. At a noise level of 0.05 the optimum sits near $0.4f_0$. Widening to $0.9f_0$ raises even the noiseless error again, because the lower edge of the passband reaches $0.1f_0$ and lets background leak in (with a flat background it does not rise). It is the same trade as the smoothing width of Post 2 and the blur width of Post 4.

## 4. How wrong does the phase get when the shift is off?

Phase shifting assumes the shifts $\delta_k$ are known exactly. In practice piezoelectric actuator non-linearity and calibration error make them drift. Taking the actual shift as $(1+\varepsilon)$ times nominal, the resulting error was measured.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig4-detuning.png" alt="Phase error of four algorithms as the shift error grows" width="700">
_Fig 5. Phase error of four algorithms when the shift is $(1+\varepsilon)$ times nominal (noiseless, 512 pixels). The 3-step (120° spacing) and 4-step (90° spacing) follow the first-order slope (grey dotted); only the Hariharan 5-step follows the second-order slope (grey dashed)._

The 3-step and 4-step errors are proportional to $\varepsilon$: tenfold in $\varepsilon$ is tenfold in phase error. The 5-step least-squares solution has the same slope with a smaller coefficient at small $\varepsilon$, but steepens as a second-order term creeps in; from $\varepsilon = 0.01$ to $0.10$ its error grows sixteenfold.

Only the algorithm of Schwider and Hariharan behaves differently. With shifts of $0, \pi/2, \pi, 3\pi/2, 2\pi$ the phase is

$$ \phi = \operatorname{atan2}\!\big(2(I_4-I_2),\ I_1-2I_3+I_5\big) $$

As $\varepsilon$ grows tenfold from 0.01 to 0.10, its error grows **a hundredfold**, from 0.000044 to 0.004379 rad. The coefficients were chosen so the first-order term cancels, which makes it 40 times more accurate than the 3-step at $\varepsilon = 0.05$. Hariharan's paper places the shifts symmetrically from $-\pi$ to $\pi$ and writes the numerator and denominator with opposite signs. The two forms differ by a constant $\pi$, which has no effect on shape.

The practical reading is this. Where the shifting hardware can be calibrated precisely, fewer steps are better, since a shorter capture leaves less room for the specimen to move. Where the calibration cannot be trusted, the answer is not more frames but an algorithm insensitive to shift error (detuning). Section 3.4 of the author's dissertation, which separately corrects amplitude and phase errors arising from retarder non-linearity and depolarization, addresses the same problem by another route.

## 5. When can wrapped phase no longer be undone?

Can the phase obtained so far be turned straight into height? No. Every value is folded into $(-\pi, \pi]$, because all the formulas above use a four-quadrant arctangent. An ordinary arctangent returns only $(-\pi/2, \pi/2)$ and loses half the phase, so the four-quadrant version, which takes the signs separately, is required. The price is that the result folds every $2\pi$.

Undoing that fold is phase unwrapping. Think of the hour hand of a clock. Where it points says nothing about how many turns it has made. Glance at it often, though, and the turns can be counted: compare with the last reading and note whether it passed twelve. If it moved more than half a turn between glances, there is no telling whether it went forwards or backwards.

One-dimensional unwrapping does exactly this: it wraps the difference between neighbours and accumulates. The method Itoh set out rests on a single assumption: **the true phase difference between adjacent pixels is smaller than $\pi$**.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig5-unwrap-limit.png" alt="Unwrapping error against phase increment per pixel, and the rate of 2-pi jumps against noise" width="750">
_Fig 6. Left: maximum error of one-dimensional unwrapping as the carrier is raised (512 pixels, 4-step phase shifting, noiseless). The horizontal axis is the maximum phase increment per pixel; values below $10^{-13}$ rad are drawn at $10^{-13}$. Past $\pi$ per pixel (red dashed) the error leaps to hundreds of radians. Right: the rate of $2\pi$ jumps when noise is added to the base fringe (maximum increment 0.51 rad per pixel), from 200 runs per noise level. Between noise 0.08 and 0.20 the rate goes from 0% to 100%._

The cliff in the left panel is that assumption. Up to a phase increment of 3.108 rad per pixel the unwrapping error stays at the $10^{-13}$ rad level, effectively exact. At 3.158 rad it leaps to 427 rad. The boundary is exactly $\pi = 3.1416$. Once broken, the result stays off by a multiple of $2\pi$ for the rest of the line.

The structure matches the Nyquist limit of section 6 of Post 1. There, patterns finer than half a period per pixel folded down in frequency; here, phase differences beyond $\pi$ per pixel cannot be undone. They are in fact the same condition. Packing fringes more densely raises phase sensitivity but eventually makes unwrapping impossible.

The right panel is the noise limit. The true phase difference of this fringe is at most 0.51 rad per pixel, leaving a margin of 2.63 rad below $\pi$, and noise eats that margin. Up to a noise level of 0.08 not one of 200 trials failed; at 0.10, 2.5% failed; at 0.14, 52%; at 0.20, all of them. The steepness of that transition matters. A small increase in noise turns a measurement that always worked into one that fails outright.

Comparing the margin with the noise on neighbour differences gives a rule of thumb. In 4-step phase shifting, the standard deviation of the phase difference between adjacent pixels is roughly the noise divided by the modulation $b$. At the edges, where $b$ is smallest (0.18), that is 0.44 rad at noise 0.08, a sixth of the margin; 0.56 rad at 0.10, one part in 4.7; and 0.78 rad at 0.14, one part in 3.4. So a margin more than five times the neighbour-difference noise is safe, and at a little over three times half the trials fail (for one line of 512 pixels). With more pixels, one large excursion becomes more likely, and the multiple needed grows.

## 6. Why does the answer depend on the path in two dimensions?

One dimension offers a single integration path, so there is nothing to choose. Two dimensions offer many paths between two points, and those paths can disagree.

The test is local. Walk a $2\times2$ loop and sum the wrapped differences. Having come back to the start, the sum ought to be zero, yet some loops give $2\pi$, like Penrose stairs that climb all the way round and arrive where they began. Such a point is called a residue, and inside that loop the unwrapped value depends on the path.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/en/fig6-residues-2d.png" alt="Wrapped phase, residues from undersampling, and residues from noise" width="750">
_Fig 7. Wrapped phase from 4-step phase shifting on a 200×200-pixel fringe. Left: the healthy case, with no residues. Centre: the shape peak raised tenfold so its flank exceeds $\pi$ per pixel; 226 residues form arcs around the peaks. Right: the same gradient as the left with the noise raised to 0.35; 4658 residues scatter over the whole image. Yellow and cyan dots are residues of sign +1 and −1, and the red and white lines are two paths joining the same two points._

The left panel is the healthy case, with no residues at all and a maximum phase gradient of 0.77 rad per pixel.

The centre panel raises the shape peak until the gradient on its flank reaches 4.59 rad per pixel, past $\pi$. Two hundred and twenty-six residues appear in pairs around the peaks. As with the short arcs of Post 5, this is not a defect of the algorithm but an absence of information in the data. Integrating between the same two points along the two drawn paths gives results differing by exactly $2\times2\pi$.

The right panel keeps the gradient at the healthy 0.77 rad per pixel and raises only the noise. Residues scatter across the whole image, 4658 of them. The two origins call for different responses. Undersampling requires sparser fringes or splitting the shape into separately captured ranges; noise can be reduced by averaging and filtering. The distribution tells them apart at a glance: clustered around peaks means undersampling, spread evenly means noise. Noise residues also pair with an opposite sign right beside them, so what separates the two is location, not pairing.

Real phase maps show the same problem. Fig 8 is an interferogram made from two Sentinel-1 radar acquisitions before and after the 2023 Türkiye earthquakes. One full cycle of colour is $2\pi$ of wrapped phase. Where the fringes are too dense or broken up for the colours to run on, neighbouring phase differences struggle to stay under $\pi$. It was in exactly this kind of satellite radar interferometry that Goldstein and colleagues set out residues and branch cuts.

<img src="/assets/img/posts/imgproc-phase-hilbert-unwrapping/ext-sentinel1-turkiye-interferogram.jpg" alt="Satellite radar interferogram near Gaziantep, Türkiye, with rainbow fringes over terrain" width="760">
_Fig 8. Interferogram of the area near Gaziantep, made from Copernicus Sentinel-1 radar acquisitions before and after the February 2023 Türkiye–Syria earthquakes. The colour bar at lower right runs from −3.14 to 3.14 rad, so one cycle of colour is 2π of wrapped phase. Places where the fringes crowd together or break up are where unwrapping is hard. (Source: European Space Agency, [Türkiye–Syria interferogram ESA24711647](https://commons.wikimedia.org/wiki/File:T%C3%BCrkiye%E2%80%93Syria_interferogram_ESA24711647.jpeg), [CC BY-SA 3.0 IGO](https://creativecommons.org/licenses/by-sa/3.0/igo/), contains modified Copernicus Sentinel data (2023), processed by ESA)_

A phase map with residues can still be forced open with branch cuts or minimum-norm methods. Those do not create the missing information, though; they choose where the error goes, onto the cuts for branch cuts and spread over the whole map for minimum-norm. In metrology the better response is to treat the appearance of residues as a failed measurement and change the conditions.

## Summary

Which method to use comes down to what can be lived with.

| Situation | Method | Reason |
|---|---|---|
| Must capture in one shot, or the specimen moves | Fourier method, starting with a sideband width near $0.4f_0$ | Solves from one frame and its bandpass rejects noise, at the cost of smearing steep shape |
| One-dimensional signal whose modulation varies far more slowly than the fringes | Hilbert analytic signal | Remove the background first and discard the ends |
| Specimen is still and the shifter can be calibrated | 4-step phase shifting | Exact without a spatial assumption |
| The shifter's calibration cannot be trusted | Hariharan 5-step | Shift error enters only at second order (40 times better than the 3-step at $\varepsilon = 0.05$) |
| Residues appear | Do not unwrap; change the conditions | Use sparser fringes or reduce the noise |

A single fringe frame carries three unknowns (background, modulation, phase) against one equation. Solving from one frame requires the spatial assumption that phase varies slowly, which leaves error even without noise; in exchange, the Fourier method's bandpass rejects noise. Solving from several frames is exact without that assumption, but requires known shifts, and when those are wrong the choice of algorithm sets the order of the resulting error. In the recovery step, for spatial unwrapping, which links neighbours, $\pi$ per pixel is a hard limit, and in two dimensions residues make the result path-dependent.

## Closing the series

One thread runs through all six posts. Every operation in image processing rests on an assumption, and the way that assumption breaks is the bias of the method. The sampling of Post 1, the smoothing width of Post 2, the periodic extension of Post 3, the grid-locked bias of Post 4, the algebraic distance of Post 5, and the spatial assumption and shift calibration of this post were all the same story. Choosing a method in metrology means choosing which assumption to live with.

This series covered tools. Two directions lead on from them. One is structured-light projection and moiré. Attaching the phase shifting of this post to a projected grating gives shape measurement, and it explains why variants such as high-precision computer-generated moiré profilometry (HCGMP), which cancels background light with complementary gratings, came about. The other is defect and contaminant detection trained on normal data alone, where the open question is how far a reconstruction error can be trusted as an indicator on metrology images.

## References

- M. Takeda, H. Ina, S. Kobayashi, "[Fourier-transform method of fringe-pattern analysis for computer-based topography and interferometry](https://doi.org/10.1364/JOSA.72.000156)," _J. Opt. Soc. Am._, 72(1), 156-160, 1982.
- P. Hariharan, B. F. Oreb, T. Eiju, "[Digital phase-shifting interferometry: a simple error-compensating phase calculation algorithm](https://doi.org/10.1364/AO.26.002504)," _Applied Optics_, 26(13), 2504-2506, 1987.
- K. Itoh, "[Analysis of the phase unwrapping algorithm](https://doi.org/10.1364/AO.21.002470)," _Applied Optics_, 21(14), 2470, 1982.
- R. M. Goldstein, H. A. Zebker, C. L. Werner, "[Satellite radar interferometry: Two-dimensional phase unwrapping](https://doi.org/10.1029/RS023i004p00713)," _Radio Science_, 23(4), 713-720, 1988 (residues and branch cuts).
- E. Bedrosian, "A product theorem for Hilbert transforms," _Proc. IEEE_, 51(5), 868-869, 1963.
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025.

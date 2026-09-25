---
title: "Metrology Image Processing 2 — Convolution and Kernels: Smoothing, Differentiation, and Borders"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-convolution-kernels/
page_id: imgproc-convolution-kernels
date: 2026-10-06 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, convolution, kernel, gaussian-filter, edge-detection, median-filter]
description: "Smoothing buys lower noise at the cost of blurring something. This post measures that trade directly, as edge position scatter against edge width."
math: true
---

[Post 1](/en/posts/imgproc-sampling-noise-model/) fixed what a single pixel value counts and how much it scatters. Now those values get mixed together. Almost every operation that combines neighbouring pixels reduces to convolution, and smoothing, differentiation and edge detection all come down to choosing one kernel.

Mixing has a price. Averaging neighbours to suppress noise blurs edges; differentiating to sharpen edges brings the noise back. This post does not settle for describing that trade in words. It measures it. The goal is a concrete answer to the question of how wide the smoothing kernel should be.

## 1. Why mix neighbours at all?

Section 5 of Post 1 said that the way to reduce noise is to capture several frames and average them. What if the specimen is moving, or only one exposure is allowed? The only remaining option is to average over space instead of time.

Spatial averaging carries an assumption: that neighbouring pixels are looking at roughly the same thing. On a flat surface the assumption holds well. On an edge it fails completely. The rest of this post is, in effect, an accounting of what that failure costs.

The operation itself is simple. Each output pixel is a weighted sum of input pixels, and the set of weights is the kernel.

$$ g(x) = \sum_{k} w(k)\, f(x-k) $$

The essential point is that the weights do not vary with position. The same kernel applies everywhere in the image, which makes the operation linear and shift-invariant, and therefore a single multiplication in the frequency domain. That is why the arguments for choosing a kernel in sections 3 and 4 all reduce to frequency response.

## 2. How much noise does a kernel remove?

Given a kernel, can the output noise be computed in advance? If the input pixels have independent noise of equal variance $\sigma_{in}^2$, the variance of a weighted sum is proportional to the sum of squared weights.

$$ \sigma_{out}^2 = \sigma_{in}^2 \sum_k w(k)^2 $$

Two conditions separate out here. Preserving the brightness of flat regions requires $\sum w = 1$; reducing noise requires $\sum w^2 < 1$. The second condition is not free. Once negative weights appear, $\sum w^2$ can exceed 1 while $\sum w = 1$ still holds. The common 3×3 sharpening kernel $[0,-1,0;\,-1,5,-1;\,0,-1,0]$ is such a case, with $\sum w^2 = 29$, amplifying noise by a factor of 5.4.

A box kernel of width $k$ has $\sum w^2 = 1/k$, so noise falls as $1/\sqrt{k}$ — the same form as the $1/\sqrt{N}$ of frame averaging in Post 1. The cost differs, though. Frame averaging spends time; spatial averaging spends resolution.

A two-dimensional Gaussian is the one-dimensional kernel applied once along each axis, so its noise ratio is the square of the one-dimensional ratio, namely $\sum w^2$.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig1-noise-gain.png" alt="Noise reduction versus smoothing width, theory compared with measurement" width="700">
_Fig 1. Noise reduction is set by the sum of squared kernel weights_

Filtering an actual white-noise image reproduces the theoretical curve to three decimal places: at $\sigma=2$, theory gives $0.1411$ and measurement $0.1407$. Noise drops to one seventh. Reaching the same reduction by frame averaging would take 50 frames, so spatial smoothing is cheap.

Box kernels matched to the same variance are plotted alongside. At the width-5 box corresponding to $\sigma=1.5$, the box gives $0.1997$ against the Gaussian's $0.1877$. On noise alone the difference is small. The reason to avoid the box lies elsewhere.

## 3. Box or Gaussian?

If several kernels remove comparable amounts of noise, what decides between them? Frequency response. The two pass entirely different bands.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig2-kernel-shape-response.png" alt="Shapes and frequency responses of box and Gaussian kernels" width="750">
_Fig 2. A box kernel cuts high frequencies poorly, leaving sidelobes_

The box response does not settle to zero past its first null but repeats in sidelobes. A width-9 box peaks at $0.227$ in its largest sidelobe and still passes $0.111$ at Nyquist. A $\sigma=2$ Gaussian passes $1.5\times10^{-5}$ there, which is zero for practical purposes.

Section 6 of Post 1 is where this difference bites. Low-pass filtering before sampling is meant to prevent aliasing, but a box lets more than 10% of the content above Nyquist survive and fold straight back down. The alternating sign of the sidelobes is a second problem: some frequency components pass with inverted phase, so a bright stripe in the original can appear as a dark one.

The Gaussian became standard for reasons beyond the absence of sidelobes. A two-dimensional Gaussian separates into two one-dimensional passes, which is cheap, and repeated Gaussian filtering yields another Gaussian whose width is simply $\sqrt{\sigma_1^2+\sigma_2^2}$.

The benefit of separability grows with kernel size.

| Kernel | Direct | Two 1-D passes | Ratio |
|---|---|---|---|
| 5×5 | 25 | 10 | 2.5 |
| 9×9 | 81 | 18 | 4.5 |
| 15×15 | 225 | 30 | 7.5 |
| 25×25 | 625 | 50 | 12.5 |

The figures are multiplications per pixel. Reducing $k \times k$ to $2k$ is often what separates feasible from infeasible in real-time processing.

## 4. Why a derivative kernel behaves differently

Finding an edge means finding where brightness changes abruptly, which means differentiating. How does this operation, the opposite of smoothing, respond to noise?

A derivative kernel has $\sum w = 0$, because it must return zero over flat regions. Its $\sum w^2$, however, is not zero — it erases signal while keeping noise. For the central difference $[-0.5, 0, 0.5]$, $\sqrt{\sum w^2} = 0.707$. It removes only 30% of the noise while eliminating the flat component entirely.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig3-derivative-kernels.png" alt="Frequency responses of derivative kernels compared with the ideal derivative" width="700">
_Fig 3. Derivative response grows with frequency — hence the pairing with smoothing_

The ideal derivative responds in direct proportion to frequency. The central difference tracks that line well at low frequencies, bends away at $f=0.5$, and falls to zero at Nyquist. Patterns near Nyquist are therefore invisible to a central difference, which is worth knowing in its own right.

The real problem is that signal and noise occupy different parts of this response curve. An edge that has passed through optics is blurred, so its gradient information sits at low frequencies — precisely where the derivative response is smallest. Noise, by contrast, is spread evenly across all frequencies and passes freely through the mid-band where the response peaks. This is how differentiation degrades the signal-to-noise ratio.

The remedy is to fuse differentiation with smoothing. A derivative-of-Gaussian (DoG) kernel at $\sigma=2$ brings $\sqrt{\sum w^2}$ down to $0.125$, a factor of 5.7 below the central difference. Fig 3 shows why: the DoG response peaks at $f=0.15$ and suppresses everything above it.

```python
# core of kernels.py (full code: _code/imgproc-convolution-kernels/)
def noise_gain(w):
    """Variance gain against white noise: output variance = input variance x this."""
    return float(np.sum(np.asarray(w, float) ** 2))

CENTRAL_DIFF = np.array([-0.5, 0.0, 0.5])
w_dog = np.gradient(gaussian_kernel(2.0))   # derivative of a Gaussian
```

The familiar Sobel kernel has the same structure: a $[-1, 0, 1]$ derivative component multiplied by a $[1, 2, 1]$ smoothing component along the orthogonal axis, shown as the green curve in Fig 3. Scharr's $[3, 10, 3]$ revises that smoothing component, with the aim of improving the rotational symmetry of the gradient direction rather than the noise.

## 5. How much smoothing to apply

Everything so far was preparation. Stronger smoothing lowers the noise and stabilizes the edge position, while simultaneously blurring the edge and making its position harder to pin down. Two effects in opposite directions imply an optimum. The following experiment locates it.

A synthetic edge with a known position was generated. The convolution of a step with a Gaussian point spread function (PSF) is an error function, so the edge profile is an erf; the blur width was set to $1.2$ pixels and the true position to $63.37$ pixels. The contrast of 1200 DN and the noise of 17.2 DN are carried over from section 5 of Post 1. For each smoothing width, noise was redrawn 400 times, and the position was estimated by interpolating the gradient maximum to the vertex of a parabola.

```python
# core of edge_experiment.py
def estimate_edge(profile, smooth_sigma):
    p = gaussian_filter1d(profile, smooth_sigma)
    g = np.abs(np.convolve(p, CENTRAL_DIFF[::-1], mode="same"))
    i = int(np.argmax(g[2:-2])) + 2
    a, b, c = g[i - 1], g[i], g[i + 1]
    return i + 0.5 * (a - c) / (a - 2.0 * b + c)   # parabola vertex
```

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig4-edge-tradeoff.png" alt="Edge position scatter traded against edge width as smoothing increases" width="700">
_Fig 4. The smoothing trade — position stabilizes, the edge blurs_

Without smoothing the position scatter is $0.094$ pixels. It improves 2.6-fold to $0.037$ pixels at $\sigma=2$, then degrades again, reaching $0.051$ pixels at $\sigma=10$. Over the same range the edge width grows monotonically from $3.2$ to $25.8$ pixels.

That the scatter worsens again is the key result. Early on, the gain from noise reduction dominates. Past a certain width the gradient peak becomes too flat for its vertex to be located reliably. A flat peak means the denominator $a-2b+c$ of the parabolic interpolation approaches zero, so small noise is amplified into large position error.

The practical rule is to set the smoothing width near the width of the optical PSF. With a PSF of $1.2$ pixels, this experiment found an optimum at $2.0$ pixels. Going far beyond the PSF spends resolution for nothing.

One more thing deserves attention: the bias of the estimated position stayed below $0.02$ pixels throughout. Smoothing with a symmetric kernel widens an edge but does not shift it to one side. Where subpixel estimation does acquire a genuine bias is treated separately in Post 4.

## 6. What happens at the borders

A kernel demands neighbouring values for every output pixel. At the very edge of an image, some of those values lie outside the array. What should fill them?

There is no correct answer. Something must be assumed about the outside, and each mode assumes something different. Zero padding assumes darkness beyond the border; replication assumes the final value persists; reflection assumes the signal folds back symmetrically.

A signal that genuinely continues past the array was constructed, and the central 96 pixels were taken as the captured image. The ground truth is the full-length signal smoothed in one piece.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig5-boundary-modes.png" alt="Error comparison of zero, replicate and reflect padding modes" width="750">
_Fig 5. Error left at the borders by each padding mode (right panel rescaled)_

Zero padding produces errors an order of magnitude apart from the others: $191$ DN RMS over the outer five pixels, peaking at $493$ DN. This is the familiar darkening of a bright image's border, and searching for edges on top of it yields a strong spurious edge tracing the frame.

The other two are two orders of magnitude smaller. Replication gives $18.7$ DN RMS and reflection $37.4$ DN. Reflection came out worse here because this scene carries a gradient across the border; folding the signal also inverts the sign of that gradient, creating a kink. Had the scene been flat near the border, the ranking could reverse.

No mode is therefore universally better. In metrology the reliable solution is not to choose the assumption well but to avoid depending on it. Shrinking the region of interest by half the kernel width keeps invented values from ever reaching the result. A little field of view is traded for the complete removal of border artifacts.

## 7. When linear filters do not work

The discussion so far assumed noise of comparable size on every pixel. Dead pixels and cosmic ray hits break that assumption: a single point takes an entirely wrong value. Is a Gaussian still the answer?

Extreme values were planted on 1.02% of the pixels and both filters were compared against the clean original.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig6-median-vs-gaussian.png" alt="Gaussian and median filters applied to an image with spike pixels" width="750">
_Fig 6. A linear filter smears a spike instead of removing it_

The numbers are unambiguous. The spiked image has an RMS error of $14.5$ DN against the original. A Gaussian raises it to $15.3$ DN — **it makes things worse**. A median filter brings it down to $9.2$ DN.

The Gaussian fails because of linearity itself. A linear filter returns a weighted sum of its input, so it cannot destroy the total error an outlier brings; it can only distribute that error among the neighbours. The grey blobs where white dots used to be, in the centre panel of Fig 6, are the result. A single contaminated point has been widened to the kernel width, so the peak error drops while the contaminated area grows.

A median filter sorts and takes the middle value. Being non-linear, it discards outliers outright. The cost is real, though. It is not a convolution, so frequency response says nothing about it, and it erases lines and corners finer than its kernel. Where fine structure is itself the measurand, a median filter cannot be applied casually.

## Summary and what comes next

Choosing a kernel fixes its noise reduction through $\sum w^2$ and its signal preservation through the frequency response. Those two numbers are the whole of smoothing kernel selection. A derivative kernel has $\sum w = 0$ with a response that grows with frequency, so it must be fused with smoothing — which is what the DoG and Sobel families are. Smoothing width has an optimum, found here at roughly 1.7 times the PSF width. Borders and spikes are the two places where the assumptions behind linear convolution break, answered respectively by shrinking the region of interest and by switching to a non-linear filter.

The next post moves to the frequency domain. Having explained kernel behaviour through frequency response, it asks what happens when filters are designed there directly. It covers the removal of periodic noise such as illumination ripple and grating patterns, which cannot be touched in the spatial domain, and goes on to why an ideal filter produces ringing, and what the extrema envelope method developed in the dissertation does about it.

## References

- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 3 (spatial filtering), ch. 10 (edge detection).
- H. Scharr, "Optimal Operators in Digital Image Processing," Ph.D. dissertation, Universität Heidelberg, 2000 (optimizing gradient kernels for rotational symmetry).
- J. Canny, "[A Computational Approach to Edge Detection](https://doi.org/10.1109/TPAMI.1986.4767851)," _IEEE Trans. Pattern Anal. Mach. Intell._, PAMI-8(6), 679-698, 1986 (the case for fusing differentiation with smoothing).
- I. T. Young, L. J. van Vliet, "[Recursive implementation of the Gaussian filter](https://doi.org/10.1016/0165-1684(95)00020-E)," _Signal Processing_, 44(2), 139-151, 1995 (separability and fast implementations).
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the standard test image in Fig 6).

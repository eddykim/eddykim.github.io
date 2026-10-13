---
title: "Metrology Image Processing 2 — Can Blurring an Image Locate Its Edges More Precisely?"
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

The operation itself is simple. Picture a small window of weights sliding across the image one pixel at a time. At each position, the pixels under the window are multiplied by the weights and summed, and that sum becomes the output pixel. Each output pixel is thus a weighted sum of input pixels, and the set of weights is the kernel.

$$ g(x) = \sum_{k} w(k)\, f(x-k) $$

Here $f$ is the input image, $g$ the output image, and $w(k)$ the kernel weight at offset $k$ from the centre. The equation defines convolution, but the kernel arrays in this post are written in the order in which they would be laid over the image, which is strictly the order of correlation. For a symmetric kernel the two coincide; for a derivative kernel only the sign changes. Only gradient magnitudes are used here, so the results are unaffected.

The essential point is that the weights do not vary with position. The same kernel applies wherever the window sits, which makes the operation linear and shift-invariant, and therefore a single multiplication in the frequency domain. That is why the arguments for choosing a kernel in sections 3 and 4 all reduce to frequency response.

## 2. How much noise does a kernel remove?

Given a kernel, can the output noise be computed in advance? If the input pixels have independent noise of equal variance $\sigma_{in}^2$, the variance of a weighted sum is proportional to the sum of squared weights.

$$ \sigma_{out}^2 = \sigma_{in}^2 \sum_k w(k)^2 $$

Two conditions separate out here. Preserving the brightness of flat regions requires $\sum w = 1$; reducing noise requires $\sum w^2 < 1$. The second condition is not free. Once negative weights appear, $\sum w^2$ can exceed 1 while $\sum w = 1$ still holds. The common 3×3 sharpening kernel $[0,-1,0;\,-1,5,-1;\,0,-1,0]$ is such a case, with $\sum w^2 = 29$, amplifying noise by a factor of 5.4.

A box kernel of width $k$ has $\sum w^2 = 1/k$, so noise falls as $1/\sqrt{k}$. This is the same form as the $1/\sqrt{N}$ of frame averaging in Post 1. The cost differs, though. Frame averaging spends time; spatial averaging spends resolution.

A two-dimensional Gaussian is the one-dimensional kernel applied once along each axis. Its ratio of standard deviations is therefore the square of the one-dimensional ratio $\sqrt{\sum w^2}$, which is simply $\sum w^2$ of the one-dimensional kernel. That is why the theory curve in Fig 1 is $\sum w^2$.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig1-noise-gain.png" alt="Noise reduction versus smoothing width, theory compared with measurement" width="700">
_Fig 1. Ratio of standard deviations measured after filtering a 768×768 white-noise image with Gaussians (circles) and with odd-width boxes of matched variance (squares). The solid line is the theoretical Σw² of the one-dimensional kernel; at σ=2 it is 0.141, so the noise drops to one seventh._

Filtering an actual white-noise image reproduces the theoretical curve to three decimal places: at $\sigma=2$, theory gives $0.1411$ and measurement $0.1407$. Noise drops to one seventh. Reaching the same reduction by frame averaging would take 50 frames, so spatial smoothing is cheap.

Box kernels matched to the same variance are plotted alongside. At the width-5 box corresponding to $\sigma=1.5$, the box gives $0.1997$ against the Gaussian's $0.1877$. The widths are rounded to odd integers, so the match is not exact. A width-5 box has a standard deviation of $1.41$ pixels, slightly below $1.5$, which tilts the comparison a little against the box. Either way, on noise alone the difference is small. The reason to avoid the box lies elsewhere.

## 3. Box or Gaussian?

If several kernels remove comparable amounts of noise, what decides between them? Frequency response. The two pass entirely different bands.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig2-kernel-shape-response.png" alt="Shapes and frequency responses of box and Gaussian kernels" width="750">
_Fig 2. Kernel shapes (left) and frequency responses (right) of a width-9 box and a σ=2 Gaussian. Past its first null the box still leaves sidelobes up to 0.23 and passes 0.11 at Nyquist, while the Gaussian is effectively zero at Nyquist._

The box response does not settle to zero past its first null but repeats in sidelobes. A width-9 box peaks at $0.227$ in its largest sidelobe and still passes $0.111$ at Nyquist. A $\sigma=2$ Gaussian passes $1.5\times10^{-5}$ there, which is zero for practical purposes.

Section 6 of Post 1 is where this difference bites. Low-pass filtering before sampling is meant to prevent aliasing. A box, however, passes up to 23% in its sidelobes beyond the first null, and all of that folds straight back down. The alternating sign of the sidelobes is a second problem: some frequency components pass with inverted phase, so a bright stripe in the original can appear as a dark one.

The Gaussian became standard for reasons beyond the absence of sidelobes. A two-dimensional Gaussian separates into two one-dimensional passes, which is cheap, and repeated Gaussian filtering yields another Gaussian whose width is simply $\sqrt{\sigma_1^2+\sigma_2^2}$.

The benefit of separability grows with kernel size.

| Kernel | Direct | Two 1-D passes | Ratio |
|---|---|---|---|
| 5×5 | 25 | 10 | 2.5 |
| 9×9 | 81 | 18 | 4.5 |
| 15×15 | 225 | 30 | 7.5 |
| 25×25 | 625 | 50 | 12.5 |

The figures are multiplications per pixel. Reducing $k \times k$ to $2k$ is often what separates feasible from infeasible in real-time processing.

## 4. Why does a derivative kernel behave differently?

Finding an edge means finding where brightness changes abruptly, which means differentiating. How does this operation, the opposite of smoothing, respond to noise?

A derivative kernel has $\sum w = 0$, because it must return zero over flat regions. Its $\sum w^2$, however, is not zero, so it erases signal while keeping noise. For the central difference $[-0.5, 0, 0.5]$, $\sqrt{\sum w^2} = 0.707$. It removes only 30% of the noise while eliminating the flat component entirely.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig3-derivative-kernels.png" alt="Frequency responses of derivative kernels compared with the ideal derivative" width="700">
_Fig 3. Frequency responses of the central difference, of a kernel combining [1,2,1] along the same axis (the diagonal section of the Sobel response), and of a σ=2 DoG, compared with the ideal derivative (dotted). All three follow the ideal derivative at low frequencies; the DoG peaks at f=0.15 and pushes everything above it nearly to zero._

The ideal derivative responds in direct proportion to frequency. The central difference tracks that line well at low frequencies, sags gradually below it, peaks at $f=0.5$ with a value of $1.0$, and falls to zero at Nyquist. Patterns near Nyquist are therefore invisible to a central difference, which is worth knowing in its own right.

The real problem is that signal and noise occupy different parts of this response curve. An edge that has passed through optics is blurred, so its gradient information sits at low frequencies. That is precisely where the derivative response is smallest. Noise, by contrast, is spread evenly across all frequencies and passes freely through the mid-band where the response peaks. This is how differentiation degrades the signal-to-noise ratio.

The remedy is to fuse differentiation with smoothing. A derivative-of-Gaussian (DoG) kernel at $\sigma=2$ brings $\sqrt{\sum w^2}$ down to $0.125$, a factor of 5.7 below the central difference. It should not be confused with the difference of Gaussians, which image processing also commonly abbreviates as DoG. In Fig 3 the DoG response peaks at $f=0.15$ and pushes the frequencies above it nearly to zero. That is where the noise reduction comes from.

```python
# excerpts from kernels.py and generate_figures.py (full code: _code/imgproc-convolution-kernels/)
def noise_gain(w):
    """Variance gain against white noise: output variance = input variance x this."""
    return float(np.sum(np.asarray(w, float) ** 2))

CENTRAL_DIFF = np.array([-0.5, 0.0, 0.5])
w_dog = np.gradient(gaussian_kernel(2.0))   # derivative of a Gaussian
```

The familiar Sobel kernel has the same structure: a $[-1, 0, 1]$ derivative component multiplied by a $[1, 2, 1]$ smoothing component along the orthogonal axis. That smoothing leaves the response along the derivative axis untouched, so along that axis Sobel responds exactly like the central difference. What it does instead is suppress noise in the orthogonal direction, lowering $\sqrt{\sum w^2}$ from $0.707$ to $0.433$. The green curve in Fig 3 is the response with the same $[1, 2, 1]$ applied along the derivative axis, which equals the diagonal cross-section of the Sobel response. Scharr's $[3, 10, 3]$ revises that smoothing component, with the aim of improving the rotational symmetry of the gradient direction rather than the noise.

## 5. How much smoothing should be applied?

Everything so far was preparation. Stronger smoothing lowers the noise and stabilizes the edge position, while simultaneously blurring the edge and making its position harder to pin down. Two effects in opposite directions imply an optimum. The following experiment measures where it lies.

A synthetic edge with a known position was generated. The convolution of a step with a Gaussian point spread function (PSF) is an error function, so the edge profile is an erf; the blur width was set to $1.2$ pixels and the true position to $63.37$ pixels. The contrast of 1200 DN and the noise of 17.2 DN are carried over from section 5 of Post 1. For each smoothing width, noise was redrawn 400 times. The position was estimated by fitting a parabola through the three samples around the gradient maximum and taking its vertex.

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
_Fig 4. Position scatter (red, left axis) and edge width (blue, right axis) of a synthetic edge with 1.2-pixel blur, from 400 noise draws at each smoothing width σ. The scatter bottoms out at 0.036–0.038 pixels over the shaded range σ=2–4, while the edge width keeps growing from 3.2 to 25.8 pixels._

Without smoothing the position scatter is $0.094$ pixels. It improves 2.6-fold to $0.037$ pixels at $\sigma=2$, stays nearly flat up to $\sigma=4$, then degrades again, reaching $0.051$ pixels at $\sigma=10$. Over the same range the edge width grows monotonically from $3.2$ to $25.8$ pixels.

That the scatter worsens again is the key result. Early on, the gain from noise reduction dominates. Past a certain width the gradient peak becomes too flat for its vertex to be located reliably. A flat peak means the denominator $a-2b+c$ of the parabolic interpolation approaches zero, so small noise is amplified into large position error.

The practical rule is to keep the smoothing width within two to three times the optical PSF width. With a PSF of $1.2$ pixels, the scatter in this experiment bottomed out at $0.036$–$0.038$ pixels for smoothing widths of $2$–$4$ pixels. The floor is flat enough that 400 trials cannot locate a minimum within it, which is why Fig 4 marks it as a band rather than a point. Going far beyond the PSF spends resolution for nothing.

One more thing deserves attention: bias. Averaged over the 400 trials, the bias stayed below $0.02$ pixels throughout. Smoothing with a symmetric kernel widens an edge but does not shift it to one side. Where the small residual bias comes from shows up once the noise is removed. On a noise-free edge the bias is $-0.023$ pixels without smoothing, shrinking to $-0.007$ at $\sigma=2$ and $-0.002$ at $\sigma=4$. The bias therefore comes not from smoothing but from locating the vertex with a three-point parabola, and smoothing actually reduces it. Where subpixel estimation acquires such bias is treated separately in Post 4.

## 6. What happens at the borders?

A kernel demands neighbouring values for every output pixel. At the very edge of an image, some of those values lie outside the array. What should fill them?

There is no correct answer. Something must be assumed about the outside, and each mode assumes something different. Zero padding assumes darkness beyond the border; replication assumes the final value persists; reflection assumes the signal folds back symmetrically. Even "reflect" means different things in different libraries. numpy's `reflect` and the OpenCV default do not repeat the edge pixel, whereas scipy.ndimage's `reflect` does. The version used here does not repeat it.

A signal that genuinely continues past the array was constructed, and the central 96 pixels were taken as the captured image. The ground truth is the full-length signal smoothed in one piece.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig5-boundary-modes.png" alt="Error comparison of zero, replicate and reflect padding modes" width="750">
_Fig 5. Error left by three padding modes when a sloped scene, of which only the central 96 pixels were captured, is filtered with a σ=3 Gaussian. Zero padding deviates by up to 493 DN at the border (left); replication and reflection leave about one tenth and one fifth of that in RMS over the outer five pixels (right, rescaled vertical axis)._

Zero padding produces errors an order of magnitude apart from the others: $191$ DN RMS over the outer five pixels, peaking at $493$ DN. This is the familiar darkening of a bright image's border, and searching for edges on top of it yields a strong spurious edge tracing the frame.

The other two are five to ten times smaller. Replication gives $18.7$ DN RMS and reflection $37.4$ DN. Reflection came out worse here because this scene carries a gradient across the border. Folding the signal also inverts the sign of that gradient, creating a kink. Had the scene been flat near the border, the ranking could reverse.

No mode is therefore universally better. In metrology the reliable solution is not to choose the assumption well but to avoid depending on it. Shrinking the region of interest (ROI) by half the kernel width keeps invented values from ever reaching the result. For a Gaussian, the half-width is the radius at which the kernel is truncated, commonly $4\sigma$; with $\sigma=3$ in this experiment, that means 12 pixels. A little field of view is traded for the complete removal of border artefacts.

## 7. Why can't a Gaussian remove a single dead pixel?

The discussion so far assumed noise of comparable size on every pixel. Dead pixels and cosmic ray hits break that assumption: a single point or a short streak takes an entirely wrong value. Fig 6 shows the traces cosmic rays left on an image taken with the lens covered. No light reached the sensor, yet bright dots and short streaks are scattered across it. Is a Gaussian still the answer?

<img src="/assets/img/posts/imgproc-convolution-kernels/ext-cosmic-ray-dark-frame.jpg" alt="Cosmic ray traces on a CCD dark frame taken with the lens covered" width="700">
_Fig 6. Part of a CCD dark frame taken with the lens covered. Cosmic ray hits appear as dots one or two pixels across, short curved worms a few pixels long, and one long straight track crossing the sensor at an angle. (Source: RadiceEdoardoLuca, [The Sound of Silence](https://commons.wikimedia.org/wiki/File:The_Sound_of_Silence.jpg), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0), cropped and brightness adjusted)_

The experiment takes the simplest case and plants only single-pixel spikes. Extreme values were planted on 1.01% of the pixels (259) and both filters were compared against the clean original.

<img src="/assets/img/posts/imgproc-convolution-kernels/en/fig6-median-vs-gaussian.png" alt="Gaussian and median filters applied to an image with spike pixels" width="750">
_Fig 7. A test image with 0 or 255 planted on 1.01% of its pixels (left), filtered with a σ=1.5 Gaussian (centre) and a 3×3 median (right). The Gaussian smears each dot into a grey blob; the median removes the dots._

The numbers are unambiguous. The spiked image has an RMS error of $14.5$ DN against the original. A Gaussian raises it to $15.3$ DN, so **it makes things worse**. A median filter brings it down to $9.2$ DN. Most of both figures, however, comes from blurring image detail rather than from the spikes. Applied to the spike-free original, the same filters already give $14.9$ DN for the Gaussian and $9.2$ DN for the median. Isolating the spike contribution, the Gaussian reduces it to $2.8$ DN and the median to $1.8$ DN. Total RMS alone cannot separate the two filters fairly; the spike share does.

The Gaussian cannot remove a spike cleanly because of linearity itself. A linear filter returns a weighted sum of its input, so it cannot destroy the total error an outlier brings; it can only distribute that error among the neighbours. The grey blobs where white dots used to be, in the centre panel of Fig 7, are the result. A single contaminated point has been widened to the kernel width, so the peak error drops while the contaminated area grows.

A median filter sorts and takes the middle value. Being non-linear, it discards outliers outright. The cost is real, though. It is not a convolution, so frequency response says nothing about it, and it erases lines narrower than half its window and rounds off sharp corners. Where fine structure is itself the measurand, a median filter cannot be applied casually.

## Summary and what comes next

Choosing a kernel fixes its noise reduction through $\sum w^2$ and its signal preservation through the frequency response. Those two numbers are the whole of smoothing kernel selection. A derivative kernel has $\sum w = 0$ with a response that grows with frequency, so it must be fused with smoothing, which is what the DoG does; the Sobel family applies its smoothing orthogonally to cut noise. Smoothing width has an optimal range, here around two to three times the PSF width. Borders and spikes are the two places where the assumptions behind linear convolution break, answered respectively by shrinking the region of interest and by switching to a non-linear filter.

The next post moves to the frequency domain. Having explained kernel behaviour through frequency response, it asks what happens when filters are designed there directly. It covers the removal of periodic noise such as illumination ripple and grating patterns, which cannot be touched in the spatial domain. It then turns to why an ideal filter produces ringing, and what the extrema envelope method introduced in the author's dissertation does about it.

## References

- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 3 (spatial filtering), ch. 10 (edge detection).
- H. Scharr, "Optimale Operatoren in der Digitalen Bildverarbeitung" (Optimal Operators in Digital Image Processing), Ph.D. dissertation, Universität Heidelberg, 2000, [doi:10.11588/heidok.00000962](https://doi.org/10.11588/heidok.00000962) (optimizing gradient kernels for rotational symmetry).
- J. Canny, "[A Computational Approach to Edge Detection](https://doi.org/10.1109/TPAMI.1986.4767851)," _IEEE Trans. Pattern Anal. Mach. Intell._, PAMI-8(6), 679-698, 1986 (the case for fusing differentiation with smoothing).
- I. T. Young, L. J. van Vliet, "[Recursive implementation of the Gaussian filter](https://doi.org/10.1016/0165-1684(95)00020-E)," _Signal Processing_, 44(2), 139-151, 1995 (separability and fast implementations).
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the test image in Fig 7, `skimage.data.camera`, photo by Lav Varshney, CC0).
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025.

---
title: "Metrology Image Processing 4 — Edges and Subpixels: The Bias That Tracks the Grid"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-edge-subpixel/
page_id: imgproc-edge-subpixel
date: 2026-10-20 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, edge-detection, subpixel, pixel-locking, bias-variance, metrology]
description: "Three ways to measure something smaller than a pixel, applied to the same edge, exposing the systematic bias each one repeats along the pixel grid."
math: true
---

[Post 2](/en/posts/imgproc-convolution-kernels/) located an edge at the vertex of a parabola fitted to the gradient peak while searching for the optimal smoothing width, and [Post 3](/en/posts/imgproc-frequency-domain-filtering/) refined extremum wavenumbers the same way. Both times the interpolation itself was deferred to Post 4. This is that post.

Two questions. How is a position finer than the pixel spacing read out at all? And why does the value so obtained carry an error that repeats along the pixel grid? The second question is the heart of this post. That error is systematic rather than random, so averaging repeated measurements does not remove it.

## 1. An edge is blurred twice

What produces the edge profile a sensor records? Two stages.

First the point spread function (PSF) of the optics blurs the step. Taking the blur as Gaussian, its convolution with a step is an error function, so the profile is an erf.

$$ I(x) = B + (A - B)\cdot\frac{1}{2}\left[1 + \mathrm{erf}\!\left(\frac{x - x_0}{\sqrt{2}\,w}\right)\right] $$

Then each pixel integrates the light over its own area. This is where the statement from section 6 of [Post 1](/en/posts/imgproc-sampling-noise-model/) — that a pixel is an area, not a point — actually gets used. Because of that second stage, a sample is not a point on the curve but the curve averaged over a pixel width.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig1-edge-formation.png" alt="An ideal step, the curve after the PSF, and the values sampled by pixels" width="700">
_Fig 1. An edge is blurred once by the optics and again by pixel integration_

Here $w$ is the blur width; note that Post 2 used the same letter for kernel weights. Blur width is not directly observable, so in practice the edge is measured by the distance over which the intensity rises from 10% to 90%. For Gaussian blur that distance is $2.563\,w$. Every result below is governed by the ratio of this edge width to the pixel spacing. What limits subpixel accuracy is that ratio, more than noise.

## 2. Three ways to read a subpixel position

Why is there more than one way to extract a single position from the same profile? Because each method assumes something different about the profile, and the differences in assumption become the differences in performance.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig2-three-families.png" alt="What the parabola, moment and area methods read from the same profile" width="750">
_Fig 2. The three families read different things from the same profile_

The interpolation family fits a parabola to the peak of the gradient magnitude and takes its vertex. Using only three samples, it is the cheapest, and it is what Posts 2 and 3 used. Its assumption is that the peak is a parabola.

The moment family inverts the first three moments of the intensity within a window to recover the position of a step. Tabatabai and Mitchell proposed it in 1984. Assuming the window holds one ideal step with the low side occupying a fraction $p$, the three moments give $p$ in closed form.

$$ p = \frac{1}{2}\left(1 + s\sqrt{\frac{1}{4 + s^2}}\right), \qquad s = \frac{m_3 - 3m_1 m_2 + 2m_1^3}{\sigma^3} $$

Here $m_k$ is the $k$-th moment of the intensities in the window, $\sigma^2 = m_2 - m_1^2$, and $s$ is the skewness.

The area family uses only the fact that a pixel value is the average of the light that pixel received. Taking the plateaus $A$ and $B$ and summing the normalized profile $(A - v)/(A - B)$ across the window gives the distance from the left end of the window to the edge. As long as the blur is symmetric about the edge, neither its shape nor its width needs to be assumed.

```python
# core of subpixel.py (full code: _code/imgproc-edge-subpixel/)
def partial_area(profile, plateau=3, window="auto"):
    lo_i, hi_i = transition_window(profile) if window == "auto" else (0, len(profile))
    seg = profile[lo_i:hi_i]
    lo, hi = seg[:plateau].mean(), seg[-plateau:].mean()
    return lo_i - 0.5 + np.sum((hi - seg) / (hi - lo))
```

## 3. A bias that repeats along the grid

What appears when the true position is shifted in small steps within a single pixel and the error of each method is measured? No noise was added, so everything seen here is systematic.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig3-pixel-locking.png" alt="Bias of the three methods against the fractional true position, at four blur widths" width="750">
_Fig 3. The bias of parabolic interpolation repeats with the pixel grid_

The error of parabolic interpolation traces a smooth S-curve with a period of one pixel. At $w = 0.6$ it reaches $0.042$ pixels, and it is exactly zero where the fractional position is $0$ or $0.5$. When the edge sits on a pixel centre or a pixel boundary, the sample points are arranged symmetrically about it and the interpolation error cancels; in between, it does not.

This is pixel locking: the estimate is pulled toward grid positions. It is dangerous in metrology because it is not random error. The same edge measured repeatedly is wrong by the same amount in the same direction, so frame averaging leaves it untouched. The worse case is measuring many edges whose grid phases happen to be spread evenly. The systematic bias then masquerades as scatter, and the result looks imprecise but accurate when it is neither.

The moment method is not immune. At $w = 0.4$ its bias reaches $0.036$ pixels, comparable in size to the parabola but opposite in phase. Pixel integration breaks its assumption of an ideal step, and the way it breaks depends on the grid phase. It does fall away quickly as blur grows.

The area curve lies almost flat at zero for $w = 0.4$ and $0.6$, since it assumes nothing beyond symmetry and is therefore indifferent to grid phase. In the $w = 2.0$ panel, however, it swings by $\pm0.011$ pixels. That is not a property of the method but of cropping the window to the transition, whose cost section 5 treats separately. The small steps visible in the curve mark where the window length changes by a whole pixel.

## 4. The blur width that minimizes bias

How does the bias depend on blur width, then? A sharper edge feels like it should be easier to locate. It is not.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig4-bias-vs-blur.png" alt="Maximum absolute bias of the three methods against PSF blur width" width="700">
_Fig 4. A peak, a U and a rise — each method has its own optimal blur_

Parabolic interpolation is worst at $w = 0.6$ and improves in both directions. With a very sharp edge only two or three samples land on the gradient peak, leaving nothing to fit; with enough blur the peak genuinely approaches a parabola and the assumption holds. The worst case lies between.

The moment method traces the opposite shape, a U with its minimum of $0.0016$ pixels at $w = 0.8$. Too sharp and pixel integration breaks the step assumption; too blurred and nothing resembling a step remains inside the window.

The area method stays below $10^{-3}$ pixels for $w \le 0.6$, but its bias grows with blur and reaches $0.011$ pixels at $w = 2.0$. This is not the symmetry assumption failing. It comes from cropping the window to the transition: with the window uncropped, the bias at the same $w = 2.0$ stays at $0.00002$ pixels. As blur grows the tails of the transition lengthen, and estimating $A$ and $B$ at window ends that have not yet reached the plateaus introduces exactly that offset.

A practical rule follows. Capturing the sharpest possible edge is not the best move. Matching the optics and the pixel spacing so that the 10-90% edge width lands at 2 to 3 pixels serves better. Section 6 of Post 1 argued for enough blur to keep content above Nyquist from folding down; subpixel bias points to the same answer. Under-blurring buys aliasing and grid bias together.

## 5. Smaller bias, larger scatter

If the area method is nearly exact and unbiased, why not simply use it? Because none of the experiments so far contained noise.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig5-noise-scatter.png" alt="Scatter of the three methods against noise, and the window-length dependence of the area method" width="750">
_Fig 5. Smaller bias comes with larger scatter, and the area window is a trade of its own_

On the left is the scatter of the estimated position at $w = 1.0$ as noise increases. The ordering is exactly the reverse of the bias ordering. At a noise level of 5 DN the moment method is the most stable at $0.0088$ pixels, followed by the area method at $0.0149$ and the parabola at $0.0217$.

The reason lies in how many pixels each method uses and how. The moment method forms an average over the whole window, so noise averages down. The parabola uses three gradient samples, and section 4 of Post 2 showed that differentiation has already amplified the noise in them. The area method sums every pixel in the window, so the noise accumulates with the signal.

The right panel shows that accumulation directly. The scatter of the area method grows as the square root of the window length $n$, from $0.0074$ at 8 pixels to $0.1015$ at 64. Shortening the window instead lets the transition overflow it, so the plateaus are misjudged and the bias explodes to $0.87$ pixels at 6 pixels. The practical window length is the 12 to 16 pixels where the bias has vanished. The automatic window, which finds the transition and adds a 3-pixel margin, chose a length of 12 at $w = 1.0$ and more than halved the scatter relative to using the full profile.

That choice carries the price seen in section 4. The automatic window grows with blur but its 3-pixel margin does not, so the longer the tails of the transition, the further the plateau estimate drifts. Halving the scatter bought a bias of $0.011$ pixels at $w = 2.0$. Scaling the margin with blur width reduces the bias and brings the scatter back.

## 6. Which method to use

Since bias and scatter move in opposite directions, only their combination answers the question. Taking the total error as $\sqrt{\text{bias}^2 + \text{scatter}^2}$, blur width and noise were swept together.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig6-total-error.png" alt="Total error of the three methods against blur width and noise" width="750">
_Fig 6. Blur width and noise change the winner_

For $w \ge 0.8$ — a 10-90% edge width beyond 2 pixels — the moment method wins at every noise level. Its bias is already down at the $0.002$ pixel level, leaving only scatter, and its scatter is the smallest of the three.

Sharper edges reverse the order. At $w = 0.3$ to $0.5$ with low noise the area method wins, because this is the range where the moment method's step assumption breaks and its bias climbs to $0.06$ pixels, favouring the unbiased option. At the same $w = 0.3$ with higher noise the parabola wins instead: the area method's scatter grows in proportion to noise, while the parabola holds a small $0.02$ pixel bias on sharp edges and uses few samples.

The practical reading is that the optics come before the method. If the edge width can be set to 2 or 3 pixels, the moment method is the simple choice, and in that regime the three methods barely differ anyway. Where the edge width cannot be arranged, noise level decides between the area method and the parabola.

Section 4.2 of the dissertation contains a case where this matters. To locate the back focal plane of an objective, the imaging system was translated along the optical axis while the focal spot was photographed, and the position of minimum spot diameter was taken as the back focal plane. The diameter was measured by image processing. If grid bias enters that diameter measurement, it distorts the diameter-versus-position curve and shifts the location of the minimum. The result was 78.9 mm from the relay lens and 6.1 mm inside the rear of the objective — numbers that deserve confidence only once the bias of the diameter estimate is known to be smaller than its scatter.

## 7. When it fails

All three methods share one premise: exactly one edge inside the window. The failure conditions are the situations that break it.

Illumination tilted to one side leaves the plateaus non-flat, and the area method falls first. It takes $A$ and $B$ from the window ends, and if those values ride a gradient the normalization is wrong from the start. The background slope must be removed first, or the window shortened to limit its influence.

Two edges inside one window derail all three. The parabola picks the taller peak, while the moment and area methods return an intermediate average of the two. This arises often when measuring fine lines or narrow grooves, and it calls for a method that models the edge pair together.

At the border of an image the window cannot be formed at all. Section 6 of Post 2 argued for shrinking the region of interest by half the kernel width rather than inventing values for the outside, and the same principle applies here.

Two practical failure checks exist. One is to shift the window by a pixel or two and watch whether the estimate moves; a healthy estimate is insensitive to window placement. The other is to run all three methods and compare. Three methods resting on different assumptions agreeing on a value is grounds for trusting it, and a wide disagreement signals that one of the premises has broken.

## Summary and what comes next

Reading a position finer than a pixel means assuming something about the profile, and when that assumption breaks differently at different grid phases, a systematic bias with a period of one pixel appears. Parabolic interpolation is worst near $w = 0.6$ at $0.042$ pixels; the moment method carries a bias of opposite phase that all but vanishes at $w = 0.8$; the area method, assuming only symmetry, is essentially unbiased with an uncropped window at the cost of scatter growing as the square root of the window length, and the reverse once the window is cropped. The ranking of the three shifts with blur and noise, but matching the optics to a 10-90% edge width of 2 to 3 pixels makes the choice itself matter less.

The next post turns to geometric fitting. Having obtained the position of one edge point to subpixel precision, the problem becomes extracting a single parameter — a line, a circle — from hundreds of such points. What gets minimized changes the answer, and on short arcs the commonly used algebraic methods carry a pronounced bias. The trade between bias and scatter seen here returns there.

## References

- A. J. Tabatabai, O. R. Mitchell, "[Edge Location to Subpixel Values in Digital Imagery](https://doi.org/10.1109/TPAMI.1984.4767502)," _IEEE Trans. Pattern Anal. Mach. Intell._, PAMI-6(2), 188-201, 1984 (moment-based subpixel estimation).
- A. Trujillo-Pino et al., "[Accurate subpixel edge location based on partial area effect](https://doi.org/10.1016/j.imavis.2012.10.005)," _Image and Vision Computing_, 31(1), 72-90, 2013 (the area-based method).
- J. Canny, "[A Computational Approach to Edge Detection](https://doi.org/10.1109/TPAMI.1986.4767851)," _IEEE Trans. Pattern Anal. Mach. Intell._, PAMI-8(6), 679-698, 1986.
- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 10 (edge detection).
- Y. Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," Ph.D. dissertation, Seoul National University, 2025, sec. 4.2 (locating the back focal plane from focal spot diameter).

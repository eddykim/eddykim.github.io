---
title: "Metrology Image Processing 4 — How Far Can You Trust an Edge Position Finer Than a Pixel?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-edge-subpixel/
page_id: imgproc-edge-subpixel
date: 2026-10-20 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, edge-detection, subpixel, pixel-locking, bias-variance, metrology]
description: "Three ways to read a position finer than a pixel, applied to the same edge, exposing the systematic bias that repeats along the pixel grid."
math: true
---

[Post 2](/en/posts/imgproc-convolution-kernels/) located an edge at the vertex of a parabola fitted to the gradient peak while searching for the optimal smoothing width, and [Post 3](/en/posts/imgproc-frequency-domain-filtering/) refined extremum wavenumbers the same way. Both times the interpolation itself was deferred to Post 4. This is that post.

Two questions. How is a position finer than the pixel spacing read out at all? And why does the value so obtained carry an error that repeats along the pixel grid? The second question is the heart of this post. That error is systematic rather than random, so averaging repeated measurements does not remove it.

A boundary found at pixel level follows the pixel corners in a staircase. A boundary found at subpixel level follows the actual outline smoothly (Fig 1).

<img src="/assets/img/posts/imgproc-edge-subpixel/ext-subpixel-edge-angiogram.png" alt="Vessel boundaries in an enlarged part of an angiogram, found at pixel level and at subpixel level" width="700">
_Fig 1. Vessel boundaries found in an enlarged part of an angiogram. On the left the boundary is found at pixel level (yellow squares); on the right the same region is processed at subpixel level (yellow dashed line). The pixel-level boundary steps along the grid, while the subpixel boundary resolves positions inside each pixel and follows the vessel wall smoothly. The full angiogram at the centre of the original has been removed, leaving the two enlargements side by side. (Source: Machanguillo1, [Subpixel edge detection](https://commons.wikimedia.org/wiki/File:Subpixel_edge_detection.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), cropped)_

## 1. What does an edge look like on the sensor?

What produces the edge profile a sensor records? The edge is blurred twice.

First the point spread function (PSF) of the optics blurs the step. Taking the blur as Gaussian, its convolution with a step is an error function, so the profile is an erf.

$$ I(x) = B + (A - B)\cdot\frac{1}{2}\left[1 + \mathrm{erf}\!\left(\frac{x - x_0}{\sqrt{2}\,w}\right)\right] $$

Here $B$ is the intensity of the dark plateau, $A$ that of the bright plateau, $x_0$ the true position and $w$ the blur width. Post 2 used the same letter for kernel weights. As $w \to 0$ the erf returns to a step, the ideal edge.

Then each pixel integrates the light over its own area. Section 6 of [Post 1](/en/posts/imgproc-sampling-noise-model/) said that a pixel is an area, not a point, and that fact is used here. Because of this second stage, a sample is not a point on the curve but the curve averaged over a pixel width. In Fig 2 the dotted step passes through the PSF to become the solid curve, and the red dots are that curve averaged pixel by pixel.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig1-edge-formation.png" alt="An ideal step, the curve after the PSF, and the values sampled by pixels" width="700">
_Fig 2. An edge rising from 100 to 1300 with blur width w=1.0 pixel and true position 11.37. The ideal step (dotted) becomes an erf curve (solid) after the PSF, and the pixels record that curve averaged over each pixel (dots). The edge is blurred once by the optics and again by pixel integration._

Blur width is not directly observable. In practice the blur is measured as the 10-90% edge width, the distance over which the intensity rises from 10% to 90%. For the Gaussian blur alone this distance is $2.563\,w$. Measured on pixel samples, which adds the pixel aperture, it is slightly wider: 2.18 pixels at $w = 0.8$ and 2.67 pixels at $w = 1.0$. Every bias below is governed by the ratio of this edge width to the pixel spacing. Once noise enters, its level matters too (sections 5 and 6), but the systematic error of subpixel estimation is set by this ratio rather than by noise.

## 2. Why is there more than one way to read a position?

Each method assumes something different about the profile, so each performs differently. The rough location of an edge is found first by a pixel-level method such as the Canny (1986) detector. Subpixel estimation then refines the position from the profile around it, and the methods fall into three families.

The interpolation family fits a parabola to the peak of the gradient magnitude and takes its vertex. Using only three samples, it is the cheapest, and it is what Post 2 used. Post 3 found the vertex by fitting a quadratic to nine samples. Its assumption is that the peak is a parabola.

The moment family inverts the first three moments of the intensity within a window to recover the position of a step. Tabatabai and Mitchell proposed it in 1984. Here the moment window is the whole 24-pixel profile. Assuming the window holds one ideal step with the low side occupying a fraction $p$, the three moments give $p$ in closed form.

$$ p = \frac{1}{2}\left(1 + s\sqrt{\frac{1}{4 + s^2}}\right), \qquad s = \frac{m_3 - 3m_1 m_2 + 2m_1^3}{\sigma^3} $$

Here $m_k$ is the $k$-th moment of the intensities in the window, $\sigma^2 = m_2 - m_1^2$, and $s$ is the skewness. The intuition lies in the intensity histogram. An ideal step leaves only two values in the window, dark and bright. The more dark pixels there are, the longer the histogram's tail towards the bright side, and the skewness measures that lean. The lean gives the dark fraction $p$, and $p$ times the window length is the distance from the left end of the window to the edge. When $s = 0$, $p = 1/2$ and the edge sits at the centre of the window.

The area family uses only the fact that a pixel value is the average of the light that pixel received. Taking the plateaus $A$ and $B$ and summing the normalized profile $(A - v)/(A - B)$ across the window gives the distance from the left end of the window to the edge. As long as the blur is symmetric about the edge, neither its shape nor its width needs to be assumed.

The area method amounts to counting dark pixels. A fully dark pixel counts as one, a half-dark pixel as one half. With symmetric blur, the darkness spilled to the right of the edge exactly cancels the brightness spilled to the left, so the count is the distance to the edge.

Fig 3 shows what each family reads from the same profile. On the left, the vertex of the parabola through the three gradient samples lies 0.028 pixels left of the true position (dotted line). In the centre, the step recovered by the moments is drawn per pixel, so it rises at the pixel boundary (11.5), 0.13 pixels from the true position (11.37). For this profile $s = -0.022$ and $p = 0.495$, so the edge lies slightly left of the window centre (11.5). On the right, the first nine pixels are almost fully dark (summing to 9.00), and the following ones (0.99, 0.91, 0.64, 0.27, 0.06, …) add 2.87, so the cumulative sum stops at 11.87. The window's left end is at −0.5, which puts the edge at 11.37, the true position.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig2-three-families.png" alt="What the parabola, moment and area methods read from the same profile" width="750">
_Fig 3. What the three families read from the profile of Fig 2 (w=1.0, true position 11.37). Left: a parabola fitted to three gradient-magnitude samples. Centre: the two-level step recovered from the intensity moments. Right: the normalized profile (bars) and its cumulative sum (line). The dotted line is the true position, and the final value of the cumulative sum, 11.87, is the distance from the window's left end (−0.5) to the edge._

```python
# core of subpixel.py (full code: _code/imgproc-edge-subpixel/)
def partial_area(profile, plateau=3, window="auto"):
    lo_i, hi_i = transition_window(profile) if window == "auto" else (0, len(profile))
    seg = profile[lo_i:hi_i]
    lo, hi = seg[:plateau].mean(), seg[-plateau:].mean()
    return lo_i - 0.5 + np.sum((hi - seg) / (hi - lo))
```

## 3. Why does the error repeat pixel after pixel?

What appears when the true position is shifted in small steps within a single pixel and the error of each method is measured? No noise was added, so everything seen here is systematic.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig3-pixel-locking.png" alt="Bias of the three methods against the fractional true position, at four blur widths" width="750">
_Fig 4. Bias of the three methods as the fractional part of the true position moves through one pixel in 120 steps, without noise. The four panels are blur widths w=0.4, 0.6, 1.0 and 2.0. The bias of parabolic interpolation (red) repeats as a one-pixel S-curve and reaches 0.042 pixels at w=0.6; the moment method (green) has the opposite phase; the area method (blue) stays at zero for small blur._

The error of parabolic interpolation traces a smooth S-curve with a period of one pixel. At $w = 0.6$ it reaches $0.042$ pixels, and it is exactly zero where the fractional position is $0$ or $0.5$. When the edge sits on a pixel centre or a pixel boundary, the sample points are arranged symmetrically about it and the interpolation error cancels; in between, it does not.

Reading between the marks of a ruler by eye, one tends to drift towards the nearest mark. Subpixel estimation has the same habit, called pixel locking: the estimate is pulled towards grid positions. It is dangerous in metrology because it is not random error. The same edge measured repeatedly is wrong by the same amount in the same direction, so frame averaging leaves it untouched. The trickier case is measuring many edges whose grid phases happen to be spread evenly. The systematic bias then looks like scatter. Each individual value is biased, yet the spread passes for ordinary imprecision.

The moment method also locks. At $w = 0.4$ its bias reaches $0.036$ pixels, comparable in size to the parabola but opposite in phase. The intermediate values produced by blur and pixel integration break its assumption of an ideal two-level step, and the way it breaks depends on the grid phase. It does fall away quickly as blur grows.

The area curve lies almost flat at zero for $w = 0.4$ and $0.6$, since it assumes nothing beyond symmetry and is therefore indifferent to grid phase. Its swing of $\pm0.011$ pixels in the $w = 2.0$ panel comes from cropping the window to the transition region. Section 4 explains the cause and section 5 the price. The small steps visible in the curve mark where the window length changes by a whole pixel.

## 4. Is a sharper edge always better?

How does the bias depend on blur width, then? A sharper edge feels like it should be easier to locate. The answer depends on the method.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig4-bias-vs-blur.png" alt="Maximum absolute bias of the three methods against PSF blur width" width="700">
_Fig 5. Maximum absolute bias of the three methods against blur width w (no noise, maximum over 120 phase steps). The top axis converts w to the 10-90% edge width as 2.563w, and values below 10⁻⁶ are drawn at that floor. The parabola forms a peak, the moment method a U and the area method a rise with blur, so the blur width with the smallest bias differs from method to method._

Parabolic interpolation is worst at $w = 0.6$ and improves in both directions. With a very sharp edge, pixel integration and the central difference make the three gradient samples shift linearly with the true position, so the vertex nearly coincides with it (zero bias as $w \to 0$). With enough blur the peak genuinely approaches a parabola and the assumption holds. The worst case lies between.

The moment method traces the opposite shape, a U with its minimum of $0.0016$ pixels at $w = 0.8$. Too sharp, and pixel integration breaks the step assumption. Too blurred, and the tails of the transition reach both ends of the 24-pixel window, and the more the edge sits off the window centre, the more asymmetrically the window clips them. With the edge exactly at the centre (11.5) this bias is zero. Widening the window to 48 pixels with the edge at its centre cuts the maximum bias at $w = 2.0$ from $0.0096$ to $0.0020$ pixels.

The area method stays below $10^{-3}$ pixels for $w \le 0.6$, but its bias grows with blur and reaches $0.011$ pixels at $w = 2.0$. This is not the symmetry assumption failing. It comes from cropping the window to the transition: with the window uncropped, the bias at the same $w = 2.0$ stays at $0.00002$ pixels. As blur grows, the tails of the transition lengthen. If the ends of the cropped window have not reached the plateaus, the $A$ and $B$ taken there are off by that much.

A practical rule follows for the moment method: do not capture the sharpest edge possible. Matching the optics and the pixel spacing so that the 10-90% edge width lands at 2 to 3 pixels serves better. Section 6 of Post 1 argued for enough blur to keep content above Nyquist from folding down. The bias of the moment method points to the same answer. Parabolic interpolation, however, is worst near $w = 0.6$ and still leaves $0.025$ to $0.038$ pixels at a 2-to-3-pixel edge width, more than on sharper edges ($0.009$ pixels at $w = 0.2$). The blur width has to be chosen together with the method.

## 5. Is an unbiased method enough?

If the area method is nearly exact and unbiased, it seems enough to use it alone. But none of the experiments so far contained noise. Here Gaussian noise is added to the profile (24 pixels, intensity 100 to 1300 DN, true position 11.37) and each estimate is repeated 600 times.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig5-noise-scatter.png" alt="Scatter of the three methods against noise, and the window-length dependence of the area method" width="750">
_Fig 6. Left: standard deviation (scatter) of the estimated position at w=1.0 as noise rises from 1 to 40 DN, 600 trials each. Right: bias (no noise) and scatter (noise 5 DN, 400 trials) of the area method against window length, with the whole profile used as the window. The dotted line is a reference proportional to the window length. The grey band marks the degenerate 6-pixel window, where the 3+3 plateau pixels cover the whole window and the estimate is fixed regardless of noise. Values below 5×10⁻⁵ are drawn at that floor. The method with smaller bias has larger scatter, and the area window is a trade of its own._

The left panel is the scatter at $w = 1.0$ as noise increases. The ordering does not follow the bias ordering. The area method, which had the smallest bias at sharp blur ($w \le 0.6$), scatters more than the moment method here. At 5 DN the moment method is the most stable at $0.0088$ pixels, followed by the area method at $0.0149$ and the parabola at $0.0217$. Beyond 10 DN the area method scatters more than even the parabola ($0.0443$ against $0.0418$ at 10 DN, $0.137$ against $0.084$ at 20 DN).

The reason lies in how many pixels each method uses and how. The moment method forms an average over the whole window, so noise averages down. The parabola uses three gradient samples, and section 4 of Post 2 showed that differentiation has already degraded their signal-to-noise ratio. The area method sums every pixel in the window and normalizes the sum by plateaus taken from a few end pixels, so the noise accumulates with the window length.

The right panel shows that accumulation directly. The scatter of the area method grows faster than the square root of the window length $n$, nearly in proportion to $n$. From $0.0216$ at 16 pixels, where the bias has vanished, it rises to $0.1015$ at 64: four times the window, 4.7 times the scatter. The plateaus $A$ and $B$ are averages of three pixels at each end, so the noise in those averages enters every term of the sum. Shortening the window instead lets the transition overflow it, the plateaus are misjudged, and at 6 pixels the bias shoots up to $0.87$ pixels. The scatter there looks like zero, but not because the estimate is precise. The 3+3 plateau pixels cover the whole window and the estimate is stuck at 2.5. The practical window length is around 16 pixels, where the bias has vanished; even at 12 pixels the bias is a small $0.010$ pixels.

The automatic window, which finds the transition and adds a 3-pixel margin, chose a length of 12 at $w = 1.0$ and more than halved the scatter relative to using the full profile. The 12-pixel bias in the figure above ($0.010$ pixels) differs from the automatic window's bias in section 4 (about $0.003$ pixels) most likely because the windows are taken differently: the figure uses an entire 12-pixel profile as the window, whereas the automatic window cuts 12 pixels out of a longer profile, positioned by the transition. The automatic window can, however, be fooled by noise: once noise crosses the gradient threshold, the transition is taken to be wider. The mean window length stays at 12 pixels up to 5 DN, then grows to 14.5 pixels at 10 DN and 21 pixels at 20 DN. This is why the area curve in the left panel steepens between 10 and 20 DN.

The price of that choice is the area-method bias seen in section 4. The automatic window grows with blur, but its 3-pixel margin does not, so the longer the tails of the transition, the further the plateau estimate drifts. Halving the scatter bought a bias of $0.011$ pixels at $w = 2.0$. Scaling the margin with blur width reduces the bias and brings the scatter back.

## 6. Which method should be used?

The short answer is that the optics come before the method. If the edge width can be set to 2 or 3 pixels, use the moment method; if not, let the noise level decide between the area method and the parabola. The basis is the total error, which combines bias and scatter. Since the two move in opposite directions, the total error was taken as $\sqrt{\text{bias}^2 + \text{scatter}^2}$ and blur width and noise were swept together.

<img src="/assets/img/posts/imgproc-edge-subpixel/en/fig6-total-error.png" alt="Total error of the three methods against blur width and noise" width="750">
_Fig 7. Total error of the three methods, as root-mean-square error (RMSE), against blur width w and noise (1, 3, 10 and 30 DN). The bias is the RMS over 120 phase steps and the scatter comes from 300 trials. Blur width and noise change the winner._

For $w \ge 0.8$, a 10-90% edge width beyond 2 pixels, the moment method wins at every noise level. Its bias is already down at the $0.002$ pixel level, leaving only scatter, and its scatter is the smallest of the three.

Sharper edges reverse the order. At $w = 0.3$ to $0.5$ with low noise the area method wins, because this is the range where the moment method's step assumption breaks and its bias climbs to $0.06$ pixels, favouring the unbiased option. At the same $w = 0.3$ with noise around 10 DN the parabola wins instead. The area method's scatter grows faster than the noise (section 5), while the parabola keeps a small $0.02$ pixel bias on sharp edges and uses few samples.

With the edge width at 2 or 3 pixels, the moment method is the simple choice. In that regime the moment and area methods both come down to about $0.002$ pixels of bias. Only parabolic interpolation still leaves a grid bias of about $0.03$ pixels.

Section 4.2 of the author's dissertation contains a case where this matters. To locate the back focal plane of an objective, the imaging system was translated along the optical axis while the focal spot was photographed, and the position of minimum spot diameter was taken as the back focal plane. The diameter was measured by image processing. If grid bias enters that diameter measurement, it distorts the diameter-versus-position curve and shifts the location of the minimum. The result was 78.9 mm from the relay lens and 6.1 mm inside the rear of the objective. Trusting these numbers first requires checking that the bias of the diameter estimate is smaller than its scatter.

## 7. When does it fail?

All three methods share one premise: exactly one edge inside the window. The failure conditions are the situations that break it.

Illumination tilted to one side leaves the plateaus non-flat, and the area method falls first. It takes $A$ and $B$ from the window ends, and if those values ride a gradient the normalization is wrong from the start. The background slope must be removed first, or the window shortened to limit its influence.

Two edges inside one window derail all three. The parabola picks the taller peak. With two edges of the same direction, the moment and area methods return a value somewhere between them. With an opposite pair, as in a fine line or a narrow groove, the two plateaus coincide, the area method diverges and the moment method returns a value outside both edges. Such pairs are common in linewidth measurement, and they call for a method that models the edge pair together. Fig 8 shows a profile across three fine bars in a real microscope image.

<img src="/assets/img/posts/imgproc-edge-subpixel/ext-foldscope-usaf-profile.png" alt="A resolution target imaged with a paper microscope, and an intensity profile across three fine bars" width="750">
_Fig 8. A USAF 1951 resolution target imaged with a paper microscope (Foldscope) (left), an enlargement of part of it (centre), and the intensity profile along the green line in the enlargement (right). In the graph the green line (Data) is the measurement and the black line (Ideal) the ideal bar profile. The bars are so fine that edges of opposite direction sit close together, and the measured profile never climbs back to the plateau between bars. For a profile like this, the premise of a single edge in the window fails from the outset. (Source: Cybulski J, Clements J, Prakash M, [Foldscope: Origami-Based Paper Microscope](https://doi.org/10.1371/journal.pone.0098781), panels E and F of the paper's figure, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0), cropped)_

At the border of an image the window cannot be formed at all. Section 6 of Post 2 argued for shrinking the region of interest by half the kernel width rather than inventing values for the outside. The same principle applies here.

Two practical failure checks exist. One is to shift the window by a pixel or two and watch whether the estimate moves; a healthy estimate is insensitive to window placement. The other is to run all three methods and compare. Three methods resting on different assumptions agreeing on a value is grounds for trusting it, and a wide disagreement signals that one of the premises has broken.

## Summary and what comes next

Reading a position finer than a pixel means assuming something about the profile. When that assumption breaks differently at different grid phases, a systematic bias with a period of one pixel appears. Parabolic interpolation is worst near $w = 0.6$, at $0.042$ pixels. The moment method carries a bias of opposite phase that all but vanishes at $w = 0.8$. The area method, assuming only symmetry, is essentially unbiased with an uncropped window. Its scatter then grows nearly in proportion to the window length, and cropping the window trades that scatter for bias. The ranking of the three shifts with blur and noise. Matching the optics to a 10-90% edge width of 2 to 3 pixels makes the choice between the moment and area methods matter less, while the grid bias of parabolic interpolation remains.

The next post turns to geometric fitting. Having obtained the position of one edge point to subpixel precision, the problem becomes extracting a single parameter, such as a line or a circle, from hundreds of such points. What gets minimized changes the answer, and on short arcs the commonly used algebraic methods carry a pronounced bias. The trade between bias and scatter seen here returns there.

## References

- A. J. Tabatabai, O. R. Mitchell, "[Edge Location to Subpixel Values in Digital Imagery](https://doi.org/10.1109/TPAMI.1984.4767502)," _IEEE Trans. Pattern Anal. Mach. Intell._, PAMI-6(2), 188-201, 1984 (moment-based subpixel estimation).
- A. Trujillo-Pino et al., "[Accurate subpixel edge location based on partial area effect](https://doi.org/10.1016/j.imavis.2012.10.005)," _Image and Vision Computing_, 31(1), 72-90, 2013 (the area-based method).
- J. Canny, "[A Computational Approach to Edge Detection](https://doi.org/10.1109/TPAMI.1986.4767851)," _IEEE Trans. Pattern Anal. Mach. Intell._, PAMI-8(6), 679-698, 1986.
- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 10 (edge detection).
- J. S. Cybulski, J. Clements, M. Prakash, "[Foldscope: Origami-Based Paper Microscope](https://doi.org/10.1371/journal.pone.0098781)," _PLoS ONE_, 9(6), e98781, 2014 (source of Fig 8).
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025.

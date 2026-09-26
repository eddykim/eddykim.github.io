---
title: "Metrology Image Processing 5 — Geometric Fitting: What Exactly Gets Minimized"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-geometric-fitting/
page_id: imgproc-geometric-fitting
date: 2026-10-28 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, circle-fitting, ellipse-fitting, ransac, algebraic-distance, metrology]
description: "Minimizing the algebraic distance to extract one radius from hundreds of edge points underestimates it by 47% on a short arc. This traces where that bias comes from."
math: true
---

[Post 4](/en/posts/imgproc-edge-subpixel/) obtained the position of a single edge point to subpixel precision. What metrology finally needs, though, is not points but one number: the diameter of a hole, the centre of a lens, the angle of a wafer notch. How is that number extracted from hundreds of points?

The answer does not end at "fit a circle by least squares." **What** gets minimized matters, and the choice can be wrong by nearly half the radius on a short arc. This post follows that fork.

## 1. What gets minimized

Suppose a circle is to be fitted to points $(x_i, y_i)$. The quantity of interest is the actual distance from each point to the circle.

$$ d_i^{geo} = \sqrt{(x_i-a)^2 + (y_i-b)^2} - R $$

This is the geometric distance. The difficulty is that it is non-linear in the parameters $(a, b, R)$, so there is no closed-form solution and an iterative method is required.

A common detour exists. Expanding the circle equation into

$$ x^2 + y^2 + Dx + Ey + F = 0 $$

makes it linear in $(D, E, F)$. Substituting each point gives

$$ d_i^{alg} = x_i^2 + y_i^2 + D x_i + E y_i + F $$

which is called the algebraic distance, and minimizing its sum of squares takes one linear solve. The method Kåsa set out in 1976 remains the most widely used.

But the algebraic distance is not a distance. For a point at distance $d$ from the centre,

$$ d^{alg} = d^2 - R^2 = (d-R)(d+R) \approx 2R \cdot d^{geo} $$

so the same geometric error counts for more when $R$ is large. That $R$ is a free parameter is what makes this a problem. Shrinking the radius shrinks every term at once, and the minimization drags the radius down with it.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig1-algebraic-vs-geometric.png" alt="Algebraic and geometric objectives compared on the same point set" width="750">
_Fig 1. The algebraic objective has its minimum below the true radius_

The right panel shows this directly. For points on a 60-degree arc, the radius was held fixed while the centre was optimized, and both objectives were plotted. The algebraic objective bottoms out near $R \approx 36$, the geometric one near $R \approx 44$. The true value is 40.

On this single realization the geometric estimate is also off, by $+3.7$, because a short arc with only 60 points carries that much scatter. What matters is that the gap between the two minima opens in one direction only, and whether that is systematic is settled in the next section by repetition.

## 2. What happens on short arcs

If the bias grows as the arc shortens, how short is short enough to matter? With a true radius of 40, 120 points and radial noise of 0.4, the arc span was reduced from 360 to 10 degrees, 400 times each.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig2-arc-span-bias.png" alt="Radius bias and scatter of four methods against arc span" width="750">
_Fig 2. On short arcs the Kåsa method systematically underestimates the radius_

Over a full circle the four methods are indistinguishable, all biased by less than 0.01. They begin to separate near 120 degrees, and past that only the Kåsa method falls away.

| Arc | Kåsa bias | Taubin bias | Geometric bias |
|---|---|---|---|
| 360° | $+0.004$ | $+0.004$ | $+0.002$ |
| 90° | $-0.364$ | $+0.037$ | $+0.036$ |
| 45° | $-5.89$ | $+0.192$ | $+0.194$ |
| 30° | $-18.87$ | $+0.633$ | $+0.642$ |
| 20° | $-32.23$ | $+2.85$ | $+2.91$ |

On a 30-degree arc the Kåsa method reports 21.1 for a radius of 40. **That is 47% low, and neither less noise nor more points removes it.** This is what a systematic bias means.

Pratt and Taubin use the same algebraic distance yet are essentially unbiased. Both changed the constraint. Seeing the difference requires opening the circle equation one place further.

$$ A(x^2+y^2) + Bx + Cy + D = 0 $$

The form $x^2+y^2+Dx+Ey+F=0$ used in section 1 is this with $A$ fixed to 1. That fixing is precisely what Kåsa does, and it nails down the weight that scales with radius. It leaves open the route of shrinking all coefficients together to lower the objective. Pratt instead imposes $B^2+C^2-4AD=1$, and Taubin normalizes the algebraic distance by the magnitude of its gradient. Both constraints close that route. The letters $A$ through $D$ here belong to this section only and differ from the conic coefficients of section 4.

Worth noting is that Taubin and the geometric fit agree to three decimal places. For circles, a well-chosen algebraic method effectively reaches the geometric optimum, which means there is no reason to iterate. Where iteration is wanted anyway, the Levenberg-Marquardt method of [Optimization 3](/en/posts/optimization-levenberg-marquardt/) converges in three or four steps from an algebraic starting point.

The right panel shows scatter, and the order reverses there. At 30 degrees Kåsa scatters by 1.28 against 3.79 for the others. The biased estimator is the steadier one — the same shape as Post 4.

## 3. Where the unbiased methods lose

Since bias and scatter move oppositely, they have to be combined.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig3-total-error.png" alt="Total radius error of four methods against arc span" width="700">
_Fig 3. Below 20 degrees the unbiased methods are the ones that fail_

Down to 20 degrees Taubin and the geometric fit win: at 30 degrees their total error is 3.85 against 18.9 for Kåsa. Below 15 degrees the order flips, with Kåsa at 36.2 while Taubin reaches 90.9 and the geometric fit 104.

This reversal must not be read as "Kåsa is better on short arcs." On a 10-degree arc the geometric fit scatters by 2760. Measuring a radius of 40 with an error of 2760 is not a measurement. Kåsa's error looks smaller only because this method always collapses toward a smaller radius, so its error cannot exceed $R$ itself. Both have failed; they merely fail differently.

The rule that follows concerns design rather than method selection. Below 90 degrees the uncertainty of the radius estimate must be computed explicitly, and below 20 degrees the radius should not be measured at all. Needing only the centre helps somewhat, but even the centre error reaches 2.9 pixels at 30 degrees. That section 4.3 of the dissertation could use a full circle when setting the reference radius in back focal plane images was a requirement, not luck.

## 4. Ellipses demand one more thing

If a circle has three parameters, a general conic has six — five up to scale, matching the centre, two semi-axes and orientation of an ellipse exactly. Fitting a conic by least squares should then suffice. One problem stands in the way.

$$ Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0 $$

This equation does not describe only ellipses. The discriminant $B^2-4AC$ is negative for an ellipse, zero for a parabola and positive for a hyperbola. An unconstrained least-squares solution gives no say in which of the three appears.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig4-ellipse-constraint.png" alt="How often a conic fit returns an ellipse, with and without the constraint" width="750">
_Fig 4. Without the constraint a conic fit does not guarantee an ellipse_

Eighty points on a 120-degree arc were fitted 300 times at increasing noise. At low noise the unconstrained fit always returns an ellipse, but the rate falls to 66% at a noise of 2 and to 15% at 5. The rest are hyperbolas or parabolas. The right panel shows one such case: the same points yield an ellipse under one method and a two-branch hyperbola under the other.

The method proposed by Fitzgibbon and colleagues in 1999 imposes $4AC - B^2 = 1$. The left side is the discriminant with its sign flipped, so setting it to 1 is the same as demanding an ellipse. Absorbing the constraint into the normalization reduces the problem to a single generalized eigenproblem, and **an ellipse comes out for any input whatsoever.** That is why the green curve sits at 100% across all noise levels.

The property cuts both ways. An ellipse is returned even when the points are not on one, so the fact that the result is an ellipse is no evidence that the fit is sound. The residual has to be checked separately.

The implementation follows the numerically stable form of Halíř and Flusser, which splits the design matrix into quadratic and linear parts and reduces a $6\times6$ eigenproblem to $3\times3$.

```python
# core of ellipse_fit.py (full code: _code/imgproc-geometric-fitting/)
D1 = np.column_stack([u * u, u * v, v * v])      # quadratic terms
D2 = np.column_stack([u, v, np.ones_like(u)])    # linear terms
S1, S2, S3 = D1.T @ D1, D1.T @ D2, D2.T @ D2
T = -np.linalg.solve(S3, S2.T)
M = S1 + S2 @ T
M = np.array([M[2] / 2.0, -M[1], M[0] / 2.0])    # constraint folded into the left side
evals, evecs = np.linalg.eig(M)
cond = 4.0 * evecs[0] * evecs[2] - evecs[1] ** 2  # ellipse condition
```

## 5. Outliers drag the least-squares fit

Every point has so far been assumed to come from the same circle. Real images do not oblige. Edges from other structures creep in, and as section 7 of Post 4 showed, a window holding two edges produces a stray point.

Least squares sums squared residuals, so a few badly misplaced points outweigh hundreds of good ones. The standard way to screen them out is random sample consensus (RANSAC).

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig5-ransac.png" alt="Least squares compared with RANSAC on a point set containing outliers" width="750">
_Fig 5. Outliers drag the least-squares fit while RANSAC holds_

The radius error was measured as the outlier fraction rose from 0 to 40%. Least squares is already at $+0.47$ with 5% outliers and reaches $+3.66$ at 40%. The RANSAC result stays at $-0.006$ even at 40%.

The logic is simple. Instead of using every point, draw the three points that minimally determine a circle, build the model, and count how many points fall within a set distance of it. Repeat a few hundred times, keep the model with the most agreement, then refit precisely using only those points.

There is one knob: the distance threshold. It is set to two or three times the scatter of the inliers, which is exactly what Post 4 measured for subpixel edge points. Too tight and inliers are discarded; too loose and outliers are embraced. In the example on the left, 30 outliers were planted and RANSAC classified 93 points as inliers — 120 minus 93 leaves 27, close to the true count.

## 6. On a real image

Whether the synthetic result carries over was checked on a real image. Canny edges were extracted from a widely used image of coins, and the largest connected component gave 226 edge points.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig6-real-image.png" alt="Circle fits to a coin contour and to a 45-degree arc of it" width="750">
_Fig 6. The same bias appears on a real image_

Using the whole contour gives a radius of 31.55 pixels, with the four methods agreeing to two decimal places. Keeping only the 32 points spanning 45 degrees gives 26.80 for Kåsa, 33.83 for Taubin and 33.85 for the geometric fit.

One caution belongs here. The true value is unknown in this case. Taking the 31.55 from the full contour as a reference makes Kåsa 15% low and Taubin 7% high, but the discrepancy on the Taubin side looks like the coin rim not being perfectly circular rather than a bias, since the value shifts with the choice of arc. This is why the synthetic experiment used a known truth and 400 repetitions. A single real image cannot separate a systematic bias from the shape error of that particular specimen.

## Summary and what comes next

Extracting one number from many points means choosing what to minimize. The algebraic distance is cheap to compute but is not a distance, so it systematically shrinks the radius as the arc gets shorter — by 47% at 30 degrees. Pratt and Taubin remove that bias by changing the constraint, and for circles they effectively reach the geometric optimum, so iteration is not required. For ellipses the constraint serves a different purpose: guaranteeing the kind of solution rather than removing a bias. And whichever method is used, a radius from an arc shorter than 20 degrees is not a measurement.

The next post turns to phase. All five posts so far have been about reading intensity. In interferometric and moiré images the information sits in the phase rather than the intensity, and phase is never measured directly. The Hilbert transform used in Post 3 gets treated head-on, along with how phase shifting cancels background and modulation depth, and the unwrapping problem that arises when the phase is recovered.

## References

- I. Kåsa, "[A circle fitting procedure and its error analysis](https://doi.org/10.1109/TIM.1976.6312298)," _IEEE Trans. Instrum. Meas._, IM-25(1), 8-14, 1976.
- V. Pratt, "[Direct least-squares fitting of algebraic surfaces](https://doi.org/10.1145/37402.37420)," _ACM SIGGRAPH Computer Graphics_, 21(4), 145-152, 1987.
- N. Chernov, C. Lesort, "[Least squares fitting of circles](https://doi.org/10.1007/s10851-005-0482-8)," _Journal of Mathematical Imaging and Vision_, 23(3), 239-252, 2005 (bias analysis of Kåsa, Pratt and Taubin, and stable implementations).
- A. Fitzgibbon, M. Pilu, R. B. Fisher, "[Direct least square fitting of ellipses](https://doi.org/10.1109/34.765658)," _IEEE Trans. Pattern Anal. Mach. Intell._, 21(5), 476-480, 1999.
- R. Halíř, J. Flusser, "Numerically stable direct least squares fitting of ellipses," _Proc. WSCG_, 125-132, 1998.
- M. A. Fischler, R. C. Bolles, "[Random sample consensus](https://doi.org/10.1145/358669.358692)," _Communications of the ACM_, 24(6), 381-395, 1981.
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the coin test image in Fig 6).
- Y. Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," Ph.D. dissertation, Seoul National University, 2025, sec. 4.3 (setting the reference radius in back focal plane images).

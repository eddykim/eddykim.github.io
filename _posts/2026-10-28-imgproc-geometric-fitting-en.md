---
title: "Metrology Image Processing 5 — Can a Radius Be Trusted When Only Part of the Circle Is Visible?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-geometric-fitting/
page_id: imgproc-geometric-fitting
date: 2026-10-28 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, circle-fitting, ellipse-fitting, ransac, algebraic-distance, metrology]
description: "Minimizing the algebraic distance to extract one radius from hundreds of edge points underestimates it by 47% on a 30-degree arc. This traces where that bias comes from."
math: true
---

[Post 4](/en/posts/imgproc-edge-subpixel/) obtained the position of a single edge point to subpixel precision. What metrology finally needs, though, is not points but one number: the diameter of a hole, the centre of a lens, the angle of a wafer notch. How is that number extracted from hundreds of points?

The answer does not end at "fit a circle by least squares." **What** gets minimized matters, and on a short arc the choice can make the radius wrong by nearly half. This post follows that fork.

## 1. What should be minimized?

Suppose a circle is to be fitted to points $(x_i, y_i)$. The quantity of interest is the actual distance from each point to the circle.

$$ d_i^{geo} = \sqrt{(x_i-a)^2 + (y_i-b)^2} - R $$

This is the geometric distance. The difficulty is that it is non-linear in the parameters $(a, b, R)$, so there is no closed-form solution and an iterative method is required.

A common detour exists. Expanding the circle equation into

$$ x^2 + y^2 + Dx + Ey + F = 0 $$

makes it linear in $(D, E, F)$. Substituting each point gives

$$ d_i^{alg} = x_i^2 + y_i^2 + D x_i + E y_i + F $$

which is called the algebraic distance, and minimizing its sum of squares takes one linear solve. The method Kåsa set out in 1976 is still widely used. Because the same method already appears in a 1972 paper by Delogne, it is also called the Delogne–Kåsa method.

But the algebraic distance is not a distance. For a point at distance $d$ from the centre,

$$ d^{alg} = d^2 - R^2 = (d-R)(d+R) \approx 2R \cdot d^{geo} $$

so the same geometric error counts for more when $R$ is large. That $R$ is a free parameter is what makes this a problem. Shrinking the radius shrinks every term at once, and the minimization drags the radius down with it. In other words, the algebraic distance is a balance that penalizes large circles.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig1-algebraic-vs-geometric.png" alt="Algebraic and geometric objectives compared on the same point set" width="750">
_Fig 1. Circle fits to 60 edge points on a 60° arc (radial noise 0.8, 2% of the radius). Left: the circles from the Kåsa method (red) and the geometric fit (blue), overlaid on the true circle (dotted). Right: the two objectives with the radius held fixed and only the centre optimized, each scaled so that its minimum is 1. The minimum of the algebraic objective lies below the true radius of 40._

The right panel shows this directly. Sixty points on a 60-degree arc were given radial noise of 0.8, and both objectives were plotted with the radius held fixed while the centre was optimized. The algebraic objective bottoms out near $R \approx 36$, the geometric one near $R \approx 44$. The true value is 40. The left panel draws the resulting circles. The Kåsa circle (red, $R=36.4$) passes through the points, yet the whole circle has shrunk inward.

On this single point set the geometric estimate is also off, by $+3.7$, because a short arc with only 60 points carries that much scatter. What matters is that the gap between the two minima opens in one direction only. Whether that is systematic is settled in the next section by repetition.

## 2. What happens on short arcs?

In real metrology, seeing only part of a circle is common. Fig 2 shows the edges of two silicon wafers. The upper wafer has a flat, a straight cut along the edge; the lower one has a notch, a V-shaped groove. To say in which direction the notch lies from the centre, the centre has to be known first, yet a single photograph captures only a short arc of a large circle.

<img src="/assets/img/posts/imgproc-geometric-fitting/ext-wafer-notch-flat.jpg" alt="The flat and notch on silicon wafer edges, photographed next to a match" width="600">
_Fig 2. The flat (top, 6-inch wafer) and notch (bottom, 200 mm wafer) of silicon wafers. Compared with the match, the edge in the photograph is only a very short part of the circumference. (Source: Kuebi = Armin Kübelbeck, [Notch and Flat 01](https://commons.wikimedia.org/wiki/File:Notch_and_Flat_01.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

Why a short arc is a problem can be guessed from a bow. Call the largest height of the bow above its string the sagitta. For an arc subtending an angle $\theta$, it is

$$ s = R\left(1-\cos\frac{\theta}{2}\right) \approx \frac{R\theta^2}{8} $$

with $\theta$ in radians. The information about the radius sits essentially in this sagitta. For the same string length, a small sagitta means a large circle and a large one a small circle. The sagitta, however, falls with the square of the arc angle. On a circle of radius 40 it is 1.36 for a 30-degree arc, 0.61 for 20 degrees and 0.34 for 15 degrees. In other words, as the arc shortens, the signal that carries the radius sinks to the size of the noise.

How short, then, is short enough to matter? With a true radius of 40, 120 points and radial noise of 0.4 (1% of the radius), the arc span was reduced from 360 to 10 degrees, and the radius bias and scatter were measured over 400 repetitions at each angle.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig2-arc-span-bias.png" alt="Radius bias and scatter of four methods against arc span" width="750">
_Fig 3. Radius bias (left) and scatter (right) of four methods, fitted 400 times at each arc span with a true radius of 40, 120 points and noise of 0.4 (1% of the radius). The left vertical axis is linear only between −1 and 1 and logarithmic outside. Pratt (orange) almost coincides with Taubin (green). On short arcs only the Kåsa method systematically underestimates the radius, while its scatter is, conversely, the smallest._

Over a full circle the four methods are indistinguishable, all biased by less than 0.01. They begin to separate near 120 degrees, and past that only the Kåsa method falls away.

| Arc | Kåsa bias | Taubin bias | Geometric bias |
|---|---|---|---|
| 360° | $+0.004$ | $+0.004$ | $+0.002$ |
| 90° | $-0.364$ | $+0.037$ | $+0.036$ |
| 45° | $-5.89$ | $+0.192$ | $+0.194$ |
| 30° | $-18.87$ | $+0.633$ | $+0.642$ |
| 20° | $-32.23$ | $+2.85$ | $+2.91$ |

On a 30-degree arc the Kåsa method reports 21.1 for a radius of 40. **That is 47% low, and more points do not remove it.** More points shrink the scatter but do not pull the mean back to the truth. That is what a systematic bias means. The bias does, however, grow with the square of the noise: at the same 30 degrees, halving the noise brings it to $-7.3$, and cutting it tenfold to $-0.3$.

Pratt and Taubin use the same algebraic distance yet are essentially unbiased. Both changed the constraint. Seeing the difference requires writing the circle equation with one more coefficient.

$$ A(x^2+y^2) + Bx + Cy + D = 0 $$

The form $x^2+y^2+Dx+Ey+F=0$ used in section 1 is this with $A$ fixed to 1. That fixing is precisely what Kåsa does, and it nails the radius-dependent weight into the objective. With $A=1$ the algebraic distance is $d^2-R^2 \approx 2R\,d^{geo}$, so the objective keeps a factor of $4R^2$. Pratt instead imposes $B^2+C^2-4AD=1$; for a circle $B^2+C^2-4AD = 4A^2R^2$, so this amounts to dividing the objective by $4R^2$. Taubin divides the sum of squared algebraic distances by the mean squared gradient magnitude, which on a circle is again close to $4R^2$. Both constraints strip out the radius-dependent weight. The letters $A$ through $D$ here belong to this section only and differ from the conic coefficients of section 4.

Worth noting is that Taubin and the geometric fit differ in bias by less than 0.01 even at 30 degrees. For circles, a well-chosen algebraic method effectively reaches the geometric optimum, so in practice iteration is usually unnecessary. Pratt and Taubin do keep a small bias proportional to the square of the noise, and with very many points the difference from the geometric fit shows up again (Al-Sharadqah and Chernov, 2009). Where iteration is wanted, the Levenberg-Marquardt method of [Optimization 3](/en/posts/optimization-levenberg-marquardt/) will do. From an algebraic starting point it converges in three or four steps on a full circle, and in about ten on a 30-degree arc.

The right panel shows scatter, and the order reverses there. At 30 degrees Kåsa scatters by 1.28 against 3.79 for the others. The less biased method can lose on scatter, as observed in Post 4, and the same happens here.

## 3. When do the less biased methods fail worse?

Since bias and scatter move oppositely, they have to be combined. As in Post 4, the total error is taken as $\sqrt{\text{bias}^2 + \text{scatter}^2}$, which equals the root-mean-square error (RMSE) about the true value.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig3-total-error.png" alt="Total radius error of four methods against arc span" width="700">
_Fig 4. Total radius error of the four methods under the conditions of Fig 3. Down to 20 degrees Taubin, Pratt and the geometric fit have a smaller total error than the Kåsa method; from 15 degrees the order flips. The shaded band is where no method is usable when the noise is 1% of the radius._

Down to 20 degrees Taubin and the geometric fit win: at 30 degrees their total error is 3.85 against 18.9 for Kåsa. Even so, at 20 degrees the Taubin total error is 11.6, about 29% of the radius, so winning does not make the result usable. From 15 degrees down the order flips, with Kåsa at 36.2 while Taubin reaches 90.9 and the geometric fit 104.

These values below 15 degrees are not numbers to take at face value. Changing the random seed five ways at 15 degrees leaves Kåsa at 36.1 to 36.2, while Taubin ranges from 50 to 135 and the geometric fit from 53 to 162. Chernov (2011) proved that under the usual noise assumptions the estimates of a circle's centre and radius have no finite mean or variance. The property stays hidden on long arcs; on short arcs it shows up as an average that does not settle even over 400 repetitions.

This reversal must not be read as "Kåsa is better on short arcs." On a 10-degree arc the geometric fit scatters by 2760. Measuring a radius of 40 with an error of 2760 is not a measurement. Kåsa's error looks smaller only because this method always collapses toward a smaller radius, so its error cannot exceed $R$ itself. Both have failed; they merely fail differently.

The rule that follows concerns the design of the measurement rather than the choice of method. With noise at 1% of the radius and 120 points, the uncertainty of the radius estimate must be computed explicitly below 90 degrees, and below 20 degrees the radius should not be measured at all. The thresholds move with the noise. At 0.1% of the radius, even a 15-degree arc gives the radius to within 4%.

The sagitta of section 2 ties the two thresholds together. When the ratio $s/\sigma$ of sagitta to noise is the same, the errors are similar even for different arc angles and noise levels. With 120 points, at 30 degrees with 1% noise ($s/\sigma = 3.4$) and at 10 degrees with 0.1% noise ($s/\sigma = 3.8$), the Taubin total error is 9.6% and 8.4% of the radius, and the Kåsa bias is $-18.9$ and $-17.1$. With 120 points and $s/\sigma$ above 3, the relative radius error of Taubin can be estimated roughly as $0.32\,\sigma/s$. Below about 1.5 this estimate breaks down and the fit diverges. In other words, before looking at the arc angle, work out how many times the noise the sagitta is.

Needing only the centre helps somewhat, but even the centre error averages 3 pixels at 30 degrees. Section 4.3 of the author's dissertation did not set the reference radius in back focal plane images by fitting a circle. It used the Brewster angle instead, located where $\tan^2\Psi$ is smallest along a diametral profile. Here $\Psi$ is the ellipsometric angle that expresses the ratio of the p and s reflection amplitudes. The reference does not depend on arc length at all.

## 4. Why might a conic fit not be an ellipse?

A round hole or a lens rim photographed at an angle appears as an ellipse in the image. If a circle has three parameters, a general conic has six. Five remain once the scale is removed, matching the centre, two semi-axes and orientation of an ellipse exactly. Fitting a conic by least squares should then suffice, but one problem stands in the way.

$$ Ax^2 + Bxy + Cy^2 + Dx + Ey + F = 0 $$

This equation does not describe only ellipses. The discriminant $B^2-4AC$ is negative for an ellipse, zero for a parabola and positive for a hyperbola. All three are sections of a cone cut by a plane, and changing only the tilt of the cutting plane turns an ellipse through a parabola into a hyperbola (Fig 5). The coefficients behave the same way: a little noise in them can flip the sign of the discriminant. An unconstrained least-squares solution gives no control over which of the three appears.

<img src="/assets/img/posts/imgproc-geometric-fitting/ext-conic-sections.jpg" alt="A double cone cut by planes at four different tilts" width="600">
_Fig 5. Sections of a double cone cut by a plane. A plane perpendicular to the axis gives a circle (1), a slight tilt an ellipse (2), a plane parallel to the side of the cone a parabola (3), and a steeper plane cutting both cones a hyperbola (4). The parabola is merely the boundary case between ellipse and hyperbola. (Source: JensVyff, [TypesOfConicSections](https://commons.wikimedia.org/wiki/File:TypesOfConicSections.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig4-ellipse-constraint.png" alt="How often a conic fit returns an ellipse, with and without the constraint" width="750">
_Fig 6. The fraction of ellipses returned when conics were fitted 300 times to 80 points on a 120° arc at increasing noise (left), and an example at noise 3 in which the unconstrained fit returns a hyperbola (right). The unconstrained fit starts to lose the guarantee at noise 1 and drops to 66% at noise 2, while the Fitzgibbon constraint yields an ellipse at every noise level._

Eighty points on a 120-degree arc were fitted 300 times at increasing noise. At low noise the unconstrained fit always returns an ellipse, but the rate falls to 66% at a noise of 2 and to 15% at 5 (it rises again to 22% at noise 8, but the guarantee is still absent). The rest are hyperbolas or parabolas. The right panel shows one such case: for the same points the Fitzgibbon fit is an ellipse and the unconstrained one a two-branch hyperbola.

The method Fitzgibbon and colleagues published in 1996, and in journal form in 1999, imposes $4AC - B^2 = 1$. The left side is the discriminant with its sign flipped, so setting it to 1 is the same as demanding an ellipse. Introducing this constraint through a Lagrange multiplier reduces the problem to a single generalized eigenproblem, and **an ellipse comes out for any input whatsoever.** That is why the green curve sits at 100% across all noise levels.

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

## 5. How are stray points screened out?

Every point has so far been assumed to come from the same circle. Real images do not oblige. Edges from other structures creep in, and as section 7 of Post 4 showed, a window holding two edges produces a stray point.

Least squares sums squared residuals, so a few badly misplaced points outweigh hundreds of good ones. Such points are called outliers, and the standard way to screen them out is random sample consensus (RANSAC).

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig5-ransac.png" alt="Least squares compared with RANSAC on a point set containing outliers" width="750">
_Fig 7. Left: an example in which 30 of 120 points on a circle of radius 40 were turned into outliers, with the least-squares circle fitted to all points (red) and the RANSAC circle (green). Right: the mean radius error over 120 trials as the outlier fraction rises from 0 to 40%. Least squares is dragged into a larger radius by the outliers, while RANSAC holds the true value even at 40%._

The radius error was measured as the outlier fraction rose from 0 to 40%. Least squares is already at $+0.47$ with 5% outliers and reaches $+3.66$ at 40%. The RANSAC result stays at $-0.006$ even at 40%.

The logic is simple. Instead of using every point, draw the three points that minimally determine a circle, build the model, and count how many points, the inliers, fall within a set distance of it. Repeat a few hundred times, keep the model with the most agreement, then refit precisely using only those points.

The knob that actually has to be chosen is the distance threshold. The number of iterations is a parameter too, but once the inlier fraction $w$ and the desired success probability $p$ are fixed it follows from $N = \log(1-p)/\log(1-w^3)$. Even with 60% inliers, 19 draws reach $p=0.99$, so a few hundred is generous. A threshold of about twice the inlier scatter is the starting point; for Gaussian noise, 1.96 times keeps 95% of the inliers. That scatter is exactly what Post 4 measured for subpixel edge points. Too tight and inliers are discarded; too loose and outliers are embraced.

The computation behind Fig 7 sets a generous threshold of 1.5, or 3.75 times the noise of 0.4, because most outliers lie far from the circle. At this threshold practically no Gaussian inlier is discarded. In the example on the left, 30 outliers were planted and RANSAC classified 93 points as inliers. That is all 90 true inliers plus 3 outliers that happened to land within the threshold.

## 6. Does the same bias appear in a real image?

Whether the synthetic result carries over was checked on a real image. A Canny detector extracted edges from a widely used image of coins, and the largest connected component gave 226 edge points.

<img src="/assets/img/posts/imgproc-geometric-fitting/en/fig6-real-image.png" alt="Circle fits to a coin contour and to a 45-degree arc of it" width="750">
_Fig 8. The contour of one coin (blue points) found by a Canny detector in the scikit-image coins test image. Left: the geometric fit to all 226 points (yellow). Right: circles from three methods fitted to only the 32 points of a 45-degree stretch at the lower right. The Taubin circle (green, dashed) almost coincides with the geometric fit (blue), and only the Kåsa circle (red) is visibly smaller._

Using the whole contour gives a radius of 31.55 pixels, with the four methods agreeing to within 0.02 pixels. Keeping only the 32 points spanning 45 degrees gives 26.80 for Kåsa, 33.83 for Taubin and 33.85 for the geometric fit.

One caution belongs here. The true value is unknown in this case. Taking the 31.55 from the full contour as a reference makes Kåsa 15% low and Taubin 7% high. The Taubin discrepancy, however, looks like the scatter of fitting a short arc with 32 points rather than a bias. Cutting eight different 45-degree pieces from the same contour, Taubin ranges from 23.5 to 45.6, while Kåsa comes out below 31.55 on all eight. This is why the synthetic experiment used a known truth and 400 repetitions. A single real image cannot separate a systematic bias from the shape error of that particular specimen.

## Summary and what comes next

The practical prescription first. For a circle, one pass of Taubin or Pratt is enough. On a short arc, work out how many times the noise the sagitta is before choosing a method. Fit ellipses with the Fitzgibbon constraint and check the residual separately. When outliers are present, select the inliers with RANSAC and refit.

Extracting one number from many points means choosing what to minimize. The algebraic distance is cheap to compute but is not a distance, so it systematically shrinks the radius as the arc gets shorter: by 47% at 30 degrees when the noise is 1% of the radius. Pratt and Taubin greatly reduce that bias by changing the constraint, and for circles they effectively reach the geometric optimum, so iteration is usually not required. For ellipses the constraint serves a different purpose, guaranteeing the kind of solution rather than removing a bias. And with 120 points and noise of about 1% of the radius, a radius from an arc shorter than 20 degrees is not a measurement, whichever method is used.

The next post turns to phase. All five posts so far have been about reading intensity. In interferometric and moiré images the information sits in the phase rather than the intensity, and phase is never measured directly. The Hilbert transform gets treated head-on, extending the frequency-domain filtering of Post 3, along with how phase shifting cancels background and modulation depth, and the unwrapping problem that arises when the phase is recovered.

## References

- I. Kåsa, "[A circle fitting procedure and its error analysis](https://doi.org/10.1109/TIM.1976.6312298)," _IEEE Trans. Instrum. Meas._, IM-25(1), 8-14, 1976.
- V. Pratt, "[Direct least-squares fitting of algebraic surfaces](https://doi.org/10.1145/37402.37420)," _ACM SIGGRAPH Computer Graphics_, 21(4), 145-152, 1987.
- G. Taubin, "[Estimation of planar curves, surfaces, and nonplanar space curves defined by implicit equations with applications to edge and range image segmentation](https://doi.org/10.1109/34.103273)," _IEEE Trans. Pattern Anal. Mach. Intell._, 13(11), 1115-1138, 1991.
- N. Chernov, C. Lesort, "[Least squares fitting of circles](https://doi.org/10.1007/s10851-005-0482-8)," _Journal of Mathematical Imaging and Vision_, 23(3), 239-252, 2005 (bias analysis of Kåsa, Pratt and Taubin, stable implementations, and the Delogne–Kåsa attribution).
- A. Al-Sharadqah, N. Chernov, "[Error analysis for circle fitting algorithms](https://doi.org/10.1214/09-EJS419)," _Electronic Journal of Statistics_, 3, 886-911, 2009 (the bias of order noise squared that remains in Pratt, Taubin and the geometric fit).
- N. Chernov, "[Fitting circles to scattered data: parameter estimates have no moments](https://doi.org/10.1007/s00184-009-0283-y)," _Metrika_, 73(3), 373-384, 2011.
- A. Fitzgibbon, M. Pilu, R. B. Fisher, "[Direct least square fitting of ellipses](https://doi.org/10.1109/34.765658)," _IEEE Trans. Pattern Anal. Mach. Intell._, 21(5), 476-480, 1999.
- R. Halíř, J. Flusser, "Numerically stable direct least squares fitting of ellipses," _Proc. WSCG_, 125-132, 1998.
- M. A. Fischler, R. C. Bolles, "[Random sample consensus: a paradigm for model fitting with applications to image analysis and automated cartography](https://doi.org/10.1145/358669.358692)," _Communications of the ACM_, 24(6), 381-395, 1981.
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the coin test image in Fig 8).
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025.

---
date: 2026-10-24 20:00:00 +0900
layout: post
title: "Ellipsometer Instrumentation 4 — Can You Trust the Angle on the Polarizer Dial?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-calibration/
page_id: ellipsometer-calibration
categories: [Optics, Instrumentation]
tags: [ellipsometry, calibration, rotating-analyzer, regression, instrumentation]
description: The angle a dial reports is not the true angle. This post follows three ways of measuring that difference, and what each of them cannot do.
math: true
---

[Post 1](/en/posts/ellipsometer-oblique-rotating-element/) sorted rotating-element ellipsometers by what they cannot measure. [Post 2](/en/posts/ellipsometer-components/) followed where each component stops working across a spectral band. [Post 3](/en/posts/ellipsometer-dual-rotating-compensator-mme/) asked why one modulating element is not enough.

All three rested on an assumption. When the polarizer dial reads $30°$, the transmission axis is taken to sit $30°$ from the plane of incidence. That assumption is false. No matter how precisely a divided circle is cut, there is no guarantee that the optical transmission axis and the mechanical zero coincide. Nor does a detector return a signal strictly proportional to the light that enters it.

The question here follows from that. **How do you measure the difference between what the dial reads and the true angle?** Measuring that difference and undoing it is called calibration. The instrument built to measure a sample must now measure itself, and the measuring tool is already misaligned. That is the one difficulty.

Three eras answered it differently, and this post follows them in order. Yet all three use one idea. Set up a relation that must hold if the instrument were ideal, then **gauge the defect by how far that relation is broken.** It is the general form of what Section 9 of Post 3 called the surplus observables watching over the instrument. The last section asks when the calibration values obtained this way have to be measured again.

## 1. Where do the errors enter?

Before looking at methods, it is worth counting what needs correcting. Azzam and Bashara divide the paths by which $\rho = \tan\Psi\,e^{i\Delta}$ picks up systematic error into seven.

Azimuth errors; optical imperfections in components and cell windows (the windows of the vacuum or reaction cells used for real-time measurement); beam deviation; stray beams; polarization and collimation errors in the source beam; polarization-dependent sensitivity of the detector; and residual mechanical errors.

Azimuth errors split further by origin. Residual calibration error of the divided circle stays constant once measured, while sample remounting error changes every time the sample is exchanged or replaced. The first is a property of the instrument; the second is generated afresh with each measurement. This distinction answers the question of Section 10, when to recalibrate.

Component imperfections are worse. Polarizing prisms have finite extinction ratios, cell windows carry stress birefringence, compensators scatter from defects. Azzam and Bashara write that these are **unavoidable**. The component specifications surveyed in Post 2 come back here as an error list.

The detector side does not end in one line either. There is dark signal, which has to be measured separately with a mechanical shutter blocking the beam path and subtracted. There is nonlinear response, in which output is not proportional to input intensity. There is polarization-dependent sensitivity, produced when residual stress in a photomultiplier tube (PMT) window acts photoelastically. The photoelastic effect that served as a modulation tool in Post 2 appears here as a defect. Array detectors add image persistence, in which signal survives readout.

The list also shifts with the configuration. A rotating-compensator ellipsometer (RCE) must first fix the compensator retardance at each wavelength. This is measured in a straight-through arrangement, with the sample removed and only the compensator between polarizer and analyzer. Because of the $\delta \propto 1/\lambda$ relation from Post 2, the retardance of an MgF$_2$ compensator rises nearly linearly with photon energy. Parameterizing that dispersion lets it be recovered in the general arrangement as well. Phase-modulation ellipsometry (PME) must control the retardance, set that value at each wavelength, and control the temperature as well, plus a higher-harmonic correction. The higher-order Bessel terms from Post 1 and the temperature sensitivity of the photoelastic modulator from Post 2 meet here.

A list this long cannot be worked through one item at a time. Practical calibration takes a different route.

## 2. What does calibration have to undo?

If the seven branches cannot be handled one by one, what should be handled instead? Putting the errors into the signal shows the answer. Written out with the errors included, the detected intensity of a rotating-analyzer ellipsometer (RAE) reads

$$I_t(t) = I_0\big[1 + \eta\,\alpha\cos 2(\omega t - A_s) + \eta\,\beta\sin 2(\omega t - A_s)\big]$$

$P_s$ and $A_s$ are the offsets of the polarizer and analyzer zero positions from the plane of incidence. When the dial reads $P$, the true angle is $P - P_s$, and likewise for the analyzer. The coefficient $\eta$ carries the detector non-ideality: it states how far the ac component is attenuated relative to the dc component. The nonlinear response and the polarization-dependent sensitivity of the previous section both enter here at once. That $\eta$ is not the coefficient of a single cause becomes a problem later.

How, then, are the measured $$(\alpha', \beta')$$ related to the true $(\alpha, \beta)$? Regrouping the expression above by the $\cos 2\omega t$ and $\sin 2\omega t$ coefficients gives

$$\begin{bmatrix}\alpha'\\ \beta'\end{bmatrix} = \eta\begin{bmatrix}\cos 2A_s & -\sin 2A_s\\ \sin 2A_s & \cos 2A_s\end{bmatrix}\begin{bmatrix}\alpha\\ \beta\end{bmatrix}$$

This is a coordinate rotation through $2A_s$ in the $(\alpha, \beta)$ plane, multiplied by a scale factor $\eta$. Knowing $A_s$ and $\eta$ is therefore enough: invert the matrix and the true values come back. The error list ran to seven branches, but what must be undone reduces to **one rotation and one scale factor**. That reduction is what makes calibration tractable.

In multichannel instruments, undoing $A_s$ mathematically matters most. A photodiode array is read out pixel by pixel, and the analyzer keeps turning during readout. So **$A_s$ differs from pixel to pixel.** No mechanical adjustment can fix that, so it has to be handled by calibration from the start. Section 8 of [Post 1](/en/posts/ellipsometer-oblique-rotating-element/) noted that RAE is achromatic and therefore measures as many wavelengths simultaneously as the array has pixels. This is the price attached to that strength.

The three errors leave different marks on the signal. A clock makes the difference easy to see. $A_s$ is a dial fitted back on at an angle, so every reading is off by the same amount. $\eta$ is a shortened hand: the direction is right but the length is misread. $P_s$ leaves the clock intact but shifts the moment at which it is read. The right panel of Figure 1 shows exactly these three.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig1-error-anatomy.png" alt="What the three errors do to the measured signal" width="780">
_Figure 1. Left, the detected intensity over half an analyzer revolution (0°–180°, one signal period). Right, the locus traced by $$(\alpha', \beta')$$ as the polarizer is scanned. The errors are exaggerated for visibility, $$P_s = 4°$$, $$A_s = 8°$$, $\eta = 0.7$ (the text works with errors around $0.1°$). The right panel uses a sample with large $\lvert\cos\Delta\rvert$ ($\Psi = 35°$, $\Delta = 50°$) so the ellipse is not flattened. In the waveform the three look alike; in the plane of the ellipse they separate into rotation ($$A_s$$), shrinkage ($\eta$), and sliding along the ellipse ($$P_s$$)._

In the waveform the three look like similar distortions. In the $$(\alpha', \beta')$$ plane they separate cleanly. Scanning the polarizer from $-90°$ to $+90°$ traces a closed ellipse with semi-axes $1$ and $\lvert\cos\Delta\rvert$; $\Psi$ does not enter. Substituting $u = \tan(P - P_s) = \tan\Psi\,\tan\theta$ gives $\alpha = \cos 2\theta$ and $\beta = \cos\Delta\,\sin 2\theta$, and $\Psi$ cancels. $A_s$ rotates this ellipse about the origin, $\eta$ shrinks it whole, and $P_s$ **leaves the ellipse untouched and only slides the operating point along it.**

Because the three act differently, they can be measured separately. Residual calibration uses that separation.

## 3. How wrong is the answer without calibration?

Before the methods, it helps to fix the size of the problem. Synthesize a signal with one error at a time, then feed the dial readings straight into the ideal formulas and read off $\Psi$ and $\Delta$. The polarizer sits at $45°$. The signal is synthesized from products of component Mueller matrices alone, and the implementation was checked separately against the closed-form expressions of Fujiwara's Eqs. 4.59–4.60.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig2-uncalibrated-error.png" alt="Uncalibrated readings fed straight into the ideal formulas" width="780">
_Figure 2. Errors propagated into $\Psi$ (left) and $\Delta$ (right) for a sample with $\Psi = 16.8°$, sweeping $\Delta$, with one error injected at a time. The polarizer is at $45°$, and the dotted line marks the $0.01°$ level. $$P_s$$ reaches only $\Psi$; $$A_s$$ disturbs $\Delta$ near $\Delta \approx 90°$; $\eta$ disturbs $\Delta$ more and more as $\Delta$ approaches $0°$ or $180°$._

The three errors go to different places. A polarizer offset $P_s = 0.1°$ shifts $\Psi$ by $0.055°$ (by $0.064°$ with the polarizer at $30°$ or $60°$) but does **not** propagate into $\Delta$ at all. In the ideal inversion $$\cos\Delta = \beta' / \sqrt{1 - \alpha'^2}$$, and the polarizer angle cancels between numerator and denominator. The $10^{-13}$ degrees left in the computation are floating-point noise.

An analyzer offset $A_s = 0.1°$ does the opposite: it barely touches $\Psi$, and near $\Delta \approx 90°$ it shifts $\Delta$ by $0.30°$, amplified threefold over the injected error. This $\Delta$ error, however, scales with $\sin\Delta$ and shrinks toward $0°$ and $180°$ ($0.016°$ at $\Delta = 3°$), while the $\Psi$ error grows to nearly $0.1°$ instead.

The detector coefficient $\eta$ leaves $\Delta$ untouched at $\Delta = 90°$ and degrades sharply as $\Delta$ approaches $0°$ or $180°$. An attenuation of one part in a thousand, $\eta = 0.999$, produces a $2.5°$ error on a sample with $\Delta = 3°$. This is the region where Post 1 put the RAE amplification factor at $1/\lvert\sin\Delta\rvert$.

A real instrument departs much further. In Table 1 of Johs's 1993 paper, $\eta$ is $0.9687$, about 3% from unity. Left uncorrected, an $\eta$ that far off shifts $\Psi$ by more than $1°$ even on a sample like SiO$_2$/Si ($\Psi = 16.8°$, $\Delta = 87°$). On samples with $\Delta$ below $10°$, the $\Delta$ error becomes as large as $\Delta$ itself, or larger.

That the vulnerable region lies at the extremes of $\Delta$ follows from the structure of RAE. The instrument obtains only $\cos\Delta$ and inverts it. Since $\cos\Delta$ is not linear in $\Delta$, its slope vanishes where $\cos\Delta = \pm1$. Change the configuration and the region moves. PME with the modulator at $0°$ measures $S_2$ and $S_3$, so on the $\Psi$ side it obtains only $\sin 2\Psi$, and the error grows at $\Psi = 45°$. Measuring twice in two configurations resolves it, but that is a poor remedy for real-time work.

The three errors so far, together with the polarization-dependent sensitivity of Section 7, collect into one table. The magnitudes are the values injected in this post's calculations.

| Error | Injected | Where it goes | How to diagnose it | Calibration | Figure showing the method |
| --- | --- | --- | --- | --- | --- |
| Polarizer offset $P_s$ | $0.1°$ | $0.055°$ in $\Psi$, none in $\Delta$ | position of the minimum of $R(P)$ | residual, regression | 3, 4 |
| Analyzer offset $A_s$ | $0.1°$ | $0.30°$ in $\Delta$ near $\Delta \approx 90°$ | direction of $$(\alpha', \beta')$$ at the minimum | residual, regression | 3 |
| Detector attenuation $\eta$ | $0.999$ | $2.5°$ in $\Delta$ at $\Delta = 3°$ | depth of the minimum of $R(P)$ | residual, regression | 3 |
| Polarization-dependent sensitivity $\epsilon$ | $0.25$ | if unmodelled, $\eta$ becomes $1.0865$ and $\Psi$ shifts by $-0.65°$ | MSE of the regression fit | include it in the model | 5 |

The conclusion is plain. An analyzer dial off by only $0.1°$ costs $0.3°$ on a sample with $\Delta \approx 90°$, and a one-part-in-a-thousand detector attenuation grows to several degrees in the unfavourable regions. Calibration is not optional.

## 4. How does residual calibration read three numbers off one curve?

The method that came out of Aspnes's 1974 paper and the 1975 paper by Aspnes and Studna rests on a compact idea. With the polarizer exactly at $0°$ there is no s component, so only p-polarized light reaches the sample, and for an isotropic sample the reflected light stays p-polarized. The moment the polarizer leaves $0°$, an s component appears and the reflected light becomes elliptically polarized. **Measure how far it departs and $P_s$ follows.**

The residual function measures how complete the polarization is.

$$R(P) = 1 - (\alpha'^2 + \beta'^2) = 1 - \eta^2(\alpha^2 + \beta^2)$$

For fully linear polarization $\alpha^2 + \beta^2 = 1$, so $R$ vanishes when $\eta = 1$. Three numbers fall out in turn.

$P_s$ comes from measuring $R(P)$ at several polarizer settings, fitting a parabola $R = c_0 + c_1 P + c_2 P^2$, and taking $P_{\min} = -c_1/2c_2$. What matters is that the expression for $R$ depends only on $P - P_s$. $P_s$ can be obtained while $A_s$ is still unknown.

$\eta$ comes from the depth of that minimum. At $P = P_s$ we have $\alpha = 1$ and $\beta = 0$, so $R(P_{\min}) = 1 - \eta^2$. **The position of the minimum gives $P_s$, and the depth of the minimum gives $\eta$.**

$A_s$ is read off at the same setting. At $P - P_s = 0$ the rotation relation of Section 2 reduces to $$(\alpha', \beta') = \eta(\cos 2A_s,\ \sin 2A_s)$$, so $$A_s = \tfrac12 \tan^{-1}(\beta'/\alpha')$$.

Three numbers from one curve. The tidiness carries a condition, however. If the polarizer is a quartz Rochon prism, the optical activity of quartz rotates the plane of polarization along the path. An optical-activity correction is then needed before the procedure holds. The optical activity flagged as a property of Rochon prisms in Post 2 returns here as a bill.

Running it gives the following.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig3-residual-calibration.png" alt="Residual calibration: three numbers from one minimum" width="780">
_Figure 3. Left, $R(P)$ synthesized with $$P_s = 0.26°$$, $$A_s = 64.3°$$, $\eta = 0.9687$, together with its quadratic fit. The position of the minimum returns $$P_s$$ exactly, but the $\eta$ read from its depth is $0.9665$, short of the true value by $2.2 \times 10^{-3}$. Right, the error in recovered $\eta$ as the scan range is narrowed; the error shrinks with the range._

With the scan taken symmetrically about the minimum, $P_s$ and $A_s$ return at machine precision. $\eta$ does not. The injected $0.9687$ comes back as $0.9665$, a systematic shortfall of $2.2 \times 10^{-3}$.

The cause is not noise; no noise was added to the synthesized signal. At first the gap looks like an implementation bug, but it shrinks as the scan range is narrowed. The cause is that **$R(P)$ is not a parabola.** Fitting the data in the left panel with a quadratic leaves a maximum absolute residual of $8.2 \times 10^{-3}$; fitting the same data with a quartic brings it to $2.3 \times 10^{-4}$. Terms of fourth order and above are genuinely present.

$P_s$ and $A_s$ come from the **position** of the minimum and are relatively insensitive to that distortion. On an interval symmetric about the minimum, even-order terms such as the quartic cannot displace the vertex; they only lift its depth. $\eta$ comes from the **value** at the minimum, so it depends directly on how accurately the fitted parabola reaches the floor.

The right panel confirms this. Narrowing the scan from $\pm5°$ to $\pm0.5°$ reduces the $\eta$ error from $2.2 \times 10^{-3}$ to $2.5 \times 10^{-7}$, a factor of roughly $8{,}800$. The narrower the interval, the better the parabolic approximation holds.

The symmetry matters too. Scanning $\pm5°$ centred on the dial zero, without knowing $P_s$, puts the minimum off to one side of the interval. The even-order terms then no longer cancel, and $P_s$ retains an error of $0.014°$. Given that a $P_s$ error of $0.1°$ shifts $\Psi$ by $0.055°$, this is not negligible.

So why not keep narrowing the interval? Because the measured $R$ values then differ by less and less and sink into the noise. **The systematic error of the parabolic approximation and the noise stand in a trade-off, and the interval must be chosen somewhere between them.** Johs notes that practice generally keeps the range under $10°$; that is the result of this trade-off.

Put together as a procedure:

1. If the polarizer is a quartz Rochon prism, correct for optical activity first. As long as the plane of polarization turns inside the prism, the direction the dial reports is wrong from the outset.
2. Scan widely around the dial zero once to find the approximate minimum. A narrow first scan is unlikely to be symmetric about the minimum, and $P_s$ is pulled.
3. Rescan symmetrically over less than $10°$, centred on that point. Symmetry keeps the even-order terms from moving $P_s$; narrowness keeps the systematic error in $\eta$ small. Too narrow, and the curve sinks into the noise.
4. Fit the parabola: read $P_s$ from the position of the minimum and $\eta$ from its depth.
5. Set the polarizer to $P_s$ and read $A_s$ from the direction of $$(\alpha', \beta')$$.
6. Check that the sample's $\Delta$ lies within $30°$–$150°$. If not, do not use this procedure. The next section explains why.

## 5. Why does residual calibration fail when the phase difference nears 0° or 180°?

Residual calibration carries a more fundamental restriction, because the $P$ dependence of $R(P)$ scales with $\sin^2\Delta$ (in Fujiwara's Eq. 4.61c, $\sin\Delta$ enters squared). Expanding that expression about $P = P_s$ shows that the curvature of the minimum scales as $(2\sin\Delta\cot\Psi)^2$. A large $\Psi$ also flattens it. Taking the curvature at $\Psi = 27°$ as 1, it is $0.26$ at $\Psi = 45°$ and $0.019$ at $\Psi = 75°$.

A sample with little absorption has $\Delta$ near $0°$ or $\pm180°$. Then $\sin\Delta \to 0$ and $R(P)$ becomes nearly independent of $P$. The minimum flattens, and the position of a flat floor cannot be read.

**Post 1 identified $\Delta \cong 0°, 180°$ as the weak region of RAE. The calibration method inherits that weakness exactly.** Where the signal being measured has vanished, that same signal cannot calibrate the instrument. This is also why Ilsin An recommends residual calibration only for samples with $\Delta$ between $30°$ and $150°$.

In a noiseless calculation the blind spot does not appear. However small the curvature, perfect data still locate the minimum exactly. The blind spot emerges when the curvature sinks below the noise, so noise must be added to the detected intensity and the trial repeated.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig4-blind-spot.png" alt="Near Delta = 0 or 180 degrees, residual calibration loses its minimum" width="780">
_Figure 4. For a sample with $\Psi = 27°$, the depth of the minimum as $\Delta$ is lowered (left) and the error in recovered $$P_s$$ (right). For the right panel, Gaussian noise of standard deviation $2\times10^{-4}$ was added to the detected intensity (dc component $0.25$) over 40 trials, and the median is taken over the successful fits only. The dashed line marks $0.01°$. The error stays below it down to $\Delta = 30°$ and crosses it at $\Delta = 10°$, where it reaches $0.047°$._

Going from $\Delta = 90°$ to $\Delta = 1°$, the curvature of the minimum collapses by more than three orders of magnitude. The scatter that noise leaves in $R$ stays near $3.6 \times 10^{-4}$. From $\Delta = 3°$ downward, the depth of the minimum is shallower than that scatter.

The right panel shows the consequence. The median error in recovered $P_s$ (successful fits only) grows from $0.0016°$ at $\Delta = 90°$ to $1.15°$ at $\Delta = 1°$, a factor of about 700. At $\Delta = 1°$, 19 of 40 fits fail outright: the fitted parabola turns concave down and has no minimum at all.

Two alternatives exist. One is zone-difference calibration, which uses a different calibration function; the procedure is similar but performs better than residual calibration for $\Delta < 30°$ and $\Delta > 150°$. The other is the regression calibration of the next section.

## 6. How does regression calibration get past those limits?

The paper Johs published in 1993 states the limits of residual calibration at the outset. One is that the usable samples and the accessible wavelength and angle-of-incidence ranges are restricted by $\Delta$. The other is that data are taken only over the narrow interval where the parabolic approximation holds, so **defects arising from misalignment or non-ideal components can pass unnoticed**.

The shift in approach can be put in one sentence. Instead of hunting for special conditions under which one parameter separates, **model the response of the optical system with exact expressions and fit every parameter at once.**

Scan the polarizer broadly from $5°$ to $85°$, measure $$(\alpha', \beta')$$ at each setting, and fit the model below with the [Levenberg-Marquardt (LM) method](/en/posts/optimization-levenberg-marquardt/).

$$\alpha_0 = \frac{\tan^2\Psi - \tan^2(P - P_s)}{\tan^2\Psi + \tan^2(P - P_s)}, \qquad \beta_0 = \frac{2\tan\Psi\cos\Delta\,\tan(P - P_s)}{\tan^2\Psi + \tan^2(P - P_s)}$$

Here $\alpha_0$ and $\beta_0$ are the true coefficients $\alpha$ and $\beta$ of Section 2, written in terms of the dial reading $P$ and the offset $P_s$. The model is this with the rotation and scale factor of Section 2 applied, and the free parameters are the five $\Psi$, $\Delta$, $P_s$, $A_s$, $\eta$. That the optical constants of the sample need not be known in advance is essential. The sample is fitted along with the instrument.

With the parabolic approximation gone there is no reason to stay in a narrow interval, and where $\sin\Delta$ goes to zero the $P$ dependence of $\alpha_0$ survives. Take the $\Delta = 1°$ sample that failed 19 times out of 40 in the previous section, add the same noise ($\sigma = 2\times10^{-4}$), and run regression calibration over this wide range. All 40 fits converge, with a median $P_s$ error of $0.0009°$, more than a thousand times smaller than the $1.15°$ of residual calibration. Regression calibration does not lose the minimum at what residual calibration calls a blind spot. How the precision fares is seen in the measurement below.

Measurement says the same. Table 1 of Johs's paper reports calibrations of an SiO$_2$/Si sample and a glass (Corning 7059) sample at $75°$ incidence and 500 nm, with the polarizer scanned over $90°$ in 90 points. The glass sample has $\Delta = 2.01°$, where the residual function becomes a flat line and residual calibration is impossible.

Regression calibration still returns $P_s$ to $\pm0.014°$ and $\eta = 0.9692 \pm 0.0004$, the same as the $0.9687$ from SiO$_2$/Si. Since $\eta$ belongs to the instrument rather than the sample, it should agree.

The 90% confidence interval on $\Delta$, however, is $\pm0.984°$, 45 times wider than the $\pm0.022°$ for SiO$_2$/Si. The blind spot does not vanish; it moves into the confidence interval. A regression result must therefore be read together with its confidence intervals. Johs regards intervals below $0.01°$ on $\Psi$, $\Delta$, $P_s$, and $A_s$ as ideal.

There is a price. Initial values are required, and when parameters are strongly correlated, their individual values cannot be separated. Above all, **a model cannot account for what is not in it.**

## 7. Where does a defect missing from the model hide?

This is where regression calibration becomes dangerous. If an effect is present in the instrument but absent from the model, the fit **absorbs it into another parameter of the model.** The fit still converges, and the result still looks plausible.

Take the polarization-dependent sensitivity of the detector from the seven branches in Section 1. Let the detector sensitivity vary with the direction of the incoming polarization as $1 + \epsilon\cos 2A$, synthesize the signal, and fit it with a model that omits that term.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig5-model-mismatch.png" alt="An unmodelled effect hides in eta and the MSE reports it" width="780">
_Figure 5. Fitted $\eta$ and the error in $\Psi$ from a model without the term, as the polarization-dependent sensitivity $\epsilon$ is increased (left), and the MSE of the same fits (right). At $\epsilon = 0.25$, $\eta$ comes out as $1.0865$ and $\Psi$ shifts by $0.65°$. The fit converges, but the MSE rises almost in proportion to $\epsilon$ and exposes the hidden effect._

At $\epsilon = 0.25$, $\eta$ comes out as $1.0865$ instead of its true value $1.0000$. $\Psi$ also shifts from $16.800°$ to $16.147°$, by $0.65°$. **The omitted effect hides inside $\eta$, and the sample parameters share the cost.** The remark in Section 2 that $\eta$ is not the coefficient of a single cause returns in this form.

How is a hidden effect detected? The mean squared error (MSE) of the fit reports it. MSE here is the residual sum of squares divided by the degrees of freedom and then square-rooted. Johs's Eqs. 6–7 weight each residual by the measurement standard deviation; the noiseless data here are left unweighted. With a model that matches exactly, the MSE is $10^{-16}$, numerically zero. With the polarization-dependent sensitivity present it rises almost linearly in $\epsilon$ for small $\epsilon$ (as $\epsilon^{0.997}$), reaching $3.0\times10^{-2}$ at $\epsilon = 0.25$.

MSE is therefore not a measure of fit quality but a **diagnostic for model mismatch**. In a real measurement weighted by the standard deviation, a correct model gives an MSE near 1; Johs's Table 1 shows 1.14 and 1.06. An MSE well above that means the data contain something the model does not. Table 2 of the same paper shows this structure.

(The polarization-dependent sensitivity used here is the simplest $\cos 2A$ form and does not reproduce the paper's numbers. What is reproduced is not a figure but the **structure**: an unmodelled effect contaminates $\eta$, and the MSE reveals it.)

## 8. How does zone averaging remove defects without a model?

Long before regression calibration, the null-ellipsometry era used a different strategy: not modelling the defects, but **cancelling them through the symmetry of the measurement arrangement.**

Place the compensator at $C = +\pi/4$ and adjust polarizer and analyzer to obtain two nulls; turn the compensator by $90°$ to $C = -\pi/4$ and obtain two more. These four arrangements are called zones, and averaging the four pairs $(P_i, A_i)$ is four-zone averaging.

The range of cancellation is unexpectedly wide. In treating two-zone averaging, Azzam and Bashara list imperfections such as azimuth errors, window birefringence and dichroism, imperfections of polarizer and analyzer, and non-ideal compensator behaviour. Four-zone averaging, they write, keeps **all** of these from reaching $\Psi$, with a single exception, the depolarization of the polarizer output. None of them reaches $\Delta$ either, with one exception, the birefringence of the cell windows. That list is not the full set of seven branches in Section 1, though. Beam deviation, stray beams, and the polarization-dependent sensitivity of the detector are not on it.

"Not reaching" here does not mean the error shrinks. To first order in small imperfections, the coupling coefficients **cancel identically** as the zones are added and subtracted. And which two survive is fixed as well. This is the content behind the single line in Section 3 of [Post 1](/en/posts/ellipsometer-oblique-rotating-element/) that four-zone averaging cancels the alignment errors of the components.

The procedure for actually locating a null rests on the same structure. Dither the polarizer slightly, read the two angles $P^+$ and $P^-$ at which the signal is equal on either side of the minimum, and take $P_n = \tfrac12(P^+ + P^-)$. This is the method of swings. It works because the detected signal is a **symmetric, parabolic function** of the departure from the null.

That is the same idea as fitting $R(P)$ with a parabola in residual calibration. It is also exactly what regression calibration discards. The three methods thread onto one line.

## 9. What do the three methods share, and where do they part?

The methods differ; the principle does not.

| Era | Method | Signal that reveals the non-ideality |
| --- | --- | --- |
| Null | Zone averaging | How far the four null angles depart from the ideal zone relations |
| Parabola | Residual calibration | How far the minimum of $R(P)$ departs from the ideal value (0) |
| Regression | MSE | How far the MSE exceeds what the noise would give |

All three set up a relation that must hold for an ideal instrument and gauge the defect by how far it breaks. Post 3 noted that DRC-MME solves fifteen normalized Mueller elements from twenty-four surviving Fourier coefficients and has nine left over. This is what those nine are for. The dual-rotating method of Post 3, which is outside the table, also gauges defects by how far five Mueller-matrix identities are violated, so it fits the same frame.

The philosophies diverge. One is **cancelling by symmetry**. Zone averaging does not model the defects. It arranges the measurement so that particular classes of defect cancel, and knows in advance which ones vanish and which remain. The other is **measuring by model**. Regression calibration sets the defects up as parameters and fits them alongside everything else. It requires a model, but it can handle arbitrary defects and **returns the defect values themselves as a bonus.**

Post 2 described the cost of a pseudo zero-order compensator and called using cheap components and making up for them in calibration the design philosophy of commercial rotating-compensator instruments, deferring the argument to this post. This is that fork. If measuring a component defect and subtracting it is cheaper than eliminating it, the ability to measure defects becomes a specification of the instrument.

## 10. When should the instrument be recalibrated?

How long do calibration values stay valid once obtained? The answer comes from separating what belongs to the instrument from what changes with each measurement. Section 1's split of azimuth errors into residual divided-circle error and sample remounting error is the starting point.

Table 1 of Johs's paper puts numbers on this. It reports two calibrations on the same instrument, one on SiO$_2$/Si and one on glass. $A_s$ comes out as $64.322°$ and $64.318°$, $0.004°$ apart, and $\eta$ as $0.9687$ and $0.9692$, equal within the confidence interval ($\pm0.0004$). $P_s$, however, comes out as $42.337°$ and $42.392°$, $0.055°$ apart, which is larger than either confidence interval ($\pm0.005°$, $\pm0.014°$). Since the plane of incidence is set by the tilt of the sample surface, this difference appears to be remounting error, which changes each time a sample is put in place.

It is not small. Reusing the previous sample's $P_s$ would shift $\Psi$ by $0.030°$ on SiO$_2$/Si and $0.045°$ on glass, three to four times the $0.01°$ Johs treats as ideal. It does not reach $\Delta$ (Section 3). $\eta$ needs attention too. An $\eta$ off by only $5\times10^{-4}$ costs just $0.005°$ in $\Delta$ on a sample with $\Delta = 87°$, but $0.74°$ on one with $\Delta = 3°$.

The decision therefore splits. If you remeasure the same sample without touching the mount, the previous calibration values can be reused. After exchanging the sample, it is safer to find at least $P_s$ again. For samples with $\Delta$ near $0°$ or $180°$, refit $\eta$ on the spot as well. Regression calibration fits all five parameters for every sample, so each measurement is itself a recalibration. And if the MSE climbs above its usual level on the same instrument and the same kind of sample, suspect the instrument before the sample (Section 7). Ilsin An reaches the same judgement: if today's sample has the same thickness as yesterday's and the mount has not been touched, yesterday's calibration can be used; if high precision is needed or there is drift, calibrate again.

Calibration, in other words, is not something done once when the instrument is installed.

## Summary

The angle a dial reports is not the true angle, and the difference can be amplified several times in the result. What has to be undone, fortunately, reduces to one rotation and one scale factor.

Residual calibration reads all three numbers off the minimum of a single curve. Because it leans on a parabolic approximation, though, it leaves a systematic error in $\eta$ and loses the minimum altogether on samples with $\Delta$ near $0°$ or $180°$. Regression calibration drops that approximation and escapes both restrictions, at the risk of absorbing unmodelled effects into other parameters. The MSE is what watches for that risk. Zone averaging removes defects by symmetry, without a model. Of the calibration values, $A_s$ and $\eta$ stay with the instrument, while $P_s$ can change when the sample is remounted, so it is safer to find it again.

If Post 1 asked what cannot be measured, Post 2 where the components stop working, and Post 3 why one element is not enough, this post asks **how to work with limitations once they are known.**

## References

- D. E. Aspnes, "Effects of component optical activity in data reduction and calibration of rotating-analyzer ellipsometers," *J. Opt. Soc. Am.* **64**, 812–819 (1974) — the origin of residual calibration.
- D. E. Aspnes and A. A. Studna, "High precision scanning ellipsometer," *Appl. Opt.* **14**, 220 (1975) — cited alongside the 1974 paper as an original source by Johs and by Fujiwara.
- B. Johs, "Regression calibration method for rotating element ellipsometers," *Thin Solid Films* **234**, 395 (1993) — the turning point at which the parabolic approximation is abandoned. Table 1 (calibration values for two samples) and Table 2 (absorption of polarization-dependent sensitivity).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, §3.8.2 and §5.4 — classification of error sources and zone averaging.
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, §4.3.3–4.3.4 and §4.4.1 — calibration procedures and the error-prone regions of each configuration.
- Ilsin An, *Ellipsometry* (in Korean, 《엘립소메트리: 타원편광분석기술의 기초와 이해》), 2nd ed., Hanyang University Press, 2017, Ch. 5 — the $\Delta$ range for residual calibration and when to reuse calibration values.
- Youngjoon Kim, "Development of snapshot angle-resolved ellipsometry using a line-scan spectrograph and back focal plane spectral interference," Ph.D. dissertation, Seoul National University, 2025.

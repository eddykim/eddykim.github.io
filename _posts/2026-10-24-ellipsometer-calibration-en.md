---
date: 2026-10-24 20:00:00 +0900
layout: post
title: "Ellipsometer Instrumentation 4 — Calibration: Between the Dial and the True Angle"
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

All three rested on an assumption. When the polarizer dial reads $30°$, the transmission axis is taken to sit $30°$ from the plane of incidence. That assumption is false. No matter how precisely a divided circle is cut, the optical transmission axis and the mechanical zero do not coincide, and a detector does not return a signal strictly proportional to the light that enters it.

The question here follows from that. **How do you measure the difference between what the dial reads and the true angle?** The instrument built to measure a sample must now measure itself, and the measuring tool is already misaligned. That is the whole difficulty.

Three eras answered it differently, and this post follows them in order. Yet all three use one idea. Set up a relation that must hold if the instrument were ideal, then **gauge the defect by how far that relation is broken.** It is the general form of what Post 3 called the surplus observables watching over the instrument.

## 1. Where the errors enter

Before looking at methods, it is worth counting what needs correcting. Azzam and Bashara divide the paths by which $\rho = \tan\Psi\,e^{i\Delta}$ picks up systematic error into seven.

Azimuth errors; optical imperfections in components and cell windows; beam deviation; stray beams; polarization and collimation errors in the source beam; polarization-dependent sensitivity of the detector; and residual mechanical errors.

Azimuth errors split further by origin. Residual calibration error of the divided circle stays constant once measured, while sample remounting error changes every time the sample is exchanged or replaced. The first is a property of the instrument; the second is generated afresh with each measurement.

Component imperfections are worse. Polarizing prisms have finite extinction ratios, cell windows carry stress birefringence, compensators scatter from defects. Azzam and Bashara write that these are **unavoidable**. The component specifications surveyed in Post 2 come back here as an error list.

The detector side does not end in one line either. There is dark signal, removed by inserting and withdrawing a mechanical shutter from the beam path. There is nonlinear response, in which output is not proportional to input intensity. There is polarization-dependent sensitivity, produced when residual stress in a photomultiplier (PMT) window acts photoelastically. The photoelastic effect that served as a modulation tool in Post 2 appears here as a defect. Array detectors add image persistence, in which signal survives readout.

The list also shifts with the configuration. A rotating-compensator ellipsometer (RCE) must first fix the compensator retardance at each wavelength, which is done in a straight-through arrangement with the sample removed and only the compensator between polarizer and analyzer. Because of the $\delta \propto 1/\lambda$ relation from Post 2, the retardance of an MgF$_2$ compensator rises nearly linearly with photon energy, and parameterizing that dispersion lets it be recovered in the general arrangement as well. Phase-modulation ellipsometry (PME) needs retardance control, wavelength-by-wavelength control, and temperature control, plus a higher-harmonic correction. The higher-order Bessel terms from Post 1 and the temperature sensitivity of the photoelastic modulator from Post 2 meet here.

A list this long cannot be worked through one item at a time. Practical calibration takes a different route.

## 2. Calibration undoes one rotation and one scale factor

Written out with the errors included, the detected intensity of a rotating-analyzer ellipsometer (RAE) reads

$$I_t(t) = I_0\big[1 + \eta\,\alpha\cos 2(\omega t - A_s) + \eta\,\beta\sin 2(\omega t - A_s)\big]$$

$P_s$ and $A_s$ are the offsets of the polarizer and analyzer zero positions from the plane of incidence. When the dial reads $P$, the true angle is $P - P_s$, and likewise for the analyzer. The coefficient $\eta$ carries the detector non-ideality: it states how far the ac component is attenuated relative to the dc component. The nonlinear response and the polarization-dependent sensitivity of the previous section both enter here at once. That $\eta$ is not the coefficient of a single cause becomes a problem later.

What matters is the relation between the measured $(\alpha', \beta')$ and the true $(\alpha, \beta)$.

$$\begin{bmatrix}\alpha'\\ \beta'\end{bmatrix} = \eta\begin{bmatrix}\cos 2A_s & -\sin 2A_s\\ \sin 2A_s & \cos 2A_s\end{bmatrix}\begin{bmatrix}\alpha\\ \beta\end{bmatrix}$$

This is a coordinate rotation by $A_s$ multiplied by a scale factor $\eta$. Knowing $A_s$ and $\eta$ is therefore enough: invert the matrix and the true values come straight back. The error list ran to seven branches, but what must be undone reduces to **one rotation and one scale factor**. That reduction is what makes calibration tractable.

Undoing $A_s$ mathematically matters most in multichannel instruments. A photodiode array is read out pixel by pixel, and the analyzer keeps turning during readout, so **$A_s$ differs from pixel to pixel.** No mechanical adjustment can fix that, so it has to be handled by calibration from the start. Post 1 noted that RAE is achromatic and therefore measures as many wavelengths simultaneously as the array has pixels; this is the price attached to that strength.

The three errors leave different marks on the signal.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig1-error-anatomy.png" alt="What the three errors do to the measured signal" width="780">
_Figure 1. Left, the detected intensity over one analyzer revolution. Right, the locus traced by $(\alpha', \beta')$ as the polarizer is scanned. The right panel uses a sample with large $\lvert\cos\Delta\rvert$ ($\Psi = 35°$, $\Delta = 50°$) so the ellipse is not flattened._

In the waveform the three look like similar distortions, but in the $(\alpha', \beta')$ plane they separate cleanly. Scanning the polarizer from $-90°$ to $+90°$ traces a closed ellipse with semi-axes $1$ and $\lvert\cos\Delta\rvert$; $\Psi$ does not enter. Substituting $u = \tan(P - P_s) = \tan\Psi\,\tan\theta$ gives $\alpha = \cos 2\theta$ and $\beta = \cos\Delta\,\sin 2\theta$, and $\Psi$ cancels. $A_s$ rotates this ellipse about the origin, $\eta$ shrinks it whole, and $P_s$ **leaves the ellipse untouched and only slides the operating point along it.**

Because the three act differently, they can be measured separately. Residual calibration uses exactly that separation.

## 3. How wrong the answer is without calibration

Before the methods, it helps to fix the size of the problem. Synthesize a signal with one error at a time, then feed the dial readings straight into the ideal formulas and read off $\Psi$ and $\Delta$.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig2-uncalibrated-error.png" alt="Uncalibrated readings fed straight into the ideal formulas" width="780">
_Figure 2. Errors propagated into $\Psi$ and $\Delta$ for a sample with $\Psi = 16.8°$, sweeping $\Delta$, with one error injected at a time. The dotted line marks the $0.01°$ level._

The three errors go to different places. A polarizer offset $P_s = 0.1°$ shifts $\Psi$ by $0.055°$ but does **not** propagate into $\Delta$ at all. In the ideal inversion $\cos\Delta = \beta' / \sqrt{1 - \alpha'^2}$, and the polarizer angle cancels between numerator and denominator. The $10^{-13}$ degrees left in the computation are floating-point noise.

An analyzer offset $A_s = 0.1°$ does the opposite. It barely touches $\Psi$ and shifts $\Delta$ by $0.30°$, amplified threefold over the injected error.

The detector coefficient $\eta$ leaves $\Delta$ untouched at $\Delta = 90°$ and degrades sharply as $\Delta$ approaches $0°$ or $180°$. An attenuation of one part in a thousand, $\eta = 0.999$, produces a $2.5°$ error on a sample with $\Delta = 3°$. This is the region where Post 1 put the RAE amplification factor at $1/\lvert\sin\Delta\rvert$.

That the vulnerable region lies at the extremes of $\Delta$ follows from the structure of RAE. The instrument obtains only $\cos\Delta$ and inverts it, and since $\cos\Delta$ is not linear in $\Delta$, its slope vanishes where $\cos\Delta = \pm1$. Change the configuration and the region moves. PME measures $S_2$ and $S_3$, so on the $\Psi$ side it obtains only $\sin 2\Psi$, and the error grows at $\Psi = 45°$. Measuring twice in two configurations resolves it, but that is a poor remedy for real-time work.

The conclusion is plain. A $0.1°$ alignment error appears as $0.3°$ in the result, and several degrees in the unfavourable regions. Calibration is not optional.

## 4. Residual calibration reads three numbers off one curve

The method Aspnes proposed in 1974 rests on a compact idea. With the polarizer exactly at $0°$ there is no s component, so only p-polarized light reaches the sample, and for an isotropic sample the reflected light stays p-polarized. The moment the polarizer leaves $0°$, an s component appears and the reflected light becomes elliptically polarized. **Measure how far it departs and $P_s$ follows.**

The residual function measures how complete the polarization is.

$$R(P) = 1 - (\alpha'^2 + \beta'^2) = 1 - \eta^2(\alpha^2 + \beta^2)$$

For fully linear polarization $\alpha^2 + \beta^2 = 1$, so $R$ vanishes when $\eta = 1$. Three numbers fall out in turn.

$P_s$ comes from measuring $R(P)$ at several polarizer settings, fitting a parabola $R = c_0 + c_1 P + c_2 P^2$, and taking $P_{\min} = -c_1/2c_2$. What matters is that the expression for $R$ depends only on $P - P_s$. $P_s$ can be obtained while $A_s$ is still unknown.

$\eta$ comes from the depth of that minimum. At $P = P_s$ we have $\alpha = 1$ and $\beta = 0$, so $R(P_{\min}) = 1 - \eta^2$. **The position of the minimum gives $P_s$, and the depth of the minimum gives $\eta$.**

$A_s$ is read off at the same setting. At $P - P_s = 0$ the rotation relation of the previous section reduces to $(\alpha', \beta') = \eta(\cos 2A_s,\ \sin 2A_s)$, so $A_s = \tfrac12 \tan^{-1}(\beta'/\alpha')$.

Three numbers from one curve. The tidiness carries a condition, however. If the polarizer is a quartz Rochon prism, the optical activity of quartz rotates the plane of polarization along the path, and an optical-activity correction is needed before the procedure holds. The optical activity flagged as a property of Rochon prisms in Post 2 returns here as a bill.

Running it gives the following.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig3-residual-calibration.png" alt="Residual calibration: three numbers from one minimum" width="780">
_Figure 3. Left, $R(P)$ synthesized with $P_s = 0.26°$, $A_s = 64.3°$, $\eta = 0.9687$, together with its quadratic fit. Right, the error in recovered $\eta$ as the scan range is narrowed._

With the scan taken symmetrically about the minimum, $P_s$ and $A_s$ return at machine precision. $\eta$ does not. The injected $0.9687$ comes back as $0.9665$, a systematic shortfall of $2.2 \times 10^{-3}$.

The cause is not noise; no noise was added to the synthesized signal. The cause is that **$R(P)$ is not a parabola.** Fitting the data in the left panel with a quadratic leaves a residual of $8.2 \times 10^{-3}$; fitting the same data with a quartic brings it to $2.3 \times 10^{-4}$. Terms of fourth order and above are genuinely present.

$P_s$ and $A_s$ come from the **position** of the minimum and are relatively insensitive to that distortion. On an interval symmetric about the minimum, even-order terms such as the quartic cannot displace the vertex; they only lift its depth. $\eta$ comes from the **value** at the minimum, so it depends entirely on how accurately the fitted parabola reaches the floor.

The right panel confirms this. Narrowing the scan from $\pm5°$ to $\pm0.5°$ reduces the $\eta$ error from $2.2 \times 10^{-3}$ to $2.5 \times 10^{-7}$, a factor of roughly $8{,}800$. The narrower the interval, the better the parabolic approximation holds.

The symmetry matters too. Scanning $\pm5°$ centred on the dial zero, without knowing $P_s$, puts the minimum off to one side of the interval; the even-order terms no longer cancel and $P_s$ retains an error of $0.014°$. Given that a $P_s$ error of $0.1°$ shifts $\Psi$ by $0.055°$, this is not negligible. Hence the practical procedure of scanning once to locate the approximate minimum, then rescanning centred on it.

So why not keep narrowing the interval? Because the measured $R$ values then differ by less and less and sink into the noise. **The systematic error of the parabolic approximation and the noise stand in a trade-off, and the interval must be chosen somewhere between them.** Johs notes that practice generally keeps the range under $10°$; that is the result of this trade-off.

## 5. Where the measurement is dark, the calibration is dark too

Residual calibration carries a more fundamental restriction, because the $P$ dependence of $R(P)$ scales with $\sin\Delta$.

A sample with little absorption has $\Delta$ near $0°$ or $\pm180°$. Then $\sin\Delta \to 0$ and $R(P)$ becomes nearly independent of $P$. The minimum flattens, and the position of a flat floor cannot be read.

**Post 1 identified $\Delta \cong 0°, 180°$ as the weak region of RAE, and the calibration method inherits that weakness exactly.** Where the signal being measured has vanished, that same signal cannot calibrate the instrument.

In a noiseless calculation the blind spot does not appear. However small the curvature, perfect data still locate the minimum exactly. The blind spot emerges when the curvature sinks below the noise, so noise must be added to the detected intensity and the trial repeated.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig4-blind-spot.png" alt="Near Delta = 0 or 180 degrees, residual calibration loses its minimum" width="780">
_Figure 4. For a sample with $\Psi = 27°$, the depth of the minimum as $\Delta$ is lowered (left), and the error in recovered $P_s$ over 40 trials with Gaussian noise of standard deviation $2\times10^{-4}$ on the detected intensity (right)._

Going from $\Delta = 90°$ to $\Delta = 1°$, the curvature of the minimum collapses by more than three orders of magnitude. The scatter that noise leaves in $R$ stays near $3.6 \times 10^{-4}$, and from $\Delta = 3°$ downward the depth of the minimum is shallower than that scatter.

The right panel shows the consequence. The median error in recovered $P_s$ grows from $0.0016°$ at $\Delta = 90°$ to $1.15°$ at $\Delta = 1°$, a factor of about 700. At $\Delta = 1°$, 19 of 40 fits fail outright: the fitted parabola turns concave down and has no minimum at all.

Two alternatives exist. One is zone-difference calibration, which uses a different calibration function; the procedure is similar but performs better than residual calibration for $\Delta < 30°$ and $\Delta > 150°$. The other is the regression calibration of the next section.

## 6. Regression calibration drops the approximation and fits the whole model

The paper Johs published in 1993 states the limits of residual calibration at the outset. One is that the usable samples and the accessible wavelength and angle-of-incidence ranges are restricted by $\Delta$. The other is that data are taken only over the narrow interval where the parabolic approximation holds, so **defects arising from misalignment or non-ideal components can pass unnoticed**.

The shift in approach summarizes as follows. Instead of hunting for special conditions under which one parameter separates, **model the response of the optical system with exact expressions and fit every parameter at once.**

Scan the polarizer broadly from $5°$ to $85°$, measure $(\alpha', \beta')$ at each setting, and fit the model below with the [Levenberg-Marquardt (LM) method](/en/posts/optimization-levenberg-marquardt/).

$$\alpha_0 = \frac{\tan^2\Psi - \tan^2(P - P_s)}{\tan^2\Psi + \tan^2(P - P_s)}, \qquad \beta_0 = \frac{2\tan\Psi\cos\Delta\,\tan(P - P_s)}{\tan^2\Psi + \tan^2(P - P_s)}$$

The model is this with the rotation and scale factor of Section 2 applied, and the free parameters are the five $\Psi$, $\Delta$, $P_s$, $A_s$, $\eta$. That the optical constants of the sample need not be known in advance is essential. The sample is fitted along with the instrument.

With the parabolic approximation gone there is no reason to stay in a narrow interval, and where $\sin\Delta$ goes to zero the $P$ dependence of $\alpha_0$ survives. Running regression calibration on the same data for the $\Delta = 1°$ sample that failed 19 times out of 40 in the previous section returns $P_s$ to a precision of $10^{-17}$. What residual calibration calls a blind spot is not one here.

There is a price. Initial values are required, and fitting five parameters simultaneously leaves correlated combinations poorly resolved. Above all, **a model cannot account for what is not in it.**

## 7. Converging is not the same as being right

This is where regression calibration becomes dangerous. If an effect is present in the instrument but absent from the model, the fit **absorbs it into another parameter of the model.** Convergence still occurs, and the result still looks plausible.

Take the polarization-dependent sensitivity of the detector from the seven branches in Section 1. Let the detector sensitivity vary with the direction of the incoming polarization as $1 + \epsilon\cos 2A$, synthesize the signal, and fit it with a model that omits that term.

<img src="/assets/img/posts/ellipsometer-calibration/en/fig5-model-mismatch.png" alt="An unmodeled effect hides in eta and the MSE reports it" width="780">
_Figure 5. Fitted $\eta$ and $\Psi$ from a model without the term, as the polarization-dependent sensitivity $\epsilon$ is increased (left), and the MSE of the same fits (right)._

At $\epsilon = 0.25$, $\eta$ comes out as $1.0865$ instead of its true value $1.0000$. $\Psi$ also shifts from $16.800°$ to $16.147°$, by $0.65°$. **The omitted effect hides inside $\eta$, and the sample parameters share the cost.** The remark in Section 2 that $\eta$ is not the coefficient of a single cause returns in this form.

How is a hidden effect detected? The MSE of the fit reports it. MSE here is the residual sum of squares divided by the degrees of freedom and then square-rooted, the definition used by spectroscopic ellipsometry software. With a model that matches exactly, noiseless data give $10^{-16}$, numerically zero. With the polarization-dependent sensitivity present it rises almost linearly in $\epsilon$ for small $\epsilon$ (as $\epsilon^{0.997}$), reaching $3.0\times10^{-2}$ at $\epsilon = 0.25$.

MSE is therefore not a measure of fit quality but a **diagnostic for model mismatch**. In a real measurement there is an expected value set by the noise level, and an MSE above it means the data contain something the model does not. This is precisely the structure Johs shows in Table 2 of that paper.

(The polarization-dependent sensitivity used here is the simplest $\cos 2A$ form and does not reproduce the paper's numbers. What is reproduced is not a figure but the **structure**: an unmodeled effect contaminates $\eta$, and the MSE reveals it.)

## 8. Zone averaging cancels by symmetry instead of modelling

Long before regression calibration, the null-ellipsometry era used a different strategy: not modelling the defects, but **cancelling them through the symmetry of the measurement arrangement.**

Place the compensator at $C = +\pi/4$ and adjust polarizer and analyzer to obtain two nulls; turn the compensator by $90°$ to $C = -\pi/4$ and obtain two more. Averaging the four pairs $(P_i, A_i)$ is four-zone averaging.

The range of cancellation is unexpectedly wide. Azzam and Bashara write that four-zone averaging frees $\Psi$ from **all** of the imperfections listed earlier, with a single exception, the depolarization of the polarizer output. $\Delta$ is likewise free of all of them, with one exception, the birefringence of the cell windows.

"Free" here does not mean the error shrinks. In adding and subtracting the zones, the coupling coefficients **cancel identically** and disappear. And the two survivors are named explicitly. This is the content behind the single line in Post 1 that four-zone averaging cancels the alignment errors of the components.

The procedure for actually locating a null rests on the same structure. Dither the polarizer slightly, read the two angles $P^+$ and $P^-$ at which the signal is equal on either side of the minimum, and take $P_n = \tfrac12(P^+ + P^-)$; this is the method of swings. It works because the detected signal is a **symmetric, parabolic function** of the departure from the null.

That is exactly the idea behind fitting $R(P)$ with a parabola in residual calibration, and exactly what regression calibration discards. The three methods thread onto one line.

## 9. One idea across three eras

The methods differ; the principle does not.

| Era | Method | Signal that reveals the non-ideality |
| --- | --- | --- |
| Null | Zone relations | How far the four null angles depart from the ideal zone relations |
| Dual rotating | Identities | How far five Mueller-matrix identities are violated |
| Regression | MSE | How far the MSE exceeds what the noise would give |

All three set up a relation that must hold for an ideal instrument and gauge the defect by how far it breaks. Post 3 noted that DRC-MME solves fifteen normalized Mueller elements from twenty-four surviving Fourier coefficients and has nine left over; this is what those nine are for.

The philosophies diverge. Zone averaging does not model the defects. It arranges the measurement so that particular classes of defect cancel, and knows in advance which ones vanish and which remain. Regression calibration sets the defects up as parameters and fits them alongside everything else. It requires a model, but it can handle arbitrary defects and **returns the defect values themselves as a bonus.**

Post 2 described the cost of a pseudo zero-order compensator and called using cheap components and making up for them in calibration the design philosophy of commercial rotating-compensator instruments, deferring the argument to this post. This is that fork. If measuring a component defect and subtracting it is cheaper than eliminating it, the ability to measure defects becomes a specification of the instrument.

## Summary

The angle a dial reports is not the true angle, and a difference of only $0.1°$ shows up as $0.3°$ in the result. What has to be undone, fortunately, reduces to one rotation and one scale factor.

Residual calibration reads all three numbers off the minimum of a single curve. Because it leans on a parabolic approximation, though, it leaves a systematic error in $\eta$ and loses the minimum altogether on samples with $\Delta$ near $0°$ or $180°$. Regression calibration drops that approximation and escapes both restrictions, at the risk of absorbing unmodeled effects into other parameters. The MSE is what watches for that risk.

If Post 1 asked what cannot be measured, Post 2 where the components stop working, and Post 3 why one element is not enough, this post asks **how to work with limitations once they are known.**

## Reproducing the calculations

Every number and figure in this post is computed from products of component Mueller matrices alone. Closed-form expressions from the literature are kept out of the synthesis and reserved as an independent path for the verification script to check against. This practice caught a transcription error from the references in each of the three preceding posts.

- `calibration.py` — error injection, residual calibration, regression calibration, ideal inversion
- `verify_calibration.py` — checks against Fujiwara's Eqs. 4.59–4.61 and Johs's Eqs. 2–5. It also confirms that the expressions reduce to the ideal RAE of Post 1 when the errors are set to zero, that the inversion returns the true values, that the locus semi-axes are $1$ and $\lvert\cos\Delta\rvert$, that injected values are recovered, and that residual and regression calibration agree
- `generate_figures.py` — Figures 1–5 in both languages

The full code sits in `_code/ellipsometer-calibration/`.

The first thing verification caught was the systematic error in $\eta$ in Section 4. The round trip would not close, which pointed at the implementation, but narrowing the scan range shrank the error with it. It was not a defect in the code; it was the price the method pays for its parabolic approximation. The right panel of Figure 3 is what came out of that check.

## References

- D. E. Aspnes, *J. Opt. Soc. Am.* **64**, 812 (1974) — the origin of residual calibration. The citation is reproduced as given in Fujiwara §4.3.3.
- B. Johs, "Regression calibration method for rotating element ellipsometers," *Thin Solid Films* **234**, 395 (1993) — the turning point at which the parabolic approximation is abandoned.
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, §3.8.2 and §5.4 — classification of error sources and zone averaging.
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, §4.3.3–4.3.4 and §4.4.1 — calibration procedures and the error-prone regions of each configuration.
- Youngjoon Kim, "Development of snapshot angle-resolved ellipsometry using a line-scan spectrograph and back focal plane spectral interference," Ph.D. dissertation, Seoul National University, 2025.

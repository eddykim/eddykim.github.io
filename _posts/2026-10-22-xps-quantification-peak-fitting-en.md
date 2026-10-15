---
title: "XPS Fundamentals 4 — Why Do Different Analysts Get Different Compositions from the Same Spectrum?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/xps-quantification-peak-fitting/
page_id: xps-quantification-peak-fitting
date: 2026-10-22 20:00:00 +0900
categories: [Surface Analysis, Spectrum Fitting]
tags: [xps, quantification, peak-fitting, shirley-background, voigt, charge-referencing]
description: "Between a peak area and a composition sit sensitivity factors, a background model, a line shape and charge referencing. This post shows where each one breaks and fits a synthetic spectrum by hand."
math: true
---

[Part 3](/en/posts/xps-spectrum-structure-satellites/) sorted out where each structure in a spectrum comes from and set criteria for what should and should not be counted as a component. Now it is time to take those criteria and extract numbers, the kind that end up in a report as "35% oxide."

That number is not fixed the moment a peak area is measured. First the energy scale has to be referenced, then one has to decide where and in what shape to draw the background, which functions to split overlapping peaks into, and which sensitivity factors to use in converting areas to atom counts. A single number in a composition table comes out only after passing through four models: charge referencing, background, line shape and sensitivity factors. The refrain of this series since Part 1, that a model sits between what is measured and what we want to know, is nowhere clearer than here.

This post goes through the four models in turn and, at the end, fits a synthetic spectrum with a known answer to see what can and cannot be recovered.

## 1. Does a Bigger Peak Mean More of That Element?

How far does the intuition "a bigger peak means more of that element" hold? The area $I_x$ of a photoelectron peak from a given core level of element $x$ can be written roughly as a product:

$$ I_x = A \, n_x \, f \, \sigma_x \, \theta \, \lambda(E_K) \, T(E_K) $$

$A$ is the analysis area, $n_x$ the atomic density, $f$ the X-ray flux, $\sigma_x$ the photoionization cross-section of that level, $\theta$ an angular factor from the X-ray and analyser geometry, $\lambda$ the inelastic mean free path from Part 2, and $T$ the transmission function of the spectrometer. Photoionization cross-sections differ by up to a few hundred times between elements and levels; the values calculated by Scofield (1976) at the Al K$\alpha$ and Mg K$\alpha$ energies are the standard.

Within one measurement, $A$ and $f$ are common, so they cancel in the ratio of two elements. The angular factor $\theta$ differs between levels, but on instruments with the X-ray source and analyser at the magic angle (54.7°) it is the same for every level and cancels as well:

$$ \frac{n_1}{n_2} = \frac{I_1 / (\sigma_1 \lambda_1 T_1)}{I_2 / (\sigma_2 \lambda_2 T_2)} $$

The $\lambda$ in this expression deserves attention. As [Part 2](/en/posts/xps-surface-sensitivity-imfp/) showed, $\lambda$ depends on kinetic energy, so each peak probes a different depth. And the expression assumes the sample is uniform within that depth. If the surface carries an adventitious carbon layer or a thin oxide covers a metal, the assumption fails, and the ratio above is not a physical composition but an apparent one weighted exponentially in depth. Dropping this assumption and using the attenuation in reverse is the thickness measurement of Part 5.

So the answer to the opening question is this. Peak sizes can be compared as they stand only for the same level of the same element, measured under the same conditions. Between different elements the cross-sections alone differ by hundreds of times, so the areas must be divided by sensitivities, and even then the ratio is a composition only if the sample is uniform in depth.

## 2. Why Is There More Than One Sensitivity Factor Table?

Why is the expression above not used directly in practice? Instead of entering $\sigma$, $\lambda$, $T$ and the angular factor one by one, analysts use a relative sensitivity factor (RSF) that lumps them together. Analysis software carries an RSF table for each element and level and computes composition by dividing peak areas by those values. To give the conclusion first: because this lumped factor differs between tables and instruments, XPS quantification is trustworthy only when comparing samples measured under the same conditions.

The trouble is that there is more than one table, built in three ways: tables determined empirically from standard samples, tables built from Scofield cross-sections multiplied by calculated $\lambda$ and transmission corrections, and tables that instrument makers ship with their analysis software. The values in the three differ. That is why quantifying the same spectrum in different software gives different compositions.

The fact that the transmission function differs between instruments is especially important. Part 1 noted that survey scans use a high pass energy and high-resolution scans a low one, and deferred the reason to here. Most XPS spectrometers use a concentric hemispherical analyser (CHA). Electrons from the sample are first retarded by a lens to a set pass energy $E_p$ and then sent between two hemispheres held at a potential difference. Only electrons with kinetic energy near $E_p$ make it through to the exit slit. With the inner hemisphere positive and the outer one negative, electrons are pulled inwards and follow an arc. An electron with more kinetic energy than $E_p$ swings round on a larger radius and lands towards the outer hemisphere; one with less lands towards the inner one. In other words, the two hemispheres act as a prism that spreads electrons out by energy (Figure 1). Scanning the retarding voltage while keeping $E_p$ fixed is called FAT (fixed analyser transmission) mode, the standard mode for XPS.

<img src="/assets/img/posts/xps-quantification-peak-fitting/ext-hemispherical-analyzer-principle.png" alt="Cutaway of a hemispherical analyser: electrons entering the slit bend between the two hemispheres and are spread out by kinetic energy on the detector" width="560">
_Figure 1. The structure of a hemispherical analyser. Electrons entering the slit at the lower right bend between the two hemispheres and reach the detector at the lower left. On the detector, position along the slit (x) is preserved, while kinetic energy ($$E_k$$) is spread out in the perpendicular direction. The three coloured bands are groups of electrons with different kinetic energies, and the width of each band along the slit direction ($$x$$) is the range of positions at the entrance. (Source: Ponor, [Hemispherical electron energy analyzer working principle](https://commons.wikimedia.org/wiki/File:Hemispherical_electron_energy_analyzer_working_principle.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

The energy resolution of this analyser is proportional to the pass energy. With entrance slit width $W$ and mean hemisphere radius $R$, the resolution is roughly $\Delta E \approx E_p \, W/(2R)$ plus a term from the angular spread of incoming electrons. Even electrons of the same energy land at different points on the exit if they enter at slightly different angles (Figure 2). The angular term appears with different coefficients in different sources, such as $\alpha^2$ and $\alpha^2/2$, depending on whether the resolution is defined as a base width or an FWHM and whether $\alpha$ is taken as a half-angle or a full angle; check the definitions of $\Delta E$ and $\alpha$ before borrowing the formula. Either way the conclusion is the same. Lowering $E_p$ improves resolution but lets fewer electrons through the analyser, weakening the signal.

In Thermo Fisher's analyser literature, the same Ag 3d peak measured at pass energies from 4 to 40 eV loses height as the linewidth narrows. The same note gives the narrowest Ag 3d₅/₂ FWHM as below 0.82 eV with a non-monochromated Mg K$\alpha$ source and below 0.46 eV with a monochromated one. The Si⁰ Gaussian width of 0.55 eV in the synthetic spectrum of Section 7 falls within this range.

In FAT mode the resolution is constant across the spectrum, but the transmission function depends on the retarding ratio $E_K/E_p$. So even on the same instrument, changing the pass energy changes the sensitivity.

<img src="/assets/img/posts/xps-quantification-peak-fitting/ext-hemispherical-analyzer-trajectories.png" alt="Calculated trajectories of electrons with five kinetic energies and five entrance angles inside a hemispherical analyser" width="680">
_Figure 2. Calculated electron trajectories in a hemispherical analyser. At the entrance slit on the lower left, five groups of electrons with different kinetic energies (yellow to red) each enter at five angles, bend between the inner hemisphere ($$V_1$$, +) and the outer one ($$V_2$$, −), and reach the exit on the right. $\alpha$ is the spread of entrance angles. Electrons landing further out have more kinetic energy. Trajectories of the same colour also spread at the exit according to their angle and overlap with neighbouring colours; this overlap is the angular term in the resolution formula. (Source: Ponor, [Hemispherical electron energy analyzer - calculated trajectories for different energies and initial angles](https://commons.wikimedia.org/wiki/File:Hemispherical_electron_energy_analyzer_-_calculated_trajectories_for_different_energies_and_initial_angles.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

The practical conclusion: XPS quantification is strong for relative comparisons and weak for absolute values. Trends among samples measured on the same instrument under the same conditions are trustworthy, but comparing decimals against literature values obtained on other instruments with other RSF tables is meaningless.

## 3. What Is the Energy Scale of an Insulating Sample Referenced To?

What is the energy axis of an insulating sample referenced to? Section 1 of [Part 1](/en/posts/xps-photoemission-binding-energy/) explained that for a conductor, the binding-energy scale is fixed because the Fermi levels of sample and spectrometer align. In an insulator, electrons cannot refill the sites left by photoelectrons, the surface charges positively, and every peak shifts bodily to higher binding energy. A charge neutralizer that floods the surface with low-energy electrons or ions prevents this, but can overcompensate and shift the peaks the other way. Either way the measured scale is itself off, so a separate reference is needed. The most common answer is the convention of setting C 1s to 284.8 eV, and that convention is still under debate. To give the conclusion first: even when this value is used, the report must state which reference was used and what value it was set to.

The convention uses adventitious carbon. As Part 2 noted, nearly every air-exposed sample is covered with a hydrocarbon layer, so the C–C/C–H component of its C 1s is set to 284.8 eV and the whole spectrum is shifted accordingly. Some, for example for polymers, use 285.0 eV.

This convention has been debated in recent years. One side regards the reference as unstable in principle. Greczynski (2024) measured 360 thin-film samples, including metals, nitrides, carbides and oxides, and reported that the C 1s position of adventitious carbon moves over a range of about 3 eV (284.1–286.9 eV) with the sample's work function. The interpretation is that adventitious carbon aligns not with the sample's Fermi level but with the vacuum level, in which case referencing to it requires measuring the sample's work function separately. Earlier, Greczynski and Hultman (2022) had examined separately the objection that these shifts might come from differential charging of the sample.

The other side responds with practical data. Following Biesinger's (2022) review of data from one multi-user facility, Morgan (2025) reviewed five years of data from a second facility. The conclusion was that regardless of spectrometer or charge-compensation method, 284.8 eV is generally suitable for electronically isolated samples, with a few exceptions where a secondary reference is preferable.

It is worth noting that the two conclusions come from different sample populations. Greczynski systematically compared films with widely different work functions, while Morgan's conclusion is restricted to electronically isolated samples. As long as the debate is unsettled, two things can be done: always report which reference was used and what value it was set to, and where possible cross-check against another reference. Alternatives include the Fermi edge of a conductor, an internal reference peak whose component is known for certain, and the modified Auger parameter from [Part 3](/en/posts/xps-spectrum-structure-satellites/). The Auger parameter is independent of charging and sidesteps the debate altogether. Part 1 noted that the binding-energy scale is not itself a measurement but a scale built on an assumption; for insulators even its zero point is the product of a choice of reference.

## 4. Where, and in What Shape, Is the Background Drawn Under a Peak?

To measure a peak area, the background must first be subtracted, yet the background is never measured. Only the sum of peak and background is. The background can only be drawn from a model, and the area depends on that model.

The simplest, a linear background, joins the two ends of the analysis window with a straight line. It is simple but ignores the stepped structure of Section 7 in Part 3. In spectra with a large step, it draws the background too high or too low under the peak.

The Shirley (1972) background accounts for that step. The background climbs like a staircase towards higher binding energy because some of the electrons leaving a peak lose energy inside the sample and pile up behind it, on the higher-binding-energy side (Section 7 of Part 3). Shirley makes the height of the step proportional to "the peak area accumulated on the lower-binding-energy side of that point." On a binding-energy axis, the background height at an energy $E$ is taken to be proportional to the peak area on its lower-binding-energy side:

$$ B(E) = I_{\text{low}} + (I_{\text{high}} - I_{\text{low}}) \frac{\int_{E_{\text{low}}}^{E} \left[ I(E') - B(E') \right] dE'}{\int_{E_{\text{low}}}^{E_{\text{high}}} \left[ I(E') - B(E') \right] dE'} $$

$I_{\text{low}}$ and $I_{\text{high}}$ are the intensities at the two ends of the window. In kinetic-energy terms, "the background is proportional to the peak area at higher kinetic energy." Because the background $B$ appears again on the right-hand side, the form is self-referential: start from a constant background and iterate until it converges. Shirley's original paper computed it only once; the present practice of repeating it to convergence was set out by Proctor and Sherwood (1982). It is simple to implement and the most widely used.

The Tougaard background goes a step closer to the physics. It writes the probability of an electron losing energy $T$ in the solid directly as an inelastic-scattering cross-section, and computes the distribution of loss electrons by convolving the spectrum with that cross-section. The most common choice is the two-constant universal cross-section:

$$ F(T) = \frac{B\,T}{(C + T^2)^2} $$

The universal values fitted for metals are $B$ = 2866 eV² and $C$ = 1643 eV², and in practice $B$ is adjusted so that the background meets the data at the end of the window. It is the most rigorous physically and can even yield depth-distribution information from the background, but it needs a wide energy window and is tricky to use. It comes back in Part 5.

<img src="/assets/img/posts/xps-quantification-peak-fitting/en/fig1-backgrounds.png" alt="Linear, Shirley and Tougaard backgrounds applied to the same synthetic Si 2p spectrum, with area errors" width="680">
_Figure 3. Peak area by background model (synthetic spectrum). The true background was generated with the Tougaard universal cross-section, which favours the Tougaard model. The vertical axis is clipped at 300 to show the background. The linear background underestimates the area by 10.2% and the Shirley background by 12.7%; only Tougaard, which has the same shape as the true background, gets within 1.0%._

Figure 3 applies the three backgrounds to the synthetic Si 2p spectrum used in Section 7. Since the true background was generated as a Tougaard-type loss, it is no surprise that the Tougaard model matches the true area to within 1.0%. The other two are what matter. The linear background underestimates the area by 10.2%, and Shirley by 12.7%. The Shirley background raises its step as soon as a peak is passed, whereas real losses accumulate slowly over tens of eV.

The size of this difference, however, depends strongly on conditions. Figure 3 sets the background step at the end of the window to 6.5% of the peak height. In the same calculation with the universal value as is (a step of 1.1% of the peak), the differences among the three models shrink to within 2%; with a step of 3.2% they range from −7% to −1% (means over 20 noise realizations; in a single realization the noise alone moves the area by about ±2%). The problem is that in a real spectrum the true background is unknown, so one cannot tell which case applies. That is why comparisons of peak areas between samples must use the same window and the same background model. Where the window ends are placed alone can move an area by several percent.

## 5. Why Do Peaks Have the Shape They Do?

Why do peaks have the shape they do? The observed shape mixes two contributions of entirely different character.

The first is the Lorentzian. A core hole is soon filled by an Auger process or fluorescence, so its lifetime $\tau$ is finite. By the uncertainty relation, the shorter the lifetime the more the energy spreads, with a width of roughly $\Gamma \sim \hbar/\tau$. The synthetic spectrum of Section 7 sets this width to 0.10 eV for Si 2p, which corresponds to a lifetime of $\hbar/\Gamma \approx 6.6$ fs. In other words, the core hole is filled within a few femtoseconds, and that short lifetime shows up as a width of 0.1 eV. The spectrum of a state decaying exponentially in time has a Lorentzian shape. This width is intrinsic to the element and level. The extra width of Ti 2p₁/₂ from the Coster–Kronig transition in Section 1 of Part 3 is also an effect of this lifetime width.

The second is the Gaussian. Instrumental contributions such as the linewidth of monochromated X-rays and analyser resolution, and sample contributions such as phonon broadening, nonuniform surface charging and fine inhomogeneity of the chemical environment, all go here. With many independent causes, their combination tends toward a Gaussian. A non-monochromated K$\alpha$ line is itself Lorentzian in shape, but by convention all instrumental width is put on the Gaussian side.

How do the two broadenings combine? Think of the instrument and sample smearing each point of the lifetime-broadened Lorentzian peak into a Gaussian again. Stamping a Gaussian at every point and adding them up, weighted by intensity, is convolution. It is like taking a blurred photograph through a blurred lens once more. In other words, because the two processes are independent, the observed line shape is neither their product nor their sum but their **convolution**, and that is the Voigt function.

<img src="/assets/img/posts/xps-quantification-peak-fitting/en/fig2-voigt-composition.png" alt="Gaussian and Lorentzian with 1 eV FWHM and their convolution, the Voigt, on linear and log axes" width="720">
_Figure 4. A Gaussian (FWHM 1.00 eV), a Lorentzian (FWHM 1.00 eV) and their convolution, the Voigt function, all normalized to the same area. (a) Linear axis; (b) log axis. Arrows mark the FWHM. The Voigt function broadens to an FWHM of 1.64 eV, and on the log axis only the Gaussian tail vanishes within 2 eV of the centre._

Figure 4(a) shows that at equal width the Gaussian rises sharply while the Lorentzian is low and broad. The Voigt function obtained by convolving them has an FWHM of 1.64 eV, less than the simple sum of the two widths (2 eV). The approximation $f_V \approx 0.5346\, f_L + \sqrt{0.2166\, f_L^2 + f_G^2}$ (Olivero and Longbothum, 1977), worth remembering, also gives 1.64 eV for two widths of 1 eV. Here $f_G$, $f_L$ and $f_V$ are the FWHMs of the Gaussian, Lorentzian and Voigt functions. The difference is clearest in the tails on the log axis of Figure 4(b). The Gaussian is effectively zero only 2 eV from the centre, while the Lorentzian and Voigt functions keep 0.7–1% of their own peak height even 6 eV out. These tails are hard to tell apart from the background, so the choices of line shape and background become entangled.

The Voigt convolution is expensive to compute, so software often uses approximations. A linear combination of a Gaussian and a Lorentzian is called a pseudo-Voigt function. Beware that notation differs between programs. In the widely used CasaXPS, `GL(m)` is the **product** of a Gaussian and a Lorentzian, and the linear combination has a separate name, `SGL(m)`. A paper that says only "fitted with GL(30)" is ambiguous unless you know which software was used.

Metals add asymmetry on top of this. They need the Doniach–Šunjić line shape of Section 4 in Part 3, which also showed that using a symmetric function on a metal peak creates a spurious oxide component to fill the tail. Knowing where each width comes from tells you what to tie together and what to leave free when setting constraints: tie the lifetime width across the same level, tie the instrumental width within one measurement, and let only the sample broadening vary between chemical states.

## 6. Does a Small Residual Mean a Correct Decomposition?

Is a decomposition correct if its residual is small? No. Adding a component adds free parameters, so the residual almost always falls. The extra component fills in noise, a misfit in the background, or anything else. The residual or $\chi^2$ alone therefore cannot determine the number of components.

What determines the number and shape of components is physics, not the data. Constraints from physics must be imposed:

- Tie the area ratio of a spin-orbit doublet to the value from $2j+1$ (p 1:2, d 2:3, f 3:4) and the splitting to the literature value. For elements whose splitting varies slightly with chemical state, or whose measured area ratio deviates a little, as seen in Part 3, use that element's literature values from Part 3.
- Give the two components of a doublet from the same species the same line shape and width. Where only one side broadens, as with Coster–Kronig in Ti 2p, set the width ratio separately.
- Restrict component positions to a range around literature values.
- For elements showing multiplet splitting, use a multiplet model fixed from standard samples.

Constraints reduce the number of free parameters and so increase the residual. In return the solution becomes stable and physically meaningful. This trade-off is the heart of fitting. Reducing the residual is not the goal in itself; the goal is to choose, among physically sensible models, one that explains the data.

Mathematically, peak fitting is the nonlinear least-squares problem of [Optimization 3](/en/posts/optimization-levenberg-marquardt/), and fitting software such as CasaXPS solves it with the Levenberg–Marquardt method or a variant. Converging to different solutions from different initial component positions is the same local-minimum problem as in [Optimization 4](/en/posts/optimization-global-heuristics/). The more components and the fewer constraints, the more nearly equally good local minima there are. Constraints are also what restore uniqueness to the solution.

## 7. Does Fitting a Spectrum with a Known Answer Give That Answer Back?

If we build a spectrum with a known answer and fit it, are the parameters recovered correctly? Let us check the discussion so far with numbers. The code is in [`_code/xps-quantification-peak-fitting/`](https://github.com/eddykim/eddykim.github.io/tree/main/_code/xps-quantification-peak-fitting), and running `xps_fit.py` prints every number in this section.

### 7.1 Building the synthetic spectrum

To lead into Part 5, the spectrum mimics the Si 2p of silicon covered with a thin SiO₂ layer. There are two components. Si⁰ has its 2p₃/₂ at 99.40 eV, a total doublet area of 1000 and a Gaussian FWHM of 0.55 eV; Si⁴⁺ of SiO₂ sits at 103.30 eV with area 600 and FWHM 1.30 eV. The true Si⁴⁺ share is 37.5%. Both are spin-orbit doublets, with a splitting of 0.61 eV and an intensity ratio of 1/2 taken from the SiO₂/Si interface analysis of Lu et al. (1995). The background is a loss background generated from the Tougaard universal cross-section of the previous section, scaled so the step at the end of the window is 6.5% of the peak height. Finally, Poisson noise was added to the counts.

### 7.2 Implementing the Shirley background

The iteration of Section 4 translates directly. With the binding-energy axis in ascending order, the cumulative area from the low-binding-energy end is the numerator of the equation.

```python
def shirley(x, y, n_end=5, tol=1e-7, max_iter=100):
    """x is binding energy, ascending. background = low end + (high - low) x cumulative area ratio."""
    lo, hi = y[:n_end].mean(), y[-n_end:].mean()
    bg = np.full_like(y, lo)
    for it in range(max_iter):
        cum = _cumarea(x, np.clip(y - bg, 0, None))   # cumulative area from the low-BE end
        new = lo + (hi - lo) * cum / cum[-1]
        if np.max(np.abs(new - bg)) < tol * hi:
            return new, it + 1
        bg = new
    return bg, max_iter
```

For this spectrum it converged in seven iterations.

### 7.3 Fitting with constraints

The constraints are imposed through parameterization. Deriving 2p₁/₂ from 2p₃/₂ instead of treating it as an independent component fixes the area ratio and splitting automatically. Each species has three parameters: position, area and width.

```python
def doublet(x, c32, area, g_fwhm):
    """Spin-orbit doublet. 2p1/2 is derived from 2p3/2, not an independent parameter."""
    a32 = area / (1 + SO_RATIO)
    return voigt(x, c32, a32, g_fwhm) + voigt(x, c32 + SO_SPLIT, a32 * SO_RATIO, g_fwhm)

def resid(p):   # p = [Si0 position, area, width, Si4+ position, area, width]
    return (doublet(x, *p[:3]) + doublet(x, *p[3:]) - y_sub) / np.sqrt(y_raw)
```

The residuals are weighted by the Poisson standard deviation. The sum of squared residuals is then $\chi^2$, and $\chi^2_{\text{red}}$, $\chi^2$ divided by the degrees of freedom, near 1 means the model explains the data down to the noise level. The minimization used `scipy.optimize.least_squares`.

<img src="/assets/img/posts/xps-quantification-peak-fitting/en/fig3-fit-result.png" alt="Synthetic Si 2p spectrum fitted with the constrained doublet model after Shirley background subtraction, normalized residual and a table of true and fitted values" width="680">
_Figure 5. Constrained fit (synthetic spectrum). Top: data, Shirley background, doublets per species (dark 2p₃/₂, light 2p₁/₂) and fit sum. Bottom: residuals divided by the Poisson standard deviation. The positions are right to within 0.01 eV, but the Si⁴⁺ area comes out low at 497 (true 600), and the residual keeps negative dips at 101–102 eV and 105–107 eV ($$\chi^2_{\text{red}}$$ 2.57)._

The result is a half success. Both positions are correct to within 0.01 eV. The areas, however, come out as 954 ± 12 for Si⁰ (true 1000) and 497 ± 11 for Si⁴⁺ (true 600), and the Si⁴⁺ width comes out narrow at 1.06 ± 0.04 eV (true 1.30). The Si⁴⁺ share is 34.3%, 3.2 percentage points below the true value. $\chi^2_{\text{red}}$ is 2.57, clearly above 1, and the residual at the bottom of Figure 5 has systematic negative dips around 101–102 eV and 105–107 eV.

The cause is the background, not the fit. Subtracting the true background exactly from the same data and fitting the same constrained model drops $\chi^2_{\text{red}}$ to 0.99 and recovers the areas within their errors, 1003.5 ± 7.6 for Si⁰ and 598.5 ± 6.8 for Si⁴⁺, with a Si⁴⁺ share of 37.4%. The fitting machinery works; the mismatch of the Shirley background seen in Section 4 turns directly into an area bias.

The most important observation here concerns the error bars. The Si⁴⁺ area from the Shirley background is 497 ± 11, but the actual error is −103. The error a fit reports is only the statistical error from noise and contains nothing of the error from choosing the wrong background model. A $\chi^2_{\text{red}}$ above 1 is the only warning sign. When a model's assumption is wrong, the number still comes out wearing a plausible error bar.

### 7.4 What happens when the constraints are released?

What happens if the spin-orbit constraint is released on the same data? Instead of doublets, the fit uses $n$ single Voigt functions with free position, area and width, choosing the best of several starting points. Then, as an analyst would, each component is labelled with the nearest chemical state in the chemical-shift table of Lu et al. (1995): relative to Si⁰, Si¹⁺ +0.98, Si²⁺ +1.82, Si³⁺ +2.65 and Si⁴⁺ +3.84 eV.

With two singlets, $$\chi^2_{\text{red}}$$ is 4.91, clearly inadequate. With three singlets, $$\chi^2_{\text{red}}$$ falls to 2.55, practically the same as the 2.57 of the correct constrained model. But the three components sit at 99.42, 100.04 and 103.49 eV. The middle one is the 2p₁/₂ of Si⁰ split off on its own; it lies +0.62 eV from Si⁰, closest to Si¹⁺ in the chemical-shift table. An analyst checking the result against the literature table would report "Si⁰ 48.7%, Si¹⁺ 17.0%, Si⁴⁺ 34.3%." A nonexistent intermediate oxide has appeared at 17%, and yet the residual cannot tell this model from the correct one.

<img src="/assets/img/posts/xps-quantification-peak-fitting/en/fig4-overfitting.png" alt="Decrease of χ² and the spurious intermediate-oxide and Si4+ shares as the number of unconstrained singlets increases" width="720">
_Figure 6. Adding components (synthetic spectra; mean and standard deviation over 20 noise realizations). (a) χ² keeps falling with the number of components and already reaches the level of the constrained model (dashed) at n = 3. (b) Meanwhile spurious components reported as Si¹⁺–Si³⁺ appear, and the spread of the results grows with n._

Figure 6 repeats this experiment 20 times with fresh noise. As the number of components goes from 2 to 7, the mean $\chi^2$ keeps falling: 1322, 709, 671, 656, 646, 640. The correct constrained model has a $\chi^2$ of 718, so with as few as three components the wrong model reaches the same level on residuals.

The residual cannot tell the correct and wrong models apart. Paired realization by realization, the three-component $\chi^2$ is smaller in 15 of 20 cases. But the difference averages 9.7, about what three extra parameters buy, and the difference in $\chi^2_{\text{red}}$, after dividing by the degrees of freedom, is −0.005 ± 0.040, indistinguishable from zero.

Statistical criteria do no better. Even the BIC (Bayesian information criterion), which penalizes $k$ parameters by $k \ln N$, picks the wrong three-component model in 6 of 20 cases. Statistical criteria do not replace physical constraints.

Adding components creates spurious species and destabilizes the solution. As components are added, the spurious intermediate-oxide share jumps from 0% to 18.1%; with more components the mean eases slightly to 15.5%, but the spread between noise realizations widens from ±0.8 to ±5.1 percentage points. With surplus components, the data cannot decide where to put them, and the solution differs from one realization to the next. That spread is what Section 6 meant by "the fit is not unique."

The Si⁴⁺ share, though, stays near 34% almost regardless of the number of components. In this synthetic sample Si⁰ and Si⁴⁺ are 3.9 eV apart, so however the components are divided, the area on the oxide side barely moves. What goes wrong without constraints is not the total amount of oxide but "which species are present." And the oxide share staying about 3 percentage points below the true 37.5% comes from the Shirley background bias of Section 7.3, not from the number of components. In short, in this synthetic sample the number of components gets "which species" wrong, and the background gets "how much" wrong.

## 8. What Should Be Reported Along with a Composition Number?

Misreading a fit is not limited to individual mistakes. Major et al. (2020) had an expert committee review 409 XPS analyses published in three journals over six months and found that about 30% were fundamentally flawed. More than 65% of the papers included peak fitting, and fitting was a common source of serious errors. That is one of the reasons behind the series of practical guides the XPS community has been issuing recently.

Following the four models covered here, the items to report alongside a composition number are:

- Charge referencing: which reference, set to what value, and whether a neutralizer was used (items ISO 19318 requires to be reported)
- Background: which model, applied over which window
- Line shape: which function (with the software name and its notation), and whether asymmetric shapes were used for metals
- Constraints: how doublet area ratios and splittings, width ties and position ranges were imposed
- Sensitivity factors: which RSF table, and how the transmission function was corrected
- Fit quality: $\chi^2_{\text{red}}$ or a residual plot

A composition number that cannot fill in this list cannot be reproduced by anyone else. As the experiment in Section 7 shows, changing the background alone moves the Si⁴⁺ share by about 3 percentage points, and releasing the constraint on the number of components alone makes an intermediate oxide that is not there appear at up to 17%.

## Summary and Next Part

This post followed the four models a peak area passes through on its way to becoming a composition number. Charge referencing rebuilds the energy scale of an insulator, and the most common reference, adventitious carbon at 284.8 eV, is still debated. Background models change areas by several to more than ten percent depending on what the true background is. The line shape is the Voigt function, the convolution of a lifetime Lorentzian and an instrumental and sample Gaussian, and software notations for its approximations differ. Sensitivity factors vary with instrument and table, making XPS stronger for relative comparisons than for absolute composition. The synthetic-spectrum experiment confirmed two things: with the wrong background model, the error bars say nothing about the actual error, and releasing the spin-orbit constraint reports 17% of a species that does not exist while the residual cannot tell the difference.

So far the sample has been assumed uniform in depth within the information depth. The final Part 5 drops that assumption. The SiO₂/Si used in Part 1 is exactly such a sample. Using the attenuation of Part 2 in reverse lets one measure the oxide thickness, but thickness is the most indirect quantity XPS yields. Part 5 covers the inverse problem of angle-resolved XPS, preferential sputtering where ion sputtering alters the sample's composition, and X-ray damage: the pitfalls in which the measurement changes what is being measured.

## References

- J. H. Scofield, "Hartree-Slater subshell photoionization cross-sections at 1254 and 1487 eV," *J. Electron Spectrosc. Relat. Phenom.* 8, 129–137 (1976). <https://doi.org/10.1016/0368-2048(76)80015-1>
- D. A. Shirley, "High-resolution X-ray photoemission spectrum of the valence bands of gold," *Phys. Rev. B* 5, 4709–4714 (1972). <https://doi.org/10.1103/PhysRevB.5.4709>
- S. Tougaard, "Universality classes of inelastic electron scattering cross-sections," *Surf. Interface Anal.* 25, 137–154 (1997). <https://doi.org/10.1002/(SICI)1096-9918(199703)25:3%3C137::AID-SIA230%3E3.0.CO;2-L>
- Z. H. Lu, S. P. Tay, T. Miller, and T.-C. Chiang, "Process dependence of the SiO₂/Si(100) interface structure," *J. Appl. Phys.* 77, 4110–4112 (1995). <https://doi.org/10.1063/1.359494>
- G. Greczynski, "Binding energy referencing in X-ray photoelectron spectroscopy: expanded data set confirms that adventitious carbon aligns to the vacuum level," *Appl. Surf. Sci.* 670, 160666 (2024). <https://arxiv.org/abs/2405.10919>
- G. Greczynski and L. Hultman, "Referencing to adventitious carbon in X-ray photoelectron spectroscopy: Can differential charging explain C 1s peak shifts?," *Appl. Surf. Sci.* 606, 154855 (2022). <https://doi.org/10.1016/j.apsusc.2022.154855>
- A. Proctor and P. M. A. Sherwood, "Data analysis techniques in X-ray photoelectron spectroscopy," *Anal. Chem.* 54, 13–19 (1982). <https://doi.org/10.1021/ac00238a008>
- J. J. Olivero and R. L. Longbothum, "Empirical fits to the Voigt line width: A brief review," *J. Quant. Spectrosc. Radiat. Transf.* 17, 233–236 (1977). <https://doi.org/10.1016/0022-4073(77)90161-3>
- D. J. Morgan, "The Utility of Adventitious Carbon for Charge Correction: A Perspective From a Second Multiuser Facility," *Surf. Interface Anal.* 57, 28–35 (2025). <https://doi.org/10.1002/sia.7360>
- M. C. Biesinger, "Accessing the robustness of adventitious carbon for charge referencing (correction) purposes in XPS analysis: Insights from a multi-user facility data review," *Appl. Surf. Sci.* 597, 153681 (2022). <https://doi.org/10.1016/j.apsusc.2022.153681>
- G. H. Major et al., "Practical guide for curve fitting in x-ray photoelectron spectroscopy," *J. Vac. Sci. Technol. A* 38, 061203 (2020). <https://doi.org/10.1116/6.0000377>
- G. H. Major et al., "Assessment of the frequency and nature of erroneous x-ray photoelectron spectroscopy analyses in the scientific literature," *J. Vac. Sci. Technol. A* 38, 061204 (2020). <https://pubs.aip.org/avs/jva/article/38/6/061204/1023649/Assessment-of-the-frequency-and-nature-of>
- Casa Software Ltd., "Peak Fitting in XPS," 2006. <https://mmrc.caltech.edu/XPS%20software/CASA/Books/peak_fitting_in_xps.pdf>
- Thermo Fisher Scientific, "Alpha110: Hemispherical Electron Energy Analyzer," application note. <https://documents.thermofisher.com/TFS-Assets/CAD/Application-Notes/D16093~.pdf>

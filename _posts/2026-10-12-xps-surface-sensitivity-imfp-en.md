---
title: "XPS Fundamentals 2 — X-rays Go Micrometres Deep, So Why Does XPS See Only the Top Few Nanometres?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/xps-surface-sensitivity-imfp/
page_id: xps-surface-sensitivity-imfp
date: 2026-10-12 20:00:00 +0900
categories: [Surface Analysis, Depth Profiling]
tags: [xps, imfp, surface-analysis, information-depth, arxps]
description: "X-rays penetrate micrometres into a sample, so why is XPS a top-10-nm technique? The inelastic mean free path and the information depth."
math: true
---

[Part 1](/en/posts/xps-photoemission-binding-energy/) got as far as the energy-conservation equation $E_K = h\nu - E_B - \phi_{sp}$, which yields the binding energy, and the fact that core-level binding energies fingerprint the elements. It did not explain why XPS is called a *surface* technique. Depth appears nowhere in Part 1's equations.

The puzzle starts with the X-rays. Al K$\alpha$ X-rays (1486.6 eV) travel several micrometres into silicon. Photoemission happens wherever the X-rays reach, so photoelectrons are generated throughout that depth. Yet the peaks in an XPS spectrum are made of signal from the top few nanometres. Surface sensitivity does not come from the X-rays. Where does it come from?

This post answers that question and takes apart a familiar sentence: "XPS sees the top 10 nm." The short version is that 10 nm is not a physical constant. It is the product of an exponential attenuation model, a conventional 95% criterion, a mean free path computed from a predictive formula, and a normal-emission geometry. Change any one of the four and the number changes.

## 1. If the X-rays Go Deep, Why Does the Signal Come Only from the Surface?

To answer this, compare how far an X-ray and a photoelectron each travel in a solid. The database of the Center for X-Ray Optics (CXRO) at Lawrence Berkeley National Laboratory gives attenuation lengths for Al K$\alpha$ photons of 7.9 µm in silicon, 4.1 µm in SiO₂ and 0.23 µm in gold. How far a photoelectron of similar energy travels without losing energy is measured by the inelastic mean free path (IMFP), defined in Section 2 and calculated in Section 3. The predictive formula of Section 3 gives about 3.1 nm for Si 2p photoelectrons in silicon and about 1.6 nm for Au 4f photoelectrons in gold. The ratios are about 2,500 in silicon and about 150 in gold, a difference of two to three orders of magnitude.

Here is what that means. Over the first 10 nm below a silicon surface, the X-ray intensity drops by only 0.13%. Within the depth from which any signal can emerge, photoelectrons are therefore generated almost uniformly. The problem is the way out. A photoelectron starting deep in the sample collides with electrons in the solid and undergoes inelastic scattering. It excites plasmons or valence electrons and loses kinetic energy bit by bit. Of the Si 2p photoelectrons that start 10 nm down, only 3.9% leave the surface without any loss. From 1 nm the figure is 72%, from 3 nm 38%, and from 5 nm 20%.

An electron that has lost energy can still leave the surface and reach the analyser. But its kinetic energy has changed, so it cannot appear at the original peak position. The stepped background on the high-binding-energy side of every peak in the survey spectrum of Part 1 is made of these electrons.

Think of a foggy day. Streetlights reach far, and light returning from distant people still reaches your eye. But light scattered many times in the fog loses all shape and becomes a grey haze. You can recognise only the faces within visibility range. In XPS the lighting is the X-ray beam, the visibility range is the IMFP, and the haze is the stepped background.

So the only electrons that contribute to a peak area are those that **escaped without a single inelastic collision**. That is what surface sensitivity is. Information from deep in the sample is not destroyed; it moves from the peak into the background. This view comes back in Part 4, when we ask what shape of background to subtract.

## 2. What Does the IMFP (Inelastic Mean Free Path) Measure?

How do we quantify "the probability of escaping without loss"? If inelastic collisions are independent events occurring with the same probability per unit path length, the probability of travelling a distance $s$ without any collision is exponential.

$$ P(s) = e^{-s/\lambda} $$

Here $\lambda$ is the IMFP introduced in Section 1. The international standard ISO 18115 defines it as the average distance an electron of a given energy travels between successive inelastic collisions. In terms of the equation, it is the distance over which the no-loss probability falls to $1/e$. Every distance $\lambda$ an electron travels, the fraction that has not yet lost energy shrinks to 37%.

If a photoelectron created at depth $z$ heads for the analyser at an angle $\theta$ from the surface normal, its path in the solid is $z/\cos\theta$. The probability that it escapes without loss is therefore

$$ P(z,\theta) = \exp\!\left(-\frac{z}{\lambda\cos\theta}\right) $$

Throughout this post, $\theta$ is always **measured from the surface normal**. Section 5 explains why that convention matters.

The equation has the same form as the Beer–Lambert law from [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/), but the thing being attenuated is different. Beer–Lambert describes photons being absorbed and disappearing. Here the electrons do not disappear; they only lose energy.

In practice three more quantities with similar names are in use, and mixing up the four is a common source of error. Following Powell (2020), they separate as follows.

| Quantity | Definition | Depends on |
|---|---|---|
| IMFP $\lambda$ | Mean distance between inelastic collisions | Material, electron energy |
| Effective attenuation length (EAL) | The value that, used in place of the IMFP in an expression that ignores elastic scattering, corrects that expression for elastic scattering | The above + geometry, application |
| Mean escape depth (MED) | Average depth from which detected electrons originate | Same as above |
| Information depth (ID) | Depth from which a specified fraction (95% or 99%) of the signal originates | Same as above |

Only the IMFP is set by the material and the electron energy alone; the other three also bring in the measurement conditions. The difference comes from elastic scattering. When an electron is deflected near a nucleus without losing energy, its actual path becomes longer than a straight line, and an electron starting at a given depth has more chances to scatter inelastically. The EAL is therefore shorter than the IMFP. In the calculations Powell summarises, the ratio EAL/IMFP over 200 eV–1.5 keV, the range used most in XPS, is 0.77–0.92 for aluminium and silicon and 0.68–0.85 for copper, silver, indium and gold. Elastic scattering can change the EAL, MED and ID by up to about 40%. The calculations in this post ignore elastic scattering and use the IMFP, so real depths are shallower than the numbers given here.

## 3. How Does the Mean Free Path Change with Electron Energy?

Which peak is more surface sensitive? The short answer is that the lower its kinetic energy, the shallower a peak sees. $\lambda$ varies with material, but it varies even more with the electron's kinetic energy. Even within SiO₂, O 1s photoelectrons (kinetic energy about 950 eV) have $\lambda$ = 2.9 nm, while Si 2p photoelectrons (about 1,380 eV) have 3.90 nm. This section first looks at the shape of that dependence, then at where the $\lambda$ values used in practice come from.

### What Shape Is the Universal Curve?

In the overlayer method, a thin film of known thickness is deposited on a substrate and the attenuation length is obtained from how much the substrate signal falls. When attenuation lengths measured this way on many materials in the late 1960s and 1970s are plotted against kinetic energy, most of them fall on a similar V-shaped curve, known as the universal curve. Seah and Dench (1979) condensed the curve for elements into a single empirical formula in units of monolayers. Those measurements assumed exponential attenuation with no elastic scattering, so strictly the compiled values are not IMFPs but what Section 2 called EALs. Today the formula serves only as a qualitative guide.

$$ \lambda_m = \frac{538}{E^2} + 0.41\,(aE)^{1/2} $$

$E$ is the electron energy above the Fermi level in eV, $a$ the monolayer thickness in nm, and $\lambda_m$ the number of monolayers. The first term describes the drop in scattering at low energy, where the electron has too little energy to excite valence electrons. The second describes high energies, where a faster electron spends less time interacting. With $a = 0.25$ nm, converted to nanometres ($a\lambda_m$), the minimum is 0.41 nm at about 41 eV. Electrons of a few tens of eV are the most surface sensitive.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/en/fig1-imfp-universal-curve.png" alt="TPP-2M IMFPs for Si, SiO2 and Au with the Seah–Dench universal curve and representative Al Kα peaks" width="640">
_Figure 1. IMFP versus kinetic energy. Solid lines: the TPP-2M predictive formula (used above 50 eV, less accurate below 200 eV). Dashed line: the Seah–Dench empirical formula (fitted to measured attenuation lengths, so it includes elastic-scattering effects). Circles mark Si 2p (Si), O 1s (SiO₂) and Au 4f (Au) photoelectrons excited by Al K$\alpha$. Near 1.4 keV, $\lambda$ ranges from 1.6 nm (gold) to 3.9 nm (SiO₂), and gold has about half the value of silicon at the same energy._

Figure 1 gives a feel for the numbers. Near 50–100 eV, $\lambda$ is 0.4–0.8 nm, the shortest part of the curve. TPP-2M is less accurate below 200 eV, though, so these values are only rough. Over the range where Al K$\alpha$ photoelectrons mostly sit, from a few hundred eV to 1.4 keV, the IMFP grows with energy and reaches 1.6 nm (gold) to about 4 nm (SiO₂) at 1.4 keV. Reaching tens of nanometres takes kinetic energies above 10 keV, the territory of HAXPES (hard X-ray photoelectron spectroscopy), which uses high-energy X-rays. At high energy the IMFP grows roughly as $E^{\,p}$. The high-energy limit of the Seah–Dench formula gives $p = 0.5$, and fits to measured attenuation lengths for various materials give $p$ between 0.54 and 0.81.

That each peak sees a different depth matters immediately for quantification. Computing a composition from the intensity ratio of two peaks in the same sample, such as O 1s and Si 2p in SiO₂, divides two averages taken over different depth ranges. This problem carries straight into quantification in Part 4 and thickness extraction in Part 5.

### Where Do the Mean Free Path Values Used in Practice Come From?

The values used in practice today are calculated, not empirical. Tanuma, Powell and Penn calculated IMFPs from energy-loss functions derived from optical constants, and built a predictive formula that reproduces the results from four material parameters: the number of valence electrons, density, molecular weight and band gap. This is TPP-2M. The solid lines in Figure 1 were computed directly from it.

```python
def tpp2m(E, Nv, rho, M, Eg=0.0):
    """TPP-2M IMFP (nm). E: kinetic energy (eV), Nv: valence electrons, rho: density (g/cm^3)."""
    Ep = 28.816 * np.sqrt(Nv * rho / M)            # free-electron plasmon energy
    U = (Ep / 28.816) ** 2
    beta = -1.0 + 9.44 / np.sqrt(Ep**2 + Eg**2) + 0.69 * rho**0.1
    gamma = 0.191 * rho**-0.5
    C, D = 19.7 - 9.1 * U, 534 - 208 * U
    return E / (Ep**2 * (beta * np.log(gamma * E) - C / E + D / E**2))
```

The coefficients are those given in Powell (2020). The relativistic correction factor $\alpha(E)$ in Powell's equation is left out, because below 1.5 keV it changes the IMFP by less than 0.4%. The full script is in [`_code/xps-surface-sensitivity-imfp/generate_figures.py`](https://github.com/eddykim/eddykim.github.io/blob/main/_code/xps-surface-sensitivity-imfp/generate_figures.py).

Section 6 of [Part 1](/en/posts/xps-photoemission-binding-energy/) said that models stand between what is measured and what we want to know. Here one of those models shows up in numbers for the first time. The $\lambda$ that goes into a thickness or composition calculation is usually not a measurement but a calculated value from TPP-2M or the NIST database (SRD 71). According to Powell, IMFPs calculated from optical data carry uncertainties of up to about 10%. TPP-2M deviates from those calculations by about 9% RMS, and experimental values from elastic-peak measurements differ from the calculations by 12–15% RMS. A thickness result is directly proportional to $\lambda$, so the uncertainty in $\lambda$ becomes the uncertainty in the thickness.

## 4. Where Does the "Top 10 nm" Information Depth Come From?

One integral gives the answer. Assume normal emission ($\theta = 0$) and integrate the contribution by depth. If photoelectrons are generated uniformly regardless of depth, the signal from depth $z$ is proportional to $e^{-z/\lambda}$. The fraction of the total signal that comes from between the surface and depth $d$ is

$$ \frac{\int_0^{d} e^{-z/\lambda}\,dz}{\int_0^{\infty} e^{-z/\lambda}\,dz} = 1 - e^{-d/\lambda} $$

Integer multiples of $\lambda$ give familiar numbers: 63.2% at $d = \lambda$, 86.5% at $2\lambda$ and 95.0% at $3\lambda$.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/en/fig2-depth-contribution.png" alt="Contribution by depth and cumulative contribution for Si 2p photoelectrons, with λ, 2λ and 3λ marked" width="720">
_Figure 2. Signal contribution by depth for Si 2p photoelectrons in silicon ($\lambda$ = 3.09 nm, TPP-2M). (a) Share of each 1 nm slice (bar height is the mean of the curve over that slice). (b) Cumulative contribution. The first $3\lambda$ = 9.3 nm give 95% of the signal, and the top nanometre alone gives 28%._

With $\lambda$ = 3.09 nm for Si 2p photoelectrons in silicon, $3\lambda$ = 9.3 nm. That is where "XPS sees the top 10 nm" comes from. Unpacking the assumptions behind the number changes the picture.

First, 95% is a convention, not something physics dictates. The ISO definition sets the information depth at a "specified percentage" of the signal and gives both 95% and 99% as examples. With 99% the depth becomes $4.6\lambda$, about 14 nm. Second, $\lambda$ depends on the peak and the material. The same calculation gives $3\lambda$ = 8.8 nm for O 1s in SiO₂ and 4.7 nm for Au 4f in gold. On gold, XPS sees only about half of the "10 nm." Third, as Section 2 showed, including elastic scattering makes the real information depth shallower. Fourth, normal emission was assumed. Section 5 deals with that.

More important is that the signal is an average weighted exponentially in depth. Figure 2(a) shows the share of each 1 nm slice: 28% from the top nanometre, 20% from the next, then 14% and 10%. The first nanometre weighs nearly twice as much as the third. For gold the top nanometre accounts for 47%. The composition XPS reports is not "the average composition of the top 10 nm" but an average weighted heavily toward the outermost layer. For a sample whose composition varies with depth, the difference changes the result.

## 5. Why Does Tilting the Detection Angle Make XPS See Shallower?

What should change if we want to look shallower? Changing the X-ray energy is hard, but tilting the sample is easy. In the equation of Section 2, $\lambda\cos\theta$ took the place of $\lambda$. Tilting the detection direction away from the normal lengthens the path an electron from a given depth must travel, so the electrons that escape without loss come from shallower depths. The mean escape depth and the information depth both shrink in proportion to $\cos\theta$.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/en/fig3-takeoff-angle.png" alt="Escape-path geometry versus detection angle, and cumulative contribution curves by angle" width="720">
_Figure 3. (a) From the same depth $z$, detection at $\theta$ = 60° doubles the escape path. (b) Cumulative contribution by angle for Si 2p in silicon ($\lambda$ = 3.09 nm). Dots mark the 95% depth $3\lambda\cos\theta$, which shrinks with tilt to 2.4 nm at 75°._

Figure 3(b) plots the cumulative contribution at several angles for Si 2p in silicon. The 95% depth is 9.3 nm at $\theta$ = 0°, 6.6 nm at 45°, 4.6 nm at 60° and 2.4 nm at 75°. The share from the top nanometre grows from 28% to 48% at 60° and 71% at 75°. Comparing spectra taken at normal and grazing emission shows which components sit toward the surface. Extending this idea to a series of angles gives angle-resolved XPS (ARXPS), the subject of Part 5. The $\cos\theta$ relation has limits, though. Powell notes that a single EAL is usually adequate up to emission angles of about 60°, but beyond that the EAL itself depends on angle. Near-grazing angles such as 75° are where elastic scattering starts to matter.

### Which Plane Is the Angle Measured From?

This is the place to flag a common practical mistake: the literature measures the detection angle from different reference planes. This post and Powell (2020) use the emission angle $\theta$ from the surface normal, so the effective depth is written $\lambda\cos\theta$. Papers that use the take-off angle measured from the surface plane write the same physics as $\lambda\sin\theta'$. This post writes the take-off angle as $\theta'$. The two are related by $\theta' = 90° - \theta$, and normal emission is $\theta$ = 0° but $\theta'$ = 90°.

A typical example is the Strohmeier equation, widely used to compute oxide thickness.

$$ d = \lambda_{ox}\sin\theta'\,\ln\!\left(\frac{N_m\lambda_m I_{ox}}{N_{ox}\lambda_{ox} I_m} + 1\right) $$

$N$ is the volume density of metal atoms, and $I$ the peak areas of the oxide and metal components. The equation is written in the take-off angle $\theta'$. At normal emission $\sin 90° = 1$ and $\lambda_{ox}$ remains unchanged, so it expresses the same physics as the normal-referenced $\cos\theta$ form. Yet many widely used practical guides that present this equation simply call $\theta'$ the "take-off angle" without saying which plane it is measured from. On an instrument whose analyser is not at 45° from the normal, plugging the normal-referenced angle into the equation gives the wrong thickness. For an analyser tilted 30° from the normal, the correct value is $\sin 60°$ = 0.87; using $\sin 30°$ = 0.5 instead gives a thickness 42% too small. Before borrowing an equation, check the reference plane of its angle.

## 6. When Is Seeing Only the Surface a Disadvantage?

Surface sensitivity is an advantage only if the information depth is filled with what we want to know about.

Nearly every air-exposed sample is covered with adventitious carbon, a layer of hydrocarbons and oxygen-containing organics picked up from the air. Its thickness ranges from under 1 nm to about 2 nm, depending on exposure. That thickness is of the same order as the IMFP. Approximating the hydrocarbon layer as polyethylene and applying TPP-2M gives $\lambda \approx$ 4.5 nm for Si 2p photoelectrons in that layer. A 1 nm adventitious layer leaves 80% of the underlying signal, and a 2 nm layer leaves 64%. Since the signal is weighted toward the surface, as Section 4 showed, the heaviest first one or two nanometres are occupied by a layer we are not interested in. Removing it by ion sputtering can alter the chemical state of what lies underneath, which is a topic for Part 5.

Bulk information is out of reach in principle. In alloys and doped materials, surface segregation, where one element accumulates at the surface, is common. When XPS then gives a composition different from the bulk, that is not a measurement error. XPS by definition measures surface composition, and reading the result as bulk composition is a misinterpretation.

The ultra-high vacuum (UHV) that XPS instruments need also follows from surface sensitivity, for two reasons. One is recontamination. The rate at which gas molecules strike a surface is proportional to pressure. At $1.3\times10^{-6}$ mbar ($10^{-6}$ Torr), enough molecules strike in about one second to cover a monolayer if every one sticks (this exposure is one langmuir). At an analysis-chamber pressure around $10^{-10}$ mbar the same exposure takes about 3.7 hours, so the surface survives the measurement. The other is gas-phase scattering of the photoelectrons. An electron that loses energy colliding with a gas molecule on its long path to the analyser drops out of the peak just as if it had lost that energy inside the solid. Figure 4 shows an XPS instrument built to meet these conditions.

<img src="/assets/img/posts/xps-surface-sensitivity-imfp/ext-xps-spectrometer.jpg" alt="Photograph of a teaching XPS instrument with a stainless-steel UHV chamber and control rack" width="520">
_Figure 4. An XPS instrument used for student experiments by the Institute for Experimental Physics of the University of Hamburg at the DESY research centre. All the analysis hardware is attached to one stainless-steel UHV chamber, with power supplies and controllers in the rack on the right. The bolted circular flanges and viewports crowding the chamber are there because, for the two reasons above, the whole chamber must be held at ultra-high vacuum. (Source: MichaelBueker, [Desyxps](https://commons.wikimedia.org/wiki/File:Desyxps.jpg), [CC BY 3.0](https://creativecommons.org/licenses/by/3.0))_

## 7. How Does This Compare with Optical Metrology?

What would an optical measurement of the same film see? In [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/) the penetration depth of light was written $\delta = \lambda/4\pi k$ (here $\lambda$ is the wavelength of light). At weakly absorbing wavelengths where $k$ is small, $\delta$ reaches hundreds of nanometres to several micrometres. Ellipsometry sees through an entire multilayer stack at once and extracts each layer's thickness and optical constants through a model. XPS sees only the top few nanometres of that stack.

What each technique does well also differs. XPS distinguishes elements and chemical states directly, but thickness is the most indirect quantity it yields. Ellipsometry gives thickness and optical constants with sub-nanometre sensitivity, but which chemical bonds make up a layer must be put into the optical-constant model in advance. The two techniques complement rather than compete with each other. For a thin oxide film, a common approach is to cross-check XPS for chemical state and approximate thickness against ellipsometry for precise thickness. Part 5 returns to this comparison when it examines the limits of XPS thickness measurement.

## Summary and Next Part

This post traced where the surface sensitivity of XPS comes from. X-rays penetrate micrometres into a solid, but only photoelectrons that escape without inelastic scattering contribute to a peak, and that probability falls as $\exp(-z/\lambda\cos\theta)$. The IMFP $\lambda$ follows the universal curve with kinetic energy and is 1.6–4 nm for shallow core levels excited by Al K$\alpha$ (kinetic energies near 1.4 keV). The $\lambda$ used in practice is usually a calculated value from a predictive formula such as TPP-2M, with an uncertainty of around 10%. The "10 nm" information depth is the product of several assumptions: the 95% convention, a silicon-like $\lambda$, normal emission, and neglect of elastic scattering. The signal is an exponentially weighted average in depth that leans toward the outermost layer, and tilting the detection angle strengthens that lean. When quoting an angle, check its reference plane first.

The next part returns to the spectrum itself. In the survey spectrum of Part 1, each element seemed to have a single peak, but at high resolution the picture changes. There are spin-orbit doublets, shake-up satellites, multiplet splitting, plasmon losses and Auger series. Telling whether each comes from initial-state or final-state effects is what makes it possible to read chemical states correctly.

## References

- C. J. Powell, "Practical guide for inelastic mean free paths, effective attenuation lengths, mean escape depths, and information depths in x-ray photoelectron spectroscopy," *J. Vac. Sci. Technol. A* 38, 023209 (2020). <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=929257>
- S. Tanuma, C. J. Powell, and D. R. Penn, "Calculations of electron inelastic mean free paths. V. Data for 14 organic compounds over the 50–2000 eV range," *Surf. Interface Anal.* 21, 165–176 (1994). <https://doi.org/10.1002/sia.740210302>
- M. P. Seah and W. A. Dench, "Quantitative electron spectroscopy of surfaces: A standard data base for electron inelastic mean free paths in solids," *Surf. Interface Anal.* 1, 2–11 (1979). <https://doi.org/10.1002/sia.740010103>
- B. R. Strohmeier, "An ESCA method for determining the oxide thickness on aluminum alloys," *Surf. Interface Anal.* 15, 51–56 (1990). <https://doi.org/10.1002/sia.740150109>
- NIST Standard Reference Database 71: Electron Inelastic-Mean-Free-Path Database. <https://www.nist.gov/srd/nist-standard-reference-database-71>
- NIST Standard Reference Database 82: Electron Effective-Attenuation-Length Database. <https://www.nist.gov/srd/nist-standard-reference-database-82>
- Center for X-Ray Optics (CXRO), X-Ray Attenuation Length. <https://henke.lbl.gov/optical_constants/atten2.html>
- Thermo Fisher Scientific, XPS Periodic Table — Carbon. <https://www.thermofisher.com/us/en/home/materials-science/learning-center/periodic-table/non-metal/carbon.html>
- D. N. G. Krishna and J. Philip, "Review on surface-characterization applications of X-ray photoelectron spectroscopy (XPS): Recent developments and challenges," *Appl. Surf. Sci. Adv.* 12, 100332 (2022). <https://doi.org/10.1016/j.apsadv.2022.100332>
- D. J. Morgan, "X-Ray Photoelectron Spectroscopy (XPS): An Introduction," Cardiff Catalysis Institute. <https://sites.cardiff.ac.uk/xpsaccess/files/2014/07/AccessXPS_Primer_Paper.pdf>

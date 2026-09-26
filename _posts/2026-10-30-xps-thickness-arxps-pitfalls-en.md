---
title: "XPS Fundamentals 5 — Thickness Measurement and Its Pitfalls"
lang: en
lang-exclusive: ["en"]
permalink: /posts/xps-thickness-arxps-pitfalls/
page_id: xps-thickness-arxps-pitfalls
date: 2026-10-30 20:00:00 +0900
categories: [Surface Analysis, Depth Profiling]
tags: [xps, arxps, film-thickness, depth-profiling, beam-damage, metrology]
description: "Attenuation itself can measure film thickness. The overlayer equation, the ARXPS inverse problem, and the pitfalls of preferential sputtering and beam damage."
math: true
---

Up to [Part 4](/en/posts/xps-quantification-peak-fitting/), the sample was assumed uniform in depth within the information depth. The $\lambda$ in the composition equation carried that assumption. Real samples, however, are nearly always layered: the native oxide on a silicon wafer, the adventitious carbon covering every surface, a deliberately deposited film. The synthetic spectrum of Part 4 was silicon covered with a thin SiO₂ layer.

Layering is not only a nuisance. The exponential attenuation of [Part 2](/en/posts/xps-surface-sensitivity-imfp/) carries depth information. The thicker the oxide, the more the signal from the substrate beneath it is reduced. Read in reverse, how much the signal was reduced tells you the thickness of the oxide. This post covers how thickness is measured on this principle, and where the method goes wrong.

The short version: thickness is the most indirect quantity XPS yields. It comes from a single intensity ratio through a model, and depends entirely on the assumption of a uniform layer and on an estimated $\lambda$. The second half of the post goes a step further. Sputtering a sample or exposing it to X-rays for a long time in order to see depth can change the very thing being measured.

## 1. The Uniform Overlayer Model

If a metal signal survives beneath an oxide, does how much of it survives tell us the thickness? Start with the simplest geometry. A uniform overlayer of thickness $d$, say SiO₂, lies on an infinitely thick substrate (Si), and the analyzer sits at an angle $\theta$ from the surface normal. As in Part 2, $\theta$ in this post is always **measured from the normal**.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/en/fig1-overlayer-geometry.png" alt="Geometry of the uniform overlayer model and intensity ratio versus thickness" width="720">
_Figure 1. (a) The uniform overlayer model. The substrate signal is attenuated by the oxide. (b) Intensity ratio versus thickness (synthetic calculation, $\lambda$ = 3 nm). The dotted line marks $3\lambda$, where the substrate signal drops to 5%._

A photoelectron from depth $z$ in the overlayer escapes without loss with probability $\exp(-z/\lambda\cos\theta)$. Integrating from the surface to $d$ gives the overlayer signal:

$$ I_o = I_o^{\infty} \left[ 1 - \exp\!\left( -\frac{d}{\lambda_o \cos\theta} \right) \right] $$

$I_o^{\infty}$ is the intensity from an infinitely thick overlayer. The substrate signal must pass through the whole overlayer of thickness $d$ and is attenuated once more:

$$ I_m = I_m^{\infty} \exp\!\left( -\frac{d}{\lambda_m \cos\theta} \right) $$

Dividing the two gives the intensity ratio $R = I_o/I_m$. When the two peaks have nearly the same kinetic energy, as for the Si 2p of SiO₂/Si, one can set $\lambda_o \approx \lambda_m = \lambda$ and the expression simplifies neatly:

$$ R = R_0 \left[ \exp\!\left( \frac{d}{\lambda\cos\theta} \right) - 1 \right], \qquad R_0 = \frac{I_o^{\infty}}{I_m^{\infty}} $$

Solving for $d$ gives the thickness from a single intensity ratio:

$$ d = \lambda \cos\theta \, \ln\!\left( 1 + \frac{R}{R_0} \right) $$

Figure 1(b) plots this relation for $\lambda$ = 3 nm. That value is the IMFP of Si 2p in SiO₂ computed with TPP-2M in Part 2 (3.88 nm) multiplied by the lower end of the EAL/IMFP ratio from Powell cited there (0.77). As thickness grows, the ratio grows exponentially, and at $3\lambda$ = 9 nm the substrate signal falls to 5% of its original value. Beyond that the substrate peak is lost in the background and thickness can no longer be measured. Tilting the detection angle to 60° makes the ratio much larger at the same thickness, increasing sensitivity to thin films but halving the measurable range.

The expression rests on five assumptions:

1. The overlayer is uniform and flat.
2. The interface between overlayer and substrate is abrupt.
3. $\lambda$ and $R_0$ are known.
4. Elastic scattering is neglected (or an effective attenuation length is used in place of $\lambda$).
5. There is no surface contamination layer, or any layer affects both signals equally.

The following sections trace where each of these assumptions breaks.

## 2. The Strohmeier Equation and the $R_0$ Trap

What form of the equation is actually used in practice? The most widely used for oxide thickness is Strohmeier's (1990):

$$ d = \lambda_o \sin\theta' \, \ln\!\left( \frac{N_m \lambda_m I_o}{N_o \lambda_o I_m} + 1 \right) $$

$N_m$ and $N_o$ are the volume densities of metal atoms in the metal and the oxide. Here $\theta'$ is the take-off angle **measured from the surface plane**. This is exactly where the angle confusion warned about in Section 5 of Part 2 arises. Compared with the expression in Section 1, $\sin\theta' = \cos\theta$, that is $\theta' = 90° - \theta$, and the $N_m\lambda_m/(N_o\lambda_o)$ inside the logarithm is a theoretical estimate of $1/R_0$. The two expressions describe the same physics. But on an instrument whose analyzer sits 30° from the normal, entering 30° directly as $\theta'$ in this equation gives, for an intensity ratio of 1.00, a thickness of 1.14 nm instead of 1.97 nm. That is why the reference plane of the angle must be checked before borrowing an equation.

A deeper trap lies in $R_0$. The Strohmeier equation computes $R_0$ from ratios of atomic densities and attenuation lengths. Seah and Spencer (2002) measured this value directly for SiO₂/Si by two routes and obtained 0.88 ± 0.03, far from the calculated 0.53 ± 0.05. Possible causes include systematic uncertainties in the IMFPs and in the density of SiO₂, and the share of intensity that shake-up satellites (Part 3) carry outside the main peak. The effect on thickness is large. For the same intensity ratio of 1.00, the measured $R_0$ gives 2.28 nm and the calculated one 3.18 nm.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/en/fig2-thickness-sensitivity.png" alt="Statistical versus systematic errors in thickness inversion: ratio error and counting noise, and the effects of λ and R0" width="720">
_Figure 2. Precision and accuracy of thickness inversion (synthetic calculation, $\lambda$ = 3 nm, $R_0$ = 0.88). (a) Statistical error. (b) Systematic error: $\lambda$ off by ±10% (green band) and using the calculated $R_0$ = 0.53 instead of the measured value (red line)._

Figure 2 sets the two kinds of error side by side. The code is in [`_code/xps-thickness-arxps-pitfalls/`](https://github.com/eddykim/eddykim.github.io/tree/main/_code/xps-thickness-arxps-pitfalls), and running `overlayer.py` prints every number in this section and in Section 4. The statistical error on the left is small. For a measurement collecting 10⁵ counts from the bare substrate, the standard deviation of thickness from counting noise is 0.004 nm for a 0.5 nm film and still only 0.049 nm for a 10 nm film. Even deliberately adding a 5% error to the intensity ratio gives a thickness error that starts at 0.023 nm and saturates near 0.15 nm. As the thickness grows, the ratio grows exponentially, so the same fractional error becomes an ever smaller share of the thickness.

The systematic error on the right is another story. Taking $\lambda$ 10% too large makes the thickness exactly 10% too large, since $d$ is proportional to $\lambda$. Using the calculated $R_0$ turns a 1 nm film into 1.51 nm (+51%), a 4 nm film into 5.19 nm (+30%) and an 8 nm film into 9.44 nm (+18%). The thinner the film, the larger the relative error.

Seah (2005) set out this contrast with measurements on SiO₂/Si. The precision of XPS thickness can be as good as 0.025 nm at one standard deviation, but in the past it was overshadowed by systematic errors of about 20% of the thickness from uncertainty in the attenuation length. With a newly calibrated attenuation length this can be reduced to about 2% at 95% confidence, but variability in procedures across ordinary laboratories often degrades it to 0.4 nm in practice. Seah and Spencer (2002) also reported that, because the crystal structure of the substrate makes intensities vary with direction, data taken at normal emission give thicknesses 18% too small. These numbers show how heavy assumptions 1 and 3 really are.

Assumption 5, the contamination layer, conversely matters very little for this method. An adventitious carbon layer attenuates both signals, but the oxide and substrate components of Si 2p differ in kinetic energy by only 4 eV. Their attenuation lengths in the hydrocarbon layer, 4.716 nm and 4.726 nm, are practically identical, so a 1 nm contamination layer changes an inferred thickness of 2 nm by less than 1 pm. That is the strength of comparing two chemical states of the same element, and this property becomes important again in Section 7.

## 3. The Thickogram and Layers of Different Elements

What if the overlayer and substrate are different elements altogether, such as a metal film on silicon or a polymer coating on a metal? Then the two peaks differ in kinetic energy, $\lambda_o \ne \lambda_m$, and the neat inversion of Section 1 no longer applies. Finding the thickness requires solving the equations numerically.

Cumpson's (2000) Thickogram turns this calculation into a single chart. The intensity ratio of the two peaks, the ratio of their sensitivity factors and the ratio of their kinetic energies are placed on the chart and the thickness read off. Because it uses intensity ratios, instrument factors cancel and the effect of contamination is reduced. Its accuracy, however, is ultimately set by the accuracy of the attenuation lengths entered, since thickness is proportional to them, as Section 2 showed.

Non-planar shapes make things harder still. In core–shell nanoparticles, where a spherical shell wraps a core, the path length of photoelectrons varies from point to point. Shard (2012) proposed a straightforward method for extracting shell thickness from XPS data of such core–shell particles. This section is kept short for one reason: to show how strong an assumption the "uniform planar layer" of Section 1 is. Change the shape and the same intensity ratio means an entirely different thickness.

## 4. ARXPS — From an Angle Series to a Depth Profile, and the Inverse Problem

Can measuring at many angles recover the layer structure in full? The expression of Section 1 assumed a single layer with an abrupt interface. When the layer structure is unknown, the principle of Section 5 of Part 2 is extended by measuring repeatedly at different detection angles. This is angle-resolved XPS (ARXPS). With $c(z)$ the oxide fraction at depth $z$, the intensity measured at angle $\theta$ is

$$ I(\theta) \propto \int_0^{\infty} c(z) \exp\!\left( -\frac{z}{\lambda\cos\theta} \right) dz $$

This is a Laplace transform of $c(z)$. Changing the angle changes the transform variable $1/(\lambda\cos\theta)$. Recovering $c(z)$ from $I(\theta)$ at several angles is therefore a problem of inverting a Laplace transform.

Integral equations with exponential kernels like this are classic ill-posed inverse problems. The exponential integrates over and smears out the fine structure of $c(z)$, so quite different distributions give nearly the same $I(\theta)$. Solving in reverse amplifies small noise in the data into large changes in the solution. Cumpson has long studied the depth-resolution limits of ARXPS; in a recent preprint (2026) he calls free-form depth-profile reconstruction "a severely ill-posed inverse Laplace-transform problem" and summarizes earlier analysis as showing that the meaningful depth resolution is governed mainly by the signal-to-noise ratio and is a substantial fraction of the depth itself.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/en/fig3-arxps.png" alt="Depth profiles of an abrupt interface and graded interfaces of different widths, with the difference in intensity ratio by angle and an uncertainty band" width="720">
_Figure 3. Can ARXPS distinguish interface structure? (synthetic calculation). (a) An abrupt interface at 2 nm and graded interfaces of different widths. (b) Difference in intensity ratio by angle. The gray band combines a ±0.5° error in detection angle and counting noise (10⁵ substrate counts)._

Figure 3 puts numbers on this ill-posedness. It compares an abrupt interface at 2 nm with graded interfaces spread out in an error-function shape. For each graded interface the center position was adjusted, along with a constant multiplying the whole intensity ratio within the measured uncertainty of $R_0$ (±0.03/0.88, about ±3.4%), to match the abrupt interface as closely as possible. The uncertainty band combines the 0.5° angle-setting error that Seah (2005) identified as a limiting factor with counting noise: ±0.65% at 0°, ±1.29% at 40° and ±5.80% at 70°.

The result: graded interfaces 0.3 nm or 0.5 nm wide stay within the band at every angle from 0° to 70°. They cannot be distinguished from the abrupt interface. Only at a width of 0.8 nm does the difference leave the band, at 21 of 29 angles, and at 1.2 nm at 25. For an interface 2 nm deep, structure of 0.5 nm or less leaves practically no trace in the angle data. This conclusion uses the full $R_0$ uncertainty. Knowing $R_0$ exactly would lower the limit somewhat, but as Section 2 showed, that value itself varied widely between sources.

So in practice one does not try to recover an arbitrary $c(z)$ freely. Instead a layer model with an assumed number and order of layers is fitted, or regularization is imposed to keep the solution smooth. It is the same idea as the way the Levenberg–Marquardt method of [Optimization 3](/en/posts/optimization-levenberg-marquardt/) adds diagonal terms to $J^T J$ to tame a poor condition number. Fitting a layer model is also prone to converging to different local minima depending on the starting values, as seen in [Optimization 4](/en/posts/optimization-global-heuristics/). When you see a result claiming to have freely recovered a complex depth profile from angle data alone, check first what assumptions and regularization went into it.

There is one more nondestructive way to see depth: the Tougaard background analysis anticipated in Section 7 of [Part 3](/en/posts/xps-spectrum-structure-satellites/) and Section 4 of Part 4. The probability of a photoelectron losing energy depends on the distance it travels in the solid, so the shape of the background under a peak encodes how the element is distributed in depth. An element only at the surface loses little energy and has a low background behind its peak; an element buried deeper has a large loss tail. Tougaard developed methods that compute the background shapes of different depth distributions from inelastic-scattering cross-sections and compare them with the measured spectrum, estimating nanoscale depth structure from a single measurement without changing the angle. This method, too, rests on models, the background and the cross-section, and shares the limits of the same exponential kernel as the angle method.

The common advantage of ARXPS and background analysis is that they do not remove material. Their common limit is that they cannot see beyond the information depth. Structures deeper than about 10 nm are out of reach of either.

## 5. Sputter Depth Profiles and Preferential Sputtering

To see deeper, why not just remove material? The most direct approach is to erode the surface step by step with an argon ion beam and repeat the XPS measurement. Converting sputter time to depth gives a composition-versus-depth curve. It can reach hundreds of nanometers beyond the information-depth limit, and it can remove the adventitious carbon layer that Part 2 wished away.

The trouble is that the ion beam changes the composition. When constituent elements are not removed in the same proportion, this is called preferential sputtering. In oxides, oxygen tends to leave first. Counsell et al. (2014) reported that eroding amorphous TiO₂ with 5 kV monatomic Ar⁺ ions preferentially removes oxygen and reduces Ti from the +4 state to +3 and +2. Ti³⁺ and Ti²⁺ that were not in the original sample appear in the spectrum. Reading this as "there is a reduced oxide at the interface" reports an artifact of the measurement as a property of the sample.

There are other artifacts too: atomic mixing, where incoming ions drive atoms deeper and blend the layers; a surface that roughens as erosion proceeds; and sputter rates that differ from layer to layer, so the conversion from sputter time to depth differs between layers. On top of this comes the information depth of XPS itself. Each step measures not just the newly exposed surface but an exponentially weighted average over a few nanometers beneath it.

<img src="/assets/img/posts/xps-thickness-arxps-pitfalls/en/fig4-sputter-artifact.png" alt="True oxygen profile and the observed profile after atomic mixing, information-depth weighting and preferential oxygen removal" width="640">
_Figure 4. How a sputter depth profile smears an interface (schematic). A 5 nm MO₂ oxide over metal M. The atomic mixing width of 1 nm, information depth of 2 nm and 30% preferential oxygen removal are arbitrary values chosen for illustration._

Figure 4 is a schematic that adds these effects one at a time. The true profile ends abruptly at 5 nm, but atomic mixing spreads the interface, the exponential weighting of the information depth smears it again, and preferential oxygen removal lowers the oxygen fraction of the whole oxide. From the observed profile alone, the oxide looks thinner than it is, its interface wider, and its composition oxygen-deficient.

There are mitigations. Lower the ion energy, and rotate the sample during erosion to reduce roughening. For organics and polymers, and when the chemical state of an oxide must be preserved, use cluster ions (Ar$_n^+$) of hundreds to thousands of argon atoms. Each atom carries little energy, so the damage does not penetrate deeply. Counsell et al. likewise reported that cluster ions greatly reduce damage to TiO₂ and ion incorporation. Still, the key message stands. A sputter depth profile is useful for composition trends, but chemical states read after sputtering must not be reported as the chemical states of the original sample.

## 6. The Measurement Changes the Sample — Beam Damage

Aren't X-rays nondestructive? Introductions to XPS often call it a nondestructive technique. Morgan (2023) addresses exactly this belief: many think of XPS as nondestructive, but the X-ray photons and the electron cascade that follows can change the analyzed area significantly. That holds even for monochromated X-rays.

Easily reduced species are especially vulnerable. The Cardiff University XPS reference pages note that Au(III) compounds are photoreduced during photoemission, with both Au(III) and Au(0) present after 10 minutes of irradiation. Biesinger et al. (2010) reported that V₂O₅ slowly turns into V(IV) under X-rays, with more than 15% converted over 24 hours with a 210 W source. Cu(II), identified in Part 3 by its satellites, is another species for which X-ray reduction must be considered. Bond scission is common in polymers.

The practical procedure is simple. Record sensitive regions twice, at the start and at the end of the measurement. The Cardiff pages recommend the same for the Au 4f region: acquire it first on its own, then again after all other regions, to assess the degree of reduction. If the two spectra differ, do not draw conclusions about chemical state from those data. Shortening the exposure or moving to fresh spots are further options.

This section goes a step beyond the thread of this series. Parts 1 through 4 looked at where the models between what is measured and what we want to know break down. Here, before any model, the object being measured changes during the measurement. It is the most fundamental form of model failure. No model will recover the original sample from one that changed while it was being measured.

## 7. What Optical Metrology and XPS Supply Each Other

The same SiO₂/Si thickness can be measured by ellipsometry or by XPS; what is the difference? Side by side:

| Item | Ellipsometry | XPS |
|---|---|---|
| What is measured | Change in polarization of reflected light ($\Psi, \Delta$) | Kinetic-energy distribution of photoelectrons |
| Depth that can be seen | Hundreds of nm to several µm when absorption is weak | About 10 nm |
| Chemical-state information | Must be put into the optical-constant model in advance | Obtained directly |
| Measurement environment | Air | Ultra-high vacuum |
| Speed | Seconds, possible in-line during processing | Minutes or more |
| Prior knowledge needed for thickness | Layer structure and optical constants of each layer | Layer structure, attenuation length, $R_0$ |

The last row is the key. Ellipsometry needs optical constants to give a thickness; XPS needs $\lambda$ and the layer structure. Both are model inversions; they differ only in the kind of prior knowledge they require.

This difference has also shown up in measurements. In a pilot study under the Consultative Committee for Amount of Substance (CCQM), Seah et al. (2004) compared 45 results from 10 methods measuring thermal oxides 1.5–8 nm thick. All methods showed excellent linearity against the reference values, but each had an inherent zero-thickness offset between 0 and 1 nm. Ellipsometry saw surface water and carbon contamination as about 1 nm of oxide, and ion-beam methods such as MEIS, NRA and RBS saw adsorbed oxygen as about 0.5 nm. Only XPS using the Si 2p peaks had zero offset, because, as Section 2 showed, the contamination layer attenuates both Si 2p components equally. XPS, on the other hand, had a large uncertainty in its scaling constant, the attenuation length. The study therefore combined the XPS zero point with the scales of the other methods and recommended multiplying the reference attenuation lengths by 0.986 ± 0.009. The two kinds of metrology filled in each other's weaknesses.

This concern is the same one first set out in [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/) on this blog. That post started from the fact that detectors measure only intensity, so the phase of the electric field is lost, and ellipsometry recovered that information by comparing polarization states. This series started from the fact that XPS measures only the kinetic energy of photoelectrons, and element, chemical state, composition and thickness were all obtained by building models on top of that. The measurement principles differ, but the structure is the same: a model sits between what is measured and what we want to know, and the assumptions of that model decide the result.

## Summary — Closing the Series

This post followed the use of attenuation in reverse to measure thickness. The uniform overlayer model gives thickness from a single intensity ratio. Its statistical precision is good, around 0.01 nm, but its accuracy is set by the model constants $\lambda$ and $R_0$, and using the calculated $R_0$ puts thin-film thicknesses off by more than 50%. Confusing angle conventions gives entirely different thicknesses from the same data. ARXPS is an inverse problem with an exponential kernel, so at a depth of 2 nm, interface structure of 0.5 nm or less leaves no trace in the angle data. Sputter depth profiles see deeper but create chemical states that were never there through preferential sputtering, and X-rays themselves change easily reduced species.

Looking back over the five parts, the same question kept returning in different forms. In [Part 1](/en/posts/xps-photoemission-binding-energy/), converting kinetic energy to binding energy required the energy-conservation equation and the assumption of Fermi-level alignment, and the measured binding energy differed from the orbital energy because of relaxation. In [Part 2](/en/posts/xps-surface-sensitivity-imfp/), "the top 10 nm" was the product of a 95% convention and a $\lambda$ computed from a predictive formula. In [Part 3](/en/posts/xps-spectrum-structure-satellites/), deciding which structures in a spectrum count as components was already interpretation, and miscounting satellites and multiplets created species that were not there. In [Part 4](/en/posts/xps-quantification-peak-fitting/), a single composition number passed through four models, charge referencing, background, line shape and sensitivity factors, and with the wrong background model the error bars said nothing about the actual error. And in this post, thickness was the most indirect quantity of all, adding the layer structure and $\lambda$ on top of all those models.

From beginning to end, what XPS actually measures is a single distribution of photoelectron kinetic energies. Element, chemical state, composition and thickness are all values computed backward by fitting models to that distribution. As long as this is remembered and each model's assumptions are written down alongside the result, XPS tells us more about the state of a surface than any other method. The moment the gap between what is measured and what we want to know is forgotten, only plausible numbers remain.

## References

- M. P. Seah and S. J. Spencer, "Ultrathin SiO₂ on Si II. Issues in quantification of the oxide thickness," *Surf. Interface Anal.* 33, 640–652 (2002). <https://doi.org/10.1002/sia.1433>
- M. P. Seah et al., "Critical review of the current status of thickness measurements for ultrathin SiO₂ on Si Part V: Results of a CCQM pilot study," *Surf. Interface Anal.* 36, 1269–1303 (2004). <https://doi.org/10.1002/sia.1909>
- M. P. Seah, "Ultrathin SiO₂ on Si. VI. Evaluation of uncertainties in thickness measurement using XPS," *Surf. Interface Anal.* 37, 300–309 (2005). <https://doi.org/10.1002/sia.2020>
- B. R. Strohmeier, "An ESCA method for determining the oxide thickness on aluminum alloys," *Surf. Interface Anal.* 15, 51–56 (1990). <https://doi.org/10.1002/sia.740150109>
- HarwellXPS, "Technical Note #1: Oxide Thickness Determination Using the Strohmeier Equation." <https://subsite.harwellxps.uk/wp-content/uploads/2018/02/HarwellXPS_TechNote_01.pdf>
- P. J. Cumpson, "The Thickogram: a method for easy film thickness measurement in XPS," *Surf. Interface Anal.* 29, 403–406 (2000). <https://doi.org/10.1002/1096-9918(200006)29:6%3C403::AID-SIA884%3E3.0.CO;2-8>
- A. G. Shard, "A straightforward method for interpreting XPS data from core–shell nanoparticles," *J. Phys. Chem. C* 116, 16806–16813 (2012). <https://doi.org/10.1021/jp305267d>
- P. J. Cumpson, "Angle-resolved XPS and AES: Depth-resolution limits and a general comparison of properties of depth-profile reconstruction methods," *J. Electron Spectrosc. Relat. Phenom.* 73, 25–52 (1995). <https://doi.org/10.1016/0368-2048(94)02270-4>
- P. J. Cumpson, "Resolution Limits in Constrained Angle-Resolved XPS Depth-Profile Reconstruction," ChemRxiv (2026). <https://doi.org/10.26434/chemrxiv.15007123/v1>
- S. Tougaard, "Accuracy of the non-destructive surface nanostructure quantification technique based on analysis of the XPS or AES peak shape," *Surf. Interface Anal.* 26, 249–269 (1998). <https://doi.org/10.1002/(SICI)1096-9918(199804)26:4%3C249::AID-SIA368%3E3.0.CO;2-A>
- J. D. P. Counsell, A. J. Roberts, W. Boxford, C. Moffitt, and K. Takahashi, "Reduced Preferential Sputtering of TiO₂ using Massive Argon Clusters," *J. Surf. Anal.* 20, 211–215 (2014). <https://www.jstage.jst.go.jp/article/jsa/20/3/20_211/_article>
- D. J. Morgan, "XPS insights: Sample degradation in X-ray photoelectron spectroscopy," *Surf. Interface Anal.* 55, 331–335 (2023). <https://doi.org/10.1002/sia.7205>
- Cardiff University XPS Access, "Gold." <https://sites.cardiff.ac.uk/xpsaccess/reference/gold/>
- M. C. Biesinger, L. W. M. Lau, A. R. Gerson, and R. St. C. Smart, "Resolving surface chemical states in XPS analysis of first row transition metals, oxides and hydroxides: Sc, Ti, V, Cu and Zn," *Appl. Surf. Sci.* 257, 887–898 (2010). <https://doi.org/10.1016/j.apsusc.2010.07.086>
- C. J. Powell, "Practical guide for inelastic mean free paths, effective attenuation lengths, mean escape depths, and information depths in x-ray photoelectron spectroscopy," *J. Vac. Sci. Technol. A* 38, 023209 (2020). <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=929257>

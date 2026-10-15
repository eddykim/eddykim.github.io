---
date: 2026-11-07 20:00:00 +0900
layout: post
title: "Ellipsometer Instrumentation 5 — Can Polarization Be Measured in One Shot Without Rotating Anything?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-channeled-spectroscopic/
page_id: ellipsometer-channeled-spectroscopic
categories: [Optics, Instrumentation]
tags: [ellipsometry, channeled-spectroscopy, snapshot, fourier-analysis, instrumentation]
description: A thick birefringent plate becomes a variable retarder once a spectrometer is placed after it. This post follows what a single shot costs.
math: true
---

Every instrument from [Post 1](/en/posts/ellipsometer-oblique-rotating-element/) through [Post 4](/en/posts/ellipsometer-calibration/) rotated or oscillated something: an analyzer, a compensator, a photoelastic modulator driven by a voltage. The axis carrying the modulation was always **time**, and the polarization was recovered by Fourier-analyzing the detected signal along that axis.

The instrument in this post has no moving parts and no active elements. It still modulates.

The starting point is a sentence Oka and Kato put at the end of their 1999 paper.

> The phase retardation of a thick birefringent plate changes appreciably with wave number. This implies that the birefringent plate can serve as a **variable retarder when it is combined with a spectroscopic device**.

Instead of rotating the retarder, make it thick. The retardance then differs from wavelength to wavelength, and the moment a spectrometer spreads the spectrum out, **a full period of the modulation appears within a single frame.** The wavelength dispersion of the spectrometer does the work the time axis used to do. A spectrum carrying such fringes is called a channeled spectrum, and measuring polarization with it is channeled spectroscopic ellipsometry. A photograph of a grating spectrometer spreading light by wavelength appears in [Post 2](/en/posts/ellipsometer-components/).

If Post 1 asked what cannot be measured, Post 2 where the components stop working, Post 3 why one modulating element is not enough, and Post 4 how to work with limitations once they are known, this post asks **what a single shot costs.**

## 1. What does a single thick retarder write into the spectrum?

Take the simplest arrangement. A thick retarder sits behind the polarizer, the light reflects off the sample, and an analyzer receives it. With the polarizer at $0°$, the retarder at $45°$ and the analyzer at $-45°$, the detected spectrum reads

$$I(\sigma) = \frac{I_{in}(\sigma)}{4}\Big[1 - a\,\cos\phi(\sigma) - b\,\sin\phi(\sigma)\Big],
\qquad a = \cos 2\Psi,\quad b = \sin 2\Psi\,\sin\Delta$$

What this arrangement actually measures is the pair $a$ and $b$; $\Psi$ and $\Delta$ are computed from them. Keeping the two apart pays off later. The sign of the analyzer azimuth matters too. With the analyzer at $+45°$ the $\sin\Delta$ term flips sign, so under this series' rotation convention $-45°$ is correct. Hagen writes the same expression with $R(45°)$ because he takes the opposite handedness of rotation.

Here $\sigma = 1/\lambda$ is the wave number. Wave number replaces wavelength because the retardance grows almost linearly with it.

$$\phi(\sigma) = 2\pi\,B_{\rm eff}\,t\,\sigma = 2\pi L\sigma$$

Strictly, a constant phase and higher-order terms are added. The constant phase is removed by calibration.

This is the $\delta = 2\pi d\,\Delta n/\lambda$ of Post 2 rewritten in wave number. The reading, however, is inverted. In Post 2 the wavelength dependence of retardance was a **defect**, which is why zero-order plates were preferred and multi-order plates avoided. Here that dependence is the **principle**. The plate is made deliberately thick to maximize it.

The slope $L = B_{\rm eff}\,t$ is the optical path difference (OPD), and this single number fixes where the channels sit. The coefficient $B_{\rm eff}$ is the first-order effective birefringence,

$$B_{\rm eff} = B(\bar\sigma) + \bar\sigma\,\frac{dB}{d\sigma}\Big\vert_{\bar\sigma}$$

This matters because it is not the birefringence $B = n_e - n_o$ itself but that value plus a dispersion term. For quartz Hagen uses about $0.00998$ over 400–1045 nm, rather than the nominal $0.009$ of Post 2. Hagen notes that for common retarder materials the value runs **about 10% above** the nominal birefringence, and indeed $0.009 \times 1.1 = 0.0099$.

Narrowing the band to 400–800 nm gives about $0.0102$ from the dispersion formula of quartz, and back-calculating from Hagen's own measured channel positions gives $0.0101$. The differences stay near 2% and do not change the scale of the thickness limits, so every calculation here uses $0.00998$. Hagen writes this coefficient $\beta$; since Posts 1 and 4 already use that symbol for the Fourier coefficients of a rotating-element instrument, it appears here as $B_{\rm eff}$.

What happens when this spectrum is Fourier transformed along the wave-number axis?

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig1-single-retarder-channels.png" alt="Spectrum and Fourier channels for a single retarder" width="780">
_Figure 1. One sample ($\Psi = 35°$, $\Delta = 50°$) seen through quartz retarders of two thicknesses: 0.4 mm ($L = 4.0\ \mu m$) on top and 3.0 mm ($L = 29.9\ \mu m$) below. Left, the detected spectrum; right, its Fourier transform along the wave-number axis. With the thin plate the $\pm L$ peaks sit against the dc peak and would mix with it under a filter window; with the thick plate the three peaks separate._

Because $\cos\phi$ and $\sin\phi$ each move to $\pm L$, peaks appear at $h = 0, \pm L$. The dc peak carries the source intensity and the sample reflectance; the $\pm L$ peaks carry $a$ and $b$. Only the latter two are wanted, so only those channels need to be cut out.

The trouble is in the upper right panel. A thin retarder gives a small $L$, and the three peaks overlap. Overlapping peaks cannot be separated. The 0.4 mm plate in the upper row is slightly thinner than the 0.46 mm lower bound derived in Section 2.

## 2. How thick must the retarder be for the channels to separate?

How large must $L$ be for the peaks to stay apart? The answer follows from the fact that the peaks have width.

Treat the source intensity as a Gaussian along the wave-number axis and take its span as six standard deviations, $6\,\mathrm{STD}$. Each peak in the Fourier domain then has a standard deviation of

$$\mathrm{STD}\big\vert_h = \frac{6}{2\pi(\sigma_{max}-\sigma_{min})}$$

The Fourier transform of a Gaussian inverts its width, and on an $h$ axis defined through the phase $2\pi h\sigma$ a factor of $2\pi$ comes along. The idea of setting a lower bound this way comes from Seung Woo Lee's dissertation (2021). Requiring the centers of the dc peak and the sidebands to sit at least $3\,\mathrm{STD}$ apart gives

$$L > \frac{36}{2\pi(\sigma_{max}-\sigma_{min})}$$

Over the visible range 400–800 nm, $\sigma$ runs from $1.25$ to $2.5\ \mu m^{-1}$, so $L > 4.6\ \mu m$ and a quartz plate must exceed **about 0.46 mm**. Post 2 put the zero-order quarter-wave plate at 17.6 μm; this is roughly 26 times thicker. That is the scale meant by a "thick" plate in channeled spectroscopy. (Lee's Eq. 3.26 writes the standard deviation as $6/(\sigma_{max}-\sigma_{min})$ and drops the $2\pi$. With the transform kernel $e^{-2\pi jh\sigma}$ defined in his own Eq. 3.9, a further division by $2\pi$ is required.) Lee's dissertation arrives at a minimum of $3200\ \mu m$, but besides the $2\pi$ it takes the birefringence of quartz as $0.0045$ and the peak separation as $1.5\,\mathrm{STD}$, so the two numbers are not directly comparable.

Real instruments are built well above this bound. Hagen's experimental instrument uses quartz plates of 3.89 mm and 7.77 mm; Oka and Kato used 13 mm and 39 mm. That is 200 to more than 2000 times a zero-order plate. Feeding Hagen's conditions (400–1035 nm, 1.01 nm resolution, thickness ratio 1:2) into this lower bound and the upper bound of Section 3 puts the base thickness between 0.37 and 5.9 mm. His 3.89 mm sits inside the upper part of that window.

Stated generally, the condition is that **the bandwidth of the spectral structure the sample produces must fit inside one channel.** Hu and colleagues write it as

$$BW_{s,0,max} < B_{\rm eff}\, t_0$$

Here $t_0$ is the base thickness, the unit shared by the retarder thicknesses ($d_1 = m t_0$ and $d_2 = n t_0$ for a ratio $m:n$). The faster $\Psi$ and $\Delta$ vary with wavelength, the wider the band, and the thicker the retarder has to be.

The essential bargain of the method appears here. **Spectral resolution on the sample is traded for polarization channels.** Oka states the same thing as a condition: the polarization state must vary slowly compared with the carrier modulation.

## 3. Why not simply make the retarder thicker?

Does that mean thicker is always better? No.

A thicker retarder packs the fringes closer together. Once they are finer than the pixel spacing of the spectrometer, the fringes cannot be recorded and the sampling theorem produces aliasing. Writing $f_0 = B_{\rm eff}\,t_0$ for the carrier frequency of the base thickness, the highest carrier frequency in an arrangement measuring the full Stokes vector is $(m+n+\tfrac12)f_0$ (the $\tfrac12$ is the half-width of the highest channel), so the condition becomes

$$\frac{1}{\Delta\sigma} \ge \big[2(m+n)+1\big]\,B_{\rm eff}\,t_0$$

where $m : n$ is the thickness ratio of the two retarders and $\Delta\sigma$ is the wave-number resolution of the spectrometer.

One caution applies when computing $\Delta\sigma$. A dispersive spectrometer samples wavelength uniformly, so **its wave-number resolution is not uniform across the band.** Since $\Delta\sigma = \Delta\lambda/\lambda^2$, the short-wavelength end is worse. The geometric mean serves as the representative value.

$$\Delta\sigma_{avg} = \frac{\Delta\lambda}{\lambda_{min}\lambda_{max}}$$

Hagen derives it as the geometric mean of the resolutions at the two ends of the band, and Hu quotes it in this compact form.

Put in numbers. Over 400–800 nm, a spectrometer with 1 nm resolution and a 1:2 thickness ratio allows a base thickness of up to about 4.6 mm, ten times the 0.46 mm lower bound of Section 2. Coarsen the resolution to 2 nm and the upper bound falls to 2.3 mm.

Combining the two inequalities leaves a window of allowed thicknesses.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig2-thickness-window.png" alt="The allowed thickness window" width="780">
_Figure 2. The range of base thickness $$t_0$$ for quartz retarders over 400–800 nm, against the spectrometer's wavelength resolution. The green region under the solid line is usable by a channeled spectropolarimeter (CSP) measuring the Stokes vector with two retarders (ratio 1:2); the darker region under the dashed line is usable by a four-retarder Mueller arrangement (ratio 2:1:7:15). The darker region closes at 1.37 nm._

The lower bound is a horizontal line. It is fixed by the source bandwidth, so it has nothing to do with spectrometer performance. The upper bound slopes: a finer spectrometer raises it.

**Going to four retarders to measure the full Mueller matrix drops the upper bound by a factor of $7.3$.** With the ratio 2:1:7:15 proposed by Hagen, the highest carrier frequency jumps from $m+n = 3$ to $m+n+p+q = 25$ times the base frequency ($p$ and $q$ are the ratios of the third and fourth retarders; both multiples are taken before adding the half-width $\tfrac12$). What matters is where the lower bound meets the dashed Mueller line. Past that point, a spectrometer coarser than $1.37$ nm leaves **no usable thickness at all.** The window closes.

In Post 3, going from one compensator to two in order to reach all sixteen Mueller elements cost speed. Moved onto the wavelength axis, the same cost appears as a **spectrometer specification**.

## 4. How are retarders arranged to measure the full Stokes vector?

One retarder yields $a$ and $b$. That suffices for $\Psi$ and $\Delta$ of an isotropic sample but not for the full Stokes vector. Oka's original arrangement uses two retarders, their fast axes crossed at $45°$, with the analyzer aligned to the fast axis of the first.

Three quasi-cosinusoidal components then appear in the detected spectrum.

$$P(\sigma) = \tfrac12 S_0 + \tfrac12 S_1\cos(2\pi L_2\sigma + \cdots) + \tfrac14\lvert S_{23}\rvert\cos(2\pi(L_1-L_2)\sigma + \cdots) - \tfrac14\lvert S_{23}\rvert\cos(2\pi(L_1+L_2)\sigma + \cdots)$$

with $S_{23} = S_2 - iS_3$. **The carrier frequencies are $L_2$, $L_1-L_2$ and $L_1+L_2$**, so the Fourier domain holds seven peaks including dc. Choose the thicknesses well and they stay apart, each one extractable by filtering.

The structure matches Post 3. There the rotation ratio of the two compensators (fast:slow), 5:3, fixed the harmonic orders; here **the thickness ratio $d_1 : d_2$ fixes the carrier frequencies.** What the time axis did is repeated on the wavelength axis.

So how is the ratio chosen? Spacing the seven peaks evenly along the OPD axis makes $1:2$ optimal. Real instruments often use $3:1$ instead. The reason lies somewhere an idealized calculation never shows.

## 5. Why do real instruments use a 3:1 thickness ratio instead of 1:2?

The palm-sized instrument Okabe and colleagues built in 2009 uses calcite retarders. The effective birefringence of calcite is about twenty times that of quartz, so the same $L$ comes from a plate about one twentieth as thick, which is what shrinks the sensing head. Okabe's calcite plates of 150 µm and 450 µm correspond to roughly 2.9 mm and 8.6 mm of quartz.

There is a price. The ordinary and extraordinary indices of calcite differ substantially (Okabe uses $n_o = 1.6548$ and $n_e = 1.4864$ at 589 nm; the standard value of $n_o$ is $1.6584$). Even though the material itself has almost no diattenuation, **the Fresnel reflectance at the entrance and exit faces differs between the axes**, so the fast and slow transmittances split. At normal incidence through both faces, the fast axis transmits $0.925$ and the slow axis $0.882$, a difference of $4.6\%$. Defining the transmittance ratio angle as $\gamma = \tan^{-1}(t_s/t_f)$, an ideal retarder gives $45°$ while calcite gives $44.32°$. Okabe writes the angle with the intensity transmittances $T_s/T_f$, but the intensity ratio gives $43.64°$; only the amplitude ratio reproduces the paper's $44.32°$ exactly.

This imbalance creates a term the signal did not have: **one whose carrier frequency is $L_1$.** The original expression held only $L_2$ and $L_1 \pm L_2$; now there are four.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig3-thickness-ratio.png" alt="Channel layout for two thickness ratios" width="780">
_Figure 3. Channel layout for one input polarization seen through two calcite retarders, with ratio 1:2 on the left and 3:1 on the right. The dark curve assumes no diattenuation; the red curve includes the Fresnel imbalance of calcite. At 1:2 the two curves coincide and the new term stays hidden; at 3:1 it rises in the empty $$L_1$$ slot._

The problem is that the two curves in the left panel coincide. The new term hides inside an existing channel, and the intensity alone cannot tell where the original signal ends. With a thickness ratio $d_1 : d_2$ of $1:2$, $L_2 = 2L_1$, so the new $L_1$ term lands at $L_2 - L_1 = L_1$, **exactly on top of** the $L_1-L_2$ channel. The ratio that was optimal in the ideal calculation, with its even spacing, becomes unusable once the new $L_1$ term is taken into account. (Okabe places the retarders before the sample, which reverses the numbering, so he describes the same situation as $2:1$.)

Switching to $3:1$ puts the peaks at $L_2, 2L_2, 3L_2, 4L_2$, nine slots spaced evenly. And the slot at $L_1 = 3L_2$ is **one that should be empty if everything were ideal.** If a signal rises there, its size reports the diattenuation and misalignment of the retarders.

In Post 3 the dual rotating compensator solved fifteen Mueller elements from twenty-four Fourier coefficients, and the nine left over watched the instrument. In Post 4 the MSE of the fit revealed effects missing from the model. Here the empty channels do the same job. **A slot that should be empty watches the instrument.** The same idea runs through three posts.

The cost is resolution. Nine slots divide what seven would have, so each channel narrows. That is the subject of the next section.

## 6. How many resolvable points does each channel keep?

A spectrometer has a finite ability to resolve. Following Hagen's worked example, a spectrometer covering 400–1035 nm at $1.01$ nm gives

$$RP_\sigma = \frac{\sigma_{max}-\sigma_{min}}{\Delta\sigma_{avg}} \approx 629$$

meaning 629 resolvable points along the wave-number axis. Hagen reports 632 for the same example, computing the wave-number range with 1045 nm and the resolution with 1035 nm; using 1035 nm throughout gives 628.7. Those 629 points are shared out among the channels.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig4-channel-budget.png" alt="Resolvable points shared among channels" width="780">
_Figure 4. Points left per Stokes component for each arrangement in Hagen's example (400–1035 nm, 1.01 nm resolution). The "channels" label on each bar gives the number of channels the arrangement creates; with four retarders, 51 channels share the 629 points and leave 12._

One retarder measuring only $\Psi$ and $\Delta$ leaves 210 points per component. Two retarders for the full Stokes vector leave 90; the $3:1$ ratio forced by calcite leaves 70; going all the way to sixteen Mueller elements leaves **12**.

This is the exact price of a single shot. A rotating-element instrument can spend all 629 points on every Stokes component. A channeled instrument takes one frame and divides the points instead.

One option Hagen offers shows the trade clearly. In a $1:2$ arrangement, **pushing the outermost channel beyond the resolvable region and using only the inner two** raises the resolution of the reconstructed spectra by 50%. In exchange, the redundant information that channel carried disappears, dynamic calibration is no longer available, and the instrument must be recalibrated periodically while thermal drift is suppressed.

Discard the redundancy and gain resolution; keep it and gain self-calibration. This is a third fork after the "cancel by symmetry" and "measure by modelling" of Post 4.

## 7. What error does demodulation leave behind?

The act of cutting out a channel, demodulation, creates error of its own. After the two-sided thickness limit and the sharing of resolution points, demodulation error is the third constraint of channeled spectroscopy. The review by Zhang and Zhao also lists Fourier-domain crosstalk among the three fundamental constraints of the method.

Applying a window in the Fourier domain to isolate one channel does two things. First comes **crosstalk**: the skirts of neighbouring channels fall inside the window and mix in. The weaker the polarization signal, the more easily it is buried by that leakage. Second comes **degraded spectral resolution**: the window acts as a low-pass filter and shaves off the high-frequency content of the sample spectrum.

Okabe wrote this as a convolution. With $w(\sigma)$ the Fourier transform of the window function, the demodulated quantities become

$$a \;\longrightarrow\; a * w, \qquad b \;\longrightarrow\; b * w$$

and $\Psi$ and $\Delta$ are computed from these. The important point is that it cannot be undone. **It is not $\Psi$ and $\Delta$ that were convolved; they are computed nonlinearly from quantities that were.**

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig5-demodulation.png" alt="Demodulation distortion and retardance drift" width="780">
_Figure 5. Calculation for a SiO₂ film (index 1.46) on silicon at 70° incidence. Left, RMS error in $\Psi$ and $\Delta$ as the channel half-width is narrowed: green solid for the 100 nm film, red dashed for the 600 nm film, dark lines for $\Psi$ and pale lines for $\Delta$; the vertical dotted line marks $L/2$. Right, the $(a, b)$ locus of a 60 nm film rotating with the retardance drift $\delta\phi$._

Two things are readable on the left. Narrowing the window raises the error, and in this calculation (SiO₂ on silicon at 70°, films of 100 and 600 nm) **the thicker film suffers more.** A thick film makes $\Psi$ and $\Delta$ oscillate rapidly with wavelength, so its band is wide and more of it is cut away. What Section 2 called trading the sample's spectral resolution for polarization channels appears here as numbers. With the window half-width set to $L/2$, a 100 nm film shows $0.39°$ of error in $\Psi$ while a 600 nm film shows $2.27°$. $\Delta$ fares worse, at $7.0°$ and $15.6°$ under the same conditions. Using $\Delta$ through a window therefore effectively requires one of the remedies below.

The remedies split three ways.

| Approach | Method | Idea |
| --- | --- | --- |
| Add hardware | time or space division to capture several spectra | rotate a retarder and record several fringe patterns, or split the beam and record several at once. The windowing step disappears |
| Put it in the model | Okabe | do not try to remove it; include $w(\sigma)$ in the forward model and fit it alongside |
| Change the tool | coherent demodulation, compressed sensing, extrema envelopes | build a demodulation that never applies a Fourier window |

Okabe's remedy touches Post 4 directly. The goal of ellipsometry is usually not $\Psi$ and $\Delta$ themselves but a thickness or a dielectric function, obtained by numerical fitting. The effect of the window can therefore be folded into the theoretical model. It is the counterpart of Post 4's regression calibration and its "a model cannot account for what is not in it." **Rather than subtracting the instrument's flaw from the data, add it to the model.**

As an example of the third category, the author's dissertation (Chapter 3) introduced the extrema envelope method (EEM) to channeled spectra. A channeled spectrum is a baseline with a single harmonic riding on it, so taking the envelopes through the local maxima and minima and averaging them cancels the harmonic and leaves the baseline. Figure 6 shows the two envelopes. The method removes at the source the Gibbs and edge errors a Fourier filter produces at the signal boundaries. Averaging two envelopes to obtain a local mean is an idea Huang and colleagues used first, in empirical mode decomposition (EMD). It has a cost: because it locates extrema directly, it does worse than a Fourier low-pass filter when noise dominates. Metrology Image Processing 3 (the frequency domain) compares the two.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/ext-signal-envelopes.png" alt="Upper and lower envelopes of an oscillating signal" width="600">
_Figure 6. Upper (red) and lower (green) envelopes drawn over an oscillating signal (blue). An example of envelopes drawn over a generic oscillating signal, not a channeled spectrum. (Source: Brews ohare, [Signal envelopes](https://commons.wikimedia.org/wiki/File:Signal_envelopes.svg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

## 8. Which dials drift in a channeled instrument?

The thesis of Post 4 was that the angle a dial reports is not the true angle. The same problem arises here, but the dial that drifts is a different one.

The first is **retardance**. A thick plate is sensitive to temperature. The real plates of Section 2 are hundreds to thousands of times thicker than a zero-order plate, so the same temperature change shifts their retardance correspondingly more. Calibrating in advance does not help if the temperature changes between calibration and measurement.

Writing the drift as $\delta\phi$, the two measured quantities move as

$$\begin{bmatrix}a'\\ b'\end{bmatrix}
= \begin{bmatrix}\cos\delta\phi & \sin\delta\phi\\ -\sin\delta\phi & \cos\delta\phi\end{bmatrix}
\begin{bmatrix}a\\ b\end{bmatrix}$$

a **rotation** by $\delta\phi$ in the $(a, b)$ plane. The right panel of Figure 5 shows it. Because the film's locus lies close to the unit circle, the rotation looks like the locus sliding along the circle: the 800 nm end moves to positions turned by $\delta\phi = 10°$ and $30°$.

Section 2 of Post 4 described how the measured pair $(\alpha, \beta)$ of a rotating-element instrument rotates by $2A_s$ because of the analyzer offset. The matrix has the same form. **In a channeled instrument the retardance drift plays that role.** What must be undone has moved from an angle to a retardance; the structure is unchanged.

The way to undo it also resembles Post 4. In the two-retarder arrangement, the three channels at different carrier frequencies pick up different phase factors from the drift. Solving that difference **recovers $\delta\phi$ from the measured values alone.** Since both retarders share a material and a temperature, $\delta\phi_2 = (d_2/d_1)\,\delta\phi_1$, so knowing one corrects both. No thermometer and no temperature controller are required. The redundant channels watch the instrument, once again.

The second dial runs deeper: **the wavelength scale itself**.

A channeled instrument writes its information into fringes along the wavelength axis. The moment the pixel-to-wavelength mapping drifts, the carrier phase is corrupted directly. A group at Jeonbuk National University identified exactly this while working on an interferometric snapshot ellipsometer. The ordinary wavelength calibration, matching known lamp lines to the pixels of the two spectrometers, was not enough. The spectrometers had to be calibrated against the dense peaks produced by the instrument's own Mach–Zehnder interferometer.

Post 4 dealt with an angular dial; Post 5 deals with a retardance and a wavelength scale. Change the instrument and the dial you must distrust changes with it.

## 9. So what was gained?

Having spent this long on the costs, the gains deserve equal space.

Okabe's instrument has a palm-sized sensing head of $220 \times 45 \times 30$ mm and an acquisition time of **20 ms**. The paper gives the reason in one line. The configuration is simple and **there are no mechanical or active components for polarization control**. The performance holds up too. Twelve SiO$_2$ films spanning 3 to 4000 nm agreed with a commercial rotating-compensator instrument, and reference samples matched their certified values. **Thickness measurements stayed within $0.11$ nm as the temperature was varied from 5 to 45 °C.**

The author's dissertation reached similar results. An instrument that adds the back-focal-plane angle resolution of Post 6 to channeled spectroscopy measured single and multilayer films of about 10–1500 nm on silicon. A single shot of under 0.05 s, with no mechanical motion, gave both thickness and refractive index; thickness repeatability was at most 0.39 nm at three standard deviations, and agreement with a commercial ellipsometer reached $R^2 = 0.9999$.

It is worth noting where the field is heading. The doubt often raised about the snapshot family is accuracy: the belief that it falls well short of commercial spectroscopic ellipsometry. Recent work shows that doubt to be wrong. The snapshot approach of Wang and colleagues, using the back focal plane, brings the thickness RMS error down to $0.2$ nm against a commercial instrument and, at $1$ ms exposure, maps 60,000 points on a 4-inch wafer in 600 s, more than a hundred times faster than a commercial spectroscopic ellipsometer.

According to Wang and colleagues, with 3D NAND at 232 tiers and CFET stacking twenty to thirty layers under 10 nm each, the number of layers whose thickness and composition must be verified has exploded. An instrument that holds its accuracy while taking one frame is what breaks that bottleneck.

## Summary

Where a rotating-element instrument scans its modulator along the time axis, a channeled instrument lets the spectrometer's wavelength dispersion do it. A thick birefringent plate becomes a variable retarder once it meets a spectroscopic device, and the polarization goes into the fringes riding on a single spectrum.

The cost is quantitative. The retarder thickness is bounded on both sides: too thin and the channels overlap, too thick and they alias. The spectrometer's resolvable points are divided among the channels, so each Stokes component keeps a smaller share, and reaching all sixteen Mueller elements turns 629 points into 12. The demodulation that cuts out a channel leaves distortion of its own.

And there is the warning Oka set down in 1999. **A single shot is not the same as a fast one.** It was written for scanning spectrometers: capturing the fine structure of the channels required many samples, and the scan took correspondingly longer. Okabe's instrument, with a multichannel spectrometer, takes 20 ms, so speed now comes too. The essential gain is still **the absence of moving parts**.

The next post lifts the restriction of viewing the sample at one angle. The back focal plane of a high-numerical-aperture objective spreads the angles of incidence out at once.

## References

- K. Oka and T. Kato, "Spectroscopic polarimetry with a channeled spectrum," *Opt. Lett.* **24**, 1475 (1999) — the origin of the method.
- N. Hagen, "Design of channeled spectropolarimeters," *Appl. Opt.* **61**, 3381 (2022) — design guidelines for thickness, tolerancing and the resolution budget.
- H. Okabe, M. Hayakawa, J. Matoba, H. Naito, K. Oka, "Error-reduced channeled spectroscopic ellipsometer with palm-size sensing head," *Rev. Sci. Instrum.* **80**, 083104 (2009) — the 3:1 thickness ratio and self-calibration.
- J. Hu, X. Chen, W. Chen, S. Yang, Y. Wang, Z. Tang, S. Liu, "Frequency properties of channeled spectropolarimetry: an information theory perspective," *Opt. Express* **32**, 3735 (2024) — the general form of the thickness window.
- B. Zhang and B. Zhao, "Channeled spectropolarimetry: A review of technological evolution, algorithmic breakthroughs, and diversifying applications," *Opt. Lasers Eng.* **196**, 109421 (2026) — twenty-six years reviewed.
- V. Dembele, M. Jin, I. Choi, W. Chegal, D. Kim, "Interferometric snapshot spectro-ellipsometry," *Opt. Express* **26**, 1333 (2018), with the follow-up on calibration in *Curr. Opt. Photonics* **4**, 345 (2020) — the interferometric branch and the wavelength-scale problem.
- J. Wang, Q. Xu, L. Peng, J. Yang, H. Zhu, J. Zhu, Y. Shi, O. Zakharov, H. Jiang, M. Xu, J. Liu, S. Liu, "Snapshot Fourier ellipsometry: Pushing to sub-nanometer accuracy for high-throughput thin film metrology," *Adv. Sci. Instrum.* **1**, 100008 (2026) — accuracy and throughput of back-focal-plane snapshot ellipsometry, and the 3D-stacking background.
- Seung Woo Lee, "Co-axial spectroscopic snapshot ellipsometry for micro-spot measurement using high-frequency modulation and selective detection of spectral signals," Ph.D. dissertation, Seoul National University, 2021.
- Young Joon Kim, "Development of snapshot angle-resolved ellipsometry using a line-scan spectrograph and back focal plane spectral interference," Ph.D. dissertation, Seoul National University, 2025, Ch. 3 — the extrema envelope method.
- Y. J. Kim, S. J. Lee, H. J. Pahk, "Snapshot angle-resolved channeled spectroscopic micro-ellipsometry for thin-film characterization," *Appl. Opt.* **64**, 3391 (2025), [doi:10.1364/AO.558576](https://doi.org/10.1364/AO.558576) — the journal version of the dissertation.
- N. E. Huang et al., "[The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis](https://doi.org/10.1098/rspa.1998.0193)," *Proc. R. Soc. Lond. A* **454**, 903 (1998) — the envelope-mean idea.

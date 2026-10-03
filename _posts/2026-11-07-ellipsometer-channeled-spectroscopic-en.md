---
date: 2026-11-07 20:00:00 +0900
layout: post
title: "Ellipsometer Instrumentation 5 — Channeled Spectroscopy: Modulating in Wavelength Instead of Rotating"
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

Instead of rotating the retarder, make it thick. The retardance then differs from wavelength to wavelength, and the moment a spectrometer spreads the spectrum out, **a full period of the modulation appears within a single frame.** The wavelength dispersion of the spectrometer does the work the time axis used to do. A spectrum carrying such fringes is called a channeled spectrum.

If Post 1 asked what cannot be measured, Post 2 where the components stop working, Post 3 why one modulating element is not enough, and Post 4 how to work with limitations once they are known, this post asks **what a single shot costs.**

## 1. Start with one retarder

Take the simplest arrangement. A thick retarder sits behind the polarizer, the light reflects off the sample, and an analyzer receives it. With the polarizer at $0°$, the retarder at $45°$ and the analyzer at $-45°$, the detected spectrum reads

$$I(\sigma) = \frac{I_{in}(\sigma)}{4}\Big[1 - a\,\cos\phi(\sigma) - b\,\sin\phi(\sigma)\Big],
\qquad a = \cos 2\Psi,\quad b = \sin 2\Psi\,\sin\Delta$$

What this arrangement actually measures is the pair $a$ and $b$; $\Psi$ and $\Delta$ are computed from them. Keeping the two apart pays off later.

Here $\sigma = 1/\lambda$ is the wave number. Wave number replaces wavelength because the retardance is proportional to it.

$$\phi(\sigma) = 2\pi\,B_{\rm eff}\,t\,\sigma = 2\pi L\sigma$$

This is the $\delta = 2\pi d\,\Delta n/\lambda$ of Post 2 rewritten in wave number. The reading, however, is inverted. In Post 2 the wavelength dependence of retardance was a **defect**, which is why zero-order plates were preferred and multi-order plates avoided. Here that dependence is the **principle**. The plate is made deliberately thick to maximize it.

The slope $L = B_{\rm eff}\,t$ is the optical path difference (OPD), and this single number fixes where the channels sit. The coefficient $B_{\rm eff}$ is the first-order effective birefringence,

$$B_{\rm eff} = B(\bar\sigma) + \bar\sigma\,\frac{dB}{d\sigma}\Big\vert_{\bar\sigma}$$

which matters because it is not the birefringence $B = n_e - n_o$ itself but that value plus a dispersion term. For quartz it comes to about $0.00998$ rather than the nominal $0.009$ used in Post 2. Hagen writes this coefficient $\beta$; since Posts 1 and 4 already use that symbol for the Fourier coefficients of a rotating-element instrument, it appears here as $B_{\rm eff}$. Hagen notes that for common retarder materials the value runs **about 10% above** the nominal birefringence, and indeed $0.009 \times 1.1 = 0.0099$.

What happens when this spectrum is Fourier transformed along the wave-number axis?

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig1-single-retarder-channels.png" alt="Spectrum and Fourier channels for a single retarder" width="780">
_Figure 1. One sample ($\Psi = 35°$, $\Delta = 50°$) seen through retarders of two thicknesses. Left, the detected spectrum; right, its Fourier transform along the wave-number axis._

Because $\cos\phi$ and $\sin\phi$ each move to $\pm L$, peaks appear at $h = 0, \pm L$. The dc peak carries the source intensity and the sample reflectance; the $\pm L$ peaks carry $a$ and $b$. Only the latter two are wanted, so only those channels need to be cut out.

The trouble is in the upper row. A thin retarder gives a small $L$, and the three peaks overlap. Overlapping peaks cannot be separated.

## 2. Thickness is squeezed from below

How large must $L$ be for the peaks to stay apart? The answer follows from the fact that the peaks have width.

Treat the source intensity as a Gaussian along the wave-number axis and take its wavelength span as $6\,\mathrm{STD}$. Each peak in the Fourier domain then has a standard deviation of

$$\mathrm{STD}\big\vert_h = \frac{6}{2\pi(\sigma_{max}-\sigma_{min})}$$

The Fourier transform of a Gaussian inverts its width, and on an $h$ axis defined through the phase $2\pi h\sigma$ a factor of $2\pi$ comes along. Requiring the dc peak and the sidebands to sit $3\,\mathrm{STD}$ apart gives

$$L > \frac{36}{2\pi(\sigma_{max}-\sigma_{min})}$$

Over the visible range 400–800 nm, $\sigma$ runs from $1.25$ to $2.5\ \mu m^{-1}$, so $L > 4.6\ \mu m$ and a quartz plate must exceed **about 0.46 mm**. Post 2 put the zero-order quarter-wave plate at 17.6 μm; this is roughly 26 times thicker. That is the scale meant by a "thick" plate in channeled spectroscopy.

Stated generally, the condition is that **the bandwidth of the spectral structure the sample produces must fit inside one channel.** Hu and colleagues write it as

$$BW_{s,0,max} < B_{\rm eff}\, t_0$$

The faster $\Psi$ and $\Delta$ vary with wavelength, the wider the band, and the thicker the retarder has to be.

The essential bargain of the method appears here. **Spectral resolution on the sample is traded for polarization channels.** Oka states the same thing as a condition: the polarization state must vary slowly compared with the carrier modulation.

## 3. And from above

Does that mean thicker is always better? No.

A thicker retarder packs the fringes closer together. Once they are finer than the pixel spacing of the spectrometer, the fringes cannot be recorded and the sampling theorem produces aliasing. In an arrangement measuring the full Stokes vector the highest carrier frequency is $(m+n+\tfrac12)f_0$, so the condition becomes

$$\frac{1}{\Delta\sigma} \ge \big[2(m+n)+1\big]\,B_{\rm eff}\,t_0$$

where $m : n$ is the thickness ratio of the two retarders and $\Delta\sigma$ is the wave-number resolution of the spectrometer.

One caution applies when computing $\Delta\sigma$. A dispersive spectrometer samples wavelength uniformly, so **its wave-number resolution is not uniform across the band.** Since $\Delta\sigma = \Delta\lambda/\lambda^2$, the short-wavelength end is worse. The geometric mean serves as the representative value.

$$\Delta\sigma_{avg} = \frac{\Delta\lambda}{\lambda_{min}\lambda_{max}}$$

Hagen derives it as the geometric mean of the resolutions at the two ends of the band, and Hu quotes it in this compact form.

Combining the two inequalities leaves a window of allowed thicknesses.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig2-thickness-window.png" alt="The allowed thickness window" width="780">
_Figure 2. Quartz retarders over 400–800 nm. The horizontal axis is the spectrometer's wavelength resolution, the vertical axis the base thickness. The shaded regions are usable._

The lower bound is a horizontal line. It is fixed by the sample bandwidth and the source width, so it has nothing to do with spectrometer performance. The upper bound slopes: a finer spectrometer raises it.

What matters is where the two meet. **Going to four retarders to measure the full Mueller matrix drops the upper bound by a factor of $7.3$,** because the highest carrier frequency jumps from $m+n$ to $m+n+p+q$. As a result, a spectrometer coarser than $1.4$ nm leaves **no usable thickness at all.** The window closes.

In Post 3, going from one compensator to two in order to reach all sixteen Mueller elements cost speed. Moved onto the wavelength axis, the same cost appears as a **spectrometer specification**.

## 4. Two retarders for the full Stokes vector

One retarder yields $a$ and $b$. That suffices for $\Psi$ and $\Delta$ of an isotropic sample but not for the full Stokes vector. Oka's original arrangement uses two retarders, their fast axes crossed at $45°$, with the analyzer aligned to the fast axis of the first.

Three quasi-cosinusoidal components then appear in the detected spectrum.

$$P(\sigma) = \tfrac12 S_0 + \tfrac12 S_1\cos(2\pi L_2\sigma + \cdots) + \tfrac14\lvert S_{23}\rvert\cos(2\pi(L_1-L_2)\sigma + \cdots) - \tfrac14\lvert S_{23}\rvert\cos(2\pi(L_1+L_2)\sigma + \cdots)$$

with $S_{23} = S_2 - iS_3$. **The carrier frequencies are $L_2$, $L_1-L_2$ and $L_1+L_2$**, so the Fourier domain holds seven peaks including dc. Choose the thicknesses well and they stay apart, each one extractable by filtering.

The structure matches Post 3. There the rotation ratio $m_1 : m_2$ fixed the harmonic orders; here **the thickness ratio $d_1 : d_2$ fixes the carrier frequencies.** What the time axis did is repeated on the wavelength axis.

So how is the ratio chosen? Spacing the seven peaks evenly along the OPD axis makes $1:2$ optimal. Real instruments often use $3:1$ instead. The reason lies somewhere an idealized calculation never shows.

## 5. The slot that should be empty

The palm-sized instrument Okabe and colleagues built in 2009 uses calcite retarders. Its large birefringence delivers the same $L$ from a much thinner plate, which is what shrinks the sensing head.

There is a price. The ordinary and extraordinary indices of calcite differ substantially ($n_o = 1.6548$, $n_e = 1.4864$). Even though the material itself has almost no diattenuation, **the Fresnel reflectance at the entrance and exit faces differs between the axes**, so the fast and slow transmittances split. At normal incidence through both faces, the fast axis transmits $0.925$ and the slow axis $0.882$, a difference of $4.6\%$. Defining the transmittance ratio angle as $\gamma = \tan^{-1}(t_s/t_f)$, an ideal retarder gives $45°$ while calcite gives $44.32°$.

This imbalance creates a term the signal did not have: **one whose carrier frequency is $L_1$.** The original expression held only $L_2$ and $L_1 \pm L_2$; now there are four.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig3-thickness-ratio.png" alt="Channel layout for two thickness ratios" width="780">
_Figure 3. Channel layout for one input polarization seen through two calcite retarders. The dark curve assumes no diattenuation; the red curve includes the Fresnel imbalance of calcite._

The left panel is the problem. With a thickness ratio $d_1 : d_2$ of $1:2$, $L_2 = 2L_1$, so the new $L_1$ term lands at $L_2 - L_1 = L_1$, **exactly on top of** the $L_1-L_2$ channel. The ratio that was optimal in the ideal calculation, with its even spacing, becomes unusable once the new $L_1$ term is taken into account. (Okabe places the retarders before the sample, which reverses the numbering, so he describes the same situation as $2:1$.)

Switching to $3:1$ puts the peaks at $L_2, 2L_2, 3L_2, 4L_2$, nine slots spaced evenly. And the slot at $L_1 = 3L_2$ is **one that should be empty if everything were ideal.** If a signal rises there, its size reports the diattenuation and misalignment of the retarders.

> Post 3 noted that a dual rotating compensator solves fifteen Mueller elements from twenty-four surviving Fourier coefficients and has nine left over. Post 4 noted that the MSE reveals effects missing from the model. Here the empty channels do the same job. **A slot that should be empty watches the instrument** — the same idea, three posts running.

The cost is resolution. Nine slots divide what seven would have, so each channel narrows. That is the subject of the next section.

## 6. The channel budget

A spectrometer has a finite ability to resolve. Following Hagen's worked example, a spectrometer covering 400–1035 nm at $1.01$ nm gives

$$RP_\sigma = \frac{\sigma_{max}-\sigma_{min}}{\Delta\sigma_{avg}} \approx 629$$

meaning 629 resolvable points along the wave-number axis. Those 629 points are shared out among the channels.

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig4-channel-budget.png" alt="Resolvable points shared among channels" width="780">
_Figure 4. Points left per Stokes component for each arrangement. The number inside each bar is how many channels that arrangement creates._

One retarder measuring only $\Psi$ and $\Delta$ leaves 210 points per component. Two retarders for the full Stokes vector leave 90; the $3:1$ ratio forced by calcite leaves 70; going all the way to sixteen Mueller elements leaves **12**.

This is the exact price of a single shot. A rotating-element instrument can spend all 629 points on every Stokes component. A channeled instrument takes one frame and divides the points instead.

One option Hagen offers shows the trade clearly. In a $1:2$ arrangement, **pushing the outermost channel beyond the resolvable region and using only the inner two** raises the resolution of the reconstructed spectra by 50%. In exchange, the redundant information that channel carried disappears, dynamic calibration is no longer available, and the instrument must be recalibrated periodically while thermal drift is suppressed.

Discard the redundancy and gain resolution; keep it and gain self-calibration. This is a third fork after the "cancel by symmetry" and "measure by modelling" of Post 4.

## 7. The distortion demodulation leaves behind

The act of cutting out a channel creates error of its own. This is the third constraint of channeled spectroscopy and the one that has taken longest to address.

Applying a window in the Fourier domain to isolate one channel does two things. First comes **crosstalk**: the skirts of neighbouring channels fall inside the window and mix in. The weaker the polarization signal, the more easily it is buried by that leakage. Second comes **degraded spectral resolution**: the window acts as a low-pass filter and shaves off the high-frequency content of the sample spectrum.

Okabe wrote this as a convolution. With $w(\sigma)$ the Fourier transform of the window function, the demodulated quantities become

$$a \;\longrightarrow\; a * w, \qquad b \;\longrightarrow\; b * w$$

and $\Psi$ and $\Delta$ are computed from these. The important point is that it cannot be undone. **It is not $\Psi$ and $\Delta$ that were convolved; they are computed nonlinearly from quantities that were.**

<img src="/assets/img/posts/ellipsometer-channeled-spectroscopic/en/fig5-demodulation.png" alt="Demodulation distortion and retardance drift" width="780">
_Figure 5. Left, RMS error in $\Psi$ and $\Delta$ as the channel window is narrowed. Right, how the measured pair moves when the retardance drifts._

Two things are readable on the left. Narrowing the window raises the error, and **a thicker film suffers more.** A thick film makes $\Psi$ and $\Delta$ oscillate rapidly with wavelength, so its band is wide and more of it is cut away. What Section 2 called trading the sample's spectral resolution for polarization channels appears here as numbers. With the window half-width set to $L/2$, a 100 nm film shows $0.39°$ of error in $\Psi$ while a 600 nm film shows $2.27°$.

The remedies split three ways.

| Approach | Method | Idea |
| --- | --- | --- |
| Add hardware | time or space division to capture several spectra | rotate a retarder and record several fringe patterns, or split the beam and record several at once. The windowing step disappears |
| Put it in the model | Okabe | do not try to remove it; include $w(\sigma)$ in the forward model and fit it alongside |
| Change the tool | coherent demodulation, compressed sensing, extrema envelopes | build a demodulation that never applies a Fourier window |

Okabe's remedy touches Post 4 directly. The goal of ellipsometry is usually not $\Psi$ and $\Delta$ themselves but a thickness or a dielectric function, obtained by numerical fitting. The effect of the window can therefore be folded into the theoretical model. It is the counterpart of Post 4's regression calibration and its "a model cannot account for what is not in it." **Rather than subtracting the instrument's flaw from the data, add it to the model.**

A recent entry in the third category is the extrema envelope method. A channeled spectrum is a baseline with a single harmonic riding on it, so taking the envelopes through the local maxima and minima and averaging them cancels the harmonic and leaves the baseline. It removes at the source the Gibbs and edge errors a Fourier filter produces at the signal boundaries.

## 8. Another dial goes out of true

The thesis of Post 4 was that the angle a dial reports is not the true angle. The same problem arises here, but the dial that drifts is a different one.

The first is **retardance**. A thick plate is sensitive to temperature. At hundreds of times the thickness of a zero-order plate, the same temperature change shifts its retardance correspondingly more. Calibrating in advance does not help if the temperature changes between calibration and measurement.

Writing the drift as $\delta\phi$, the two measured quantities move as

$$\begin{bmatrix}a'\\ b'\end{bmatrix}
= \begin{bmatrix}\cos\delta\phi & \sin\delta\phi\\ -\sin\delta\phi & \cos\delta\phi\end{bmatrix}
\begin{bmatrix}a\\ b\end{bmatrix}$$

a **rotation** by $\delta\phi$ in the $(a, b)$ plane. The right panel of Figure 5 shows it.

Section 2 of Post 4 described how the measured pair $(\alpha, \beta)$ of a rotating-element instrument rotates by $2A_s$ because of the analyzer offset. The matrix has the same form. **In a channeled instrument the retardance drift plays that role.** What must be undone has moved from an angle to a retardance; the structure is unchanged.

The way to undo it also resembles Post 4. The three channels at different carrier frequencies pick up different phase factors from the drift. Solving that difference **recovers $\delta\phi$ from the measured values alone.** Since both retarders share a material and a temperature, $\delta\phi_2 = (d_2/d_1)\,\delta\phi_1$, so knowing one corrects both. No thermometer and no temperature controller are required. The redundant channels watch the instrument, once again.

The second dial runs deeper: **the wavelength scale itself**.

A channeled instrument writes its information into fringes along the wavelength axis. The moment the pixel-to-wavelength mapping drifts, the carrier phase is corrupted directly. A group at Jeonbuk National University identified exactly this while working on an interferometric snapshot ellipsometer. Ordinary wavelength calibration — matching the known lines of a calibration lamp to the pixels of its two spectrometers — was not enough; the spectrometers had to be calibrated against the dense peaks produced by the instrument's own Mach–Zehnder interferometer.

Post 4 dealt with an angular dial; Post 5 deals with a retardance and a wavelength scale. Change the instrument and the dial you must distrust changes with it.

## 9. What was gained

Having spent this long on the costs, the gains deserve equal space.

Okabe's instrument has a sensing head of $220 \times 45 \times 30$ mm — palm sized — and an acquisition time of **20 ms**. The paper gives the reason in one line: the configuration is simple and **there are no mechanical or active components for polarization control**. The performance holds up too. Twelve SiO$_2$ films spanning 3 to 4000 nm agreed with a commercial rotating-compensator instrument, and reference samples matched their certified values. **Thickness measurements stayed within $0.11$ nm as the temperature was varied from 5 to 45 °C.**

It is worth noting where the field is heading. The doubt the snapshot family has long carried is accuracy — the belief that it falls well short of commercial spectroscopic ellipsometry. Recent work confronts that directly. The snapshot approach of Wang and colleagues, using the back focal plane, brings the thickness RMS error down to $0.2$ nm against a commercial instrument and, at $1$ ms exposure, maps 60,000 points on a 4-inch wafer in 600 s, more than a hundred times faster than a commercial spectroscopic ellipsometer.

With 3D NAND at 232 tiers and CFET stacking twenty to thirty layers under 10 nm each, the number of layers whose thickness and composition must be verified has exploded. An instrument that holds its accuracy while taking one frame is what breaks that bottleneck.

## Summary

Where a rotating-element instrument scans its modulator along the time axis, a channeled instrument lets the spectrometer's wavelength dispersion do it. A thick birefringent plate becomes a variable retarder once it meets a spectroscopic device, and the polarization goes into the fringes riding on a single spectrum.

The cost is quantitative. Thickness is squeezed from both sides — too thin and the channels overlap, too thick and they alias. The spectrometer's resolvable points are divided among the channels, so each Stokes component keeps a smaller share, and reaching all sixteen Mueller elements turns 629 points into 12. The demodulation that cuts out a channel leaves distortion of its own.

And there is the warning Oka set down in 1999. **A single shot is not the same as a fast one.** Capturing the fine structure of the channels requires a large number of samples, and the spectrometer scan takes correspondingly longer. What is gained is not speed but **the absence of moving parts**.

The next post lifts the restriction of viewing the sample at one angle. The back focal plane of a high-numerical-aperture objective spreads the angles of incidence out at once.

## Reproducing the calculations

Every number and figure in this post is computed from products of component Mueller matrices alone. Closed-form expressions from the literature are kept out of the synthesis and reserved as an independent path for the verification script.

- `channeled.py` — channeled spectrum synthesis, Fourier demodulation, thickness limits, a single-film model
- `verify_channeled.py` — checks against Oka's Eq. 4, Hagen's Eq. 1, Okabe's $44.320°$ and Hu's thickness inequality. All eleven groups pass
- `generate_figures.py` — Figures 1–5 in both languages

The full code sits in `_code/ellipsometer-channeled-spectroscopic/`.

The cross-check again surfaced problems on the literature side.

Okabe writes the transmittance ratio angle as $\gamma = \tan^{-1}(T_s/T_f)$, but substituting intensity transmittances gives $43.640°$. **Only the amplitude ratio reproduces the paper's $44.320°$ exactly.** Hagen's resolution budget is inconsistent within the paper: the wave-number range is computed with $\lambda_{max} = 1.045\ \mu m$ while the resolution uses $1035$ nm, giving $632$. Using $1.035\ \mu m$ throughout gives $628.7$. The $B_{\rm eff} = 0.00998$ quoted in that same paper works out to $0.01010$, $1.3\%$ larger, when back-calculated from its own measured channel positions.

Two convention traps turned up as well. In the single-retarder arrangement, putting the analyzer at $+45°$ flips the sign of the $\sin\Delta$ term; sweeping the azimuths confirmed that $-45°$ is correct. Hagen's Eq. 1 is written with $R(45°)$, but under this series' rotation convention his expression matches $-45°$ — the handedness is opposite. The convention of Posts 1 through 4 was left alone and matched only for the comparison.

The thickness lower bound derived in Lee's dissertation was reproduced in spirit but recomputed. The dissertation takes the birefringence of quartz as $0.0045$ and obtains a minimum thickness of $3200\ \mu m$, whereas the nominal value ($0.009$), Hagen's coefficient ($0.00998$) and the value back-calculated from Hagen's measurements ($0.0101$) are all roughly double that. This post uses $B_{\rm eff} = 0.00998$ and also takes the peak separation as $3\,\mathrm{STD}$ rather than the dissertation's $1.5\,\mathrm{STD}$. Both were changed, so the number here is not directly comparable with the one in that work. The dissertation's Eq. 3.26 writes the $h$-domain standard deviation as $6/(\sigma_{max}-\sigma_{min})$, but with the transform kernel $e^{-2\pi jh\sigma}$ of its own Eq. 3.9 it must be divided by a further $2\pi$; that is corrected here as well.

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

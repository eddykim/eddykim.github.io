---
title: "Metrology Image Processing 3 — The Frequency Domain: Periodic Noise, Ringing, and the Extrema Envelope"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-frequency-domain-filtering/
page_id: imgproc-frequency-domain-filtering
date: 2026-10-14 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, fourier-transform, notch-filter, gibbs-ringing, spectral-leakage, envelope]
description: "The frequency domain erases periodic noise the spatial domain cannot touch, at the price of a periodicity assumption and ringing. The dissertation's extrema envelope method avoids paying it."
math: true
---

[Post 2](/en/posts/imgproc-convolution-kernels/) showed that the arguments for choosing a kernel all reduce to frequency response. Why not design the filter in the frequency domain directly, then? Two things are genuinely impossible in the spatial domain and become easy there: erasing periodic noise spread across the whole image, such as illumination ripple or reflections off a nearby grating, and specifying exactly which band to pass.

It is not free. The Fourier transform assumes the signal repeats periodically outside its own extent, and for real images that assumption is almost always false. This post puts what the frequency domain offers next to what it demands in return, and ends with a method that avoids paying the price at all, drawn from the dissertation.

## 1. Reading a spectrum

What does a single point in a two-dimensional spectrum mean? Its direction from the origin gives the orientation of a fringe pattern in the image; its distance from the origin gives how closely that pattern is spaced. One fringe pattern becomes one point.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig1-reading-the-spectrum.png" alt="Amplitude spectra of a synthetic grating and a natural image" width="700">
_Fig 1. Orientation and period appear as a single location in the spectrum_

The synthetic grating on the top row superimposes a tilted fringe pattern on a horizontal one. Its spectrum carries two pairs of points facing each other across the origin. The spectrum of a real-valued image is symmetric about the origin, so one fringe pattern yields two points.

The natural image below behaves differently. Energy concentrates at low frequencies and the rest spreads out broadly, because natural images have no dominant period. That contrast is what the next section rests on: periodic noise gathers into narrow points while image content spreads out, so the frequency domain can separate them.

## 2. Periodic noise — what the spatial domain cannot do

Suppose ripple from the illumination supply, or a reflection from an adjacent grating, has laid an even fringe pattern across the whole image. Can a kernel from Post 2 remove it? No. A smoothing kernel destroys image detail along with the fringes, and if the fringe period exceeds the kernel width it does nothing at all. Fringes span the entire image, so they are not a target for a local operation.

In the spectrum those fringes are two points. Zeroing just those locations and transforming back is all it takes. This is a notch filter.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig2-notch-filter.png" alt="Notch removal of periodic noise, with the frequency on and off a Fourier bin" width="750">
_Fig 2. Stripes untouchable in the spatial domain are two points in the spectrum_

Fringes of amplitude 38 DN were added, bringing the contamination to 26.87 DN RMS, and a notch was applied. The outcome splits in two.

The third panel is the case where the fringe frequency falls exactly on a Fourier bin ($44/512$, $72/512$). The RMS after notching is 0.76 DN — and applying the same notch to the uncontaminated original also costs exactly 0.76 DN. The fringes vanished without trace, leaving only what the notch itself removed. That is a 35-fold improvement.

The fourth panel is the case where the frequency lands between bins ($43.52/512$, $71.68/512$). The RMS stalls at 5.48 DN and fringes remain visible in the corners: seven times the collateral damage survives as residue. The improvement is only 4.9-fold.

The fringe strength and the notch size are identical, yet the results differ sevenfold. What separated them? Section 4 answers that.

## 3. Cutting frequency with a knife leaves rings in space

If the appeal of the frequency domain is exact control over the passband, then an ideal low-pass filter that zeros everything above a cutoff ought to be optimal. Is it?

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig3-ideal-vs-gaussian.png" alt="Ideal low-pass filter compared with a Gaussian low-pass filter" width="750">
_Fig 3. The negative lobes of the ideal filter print a ring beside every edge_

The top row is the ideal filter. A clean disc in the frequency domain becomes, back in space, a central peak surrounded by rings of alternating sign. The first negative trough reaches $-13.2\%$ of the peak height. That ring is copied beside every edge in the image, so the result carries alternating bright and dark borders — visible around the tripod and the figure's outline.

The numbers say the same. The original spans $[0, 255]$ in pixel value; the ideally filtered result spans $[-29.2, 272.0]$. Values that never existed have appeared. For metrology that reads brightness directly, that alone is an error.

The Gaussian filter on the bottom row is Gaussian in space as well, so it has no negative trough at all ($-1.0\times10^{-18}$, at the level of numerical error). Its output spans $[3.5, 244.0]$, inside the original range.

Post 2 noted that a box kernel's frequency response leaves sidelobes. What appears here is the same phenomenon seen from the other side. A box in one domain is a sinc in the other, and the tails of that sinc are the sidelobes and the rings alike. That the Gaussian is the only function smooth in both domains produces the same conclusion across two posts. This oscillation is the Gibbs phenomenon.

## 4. The assumption the Fourier transform carries

What is the bright cross running through the middle of the spectrum? Does the specimen really contain horizontal and vertical structure that strong?

It does not. The cross comes from the computation, not the specimen. The Fourier transform treats the image as repeating periodically outside its extent, which means the right edge is joined to the left edge. If the two edges differ in brightness, a step appears there that exists nowhere in reality. That step is the cross.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig4-periodic-extension.png" alt="The spectral cross from periodic extension and its suppression by a Hann window" width="750">
_Fig 4. The cross is made by the periodicity assumption, not by the specimen_

The top and bottom edges of this image differ by 91.1 DN RMS. Multiplying by a Hann window, which presses the borders to zero and joins the two edges smoothly, removes the step. The brightness of the cross relative to its surroundings falls from 3.35 to 1.00 — the cross disappears entirely in a statistical sense. The profile on the right shows the same thing: the windowed curve in blue sits lower in the tails.

The same mechanism solves the puzzle from section 2. When a fringe frequency falls between bins, the fringe pattern does not fit a whole number of times across the image width either. Periodic extension cuts it, and that discontinuity scatters energy into neighbouring bins. This is spectral leakage. What should have been one point is smeared, so a small notch cannot contain it. Enlarging the notch reduces the residue but carves away image content, trading one error for another.

Windowing suggests itself as the remedy, but windows are not free either. Pressing the borders to zero discards the information there. That is tolerable when only the centre of the field matters. It is fatal when the entire band must be used — which is exactly the situation in the next section.

## 5. The problem of extracting a baseline

Consider a concrete case of how these defects bite in real metrology: the signal of the channeled spectroscopic ellipsometry treated in the dissertation.

Light passing a multi-order retarder acquires a phase retardance nearly linear in wavenumber. The detected spectrum is therefore a slowly varying baseline signal (BLS) with a single harmonic riding on top. Below, $\sigma$ denotes wavenumber, the reciprocal of wavelength $1/\lambda$. This follows the convention of spectroscopy and the notation of the dissertation; note that Post 2 used the same symbol for a smoothing width.

$$ I(\sigma) = I_{BLS}(\sigma) + A(\sigma)\cos\!\left(2\pi\,\Delta n d\,\sigma - B(\sigma)\right) $$

The baseline carries the source spectrum and the sample reflectance; the modulated signal (MDS) carries the polarization information. Separating the two accurately is where analysis begins. Since the baseline normalizes the modulated signal, an error in the baseline propagates directly into $\Psi$ and $\Delta$.

Conventional channeled spectral polarimetry performs this separation with a Fourier transform along the wavenumber axis followed by a low-pass filter: cut above the modulation frequency and the baseline remains. The implementation is simple, and it inherits every problem from the three preceding sections. The ringing of section 3 appears, and the periodic-extension assumption of section 4 breaks at both ends of the band. A spectrum has no reason to show equal reflectance at 400 nm and 770 nm, so the two ends almost always differ. On top of that, whatever content the baseline itself held above the cutoff is erased with the modulation.

## 6. The extrema envelope method

Can the two components be separated without entering the frequency domain? The extrema envelope method (EEM) proposed in the dissertation uses a single property: the modulated signal oscillates symmetrically about the baseline. At the maxima and minima,

$$ I_{max} = I_{BLS} + A, \qquad I_{min} = I_{BLS} - A $$

The source intensity is folded into the baseline here. The dissertation factors it out as $I_{in}(I_{BLS} \pm A)$, but it divides out during normalization, so the conclusion is unchanged. The mean of the two envelopes is therefore the baseline.

$$ I_{BLS}(\sigma) = \frac{I_{max}(\sigma) + I_{min}(\sigma)}{2} $$

No frequency domain means no periodicity assumption and no cutoff to choose.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig5-extrema-envelope.png" alt="The extrema envelope method in operation, with a local quadratic fit zoomed" width="750">
_Fig 5. The mean of the two envelopes is the baseline_

The procedure has three steps. The signal is differentiated numerically and the sign changes locate the extrema. A quadratic is then fitted around each one to refine its wavenumber and intensity to subpixel precision. Finally the maxima and the minima are each joined by a spline to form the upper and lower envelopes, which are averaged.

The zoomed panel explains why the middle step is needed. Taking the sample point as the extremum locks it to the pixel grid, and that grid error rides the envelope into the baseline. It is the same device used in section 5 of Post 2, where an edge position was interpolated to the vertex of a parabola, and it is needed for the same reason.

One principle was kept in the implementation: the envelopes are never extrapolated. Before the first extremum and after the last there is no basis for an envelope, so no value is returned there. Of 1400 points, 1293 remain valid. Marking an unusable region as unusable is safer in metrology than inventing values to fill the band to its edges.

```python
# core of baseline_extraction.py (full code: _code/imgproc-frequency-domain-filtering/)
def eem_baseline(x, y, half=3, detect_sigma=0.0):
    maxima, minima = find_extrema(x, y, half, detect_sigma)
    up = CubicSpline(maxima[:, 0], maxima[:, 1])
    lo = CubicSpline(minima[:, 0], minima[:, 1])
    valid = (x >= max(maxima[0, 0], minima[0, 0])) & (x <= min(maxima[-1, 0], minima[-1, 0]))
    out = np.full_like(x, np.nan)
    out[valid] = 0.5 * (up(x[valid]) + lo(x[valid]))   # mean of the two envelopes
    return out, ...
```

## 7. Which one is better

Both methods were applied to the same signal and compared against the true baseline.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig6-lpf-vs-eem.png" alt="Baseline extraction error of FT+LPF and EEM, and their dependence on noise" width="750">
_Fig 6. EEM wins at the borders but loses once noise grows_

The left panel is the error profile without noise. The Fourier low-pass error spikes at both ends of the band and oscillates inward as it decays — the Gibbs ringing from the periodic-extension discontinuity promised in section 4. Its RMS is 0.00654 overall but 0.01828 over the outer 5%, nearly three times higher. The extrema envelope error hugs zero and is indistinguishable in the same plot: 0.00003 overall and 0.00008 at the borders, more than 200 times smaller.

On this evidence EEM looks decisively better, but the right panel overturns that. Repeating the comparison at increasing noise makes the two curves cross. Fourier low-pass error fails to double while noise grows sevenfold, whereas the extrema envelope degrades roughly in proportion to noise. The crossover sits just under 10% of the modulation amplitude, and beyond it the Fourier low-pass is the better choice.

The reason is clear once you ask what EEM relies on. It trusts the position and height of each extremum outright, and noise is exactly what disturbs them. Worse, noise flips the sign of the derivative repeatedly, so spurious extrema pour in. With no protection at the detection step, a noise level of 0.005 turns 25 extrema into 43 and drives the RMS error to 0.128 — twenty times worse than the Fourier low-pass.

Smoothing was therefore applied at the detection step alone, the smoothing of Post 2. A fringe period spans 56 pixels, so a smoothing width of 6 pixels leaves the extrema essentially unflattened while eliminating the spurious ones. The position and intensity of each extremum are then measured on the original signal, not the smoothed one. Separating detection from measurement is what makes EEM practical.

The selection rule reduces to one statement. Where the signal-to-noise ratio is adequate and both ends of the band must be used, the extrema envelope is better. Where noise dominates, or the modulation frequency is cleanly separated, the Fourier low-pass is safer. The dissertation chose EEM because spectroscopic ellipsometry is the former case: the short-wavelength end of the band could not be discarded, and the optics delivered enough signal-to-noise in a single shot.

## Summary and what comes next

The frequency domain does two things the spatial domain cannot. It erases periodic noise spread across an entire image with a few points, and it specifies a passband exactly. The price comes in three parts. Cutting a band with a knife leaves rings in space; the periodic extension the Fourier transform assumes is almost always violated by real signals; and the leakage from that discontinuity means even a notch fails when the frequency does not land on a bin. The extrema envelope method avoids paying any of it, at the cost of fragility against noise. That neither option is free is the conclusion of this post.

The next post takes up edges and subpixel estimation. Post 2 interpolated an edge position to the vertex of a parabola while locating the optimal smoothing width, and this post refined extremum positions the same way. Post 4 confronts that interpolation directly: how something smaller than a pixel is measured, and why such estimates acquire a periodic bias that tracks the pixel grid.

## References

- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 4 (frequency-domain filtering, notch filters, ringing).
- F. J. Harris, "[On the use of windows for harmonic analysis with the discrete Fourier transform](https://doi.org/10.1109/PROC.1978.10837)," _Proceedings of the IEEE_, 66(1), 51-83, 1978 (the standard reference on window functions and spectral leakage).
- K. Oka, T. Kato, "[Spectroscopic polarimetry with a channeled spectrum](https://doi.org/10.1364/OL.24.001475)," _Optics Letters_, 24(21), 1475-1477, 1999 (channeled spectral polarimetry and Fourier-based component separation).
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the standard test image in Figs 1-4).
- Y. Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," Ph.D. dissertation, Seoul National University, 2025, secs. 3.2 and 3.3 (signal model and the extrema envelope method).

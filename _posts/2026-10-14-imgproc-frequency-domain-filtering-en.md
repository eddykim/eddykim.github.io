---
title: "Metrology Image Processing 3 — What Is Lost When the Fourier Transform Removes Stripes from an Image?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-frequency-domain-filtering/
page_id: imgproc-frequency-domain-filtering
date: 2026-10-14 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, fourier-transform, notch-filter, gibbs-ringing, spectral-leakage, envelope]
description: "What removing stripe noise in the frequency domain costs (ringing, periodic extension), compared with the extrema envelope method that avoids it."
math: true
---

[Post 2](/en/posts/imgproc-convolution-kernels/) showed that the arguments for choosing a kernel all reduce to frequency response. Why not design the filter in the frequency domain directly, then? Two things are genuinely impractical with a small spatial kernel and become easy there. One is erasing periodic noise spread across the whole image, such as illumination ripple or reflections off a nearby grating. The other is specifying exactly which band to pass.

It is not free. The Fourier transform assumes the signal repeats periodically outside its own extent, and for real images that assumption is almost always false. This post puts what the frequency domain offers next to what it demands in return. It ends with a method that avoids paying the price at all, drawn from the author's doctoral dissertation.

## 1. What does a single point in a spectrum point to?

What does a single point in a two-dimensional spectrum mean? Its direction from the origin is perpendicular to the fringes, the direction in which brightness varies; its distance from the origin gives how closely they are spaced. One fringe pattern becomes a pair of points facing each other across the origin.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig1-reading-the-spectrum.png" alt="Amplitude spectra of a synthetic grating and a natural image" width="700">
_Fig 1. Amplitude spectra, on a log scale, of a synthetic grating (top) and a natural image (bottom). The grating spectrum is cropped to the central ±64 bins, and its one-bin points are enlarged so that they can be seen. Each fringe pattern produces a pair of points across the origin, and the centre point is the mean brightness. In the natural image, energy gathers at low frequencies and the rest spreads out._

The synthetic grating on the top row superimposes a tilted fringe pattern on a horizontal one. Its spectrum carries two pairs of points facing each other across the origin. The points come in pairs because the spectrum of a real-valued image is symmetric about the origin. The centre point is the mean brightness of the image, its DC component.

Each cell of the spectrum is a frequency bin. In an image 512 pixels wide, a point moves one bin further from the origin each time one more fringe fits across the width. The grating was built so that the tilted fringes fit exactly 31 times across and 13 times down, and the horizontal fringes 51 times down. Each point therefore lands cleanly on a single bin. Section 4 shows what happens when that condition fails.

The natural image below behaves differently. Energy concentrates at low frequencies and the rest spreads out broadly, and the spread exists because a natural image has no period that gathers into a single point. That contrast is what the next section rests on. Periodic noise gathers into narrow points while image content spreads out, so the frequency domain can separate them.

## 2. How do you remove stripes spread across a whole image?

Suppose ripple from the illumination supply, or a reflection from an adjacent grating, has laid an even fringe pattern across the whole image. Can a kernel from Post 2 remove it? Not with a small kernel. A smoothing kernel destroys image detail along with the fringes, and if the fringe period exceeds the kernel width it does nothing at all. Fringes span the entire image, so they are not a target for a local operation.

In the spectrum those fringes are two points. Zeroing just those locations and transforming back is all it takes. This is a notch filter.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig2-notch-filter.png" alt="Notch removal of periodic noise, with the frequency on and off a Fourier bin" width="750">
_Fig 2. Periodic noise removed with a notch filter. From left: the image with 38 DN fringes added, its spectrum (circles: noise peaks), and the notch result when the fringe frequency falls on a bin ($44/512$, $72/512$) and between bins ($43.52/512$, $71.68/512$). With the same notch, on-bin fringes vanish (0.76 DN RMS), while between-bin fringes survive at the borders (5.48 DN RMS)._

Fringes with an amplitude of 38 digital numbers (DN) were added, bringing the contamination to 26.87 DN root mean square (RMS), and a notch was applied. The two circled points in the second panel are those fringes; each occupies a single bin and is hard to see without the marks. The outcome then splits in two.

The third panel is the case where the fringe frequency falls exactly on a Fourier bin ($44/512$, $72/512$). The RMS after notching is 0.76 DN, and applying the same notch to the uncontaminated original also costs exactly 0.76 DN. The fringes vanished without trace, leaving only what the notch itself removed. That is a 35-fold improvement.

The fourth panel is the case where the frequency lands between bins ($43.52/512$, $71.68/512$). The RMS stalls at 5.48 DN and fringes remain visible in the corners: seven times the collateral damage survives as residue. The improvement is only 4.9-fold.

The fringe strength and the notch size are identical, yet the results differ sevenfold. What separated them? Section 4 answers that.

## 3. Why does a sharp frequency cut leave rings beside edges?

If the appeal of the frequency domain is exact control over the passband, then an ideal low-pass filter that zeros everything above a cutoff frequency ought to be optimal. Is it?

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig3-ideal-vs-gaussian.png" alt="Ideal low-pass filter compared with a Gaussian low-pass filter" width="750">
_Fig 3. An ideal low-pass filter (top) compared with a Gaussian low-pass filter (bottom). The ideal filter cuts off at 0.06 cycles per pixel, and the Gaussian uses the same value as its standard deviation. From left: the filter in the frequency domain, a central cut through its impulse response, and the filtered image. The ideal filter's impulse response has a negative trough reaching 13.2% of the peak, which leaves a ring at every edge; the Gaussian has none._

The top row is the ideal filter. A clean disc in the frequency domain becomes, back in space, an impulse response with a central peak surrounded by rings of alternating sign. The first negative trough reaches $-13.2\%$ of the peak height. That ring is copied beside every edge in the image, so the result carries alternating bright and dark borders. This is called ringing, and it is plain around the tripod and the photographer's outline.

The numbers say the same. The original spans $[0, 255]$ in pixel value; the ideally filtered result spans $[-29.2, 272.0]$. Values that never existed have appeared. For metrology that reads brightness directly, that alone is an error.

The Gaussian filter on the bottom row is Gaussian in space as well, so it has no negative trough at all ($-1.0\times10^{-18}$, at the level of numerical error). Its output spans $[3.5, 244.0]$, inside the original range.

Post 2 noted that a box kernel's frequency response leaves sidelobes. What appears here is the same phenomenon seen from the other side. A box in one domain is a sinc in the other (for a circular filter, its two-dimensional counterpart $2J_1(x)/x$), and those tails are the sidelobes and the rings alike. This oscillation at a sharply cut boundary is the Gibbs phenomenon. Fig 4 shows it at a single one-dimensional step.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/ext-gibbs-phenomenon.png" alt="A square wave synthesized from 25 Fourier terms, with Gibbs overshoot beside each step" width="650">
_Fig 4. A square wave rebuilt from the first 25 terms of its Fourier series. Every step is flanked by an overshoot and a train of ripples. Adding more terms narrows the overshoot but leaves its height at about 9% of the step. (Source: Johannes Rössel, [Gibbs phenomenon 50](https://commons.wikimedia.org/wiki/File:Gibbs_phenomenon_50.svg), public domain)_

Because a Gaussian transforms into a Gaussian, it takes no negative values in either domain, and that property produces the same conclusion across two posts. For metrology that uses brightness values directly, the edge of the passband should therefore be tapered smoothly, as a Gaussian does, rather than cut sharply.

## 4. Where does the cross in the middle of the spectrum come from?

What is the bright cross running through the middle of the spectrum? Does the specimen really contain horizontal and vertical structure that strong?

It does not. The cross comes from the computation, not the specimen. The Fourier transform treats the image as repeating periodically outside its extent; this is called periodic extension. The right edge is joined to the left edge, and the bottom edge to the top. If two joined edges differ in brightness, a step appears there that exists nowhere in reality. That step is the cross.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig4-periodic-extension.png" alt="The spectral cross from periodic extension and its suppression by a Hann window" width="750">
_Fig 5. The spectrum of a natural image computed without a window (left) and after a Hann window (centre), with a vertical cut through each (right). The cuts are normalized to each spectrum's DC component and drawn on a log scale. Without a window a bright cross runs through the centre; with the Hann window it is greatly weakened and the tails of the cut (blue) sit lower. The cross is made by the periodicity assumption, not by the specimen._

The top and bottom edges of this image differ by 91.1 DN RMS. A Hann window, one kind of window function, presses the borders to zero and joins the two edges smoothly, which removes the step. The brightness of the cross relative to its surroundings falls from 3.35 to 1.00, so the cross disappears, although a faint vertical band remains in the centre panel. The profile on the right shows the same thing: the windowed curve in blue sits lower in the tails.

The same mechanism solves the puzzle from section 2. When a fringe frequency falls between bins, the fringe pattern does not fit a whole number of times across the image width either. Periodic extension cuts it, and that discontinuity scatters energy into neighbouring bins. This is spectral leakage. What should have been one point is smeared, so a small notch cannot contain it.

Would a larger notch help? The notch in Fig 2 has a radius of 0.012 cycles per pixel, about 6 bins. For the between-bin fringes, doubling the radius lowers the residue from 5.48 DN to 4.12 DN, but quadrupling it raises the residue again to 4.82 DN. At that size the notch is already carving 4.00 DN out of the image content by itself. Leakage, in other words, is not cured by a bigger notch; the 0.76 DN of the on-bin case is never reached.

Windowing suggests itself as the remedy, but windows are not free either. Pressing the borders to zero discards the information there. That is tolerable when only the centre of the field matters. It is fatal when the entire band must be used, which is exactly the situation in the next section.

## 5. What goes wrong when a baseline is extracted from a spectral signal?

Why must the baseline of a spectral signal be extracted accurately? Because any error in it passes straight into the polarization information being measured. The signal of channeled spectroscopic ellipsometry, treated in the author's dissertation, shows how.

Light passing a multi-order retarder acquires a phase retardance nearly linear in wavenumber. The detected spectrum is therefore a slowly varying baseline signal (BLS) with a single harmonic riding on top. Below, $\sigma$ denotes wavenumber, the reciprocal of wavelength $1/\lambda$. This follows the convention of spectroscopy and the notation of the dissertation; note that Post 2 used the same symbol for a smoothing width.

$$ I(\sigma) = I_{BLS}(\sigma) + A(\sigma)\cos\!\left(2\pi\,\Delta n d\,\sigma - B(\sigma)\right) $$

Here $A$ is the modulation amplitude, $B$ the phase carried by the sample, and $\Delta n d$ the product of the retarder's birefringence and thickness. The baseline carries the source spectrum and the sample reflectance; the modulated signal carries the polarization information. The split is this clean only in the azimuth configuration the dissertation chose (polarizer 0°, retarder 45°, analyser 45°); in a general configuration the baseline also contains polarization terms. Separating the two accurately is where analysis begins. Since the baseline normalizes the modulated signal, an error in the baseline propagates directly into the ellipsometric angles $\Psi$ and $\Delta$.

Conventional channeled spectral polarimetry performs this separation with a Fourier transform along the wavenumber axis followed by a low-pass filter, written FT+LPF from here on. Cut above the modulation frequency and the baseline remains. The implementation is simple, and it inherits every problem from the three preceding sections. The ringing of section 3 appears, and the periodic-extension assumption of section 4 breaks at both ends of the band. A spectrum has no reason to show equal reflectance at 400 nm and 770 nm, so the two ends almost always differ. On top of that, whatever content the baseline itself held above the cutoff is erased with the modulation.

## 6. Can a baseline be extracted without the Fourier transform?

Can the two components be separated without entering the frequency domain? Waves suggest a way. Halfway between the line through the crests and the line through the troughs lies mean sea level. If an oscillation is symmetric about its baseline, the mean of its upper and lower envelopes is that baseline.

The author's dissertation (Chapter 3) introduced the extrema envelope method (EEM) to channeled spectra on this idea. It uses a single property: the modulated signal oscillates symmetrically about the baseline. Taking the mean of the upper and lower envelopes as the local mean appeared earlier, in Huang et al.'s empirical mode decomposition (EMD) of 1998. At the maxima and minima,

$$ I_{max} = I_{BLS} + A, \qquad I_{min} = I_{BLS} - A $$

so the mean of the two envelopes is the baseline.

$$ I_{BLS}(\sigma) = \frac{I_{max}(\sigma) + I_{min}(\sigma)}{2} $$

The source intensity is folded into the baseline here. The dissertation factors it out as $$I_{in}\!\left(I_{BLS} \pm A\right)$$, but it divides out during normalization, so the conclusion is unchanged. The first equation is also an approximation that holds when the baseline and amplitude vary slowly within one fringe period. For the signal in this post its error is at most $2.4\times10^{-4}$, about 0.1% of the modulation amplitude.

No frequency domain means no periodicity assumption and no cutoff to choose. Two settings remain instead: the smoothing width used to detect extrema and the fitting window around each one. As section 7 shows, both are chosen relative to one fringe period.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig5-extrema-envelope.png" alt="The extrema envelope method in operation, with a local quadratic fit zoomed" width="750">
_Fig 6. The extrema envelope method applied to a synthetic channeled spectrum (grey, 1400 points, 25 fringes). Left: the upper and lower envelopes, splines through the maxima (red ×) and minima (blue ×), and their mean, the EEM baseline (green), which overlaps the true baseline (black dashed) indistinguishably. Right: one extremum zoomed; the dashed curve is a quadratic fitted to nine samples and × marks its vertex._

The procedure has three steps. The signal is differentiated numerically and the sign changes locate the extrema. A quadratic is then fitted to the nine samples around each one to refine its wavenumber and intensity to subpixel precision. Finally the maxima and the minima are each joined by a spline to form the upper and lower envelopes, which are averaged.

Why is the middle step needed? Taking the sample point as the extremum locks it to the pixel grid, and that grid error rides the envelope into the baseline. Without refinement the baseline RMS error for this signal is 0.00053; with it, 0.00003, about eighteen times smaller. The samples are dense enough, as the zoomed panel shows, that the difference is invisible by eye, yet this step removes most of the baseline error. It is the same device used in section 5 of Post 2, where an edge position was interpolated to the vertex of a parabola, and it is needed for the same reason.

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

## 7. Which should you use, the Fourier low-pass or the extrema envelope?

The answer depends on the noise. For this synthetic signal, when the noise is below about 9% of the mean modulation amplitude the extrema envelope is better; above it, the Fourier low-pass is. The grounds are a comparison of both methods on the same signal against the true baseline.

<img src="/assets/img/posts/imgproc-frequency-domain-filtering/en/fig6-lpf-vs-eem.png" alt="Baseline extraction error of FT+LPF and EEM, and their dependence on noise" width="750">
_Fig 7. Baseline extraction error of the Fourier low-pass (FT+LPF, red) and the extrema envelope method (EEM, green). Left: error against wavenumber without noise; the low-pass, cut at 17 bins, spikes at both ends of the band and oscillates inward as it decays. Right: RMS error against noise standard deviation (log scale), with the low-pass cutoff chosen as the best bin at each noise level. Once the noise exceeds about 9% of the mean modulation amplitude the two curves swap (the dotted line marks 0.02, the first level where they do)._

The left panel is the error profile without noise. The Fourier low-pass error spikes at both ends of the band and oscillates inward as it decays: the Gibbs ringing from the periodic-extension discontinuity promised in section 4. Over the range where EEM returns values (1293 of 1400 points) its RMS is 0.00654, but over the points of that range within the outer 5% of the band it is 0.01828, nearly three times higher. Including the band ends where EEM returns nothing, the low-pass error is 0.01737 overall and 0.05196 over the outer 5%. The extrema envelope error hugs zero and is indistinguishable in the same plot: 0.00003 over the same range and 0.00008 at its borders, more than 200 times smaller.

On this evidence EEM looks decisively better, but the right panel overturns that. Repeating the comparison at increasing noise makes the two curves cross. Fourier low-pass error fails to double while the noise is raised sevenfold, from 0.01 to 0.07, whereas the extrema envelope degrades roughly in proportion to noise. The crossover sits at a noise standard deviation of about 0.0196, roughly 9% of the mean modulation amplitude of 0.229. Beyond it the Fourier low-pass is the better choice. In this comparison, however, the low-pass cutoff was picked at each noise level by looking at the true baseline, and a single noise seed was used. In a real measurement, where the true value is unknown, the low-pass fares somewhat worse than this.

The reason is clear once you ask what EEM relies on. It trusts the position and height of each extremum outright, and noise is exactly what disturbs them. Worse, noise flips the sign of the derivative repeatedly, so spurious extrema pour in. With no protection at the detection step, a noise level of 0.005 turns roughly 25 maxima and 25 minima into 43 of each and drives the RMS error to 0.128, twenty times worse than the Fourier low-pass.

Smoothing was therefore applied at the detection step alone, the smoothing of Post 2. A fringe period spans 56 pixels, so the smoothing width was set to 6 pixels. This smoothing scales the fringe amplitude by 0.80, but it keeps the position and number of the extrema while eliminating the spurious ones. Each extremum's position and intensity are measured on the original signal, not the smoothed one, so the reduced amplitude never enters the result. Separating detection from measurement is what makes EEM practical.

Carried over to practice, the comparison on this synthetic signal says: where the signal-to-noise ratio (SNR) is adequate and both ends of the band must be used, the extrema envelope is better. Where noise dominates, or the modulation frequency is cleanly separated, the Fourier low-pass is safer. The author's dissertation chose EEM because of errors at the signal boundaries. It set out to avoid the Gibbs and edge errors that a Fourier low-pass produces at both ends of the band, and the loss of the baseline's own high-frequency content (sec. 3.3). It did not examine sensitivity to noise. The comparison in this section adds the caveat that the choice pays off only when the signal-to-noise ratio is adequate.

## Summary and what comes next

The frequency domain does two things that are impractical with a small spatial kernel. It erases periodic noise spread across an entire image with a few points, and it specifies a passband exactly. The price comes in three parts. Cutting a band sharply leaves rings in space. The periodic extension the Fourier transform assumes is almost always violated by real signals. And the leakage from that discontinuity means even a notch fails when the frequency does not land on a bin. The extrema envelope method avoids paying any of it, at the cost of fragility against noise. That neither option is free is the conclusion of this post.

The next post takes up edges and subpixel estimation. Post 2 interpolated an edge position to the vertex of a parabola while locating the optimal smoothing width, and this post refined extremum positions the same way. Post 4 confronts that interpolation directly: how something smaller than a pixel is measured, and why such estimates acquire a periodic bias that tracks the pixel grid.

## References

- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 4 (frequency-domain filtering, notch filters, ringing).
- F. J. Harris, "[On the use of windows for harmonic analysis with the discrete Fourier transform](https://doi.org/10.1109/PROC.1978.10837)," _Proceedings of the IEEE_, 66(1), 51-83, 1978 (the standard reference on window functions and spectral leakage).
- K. Oka, T. Kato, "[Spectroscopic polarimetry with a channeled spectrum](https://doi.org/10.1364/OL.24.001475)," _Optics Letters_, 24(21), 1475-1477, 1999 (channeled spectral polarimetry and Fourier-based component separation).
- N. E. Huang et al., "[The empirical mode decomposition and the Hilbert spectrum for nonlinear and non-stationary time series analysis](https://doi.org/10.1098/rspa.1998.0193)," _Proc. R. Soc. Lond. A_, 454(1971), 903-995, 1998 (the envelope mean and the limits of spline envelopes).
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the standard test image in Figs 1-3 and 5).
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025, secs. 3.2 and 3.3 (signal model and the extrema envelope method).

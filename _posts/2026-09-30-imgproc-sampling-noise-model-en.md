---
title: "Metrology Image Processing 1 — What a Pixel Value Actually Counts"
lang: en
lang-exclusive: ["en"]
permalink: /posts/imgproc-sampling-noise-model/
page_id: imgproc-sampling-noise-model
date: 2026-09-30 20:00:00 +0900
categories: [Computation, Image Processing]
tags: [image-processing, noise, shot-noise, photon-transfer-curve, sampling, aliasing]
description: Before extracting any signal from an image, fix what a single pixel value counts and how much it scatters. A photon transfer curve and an aliasing experiment settle both.
math: true
---

Texts on image processing usually open with the sentence "an image is a two-dimensional array." For pictures meant to be looked at, that is enough. For measurement it is not. If you do not know what the number in that array counts, and how large its uncertainty is, then no filter and no fit downstream can be assigned a precision. Claiming a thickness repeatability of 0.39 nm requires tracing that number all the way back to its origin.

This series therefore starts at the sensor rather than at kernels or edge detection. Two quantities are fixed in this post: what a pixel value is proportional to, and how much that value scatters between repeated exposures. Once both are known, every operation introduced in later posts can be evaluated for how much precision it gains or destroys.

## 1. What is a pixel value proportional to?

Suppose a pixel in a 12-bit image reads 1500. Fifteen hundred of what? Answering that means retracing the path from light to number.

<img src="/assets/img/posts/imgproc-sampling-noise-model/en/fig1-photon-to-dn-chain.png" alt="Conversion chain from photons to pixel values and the noise entering at each stage" width="750">
_Fig 1. From photons to pixel values, and where each noise source enters_

During the exposure, an average of $\mu_p$ photons reaches one pixel. Some fraction of them liberates electrons; that fraction is the quantum efficiency $\eta$. Dark current adds $\mu_d$ electrons that have nothing to do with light. The accumulated charge is converted to a voltage and then digitized by an analog-to-digital converter, and the final conversion ratio is the system gain $K$. The digitized value is called a digital number (DN), so the gain is expressed in DN per electron.

$$ \mu_y = K(\eta\,\mu_p + \mu_d) + y_0 $$

Here $y_0$ is an offset deliberately added so that dark regions do not clip at zero. The equation says something simple. A pixel value is proportional to the photon count, the constant of proportionality is $K\eta$, and one constant sits on top. The number 1500 means nothing on its own; only $K$ and $y_0$ translate it into a count of electrons.

Rather than a real camera, this post uses a synthetic sensor that implements the chain above directly in code. Estimation methods can only be checked against a known truth. The parameters were set to a quantum efficiency of 0.60, a gain of 0.25 DN/e$^-$, a read noise of 3.0 e$^-$, a full well of 15000 e$^-$, 12 bits, and an offset of 100 DN.

## 2. How much does that value scatter?

Photograph the same scene twice without changing anything, and the pixel values differ. By how much? The answer underpins the whole series.

The scatter has three origins. First, photon arrival is itself random. The electron count follows a Poisson distribution, so its variance equals its mean $\mu_e$. This is shot noise, and it cannot be removed as long as there is light. Second, the readout circuit contributes read noise, a Gaussian independent of signal level. Third, truncation to integers produces quantization noise. The three sources are mutually independent, so their variances simply add.

$$ \sigma_y^2 = K^2\sigma_d^2 + \sigma_q^2 + K(\mu_y - y_0) $$

Here $\sigma_d$ is the noise measured in dark regions, combining read noise with dark-current shot noise. A short exposure leaves no time for dark current to accumulate, so the second contribution effectively vanishes. At the 10 ms used here, dark current amounts to only 0.05 e$^-$, and $\sigma_d$ may be treated as the read noise from this point on.

What deserves attention is that the last term is not constant but **proportional** to the signal. The brighter the exposure, the larger the absolute noise. Brightness still pays off, because the signal grows faster than the noise. Writing out the signal-to-noise ratio (SNR) makes this explicit.

$$ \mathrm{SNR} = \frac{\mu_e}{\sqrt{\sigma_d^2 + \mu_e}} $$

Once the signal is large enough, $\sigma_d^2$ in the denominator becomes negligible and the SNR approaches $\sqrt{\mu_e}$. Quadrupling the signal only doubles the SNR. This square-root law is the practical criterion for setting exposure in metrology imaging.

## 3. Measuring the camera directly

Datasheets sometimes list $K$ and $\sigma_d$, but often they do not, or the listed values disagree with the hardware. Can both constants be recovered from images alone? The variance equation above already contains the answer. Plot $\sigma_y^2$ against $\mu_y - y_0$ and the result is a straight line whose slope is $K$ and whose intercept is $K^2\sigma_d^2 + \sigma_q^2$. Collecting means and variances under uniform illumination at varying exposure traces that line. This is the photon transfer curve.

The variance must not be computed as the spread across pixels in a single frame. Per-pixel sensitivity deviations would be folded in and inflate it. Instead, two frames are captured under identical conditions and subtracted. Whatever is fixed appears identically in both and cancels, leaving only the temporally varying noise. The variance of the difference is twice the original, so it is halved.

```python
# core of sensor_model.py (full code: _code/imgproc-sampling-noise-model/)
for p in photon_levels:
    a = sensor.capture(p, exposure_s, rng)
    b = sensor.capture(p, exposure_s, rng)
    means.append(0.5 * (a.mean() + b.mean()))
    variances.append(0.5 * (a - b).var())   # the fixed pattern cancels

slope, intercept = np.polyfit(signal_dn[m], var_dn[m], 1)
K = slope
sigma_dark_e = np.sqrt(intercept - sigma_q_sq) / K   # remove the quantization part
```

<img src="/assets/img/posts/imgproc-sampling-noise-model/en/fig2-photon-transfer-curve.png" alt="Photon transfer curve with a linear-range fit and a low-signal inset" width="700">
_Fig 2. Photon transfer curve — the slope gives the gain $K$, the intercept gives the dark noise_

The slope of the fit yields $K = 0.2511$ DN/e$^-$, recovering the true value of $0.2500$ DN/e$^-$ to within 0.4%. Obtaining a sensor constant from a single least-squares line is the backbone of the camera characterization procedure specified by EMVA Standard 1288. Fitting the line is itself the simplest instance of the least-squares problem treated in [Optimization 1](/en/posts/optimization-gradient-descent/).

The intercept requires one further step. The fitted intercept is $0.627$ DN$^2$, and it holds not only the dark noise but also the quantization noise $\sigma_q^2 = 1/12$ DN$^2$. Converting the raw intercept gives $3.15$ e$^-$, overshooting the true $3.01$ e$^-$ by 5%. Removing the quantization term first gives $2.94$ e$^-$, cutting the error to 2%. For a low-noise sensor whose intercept is below 1 DN$^2$, $1/12$ is not a negligible quantity.

The fitting range was deliberately restricted to low signals. The intercept is roughly $0.6$ DN$^2$, while the slope term climbs to 900 DN$^2$ at the bright end. Across three orders of magnitude, a fit spanning the bright data would bury the intercept in scatter. The inset of Fig 2 shows the region where the intercept actually resolves.

The sharp downturn at the right end of the curve marks saturation. It occurs at a signal of 3602 DN, or 14410 e$^-$, consistent with the full well of 15000 e$^-$. Charge stops accumulating, so the fluctuation stops as well. A pixel past this point cannot be used as a measurement, however smooth its value may look.

## 4. Does adding bits add precision?

A 12-bit camera clearly beats an 8-bit one. Does a 16-bit camera beat the 12-bit one? Each added bit halves the quantization step $\Delta$, and under a uniform-distribution assumption the quantization noise falls as $\sigma_q = \Delta/\sqrt{12}$. The difficulty is that this term is only one part of the total noise.

<img src="/assets/img/posts/imgproc-sampling-noise-model/en/fig3-bit-depth-vs-noise.png" alt="Quantization noise compared with read noise as a function of bit depth" width="700">
_Fig 3. More bits never take you below the read noise_

Dividing the full well of 15000 e$^-$ by $2^B$ gives the quantization step in electrons. Quantization noise equals the 3 e$^-$ read noise at 10.5 bits. Beyond that, the combined noise $\sqrt{\sigma_{read}^2 + \sigma_q^2}$ flattens toward the read noise. At 12 bits the combined noise is 3.18 e$^-$, only 6% above the read noise; at 14 bits it is 3.01 e$^-$ and the difference disappears.

Going past 12 bits on this sensor therefore doubles the data volume without adding information. Conversely, a low-noise sensor whose read noise is far below the quantization step gains the full benefit of every added bit. Whether a bit depth suffices is settled by comparison against the read noise, not by the camera specification sheet.

## 5. How many frames should be averaged?

The cheapest way to reduce noise is to capture several frames and average them. Averaging $N$ frames reduces random noise by $1/\sqrt{N}$. Does the result keep improving with longer acquisition?

<img src="/assets/img/posts/imgproc-sampling-noise-model/en/fig4-frame-averaging-fpn-floor.png" alt="Residual spread against frame count, with and without fixed-pattern noise" width="700">
_Fig 4. Frame averaging stops at the fixed-pattern floor_

Without a fixed pattern (blue), the residual spread follows the $1/\sqrt{N}$ reference exactly. Averaging 512 frames takes a single-frame spread of 17.20 DN down to 0.77 DN. With a fixed pattern (red), the curve starts at 21.06 DN, settles near 12 DN, and stops there no matter how many frames are added.

The reason lies in the nature of fixed-pattern noise. Sensitivity and dark current differ slightly from pixel to pixel, and those deviations do not change over time. The same value is added every frame, so averaging leaves it intact. The simulation used a 1% pixel-to-pixel sensitivity deviation, and at a signal of 1200 DN the floor came out at 11.99 DN — exactly 1.00%.

The practical conclusion is unambiguous. Averaging is free down to the floor, but going below it requires a different approach. Sensitivity deviations are removed by flat-field correction, dividing by an image of a uniform surface; dark-current deviations are removed by dark correction, subtracting an image taken with the light blocked. Neither reduces noise. Both identify a fixed value and subtract it. Increasing the frame count without knowing where the noise floor sits is a common waste.

## 6. The limit set by pixel spacing

So far only the magnitude and scatter of pixel values have been considered. What constraint follows from the fact that pixels are arranged in space? Replacing the lens with a better one changes nothing about it, so what exactly is lost when pixels are large?

A sensor with pixel spacing $p$ can reproduce spatial frequencies up to $f_{Nyq} = 1/(2p)$. A period must be sampled at least twice to be recognized as a period, which is the Nyquist limit. The trouble is that components above this limit do not simply disappear.

<img src="/assets/img/posts/imgproc-sampling-noise-model/en/fig5-frequency-folding.png" alt="Spatial frequencies above Nyquist folding back to lower frequencies" width="680">
_Fig 5. Frequency folding — why aliasing is dangerous in metrology_

Frequencies above Nyquist fold back into lower frequencies in a triangular pattern. A pattern at 1.2 times Nyquist, for instance, appears as one at 0.8 times. Coarse fringes that were never in the original show up in the image. This is aliasing.

<img src="/assets/img/posts/imgproc-sampling-noise-model/en/fig6-aliasing.png" alt="Aliasing in a zone plate and a natural image, before and after low-pass filtering" width="750">
_Fig 6. Content above Nyquist does not vanish — it folds down in frequency_

The top row shows a zone plate, whose frequency rises with distance from the center. In the original the rings merely grow denser outward, but the centre panel, sampled four times more coarsely, has acquired large rings and checkered patches that were never there. The bottom row applies the same treatment to a widely used test image, where coarse speckle appears across the grass and around the tripod.

The right column was low-pass filtered with a Gaussian before sampling. Components above Nyquist were erased in advance, so nothing remains to fold. Detail was lost, but no spurious pattern was created. Which is preferable — losing information or gaining falsehood — depends on the application, and in metrology it is almost always the former. Blurred detail is visible as blur, whereas a pattern produced by aliasing is indistinguishable from real signal. This is how a periodic structure absent from a surface ends up being read out as topography.

One addition: a pixel is an area rather than a point, and it integrates the light falling within it. That integration acts as a $\mathrm{sinc}(fp)$ attenuation and suppresses high frequencies to some degree. Its first zero lies at twice the Nyquist frequency, however, so roughly 64% still passes near Nyquist itself. Pixel aperture alone is no substitute for optical low-pass filtering.

## Summary and what comes next

Two numbers were fixed in this post. A pixel value is proportional to the photon count as $\mu_y = K(\eta\mu_p + \mu_d) + y_0$, and its variance grows with signal as $\sigma_y^2 = K^2\sigma_d^2 + \sigma_q^2 + K(\mu_y - y_0)$. A photon transfer curve recovers $K$ and $\sigma_d$ from images alone, provided the quantization term is removed from the intercept first. Once recovered, questions such as whether the bit depth suffices, where frame averaging stops, and how much exposure to give all become calculable. Spatially, Nyquist sets the reproduction limit, and content beyond it folds down into false signal instead of vanishing.

The next post takes up convolution. It follows the variance obtained here through a kernel, and puts a number on the trade in which reducing noise blurs edges. Widening the smoothing kernel shrinks the scatter of the edge position but widens the edge itself, and where those two curves meet sets the subpixel precision limit treated in Post 4.

## References

- EMVA, "[EMVA Standard 1288: Standard for Characterization of Image Sensors and Cameras, Release 3.0](https://www.emva.org/wp-content/uploads/EMVA1288-3.0.pdf)," European Machine Vision Association, 2010.
- J. R. Janesick, _Photon Transfer: DN → λ_, SPIE Press, 2007 (photon transfer curve and extraction of sensor constants).
- R. C. Gonzalez, R. E. Woods, _Digital Image Processing_, 4th ed., Pearson, 2018, ch. 2 (sampling and quantization), ch. 4 (aliasing).
- S. van der Walt et al., "[scikit-image: image processing in Python](https://doi.org/10.7717/peerj.453)," _PeerJ_, 2, e453, 2014 (the standard test image in Fig 6).
- Y. Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," Ph.D. dissertation, Seoul National University, 2025, sec. 4.2 (locating the back focal plane by measuring focal spot diameter with image processing).

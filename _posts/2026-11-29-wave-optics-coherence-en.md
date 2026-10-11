---
layout: post
date: 2026-11-29 20:00:00 +0900
title: "Wave Optics Foundations 2 — When Does Light Interfere? Temporal and Spatial Coherence"
lang: en
lang-exclusive: ["en"]
permalink: /posts/wave-optics-coherence/
page_id: wave-optics-coherence
categories: [Optics, Wave Optics]
tags: [wave-optics, coherence, interferometry, fourier-transform-spectroscopy, speckle]
description: "Coherence measures how long and over how wide a region the phase of light stays orderly. This post follows how the spectrum and the source size set it, and how it is used in measurement."
math: true
---

Shine a laser pointer on a white wall and the spot is not evenly bright. Up close it is a dense tangle of bright and dark grains, and the grains wriggle as you move your head (Figure 1). Light the same wall with a flashlight and no grains appear. The wall is the same; only the light has changed.

<img src="/assets/img/posts/wave-optics-coherence/ext-laser-speckle.jpg" alt="Speckle pattern from laser light scattered by a rough surface" width="420">
_Figure 1. Speckle pattern on a wall from laser light scattered off a plastic surface. (Source: Epzcaw, [Objective speckle](https://commons.wikimedia.org/wiki/File:Objective_speckle.jpg), public domain)_

[Part 1](/en/posts/wave-optics-interference/) left two puzzles. In Young's experiment with sunlight, only a few central fringes were visible before they faded. Two lamps in a room produced no fringes at all. The grains of Figure 1 are a third puzzle. Laser light interferes with itself even after scattering from a rough wall and forms a pattern; flashlight light does not.

All three have one answer. Part 1 treated light as a perfect sine wave, a wave whose phase stays orderly forever and everywhere. Real light is not like that. The time and the region over which its phase stays orderly are limited, and those limits decide when interference can be seen. The degree of this order is called coherence [1, 2, 3]. This post splits coherence into its temporal and spatial parts, follows how the spectrum and the size of the source set each, and shows how short coherence turns out to be useful in measurement.

## 1. How long is the phase remembered?

In an ordinary source, light is emitted by many atoms independently. A single atom emits for something like $10^{-8}$ to $10^{-9}$ s, and thermal motion and collisions smear out the frequency [1]. Watched at one point, the light is therefore a string of short sine-wave pieces joined with random phases [1, 3]. Within a piece the phase is predictable; when the piece changes, the phase is reset.

This "memory of phase" is measured by the autocorrelation: how much the field at time $t$ and the field a time $\tau$ later move together, on average. Normalized by the intensity, it is the complex degree of coherence $g(\tau)$ [3]:

$$ g(\tau) = \frac{\langle U^*(t)\,U(t+\tau)\rangle}{\langle \lvert U(t)\rvert^2\rangle} $$

At $\tau = 0$, $\lvert g\rvert = 1$; as $\tau$ grows the phases at the two times become unrelated and $\lvert g\rvert$ falls to zero [3]. The width of the delay over which correlation survives is the coherence time $\tau_c$, and the distance light travels in that time, $l_c = c\tau_c$, is the coherence length [1, 3]. The coherence length is the length over which the wave can be trusted to be a sinusoid [1].

<img src="/assets/img/posts/wave-optics-coherence/en/fig2-random-light-coherence-time.png" alt="Field, degree of coherence and spectrum of two random light waves with different spectral widths" width="780">
_Figure 2. Two random light waves that differ only in spectral width. (a), (d) Field (colour) and envelope (black). (b), (e) Degree of coherence versus delay. Lines are measured directly on the generated waves; circles are computed by Fourier transforming the spectrum. (c), (f) Power spectra. A spectrum ten times broader gives a coherence time ten times shorter._

Figure 2 makes this visible. The top row is light with a narrow spectrum, the bottom row light with a spectrum ten times broader. The narrow-band light has a slowly undulating envelope and keeps its sinusoidal shape for tens of periods. The broad-band light changes amplitude and phase every few periods. The widths of the two degrees of coherence differ by the same factor of ten.

## 2. Spectrum and coherence are a Fourier pair

The degree of coherence in Figure 2 was obtained in two ways: by measuring the autocorrelation of the generated random wave (lines) and by Fourier transforming the spectrum (circles). The two coincide, and not by accident. The autocorrelation of light and its power spectrum are Fourier transforms of each other; this is the Wiener–Khinchin theorem [1, 3]. Since the coherence time and the spectral width are the widths of a Fourier pair, they are inversely proportional [3]:

$$ \tau_c \approx \frac{1}{\Delta\nu}, \qquad l_c = c\,\tau_c \approx \frac{\lambda^2}{\Delta\lambda} $$

[2, 3]. The exact constant ranges from about 0.3 to 1 depending on the spectral shape and on how the widths are defined [3]. Defining the coherence time as $\int\lvert g\rvert^2 d\tau$, for example, gives $\tau_c = 0.66/\Delta\nu$ for a Gaussian spectrum of full width $\Delta\nu$ and $0.32/\Delta\nu$ for a Lorentzian [3].

The intuition is this. Sine waves of slightly different frequencies all start in phase but drift apart with time. If their frequencies differ by $\Delta\nu$, after roughly $1/\Delta\nu$ they have slipped by a full cycle, and their sum no longer resembles its starting shape. The wider the spread of frequencies, the sooner this happens.

This relation gives the coherence length of each source [3].

| Source | Coherence time | Coherence length |
|---|---|---|
| Filtered sunlight (0.4–0.8 µm) | 2.7 fs | 0.8 µm |
| LED (1 µm, 50 nm width) | 67 fs | 20 µm |
| Low-pressure sodium lamp | 2 ps | 0.6 mm |
| Multimode helium–neon laser | 0.67 ns | 20 cm |
| Single-mode helium–neon laser | 1 µs | 300 m |

The coherence length of sunlight is 0.8 µm, comparable to its wavelength. That is why the sunlight fringes in Figure 4 of Part 1 faded after a few fringes. Once the path difference exceeds the coherence length the fringes wash out, and the number of visible fringes is roughly $l_c/\lambda = \nu_0/\Delta\nu$ [3]. A single-mode laser, by contrast, has a coherence length of hundreds of metres; in a laboratory-sized interferometer it is effectively a perfect sine wave.

Increasing the coherence time is simple: trim the spectrum with a narrow-band filter. The price is intensity. Every frequency thrown away takes its light with it [3].

## 3. Measuring coherence with a Michelson interferometer

The Michelson interferometer of Part 1 splits light in two, delays one half, and recombines them. It overlaps light with a delayed copy of itself, which is exactly the definition of the degree of coherence. So when the mirror is moved and the fringe visibility is measured, that visibility is $\lvert g\rvert$ at the corresponding path difference [3]. In Part 1 the visibility was set by the intensity ratio of the two beams; for partially coherent light it is multiplied by $\lvert g\rvert$ [3]:

$$ V = \frac{2\sqrt{I_1 I_2}}{I_1 + I_2}\,\lvert g\rvert $$

<img src="/assets/img/posts/wave-optics-coherence/en/fig3-michelson-interferograms.png" alt="Michelson interferograms of an LED, the sodium doublet and white light, and the spectrum recovered by Fourier transform" width="780">
_Figure 3. Computed output of a Michelson interferometer versus path difference. (a) LED: fringes only within about ±10 µm of path difference. (b) Sodium D doublet: the fringes are so dense they appear as a filled band, and its edge (black) — the visibility — drops to zero and returns every 0.58 mm. (c) White light: only a few central fringes. (d) Spectrum recovered by Fourier transforming the interferogram of (b). Scanning to 1.8 mm separates the 589.0 nm and 589.6 nm lines; scanning only to 0.3 mm merges them into one._

In Figure 3(a) the LED fringes exist only within about ±10 µm of path difference. The white light of (c) shows only two or three central fringes. The distance the mirror moves before the fringes vanish is the coherence length.

The sodium lamp of (b) is more interesting. The yellow light of sodium consists of two lines, at 589.0 nm and 589.6 nm. Each is nearly a pure sine wave, but their wavelengths differ slightly, so as the path difference grows their fringe patterns slide apart. At a path difference of $\lambda^2/(2\Delta\lambda) = 0.29$ mm the bright fringes of one line fall on the dark fringes of the other and the pattern disappears; at 0.58 mm they line up again and the fringes return. It is the same as the beating of two musical notes with slightly different pitch. Since the mirror moves half the path difference, the fringes vanish and reappear every 0.29 mm of mirror travel. This is how Michelson read the structure of spectral lines from fringe visibility [2].

### From interferogram to spectrum

Read the Wiener–Khinchin theorem backwards and the spectrum can be obtained by Fourier transforming the interferogram recorded while moving the mirror [3]. This is Fourier transform spectroscopy. Figure 3(d) transforms the sodium interferogram of (b). Recording out to 1.8 mm of path difference separates the two lines cleanly; recording only to 0.3 mm merges them. The smallest resolvable wavenumber difference is roughly the inverse of the scanned path difference. When 1/(0.3 mm) = 33 cm⁻¹ exceeds the separation of the two lines (17 cm⁻¹), they are not resolved. Scanning further increases the resolution.

The Fourier transform infrared spectrometer (FTIR), widely used in infrared spectroscopy, works on this principle. Bell summarizes the advantages of interferometric spectroscopy in two points [4]. A grating spectrometer measures only the narrow band of wavelengths passing its exit slit at any moment, whereas an interferometer measures all wavelengths at every moment. This is the Fellgett, or multiplex, advantage. An interferometer can also accept light through a large circular aperture rather than a narrow slit, collecting more light at the same resolution. This is the Jacquinot, or throughput (étendue), advantage [4]. The étendue of [Geometrical Optics Foundations 3](/en/posts/geometric-optics-stops-pupils-etendue/) appears again here.

## 4. Short coherence is useful

Short coherence looks like a defect, but in measurement it turns out to be useful.

### Film thickness and coherence

The oil and soap films of Part 1 showed colours only when thin. A thick window pane reflects light from both its faces, yet shows no coloured bands. Once the path difference $2nt$ between the two reflections exceeds the coherence length, the two beams no longer interfere and their intensities simply add. The coherence length of sunlight as seen by the eye is under 1 µm, so films thicker than a few micrometres lose their fringes.

<img src="/assets/img/posts/wave-optics-coherence/en/fig4-thin-film-coherence.png" alt="Reflectance of a water film seen by the eye and by a spectrometer, and reflection spectra of thicker films" width="780">
_Figure 4. (a) Reflectance of a water film in air versus thickness, viewed through a broad band like the eye (150 nm FWHM, red) and a narrow band like one spectrometer pixel (2 nm FWHM, blue). For the eye the fringes vanish beyond about 1 µm; for the spectrometer they remain at 30 µm. (b) Reflection spectra of films 5 µm and 30 µm thick. The thicker the film, the denser the fringes._

Look at the same film with a spectrometer, however, and the fringes return (Figure 4). Each spectrometer pixel receives only a narrow band of wavelengths, so the light it sees has a long coherence length. A 2 nm band has a coherence length of about 0.15 mm, enough for the two reflections from a film tens of micrometres thick to interfere. The fringes now appear along the wavelength axis rather than the thickness axis (Figure 4(b)). The thicker the film, the denser the fringes in wavelength, and their spacing gives the thickness. This is how the spectral reflectometers that measure film thickness in semiconductor processing work, as do [ellipsometers](/en/posts/ellipsometry-thin-film-multilayer-reflectance/). Conversely, a very thick layer such as a substrate washes out the fringes even within one spectrometer pixel, so multilayer calculations add the two reflections from such a layer by intensity, not by amplitude [5].

### White-light interferometry and OCT

The white-light fringes of Figure 3(c) exist only near zero path difference. Turned around, this locates "zero path difference" very sharply. A white-light interferometer uses the sample surface as one arm of a Michelson and moves the reference mirror; for each pixel it finds the mirror position at which the fringes are strongest and takes that as the height of the pixel. A laser interferometer repeats its fringes every half wavelength and cannot tell apart steps that differ by more than half a wavelength, but white light has only one fringe packet and no such ambiguity.

Applying the same principle inside tissue gives optical coherence tomography (OCT) [3]. As the reference mirror moves, the reference light interferes with light returning from boundaries inside the sample, producing a fringe packet about one coherence length wide for each boundary (Figure 5) [3]. The position of the packet gives the depth of the boundary, and its size the strength of the reflection. Since the depth resolution is set by the coherence length, a broader spectrum resolves depth more finely. For a Gaussian spectrum the full width at half maximum of the fringe packet is

$$ \Delta z = \frac{2\ln 2}{\pi}\,\frac{\lambda_0^2}{\Delta\lambda} $$

For a source centred at 840 nm with a 50 nm FWHM, this is 6.2 µm, which is why two boundaries 8 µm apart are resolved in Figure 5. OCT is widely used in clinical medicine, and frequency-domain schemes, which sweep the wavelength of the source or measure the spectrum with a spectrometer instead of moving the reference mirror, offer better sensitivity and speed [3]. It is the same idea as the Fourier transform spectroscopy of Section 3.

<img src="/assets/img/posts/wave-optics-coherence/en/fig5-oct-ascan.png" alt="OCT interference signal and envelope" width="780">
_Figure 5. Computed time-domain OCT (840 nm, 50 nm FWHM) of a sample with boundaries at depths of 100, 130, 138 and 200 µm (dotted). (a) Interference signal as the reference mirror moves. (b) Its envelope. Each boundary gives a peak 6.2 µm wide, and the two boundaries 8 µm apart are resolved._

## 5. Spatial coherence: the size of the source

So far we have asked how orderly the phase stays in time at one point. The same question can be asked in space: at the same instant, how closely do the fields at two laterally separated points move together? Young's experiment measures exactly this. For the light from two holes to interfere, the light at the two holes must be coherent with each other. This is spatial coherence [1, 3].

In Part 1, Young placed one more pinhole in front of the two holes, and spatial coherence is the reason. An extended source is a collection of many mutually independent point sources [2]. Each point source makes its own perfect Young's pattern, but a point at a different position shifts its pattern sideways. Patterns from independent points add in intensity, so the shifted patterns overlap and the contrast blurs [1, 3]. The larger the source and the farther apart the holes, the more the patterns shift and the more they blur (Figure 6(a)). The pinhole shrinks the source to prevent this blurring. A pinhole imparts spatial coherence but not temporal coherence [3].

<img src="/assets/img/posts/wave-optics-coherence/en/fig6-spatial-coherence.png" alt="Young's fringes for different source sizes and the spatial degree of coherence of sunlight" width="780">
_Figure 6. (a) Fringes on a screen 1 m behind two holes 50 µm apart lit at 550 nm by a uniform disc source of various angular diameters. (b) Visibility versus hole separation for sunlight (angular diameter 0.53°). Circles are measured by adding in intensity the patterns of the points of the source; the dashed line is the result of the van Cittert–Zernike theorem. The visibility first reaches zero at a separation of 73 µm._

Figure 6(b) shows the visibility measured with sunlight as the hole separation varies. For a uniformly bright disc source the visibility is

$$ \lvert g\rvert = \left\lvert \frac{2J_1(x)}{x} \right\rvert, \qquad x = \frac{\pi d\,\theta_s}{\lambda} $$

and it first reaches zero at a separation $d = 1.22\,\lambda/\theta_s$ [1, 3], where $\theta_s$ is the angular diameter of the source and $J_1$ the first-order Bessel function. For sunlight at 550 nm this separation is 73 µm. To do Young's experiment with bare sunlight, the two holes must be within tens of micrometres of each other [3]. Light coming directly from an extended surface, such as just in front of a hot filament, has a coherence area of the order of the wavelength squared and is effectively incoherent [3].

The shape of this formula will return in the next part. The amplitude of the diffraction pattern of a circular aperture — the Airy pattern — is $2J_1(x)/x$. This is no coincidence. As van Cittert showed in 1934 and Zernike in 1938 [2], the spatial degree of coherence produced far away by a source made of independent points is the Fourier transform of the source's intensity distribution. Mathematically it is the same as the diffraction amplitude of coherent light through an aperture of the same shape. This is the van Cittert–Zernike theorem. Even an incoherent source gains coherence as its light spreads out, which is why light from very distant stars is coherent over a wide region on the ground [3].

### Measuring the size of a star

Turned around, this relation measures the size of a source. In 1890, following a suggestion by Fizeau, Michelson proposed a stellar interferometer: two mirrors spread apart in front of a telescope, with their separation varied until the fringes vanish [1]. At the separation $h$ where the fringes first vanish, the angular diameter of the star is $\theta_s = 1.22\lambda/h$. In December 1920, with the instrument mounted on the 2.5 m telescope at Mount Wilson, the fringes of Betelgeuse vanished at a mirror separation of 121 inches (about 3.1 m). With a wavelength of 570 nm this gives an angular diameter of 0.047 arcseconds, and since the distance was known, a diameter about 280 times that of the Sun [1]. The telescope alone could not resolve this angle; measuring coherence did.

<img src="/assets/img/posts/wave-optics-coherence/ext-alma-antennas.jpg" alt="ALMA radio telescope antennas on the Chajnantor plateau in Chile" width="680">
_Figure 7. ALMA's 12 m antennas on the Chajnantor plateau in Chile (about 5000 m altitude). The coherence of the signals received by each pair of antennas is measured, and the brightness distribution of the sky is reconstructed. (Source: Iztok Bončina/ESO, [ALMA Antennas on Chajnantor](https://commons.wikimedia.org/wiki/File:ALMA_Antennas_on_Chajnantor.jpg), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0))_

Radio astronomy developed this principle into an imaging technique [1]. A radio interferometer spreads many antennas over a wide area (Figure 7) and measures the coherence of the signals received by each pair. By the van Cittert–Zernike theorem, each such value is one point of the Fourier transform of the sky's brightness distribution. The more varied the separations and orientations of the antenna pairs, the more of the Fourier plane is filled, and inverting the transform gives an image.

## 6. Speckle: noise made by coherence

Back to Figure 1. The wall is rough on the scale of a wavelength. When laser light scatters from it, light from many small patches of the wall reaches each point of the eye along different optical paths. Because the laser's coherence is long, all of this light interferes. The paths are random, so the phases are random: in the phasor picture of Part 1, many arrows of similar length pointing in random directions, joined head to tail. In some places the sum is long; in others it is nearly zero. This produces the bright and dark grains called speckle.

<img src="/assets/img/posts/wave-optics-coherence/en/fig8-speckle.png" alt="Computed speckle pattern, average of several patterns, and intensity distribution" width="780">
_Figure 8. Computed image, through a lens with a circular pupil, of coherent light scattered by a rough surface. (a) A single speckle pattern: the contrast (standard deviation over mean) is 1. (b) Average of 16 independent speckle patterns: the contrast drops to 1/4. (c) Intensity distribution. A single speckle pattern follows an exponential distribution (dashed), with darkness the most common outcome._

A sum of random arrows has a statistically fixed form. The intensity follows an exponential distribution, the most probable intensity is zero, and the standard deviation equals the mean, so the contrast is 1 (Figure 8(a), (c)). In imaging terms, signal and noise are of the same size. Flashlight light shows no speckle because its coherence is short (in time) and its source is large (in space). Components of slightly different wavelength, and components arriving from different directions, each make a different speckle pattern, and these add in intensity and average out. Averaging $N$ independent speckle patterns reduces the contrast to $1/\sqrt{N}$ (Figure 8(b)).

That is why instruments that image with lasers often reduce coherence deliberately: they rotate a diffuser, or illuminate with several wavelengths or angles, and average independent speckle patterns. Conversely, some measurements use the motion of speckle itself as the signal. Since the speckle pattern changes when a surface moves even slightly, comparing speckle at two times measures deformation or vibration. Coherence is a measuring tool in the form of interference fringes on one side, and noise to be removed on the other.

## Summary

- Real light is a string of short sine-wave pieces joined with random phases. The time over which phase is remembered is the coherence time $\tau_c$, and the corresponding distance $l_c = c\tau_c$ is the coherence length.
- The degree of coherence and the spectrum are a Fourier pair (the Wiener–Khinchin theorem), so $l_c \approx \lambda^2/\Delta\lambda$: under 1 µm for sunlight, hundreds of metres for a single-mode laser.
- The fringe visibility of a Michelson interferometer is the degree of coherence. The sodium doublet beats every 0.29 mm of mirror travel, and Fourier transforming the interferogram gives the spectrum (Fourier transform spectroscopy).
- Films thicker than the coherence length lose their fringes, and narrowing the band with a spectrometer brings the fringes back along the wavelength axis. Short coherence becomes the depth resolution of white-light interferometry and OCT.
- An extended source reduces spatial coherence. The degree of coherence is the Fourier transform of the source's brightness distribution (van Cittert–Zernike), and for a disc source it first reaches zero at $1.22\lambda/\theta_s$. This measures the sizes of stars and builds images in radio interferometers.
- Coherent light scattered by a rough surface forms speckle with contrast 1. Averaging $N$ independent speckle patterns reduces the contrast to $1/\sqrt{N}$.

Section 5 showed that the spatial degree of coherence has the same shape as the diffraction pattern of a circular aperture. Why does light passing through an aperture draw such a pattern? Why a lens images a star as a small disc rather than a point, and why the diffraction images of Geometrical Optics Foundations 4 and the degree of coherence of this post come from the same mathematics, is the subject of the next part: diffraction.

## References

1. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Section 7.4.3 — coherence time and coherence length; Section 9.2.1 — temporal and spatial coherence; Section 11.3.5 — correlation and the Wiener–Khinchin theorem; Chapter 12 — visibility and the mutual coherence function, the Michelson stellar interferometer and the measurement of Betelgeuse, radio interferometry.
2. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Section 7.5.8 — width of spectral lines and fringe visibility; Chapter 10 — interference and diffraction with partially coherent light, history of coherence theory (Verdet, Michelson, van Cittert, Zernike).
3. B. E. A. Saleh and M. C. Teich, *Fundamentals of Photonics*, 3rd ed. (Wiley, 2019), Sections 12.1–12.3 — statistics of random light, the degree of coherence and coherence time, the Wiener–Khinchin theorem, coherence lengths of sources (Table 12.1-2), interference of partially coherent light, Fourier transform spectroscopy, OCT, Young's experiment and spatial coherence, source size and coherence distance.
4. R. J. Bell, *Introductory Fourier Transform Spectroscopy* (Academic Press, 1972), Chapter 1 — the Fellgett (multiplex) and Jacquinot (throughput) advantages; Chapter 4 — interferogram and spectrum of the sodium doublet.
5. C. C. Katsidis and D. I. Siapkas, "General transfer-matrix method for optical multilayer systems with coherent, partially coherent, and incoherent interference," *Applied Optics* 41, 3978–3987 (2002). [doi:10.1364/AO.41.003978](https://doi.org/10.1364/AO.41.003978)

---
title: "Electron Microscopy Foundations 6 — Phase Contrast and the CTF: Is That Bright Spot an Atom?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-tem-phase-contrast-ctf/
date: 2026-11-01 20:00:00 +0900
page_id: electron-microscopy-tem-phase-contrast-ctf
categories: [Electron Microscopy, TEM Imaging]
tags: [electron-microscopy, tem, hrtem, contrast-transfer-function, scherzer-defocus]
description: A bright spot in a high-resolution image is not guaranteed to be an atom. Change the defocus alone and the contrast inverts.
math: true
---

[Post 5](/en/posts/electron-microscopy-tem-diffraction-contrast/) formed images by selecting either the transmitted beam or a diffracted beam with the objective aperture. Selecting means discarding, and whatever information the discarded beams carried is lost with them.

This post discards nothing. The aperture opens wide enough to pass the transmitted beam together with several diffracted beams, and they are allowed to interfere in the image plane. High-resolution TEM (HRTEM) images, in which the crystal lattice periodicity appears directly, form this way.

There is a trap in it. Seeing a lattice pattern does not mean the bright spots sit where the atoms are. Photograph the same specimen under the same conditions, shift the focus slightly, and the contrast inverts completely. Why that happens, and how far such an image can be trusted, is the subject here.

## 1. A specimen changes phase, not amplitude

What happens to an electron wave that crosses a thin specimen?

If it is thin enough, the electron is neither absorbed nor scattered away through a large angle. Only the phase changes. Where an atom sits the potential runs deeper, the electron wavelength shortens slightly, and the phase advances by that much. The wave just past the specimen can therefore be written

$$ \psi_{exit}(x,y) = \exp\!\left[\,i\,\sigma_e V_p(x,y)\right] $$

with $V_p$ the potential projected along the thickness and $\sigma_e$ the interaction constant. An object that alters only phase is a phase object, and for $\sigma_e V_p \ll 1$ it expands under the weak phase object approximation (WPOA).

$$ \psi_{exit} \approx 1 + i\,\sigma_e V_p $$

The problem surfaces here. A detector measures intensity, never phase. Computing the intensity of that wave gives

$$ \lvert \psi_{exit} \rvert^2 = 1 + \sigma_e^2 V_p^2 \approx 1 $$

and the first-order term in $V_p$ has vanished. **A perfect lens forming a perfect image would show nothing at all.** Where the atoms are is recorded in the phase, and measuring intensity strips that record away.

So what converts phase into intensity? The answer is unexpected: the aberrations of the lens.

## 2. The aberration function and the contrast transfer function

In post 2 spherical aberration was the chief obstacle to resolution. Its role reverses here.

Lens imperfections add a phase to the wave that depends on spatial frequency $k$. That added phase is the aberration function $\chi(k)$. The diffraction vector $g$ of post 5 is one particular value of $k$, fixed by the crystal; the argument here holds with or without a crystal, so $k$ runs continuously.

$$ \chi(k) = \pi\,\Delta f\,\lambda\,k^2 + \frac{1}{2}\,\pi\,C_s\,\lambda^3 k^4 $$

Defocus $\Delta f$ enters as $k^2$ and spherical aberration $C_s$ as $k^4$. By convention $\Delta f$ is negative at underfocus. That the two terms can carry opposite signs proves decisive shortly.

Applying this phase to a weak phase object and recomputing the intensity brings the vanished first-order term back.

$$ I(x,y) \approx 1 + 2\,\sigma_e V_p \otimes \mathcal{F}^{-1}\!\left[\sin\chi(k)\right] $$

Here $\sin\chi(k)$ is the contrast transfer function (CTF). It states, frequency by frequency, how much phase information reaches the intensity and with what sign.

Where $\chi$ is zero, $\sin\chi$ is zero as well. A perfect lens with no aberration shows no phase object, which is section 1's conclusion restated. Phase contrast is obtained not by removing aberration but by leaving the right amount of it.

<img src="/assets/img/posts/electron-microscopy-tem-phase-contrast-ctf/en/fig1-ctf-curves.png" alt="Contrast transfer function curves compared across defocus values" width="820">
_Fig 1. $\sin\chi(k)$ for three defocus values_

The three curves differ only in defocus. The solid red one is the Scherzer defocus explained below; the others are 20 nm under and 20 nm over. Same lens, same specimen, yet shifting the focus alone completely changes which spatial frequencies get through, how strongly, and with which sign.

Note especially how often the curves cross zero. Information at a frequency where $\sin\chi = 0$ never reaches the image at all, and information in a band where the sign has flipped arrives with its contrast reversed.

## 3. Scherzer defocus — cancelling aberration with aberration

If the curve oscillates like that, when does a usable image appear?

When $\sin\chi$ stays flat near $-1$ over a wide band. Every spatial frequency in such a band transfers with the same sign and comparable strength, so the image faithfully resembles the object. The opposite signs available to the two terms of $\chi$ serve here: taking $\Delta f$ negative at underfocus makes the $k^2$ term negative, cancelling the positive $k^4$ term across a broad range.

The defocus that cancels best is the Scherzer defocus.

$$ \Delta f_{Sch} = -1.2\sqrt{C_s\lambda} $$

Substituting the $C_s = 0.5$ mm and $\lambda = 2.51$ pm of post 2 gives $-42.5$ nm. The solid red curve of Figure 1 is that condition, and $\sin\chi$ holds negative out to $k = 5.2$ nm$^{-1}$. This band is the Scherzer passband, and its edge is the point resolution.

$$ \frac{1}{5.2\ \mathrm{nm^{-1}}} = 0.192\ \mathrm{nm} $$

Post 2 obtained 0.197 nm from the Scherzer formula $0.66(C_s\lambda^3)^{1/4}$. The two differ by about 2 %, since the formula approximates where the passband ends, but they are the same quantity computed two ways. The 0.2 nm quoted in post 1 for uncorrected HRTEM is ultimately this number.

## 4. Damping envelopes and the information limit

Does information beyond the point resolution disappear entirely? No — it merely becomes unreadable.

Mathematically $\sin\chi$ never stops oscillating however large $k$ grows. In practice two effects erase the high frequencies.

One is temporal coherence. The gun energy spread of post 1 and instabilities in the lens supply make the focal position spread slightly. Each focus within that spread carries a different $\chi$, and the contributions cancel, producing an envelope that falls exponentially in $k^4$.

The other is spatial coherence. The illuminating beam is not a perfect plane wave but carries some convergence, so components arriving at slightly different angles experience different $\chi$.

<img src="/assets/img/posts/electron-microscopy-tem-phase-contrast-ctf/en/fig2-envelopes.png" alt="Damping envelopes and the relation between point resolution and information limit" width="820">
_Fig 2. The effective transfer function with damping envelopes applied_

The thin grey line is the undamped $\sin\chi$, oscillating at full amplitude to the right edge. The blue dashed and green dotted curves are the two envelopes, and the heavy red curve is what actually transfers. The frequency at which the envelope falls to $1/e^2$ is the information limit, which under these conditions is $k = 7.1$ nm$^{-1}$, or 0.140 nm.

Between the point resolution at 0.192 nm and the information limit at 0.140 nm lies the shaded band. Information there does not disappear. Because $\sin\chi$ changes sign repeatedly across it, however, some frequencies arrive faithfully and others inverted, mixed together. No amount of looking at the image distinguishes them.

Using that band requires computation. Reconstructing the original exit wave from a series of images taken at stepped defocus — focal series reconstruction — does it, as does the electron ptychography of post 7. It is territory read by calculation rather than by eye.

## 5. Is that bright spot an atom?

The opening question can now be answered by applying the preceding calculation to an actual lattice.

A phase object was built by placing Gaussian projected potentials on a square lattice of 0.25 nm spacing, and images were synthesized by passing it through the transfer function above. The exit wave was propagated directly rather than expanded under the weak phase object approximation.

```python
# core of generate_figures.py (full code: _code/electron-microscopy-tem-phase-contrast-ctf/)
psi = np.exp(1j * v)                       # phase object: it alters only phase
transfer = np.exp(1j * chi(kr, df)) * envelope_temporal(kr) \
    * envelope_spatial(kr, df)
img = np.abs(np.fft.ifft2(np.fft.fft2(psi) * transfer)) ** 2
```

<img src="/assets/img/posts/electron-microscopy-tem-phase-contrast-ctf/en/fig3-lattice-images.png" alt="Three lattice images synthesized at different defocus compared with the true atom positions" width="900">
_Fig 3. The same phase object imaged at three defocus values. All three images share one intensity scale_

Leftmost is the truth — the projected potential, showing where the atoms are. The other three differ only in defocus.

At the Scherzer defocus of $-42.5$ nm, $\sin\chi = -0.82$ at the lattice frequency of 4 nm$^{-1}$, and the atom sites appear **dark**, exactly the sign the weak phase object approximation predicts.

At $-50.1$ nm, $\sin\chi$ passes through zero. The lattice pattern effectively disappears, its amplitude falling by a factor of thirty. The specimen is unchanged, yet the lattice is gone from the image.

At $-62.5$ nm, $\sin\chi$ reaches $+1$ and the atom sites appear **bright**. The contrast has inverted completely.

All three came from the same specimen, the same lens and the same voltage, with only the focus changed, and none of them is a faulty exposure. Reading bright spots in an HRTEM image as atoms is therefore unsafe. Interpretation requires knowing the defocus and the thickness, which in practice means matching against simulated images computed across a range of conditions.

A thicker specimen makes matters worse. The phase object approximation itself breaks down as electrons scatter repeatedly inside the specimen, and then contrast inverts with thickness as well. Hence the need to match thickness and defocus together when interpreting HRTEM images.

## 6. What aberration correctors changed

What happens to the CTF when $C_s$ shrinks?

Post 2 noted that correctors bring $C_s$ down to the micrometre range. The $k^4$ term of $\chi$ shrinks accordingly, widening the passband and improving the point resolution. Push the first zero past the information limit and point resolution and information limit coincide, at which moment the ambiguous shaded band disappears. The range that can be read by eye extends all the way to the information limit.

Something interesting follows. Drive $C_s$ to zero and $\chi$ approaches zero too, returning to the problem of section 1: phase contrast vanishes. Corrected instruments therefore sometimes tune $C_s$ to a small negative value and work at overfocus. Atoms then appear bright against a dark background, which is often easier to interpret. This is negative $C_s$ imaging.

The goal has shifted from eliminating aberration to leaving exactly as much as wanted, with the sign chosen deliberately. Scherzer's theorem in post 2 forbade negative aberration in a rotationally symmetric lens; once correctors circumvented that prohibition, aberration became a design variable.

## Summary and what comes next

A thin specimen alters the phase of the electron wave rather than its amplitude, so a perfectly formed image carries nothing in its intensity. What converts phase into intensity is lens aberration, and the contrast transfer function $\sin\chi(k)$ describes that conversion frequency by frequency.

The defocus and spherical aberration terms of $\chi$ can oppose one another. At the Scherzer defocus of $-42.5$ nm the cancellation opens a passband out to $k = 5.2$ nm$^{-1}$, whose edge is the point resolution of 0.192 nm — the same quantity as post 2's 0.197 nm. Beyond it information still transfers, but $\sin\chi$ mixes signs so that it cannot be read without computation, and the coherence envelopes cut the information off entirely at 0.140 nm.

Figure 3 carries the central conclusion. Photograph one specimen at three focus settings and the atoms come out dark, invisible, or bright. A high-resolution image is not a photograph of an atomic arrangement but the output of a transfer function, and interpreting it requires knowing the conditions that produced it.

The next post closes the series. Illumination has been a broad beam until now; post 7 converges it back to a point and scans. It covers why reciprocity holds in scanning transmission electron microscopy, and how a high-angle annular dark field detector produces contrast proportional to atomic number, free of the sign-reversal problem met here. EDS and EELS for reading composition follow, along with 4D-STEM, which stores whole diffraction patterns and recovers phase by computation, and finally FIB specimen preparation and semiconductor metrology applications.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 27–31 (phase contrast, the CTF, image simulation).
- J. C. H. Spence, *High-Resolution Electron Microscopy*, 4th ed., Oxford University Press, 2013 (standard treatment of the weak phase object approximation and transfer theory).
- O. Scherzer, "The theoretical resolution limit of the electron microscope," *Journal of Applied Physics*, vol. 20, pp. 20–29, 1949 (derivation of the optimum defocus).
- C. L. Jia, M. Lentzen, and K. Urban, "Atomic-resolution imaging of oxygen in perovskite ceramics," *Science*, vol. 299, pp. 870–873, 2003 (negative $C_s$ imaging).
- [exCTF simulator, *Journal of Analytical Science and Technology*](https://link.springer.com/article/10.1186/s40543-020-00231-9) (CTF before and after aberration correction).

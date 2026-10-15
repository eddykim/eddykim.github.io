---
title: "Electron Microscopy Foundations 1 — If the Electron's Wavelength Is Far Smaller Than an Atom, Why Can We Barely See One?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-why-electrons/
page_id: electron-microscopy-why-electrons
date: 2026-09-28 20:00:00 +0900
categories: [Electron Microscopy, Electron Optics]
tags: [electron-microscopy, de-broglie, resolution, electron-gun, vacuum]
description: A 2.5 pm electron still resolves only 0.2 nm. This post traces the gap to the lens, then explains why electron guns and vacuum matter.
math: true
---

Anyone who works with optical metrology long enough keeps running into the same wall: wavelength. Measuring a film far thinner than the wavelength is possible, as [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/) showed, but converting thickness into a phase and reading it back is not the same as seeing the layer. For resolution, the ability to tell two neighbouring structures apart as two, visible light stops near 200 nm.

So one looks for a shorter wavelength, and the electron presents itself. Treated as a wave, an electron has a wavelength five orders of magnitude shorter than visible light. On paper it should pick out individual atoms with room to spare. Yet a real transmission electron microscope (TEM) delivers about 0.2 nm, and even an aberration-corrected instrument reaches only 0.05 nm. That is dozens of times worse than the wavelength.

Electron microscopes fall into three broad kinds, depending on how they use the electrons. A TEM floods a very thin specimen with electrons and forms an image from those that pass through. A scanning electron microscope (SEM) sweeps a fine electron probe across the surface of a thick specimen, line by line, and counts the electrons thrown back at each point. A scanning transmission electron microscope (STEM) combines the two: it scans a probe over a thin specimen and counts the electrons that come out the other side.

This series starts from the gap between wavelength and resolution. The history of electron microscopy is less a story about obtaining a short wavelength than a record of fighting everything that prevents that wavelength from becoming resolution. This post measures the gap and sketches the culprit. Its identity, and what can be done about it, belong to post 2 onward. Put in a reader's own words, the questions the series sets out to answer are these, with the place each is taken up in brackets.

- How short is an electron's wavelength, and what sets it? (sections 2 and 3)
- If the wavelength is eighty times smaller than the spacing between atoms, why can we barely see atoms? (sections 4 and 5)
- Why are there several kinds of electron gun, and how does one choose? (section 6)
- Why must the specimen sit in a vacuum? (section 7)
- Why can a lens made of magnetic field never be concave? ([post 2](/en/posts/electron-microscopy-electron-optics/))
- If the beam is 1 nm wide, can you see 1 nm detail? ([post 3](/en/posts/electron-microscopy-beam-specimen-interaction/))
- Why do electron micrographs look three-dimensional in black and white? (post 4)
- Why do defects in a crystal show up as dark lines? (post 5)
- Is a bright spot in an atomic-resolution image really an atom? (post 6)
- Can an electron microscope tell which atom is which? (post 7)

## 1. Why does the optical microscope stop at 200 nm?

Fix the baseline first. Where does the claim that visible light stops at 200 nm come from?

Abbe's limit on resolution reads

$$ d = \frac{\lambda}{2\,\mathrm{NA}} $$

where $\lambda$ is the wavelength and $\mathrm{NA} = n\sin\alpha$ is the numerical aperture (NA), a measure of how wide a cone of light the lens collects from the specimen. Substituting the middle of the visible range, $\lambda = 550$ nm, together with $\mathrm{NA} = 1.4$ for an oil-immersion objective gives $d \approx 196$ nm.

The equation says something simple. To improve resolution, shorten the wavelength or enlarge the numerical aperture. But $\mathrm{NA} = n\sin\alpha$ tops out at 1.4–1.5 for ordinary oil-immersion objectives, and even specialist objectives that pair high-index immersion fluid with high-index cover glass stop near 1.7. Wavelength is the only remaining route, and visible light closes that route too.

## 2. If an electron is a wave, what is its wavelength?

Shortening the wavelength means abandoning light. Using electrons became thinkable because de Broglie assigned a wavelength to matter.

$$ \lambda = \frac{h}{p} $$

Here $h$ is Planck's constant and $p$ the momentum. The larger the momentum, the shorter the wavelength. But does an electron really behave like a wave? In 1927 Davisson and Germer, and independently G. P. Thomson, fired electrons at crystals and recorded diffraction patterns just like those of light. Electron microscopy rests on that property.

Accelerating an electron through a potential $V$ gives it kinetic energy $eV$, so in the non-relativistic regime the momentum is $p = \sqrt{2m_0 eV}$ and the wavelength becomes

$$ \lambda = \frac{h}{\sqrt{2m_0 eV}} $$

with $m_0$ the electron rest mass and $e$ the elementary charge. The wavelength falls as the inverse square root of the voltage: quadruple the voltage and the wavelength halves. A single knob, the accelerating voltage, sets the wavelength. Even 100 V yields 123 pm, about one four-thousandth of visible light. This is the basis on which electron microscopy stands.

## 3. Why does relativity matter at high voltage?

Real TEMs run at 100 kV to 300 kV. Does the expression above still hold there?

The rest-mass energy of an electron, $m_0c^2$, is 511 keV. At 200 kV the electron gains nearly 40 % of that and travels at 0.70 times the speed of light. Mass increase can no longer be ignored. Solving the relativistic relation $E^2 = (pc)^2 + (m_0c^2)^2$ for momentum attaches a correction term.

$$ \lambda = \frac{h}{\sqrt{2m_0 eV\left(1 + \dfrac{eV}{2m_0c^2}\right)}} $$

The bracketed factor is the relativistic correction. At low voltage it approaches unity and recovers the previous expression; at high voltage it pulls the wavelength further down.

```python
# core of generate_figures.py (full code: _code/electron-microscopy-why-electrons/)
H, M0, QE, C = 6.62607015e-34, 9.1093837015e-31, 1.602176634e-19, 2.99792458e8

def wavelength(volts, relativistic=True):
    """De Broglie wavelength [m] of an electron accelerated through V volts."""
    p2 = 2.0 * M0 * QE * np.asarray(volts, dtype=float)
    if relativistic:
        p2 = p2 * (1.0 + QE * np.asarray(volts) / (2.0 * M0 * C**2))
    return H / np.sqrt(p2)
```

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig1-wavelength-vs-voltage.png" alt="Electron de Broglie wavelength against accelerating voltage, and the error of the non-relativistic approximation" width="750">
_Fig 1. De Broglie wavelength of an electron against accelerating voltage (left), and how far the non-relativistic approximation overestimates it (right). The solid line is the relativistic expression and the dashed line the non-relativistic one. At 200 kV they differ by 9.3 %._

The two curves on the left nearly coincide below 30 kV and separate above it. The right-hand panel expresses that separation as a percentage: an error of 1.5 % at 30 kV grows to 9.3 % at 200 kV and 13.7 % at 300 kV. For the 1–30 kV range typical of an SEM, the non-relativistic form causes little harm. In the TEM range (100–300 kV), dropping the correction misstates the wavelength by 5–14 %.

At the standard TEM condition of 200 kV, then, the electron wavelength is 2.51 pm. Interatomic spacings run around 0.2 nm, or 200 pm, so the wavelength is one eightieth of the spacing between atoms.

## 4. Why does resolution stall at 0.2 nm when the wavelength is picometres?

Now put that wavelength into Abbe's equation. What happens?

Even at $\mathrm{NA} = 1$ the result is $d = \lambda/2 = 1.3$ pm, more than 150 times smaller than the spacing between atoms. What was actually achieved looks nothing like it.

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig2-resolution-gap.png" alt="Electron de Broglie wavelength compared with electron microscope resolution achieved in each era" width="750">
_Fig 2. Electron microscope resolution achieved over time (bars, log scale) against the 2.5 pm wavelength of a 200 kV electron (blue dashed line). Hatched bars are commonly quoted estimates that could not be traced to a primary source. No bar reaches the wavelength, and uncorrected HRTEM (0.2 nm) sits about eighty times above it._

The first electron microscope was a two-stage instrument built by Knoll and Ruska in 1931. Ruska's 1933 TEM was the first to surpass the optical microscope. Its resolution is often quoted at around 50 nm, although that figure does not appear in the primary sources. The 1970s are often said to have reached about 0.3 nm. Before aberration correctors became commercial, high-resolution TEM (HRTEM) sat near 0.2 nm for years. Aberration-corrected STEM opened up 0.05 nm, and recent electron ptychography reports about 0.02 nm (15–20 pm, depending on how resolution is measured).

Not one bar reaches down to the dashed line at 2.5 pm. Measured against uncorrected HRTEM, the ratio of resolution to wavelength is about eighty. An optical microscope operates at roughly 0.4 times its wavelength; the electron microscope, by comparison, makes little use of its own.

Fig 3 shows what a commercial electron microscope looked like in 1939, soon after the first one overtook the optical microscope. Electrons leave the gun at the top and travel down a vertical stack of lenses, the column, to form an image near the bottom.

<img src="/assets/img/posts/electron-microscopy-why-electrons/ext-siemens-1939-microscope.jpg" alt="The 1939 Siemens commercial electron microscope at the Deutsches Museum" width="400">
_Fig 3. Siemens' first commercial electron microscope (1939, designed by Ruska), on display at the Deutsches Museum. It was the first electron microscope sold with a resolving power better than a light microscope; the lenses are stacked inside the vertical column. (Source: J Brew, [Ernst Ruska Electron Microscope - Deutsches Museum - Munich-edit](https://commons.wikimedia.org/wiki/File:Ernst_Ruska_Electron_Microscope_-_Deutsches_Museum_-_Munich-edit.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

## 5. Where does the gap come from?

What opens that factor of eighty? Abbe's equation already holds the clue.

Resolution depends not on wavelength alone but on the ratio of wavelength to numerical aperture. Optical microscopy pushed $\mathrm{NA}$ to 1.4. What does an electron objective manage? The aperture semi-angle $\alpha$ in practical use is about 10 mrad, or 0.01 rad. At such small angles $\sin\alpha \approx \alpha$, so the numerical aperture is 0.01, less than a hundredth of the optical 1.4.

Writing the diffraction-limited resolution in terms of the Airy disk gives

$$ d_d = 0.61\,\frac{\lambda}{\alpha} $$

Substituting $\lambda = 2.51$ pm and $\alpha = 10$ mrad returns $d_d \approx 0.15$ nm. This counts diffraction alone; spherical aberration makes it worse. Even so, it is the same order as the 0.2 nm that HRTEM actually delivers. The eighty-fold gap therefore does not arise because diffraction grew with wavelength. It arises because the numerical aperture is less than a hundredth of an optical lens.

Why not simply open the aperture? Here electron optics parts ways with the optics of light. [Ray Tracing 1](/en/posts/raytracing-spherical-lens-refraction/) traced spherical aberration in a spherical lens, where rays further from the axis cross the axis earlier. Light optics cancels that aberration by combining elements of differing index and curvature, so that a concave element undoes what a convex one introduced.

Electron optics has no such concave element. Scherzer proved in 1936 that an electron lens which is rotationally symmetric, static, free of space charge on the axis and forming a real image always has a positive spherical aberration coefficient $C_s$. Opening the aperture therefore inflates the aberration disc $C_s\alpha^3$ as the cube of the angle, cancelling whatever diffraction gave back. The figure of about 10 mrad comes from the place where diffraction and spherical aberration pull against each other ([post 2](/en/posts/electron-microscopy-electron-optics/) works it out as 6.5 mrad). Instrument builders did not pick it arbitrarily.

Post 2 sets up that balance, derives the optimum aperture angle and the resolution it implies, goes through Scherzer's four conditions one by one, and explains what an aberration corrector does to sidestep the theorem.

## 6. How does one choose an electron gun?

A short wavelength is not enough. Talking only about resolution misses something: forming an image requires enough electrons to arrive.

A scanning instrument focuses electrons onto a single point on the specimen and counts the signal leaving that point. The smaller the point, the fewer electrons inside it, and a weak signal drowns in noise. A TEM, which floods the whole specimen at once, is no better off: the more the image is magnified, the fewer electrons land on each detector pixel. Either way, resolution without signal produces no image. This is why the quantity that matters for an electron source is not total current but brightness.

$$ \beta = \frac{I}{A\,\Omega} $$

Current per unit area per unit solid angle. At a fixed accelerating voltage, and ignoring lens aberrations and Coulomb repulsion between electrons, brightness is conserved as lenses demagnify the beam and apertures trim it. How much current can finally be packed into a small probe is therefore decided at the source.

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig3-electron-guns.png" alt="Brightness, energy spread and source size for tungsten and LaB6 thermionic guns and Schottky and cold field emission guns" width="800">
_Fig 4. Brightness (left), energy spread (centre) and source size (right) of four electron guns. Brightness and energy spread are the 100 kV values of Williams & Carter, Table 5.1. A cold field emission gun is three orders of magnitude brighter than a tungsten thermionic gun, with a tenth of its energy spread._

Electrons are extracted in two broad ways. Thermionic emission heats a metal until thermal energy carries electrons over the work function; field emission applies a strong field that thins the barrier until they tunnel through. The four guns in Fig 4 sit between these two. A tungsten filament and a LaB₆ crystal are thermionic sources. LaB₆ has a lower work function, so it emits readily at a lower temperature and is tens of times brighter than tungsten. A cold field emission gun (CFEG) pulls electrons out of a sharp tungsten tip at room temperature by field alone. A Schottky field emission gun lies in between: it keeps the tip hot and uses the field to lower the barrier (Fig 5).

<img src="/assets/img/posts/electron-microscopy-why-electrons/ext-schottky-emitter.jpg" alt="A Schottky field emitter viewed at increasing magnification" width="760">
_Fig 5. A Schottky emitter at increasing magnification. On the far left a V-shaped tungsten filament glows at about 1700 K with 2.2 A flowing; the three electron micrographs to the right zoom in on the single-crystal tungsten tip welded to it (scale bars 1 mm, 100 µm, 1 µm). The ZrO₂ reservoir partway up the tip lowers the work function of its surface, and W{100} marks the tungsten {100} facet at the apex. (Source: ErwinMeier, [Schottky-Emitter 01](https://commons.wikimedia.org/wiki/File:Schottky-Emitter_01.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

The left panel of Fig 4 puts a tungsten thermionic gun near $10^6$ A/cm²·sr against $10^9$ for a cold field emission gun, both at 100 kV. That is three orders of magnitude. The energy spread $\Delta E$ in the middle panel follows the same trend, 3 eV for tungsten against 0.3 eV for cold field emission. Since energy spread feeds directly into chromatic aberration, that difference reaches resolution as well. Source size on the right shrinks from 50 µm to 5 nm, four orders of magnitude.

Cold field emission for everything, then? There is a price. A thermionic gun operates around $10^{-3}$ Pa, whereas a field emission gun demands ultra-high vacuum near $10^{-8}$ Pa. A single adsorbed monolayer on the emitter tip alters the work function and makes the current drift. This is why a CFEG must periodically heat its tip to drive the adsorbed layer off.

It is also worth separating brightness from sheer output. Brightness is current per unit area per unit solid angle, so a cold field emitter, whose emitting area is extremely small, ranks highest in brightness while falling short of a thermionic gun in total current. That is why Schottky field emission is preferred for analytical work in the SEM, where a wide field must be scanned quickly or enough X-rays generated. The choice trades brightness against total current, stability and cost.

## 7. Why must an electron microscope sit in a vacuum?

The previous section noted that a field emission gun needs ultra-high vacuum. Yet the whole column of an electron microscope is evacuated, not just the gun. Why?

Being charged, electrons scatter readily off gas molecules. The mean free path of a gas, the average distance a molecule travels between collisions with other molecules, is

$$ \ell = \frac{k_B T}{\sqrt{2}\,\pi d_m^2 P} $$

with $k_B$ Boltzmann's constant, $T$ the temperature, $d_m$ the molecular collision diameter and $P$ the pressure. A fast electron is effectively pointlike and far faster than the molecules, so the geometric cross-section $\pi d_m^2$ gives way to an energy-dependent electron–gas scattering cross-section $\sigma_T$ (elastic plus inelastic):

$$ \ell_e = \frac{k_B T}{\sigma_T P} $$

Either way, the inverse dependence on pressure is the point.

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig4-mean-free-path.png" alt="Mean free path against pressure, compared with the length of an electron microscope column" width="750">
_Fig 6. Mean free path of air molecules (grey dashed) and of 200 keV electrons (blue solid) against pressure; the horizontal axis runs toward higher vacuum on the right. The red dashed line marks a 1 m column; the electron mean free path crosses it at about 10 Pa._

At atmospheric pressure the mean free path of air molecules is 66 nm. A 200 keV electron sees a cross-section about a thousand times smaller, so its own mean free path is around 0.1 mm; over a column a little over a metre long it would scatter some ten thousand times before reaching the specimen. There would be no image, and no beam to speak of.

The red line in Fig 6 marks a column length of 1 m. The electron mean free path crosses it from roughly 10 Pa, and about 0.1 Pa is enough for scattering to be negligible. At 0.1 Pa the electron mean free path is about 100 m, so only about 1 % of electrons scatter on their way down a 1 m column. That real instruments go far beyond this, down to $10^{-4}$ Pa and below, has less to do with letting electrons through than with specimen contamination and gun lifetime. Residual hydrocarbon molecules cracked by the beam deposit as a carbon film on the specimen surface, darkening it as the observation proceeds. Much of the electron microscope's reputation as a demanding instrument traces back to these vacuum requirements.

## Summary

To collect what this post established: accelerated through 200 kV, an electron has a de Broglie wavelength of 2.51 pm, and omitting the relativistic correction in that range misstates it by more than 9 %. Feeding that wavelength into Abbe's equation predicts picometre resolution, yet what was achieved is around 0.2 nm, a gap of roughly eighty.

The cause turned out to be aperture angle rather than wavelength. Electron microscopes work at an aperture semi-angle of about 10 mrad, which as a numerical aperture is less than a hundredth of an optical objective. The reason it cannot be opened is Scherzer's theorem on the unavoidable spherical aberration of rotationally symmetric electron lenses. Beyond that, forming an image requires enough current in a small probe, which makes source brightness a second constraint, and electrons need vacuum to cross the column without scattering.

[Post 2](/en/posts/electron-microscopy-electron-optics/) enters the electron optics behind this gap. It covers how a lens that focuses electrons with a magnetic field comes to obey the same imaging equation as a glass lens, why nothing equivalent to a concave element can be built, and how minimizing the combined effect of diffraction and spherical aberration yields the optimum aperture angle and a resolution $d_{min} \propto C_s^{1/4}\lambda^{3/4}$. It then turns to how an aberration corrector works around Scherzer's theorem.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 1, 5, 6, 8 (electron wavelength, sources, lenses and resolution, vacuum).
- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 5 (electron guns and brightness).
- C. Davisson and L. H. Germer, "Diffraction of electrons by a crystal of nickel," *Physical Review*, vol. 30, pp. 705–740, 1927 (experimental confirmation of electron diffraction).
- O. Scherzer, "Über einige Fehler von Elektronenlinsen," *Zeitschrift für Physik*, vol. 101, pp. 593–603, 1936 (spherical aberration of rotationally symmetric electron lenses).
- P. E. Batson, N. Dellby, and O. L. Krivanek, "Sub-ångstrom resolution using aberration corrected electron optics," *Nature*, vol. 418, pp. 617–620, 2002.
- C. Hetherington, "[Aberration correction for TEM](https://www.sciencedirect.com/science/article/pii/S1369702104005711)," *Materials Today*, vol. 7, no. 12, pp. 50–55, 2004 (resolution before and after correction).
- [FE electron gun, JEOL glossary](https://www.jeol.com/words/semterms/20121024.062458.php) (brightness and energy spread of field emission guns).
- [Scanning Electron Microscopy, Thermo Fisher Scientific](https://www.thermofisher.com/us/en/home/materials-science/learning-center/applications/scanning-electron-microscope-sem-electron-column.html) (electron column overview).

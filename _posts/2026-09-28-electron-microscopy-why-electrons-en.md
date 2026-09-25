---
title: "Electron Microscopy Foundations 1 — Why Electrons: a 2.5 pm Wavelength and a 0.2 nm Resolution"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-why-electrons/
page_id: electron-microscopy-why-electrons
date: 2026-09-28 20:00:00 +0900
categories: [Electron Microscopy, Electron Optics]
tags: [electron-microscopy, de-broglie, resolution, electron-gun, vacuum]
description: An electron's wavelength is 2.5 pm, yet resolution stalls at 0.2 nm. This post measures that eighty-fold gap and identifies what opens it.
math: true
---

Anyone who works with optical metrology long enough keeps running into the same wall: wavelength. Measuring a film far thinner than the wavelength is possible — [Ellipsometry 1](/en/posts/ellipsometry-electromagnetic-fresnel/) covered exactly that — but converting thickness into a phase and reading it back is not the same as seeing the layer. For resolution, the ability to tell two neighbouring structures apart as two, visible light stops near 200 nm.

So one looks for a shorter wavelength, and the electron presents itself. Treated as a wave, an electron has a wavelength five orders of magnitude shorter than visible light. On paper, atomic nuclei should be within reach. Yet a real transmission electron microscope (TEM) delivers about 0.2 nm, and even an aberration-corrected instrument reaches only 0.05 nm. That is dozens of times worse than the wavelength.

This series starts from that gap. The history of electron microscopy is less a story about obtaining a short wavelength than a record of fighting everything that prevents that wavelength from becoming resolution. This post measures the gap precisely and sketches the outline of the culprit. Its identity, and what can be done about it, belong to post 2 onward.

## 1. Where the optical microscope stops

Fix the baseline first. Where does the claim that visible light stops at 200 nm come from?

Abbe's limit on resolution reads

$$ d = \frac{\lambda}{2\,\mathrm{NA}} $$

where $\lambda$ is the wavelength and $\mathrm{NA} = n\sin\alpha$ is the numerical aperture, a measure of how wide a cone of light the lens collects from the specimen. Substituting the middle of the visible range, $\lambda = 550$ nm, together with $\mathrm{NA} = 1.4$ for an oil-immersion objective gives $d \approx 196$ nm.

The equation says something simple. To improve resolution, shorten the wavelength or enlarge the numerical aperture. But $\mathrm{NA} = n\sin\alpha$ caps out near 1.5 no matter how high the index of the immersion fluid, since $\sin\alpha$ cannot exceed one. Wavelength is the only remaining route, and visible light closes that route too.

## 2. Treating an electron as a wave

Shortening the wavelength means abandoning light. Using electrons became thinkable because de Broglie assigned a wavelength to matter.

$$ \lambda = \frac{h}{p} $$

The larger the momentum $p$, the shorter the wavelength. Accelerating an electron through a potential $V$ gives it kinetic energy $eV$, so in the non-relativistic regime the momentum is $p = \sqrt{2m_0 eV}$ and the wavelength becomes

$$ \lambda = \frac{h}{\sqrt{2m_0 eV}} $$

A single knob, the accelerating voltage, sets the wavelength. Even 100 V yields 123 pm, about one four-thousandth of visible light. This is the basis on which electron microscopy stands.

## 3. Relativity enters as the voltage rises

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
_Fig 1. De Broglie wavelength against accelerating voltage (left), and how far the non-relativistic approximation overestimates it (right)_

The two curves on the left nearly coincide below 30 kV and separate above it. The right-hand panel expresses that separation as a percentage: an error of 1.5 % at 30 kV grows to 9.3 % at 200 kV and 13.7 % at 300 kV. For the 1–30 kV range typical of a scanning electron microscope (SEM), the non-relativistic form causes little harm. In the TEM range, dropping the correction misstates the wavelength by more than ten percent.

At the standard TEM condition of 200 kV, then, the electron wavelength is 2.51 pm. Interatomic spacings run around 0.2 nm, or 200 pm, so the wavelength is one eightieth of the spacing between atoms.

## 4. Picometres on paper, nanometres in practice

Now put that wavelength into Abbe's equation. What happens?

Even at $\mathrm{NA} = 1$ the result is $d = \lambda/2 = 1.3$ pm, a figure approaching nuclear dimensions. What was actually achieved looks nothing like it.

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig2-resolution-gap.png" alt="Electron de Broglie wavelength compared with electron microscope resolution achieved in each era" width="750">
_Fig 2. The wavelength of a 200 kV electron (dashed line) against resolution actually achieved_

Ruska's first TEM, built in 1933, managed around 50 nm. The 1970s brought 0.3 nm. High-resolution TEM (HRTEM) sat near 0.2 nm for years before aberration correctors became commercial. Aberration-corrected scanning TEM (STEM) opened up 0.05 nm, and recent electron ptychography reports the 0.02 nm range.

Not one bar reaches down to the dashed line at 2.5 pm. Measured against uncorrected HRTEM, the ratio of resolution to wavelength is about eighty. An optical microscope operates at roughly 0.4 times its wavelength; the electron microscope, by comparison, barely exploits its own.

## 5. Where the gap comes from

What opens that factor of eighty? Abbe's equation already holds the clue.

Resolution depends not on wavelength alone but on the ratio of wavelength to numerical aperture. Optical microscopy pushed $\mathrm{NA}$ to 1.4. What does an electron objective manage? The aperture semi-angle $\alpha$ in practical use is about 10 mrad, or 0.01 rad. At such small angles $\sin\alpha \approx \alpha$, so the numerical aperture is 0.01 — smaller than an optical objective by a factor of 140.

Writing the diffraction-limited resolution in terms of the Airy disc gives

$$ d_d = 0.61\,\frac{\lambda}{\alpha} $$

Substituting $\lambda = 2.51$ pm and $\alpha = 10$ mrad returns $d_d \approx 0.15$ nm, the same order as the 0.2 nm that HRTEM actually delivers. The eighty-fold gap therefore does not arise because diffraction grew with wavelength. It arises because the aperture angle is stuck at one hundredth of an optical lens.

Why not simply open the aperture? Here electron optics parts ways with the optics of light. [Geometrical Optics 1](/en/posts/raytracing-spherical-lens-refraction/) traced spherical aberration in a spherical lens, where rays further from the axis cross the axis earlier. Light optics cancels that aberration by combining elements of differing index and curvature — a concave element undoing what a convex one introduced.

Electron optics has no such concave element. Scherzer proved in 1936 that any electron lens that is rotationally symmetric, static and free of space charge has a strictly positive spherical aberration coefficient $C_s$. Opening the aperture therefore inflates the aberration disc $C_s\alpha^3$ as the cube of the angle, consuming whatever diffraction gave back. The value of 10 mrad marks where diffraction and spherical aberration balance; it is not a number instrument builders settled for out of laziness.

Setting up that balance, deriving the optimum aperture angle and the resolution it implies, and explaining what an aberration corrector does to sidestep Scherzer's theorem, is the content of post 2.

## 6. A short wavelength is not enough — the electron gun

Talking only about resolution misses something. Forming an image requires enough electrons to arrive.

A scanning instrument focuses electrons onto a single point on the specimen and counts the signal leaving that point. The smaller the point, the fewer electrons inside it, and a weak signal drowns in noise. A TEM, which floods the whole specimen at once, is no better off: the more the image is magnified, the fewer electrons land on each detector pixel. Either way, resolution without signal produces no image. This is why the quantity that matters for an electron source is not total current but brightness.

$$ \beta = \frac{I}{A\,\Omega} $$

Current per unit area per unit solid angle. Brightness is conserved as lenses demagnify the beam and apertures trim it, so how much current can finally be packed into a small probe is decided at the source.

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig3-electron-guns.png" alt="Brightness, energy spread and source size for tungsten and LaB6 thermionic guns and Schottky and cold field emission guns" width="800">
_Fig 3. Brightness, energy spread and source size of four electron guns_

Electrons are extracted in two broad ways. Thermionic emission heats a metal until thermal energy carries electrons over the work function; field emission applies a strong field that thins the barrier until they tunnel through.

The left panel of Figure 3 puts a tungsten thermionic gun near $10^5$ A/cm²·sr against $10^9$ for a cold field emission gun (CFEG) — four orders of magnitude. The energy spread $\Delta E$ in the middle panel follows the same trend, 2.3 eV for tungsten against 0.3 eV for cold field emission. Since energy spread feeds directly into chromatic aberration, that difference reaches resolution as well. Source size on the right shrinks from 50 µm to 5 nm, another four orders.

Cold field emission for everything, then? There is a price. A thermionic gun operates around $10^{-3}$ Pa, whereas a field emission gun demands ultra-high vacuum near $10^{-8}$ Pa. A single adsorbed monolayer on the emitter tip alters the work function and makes the current drift. This is why a CFEG must periodically heat its tip to drive the adsorbed layer off.

It is also worth separating brightness from sheer output. Brightness is current per unit area per unit solid angle, so a cold field emitter — whose emitting area is extremely small — ranks highest in brightness while falling short of a thermionic gun in total current. That is why Schottky field emission is preferred for analytical work, where a wide field must be scanned quickly or enough X-rays generated. The choice trades brightness against total current, stability and cost.

## 7. An electron travels centimetres in air

Since vacuum has come up, it deserves its own note. Electron microscopes are vacuum instruments for more than the sake of gun lifetime.

Being charged, electrons scatter readily off gas molecules. The mean free path, the average distance travelled between collisions, is

$$ \ell = \frac{k_B T}{\sqrt{2}\,\pi d_m^2 P} $$

with $d_m$ the molecular collision diameter and $P$ the pressure. The inverse dependence on pressure is the point.

<img src="/assets/img/posts/electron-microscopy-why-electrons/en/fig4-mean-free-path.png" alt="Mean free path against pressure, compared with the length of an electron microscope column" width="750">
_Fig 4. Electron mean free path against vacuum level_

At atmospheric pressure the mean free path is 66 nm. An electron microscope column runs a little over a metre, so in air an electron would scatter more than ten million times before reaching the specimen. There would be no image, and no beam to speak of.

The red line in Figure 4 marks a column length of 1 m. Crossing it requires $6.7 \times 10^{-3}$ Pa, roughly the $10^{-2}$ Pa level. That real instruments go far beyond this, down to $10^{-4}$ Pa and below, has less to do with letting electrons through than with specimen contamination and gun lifetime. Residual hydrocarbon molecules cracked by the beam deposit as a carbon film on the specimen surface, darkening it as the observation proceeds. Much of the electron microscope's reputation as a demanding instrument traces back to these vacuum requirements.

## Summary and what comes next

To collect what this post established: accelerated through 200 kV, an electron has a de Broglie wavelength of 2.51 pm, and omitting the relativistic correction in that range misstates it by more than 9 %. Feeding that wavelength into Abbe's equation predicts picometre resolution, yet what was achieved is around 0.2 nm — a gap of roughly eighty.

The cause turned out to be aperture angle rather than wavelength. Electron microscopes work at an aperture semi-angle of 10 mrad, one hundredth of an optical objective, and the reason it cannot be opened is Scherzer's theorem on the unavoidable spherical aberration of rotationally symmetric electron lenses. Beyond that, forming an image requires enough current in a small probe, which makes source brightness a second constraint, and electrons need vacuum to cross the column without scattering.

Post 2 enters the electron optics behind this gap. It covers how a lens that focuses electrons with a magnetic field comes to obey the same imaging equation as a glass lens, why nothing equivalent to a concave element can be built, and how minimizing the combined effect of diffraction and spherical aberration yields the optimum aperture angle and a resolution $d_{min} \propto C_s^{1/4}\lambda^{3/4}$. It then turns to how an aberration corrector works around Scherzer's theorem.

## References

- D. B. Williams and C. B. Carter, *Transmission Electron Microscopy: A Textbook for Materials Science*, 2nd ed., Springer, 2009, ch. 1–6 (electron wavelength, sources, vacuum).
- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 2 (electron guns and brightness).
- O. Scherzer, "Über einige Fehler von Elektronenlinsen," *Zeitschrift für Physik*, vol. 101, pp. 593–603, 1936 (spherical aberration of rotationally symmetric electron lenses).
- P. E. Batson, N. Dellby, and O. L. Krivanek, "Sub-ångstrom resolution using aberration corrected electron optics," *Nature*, vol. 418, pp. 617–620, 2002.
- [Aberration correction for TEM, *Materials Today*](https://www.sciencedirect.com/science/article/pii/S1369702104005711) (resolution before and after correction).
- [FE electron gun, JEOL glossary](https://www.jeol.com/words/semterms/20121024.062458.php) (brightness and energy spread of field emission guns).
- [Scanning Electron Microscopy, Thermo Fisher Scientific](https://www.thermofisher.com/us/en/home/materials-science/learning-center/applications/scanning-electron-microscope-sem-electron-column.html) (electron column overview).

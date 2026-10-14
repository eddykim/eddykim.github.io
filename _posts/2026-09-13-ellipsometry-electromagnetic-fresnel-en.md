---
title: "Ellipsometry Foundations 1 — If a Detector Measures Only Intensity, How Do We Learn the Phase of Light?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometry-electromagnetic-fresnel/
page_id: ellipsometry-electromagnetic-fresnel
date: 2026-09-13 20:00:00 +0900
categories: [Optics, Electromagnetics]
tags: [ellipsometry, electromagnetics, maxwell-equations, fresnel-equations, polarization, thin-film]
description: The mathematical chain from Maxwell's equations to the Fresnel reflection coefficients, and why the fact that a detector measures only intensity gives rise to ellipsometry as a measurement technique.
math: true
---

Ellipsometry determines the thickness and refractive index of a thin film by measuring how the polarization state of light changes on reflection from a sample surface. Understanding it means returning to the most basic question in electromagnetics: how light reflects and refracts at the boundary between two media. This series covers that background over three posts. This first one follows the mathematical chain from Maxwell's equations to the Fresnel equations. [Post 2](/en/posts/ellipsometry-polarization-mueller-matrix/) treats the Stokes vector and Mueller matrix that describe a polarization state, and [Post 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/) covers reflection from a multilayer thin-film stack and the ellipsometric parameters $\Psi$ and $\Delta$.

There is one point in particular worth watching for here. The Fresnel equations themselves are a standard result found in any optics textbook. But answering the question "why is a technique that compares polarization states needed at all?" requires pausing at the Poynting vector, midway through their derivation. A detector does not measure the electric field itself, and what is lost as a result is the motivation running through this entire series.

## 1. How do Maxwell's equations produce a wave?

Where does the statement that light is an electromagnetic wave come from? The starting point is Maxwell's equations. In vacuum, free of charge and current, they take the following form.

$$ \nabla \cdot \mathbf{E} = 0, \qquad \nabla \cdot \mathbf{B} = 0 $$

$$ \nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t}, \qquad \nabla \times \mathbf{B} = \mu_0 \epsilon_0 \frac{\partial \mathbf{E}}{\partial t} $$

Here $\mathbf{E}$ is the electric field, $\mathbf{B}$ the magnetic field, and $\mu_0$ and $\epsilon_0$ the permeability and permittivity of free space. Combining these four gives a wave equation satisfied separately by the electric and magnetic fields. Take the curl of each curl equation, then apply the identity $\nabla\times(\nabla\times\mathbf{A}) = \nabla(\nabla\cdot\mathbf{A}) - \nabla^2\mathbf{A}$ together with the divergence conditions $\nabla\cdot\mathbf{E}=\nabla\cdot\mathbf{B}=0$.

$$ \nabla^2 \mathbf{E} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2} = 0, \qquad \nabla^2 \mathbf{B} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{B}}{\partial t^2} = 0 $$

In other words, the electric and magnetic fields are waves travelling at a speed fixed by two vacuum constants alone, $c = 1/\sqrt{\mu_0\epsilon_0} \approx 3.00\times10^{8}\,\mathrm{m/s}$. That number is the speed of light.

The solutions to this wave equation take the form of plane waves.

$$ \mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0\, e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}, \qquad \mathbf{B}(\mathbf{r}, t) = \mathbf{B}_0\, e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})} $$

Here $\mathbf{k}$ is the wave vector and $\omega$ the angular frequency, related by $\lvert \mathbf{k} \rvert = \omega/c$. Taking the time factor as $e^{j\omega t}$ follows the convention widely used in the ellipsometry literature since the 1968 Nebraska conference (Azzam & Bashara, Fujiwara). The literature is not uniform on this, however. The sign convention matters again in section 3, when the complex refractive index is defined, so it is fixed here in advance.

The exponential form is used because it makes calculation easier. The physical quantity is the real part of this complex expression, but differentiation and phase arithmetic are far simpler in exponential form, so the notation is carried all the way to the Fresnel equations. One important fact is already visible at this point. The electric field is a complex quantity carrying not only the amplitude $\mathbf{E}_0$ but also the phase $\omega t - \mathbf{k}\cdot\mathbf{r}$.

## 2. What does a detector see of the electric field?

If the field carries both amplitude and phase, which of the two does a detector read? To answer, start from the energy the wave carries. That energy flow is given by the Poynting vector $\mathbf{S}$.

$$ \mathbf{S} = \frac{1}{\mu_0}(\mathbf{E} \times \mathbf{B}) $$

The angular frequency of visible light is on the order of $10^{15}\,\mathrm{rad/s}$; at 633 nm one period lasts about 2.1 fs. No detector, whether photodiode, CCD or CMOS, can follow a field oscillating at that rate in real time. What a detector actually responds to is the energy flow integrated over the exposure, that is, its time average.

$$ \langle \mathbf{S} \rangle = \frac{1}{2}\,\mathrm{Re}(\mathbf{E} \times \mathbf{H}^*) $$

Here $\mathbf{H} = \mathbf{B}/\mu_0$, and $\mathbf{E}$, $\mathbf{H}$ are complex amplitudes. Expanding this for a plane wave in vacuum, the intensity $I$ the detector measures reduces to a scalar proportional to the squared magnitude of the field amplitude.

$$ I = \lvert \langle \mathbf{S} \rangle \rvert = \frac{1}{2}\sqrt{\frac{\epsilon_0}{\mu_0}}\, \lvert \mathbf{E}_0 \rvert^2 $$

Inside a medium of refractive index $n$ the right-hand side gains a factor $n$, but the intensity is still proportional to $\lvert \mathbf{E}_0 \rvert^2$.

This expression is the starting point of the whole series. The intensity $I$ keeps only the amplitude $\lvert \mathbf{E}_0 \rvert$ and erases the phase completely in the squaring and time-averaging. Two waves of equal amplitude but different phase therefore look identical to a detector. Measurements that record only intensity, such as reflectance and transmittance, are usually called photometry. From the intensity at a single wavelength they have no access to the phase. A Kramers–Kronig analysis of a broadband reflectance spectrum can recover it, but only with extrapolation beyond the measured range and for a bulk sample. Section 7 returns to this point.

## 3. Why does light bend at an interface?

Now consider what happens when this wave meets a boundary. When a plane wave meets the interface between two media of different refractive index, part of it is reflected and part refracted into the second medium. This is not a new physical law. It follows directly from the Maxwell boundary condition that the tangential components of $\mathbf{E}$ and $\mathbf{H}$ be continuous across the interface. For the non-magnetic media treated here, continuity of the tangential $\mathbf{B}$ says the same thing.

The plane containing the propagation direction of the incident light and the surface normal is called the plane of incidence. Every polarization component is decomposed into s-polarization, perpendicular to that plane, and p-polarization, parallel to it. The angles of incidence $\theta_i$, reflection $\theta_r$ and refraction $\theta_t$ are measured from the normal. Fig 1 collects these definitions. The ⊙ symbol marks oscillation out of the page (s), and the double arrow marks oscillation in the page, perpendicular to the ray (p).

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/en/fig1-concept-diagram.png" alt="s/p polarization and the incident, reflected and refracted rays on the plane of incidence" width="600">
_Fig 1. Incident, reflected and refracted rays on the plane of incidence (the page). s-polarized light oscillates perpendicular to the page (⊙); p-polarized light oscillates in the plane of incidence, perpendicular to the ray (↔). Medium 2 has the higher index, so the refracted ray bends toward the normal (schematic angles $$\theta_i = 35°$$, $$\theta_t = 22°$$)._

The reflection angle always equals the incidence angle, and the refraction angle is set by Snell's law.

$$ \theta_i = \theta_r, \qquad N_1 \sin\theta_i = N_2 \sin\theta_t $$

Here $N_1$ and $N_2$ are the refractive indices of the two media. A real absorbing material is described by a complex refractive index $N = n - jk$. Its real part $n$ is the ordinary refractive index, and its imaginary part $k$, the extinction coefficient, expresses how strongly the material absorbs light. Note that this $k$ shares its letter with the bold wave vector $\mathbf{k}$. The minus sign pairs with the $e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}$ convention fixed in section 1: substituting the wavenumber $\lvert \mathbf{k} \rvert = N\omega/c$ under that convention yields a decay $e^{-k\omega z/c}$ along the direction of propagation. Under the physics convention, with a time factor $e^{-j\omega t}$, the same decay requires $N = n + jk$ instead. Either is valid, but mixing the two turns absorption into gain, so one convention is used throughout.

This is the same Snell's law that the lens refraction calculation of [Ray Tracing 1](/en/posts/raytracing-spherical-lens-refraction/) took as given without asking why it holds. The answer is here. At any point on the interface, the wave above and the wave below cannot have different phases. So the spacing of wave crests, measured along the interface, must be the same on both sides. In the higher-index medium 2 the wavelength is shorter, and the wavefronts can keep the same spacing only by tilting. That tilt is the refraction angle.

In equations: the plane-wave solution $e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}$ must join in phase across the incident, reflected and refracted waves at every point on the interface and at every instant. That requires the projections of their wave vectors along the interface to agree. This condition is called phase matching. Written in terms of angles it reads $\lvert \mathbf{k}_1 \rvert \sin\theta_i = \lvert \mathbf{k}_2 \rvert \sin\theta_t$, and with $\lvert \mathbf{k} \rvert = N\omega/c$ it becomes exactly $N_1\sin\theta_i = N_2\sin\theta_t$. Phase matching *is* Snell's law. Fig 2 shows the process in terms of wavefronts.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/ext-wavefront-refraction.png" alt="Wavefronts of a plane wave bending as it crosses an interface" width="450">
_Fig 2. A plane wave passing from a medium of index $$n_1$$ into one of index $$n_2$$. The thick grey lines are the two edges of the beam, and the thin lines (blue above, green below) are wavefronts. The wavefronts continue unbroken at the yellow points on the interface, and because their spacing narrows in the lower medium, the direction of travel (dashed arrow) bends toward the normal. (Source: 老陳, [Wavefront Refraction](https://commons.wikimedia.org/wiki/File:Wavefront_Refraction.svg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

## 4. How much light does an interface reflect?

Snell's law says only at what angle refraction occurs, not how much light is reflected and transmitted. The Fresnel equations set those proportions. The coefficients they give carry both an amplitude ratio and a phase shift, which puts them at the centre of this series. Solving the tangential continuity condition of section 3 separately for s- and p-polarization gives the reflection coefficient $r$, defined as the ratio of the reflected to the incident field amplitude.

$$ r_p = \frac{N_2\cos\theta_i - N_1\cos\theta_t}{N_2\cos\theta_i + N_1\cos\theta_t}, \qquad r_s = \frac{N_1\cos\theta_i - N_2\cos\theta_t}{N_1\cos\theta_i + N_2\cos\theta_t} $$

The transmission coefficients $t_p$ and $t_s$ follow the same way and reappear in Post 3 when calculating multiple reflections inside a film. The point to notice is that $N_1$ and $N_2$ are complex in general. For an absorbing medium $N=n-jk$, so $r_p$ and $r_s$ are complex too. A complex reflection coefficient carries two pieces of information at once: its magnitude, the amplitude ratio, and its argument, the phase shift. This is the first opening on the problem raised in section 2. The phase is alive in the reflection coefficient $r$ all along.

Start with the limiting case. At normal incidence ($\theta_i = 0$), $\cos\theta_i = \cos\theta_t = 1$, so $r_p = (N_2-N_1)/(N_2+N_1)$ and $r_s = -(N_2-N_1)/(N_2+N_1)$. From air ($N_1=1.0$) into glass ($N_2=1.5$), $r_p = +0.2$ and $r_s = -0.2$, and both reflectances are $0.2^2 = 0.04$. At normal incidence nothing distinguishes s from p, so the two polarizations are physically reflected alike. The sign difference $r_p = -r_s$ comes from the convention for the p reference direction. This post uses the same convention as Born & Wolf and Fujiwara, as do Posts 2 and 3. Under it, the 180° phase difference between the two reflected waves at normal incidence serves as the reference value.

## 5. Why does p-reflectance vanish at one angle?

The reflection coefficient carries phase. How, then, does the reflectance, which is what intensity shows, change with angle? The measurable reflectance $R$ is the squared magnitude of the reflection coefficient.

$$ R_p = \lvert r_p \rvert^2, \qquad R_s = \lvert r_s \rvert^2 $$

Computing $R_p$ and $R_s$ against angle of incidence for light passing from air ($N_1=1.0$) into glass ($N_2=1.5$) gives Fig 3.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/en/fig2-fresnel-reflectance.png" alt="Fresnel reflectance R_p and R_s against angle of incidence, with the Brewster angle marked" width="600">
_Fig 3. Fresnel reflectance for incidence from air ($$N_1=1.0$$) onto glass ($$N_2=1.5$$). $$R_s$$ (purple) rises monotonically with angle, while $$R_p$$ (green) falls to zero at the Brewster angle of 56.3° (dashed line). Both equal 0.04 at normal incidence._

The two curves start together at 0.04. $R_s$ increases monotonically with angle, but $R_p$ falls, touches zero at one angle, and rises again. That angle is the Brewster angle. Solving the condition that the numerator of $r_p$ vanish, $N_2\cos\theta_i = N_1\cos\theta_t$, together with Snell's law gives $\theta_B = \arctan(N_2/N_1)$, which is 56.3° for the air-glass pair. At the Brewster angle no p-polarized light is reflected, and the reflected light is purely s-polarized. The reflected and refracted rays are then exactly perpendicular (the refraction angle in glass is 33.7°, and 56.3° + 33.7° = 90°). This property is used in practice to calibrate the azimuth reference of a polarizer, that is, the direction of the plane of incidence. Fig 4 shows the situation.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/ext-brewster-angle.png" alt="Reflection and refraction at the Brewster angle; the reflected light is purely s-polarized" width="500">
_Fig 4. Unpolarized light incident at the Brewster angle $$\theta_B$$, with its reflected and refracted rays. Dots mark oscillation perpendicular to the page (s) and double arrows oscillation in the plane of incidence (p). Only dots remain on the reflected ray, so it is s-polarized; the refracted ray keeps both components and is only slightly polarized, and the two rays meet at a right angle. (Source: Pajs, [Brewsters-angle](https://commons.wikimedia.org/wiki/File:Brewsters-angle.svg), public domain)_

There is one thing Fig 3 cannot show: the sign of $r_p$. In glass, $r_s$ grows more negative with angle (−0.303 at 45°, −0.547 at 70°). $r_p$, by contrast, shrinks and passes through zero at the Brewster angle into negative values (+0.092 at 45°, −0.206 at 70°). Crossing the Brewster angle therefore flips the phase difference between the two reflected waves from 180° to 0°. Reflectance does not see this. $R_p$ is 0.0085 at both 45° and 63.6°, yet $r_p$ is +0.092 at one and −0.092 at the other. The conclusion of section 2, that intensity erases phase, appears here as numbers.

An absorbing sample behaves a little differently. When $N$ is complex the numerator of $r_p$ never vanishes exactly, so $R_p$ reaches only a minimum. That angle is the pseudo-Brewster angle. For Si at 633 nm ($N \approx 3.88 - 0.016j$) the minimum is $R_p \approx 4\times10^{-6}$ at about 75.6°. Near that angle the phase difference changes fastest: 179.4° at 70° and 173.9° at 75°, it passes 90° near 75.6° and reaches 0.6° by 80°. The neighbourhood of the pseudo-Brewster angle is where a small change in the sample shows up most strongly in the phase difference. That is why ellipsometers measuring samples on Si substrates set the angle of incidence in the 70s. Why oblique incidence is needed, and how an instrument uses this angle, is covered in section 1 of [Ellipsometer Instrumentation 1](/en/posts/ellipsometer-oblique-rotating-element/).

## 6. Why are absorbing substrates measured in reflection?

So far the light has split at the interface. Now follow the light that enters the medium. In an absorbing medium the intensity decays exponentially, at a rate set by the extinction coefficient $k$ (the Beer–Lambert law).

$$ I(z) = I_0\, e^{-\alpha_{abs} z}, \qquad \alpha_{abs} = \frac{4\pi k}{\lambda} $$

Here $\lambda$ is the vacuum wavelength. The depth at which the incident intensity has fallen to $1/e$ is the penetration depth $\delta$.

$$ \delta = \frac{1}{\alpha_{abs}} = \frac{\lambda}{4\pi k} $$

This explains why thin-film metrology is mostly done in reflection. Where the substrate absorbs, $\delta$ is far shorter than the substrate thickness: about 3 µm for silicon at 633 nm, against a wafer several hundred µm thick. The light does not pass through the sample, and transmission measurement becomes impossible. This is one reason ellipsometry of films on silicon and other absorbing substrates is built around a reflection geometry.

## 7. What do we gain by measuring polarization instead of intensity?

The thread left hanging in section 2 can now be picked up. The intensity a detector records, $I \propto \lvert E_0 \rvert^2$, carries only the magnitude $\lvert r \rvert$ of the reflection coefficient; the phase of $r$ vanishes in the measurement. However precisely the reflectance $R=\lvert r \rvert^2$ is measured at a single wavelength, it gives no direct access to the phase difference between the s- and p-polarized reflected waves. Yet, as section 5 showed, that phase difference is exactly the quantity that responds sensitively to film thickness and optical constants.

This is why ellipsometry measures the change in polarization state rather than intensity. Preparing the incident light in a known polarization state (usually linear) and measuring how that state has changed after reflection yields the ratio of the two reflection coefficients:

$$ \rho = \frac{r_p}{r_s} = \tan\Psi \, e^{j\Delta} $$

Here $\tan\Psi$ is the amplitude ratio of the two coefficients and $\Delta$ their phase difference. The phase difference followed through sections 4 and 5 is this $\Delta$. A reflectance measurement can at best obtain $\lvert r_p \rvert$ and $\lvert r_s \rvert$ separately. A change in polarization state, by contrast, is a relative quantity that compares the two components against each other, so it yields the magnitude ratio and the phase difference together. The inaccessibility of phase found in section 2 is circumvented by comparing the two polarization components. Ellipsometry gives up absolute intensity and takes the ratio of two polarizations instead.

This pair $\Psi, \Delta$ is what ellipsometry actually delivers. Recovering thickness and refractive index from it takes two steps. First, a thin-film reflection model is built that computes $\Psi, \Delta$ from a given thickness and index; that is done in [Post 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/). Second, the model is fitted to the measured values by nonlinear least squares, which is the subject of the Levenberg-Marquardt method in [Optimization 3](/en/posts/optimization-levenberg-marquardt/).

## Summary and what comes next

This post followed a single mathematical chain from Maxwell's equations in vacuum through the wave equation, the plane-wave solution and the Poynting vector, to Snell's law and the Fresnel equations at an interface. From a metrology standpoint the crucial link was the time average of the Poynting vector: intensity erases the phase of the electric field. The Fresnel coefficients $r_p$ and $r_s$ are complex and do carry phase. At an air-glass boundary the phase difference between the two reflected waves flips from 180° to 0° across the Brewster angle of 56.3°, and on Si it changes fastest near the pseudo-Brewster angle of 75.6°. Ellipsometry recovers this phase difference by measuring the ratio of the two polarization components instead of intensity.

The next post covers the tools that describe a polarization state quantitatively: the Stokes vector and the Mueller matrix. The change in polarization state, noted here only as $\rho = \tan\Psi\, e^{j\Delta}$, is taken up there as the mathematical object it actually is. The instrument-side question, how $\Psi$ and $\Delta$ can be extracted with a detector that reports only intensity, is picked up by [Ellipsometer Instrumentation 1](/en/posts/ellipsometer-oblique-rotating-element/). This Foundations series sets up the calculation; the Instrumentation series covers the hardware that turns it into a measurement.

## References

- Youngjoon Kim, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발" [Development of snapshot angle-resolved ellipsometry using a line-scan spectrometer and back-focal-plane spectral interference], Ph.D. dissertation, Seoul National University, 2025 (in Korean), sections 2.1–2.2.
- E. Hecht, *Optics*, 5th ed., Pearson, 2017, ch. 3–4 (electromagnetic waves and boundary conditions).
- M. Born and E. Wolf, *Principles of Optics*, 7th ed., Cambridge University Press, 1999, ch. 1 (the standard derivation of the Fresnel equations).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, §1.4 (time-factor convention) and §5.3.2 (polarizer azimuth calibration at the Brewster angle).
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 2 and 4 (the Fresnel equations in §2.3; the Brewster and pseudo-Brewster angles in §2.3.4; their link to the ellipsometric parameters in §4.1).
- M. A. Green, "Self-consistent optical parameters of intrinsic silicon at 300 K including temperature coefficients," *Solar Energy Materials and Solar Cells* 92, 1305–1310 (2008) (optical constants of Si at 633 nm).

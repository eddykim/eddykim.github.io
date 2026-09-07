---
title: "Ellipsometry Foundations 1 — From Electromagnetic Waves to the Fresnel Equations"
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

Ellipsometry determines the thickness and refractive index of a thin film by measuring how the polarization state of light changes on reflection from a sample surface. Understanding it means returning to the most basic question in electromagnetics: how light reflects and refracts at the boundary between two media. This series covers that background over three posts. This first one follows the mathematical chain from Maxwell's equations to the Fresnel equations; the second treats the Stokes vector and Mueller matrix that describe a polarization state; the third covers reflection from a multilayer thin-film stack and the ellipsometric parameters $\Psi$ and $\Delta$.

There is one point in particular worth watching for here. The Fresnel equations themselves are a standard result found in any optics textbook, but answering the question "why is a technique that compares polarization states needed at all?" requires pausing at the Poynting vector, midway through their derivation. That the quantity a detector actually measures is not the electric field itself — and what information is lost as a result — is the motivation running through this entire series.

## 1. From Maxwell's equations to a plane wave

The behaviour of electromagnetic waves is described by Maxwell's equations. In a medium free of charge and current — vacuum, or a dielectric — they take the following form.

$$ \nabla \cdot \mathbf{E} = 0, \qquad \nabla \cdot \mathbf{B} = 0 $$

$$ \nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}}{\partial t}, \qquad \nabla \times \mathbf{B} = \mu_0 \epsilon_0 \frac{\partial \mathbf{E}}{\partial t} $$

Here $\mathbf{E}$ is the electric field, $\mathbf{B}$ the magnetic field, and $\mu_0$ and $\epsilon_0$ the permeability and permittivity of free space. Combining these four gives a wave equation satisfied separately by the electric and magnetic fields: take the curl of each curl equation, then apply the identity $\nabla\times(\nabla\times\mathbf{A}) = \nabla(\nabla\cdot\mathbf{A}) - \nabla^2\mathbf{A}$ together with the divergence conditions $\nabla\cdot\mathbf{E}=\nabla\cdot\mathbf{B}=0$.

$$ \nabla^2 \mathbf{E} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{E}}{\partial t^2} = 0, \qquad \nabla^2 \mathbf{B} - \mu_0 \epsilon_0 \frac{\partial^2 \mathbf{B}}{\partial t^2} = 0 $$

The solutions to this wave equation take the form of plane waves.

$$ \mathbf{E}(\mathbf{r}, t) = \mathbf{E}_0\, e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}, \qquad \mathbf{B}(\mathbf{r}, t) = \mathbf{B}_0\, e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})} $$

Here $\mathbf{k}$ is the wave vector and $\omega$ the angular frequency, related by $\lvert \mathbf{k} \rvert = \omega/c$ with $c = 1/\sqrt{\mu_0\epsilon_0}$. Taking the time factor as $e^{j\omega t}$ follows the convention of the ellipsometry literature. That sign convention matters again in section 3, when the complex refractive index is defined, so it is fixed here in advance. The exponential form is used for computational convenience: the physical quantity is the real part of this complex expression, but differentiation and phase arithmetic are far simpler in exponential form, so the notation is carried all the way through to the Fresnel equations. One important fact is already visible at this point. The electric field is a complex quantity carrying not only the amplitude $\mathbf{E}_0$ but also the phase $\omega t - \mathbf{k}\cdot\mathbf{r}$.

## 2. The Poynting vector, and what a detector actually sees

The energy flow carried by an electromagnetic wave is given by the Poynting vector $\mathbf{S}$.

$$ \mathbf{S} = \frac{1}{\mu_0}(\mathbf{E} \times \mathbf{B}) $$

The angular frequency of visible light is on the order of $10^{15}\,\mathrm{rad/s}$. No detector — photodiode, CCD, CMOS, whatever it may be — can follow an electric field oscillating at that rate in real time. What a detector actually responds to is the energy flow integrated over the exposure, that is, its time average.

$$ \langle \mathbf{S} \rangle = \frac{1}{2}\,\mathrm{Re}(\mathbf{E} \times \mathbf{H}^*) $$

Expanding this for a plane wave in an isotropic medium, the intensity $I$ the detector measures reduces to a scalar proportional to the squared magnitude of the field amplitude.

$$ I = \lvert \langle \mathbf{S} \rangle \rvert = \frac{1}{2}\sqrt{\frac{\epsilon_0}{\mu_0}}\, \lvert \mathbf{E}_0 \rvert^2 $$

This expression is the starting point of the whole series. Section 1 ended by noting that the electric field carries both amplitude and phase; the intensity $I$ keeps only the amplitude, $\lvert \mathbf{E}_0 \rvert$, and erases the phase completely in the squaring and time-averaging. Measurements that record only intensity — reflectance, transmittance, what is usually called photometry — have no access in principle to that phase information. This point is taken up again later.

## 3. At the interface: refraction and Snell's law

When a plane wave meets the boundary between two media of different refractive index, part of it is reflected and part refracted into the second medium. This is not a new physical law: it follows directly from the Maxwell boundary condition that the tangential components of the electric and magnetic fields must be continuous across the interface.

The plane containing the propagation direction of the incident light and the surface normal is called the plane of incidence, and every polarization component is decomposed into s-polarization, perpendicular to that plane, and p-polarization, parallel to it. Defining the angles of incidence $\theta_i$, reflection $\theta_r$ and refraction $\theta_t$ with respect to the normal, the reflection angle always equals the incidence angle and the refraction angle is set by Snell's law.

$$ \theta_i = \theta_r, \qquad N_1 \sin\theta_i = N_2 \sin\theta_t $$

Here $N_1$ and $N_2$ are the refractive indices of the two media. A real absorbing material is described by a complex refractive index $N = n - jk$, whose real part $n$ is the ordinary refractive index and whose imaginary part $k$, the extinction coefficient, expresses how strongly the material absorbs light. The minus sign pairs with the $e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}$ convention fixed in section 1: substituting $\mathbf{k} = N\omega/c$ under that convention yields a decay $e^{-k\omega z/c}$ along the direction of propagation. Under the physics convention, with a time factor $e^{-j\omega t}$, the same decay requires writing $N = n + jk$ instead. Either is valid, but mixing the two turns absorption into gain, so one convention is used throughout.

Snell's law here is the very expression carried over into the lens refraction calculation of [Geometrical Optics 1](/en/posts/raytracing-spherical-lens-refraction/). That post used it as a given result without addressing why it holds; the answer is here. The condition that the tangential component of the wave vector match on both sides of the interface — phase matching — *is* Snell's law. For the plane-wave solution $e^{j(\omega t - \mathbf{k}\cdot\mathbf{r})}$ to join continuously in phase across the incident, reflected and refracted waves at every point on the interface and at every instant, the projections of their wave vectors along the interface must agree; writing that condition in terms of angles gives exactly $N_1\sin\theta_i = N_2\sin\theta_t$.

## 4. The Fresnel equations — coefficients carrying both amplitude ratio and phase difference

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/en/fig1-concept-diagram.png" alt="s/p polarization and the reflected and refracted vectors on the plane of incidence" width="600">
_Fig 1. s/p polarization and the incident, reflected and refracted rays defined on the plane of incidence_

Snell's law says only at what angle refraction occurs, not how much light is reflected and transmitted at the interface. The Fresnel equations set those proportions. Solving the continuity condition on the tangential field components separately for s- and p-polarization gives the reflection coefficient $r$, defined as the ratio of the reflected to the incident field amplitude.

$$ r_p = \frac{N_2\cos\theta_i - N_1\cos\theta_t}{N_2\cos\theta_i + N_1\cos\theta_t}, \qquad r_s = \frac{N_1\cos\theta_i - N_2\cos\theta_t}{N_1\cos\theta_i + N_2\cos\theta_t} $$

The transmission coefficients $t_p$ and $t_s$ follow the same way, and these coefficients reappear in the third post when calculating multiple reflections inside a film. The point to notice is that $N_1$ and $N_2$ are complex in general: for an absorbing medium $N=n-jk$, so $r_p$ and $r_s$ are complex too. A complex reflection coefficient carries two pieces of information at once — its magnitude, the amplitude ratio, and its argument, the phase shift. This is the first opening on the problem raised in section 2, that intensity measurement loses the phase: the phase is alive in the reflection coefficient $r$ itself all along.

## 5. Reflectance and the Brewster angle

The measurable reflectance $R$ is the squared magnitude of the reflection coefficient.

$$ R_p = \lvert r_p \rvert^2, \qquad R_s = \lvert r_s \rvert^2 $$

Computing $R_p$ and $R_s$ against angle of incidence for the index pair $N_1=1.0$ (air) and $N_2=1.5$ (glass, as an example) gives the following.

<img src="/assets/img/posts/ellipsometry-electromagnetic-fresnel/en/fig2-fresnel-reflectance.png" alt="Fresnel reflectance R_p and R_s against angle of incidence, with the Brewster angle marked" width="600">
_Fig 2. Fresnel reflectance $R_p$ and $R_s$ (N1=1.0, N2=1.5), with the Brewster angle marked_

$R_s$ increases monotonically with angle, but $R_p$ passes through a point where it is exactly zero. That angle is the Brewster angle; solving the condition that the numerator of $r_p$ vanish, $N_2\cos\theta_i = N_1\cos\theta_t$, together with Snell's law gives $\theta_B = \arctan(N_2/N_1)$ — 56.3° for the air-glass pair. At the Brewster angle no p-polarized component is reflected at all and the reflected light is purely s-polarized, a property used in practice to produce reference light of known polarization for calibrating a source or a detector chain. Pursuing the Brewster angle further is outside the scope of this series; it appears here only as an illustration of the Fresnel equations in use.

## 6. Absorption: the Beer-Lambert law and penetration depth

As light propagates into an absorbing medium its intensity decays exponentially, governed by the extinction coefficient $k$ (the Beer-Lambert law).

$$ I(z) = I_0\, e^{-\alpha_{abs} z}, \qquad \alpha_{abs} = \frac{4\pi k}{\lambda} $$

The depth at which the incident intensity has fallen to $1/e$ is the penetration depth $\delta$.

$$ \delta = \frac{1}{\alpha_{abs}} = \frac{\lambda}{4\pi k} $$

This explains why thin-film metrology is mostly done in reflection. In a wavelength region of strong absorption, $\delta$ falls below the film thickness itself, the light does not pass through the sample, and transmission measurement becomes impossible. That is why ellipsometry of films on silicon and other absorbing substrates is built around a reflection geometry.

## 7. Why ellipsometry

The thread left hanging in section 2 can now be picked up. The intensity a detector records, $I \propto \lvert E_0 \rvert^2$, carries only the magnitude $\lvert r \rvert$ of the reflection coefficient; the phase of $r$ vanishes in the measurement. However precisely the reflectance $R=\lvert r \rvert^2$ is measured, it gives no access to the phase difference between the s- and p-polarized reflected waves — precisely the quantity that responds sensitively to film thickness and optical constants.

This is where ellipsometry's choice to measure the change in polarization state, rather than intensity, comes from. Preparing the incident light in a known polarization state (usually linear) and measuring how that state has changed after reflection yields the ratio of the two reflection coefficients:

$$ \rho = \frac{r_p}{r_s} = \tan\Psi \, e^{j\Delta} $$

Here $\tan\Psi$ is the amplitude ratio of the two coefficients and $\Delta$ their phase difference. A reflectance measurement, recording only intensity — the time-averaged Poynting vector — can at best obtain the magnitudes $\lvert r_p \rvert$ and $\lvert r_s \rvert$ separately, whereas measuring the change in polarization state, a relative quantity, yields the entire ratio of the two coefficients at once: magnitude ratio and phase difference together. The inaccessibility of phase established in section 2 is circumvented not by measuring an absolute quantity but by comparing two polarization components against each other. This pair $\Psi, \Delta$ is what ellipsometry actually measures, and recovering thickness and refractive index from it is the subject of the third post.

## Summary and what comes next

This post followed a single mathematical chain from Maxwell's equations in vacuum through the wave equation, the plane-wave solution, the Poynting vector, and on to Snell's law and the Fresnel equations at an interface. From a metrology standpoint the crucial link in that chain was the time average of the Poynting vector: intensity erases the phase of the electric field. The Fresnel coefficients $r_p$ and $r_s$ are complex and do carry phase, but that information stays out of reach as long as only intensity is measured. Ellipsometry circumvents the limit by measuring the relative ratio of two polarization components — amplitude ratio plus phase difference — instead of an absolute intensity.

The next post covers the tools that describe a polarization state quantitatively: the Stokes vector and the Mueller matrix. The change in polarization state, noted here only as $\rho = \tan\Psi\, e^{j\Delta}$, is taken up there as the mathematical object it actually is, and as something that is actually measured.

## References

- Youngjoon Kim, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발" [Development of snapshot angle-resolved ellipsometry using a line-scan spectrometer and back-focal-plane spectral interference], Ph.D. dissertation, Seoul National University, 2025 (in Korean), sections 2.1–2.2.
- E. Hecht, *Optics*, 5th ed., Pearson, 2016, ch. 3–4 (electromagnetic waves and boundary conditions).
- M. Born and E. Wolf, *Principles of Optics*, 7th ed., Cambridge University Press, 1999, ch. 1 (the standard derivation of the Fresnel equations).
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 2 (linking the ellipsometric parameters to the Fresnel equations).

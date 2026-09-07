---
title: "Ellipsometry Foundations 2 — Polarization, Stokes Vectors and Mueller Matrices"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometry-polarization-mueller-matrix/
page_id: ellipsometry-polarization-mueller-matrix
date: 2026-09-23 20:00:00 +0900
categories: [Optics, Polarization]
tags: [ellipsometry, polarization, jones-vector, stokes-vector, mueller-matrix]
description: What a Jones vector can and cannot express compared with a Stokes vector, and the Mueller matrices of polarizers and retarders mapped onto the hardware they describe.
math: true
---

[Post 1](/en/posts/ellipsometry-electromagnetic-fresnel/) ended by arriving at the ratio of Fresnel reflection coefficients, $\rho = r_p/r_s = \tan\Psi\, e^{j\Delta}$. Taken alone, that single value looks like everything ellipsometry needs. But $\rho$ compresses the amplitude ratio and phase difference of two components, s and p, into one scalar, and it does so on the assumption that the light is always in an ideal, fully polarized state. In a real measurement there is always some degree of depolarization — part of the polarization state randomized by surface roughness on the sample, multiple scattering, or non-ideal elements inside the optics. A scalar $\rho$ has no way of expressing that situation at all.

This post compares the two mathematical languages used to describe a polarization state, the Jones vector and the Stokes vector, examines the difference in what each can express, and explains why a real ellipsometer standardizes on Mueller matrices rather than Jones matrices. It then maps the Mueller matrices of the two core components — the polarizer and the retarder — onto the actual hardware one at a time.

Before that, a picture of where those components sit. An ellipsometer consists of a Polarization State Generator (PSG), which puts the light from the source into a chosen polarization state; the sample it reflects from; and a Polarization State Analyzer (PSA), which reads out how that state has changed.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig3-psg-psa-diagram.png" alt="The PSG-sample-PSA chain as a product of Mueller matrices" width="700">
_Fig 1. The PSG-sample-PSA chain as a sequential product of Mueller matrices_

The PSG and PSA are each built from a polarizer and a retarder. The polarization state changes as the light passes through these elements from left to right, and the whole change equals the product of each element's Mueller matrix taken in order. The aim of this post is to understand what $S_{in}$, $M_P$, $M_C$, $M_{sample}$ and $S_{out}$ in that diagram each are, and why they take the forms they do.

## 1. Defining a polarization state — the trajectory of the field

From the expression in post 1, the electric field of a plane wave travelling along z decomposes into oscillations along x and y.

$$ E(z,t) = \lvert E_x \rvert \cos(\omega t - kz + \phi_{e,x})\, \hat{x} + \lvert E_y \rvert \cos(\omega t - kz + \phi_{e,y})\, \hat{y} $$

At a fixed position — $z=0$, say — the trajectory traced by $(E_x, E_y)$ over time is precisely what defines the polarization state. Only two numbers shape that trajectory: the amplitude ratio $\lvert E_y \rvert / \lvert E_x \rvert$ and the phase difference $\phi_{e,y}-\phi_{e,x}$. A phase difference of 0 or $\pi$ gives a straight line through the origin regardless of the amplitude ratio — linear polarization — with the amplitude ratio setting only the slope of that line. Equal amplitudes ($\lvert E_x\rvert = \lvert E_y \rvert$) together with a phase difference of $\pi/2$ or $3\pi/2$ give a circle: circular polarization. Every other combination gives an ellipse, that is, elliptical polarization. Linear and circular are in fact two limiting special cases of the elliptical.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig1-polarization-trajectories.png" alt="Field trajectories for linear, circular and elliptical polarization" width="700">
_Fig 2. Trajectories of $(E_x, E_y)$ by amplitude ratio and phase difference_

All three trajectories in Figure 2 come from the same expression with only those two numbers changed. Describing a polarization state therefore reduces to the question of what mathematical object should represent this trajectory.

## 2. The Jones vector — the language of full polarization

The most direct approach is to collect the amplitudes and phases that produced the trajectory into a vector. That is the Jones vector.

$$ \mathbf{E} = \begin{bmatrix} \lvert E_x \rvert\, e^{j\phi_{e,x}} \\ \lvert E_y \rvert\, e^{j\phi_{e,y}} \end{bmatrix} $$

The Jones vector expresses an entire polarization state with two complex numbers — four reals: two amplitudes and two phases, of which the overall phase is irrelevant to measurement, leaving three effective degrees of freedom. That compactness also simplifies propagation through an optical system: represent each element by a 2×2 Jones matrix, and the state after several elements is just the ordered product of those matrices (Jones calculus).

The difficulty is the premise on which the representation rests. The two components of a Jones vector presuppose perfectly coherent oscillation holding a fixed phase relation — in other words, that the single exact trajectory of Figure 2 persists without wavering. Real light maintains that relation only within a finite coherence time; measured over any longer interval, the phase difference and amplitude ratio typically fluctuate slightly. The larger those fluctuations, the further the state drifts from "one fully polarized ellipse," and a Jones vector has nowhere to put that drift. An unpolarized state — natural light, in which the phase relation between $E_x$ and $E_y$ is entirely random — cannot be written as a Jones vector at all.

## 3. The Stokes vector — polarization defined by measurable intensities

The Stokes vector takes a different approach. Instead of handling the field directly, it defines the polarization state through combinations of the intensities a detector measures behind four polarization filters.

$$ S_0 = \lvert E_x \rvert^2 + \lvert E_y \rvert^2, \qquad S_1 = \lvert E_x \rvert^2 - \lvert E_y \rvert^2 $$

$$ S_2 = 2\lvert E_x \rvert \lvert E_y \rvert \cos(\phi_{e,x}-\phi_{e,y}), \qquad S_3 = 2\lvert E_x \rvert \lvert E_y \rvert \sin(\phi_{e,x}-\phi_{e,y}) $$

$S_0$ is the total intensity, $S_1$ the intensity difference between horizontal (0°) and vertical (90°) components, $S_2$ that between $+45°$ and $135°$, and $S_3$ that between left- and right-circular components. What matters is that all four are combinations of "intensity transmitted through a particular polarization filter" — quantities a detector measures directly. The components of a Jones vector (the absolute phase of the field) cannot be measured directly, whereas the components of a Stokes vector are in principle determined by four intensity measurements.

The more consequential difference is the range of states each can express. For a fully polarized state $S_0^2 = S_1^2+S_2^2+S_3^2$ holds exactly, but real light generally satisfies only the inequality:

$$ S_0^2 \ge S_1^2 + S_2^2 + S_3^2 $$

The slack in that inequality is exactly the territory a Jones vector cannot reach. Defining the degree of polarization $P = \sqrt{S_1^2+S_2^2+S_3^2}/S_0$, $P=1$ is fully polarized, $0<P<1$ partially polarized, and $P=0$ unpolarized. Figure 3 renders all three in the same Stokes language.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig2-partial-polarization.png" alt="Ensembles of field trajectories for full, partial and zero polarization with the measured degree of polarization" width="700">
_Fig 3. Fully polarized to partially polarized to unpolarized: a Jones vector covers only the left panel, a Stokes vector all three_

Each panel can be read as the result of observing, repeatedly and over a period longer than the coherence time, a state that is fully polarized at any instant — a single ellipse like the one in Figure 2. On the left the same ellipse recurs exactly every time (P=1.00, the ellipse of Figure 2). On the right a completely random ellipse appears each time, so that averaging leaves no preferred direction at all (P=0.02). The middle sits between them, with signal and noise mixed half and half (P=0.50). A Jones vector can express only the "single fixed frame" of the left panel. Measuring a rough surface or a thick film with a real ellipsometer routinely produces something closer to the middle panel, and Jones calculus is then a tool that cannot be applied in principle. This is why ellipsometry takes the Stokes vector as its basic representation.

## 4. The Mueller matrix — the linear transformation through a system

Passing through an optical system changes the polarization state. Expressing that change on Stokes vectors requires a 4×4 real matrix: the Mueller matrix.

$$ S_{out} = M\, S_{in} $$

The entire PSG-sample-PSA chain of Figure 1 equals the product of the individual Mueller matrices taken in order from the incident side ($S_{out} = M_A M_C' M_{sample} M_C M_P S_{in}$). For a system of ideal elements only, Jones matrices could in fact do the same arithmetic. The real reason for Mueller matrices is the depolarization seen in section 3. A sample with surface roughness, for example, scatters part of the incident fully polarized light in random directions and reflects a partially polarized state — a transformation no 2×2 Jones matrix can reproduce, since the output is not expressible as a Jones vector to begin with. A Mueller matrix, by contrast, is a linear map taking any physically valid Stokes vector to another, fully polarized or partially polarized or unpolarized alike, so every real optical phenomenon including depolarization fits in one framework.

## 5. The rotation matrix — changing a component's alignment angle

Before writing the Mueller matrices of a polarizer and a retarder, the rotation matrix has to be defined. It expresses, as a rotation of the reference axes of the Stokes vector, exactly the operation of turning a polarizer by hand on the bench to change its alignment angle.

$$ M_R(\Omega) = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos 2\Omega & \sin 2\Omega & 0 \\ 0 & -\sin 2\Omega & \cos 2\Omega & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix} $$

The angle enters as $2\Omega$, not $\Omega$. This follows from the geometric symmetry of polarization: a linear state returns to itself after a rotation of only $180°$ in azimuth ($E$ and $-E$ are the same line), so the Stokes parameters $(S_1,S_2)$ vary with a period of twice the physical azimuth. Aligning a polarizer or retarder at an arbitrary angle $\Omega$ below always appears as sandwiching its reference-orientation Mueller matrix: $M_R(-\Omega)\,(\cdot)\,M_R(\Omega)$.

## 6. Polarizer and retarder — the expressions and the hardware they correspond to

An ideal polarizer transmits the polarization component along one axis and blocks the rest. At an alignment angle $\Omega_P$ its Mueller matrix is:

$$ M_P(\Omega_P) = M_R(-\Omega_P)\, \frac{1}{2}\begin{bmatrix} 1 & 1 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}\, M_R(\Omega_P) $$

Taking apart what the central matrix (reference orientation, $\Omega_P=0$) does: only $S_1$ survives alongside $S_0$, while $S_2$ and $S_3$ are driven to zero — it transmits the horizontal component and erases every other piece of polarization information. In the laboratory a wire-grid polarizer or a Glan-Taylor prism plays this role. The angle at which the part is turned in its mount is $\Omega_P$ itself.

A retarder, rather than blocking a component, delays the phase of the component along one axis. For an alignment angle $\Omega_C$ and retardance $\phi$:

$$ M_C(\Omega_C, \phi) = M_R(-\Omega_C) \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & \cos\phi & \sin\phi \\ 0 & 0 & -\sin\phi & \cos\phi \end{bmatrix} M_R(\Omega_C) $$

The central matrix in the reference orientation leaves $S_1$ alone and rotates the $(S_2, S_3)$ plane by $\phi$. As the simulation below shows directly, that rotation is the mechanism turning linear polarization into circular. Physically this is a wave plate, a Babinet-Soleil compensator, or a liquid crystal variable retarder.

The retardance $\phi$ is set by the birefringence $\Delta n$ of the element, its physical thickness $d$, and the wavenumber $\sigma=1/\lambda$.

$$ \phi(\sigma) = 2\pi \sigma \Delta n\, d $$

This shows that $\phi$ depends on wavelength. For a monochromatic ellipsometer it hardly matters, but in spectroscopic ellipsometry, which handles a broad band at once, it means each wavelength component emerges from the retarder in a different polarization state. Turning that property from an error source into a tool for wavelength discrimination is the central idea behind channelled spectroscopic ellipsometry, but that story lies beyond this background series and is noted here only in passing.

## 7. Simulation — how polarization changes with retardance

Multiplying the expressions above out directly: an unpolarized source ($S_{in}=[1,0,0,0]^T$) passes a polarizer aligned at $0°$ and then a retarder aligned at $45°$. Sweeping only the retardance $\phi$ from 0 to $\pi$, the output Stokes vector $S_{out} = M_C(45°,\phi)\, M_P(0°)\, S_{in}$ was computed.

```python
# core of generate_figures.py (full code: _code/ellipsometry-polarization-mueller-matrix/)
def mueller_rotation(omega):
    c, s = np.cos(2 * omega), np.sin(2 * omega)
    return np.array([[1, 0, 0, 0], [0, c, s, 0], [0, -s, c, 0], [0, 0, 0, 1]])

def mueller_retarder(omega_c, phi):
    c, s = np.cos(phi), np.sin(phi)
    base = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, c, s], [0, 0, -s, c]])
    return mueller_rotation(-omega_c) @ base @ mueller_rotation(omega_c)

S_out = mueller_retarder(np.pi/4, phi) @ mueller_polarizer(0.0) @ S_in
```

Carrying out the matrix product gives a closed form for the output.

$$ S_{out} = \tfrac{1}{2}\begin{bmatrix} 1 & \cos\phi & 0 & \sin\phi \end{bmatrix}^T $$

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig4-retarder-phase-sweep.png" alt="Output polarization against retardance after a polarizer and a retarder" width="750">
_Fig 4. Output polarization state against retardance $\phi$, after a polarizer (0°) and a retarder (45°)_

The five ellipses along the top show directly how the state evolves as $\phi$ grows. At $\phi=0$ the retarder does nothing and the horizontal linear polarization from the polarizer emerges unchanged. As $\phi$ increases the trajectory opens out, becoming a full circle at $\phi=\pi/2$ — exactly a quarter wave — the standard result when the retarder is aligned $45°$ from the axis of the incident linear polarization. Increasing $\phi$ further narrows the circle back into an ellipse, until at $\phi=\pi$ (half wave) it has flipped entirely into linear polarization perpendicular to the original.

The plot along the bottom shows the same result as the Stokes components $S_1/S_0$, $S_2/S_0$ and $S_3/S_0$. As the closed form says, $S_2$ stays exactly zero for every $\phi$. Setting the alignment to $45°$ makes the retarder's Mueller matrix hold the $S_2$ axis fixed while rotating only the $(S_1, S_3)$ plane by $\phi$, and the input leaving the polarizer already has $S_2 = 0$. What moves instead are $S_1$ and $S_3$, trading places as $\cos\phi$ and $\sin\phi$ — all reproduced by a single Mueller matrix product. That two components alone, with three numbers between them — two alignment angles and a retardance — can move freely among linear, elliptical and circular states is exactly why these two make up the core of a PSG and PSA.

## Summary and what comes next

This post examined the difference in expressive power between the two languages for polarization. The Jones vector can represent only fully polarized states but keeps the arithmetic compact; the Stokes vector, defined through measurable intensities, extends to partial polarization and to unpolarized light. Real ellipsometers standardize on Mueller rather than Jones calculus because situations that break the fully-polarized assumption — surface roughness, depolarization — genuinely occur. On that foundation, the Mueller matrices of the polarizer and the retarder were defined, with the rotation matrix expressing alignment angle and the retardance formula expressing wavelength dependence.

The next post applies this Mueller formalism to reflection from a thin-film sample. The ellipsometric parameters $\Psi$ and $\Delta$, which emerge when the Fresnel coefficients derived in post 1 combine with multiple reflection inside a film, are treated there in terms of how they are measured and interpreted in the Stokes-Mueller language defined here.

## References

- Youngjoon Kim, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발" [Development of snapshot angle-resolved ellipsometry using a line-scan spectrometer and back-focal-plane spectral interference], Ph.D. dissertation, Seoul National University, 2025 (in Korean), sections 1.2 and 2.3.
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 2 (Jones and Stokes vectors, Mueller matrices).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1987, ch. 1–2 (mathematical representations of polarization).
- E. Collett, *Field Guide to Polarization*, SPIE Press, 2005 (Stokes vectors, degree of polarization, Mueller calculus).

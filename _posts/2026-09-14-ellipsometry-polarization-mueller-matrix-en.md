---
title: "Ellipsometry Foundations 2 — How Many Numbers Does It Take to Describe the Polarization of Light?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometry-polarization-mueller-matrix/
page_id: ellipsometry-polarization-mueller-matrix
date: 2026-09-14 20:00:00 +0900
categories: [Optics, Polarization]
tags: [ellipsometry, polarization, jones-vector, stokes-vector, mueller-matrix]
description: What a Jones vector can and cannot express compared with a Stokes vector, and the Mueller matrices of polarizers and retarders mapped onto the hardware they describe.
math: true
---

[Post 1](/en/posts/ellipsometry-electromagnetic-fresnel/) ended by arriving at the ratio of Fresnel reflection coefficients, $\rho = r_p/r_s = \tan\Psi\, e^{j\Delta}$. Taken alone, that single value looks like everything ellipsometry needs. But $\rho$ compresses the amplitude ratio and phase difference of two components, s and p, into one scalar, and it rests on the assumption that the light is always in an ideal, fully polarized state. In a real measurement, surface roughness on the sample, multiple scattering and non-ideal optical elements usually randomize part of the polarization state. This is called depolarization. A scalar $\rho$ has no way of expressing it.

This post compares the two mathematical languages used to describe a polarization state, the Jones vector and the Stokes vector, and examines the difference in what each can express. It then explains why Mueller matrices rather than Jones matrices are needed once depolarization has to be handled. Finally, it maps the Mueller matrices of the two core components of an ellipsometer, the polarizer and the retarder, onto the actual hardware one at a time.

First, a picture of where those components sit. An ellipsometer consists of a Polarization State Generator (PSG), which puts the light from the source into a chosen polarization state; the sample it reflects from; and a Polarization State Analyzer (PSA), which reads out how that state has changed.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig3-psg-psa-diagram.png" alt="The PSG-sample-PSA chain as a product of Mueller matrices" width="700">
_Fig 1. The PSG-sample-PSA structure of an ellipsometer. Light travels from left to right, but Mueller matrices multiply from the right, so the element on the incident side, $$M_P$$, stands rightmost in the expression._

The PSG and PSA are each built from a polarizer and a retarder. The polarization state changes at each element the light passes, and the whole change equals the product of the elements' Mueller matrices taken in order. The aim of this post is to understand what $S_{in}$, $M_P$, $M_C$, $M_{sample}$ and $S_{out}$ in that diagram each are, and why they take the forms they do. The $$M_C'$$ and $M_A$ on the PSA side are the same kinds of elements as $M_C$ and $M_P$ in the PSG.

## 1. What determines a polarization state?

What exactly does it mean for light to be "polarized"? From the expression in post 1, the electric field of a plane wave travelling along z decomposes into oscillations along x and y.

$$ E(z,t) = \lvert E_x \rvert \cos(\omega t - kz + \phi_{e,x})\, \hat{x} + \lvert E_y \rvert \cos(\omega t - kz + \phi_{e,y})\, \hat{y} $$

At a fixed position, $z=0$ say, the trajectory traced by $(E_x, E_y)$ over time is the polarization state. Only two numbers shape that trajectory: the amplitude ratio $\lvert E_y \rvert / \lvert E_x \rvert$ and the phase difference $\phi_{e,y}-\phi_{e,x}$. A phase difference of 0 or $\pi$ gives a straight line through the origin regardless of the amplitude ratio. That is linear polarization, and the amplitude ratio sets only the slope of the line. If one amplitude is zero, the light is also linearly polarized whatever the phase difference. Equal amplitudes ($\lvert E_x\rvert = \lvert E_y \rvert$) with a phase difference of $\pi/2$ or $3\pi/2$ give a circle: circular polarization. Every other combination gives an ellipse, that is, elliptical polarization. Linear and circular polarization are the two limiting special cases of the elliptical.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig1-polarization-trajectories.png" alt="Field trajectories for linear, circular and elliptical polarization" width="700">
_Fig 2. Three trajectories from the same expression with only the amplitude ratio and phase difference changed. From left, (amplitude ratio, phase difference) is (1, 0), (1, π/2) and (0.55, π/4). The black dot marks $t=0$ and the red arrow the sense of rotation. The circle in the middle turns clockwise as seen looking into the source, so it is right-circular._

All three trajectories in Figure 2 come from the same expression with only those two numbers changed. Post 1's $\Psi$ and $\Delta$ are exactly these two numbers. With light linearly polarized at 45° on the sample, the p and s components start with equal amplitude and phase, so the reflected ellipse has amplitude ratio $\lvert E_p \rvert/\lvert E_s \rvert = \tan\Psi$ and phase difference $\phi_p-\phi_s = \Delta$. In short, an ellipsometer measures the shape of the reflected ellipse. Describing a polarization state therefore reduces to the question of what mathematical object should represent this trajectory.

Before going further, the sign conventions of this post are collected in one place. Textbooks sometimes take the handedness of circular light and the sense of angles in opposite ways. Without a stated convention, the same matrix gives answers of opposite sign.

| Item | Convention in this post |
| --- | --- |
| Time factor | $e^{j\omega t}$, as in post 1. The physical field has the form $\cos(\omega t - kz + \phi)$ above |
| Phase difference | $\delta = \phi_{e,y}-\phi_{e,x}$ |
| Viewing direction | Looking into the source: x points right, y points up, and the light comes toward the viewer |
| Handedness | Clockwise rotation in this view is right-circular. $\delta=+\pi/2$ is right-circular and corresponds to $S_3>0$ in section 3 |
| Azimuth | Measured from the x axis toward the y axis, i.e. counterclockwise in this view (section 5) |
| Retarder angle | The azimuth of the fast axis (section 6) |
| Link to p and s | x is the p direction and y the s direction. The sample Mueller matrix in [post 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/) uses the same assignment |

This handedness convention is the one used in optics and ellipsometry textbooks such as Born & Wolf, Azzam & Bashara and Fujiwara. The IEEE convention, which judges handedness looking along the direction of travel, swaps the names right and left.

## 2. What can a Jones vector not describe?

The most direct way to turn the trajectory into a mathematical object is to collect the amplitudes and phases that produced it into a vector. That is the Jones vector.

$$ \mathbf{E} = \begin{bmatrix} \lvert E_x \rvert\, e^{j\phi_{e,x}} \\ \lvert E_y \rvert\, e^{j\phi_{e,y}} \end{bmatrix} $$

The Jones vector expresses a polarization state with two complex numbers, four reals in all: two amplitudes and two phases. The overall phase does not show up in measurement, so three degrees of freedom remain. That compactness also simplifies propagation through an optical system. Represent each element by a 2×2 Jones matrix, and the state after several elements is just the ordered product of those matrices (Jones calculus).

The difficulty lies in the premise. A Jones vector assumes that the two components hold a fixed phase relation, that is, a perfectly coherent oscillation. It treats the single trajectory of Figure 2 as persisting without wavering. Light from a real source often keeps its amplitude ratio and phase difference fixed only over intervals of the order of the coherence time and lets them wander over longer ones. Conversely, light that keeps them fixed, such as light fresh from a polarizer, is fully polarized however short its coherence time. The larger the wandering, the further the state drifts from "one fully polarized ellipse." A Jones vector has nowhere to put that wandering. An unpolarized state, such as natural light, in which the phase relation between $E_x$ and $E_y$ is entirely random, cannot be written as a Jones vector at all.

## 3. Why is the Stokes vector defined by intensities?

If the wandering phase cannot be stored, what should polarization be written in terms of? The Stokes vector uses the quantity a detector actually measures: intensity. Instead of handling the field directly, it defines the polarization state through combinations of the intensities transmitted by four polarization filters.

$$ S_0 = \lvert E_x \rvert^2 + \lvert E_y \rvert^2, \qquad S_1 = \lvert E_x \rvert^2 - \lvert E_y \rvert^2 $$

$$ S_2 = 2\lvert E_x \rvert \lvert E_y \rvert \cos(\phi_{e,y}-\phi_{e,x}), \qquad S_3 = 2\lvert E_x \rvert \lvert E_y \rvert \sin(\phi_{e,y}-\phi_{e,x}) $$

$S_0$ is the total intensity, $S_1$ the intensity difference between horizontal (0°) and vertical (90°) components, $S_2$ that between $+45°$ and $135°$, and $S_3$ that between right- and left-circular components (handedness as in section 1). What matters is that all four are combinations of "intensity transmitted through a particular polarization filter," quantities a detector measures directly. The components of a Jones vector (the absolute phase of the field) cannot be measured directly. The components of a Stokes vector, in principle, follow from four intensity measurements.

Consider the reflected light of section 1 in this language. With x as p and y as s, and light linearly polarized at 45° on the sample, the normalized Stokes vector of the reflected light is:

$$ \frac{S}{S_0} = \begin{bmatrix} 1 & -\cos2\Psi & \sin2\Psi\cos\Delta & -\sin2\Psi\sin\Delta \end{bmatrix}^T $$

Here $S_1$ depends on $\Psi$ alone, while $S_2$ and $S_3$ depend on both $\Psi$ and $\Delta$. The minus sign on $S_3$ arises because, under this assignment, the phase difference $\phi_{e,y}-\phi_{e,x}$ is $\phi_s-\phi_p = -\Delta$. For example, $\Psi=20°$ and $\Delta=100°$ give $[1, -0.766, -0.112, -0.633]^T$. In short, the $(\Psi, \Delta)$ an ellipsometer measures are carried directly by three Stokes components.

The more consequential difference is the range of states each can express. For a fully polarized state $S_0^2 = S_1^2+S_2^2+S_3^2$ holds exactly, but real light generally satisfies only the inequality:

$$ S_0^2 \ge S_1^2 + S_2^2 + S_3^2 $$

The slack in that inequality is exactly the territory a Jones vector cannot reach. Defining the degree of polarization $P = \sqrt{S_1^2+S_2^2+S_3^2}/S_0$, $P=1$ is fully polarized, $0<P<1$ partially polarized, and $P=0$ unpolarized. Figure 3 tells these three cases apart by the degree of polarization $P$.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/en/fig2-partial-polarization.png" alt="Ensembles of field trajectories for full, partial and zero polarization with the measured degree of polarization" width="700">
_Fig 3. Ensembles of ellipses, each fully polarized at its instant. The share of one fixed ellipse (amplitude ratio 0.55 and phase difference π/4, as on the right of Figure 2) is 1, 0.5 and 0 from left to right, with random ellipses filling the rest. The degree of polarization computed from 2,000 ellipses is 1.00, 0.50 and 0.02; the middle and right panels draw 400 of them. A Jones vector can represent only the left panel._

Each panel can be read as the result of observing, repeatedly and over a period longer than the coherence time, a state that is fully polarized at any instant: a single ellipse like those of Figure 2. On the left the same ellipse recurs exactly every time (P=1.00, the same shape as the ellipse on the right of Figure 2). On the right a completely random ellipse appears each time, so that averaging leaves no preferred direction (P=0.02). The middle sits between them, with signal and noise mixed half and half (P=0.50). A Jones vector can express only the "single fixed frame" of the left panel. Measuring a rough surface, a film of non-uniform thickness, or a transparent substrate with backside reflection yields partially polarized light noticeably below P=1 (around 0.7–0.9 in published examples). Jones calculus cannot describe this case at all. This is why ellipsometry turns to the Stokes vector wherever depolarization has to be handled.

## 4. What matrix handles depolarization?

A language for partial polarization is now in hand. What matrix, then, describes a sample or element that scrambles polarization? A transformation acting on Stokes vectors is written as a 4×4 real matrix. That is the Mueller matrix.

$$ S_{out} = M\, S_{in} $$

Here each element of $M$ says how much of an input Stokes component is carried into an output Stokes component. The simplest example is an ideal depolarizer, $M = \mathrm{diag}(1, a, a, a)$ with $0 \le a \le 1$. It leaves $S_0$ alone and scales $S_1, S_2, S_3$ by $a$, so fully polarized input ($P=1$) leaves as partially polarized light with $P=a$. With $a=0.8$ the result falls within the published 0.7–0.9 quoted in section 3. In short, the "transformation that lowers P," which a Jones matrix cannot write, a Mueller matrix needs only a few diagonal entries to express it.

The entire PSG-sample-PSA chain of Figure 1 equals the product of the individual Mueller matrices, ordered so that the incident-side element stands on the right ($$S_{out} = M_A M_C' M_{sample} M_C M_P S_{in}$$). For a system of ideal elements only, Jones matrices could do the same arithmetic. Mueller matrices become necessary because of the depolarization seen in section 3. A sample whose surface is rough on the scale of the wavelength, for example, scatters part of the incident fully polarized light in random directions and reflects a partially polarized state. No 2×2 Jones matrix can reproduce that transformation, because the output cannot be written as a Jones vector. Smaller roughness hardly depolarizes and is often modelled as an effective medium approximation (EMA) layer.

A Mueller matrix, by contrast, is a linear map taking any physically valid Stokes vector, fully polarized, partially polarized or unpolarized, to another. Polarization transformations in general, depolarization included, therefore fit in one framework. Mueller calculus is a linear combination of intensities, though, so calculations that need absolute phase, such as interference of coherent beams, remain the job of the Jones side. In short, Jones is the language of calculation and Mueller the language of measurement.

## 5. How does the matrix change when a component is rotated?

Turn a polarizer in its mount by 30° on the bench: what changes in the equations? There is no need to memorize a separate matrix for every angle. Rotating the reference axes of the Stokes vector does the same job as rotating the part, and that operation is the rotation matrix.

$$ M_R(\Omega) = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & \cos 2\Omega & \sin 2\Omega & 0 \\ 0 & -\sin 2\Omega & \cos 2\Omega & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix} $$

The angle enters as $2\Omega$, not $\Omega$. The reason is that a linear state returns to itself after only $180°$ of azimuth ($E$ and $-E$ are the same line). The Stokes components $(S_1,S_2)$ therefore turn through twice the azimuth, $2\Omega$, and $180°$ of azimuth makes one full period. Below, placing a polarizer or retarder at an arbitrary azimuth $\Omega$ always appears as sandwiching its reference-orientation Mueller matrix: $M_R(-\Omega)\,(\cdot)\,M_R(\Omega)$.

One number confirms the sense of rotation. Set the polarizer of section 6 at $\Omega_P=45°$ and send in unpolarized light: the output is $\tfrac{1}{2}[1, 0, 1, 0]^T$, linear polarization at $+45°$ with $S_2>0$. The azimuth is thus measured from the x axis toward the y axis, as the convention of section 1 states.

## 6. What do the polarizer and retarder matrices correspond to in hardware?

Now the two components of the PSG and PSA can be written as matrices. An ideal polarizer transmits the polarization component along one axis and blocks the rest. At azimuth $\Omega_P$ its Mueller matrix is:

$$ M_P(\Omega_P) = M_R(-\Omega_P)\, \frac{1}{2}\begin{bmatrix} 1 & 1 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}\, M_R(\Omega_P) $$

Taking apart what the central matrix (reference orientation, $\Omega_P=0$) does: only $S_1$ survives alongside $S_0$, while $S_2$ and $S_3$ become zero. In other words, it transmits the horizontal component and erases all other polarization information. The leading factor $\tfrac{1}{2}$ means that unpolarized input loses half its intensity. Send in linear polarization at azimuth $\theta$, $[1, \cos2\theta, \sin2\theta, 0]^T$, and the output intensity is $$S_0' = (1+\cos2\theta)/2 = \cos^2\theta$$: 0.75 at $\theta=30°$, 0.25 at $60°$, and 0 at $90°$. In short, this matrix contains Malus's law. In the laboratory a wire-grid polarizer or a Glan-Taylor prism plays this role. The angle at which the part is turned in its mount is $\Omega_P$ itself.

A retarder, rather than blocking a component, introduces a phase difference between the components along its fast and slow axes. For a fast axis at azimuth $\Omega_C$ and retardance $\phi$:

$$ M_C(\Omega_C, \phi) = M_R(-\Omega_C) \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & \cos\phi & \sin\phi \\ 0 & 0 & -\sin\phi & \cos\phi \end{bmatrix} M_R(\Omega_C) $$

The central matrix in the reference orientation leaves $S_1$ alone and rotates the $(S_2, S_3)$ plane by $\phi$. For example, $+45°$ linear polarization $[1, 0, 1, 0]^T$ given $\phi=\pi/2$ becomes $[1, 0, 0, -1]^T$, left-circular light. The slow-axis (y) component lags the fast-axis (x) component, so the phase difference $\delta$ becomes $-\pi/2$. As the simulation in section 7 shows, this rotation is the very mechanism that turns linear polarization into circular. Physical realizations include a wave plate, a Babinet-Soleil compensator, or a liquid crystal variable retarder. Figure 4 shows a wave plate in a rotation mount.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/ext-waveplate-rotation-mount.jpg" alt="A one-inch wave plate in a rotation mount" width="450">
_Fig 4. A one-inch wave plate in a rotation mount. The degree scale around the rim reads how far the plate has been turned, and that turn changes the fast-axis azimuth $$\Omega_C$$. (Source: Geek3, [Optical-waveplate-1inch](https://commons.wikimedia.org/wiki/File:Optical-waveplate-1inch.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

What, then, sets the retardance $\phi$? Retarders are made of birefringent material, whose refractive index depends on the direction of polarization. The letter under the calcite crystal in Figure 5 appears twice for the same reason: two mutually perpendicular polarizations see different indices and take different paths.

<img src="/assets/img/posts/ellipsometry-polarization-mueller-matrix/ext-calcite-double-refraction.jpg" alt="The letter A seen doubled through a calcite crystal" width="450">
_Fig 5. The letter "A" seen through a calcite crystal. Inside the crystal, two mutually perpendicular polarizations see different refractive indices and take different paths, so the letter appears doubled. (Source: Twyla Baker, [Calcite and double refraction (GeoDIL number - 670)](https://commons.wikimedia.org/wiki/File:Calcite_and_double_refraction_%28GeoDIL_number_-_670%29.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

A wave plate is cut with the optic axis of the crystal parallel to its face, so the two polarizations of normally incident light travel the same path. The image then does not split; one polarization simply lags. That lag is the retardance, set by the birefringence $\Delta n$ (the difference between the two indices), the physical thickness $d$, and the wavenumber $\sigma=1/\lambda$.

$$ \phi(\sigma) = 2\pi \sigma \Delta n\, d $$

What matters here is that $\phi$ is proportional to $1/\lambda$. In short, a retarder's retardance differs from wavelength to wavelength. Even ignoring the wavelength dependence of $\Delta n$, a plate that is exactly a quarter wave ($90°$) at 633 nm gives $142°$ at 400 nm, $114°$ at 500 nm and $71°$ at 800 nm. In practice $\Delta n$ also varies with wavelength, so the values differ somewhat. For a monochromatic ellipsometer this hardly matters. In spectroscopic ellipsometry, which handles a broad band at once, each wavelength component leaves the retarder in a different polarization state. This is why a spectroscopic ellipsometer must calibrate the retarder's retardance wavelength by wavelength. Using that property as a tool for wavelength discrimination rather than as an error source is the idea of channelled spectroscopic ellipsometry, covered in part 5 of the ellipsometer series (channelled spectroscopy).

## 7. How does the polarization change with retardance?

How does the polarization state move when these matrices are actually multiplied? An unpolarized source ($S_{in}=[1,0,0,0]^T$) passes a polarizer at azimuth $0°$ and then a retarder with its fast axis at $45°$. Sweeping only the retardance $\phi$ from 0 to $\pi$, the output Stokes vector $S_{out} = M_C(45°,\phi)\, M_P(0°)\, S_{in}$ was computed.

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
_Fig 6. Polarization state of light after a polarizer (0°) and a retarder (fast axis at 45°). The top row shows the trajectories at $\phi$ = 0, π/4, π/2, 3π/4 and π (red arrows give the sense of rotation); the bottom row shows the normalized Stokes components (dashed lines mark the $\phi$ of the top row). $$S_2$$ stays at zero; the light is right-circular at $\phi=\pi/2$ ($$S_3/S_0=1$$) and vertically linear at $\phi=\pi$._

The five trajectories along the top show how the state evolves as $\phi$ grows. At $\phi=0$ the retarder does nothing, and the horizontal linear polarization from the polarizer emerges unchanged. As $\phi$ increases the trajectory opens out and becomes a full circle at $\phi=\pi/2$ (a quarter wave). This is the standard result when the retarder is set $45°$ from the axis of the incident linear polarization. The arrow turns clockwise, so the circle is right-circular, in agreement with $S_3/S_0=1$ in the bottom row. Increasing $\phi$ further narrows the circle back into an ellipse, until at $\phi=\pi$ (a half wave) it has flipped into linear polarization perpendicular to the original.

The plot along the bottom shows the same result as the Stokes components $S_1/S_0$, $S_2/S_0$ and $S_3/S_0$. As the closed form says, $S_2$ stays exactly zero for every $\phi$, for two reasons. With the azimuth at $45°$, the retarder's Mueller matrix holds the $S_2$ axis fixed and rotates only the $(S_1, S_3)$ plane by $\phi$. And the input leaving the polarizer already has $S_2 = 0$. What moves instead are $S_1$ and $S_3$, trading places as $\cos\phi$ and $\sin\phi$. Two components (a polarizer and a retarder) and three numbers (two azimuths and one retardance) are enough to produce linear, elliptical and circular polarization. That is why the PSG and PSA are built from these two components.

## Summary

This post examined the difference in expressive power between the two languages for polarization. The Jones vector can represent only fully polarized states but keeps the arithmetic compact. The Stokes vector, defined through measurable intensities, extends to partial polarization and to unpolarized light. Ellipsometer equations are usually set up in Jones calculus, but situations that break the fully polarized assumption, such as surface roughness and depolarization, genuinely occur, and there Mueller calculus becomes necessary. On that foundation, the Mueller matrices of the polarizer and the retarder were defined, with the rotation matrix expressing azimuth and the retardance formula expressing wavelength dependence. The polarizer matrix contains Malus's law, and the retarder matrix a rotation of the $(S_2, S_3)$ plane. Signs follow the convention of section 1: light that turns clockwise as seen looking into the source is right-circular, with $S_3>0$.

[Post 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/) computes the reflection coefficients $r_p, r_s$ of thin-film multilayers. It treats how the Fresnel coefficients of post 1 combine with multiple reflection inside a film, and uses the resulting $r_p, r_s$ to build the sample Mueller matrix $M_{sample}$ at the centre of this post's chain (post 3, section 7).

## References

- Youngjoon Kim, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발" [Development of snapshot angle-resolved ellipsometry using a line-scan spectrometer and back-focal-plane spectral interference], Ph.D. dissertation, Seoul National University, 2025 (in Korean), sections 1.2 and 2.3.
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 3 (Jones and Stokes vectors, Mueller matrices).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, ch. 1–2 (mathematical representations of polarization).
- E. Collett, *Field Guide to Polarization*, SPIE Press, 2005 (Stokes vectors, degree of polarization, Mueller calculus).

---
layout: post
date: 2026-12-13 20:00:00 +0900
title: "Wave Optics Foundations 5 — Why Does Laser Light Seem Not to Spread, Even Far Away?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/wave-optics-gaussian-beams/
page_id: wave-optics-gaussian-beams
categories: [Optics, Wave Optics]
tags: [wave-optics, gaussian-beam, laser-resonator, gouy-phase, beam-quality]
description: "How thin and how far a laser beam can go, followed through the Gaussian beam's waist and Rayleigh range, the Gouy phase, the ABCD law, resonator modes and M²."
math: true
---

Point a laser pointer at a wall and the light is still a small dot on a wall 10 m away. A flashlight beam spreads widely after only a few steps. So a laser beam looks like light that does not spread. Yet send a laser beam to the Moon, as in Figure 1, and depending on how thick the beam was to begin with, its diameter at the Moon grows to anywhere from hundreds of metres to thousands of kilometres.

<img src="/assets/img/posts/wave-optics-gaussian-beams/ext-lunar-laser-ranging.jpg" alt="Green laser beams sent into the night sky from NASA Goddard's laser ranging facility" width="640">
_Figure 1. The laser ranging facility at NASA's Goddard Space Flight Center. Two green beams head into the sky; both aim at the Lunar Reconnaissance Orbiter, a spacecraft orbiting the Moon. (Source: NASA, [Goddard Spaceflight Center Laser Ranging Facility](https://commons.wikimedia.org/wiki/File:Goddard_Spaceflight_Center_Laser_Ranging_Facility.jpg), public domain)_

[Part 3](/en/posts/wave-optics-diffraction/) showed that light passing an aperture of width $D$ spreads by an angle of about $\lambda/D$. Laser beams are no exception. What sets them apart is that they take the shape that spreads least within that limit. That shape is a Gaussian function, and such a beam is called a Gaussian beam. This post answers the following questions.

- Why does a laser beam at first seem not to spread, and then spread after all? Where is the boundary?
- How small can a lens focus a beam, and how long does the beam keep that size?
- Why is the light from a laser Gaussian in shape? What sets the size of the beam?
- As previewed at the end of [Part 4](/en/posts/wave-optics-fourier-imaging/), why does the phase of light passing through focus shift by π?
- How far does a real laser beam depart from the ideal Gaussian, and how is that measured?

## 1. Which beam shape spreads the least?

What is the difference between a flashlight and a laser pointer? In the language of [Part 2](/en/posts/wave-optics-coherence/), it is spatial coherence. In a flashlight, the many points of an extended source emit independently, so light from many directions is mixed from the start. A laser beam has its phase ordered across the whole cross-section and moves as a single wave. Even such a wave cannot avoid the diffraction of Part 3. So what cross-sectional shape spreads the least?

The answer is the Gaussian. Figure 2 shows a photograph of a commercial green laser pointer beam projected onto a screen, along with its intensity profile.

<img src="/assets/img/posts/wave-optics-gaussian-beams/en/fig2-tem00-profile.png" alt="Photograph of a green laser pointer beam cross-section and a Gaussian fit to its intensity profile" width="780">
_Figure 2. (a) Beam cross-section of a 5 mW green laser pointer. (b) The green channel of the photograph read along the dashed line (green) and a Gaussian fit (black dashes). The standard deviation of the residual is 0.8 % of the peak, so the shape is nearly identical to a Gaussian. The camera's gamma correction raises intensity to a power, which keeps the Gaussian shape but may change the width. (Photo source: Zaereth, [Green laser pointer TEM00 profile](https://commons.wikimedia.org/wiki/File:Green_laser_pointer_TEM00_profile.JPG), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

Where does the Gaussian beam come from? Light that travels almost in one direction, like a beam, obeys the paraxial wave equation. The simplest solution of that equation is a spherical wave spreading from a point (a paraboloidal wave in the paraxial case). Saleh and Teich use a complex source point here [1]. Shifting the source position $z$ by an imaginary amount still gives a solution. But placing the source at an imaginary position turns the quadratic phase $e^{ikr^2/2z}$ of the paraboloidal wave into a Gaussian amplitude $e^{-r^2/W^2}$ across the cross-section. In other words, a Gaussian beam is "a spherical wave whose source sits at an imaginary distance". All the formulas that follow come from this picture. The field at a distance $z$ from the waist (the narrowest point) is [1, 2]

$$
U(r, z) = \frac{W_0}{W(z)}\exp\!\left[-\frac{r^2}{W(z)^2}\right]\exp\!\left[i\left(kz + \frac{kr^2}{2R(z)} - \zeta(z)\right)\right]
$$

Here $W_0$ is the waist radius, $W(z)$ is the radius at which the intensity falls to $1/e^2$ of the centre (the beam radius), $R(z)$ is the radius of curvature of the wavefront, and $\zeta(z)$ is the Gouy phase discussed in Section 2. All three are set by a single distance $z_0$:

$$
W(z) = W_0\sqrt{1 + \left(\frac{z}{z_0}\right)^2},\qquad R(z) = z\left[1 + \left(\frac{z_0}{z}\right)^2\right],\qquad z_0 = \frac{\pi W_0^2}{\lambda}
$$

$z_0$ is called the Rayleigh range. In other words, out to $z_0$ from the waist the beam radius grows only by a factor of $\sqrt{2}$; beyond that it grows in proportion to distance. $z_0$ is exactly the boundary where the beam goes from seeming not to spread to spreading. Far away, the beam spreads along a cone of half-angle

$$
\theta_0 = \frac{\lambda}{\pi W_0}
$$

[1]. Since $W_0\theta_0 = \lambda/\pi$, narrowing the waist increases the spread by the same factor. This product is a lower bound that cannot be reduced. As Section 5 shows, real beams have a larger product.

<img src="/assets/img/posts/wave-optics-gaussian-beams/en/fig3-gaussian-beam-shape.png" alt="Side view of a Gaussian beam's intensity, and how the beam width and wavefront curvature change with distance" width="780">
_Figure 3. (a) Side view of the intensity of a Gaussian beam with waist radius $$W_0 = 10\lambda$$, propagated by the angular spectrum method (Part 3). The cyan curve is the formula's $W(z)$, the white dots are the far-field cone, and the vertical dashed lines are $$\pm z_0$$. (b) The beam radius (blue, $$W/W_0$$) and the wavefront radius of curvature (red, plotted as $$R/4z_0$$ so that its minimum is 0.5). Circles are the radius measured from wave propagation. The radius of curvature is infinite at the waist (a flat wavefront), reaches its minimum of $$2z_0$$ at $$\pm z_0$$, and far away becomes equal to the distance $z$._

Let us read Figure 3. Near the waist the beam looks like a parallel bundle and the wavefront is nearly flat. Past $z_0$ the beam approaches the dotted cone, and the wavefront becomes like a sphere centred on the waist. In other words, a Gaussian beam behaves like a plane wave up close and like a spherical wave far away, and $z_0$ connects the two. The beam radius measured from wave propagation agrees with the formula to within 0.05 %, and so does the radius of curvature. The Gouy phase of the next section agrees to within 0.0002 rad.

Now let us plug in numbers. A 532 nm green beam 1 mm in diameter (beam radius 0.5 mm) has $z_0$ = 1.48 m. On a wall 10 m away its diameter is 6.8 mm, still a small dot. The same beam reaching the Moon (380,000 km) widens to 260 km across. It seems not to spread over the first few metres because $z_0$ is 1.5 m. In Saleh and Teich's example, a He-Ne laser (633 nm) beam with a 0.1 mm waist diameter has a Rayleigh range of only 12 mm and a full divergence of 8 mrad [1]. When this beam reaches the Moon, its diameter exceeds 3000 km. Conversely, a wider beam spreads less. Expand the green beam to 1 m in diameter before firing, and diffraction alone would shrink its diameter at the Moon to 260 m. In practice, turbulence in the atmosphere spreads the beam much more than this. That is why a beam meant to travel far is expanded with a telescope before it is fired (the beam expander of Section 3).

Why does the Gaussian shape spread the least? Cutting light with an aperture of the same diameter produces the Airy pattern of Part 3, which spreads by a full angle of $2.44\lambda/D$. A Gaussian beam with waist diameter $D = 2W_0$ spreads by a full angle of $1.27\lambda/D$, a little over half as much [2]. Sharp edges create high spatial frequencies, which scatter to wide angles. A Gaussian fades out smoothly at its edges and has no such high frequencies. The Fourier transform of a Gaussian is again a Gaussian. The product of spatial width and spread angle mentioned in Part 3 is smallest for the Gaussian when measured by second moments. Inside the beam radius $W$ lies 86 % of the beam power [1].

## 2. Why does the phase of light lag as it passes through focus?

Why does the phase of light passing through focus drift away from a plane wave's? Besides $kz$ and the curvature term, the formula above contains one more phase:

$$
\zeta(z) = \tan^{-1}\frac{z}{z_0}
$$

Here $z_0$ is the Rayleigh range. $\zeta$ is called the Gouy phase. On the axis, the phase of a Gaussian beam lags that of a plane wave of the same frequency by $\zeta$. $\zeta$ goes from $-\pi/2$ far before the waist, through 0 at the waist, to $+\pi/2$ far after the waist. In other words, passing through focus the beam falls behind the plane wave by a total of $\pi$ across the waist [1]. L. G. Gouy first showed this in 1890 [4]. The "phase anomaly" near focus introduced in Part 3, Section 5, following Born and Wolf, is exactly this.

Why does it fall behind? Confining a beam sideways mixes in oblique components in the plane-wave decomposition of Part 3, Section 7. An obliquely travelling plane wave accumulates phase along the axis with wavenumber $k_z = \sqrt{k^2 - k_x^2 - k_y^2}$, more slowly than $k$. So the average phase of the confined beam falls behind a plane wave travelling straight along the axis. Near the waist the beam is narrowest, so the phase falls behind fastest; far away the beam is wide and its phase hardly falls any further behind [1, 5].

<img src="/assets/img/posts/wave-optics-gaussian-beams/en/fig4-gouy-phase.png" alt="Gouy phase of a Gaussian beam and higher-order modes, and the resonance frequencies of a resonator" width="780">
_Figure 4. (a) How far the on-axis phase lags behind a plane wave. Lines are the formula $$(N+1)\tan^{-1}(z/z_0)$$ and circles are values measured from wave propagation. HG00 (the Gaussian) falls behind by π in total, HG10 by 2π, HG11 by 3π. (b) Resonance frequencies of a resonator 30 cm long with mirrors of 60 cm radius of curvature. Longitudinal modes are 500 MHz apart, and transverse modes are shifted by 167 MHz because of the Gouy phase. The l + m = 3 mode coincides with the next longitudinal mode. The heights of the lines only tell the orders apart and carry no meaning._

Why does this phase matter? First, the more complex the shape of a higher-order mode (Section 5), the more it falls behind. A mode of order $N$ lags by $(N+1)\zeta$ (Figure 4(a)) [1]. Second, in the laser resonator we will see in Section 4, the Gouy phase is added to the phase accumulated in one round trip, so the resonance frequency differs slightly from mode to mode (Figure 4(b)). Within each longitudinal-mode spacing $c/2L$, set by the resonator length, the transverse modes sit at positions shifted by the Gouy phase [1]. In the resonator of Figure 4(b) the shift is exactly 1/3 of the spacing, because the resonator's Gouy phase difference is $\pi/3$.

## 3. Where does a lens put the waist of a beam, and how small is it?

When a lens focuses a beam, where does the waist form, and how small is it? Rewriting the "source at an imaginary distance" of Section 1 makes the calculation simple. Bundle the state of the beam into a single complex number:

$$
q = z - i z_0,\qquad \frac{1}{q} = \frac{1}{R} + i\frac{\lambda}{\pi W^2}
$$

Here $q$ is called the complex radius of curvature (the q-parameter) of the beam. Its real part is the distance from the waist and the magnitude of its imaginary part the Rayleigh range. The real part of its reciprocal holds the wavefront curvature and the imaginary part holds the beam radius. Some textbooks use the opposite sign convention [1]. When the beam travels a distance $d$, $d$ is added to $q$; when it passes a lens of focal length $f$, $1/f$ is subtracted from $1/q$. Combining these two rules into one gives [1]

$$
q_2 = \frac{A q_1 + B}{C q_1 + D}
$$

Here $A, B, C, D$ are exactly the elements of the ray transfer matrix of [Geometrical Optics Foundations 2](/en/posts/geometric-optics-paraxial-abcd/). In other words, the matrices that traced rays trace Gaussian beams as they are, in media of the same refractive index. Multiply the matrices of several lenses, mirrors, even graded-index media into one, and plug it into this formula. This is called the ABCD law.

<img src="/assets/img/posts/wave-optics-gaussian-beams/en/fig5-lens-focusing.png" alt="Propagation of a beam after a lens and the ABCD law, output waist position, and focal spot size and depth versus f-number" width="780">
_Figure 5. (a) A beam with a 20 µm waist radius (Rayleigh range 2.5 mm), placed with its waist 2.5 mm in front of a lens of 1 mm focal length, is focused by the lens. The Rayleigh range and the distance in front of the lens are separate quantities that happen to share the value 2.5 mm. The colour is the intensity calculated by wave propagation (angular spectrum and lens phase), and the cyan curve is the beam radius from the ABCD law. The waist (white dashes) forms 0.18 mm behind the focal point (green dots). (b) Output waist position versus input waist position. The black dashes are the geometrical-optics image position, and the coloured curves are for input-beam Rayleigh ranges of 0.25f, 0.5f and f. (c) Diameter and depth of the focus when a collimated beam (532 nm) is focused at f-number $$F_\# = f/D$$. Below f-number 1 the paraxial approximation breaks down and the formula no longer holds (Part 6)._

In Figure 5(a) the calculation and the ABCD law overlap. In this post's calculation, for three different distances in front of the lens, the waist position and size found by wave propagation agree with the ABCD law to within 1 µm. Yet the waist does not form at the focal point. In geometrical optics a point 2.5 mm in front of the lens would be imaged 1.67 mm behind it, but the waist forms at 1.18 mm. Figure 5(b) shows why. The geometrical-optics image position (dashes) runs off to infinity as the input approaches the focal point, but the output waist of a Gaussian beam moves only within a finite range. The thicker the input beam (the longer its Rayleigh range), the more it departs from geometrical optics. A beam with a Rayleigh range near zero, like a point source, approaches the lens formula. If the input waist is at the front focal point, the output waist forms at the back focal point regardless of the beam's thickness. That is why all the curves pass through (1, 1). In other words, when dealing with beams, compute $q$ instead of using the lens formula.

The most common case in practice is focusing a collimated beam with a lens. Focusing a collimated beam of diameter $D = 2W_0$ with focal length $f$ gives a focal spot of diameter [1]

$$
2W_0' \approx \frac{4}{\pi}\lambda F_\#,\qquad F_\# = \frac{f}{D}
$$

In other words, the spot size is proportional to the f-number, and the depth of focus $$2z_0' = 2\pi W_0'^2/\lambda$$ is proportional to the square of the f-number (Figure 5(c)). Focusing at 532 nm with an f-number of 2 gives a spot 1.4 µm across and 5.4 µm deep. At f-number 10 it is 6.8 µm and 135 µm. This trade-off, that a small spot lasts only a short distance, has to be weighed when designing laser machining, optical recording and laser-scanning microscopes. To keep a small spot over a long depth, one of the two must be given up. Conversely, widening the beam with a telescope used backwards (a beam expander) reduces the divergence by as much as the waist grows, exactly as in the lunar ranging of Section 1.

## 4. Why does a laser emit a Gaussian beam?

Why is the beam from a laser Gaussian of all shapes? Light in a laser is amplified as it bounces back and forth between two mirrors countless times. A shape survives only if the light returns in the same shape after one round trip. The wavefront of a Gaussian beam is spherical, so if a mirror's curvature equals the wavefront curvature at that point, the beam retraces its own path after reflection [1, 3]. In other words, the mirror curvatures and spacing pick out the one beam that satisfies the $R(z)$ condition, and that beam's $z_0$ and $W_0$ set the size of the laser beam. Fowles puts it the other way round: pick any two wavefronts of a mode and place mirrors of the same curvature there, and you always get a resonator with the same mode [3].

If the two mirrors have the same curvature (symmetric) and a radius of curvature equal to their spacing ($R = L$), their focal points meet in the middle. This is called a confocal resonator. Its waist is in the middle, $W_0 = \sqrt{\lambda L/2\pi}$, and the beam radius on the mirrors is $\sqrt{2}$ times that [1, 2]. In Hecht's example, a confocal He-Ne laser 28 cm long has a waist radius of 0.168 mm and a full divergence of 2.40 mrad [2]. This post's calculation gives the same values.

Not every mirror arrangement has such a beam. If the mirrors' curvature is too small or their spacing too long, light leaks out sideways with every round trip. [Geometrical Optics Foundations 2](/en/posts/geometric-optics-paraxial-abcd/) showed the condition for rays to stay confined in a row of lenses, and a resonator is such a row of lenses folded up. A mirror is equivalent to a lens whose focal length is half its radius of curvature. With $g_i = 1 - L/R_i$, the condition for light to stay confined is [1, 3]

$$
0 \le g_1 g_2 \le 1
$$

In other words, if $g_1 g_2$ lies between 0 and 1 the resonator is stable; outside that range it is unstable and light leaks out (Figure 6(c)). Two flat mirrors ($g_1 = g_2 = 1$), the confocal resonator ($g_1 = g_2 = 0$) and the concentric resonator ($g_1 = g_2 = -1$) all sit on the boundary. A resonator on the boundary becomes unstable with the slightest misalignment. According to Fowles, a flat-mirror resonator must be aligned to about one arcsecond (1/3600 of a degree), whereas a confocal resonator tolerates misalignment of about a quarter of a degree [3]. Unstable resonators have their uses too. High-gain, high-power lasers (CO₂ lasers and the like) sometimes use unstable resonators deliberately [2]. In this post's calculation, sending rays around 382 randomly chosen resonators, whether they stay confined agrees with the inequality above in every case.

But light does not start out Gaussian. How does a resonator pick out the Gaussian shape? In 1961 A. G. Fox and T. Li solved this question on a computer [3, 6]. Bouncing between mirrors is the same as passing through a passage lined with identical apertures at regular intervals. Feed in uniform light and compute the round trips one after another with the diffraction integral, and depending on the resonator, within a few dozen to a few hundred passes the shape settles; after that only its amplitude decays, by a constant fraction per pass. That shape is the resonator's mode. In the same year G. D. Boyd and J. P. Gordon solved the modes of the confocal resonator in terms of spheroidal wave functions and showed that, for a large Fresnel number, they are approximated by products of a Gaussian and Hermite polynomials [3]. The 1966 paper in which H. Kogelnik and Li collected these results is still the basic reference for laser-beam calculations [1, 7].

<img src="/assets/img/posts/wave-optics-gaussian-beams/en/fig6-resonator-fox-li.png" alt="How a resonator mode forms under Fox–Li iteration, loss per pass, and the stable region" width="780">
_Figure 6. (a) Intensity on the mirror, computed by feeding uniform light into a symmetric resonator (633 nm) 30 cm long with mirrors of 60 cm radius of curvature and 0.6 mm aperture radius, and iterating round trips. After 20 passes a rounded central peak has formed but the width keeps growing, and only after 1500 passes does it overlap the Gaussian calculated from the self-consistent q (black dashes). (b) The fraction of light lost per pass. With a small aperture (0.45 mm) it converges to 0.7 % within a few dozen passes; with a large one (1 mm) higher modes linger and the loss oscillates. (c) The stable region $$0 \le g_1g_2 \le 1$$ (blue) and representative resonators._

Figure 6(a) shows uniform light turning into a Gaussian shape. After one pass, diffraction fringes from the aperture edge appear; after a few more, the edges are shaved off and the middle rises into a rounded peak. After 1500 passes the beam radius on the mirror is 265 µm, the same as the 264 µm calculated by matching the mirror curvature. Figure 6(b) carries a practical lesson. Higher-order modes are wider and are cut off more at the mirror aperture. With a small aperture the higher modes die quickly and only the Gaussian mode remains, but the Gaussian mode also loses a bit more. With a large aperture the loss is small, but higher modes survive for a long time. That is why single-transverse-mode lasers put an aperture inside the resonator to filter out the higher modes. The criterion for this judgement is the Fresnel number $N_F = a^2/\lambda L$, the square of the aperture radius measured in units of $\sqrt{\lambda L}$ (the very number from Part 3) [3]. For the resonator of Figure 6, at $N_F = 1.1$ the higher modes die quickly and only the Gaussian remains, while at $N_F = 1.9$ the loss is still falling after 600 passes. In this resonator, mode selection is reliable when $N_F$ is close to 1.

A brief note on history. In July 1960 T. H. Maiman announced the ruby laser, and in December 1960 A. Javan, W. R. Bennett and D. R. Herriott built the first continuously operating He-Ne laser, which they reported in February 1961 [2]. The papers of Fox and Li and of Boyd and Gordon appeared the same year. The laser and the theory of its beam were born almost at the same time.

## 5. How far does a real laser beam depart from the Gaussian?

Is the Gaussian the only thing a resonator picks out? No. A Gaussian multiplied by Hermite polynomials is also a solution of the paraxial wave equation, and it keeps its shape as it spreads. These are called Hermite–Gaussian (HG) modes, labelled $$\mathrm{HG}_{lm}$$ by the numbers of nodes $l$, $m$ in the $x$ and $y$ directions. Written in cylindrical coordinates they become Laguerre–Gaussian (LG) modes (Figure 7(a)) [1]. The $\mathrm{TEM}_{00}$ in laser specifications is exactly the Gaussian mode. The higher the mode order, the larger the Gouy phase, $(l + m + 1)\zeta$ (Section 2). In this post's calculation, propagating four HG modes and three LG modes, the phase lag at $z = z_0$ agrees with the formula to within 0.002 rad.

<img src="/assets/img/posts/wave-optics-gaussian-beams/en/fig7-higher-modes-m2.png" alt="Intensity and phase of Hermite–Gaussian and Laguerre–Gaussian modes, and beam width versus beam quality M²" width="780">
_Figure 7. (a) Intensity of higher-order modes seen at the waist. $$\mathrm{LG}_0^1$$ is a doughnut with an empty centre, and its phase winds by $2\pi$ as it goes once around the axis, a helix (last panel). Each intensity panel is normalised to its own maximum. (b) Beam quality. Second-moment width in $x$ of a Gaussian ($M^2 = 1$) and of a beam that mixes $$\mathrm{HG}_{00}$$ and $$\mathrm{HG}_{10}$$ incoherently with equal intensity ($M^2 = 2$). Lines are the hyperbolic formula and circles are values measured from wave propagation. The mixed beam is $\sqrt{2}$ times wider both at the waist and far away._

Among LG modes, a beam whose phase winds in a helix, like $\mathrm{LG}_0^1$, has an undefined phase at its centre and zero intensity there. Such a beam is called an optical vortex, and it is used as an example showing that light carries orbital angular momentum [1, 2].

Real laser beams are often not pure $\mathrm{TEM}_{00}$. Higher modes mix in, or, as with semiconductor lasers, the shape is distorted from the start. How is such a beam evaluated? Define the beam width by the second moment of the intensity, and the width of any beam traces the same hyperbola with distance as a Gaussian. Only the product of the waist radius $$W_{0m}$$ and divergence angle $\theta_m$ is larger than for a Gaussian. The ratio

$$
M^2 = \frac{\pi W_{0m}\theta_m}{\lambda}
$$

is called the beam quality factor [1]. In other words, $M^2$ is "how many times more it spreads than a Gaussian with the same waist", and it is always at least 1. The mixed beam in Figure 7(b) has $M^2 = 2$, and the calculation gives 2.00 as well. As typical values, Saleh and Teich give below 1.1 for He-Ne lasers, 1.1–1.3 for ion lasers, 1.1–1.7 for single-transverse-mode semiconductor lasers and 3–4 for multimode lasers [1]. The practical consequence is simple: the spot-size formula of Section 3 must be multiplied by $M^2$. A beam with $M^2 = 1.5$ focused by the same lens gives a spot 1.5 times larger. To measure it, focus the beam with a lens, measure its width before and after the focus with a camera, and fit a hyperbola to read off the waist and the divergence together [1].

Is there really no beam that spreads less than a Gaussian? The Bessel beam glimpsed in Part 3 keeps the bright core near its axis the same size over long distances. But its outer rings go on without end, so its energy is spread widely and its second-moment width is infinite [1, 2]. Among beams that hold finite energy in a small region, the one that spreads least is still the Gaussian, because, as seen in Section 1, it has the smallest product of spatial width and spread angle when measured by second moments.

## Summary

- A laser beam cannot escape diffraction, and the Gaussian shape spreads least within that limit. Waist radius $W_0$, Rayleigh range $z_0 = \pi W_0^2/\lambda$, divergence angle $\theta_0 = \lambda/\pi W_0$, and $W_0\theta_0 = \lambda/\pi$. Within $z_0$ the beam is nearly parallel; beyond it, it spreads in a cone.
- A green beam 1 mm across has $z_0$ of 1.5 m and becomes 260 km across at the Moon. Widening the beam makes it spread less.
- A beam passing through focus falls behind a plane wave by π in phase (the Gouy phase). Light confined sideways accumulates phase more slowly along the axis because of its oblique components. Higher modes fall behind $(N+1)$ times as much, and a resonator's resonance frequencies split by mode.
- A beam is described by a single complex number $q$, which passes through an optical system as $q_2 = (Aq_1 + B)/(Cq_1 + D)$ with the ABCD matrices of geometrical optics. The output waist does not go to the geometrical image position. The focal spot of a collimated beam has diameter $$(4/\pi)\lambda F_\#$$, and its depth scales as $$F_\#^2$$.
- A resonator picks out the beam whose wavefront curvature matches the mirror curvature. The stability condition is $0 \le g_1g_2 \le 1$. Fox–Li iteration grows a Gaussian mode from uniform light, and an aperture filters out higher modes.
- The quality of a real beam is measured by $M^2$: how many times the product of waist and divergence exceeds a Gaussian's. The focal spot grows by the same factor.

Every formula in this post rests on the paraxial approximation, the assumption that the beam travels almost in one direction. Focus light with a lens whose numerical aperture (NA) approaches 1, like a microscope objective, and that assumption breaks. When light converges on the focus from all sides, the direction of the electric field, that is, the polarization, changes the shape of the focus. The final Part 6 is high-NA vector focusing.

## References

1. B. E. A. Saleh and M. C. Teich, *Fundamentals of Photonics*, 3rd ed. (Wiley, 2019), Section 3.1 — the Gaussian beam (complex source point, beam radius, curvature and Gouy phase, Rayleigh range and depth of focus, the 86 % criterion, beam quality $M^2$ and typical values). Section 3.2 — beams through lenses and mirrors, focal spot size, the ABCD law. Sections 3.3–3.5 — Hermite–Gaussian and Laguerre–Gaussian modes, Bessel beams. Section 11.2 — stability of spherical-mirror resonators, Gaussian modes, resonance frequencies.
2. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Section 13.1 — principles and history of the laser, resonator types and stability, waist and divergence of the $\mathrm{TEM}_{00}$ beam (Example 13.6), transverse and longitudinal modes. Section 10.2.7 — Bessel beams.
3. G. R. Fowles, *Introduction to Modern Optics*, 2nd ed. (Holt, Rinehart and Winston, 1975), Section 9.6 — Fox and Li's computation of resonator modes, diffraction loss and the Fresnel number, size and wavefronts of the Gaussian beam, alignment tolerance of resonators. Section 10.5 — stability condition of resonators.
4. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Section 8.8.4 — the phase anomaly near focus (Gouy phase).
5. S. Feng and H. G. Winful, "Physical origin of the Gouy phase shift," *Optics Letters* 26, 485–487 (2001). [doi:10.1364/OL.26.000485](https://doi.org/10.1364/OL.26.000485)
6. A. G. Fox and T. Li, "Resonant modes in a maser interferometer," *Bell System Technical Journal* 40, 453–488 (1961).
7. H. Kogelnik and T. Li, "Laser beams and resonators," *Applied Optics* 5, 1550–1567 (1966). [doi:10.1364/AO.5.001550](https://doi.org/10.1364/AO.5.001550)

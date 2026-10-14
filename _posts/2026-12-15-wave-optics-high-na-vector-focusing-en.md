---
layout: post
date: 2026-12-15 20:00:00 +0900
title: "Wave Optics Foundations 6 — Why Is the Spot of Light from a Microscope Lens Not Round?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/wave-optics-high-na-vector-focusing/
page_id: wave-optics-high-na-vector-focusing
categories: [Optics, Wave Optics]
tags: [wave-optics, vector-diffraction, high-na, radial-polarization, richards-wolf]
description: "The focus of a lens with numerical aperture near 1 is computed by vector diffraction including polarization, tracing where the elongated spot, the longitudinal field, radial polarization and the Maltese cross come from."
math: true
---

Numbers are engraved on the side of a microscope objective. The lens in Figure 1 reads "100×/1,25 Oil": magnification 100, numerical aperture (NA) 1.25, used immersed in oil. In oil with refractive index 1.515, NA 1.25 corresponds to a half-angle of 56°. Light passing the edge of the lens is bent into the focus at 56° to the optical axis.

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/ext-oil-objective-100x.jpg" alt="Side marking of a Zeiss CP-Achromat 100×/1.25 oil-immersion objective" width="600">
_Figure 1. An oil-immersion objective (Zeiss CP-Achromat). "100×/1,25 Oil" gives the magnification, numerical aperture and immersion medium, and "∞/0,17" means infinity-corrected optics and a 0.17 mm cover glass. (Source: QuodScripsiScripsi, [Objective zeiss 100x](https://commons.wikimedia.org/wiki/File:Objective_zeiss_100x.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

[Wave Optics Foundations 3](/en/posts/wave-optics-diffraction/) and [Part 4](/en/posts/wave-optics-fourier-imaging/) computed the focus of such a lens as an Airy disk. That is scalar theory, which treats light as a single amplitude. But compute the focus of linearly polarized light gathered by an NA 0.9 lens, and the focus is not round. It is 0.73λ along the polarization and 0.54λ across it, about 35 % elongated. Scalar theory's answer is a circle 0.57λ across. On top of that, the focus acquires an electric-field component pointing along the direction of travel. This post traces where these differences come from.

- Up to what NA does scalar theory hold? Why does polarization matter at high NA?
- Why is the focus elongated along the polarization? How does an electric field along the direction of travel arise?
- Can changing the polarization make the focus smaller?
- What traces does this effect leave on microscopy and metrology?

This is the last part of Wave Optics Foundations. [Wave Optics Foundations 5](/en/posts/wave-optics-gaussian-beams/) rested on the paraxial approximation, that light travels almost in one direction. This post drops that assumption.

## 1. Up to what NA does scalar theory hold?

Light is an electromagnetic wave, and its electric field is a vector with a direction. Yet Wave Optics Foundations 3–5 treated light as a single complex amplitude. How could they? Born and Wolf set out the justification as follows [1]. If rays make angles with the optical axis of no more than 10–15°, reflection losses at lens surfaces are a few percent and the direction of polarization rotates by only a few degrees. The direction of the electric field can then be taken as the same everywhere on the reference sphere, and the field at the focus has no component along the direction of travel ($E_z = 0$). In other words, scalar theory holds only for systems of "moderate numerical aperture" [1]. Goodman likewise notes that scalar theory holds only when the aperture is much larger than the wavelength and the field is observed well away from it [2].

What changes at high NA? Take the simplest scene. Two plane waves of equal intensity cross at $\pm\theta$ to the optical axis on their way to the focus (Figure 2(a)). If both electric fields point out of the page (s polarization), the two fields are always parallel and the interference fringes have visibility 1. If the fields lie in the plane of the page (p polarization), each field must be perpendicular to its own direction of travel, so they are tilted by $2\theta$ relative to each other. The dot product of the two vectors is $\cos 2\theta$, and the fringe visibility drops by that much.

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/en/fig2-two-beam-polarization.png" alt="Polarization of two crossing plane waves and the visibility of their interference fringes" width="780">
_Figure 2. (a) Two plane waves crossing at $\pm\theta$ toward the focus. With s polarization (⊙) the two fields are parallel; with p polarization (↔) each is perpendicular to its own direction of travel, so they are tilted $2\theta$ to each other. (b) Fringe visibility. It is always 1 for s polarization and $\lvert\cos 2\theta\rvert$ for p polarization. Circles are values computed by adding the two field vectors directly. The vertical dotted lines are the angles corresponding to NA 0.3, 0.71 and 0.9 in air._

In Figure 2(b), at NA 0.3 (17°) the p-polarized visibility is 0.82. At NA 0.71 (45°) the two fields are at right angles and the fringes vanish completely. The Fresnel–Arago law of [Wave Optics Foundations 1](/en/posts/wave-optics-interference/), that orthogonal polarizations do not interfere, is exactly what happens inside the focus. At NA 0.9 (64°) the visibility climbs back to 0.62, but the bright and dark places swap. Where the two waves are in phase, the transverse components reinforce and the longitudinal ($z$) components cancel. Where they are out of phase, the longitudinal components reinforce instead. Beyond 45° each wave's longitudinal component is larger than its transverse one, so the places where the longitudinal components reinforce become the brighter ones. In other words, at a high-NA focus, light of the same polarization adds differently depending on the direction it comes from. This is what scalar theory misses. The focus is a sum of plane waves arriving from all sides, so this effect changes the shape of the focus.

## 2. How is the focus of a real lens computed?

So how is the focus of a real lens computed? In 1959 Richards (B. Richards) and Wolf (E. Wolf) solved this problem with electromagnetic theory [1, 3]. Their picture is simple. When collimated light passes an aberration-free lens, the lens transfers each ray onto a reference sphere of radius $f$ centred on the focus. From each point of that sphere a plane wave heads for the focus, and the field at the focus is the sum of those plane waves. It is the plane-wave decomposition of [Wave Optics Foundations 3, Section 7](/en/posts/wave-optics-diffraction/) used in reverse. This time, however, each plane wave carries its own electric-field vector. Figure 3(a) shows the scene.

To turn this picture into a calculation, two things must be settled. First, where on the sphere does a ray entering the lens at height $h$ land? An aplanatic lens is one that satisfies the sine condition, and then $h = f\sin\theta$ [1]. This is the very sine condition that appeared in [Geometrical Optics Foundations 4](/en/posts/geometric-optics-aberrations/) as the condition for zero coma. Second, what is that ray's amplitude? The same energy must flow through an area element in front of the lens and the corresponding area element on the sphere. This condition puts the following factor on the amplitude on the sphere [3]:

$$
A(\theta) = A_0\sqrt{\cos\theta}
$$

In other words, light that entered with uniform intensity is weaker at the edge of the sphere. A ring at the edge of the lens is spread over a band on the sphere that is wider than the ring itself.

What happens to the polarization? Split the electric field in front of the lens into radial and azimuthal components, and the answer is simple [3]. The azimuthal component is perpendicular to the plane in which the ray bends, so its direction is unchanged when the ray bends. The radial component tilts with the ray by $\theta$ and becomes the direction $\hat{\theta}$ perpendicular to the ray.

$$
\hat{\rho} \;\to\; \hat{\theta} = \cos\theta\,(\cos\varphi\,\hat{x} + \sin\varphi\,\hat{y}) - \sin\theta\,\hat{z},\qquad \hat{\varphi} \;\to\; \hat{\varphi}
$$

Here $\varphi$ is the azimuth over the pupil. In other words, radial polarization tilts and gains a $z$ component $-\sin\theta$. In front of the lens, $x$-polarized light splits by azimuth into a radial component $\cos\varphi$ and an azimuthal component $-\sin\varphi$, so on the sphere it has all three components (Figure 3(b)–(d)).

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/en/fig3-reference-sphere-polarization.png" alt="Reference sphere of an aplanatic lens and the three components of the electric field on it" width="780">
_Figure 3. (a) An aplanatic lens. Collimated light entering at height $h$ passes through the point $h = f\sin\theta$ on a reference sphere of radius $f$ on its way to the focus. The electric field (↕) tilts to stay perpendicular to the ray. (b)–(d) Magnitude of each component of the field on the reference sphere, plotted in pupil coordinates, for $x$-polarized light entering an NA 0.9 lens. The $x$ component weakens by $\sqrt{\cos\theta}$ toward the edge, a small $y$ component appears along the diagonals, and the $z$ component is large at both edges along the polarization direction ($x$)._

Figure 3(d) is the key. The $z$ component is large at the left and right edges of the pupil, with opposite signs. The fields of light arriving from the left and from the right tilt in opposite directions, like the p polarization of Figure 2. These two $z$ components cancel at the exact centre of the focus, but slightly off centre in the $x$ direction they fail to cancel and remain.

Now add up the plane waves. The integral over azimuth can be done in advance with Bessel functions, so the focal field for $x$-polarized input can be written with three one-dimensional integrals [3]:

$$
\begin{aligned}
E_x &\propto I_0 + I_2\cos 2\varphi, \qquad E_y \propto I_2\sin 2\varphi, \qquad E_z \propto -2i\,I_1\cos\varphi \\
I_m &= \int_0^{\alpha}\sqrt{\cos\theta}\,\sin\theta\; w_m(\theta)\, J_m(k\rho\sin\theta)\, e^{ikz\cos\theta}\, d\theta
\end{aligned}
$$

Here $\alpha$ is the lens half-angle ($\mathrm{NA} = n\sin\alpha$), $(\rho, \varphi, z)$ are cylindrical coordinates near the focus, and the weights are $w_0 = 1+\cos\theta$, $w_1 = \sin\theta$, $w_2 = 1-\cos\theta$. In other words, $I_0$ is the main term, nearly the Airy pattern of scalar theory, and $I_1$ and $I_2$ are corrections produced by the tilted polarization. When $\theta$ is small, $w_1$ and $w_2$ go to 0 and scalar theory returns. This post's calculation compared these integrals with a direct sum of tens of thousands of vector plane waves on the reference sphere and obtained the same field to within $10^{-6}$. At NA 0.1 the focal intensity agrees with the scalar Airy pattern $[2J_1(v)/v]^2$ to within 0.2 %.

## 3. What shape is the focus of linearly polarized light?

Now the opening question can be answered. Figure 4 shows the intensity of each component in the focal plane when $x$-polarized light is gathered by an NA 0.9 lens (in air).

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/en/fig4-linear-focus-components.png" alt="Ex, Ey and Ez components and total intensity of a linearly polarized focus at NA 0.9" width="780">
_Figure 4. Focal-plane intensity for $x$-polarized light gathered by an NA 0.9 lens (maximum of the total intensity set to 1). (a) $\lvert E_x\rvert^2$. (b) $\lvert E_y\rvert^2$: up to 0.005 at four places along the diagonals. (c) $\lvert E_z\rvert^2$: two lumps up to 0.15 on either side along the polarization. (d) Total intensity. The cyan curve is the half-maximum contour and the white dashed circle is the half-maximum circle of the scalar Airy pattern. The unit λ is the wavelength in the medium._

Let us read Figure 4. $E_x$ (a) is nearly the Airy disk of scalar theory. $E_y$ (b) is small, at 0.5 %. The problem is $E_z$ (c). As seen in Section 2, it is zero at the centre and reaches 15 % of the central intensity at two places about 0.4λ off centre along the polarization. Adding these two lumps to $E_x$ stretches the total intensity along the polarization (d). The full width at half maximum is 0.73λ along the polarization and 0.54λ across it. Compared with scalar theory's 0.57λ, it is 28 % wider along the polarization and actually 6 % narrower across it.

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/en/fig5-na-dependence.png" alt="Full width at half maximum of the focus and strength of the longitudinal component versus NA" width="780">
_Figure 5. (a) FWHM of the linearly polarized focus multiplied by NA. The scalar Airy pattern gives 0.514 regardless of NA (dashes). As NA grows, the width along the polarization (red) widens and the width across it (blue) narrows. (b) Maximum of the longitudinal component $\lvert E_z\rvert^2$ as a percentage of the central intensity. At low NA it is proportional to NA²._

Figure 5 shows the size of this effect by NA. At NA 0.3 the difference between the two widths is a little over 2 % and $E_z$ is a little over 1 %. Scalar theory is enough. At NA 0.6 the width difference is 10 % and $E_z$ exceeds 5 %. At NA 0.9 and above, the difference exceeds 30 %. At low NA the strength of $E_z$ is proportional to NA². In this post's calculation, doubling NA multiplies it by exactly four. A rule of thumb: below NA 0.5 scalar calculation is fine, and above 0.7 vector calculation is needed.

One point is worth pausing on. How does the NA 1.4 of an oil-immersion lens compare with NA 0.9 in air? What sets the vector effect is the ray angle, that is, $\sin\alpha = \mathrm{NA}/n$. An NA 1.4 oil lens has $\sin\alpha = 1.4/1.515 = 0.92$, so its focus is distorted about as much as with NA 0.92 in air. The unit of length, however, is the wavelength in the medium, $\lambda/n$. So oil immersion makes the focus $1/n$ times smaller, but it does not reduce the proportion of the distortion.

## 4. Can changing the polarization make a smaller focus?

If polarization changes the shape of the focus, can a well-chosen polarization make a smaller focus? Let us compare three polarizations (Figure 6).

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/en/fig6-cylindrical-polarizations.png" alt="NA 0.9 foci of circular, radial and azimuthal polarization, and spot area versus annular pupil" width="780">
_Figure 6. Foci formed by an NA 0.9 lens (each panel normalized to its own maximum; cyan is the half-maximum contour). (a) Circular polarization: round. (b) Total intensity of radial polarization. (c) Only the $\lvert E_z\rvert^2$ component of radial polarization. (d) Azimuthal polarization: a doughnut with an empty centre. (e) Half-maximum area with an annular pupil whose centre is blocked. The horizontal axis is the inner NA of the annulus (outer NA 0.9). Stars are the values reported by Dorn et al._

Circular polarization (a) gives a round focus. It is the sum of the $x$- and $y$-polarized foci with a 90° phase difference, so their elongations point in perpendicular directions and average out. Its half-maximum area is 0.30λ², nearly the same as linear polarization (0.31λ²).

Radial polarization is a beam whose electric field points outward everywhere in the pupil. In the picture of Section 2, the whole input is radial, so all of it tilts on the sphere and gains a $z$ component. And this $z$ component has the same sign from every direction, so it reinforces rather than cancels at the centre of the focus. A sharp $E_z$ peak therefore forms at the centre (c). The transverse component $E_\rho$, on the other hand, is zero at the centre and forms a ring around it. The total intensity (b) is the sum of the two.

So does radial polarization make a smaller focus? Figure 6(e) gives the answer. The conclusion first: it depends on the conditions. Using the full pupil (horizontal axis 0), radial polarization has a half-maximum area of 0.47λ², actually wider than linear polarization's 0.31λ². The ring of the transverse component widens the spot. But block the centre of the pupil and use only the large angles at its edge, and things change. At large angles the share of the $z$ component grows in proportion to $\sin\theta$, and the transverse component shrinks. Blocking out to an inner NA beyond about 0.6, radial polarization overtakes linear polarization. With a narrow annulus of NA 0.85–0.9, radial polarization gives 0.16λ² and linear polarization under the same conditions gives 0.27λ². This agrees with the 0.16λ² and 0.26λ² that Dorn (R. Dorn), Quabis (S. Quabis) and Leuchs (G. Leuchs) reported from measurements in 2003 [3, 4]. In other words, "radial polarization makes a smaller focus" is a statement about high NA with only the edge of the pupil in use. The price is discarding a substantial part of the light.

With azimuthal polarization (d), the electric field points around the pupil. The azimuthal component does not tilt, so there is no $z$ component, and the transverse components cancel at the centre, leaving a doughnut with an empty centre. This doughnut has many uses. Examples are the STED microscope, which shines light that switches off fluorescent molecules in a doughnut shape and leaves only the centre to raise resolution, and optical tweezers that trap particles in the dark centre [3]. The strong longitudinal component of radial polarization is used to read the orientation of molecules or to concentrate light at the tip of a metal probe [3].

## 5. What does high NA leave on measurements?

The calculations so far treated the lens as a perfect aplanatic lens. A real lens is made of many glass surfaces, and at each surface the transmittances for p and s polarization differ. The Fresnel coefficients of [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/) apply as they are. The larger the angle at which a ray passes, the larger the difference between p and s, and the original polarization rotates a little.

<img src="/assets/img/posts/wave-optics-high-na-vector-focusing/en/fig7-maltese-cross.png" alt="Maltese-cross pattern between crossed polarizers and extinction ratio versus NA" width="780">
_Figure 7. (a) Fraction of light leaking through the pupil (log scale) when a cone of NA 0.9 passes a glass plate of refractive index 1.5 (both faces) placed between an $x$ polarizer and a $y$ polarizer. No light leaks on the axes along and perpendicular to the polarization; the most leaks along the diagonals. (b) Extinction ratio over the whole pupil (transmitted light / leaked light), calculated with only 2 and 8 uncoated surfaces._

Figure 7(a) is a view of the pupil between crossed polarizers. On the axes along the polarization ($x$) and perpendicular to it ($y$), the ray's plane of incidence is parallel or perpendicular to the polarization, so the polarization does not rotate. Along the diagonals the p and s components are mixed, the polarization rotates, and that much light leaks through the analyzer. The result is a dark cross and four bright quadrants. This is the Maltese cross that appears when looking at the back focal plane (BFP) of the objective in a polarizing microscope [5, 6]. This leaked light lowers the extinction ratio of the polarizing microscope. In Figure 7(b), just two uncoated air–glass surfaces drop the extinction ratio to about 1000 at NA 0.9, and eight surfaces drop it to about 250.

Measured values for real objectives are of the same order. The extinction ratios measured by Shribak (M. Shribak), Inoué (S. Inoué) and Oldenbourg (R. Oldenbourg) were 600 for an NA 1.25 oil lens, 430 for an NA 1.40 oil lens and 600 for an NA 0.7 dry lens [5]. Adding the rectifier devised by Inoué and Hyde (1957), a zero-power meniscus lens and a half-wave plate behind the lens, raised the extinction ratio of the NA 1.25 lens to 5100 [5]. As a practical remedy, Murphy and Davidson note that closing the condenser diaphragm cuts out the large-angle rays and improves extinction, but costs resolution [6]. The condenser-diaphragm story of [Wave Optics Foundations 4](/en/posts/wave-optics-fourier-imaging/) comes back here.

In metrology this effect shows up as an error. An ellipsometer that measures thin films with focused light uses rays up to large angles, and the objective itself changes the polarization of those rays. Chen (X. Chen) et al. measured the Mueller matrix of an NA 0.95 objective across the pupil. The centre of the pupil barely touched the polarization, but at the edge retardance approached 30° [7]. For an NA 0.55 lens it was below 1.5°, negligible. So high-NA metrology first measures the polarization properties of the lens and puts them into the model. This is the same story told in [Ellipsometer Instrumentation 6](/en/posts/ellipsometer-back-focal-plane-micro/), Section 7, which explained that focusing at high NA breaks the plane-wave assumption. The polarization changes of the objective and the vector field in the focus must be handled together. The longitudinal component $E_z$ in the focus is invisible to an ordinary microscope, but there is work that imaged it directly by converting the polarization [8].

## Summary

- Scalar theory holds when ray angles are small (10–15° by Born and Wolf's criterion). At high NA, the electric fields of plane waves converging from all sides tilt relative to one another, and light of the same polarization adds differently depending on direction. When two plane waves converge at 45°, p polarization does not interfere at all.
- Richards–Wolf theory sees an aplanatic lens as transferring rays onto a reference sphere. The amplitude on the sphere carries a factor $\sqrt{\cos\theta}$, radial polarization tilts and gains a $z$ component, and azimuthal polarization is unchanged. The focal field is the sum of these vector plane waves.
- The focus of linearly polarized light is elongated along the polarization. At NA 0.9 it is 0.73λ × 0.54λ, and a longitudinal component $E_z$ appears on either side along the polarization at 15 % of the central intensity. Below NA 0.5 scalar calculation is fine; above 0.7 vector calculation is needed.
- Circular polarization gives a round focus. Radial polarization gathers a strong $E_z$ at the centre, but a smaller focus is obtained only when using just the edge of the pupil (0.16λ² with an NA 0.85–0.9 annulus versus 0.27λ² for linear polarization). Azimuthal polarization gives a doughnut.
- Differences in Fresnel transmission at real lens surfaces rotate the polarization, produce the Maltese cross and drop the extinction ratio to the hundreds. High-NA polarization measurements must measure the lens's polarization properties and put them into the model.

This concludes Wave Optics Foundations. In [Geometrical Optics Foundations 1](/en/posts/geometric-optics-fermat-eikonal/), light was rays in the limit of zero wavelength. At the end of Geometrical Optics Foundations 4, once aberrations became as small as a wavelength, the ray picture gave way to diffraction. Wave Optics Foundations 1–2 set out the rule for adding light (interference) and the condition under which it can be added (coherence). Part 3 used that rule to compute the light behind an aperture, and Part 4 rewrote lens imaging as a filter on spatial frequencies. Part 5 saw the laser beam as the solution that spreads least. And in this post, light became once again an electromagnetic wave with a direction. Rays, scalar waves and vector waves are not different theories but different levels of detail in looking at the same light. Wavelength and angle decide which picture is enough.

## References

1. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Section 4.5.1 — the sine condition and aplanatic pairs. Section 8.4 — conditions for the transition to a scalar theory (ray angles of 10–15°, $E_z = 0$ at the focus). Footnote in Section 8.5.2 — literature on vector theories of focusing (Ignatowsky, Richards and Wolf, and others). Section 8.8.1 — the scalar Debye integral.
2. J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed. (McGraw-Hill, 1996), Section 3.2 — conditions under which scalar theory holds.
3. Q. Zhan, "Cylindrical vector beams: from mathematical concepts to applications," *Advances in Optics and Photonics* 1, 1–57 (2009), Section 5 — the Richards–Wolf integrals (sine condition and $\sqrt{\cos\theta}$ weighting, rotation of polarization on the reference sphere, focal fields of linear, radial and azimuthal polarization), the smaller focus of radial polarization, applications. [doi:10.1364/AOP.1.000001](https://doi.org/10.1364/AOP.1.000001)
4. R. Dorn, S. Quabis, and G. Leuchs, "Sharper focus for a radially polarized light beam," *Physical Review Letters* 91, 233901 (2003). (Cited from Reference 3.)
5. M. Shribak, S. Inoué, and R. Oldenbourg, "Polarization aberrations caused by differential transmission and phase shift in high-numerical-aperture lenses: theory, measurement, and rectification," *Optical Engineering* 41, 943–954 (2002). [doi:10.1117/1.1467669](https://doi.org/10.1117/1.1467669)
6. D. B. Murphy and M. W. Davidson, *Fundamentals of Light Microscopy and Electronic Imaging*, 2nd ed. (Wiley-Blackwell, 2013), Chapter 6 — numerical aperture and immersion. Chapter 9 — the polarization cross in the objective's back aperture in a polarizing microscope.
7. X. Chen *et al.*, "Calibration of polarization effect of a high-numerical-aperture objective lens with Mueller matrix polarimetry," *Measurement Science and Technology* 30, 025201 (2019). [doi:10.1088/1361-6501/aaf4d0](https://doi.org/10.1088/1361-6501/aaf4d0)
8. Z. Gao *et al.*, "Polarization-conversion microscopy for imaging the vectorial polarization distribution in focused light," *Optica* 8, 984–994 (2021). [doi:10.1364/OPTICA.422836](https://doi.org/10.1364/OPTICA.422836)

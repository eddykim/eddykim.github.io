---
layout: post
date: 2026-12-05 20:00:00 +0900
title: "Wave Optics Foundations 3 — Why Is the Centre of a Shadow Bright? Fresnel and Fraunhofer Diffraction"
lang: en
lang-exclusive: ["en"]
permalink: /posts/wave-optics-diffraction/
page_id: wave-optics-diffraction
categories: [Optics, Wave Optics]
tags: [wave-optics, diffraction, fresnel-diffraction, fraunhofer-diffraction, airy-disk]
description: "How light spreads past apertures and edges, from the Huygens–Fresnel principle through Fresnel zones, the Airy disk and gratings to the plane-wave decomposition."
math: true
---

Hang a coin-sized disk on a thread and light it with a small point source, and it casts a round shadow on the wall. If light travelled only in straight lines, the whole shadow would be dark. Yet right at its centre there is a bright spot (Figure 1).

<img src="/assets/img/posts/wave-optics-diffraction/ext-arago-spot.jpg" alt="Photograph of the bright spot (Arago spot) at the centre of a disk's shadow" width="460">
_Figure 1. The shadow of a disk 5.8 mm in diameter. Sunlight was narrowed through a hole about 0.5 mm across to act as a point source, with 1.53 m from source to disk and 1.83 m from disk to the paper screen. The bright spot at the centre of the shadow is the Arago spot, and rainbow-coloured rings surround the shadow's edge. (Source: Aleksandr Berdnikov, [A photograph of the Arago spot](https://commons.wikimedia.org/wiki/File:A_photograph_of_the_Arago_spot.png), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

The spot comes with a story. In 1818 the French Academy of Sciences ran a prize competition on diffraction. Fresnel (A. Fresnel) submitted a memoir that treated light as a wave and calculated diffraction. One of the judges, Poisson (S. D. Poisson), a supporter of the particle theory of light, drew a strange conclusion from Fresnel's theory. If the theory were right, he argued, the centre of a disk's shadow would have to be bright, which was absurd. The chair of the committee, Arago (F. Arago), did the experiment at once, and the spot was there. Fresnel won the prize, and the spot is now called the Arago spot or Poisson spot. Maraldi (G. F. Maraldi) had seen it about a century earlier, but the observation was forgotten [1, 2].

This post starts from that scene and answers three questions.

- If light travels in straight lines, how can the inside of a shadow be bright? And why do we rarely see such patterns in daily life?
- How does the light behind an aperture change shape with distance? How is that shape related to the shape of the aperture?
- Even a lens with no aberrations at all images a star as a small disk, not a point. What sets the size of that disk, and how does it bound the performance of telescopes, microscopes and spectrometers?

[Part 1](/en/posts/wave-optics-interference/) covered interference between two beams, and [Part 2](/en/posts/wave-optics-coherence/) covered coherence, the condition for that interference to be visible. Diffraction is interference among countless beams. Section 5 of Part 2 showed that the spatial degree of coherence takes the form $2J_1(x)/x$, and [Geometrical Optics Foundations 4](/en/posts/geometric-optics-aberrations/) ended by saying that once aberrations shrink to a fraction of a wavelength, the image becomes an Airy disk. Sections 4 and 5 of this post show that both results come from the same mathematics.

## 1. Why does light bend at edges?

You can hear sound from behind a wall, but you cannot see light from behind it. Both are waves, so what is different? The answer is the size of the wavelength. Sound at 500 Hz has a wavelength of about 68 cm, comparable to walls and trees, so it spreads widely around obstacles. Visible light has a wavelength of about 0.5 µm, tens of thousands of times smaller than everyday objects [1]. Water waves behave the same way. When a gap is narrower than the wavelength, the waves spread out in semicircles beyond it. When the gap is much wider than the wavelength, they split into a straight bundle of waves and a calm shadow beside it. Straight-line propagation in geometrical optics is the limit of zero wavelength.

So how do we calculate the light behind an aperture when the wavelength is finite? In 1678 Huygens (C. Huygens) regarded each point of a wavefront as a small source emitting a new spherical wave. Fresnel added the interference of Part 1 to this. Add up all the small spherical waves from the points of the aperture, keeping their phases, and you obtain the light at any point behind it. This is the Huygens–Fresnel principle, and in its modern form it reads [3]:

$$
U(P) = \frac{1}{i\lambda}\iint_{\text{aperture}} U(Q)\,\frac{e^{ikr}}{r}\,\cos\theta\; dS
$$

Here $U(Q)$ is the complex amplitude at a point $Q$ of the aperture, $r$ is the distance from $Q$ to the observation point $P$, $k = 2\pi/\lambda$, and $\theta$ is the angle between $QP$ and the aperture normal. In other words, each patch of the aperture sends out a small spherical wave $e^{ikr}/r$, and at the observation point all those arrows (the phasors of Part 1) are added up.

The formula carries three things that are absent from Huygens' picture. First, $1/\lambda$: the strength of the secondary sources is inversely proportional to the wavelength. Second, $1/i$: the secondary sources lead the incident light in phase by 90°. Third, $\cos\theta$: the secondary sources radiate most strongly forward and weakly to the side. Fresnel put these in almost by hand to make his results fit. In 1882 Kirchhoff (G. Kirchhoff) derived them from the wave equation. The first two factors can be understood through one fact: the field seen from afar is proportional to the time derivative of the field in the aperture [3]. Differentiating a sine wave multiplies its amplitude by the frequency (hence $1/\lambda$) and advances its phase by 90°.

One point is worth pausing on. Nothing in the aperture actually emits light. The real cause is the charges oscillating at the edge of the screen. The Huygens–Fresnel principle is a convenient device that replaces their effect with fictitious sources over the aperture. Kirchhoff's formula has a weakness, its boundary conditions contradict each other, and the Rayleigh–Sommerfeld formula fixes it. The two differ only in the obliquity factor, $(1+\cos\theta)/2$ versus $\cos\theta$, and give the same answer at small angles [3]. For this theory, which treats light as a single scalar rather than a vector, to hold, the aperture must be much larger than the wavelength and we must look more than a few wavelengths away from it [3]. This post stays within that range.

What do we get if we actually carry out the integral? Figure 2 shows a 100 µm wide slit illuminated by a plane wave, with the intensity behind it calculated as a function of distance.

<img src="/assets/img/posts/wave-optics-diffraction/en/fig2-slit-near-to-far.png" alt="Light spreading behind a 100 µm slit with distance, and intensity cross-sections at four distances" width="780">
_Figure 2. Intensity behind a slit 100 µm wide (half-width a = 50 µm) illuminated by a 0.5 µm plane wave. (a) Side view out to 60 mm behind the slit; the white dashes are the geometric shadow boundary. (b)–(e) Cross-sections at distances where the Fresnel number $N_F = a^2/\lambda z$ is 10, 1, 0.5 and 0.1. Grey is the geometrical-optics answer, red is the diffraction calculation and the black dashes are the Fraunhofer approximation. Just behind the slit, ripples sit on the shape of the shadow; further away the shape turns into a sinc² pattern and widens._

Figure 2 reads in three stages. Just behind the slit (b), the light keeps nearly the shape of the slit, with ripples only near the edges. A little further (c), the ripples from the edges spread to the middle and change the shape substantially. Far away (e), the light spreads much wider than the slit, and from then on it only grows in size without changing shape. The first two stages are called Fresnel diffraction and the last Fraunhofer diffraction. They are not different phenomena but the same integral seen at different distances [4].

A single number decides which stage applies.

$$
N_F = \frac{a^2}{\lambda z}
$$

Here $a$ is the half-width of the aperture (the radius, for a circle) and $z$ is the distance from the aperture to the observation plane. This number is called the Fresnel number. In other words, it measures the area of the aperture in units of the area $\lambda z$. If $N_F$ is much larger than 1 we see a shadow, near 1 a Fresnel pattern, and much smaller than 1 a Fraunhofer pattern. The $N_F = 0.5$ of Figure 2(d) is roughly the boundary. At that distance, the far-field spread $\lambda z/2a$ equals the slit's half-width $a$.

We can now answer half of the first question. Think of the shadow of a hand on a wall 1 m away. Taking the half-width of the hand as 2.5 cm gives $N_F = 1250$. The shadow follows geometrical optics, and fringes appear only in a band about $\sqrt{\lambda z} \approx 0.7$ mm wide at the edge. Moreover, a lamp is an extended source, so its spatial coherence, discussed in Part 2, is short. The patterns made by different points of the source are slightly offset, and as they overlap even that narrow band is washed out. That is why the photograph in Figure 1 narrowed the sunlight through a 0.5 mm hole.

## 2. A bright spot in the shadow: Fresnel zones

Why is the centre of a disk's shadow bright? The answer appears when the integral above is split up cleverly. Stand at a point $P$ on the axis behind a circular aperture and look back at the aperture. The distance from the centre of the aperture to $P$ is $z$, and the distance to $P$ grows slightly as you move away from the centre. Draw a circle each time that distance grows by half a wavelength, and the aperture is divided into annular bands. These bands are called Fresnel zones. The outer radius of the $m$-th zone is

$$
r_m \approx \sqrt{m\lambda z}
$$

In other words, the zone radius grows as the square root of the zone number, and every zone has the same area $\pi\lambda z$. For collimated He-Ne laser light (633 nm) viewed from a point 1.58 m away, the first zone has a radius of 1 mm [1]. The Fresnel number defined above is nothing other than the number of zones inside the aperture.

Neighbouring zones differ on average by half a wavelength in their distance to $P$, so they nearly cancel each other's contributions. Drawn with phasors, this situation becomes Figure 3(a). Within one zone the phase turns from 0 to 180°, so chopping a zone's contribution into small pieces and chaining them gives a semicircle. The next zone draws the other semicircle. As the zones move outward, the obliquity factor $\cos\theta$ and the distance $1/r$ weaken their contributions slightly, so the circle does not close but winds inward into a spiral.

<img src="/assets/img/posts/wave-optics-diffraction/en/fig3-fresnel-zones-arago.png" alt="Phasor spiral of Fresnel zones, on-axis intensity behind a circular hole, cross-section of a disk's shadow, on-axis intensity of a zone plate" width="780">
_Figure 3. (a) Phasors chaining the contributions of 60 zones seen from an on-axis point (20λ from the aperture). Red and blue alternate, one zone each. The black arrow is the sum with only the first zone open, the green arrow the sum with nothing blocked. (b) On-axis intensity behind a circular hole of radius 1 mm. The horizontal axis is the Fresnel number, and the numbers are the counts of open zones. (c) Intensity on the screen calculated for the geometry of Figure 1 (disk diameter 5.8 mm, 550 nm point source). The grey band is the geometric shadow, and the inset magnifies the centre. (d) On-axis intensity of a zone plate with only the 20 odd-numbered zones open (log scale). At the focal length f = 10 cm the intensity rises 1600-fold, and foci also appear at f/3, f/5 and f/7._

Several things can be read straight off Figure 3(a).

- **With nothing blocked (green)**, the sum is the arrow from the start of the spiral to its centre. It is half as long as the first zone alone (black). In other words, an unobstructed plane wave equals half the contribution of the first zone.
- **With only the first zone open**, the amplitude is twice the unobstructed value and the intensity four times. Blocking most of the light with an aperture actually makes the axis four times brighter.
- **With two zones open**, the two semicircles almost close the circle and the sum is close to zero. Widening the aperture darkens the axis.

Figure 3(b) shows this as a function of distance. Changing the distance behind a 1 mm radius hole changes the number of open zones $N_F$, and the on-axis intensity swings between 4 for odd $N_F$ and 0 for even $N_F$. Moving a screen back and forth behind the aperture makes the centre of the pattern alternately bright and dark.

Now consider the disk. The disk blocks the first few zones. The light on the axis comes from the zones outside the disk, and in phasor terms it is the arrow from some point in the middle of the spiral to its centre. Its length is half the contribution of the zone just outside the blocked ones. The zone contributions weaken very slowly, so this value is almost the same as with nothing blocked. In other words, the centre of the disk's shadow is as bright as it would be without the disk, however many zones the disk covers. The exact on-axis intensity is $z^2/(z^2+a^2)$, close to 1 when the disk is far away. This holds only at the one point on the axis; move slightly to the side and the contributions of the outer zones fall out of step, so it quickly turns dark. The conclusion Poisson called absurd follows from a single zone picture.

Figure 3(c) calculates the same geometry as Figure 1. The disk covers about 18 zones ($N_F \approx 18$). The intensity at the centre of the shadow is 1.00 times the value without the disk, and the first dark ring of the spot has a radius of 0.13 mm on the screen. The spot in the photograph looks larger because the source hole (0.5 mm) is imaged as a 0.6 mm blur on the screen and sunlight contains many wavelengths. The bright fringes outside the shadow also come out of the calculation.

Knowing the zones also lets us focus light without a lens. Cover all the even-numbered zones and leave only the odd ones open, and the contributions of the remaining zones all add in the same direction. This is called a zone plate, and Rayleigh (Lord Rayleigh) devised it in 1871 [1]. Cut the zone radii as $r_m = \sqrt{m\lambda f}$ and a focus appears at distance $f = r_1^2/\lambda$. Each open zone adds an amplitude of 2, so opening 20 odd zones gives an amplitude of 40 and an intensity 1600 times larger (Figure 3(d)) [1]. The same plate also makes foci at $f/3$, $f/5$, $f/7$, …. Seen from distance $f/3$, each zone of the plate covers three zones at that distance; two of the three cancel and one remains. The on-axis intensity is therefore the same as at $f$, but the focus is narrower by the same factor and collects less energy. Because the focal length is inversely proportional to the wavelength, the chromatic aberration is also large. Even so, zone plates replace lenses in the X-ray region. At X-ray wavelengths the refractive index of matter is very close to 1, which makes refractive lenses hard to build.

## 3. Fringes at an edge and Babinet's principle

When a knife edge blocks half the light, where exactly is the boundary of its shadow? For anyone measuring the edge of an object in an image taken with coherent light, this is a practical question.

A straight edge has no circular symmetry, so the zones are strips parallel to the edge rather than concentric circles. Chaining the strips' contributions gives the Cornu spiral of Figure 4(c), an S-shaped curve joining two eyes. The amplitude at an observation point is the length of the arrow from a point on the curve to the upper eye. As the observation point moves into the shadow, the starting point moves toward the upper eye and the arrow shortens. As it moves into the lit side, the starting point travels along the lower spiral and the arrow length rises and falls.

<img src="/assets/img/posts/wave-optics-diffraction/en/fig4-straight-edge.png" alt="Intensity behind a knife edge, the curve rescaled by √(λz), and the Cornu spiral" width="780">
_Figure 4. (a) Intensity at three distances behind a straight edge illuminated by a plane wave (0.5 µm). x = 0 is the geometric shadow boundary, and the negative side is the shadow. (b) Rescaling the horizontal axis to $w$, proportional to $\sqrt{\lambda z}$, makes the three curves collapse onto one. The intensity is 1/4 at the boundary and 1.37 at the first maximum. (c) The Cornu spiral. The amplitude at an observation point is the length of the arrow from the marked point to the upper eye (blue: boundary, red: first maximum, green: in the shadow, w = −1)._

Three things can be read from Figure 4. First, the intensity at the geometric boundary (x = 0) is exactly 1/4. At the boundary the starting point sits at the very middle of the spiral, so the arrow is half as long as the unobstructed one, which is 1/4 in intensity. Second, the intensity falls off monotonically into the shadow without fringes, while on the lit side it rises to a first maximum of 1.37 and then settles in decaying ripples [2, 5]. Third, the width of the fringes is proportional to $\sqrt{\lambda z}$ (Figure 4(b)). Increasing the distance 25-fold widens the fringes only 5-fold.

The first fact has a direct practical use. When viewing an edge with coherent light, taking the edge where the intensity is half (50 %) shifts it away from the shadow. With coherent light the edge is where the amplitude is half, which is where the intensity is 1/4. The same happens when a lens forms the image with coherent illumination, whereas with incoherent illumination 50 % is correct. This difference returns in Part 4 on Fourier imaging. Methods for reading an edge in an image more precisely than a pixel are covered in [Metrology Image Processing 4](/en/posts/imgproc-edge-subpixel/). When applying them to a coherent image, the reference level has to be set again first.

There is one more principle that pairs with the edge. Babinet (J. Babinet) pointed it out in 1837, and it concerns two screens that are each other's inverse, such as an aperture and an opaque screen of the same shape [2]. The diffraction integral is linear, so adding the fields made by each screen gives the unobstructed plane wave. In the Fraunhofer region the unobstructed plane wave gathers only at a single point (the geometrical image point) and is zero in all other directions. Outside that point, therefore, the fields of the two screens have the same magnitude and opposite sign. In other words, their intensity patterns are identical. This is Babinet's principle. In the Fresnel region only the sum of the fields holds, and the intensity patterns are not the same [1].

Thanks to this principle, the diameter of a hair or a thin wire can be measured with a single laser. Since it makes the same pattern as a slit of the same width, measuring the spacing of the zeros with the slit formula of the next section gives the diameter.

## 4. Far away it becomes a Fourier transform: Fraunhofer diffraction

Light that travelled far, as in Figure 2(e), stopped changing shape. How is that shape related to the aperture? Expanding the distance from an aperture point $(\xi, \eta)$ to an observation point $(x, y, z)$ for large $z$ gives

$$
r \approx z + \frac{x^2+y^2}{2z} - \frac{x\xi + y\eta}{z} + \frac{\xi^2+\eta^2}{2z}
$$

The last term is a phase that varies only across the aperture. Keeping it is the Fresnel approximation, and the remaining integral becomes the Fourier transform of the aperture field multiplied by a quadratic phase [3]. If the phase $k(\xi^2+\eta^2)/2z$ that this term produces is much smaller than 1 radian over the whole aperture, that is, if $N_F \ll 1$, the term can be dropped. What remains is

$$
U(x, y) \propto \iint U(\xi, \eta)\, \exp\!\left[-i 2\pi\left(\frac{x}{\lambda z}\xi + \frac{y}{\lambda z}\eta\right)\right] d\xi\, d\eta
$$

Here the proportionality constant contains a phase that depends on $x, y$, but it does not affect the intensity. In other words, the far pattern is the Fourier transform of the aperture function, and a position $x$ on the screen corresponds to the spatial frequency $x/\lambda z$. This is Fraunhofer diffraction [3, 4]. Section 7 returns to the correspondence between spatial frequency and direction.

The condition $z \gg a^2/\lambda$ is stricter than it looks. For a 1 mm diameter hole, $a^2/\lambda$ is 0.5 m, so a few metres suffice. But illuminate a 2.5 cm aperture with 0.6 µm light and you must go much further than 1.6 km [3]. That is why the laboratory uses lenses. A convex lens pulls the pattern at infinity into its focal plane. This is the starting point of Part 4 on Fourier imaging. The fact that the far-field pattern is a Fourier transform is also what the phase-retrieval problem of [Ptychography 2](/en/posts/ptycho-phase-problem-iterative-retrieval/) rests on. The detector measures only the intensity of the Fourier transform and loses the phase.

The simplest example is a slit of width $b$. The Fourier transform of a rectangle function is a sinc function, so the intensity is

$$
I(\theta) = I_0\left(\frac{\sin\beta}{\beta}\right)^2,\qquad \beta = \frac{\pi b\sin\theta}{\lambda}
$$

Here $\theta$ is the angle between the observation direction and the axis. The first zero lies at $\sin\theta = \lambda/b$. At this angle, pairing the light from the upper half of the slit with that from the lower half gives each pair a path difference of half a wavelength, so they all cancel [1]. The side lobes fall off quickly, to 4.7 %, 1.6 % and 0.8 % of the centre [2]. More important than the numbers is the shape: the narrower the slit, the wider the pattern. The product of width and spread angle, $b\sin\theta \approx \lambda$, is constant. This is the same Fourier-pair property as the statement in Part 2 that shorter pulses have broader spectra.

Figure 5 calculates four aperture shapes and their far-field patterns.

<img src="/assets/img/posts/wave-optics-diffraction/en/fig5-aperture-far-field.png" alt="A narrow rectangle, a square, a circle and a segmented hexagonal mirror aperture, with their far-field diffraction patterns" width="780">
_Figure 5. Top row: apertures; bottom row: their far-field intensity (log scale, six decades). (a) A rectangle narrow in the horizontal direction spreads widely in the horizontal direction. (b) A square spreads in a cross perpendicular to its two pairs of sides. (c) A circle makes concentric rings. (d) An aperture made of 18 hexagonal segments and 3 secondary-mirror struts, like the Webb Space Telescope. Six spikes stand perpendicular to the three directions in which straight edges are most common, and the vertical strut adds a horizontal spike._

Two rules can be read from Figure 5. First, the pattern is wide in the direction in which the aperture is narrow (a). Second, a straight edge sends out a spike perpendicular to itself (b, d). The secondary sources lined up along a straight edge cancel each other in directions parallel to the edge, but not in the direction perpendicular to it. That is also why street lamps shine as stars in night photographs: the straight edges of the camera's aperture blades make the spikes.

<img src="/assets/img/posts/wave-optics-diffraction/ext-webb-focused-star.jpg" alt="A star imaged by the Webb Space Telescope after mirror alignment, with eight diffraction spikes" width="640">
_Figure 6. A star imaged by the Webb Space Telescope after it completed primary-mirror alignment (fine phasing) on 11 March 2022. The six thick spikes come from the straight edges of the hexagonal mirror segments, and the two thin horizontal spikes come from the vertical strut. (Source: NASA/STScI, [Webb's First Image of Focused Star](https://commons.wikimedia.org/wiki/File:Webb%27s_First_Image_of_Focused_Star_%28potm2203a%29.tiff), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0))_

Figure 6 is the real thing. Webb's eight spikes are made by the hexagonal mirror segments and the struts [6]. Compared with the calculation in Figure 5(d), the number and directions of the spikes match exactly: six perpendicular to the hexagon edges, and two horizontal ones perpendicular to the vertical strut. The two lower struts run parallel to hexagon edges, so they make no spikes of their own; theirs overlap the spikes from the hexagon edges. For astronomers, the shape of these spikes is a fingerprint of the telescope's aperture.

The link to Part 2 also appears here. The van Cittert–Zernike theorem said that the spatial degree of coherence is the Fourier transform of the source's brightness distribution. This section says that the far field is the Fourier transform of the aperture function. The same mathematics sets coherence on one side and diffraction patterns on the other. That is why the degree of coherence of a round star and the diffraction pattern of a round aperture are both $2J_1(x)/x$. The next section is about that round aperture.

## 5. The circular aperture: the Airy disk and resolution

The apertures of telescopes and microscopes are usually round. What does the far-field pattern of a round aperture look like, and what does its size determine? Computing the Fourier transform of a circle gives the first-order Bessel function $J_1$. Airy (G. B. Airy) first obtained the result in 1835 [2].

$$
I(\theta) = I_0\left[\frac{2J_1(x)}{x}\right]^2,\qquad x = \frac{\pi D\sin\theta}{\lambda}
$$

Here $D$ is the aperture diameter. The bright central disk is called the Airy disk, and the whole pattern including the surrounding rings the Airy pattern. The first dark ring lies at $x = 3.832$, that is, at the angle

$$
\sin\theta_1 = 1.22\,\frac{\lambda}{D}
$$

In other words, it is 22 % larger than the slit's $\lambda/b$. A circle narrows toward its edge, so it behaves like an effectively narrower aperture than a slit of the same width. When a lens forms this pattern in its focal plane, the angle becomes a distance and the radius of the first ring becomes $1.22\lambda f/D$. Figure 7 is a real pattern.

<img src="/assets/img/posts/wave-optics-diffraction/ext-pinhole-airy.jpg" alt="Photograph of the Airy pattern of a 0.09 mm pinhole illuminated by a red laser" width="600">
_Figure 7. The Airy pattern of a 0.09 mm diameter pinhole illuminated by a red laser, photographed 65 mm behind it. The exposure saturates the central disk so that dozens of outer rings show. The photograph is 17.3 mm wide. (Source: Bautsch, [Beugungsscheibchen.k.720](https://commons.wikimedia.org/wiki/File:Beugungsscheibchen.k.720.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

This photograph lets us check the formula. The Fresnel number of this setup is 0.05, well in the Fraunhofer region. Calibrating the scale with the width of the photograph and measuring the dark rings gives radii of 0.52, 0.96 and 1.39 mm for the first three. Their ratio is 1 : 1.86 : 2.69, which matches the ratio 1 : 1.83 : 2.66 of the Airy zeros 3.832, 7.016 and 10.174. The absolute size is 7–9 % smaller than the $1.22\lambda z/D$ = 0.56–0.57 mm obtained by taking the red laser's wavelength as 635–650 nm. The stated pinhole diameter or distance is probably off by that much. The ratio of ring sizes agrees regardless of such scale errors.

<img src="/assets/img/posts/wave-optics-diffraction/en/fig8-airy-resolution-focus.png" alt="Cross-sections of the Airy and slit patterns, encircled energy, two points at the Rayleigh spacing, and the intensity near focus" width="780">
_Figure 8. (a) Far-field intensity of a circle of diameter D and a slit of width D (log scale). The dotted lines are the first three zeros of the Airy pattern. (b) Energy inside a given radius of the Airy pattern. The first, second and third dark rings enclose 83.8 %, 91.0 % and 93.8 %. (c) The image of two mutually incoherent points separated by the radius of the first dark ring (the Rayleigh criterion). The dip in the horizontal cut is 0.735 of the peak. (d) Intensity near the focus of an aberration-free lens (log scale, optical coordinates). The cyan dashes are the geometrical cone of light, and the white circles mark the first zeros on the axis._

Figures 8(a) and (b) put numbers on the shape of the Airy pattern. The rings have 1.75 %, 0.42 % and 0.16 % of the central brightness, much dimmer than the side lobes of a slit [2]. Of the energy, 84 % lies in the central disk and 91 % inside the second dark ring (Rayleigh's result of 1881) [2]. So when talking about the sharpness of an image, looking only at the size of the disk is not far wrong.

When two stars sit close together, their Airy disks overlap. How far apart must they be to look like two? Rayleigh proposed taking as the boundary the case where the centre of one disk falls on the first dark ring of the other. The dip between the two peaks is then 0.735 of the peak (Figure 8(c)), and the eye can tell the two apart [2]. This is the Rayleigh criterion.

$$
\Delta\theta_{\min} = 1.22\,\frac{\lambda}{D}
$$

Rayleigh himself called this criterion a convenient convention [1]. The Sparrow criterion, which takes as the boundary the spacing at which the dip just disappears, is also used [1]. Either way, the core is the same: resolution is set by $\lambda/D$, and to see smaller things you must enlarge the aperture or shorten the wavelength. Here are typical values.

| Instrument | Diameter D | Wavelength | Smallest resolvable angle |
|---|---|---|---|
| Human eye (pupil) | 2 mm | 550 nm | about 1′ |
| Hubble Space Telescope | 2.4 m | 550 nm | 0.058″ |
| Palomar telescope | 5 m | 550 nm | 0.028″ |

The values for the eye and the Palomar telescope come from the textbooks [1, 2]; Hubble's was calculated with the same formula. For microscopes we speak of distance on the specimen instead of angle. The smallest separation at which two incoherently emitting points can be resolved is $0.61\lambda/\mathrm{NA}$ [2]. Here $\mathrm{NA}$ is the numerical aperture of the objective. With green light and an NA 0.95 objective, that is about 0.35 µm. If the specimen is illuminated coherently, the coefficient changes [2]. That difference is also covered in Part 4.

The Strehl ratio and Maréchal criterion of [Geometrical Optics Foundations 4](/en/posts/geometric-optics-aberrations/) take this pattern as their reference. They set the central intensity of the aberration-free Airy disk to 1 and measured how much aberrations cut it. "Diffraction-limited" means that the size of the image is set by this disk, not by aberrations.

The Airy disk is a cross-section in the focal plane. How does the light spread in front of and behind the focus? Figure 8(d) calculates the region near the focus of an aberration-free lens, seen from the side. Following Born and Wolf, the axial distance $z$ and lateral distance $r$ were converted into the optical coordinates $u = k\,\mathrm{NA}^2 z$ and $v = k\,\mathrm{NA}\,r$ [2]. In geometrical optics the light would follow the cone (dashed lines), gather at one point and spread again. Real light forms a bright region near the focus that is shaped like a cylinder, thicker than the cone. The on-axis intensity is $[\sin(u/4)/(u/4)]^2$, and the first zero is at $u = \pm 4\pi$, that is, $z = \pm 2\lambda/\mathrm{NA}^2$. Accepting a 20 % drop in intensity, the depth of focus is [2]

$$
\Delta z \approx \pm 2\lambda N^2
$$

Here $N = f/D$ is the f-number. In other words, the depth of focus is proportional to the square of the f-number. An f/10 lens with 500 nm light gives ±0.1 mm; opening it to f/2 gives ±4 µm, 25 times shallower. The photographer's rule of thumb that stopping down increases depth of field comes out of wave optics this way. Also hidden in this figure is the phenomenon in which the phase passing through focus differs by π from that of an ideal spherical wave (the Gouy phase) [2]. That story is told in Part 5 on Gaussian beams.

## 6. Many slits: the diffraction grating

Hold a CD in sunlight and you see a rainbow (Figure 9). The information-bearing track on a CD's surface is wound in a tight spiral. Neighbouring turns of the track act like many slits placed side by side. What changes when there are many slits?

<img src="/assets/img/posts/wave-optics-diffraction/ext-cd-rainbow.jpg" alt="Rainbow-coloured rings made by light reflected from a CD and split by wavelength" width="460">
_Figure 9. Diffraction from the surface of a CD. The closely spaced track acts as a reflection grating and sends white light to different angles by wavelength. (Source: Hapecko, [CD rainbow](https://commons.wikimedia.org/wiki/File:CD_rainbow.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

Suppose $N$ slits of width $b$ are placed at spacing $d$. The far field is the field of a single slit multiplied by the sum of the phase factors coming from the positions of the $N$ slits. In terms of intensity [1]:

$$
I(\theta) = I_0\left(\frac{\sin\beta}{\beta}\right)^2\left(\frac{\sin N\alpha}{N\sin\alpha}\right)^2,\qquad \alpha = \frac{\pi d\sin\theta}{\lambda}
$$

Here $\beta$ is the single-slit phase of Section 4. The first factor is the single-slit pattern and sets the overall envelope. The second factor has the same form as the interference of $M$ equal-intensity beams in [Part 1, Section 6](/en/posts/wave-optics-interference/). In other words, the pattern of a grating is "the diffraction of one slit × the interference of $N$ slits". There is no physical boundary between diffraction and interference; the names are only a convention depending on what and how many things are added [1].

<img src="/assets/img/posts/wave-optics-diffraction/en/fig10-grating.png" alt="Far-field intensity of 2, 5 and 20 slits, and the sodium D lines resolved by two gratings" width="780">
_Figure 10. (a) Far-field intensity of N slits of width 0.3d, drawn so that the principal peaks equal 1. The dashes are the single-slit envelope. The more slits, the thinner the principal peaks at the same positions, with N − 2 small peaks in between. (b) The sodium D lines (589.0 and 589.6 nm) entering the second diffraction order of a 600 lines/mm grating. The dotted lines are each line and the solid line their sum. With 400 lines they appear as one lump; with 600 lines they split in two._

In Figure 10(a) the positions of the principal peaks do not depend on the number of slits. They occur where the path differences of the light from all slits are integer multiples of the wavelength, so, including obliquely incident light, the condition is

$$
d(\sin\theta_m - \sin\theta_i) = m\lambda
$$

Here $\theta_i$ is the angle of incidence and $\theta_m$ the angle of the $m$-th diffraction order. This is the grating equation. The intensity of the principal peaks grows as $N^2$ and their width narrows as $1/N$. In a direction off by 1/100 of a wavelength per line, the 50th line is off by half a wavelength and cancels the first line's light [1]. The more lines, the smaller the deviation that causes cancellation, so the peaks are sharp.

Sharp peaks become the power to separate wavelengths. Two lights of slightly different wavelength are bent at slightly different angles. Taking as the Rayleigh criterion the case where one wavelength's peak falls on the other's first zero, the wavelength difference a grating can separate is [1, 2]

$$
\frac{\lambda}{\Delta\lambda} = mN
$$

In other words, the resolving power is the product of the diffraction order and the number of illuminated lines. The sodium D lines are 0.597 nm apart, so $\lambda/\Delta\lambda \approx 987$. In the second order, 494 lines are enough [1]. Figure 10(b) shows this. With 400 lines ($mN$ = 800) the two lines form one lump; with 600 lines ($mN$ = 1200) they split in two. A grating 15 cm long with 90,000 lines has $mN$ = 180,000 in the second order and resolves two lines 0.003 nm apart near 540 nm [1]. Where the Fabry–Pérot of Part 1 overlapped beams bouncing hundreds of times off highly reflective mirrors, the grating does the same with tens of thousands of slits side by side. The Fourier transform spectroscopy of Part 2 obtained its resolution by scanning a long path difference. All three are alike in that the resolution is set by how long a path difference the interfering light spans. For a grating that path difference is the one between light from its two ends, $Nd(\sin\theta_m - \sin\theta_i)$. On the way in and on the way out it cannot exceed the grating width $W$, so this value cannot be longer than $2W$, and the resolving power cannot exceed $2W/\lambda$ [2].

Gratings have their price too. Higher orders give better resolution, but the $m$-th order of wavelength $\lambda$ and the $(m+1)$-th order of wavelength $\lambda m/(m+1)$ overlap at the same angle. The wavelength range free of overlap (the free spectral range) is $\lambda/m$ [1]. And most of the light goes into the useless zeroth order. Calculating the efficiency of thin gratings, an amplitude grating whose transmittance varies sinusoidally sends at most 6.25 % into the first order, while a phase grating that varies only the phase sinusoidally sends at most 33.8 % [3]. Practical gratings therefore cut the grooves into a sawtooth profile to steer light into the desired order. This is called blazing. The idea of pulling energy out of the zeroth order came from Rayleigh in 1888, and Wood (R. W. Wood) actually ruled such gratings in 1910 [1].

## 7. How to compute diffraction: split light into plane waves

So far we have viewed diffraction as "the sum of spherical waves from each point of the aperture". The same calculation can be viewed the other way round: split the field just behind the aperture into a sum of plane waves travelling in different directions. Why is this view needed? It shows at a glance why patterns smaller than a wavelength cannot travel far, and the standard way to compute diffraction on a computer also comes from it.

Consider a plane wave passing a pattern of period $\Lambda$. In the picture of Section 1, light goes in the direction in which the light from two points one period apart differs in path by one wavelength. That direction is set by

$$
\sin\theta = \frac{\lambda}{\Lambda}
$$

In other words, a component of spatial frequency $1/\Lambda$ is a plane wave travelling at angle $\theta$. This is the grating equation of Section 6 with $m = 1$. Any field can be written by Fourier transform as a sum of spatial-frequency components, so any field is a sum of plane waves in different directions. This is called the angular spectrum [3, 4]. Free space spreads these plane waves out by direction. It is like a prism spreading light by colour, and Saleh and Teich liken it to a "spatial prism" [4].

But if $\Lambda < \lambda$, then $\sin\theta > 1$ and no angle fits. Where does such a component go? Figure 11 gives the answer.

<img src="/assets/img/posts/wave-optics-diffraction/en/fig11-angular-spectrum.png" alt="A tilted plane wave behind a pattern of period 1.5λ, an evanescent wave behind a pattern of period 0.8λ, and the transfer function of free space" width="780">
_Figure 11. (a) The field (real part) of light after a phase pattern of period 1.5λ. It becomes a single plane wave travelling at 41.8°. (b) Light after a pattern of period 0.8λ. It cannot travel and decays by 1/e every 0.21λ. (c) The factor by which the amplitude of each spatial-frequency component is multiplied after travelling a distance z through free space. Inside the grey band ($\lvert f\lambda\rvert < 1$) components pass unchanged with magnitude 1; outside, they vanish within a fraction of a wavelength._

A component whose period is shorter than the wavelength becomes an evanescent wave that stays attached to the pattern and decays exponentially (Figure 11(b)). Travelling through free space amounts to multiplying each component by the transfer function [3]

$$
H(f_x, f_y) = \exp\!\left(i 2\pi z\sqrt{\frac{1}{\lambda^2} - f_x^2 - f_y^2}\right)
$$

Here $f_x, f_y$ are spatial frequencies. Components with a positive quantity under the square root change only in phase and keep their magnitude. For a negative quantity the exponent becomes real and the magnitude shrinks (Figure 11(c)). In other words, free space is a filter that passes only spatial frequencies up to $1/\lambda$. Even if an object has patterns with periods shorter than the wavelength, that information vanishes before travelling even one wavelength from the object. That is why a lens far from the object cannot see periods shorter than the wavelength. This is the fundamental limit of resolution [3]. The $\lambda/D$ limit of Section 5 arises because the lens accepts only some of these plane waves, and this limit remains even if the lens accepts every direction. Getting past it requires bringing a probe closer to the object than a wavelength to pick up the evanescent waves directly. Near-field microscopes do exactly that [3, 4].

Treating $f$ as small in the transfer function's square root and keeping the binomial expansion to the second term gives $e^{ikz}\exp[-i\pi\lambda z(f_x^2+f_y^2)]$. This is the same as the Fresnel approximation of Section 4 [3]. In other words, the Fresnel approximation is a paraxial approximation that handles only small angles. So as long as the angles are small, it works well even close to the aperture, regardless of how large $N_F$ is. For a 1 cm radius aperture, a 1 m distance and 0.5 µm light, $N_F$ = 200, yet the error of the Fresnel approximation is negligible [4].

This view is itself the computational method. Split the aperture field with an FFT, multiply by $H$, and transform back with an inverse FFT. Figures 2, 4 and 11 of this post and the probe of [Ptychography 3](/en/posts/ptycho-overlap-pie-epie/) were all computed this way. The practical catch is sampling. Two points are worth noting.

- **Grid spacing.** To hold all the propagating plane waves, the pixel spacing must be $\lambda/2$ or less. With a coarser grid, the large-angle components are cut off.
- **Propagation distance.** On a grid of $N$ pixels with spacing $\Delta x$, the phase of the transfer function must not jump by more than π between neighbouring samples, which requires roughly $z \le N\Delta x^2/\lambda$. Propagating further makes light that crossed the periodic boundary come back from the other side and create spurious patterns. In that case widen the grid, or switch to a far-field method (the single-FFT Fresnel approximation, with output spacing $\lambda z/N\Delta x$).

In summary: near and at large angles, use the angular spectrum method; far and at small angles, the Fresnel FFT; very far, the Fraunhofer Fourier transform. Two numbers, the Fresnel number $N_F$ and the largest angle, are enough to choose.

## Summary

- Diffraction is the interference of countless secondary sources. The Huygens–Fresnel principle adds the spherical waves from each point of the aperture, keeping their phases. The scalar theory holds when the aperture is much larger than the wavelength.
- A single number, the Fresnel number $N_F = a^2/\lambda z$, sets the regime. Large means a shadow, near 1 a Fresnel pattern, small a Fraunhofer pattern. Diffraction is rarely seen in daily life because $N_F$ is large and light sources are extended.
- Seen through Fresnel zones, the axis behind a circular hole is four times brighter when an odd number of zones is open and dark when an even number is. The centre of a disk's shadow is as bright as without the disk (the Arago spot). Keeping only the odd zones makes a lens called a zone plate.
- At a straight edge in coherent light the intensity at the geometric boundary is 1/4, and the fringe width is proportional to $\sqrt{\lambda z}$. In a coherent image the edge is taken at 25 %, not 50 %.
- The far pattern is the Fourier transform of the aperture function. It spreads widely in the narrow direction, and straight edges send spikes perpendicular to themselves.
- A circular aperture makes the Airy pattern. The first dark ring is at $1.22\lambda/D$ and encloses 84 % of the energy. This is the yardstick for the resolution of telescopes and microscopes (the Rayleigh criterion). The depth of focus is $\pm 2\lambda N^2$.
- A grating is the product of one slit's diffraction and $N$ slits' interference. Its resolving power is $mN$, and blazing gathers light into the desired order.
- Every field is a sum of plane waves. A component of period $\Lambda$ travels at $\sin\theta = \lambda/\Lambda$, and components with $\Lambda < \lambda$ die out as evanescent waves. Free space is a filter with bandwidth $1/\lambda$.

Section 4 said that a convex lens pulls the pattern at infinity into its focal plane. Then the focal plane of a lens holds the Fourier transform of the object. If image formation by a lens is viewed as "splitting the object into plane waves, filtering some of them through the lens aperture and gathering them again", the whole of imaging can be written as a single filter. Why coherent and incoherent illumination produce different images, and what the MTF measures, is the subject of the next part: Fourier imaging.

## References

1. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Chapter 10 — diffraction. Section 10.1: the water-wave and sound analogies, the Fraunhofer condition, the Fresnel number. Section 10.2: the slit, many slits, the circular aperture and the Airy pattern, resolution (Rayleigh and Sparrow criteria, the eye and telescopes), the diffraction grating (resolving power, free spectral range, blazing, the sodium D-line example). Section 10.3: Fresnel zones, the Arago spot and the 1818 prize competition, the zone plate, Fresnel integrals and the straight edge, Babinet's principle.
2. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Chapter 8 — elements of the theory of diffraction. Section 8.1: history. Section 8.2: the Huygens–Fresnel principle and Fresnel zones. Section 8.3: Kirchhoff's theory, Babinet's principle. Section 8.5: Fraunhofer patterns (Table 8.1 slit, Table 8.2 circular aperture, encircled energy). Section 8.6: resolving power of gratings, resolution of telescopes and microscopes. Section 8.7: the straight edge. Section 8.8: the three-dimensional intensity distribution near focus and the depth of focus.
3. J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed. (McGraw-Hill, 1996), Chapter 3 — scalar diffraction theory (Kirchhoff and Rayleigh–Sommerfeld, interpretation of the Huygens–Fresnel principle, the angular spectrum and evanescent waves). Chapter 4 — the Fresnel and Fraunhofer approximations, the circular aperture, diffraction efficiencies of thin sinusoidal amplitude and phase gratings.
4. B. E. A. Saleh and M. C. Teich, *Fundamentals of Photonics*, 3rd ed. (Wiley, 2019), Sections 4.1–4.3 — spatial frequencies and plane waves, the transfer function of free space and the validity of the Fresnel approximation, far-field patterns and the optical Fourier transform, diffraction patterns as a function of the Fresnel number.
5. G. R. Fowles, *Introduction to Modern Optics*, 2nd ed. (Holt, Rinehart and Winston, 1975), Chapter 5 — the Fresnel–Kirchhoff formula, the zone plate, the Cornu spiral, resolving power of gratings.
6. NASA, ESA, CSA, L. Hustak and J. DePasquale (STScI), "Webb's Diffraction Spikes" (2022). [https://www.flickr.com/photos/nasawebbtelescope/52201717519/](https://www.flickr.com/photos/nasawebbtelescope/52201717519/)

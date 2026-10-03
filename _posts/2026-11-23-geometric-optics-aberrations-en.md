---
layout: post
date: 2026-11-23 20:00:00 +0900
title: "Geometrical Optics Foundations 4 — Why Can't a Lens Image a Point as a Point? Seidel and Wave Aberrations"
lang: en
lang-exclusive: ["en"]
permalink: /posts/geometric-optics-aberrations/
page_id: geometric-optics-aberrations
categories: [Optics, Geometrical Optics]
tags: [geometric-optics, aberration, seidel, zernike, strehl-ratio]
description: "Aberrations are the higher-order terms paraxial optics throws away. This post views the five Seidel aberrations and chromatic aberration as wavefront errors, and asks how small they must be."
math: true
---

In April 1990 the Hubble Space Telescope reached orbit. Its 2.4 m primary mirror would photograph stars without the blurring of the atmosphere. The first images were blurred anyway. A star did form a bright disk of the expected size, but that disk held only 12% of the light, against an expected 70%. About 70% of the light was spread into a hazy halo some 1.5 arcseconds across [1]. The cause was the shape of the mirror. A lens spacing in the instrument used to test the mirror's figure was off by 1.3 mm, and the mirror polished to satisfy that instrument came out too flat near its edge. The error was ten times the design tolerance [6]. Light reflected from the edge came to focus 38 mm beyond light reflected from the centre [1]. The left panel of Figure 1 is the galaxy M100 as imaged at the time.

<img src="/assets/img/posts/geometric-optics-aberrations/ext-hubble-m100.jpg" alt="The core of galaxy M100 imaged by three generations of Hubble cameras" width="760">
_Figure 1. The core of galaxy M100 imaged by the Hubble Space Telescope. Left: the original camera (WFPC1), 1993. Centre: WFPC2, installed in December 1993 with optics that cancel the mirror's flaw. Right: WFC3, installed in 2009, image taken in 2018. The mirror is the same in all three; only the optics inside the camera changed. (Source: NASA, ESA, STScI and Judy Schmidt, [PIA22913](https://commons.wikimedia.org/wiki/File:PIA22913-HubbleSpaceTelescope-ComparisonOfCameraImages-20181204.jpg), public domain)_

This episode holds every question of this post. Why did a mirror that was only "slightly" wrong ruin the image so badly? Why did the damage take the form of a disk with a halo? How could the telescope be repaired with a few small mirrors inside a camera, without replacing the primary? And conversely, how small must an error be before it no longer counts as an error?

So far this series has dealt with ideal lenses. The matrices of [Part 2](/en/posts/geometric-optics-paraxial-abcd/) were a paraxial approximation valid near the axis, and [Part 3](/en/posts/geometric-optics-stops-pupils-etendue/) added stops and pupils on top. This post restores the terms the paraxial approximation dropped. Those terms are aberrations. The marginal and chief rays defined in Part 3 become the two reference rays for computing them, and aberrations are redefined as the difference between the ideal sphere and the actual wavefront at the exit pupil.

## 1. What paraxial optics throws away

Part 2 began by straightening Snell's law: the sine of an angle was replaced by the angle itself, $\sin\theta \approx \theta$. Keeping one more term of the series,

$$ \sin\theta = \theta - \frac{\theta^3}{6} + \frac{\theta^5}{120} - \cdots $$

gives what is called third-order theory [1]. In the 1850s L. von Seidel showed that the difference between paraxial and third-order theory falls into five types. These five are called the Seidel aberrations or the primary aberrations [1, 2]. Whatever remains between an exact ray trace and third-order theory is shared among aberrations of higher order [1].

Aberrations became a practical problem with photography. Before Daguerre's process of 1839, optical design meant mostly telescope objectives with narrow fields. A camera needs a wide field and a large aperture at the same time. Around 1840 J. Petzval, by adding higher-order terms to the paraxial formulas, built a portrait lens faster than any lens of its day [2]. Modern design software does essentially the same thing. It adjusts shapes, thicknesses, glasses, spacings and the stop position so that the aberration introduced at one surface is cancelled at another [1].

Aberrations divide into two kinds [1]. Monochromatic aberrations appear even with light of a single wavelength; chromatic aberrations arise because the refractive index depends on wavelength. Of the five monochromatic aberrations, spherical aberration, coma and astigmatism blur the image, while field curvature and distortion keep it sharp but warp its shape or position [1, 2]. Sections 2 and 3 build a framework for all five at once, Sections 4–7 take them one by one, and Section 8 treats chromatic aberration.

## 2. Ray aberrations and wave aberrations

There are two ways to measure an aberration. One is to measure how far a ray misses in the image plane. The distance from the paraxial image point to where the actual ray lands is the ray aberration [2]. The longitudinal spherical aberration of Part 1 is an example: the edge ray of a 200 mm diameter lens crossed the axis 15.84 mm ahead of the paraxial focus.

The other is to look at the wavefront. In an aberration-free system, light from an object point leaves the exit pupil as a perfect spherical wave centred on the image point. As Part 3 showed, the exit pupil is the window through which light enters image space. The sphere centred on the image point and passing through the centre of the exit pupil is called the Gaussian reference sphere [2, 3]. The optical path length by which the actual wavefront departs from this sphere is the wave aberration $W$ [1, 2, 3].

A finish line makes a good picture. The image point is the finish, and the wavefront is the line of runners who all started at the same moment. For every runner to arrive together, all must now be the same distance from the finish, on a circle centred on it. That circle is the reference sphere. A runner on the edge who is ahead of the circle arrives early, or arrives somewhere other than the finish. $W$ is how far ahead of or behind the circle each runner stands, measured as optical path. In this post $W > 0$ when the actual wavefront is ahead of the reference sphere, towards the image. Born and Wolf use the opposite sign [2]. Conventions differ between books, so check before copying a formula.

The two descriptions are two faces of one thing. Part 1 showed that rays are normals to the wavefront. Wherever the wavefront is tilted relative to the reference sphere, its normal — the ray — misses the image point (Figure 2(a)). The miss is proportional to the slope of the wavefront. In normalized pupil coordinates $\rho$ (equal to 1 at the rim of the exit pupil),

$$ \varepsilon_y \approx -\frac{1}{n'\,\mathrm{NA}}\,\frac{\partial W}{\partial \rho_y} $$

where $n'$ is the image-space index and NA the image-side numerical aperture [2]. Born and Wolf give the exact relation in terms of the radius of the reference sphere; this is that relation with the radius held constant [2].

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig2-wave-and-ray-aberration.png" alt="Relation between wave and ray aberration: a sketch and a real-ray calculation for a plano-convex lens" width="780">
_Figure 2. (a) Sketch (exaggerated). At the edge the actual wavefront (blue) is ahead of the reference sphere (grey dashed) towards the image. Its normals — the rays — cross the axis ahead of the image point P*, more so towards the edge. (b) Real-ray trace of an on-axis point through an f = 100 mm, F/5 N-BK7 plano-convex lens. Blue: wave aberration. Red: the height by which each ray misses at the Gaussian image plane. Black circles: the miss computed by differentiating the wave aberration._

Figure 2(b) checks this on a real lens. For an f = 100 mm N-BK7 plano-convex lens at F/5, the wavefront at the rim of the exit pupil is 4.7 waves ahead of the reference sphere (wavelength 587.6 nm). The marginal ray misses the Gaussian image point by 112 µm. The miss obtained by differentiating the wave aberration (black circles) cannot be told apart from the real-ray result (red line).

This relation also shows how an aberration grows with aperture. The wave aberration of spherical aberration goes as $\rho^4$. Differentiating, the transverse ray aberration goes as $\rho^3$, and the longitudinal aberration — the shift of the axis crossing — is that divided by the ray slope ($\propto\rho$), so it goes as $\rho^2$. The statement in Part 1 that longitudinal spherical aberration scales with the square of the ray height (half the height, a quarter of the aberration) is the end of this chain. The same aberration grows as the fourth power in the wavefront, the third power in the image plane, and the second power along the axis.

Why use the wavefront as the reference, when ray aberrations are easier to see? There are two reasons. First, the wave aberration holds every ray in a single function $W(\rho)$; ray aberrations follow by differentiation. The information sits in one place. Second, once aberrations shrink to a fraction of a wavelength, the ray picture breaks down and diffraction determines the image. The starting point for computing that image is the wavefront in the exit pupil [2, 3]. Goodman describes this as inserting a phase plate in the exit pupil: an aberration-free spherical wave with a phase $kW$ added is called the generalized pupil function, and diffraction calculations start from it [3]. Section 10 follows this route to answer "how small is small enough".

## 3. Why exactly five?

What is $W$ a function of? Fix an off-axis object point; its rays are then distinguished by the point $(\rho, \theta)$ where they cross the exit pupil. Writing the normalized distance of the object point from the axis as $h$, we have $W = W(h, \rho, \theta)$, with $\theta$ measured in the pupil from the direction of the object point (the meridional plane).

Rotational symmetry of the lens restricts this function strongly [2]. Turning the whole lens about its axis changes nothing, so $W$ can only be built from quantities that do not change when both the object-point vector and the pupil-point vector rotate. Two vectors give exactly three such quantities: the squared length of the object vector $h^2$, the squared length of the pupil vector $\rho^2$, and their dot product $h\rho\cos\theta$ [2]. All a lens "knows" is how far the object point is from the axis, where the ray crosses the pupil, and the direction between the two.

Expand $W$ in these three quantities. The constant term is zero by choice of reference. The second-order terms, $\rho^2$ and $h\rho\cos\theta$, shift the focus and shift the image sideways; they belong to paraxial theory. Taking the paraxial image point as reference makes them zero [2]. The lowest surviving order is the fourth, and there are exactly five fourth-order combinations that involve pupil coordinates [2]. A pure $h^4$ term, with no pupil coordinate, adds the same value to every ray and has no effect on the image [2].

| Combination | Wave aberration term | Name | In the pupil | In the field |
|---|---|---|---|---|
| $(\rho^2)^2$ | $W_{040}\,\rho^4$ | spherical aberration | $\rho^4$ | independent of field |
| $\rho^2 \cdot h\rho\cos\theta$ | $W_{131}\,h\rho^3\cos\theta$ | coma | $\rho^3$ | $h$ |
| $(h\rho\cos\theta)^2$ | $W_{222}\,h^2\rho^2\cos^2\theta$ | astigmatism | $\rho^2$ | $h^2$ |
| $h^2 \cdot \rho^2$ | $W_{220}\,h^2\rho^2$ | field curvature | $\rho^2$ | $h^2$ |
| $h^2 \cdot h\rho\cos\theta$ | $W_{311}\,h^3\rho\cos\theta$ | distortion | $\rho$ | $h^3$ |

The three subscripts are the powers of $h$, $\rho$ and $\cos\theta$. Born and Wolf write the same five terms with coefficients $B$, $C$, $D$, $E$, $F$ [2]. This one table already gives the character of all five. Only spherical aberration lacks $h$, so only it appears for an on-axis object point. Coma grows linearly with the field, astigmatism and field curvature quadratically, distortion as the cube. With pupil diameter the order reverses: spherical aberration grows fastest and distortion slowest. That is why stopping down shrinks spherical aberration dramatically but leaves distortion almost untouched.

Figure 3 computes the wavefront and image for each of the five terms with its coefficient set to one wave. The middle row shows where rays land, using the derivative relation of Section 2; the bottom row is the diffraction image of the same wavefront.

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig3-seidel-gallery.png" alt="Wavefront, ray hits and diffraction image for the five Seidel aberrations" width="780">
_Figure 3. The five Seidel aberrations, each with a coefficient of one wave. Top: wave aberration in the exit pupil (red ahead, blue behind). Middle: ray hits in the Gaussian image plane. Bottom: diffraction image of the same wavefront (square root of intensity). The field point is upward, the cross marks the paraxial image point, and the number is the peak intensity relative to the aberration-free peak. Axes are in units of λ/NA._

It is worth memorizing the five shapes. Spherical aberration surrounds a central spot with a hazy halo — the left panel of Figure 1. Coma is a comet trailing a tail to one side. Astigmatism stretches the point into a line. Field curvature gives a defocused disk in this plane, meaning that this field point comes to focus in another plane. Distortion does not blur at all; it only moves the image, and the peak intensity stays at 1.00.

## 4. Spherical aberration

Spherical aberration means that rays through the edge of a lens focus somewhere other than rays through its centre [1]. In a positive lens the edge rays are bent too strongly and focus ahead of the paraxial focus; this is called undercorrected spherical aberration. A negative lens does the opposite [1]. Part 1 showed that the envelope of this ray bundle is a caustic.

### Where to focus

With spherical aberration there is no single "focus". Figure 4(a) shows the region around focus of a plano-convex lens opened up to F/2.5. The paraxial focus and the marginal focus are 4.57 mm apart. Where in between should the screen go for the best image?

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig4-spherical-focus.png" alt="Rays near focus of a lens with spherical aberration, and blur and Strehl ratio versus focus position" width="780">
_Figure 4. (a) Real rays near the focus of an f = 100 mm, F/2.5 plano-convex lens (colour = entrance height, yellow at the edge). The further out a ray enters, the earlier it crosses the axis. The circle of least confusion, where the blur diameter is smallest, lies three quarters of the way from the paraxial towards the marginal focus. (b) For pure third-order spherical aberration ($W_{040} = 0.5\lambda$): geometric blur radius (orange) and central image intensity (blue) versus focus position. The intensity peaks halfway._

The geometrical answer is the plane of smallest blur, the circle of least confusion [1]. With pure third-order spherical aberration it lies exactly three quarters of the way from the paraxial towards the marginal focus. The real lens of Figure 4(a) puts it at 0.76.

Diffraction changes the answer. Born and Wolf showed that for small spherical aberration the point of maximum central intensity — the diffraction focus — lies exactly midway between the paraxial and marginal foci [2]. Figure 4(b) shows the difference. With $W_{040} = 0.5\lambda$ the central intensity peaks at the halfway point (0.95) and drops to 0.77 at the geometrical circle of least confusion. Making the edge of the ray bundle as small as possible and putting as much light as possible into the centre are different goals. When aberrations are of the order of a wavelength, the second goal is the right one.

The halfway point has a simple origin. Shifting focus by $\delta z$ adds a $\rho^2$ term to the wavefront. The combination $\rho^4 - \rho^2$ is spread more evenly around zero across the pupil than $\rho^4$, and its mean square is a quarter as large. Section 10 shows that this is the idea behind the Zernike polynomials.

### Reducing it by bending the lens

A given focal length can be realized with many lens shapes. With radii $R_1$ and $R_2$, define the shape factor $q = (R_2 + R_1)/(R_2 - R_1)$: $q = 0$ is a symmetric biconvex lens, $q = +1$ a plano-convex lens with its curved side towards the object, and $q = -1$ the same lens with its flat side towards the object. Changing $q$ while keeping the focal length is called lens bending.

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig5-lens-bending.png" alt="Spherical aberration and coma versus lens shape factor" width="700">
_Figure 5. Spherical aberration (middle) and coma at a 5° field (bottom) of an f = 100 mm, F/5 N-BK7 lens bent into different shapes. Lines: third-order theory from two paraxial rays (Section 9). Squares: longitudinal spherical aberration measured with real rays. The stop is at the lens and the object is at infinity._

Hecht points out that simply turning a plano-convex lens around reduces its spherical aberration considerably [1]. Figure 5 puts a number on it. With the flat side towards the object ($q = -1$) the longitudinal spherical aberration is 4.38 mm; turned around so the curved side faces the object ($q = +1$) it is 1.10 mm, about a quarter. The intuition comes from prisms. The edge of a lens bends rays like a small prism, and a prism deviates least when the entering and leaving angles are similar [1]. A collimated beam meeting the curved side first shares the bending between the two surfaces; meeting the flat side first, the rear surface must do all of it.

Spherical aberration is smallest near $q \approx 0.7$: curved side towards the object, rear side nearly flat [1]. The thin-lens formula gives $q = 2(n^2 - 1)/(n + 2) = 0.74$; for the real lens, 4 mm thick, it is 0.70. Even the minimum is not zero. A singlet with two spherical surfaces cannot remove the spherical aberration for an object at infinity. Removing it takes a pair of positive and negative lenses [1], or an aspheric surface.

### Points without spherical aberration

There is an exception. Even a single spherical surface has a pair of object and image points that it images free of spherical aberration. These are called aplanatic points [1]. Oil-immersion objectives exploit them. The specimen sits in oil with the same index as glass, which optically joins it to the first lens, and the curved surface of that lens is made so that the specimen sits at its aplanatic point. Light from the specimen leaves the first lens as a narrower bundle with no spherical aberration. A meniscus lens working on the same principle follows and narrows the bundle once more [1]. It is the classical way to accept a large numerical aperture without spherical aberration.

### The microscope cover slip

Spherical aberration is not only a property of lenses. Murphy and Davidson describe how microscope users create it without noticing [5]. High-power dry objectives are corrected for viewing through a 0.17 mm cover slip (grade #1.5). A cover slip is a flat glass plate in a converging beam, and a flat plate lengthens the optical path more for oblique rays, which produces spherical aberration. The objective carries the opposite of the aberration of a 0.17 mm plate, so any difference in thickness is left uncorrected.

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig6-coverslip.png" alt="Wave aberration and Strehl ratio from a cover slip thickness error versus objective numerical aperture" width="760">
_Figure 6. Spherical aberration left in a dry objective when the cover slip (n = 1.523) departs from its design thickness (wavelength 550 nm), assuming an objective that obeys the sine condition with a uniformly filled pupil, after refocusing. (a) Root-mean-square (RMS) wave aberration. (b) Central image intensity (Strehl ratio, Section 10). The dotted vertical line marks NA 0.4._

Figure 6 computes the size of the effect. A 10 µm thickness error leaves only 0.002 waves (RMS) in an NA 0.4 objective — negligible, consistent with Murphy and Davidson's remark that cover slip thickness hardly matters below NA 0.4 [5]. At NA 0.85, however, the same 10 µm lowers the central intensity to 0.80, and at NA 0.95 to 0.44, because the wave aberration grows faster than the fourth power of the NA. That is why high-NA dry objectives carry a correction collar, which changes a lens spacing as it is turned [5]. Imaging deep inside a specimen produces spherical aberration in the same way.

## 5. Coma

Coma appears for object points even slightly off the axis [1]. The name comes from the Latin for a comet's tail: a star is recorded as a comet trailing a tail (Figure 7). Its asymmetry makes it perhaps the most objectionable aberration [1], and in telescopes it must be removed because it spoils accurate measurement of star positions [2].

<img src="/assets/img/posts/geometric-optics-aberrations/ext-coma-corrector.jpg" alt="A star near the edge of the field before and after a coma corrector" width="500">
_Figure 7. A star near the edge of the field of an F/3.9 Newtonian reflector. Left: without a corrector, showing the comatic tail. Right: with a coma corrector. (Source: Rawastrodata, [Baader Rowe Coma Corrector Comparison](https://commons.wikimedia.org/wiki/File:Baader_Rowe_Coma_Corrector_Comparison.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

### Each ring of the aperture has its own magnification

The easiest route to coma is the principal surface of Part 2. There we saw that the principal "plane" is really a curved surface: in the lens of Part 1, a ray at a height of 100 mm bent at a point 1.1 mm away from the paraxial principal plane. Hecht identifies this curvature as the origin of coma [1]. With a curved principal surface, the centre and the edge of a lens behave like lenses of different focal length, and different focal lengths mean different magnifications. For an on-axis object point this does no harm, since every image point lies on the axis. For an off-axis point, the image formed by the central ring and the image formed by the outer ring land at different heights [1].

More precisely, rays crossing the pupil on one ring (radius $\rho$) trace out a circle in the image plane. Its radius is proportional to $\rho^2$, and its centre is displaced from the paraxial image point by twice that radius. As a ray goes once round the ring, its image goes twice round the circle [2]. Larger rings give larger circles further out, so all the circles touch two lines through the paraxial image point that open at 60° (Figure 8(a)) [1, 2]. The comet's head is the paraxial image point, where the central rings collect, and its tail is the fan of overlapping circles from the outer rings. A little over half of the light falls in the roughly triangular region near the apex [1].

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig8-coma.png" alt="Ray hits for coma: pure third-order coma and a real-ray paraboloid" width="760">
_Figure 8. (a) Pure third-order coma ($W_{131} = 1\lambda$). Rays through one ring of the pupil go twice round one circle in the image plane. The three dots are rays at 0°, 45° and 90° on the ring. The circles touch two lines (dashed) opening at 60°. (b) Real-ray trace of light entering an f = 300 mm, F/3 paraboloidal mirror at a field of 0.2°. Colour indicates the ring radius in the pupil._

A paraboloidal mirror brings light parallel to its axis to a point with no spherical aberration [1]. Tilt the light slightly and coma appears. Figure 8(b) shows an F/3 paraboloid with light arriving 0.2° off axis. The tail length of 22.4 µm agrees within 3% with the third-order value of 21.8 µm. The Newtonian telescope of Figure 7 uses exactly such a primary mirror. Coma grows linearly with field angle and as the inverse square of the f-number, so a fast telescope such as an F/3.9 suffers badly at the edge of its field. Hence the separate corrector lens.

### The sine condition

If coma comes from rings of the aperture having different magnifications, the cure is to make every ring's magnification the same. The optical sine theorem, found independently by Abbe and Helmholtz in 1873, gives that magnification [1]. With object- and image-space indices $n_o, n_i$, object and image heights $y_o, y_i$, and ray angles to the axis $\alpha_o, \alpha_i$,

$$ n_o\,y_o \sin\alpha_o = n_i\,y_i \sin\alpha_i $$

For the magnification $y_i / y_o$ to be the same for every ray, $\sin\alpha_o / \sin\alpha_i$ must be constant over the aperture. This is the sine condition [1]. If it holds, primary coma is absent [2]. A system free of both spherical aberration and coma is called an aplanat [1]; it is the basic design condition for microscope objectives, which accept light over large angles. The Lagrange invariant $nyu$ of Part 2 is the sine theorem in the paraxial region.

### Removing it by moving the stop

The lower panel of Figure 5 shows that lens bending can also bring coma to zero. At a 5° field coma vanishes for $q = 0.81$, close to the shape of minimum spherical aberration ($q = 0.70$) [1]. Born and Wolf give this shape for a thin lens of index 1.5 as $R_1 = 5f/9$, $R_2 = -5f$ [2]. A plano-convex lens with its curved side towards the object is therefore a good compromise that nearly minimizes spherical aberration and coma together.

The other remedy is to move the stop, a method found by W. H. Wollaston in 1812 [1]. Moving the stop makes the chief ray pass through a different part of the lens, and the bundle centred on it uses a different part of the lens. If the lens has spherical aberration, a bundle displaced to one side loses its internal symmetry, and coma appears or disappears [1]. An important rule follows. Moving the stop changes an aberration only if one of the aberrations before it in the list is present. The order is spherical aberration, coma, astigmatism, Petzval curvature, distortion [1]. Spherical aberration and Petzval curvature do not depend on the stop position, and coma depends on it only when spherical aberration is present. The calculation of Section 9 states this rule as a formula.

## 6. Astigmatism and field curvature

### Two focal lines

Further off axis, astigmatism appears [1]. Split the ray bundle from an off-axis point into two sections. The plane containing the object point and the axis is the meridional (or tangential) plane; the plane containing the chief ray and perpendicular to the meridional plane is the sagittal plane [1, 2]. For an on-axis point the two are indistinguishable, but an off-axis bundle crosses the lens obliquely, and the lens appears curved differently in the two planes — just as a round plate seen at an angle looks elliptical. Meridional rays meet the lens more obliquely, are bent more strongly, and focus closer [1].

So there are two foci. Figure 9(a) follows the image of a 10° field point of an F/10 plano-convex lens through focus. Where the meridional rays meet first, the image is a horizontal line, because the tangential focal line is perpendicular to the meridional plane [2]. At the sagittal focus, 2.9 mm further on, it becomes a vertical line. Halfway between, the image is the roundest blur, the circle of least confusion [1]. Shine light obliquely through a lens onto a piece of wire mesh and you can see this directly: depending on the screen position, only the horizontal or only the vertical wires are sharp [1].

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig9-astigmatism-field-curvature.png" alt="The spot of an astigmatic image through focus, and tangential, sagittal and Petzval image surfaces versus field" width="780">
_Figure 9. (a) The 10° image of an f = 100 mm, F/10 plano-convex lens (curved side first, stop at the lens), followed from the tangential to the sagittal focus. Up is away from the axis. (b) Tangential (T) and sagittal (S) focus positions and the Petzval surface versus field angle for the same lens. Lines: third-order theory; symbols: real rays. (c) With the lens reversed and the stop at the point where the centre of curvature of the rear surface appears through the flat face, T and S meet on the Petzval surface._

This astigmatism shares its name with the astigmatism of the eye but not its cause. The eye's astigmatism comes from a refracting surface curved differently in different directions; a lens's astigmatism arises even when the lens is perfectly rotationally symmetric [1].

### The image is not formed on a plane

Even without astigmatism one problem remains. The surface on which a lens sharply images a flat object is not a plane but a curved surface. This is Petzval field curvature [1]. For a system of thin lenses, the departure of this Petzval surface from the paraxial image plane is

$$ \Delta x = \frac{y_i^2}{2}\sum_j \frac{1}{n_j f_j} $$

where $y_i$ is the image height and $n_j$, $f_j$ are the index and focal length of each lens [1]. The formula contains neither lens shapes, nor spacings, nor the stop position [1, 2]. Once the powers and glasses are fixed, the Petzval surface cannot be changed. For a single positive lens it is a bowl curving towards the lens. In Figure 9(b) the Petzval surface of the 100 mm lens lies 1.0 mm in front of the paraxial image plane at a 10° field.

With astigmatism, the two focal lines trace two surfaces, the tangential and sagittal image surfaces. These are related to the Petzval surface in a striking way: the tangential surface is always three times as far from the Petzval surface as the sagittal surface, on the same side [1, 2]. At 10° in Figure 9(b) the tangential focus is 4.36 mm and the sagittal focus 1.45 mm in front of the Petzval surface. Remove the astigmatism and both surfaces collapse onto the Petzval surface [1].

### Schmidt's idea

Astigmatism can be removed through the stop position. Figure 9(c) is a dramatic example. Reverse the plano-convex lens so its flat side faces the object: collimated light then enters the flat face without bending and is refracted only at the curved rear surface. Now put the stop at the centre of curvature of that surface — at the point where it appears through the flat face, 31.4 mm in front of the lens. Every chief ray then passes through the centre of curvature and meets the curved surface at normal incidence. For each bundle, its own chief ray is effectively an optical axis. With no off-axis points, there is neither coma nor astigmatism. The two coefficients come out exactly zero, and the real-ray T and S foci meet on the Petzval surface.

Hecht explains that this is how the Schmidt camera works [1]. A stop at the centre of curvature of a spherical mirror removes coma and astigmatism, and the remaining spherical aberration is removed by a thin corrector plate at the stop. Petzval curvature, however, cannot be removed by any stop. Instead of flattening the image surface, the Schmidt camera bent the film to match it [1].

<img src="/assets/img/posts/geometric-optics-aberrations/ext-kepler-focal-plane.jpg" alt="The Kepler space telescope focal plane, with CCDs arranged on a curved surface" width="640">
_Figure 10. The focal plane of the Kepler space telescope. The CCD modules are arranged on a curved surface following the curved image surface of the 0.95 m Schmidt telescope, and a field flattener on each module maps the spherical image surface onto the flat CCDs [7]. (Source: NASA and Ball Aerospace, [Kepler CCD matrix](https://commons.wikimedia.org/wiki/File:Kepler_CCD_matrix.jpg), public domain)_

The Kepler space telescope, which found thousands of exoplanets, is a modern Schmidt camera (Figure 10). To cover the wide field of its 0.95 m Schmidt telescope, its CCD modules were laid out along the curved image surface, each with its own field-flattening lens [7]. Camera lenses with flat sensors instead pair positive and negative powers so that the Petzval sum $\sum 1/(n_j f_j)$ vanishes, the Petzval condition [1]. Any remaining curvature is removed by a negative lens placed close to the image plane, a field flattener; a lens right at the image plane barely affects the other aberrations [1]. The "plan" in a microscope objective's name indicates that this field curvature has been corrected [5].

## 7. Distortion

Distortion does not blur the image. Every point is imaged sharply, but in the wrong place [1, 2], because the magnification changes from the centre of the image outwards. If magnification increases outwards, the corners of a grid are pulled out into pincushion distortion; if it decreases, they are pulled in, giving barrel distortion [1]. Straight lines through the axis stay straight; only lines that miss the axis are curved [2].

In the table of Section 3, distortion enters the pupil coordinate only to the first power. That term gives every ray in the pupil the same tilt, which moves the whole image sideways. So all that matters is where a single ray — the chief ray, representing the bundle — lands. Stopping down leaves the chief ray as it was, so it does not reduce distortion [1].

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig11-distortion.png" alt="Distortion of the image of a grid for different stop positions" width="780">
_Figure 11. The image of a ±20° grid at infinity formed by an f = 100 mm biconvex lens. Blue: real rays. Grey dashed: paraxial image. The lens is the same in all three panels; only the stop moves: (a) in front, (b) at the centre of the lens, (c) behind._

In Figure 11 the lens is the same in all three cases; only the stop moves. With the stop in front, there is −4.5% barrel distortion at 20°; with the stop behind, +12.6% pincushion; with the stop at the centre of the lens, none. Hecht's explanation is intuitive [1]. With a front stop, the chief ray of an off-axis bundle passes towards the edge of the lens, the object distance measured along it becomes longer, and the magnification drops. A rear stop does the opposite. A stop at the centre of the lens makes the chief ray pass through the lens centre undeviated, so there is no distortion. A pinhole camera is distortion-free for the same reason [1].

This is why lenses are arranged symmetrically about the stop. The barrel distortion of the front group cancels the pincushion of the rear group. In a fully symmetric lens used at unit magnification, distortion, coma and lateral colour all vanish exactly [1]. Photographic lens designs are commonly arranged close to symmetric about the stop for this reason [1].

Born and Wolf class distortion and field curvature as relatively harmless: the image stays sharp and only its geometry is wrong, so the error can be removed by calculation [2]. Image-based metrology does exactly this. Camera calibration photographs a checkerboard, extracts distortion coefficients, and maps measured coordinates back through them. An image blurred by spherical aberration or coma cannot be undone in such a simple way.

## 8. Chromatic aberration

So far we have used light of one wavelength. The refractive index of glass changes with wavelength, generally increasing towards shorter wavelengths, so the same lens bends blue light more. A positive lens brings blue to focus nearest and red farthest. This difference in focus position is axial chromatic aberration [1]. Since the focal length varies with wavelength, so does the magnification, and towards the edge of the image each colour forms an image of a different size. That is lateral chromatic aberration [1]. The colour fringes along the dark edges in the lower photograph of Figure 12 are an example.

<img src="/assets/img/posts/geometric-optics-aberrations/ext-chromatic-aberration.jpg" alt="Photographs without and with colour fringing from chromatic aberration" width="423">
_Figure 12. Corner crops of photographs taken with a digital camera's built-in lens (top) and with an added wide-angle converter (bottom). The bottom image shows colour fringes along dark edges. (Source: Stan Zurek, [Chromatic aberration (comparison)](https://commons.wikimedia.org/wiki/File:Chromatic_aberration_%28comparison%29.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/))_

Chromatic aberration is much larger than the monochromatic kind. Hecht notes that the change of the monochromatic aberrations with wavelength is inconsequential unless the system is very well corrected, but that chromatic aberration itself is far more significant [1]. The red curve in Figure 13(a) shows how the focus of an f = 100 mm N-BK7 singlet moves with wavelength. The foci for the blue F line (486.1 nm) and the red C line (656.3 nm) are 1.5 mm apart — more than the 1.1 mm of longitudinal spherical aberration of the F/5 lens in Section 2. After a few experiments Newton wrongly concluded that chromatic aberration could not be removed from refracting telescopes, and turned to designing reflectors [1].

### The achromatic doublet

The way out is to pair a positive and a negative lens. The point is to keep the overall bending while cancelling the part of it that depends on colour [1]. The dispersion of a glass is expressed by its Abbe number,

$$ V_d = \frac{n_d - 1}{n_F - n_C} $$

where $n_d$, $n_F$ and $n_C$ are the indices at the Fraunhofer d line (587.6 nm), F line and C line [1]. The numerator measures bending power and the denominator its variation with colour, so a high Abbe number means low dispersion. Crown glasses such as N-BK7 have Abbe numbers around 64; flint glasses such as N-F2 around 36. For two thin lenses in contact, the F and C foci coincide when

$$ f_1 V_1 + f_2 V_2 = 0 $$

[1]. Abbe numbers are positive, so one lens must be positive and the other negative. The positive lens is made of high-Abbe crown and a weaker negative lens of low-Abbe flint. Although weak, the flint lens is highly dispersive and cancels the colour spread of the crown lens, leaving a net positive power. For a focal length of 100 mm that means a crown positive lens of $f_1 = 43$ mm and a flint negative lens of $f_2 = -76$ mm.

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig13-chromatic-focus.png" alt="Focus position versus wavelength for a singlet and an achromatic doublet" width="760">
_Figure 13. (a) Paraxial focus position versus wavelength for an N-BK7 singlet (red) and a cemented N-BK7 + N-F2 achromat (blue), both of focal length about 100 mm, relative to the d-line focus. (b) The achromat alone, enlarged. The F and C foci coincide, but wavelengths in between focus in front and wavelengths outside focus behind. The remainder is the secondary spectrum._

Figure 13(b) shows the resulting achromat. Splitting the power with the thin-lens formula and then adding thickness shifts the two lines slightly apart, so only the curvature of the last surface was readjusted. The F and C foci now coincide exactly, but the d-line focus lies 52 µm in front of them. Only two wavelengths were matched, so those in between still miss slightly. This is the secondary spectrum, here about 1/2000 of the focal length. The singlet's 1.5 mm has shrunk to 52 µm, an improvement of nearly thirty times.

Microscope objective classes follow these steps [5]. Achromats bring red (656 nm) and blue (486 nm) to a common focus and correct spherical aberration in the green (540 nm). Fluorite objectives, or semi-apochromats, include materials of unusual dispersion such as fluorite to reduce the secondary spectrum. Apochromats bring four or five colours to a common focus and correct spherical aberration at several colours; they can contain more than ten lens elements [5]. The prefix "plan" adds correction of the field curvature of Section 6.

## 9. Computing aberrations from two rays

So far we have looked at what aberrations look like. A designer needs numbers: how large the aberrations of a ten-element system are, and which surface produces them. Remarkably, third-order aberrations can be computed from just the two rays of Part 3, the paraxial marginal and chief rays [2]. Not a single real ray is traced.

The method is bookkeeping, surface by surface. Knowing the heights and angles of the two rays at a surface gives the five aberrations that surface contributes, and the aberrations of the whole system are their sums. Born and Wolf put Seidel's results in terms of two paraxial rays [2]. The form common in lens design runs as follows. At each surface, take the paraxial marginal-ray height $y$ and incoming slope $u$, the chief-ray $\bar{y}$ and $\bar{u}$, the surface curvature $c$, and the indices $n, n'$ before and after. Form the refraction invariants $A = n(u + yc)$ and $\bar{A} = n(\bar{u} + \bar{y}c)$. $A$ is the index times the angle of incidence, so by Snell's law it is the same on both sides of the surface. Then

$$ S_I = -\sum A^2 y\,\Delta\!\left(\frac{u}{n}\right),\quad S_{II} = -\sum A\bar{A} y\,\Delta\!\left(\frac{u}{n}\right),\quad S_{III} = -\sum \bar{A}^2 y\,\Delta\!\left(\frac{u}{n}\right) $$

$$ S_{IV} = -\sum H^2 c\,\Delta\!\left(\frac{1}{n}\right),\qquad S_V = \sum \frac{\bar{A}}{A}\left(S_{III} + S_{IV}\right)_{\text{surface}} $$

Here $\Delta(\cdot)$ is the value after the surface minus the value before, and $H$ is the Lagrange invariant of Part 3. The coefficients of the table in Section 3 follow from the five sums as $W_{040} = S_I/8$, $W_{131} = S_{II}/2$, $W_{222} = S_{III}/2$, $W_{220} = (S_{III} + S_{IV})/4$ and $W_{311} = S_V/2$.

There is no need to memorize these formulas; their structure is what matters. Spherical aberration depends only on the marginal ray's angle of incidence ($A$), astigmatism only on the chief ray's ($\bar{A}$), and coma on their product. The more strongly rays bend at a surface, the more aberration it produces. And $S_{IV}$ contains no ray data at all. It is the Petzval sum, fixed by curvatures and indices alone, which is why Section 6 found Petzval curvature independent of lens shape and stop position. Moving the stop changes only the chief ray's $\bar{A}$. Then $S_I$ is unchanged, $S_{II}$ changes only if $S_I$ is present, and $S_{III}$ only if $S_I$ or $S_{II}$ is present. This is the ordering rule Hecht describes, quoted in Section 5.

Figure 9(c) also follows from these formulas. With the stop at the centre of curvature of the surface, the chief ray meets the surface at normal incidence, so $\bar{A} = 0$, and $S_{II}$ and $S_{III}$ vanish together.

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig14-seidel-budget.png" alt="Ray paths of a Cooke triplet and the Seidel contributions of each surface" width="780">
_Figure 14. (a) A Cooke triplet (focal length 50 mm, F/5, 20° half-field, N-SK16 and N-F2 glasses; prescription from the examples of the open-source optical design library Optiland). Orange: marginal rays from the axial object point. Green: the 20° bundle. (b) Seidel contributions of each surface computed from two paraxial rays (colours) and their sums (black)._

Figure 14 is that ledger. The Cooke triplet, published by H. D. Taylor in 1893, is a photographic lens with the smallest number of elements in which all seven third-order aberrations — the five monochromatic and two chromatic — can essentially be made to vanish [1]. The contributions of individual surfaces are large. For distortion $S_V$, surface 2 contributes −0.35, surface 3 +0.29, surface 4 −0.26 and surface 5 +0.31 (in mm), yet the sum is −0.002, 0.5% of the largest contribution. Coma also sums to 1.5% of its largest term. The central negative lens cancels the aberrations of the positive lenses on either side with the opposite sign. If the matrices of Part 2 reduced ten lenses to a single 2×2 matrix, Seidel sums reduce their aberrations to five numbers.

Although computed from paraxial data alone, these sums agree well with real rays. Fitting a polynomial to a real-ray wavefront at small aperture and field gives third-order coefficients within about 1% of the sums. At larger aperture and field, fifth- and higher-order terms grow and the difference increases. Practical design starts from the Seidel sums and finishes with real-ray tracing [2].

## 10. How small is small enough?

Back to the opening question: how small must an aberration be to count as none? Born and Wolf note that ordinary instruments may have wave aberrations of forty or fifty wavelengths, but precision instruments such as astronomical telescopes and microscopes must reduce them to a fraction of a wavelength [2]. At that level the image is set by diffraction, not by the scatter of rays. Even an aberration-free lens cannot image a point as a point; it forms a small pattern called the Airy disk [1]. The wave optics series computes that pattern. Aberrations spoil it.

### Rayleigh's quarter wave

Rayleigh proposed a practical criterion. When the wavefront in the exit pupil departs from the reference sphere by more than a quarter wavelength (peak to valley), the image degrades noticeably [1, 2]. He showed that a quarter wave of spherical aberration reduces the central intensity of the image by about 20% [1, 2]. The criterion fits many aberrations reasonably well and is widely used as Rayleigh's quarter-wavelength rule [2].

### The Strehl ratio and the Maréchal criterion

A more quantitative measure is the central intensity of the image. The ratio of the peak intensity with aberrations to the peak intensity without them, for the same aperture, is the Strehl ratio [2, 4]. For small aberrations this ratio depends not on the type of aberration but only on the root-mean-square deviation $\sigma$ of the wavefront (the RMS wave aberration) [2, 4]:

$$ S \approx 1 - \left(\frac{2\pi\sigma}{\lambda}\right)^2 \approx \exp\!\left[-\left(\frac{2\pi\sigma}{\lambda}\right)^2\right] $$

Maréchal proposed to call a system well corrected when its Strehl ratio is at least 0.8, which is equivalent to an RMS wave aberration of at most $\lambda/14$ [2]. This is usually what "diffraction-limited" means.

<img src="/assets/img/posts/geometric-optics-aberrations/en/fig15-strehl-rms.png" alt="Balancing spherical aberration, and Strehl ratio versus RMS for different aberrations" width="780">
_Figure 15. (a) Adding a focus term $-\rho^2$ (blue) to spherical aberration $\rho^4$ (red) cuts the RMS to a quarter; the result has the shape of the Zernike spherical-aberration polynomial (dotted). (b) Exact Strehl ratio (coloured lines) for spherical aberration, coma and astigmatism, each with tilt and focus optimized, compared with Maréchal's approximation (black dashed). The dots are the tolerances Born and Wolf give for a Strehl ratio of 0.8._

Figure 15(b) checks this. Spherical aberration, coma and astigmatism produce entirely different images, yet at equal RMS their Strehl ratios fall on almost the same curve. Born and Wolf's table of tolerances gives the coefficients that yield a Strehl ratio of 0.8 as 0.94 waves for spherical aberration, 0.60 waves for coma and 0.35 waves for astigmatism [2]. Computing them, all three give an RMS of 0.070–0.071 waves, that is $\lambda/14$, and a Strehl ratio of 0.82. The tolerances only look different because the same RMS corresponds to different coefficients for different aberrations.

Return now to the cover slip of Figure 6. That 10 µm produces a Strehl ratio of 0.8 at NA 0.85 means the objective sits right on the diffraction-limited boundary for that thickness error.

### Pre-balanced aberrations: the Zernike polynomials

The calculation behind Figure 15(b) hides a condition. Rather than stopping at a coefficient, tilt and focus were adjusted to find the best position. As Section 4 showed, moving the focus halfway between the paraxial and marginal foci cuts the RMS of spherical aberration to a quarter (Figure 15(a)). Born and Wolf call this balancing of aberrations: mixing suitable amounts of lower-order aberration into a higher-order one maximizes the central intensity [2]. Strictly speaking, tilt and focus are not aberrations at all. They move the whole intensity distribution sideways or along the axis without changing its shape [2, 4].

F. Zernike's circle polynomials build this balancing in from the start [2, 4]. Zernike spherical aberration is $6\rho^4 - 6\rho^2 + 1$: the $\rho^4$ term premixed with the focus term ($\rho^2$) that best cancels it, and a piston term. Zernike coma, $3\rho^3 - 2\rho$, has tilt mixed in [4]. These functions are orthogonal over the disk. So when a wavefront is expanded in Zernike polynomials, the sum of the squared coefficients is the squared RMS [4], and each coefficient stands for the share of the image degradation due to that aberration, unmixed with the others. This is why interferometers and wavefront sensors report their measurements as Zernike coefficients, and why adaptive optics removes those coefficients in real time [4].

### How Hubble was fixed

Back to Hubble. A primary mirror too flat at the edge means a wavefront in the exit pupil carrying spherical aberration of the $\rho^4$ form. Wave aberrations simply add in the exit pupil, so placing in the beam a mirror that produces $\rho^4$ of the opposite sign cancels it. The new camera WFPC2 and the corrective package COSTAR for the other instruments, installed by astronauts in December 1993, contained exactly such optics [1, 6]. After the repair, more than 70% of the starlight fell in the central disk [1]. The mirror itself still has its original shape. The centre panel of Figure 1 is the result.

## Summary

- Aberrations are the terms the paraxial approximation drops. They are defined as the optical path $W$ by which the actual wavefront departs from the reference sphere at the exit pupil, and the distance by which a ray misses the image point is the slope of $W$.
- By rotational symmetry, $W$ is a function of $h^2$, $\rho^2$ and $h\rho\cos\theta$. Its five fourth-order combinations are spherical aberration, coma, astigmatism, field curvature and distortion. The first three blur the image; the last two change its shape and position.
- Spherical aberration is reduced by lens bending and removed with positive–negative pairs or aspheres. The best geometrical focus (three quarters of the way) and the best diffraction focus (halfway) differ.
- Coma is a magnification that varies over the rings of the aperture, and it vanishes when the sine condition holds. It can also be controlled by the stop position.
- Astigmatism produces two focal lines, and the tangential image surface is three times as far from the Petzval surface as the sagittal one. Petzval curvature is fixed by the powers and glasses alone and cannot be changed by the stop.
- Distortion depends only on the chief ray, so the stop position sets it, and symmetric arrangements cancel it. In metrology it is removed by calibration.
- Chromatic aberration is removed with two glasses of different Abbe numbers ($f_1V_1 + f_2V_2 = 0$). Since only two colours are matched, a secondary spectrum remains.
- Third-order aberrations can be computed surface by surface from just the paraxial marginal and chief rays, and summed (the Seidel sums).
- Small aberrations are judged by the RMS wave aberration. A Strehl ratio of 0.8, or an RMS of $\lambda/14$, marks the diffraction limit, and Zernike polynomials are aberrations pre-balanced for that judgement.

This concludes the geometrical optics foundations. The rays of Part 1 came from the limit of zero wavelength; at the end of this post, as aberrations shrank to the size of a wavelength, the ray picture gave way to diffraction. The next series, wave optics foundations, starts from there. Its first part is interference, what happens when two waves overlap. What the bottom row of Figure 3 and the Strehl ratios of Figure 15 actually compute will become clear step by step in that series.

## References

1. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Section 5.7 — the Cooke triplet, the Schmidt camera. Section 6.3 — third-order theory and the Seidel aberrations, spherical aberration (circle of least confusion, Rayleigh's quarter wave, the Hubble Space Telescope, lens orientation, aplanatic points and oil-immersion objectives), coma (curved principal surfaces, comatic circles, the sine condition, Wollaston's stop), astigmatism, Petzval field curvature and field flatteners, distortion and stop position, chromatic aberration and achromatic doublets. Chapter 1 — Newton and chromatic aberration.
2. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Sections 5.1–5.6 — wave and ray aberrations, the symmetry argument and primary aberrations, comatic circles, tangential and sagittal image surfaces, distortion, Seidel formulae in terms of two paraxial rays, Petzval's theorem, primary aberrations of a thin lens. Sections 9.1–9.3 — diffraction theory of aberrations, Zernike circle polynomials and balancing of aberrations, the Maréchal criterion and tolerances for primary aberrations.
3. J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed. (McGraw-Hill, 1996), Section 6.4 — the Gaussian reference sphere, the generalized pupil function, effects of aberrations on imaging.
4. R. K. Tyson and B. W. Frazier, *Principles of Adaptive Optics*, 5th ed. (CRC Press, 2022), Section 1.3 — the Strehl ratio and RMS wavefront error, representation of Seidel terms, the Zernike series and its properties.
5. D. B. Murphy and M. W. Davidson, *Fundamentals of Light Microscopy and Electronic Imaging*, 2nd ed. (Wiley-Blackwell, 2013), Chapter 4 — principal aberrations of lenses, objective classes (achromat, fluorite, apochromat, plan), cover slip thickness and correction collars.
6. NASA, [Hubble's Mirror Flaw](https://science.nasa.gov/mission/hubble/observatory/design/optics/hubbles-mirror-flaw/) — the cause of the primary mirror flaw and the 1993 repair.
7. Kepler & K2 Science Center, [Characteristics of the Kepler Space Telescope](https://keplergo.github.io/KeplerScienceWebsite/the-kepler-space-telescope.html) — the Schmidt telescope, curved focal plane and per-module field flatteners.

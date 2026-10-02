---
layout: post
date: 2026-11-17 20:00:00 +0900
title: "Geometrical Optics Foundations 3 — Stops and Pupils: Which Rays Get Through, and How Bright Is the Image?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/geometric-optics-stops-pupils-etendue/
page_id: geometric-optics-stops-pupils-etendue
categories: [Optics, Geometrical Optics]
tags: [geometric-optics, aperture-stop, pupil, etendue, telecentric]
description: "Starting from why stopping down never crops a photo, this post follows pupils, chief rays, f-numbers, vignetting, étendue and telecentricity."
math: true
---

[Part 2](/en/posts/geometric-optics-paraxial-abcd/) described how rays pass through an optical system using 2×2 matrices. But a matrix knows nothing about the size of a lens. Taken literally, it would let a ray one metre off axis pass straight through. Real lenses have rims, and inside a camera lens there is an adjustable diaphragm. This post deals with which rays those openings let through, and how bright the image turns out as a result.

Figure 1 shows the same camera lens seen from the front. On the left the diaphragm is wide open at f/1.4; on the right it is closed down to f/16. The bright disc inside looks like the diaphragm itself. It is not: it is the image of the diaphragm seen through the lens elements in front of it. What this disc is will be the subject of Section 2.

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/ext-minolta-aperture.jpg" alt="A lens viewed from the front with its diaphragm open and closed" width="640">
_Figure 1. A Minolta 50 mm f/1.4 lens viewed from the front, at f/1.4 (left) and f/16 (right). The bright disc inside is the image of the diaphragm formed by the lens elements in front of it — the entrance pupil. (Source: Leonrw, [16 minolta 50mm](https://commons.wikimedia.org/wiki/File:16_minolta_50mm.jpg), [CC BY 3.0](https://creativecommons.org/licenses/by/3.0))_

Anyone who has taken photographs knows one thing. Stopping down the lens makes the picture darker, but it never crops the edges. The edge of the film or sensor, on the other hand, cuts the image off. Both are openings that block light, so why do they behave so differently? The answer lies in where the opening sits.

## 1. Two Kinds of Stop

Openings in an optical system fall into two kinds according to what they do [1, 2]. One limits how wide the bundle of rays from an axial object point can spread; it is called the aperture stop. The adjustable blades inside a camera lens are an aperture stop. The other limits how much of the object ends up in the image; it is called the field stop. In a camera, the edge of the sensor plays this role. The aperture stop sets how much light reaches each image point, and the field stop sets how large the image is [1].

The phase-space picture of Part 2 makes the difference clear. Plot the rays arriving at the image plane by their height $y$ and arrival angle $u$. For a camera photographing a distant scene, the rays converging on an image point $y$ have crossed the lens at different heights, so their arrival angles are spread over a certain range. That range is set by the size of the opening at the lens. All the rays at the image plane therefore form the grey band in Figure 2.

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/en/fig2-phase-space-stops.png" alt="How the aperture stop and the field stop cut a ray bundle in phase space" width="760">
_Figure 2. Height (horizontal) and angle (vertical) of the rays arriving at the image plane of a 50 mm lens imaging a distant scene. (a) Closing the aperture stop thins the band vertically: every height survives, with fewer rays arriving at each. (b) Narrowing the field stop cuts off the ends of the band: the remaining heights still receive all their rays, while the outer heights disappear entirely._

Closing the aperture stop thins this band vertically (Figure 2(a)). Every image height survives; only the number of rays arriving at each height falls. The whole picture therefore darkens evenly. Narrowing the field stop cuts the band horizontally (Figure 2(b)). The remaining heights still receive all their rays, but the outer heights vanish altogether. The image is cropped.

Why the same kind of opening cuts angles in one place and heights in another was shown in Section 6 of Part 2: the matrix between focal planes swaps position and angle. An opening in a plane conjugate to the image — a plane where an image is formed — cuts heights. An opening in a plane where position and angle have been swapped relative to the image cuts angles. An optical system contains these two kinds of plane, alternating along its length.

The microscope shows this structure most clearly. Murphy and Davidson arrange the planes of a microscope into two sets [4]. One set, the field planes, are where images form: the field diaphragm of the illuminator, the specimen, the intermediate image inside the eyepiece, and the retina. The other set, the aperture planes, are the lamp filament, the condenser's front aperture, the back focal plane of the objective, and the pupil of the eye. The two sets interleave (Figure 3).

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/ext-kohler-illumination.png" alt="Illumination and imaging paths of a Köhler-illuminated microscope, with aperture planes 1–4 and field planes A–D" width="760">
_Figure 3. A transmitted-light microscope with Köhler illumination. Top: the illumination path. Bottom: the imaging path. Planes 1–4 are aperture planes (lamp filament, condenser diaphragm, objective back focal plane, eye pupil); A–D are field planes (field diaphragm, specimen, intermediate image, retina). The points where rays come to a focus mark the conjugate planes of each set. (Source: Mrmw, Dietzel65, Zephyris (original), [Kohler Illumination en](https://commons.wikimedia.org/wiki/File:Kohler_Illumination_en.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

The distinction is practical. If a speck of dirt appears sharp in the microscope image, it lies in one of the four field planes. If it rotates when you turn the eyepiece, it is near the intermediate image in the eyepiece; if it moves with the slide, it is on the specimen [4]. Dirt in an aperture plane is not imaged sharply; it only changes the brightness a little. Köhler illumination puts this principle to work. When the filament is imaged onto the condenser diaphragm, an aperture plane, every point of the filament illuminates the whole specimen evenly, so the specimen is uniformly lit even though the filament itself is irregular [4].

## 2. Entrance and Exit Pupils

Back to the bright disc of Figure 1. The diaphragm is buried in the middle of the lens. What our eye sees when it looks in from the front is the image of the diaphragm formed by the lenses in front of it. This is the entrance pupil [1, 2]. Likewise, the image of the diaphragm seen from behind is the exit pupil. A bundle of rays from an object point aims at the entrance pupil on the way in, and converges on the image point as if it came from the exit pupil [1]. If the aperture stop sets the width of the bundle inside the system, the pupils translate that limit into the language of object space and image space.

Pupils are practical quantities. When you look through binoculars, your eye must sit at the exit pupil. The human pupil ranges from 2 to 8 mm in diameter with lighting conditions, so binoculars for night use need an exit pupil of at least 8 mm, while 3–4 mm is enough for daytime [1].

When a system has several openings, which one is the aperture stop? It need not be the smallest. The method given by Born and Wolf and by Hecht is this. Image each opening into object space through the lenses in front of it, then look at those images from the axial object point. The image that subtends the smallest angle is the entrance pupil, and the opening it came from is the aperture stop [1, 2].

Hecht's Example 5.6 is a good exercise [1]. A lens 140 mm in diameter with a focal length of 100 mm has a 40 mm hole 80 mm behind it, and the object point is 200 mm in front of the lens. The lens itself has nothing in front of it, so it is seen as it is. The hole is seen through the lens, and since it lies inside the focal length, its image is magnified and virtual. Using the matrices of Part 2, that image is 400 mm to the right of the lens and 200 mm in diameter. From the object point, the lens rim subtends a half-angle of 19.3° and the image of the hole only 9.5° (Figure 4). Because the image of the hole subtends the smaller angle, this hole — smaller than the lens and hidden behind it — is the aperture stop.

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/en/fig4-hecht-entrance-pupil.png" alt="Hecht's Example 5.6: the image of the hole becomes the entrance pupil" width="760">
_Figure 4. Hecht's Example 5.6. Seen through the lens, the 40 mm hole behind it forms a 200 mm virtual image 400 mm to the right. From the object point this image subtends a half-angle of 9.5°, smaller than the lens rim (19.3°), so the hole is the aperture stop and its image is the entrance pupil. The blue lines are the outermost rays that actually get through; the image point forms 200 mm behind the lens._

One point is worth noting. Which opening acts as the aperture stop can depend on where the object is [1]. In this example, moving the object point further away shrinks the angle subtended by the lens rim quickly but that of the hole's image only slowly. Once the object is more than 933 mm in front of the lens, the lens rim subtends the smaller angle, and the lens itself becomes the aperture stop.

The exit pupil is also the door to wave optics. As Goodman stresses, the wavefront converging on the image is truncated at the exit pupil on its way there [3]. Without aberrations, the wavefront leaving the exit pupil is a perfect sphere centred on the image point, and how much the image spreads is determined by how that sphere is cut off at the exit pupil [3]. When the wave optics series computes the point-spread function, the exit pupil is where it starts.

## 3. Chief and Marginal Rays

Designing an optical system does not require tracing every ray. Two rays determine almost all the paraxial properties of a system [1, 2]. One is the marginal ray, which leaves the axial object point and passes the edge of the entrance pupil. Where it crosses the axis again is the image position, and its angle to the axis fixes the width of the bundle — the aperture. The other is the chief ray, which leaves the edge of the object and passes through the centre of the entrance pupil [1, 2]. It passes through the centre of the aperture stop, travels as if it emerged from the centre of the exit pupil, and the height at which it meets the image plane is the image size. Think of it as the ray that represents the middle of the bundle from an off-axis object point.

Figure 5 draws these two rays for a system of two convex lenses with a stop between them. The object is 150 mm in front of the first lens and 10 mm tall. The chief ray crosses the axis exactly at the centre of the stop. The entrance pupil is the image of the stop; since the stop lies inside the focal length of the first lens, the image is magnified and virtual, and it falls behind the stop (55 mm from the front of the first lens, radius 9.4 mm). Extending the two object-space rays straight ahead, they aim at the centre and the edge of exactly this entrance pupil. The exit pupil, conversely, lies in front of the stop (at 24 mm, radius 7.2 mm), and extending the two image-space rays backwards, they appear to come from it. The image forms 88 mm from the front of the first lens — 27 mm behind the second lens — at a height of −4.0 mm.

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/en/fig5-marginal-chief.png" alt="Marginal ray, chief ray, entrance pupil and exit pupil in a two-lens system with a stop" width="780">
_Figure 5. Two N-BK7 convex lenses with a stop (radius 5 mm) between them. Blue: the marginal ray from the axial object point to the edge of the entrance pupil. Red: the chief ray from the top of the object to the centre of the entrance pupil. Orange marks the entrance pupil and green the exit pupil; both are virtual. Dotted lines extend the object-space rays straight to the entrance pupil and the image-space rays straight back to the exit pupil._

Used together, the two rays give the invariant of Part 2, Section 4 a new face. With the chief ray's height and angle written $(\bar{y}, \bar{u})$ and the marginal ray's $(y, u)$, the quantity $H = n(\bar{y}u - y\bar{u})$ is the same at every surface [2]. In this system $H = 0.458$ mm·rad everywhere. At the object plane the marginal ray has zero height, so $H$ is the object height times the aperture angle; at the image plane it is the image height times the image-side aperture angle. In other words, the product "size of the field × width of the bundle" does not change through the system. Enlarge the image and the bundle narrows; widen the bundle and the image shrinks. Section 6 will show that this product limits how much light a system can carry.

## 4. f-Number and Numerical Aperture

Two numbers are used to express the width of the bundle. For cameras and telescopes looking at distant objects, the f-number is the focal length divided by the entrance pupil diameter [1, 2]:

$$ F = \frac{f}{D} $$

A lens with a 50 mm focal length and a 25 mm entrance pupil has an f-number of 2, written f/2 [1]. For microscopes looking at nearby specimens, the numerical aperture is defined from the object-side half-angle $\theta$ and the refractive index of the medium on that side [2, 4]:

$$ \mathrm{NA} = n\sin\theta $$

Why a smaller f-number gives a brighter image is simple. The light entering the lens from one part of a distant object is proportional to the area of the entrance pupil, that is, to $D^2$. The area of the image over which that light spreads is proportional to the square of the image size, that is, to $f^2$. The irradiance of the image is therefore proportional to $(D/f)^2 = 1/F^2$ [1]. The f-stop markings 1, 1.4, 2, 2.8, 4, 5.6 grow by a factor of $\sqrt{2}$ so that each step halves the light exactly. That is why f/1.4 at 1/500 s, f/2 at 1/250 s and f/2.8 at 1/125 s all deliver the same amount of light [1].

The two numbers are related. For a lens imaging a distant object with an image-side aperture angle $\theta'$, paraxially $\sin\theta' \approx D/2f$, so $F \approx 1/(2\,\mathrm{NA}')$. For a lens obeying the sine condition of Part 2, the relation becomes exact.

The role of the refractive index $n$ in the numerical aperture shows up in oil-immersion objectives. Filling the space between specimen and lens with oil of index 1.515 instead of air lets the lens collect steep rays that would otherwise be totally reflected at the coverslip and never reach it [4]. Abbe's NA 1.3 oil-immersion objective, advertised in a Zeiss catalogue of 1888, is believed to have resolved down to the theoretical limit of 0.26 µm in green light [4]. Microscope image brightness scales as $(\mathrm{NA}/M)^2$ in transmitted light, and as $\mathrm{NA}^4/M^2$ when the objective both illuminates and collects, as in fluorescence [4]. A wider aperture gathers more light, and a higher magnification spreads that light over a larger image. For the same reason a 40×/1.3 objective gives a brighter image than a 60×/1.4 [4].

## 5. Vignetting and the cos⁴ Law

So far we have reasoned from an axial object point. What about the oblique bundle from the edge of the object? On axis only the aperture stop cuts the bundle, but an oblique bundle can also catch on other openings such as lens rims [1, 2]. Some rays that could have passed the aperture stop are then blocked elsewhere, and the image gradually darkens towards its edges. This is vignetting.

For the test system of Figure 5, the fraction of rays from each object height that get through the entrance pupil is the blue line in Figure 6(b). All rays pass up to an object height of 30 mm. At 35 mm the rim of the first lens's front surface starts to clip the bundle and 94.4% get through; this drops to 84% at 40 mm and 60% at 50 mm, and the light is cut off completely near 62 mm. At 35 mm only that one rim clips, so the same fraction can be computed exactly as the overlap area of two circles: the entrance pupil, and the rim projected from the object point onto the entrance-pupil plane. It agrees with the grid count to within $10^{-5}$. Born and Wolf note that designers sometimes use vignetting deliberately, to cut off strongly aberrated rays near the edge of the field [2].

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/en/fig6-cos4-vignetting.png" alt="Irradiance from a disc source compared with the cos⁴ law, and vignetting and relative illumination of the test system" width="780">
_Figure 6. (a) Irradiance produced by a uniformly bright disc (the pupil), integrated directly as a function of the chief-ray angle. A small pupil follows the cos⁴ law; a large pupil falls off more slowly. (b) For the test system of Figure 5, the fraction of rays passing the pupil against object height (blue) and that fraction multiplied by cos⁴ — the relative illumination (orange dashed)._

Even with no vignetting at all, the edges of the image grow darker. Born and Wolf show that the image irradiance falls as the fourth power of the cosine of the chief-ray angle $\varphi$ [2]:

$$ E(\varphi) = E(0)\cos^4\varphi $$

Why the fourth power? Picture looking at the exit pupil from a point in the image plane, and four cosines appear one by one. Seen obliquely, the exit pupil looks like a flattened ellipse: one. The light strikes the image plane at a slant and spreads over a larger area: another. The distance to the exit pupil grows by $1/\cos\varphi$, and the inverse-square dependence on distance contributes the last two. At 30° the irradiance drops to 56%, at 45° to 25%. Figure 6(a) treats the exit pupil as a uniformly bright disc and integrates the irradiance at the image plane directly. When the pupil radius is 1% of the distance, the result matches cos⁴ to within $3\times10^{-5}$. When the pupil radius is half the distance, however, the irradiance falls more slowly, reaching 1.28 times cos⁴ at 40°. The cos⁴ law is an approximation that assumes a small exit pupil [2].

## 6. Étendue and the Brightness Theorem

Part 2, Section 4 showed that the area a bundle of rays occupies in phase space is unchanged by any optical system. Real rays have two height directions and two angle directions, so in three dimensions this area becomes "area of the surface the bundle crosses × solid angle it spreads over × $n^2$". This quantity is called étendue, and in a lossless system it is conserved. It is the Lagrange invariant of Section 3 taken in both directions — a quantity proportional to its square.

A remarkable conclusion follows. Born and Wolf show that, without losses, the brightness (radiance) of an image cannot exceed that of the object [2]. When object and image are in the same medium, equality holds only if there are no losses at all. A lens seems to concentrate light and make the image brighter, but all it really does is gather light of the same brightness from a wider range of angles. The irradiance at a point of the image is $E = \pi B\sin^2\theta'$: the radiance $B$ of the object times a factor set by the image-side aperture [2].

Think of burning a piece of paper with a magnifying glass. For an f/1 magnifier, with a diameter equal to its focal length, the irradiance of the Sun's image at the focus is $\pi B\sin^2\theta' \approx \pi B/4$ — a quarter of the light leaving the Sun's surface. Even with a bigger lens and a shorter focal length that push $\sin\theta'$ towards 1, the irradiance cannot exceed $\pi B$, the emittance of the Sun's surface itself. No arrangement of lenses and mirrors can concentrate sunlight into a spot hotter than the surface of the Sun. If it could, heat would flow on its own from a colder body to a hotter one, violating the second law of thermodynamics. The area conservation of geometrical optics meshes exactly with thermodynamics.

The microscope brightness law $(\mathrm{NA}/M)^2$ of Section 4 tells the same story. The étendue an objective can collect is proportional to the specimen area times NA², and that light is spread over an image magnified by $M$, an area $M^2$ times larger.

## 7. Telecentricity

Where the aperture stop sits affects not only brightness but also measurement accuracy. Born and Wolf single out the case where the stop is placed in the back focal plane of the lenses in front of it [2]. The entrance pupil then lies at infinity, and every chief ray in object space runs parallel to the axis. Such a system is called telecentric on the object side. Placing the stop in the front focal plane of the lenses behind it makes the image-side chief rays parallel: telecentric on the image side [2].

Why does this matter? When the object moves slightly out of focus, its image blurs. The centre of a blurred spot is where the chief ray lands, so measurements of size are taken from these centres. In an ordinary lens the chief ray enters at a slant, heading for the centre of the lens, so as the object moves back and forth, the height at which the chief ray meets the image plane changes, and so does the measured size. In a telecentric lens the chief ray is parallel to the axis, so it enters at the same height whether the object moves forwards or backwards, and the measured size does not change. This is why Born and Wolf describe the arrangement as useful for measuring the size of objects [2].

Figure 7 computes this for a 50 mm focal-length lens imaging an object at unit magnification. With the stop at the lens, the measured size shifts by about 1% for every millimetre of object displacement. Move the stop to the back focal plane, and the measured size stays put across ±5 mm. This is why optical metrology — measuring the line widths of semiconductor patterns or the dimensions of machined parts — relies on telecentric lenses.

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/en/fig7-telecentric.png" alt="Error in measured size against focus error for two stop positions" width="640">
_Figure 7. A 50 mm thin lens images an object at unit magnification onto a fixed sensor. When the object moves by δ, how much does the measured size (the centre of the blurred spot) change? With the stop at the lens (red) it changes by about 1% per millimetre; with the stop at the back focal plane (blue) it does not change._

There is a price. For the chief rays to be parallel to the axis, the front lens must be as large as the field to be measured, so that it can accept rays from the edge of the object. The bi-telecentric lens of Figure 8 has a front element 208 mm in diameter, far larger than its camera mount at the back.

<img src="/assets/img/posts/geometric-optics-stops-pupils-etendue/ext-bi-telecentric-lens.jpg" alt="A bi-telecentric lens with a very large front element" width="560">
_Figure 8. A bi-telecentric lens. Its front element is 208 mm in diameter, as large as the field it is meant to measure. The small cylinder at the back is the camera mount. (Source: Laserlicht, [Bi telecentric lens](https://commons.wikimedia.org/wiki/File:Bi_telecentric_lens.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

## Summary

- Openings in an optical system come in two kinds. The aperture stop cuts angles and darkens the whole image evenly; the field stop cuts heights and sets the extent of the image. Which role an opening plays depends on whether it sits in a plane conjugate to the image or in a plane where position and angle are swapped.
- The entrance and exit pupils are the images of the aperture stop seen from the front and from behind. When there are several openings, the image that subtends the smallest angle from the object point is the entrance pupil, and the opening it came from is the aperture stop.
- Two rays, the marginal ray and the chief ray, fix the image position, the image size and the pupils, and the Lagrange invariant built from them (field × aperture) is unchanged through the system.
- Image irradiance scales as $1/F^2$ and falls towards the edges because of vignetting and the cos⁴ law. The cos⁴ law is an approximation for a small exit pupil.
- Because étendue is conserved, an image cannot be brighter than its object. That is why a lens cannot make a spot hotter than the surface of the Sun, and it is where geometrical optics meets thermodynamics.
- Placing the stop at a focal plane makes a system telecentric, so the measured size does not change when the focus is off.

Up to now we have assumed that rays stay in the paraxial regime. In reality the marginal ray crosses the axis ahead of the paraxial prediction, and the chief ray lands at a different height. The next post deals with these departures — aberrations. The marginal and chief rays defined here become the two reference rays for computing them, and aberrations will be redefined as the difference between the ideal sphere and the actual wavefront at the exit pupil.

## References

1. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Section 5.3 — aperture and field stops, entrance and exit pupils, chief and marginal rays, vignetting, f-number and diaphragm scales, Example 5.6.
2. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Sections 4.8.2–4.8.3 — stops and pupils (the theory of Abbe and von Rohr), entrance windows and field angle, telecentricity, numerical aperture and f-number, vignetting, brightness and illumination of images, the cos⁴ law; Section 4.4.5 — the Smith–Helmholtz formula.
3. J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed. (McGraw-Hill, 1996), Appendix B.5 — entrance and exit pupils, and diffraction calculations starting from the exit pupil.
4. D. B. Murphy and M. W. Davidson, *Fundamentals of Light Microscopy and Electronic Imaging*, 2nd ed. (Wiley-Blackwell, 2013), Chapter 1 — field and aperture planes, Köhler illumination; Chapter 4 — image brightness of objectives; Chapter 6 — numerical aperture and oil-immersion objectives.

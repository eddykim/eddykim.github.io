---
layout: post
date: 2026-11-11 20:00:00 +0900
title: "Geometrical Optics Foundations 1 — Where Do Rays Come From? Fermat's Principle and the Eikonal"
lang: en
lang-exclusive: ["en"]
permalink: /posts/geometric-optics-fermat-eikonal/
page_id: geometric-optics-fermat-eikonal
categories: [Optics, Geometrical Optics]
tags: [geometric-optics, fermat-principle, eikonal, mirage, caustic]
description: "Send the wavelength to zero and rays fall out of the wave equation. Mirages and caustics show what that limit means and where it fails."
math: true
---

Look far down a straight asphalt road on a summer afternoon and you will see a pool of water lying across it. Walk towards it and it retreats at the same pace; you never reach it (Figure 1). What the pool reflects is the sky.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/ext-hot-road-mirage.jpg" alt="An inferior mirage that looks like a puddle on a hot road" width="600">
_Figure 1. An inferior mirage over a hot road. The distant cars appear reflected upside down in the road surface, as if the road were wet. (Source: Brocken Inaglory, [Hot road mirage](https://commons.wikimedia.org/wiki/File:Hot_road_mirage.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/))_

Textbooks explain this with geometrical optics. The air near the ground is hot, its refractive index is lower, and the index gradient bends rays upwards. In 1959, however, C. V. Raman, who won the Nobel Prize for the scattering effect that bears his name, called this explanation "a kind of make-believe" [5]. His reasoning went as follows. In a layered medium, Snell's law tilts a ray closer and closer to the horizontal until it runs level, and from then on it simply travels parallel to the layers. With no interface, there can be no reflection. Raman therefore reworked the mirage with wave optics.

Raman's objection was wrong. Saying exactly why, however, requires first asking what a ray is. As [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/) showed, light is an electromagnetic wave governed by Maxwell's equations, and a wave has no lines in it. The [ray-tracing series (Geometrical Optics 1)](/en/posts/raytracing-spherical-lens-refraction/) drew light as straight lines that bend only at interfaces, but it never asked why that picture holds.

This series builds the theoretical frame of geometrical optics in four parts. Part 1 deals with the concept of a ray itself; Part 2 with paraxial optics and ray transfer matrices; Part 3 with stops and fields; Part 4 with aberrations. This post derives rays from the limit of vanishing wavelength and uses them to compute the mirage, which answers Raman. It then shows why Fermat's principle is about a *stationary* path rather than a *least* one, and ends at the caustic, where the approximation breaks down.

## 1. The Ray Hypothesis

In his *Opticks*, Newton asked whether rays of light are streams of very small particles emitted by shining bodies. His contemporary Huygens saw light as a wave spreading out in all directions [4]. We now know that light is a wave. Why, then, may we draw it as lines?

Compare the sizes involved. Visible light has a wavelength of about 0.5 µm, while a lens is tens of millimetres across. The ratio is below $10^{-4}$. Pass light through a small hole in an opaque screen and you get a thin pencil of light. Its edge looks sharp at first glance, but on closer inspection it carries bright and dark diffraction fringes over a width of about one wavelength [1]. When the hole is much larger than the wavelength, the fringes can be ignored. Shrinking the hole makes the pencil thinner, but once the hole approaches the wavelength in size, the light spreads out instead.

An infinitely thin pencil, that is, a ray, therefore exists only in the limit of zero wavelength. The optics of that limit is geometrical optics. The task is clear: take the wave equation, send the wavelength to zero, and see what survives.

## 2. The Eikonal Equation

A plane wave in a homogeneous medium is written $U = A\,e^{ik_0 n\,\mathbf{s}\cdot\mathbf{r}}$, where $k_0 = 2\pi/\lambda_0$ is the vacuum wavenumber and $\mathbf{s}$ the direction of travel. In a medium whose index varies from place to place, this form no longer holds as it stands. Over a region only a few wavelengths across, though, the wave should still look like a plane wave. Written as an equation, that expectation reads

$$
U(\mathbf{r}) = A(\mathbf{r})\,e^{ik_0 S(\mathbf{r})}
$$

Here $A$ is a slowly varying amplitude and $S$ is the phase expressed as a length. For a plane wave, $S = n\,\mathbf{s}\cdot\mathbf{r}$. Substituting this form into the scalar Helmholtz equation $\nabla^2 U + k_0^2 n^2 U = 0$ and dividing by $k_0^2$ gives

$$
\left(n^2 - \lvert\nabla S\rvert^2\right)A
+ \frac{i}{k_0}\left(2\nabla A\cdot\nabla S + A\,\nabla^2 S\right)
+ \frac{1}{k_0^2}\nabla^2 A = 0
$$

As the wavelength goes to zero, $k_0$ goes to infinity. The second and third terms vanish and only the first survives:

$$
\lvert\nabla S\rvert^2 = n^2(\mathbf{r})
$$

This is the eikonal equation. The name comes from the Greek word for "image" and was coined by H. Bruns in 1895 [1].

We obtained it from a scalar wave that ignores polarization, but starting from the vector electromagnetic field leads to the same equation. Different textbooks take different routes there. Born and Wolf substitute the same form into Maxwell's equations and obtain the eikonal equation as the condition for the resulting simultaneous equations in the electric and magnetic fields to have a non-trivial solution [1]. Orfanidis, from the same substitution, requires the electric field, the magnetic field and the direction of travel to form a mutually orthogonal right-handed triad; that requirement can only be met if the direction vector $\nabla S/n$ has unit length, which is the eikonal equation again [3]. All three routes reach the same conclusion. At short wavelengths the electromagnetic field behaves locally like a plane wave, and $\nabla S$ sets that plane wave's direction and wavenumber.

What does the equation say? Surfaces of constant $S$ are wavefronts. The distance between two neighbouring wavefronts $S$ and $S + dS$ is $dS/n$. Where the index is high, light travels slowly, so it takes less distance to accumulate the same phase and the wavefronts crowd together. Figure 2 shows computed wavefronts from a point source in a medium whose index grows upwards ($n^2 = 1 + 1.5y$). The red wavefronts get closer together towards the top.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/en/fig2-wavefronts-rays.png" alt="Wavefronts and rays from a point source in a graded-index medium" width="700">
_Figure 2. A point source in a medium whose index increases upwards. Red curves are wavefronts, joining points of equal optical path length; blue curves are rays. The two cross at right angles everywhere, and the rays bend upwards towards the higher index._

The two discarded terms come at a price. Dropping them requires the amplitude and the index to change very little over the distance of one wavelength. Section 7 returns to the places where this condition fails.

## 3. Rays and Optical Path Length

The eikonal equation gives wavefronts, not lines. Where are the rays? A ray is defined as a curve perpendicular to the wavefronts. With arc length $s$,

$$
n\frac{d\mathbf{r}}{ds} = \nabla S
$$

The left-hand side, $n\,d\mathbf{r}/ds = n\mathbf{s}$, is called the ray vector. The definition is not arbitrary. In this approximation the time-averaged electric and magnetic energies are equal, the time-averaged Poynting vector points along the ray, and energy is carried in that direction at speed $v = c/n$ [1, 3]. A ray is the path along which energy flows. The identification holds only in isotropic media, however. Inside a birefringent crystal, the wavefront normal and the energy flow part ways.

Integrating along a ray gives the optical path length (OPL):

$$
[P_1P_2] = \int_{P_1}^{P_2} n\,ds = S(P_2) - S(P_1) = c\int_{P_1}^{P_2} dt
$$

Since $n\,ds = (c/v)\,ds = c\,dt$, the optical path length is the travel time multiplied by the vacuum speed of light. It is the distance light would have covered in vacuum in the same time, and also the accumulated phase counted in wavelengths [2]. The wavefronts in Figure 2 were drawn by accumulating optical path length along 26 rays from the source and joining the points with equal values. That wavefronts and rays meet at right angles then follows on its own in the resulting curves.

What about intensity? Setting the $1/k_0$ term from before to zero gives the transport equation $2\nabla A\cdot\nabla S + A\nabla^2 S = 0$. Rearranged, it says that along a thin tube formed by a bundle of rays, the product of intensity and cross-sectional area stays constant [1]. Light brightens where rays converge and dims where they spread. For a point source, the cross-section grows with the square of the distance, which gives the inverse-square law. This law hides a problem, though. If rays converge to a point or a line, the cross-section goes to zero and the intensity becomes infinite. Section 7 deals with this.

## 4. The Ray Equation and the Mirage

Can rays be drawn from the index distribution alone, without first finding the wavefronts? Differentiating the ray definition once more with respect to $s$ and using the eikonal equation eliminates $S$ [1, 3]:

$$
\frac{d}{ds}\left(n\frac{d\mathbf{r}}{ds}\right) = \nabla n
$$

This is the ray equation. In a uniform medium the right-hand side is zero, so rays are straight lines. Where the index varies, the curvature $1/\rho$ of a ray equals $\boldsymbol{\nu}\cdot\nabla\ln n$ [1], where $\boldsymbol{\nu}$ is the unit vector pointing to the side the ray curves towards. Curvature cannot be negative, so a ray **always bends towards the higher index**. Wavefronts make this intuitive. The end of a wavefront in the higher index moves slowly and the other end moves quickly, so the whole wavefront swings towards the higher index.

Raman's objection falls apart on this curvature formula. Even at the moment a ray becomes level, the index gradient still points vertically, and by the curvature formula the ray keeps bending upwards. There is no reason for it to carry on horizontally. S. Vince and W. H. Wollaston had already pointed this out in 1799 and 1800, and James Thomson, Lord Kelvin's elder brother, took up the same question in 1873 [5]. Raman's wave-optical calculation was itself correct. What was wrong was the premise that geometrical optics falls short.

In a stratified medium, where the index depends only on the height $y$, the horizontal component of the ray vector, $n\cos\theta$, is conserved ($\theta$ measured from the horizontal). Rewritten with the angle from the layer normal, this says that $n\sin\theta_\perp$ is constant. It is Snell's law for an infinite stack of infinitely thin layers.

The mirage can now be computed. In the visible, air has an index of $n - 1 \approx 2.8\times10^{-4}$ (15 °C, 1 atm), very close to 1. The Gladstone–Dale relation makes $n - 1$ proportional to the gas density [2], and at constant pressure the density is inversely proportional to absolute temperature. Hotter air therefore has a lower index. Take the air just above the ground as 50 °C and the air above as 30 °C, joined by an exponential boundary layer 10 cm thick. The index difference between the ground and eye level (1.5 m) is then $\Delta n = 1.63\times10^{-5}$.

How can such a small difference turn a ray around? Trace backwards a ray leaving the eye at an angle $\theta$ below the horizontal. The conserved quantity gives $n(y)\cos\theta(y) = n_{\text{eye}}\cos\theta$. The ray levels off and starts back up at the height where $n(y) = n_{\text{eye}}\cos\theta$, and such a height exists only if the index at the ground is lower still. The condition for the ray to return is therefore

$$
\theta < \theta_c = \arccos\frac{n_{\text{ground}}}{n_{\text{eye}}} \approx \sqrt{2\Delta n} = 5.71\ \text{mrad} = 0.327^\circ
$$

Figure 3(b) shows the ray equation integrated numerically. Rays looking down at less than 0.327° curve upwards before reaching the ground and head off into the sky. An eye looking in those directions receives sky instead of road. Because the eye assumes light travels straight, the sky appears to lie below the road surface, as if reflected. That is the pool. Rays looking down more steeply hit the road, so the road is what you see. In this model the boundary between the two, where the pool begins, lies about 290 m ahead. What defines the pool is the viewing angle, not the distance, so as you approach, the pool retreats with you.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/en/fig3-road-mirage-rays.png" alt="Temperature and index profiles above a road, and rays traced back from the eye" width="760">
_Figure 3. (a) Model temperature (solid) and refractive index (dashed). (b) Rays traced back from an eye 1.5 m high at 0.05°–0.65° below the horizontal, with the vertical scale greatly exaggerated. Rays shallower than the 0.327° critical angle (blue) turn back up and see the sky; steeper rays (grey) hit the road. The orange dotted line is the apparent direction the eye assumes._

Because the critical angle scales as $\sqrt{2\Delta n}$, an index difference of order $10^{-5}$ gives a critical angle of only about 0.3°. In A. T. Young's words, mirages are confined to a narrow strip, a fraction of a finger's width at arm's length [7]. Figure 3(b) shows one more thing worth noticing: the blue rays that turn back up cross one another. Light from a single distant point can therefore reach the eye by two paths, a straight one and one that has curved back up, which places an inverted image beneath the upright one. Such crossings require the ray curvature to change with height. If the index gradient were the same at every height, all rays would bend alike and never cross, and there would be no inverted image [7]. A boundary layer with the gradient concentrated near the ground is what makes the mirage. Conversely, when the lower air is colder, as over a cold sea, rays bend downwards and objects beyond the horizon appear lifted: a superior mirage (Figure 4).

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/ext-farallon-mirages.jpg" alt="Inferior mirage, no mirage and superior mirage of the Farallon Islands" width="600">
_Figure 4. A composite of the same islands (the Farallon Islands) photographed on different days. The top-right panel is an inferior mirage: the island floats above the sea with an inverted image hanging beneath it. The second panel shows a day without a mirage. The two lower panels and the background photograph are superior mirages, in which the island is stretched upwards into towers and stacked in several layers. (Source: Brocken Inaglory, [Farallon Islands at inferior mirage no mirage and superior mirage](https://commons.wikimedia.org/wiki/File:Farallon_Islands_at_inferior_mirage_no_mirage_and_superior_mirage.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

On the scale of the whole atmosphere, the index varies with the distance $r$ from the centre of the Earth. In such a spherically symmetric medium, $n r\sin\phi$ is constant along a ray, where $\phi$ is the angle between the ray and the radial direction [1]. This is Bouguer's formula, the optical counterpart of conservation of angular momentum in mechanics. It is why the Sun remains visible for a while after it has dropped below the horizon [2].

## 5. Snell's Law Returns at an Interface

What happens where the index jumps discontinuously, as at the surface of a lens? The ray equation uses the gradient $\nabla n$, so it cannot be applied across a discontinuity as it stands. Instead, go back to the definition of a ray. The ray vector $n\mathbf{s}$ is the gradient of $S$, so its curl is zero:

$$
\nabla\times(n\mathbf{s}) = 0
$$

Take a long, thin loop straddling the interface and apply Stokes' theorem. As the width of the loop goes to zero, the tangential component of the ray vector must be continuous across the interface [1]. With angles of incidence and refraction $\theta_1$ and $\theta_2$, this condition is $n_1\sin\theta_1 = n_2\sin\theta_2$, and it also places the refracted ray in the plane of incidence. Applied to a ray returning into the same medium, it gives the law of reflection.

Ellipsometry Foundations 1 obtained the same law from phase matching of plane waves at a flat interface. That derivation holds at any wavelength but needs a flat boundary. This one works for curved surfaces, but only when the radius of curvature is much larger than the wavelength. Each derivation covers the other's gap. It is this result that let the ray tracer in Geometrical Optics 1 apply Snell's law at each point of a spherical lens, measured from the local surface normal.

The argument can also be run in reverse: a continuous medium can be approximated by a stack of thin discrete layers. When the index steps between layers are small, the reflection at each interface is negligible, but the ray that reaches the lowest layer is totally reflected at some interface. As the layers get thinner, this total reflection turns into the continuous bending of the smooth medium. M. V. Berry showed that the limit is subtler than it looks [5]. Near an almost-level ray there is a narrow interval in which the continuum approximation breaks down, and however thin the layers are made, the interval only narrows; it never disappears. The "level ray" that Raman misunderstood lives in exactly this interval.

## 6. Fermat's Principle

The ray equation and Snell's law follow from a single principle. Hero of Alexandria held that reflected light takes the shortest path; in 1657 Fermat proposed extending it to refraction [2], and in 1662 he derived the law of refraction from the statement that light takes the path of least time [8]. Since optical path length is $c$ times travel time, least time means least optical path length.

The ray equation can be derived back from this principle [3]. Parametrize a path as $\mathbf{r}(\tau)$; its optical path length is $\int n(\mathbf{r})\lvert\dot{\mathbf{r}}\rvert\,d\tau$. Treating the integrand $L = n\lvert\dot{\mathbf{r}}\rvert$ like a Lagrangian in mechanics, the Euler–Lagrange equation reads

$$
\frac{d}{d\tau}\left(n\frac{\dot{\mathbf{r}}}{\lvert\dot{\mathbf{r}}\rvert}\right) = \lvert\dot{\mathbf{r}}\rvert\,\nabla n
$$

Changing variables to $ds = \lvert\dot{\mathbf{r}}\rvert\,d\tau$ recovers the ray equation of Section 4 exactly. The route from the eikonal and the route from the variational principle meet at the same equation.

A flat interface makes a simple check. Let a path run from $S = (0, 1)$ in air to $P = (2, -1)$ in glass ($n = 1.5$), crossing the interface $y = 0$ at $x$. Plotted against $x$, the optical path length has a single minimum (Figure 5(a)), and at that point $n_1\sin\theta_1 = n_2\sin\theta_2$ holds to within $10^{-8}$.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/en/fig5-fermat-opl.png" alt="Minimum optical path length for refraction, and reflected optical path lengths for three mirrors" width="780">
_Figure 5. (a) Optical path length of a path from air into glass. The minimum sits at the point given by Snell's law. (b) Optical path length of reflected paths joining the two foci S and P of an ellipse, plotted against the position of the reflection point. (c) The three mirrors share the same tangent at the vertex Q. For the plane mirror (blue) Q is a minimum; for the elliptical mirror (black) every point gives the same length; for the mirror more curved than the ellipse (red) Q is a maximum._

"Least", however, is not quite right. Figures 5(b) and (c) show reflected paths between the two foci $S$ and $P$ of an ellipse [2]. On an elliptical mirror every reflection point gives the same optical path length, equal to the major axis. A plane mirror tangent to the ellipse at Q makes the path through Q a minimum; a mirror more curved than the ellipse makes it a maximum. In all three cases the light reflects at Q, and the only thing they share is that the first-order change in optical path length vanishes. That is why the modern statement of Fermat's principle says the optical path length is **stationary**:

$$
\delta\int_{P_1}^{P_2} n\,ds = 0
$$

Here the textbooks disagree. Hecht calls the third case a maximum [2], whereas Born and Wolf, citing Carathéodory, note that the stationary value is never a true maximum [1]. They are comparing against different sets of paths. Among paths that only move the single reflection point along the mirror, Q is a maximum. Once arbitrary wiggly curves are allowed, a path can always be lengthened by jiggling it slightly, so Q becomes a saddle point.

When, then, is a ray a genuine minimum? Born and Wolf prove that it is one within a region where the rays from a point do not cross each other, and state that the minimum is lost once the ray passes a conjugate point, where neighbouring rays meet again [1]. We checked this numerically. A parabolic medium with $n^2 = n_0^2(1 - \alpha^2 y^2)$ ($n_0 = 1.5$, $\alpha = 0.5$) models a graded-index (GRIN) lens or fibre. We found the rays that start at the origin and pass through a given end point, divided each path into an 80-segment polyline, and computed the Hessian (second-derivative matrix) of its optical path length. Separately, we counted conjugate points as the number of times a neighbouring ray, launched at a slightly different angle, crosses the original ray.

| End point | Launch angle | Conjugate points passed | Negative Hessian eigenvalues |
|---|---|---|---|
| (4, 0.5) | 16.72° | 0 | 0 |
| (4, 0.5) | 43.83° | 1 | 1 |
| (9, 0.5) | −14.51° | 1 | 1 |
| (9, 0.5) | −40.12° | 2 | 2 |

Each conjugate point passed adds one direction along which the optical path length decreases. Only the ray in the first row is a genuine minimum. Minimizing the optical path length directly, starting from the straight line between the end points, converges onto that ray (height difference $1.6\times10^{-5}$). The other three are stationary but not minima. They are nonetheless all paths that light really takes.

Why stationary, of all things? Waves give the answer [2]. Place $S$ and $P$ 10 mm above a plane mirror and 20 mm apart. The wave from $S$ reaches $P$ by way of every point on the mirror. The wave arriving via point $x$ has phase $k\cdot\text{OPL}(x)$, and the field at $P$ is the sum of small arrows, or phasors, with these phases. Figure 6(a) shows that near the stationary point $x = 0$ the optical path length barely changes, so the phase turns slowly. Farther away, the phase spins rapidly.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/en/fig6-phasor-sum.png" alt="Phase of the wave via each mirror point, and the spiral traced by the running phasor sum" width="760">
_Figure 6. Reflection from a plane mirror, wavelength 0.5 µm. (a) Real part of the phase against the reflection point. It varies slowly only near the stationary point (orange band). (b) The running sum of the phasors from the left end (blue) to the right end (red). The region near the stationary point (thick orange line) produces most of the sum, while both ends curl into small spirals that cancel out._

Figure 6(b) traces the phasors added one after another. Far from the stationary point, the phasors point every which way and circle in place. Only within ±0.08 mm of the stationary point do they line up and build the sum. The magnitude of the total agrees with the stationary-phase estimate $\sqrt{2\pi/(k\,\text{OPL}'')}$ to within 4%. A ray is the path that interferes constructively with its neighbours, and that condition is precisely that the first-order change in optical path length vanishes. Whether it is a minimum or a maximum does not matter. This picture leads on to the Huygens–Fresnel principle, which the wave optics series takes up properly in its post on diffraction.

## 7. Where Geometrical Optics Breaks Down

The terms dropped in Section 2 could be neglected only while the amplitude and the index change very little over one wavelength. Two typical places violate this [1]. The first is the edge of a shadow, where light meets darkness and the amplitude changes abruptly. The second is where rays converge: foci and caustics.

A caustic is the envelope of a family of rays [2]: the curve or surface joining the points where neighbouring rays meet. Born and Wolf call it the focal surface of the ray congruence [1]. By the intensity law of Section 3, the cross-section of a ray tube shrinks to zero there, so the intensity should be infinite. Figure 7 shows a caustic formed by sunlight reflecting off the inside of a gold ring onto paper. It is the same kind of cusped bright curve often seen inside a teacup.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/ext-ring-caustic.jpg" alt="A bright curve formed by sunlight reflected inside a gold ring" width="520">
_Figure 7. A reflection caustic formed by the inner surface of a gold ring. Sunlight entering from one side reflects off the curved wall and traces a bright curve with a sharp point, the cusp. (Source: Bautsch, [Katakaustik.Goldring](https://commons.wikimedia.org/wiki/File:Katakaustik.Goldring.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

Computing the same situation gives Figure 8(a). Parallel light reflects once off the inside of a circular mirror of radius $R$. The numerically computed envelope of the reflected rays agrees with the closed form of the nephroid, $x = \tfrac{R}{4}(3\cos t - \cos 3t)$, $y = \tfrac{R}{4}(3\sin t - \sin 3t)$, to about $10^{-5}R$. The cusp lies $R/2$ from the centre, at the paraxial focus of a spherical mirror. Notice too that the caustic floats in the space inside the ring rather than lying on the wall where the reflection happens. Berry makes the same point to argue that the caustics of a mirage generally do not coincide with the heights at which the rays turn back [5]. The boundary that separates the heights in a distant scene that appear in the mirage from those that do not, which Berry calls the "vanishing line", is set by the points whose caustic passes through the eye.

<img src="/assets/img/posts/geometric-optics-fermat-eikonal/en/fig8-caustics.png" alt="Reflection caustic of a circular mirror, refraction caustic of a biconvex lens, and ray density across the caustic" width="780">
_Figure 8. (a) Reflection caustic inside a circular mirror. The computed envelope (red) coincides with the closed-form nephroid (dashed). (b) Refraction caustic behind an N-BK7 biconvex lens (R = ±1000 mm, thickness 100 mm, diameter 200 mm, wavelength 750 nm). Marginal rays cross the axis at 1044.2 mm, near-axis rays at 1060.1 mm. The black bar marks the circle of least confusion, where the blur is smallest. (c) Rays counted across the plane z = 1053 mm. The count spikes where the caustic passes (±0.2 mm)._

A lens behaves the same way. Figure 8(b) sends parallel light into an N-BK7 biconvex lens 200 mm in diameter. The glass index is the catalogue value, which is defined relative to air, so the air around the lens is given an index of exactly 1. Near-axis rays cross the axis at 1060.05 mm, which matches the thick-lens formula [4] to $10^{-11}$ mm. Marginal rays cross at 1044.21 mm, 15.8 mm earlier. This is spherical aberration. The focal shift grows roughly as the square of the ray height [4]; for this lens, the shift divided by the height squared stays constant to within 1.4% from 10 mm to 100 mm. The envelope of the rays forms the caustic, and its cusp sits at the paraxial focus [1]. At 1048.15 mm, where the marginal ray meets the opposite branch of the caustic, the cross-section of the bundle is smallest, 0.81 mm across [2]. The position of this circle of least confusion agrees with a direct minimization of the blur to within $10^{-5}$ mm.

Figure 8(c) counts the heights at which rays cross the plane z = 1053 mm, between the two foci. The ray density spikes at ±0.2 mm, where the caustic passes. According to geometrical optics, the intensity there is infinite. Real light is not. Studying the rainbow in 1838, Airy showed that near a caustic the wave oscillates on one side and dies away exponentially on the other [6]; the function describing this is now called the Airy function. The intensity on the caustic is finite but grows as the wavenumber $k$ to the power 1/3, while the fringe spacing shrinks as $k^{-2/3}$ [6]. As the wavelength goes to zero, the intensity goes to infinity and the fringe spacing to zero, recovering the geometrical picture. A perfect focus, where a spherical wave collapses to a point, is more extreme still: there the intensity grows as $k^2$ [6]. As Hecht puts it, nature abhors infinities [2]. Geometrical optics fails exactly where it predicts the brightest light, and that is why the shape of a focus calls for wave optics.

## Summary

- Sending the wavelength to zero in the wave equation leaves the eikonal equation $\lvert\nabla S\rvert = n$. The scalar wave, Born and Wolf's vector derivation and Orfanidis's orthogonality condition all arrive at it. Surfaces of constant $S$ are wavefronts, the curves perpendicular to them are rays, and rays are the paths along which energy flows.
- By the ray equation, a ray always bends towards the higher index. A level ray is no exception, so Raman's objection does not stand. The index difference in the air above a road is only $1.6\times10^{-5}$, yet for lines of sight grazing the road at less than 0.327° the rays turn back up and produce a sky-reflecting pool about 290 m ahead.
- Snell's law and the law of reflection at an interface re-emerge from the vanishing curl of the ray vector, with the proviso that the radius of curvature be much larger than the wavelength.
- The precise form of Fermat's principle is that the optical path length is stationary, and its Euler–Lagrange equation is the ray equation. A ray is a minimum only until it passes a conjugate point. Stationarity is required because it is the condition for constructive interference with neighbouring paths.
- Geometrical optics breaks down at shadow edges, foci and caustics. In wave terms, the intensity on a caustic is finite, scaling as the wavelength to the power −1/3, and it is decorated by Airy-function fringes.

In a homogeneous medium a ray is a straight line. Everything interesting happens at interfaces, and lens design starts by linearizing Snell's law at those interfaces with $\sin\theta \approx \theta$. The next post covers paraxial optics and ray transfer (ABCD) matrices. We will see the paraxial focus of 1060.05 mm, found here by real-ray tracing, come straight out of a product of a few 2×2 matrices.

## References

1. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Chapter 3 — vector derivation of the eikonal equation, the intensity law, limits of validity of geometrical optics, the ray equation and ray curvature, interface laws, and the minimum property of Fermat's principle with conjugate points.
2. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Sections 4.5 and 6.3 — Hero and Fermat, optical path length in layered media, mirages and the Gladstone–Dale relation, stationary paths seen through phasors, the three cases of the elliptical mirror, spherical aberration and the caustic, and the circle of least confusion.
3. S. J. Orfanidis, *Electromagnetic Waves and Antennas* (Rutgers University, online textbook, 2004), Sections 6.9–6.10 — the eikonal equation from the orthogonality condition, the energy transport velocity, and the ray equation from the Euler–Lagrange equations.
4. G. R. Fowles, *Introduction to Modern Optics*, 2nd ed. (Holt, Rinehart and Winston, 1975), Sections 1.1 and 10.2 — Newton's and Huygens's views of light, the thick-lens formula, and spherical aberration of a single lens.
5. M. V. Berry, "Raman and the mirage revisited: confusions and a rediscovery," *Eur. J. Phys.* **34**, 1423–1437 (2013) — Raman's level-ray objection and its error, the singularity in the discrete-layer limit, mirage caustics and the vanishing line. [Author's copy](https://michaelberryphysics.wordpress.com/wp-content/uploads/2013/06/berry465.pdf)
6. M. V. Berry and C. Upstill, "Catastrophe optics: morphologies of caustics and their diffraction patterns," *Progress in Optics* **18**, 257–346 (1980) — classification of caustics, the Airy function at a fold caustic, and wavelength scaling of intensity and fringe spacing. [Author's copy](https://michaelberryphysics.wordpress.com/wp-content/uploads/2022/02/berry089-1.pdf)
7. A. T. Young, "An Introduction to Mirages" — the angular size of mirages and why the index gradient must change with height. [Web page](https://aty.sdsu.edu/mirages/mirintro.html)
8. P. de Fermat, *Œuvres de Fermat*, t. II, P. Tannery and C. Henry, eds. (Gauthier-Villars, 1894), Letters LXXXVI (1657) and CXII (1662).

---
layout: post
date: 2026-11-13 20:00:00 +0900
title: "Geometrical Optics Foundations 2 — Paraxial Optics and Ray-Transfer Matrices: A Lens as One 2×2 Matrix"
lang: en
lang-exclusive: ["en"]
permalink: /posts/geometric-optics-paraxial-abcd/
page_id: geometric-optics-paraxial-abcd
categories: [Optics, Geometrical Optics]
tags: [geometric-optics, paraxial-optics, ray-transfer-matrix, cardinal-points, grin]
description: "Linearize Snell's law and a lens becomes a 2×2 matrix. This post follows why that works and what the matrix holds."
math: true
---

[Part 1](/en/posts/geometric-optics-fermat-eikonal/) ended with a promise. The paraxial focus of 1060.05 mm, found by tracing real rays through a 200 mm N-BK7 biconvex lens, would come straight out of a product of a few 2×2 matrices. Keeping that promise is where this post begins.

Real optical systems rarely stop at a single lens. Figure 1 shows a camera lens with a thin sheet of 405 nm laser light sent straight through it, so that the section the beam crosses lights up in fluorescence. Each piece of glass fluoresces in a different colour, which reveals that one "lens" is really a stack of several elements.

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/ext-pentax-lens-fluorescence.jpg" alt="A camera lens whose glass elements fluoresce along a laser line sent through it" width="640">
_Figure 1. A Pentax 55 mm F1.4 lens with a 405 nm laser line sent through it. The lens is not cut; only the section along the beam glows. Each element fluoresces in a different colour, exposing the multi-element design. The author stretched the image horizontally to correct the perspective. (Source: yellowcloud from Germany, [Laser Fluorescence Cutaway of Pentax Lens - Stretched](https://commons.wikimedia.org/wiki/File:Laser_Fluorescence_Cutaway_of_Pentax_Lens_-_Stretched.jpg), [CC BY 2.0](https://creativecommons.org/licenses/by/2.0))_

Designing such a system means knowing how rays cross dozens of surfaces. But what if, however many surfaces there are, the relation between the entrance and the exit could be summed up in four numbers? If the focus, the image position and the magnification all followed directly from those four numbers, design would become far lighter. The regime where this summary works is paraxial optics, and the four numbers form the ray-transfer matrix, usually called the ABCD matrix. This post follows, in turn, why the matrix works, what it contains and where it breaks down.

## 1. There Is No Perfect Imaging

What we want from a lens is an image. Every ray from one point of the object should meet at one point of the image; a straight line on the object should stay straight in the image; a plane of the object should map onto a plane of the image. Born and Wolf show that these three demands alone fix the form of ideal imaging. A transformation that sends points to points, lines to lines and planes to planes, the way straight lines stay straight in a perspective drawing, is called a projective transformation (collineation). Add symmetry about the optical axis, and special points such as the foci and principal points follow from the properties of the transformation alone, before anything is known about the lens [3].

The question is whether a real lens can live up to this ideal. For an object of finite size, it cannot. The easiest way to see why is a fact we will compute in Section 6: a lens that shrinks width and height by a factor $m$ shrinks depth by $m^2$. A lens that halves the width quarters the depth. The image of a cube is then not a cube but a flattened box, and no imaging can carry the whole of three-dimensional space into a similar copy. Landau and Lifshitz prove this in general. The only optical systems that image a finite object perfectly at every point are trivial ones, like a plane mirror, which only moves or flips the object [4].

There is still a regime where imaging is almost perfect: narrow bundles of rays travelling close to the axis of an axially symmetric system [4]. Near the axis, a lens surface looks like a small bowl curved the same way in every direction, so a narrow bundle leaving one point comes back together at one point. Paraxial optics is the study of this regime. It is therefore not an approximation picked for convenience; it is the only regime in which ideal imaging actually holds. How far a system strays from it is what we call aberration, the subject of Part 4.

## 2. The Paraxial Approximation

What does "close to the axis" mean in an equation? It means the angle between a ray and the axis is small. For small angles an arc and its chord are nearly equal, so $\sin\theta \approx \tan\theta \approx \theta$ [1, 2]. Written out in full, $\sin\theta = \theta - \theta^3/6 + \cdots$, and the paraxial approximation drops everything from $\theta^3$ on. Snell's law $n_1\sin\theta_1 = n_2\sin\theta_2$ then becomes $n_1\theta_1 = n_2\theta_2$. The angles are simply proportional to each other; the relation has become linear. That linearity is the key to everything in this post.

How good is the approximation? Figure 2 plots the relative error of $\sin\theta \approx \theta$ against angle, together with the paraxial error in the refraction angle for light entering glass ($n = 1.5$) from air. $\sin\theta \approx \theta$ is off by 1% at 14° and by 10% at 43°. The refraction angle is a little more forgiving: 1% at 18.6°.

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/en/fig2-paraxial-error.png" alt="Relative error of the paraxial approximation against angle" width="640">
_Figure 2. Relative error of the paraxial approximation. Blue: sin θ ≈ θ. Red: the refraction angle for light entering glass from air. The dashed line marks 5.7°, the angle at which the marginal ray of the Part 1 lens meets its front surface._

In the Part 1 lens, the marginal ray meets the front surface at only 5.7°, where $\sin\theta \approx \theta$ is wrong by just 0.17%. Yet the marginal ray came to a focus 15.8 mm ahead of the paraxial focus. How does such a small error grow so large? Errors from the two surfaces add up, and the marginal ray leaves the lens 1.5% steeper than the paraxial prediction. The distance a ray travels before reaching the axis is inversely proportional to its slope, so on slope alone a ray that should have travelled 960 mm would cross the axis 1.5%, about 14 mm, early. It also leaves 0.3 mm higher than predicted, which gives back about 3 mm, so it actually crosses the axis about 11 mm early. On top of that, the curved back surface places the marginal ray's exit point 4.7 mm further forward, and together these make 15.8 mm. A small angular error is magnified by the long distance it travels. The systematic account of such deviations, which come from the dropped $\theta^3$ terms, is Seidel aberration theory.

## 3. Two Matrices: Travel and Refraction

In the paraxial regime a ray is written as two numbers: its height $y$ at a reference plane and its angle $u$ to the axis, much as a car on a road is described by its position and heading. A ray also meets only two kinds of events in an optical system: straight-line travel through a uniform medium, and refraction at a boundary.

Take travel first. Over a distance $d$, the direction stays the same and the height changes by (distance) × (slope), just as a car climbing a slope rises in proportion to the distance it covers.

Refraction is the opposite. At the instant a ray crosses a lens surface, its height stays the same and only its direction changes. By how much? A spherical surface tilts more the further you go from the axis. At height $y$ the surface normal makes an angle of about $y/R$ with the axis, so applying the paraxial Snell's law gives a bend proportional to height. Writing Part 1's condition — the tangential component of the ray vector is continuous — in paraxial form gives [5]

$$ n_2u_2 = n_1u_1 - \Phi\,y, \qquad \Phi = \frac{n_2 - n_1}{R} $$

where $\Phi$ is the power of the surface [3]. This equation is the heart of what a lens does: **the further a ray is from the axis, the harder it is bent towards the axis.** A ray entering parallel to the axis ($u_1 = 0$) at height $y$ leaves at $u_2 = -\Phi y/n_2$, so it must travel $y/(\Phi y/n_2) = n_2/\Phi$ to reach the axis. That distance does not contain $y$. Whatever height a ray enters at, it meets the axis at the same distance, and that point is the focus. The reason a focus exists at all is packed into one line: the bend is proportional to height.

Both events only multiply the two incoming numbers by constants and add them, so each can be written as a single 2×2 matrix. Write the ray as $(y, v)$, where $v = nu$ is the reduced angle, the angle multiplied by the refractive index [2]. Why the reduced angle? Because of the paraxial Snell's law from Section 2, $n_1u_1 = n_2u_2$. The quantity that survives a plane boundary unchanged is $nu$, not $u$, so using it as the variable spares us a conversion every time the medium changes.

$$ T(d) = \begin{bmatrix} 1 & d/n \\ 0 & 1 \end{bmatrix}, \qquad R(\Phi) = \begin{bmatrix} 1 & 0 \\ -\Phi & 1 \end{bmatrix} $$

Each of the four elements has a meaning.

| Element | Meaning |
|---|---|
| $A$ | how much of the incoming height survives as outgoing height |
| $B$ | how much the incoming angle changes the outgoing height |
| $C$ | how much the incoming height changes the outgoing angle (minus the power) |
| $D$ | how much of the incoming angle survives as outgoing angle |

The travel matrix has only $B$ (angle changes height); the refraction matrix has only $C$ (height changes angle). A whole optical system is the product of these matrices in the order the ray meets them. The product is taken from the right, because the matrix of the first surface must act on the ray first [1] — just as you put on socks before shoes and take off the shoes first. For the Part 1 lens there are three matrices: refraction at the front surface, 100 mm of travel through glass, and refraction at the back surface.

$$ M = R(\Phi_2)\,T(t)\,R(\Phi_1) = \begin{bmatrix} 0.96614 & 66.145 \\ -0.0010063 & 0.96614 \end{bmatrix} $$

These four numbers are the entire Part 1 lens. A ray entering parallel to the axis, $(y, 0)$, leaves as $(Ay, Cy)$. The air outside the lens has index 1, so the reduced angle $Cy$ is the actual angle, and this ray crosses the axis $-A/C = 960.05$ mm behind the back surface, or 1060.05 mm from the front surface. That agrees to $10^{-11}$ mm with the value Part 1 obtained by narrowing real rays towards the axis. There we followed rays one by one and solved for their intersections; here we simply multiplied three 2×2 matrices.

Different textbooks write the same matrix in slightly different ways. Saleh and Teich, and Fowles, write a ray as height and ordinary angle $(y, u)$ [1, 6]; the matrix between two media of different index then has determinant $n_1/n_2$. Goodman uses reduced angles $(y, nu)$, as this post does, so that every matrix has determinant 1 [2]. Hecht reverses the order to $(nu, y)$ [5]. According to Hecht, T. Smith formulated a way of handling the ray-tracing equations with matrices in the early 1930s, but it went largely unnoticed for some thirty years before being taken up again in the 1960s [5]. The idea of composing a whole system by multiplying 2×2 arrays surface by surface goes back further, to a 1913 paper by R. A. Sampson [7]. The physics is the same in every convention. In the reduced-angle convention, though, the determinant is always 1, and as the next section shows, that 1 means something.

## 4. Why the Determinant Is 1

The travel and refraction matrices both have determinant 1, so their product — any optical system — does too. Is this an accident of the algebra, or does it tell us something?

A picture makes the meaning clear. Plot each ray as a point on a plane whose horizontal axis is the height $y$ and whose vertical axis is the angle $u$. This plane is called phase space, and a bundle of rays becomes a region in it. The grey square in Figure 3 is a bundle containing every ray with height within ±1 mm and angle within ±20 mrad.

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/en/fig3-phase-space.png" alt="A ray bundle in phase space changing shape under travel, a lens, and focal-plane-to-focal-plane transfer" width="780">
_Figure 3. A ray bundle seen in phase space (horizontal: height, vertical: angle). Grey is before, blue is after, and dots of the same colour are the same ray. (a) After 50 mm of travel, upward rays shift right and downward rays shift left. (b) A thin lens of 50 mm focal length bends high rays down and low rays up. (c) From the front focal plane to the back focal plane, the square turns through 90°. In all three cases the area stays at 80.0 mm·mrad._

Figure 3(a) shows the bundle after 50 mm of travel. Upward-pointing rays ($u > 0$) gain height and downward-pointing ones lose it, so the square is sheared sideways into a parallelogram. Push a deck of cards sideways and its shape skews, but the number of cards — the area — stays the same. Figure 3(b) shows the bundle after a thin lens. This time the angle changes according to height, so the shear is vertical, and again the area is unchanged. The determinant is exactly the factor by which a matrix scales this area. Since travel and refraction are both area-preserving shears, any number of them preserves the area, and the determinant is 1.

This conservation has a practical consequence: the width of a bundle and its spread of angles cannot both be reduced at once. When you focus sunlight with a magnifying glass, the smaller you make the Sun's image, the wider the cone of light converging on it. With the area fixed, shrinking one side forces the other to grow. In Part 3 this area appears under the name étendue, the quantity that limits how much light an optical system can carry.

Writing the area as a formula recovers an old invariant. The area of the parallelogram spanned by two rays $(y_1, v_1)$ and $(y_2, v_2)$, namely $y_1v_2 - y_2v_1$, is the same at every surface. Take one ray from the centre of the object and the other from its edge, and this becomes $n\,y\,u = n'\,y'\,u'$, the Smith–Helmholtz formula [3]. It is also called the Lagrange invariant, and in restricted forms it goes back to Huygens and Cotes [3].

Why is the area conserved? Landau and Lifshitz's explanation returns to the optical path length of Part 1. Treat the optical path length between two points as a function of their positions. Differentiating it with respect to the entry point gives the incoming direction, and with respect to the exit point the outgoing direction [4]. The input and output of a ray are thus tied together by a single function, the optical path length. In mechanics, the position and momentum of a particle are tied together in the same way, and that is why phase-space area is conserved there too (Liouville's theorem). The root is the same as that of $\nabla\times(n\mathbf{s}) = 0$, which Part 1 used to derive the boundary conditions; in integral form it is Lagrange's integral invariant $\oint n\mathbf{s}\cdot d\mathbf{r} = 0$.

Goodman reads the same fact in the language of waves. Dividing the reduced angle by the vacuum wavelength gives the local spatial frequency at that point, so height and reduced angle form a Fourier pair, like position and spatial frequency [2]. Just as enlarging a photograph makes its patterns coarser, stretching one member of the pair shrinks the other, and so the lateral magnification $m_t$ and the angular magnification $m_a$ between two conjugate planes multiply to 1 [2]. Extending this relation beyond the paraxial regime gives Abbe's sine condition [2], which is exactly the condition used in [Ellipsometer Instrumentation 6](/en/posts/ellipsometer-back-focal-plane-micro/) to convert back-focal-plane radius into angle of incidence.

## 5. The System Matrix and the Cardinal Points

Can a thick lens be treated as a single thin lens? A thin lens is the model in which a ray bends once, at one place. In a thick lens a ray bends twice, at the front and back surfaces, and its height changes in between, so the thin-lens model does not fit as it stands. Two imaginary planes fix this. When a ray reaches the front principal plane H, think of it as teleported at the same height to the back principal plane H′, where it bends once as if at a thin lens [2]. The space between the two principal planes effectively does not exist for the ray. A thin lens is just the special case where the two planes coincide.

The principal planes and foci can be found by drawing rays. The point where rays entering parallel to the axis converge is the back focal point F′, and the plane where the extensions of the incoming and outgoing rays meet is the back principal plane H′. The front side works the same way: a ray leaving the front focal point F emerges parallel and locates the front principal plane H. With the matrix in hand, no drawing is needed. With input-side index $n$ and output-side index $n'$, the back focal point lies $-n'A/C$ from the output plane and the front focal point $nD/C$ from the input plane, and each focal length is the reciprocal of the power $-C$ multiplied by the index on that side [1].

For the Part 1 lens, the focal length is 993.70 mm and the front focal point lies 960.05 mm in front of the front vertex. The two principal planes H and H′ each lie 33.64 mm inside the lens from their respective vertices (Figure 4(a)). The lens is equivalent to a thin lens of focal length 993.70 mm with a 32.7 mm "teleport gap" between the two principal planes. Born and Wolf's thick-lens formulas also give the principal-plane positions in closed form [3], and they match the values read from the matrix to $10^{-13}$ mm.

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/en/fig4-cardinal-points.png" alt="Focal points and principal planes of the Part 1 lens, and how the principal plane curves with ray height" width="780">
_Figure 4. (a) Cardinal points of the Part 1 lens. The blue rays enter parallel and converge at F′; the red rays leave F and emerge parallel. The axial scale is compressed. (b) Where the extensions of a real ray entering at height h meet. The higher the ray, the further in front of the paraxial principal plane H′ they meet._

So is the principal plane really a plane? Figure 4(b) plots where the extensions of a real ray entering at height $h$ meet. At $h = 30$ mm the meeting point is only 0.1 mm from the paraxial principal plane, but at $h = 100$ mm it has moved 1.1 mm forward. The places where rays appear to bend lie not on a plane but on a curved surface. As Goodman notes, for large-aperture lenses the principal plane becomes a curved surface [2]. How the shape of that surface connects to the sine condition is left for Part 4.

The position of the principal planes moves a great deal with the shape of the lens (Figure 5). The more the lens is bent to one side, the further the principal planes follow in that direction, and in a crescent-shaped meniscus lens they can even move outside the glass [3].

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/ext-lens-shapes-principal-planes.png" alt="Positions of the principal planes H and H′ for various lens shapes" width="560">
_Figure 5. Lenses of various shapes and the positions of their principal planes H and H′. Lenses 1–4 are converging and 5–8 diverging. In the meniscus lenses 4 and 8, both principal planes lie outside the lens. (Source: OwlTheCat, [Lens shapes 2](https://commons.wikimedia.org/wiki/File:Lens_shapes_2.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

A principal plane outside the lens is useful. Place a concave lens of focal length −50 mm 60 mm behind a convex lens of focal length 100 mm. Both the two-thin-lens formula [3] and the matrix product give a focal length of 500 mm. Yet the actual distance from the first lens to the focus is only 260 mm, because the back principal plane H′ lies 240 mm in front of the first lens, out in empty space. A long focal length fits into a short barrel, and that is the principle of the telephoto lens.

There is an extreme in the other direction as well. A biconvex lens with equal radii has $C = 0$, and hence a focus at infinity, when its thickness reaches $2nr/(n - 1)$ [3]. Parallel light then leaves as parallel light, and the lens has become a telescope. For the Part 1 lens, that thickness would be 5.9 m.

## 6. The Imaging Condition and Conjugate Planes

How is the condition for forming an image written as a matrix? Consider the full matrix from the object plane to the image plane. A ray leaving one point of the object must arrive at the same point of the image plane, whatever angle it left at. The height at the image plane must therefore not depend on the starting angle — in the terms of the table, $B$, "how much the incoming angle changes the outgoing height", must be zero [1, 2]. Two planes paired this way are called conjugate planes. Then $A$ becomes the lateral magnification, the ratio of image size to object size, and $D$ the angular magnification, the ratio of ray angles; since the determinant is 1, $AD = 1$ [2]. Halve the image and the ray angles double. The area conservation of Section 4 shows up again.

Put an object 3000 mm in front of the Part 1 lens, and the condition $B = 0$ places the image 1444.10 mm behind the back surface, with a lateral magnification of −0.487 (the minus sign means the image is inverted). The object and image distances $s$ and $s'$ measured from the principal planes satisfy the Gaussian lens equation $1/s + 1/s' = 1/f$, and the distances $z$ and $z'$ measured from the foci satisfy Newton's equation $zz' = ff'$, both to within numerical precision [1, 3]. Every textbook imaging formula comes out of a single matrix.

Now for the calculation deferred from Section 1. If the object moves a little along the axis, how far does the image move? This ratio is the longitudinal magnification. Moving the object from 3000 mm to 1 mm further away brings the image 0.237 mm closer. The longitudinal magnification 0.237 is the square of the lateral magnification 0.487: width and height shrink by 0.49, depth by 0.24. As Landau and Lifshitz point out, longitudinal and lateral magnification generally differ [4], and unless the magnification has magnitude 1, the image of a cube is a flattened or elongated box. This confirms in numbers, within the paraxial regime itself, the conclusion of Section 1 that there is no perfect three-dimensional imaging.

Among conjugate planes there is one more special pair. The matrix from the front focal plane to the back focal plane has $A = D = 0$ [2]. Here $\Phi$ is the power of the whole lens, $-C$, not of a single surface:

$$ M_{FF'} = \begin{bmatrix} 0 & 1/\Phi \\ -\Phi & 0 \end{bmatrix} $$

This matrix swaps height and angle. Light from one point in the front focal plane leaves the lens as a parallel beam travelling in one direction, and a parallel beam arriving from one direction converges to one point in the back focal plane. It is why, in a camera focused at infinity, each star becomes one point on the sensor. The 90° turn of the square in Figure 3(c) is the same thing. This is why angles of incidence spread out as radii in the back focal plane, and in wave optics it becomes the Fourier-transforming property of a lens.

## 7. Graded-Index Media and Periodic Systems

So far refraction has happened only at boundaries. Can a medium whose index varies continuously, like the air of Part 1's mirage, also be written as a matrix? In the paraxial regime the ray equation reduces to $d^2y/dz^2 = (1/n)\,dn/dy$ [1]. The parabolic medium of Part 1, $n^2 = n_0^2(1 - \alpha^2y^2)$, has its highest index at the centre and lower index further out. Rays bend towards higher index, so a ray that strays from the centre is pulled back towards it — like a ball rolling in the bottom of a bowl, or a mass on a spring, pulled back by a force proportional to its displacement. The equation becomes $d^2y/dz^2 = -\alpha^2y$, the equation of simple harmonic motion, and its solution written as a matrix is

$$ G(z) = \begin{bmatrix} \cos\alpha z & \sin\alpha z/(n_0\alpha) \\ -n_0\alpha\sin\alpha z & \cos\alpha z \end{bmatrix} $$

The hallmark of simple harmonic motion is that the period does not depend on how far the motion swings — the same isochronism that makes the period of a gently swinging pendulum independent of its amplitude. All paraxial rays therefore oscillate with the same period $2\pi/\alpha$, and rays leaving one point come back together every half period. This is the principle of the graded-index (GRIN) rod lens known by the trade name SELFOC [1]. A rod of length $d$ acts as a lens of focal length $1/(n_0\alpha\sin\alpha d)$ [1]. The same matrix also results from slicing the medium into very thin plates, each with a weak lens: cutting it into 4000 slices and multiplying reproduces the continuous matrix to $2\times10^{-8}$.

Then why did the conjugate point in Part 1 appear earlier than the paraxial half period $\pi/\alpha$? A pendulum swung widely has a period that depends on amplitude, and rays are no different. Part 1 solved the ray equation without approximation, and Figure 6(a) compares rays leaving one point of the same medium as given by the ray equation (solid) and the paraxial matrix (dashed). At 5° the two nearly coincide; at 45° they separate widely. The distance at which a real ray returns to the axis is $(\pi/\alpha)\cos\theta_0$, shorter than the paraxial value $\pi/\alpha$ (Figure 6(b)). A pendulum's period grows with amplitude, whereas rays in this medium return sooner. The shortfall is 1.5% at a launch angle of 10°, 13% at 30° and 29% at 45°. The difference between matrix and ray equation scales as the square of the amplitude: reduce the amplitude tenfold and the relative difference drops almost exactly a hundredfold. It is the error introduced by the terms dropped in Section 2.

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/en/fig6-grin-rod.png" alt="Rays in a parabolic graded-index medium and the distance at which they return to the axis" width="780">
_Figure 6. The medium n² = n₀²(1 − α²y²) (n₀ = 1.5, α = 0.5). (a) Rays leaving the origin at 5°, 15°, 30° and 45°. Solid: ray equation; dashed: paraxial matrix. (b) Distance at which a ray first returns to the axis. Dots are integrations of the ray equation, the orange line is (π/α)cos θ₀, and the red dashed line is the paraxial π/α._

If a continuous medium can be seen as a repetition of thin lenses, we can turn the idea around and line up identical lenses at equal spacing. With lenses of focal length $f$ spaced $d$ apart, each lens bends the ray back towards the axis once. With a suitable spacing the ray rises and falls about the axis and stays trapped. With too large a spacing, things change. The ray overshoots the axis and travels far to the other side before reaching the next lens, which, meeting it at a greater height, bends it harder and throws it even further to the opposite side. Repeated, this sends the ray away.

Where is the boundary? Saleh and Teich form $b = (A + D)/2$ from the matrix of a single stage and show that the ray oscillates boundedly when $\lvert b\rvert \le 1$ and grows exponentially otherwise [1]. When it oscillates, its phase advances by $\cos^{-1}b$ per stage. For the lens chain, $b = 1 - d/2f$, so the condition is $0 \le d \le 4f$ [1]. At the boundary $d = 4f$, however, the ray grows linearly from stage to stage, so strictly the ray is confined only for $0 < d < 4f$ [6].

<img src="/assets/img/posts/geometric-optics-paraxial-abcd/en/fig7-periodic-guide.png" alt="Stability condition of a periodic lens chain and the ray height at each stage" width="780">
_Figure 7. Identical lenses of focal length f spaced d apart. (a) b = 1 − d/2f and the stable region (green). (b) Height at each stage of a ray starting at height 1 parallel to the axis. For d = f and d = 2f the ray returns to its start every 6 and 4 stages; for d = 3.5f it stays bounded but never repeats; for d = 4.2f it diverges at once._

If the phase advance per stage, $\cos^{-1}b$, is a rational fraction of a full turn $2\pi$, the ray returns exactly to where it started after a certain number of stages [1]. It returns every 6 stages for $d = f$ and every 4 for $d = 2f$, and raising the matrix to that power gives the identity (Figure 7(b)). The same calculation applies to a ray bouncing between two mirrors, so the same condition decides whether a laser resonator can trap light [6]. When the wave optics series reaches Gaussian beams, the same matrix returns as the tool that carries the beam size and wavefront curvature through a system.

## Summary

- No optical system images a finite object perfectly except in trivial cases, because shrinking width by $m$ shrinks depth by $m^2$. Paraxial optics is the only regime in which ideal imaging holds; departures from it are aberrations.
- In the paraxial regime a lens is a linear device that "bends rays harder the further they are from the axis", and an optical system is a product of two kinds of 2×2 matrices, travel and refraction. The Part 1 lens reduces to a product of three matrices, whose focus of 1060.05 mm matches real-ray tracing to $10^{-11}$ mm.
- In the reduced-angle convention every matrix has determinant 1, because travel and refraction are area-preserving shears in phase space. That is why the width and angular spread of a bundle cannot both be reduced, and it is the origin of the Smith–Helmholtz (Lagrange) invariant and of $m_t m_a = 1$ between conjugate planes.
- The foci, principal planes, imaging condition and magnifications are all combinations of $A, B, C, D$. Pushing the rear principal plane out in front of the lens makes a telephoto lens; at large apertures the principal plane becomes curved.
- Continuous graded-index media and periodic lens chains are handled with the same matrices. The difference between the paraxial matrix and real rays grows as the square of the amplitude.

This post covered only how rays pass through a system. Real lenses have rims and stops, and some rays never get through at all. The next post covers stops and fields, which decide which rays pass. Its subjects are how the entrance and exit pupils, the chief and marginal rays, the numerical aperture and the f-number are defined on top of the matrix, and how the phase-space area met here becomes étendue and limits how much light a system can carry.

## References

1. B. E. A. Saleh and M. C. Teich, *Fundamentals of Photonics*, 3rd ed. (Wiley, 2019), Sections 1.3–1.4 — the paraxial ray equation, parabolic graded-index slabs and GRIN lenses, ray-transfer matrices, formulas for the cardinal points, stability and periodic trajectories of periodic systems.
2. J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed. (McGraw-Hill, 1996), Appendix B — the reduced-angle convention with unit determinant, reduced angle and local spatial frequency, magnification relations between conjugate planes and the sine condition, the focal-plane-to-focal-plane matrix, principal planes and their curvature.
3. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Sections 4.3–4.4 — ideal imaging as a projective transformation, Newton's equation, surface power, thick-lens cardinal points and the telescopic condition, combining two thin lenses, the Smith–Helmholtz formula.
4. L. D. Landau and E. M. Lifshitz, *The Classical Theory of Fields*, 3rd rev. English ed. (Pergamon, 1971), §§55–56 — optical path length between two points and the angular eikonal, the impossibility of perfect imaging, imaging with narrow bundles, longitudinal and lateral magnification.
5. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Sections 6.1–6.2 — thick lenses and principal planes, the paraxial refraction equation, the history of the matrix method and ray-vector conventions.
6. G. R. Fowles, *Introduction to Modern Optics*, 2nd ed. (Holt, Rinehart and Winston, 1975), Sections 10.3–10.5 — ray equations and ray matrices, periodic lens waveguides and optical resonators.
7. R. A. Sampson, "A new treatment of optical aberrations," *Phil. Trans. R. Soc. Lond. A* **212**, 149–185 (1913). <https://doi.org/10.1098/rsta.1913.0005> — an early case of writing surfaces and gaps as 2×2 arrays and composing a system from their product.

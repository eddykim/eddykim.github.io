---
title: "Ray Tracing 3 — Why One Formula Is Enough for Reflection"
lang: en
lang-exclusive: ["en"]
permalink: /posts/raytracing-reflection-family/
page_id: raytracing-reflection-family
date: 2026-09-09 20:00:00 +0900
categories: [Optics, Geometrical Optics]
tags: [ray-tracing, geometric-optics, reflection, python, matlab]
description: Why a single reflection vector formula serves the flat mirror, spherical mirror, beam splitter and polygon mirror alike, and how far the reflected ray drifts when a vertex normal is computed wrongly.
math: true
---

[Post 2](/en/posts/raytracing-total-internal-reflection/) closed by noting that the reflection vector formula appears six times across the library: in `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter`, and in the total internal reflection (TIR) branches of `Prism` and `SphericalLens`. This post asks why one formula suffices, and looks at how it is used in the four reflective components: the flat mirror, the spherical mirror, the beam splitter and the arbitrary mirror.

## 1. Why does no refractive index enter reflection?

Refraction takes the refractive indices of both media as inputs; Snell's law itself is $n_1 \sin\theta_1 = n_2 \sin\theta_2$. On top of that, the lens, the prism and the TIR branch of this library each computed the refraction angle in their own way.

Reflection is different. It says only that the angle of incidence equals the angle of reflection, and no refractive index appears. Split the incident vector $\vec v$ into a component along the normal $\hat n$ and one along the tangent:

$$ \vec v = (\vec v \cdot \hat n)\hat n + \big[\vec v - (\vec v \cdot \hat n)\hat n\big] $$

Reflection keeps the tangential component and flips the sign of the normal one.

$$ \vec v' = -(\vec v \cdot \hat n)\hat n + \big[\vec v - (\vec v \cdot \hat n)\hat n\big] = \vec v - 2(\vec v \cdot \hat n)\hat n $$

Here $\vec v$ is the incident direction and $\hat n$ the unit normal of the surface. A ball bouncing off a wall behaves the same way: the velocity along the wall is unchanged, and only the velocity into the wall reverses. The factor 2 is there because the normal component is subtracted once to cancel it and once more to point it the other way. The limits behave as expected. At normal incidence ($\vec v=-\hat n$) the result is $\vec v'=\hat n$, so the ray goes back the way it came; at grazing incidence ($\vec v\cdot\hat n=0$) the direction does not change.

For the same reason an error in the normal is doubled. Rotate the normal by $\delta$ with the incident ray fixed, and the reflected ray turns by $2\delta$. Turning the normal by 0.5° and 1° turns the reflected ray by 1.000° and 2.000°. Section 3 shows this doubling in practice.

No property of the medium appears anywhere. A mirror, total internal reflection inside a prism, total internal reflection at the exit face of a lens: this one expression covers them all.

```python
def reflect_vector(v, n):
    """Perfect reflection about the normal n. v and n are (vx, vy); n is a unit vector.

    Passing a sign-flipped n gives the same result -- (v.n) flips with it and
    the two sign changes cancel.
    """
    v = np.asarray(v, dtype=float)
    n = np.asarray(n, dtype=float)
    return v - 2 * np.dot(v, n) * n
```

The sign of the normal does not matter either, so all four reflective components take the direction straight from the boundary geometry without orienting it. The only routine that does orient its normal is the prism of [post 2](/en/posts/raytracing-total-internal-reflection/), and not because of reflection. The prism first uses that same normal to get the Snell's law angle of incidence $\theta_1 = \arccos(\vec v \cdot \hat n)$, which becomes its supplement if the normal points the other way. It is refraction that needs a signed normal, not reflection.

Because there is only one formula, stacking reflections stays simple. Two mirrors at right angles flip the two components one after the other, giving $\vec v' = -\vec v$: the ray returns exactly the way it came, whatever its angle of incidence. A 90° V-shaped `ArbitraryMirror` hit by rays tilted 0°, 10° and 25° returns each of them after two reflections with $\vec v_{out}\cdot\vec v_{in} = -1.000000000000$. Open the V to 100° and the returning ray is off by 20°, twice the 10° departure from a right angle. In three dimensions, three mutually perpendicular faces do the same job; that is the corner-cube retroreflector.

## 2. How does each mirror find its normal?

With the formula fixed, the only thing that differs between routines is how the normal is found.

- `FlatMirror`: the mirror is a straight line, so the normal is a single constant: the line coefficients `(a_e, b_e)` scaled to unit length.
- `SphericalMirror`: the mirror is a circular arc, so the normal is the radial direction from the centre of curvature `(cx, cy)` to the intersection point.
- `BeamSplitter`: the normal is computed exactly as for `FlatMirror`. Even the constructor is shared: `make_beam_splitter` calls `make_flat_mirror` to build the same planar shape. That is reasonable, since a plate beam splitter is physically a semi-transparent flat mirror. The difference is that it returns a transmitted ray, continuing in the incident direction, alongside the reflected one.
- `ArbitraryMirror`: the boundary is a polygon, a chain of segments, so a hit on an edge takes the direction perpendicular to that edge. A hit exactly on a vertex has two candidate edges and needs separate handling (section 3).

All four were checked for the law of reflection about the normal. The normal used in the check is computed afresh from the boundary geometry rather than taken from inside the component, so that a component using the wrong normal would be caught. The spherical mirror is the exception: its check uses the component's centre of curvature, so it gets a second, independent test below.

```python
def ang_err(vin, vout, n):
    ang_in = np.rad2deg(np.arccos(abs(np.dot(vin, n))))
    ang_out = np.rad2deg(np.arccos(abs(np.dot(vout, n))))
    return ang_in, ang_out, ang_in - ang_out
```

```text
FlatMirror             incidence=20.000000°  reflection=20.000000°  diff=0.0e+00°
SphericalMirror        incidence=7.662256°  reflection=7.662256°  diff=-1.9e-13°
BeamSplitter (refl.)   incidence=45.000000°  reflection=45.000000°  diff=3.6e-14°
ArbitraryMirror (edge) incidence=26.565051°  reflection=26.565051°  diff=-1.4e-14°
```

The incidence angles can be worked out by hand. The flat mirror is tilted 20° from vertical, so a horizontal ray meets it at 20°. On the spherical mirror a horizontal ray at height 40 mm strikes a circle of radius 300 mm, giving $\arcsin(40/300) = 7.662256°$. The beam splitter is a 45° plate, and the polygon edge has slope 1/2, giving $\arctan(1/2) = 26.565051°$. In every case the reflection angle matches to within floating-point error. Fig 1 shows the four paths.

<img src="/assets/img/posts/raytracing-reflection-family/en/fig1-reflection-family.png" alt="Ray paths through the four reflective components" width="700">
_Fig 1. One ray traced through each reflective component. Top left: a flat mirror tilted 20° from vertical. Top right: a concave spherical mirror of radius 300 mm and diameter 200 mm, hit at a height of 40 mm. Bottom left: a 45° beam splitter (red reflected, blue transmitted). Bottom right: a ray falling vertically onto an edge of a zigzag polygon mirror. All four use the same reflection formula; only the way of finding the normal differs._

Equal angles are a necessary condition, not a sufficient one. `ang_err` compares magnitudes only, so a ray that passed straight through without reflecting would also score zero. The spherical mirror was therefore checked against theory as well. Parallel rays close to the axis of a concave spherical mirror should cross it at the paraxial focus, a distance $R/2$ from the vertex: 150 mm for a 300 mm radius (`verify_reflection.py`).

```text
spherical mirror R=300mm, paraxial focus R/2=150.0mm
height h=   1mm: crosses axis at 149.999mm  closed form=149.999mm
height h=  10mm: crosses axis at 149.917mm  closed form=149.917mm
height h=  40mm: crosses axis at 148.649mm  closed form=148.649mm
height h=  95mm: crosses axis at 141.862mm  closed form=141.862mm
```

A ray 1 mm off axis crosses at 149.999 mm, in agreement with $R/2$. Higher rays cross closer to the mirror; at the rim (95 mm) the crossing is at 141.862 mm, 8.1 mm short. This is spherical aberration. With $\sin\theta = h/R$, a ray at height $h$ crosses the axis at $R - R/(2\cos\theta)$, and the simulation matches this closed form to the third decimal place. Marginal rays focusing short is the same behaviour seen for the biconvex lens in [post 1](/en/posts/raytracing-spherical-lens-refraction/).

## 3. How is the normal chosen at a polygon vertex?

`ArbitraryMirror` represents even a smooth curved mirror as a polygon, which raises a practical question first: how many edges are enough? A ray hitting an edge uses that edge's perpendicular, which differs from the true surface normal by at most half the arc spanned by one edge. Reflection doubles that, so the reflected direction is off by at most one edge's arc. The spherical mirror of Fig 1 (total arc 38.942°) was cut into $N$ edges, and 381 parallel rays were compared with the smooth `SphericalMirror` (`verify_rules_of_thumb.py`).

```text
N= 10 edges: arc per edge=3.894°  max reflected-direction error=3.871°
N= 20 edges: arc per edge=1.947°  max reflected-direction error=1.935°
N= 39 edges: arc per edge=0.999°  max reflected-direction error=0.997°
N=100 edges: arc per edge=0.389°  max reflected-direction error=0.388°
N=390 edges: arc per edge=0.100°  max reflected-direction error=0.100°
```

To keep the reflected direction within $\varepsilon$, then, use at least (total arc)/$\varepsilon$ edges: 39 for 1° on this mirror, 390 for 0.1°.

The vertex is the awkward case. Each edge has its own perpendicular, so at a vertex it is unclear which one applies. The routine therefore fits a circle through three adjacent points (previous point, vertex, next point) and takes the radial direction joining the vertex and the circle's centre as the normal; the code points it from the vertex towards the centre. Three non-collinear points determine a circle uniquely, so local curvature is approximated without any extra parameter. If the points were sampled from a circle, the fitted circle is that circle: on a 20-edge polygon of the mirror above, the vertex normal agrees with the true radius to 0.0°.

```python
def _vertex_normal(x0, y0, x1, y1, x2, y2, xc, yc, dup_point):
    """Local-curvature normal at the vertex (xc,yc), from the circle through 3 adjacent points."""
    M = np.array([[x0, y0, 1], [x1, y1, 1], [x2, y2, 1]], dtype=float)
    if dup_point:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x1**2 + y1**2], dtype=float)
    else:
        P = -np.array([x0**2 + y0**2, x1**2 + y1**2, x2**2 + y2**2], dtype=float)
    J = np.linalg.solve(M, P)
    acx, acy = -J[0] / 2, -J[1] / 2
    n = np.array([acx - xc, acy - yc])
    return n / np.linalg.norm(n)
```

Writing the circle as $x^2+y^2+J_0 x+J_1 y+J_2=0$ and substituting the three points gives a linear system: each row of `M` holds one point's coordinates and the matching entry of `P` holds that point's $-(x^2+y^2)$. Transcribing a system like this invites a familiar mistake: copy a line, paste it, and forget to change one index. The `dup_point=True` branch reproduces that mistake. It puts the second point `(x1, y1)`, the vertex itself, into the third entry of `P` instead of the third point `(x2, y2)`. How much can one wrong index move the reflected ray?

The system does not become degenerate. The third row of `M` is still `(x2, y2)`, so the matrix stays non-singular and `solve` returns an answer without a warning. The mistake therefore produces plausible numbers and no error message. What breaks is the premise that both sides of the third equation refer to the same point. The resulting circle passes through the first two points, the previous point and the vertex, and misses only the third. The two versions were compared while making the adjoining edges increasingly unequal in length and slope (`verify_vertex_normal.py`).

```text
symmetric vertex : n_dup=[-0.210, -0.978]  n_uniq=[0.000, -1.000]  diff=12.095deg  third point missed by 35.0mm
asymmetric vertex: n_dup=[-0.057, -0.998]  n_uniq=[0.310, -0.951]  diff=21.319deg  third point missed by 76.1mm
more asymmetric  : n_dup=[0.161, -0.987]  n_uniq=[0.685, -0.728]  diff=33.983deg  third point missed by 143.4mm
```

In these coordinates the discrepancy grows with asymmetry, 12° to 21° to 34°, and the distance by which the corrupted circle misses the third point grows with it, 35 mm to 76 mm to 143 mm. The ordering, however, is not set by the shape of the vertex. With the right-hand sides of the second and third equations now equal, subtracting them forces the centre $\vec c$ to satisfy $\vec c\cdot(\vec p_2-\vec p_1)=0$. That pins the centre to a line through the coordinate origin, so the same mirror gives a different error depending on where it sits in the coordinate system. Shifting the three vertices by 500 mm in $y$ reverses the trend: 17.5° to 13.6° to 9.1°. A ray sent perpendicularly onto a symmetric vertex ought to retrace its own path (Fig 2). With the normal from the reused third point it leaves 24° off course instead, twice the 12° error in the normal, the doubling of section 1.

<img src="/assets/img/posts/raytracing-reflection-family/en/fig2-vertex-normal-typo.png" alt="Effect of the vertex normal calculation" width="600">
_Fig 2. A ray falls vertically onto the symmetric vertex (0, 50), whose neighbours are (−100, 0) and (100, 0). Reusing the third point (red solid) tilts the normal by 12.1°, which turns the reflected ray by twice that, 24.2°. Using all three points (blue dashed) retraces the incident path, so the blue dashed line hides the grey incident ray._

## 4. What does this implementation leave out?

The reflection formula is independent of the medium, and one copy of it suffices, but the accuracy rests entirely on how the normal is found. As `ArbitraryMirror` shows, an edge hit and a vertex hit compute the normal in quite different ways, and the result is correspondingly sensitive to how that is done.

Reflectance is not modelled either. Even `BeamSplitter` emits the reflected and transmitted rays as two full-strength rays. A real beam splitter divides the intensity according to its reflection-to-transmission ratio, such as 50:50. The one in Fig 3 is aluminium-coated to reflect 80%.

<img src="/assets/img/posts/raytracing-reflection-family/ext-plate-beamsplitter.jpg" alt="A green laser beam splitting at an aluminium-coated glass plate" width="600">
_Fig 3. A green laser on a plate beam splitter aluminium-coated to reflect 80%. The beam entering from the left splits at the single point where it meets the plate. The plate has a visible thickness, whereas the code treats it as a line segment of zero thickness. (Source: Zaereth, [Flat metal-coated beamsplitter](https://commons.wikimedia.org/wiki/File:Flat_metal-coated_beamsplitter.png), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

Ignoring that thickness has a cost of its own: a plate beam splitter displaces the transmitted ray sideways. For glass of index 1.5168 (N-BK7 at 587.6 nm) at 45° incidence, the shift is 0.334 mm for a 1 mm plate and 1.003 mm for a 3 mm plate (`verify_rules_of_thumb.py`). The code sets this shift to zero and continues the transmitted ray along the incident direction. A cube beam splitter, two right-angle prisms cemented together, avoids the shift because light enters its faces at normal incidence.

Parts of the library remain untouched by this series. `Aperture`, `Stops` and `FlatSensor` neither reflect nor refract; they block rays or merely record intensity, which puts them in a different category from the reflective family. `CompoundSphericalLens`, the N-surface lens, generalizes `SphericalLens` and repeats principles already covered.

## Summary

The reflection formula $\vec v' = \vec v - 2(\vec v\cdot\hat n)\hat n$ contains neither the medium nor the sign of the normal. The flat mirror, spherical mirror, beam splitter and polygon mirror therefore share it, and differ only in how they find the normal. The price is that an error of $\delta$ in the normal becomes $2\delta$ in the reflected ray. When a curved mirror is entered as a polygon, the arc of one edge bounds the error in the reflected direction; one wrong index in the vertex calculation cost 24° in this example. The accuracy of a reflection is the accuracy of its normal.

Each post in this series checked one calculation by two independent routes. [Post 1](/en/posts/raytracing-spherical-lens-refraction/) checked refraction at a spherical surface against optical reversibility and the thick-lens formula; [post 2](/en/posts/raytracing-total-internal-reflection/) checked total internal reflection against the law of reflection and the critical angle; this post checked reflection by comparing angles and by the paraxial focus $R/2$ of a spherical mirror. Why light travels in straight lines and bends only at boundaries, where the law of reflection comes from, and how stops and aberrations shape an image are the subject of the Geometrical Optics Background series.

## References

- E. Hecht, *Optics*, 5th ed. (Pearson, 2017), §4.3 (law of reflection) and §5.4 (moving mirrors: a mirror turned by α turns the reflected ray by 2α).
- A. S. Glassner (ed.), *An Introduction to Ray Tracing* (Academic Press, 1989): the reflection formula in vector form. [Free PDF](https://www.realtimerendering.com/raytracing/An-Introduction-to-Ray-Tracing-The-Morgan-Kaufmann-Series-in-Computer-Graphics-.pdf)
- M. Pharr, W. Jakob, G. Humphreys, *Physically Based Rendering*, 4th ed., "Specular Reflection and Transmission". [pbr-book.org](https://pbr-book.org/4ed/Reflection_Models/Specular_Reflection_and_Transmission)

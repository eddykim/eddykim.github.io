---
title: "Geometrical Optics 3 — Why One Formula Is Enough for Reflection"
lang: en
lang-exclusive: ["en"]
permalink: /posts/raytracing-reflection-family/
page_id: raytracing-reflection-family
date: 2026-09-11 20:00:00 +0900
categories: [Optics, Geometrical Optics]
tags: [ray-tracing, geometric-optics, reflection, python, matlab]
description: Why a single reflection vector formula suffices across the whole library, how it is used in the flat mirror, spherical mirror, beamsplitter and arbitrary mirror, and how far the result drifts when the vertex normal of the arbitrary mirror is computed wrongly.
math: true
---

[Post 2](/en/posts/raytracing-total-internal-reflection/) closed by noting that the reflection vector formula appears six times across the library: in `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter`, and in the TIR branches of `Prism` and `SphericalLens`. This post examines why one formula suffices, and how it is actually used in the four reflective components — flat mirror, spherical mirror, beamsplitter and arbitrary mirror.

## 1. Why there is only one reflection formula

Refraction needed a different angular relation every time the medium changed. Snell's law itself, $n_1 \sin\theta_1 = n_2 \sin\theta_2$, takes the refractive indices of both media as inputs. That is why the angle expressions differed slightly for the lens, the prism and the TIR branch.

Reflection is different. It states only that the angle of incidence equals the angle of reflection, and no refractive index enters at all. Decomposing the incident vector $\vec v$ into a component along the normal $\hat n$ and a component along the tangent:

$$ \vec v = (\vec v \cdot \hat n)\hat n + \big[\vec v - (\vec v \cdot \hat n)\hat n\big] $$

Reflection leaves the tangential component alone and flips the sign of the normal component.

$$ \vec v' = -(\vec v \cdot \hat n)\hat n + \big[\vec v - (\vec v \cdot \hat n)\hat n\big] = \vec v - 2(\vec v \cdot \hat n)\hat n $$

No property of the medium appears anywhere. A mirror, total internal reflection inside a prism, TIR at the exit surface of a lens — this one expression covers them all.

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

The sign of the normal does not matter either. All four reflective components therefore take the direction straight from the boundary geometry, with no step to orient it. The one routine that does orient its normal is the prism of [post 2](/en/posts/raytracing-total-internal-reflection/), and not because of reflection: the prism first uses that same normal to obtain the Snell's law angle of incidence $\theta_1 = \arccos(\vec v \cdot \hat n)$, which turns into its supplement if the normal points the other way. It is the refraction that demands a signed normal, not the reflection.

## 2. Four ways to find the normal

If the formula is fixed, the only thing that differs between routines is how the normal is obtained.

- `FlatMirror`: the mirror is a straight line, so the normal is a single constant. Normalizing the line coefficients `(a_e, b_e)` is the whole calculation.
- `SphericalMirror`: the mirror is a circular arc, so the normal is the radial direction from the centre of curvature `(cx, cy)` to the intersection point.
- `BeamSplitter`: the normal is computed exactly as for `FlatMirror`. Even the constructor reuses it — `make_beam_splitter` calls `make_flat_mirror` to build the same planar shape, which is fitting, since a beamsplitter is physically a semi-transparent flat mirror. The difference is that it returns not only the reflected ray but also a transmitted ray continuing in the incident direction.
- `ArbitraryMirror`: the boundary is a polygon, a chain of segments, so a hit on an edge takes the direction perpendicular to that edge. A hit exactly on a vertex is another matter, treated in section 3.

All four were checked for the law of reflection about the normal. The normal used in the check is computed afresh from the boundary geometry rather than taken from inside the component, so that a component using the wrong normal would be caught as well.

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

All four satisfy the law of reflection to within floating-point error.

<img src="/assets/img/posts/raytracing-reflection-family/en/fig1-reflection-family.png" alt="Ray paths through the four reflective components" width="700">
_Fig 1. Ray paths for the four reflective components — one reflection formula, four ways of finding the normal_

## 3. Approximating local curvature at a polygon vertex, and its error

`ArbitraryMirror` represents even a smooth curved mirror by subdividing it into a polygon. The awkward case is the vertex. The perpendicular normal differs from edge to edge, so at a vertex it is unclear which edge's normal applies. The routine therefore fits a circle through the three adjacent points — previous point, vertex, next point — and takes the direction from that circle's centre to the vertex as the normal. Three points determine a circle uniquely, so the local curvature is approximated without any extra parameter.

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

Writing the circle as $x^2+y^2+J_0 x+J_1 y+J_2=0$ and substituting the three points gives a linear system in which each row of the matrix `M` holds one point's coordinates and the matching entry of `P` holds that point's $-(x^2+y^2)$. What `dup_point=True` does is put the second point `(x1, y1)` — the vertex itself — into the third entry of `P`, in place of the third point `(x2, y2)`.

This does not collapse the three points into a degenerate case. The third row of `M` is still `(x2, y2)`, so the matrix stays non-singular and `solve` succeeds without a warning. What breaks is the premise that both sides of the third equation refer to the same point. The circle that comes out passes through the first two points — the previous point and the vertex — and misses only the third. The two cases were compared while varying how unequal the adjoining edges are in length and slope, from symmetric to increasingly asymmetric (`verify_vertex_normal.py`).

```text
symmetric vertex : n_dup=[-0.210, -0.978]  n_uniq=[0.000, -1.000]  diff=12.095deg  third point missed by 35.0mm
asymmetric vertex: n_dup=[-0.057, -0.998]  n_uniq=[0.310, -0.951]  diff=21.319deg  third point missed by 76.1mm
more asymmetric  : n_dup=[0.161, -0.987]  n_uniq=[0.685, -0.728]  diff=33.983deg  third point missed by 143.4mm
```

The result runs opposite to expectation. The symmetric vertex looks like the case that should go furthest astray, yet it deviates least, and the discrepancy grows with asymmetry: 12° to 21° to 34°. The miss distance printed alongside explains why directly. The distance by which the corrupted circle misses the third point grows from 35 mm to 76 mm to 143 mm, and the angular error of the normal tracks that growth. A ray sent perpendicularly onto a symmetric vertex ought to retrace its own path (Fig. 2); computed with the normal from the reused third point, it leaves 12° off course instead — and at an asymmetric vertex the distortion is larger still.

<img src="/assets/img/posts/raytracing-reflection-family/en/fig2-vertex-normal-typo.png" alt="Effect of the vertex normal calculation" width="600">
_Fig 2. Reusing the third point (red solid) turns the reflected ray by 12.1°; using all three points (blue dashed) retraces the incident path exactly_

## 4. Limitations and what remains

The reflection formula itself is independent of the medium and one instance of it suffices, but the accuracy rests entirely on the logic that finds the normal. As `ArbitraryMirror` shows, hitting an edge and hitting a vertex require altogether different calculations for the normal, and the result is correspondingly sensitive to how that calculation is done.

Reflectance is not modelled either. Even `BeamSplitter` produces the reflected and transmitted rays as two full-strength rays with no division of intensity — a real beamsplitter would split the intensity according to its reflection-to-transmission ratio, 50:50 for instance.

Parts of the library remain untouched by this series. `Aperture`, `Stops` and `FlatSensor` neither reflect nor refract; they block rays or merely record intensity, which puts them in a different category from the reflective family. `CompoundSphericalLens`, the N-surface lens, is a generalization of `SphericalLens` and repeats principles already covered.

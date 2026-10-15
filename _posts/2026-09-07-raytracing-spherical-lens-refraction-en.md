---
title: "Ray Tracing 1 — How Do You Check That Your Own Lens Ray Calculation Is Right?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/raytracing-spherical-lens-refraction/
page_id: raytracing-spherical-lens-refraction
date: 2026-09-07 20:00:00 +0900
categories: [Optics, Geometrical Optics]
tags: [ray-tracing, geometric-optics, snells-law, python, matlab]
description: Rebuilding the spherical-lens refraction routine of a MATLAB 2D ray tracer in Python, from finding the ray-sphere intersection to checking the result against optical reversibility and the thick-lens formula.
math: true
---

Universities and research groups often have no ray-tracing or optical design software available. That was the case in the author's lab, where the need was not pressing either. The systems being designed were built around infinity-corrected objectives and assembled from catalogue components specified over the visible band actually in use. In an infinity-corrected system the space between objective and tube lens is collimated, so its length does not affect imaging, and filters or beam splitters can be inserted or removed without shifting focus. As the measurements demanded higher precision, however, and as observing the back focal plane (BFP) of high-magnification objectives became routine, the lab implemented its own ray tracer for designing jigs and pinning down focal distances.

The starting point was a MATLAB code, which the author ported to Python and cleaned up. With more than sixty functions and components ranging over lenses, mirrors, prisms, beam splitters and apertures, the question left after the port was how to confirm that the translated code actually worked. A plausible-looking ray diagram proves nothing: a refraction angle that is slightly wrong, or a wrong intersection point chosen near the rim of a lens, is invisible to the eye. This series works through the implementation and verification of one component at a time. The underlying theory is covered separately in the Geometrical Optics Foundations series.

The subject of this first post is the most elementary calculation of all: a ray refracting at a spherical surface. The principle is rebuilt from scratch and the result is then verified. The code shown here is a reduced reconstruction containing only what the verification needs (`make_spherical_lens`, `refract_through_lens` and so on; full code in `_code/raytracing-spherical-lens-refraction/`).

## 1. Where does the ray meet the lens surface?

Tracing a ray through a single lens begins with one question: which surface does the ray strike, and where? In the code a ray is four numbers, `[x, y, vx, vy]`. The pair $(x_0, y_0)=(x, y)$ is the position from which the ray sets out toward the surface currently under consideration, and $(v_x, v_y)$ is its direction vector. For the first surface of the lens (the one the ray meets first), $(x_0, y_0)$ is the source position; for the second surface it is the point where the ray has just been refracted at the first. Writing that starting point and direction as a general-form line simplifies everything that follows.

$$ ax + by = c, \qquad a = -v_y,\ \ b = v_x,\ \ c = a x_0 + b y_0 $$

Here $(a,b)=(-v_y, v_x)$ is perpendicular to the direction vector. The line is described by its normal rather than its slope, so vertical and horizontal rays alike fit this single form without special cases. One face of the lens is an arc of a circle with centre $(c_x, c_y)$ and radius of curvature $R$. Parameterizing points on that circle as $(R\cos w + c_x,\ R\sin w + c_y)$ and substituting into the line equation yields a trigonometric equation in the angle $w$. A line meets a circle in at most two points, so the equation has two solutions in general.

The difficulty is that only one of them is where the ray physically strikes the lens surface. The other is either the point met by extending the ray backwards, or a point outside the region the lens surface occupies, somewhere on the rest of the sphere. The two are therefore treated as candidates rather than answers, and the true intersection is the one that falls within the lens boundary. The code obtains the two candidates from two different trigonometric identities, one using `arccos` and one using `arcsin`. Figure 1(a) shows the case of a ray arriving at height $y=100$ mm. The `arcsin` candidate lies on the lens surface ($x=44.99$ mm); the `arccos` candidate lies on the far side of the same circle ($x=-1944.99$ mm). Only the first falls inside the extent of the lens surface (green box).

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig4-candidates-sag.png" alt="The two intersections of an extended ray with a spherical surface, and the sag measured from the vertex tangent plane" width="760">
_Fig 1. (a) The line through the $y=100$ mm ray meets the full circle of surface R2 (grey dashed) twice. Only the `arcsin` candidate (filled circle), which lies inside the extent of the lens surface (green box), is kept; the `arccos` candidate on the far side of the circle (open circle) is rejected. (b) Close-up of surface R2, stretched horizontally. Measured from the vertex tangent plane (black dashed), the sag is 31.75 mm at the rim and 5.01 mm at $y=100$ mm, where the ray lands._

```python
def find_entry_surface(ray, lens):
    """Find which of the lens's two surfaces the ray reaches first.

    The ray is treated as the line ax+by=c and each surface boundary as a
    chain of segments. Returns (surface index (1 or 2), distance to the
    intersection), or None if the ray misses the lens.
    """
    x, y, vx, vy = ray
    a_r, b_r = -vy, vx
    c_r = a_r * x + b_r * y

    best = None
    for surf_idx, boundary in enumerate(lens["BOUNDARY"], start=1):
        for i in range(len(boundary) - 1):
            x1, y1 = boundary[i]
            x2, y2 = boundary[i + 1]
            a_e, b_e = y2 - y1, x1 - x2
            c_e = a_e * x1 + b_e * y1

            M = np.array([[a_r, b_r], [a_e, b_e]])
            if abs(np.linalg.det(M)) < 1e-8:
                continue
            J = np.linalg.solve(M, np.array([c_r, c_e]))
            ...
```

At this stage the lens boundary is approximated by a chain of straight segments, a polygon standing in for the arc, purely to decide which surface is reached first. The refraction itself (`refract_through_lens`) solves the circle equation exactly, obtains the two candidate angles $w_c$ (from `arccos`) and $w_s$ (from `arcsin`), and keeps whichever falls within the region the lens surface actually occupies (the horizontal and vertical extent of its boundary points). If both do, it takes the one nearer the starting point.

```python
def _pick_candidate(boundary, offset, wc, xc, yc, ws, xs, ys, ref_x, ref_y):
    """Of the arccos/arcsin candidates, keep the one inside the lens boundary.

    If both lie inside, which is the usual case, keep the one nearer the
    reference point.
    """
    def inside(px, py):
        return (
            px <= np.max(boundary[:, 0]) + offset and px >= np.min(boundary[:, 0]) - offset
            and py <= np.max(boundary[:, 1]) + offset and py >= np.min(boundary[:, 1]) - offset
        )

    in_c, in_s = inside(xc, yc), inside(xs, ys)
    if in_c and in_s:
        if np.hypot(xc - ref_x, yc - ref_y) < np.hypot(xs - ref_x, ys - ref_y):
            return wc, xc, yc
        return ws, xs, ys
    elif in_c:
        return wc, xc, yc
    else:
        return ws, xs, ys
```

To see what `offset` is, one more property of the surface is needed. A sphere is not a plane. Measured against the plane tangent to the lens at its vertex, the surface falls away from that plane the further one moves from the optical axis. That departure is the sag, and at the rim of the lens ($D/2$) it equals $\mathrm{sag} = R - \sqrt{R^2 - (D/2)^2}$. For the lens in this post ($R=1000$ mm, $D=500$ mm) the sag at the rim is 31.75 mm (Figure 1(b)), so a lens 100 mm thick at the centre is only $100 - 2 \times 31.75 = 36.5$ mm thick at the edge. The `offset` is 1% of the sag, 0.32 mm. It keeps a point sitting exactly on the lens boundary from being judged to lie outside it because of floating-point error.

## 2. Which way does the ray bend at that point?

With the intersection found, the next question is how much the ray bends there. Snell's law answers it.

$$ n_1 \sin\theta_1 = n_2 \sin\theta_2 $$

Here $\theta_1$ and $\theta_2$ are the angles that the incident and refracted rays make with the surface normal at that point, and $n_1$, $n_2$ are the refractive indices on either side. Passing into the higher-index medium, the ray bends toward the normal. The code takes the normal direction from the circle's parameter angle $w$, forms $\theta_1$ from the angle relative to the incident direction vector, solves Snell's law for $\theta_2$, and converts the result back into a direction vector.

```python
theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
vx1, vy1 = np.cos(theta_t1 - (np.pi - w1)), np.sin(theta_t1 - (np.pi - w1))
```

The detailed derivation of these angle terms, which depends on the chosen parameterization of the circle, is not reworked here. Whether the expression really obeys Snell's law is checked numerically in section 3, against optical reversibility and against theory.

Sending a collimated beam through this lens gives the following paths.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig1-lens-focusing.png" alt="Parallel rays converging to a focus after a lens" width="600">
_Fig 2. Eleven parallel rays pass through a biconvex lens (R1=1000mm, R2=-1000mm, t=100mm, D=500mm, N-BK7, 750 nm) and converge near $x \approx -1000$ mm. Spherical aberration is invisible at this scale._

The lens sits at the origin, and the source is a bundle of eleven parallel rays 200 mm across, starting at $x=1000$ mm and travelling in the $-x$ direction. The rays appear to meet at a single point near $x \approx -1000$ mm behind the lens. Zoomed in, the outer rays cross the axis slightly earlier; section 3 measures by how much. A real lens behaves the same way. Figure 3 shows five parallel laser beams sent into a thin biconvex lens, which is exactly the situation a 2D ray tracer computes.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/ext-laser-biconvex-lens.jpg" alt="Five parallel laser beams entering a thin biconvex lens and crossing near one spot" width="700">
_Fig 3. Parallel beams from a five-line laser source (left) pass through a thin biconvex lens, cross near one spot and spread out again. The beams converging toward a point behind the lens is the same behaviour as in Figure 2. (Source: Карас Студио, [28. Сферна аберација со ласери](https://commons.wikimedia.org/wiki/File:28._%D0%A1%D1%84%D0%B5%D1%80%D0%BD%D0%B0_%D0%B0%D0%B1%D0%B5%D1%80%D0%B0%D1%86%D0%B8%D1%98%D0%B0_%D1%81%D0%BE_%D0%BB%D0%B0%D1%81%D0%B5%D1%80%D0%B8.ogv), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0), frame from the video, cropped)_

One point deserves care here. Since the rays travel from $+x$ toward $-x$, it is tempting to assume they meet surface R1, the one defined first, but `find_entry_surface` reports surface 2 (R2) as the entry surface. R1 and R2 are merely the order of the parameters in the lens definition; they do not correspond to physical left and right as seen from the source. This lens has R1=1000 (convex) and R2=-1000 (convex), and computing the vertex coordinates of the two surfaces directly places the R1 surface at $x=-t/2$ and the R2 surface at $x=+t/2$. A ray arriving from the $+x$ side therefore meets R2 first. This sign convention returns in the reversibility check of section 3.

## 3. Is the calculation right?

In short, both checks pass. A reversed ray returns along its own path to floating-point round-off (1e-16), and the paraxial focus agrees with the thick-lens formula to within 0.062%. The rest of this section shows how each was established.

### Does a reversed ray retrace its path?

Optical reversibility is the principle that a ray can retrace its own path exactly in the opposite direction. Snell's law keeps its form when time is run backwards, so a refracted ray fired back along its reversed direction must return along the path it came in on. This makes a useful cross-check on the two branches of `refract_through_lens` (the `entry_surface==1` case and its counterpart), because forward entry exercises one branch and the retrace exercises the other. If the two branches implement the same physics without contradiction, the retrace must close.

```python
ray_fwd = np.array([1000.0, 100.0, -1.0, 0.0])
entry_fwd = find_entry_surface(ray_fwd, lens)[0]
path_fwd = refract_through_lens(ray_fwd, entry_fwd, lens, n_lens, N_AIR)

# advance the exit ray 500 mm beyond the lens, then reverse its direction only
x2, y2, vx2, vy2 = path_fwd[-1]
far_point = np.array([x2 + 500 * vx2, y2 + 500 * vy2])
ray_retrace = np.array([far_point[0], far_point[1], -vx2, -vy2])
entry_retrace = find_entry_surface(ray_retrace, lens)[0]
path_retrace = refract_through_lens(ray_retrace, entry_retrace, lens, n_lens, N_AIR)
```

Running the forward trace and the retrace side by side gives:

```text
forward entry surface=2, retrace entry surface=1 (must be the opposite face)
original incident ray start: (1000.000, 100.000), direction (-1.000, 0.000)
retrace arrival (first refraction point): (44.987, 100.000), direction (1.000, -0.000)
position error: 0.00e+00 mm
```

The entry surface comes back exactly reversed (2 → 1). The retraced ray returns to the original first refraction point (44.987, 100 mm) with zero position error, and recovers its direction, reversed in sign, to within floating-point round-off (1e-16). In the upper panel of Figure 4 the two paths cannot be told apart. The lower panel shows the height difference between them at the same $x$: zero on the exit side, and no more than $3 \times 10^{-13}$ mm on the incident side even at $x=1000$ mm.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig2-reversibility.png" alt="Reversibility check: forward and retraced paths, and the height difference between them" width="600">
_Fig 4. Top: the retrace of the reversed exit ray (red dashed) overlaps the forward path (blue solid) all the way from the lens to x=1000 mm (vertical scale stretched). Bottom: height difference between the two paths at the same $x$ (log scale). Where the difference is exactly zero it is drawn at the floor ($10^{-17}$). The largest difference, $2.8 \times 10^{-13}$ mm, occurs at $x=1000$ mm._

### Does the focus match the thick-lens formula?

The reversibility check shows that the two branches of `refract_through_lens` behave consistently with each other, but it says nothing about whether the result matches the physical lens: if both branches were wrong in the same way, reversibility would still hold. The result is therefore compared against a theoretical expression derived independently of this code. The verdict, stated up front: paraxial rays agree almost exactly, and only the outer rays depart. The thick-lens formulae give the focal length and the back focal distance (BFD) with the lens thickness taken into account.

$$ \frac{1}{f} = (n-1)\left[\frac{1}{R_1} - \frac{1}{R_2} + \frac{(n-1)t}{n R_1 R_2}\right], \qquad \mathrm{BFD} = f\left[1 - \frac{(n-1)t}{n R_1}\right] $$

Here $n$ is the index of the glass relative to air and $t$ is the centre thickness. The focal length $f$ is measured from the rear principal plane, while the BFD is measured from the rear vertex of the lens to the focus. As $t \to 0$ the third term vanishes and the thin-lens formula remains. Treated as a thin lens, this lens would have $f$ = 976.9 mm, 16.8 mm (1.7%) shorter than the 993.7 mm obtained with the 100 mm thickness included. The BFD (960.1 mm) is in turn 33.6 mm shorter than $f$; that difference is how far the rear principal plane sits inside the lens from the rear vertex. In practice, a ruler laid from the back of the lens to the focus measures the BFD, not $f$.

In the BFD expression, $R_1$ is the radius of the surface the light reaches first. As established in section 2, in this arrangement that surface is R2, but the lens is symmetric, $\lvert R_1 \rvert = \lvert R_2 \rvert$, so the BFD is the same from either side.

These formulae rest on the paraxial approximation, valid for rays very close to the optical axis. The ray trace computes the true geometry with no such approximation, so there is no reason for a ray well off axis, least of all the marginal ray at the edge of the beam, to agree with them exactly. Exact agreement would in fact be the suspicious outcome. The discrepancy is spherical aberration.

```python
n = n_lens  # N-BK7, 750nm
inv_f = (n - 1) * (1 / R1 - 1 / R2 + (n - 1) * T / (n * R1 * R2))
f = 1 / inv_f
bfd = f * (1 - (n - 1) * T / (n * R1))
expected_focus_x = -T / 2 - bfd
```

For this lens, theory and simulation separate as follows.

```text
n(N-BK7, 750nm) = 1.511835
thick-lens theory (paraxial): f=993.698mm, BFD=960.056mm, focus x=-1010.056
simulated paraxial focus (|y0|<=20mm) x=-1009.431, vs theory +0.625mm (+0.062%)
simulated marginal focus (|y0|=100mm) x=-994.211, vs paraxial +15.220mm (spherical aberration)
```

Rays entering within 20 mm of the axis, 8% of the 250 mm lens radius, agree with theory to 0.062%, and even that small offset lies toward the lens, the direction of spherical aberration. Rays entering at the edge of the beam, 100 mm off axis, focus closer to the lens than the paraxial focus. The 15.220 mm in the output is measured from the focus of the $y_0=20$ mm ray; measured from the paraxial theoretical focus it is 15.845 mm. That is the direction of undercorrected spherical aberration in a simple biconvex lens, where the focus moves toward the lens as the ray moves outward (Hecht, section 6.3).

How does the shift depend on ray height? Tracing rays at heights $y_0$ in 5 mm steps gives Figure 5.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig3-focus-vs-theory.png" alt="Focus shift versus ray height" width="600">
_Fig 5. Shift of the focus toward the lens, measured from the paraxial theoretical focus, as a function of ray height $$y_0$$ (red dots: ray trace). The blue line is $$k y_0^2$$ with the coefficient $k$ taken from the $$y_0=1$$ mm ray. The points follow the parabola; even at $$y_0=100$$ mm they depart from it by only 1.4%._

The points follow a parabola: the shift grows as the square of the ray height. Dividing the shift by $y_0^2$ gives about $1.563 \times 10^{-3}$ mm$^{-1}$ over $y_0$ = 5–100 mm, with a spread under 1.4%. The 0.625 mm at $y_0=20$ mm thus becomes about 25 times larger, 15.8 mm, at $y_0=100$ mm. The 0.625 mm offset of the paraxial focus toward the lens, noted earlier, is the same law evaluated at $$y_0=20$$ mm. As a working rule for this lens over $$y_0\le100$$ mm, halving the beam diameter cuts the spherical aberration to a quarter.

What this aberration looks like in a real lens is shown in Figure 6: a laser beam passing through a plano-convex lens in a scattering medium.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/ext-lens-spherical-aberration.jpg" alt="A laser beam passing a plano-convex lens and converging into a caustic rather than a point" width="600">
_Fig 6. A laser beam entering from the left leaves a plano-convex lens through its convex face (the rounded part in the middle) and becomes visible in a scattering medium. It does not converge to a point: light near the edge comes together first, close to the lens, tracing a sharp caustic. This is the tendency computed in Figure 5, that rays farther from the axis focus closer to the lens, seen with the eye in a different lens. (Source: Adam Gabrys, [Lens spherical aberration RL](https://commons.wikimedia.org/wiki/File:Lens_spherical_aberration_RL.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

## 4. Summary and limitations

This post rebuilt in Python the calculation of a ray refracting at a spherical surface: finding the entry point as a line-circle intersection and obtaining the refracted direction there by Snell's law (`refract_through_lens`). The result was checked in two ways. The reversibility check reproduced the incident path to floating-point round-off. Against the thick-lens formulae, the paraxial focus agreed within 0.062%, and the outer rays departed only by a spherical aberration that grows as the square of the ray height, in the expected direction. The reversibility check passes even if both branches are wrong in the same way, so the independent check is the thick-lens comparison. Both methods reaching the same conclusion is good reason to trust the spherical refraction in `refract_through_lens` under the conditions of this post: a symmetric biconvex lens with R=±1000 mm, one wavelength, and a single $$y_0=100$$ mm ray for reversibility. Asymmetric and meniscus lenses were not checked.

What this function does not do is test for total internal reflection (TIR). When $\sin\theta_2 = (n_1/n_2)\sin\theta_1$ exceeds 1, `arcsin` leaves its domain. The function does not intercept that case and passes the `NaN` returned by NumPy straight into the next calculation. The original MATLAB code behaved the same way in its spherical-lens routine; total internal reflection was tested explicitly only in the prism routines. At the first surface, where light passes from air into the denser glass, total internal reflection cannot occur at all. At the second surface, glass to air, an ordinary lens keeps the internal angle well below the critical angle (about 41.4° for N-BK7 at 750 nm). For the lens in this post the internal angle at the second surface is 7.5° for the $y_0=100$ mm ray and only 19.2° for a ray just inside the rim ($y_0=249$ mm). Even with light across the full 500 mm diameter, it stays below half the critical angle. A combination of strong curvature and a wide aperture, however, can push the internal angle past the critical angle. [Post 2](/en/posts/raytracing-total-internal-reflection/) covers the total-internal-reflection test used for prisms, and what happens when it is applied to a lens instead.

## References

- E. Hecht, *Optics*, 5th ed., Pearson, 2017. Section 4.7 (total internal reflection), section 6.1 (thick lenses and back focal distance), section 6.3 (spherical aberration).
- SCHOTT, "TIE-29: Refractive Index and Dispersion," Technical Information. Basis for treating catalogue indices as relative to air. [PDF](https://media.schott.com/api/public/content/aaa572afd854434fb7b3faa4bc46103f)
- N-BK7 refractive index data (SCHOTT catalogue Sellmeier coefficients), refractiveindex.info. [N-BK7](https://refractiveindex.info/?shelf=specs&book=SCHOTT-optical&page=N-BK7)

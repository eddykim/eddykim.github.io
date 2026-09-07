---
title: "Geometrical Optics 1 — Refraction at a Spherical Surface, and How to Verify It"
lang: en
lang-exclusive: ["en"]
permalink: /posts/raytracing-spherical-lens-refraction/
page_id: raytracing-spherical-lens-refraction
date: 2026-09-07 20:00:00 +0900
categories: [Optics, Geometrical Optics]
tags: [ray-tracing, geometric-optics, snells-law, python, matlab]
description: Rebuilding the spherical-lens refraction routine of a MATLAB 2D ray tracer in Python — from finding the ray-sphere intersection to verifying the result against optical reversibility and the thick-lens formula.
math: true
---

Universities and research groups often have no ray-tracing or optical design software available. That was the case in the lab this work comes from, where the need was not pressing either. The systems being designed were built around infinity-corrected objectives and assembled from catalogue components specified over the visible band actually in use. In an infinity-corrected system the space between objective and tube lens is collimated, so the length of that space does not affect imaging; filters and beamsplitters can be inserted or removed without shifting focus. As the measurements demanded higher precision, however, and as observing the back focal plane (BFP) of high-magnification objectives became routine, the lab implemented its own ray tracer for designing jigs and pinning down focal distances.

The starting point was an open-source MATLAB implementation, which was ported to Python and cleaned up. With more than sixty functions and components ranging over lenses, mirrors, prisms, beamsplitters and apertures, the question left after the port was how to confirm that the translated code actually worked. A plausible-looking ray diagram proves nothing: a refraction angle that is slightly wrong, or a wrong intersection point chosen near the rim of a lens, is invisible to the eye. This series introduces the geometrical optics involved and works through the implementation and verification of one component at a time.

The subject of this first post is the most elementary calculation of all — a ray refracting at a spherical surface. The principle is rebuilt from scratch and the result is then verified. The code shown here is a reduced reconstruction containing only what the verification needs (`make_spherical_lens`, `refract_through_lens` and so on; full code in `_code/raytracing-spherical-lens-refraction/`).

## 1. Ray-sphere intersection — two candidate solutions

Tracing a ray through a single lens begins with one question: which surface does the ray strike, and where? In the code a ray is four numbers, `[x, y, vx, vy]`. The pair $(x_0, y_0)=(x, y)$ is the position from which the ray sets out toward the surface currently under consideration, and $(v_x, v_y)$ is its direction vector. For the first surface of the lens — the one the ray meets first — $(x_0, y_0)$ is the source position, but for the second surface it is the point where the ray has just been refracted out of the first surface. In other words it is not the source itself, but wherever the ray departs from on its way to this particular surface. Writing that starting point and direction as a general-form line simplifies everything that follows.

$$ ax + by = c, \qquad a = -v_y,\ \ b = v_x,\ \ c = a x_0 + b y_0 $$

Taking the coefficients as $(a,b)=(-v_y, v_x)$, perpendicular to the direction vector, lets vertical and horizontal rays alike be expressed in this single form without special cases. One face of the lens is an arc of a circle with centre $(c_x, c_y)$ and radius of curvature $R$, so parameterizing points on that circle as $(R\cos w + c_x,\ R\sin w + c_y)$ and substituting into the line equation yields a trigonometric equation in the angle $w$. A line meets a circle in at most two points, and the equation accordingly has two solutions in general.

The difficulty is that only one of those two solutions is where the ray physically strikes the lens surface. The other is either the point met by extending the ray backwards, or a point lying outside the physical diameter $D$ of the lens — on the extension of the sphere, where no glass exists. The two are therefore treated as candidates rather than answers, and the true intersection is selected afterwards by testing whether the point falls within the lens boundary. The code obtains the two candidates from two different trigonometric identities, one using `arccos` and one using `arcsin`.

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

At this stage the lens boundary is approximated by a chain of straight segments — a polygonal approximation of the arc — purely to decide which surface is reached first. The refraction itself (`refract_through_lens`) solves the circle equation exactly, obtains the two candidate angles $w_c$ (from `arccos`) and $w_s$ (from `arcsin`), and keeps whichever falls inside the physical diameter $D$.

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

A sphere is not a plane. Measured against the plane tangent to the lens at its vertex, the surface departs further from that plane the further one moves from the optical axis. That departure is the sag, and at the rim of the lens ($D/2$) it equals $\mathrm{sag} = R - \sqrt{R^2 - (D/2)^2}$. The `offset` is this sag scaled by a small factor (0.01). It prevents a point sitting exactly on the lens boundary from being misjudged as lying outside it because of floating-point error.

Once the intersection is known, Snell's law is applied there.

$$ n_1 \sin\theta_1 = n_2 \sin\theta_2 $$

Here $\theta_1$ and $\theta_2$ are the angles that the incident and refracted rays make with the surface normal at that point. The code takes the normal direction from the circle's parameter angle $w$, forms $\theta_1$ from the angle relative to the incident direction vector, solves Snell's law for $\theta_2$, and converts the result back into a direction vector.

```python
theta_t1 = np.arcsin(n_air * np.sin(np.pi - w1 + np.arctan2(vy, vx)) / n_lens)
vx1, vy1 = np.cos(theta_t1 - (np.pi - w1)), np.sin(theta_t1 - (np.pi - w1))
```

The detailed derivation of these angle terms — the part that depends on which parameterization of the circle was chosen — is not reworked here. What matters is whether the expression actually produces results consistent with Snell's law, and that is checked numerically in section 2 against optical reversibility and against theory.

Sending a collimated beam through this lens gives the following paths.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig1-lens-focusing.png" alt="Parallel rays converging to a focus after a lens" width="600">
_Fig 1. Eleven parallel rays pass through a biconvex lens (R1=1000mm, R2=-1000mm, t=100mm, D=500mm, N-BK7) and converge to a point_

The lens sits at the origin, and the source is a bundle of eleven parallel rays 200 mm across, starting at $x=1000$ mm and travelling in the $-x$ direction. One point deserves care here. Since the rays travel from $+x$ toward $-x$, it is tempting to assume they meet surface R1 — the one defined first — but `find_entry_surface` reports surface 2 (R2) as the entry surface. R1 and R2 are merely the order of the parameters in the lens definition; they do not correspond to physical left and right as seen from the source. This lens has R1=1000 (convex) and R2=-1000 (convex), and computing the vertex coordinates of the two surfaces directly places the R1 surface at $x=-t/2$ and the R2 surface at $x=+t/2$. A ray arriving from the $+x$ side therefore meets R2 first. This sign convention returns in the reversibility check of section 2.

## 2. Verification — reversibility and comparison against theory

### Optical reversibility

Optical reversibility is the principle that a ray can retrace its own path exactly in the opposite direction. It holds because Snell's law is symmetric under time reversal: reverse the direction of a refracted ray and fire it back, and it must retrace precisely the path it came in on. This makes a useful cross-check on the two branches of `refract_through_lens` (the `entry_surface==1` case and its counterpart), because forward entry exercises one branch and the retrace exercises the other. If the two branches implement the same physics without contradiction, the retrace must close.

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

The entry surface comes back exactly reversed (2 → 1), and the retraced ray recovers the starting $y$ coordinate of the original incident ray (100 mm) and its direction — same magnitude, opposite sign — with no floating-point error at all.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig2-reversibility.png" alt="Reversibility check" width="600">
_Fig 2. Firing the exit ray backwards retraces the incident path exactly; the blue solid line and the red dashed line coincide completely between the lens and x=1000_

### Comparison with the thick-lens formula — near-exact on axis, offset at the rim

The reversibility check shows that the two branches of `refract_through_lens` behave consistently with each other, but it says nothing about whether the result matches the physical lens: if both branches were wrong in the same way, reversibility would still hold. The result is therefore compared against an independent theoretical expression, derived with no reference whatsoever to this ray-tracing code. The thick-lens formulae give the focal length and the back focal distance (BFD) with the lens thickness taken into account.

$$ \frac{1}{f} = (n-1)\left[\frac{1}{R_1} - \frac{1}{R_2} + \frac{(n-1)t}{n R_1 R_2}\right], \qquad \mathrm{BFD} = f\left[1 - \frac{(n-1)t}{n R_1}\right] $$

In the BFD expression, $R_1$ is the radius of the surface the light reaches first. As established in section 1, in this arrangement that surface is R2 — but the lens is symmetric, $\lvert R_1 \rvert = \lvert R_2 \rvert$, so the BFD is the same from either side.

These formulae rest on the paraxial approximation, valid for rays very close to the optical axis. The ray trace computes the true geometry with no such approximation, so there is no reason for a marginal ray — one any distance off axis — to agree with them exactly. Exact agreement would in fact be the suspicious outcome. That discrepancy is spherical aberration.

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
simulated paraxial focus (|y0|<=20mm) x=-1010.229, vs theory -0.172mm (-0.017%)
simulated marginal focus (|y0|=100mm) x=-995.010, vs paraxial +15.219mm (spherical aberration)
```

Rays entering within 20 mm of the axis — 8% of the 250 mm lens radius — agree with theory to 0.017%, which is exact for practical purposes. Rays entering at the edge of the aperture, 100 mm off axis, focus 15.2 mm closer to the lens than the paraxial focus. That is the direction expected of undercorrected spherical aberration in a simple biconvex lens, where the focal distance shortens as the ray moves outward.

<img src="/assets/img/posts/raytracing-spherical-lens-refraction/en/fig3-focus-vs-theory.png" alt="Comparison of focus positions" width="600">
_Fig 3. The paraxial theoretical focus (black dashed) and the simulated paraxial focus (blue dotted) nearly coincide; only the marginal-ray focus (red dotted) is displaced, by the amount of the spherical aberration_

## 3. Summary and limitations

This post rebuilt in Python the calculation of a ray refracting at a spherical surface — finding the entry point by line-circle intersection and obtaining the refracted direction there by Snell's law (`refract_through_lens`) — and checked the result two independent ways, by optical reversibility and against the thick-lens formulae. The reversibility check reproduced the incident path with no floating-point error, and the comparison with theory placed the paraxial focus within 0.017% while the marginal ray departed only by the spherical aberration expected, in the expected direction. Two independent methods reaching the same conclusion is good reason to trust the spherical refraction in `refract_through_lens`.

What this function does not do is test for total internal reflection. When $\sin\theta_2 = (n_1/n_2)\sin\theta_1$ exceeds 1, `arcsin` leaves its domain; the function does not intercept that case and passes the `NaN` returned by NumPy straight into the next calculation. The original MATLAB behaved the same way in its spherical-lens routine — total internal reflection was tested explicitly only in the prism routines. Rays entering a lens from air, into a denser medium, meet the surface at shallow enough angles that the condition is rarely reached, but a combination of strong curvature and a wide aperture can reach it. The next post covers the total-internal-reflection logic used for prisms, and what happens when it is applied to a lens instead.

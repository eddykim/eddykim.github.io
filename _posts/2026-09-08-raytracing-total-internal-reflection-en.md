---
title: "Geometrical Optics 2 — Detecting Total Internal Reflection, and Applying It to a Lens"
lang: en
lang-exclusive: ["en"]
permalink: /posts/raytracing-total-internal-reflection/
page_id: raytracing-total-internal-reflection
date: 2026-09-08 20:00:00 +0900
categories: [Optics, Geometrical Optics]
tags: [ray-tracing, geometric-optics, snells-law, python, matlab]
description: Verifying the total internal reflection test in the prism routine, then applying the same logic to the exit surface of the spherical lens refraction function — and counting how often the reflection formula recurs across the library.
math: true
---

[Post 1](/en/posts/raytracing-spherical-lens-refraction/) documented the implementation of the lens refraction function and closed with a stated limitation: the function does not test for total internal reflection (TIR). When $\sin\theta_2 = (n_1/n_2)\sin\theta_1$ exceeds 1, `arcsin` leaves its domain, and the lens routine passes the resulting `NaN` straight through instead of intercepting it. The prism routine (`CalculateRayPath_Prism`), by contrast, was written to handle that case from the outset. This post verifies the prism's TIR test first, then transplants the same logic into the exit surface of the lens function.

## 1. TIR in the prism — the law of reflection on a polygonal boundary

The lens and prism routines differ from the very first step, in how they find the ray-surface intersection. A lens surface is a circular arc, so the intersection angle is solved directly with `arccos` and `arcsin` (see [post 1](/en/posts/raytracing-spherical-lens-refraction/)); a prism boundary is a chain of segments joining vertices, so line-line intersection and vector arithmetic suffice. The lens routine also terminates after one entry and one exit, whereas the prism routine loops with `while` until the ray no longer meets any boundary — meaning any number of internal reflections can occur inside it.

```python
def trace_ray_through_prism(ray, prism, n_lens, n_air, max_bounces=20, verbose=False):
    """Trace the ray until it no longer meets any prism boundary."""
    xs, ys, vx, vy = map(float, ray)
    boundary = prism["BOUNDARY"]
    air2glass = True
    path = [[xs, ys, vx, vy]]

    for _ in range(max_bounces):
        # ... (finding the nearest boundary intersection (xc, yc), omitted) ...

        # surface normal -- perpendicular to the segment, signed to face the incoming ray
        nx, ny = -(y2 - y1), (x2 - x1)
        nnorm = np.hypot(nx, ny)
        nx, ny = nx / nnorm, ny / nnorm
        if vx * nx + vy * ny < 0:
            nx, ny = -nx, -ny

        theta_in = np.arccos(np.clip((vx * nx + vy * ny) / np.hypot(vx, vy), -1, 1))
        n1, n2 = (n_air, n_lens) if air2glass else (n_lens, n_air)
        sin_out = n1 * np.sin(theta_in) / n2

        if abs(sin_out) <= 1:
            theta_out = np.arcsin(sin_out)
            # ... refracted direction from Snell's law ...
            air2glass = not air2glass
        else:
            vout = np.array([vx, vy]) - 2 * (vx * nx + vy * ny) * np.array([nx, ny])

        vx, vy = vout / np.linalg.norm(vout)
        xs, ys = xc + 0.5 * vx, yc + 0.5 * vy  # avoid re-hitting the same face
        path.append([xc, yc, vx, vy])
    return np.array(path)
```

The domain is checked first with `abs(sin_out) <= 1`; if it is exceeded, the ray is reflected about the surface normal $\hat n$ instead of refracted.

$$ \vec v' = \vec v - 2(\vec v \cdot \hat n)\hat n $$

This code was given a ray at normal incidence on an equilateral N-BK7 prism of 600 mm side.

<img src="/assets/img/posts/raytracing-total-internal-reflection/en/fig1-prism-tir-path.png" alt="Total internal reflection path inside an equilateral prism" width="600">
_Fig 1. A normally incident ray undergoes one total internal reflection inside an equilateral prism_

A ray entering along the normal looks as though it should pass straight through, but running the code shows it reflecting off the inner face first (internal angle 60°) and only then leaving through another face without deviation (internal angle 0° — it happens to strike that face perpendicularly).

```text
face=(0.0,-300.0)-(0.0,300.0)  theta_in=0.000deg  transmitted  theta_out=0.000deg
face=(0.0,300.0)-(519.6,0.0)  theta_in=60.000deg  total internal reflection (TIR)
face=(519.6,0.0)-(0.0,-300.0)  theta_in=0.000deg  transmitted  theta_out=0.000deg
```

The reason is straightforward. The critical angle of N-BK7 is $\theta_c = \arcsin(1/n) = 41.41°$ (with $n=1.511835$ at 750 nm), while the apex angle of an equilateral triangle is 60°. A ray entering normal to the first face continues undeviated, so it strikes the next face at exactly $60° - 0° = 60°$ internally. Since $60° > 41.41°$, total internal reflection is unavoidable. Transmitting the prism, in other words, requires a sufficiently *large* angle of incidence. With $\theta_1$ the angle of incidence, $\theta_1'$ the refraction angle at the first face, $\theta_2'$ the internal angle at the second face and $A$ the apex angle, the relation $\theta_1' + \theta_2' = A$ holds, so avoiding the critical angle at the second face requires

$$ \theta_1' > A - \theta_c, \qquad \theta_1 > \arcsin\big(n \sin(A - \theta_c)\big) $$

Substituting the values gives $\theta_1 > 28.8131°$. Sweeping the angle of incidence from 0° by bisection to locate the TIR-to-transmission boundary the code actually observes (`verify_tir.py`) returns 28.813100°, agreeing with theory to within 7.6e-10°.

<img src="/assets/img/posts/raytracing-total-internal-reflection/en/fig2-prism-angle-sweep.png" alt="Total internal reflection and transmission across an incidence-angle sweep" width="600">
_Fig 2. Once the angle of incidence exceeds the threshold (28.8°), the ray transmits instead of reflecting_

This is not an exotic effect. The Porro prism — the one in binoculars, which folds the light path by total internal reflection alone, with no mirror coating — works on exactly this principle. Light entering normal to the entrance face strikes the two faces inclined at 45° with an internal angle of 45° each, and since 45° > 41.4° (the critical angle of BK7) it must reflect totally, so those faces act as mirrors without any coating.

## 2. The same logic in the lens — handling TIR at the exit surface

The exit-surface refraction in the lens routine reads as follows (carried over unchanged from [post 1](/en/posts/raytracing-spherical-lens-refraction/)).

```python
phi2 = np.arctan2(vy1, vx1)
theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)  # no domain check
vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
```

When `n_lens * sin(w2 - phi2) / n_air` exceeds 1 — that is, when the internal angle exceeds the critical angle — this `arcsin` has no real solution at all. The original code does not intercept that case, so rays past the critical angle quietly produce `NaN`.

This function does not decide on its own which surface the ray enters through; the caller must supply `entry_surface`. Since R1 and R2 are merely the order of the parameters in the lens definition and bear no relation to which surface faces the source (see [post 1](/en/posts/raytracing-spherical-lens-refraction/)), the value passed must always come from `find_entry_surface` first. Otherwise `NaN` can appear somewhere unrelated to total internal reflection — in the exit-surface intersection itself.

Checking this with the y0 = 36 mm ray of a plano-convex lens (R1 set to $10^6$ mm, effectively flat; R2 = -50 mm, T = 20 mm, D = 78 mm):

```text
find_entry_surface: surface 1

before the fix:
[[-9.99935200e+00  3.60000000e+01  1.00000000e+00 -1.21878769e-05]
 [-5.30123745e+00  3.59999427e+01             nan             nan]]
-> RuntimeWarning: invalid value encountered in arcsin
```

The fix applies the same pattern as the prism: outside the domain, reflect about the exit-surface normal $(\cos w_2, \sin w_2)$ rather than refract. At a point on a sphere the normal is the radial direction from the centre through that point, so the circle's parameter angle $w_2$ can be used directly.

```python
sin_t2 = n_lens * np.sin(w2 - phi2) / n_air
if abs(sin_t2) <= 1:
    theta_t2 = np.arcsin(sin_t2)
    vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
else:
    n2x, n2y = np.cos(w2), np.sin(w2)
    dot2 = vx1 * n2x + vy1 * n2y
    vx2, vy2 = vx1 - 2 * dot2 * n2x, vy1 - 2 * dot2 * n2y
```

```text
after the fix:
[[-9.99935200e+00  3.60000000e+01  1.00000000e+00 -1.21878769e-05]
 [-5.30123745e+00  3.59999427e+01  3.68088814e-02 -9.99322324e-01]]
```

Whether this is genuinely a reflection was verified two ways (`verify_tir.py`). First, that the law of reflection holds about the normal, angle of incidence against angle of reflection: 46.05508420° versus 46.05508420°, differing by 7.1e-15°. Second, the onset of TIR was located by bisecting over y0, and the internal angle at that boundary was compared against the theoretical critical angle. The boundary the code finds is y0 = 33.072057 mm, where the internal angle is 41.410389° against a theoretical $\theta_c = 41.410389°$ — agreement to within 2.7e-10°.

<img src="/assets/img/posts/raytracing-total-internal-reflection/en/fig3-lens-marginal-tir.png" alt="Total internal reflection of the outermost rays in a plano-convex lens" width="600">
_Fig 3. Of 39 rays entering a plano-convex lens (R2=-50mm), the 6 outermost undergo total internal reflection at the exit surface_

Over the half aperture of 39 mm (D = 78 mm), TIR sets in only from 33 mm outward — roughly the outer 15% of the aperture. For a gently curved lens that band either does not exist or lies well outside the aperture, so it is rarely met in practice; but with curvature as strong as this one, where R2 approaches D/2 and the surface is nearly hemispherical, it cannot be ignored.

## 3. Limitations, and a reflection formula that keeps recurring

Unlike the prism routine, the corrected lens function still models exactly one entry and one exit. The direction of a ray totally reflected at the exit surface is now computed correctly, but where that ray goes next inside the lens is not implemented: there is no loop running until the ray stops meeting boundaries, as there is for the prism. Tracing the reflected ray onward would mean rewriting the lens function along the lines of the prism one.

One last thing was worth checking. The formula used twice in this post for reflection — once for the prism's internal reflection, once for the lens TIR — $\vec v' = \vec v - 2(\vec v\cdot\hat n)\hat n$, appears in exactly this form six times across the library: in `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter` (the reflected component), and now in the TIR branches of `Prism` and `SphericalLens`. All that differs between them is how the intersection and the normal are found — plane algebra, candidate roots on a circle, or local curvature along a polygon — and whether the routine reflects once and stops or loops until the ray escapes. Refraction demanded a different angle expression every time; reflection, whether at a mirror, a prism or a lens, reduces to this one vector formula. Post 3 checks whether the formula really is implemented consistently across the code, verifying `FlatMirror`, `SphericalMirror`, `ArbitraryMirror` and `BeamSplitter` together.

---
title: "Ray Tracing 2 — Why Can Light Entering the Edge of a Lens Fail to Get Out?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/raytracing-total-internal-reflection/
page_id: raytracing-total-internal-reflection
date: 2026-09-08 20:00:00 +0900
categories: [Optics, Geometrical Optics]
tags: [ray-tracing, geometric-optics, snells-law, python, matlab]
description: Verifying the total internal reflection test in the prism routine, adding the same test to the exit surface of the spherical-lens routine, and following where the totally reflected light goes.
math: true
---

[Post 1](/en/posts/raytracing-spherical-lens-refraction/) documented the lens refraction function and closed with a stated limitation: the function does not test for total internal reflection (TIR). When $\sin\theta_2 = (n_1/n_2)\sin\theta_1$ exceeds 1, `arcsin` leaves its domain, and the lens routine passes the resulting `NaN` on to the next step instead of intercepting it. The prism routine (`CalculateRayPath_Prism`), by contrast, was written to handle that case from the outset. This post verifies the prism's TIR test first, then adds the same test to the exit surface of the lens function.

## 1. How does the prism routine detect total internal reflection?

How does the prism routine notice that `arcsin` is about to leave its domain? Before computing a refraction angle it checks whether the sine exceeds 1, and if it does, it reflects the ray off that face instead. In effect it applies the law of reflection on each segment of a polygonal boundary.

The lens and prism routines already differ in how they find the ray-surface intersection. A lens surface is a circular arc, so the intersection angle is solved directly with `arccos` and `arcsin` (see [post 1](/en/posts/raytracing-spherical-lens-refraction/)). A prism boundary is a chain of segments joining vertices, so line-line intersection and vector arithmetic suffice. The lens routine also stops after one entry and one exit, whereas the prism routine keeps looping until the ray meets no further boundary, so any number of internal reflections can happen inside it. The excerpt below writes that loop as a `for` loop capped at `max_bounces=20`.

```python
def trace_ray_through_prism(ray, prism, n_lens, n_air, max_bounces=20, verbose=False):
    """Trace the ray until it no longer meets any prism boundary."""
    xs, ys, vx, vy = map(float, ray)
    boundary = prism["BOUNDARY"]
    air2glass = True
    path = [[xs, ys, vx, vy]]

    for _ in range(max_bounces):
        # ... (finding the nearest boundary intersection (xc, yc), omitted) ...

        # surface normal -- perpendicular to the segment, signed to point along the ray (v·n >= 0)
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

The domain is checked first with `abs(sin_out) <= 1`; outside it, the ray is reflected about the surface normal $\hat n$ instead of refracted.

$$ \vec v' = \vec v - 2(\vec v \cdot \hat n)\hat n $$

Here $\vec v$ is the incident direction, $\vec v'$ the reflected direction and $\hat n$ the unit normal of the face. The formula flips the component along the normal and leaves the component along the face untouched. Replacing $\hat n$ with $-\hat n$ changes nothing, because the normal appears twice, so for reflection it does not matter which way the normal points ([post 3](/en/posts/raytracing-reflection-family/) returns to this).

The code was run on an equilateral N-BK7 prism with 600 mm sides and a ray at normal incidence.

<img src="/assets/img/posts/raytracing-total-internal-reflection/en/fig1-prism-tir-path.png" alt="Total internal reflection path inside an equilateral prism" width="600">
_Fig 1. A ray entering one face of an equilateral N-BK7 prism (600 mm sides, 750 nm) at normal incidence. It meets the second face at an internal angle of 60° and is totally reflected, then strikes the third face perpendicularly and leaves undeviated._

A ray entering along the normal looks as though it should pass straight through. Running the code shows it reflecting totally off the second face first, then striking the third face perpendicularly and leaving without deviation. The internal angle here is the angle of incidence measured inside the glass from the face normal (`theta_in` in the output below).

```text
face=(0.0,-300.0)-(0.0,300.0)  theta_in=0.000deg  transmitted  theta_out=0.000deg
face=(0.0,300.0)-(519.6,0.0)  theta_in=60.000deg  total internal reflection (TIR)
face=(519.6,0.0)-(0.0,-300.0)  theta_in=0.000deg  transmitted  theta_out=0.000deg
```

The reason is simple. The critical angle of N-BK7 is $\theta_c = \arcsin(1/n) = 41.41°$ (with $n=1.511835$ at 750 nm), while the apex angle of an equilateral prism is 60°. A ray entering normal to the first face continues undeviated, so it strikes the next face at exactly $60° - 0° = 60°$. Since $60° > 41.41°$, total internal reflection is unavoidable.

What the critical angle means is easiest to see in a picture. When light passes from a denser to a less dense medium, the refraction angle is larger than the angle of incidence. As the angle of incidence grows, the refracted ray tilts towards the surface, and when the refraction angle reaches 90° it grazes along it (Fig 2). The critical angle is the angle of incidence at which the refracted ray grazes the surface. Beyond it there is no refracted ray left to carry the light away, so all of it comes back. In the code, that is the case `abs(sin_out) > 1`.

<img src="/assets/img/posts/raytracing-total-internal-reflection/ext-critical-angle-diagram.png" alt="Rays at increasing angles of incidence, from refraction through grazing exit to total reflection" width="700">
_Fig 2. Rays leaving a denser medium (below, n1) for a less dense one (above, n2) at increasing angles of incidence. At θ1 the ray refracts out at θ2; at the critical angle θc the refracted ray grazes the surface at 90°; beyond it the light is totally reflected. The faint lines are the partially reflected light. (Source: Jfmelero (adapted by Gavin R Putland), [ReflexionTotal en](https://commons.wikimedia.org/wiki/File:ReflexionTotal_en.svg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

The higher the index, the smaller the critical angle. At 589 nm it is 48.6° for water (n = 1.333), 43.3° for fused quartz (1.4584), 41.2° for N-BK7 (1.5168 at the d line) and 24.4° for diamond (2.417), using the indices in Hecht's Table 4.1 and the SCHOTT datasheet.

How, then, can light get through this prism without total internal reflection? By entering at a large enough angle of incidence. With $\theta_1$ the angle of incidence, $\theta_1'$ the refraction angle at the first face, $\theta_2'$ the internal angle at the second face and $A$ the apex angle, the relation $\theta_1' + \theta_2' = A$ holds. Staying below the critical angle at the second face therefore requires

$$ \theta_1' > A - \theta_c, \qquad \theta_1 > \arcsin\big(n \sin(A - \theta_c)\big) $$

The ray must be bent far enough at the first face for its angle at the second face to drop below the critical angle. This holds for rays arriving on the side of the first-face normal away from the apex (from below, as in Fig 3). A ray arriving from the apex side has $\theta_2' = A + \theta_1'$ and is totally reflected whatever its angle of incidence. Substituting the values gives $\theta_1 > 28.8131°$. Bisecting over the angle of incidence from 0° to find the TIR-to-transmission boundary the code actually produces (`verify_tir.py`) returns 28.813100°, within $7.6\times10^{-10}$° of the theoretical value.

<img src="/assets/img/posts/raytracing-total-internal-reflection/en/fig2-prism-angle-sweep.png" alt="Total internal reflection and transmission across a sweep of incidence angles" width="600">
_Fig 3. Six rays entering the same prism at angles of incidence from 0° to 45°. Up to 28.813°, just below the threshold of 28.8131°, they are totally reflected at the second face; at 35° and 45° they pass through it._

The exit angles of the transmitted rays show what the threshold means. Measured from the normal of the second face, the exit angle is 53.5° for an angle of incidence of 45°, 67.6° for 35° and 80.3° for 30°, and 89.3° at 28.82°, just above the threshold. The closer the ray comes to the threshold, the more nearly it grazes the face on the way out. The threshold is where the exit angle reaches 90°, which is the definition of the critical angle from Fig 2.

This is not an exotic effect. The Porro prism used in binoculars folds the light path by total internal reflection alone, with no mirror coating. Light entering a right-angle prism normal to its hypotenuse meets each of the two short faces, which sit at 45° to the hypotenuse, at an internal angle of 45°. Since 45° > 41.4° (the N-BK7 critical angle at 750 nm found above; 41.2° at the d line, 587.6 nm), the light must be totally reflected, and the faces act as mirrors without any coating. After the two reflections the light leaves in the direction opposite to the one it came from (Fig 4).

<img src="/assets/img/posts/raytracing-total-internal-reflection/ext-porro-prism-photo.jpg" alt="A light beam reflected twice by total internal reflection inside a right-angle prism" width="600">
_Fig 4. A light beam entering a right-angle prism normal to its hypotenuse. It is totally reflected at each short face in turn and leaves travelling back the way it came; no light is seen escaping through the short faces. (Source: Zátonyi Sándor, (ifj) Fizped, [Visszatükrözés prizmán 2](https://commons.wikimedia.org/wiki/File:Visszat%C3%BCkr%C3%B6z%C3%A9s_prizm%C3%A1n_2.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

For a 45° face to reflect totally, $\sin 45° > 1/n$, that is $n > \sqrt{2} \approx 1.414$. Fused quartz (1.4584) and N-BK7 both clear this; water (1.333) does not. For the same reason a wet reflecting face stops working. At a glass-water boundary the critical angle rises to $\arcsin(1.333/1.5168) = 61.5°$, so light arriving at 45° escapes into the water wherever the face is wet. A TIR face needs no coating, but it does need to be kept clean.

## 2. Which rays are totally reflected at the lens exit surface?

What breaks, and for which rays, when the lens lacks the same test? The exit-surface refraction in the lens routine reads as follows (carried over unchanged from [post 1](/en/posts/raytracing-spherical-lens-refraction/)).

```python
phi2 = np.arctan2(vy1, vx1)
theta_t2 = np.arcsin(n_lens * np.sin(w2 - phi2) / n_air)  # no domain check
vx2, vy2 = -np.cos(theta_t2 + w2), -np.sin(theta_t2 + w2)
```

When `n_lens * sin(w2 - phi2) / n_air` exceeds 1, that is, when the internal angle exceeds the critical angle, this `arcsin` has no real solution at all. The original code does not intercept the case, so a ray past the critical angle produces `NaN` with no exception, only a single `RuntimeWarning`.

The function does not decide on its own which surface the ray enters through, so the caller has to pass the entry surface (`entry_surface`) found by `find_entry_surface` (see [post 1](/en/posts/raytracing-spherical-lens-refraction/)). Passing the wrong one gives a wrong path or a `NaN` that has nothing to do with total internal reflection.

The test case is a plano-convex lens: R1 = $10^6$ mm, effectively flat, R2 = -50 mm, T = 20 mm, D = 78 mm. A collimated ray enters the flat side at a height of y0 = 36 mm.

```text
find_entry_surface: surface 1

before the fix:
[[-9.99935200e+00  3.60000000e+01  1.00000000e+00 -1.21878769e-05]
 [-5.30123745e+00  3.59999427e+01             nan             nan]]
-> RuntimeWarning: invalid value encountered in arcsin
```

The last line is that warning. The fix follows the prism: outside the domain, the ray is reflected about the exit-surface normal $(\cos w_2, \sin w_2)$ rather than refracted. The normal at a point on a sphere lies along the radius, so the circle's parameter angle $w_2$ can be used directly. With R2 negative, as here, this vector points from the surface towards the centre, but as shown above the reflection formula does not care about the sign of the normal.

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

Two checks confirm that this is a genuine reflection (`verify_tir.py`). First, the law of reflection holds about the normal: the angles of incidence and reflection are both 46.05508420°, differing by $7.1\times10^{-15}$°. Second, the ray height at which TIR sets in was found by bisection, and the internal angle at that height was compared with the theoretical critical angle. The code puts the boundary at y0 = 33.072057 mm, where the internal angle is 41.410389°, matching $\theta_c = 41.410389°$ to within $2.7\times10^{-10}$°.

That boundary also follows directly from geometry. A collimated ray entering the flat side is not bent there, so a ray at height $h$ meets the curved surface at an angle to the normal with $\sin\theta = h/\lvert R_2\rvert$. TIR begins when this exceeds $1/n$, at $h = \lvert R_2\rvert/n = 33.07$ mm.

<img src="/assets/img/posts/raytracing-total-internal-reflection/en/fig3-lens-marginal-tir.png" alt="Total internal reflection of the outermost rays in a plano-convex lens, and where they go next" width="600">
_Fig 5. Thirty-nine collimated rays, 2 mm apart, entering the flat side of a plano-convex lens (R2 = -50 mm, T = 20 mm, N-BK7, 750 nm). The six outside a height of 33.07 mm (red) are totally reflected at the exit surface. The lens function does not compute what happens after the reflection, so the red paths beyond that point come from a separate loop written for this figure, with the lens rim treated as a smooth surface._

Over the semi-aperture of 39 mm (D = 78 mm), TIR occurs from a height of 33 mm out to the edge: the outer 15% of the radius, or 28% of the area of a circular aperture. Fig 5 also shows where the reflected light ends up. The four rays at 36 mm and 38 mm meet the opposite side of the curved surface at the same angle as before (46.06° and 49.46°), are totally reflected a second time and leave through the flat face back towards the source; the two ends of a chord make equal angles with the radii of a circle. The two rays at 34 mm are reflected once and leave through the rim of the lens. None of them reach the focus, so the totally reflected light contributes nothing to the image.

When does this need worrying about? TIR enters the aperture when the boundary height $\lvert R_2\rvert/n$ is smaller than the semi-aperture D/2. The focal length of a plano-convex lens is $f = \lvert R_2\rvert/(n-1)$, so in terms of the f-number the condition reads $f/D < n/\big(2(n-1)\big)$, about 1.48 for N-BK7 at 750 nm. The lens here, with f = 97.7 mm, is F/1.25 and well inside that range. A fair rule of thumb, then: an N-BK7 plano-convex lens taking collimated light on its flat side shows edge TIR when it is faster than about F/1.48. That is why the effect is rarely seen with gently curved lenses. The opposite orientation, with the curved side facing the source, was not computed here.

## 3. What does the corrected lens function still leave out?

Unlike the prism routine, the corrected lens function still handles exactly one entry and one exit. The direction of a ray totally reflected at the exit surface is now correct, but the function does not compute where that ray goes next inside the lens, because it has no loop that runs until the ray stops meeting boundaries. To draw the paths in Fig 5, such a loop was added to the figure code, and its first intersection and direction agree with the lens function's exit row to within $2.5\times10^{-14}$ in every coordinate and direction component. Tracing reflected rays onward inside the library itself would mean rewriting the lens function along the lines of the prism routine.

## Summary

Detecting total internal reflection comes down to one test: whether the value passed to `arcsin` exceeds 1. If it does, the ray is reflected with $\vec v' = \vec v - 2(\vec v\cdot\hat n)\hat n$ instead of refracted. An equilateral N-BK7 prism at 750 nm always reflects totally at its second face for angles of incidence below 28.8°. A plano-convex lens taking collimated light on its flat side shows exit-surface TIR beyond a height of $\lvert R_2\rvert/n$, 33.07 mm for the lens here, and that light leaves towards the source or through the rim instead of reaching the focus.

The formula used twice in this post for reflection, once for the prism's internal reflection and once for the lens TIR, appears in exactly this form six times across the library: in `FlatMirror`, `SphericalMirror`, `ArbitraryMirror`, `BeamSplitter` (the reflected component), and now in the TIR branches of `Prism` and `SphericalLens`. All that differs between them is how the intersection and the normal are found (plane algebra, candidate roots on a circle, or local curvature along a polygon) and whether the routine reflects once and stops or loops until the ray escapes. Refraction needed a different angle expression every time; reflection, whether at a mirror, a prism or a lens, reduces to this one vector formula. [Post 3](/en/posts/raytracing-reflection-family/) checks whether the formula really is implemented consistently across the code, verifying `FlatMirror`, `SphericalMirror`, `ArbitraryMirror` and `BeamSplitter` together.

## References

- E. Hecht, *Optics*, 5th ed. (Pearson, 2017), §4.7 (total internal reflection and the critical angle), Table 4.1 (refractive indices at 589 nm), §5.5 (prisms: the apex-angle relation and the Porro prism).
- SCHOTT, N-BK7 optical glass datasheet (d-line index $n_d = 1.5168$). [SCHOTT optical glass](https://www.schott.com/en-us/products/optical-glass-p1000267)

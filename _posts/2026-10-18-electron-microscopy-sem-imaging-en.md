---
title: "Electron Microscopy Foundations 4 — SEM Imaging: the Detector Makes the Contrast"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-sem-imaging/
date: 2026-10-18 20:00:00 +0900
page_id: electron-microscopy-sem-imaging
categories: [Electron Microscopy, SEM Imaging]
tags: [electron-microscopy, sem, detector, edge-effect, charging]
description: An SEM has no lens that forms the image. What you end up seeing is decided by where the detector sits.
math: true
---

[Post 3](/en/posts/electron-microscopy-beam-specimen-interaction/) established that one beam striking one point yields signals from very different depths. Characteristic X-rays came from nearly the whole interaction volume, backscattered electrons from the upper third, secondary electrons from within 5 nm of the surface. So how are those signals actually collected and turned into an image?

Here the SEM parts ways with the optical microscope and the TEM. An optical microscope has a lens that maps each point of the specimen onto a point of the image plane. An SEM has no such lens. Its objective only converges the beam to a point; it does not spatially resolve the signal leaving that point. The detector amounts to a current meter with no sense of direction.

An image appears nonetheless. This post starts by building up how, then shows that where a detector sits, and over what angle it collects, is what defines contrast. It then computes the $\sec\theta$ law behind topographic contrast and the edge effect that follows from it, and closes with the charging problem on insulators and how low-voltage SEM answers it.

## 1. Time, not a lens, forms the image

If the detector reports only signal strength, what tells us which point on the specimen that signal came from?

Scanning does. Deflection coils sweep the beam across the specimen in a raster while the detector reads current continuously. Plotting the value read at the moment the beam sat at $(x_i, y_j)$ onto the same coordinate on screen builds the image. Position is held by the scan circuit that knows where the beam is; the detector supplies only brightness.

Several things follow from this arrangement.

No lens sets the magnification. Magnification is the screen size divided by the scanned width on the specimen. Shrinking the scan from 100 µm to 10 µm multiplies it by ten. Lens current only focuses, so the image rotation seen in post 2 stays decoupled from changing magnification in an SEM.

Depth of field is large. The aperture semi-angle is a few milliradians, so the beam diverges slowly on either side of focus. Unlike an optical microscope, where only a thin layer stays sharp at high magnification, an SEM renders structures of widely differing height sharply throughout. Half of why SEM photographs look three-dimensional sits here.

Time builds the image. The brightness of one pixel is the number of electrons gathered during the dwell time the beam spent there. Longer dwell improves signal-to-noise but stretches the exposure, and specimen drift over that period smears the result. Noise, exposure time and drift trade against one another, and negotiating them is the daily work of running an SEM.

## 2. Where the detector sits decides what it sees

If scanning handles position, what does the detector decide? Which signal gets collected. That choice is made by its position and its collection angle.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/en/fig1-detector-layout.png" alt="Layout of the ET detector, in-lens detector and annular BSE detector with secondary and backscattered electron paths" width="820">
_Fig 1. Layout of three detector types and the paths each signal takes_

Energy is what separates them. As post 3 showed, secondary electrons sit below 50 eV while backscattered electrons retain nearly the incident energy. A few tens of eV bends easily in a field of a few hundred volts; a few tens of keV barely notices the same field.

The Everhart-Thornley (ET) detector exploits that difference. Mounted off to the side with a grid held near $+200$ V, it draws in secondary electrons regardless of the direction they left the specimen, as the orange curves in Figure 1 show. Collection efficiency is high as a result, and even electrons emitted away from the detector are captured. Backscattered electrons, travelling straight, reach it only if they happen to head that way. An ET signal is therefore mostly secondary with a small backscattered admixture.

The annular BSE detector is designed the opposite way. It sits directly below the objective as a ring around the optical axis and carries no bias. Only electrons arriving straight along a line of sight register, which selects backscattered electrons exclusively, and the ring's axial symmetry keeps the signal free of directional bias.

The in-lens detector sits inside the objective, above the bore in the pole piece. It exploits the way the lens field spirals secondary electrons toward the axis, catching those drawn up through the bore. A shorter working distance strengthens the effect and raises collection efficiency, which helps particularly when signal runs short at low voltage. The cost is that the specimen must sit close to the lens, which rules out large or heavily tilted specimens.

Strike the same point with the same beam and change only which detector reads it, and the image changes. That is the practical meaning of saying no lens forms the image.

## 3. Topographic contrast — why tilted surfaces look bright

A secondary electron image from an ET detector shows surface relief because secondary electron yield depends on surface tilt. Why should it?

Post 3 established that only secondary electrons generated within about 5 nm of the surface escape. A beam entering perpendicular traverses that 5 nm layer over a path of 5 nm. If the surface is tilted by $\theta$, the beam crosses the same layer obliquely and the path becomes $5/\cos\theta$ nm. More energy is deposited within the escapable region, and proportionally more secondary electrons are produced.

$$ \frac{\delta(\theta)}{\delta(0)} = \sec\theta $$

<img src="/assets/img/posts/electron-microscopy-sem-imaging/en/fig2-edge-effect.png" alt="Secondary electron yield against tilt, and the edge effect appearing in a line scan across a trapezoidal feature" width="820">
_Fig 2. Secondary electron yield against tilt (left) and a line scan reproduced from that law (right)_

The left panel is that law: 1.41 at $45°$, 2 at $60°$, 3.86 at $75°$. The brightness difference between a flat floor and a steep sidewall comes from this curve.

As $\theta$ approaches $90°$ the secant diverges, but real surfaces do not. At steep incidence electrons pass out through the side face rather than remaining in the specimen, and no surface is perfectly flat at atomic scale. Hence the cap at six in the figure.

The right panel reproduces a line scan from this law alone. Sweeping across a 200 nm trapezoidal feature, the local surface tilt at each point gives $\sec\theta$ as the signal. The flat top and the floor return the same level; only the sidewalls rise, by nearly a factor of 2.5. This is the edge effect that makes outlines glow.

The effect is stronger at thin protrusions and sharp corners. There the interaction volume overlaps two or more free surfaces, so escape routes open in several directions at once. Particle edges burning out white in SEM photographs are this. It flatters the eye but troubles anyone measuring a linewidth, since the boundary becomes hard to place — which is why CD-SEM systems in semiconductor metrology carry dedicated algorithms for deciding a linewidth from the shape of this waveform.

## 4. Compositional contrast, and separating the two

What does a backscattered electron image show? The backscatter coefficient $\eta$ computed in post 3 rises monotonically with atomic number, and that alone becomes contrast. Heavier elements appear brighter, and a difference of one or two in atomic number is detectable.

The difficulty is that the two contrasts mix. Backscattered electrons also respond to surface tilt, and a secondary electron image carries a backscattered admixture. For clean compositional contrast, an annular BSE detector used symmetrically about the axis works better: summing the whole ring cancels directional dependence and leaves only the atomic number term.

Split the ring in half and take the difference instead, and the compositional term cancels while the tilt term survives. Choosing the sum or the difference from one detector to obtain a compositional image or a topographic image separately is the standard way of reading BSE signals.

## 5. Trading off voltage, working distance and aperture

If the preceding sections covered what becomes visible, sitting down at the instrument means choosing a few numbers. They interlock, so improving one degrades another.

Raising the accelerating voltage shrinks the beam diameter, since resolution scales as $\lambda^{3/4}$ as post 2 showed and the wavelength falls. But the interaction volume grows, as post 3 showed, so signal emerges from a wider region and fine surface detail is buried instead. Seeing the surface argues for lowering the voltage.

Shortening the working distance shortens the objective's focal length, lowering the spherical aberration coefficient and improving resolution. In-lens collection efficiency improves too. Depth of field suffers in exchange, and the specimen sits close enough to risk striking the lens.

Closing the aperture shifts the balance of post 2 away from spherical aberration, shrinking the beam and deepening the depth of field. The price is beam current: signal weakens, so holding the same signal-to-noise ratio means a longer dwell and a slower exposure.

All three are the same trade among resolution, signal and depth of field. High-resolution surface imaging pushes toward low voltage, short working distance, small aperture and long dwell; a quick overview of a rough specimen pushes the other way.

## 6. Charging — the assumption that fails on insulators

A quiet assumption has run underneath everything so far: that electrons entering the specimen leave it somewhere. In a metal they drain to ground and cause no trouble. In an insulator?

They accumulate at the surface. The field from that accumulated negative charge repels and deflects the beam that follows. Images bloom, streak, or shift bodily as the beam is thrown off. It is the first wall anyone imaging an insulator meets.

The traditional answer is a few nanometres of sputtered gold or platinum to provide a conductive path. Simple, but it damages the specimen, and the coating itself obscures structure at the nanometre scale, ruling it out for high-resolution work.

Another route exists. Balance the number of electrons arriving against the number leaving. The outgoing side is the sum of the secondary electron yield $\delta$ and the backscatter coefficient $\eta$, so a total yield $\sigma_T = \delta + \eta$ of exactly one leaves zero net accumulated charge.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/en/fig3-charging-crossover.png" alt="Total electron yield against accelerating voltage and the two crossover points where charging vanishes" width="800">
_Fig 3. Total yield curves for two insulators_

The total yield rises at low voltage, peaks at a few hundred eV, and falls again. At low voltage the electron dumps its energy just below the surface and generates many escapable secondaries; at high voltage it penetrates deep and the secondaries it generates cannot get out. The escape-depth argument of post 3 carries straight over.

The curve therefore crosses one twice, at a lower crossover $E_1$ and an upper crossover $E_2$. Between them $\sigma_T > 1$, more electrons leave than arrive, and the specimen charges positive; outside them $\sigma_T < 1$ and it charges negative. For $\mathrm{SiO_2}$ the upper crossover sits at 2.72 kV, for PMMA near 1.54 kV.

What makes $E_2$ special is that it is stable. Drop slightly below it and the specimen charges positive, raising the surface potential, which recaptures secondary electrons and reduces the outgoing count until $\sigma_T$ returns to one. Rise slightly above and negative charging decelerates the incoming electrons, pulling their effective landing energy back toward $E_2$. Feedback acts in both directions, converging on $E_2$.

This is what low-voltage SEM does. Working in the 1–2 kV range images insulators without any coating, which opens the way to viewing biological specimens, polymers and semiconductor photoresist in their original state. The cost is real. The electron wavelength lengthens and the diffraction limit worsens, chromatic aberration grows at low energy, and beam current falls with it. The aberration correctors of post 2 and the narrow energy spread of the cold field emission gun in post 1 earn their keep here in particular.

## Summary and what comes next

An SEM has no lens mapping specimen points onto image points. A scan circuit that pairs time with position does that job while the detector supplies only brightness. Magnification therefore follows the scanned width rather than any lens, and the small aperture angle grants a large depth of field.

What becomes visible is decided by where the detector sits. Placed to the side with a field applied, it draws in low-energy secondary electrons and shows topography; placed as a ring around the axis, it filters for straight-travelling backscattered electrons and shows composition. Topographic contrast rests on the path through the escape layer lengthening as $\sec\theta$, and that alone reproduces the edge effect that brightens sidewalls.

On insulators the assumption that charge drains away fails. Imaging near $E_2$, where the total yield passes one, balances electrons arriving against electrons leaving and permits imaging without a coating, and the feedback stability of that point is what makes low-voltage SEM practical.

That covers reading signals emitted from a specimen surface. From the next post onward the specimen is thinned until electrons pass through it. Post 5 constructs the diffraction pattern a thin crystal produces using the Ewald sphere, and shows how selecting the transmitted beam or a diffracted beam with the objective aperture separates bright field from dark field. How defects such as dislocations and grain boundaries appear in a TEM image is decided there.

## References

- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 6–10 (image formation, detectors, contrast mechanisms).
- L. Reimer, *Scanning Electron Microscopy: Physics of Image Formation and Microanalysis*, 2nd ed., Springer, 1998, ch. 4–6 (angular dependence of secondary electron yield, detector efficiency).
- D. C. Joy and C. S. Joy, "Low voltage scanning electron microscopy," *Micron*, vol. 27, pp. 247–263, 1996 (total yield curves and the $E_2$ crossover).
- H. Seiler, "Secondary electron emission in the scanning electron microscope," *Journal of Applied Physics*, vol. 54, pp. R1–R18, 1983 (standard account of secondary electron emission).
- [SEM Detectors: Secondary vs Backscattered Electrons, Tescan](https://tescan.com/news/sem-detectors-secondary-backscattered-electrons) (comparison of detector types).
- [Different Types of SEM Imaging — BSE and Secondary Electron Imaging, AZoM](https://www.azom.com/article.aspx?ArticleID=14309) (worked examples of topographic and compositional contrast).

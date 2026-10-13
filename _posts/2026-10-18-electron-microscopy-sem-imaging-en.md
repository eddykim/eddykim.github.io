---
title: "Electron Microscopy Foundations 4 — Why Do Electron Micrographs Look Three-Dimensional in Black and White?"
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

[Post 3](/en/posts/electron-microscopy-beam-specimen-interaction/) established that one beam striking one point yields signals from very different depths. Characteristic X-rays came from nearly the whole interaction volume, backscattered electrons (BSE) from the upper third, and secondary electrons (SE) from within 5 nm of the surface (a value for metals). So how are those signals actually collected and turned into an image?

Here the scanning electron microscope (SEM) parts ways with the optical microscope and the transmission electron microscope (TEM). An optical microscope has a lens that maps each point of the specimen onto a point of the image plane. An SEM has no such lens. Its objective only converges the beam to a point; it does not spatially resolve the signal leaving that point. The detector amounts to a current meter with no sense of direction.

An image appears nonetheless. This post starts by building up how, then shows that where a detector sits, and over what angle it collects, decides the contrast. It then computes the $\sec\theta$ law behind topographic contrast and the tilt contrast that follows from it. It closes with the charging problem on insulators and how low-voltage SEM solves it.

## 1. How does a picture appear on screen when no lens forms an image?

If the detector reports only signal strength, what tells us which point on the specimen that signal came from?

Scanning does. Deflection coils sweep the beam across the specimen in a raster while the detector reads current continuously. Plotting the value read at the moment the beam sat at $(x_i, y_j)$ onto the same coordinate on screen builds the image. Position is held by the scan circuit that knows where the beam is; the detector supplies only brightness. Each pixel ends up with a single number, so an SEM image is inherently black and white. The colours seen in magazines are added afterwards.

Several things follow from this arrangement.

No lens sets the magnification. Magnification is the screen size divided by the scanned width on the specimen. Shrinking the scan from 100 µm to 10 µm multiplies it by ten. Lens current only focuses, so the image rotation seen in [post 2](/en/posts/electron-microscopy-electron-optics/) stays decoupled from changing magnification in an SEM.

Depth of field is large. The aperture semi-angle is a few milliradians, so the beam diverges slowly on either side of focus. A Korean SEM textbook gives the estimate $D \approx 0.2/(\alpha M)$ mm, where $\alpha$ is the aperture semi-angle in radians and $M$ the magnification. With a 50 µm aperture radius and a 10 mm working distance, $\alpha$ = 5 mrad, and the depth of field is 40 µm at 1,000× and 4 µm at 10,000×. For an optical objective of numerical aperture 0.9, the wave-optical term $\lambda/\mathrm{NA}^2$ alone is about 0.7 µm at 550 nm. At the same 1,000×, then, the SEM sees tens of times deeper, and structures with large height differences stay sharp throughout. That is half of why a black-and-white SEM image looks three-dimensional. The other half, the shading, is the subject of section 3.

Time builds the image. The brightness of one pixel is the number of electrons gathered during the dwell time the beam spent there. A 1024×768 frame takes 7.9 s at 10 µs per pixel and 0.8 s at 1 µs. Longer dwell improves the signal-to-noise ratio (SNR) but stretches the exposure, and specimen drift over that period smears the result. Noise, exposure time and drift trade against one another, and negotiating them is the daily work of running an SEM.

## 2. Why does the same spot look different on each detector?

If scanning handles position, what does the detector decide? Which signal gets collected. That choice is made by its position and its collection angle.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/en/fig1-detector-layout.png" alt="Layout of the ET detector, in-lens detector and annular BSE detector with secondary and backscattered electron paths" width="820">
_Fig 1. Schematic of three detector types in an SEM chamber (not to scale). Secondary electrons (orange) are bent in by the +200 V grid of the ET detector and by the objective field; backscattered electrons (red) travel straight and reach the annular detector only along a line of sight._

Energy is what separates them. As post 3 showed, secondary electrons sit below 50 eV, while most backscattered electrons keep an energy close to the incident one. A few tens of eV bends easily in a field of a few hundred volts; a few tens of keV barely notices the same field.

The Everhart–Thornley (ET) detector exploits that difference. Mounted off to the side with a grid held near $+200$ V, it draws in secondary electrons regardless of the direction they left the specimen, as the orange curves in Figure 1 show. Collection efficiency is high as a result, and even electrons emitted away from the detector are captured. Backscattered electrons, travelling straight, reach it only if they happen to head that way. Direct backscattered electrons are only a small share. Secondaries generated where backscattered electrons leave the specimen or strike the pole piece and chamber walls are drawn in too, however, so an ET signal is a mixture carrying substantial backscatter information.

How do the collected electrons become a current signal? A few electrons of a few eV are far too weak to measure directly, so the ET detector converts them to light and amplifies that. The front grid, a Faraday cage, can be biased anywhere from −50 to +250 V. Electrons passing it are accelerated again by the +10 to 12 kV on the scintillator and strike it, producing light. The light travels down a light guide to a photomultiplier tube (PMT), where it becomes a large current (Figure 2). The ET detector, in short, amplifies the signal by turning electrons into light and light back into electrons.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/ext-everhart-thornley-detector.jpg" alt="Cross-section schematic of an ET detector with Faraday cage, scintillator, light guide and photomultiplier" width="700">
_Fig 2. Cross-section schematic of an ET detector. Secondary electrons from the sample at lower right pass the Faraday cage (+200 V) and strike the scintillator (+10 kV); the light (hν) runs along the optical waveguide to the photocathode and dynodes of the photomultiplier at left, where it becomes a current. (Source: Emmanuel de Chambost (France), [Everhart-Thornley detector](https://commons.wikimedia.org/wiki/File:Everhart-Thornley_detector.JPG), public domain)_

The arrangement offers a useful trick. Bias the cage negative (−50 V) and secondary electrons are repelled, leaving only backscattered electrons that travel straight at the detector. The same ET detector then gives a backscatter image.

The annular BSE detector is designed the opposite way. It sits directly below the objective as a ring around the optical axis and applies no field to draw electrons in, so only electrons arriving straight along a line of sight reach it. A semiconductor detector also responds only to electrons that penetrate its entrance electrode, a threshold of roughly 1–3 keV. Secondary electrons of a few eV therefore leave no signal, and only backscattered electrons register. The ring's axial symmetry keeps the signal free of directional bias.

The in-lens detector sits inside the objective, above the bore in the pole piece. It exploits the way the lens field spirals secondary electrons toward the axis, catching those drawn up through the bore. A shorter working distance strengthens the effect and raises collection efficiency, which helps particularly when signal runs short at low voltage. The cost is that the specimen must sit close to the lens, which makes large or heavily tilted specimens difficult.

Strike the same point with the same beam and change only which detector reads it, and the image changes. What an SEM shows is chosen by the detector, not by a lens.

## 3. Why do tilted surfaces look bright?

A secondary electron image from an ET detector shows surface relief because secondary electron yield depends on surface tilt. Why should it?

Post 3 established that only secondary electrons generated within about 5 nm of the surface escape. That is a value for metals; in insulators the depth reaches tens of nanometres. A beam entering perpendicular traverses that 5 nm layer over a path of 5 nm. If the surface is tilted by $\theta$, the beam crosses the same layer obliquely and the path becomes $5/\cos\theta$ nm. More energy is deposited within the escapable region, and proportionally more secondary electrons are produced.

$$ \frac{\delta(\theta)}{\delta(0)} = \sec\theta $$

Here $\delta(\theta)$ is the secondary electron yield of a face tilted by $\theta$, the number of secondaries emitted per incident electron, and $\delta(0)$ is the value for a flat face. The steeper the face, the more energy lands in the escape layer and the brighter that face appears.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/en/fig2-edge-effect.png" alt="Secondary electron yield against tilt, and the tilt contrast appearing in a line scan across a trapezoidal feature" width="820">
_Fig 3. Left: the $\sec\theta$ law from the path length through the escape layer (capped at six above 80°). Right: the signal from scanning a 200 nm high trapezoid with that law alone; it peaks at 2.44 times the flat level at the middle of each sidewall (±300 nm, tilt 65.8°). The edge effect at corners is not included._

The left panel is that law: 1.41 at $45°$, 2 at $60°$, 3.86 at $75°$. The brightness difference between a flat floor and a steep sidewall comes from this curve.

As $\theta$ approaches $90°$ the secant diverges, but the measured yield does not. At steep incidence electrons pass out through the side face rather than remaining in the specimen, and no surface is perfectly flat at atomic scale. Hence the cap at six in the figure. Measured yields also rise more slowly than the secant.

The right panel reproduces a line scan from this law alone. Sweeping across a 200 nm trapezoidal feature, the local surface tilt at each point gives $\sec\theta$ as the signal. The flat top and the floor return the same level; only the sidewalls rise, up to 2.44 times. The peak sits at the middle of the sidewall, where the tilt is greatest. This tilt (inclination) contrast is the first reason outlines glow. The calculation contains only the tilt; the edge effect at corners adds on top of it.

The edge effect proper appears at thin protrusions and sharp corners. There the interaction volume overlaps two or more free surfaces, so escape routes open in several directions at once. Particle edges burning out white in SEM photographs are this. It flatters the eye but troubles anyone measuring a linewidth, since the boundary becomes hard to place. That is why critical-dimension SEMs (CD-SEM) in semiconductor metrology carry dedicated algorithms for deciding a linewidth from the waveform.

Most of the shading that makes a black-and-white micrograph look solid, then, comes from tilt contrast and the edge effect.

## 4. Why do heavy elements look bright, and how are topography and composition separated?

What does a backscattered electron image show? The backscatter coefficient $\eta$ computed in post 3 rises monotonically with atomic number, and that alone becomes contrast. Heavier elements appear brighter; this is compositional contrast.

How close can two elements be and still be told apart? Among light elements, a difference of one or two in atomic number is detectable. It gets harder further up the table. Neighbour-element contrast is 7.3% for Al/Si but falls to 2.6% for Cu/Zn and 0.62% for Pt/Au. Separating heavy neighbours takes more beam current or a longer frame.

The difficulty is that the two contrasts mix. Backscattered electrons also respond to surface tilt, and a secondary electron image carries backscatter information through the secondaries that backscattered electrons generate. For clean compositional contrast, an annular BSE detector used symmetrically about the axis works better. Summing the whole ring greatly reduces the directional dependence and leaves the atomic number term dominant, although strongly tilted faces still appear dark in the sum.

Split the ring in half and take the difference instead, and the compositional term cancels while the tilt term survives. Taking the sum and the difference from one detector to obtain separate compositional and topographic images is the standard way of reading BSE signals.

In a real image the difference between the two contrasts is plain to see (Figure 4).

<img src="/assets/img/posts/electron-microscopy-sem-imaging/ext-se-bse-comparison.jpg" alt="Secondary electron and backscattered electron images of the same area of a tin oxide varistor ceramic" width="760">
_Fig 4. Secondary electron (left) and backscattered electron (right) images of the same area of a SnO₂ varistor ceramic. In the SE image the grain-boundary edges stand out bright; in the BSE image the same boundaries remain only as thin lines, and the triangular precipitate at upper left appears dark because its mean atomic number is lower than its surroundings. The slight brightness differences between grains appear to come from crystal orientation, which this post does not cover. (Source: Mateus Gallucci Masteghin & Marcelo Ornaghi Orlandi, [1516-1439-mr-1980-5373-MR-2016-0210-gf02](https://commons.wikimedia.org/wiki/File:1516-1439-mr-1980-5373-MR-2016-0210-gf02.jpg), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0))_

## 5. What do voltage, working distance and aperture each buy, and what do they cost?

The preceding sections covered what becomes visible. Sitting down at the instrument means choosing a few numbers. They interlock, so improving one degrades another. On what basis are they chosen?

Raising the accelerating voltage shrinks the beam diameter, since resolution scales as $\lambda^{3/4}$ as post 2 showed and the wavelength falls. Gun brightness also rises in proportion to the voltage, putting more current into the same spot, and the relative energy spread $\Delta E/E$ falls, reducing chromatic aberration; both count for a good deal. But the interaction volume grows, as post 3 showed, so signal emerges from a wider region and fine surface detail is buried instead. Seeing the surface argues for lowering the voltage.

Shortening the working distance shortens the objective's focal length, lowering the spherical aberration coefficient and improving resolution. In-lens collection efficiency improves too. In exchange the aperture semi-angle grows, so depth of field suffers, and the specimen sits close enough to risk striking the lens.

Closing the aperture shifts the balance of post 2 away from spherical aberration, shrinking the beam and deepening the depth of field. That holds only while the aperture is larger than its optimum; close it further and diffraction thickens the beam again. Beam current also falls, so holding the same signal-to-noise ratio means a longer dwell and a slower exposure.

| Change | Resolution | Signal | Depth of field | Reason |
|---|---|---|---|---|
| Raise the voltage | Smaller beam, but fine surface detail is buried | Larger | Nearly unchanged | Wavelength, gun brightness, $\Delta E/E$ / larger interaction volume |
| Shorten the working distance | Better | Higher in-lens collection | Shallower | Shorter focal length and smaller spherical aberration coefficient / larger aperture semi-angle |
| Close the aperture (down to the optimum) | Better | Smaller | Deeper | Less spherical aberration / less beam current / smaller aperture semi-angle |

All three are the same trade among resolution, signal and depth of field. High-resolution surface imaging pushes toward low voltage, short working distance, small aperture and long dwell. A quick overview of a rough specimen pushes the other way.

## 6. Why do insulators bloom and jump under the beam?

A quiet assumption has run underneath everything so far: that electrons entering the specimen leave it somewhere. In a metal they drain to ground and cause no trouble. In an insulator?

They accumulate at the surface. The field from that accumulated negative charge repels and deflects the beam that follows. Images bloom, streak, or shift bodily as the beam is thrown off (Figure 5). It is the first wall anyone imaging an insulator meets.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/ext-charging-uncoated.jpg" alt="Charging artefacts in an SEM image of an uncoated red alga gemma" width="600">
_Fig 5. SEM image of an uncoated gemma of the red alga Hildenbrandia rivularis. The horizontal streaks, the displaced scan lines and the washed-out bright patches are all traces of accumulated charge disturbing the beam and the signal. (Source: Smith609, [Hildenbrandia rivularis gemma SEM](https://commons.wikimedia.org/wiki/File:Hildenbrandia_rivularis_gemma_SEM.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

The traditional answer is a few nanometres of sputtered gold or platinum to provide a conductive path. It is simple, but it covers the surface irreversibly and can harm heat-sensitive specimens during deposition, and a thick or island-forming coating such as pure gold hides nanometre-scale structure. High-resolution work therefore uses only 1–2 nm of platinum, chromium or iridium.

Another route exists. Balance the number of electrons arriving against the number leaving. The outgoing side is the sum of the secondary electron yield $\delta$ and the backscatter coefficient $\eta$. Call that sum the total yield $\sigma_T = \delta + \eta$, the number of electrons leaving per incident electron. Each incident electron then leaves a charge of $1-\sigma_T$ behind, and a total yield of exactly one leaves no net charge.

<img src="/assets/img/posts/electron-microscopy-sem-imaging/en/fig3-charging-crossover.png" alt="Total electron yield against accelerating voltage and the two crossover points where charging vanishes" width="800">
_Fig 6. Total yield of SiO₂ (solid) and PMMA (dashed) from the universal yield curve model with representative parameters. SiO₂ crosses one at 72 eV and 2.72 kV, PMMA at 66 eV and 1.54 kV. The shading marks where SiO₂ charges positive._

The total yield rises at low voltage, peaks at a few hundred eV, and falls again. At low voltage the electron dumps its energy just below the surface and generates many escapable secondaries. At high voltage it penetrates deep, and the secondaries generated there cannot reach the surface. The escape-depth argument of post 3 carries straight over.

The curve therefore crosses one twice, at a lower crossover $E_1$ and an upper crossover $E_2$. Between them $\sigma_T > 1$, more electrons leave than arrive, and the specimen charges positive (the shading in Figure 6); outside them $\sigma_T < 1$ and it charges negative. In the model curves of Figure 6 (the Lin–Joy universal yield curve with representative parameters), the upper crossover comes out at 2.72 kV for $\mathrm{SiO_2}$ and 1.54 kV for PMMA. Reported measurements are about 3.0 kV for quartz and about 1.6 kV for PMMA, with wide scatter between methods and specimens.

Where does the usual 30 kV fall? At the arrow in Figure 6. For the $\mathrm{SiO_2}$ model, $\sigma_T$ is 0.32 there, so 68% of the incoming electrons stay in the specimen. Even at 5 kV, 28% stay. At 2 kV, by contrast, $\sigma_T$ is 1.20 and the specimen charges positive. That is why an insulator viewed at 30 kV blooms so quickly.

What makes $E_2$ special is that it is stable. Drop slightly below it and the specimen charges positive, raising the surface potential, which recaptures secondary electrons and reduces the outgoing count until $\sigma_T$ returns to one. Rise slightly above and negative charging decelerates the incoming electrons, pulling their effective landing energy back toward $E_2$. Feedback acts in both directions, converging on $E_2$.

This is what low-voltage SEM does. Working in the 1–2 kV range images insulators without any coating, which opens the way to viewing biological specimens, polymers and semiconductor photoresist in their original state. The cost is real. The electron wavelength lengthens and the diffraction limit worsens, chromatic aberration grows at low energy, and beam current falls with it. Correctors that also handle chromatic aberration, among the correctors of post 2, and the narrow energy spread of the cold field emission gun in [post 1](/en/posts/electron-microscopy-why-electrons/) earn their keep here in particular.

## Summary and what comes next

An SEM has no lens mapping specimen points onto image points. A scan circuit pairs time with position, and the detector gives each pixel a single brightness. The image is therefore black and white, and its magnification follows the scanned width. With an aperture semi-angle of a few milliradians, depth of field reaches tens of micrometres even at 1,000×.

What becomes visible is decided by where the detector sits and how it is biased. An ET detector to the side with a positive grid draws in low-energy secondary electrons and shows topography; an annular detector around the axis takes only straight-travelling, high-energy backscattered electrons and shows composition. Topographic contrast rests on the path through the escape layer lengthening as $\sec\theta$. That alone makes the sidewalls of a 200 nm trapezoid 2.44 times brighter than the flat, and the edge effect at corners adds on top. The solid look of a black-and-white micrograph comes from this deep focus and this shading.

On insulators the assumption that charge drains away fails. In the $\mathrm{SiO_2}$ model, 68% of the incoming electrons stay at 30 kV. Lowering the voltage to near $E_2$ (2.72 kV in the model, about 3.0 kV measured), where the total yield passes one, balances electrons arriving against electrons leaving and permits imaging without a coating. That point is also stable under feedback. The detector chooses the contrast, and on an insulator the balance of charge chooses the voltage.

That covers reading signals emitted from a specimen surface. From the next post onward the specimen is thinned until electrons pass through it. Post 5 constructs the diffraction pattern a thin crystal produces using the Ewald sphere, and shows how selecting the transmitted beam or a diffracted beam with the objective aperture separates bright field from dark field. How defects such as dislocations and grain boundaries appear in a TEM image is decided there.

## References

- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 5–11 (instrumentation and detectors, image formation, contrast mechanisms, neighbour-element backscatter contrast, charging, low-voltage SEM).
- L. Reimer, *Scanning Electron Microscopy: Physics of Image Formation and Microanalysis*, 2nd ed., Springer, 1998, ch. 4–6 (angular dependence of secondary electron yield, detector efficiency).
- J. Yoon, C. Yang, J. Kim and S. Lee, *Scanning Electron Microscopy and X-ray Microanalysis* (in Korean), 2005, ch. 4 (p. 86 depth-of-field estimate, pp. 87–88 structure and bias voltages of the E-T detector).
- D. C. Joy and C. S. Joy, "Low voltage scanning electron microscopy," *Micron*, vol. 27, pp. 247–263, 1996 (total yield curves and the $E_2$ crossover).
- Y. Lin and D. C. Joy, "A new examination of secondary electron yield data," *Surface and Interface Analysis*, vol. 37, pp. 895–900, 2005 (the universal yield curve used in Figure 6).
- H. Seiler, "Secondary electron emission in the scanning electron microscope," *Journal of Applied Physics*, vol. 54, pp. R1–R18, 1983 (standard account of secondary electron emission).
- M. G. Masteghin and M. O. Orlandi, "Grain-Boundary Resistance and Nonlinear Coefficient Correlation for SnO₂-Based Varistors," *Materials Research*, CC BY 4.0. [SciELO article](https://www.scielo.br/j/mr/a/FQ79NmGxdhMqvxcWkb85F6h/?lang=en) (the SE and BSE images in Figure 4).
- [SEM Detectors: Secondary vs Backscattered Electrons, Tescan](https://tescan.com/news/sem-detectors-secondary-backscattered-electrons) (comparison of detector types).
- [Different Types of SEM Imaging — BSE and Secondary Electron Imaging, AZoM](https://www.azom.com/article.aspx?ArticleID=14309) (worked examples of topographic and compositional contrast; article supplied by Thermo Fisher Scientific Phenom-World).

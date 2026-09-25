---
title: "Electron Microscopy Foundations 3 — After the Beam Lands: Interaction Volume and the Origin of Each Signal"
lang: en
lang-exclusive: ["en"]
permalink: /posts/electron-microscopy-beam-specimen-interaction/
date: 2026-10-08 20:00:00 +0900
page_id: electron-microscopy-beam-specimen-interaction
categories: [Electron Microscopy, Electron Scattering]
tags: [electron-microscopy, monte-carlo, scattering, interaction-volume, backscattered-electron]
description: Narrow the beam to 1 nm and the signal still comes from a micrometre. This post follows that gap from cross sections to a Monte Carlo calculation.
math: true
---

Everything up to [post 2](/en/posts/electron-microscopy-electron-optics/) was about collecting electrons onto a single point on the specimen. Scherzer's theorem pins the aperture angle near 6.5 mrad, and resolution ends up trapped at a value proportional to $C_s^{1/4}\lambda^{3/4}$. Suppose all that effort narrows the beam to 1 nm. Does a 1 nm image follow?

It does not. Electrons do not stop at the surface. They travel inward, collide with atoms hundreds of times, spread sideways, and generate signal throughout the region they spread into. Fire 30 kV electrons at silicon and that region runs several micrometres in both depth and width. No amount of narrowing the beam removes this volume from the resolution budget.

Yet scanning electron micrographs do resolve structures at the nanometre scale. How do the two statements coexist? The answer lies in having lumped every "signal" together. This post builds up what happens to an electron inside the specimen, starting from the scattering cross section, computes the interaction volume directly with a Monte Carlo calculation, and then separates the depths from which secondary electrons, backscattered electrons and characteristic X-rays each emerge. Every contrast mechanism in post 4 grows out of that separation.

## 1. What does the electron collide with?

A specimen contains nuclei and orbital electrons. The incident electron interacts with both, but the outcomes differ.

Collisions with a nucleus are elastic. A nucleus outweighs an electron by thousands, so the electron loses almost no energy and simply changes direction, much as a billiard ball rebounds off a cushion. This scattering spreads electrons sideways and, applied often enough, turns them around entirely.

Collisions with orbital electrons are inelastic. The masses match, so energy transfers readily while the direction barely changes. The transferred energy generates secondary electrons inside the specimen, and knocking out an inner-shell electron produces a characteristic X-ray as the vacancy fills.

Elastic scattering changes direction; inelastic scattering removes energy. A trajectory is the interleaving of the two, and the calculation below is built from exactly these two pieces.

## 2. How often, and how sharply?

Quantifying elastic scattering requires a probability. Expressing as an area the likelihood that one atom deflects a passing electron gives the scattering cross section. The screened Rutherford model supplies it.

$$ \sigma \;\propto\; \frac{Z^2}{E^2} \cdot \frac{1}{\alpha_s(1+\alpha_s)}, \qquad \alpha_s = \frac{3.4\times10^{-3}\,Z^{0.67}}{E} $$

$Z$ is the atomic number and $E$ the electron energy in keV. The screening parameter $\alpha_s$ measures how far orbital electrons shield the nuclear charge. Nearly every conclusion in this post sits in the first two factors.

The $Z^2$ in the numerator says heavier elements deflect far more often and far more sharply. Tungsten ($Z=74$) exceeds silicon ($Z=14$) by a factor of 28 in $Z^2$. The $E^2$ in the denominator says slower electrons deflect more readily, so lowering the accelerating voltage costs an electron its direction before it penetrates far.

The cross section gives a mean free path. With $N$ atoms per unit volume,

$$ \ell_{el} = \frac{1}{N\sigma} = \frac{A}{N_A\,\rho\,\sigma} $$

and the distance to the next scattering event is drawn from an exponential distribution with this mean.

Energy loss comes from Bethe's stopping power.

$$ \frac{dE}{ds} = -78500\,\frac{\rho\,Z}{A\,E}\,\ln\!\left(\frac{1.166\,E}{J}\right) \quad [\mathrm{keV/cm}] $$

$J$ is the mean ionization potential. As $E$ falls the logarithm's argument drops below one and the sign flips, so low energies call for the Joy-Luo modification, which replaces the argument with $1.166(E+0.85J)/J$. The calculations below use it.

## 3. Monte Carlo — trajectories from a four-line rule

Two equations suffice to generate one electron's trajectory. No analytic solution is needed; random numbers and step-by-step tracking do the work.

```python
# core of monte_carlo.py (full code: _code/electron-microscopy-beam-specimen-interaction/)
lam = mean_free_path(mat, e)                     # 1. mean free path
step = -lam * np.log(rng.random())               # 2. step from an exponential draw
loss = -stopping_power(mat, e) * step            # 3. Bethe energy loss
e -= loss
a = screening(mat, e)                            # 4. angle from screened Rutherford
cos_t = 1.0 - 2.0 * a * r / (1.0 + a - r)        #    (r is a uniform random number)
direction = _scatter_direction(direction, cos_t, 2 * np.pi * rng.random())
```

Four lines. The loop repeats until the electron crosses back out of the surface or its energy falls below 0.5 keV. Crossing back out counts it as a backscattered electron (BSE).

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/en/fig1-interaction-volume.png" alt="Monte Carlo trajectories and interaction volumes for 5 kV and 30 kV electrons in silicon and tungsten" width="820">
_Fig 1. Interaction volume by material and accelerating voltage. Blue trajectories stop inside the specimen; red ones return through the surface_

The four panels display the two factors from section 2 directly.

Raising the accelerating voltage from 5 kV to 30 kV (left to right) inflates the whole volume. In silicon the penetration depth grows from 0.47 µm to 9.31 µm, a factor of 20. Higher energy means a smaller cross section, less deflection, and more energy to spend before stopping.

Raising the atomic number (top to bottom) shrinks the volume and reshapes it. At 30 kV silicon reaches 9.31 µm against 1.68 µm for tungsten, a fifth as far. The shape matters more. Silicon's volume is a teardrop that swells sideways below the surface, while tungsten's is closer to a hemisphere pressed against it. In silicon an electron travels some distance before losing its direction; in tungsten it deflects sharply on entry. That is the $Z^2$ at work.

The dashed semicircle marks the Kanaya-Okayama penetration depth.

$$ R_{KO} = \frac{0.0276\,A\,E_0^{1.67}}{Z^{0.89}\,\rho} \quad [\mu\mathrm{m}] $$

Note that most trajectories crowd into the inner half of that arc. $R_{KO}$ estimates the farthest distance reached by the few electrons that travel furthest; it is not a typical value. Most of the signal originates well above it.

## 4. Look at depth, not volume

So the interaction volume spans micrometres. How then does an SEM resolve nanometres?

Because not all of the signal generated in that volume is used. Each signal escapes from a different depth, and the differences span three orders of magnitude.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/en/fig2-depth-distribution.png" alt="Energy deposition, X-ray generation and backscattered electron depth distributions in silicon at 15 kV, with the secondary electron escape depth" width="820">
_Fig 2. Depth distributions by signal in silicon at 15 kV. The right panel magnifies the top 20 nm_

The blue solid curve on the left is energy deposited against depth. Secondary electrons (SE) arise from the energy handed over in inelastic scattering, so their generation follows this curve. It peaks near 0.5 µm and extends to 2.5 µm.

The green dashed curve keeps only those segments where the electron still carries more than 1.84 keV, the ionization energy of the silicon K shell, and can therefore produce a characteristic X-ray. It nearly coincides with the blue curve. An electron entering at 15 keV holds far more than 1.84 keV over most of its path, so X-rays emerge from essentially the entire interaction volume. This is why energy dispersive X-ray spectroscopy (EDS) is limited to micrometre spatial resolution when measuring composition.

The red dotted curve shows the maximum depth reached by electrons that eventually returned. It peaks near 0.3 µm and rarely exceeds 1 µm. Backscattered electrons sample only the upper third or so of the interaction volume.

The orange band marks the escape depth of secondary electrons. Their energy is very low, below 50 eV, so they scatter again and stop almost immediately inside the specimen. Only those generated within a few nanometres of the surface get out. The representative value of 5 nm amounts to **0.171 %** of the 2.93 µm penetration depth under these conditions. That is why the band appears as a single vertical line on the left, and why the top 20 nm is magnified separately on the right.

One more signal is missing from Figure 2: the Auger electron. When an incident electron ejects an inner-shell electron, an outer-shell electron drops in to fill the vacancy, and the leftover energy sometimes kicks out another orbital electron instead of leaving as an X-ray. That ejected electron is the Auger electron, and its energy is element-specific, so it carries compositional information. Its escape depth, 1–2 nm, is shallower still than that of secondary electrons. Auger electron spectroscopy rests on exactly this, which is also why the slightest surface contamination buries the signal. Attacking the same class of problem with X-rays instead of electrons gives [XPS](/en/posts/xps-photoemission-binding-energy/).

That resolves the opening contradiction. The interaction volume spans micrometres, but the information secondary electrons carry comes from a few nanometres of surface. Only secondary electrons generated in the narrow column directly beneath the entry point — before the beam has spread sideways — ever escape, so the resolution of an SE image approaches the beam diameter. An SEM shows topography at nanometre scale not by defeating the interaction volume but by viewing only a very thin shell of it.

## 5. The backscatter coefficient and compositional contrast

If secondary electrons show surface shape, what do backscattered electrons show? The $Z^2$ from section 2 returns here.

The fraction of incident electrons that return is the backscatter coefficient $\eta$. Since scattering scales with $Z^2$, $\eta$ should climb with atomic number. The Monte Carlo can check that.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/en/fig3-backscatter-vs-z.png" alt="Backscatter coefficient against atomic number from Monte Carlo compared with the Reuter empirical fit" width="750">
_Fig 3. Backscatter coefficient against atomic number. Each point tracks 3,000 electrons_

The values rise monotonically from 0.06 for carbon ($Z=6$) to 0.54 for gold ($Z=79$), steeply among light elements and more gently toward heavy ones. The grey dashed line is Reuter's empirical fit to measurements.

This curve underwrites compositional contrast. Place two phases of differing atomic number side by side and the number of returning electrons differs, so a BSE detector renders the heavier one brighter. The slope is what makes a difference of one or two in atomic number detectable.

Figure 3 also holds a discrepancy worth naming. At high atomic number the Monte Carlo points sit systematically above the empirical fit. Silicon agrees well, 0.168 against 0.164, but tungsten runs about 12 % high, 0.532 against 0.476. This is not statistical fluctuation but a limitation of the model. The screened Rutherford cross section treats the nucleus as a point charge and compresses orbital screening into a single constant, and the approximation degrades as atomic number grows. Accurate values require tabulated Mott cross sections from relativistic wave calculations.

Leaving the discrepancy visible seems the more honest choice. A Monte Carlo of this kind demonstrates trends, not measurements. It is enough to understand why the interaction volume takes the shape it does and which depth each signal comes from, but using it as a correction factor for quantitative analysis would mean replacing the cross section first.

## Summary and what comes next

What happens to an electron inside a specimen splits two ways. Elastic scattering off nuclei changes direction; inelastic scattering off orbital electrons removes energy. Almost every property of the interaction volume follows from the single fact that the elastic cross section scales as $Z^2/E^2$. Raising the accelerating voltage enlarges the volume; raising the atomic number shrinks it and turns the teardrop into a hemisphere.

That the volume spans micrometres while an SEM resolves nanometres comes down to escape depth. Characteristic X-rays emerge from nearly the whole volume, backscattered electrons from the upper third, secondary electrons from within 5 nm of the surface, and Auger electrons from a shallower 1–2 nm. The 5 nm for secondary electrons is 0.171 % of the penetration depth. The same beam striking the same point is read at depths differing by three orders of magnitude depending on which detector receives it.

The next post takes up those detectors. It builds up how scanning forms an image, then shows how the position and collection angle of secondary and backscattered electron detectors inside the column produce topographic and compositional contrast respectively. It also covers the edge effect, where secondary electron yield rises on tilted surfaces and brightens outlines, along with the charging problem encountered on insulators and how low-voltage SEM answers it.

## References

- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 1–5 (beam-specimen interaction, interaction volume, backscatter coefficient).
- D. C. Joy, *Monte Carlo Modeling for Electron Microscopy and Microanalysis*, Oxford University Press, 1995 (the standard single-scattering implementation).
- D. C. Joy and S. Luo, "An empirical stopping power relationship for low-energy electrons," *Scanning*, vol. 11, pp. 176–180, 1989 (low-energy Bethe modification).
- K. Kanaya and S. Okayama, "Penetration and energy-loss theory of electrons in solid targets," *Journal of Physics D*, vol. 5, pp. 43–58, 1972 (penetration depth fit).
- L. Reimer, *Scanning Electron Microscopy: Physics of Image Formation and Microanalysis*, 2nd ed., Springer, 1998, ch. 3–4 (scattering theory and signal generation depth).
- [Scanning Electron Microscopy, Carleton SERC](https://serc.carleton.edu/research_education/geochemsheets/techniques/SEM.html) (overview of escape depths by signal).

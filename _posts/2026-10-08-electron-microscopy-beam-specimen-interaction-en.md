---
title: "Electron Microscopy Foundations 3 — If the Beam Is 1 nm Wide, Can You See 1 nm Detail?"
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

Everything up to [post 2](/en/posts/electron-microscopy-electron-optics/) was about collecting electrons onto a single point on the specimen. Scherzer's theorem makes spherical aberration unavoidable, so the aperture angle is pinned at its balance with diffraction (6.5 mrad under the conditions of post 2), and resolution ends up trapped at a value proportional to $C_s^{1/4}\lambda^{3/4}$. Suppose all that effort narrows the beam to 1 nm. Does a 1 nm image follow?

It does not. Electrons do not stop at the surface. They travel inward, collide with atoms hundreds of times, spread sideways, and generate signal throughout the region they spread into. Fire 30 kV electrons at silicon and that region runs several micrometres in both depth and width. This region is called the interaction volume, and no amount of narrowing the beam removes it from the resolution budget.

Yet scanning electron microscope (SEM) images do resolve structures at the nanometre scale. How do the two statements coexist? The answer lies in having lumped every "signal" together. This post builds up what happens to an electron inside the specimen, starting from the scattering cross section, and computes the interaction volume directly with a Monte Carlo calculation. It then separates the depths from which secondary electrons (SE), backscattered electrons (BSE) and characteristic X-rays each emerge. Every contrast mechanism in post 4 grows out of that separation.

## 1. What does the electron collide with?

A specimen contains nuclei and orbital electrons. The incident electron interacts with both, but the outcomes differ.

Collisions with a nucleus are elastic. A nucleus outweighs an electron by tens to hundreds of thousands of times (about 50,000 for silicon), so the electron loses almost no energy and simply changes direction. A 30 keV electron can hand a silicon nucleus at most 2.4 eV in one collision, less than a ten-thousandth of its energy. It is much like a billiard ball rebounding off a cushion. This scattering spreads electrons sideways and, applied often enough, turns them around entirely.

Collisions with orbital electrons are inelastic. The masses match, so energy transfers readily while the direction barely changes. The transferred energy generates secondary electrons inside the specimen, and knocking out an inner-shell electron produces a characteristic X-ray as the vacancy fills. Fig 1 sets the three signals side by side on a single atom.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/ext-electron-emission-mechanisms.png" alt="Three ways an incident electron produces a secondary electron, a backscattered electron and a characteristic X-ray" width="380">
_Fig 1. Three ways a primary electron (PE) produces a signal from one atom. Top: it knocks out an outer orbital electron as a secondary electron (SE). Middle: it swings sharply around the nucleus and leaves as a backscattered electron (BSE). Bottom: an outer electron fills an inner-shell vacancy and emits a characteristic X-ray. (Source: Rob Hurt, [Electron emission mechanisms](https://commons.wikimedia.org/wiki/File:Electron_emission_mechanisms.svg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

Elastic scattering changes direction; inelastic scattering removes energy. A trajectory is the interleaving of the two, and the calculation below is built from exactly these two pieces.

## 2. How often, and how sharply?

Quantifying elastic scattering requires a probability. The scattering cross section is the probability that one atom deflects a passing electron, expressed as an area. The screened Rutherford model supplies it.

$$ \sigma \;\propto\; \frac{Z^2}{E^2} \cdot \frac{1}{\alpha_s(1+\alpha_s)}, \qquad \alpha_s = \frac{3.4\times10^{-3}\,Z^{0.67}}{E} $$

$Z$ is the atomic number and $E$ the electron energy in keV. The screening parameter $\alpha_s$ measures how far orbital electrons shield the nuclear charge. Most of the qualitative conclusions in this post sit in the leading $Z^2/E^2$ term.

The two parts of the formula, $Z^2/E^2$ and the $\alpha_s$ terms, answer the two questions in the heading separately. First, how sharply. The cross section for large-angle scattering follows the leading factor, $Z^2/E^2$, so for large-angle deflection tungsten ($Z=74$) exceeds silicon ($Z=14$) by a factor of 28 (27.7 counting only deflections beyond 90°).

Second, how often. That is a question for the total cross section, which counts small angles too. Since $\alpha_s$ itself carries $Z^{0.67}/E$, part of the leading factor cancels in this energy range, where $\alpha_s \ll 1$. The total cross section varies roughly as $Z^{4/3}/E$, and per atom tungsten scatters about ten times as often as silicon. In other words, heavy elements deflect electrons more often (about tenfold here) but, above all, far more sharply. The $E^2$ in the denominator says slower electrons deflect more readily, so lowering the accelerating voltage costs an electron its direction before it penetrates far.

The cross section gives a mean free path. With $N$ atoms per unit volume,

$$ \ell_{el} = \frac{1}{N\sigma} = \frac{A}{N_A\,\rho\,\sigma} $$

and the distance to the next scattering event is drawn from an exponential distribution with this mean. Here $A$ is the atomic weight, $N_A$ Avogadro's number and $\rho$ the density. A 15 keV electron in silicon is deflected on average every 18 nm; in tungsten, every 1.6 nm.

Energy loss comes from Bethe's stopping power, the energy an electron loses per unit path length.

$$ \frac{dE}{ds} = -78500\,\frac{\rho\,Z}{A\,E}\,\ln\!\left(\frac{1.166\,E}{J}\right) \quad [\mathrm{keV/cm}] $$

$J$ is the mean ionization potential, 172 eV for silicon. As $E$ falls the logarithm's argument drops below one and the sign flips, so low energies call for the Joy-Luo modification, which replaces the argument with $1.166(E+0.85J)/J$. The calculations below use it. A 15 keV electron then loses about 2.8 keV per micrometre of silicon. Spending all 15 keV takes a path of micrometres, and along the way the electron is deflected about 320 times on average.

## 3. How do random numbers build an electron trajectory?

Two equations suffice to generate one electron's trajectory. No analytic solution is needed; random numbers and step-by-step tracking do the work. The rule has only four steps.

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

That is all four steps. The loop repeats until the electron crosses back out of the surface or its energy falls below 0.5 keV. Crossing back out counts it as a backscattered electron.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/en/fig1-interaction-volume.png" alt="Monte Carlo trajectories and interaction volumes for 5 kV and 30 kV electrons in silicon and tungsten" width="820">
_Fig 2. Interaction volume by material and accelerating voltage. Each panel shows 220 trajectories: blue ones stop inside the specimen, red ones return through the surface as backscattered electrons. The dashed semicircle is the Kanaya-Okayama range, and η is the backscatter coefficient from a separate run of 3,000 electrons. The panels use different scales; raising the voltage enlarges the volume, and in tungsten it becomes a hemisphere pressed against the surface._

The four panels display the two factors from section 2 directly.

Raising the accelerating voltage from 5 kV to 30 kV (left to right) inflates the whole volume. In silicon the penetration depth grows from 0.47 µm to 9.31 µm, a factor of 20. Higher energy means a smaller cross section, less deflection, and more energy to spend before stopping.

Raising the atomic number (top to bottom) shrinks the volume and reshapes it. At 30 kV silicon reaches 9.31 µm against 1.68 µm for tungsten, a fifth as far. Most of that difference, though, comes from density (2.33 against 19.25 g/cm³, a factor of 8.3). Density enters the mean free path and the stopping power alike and so shrinks every length as $1/\rho$; give tungsten silicon's density and its maximum depth is only 16 % shallower than silicon's. What atomic number itself changes is the shape more than the size. Silicon's volume is a teardrop that swells sideways below the surface, while tungsten's is closer to a hemisphere pressed against it. In silicon an electron travels some distance before losing its direction; in tungsten it deflects sharply on entry. That is the $Z^2$ at work.

The dashed semicircle marks the Kanaya-Okayama penetration depth.

$$ R_{KO} = \frac{0.0276\,A\,E_0^{1.67}}{Z^{0.89}\,\rho} \quad [\mu\mathrm{m}] $$

Here $E_0$ is the incident energy in keV. Most trajectories crowd into the inner half of that arc. $R_{KO}$ estimates the farthest distance reached by the few electrons that travel furthest; it is not a typical value. In silicon at 15 kV, 90 % of the energy is deposited within 1.55 µm, 53 % of $R_{KO}$ (2.93 µm); at 30 kV the figure is 5.45 µm, or 59 %. So 90 % of the energy, and most of the signal, comes from little more than half of $R_{KO}$.

## 4. From what depth does each signal come?

So the interaction volume spans micrometres. How then does an SEM resolve nanometres?

Because not all of the signal generated in that volume is used. Each signal escapes from a different depth, and the differences span three orders of magnitude.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/en/fig2-depth-distribution.png" alt="Energy deposition, X-ray generation and backscattered electron depth distributions in silicon at 15 kV, with the secondary electron escape depth" width="820">
_Fig 3. Depth distributions by signal in silicon at 15 kV, from 1,200 electrons. On the left, the blue solid curve is energy deposition, the green dashed curve characteristic X-ray generation, the red dotted curve the maximum depth reached by backscattered electrons before turning back, and the orange band the secondary electron escape depth (5 nm). The right panel magnifies the top 20 nm; at the left panel's scale the escape layer is a single line._

The blue solid curve on the left is energy deposited against depth. Secondary electrons arise from the energy handed over in inelastic scattering, so their generation follows this curve. It peaks near 0.5 µm and extends to 2.5 µm.

The green dashed curve keeps only those segments where the electron still carries more than 1.84 keV, the ionization energy of the silicon K shell, and can therefore produce a characteristic X-ray. It nearly coincides with the blue curve. An electron entering at 15 keV holds far more than 1.84 keV over most of its path, so X-rays emerge from essentially the entire interaction volume. This is why energy dispersive X-ray spectroscopy (EDS) is limited to micrometre spatial resolution when measuring composition. Turned around, the way to sharpen EDS is to lower the voltage. In silicon, 90 % of the energy deposited by electrons above 1.84 keV lies within 1.51 µm at 15 kV, 0.74 µm at 10 kV and 0.20 µm at 5 kV. There is a price, of course: as the voltage approaches the ionization energy, the stretch of path that can make X-rays shortens and the signal weakens.

The red dotted curve shows the maximum depth reached by electrons that eventually returned. It peaks between 0.3 and 0.5 µm and rarely exceeds 1 µm. Backscattered electrons sample only the upper third or so of the interaction volume.

The orange band marks the escape depth of secondary electrons. Their energy is very low, below 50 eV, so they scatter again and stop almost immediately inside the specimen. In metals and semiconductors only those generated within a few nanometres of the surface get out (in insulators, a few tens of nanometres). The representative value of 5 nm is about 0.2 % of the 2.93 µm penetration depth under these conditions, roughly one part in 600. That is why the band appears as a single vertical line on the left, and why the top 20 nm is magnified separately on the right.

One more signal is missing from Fig 3: the Auger electron. When an incident electron ejects an inner-shell electron, an outer-shell electron drops in to fill the vacancy. The leftover energy sometimes kicks out another orbital electron instead of leaving as an X-ray (Fig 4). That ejected electron is the Auger electron, and its energy is element-specific, so it carries compositional information. Its escape depth is a few nanometres or less, and around 1 nm for low-energy Auger lines (the mean free path of electrons of tens to hundreds of eV can be read from Fig 1 of [XPS Fundamentals 2](/en/posts/xps-surface-sensitivity-imfp/)). Auger electron spectroscopy (AES) rests on exactly this, which is also why the slightest surface contamination buries the signal. Attacking the same class of problem with X-rays instead of electrons gives [X-ray photoelectron spectroscopy (XPS)](/en/posts/xps-photoemission-binding-energy/).
<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/ext-auger-process.png" alt="Schematic of the Auger process in terms of atomic orbitals and energy levels" width="640">
_Fig 4. The Auger process. (a) An incident electron (red arrow) ejects a 1s electron; a 2s electron drops into the vacancy, and the energy released ejects a 2p electron as the Auger electron. (b) The same process on an energy-level diagram: an L₁ electron fills the K-shell vacancy and an L₂,₃ electron leaves above the vacuum level (Vac). (Source: A. Carlson, [Auger Process](https://commons.wikimedia.org/wiki/File:Auger_Process.svg), public domain)_

That resolves the opening contradiction. The interaction volume spans micrometres, but the information secondary electrons carry comes from a few nanometres of surface. Secondary electrons generated in the narrow column directly beneath the entry point, before the beam spreads sideways (SE1), give a resolution close to the beam diameter. Backscattered electrons leaving the surface also generate secondaries in the same shallow layer (SE2), but these are spread over the wide area the backscattered electrons exit from and form a background without resolution. Tracking 3,000 electrons in silicon at 15 kV with the same Monte Carlo as Fig 5, backscattered electrons leave the surface at a median distance of about 0.8 µm from the entry point, and 90 % of them within 1.4 µm. SE1 thus comes from a spot about the size of the beam, SE2 from a disc about 3 µm across. In yield, SE2 is about $\beta\eta$ times SE1; Seiler gives $\beta$ = 2–3, and $\eta$ for silicon is 0.17, so SE2 is 30–50 % of SE1. An SEM shows topography at nanometre scale not by defeating the interaction volume but by viewing only a very thin shell of it.

## 5. Why do heavy elements look bright in a backscattered electron image?

If secondary electrons show surface shape, what do backscattered electrons show? The $Z^2$ from section 2 returns here.

The fraction of incident electrons that return is the backscatter coefficient $\eta$. Since large-angle scattering scales with $Z^2$, $\eta$ should climb with atomic number. A Monte Carlo run on six elements at 20 kV can check that.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/en/fig3-backscatter-vs-z.png" alt="Backscatter coefficient against atomic number from Monte Carlo compared with the Reuter empirical fit" width="750">
_Fig 5. Backscatter coefficient against atomic number. Points are Monte Carlo results at 20 kV with 3,000 electrons per element; the grey dashed line is Reuter's empirical fit. The two agree for C and Si; from Cu upward the Monte Carlo points rise above the fit, reaching about 12 % high at W._

The values rise monotonically from 0.06 for carbon ($Z=6$) to 0.54 for gold ($Z=79$), steeply among light elements and more gently toward heavy ones. The grey dashed line is Reuter's empirical fit to measurements.

This curve underwrites compositional contrast. Place two phases of differing atomic number side by side and the number of returning electrons differs, so a BSE detector renders the heavier one brighter (Fig 6 is a colourised image, so the phases separate by colour rather than brightness). Among light and medium elements the slope is steep, so even a small difference in atomic number gives comparatively large contrast. Where the curve flattens, among heavy elements, the same difference yields much weaker contrast.

<img src="/assets/img/posts/electron-microscopy-beam-specimen-interaction/ext-gold-arsenopyrite-bse.jpg" alt="Backscattered electron image of native gold grown on an arsenopyrite crystal" width="600">
_Fig 6. Backscattered electron image of native gold grown on an arsenopyrite (FeAsS) crystal; the yellow areas are gold, the blue-grey ones arsenopyrite. Gold ($Z=79$) and arsenopyrite, whose mean atomic number is about 27, differ strongly in backscatter coefficient, so the two minerals separate cleanly. The colours were added by the image's creator. (Source: Boris Lobastov, [Gold on arsenopyrite SEM image](https://commons.wikimedia.org/wiki/File:Gold_on_arsenopyrite_SEM_image.png), [CC BY 4.0](https://creativecommons.org/licenses/by/4.0))_

Fig 5 also holds a discrepancy worth naming. At high atomic number the Monte Carlo points sit systematically above the empirical fit. Silicon agrees well, 0.168 against 0.164, but tungsten runs about 12 % high, 0.532 against 0.476. This is not statistical fluctuation: the difference (0.056) is six times the standard deviation (about 0.009). The first suspect is the cross section. The screened Rutherford cross section is a first Born approximation, which treats the nuclear charge as a weak perturbation, with orbital screening added as a single empirical constant, and the approximation degrades as the nuclear charge grows. This calculation alone, however, cannot tell whether the cause is the cross section or another piece, such as the stopping-power formula or stopping the tracking at 0.5 keV. Checking the cross section would mean rerunning with tabulated Mott cross sections from relativistic wave calculations.

Leaving the discrepancy visible seems the more honest choice. A Monte Carlo of this kind demonstrates trends, not measurements. It is enough to understand why the interaction volume takes the shape it does and which depth each signal comes from, but using it as a correction factor for quantitative analysis would mean validating each piece of the model, starting with the cross section.

## Summary and what comes next

What happens to an electron inside a specimen splits two ways. Elastic scattering off nuclei changes direction; inelastic scattering off orbital electrons removes energy. Almost every property of the interaction volume's shape follows from the fact that elastic scattering, large-angle scattering in particular, goes as $Z^2/E^2$. Raising the accelerating voltage enlarges the volume; elements of high atomic number are usually dense as well, which shrinks it, and the atomic number itself turns the teardrop into a hemisphere.

That the volume spans micrometres while an SEM resolves nanometres comes down to escape depth. Characteristic X-rays emerge from nearly the whole volume, backscattered electrons from the upper third, secondary electrons from within 5 nm of the surface (a few tens of nanometres in insulators), and Auger electrons from a few nanometres or less. The 5 nm for secondary electrons is about 0.2 % of the penetration depth. The same beam striking the same point is read at depths differing by three orders of magnitude depending on which detector receives it.

The next post takes up those detectors. It builds up how scanning forms an image, then shows how the position and collection angle of secondary and backscattered electron detectors inside the column produce topographic and compositional contrast respectively. It also covers the edge effect, where secondary electron yield rises on tilted surfaces and brightens outlines, along with the charging problem encountered on insulators and how low-voltage SEM answers it.

## References

- J. I. Goldstein et al., *Scanning Electron Microscopy and X-Ray Microanalysis*, 4th ed., Springer, 2018, ch. 1–5 (beam-specimen interaction, interaction volume, backscatter coefficient).
- D. C. Joy, *Monte Carlo Modeling for Electron Microscopy and Microanalysis*, Oxford University Press, 1995 (the standard single-scattering implementation).
- D. C. Joy and S. Luo, "An empirical stopping power relationship for low-energy electrons," *Scanning*, vol. 11, pp. 176–180, 1989 (low-energy Bethe modification).
- K. Kanaya and S. Okayama, "Penetration and energy-loss theory of electrons in solid targets," *Journal of Physics D*, vol. 5, pp. 43–58, 1972 (penetration depth fit).
- W. Reuter, in *Proc. 6th Int. Conf. on X-Ray Optics and Microanalysis*, G. Shinoda, K. Kohra, T. Ichinokawa (eds.), University of Tokyo Press, 1972, p. 121 (backscatter coefficient fit, Goldstein 4th ed. eq. 2.2).
- L. Reimer, *Scanning Electron Microscopy: Physics of Image Formation and Microanalysis*, 2nd ed., Springer, 1998, ch. 3–4 (scattering theory and signal generation depth).
- H. Seiler, "Secondary electron emission," in *Electron Beam Interactions with Solids*, SEM Inc., 1984, pp. 33–42 (50 eV secondary-electron convention, escape depth T ≈ 5λ).
- [Scanning Electron Microscopy, Carleton SERC](https://serc.carleton.edu/research_education/geochemsheets/techniques/SEM.html) (qualitative overview of SEM signals).

---
title: "XPS Fundamentals 3 — Spectral Structure: Spin-Orbit, Satellites, Auger"
lang: en
lang-exclusive: ["en"]
permalink: /posts/xps-spectrum-structure-satellites/
page_id: xps-spectrum-structure-satellites
date: 2026-10-16 20:00:00 +0900
categories: [Surface Analysis, Photoelectron Spectroscopy]
tags: [xps, spin-orbit, shake-up, auger-parameter, lineshape]
description: "One element does not make just one peak. Where each structure in an XPS spectrum comes from, from spin-orbit splitting to satellites and Auger peaks."
math: true
---

[Part 2](/en/posts/xps-surface-sensitivity-imfp/) started from the fact that only electrons escaping without loss contribute to a peak area, and worked out the information depth from there. Along the way it treated the spectrum itself simply. So did the picture set up in [Part 1](/en/posts/xps-photoemission-binding-energy/): one core level gives one peak, and a change of chemical state shifts that peak slightly.

Real high-resolution spectra do not look like that. A single element produces several peaks, and many of them have nothing to do with chemical state. Some are pairs dictated by atomic quantum mechanics. Some are traces of another electron being excited at the moment of photoionization. Some are not photoelectrons at all but electrons from a different process.

The distinction matters because of peak fitting in Part 4. Fitting splits a spectrum into components, but the data do not tell you what counts as a "component." Count a satellite or a multiplet as a separate chemical species and a compound that does not exist lands in the composition table, and all the quantification after that goes wrong. This post sorts out where each structure in a spectrum comes from, to give criteria for what should and should not be counted as a component.

## 1. Spin-Orbit Splitting — A Split Unrelated to Chemistry

Why does a Ti 2p spectrum always have two peaks? Not because titanium is in two chemical states. Whether it is metallic titanium or TiO₂, the 2p level always shows up as two peaks.

The cause is spin-orbit coupling of the core hole. When an electron is removed from a level with orbital angular momentum $l$, the spin and orbital angular momentum of the remaining hole couple into two states with total angular momentum $j = l \pm 1/2$. The two states have different energies, so the peak splits in two. This is spin-orbit splitting. An s level ($l = 0$) has no orbital angular momentum and does not split, which is why O 1s and C 1s show no spin-orbit doublet.

The area ratio of the two peaks comes from the degeneracy $2j+1$ of each state. For a p level ($l = 1$), $j = 1/2$ and $3/2$ have degeneracies 2 and 4, giving an area ratio of 1:2. For d levels, $j = 3/2 : 5/2$ = 4:6 = 2:3; for f levels, $j = 5/2 : 7/2$ = 6:8 = 3:4. These are ratios of **areas**, not peak heights. If the two components have different widths, the height ratio departs from these values.

<img src="/assets/img/posts/xps-spectrum-structure-satellites/en/fig1-spin-orbit-doublet.png" alt="Synthetic TiO2-like Ti 2p doublet and d-level doublet with components and area ratios" width="720">
_Figure 1. Spin-orbit doublets (synthetic spectra). (a) TiO₂-like Ti 2p: 2p₃/₂ at 458.5 eV, splitting 5.7 eV, with 2p₁/₂ made broader. (b) A d-level doublet with an arbitrary splitting._

The splitting grows with atomic number. According to Thermo Fisher's XPS reference data, Al 2p splits by only 0.44 eV, enough to make the metal peak look slightly asymmetric but not to separate it; Ti 2p splits by about 6 eV, Cu 2p by 19.75 eV and Zn 2p by about 23 eV. Within a given element, levels closer to the nucleus split more.

In practice the area ratio and the splitting are imposed as constraints in fitting. Leave both free and the fit matches well, but often at the cost of a physically meaningless solution. Still, these values should not be treated as immutable constants. Three caveats apply.

First, the splitting varies slightly with chemical state. In the same reference data, the Ti 2p splitting is 6.1 eV for the metal and 5.7 eV for TiO₂. Second, the two components need not have the same width. A Ti 2p₁/₂ hole decays quickly through a Coster–Kronig transition, in which the hole moves to 2p₃/₂ within the same shell and the energy is released. A shorter lifetime means a broader natural linewidth, so 2p₁/₂ is noticeably broader than 2p₃/₂. Figure 1(a) reflects this by giving 2p₁/₂ four times the Lorentzian width. Tying the two widths together is wrong here. Third, measured area ratios also deviate slightly from $2j+1$. Thermo Fisher gives the intensity ratio of the two Cu 2p components as 0.508 rather than 0.5. The two numbers in Figure 1(a) illustrate a related problem. The area ratio was set to 0.50, but integrating within the plotted energy window gives 0.48, because the Lorentzian tail of the broad 2p₁/₂ extends beyond the window. An area depends on where the integration range is cut, a problem that returns in Section 4 and in Part 4.

## 2. Initial State or Final State — A Criterion for Classifying Structure

Which structures in a spectrum tell us about the sample's chemistry, and which are created by the measurement process? Answering that means going back to the equation in Section 3 of Part 1.

$$ E_B = E_f(N-1) - E_i(N) $$

The binding energy is the energy difference between the $N$-electron system before photoionization and the $(N-1)$-electron system after. Part 1 drew two strands from this equation. Initial-state effects change $E_i(N)$, the charge state of the atom before photoionization, and are the main cause of chemical shifts. Final-state effects come from $E_f(N-1)$, the state the system is left in after the core hole forms.

The key fact is that there need not be only one final state. Right after the core hole forms, the system does not only go to its lowest-energy state. Several final states with slightly different energies are possible: one in which a valence electron is excited at the same time, others in which the hole's spin pairs with valence spins in different ways. If $E_f(N-1)$ can take several values, photoelectrons from the same core level emerge with several kinetic energies and the peak splits into several. That is the origin of satellite structure.

In this framework, the spin-orbit splitting of Section 1 is also a final-state effect: the hole has two final states, $j = 1/2$ and $3/2$. Being fixed for each atom, though, it carries almost no chemical information. The structures in Sections 3 to 5 are all variations on the same principle: one initial state with several final states gives several peaks.

## 3. Satellite Structure — Shake-up, Shake-off, Multiplet Splitting

What is the bump beside the main peak? The most common one is a shake-up satellite. At the moment a core electron leaves, the atom's potential changes abruptly, and that jolt can excite a valence electron into an empty higher level. The photoelectron's kinetic energy is reduced by the energy spent on that excitation, so a shake-up satellite appears on the **high-binding-energy side** of the main peak.

Two examples are typical. One is the $\pi \to \pi^*$ satellite in the C 1s of aromatic polymers: π electrons of the benzene ring are excited into π* orbitals, producing a small peak a few eV from the main line. The other is copper, where the satellite becomes the decisive evidence for the chemical state.

<img src="/assets/img/posts/xps-spectrum-structure-satellites/en/fig2-satellites.png" alt="Synthetic Cu 2p3/2 spectra for Cu(0)-like and Cu(II)-like cases; only Cu(II) has shake-up satellites" width="720">
_Figure 2. Presence of satellites in the Cu 2p₃/₂ region (synthetic spectra). Main-peak positions follow the measurements of Biesinger et al. (2010), and the satellite position follows Thermo Fisher's data (around 943 eV). Peak widths are approximate, and the satellite shape and intensity are simulated._

The values compiled by Biesinger and co-workers put Cu metal 2p₃/₂ at 932.63 eV and Cu₂O at 932.18 eV, less than 0.5 eV apart. That is on the order of the instrumental resolution, so binding energy alone can hardly separate them. A Cu(II) compound such as CuO, however, not only moves its main peak up to 933.76 eV and broadens it to 3 eV, but also carries a distinct shake-up satellite around 940–945 eV. Cu(0) and Cu(I) have no such satellite, or only a very weak one. The difference comes from the electron configuration. Cu(II) is d⁹, with a vacancy in the 3d shell that leaves room for a valence excitation; Cu(0) and Cu(I) are d¹⁰, with a full 3d shell. For copper, **whether a satellite is present is itself evidence of the oxidation state**.

The satellite is not a separate chemical species. It is simply where part of the Cu(II) photoelectrons went after losing energy. When counting Cu(II), the satellite area must therefore be added to the main peak; counting only the main peak undercounts Cu(II). Biesinger and co-workers describe a method that takes the satellite-to-main-peak area ratio from a Cu(II) standard and uses it to compute the Cu(II) fraction in samples mixing Cu(0) or Cu(I) with Cu(II).

Shake-up has two relatives. In shake-off, the valence electron is not just excited but leaves the atom altogether. Since that electron carries off a continuous range of kinetic energy, the result is a broad background rather than a distinct peak. In shake-down, the satellite appears on the **lower**-binding-energy side of the main peak. In the 3d spectra of lanthanum and cerium compounds, a final state is possible in which electrons from the surrounding ligands move into the empty 4f level and screen the core hole better. That state is lower in energy and gives a peak at lower binding energy.

Multiplet splitting is somewhat different. When a core hole forms in an atom with unpaired electrons in its valence shell, the hole's spin and the unpaired valence spins can couple in several ways. Each coupling has a different energy, so one peak splits into several closely spaced components. The effect is prominent in the 2p levels of first-row transition metals with many unpaired d electrons (Cr, Mn, Fe, Co, Ni) and in rare earths. Gupta and Sen (1974) calculated this structure systematically from atomic theory. Conversely, a state with no unpaired electrons, such as V(V), shows no multiplet splitting.

Multiplets are what practitioners most need to watch for. Fit Fe 2p or Mn 2p without knowing about multiplet splitting and each split component gets read as a different chemical species. A sample with a single oxidation state then yields three or four "species." For such elements one should use a model in which the positions and ratios of the multiplet components have been fixed from standard samples, as Biesinger et al. (2011) did.

## 4. Plasmon Losses and the Asymmetric Line Shape of Metals

Why do metal peaks have a long tail on one side? A metal holds a sea of conduction electrons that are not bound to any atom. A photoelectron crossing that sea loses energy in two ways.

The first is plasmon loss. A plasmon is a collective oscillation of the whole conduction-electron gas, and its energy is fixed by the electron density. Each time a photoelectron excites a plasmon, it loses exactly that amount of energy, so loss peaks line up at regular intervals on the high-binding-energy side of the main peak. They are especially clear in aluminum, which is close to a free-electron metal. In Al 2p the bulk plasmon appears 15.3 eV from the main peak and the surface plasmon 10.4 eV away, with several orders of the bulk plasmon lined up at integer multiples. Thermo Fisher's data warn that these loss peaks are strong enough to obscure other peaks at higher binding energy, such as Si 2p and Si 2s.

The second is an asymmetry of the line shape itself. In a metal, lifting an electron from just below the Fermi level to an empty level just above it costs almost no energy. Such small excitations can occur continuously when the core hole forms, so the photoelectron can lose any amount of energy starting from nearly zero. As a result the main peak itself becomes asymmetric, trailing toward higher binding energy. Doniach and Šunjić (1970) derived this shape theoretically, and it is usually called the DS line shape.

```python
def doniach_sunjic(be, e0, alpha, f):
    """Doniach–Šunjić line shape. alpha: asymmetry, f: lifetime width (HWHM)."""
    u = -(be - e0)          # the original form is in kinetic energy, so flip for the BE axis
    return (np.cos(np.pi * alpha / 2 + (1 - alpha) * np.arctan(u / f))
            / (f**2 + u**2) ** ((1 - alpha) / 2))
```

With $\alpha = 0$ this reduces to a Lorentzian; larger $\alpha$ gives a longer tail. The observed line shape is this convolved with a Gaussian representing instrumental resolution. Figure 3 shows a metal peak built by convolving a DS line shape with $\alpha$ = 0.12 with a Gaussian. For this peak, 68% of the area lies on the high-binding-energy side of the maximum; a symmetric peak would have 50%.

<img src="/assets/img/posts/xps-spectrum-structure-satellites/en/fig3-asymmetric-lineshape.png" alt="Asymmetric metal peak compared with a symmetric Voigt, and residuals and a spurious component when fitting with symmetric functions" width="760">
_Figure 3. Fitting a metal's asymmetric line with symmetric functions (synthetic spectra). (a) DS⊗Gauss compared with a symmetric Voigt. (b) Fit with one symmetric Voigt, with residual. (c) Fitting with two symmetric Voigts creates a spurious component at higher binding energy._

Insulator and oxide peaks, lacking conduction electrons, are close to symmetric; metal peaks are asymmetric. This difference produces one of the most common failures in the peak fitting of Part 4. Fit an asymmetric metal peak with a single symmetric Voigt, as in Figure 3(b), and a clear structure remains in the residual. Add one more component to remove it and you get Figure 3(c): the fitting program fills the tail with a second component 1.4 eV above the main peak holding 43% of the area. A shift of 1.4 eV sits within the typical range of chemical shifts seen in Part 1 (0.5–4 eV), exactly where one would be tempted to read an oxide. A pure metal peak that was never oxidized ends up reporting "43% oxide." Moreover, even with two components the root-mean-square residual does not drop to the level of the added noise. A function with the wrong shape does not fit no matter how many copies are added.

The DS line shape has one more awkward property. For $\alpha > 0$ the tail falls off slowly enough that the area diverges. Because the area depends on where the integration range is cut, quantification with asymmetric line shapes requires fixing the range and the background together. The out-of-window tail problem seen for Ti 2p₁/₂ in Section 1 is worse here.

Plasmon losses, asymmetric tails and the stepped background from Part 2 are not separate phenomena. All of them are photoelectrons losing energy in the solid. A fixed amount of loss (a plasmon) gives a peak; a continuous loss starting from zero gives the tail of a peak; many large losses give the background.

## 5. Why Auger Peaks Appear in an XPS Spectrum

Why are Auger peaks mixed into a photoelectron spectrum? The survey spectrum in Part 1 pointed out O KLL and C KLL in passing; now it is time to see what they are.

The core hole created by photoionization does not last long. An electron from an outer level soon drops into it, and the energy difference between the two levels has to go somewhere. There are two routes. One is fluorescence, emission of an X-ray photon. The other is the Auger process, which hands that energy to another electron in the atom and ejects it. The lighter the element, the more the Auger process dominates, and fluorescence gains share as atomic number increases. For the elements XPS usually deals with, many Auger electrons are emitted, so Auger peaks appear in the spectrum alongside the photoelectron peaks.

Auger peaks are named after the three levels involved, in X-ray notation (K = 1s, L = 2s and 2p, M = 3s, 3p and 3d, and so on). KLL is a process with a hole in the K shell, an L-shell electron filling it and another L-shell electron being ejected. Subshells are sometimes included, as in L₃M₄₅M₄₅ for copper. With binding energies $E_1, E_2, E_3$ for the three levels, the kinetic energy of the Auger electron is roughly

$$ E_K(\text{Auger}) \approx E_1 - E_2 - E_3 $$

In reality the interaction between the two holes and relaxation shift it from this value, but the key point is that the photon energy $h\nu$ does not appear. The Auger kinetic energy is set only by level spacings within the atom, regardless of which source created the initial hole. The photoelectron kinetic energy $E_K = h\nu - E_B - \phi_{sp}$, by contrast, moves with $h\nu$.

This difference makes the two separable. Plotted on a binding-energy axis, photoelectron peaks stay in place when the source changes. Auger peaks have fixed kinetic energies, so their apparent binding-energy position moves by the change in $h\nu$. Switching from Al K$\alpha$ (1486.6 eV) to Mg K$\alpha$ (1253.6 eV) shifts only the Auger peaks, by 233.0 eV on the binding-energy axis. "The Auger kinetic energy does not depend on $h\nu$" and "Auger peaks move on the binding-energy axis when the source changes" are the same fact seen on two different axes. This is the use of the dual anode anticipated in Section 2 of Part 1. When a photoelectron peak and an Auger peak overlap, taking a second spectrum with the other anode moves only the Auger peak out of the way.

## 6. The Auger Parameter and the Wagner Plot

How do you distinguish chemical states for an element that barely shows a chemical shift? Zinc is the classic example. In the literature averages compiled by Biesinger and co-workers, Zn 2p₃/₂ sits at 1021.62 eV for the metal and 1021.96 eV for ZnO, only 0.34 eV apart. They are close enough that different sources report the metal and oxide in opposite order. Binding energy alone can hardly say whether the zinc is oxidized.

The Auger peak is another matter. The zinc L₃M₄₅M₄₅ transition involves the 3d level, close to the valence region, twice. With two holes left in the final state, it responds much more sensitively to how well the surroundings screen those holes. So a quantity combining the photoelectron and Auger peaks is defined, the modified Auger parameter:

$$ \alpha' = E_B(\text{photoelectron}) + E_K(\text{Auger}) $$

The original Auger parameter $\alpha$ was defined as the kinetic-energy difference between the two peaks, and $\alpha' = \alpha + h\nu$. $\alpha'$ is powerful in two ways. First, it is independent of the source, since neither $E_B$ nor the Auger $E_K$ depends on $h\nu$. More importantly, it is independent of charging. If an insulating surface charges to a potential $V$, the apparent binding energy of every photoelectron peak rises by $V$ and the kinetic energy of every Auger peak falls by $V$. The two shifts cancel exactly in the sum. Part 1 noted that the binding-energy scale of an insulator is unsettled by charging; $\alpha'$ is free of that. The debate over charge correction is taken up in Part 4.

Applied to zinc, the difference is plain. $\alpha'$ is 2013.85 eV for the metal and 2010.14 eV for ZnO, 3.7 eV apart, more than ten times the binding-energy difference. Biesinger and co-workers likewise recommend $\alpha'$ for identifying zinc chemical states.

<img src="/assets/img/posts/xps-spectrum-structure-satellites/en/fig4-wagner-plot.png" alt="Wagner plots for Cu and Zn with literature averages and α' contours" width="720">
_Figure 4. Wagner plot. Points are the literature averages from Tables 7 and 10 of Biesinger et al. (2010). The horizontal (binding-energy) axis increases to the left by convention, and the gray diagonals are contours of constant $\alpha'$._

Spreading these two quantities over one plane gives the Wagner plot, also called a chemical state plot. The photoelectron binding energy goes on the horizontal axis and the Auger kinetic energy on the vertical axis, with one point per compound. By convention the horizontal axis runs with binding energy increasing to the left. Points of constant $\alpha'$ lie on $E_K = \alpha' - E_B$, which in this arrangement is a diagonal of slope +1. With the horizontal axis in the ordinary direction the same lines would appear with slope −1, so check the axis direction first when reading someone else's plot.

Copper in Figure 4(a) shows why both axes are needed. Cu₂O nearly coincides with the metal in binding energy but has an $\alpha'$ 2 eV lower, dropping two diagonals down. CuO, conversely, differs from the metal in $\alpha'$ by only 0.3 eV but has a binding energy 1 eV higher, moving it sideways. On either axis alone, two of the three chemical states overlap; on both axes together, all three separate.

The change in $\alpha'$ also has a physical interpretation. A photoelectron process leaves one hole, and an Auger process leaves two holes in the final state. If $R$ is the extra-atomic relaxation energy, the energy lowered by neighboring atoms screening a hole, a two-hole state gains roughly twice that. Add the approximation that the core levels involved all shift by the same amount through initial-state effects, and the initial-state contribution cancels in $\alpha'$, leaving

$$ \Delta\alpha' \approx 2\,\Delta R $$

So a shift in $\alpha'$ mainly measures a difference in final-state relaxation. Section 4 of Part 1 said that a binding-energy shift $\Delta E_B$ mixes initial- and final-state effects; separating out the final-state part with $\Delta\alpha'$ lets the remainder be read as the initial-state contribution. This separation holds only as far as the approximation above holds.

## 7. So Why Is the Background Stepped?

Why is the background always higher to the left of a peak, on the high-binding-energy side? All the ingredients for an answer are now in place. As Part 2 showed, photoelectrons that lose energy through inelastic scattering in the solid appear at lower kinetic energy, that is, at higher binding energy. As Section 4 showed, the energy lost varies continuously from small to large amounts. So every peak spreads energy-loss electrons across the whole range above it in binding energy, and the background steps up one level each time a peak is passed.

This picture is the physical basis for the background models of Part 4. The widely used Shirley background assumes that the background height at a given energy is proportional to the peak area on its lower-binding-energy side. It turns this section's logic, "the larger the peaks already passed, the more loss electrons lie behind them," into an equation.

Going one step further, the background can be treated as information rather than something to discard. The probability of losing energy depends on the distance a photoelectron has traveled in the solid; electrons from deeper down lose more. The shape of the background under a peak therefore encodes how the element is distributed in depth. Tougaard developed methods from this viewpoint that compute the background from inelastic-scattering cross sections and, conversely, estimate depth distributions from the background. That approach returns in Part 5 on thickness.

## Summary and Next Part

This post traced why one element produces several peaks. The criterion was Part 1's $E_B = E_f(N-1) - E_i(N)$: one initial state with several final states gives several peaks. Spin-orbit splitting is a pair of peaks from the two $j$ states of the hole, with an area ratio set by $2j+1$, though the splitting, the linewidths and the measured area ratio all vary a little. Shake-up satellites are final states with a valence electron excited at the same time, appear at higher binding energy, and in copper the presence or absence of a satellite settles the oxidation state. Multiplet splitting divides the peaks of elements with unpaired valence electrons into several components and readily invents species that are not there. Fitting a metal's asymmetric line with symmetric functions creates a spurious oxide component. Auger peaks can be told apart by the source independence of their kinetic energies, and the modified Auger parameter $\alpha'$, also independent of charging, separates chemical states that binding energy cannot. The step in the background is another face of the same physics of photoelectrons losing energy.

The next part takes this knowledge into quantification. Extracting composition from peak areas needs sensitivity factors; measuring an area requires deciding where to draw the background; splitting overlapping peaks requires choosing line shapes and constraints. For an insulator, the energy scale itself has to be rebuilt with charge correction first. Part 4 takes apart, one by one, how a single number in a composition table passes through all these models.

## References

- M. C. Biesinger, L. W. M. Lau, A. R. Gerson, and R. St. C. Smart, "Resolving surface chemical states in XPS analysis of first row transition metals, oxides and hydroxides: Sc, Ti, V, Cu and Zn," *Appl. Surf. Sci.* 257, 887–898 (2010). <https://doi.org/10.1016/j.apsusc.2010.07.086>
- M. C. Biesinger et al., "Resolving surface chemical states in XPS analysis of first row transition metals, oxides and hydroxides: Cr, Mn, Fe, Co and Ni," *Appl. Surf. Sci.* 257, 2717–2730 (2011). <https://doi.org/10.1016/j.apsusc.2010.10.051>
- R. P. Gupta and S. K. Sen, "Calculation of multiplet structure of core p-vacancy levels," *Phys. Rev. B* 10, 71–77 (1974). <https://doi.org/10.1103/PhysRevB.10.71>
- S. Doniach and M. Šunjić, "Many-electron singularity in X-ray photoemission and X-ray line spectra from metals," *J. Phys. C* 3, 285–291 (1970). <https://doi.org/10.1088/0022-3719/3/2/010>
- C. D. Wagner, "Chemical shifts of Auger lines, and the Auger parameter," *Faraday Discuss. Chem. Soc.* 60, 291 (1975). <https://doi.org/10.1039/DC9756000291>
- G. Moretti, "Auger parameter and Wagner plot in the characterization of chemical states by X-ray photoelectron spectroscopy: a review," *J. Electron Spectrosc. Relat. Phenom.* 95, 95–144 (1998). <https://doi.org/10.1016/S0368-2048(98)00249-7>
- M. Balal et al., "Intrinsic and extrinsic plasmons in the hard x-ray photoelectron spectra of nearly free electron metals," *Phys. Rev. B* 109, 205419 (2024). <https://arxiv.org/abs/2404.16620>
- Casa Software Ltd., "Peak Fitting in XPS," 2006. <https://mmrc.caltech.edu/XPS%20software/CASA/Books/peak_fitting_in_xps.pdf>
- Thermo Fisher Scientific, XPS Periodic Table (Cu, Ti, Al, Zn). <https://www.thermofisher.com/us/en/home/materials-science/learning-center/periodic-table.html>
- "Auger Peaks and the Auger Parameter," XPS Reference Pages. <https://www.xpsfitting.com/2012/08/auger-peaks-and-auger-parameter.html>
- "Workshop Exercises: The Auger Parameter and Wagner Plots," XPS Reference Pages. <https://www.xpsfitting.com/2016/09/workshop-exercises-auger-parameter-and.html>
- "Initial vs Final State Effects," HarwellXPS. <https://www.harwellxps.guru/xpskb/initial-vs-final-state-effects/>
- S. Hüfner, *Photoelectron Spectroscopy: Principles and Applications*, 3rd ed., Springer, 2003.
- D. J. Morgan, "X-Ray Photoelectron Spectroscopy (XPS): An Introduction," Cardiff Catalysis Institute. <https://sites.cardiff.ac.uk/xpsaccess/files/2014/07/AccessXPS_Primer_Paper.pdf>

---
title: "XPS Fundamentals 1 — If Only the Energy of Ejected Electrons Is Measured, How Do We Tell the Elements and Their Oxidation States?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/xps-photoemission-binding-energy/
page_id: xps-photoemission-binding-energy
date: 2026-10-02 20:00:00 +0900
categories: [Surface Analysis, Photoelectron Spectroscopy]
tags: [xps, photoemission, binding-energy, chemical-shift, surface-analysis]
description: "XPS measures only the kinetic energy of photoelectrons. How elements and chemical states are read from that, and why binding energy is not an orbital energy."
math: true
---

X-ray photoelectron spectroscopy (XPS), also known as ESCA (electron spectroscopy for chemical analysis), is one of the most widely used techniques in surface analysis. What a process engineer wants from it usually fits in one sentence: which elements are on this surface, in what proportions, and in what chemical states. Did the oxide grow as intended? How much carbon is left after cleaning? Has an unwanted compound formed at the interface? All of these are questions of that form.

What the instrument physically does is much simpler. It shines X-rays on the sample and counts the electrons that come out, sorted by energy (Figure 1). The output is a single curve of counts versus electron kinetic energy. Nothing on that curve says "30% oxygen" or "Si⁴⁺". Element, chemical state, composition and thickness are all obtained by fitting a model to the curve and working backwards. This series follows the gap between what is measured and what we want to know. The things to watch are where each model enters and where its assumptions break.

<img src="/assets/img/posts/xps-photoemission-binding-energy/ext-xps-system-schematic.gif" alt="Schematic of an XPS system: X-rays strike the sample and photoelectrons pass through a collection lens and energy analyser to a detector" width="700">
_Figure 1. The basic layout of an XPS instrument. The red cone on the left is the X-ray beam focused on an SiO₂/Si sample, and the blue arrows are the emitted photoelectrons. The electrons pass through a collection lens into an energy analyser that measures their kinetic energy, and a detector counts them to build the Si 2p spectrum at lower right. The depth (70–110 Å) and vacuum (below 10⁻⁸ torr) quoted in the drawing are the subject of [Part 2](/en/posts/xps-surface-sensitivity-imfp/); the energies labelled kV should be read as keV. (Source: Bvcrist, [System2](https://commons.wikimedia.org/wiki/File:System2.gif), public domain)_

The series has five parts, each answering one question that comes up in front of the instrument.

- Part 1 (this post): How does the energy of an ejected electron become a binding energy? Why does a peak of the same element move with its chemical state?
- [Part 2](/en/posts/xps-surface-sensitivity-imfp/): X-rays go micrometres deep, so why does XPS see only the top few nanometres? The key is the inelastic mean free path (IMFP) and the information depth.
- Part 3: Why does one element give several peaks? It covers spin-orbit splitting, satellite structure and Auger peaks.
- Part 4: Why do different analysts get different compositions from the same spectrum? It covers quantification from peak areas and peak fitting.
- Part 5: How far can you trust an oxide thickness measured by XPS? It covers thickness measurement, angle-resolved XPS, and the ways a measurement can alter the sample it is measuring.

## 1. What Can the Energy of an Ejected Electron Tell Us?

The instrument measures one thing: the kinetic energy of an electron. What can be drawn from that number? The starting point is the photoelectric effect. When an atom absorbs a photon of energy $h\nu$, one electron takes that energy and leaves the atom. Part of the photon energy pays for releasing the electron from its bound state, the binding energy $E_B$, and the remainder becomes kinetic energy. In a solid there is one more cost: the work function, the energy needed to carry the electron across the surface into vacuum. The basic equation of XPS is therefore

$$ E_K = h\nu - E_B - \phi_{sp} $$

The source fixes $h\nu$ and the spectrometer measures $E_K$. If $\phi_{sp}$ is known, $E_B$ follows in one line. Take a C 1s electron with a binding energy of about 284.8 eV. Ejected by Al K$\alpha$ (1486.6 eV), it arrives with $E_K + \phi_{sp}$ = 1201.8 eV; ejected by Mg K$\alpha$ (1253.6 eV), with 968.8 eV. The kinetic energy changes by 233.0 eV with the source, but the binding energy does not. The subscript deserves attention. $\phi_{sp}$ is the work function of the **spectrometer**, not of the sample.

<img src="/assets/img/posts/xps-photoemission-binding-energy/en/fig1-energy-diagram.png" alt="Energy-level diagram of sample and spectrometer, with aligned Fermi levels cancelling the sample work function" width="640">
_Figure 2. Energy levels in photoemission. The sample and spectrometer are in electrical contact, so their Fermi levels coincide. The spectrometer measures kinetic energy from its own vacuum level, so only the spectrometer work function, not the sample's, remains in the equation._

Figure 2 shows why the sample work function drops out. A conducting sample is electrically grounded to the spectrometer. When two metals are connected, electrons flow between them until their Fermi levels $E_F$ sit at the same height. A photoelectron leaving the sample has kinetic energy $h\nu - E_B - \phi_s$ relative to the sample's vacuum level. By the time it reaches the spectrometer entrance, the relevant reference is the spectrometer's vacuum level instead. The two vacuum levels differ by the contact potential $\phi_s - \phi_{sp}$, and the electron is accelerated or retarded by that amount. The kinetic energy the spectrometer sees is

$$ E_K = (h\nu - E_B - \phi_s) + (\phi_s - \phi_{sp}) = h\nu - E_B - \phi_{sp} $$

$\phi_s$ cancels exactly. There is no need to know the work function of every sample; only the fixed $\phi_{sp}$ of the instrument has to be calibrated. In practice, the energy axis is set by measuring a reference metal such as gold or silver whose binding energies are well known. The same equation also fixes the zero of $E_B$ at the Fermi level. The binding energy XPS reports is referenced to the Fermi level of the solid, not to the vacuum level as an atomic ionisation energy would be.

This cancellation holds only if the sample and spectrometer share a Fermi level. An insulator breaks that assumption. If electrons cannot flow in to replace the photoelectrons that left, positive charge builds up at the surface and every peak shifts to higher binding energy by the resulting potential. This is the charging problem, taken up again with charge correction in Part 4. In short, the XPS binding-energy scale is not a measurement. It is a scale built from one energy-conservation equation and one assumption about Fermi-level alignment.

## 2. Ultraviolet Light Also Ejects Electrons, So Why Use X-rays?

Ultraviolet light also produces photoelectrons, and ultraviolet photoelectron spectroscopy (UPS) is an established technique. Why go to X-rays when the goal is composition?

It depends on which electrons we want to remove. The electrons of an atom fall into two groups. The outer valence electrons take part in chemical bonding directly, and their energy levels shift a great deal depending on the bonding partner. In a solid they smear out into broad bands. The core-level electrons further in are tightly bound near the nucleus and barely participate in bonding. Their binding energies stay close to values set by the atomic number, which makes them a fingerprint for identifying elements. The difficulty is that these binding energies range from tens of eV to well over a thousand. The ultraviolet light from a helium discharge lamp carries only a few tens of eV, enough to reach the valence electrons and little else. Reaching the core levels takes photons in the keV range, which means X-rays.

A laboratory XPS source is an X-ray tube (Figure 3). An electron beam (typically accelerated through 15 kV) strikes a metal anode, knocks out 1s electrons of the anode atoms, and 2p electrons dropping into those vacancies emit characteristic X-rays. Aluminium and magnesium are the standard anode materials.

- Al K$\alpha$: $h\nu = 1486.6$ eV
- Mg K$\alpha$: $h\nu = 1253.6$ eV

<img src="/assets/img/posts/xps-photoemission-binding-energy/ext-xray-source-silver-target.jpg" alt="Inside an ultra-high-vacuum chamber, the nose of an X-ray source points down at a silver foil sample on a copper mount" width="600">
_Figure 3. The X-ray source and sample inside an ultra-high-vacuum XPS chamber. The cylinder at the top is the end of a water-cooled aluminium/magnesium dual-anode X-ray gun, and a silver foil sample is screwed onto the copper mount below. Here the source illuminates the sample directly, with no monochromator in between. (Source: MichaelBueker, [Silver Target in XPS Spectrometer cropped](https://commons.wikimedia.org/wiki/File:Silver_Target_in_XPS_Spectrometer_cropped.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

These lines were chosen because they combine sufficient energy with a narrow linewidth. Even so, the linewidth of unmonochromated Al K$\alpha$ (the FWHM of the unresolved K$\alpha_1$/K$\alpha_2$ doublet) is about 0.85 eV, too broad to resolve many chemical shifts. The tube also emits weak K$\alpha_{3,4}$ satellite lines next to the main line, together with bremsstrahlung, the continuous spectrum radiated by decelerating electrons. The satellites create spurious peaks, and the bremsstrahlung raises the background.

Most modern instruments therefore use a monochromator. It relies on Bragg diffraction from a quartz crystal to select a narrow slice of Al K$\alpha$ and focus it onto the sample. After monochromatization the linewidth drops to roughly 0.2–0.3 eV, and the satellites and bremsstrahlung are filtered out. The price is intensity. The crystal reflects only a narrow energy band, and it collects only a small solid angle of the X-rays leaving the tube, so most of the photons are discarded. In an early monochromated instrument (Kelly and Tyler, 1973), the K$\alpha$ flux reaching the sample was about one hundredth of the unmonochromated value. Cutting the linewidth to roughly a third cost two orders of magnitude in intensity. Modern instruments may differ in the ratio, but the trade between resolution and signal is the same, and it will come up repeatedly in later parts.

A dual anode with both aluminium and magnesium has no monochromator, but it is useful in another way. When the source energy changes, photoelectron peaks stay put on the binding-energy axis while Auger peaks move. Part 3 shows how this separates overlapping peaks. As a rule, use a monochromated source when the aim is to resolve chemical states, and switch dual-anode sources when an Auger peak overlaps a photoelectron peak.

## 3. Is the Binding Energy the Energy of the Orbital the Electron Came From?

At first sight the $E_B$ from Section 1 looks like the energy of the orbital the electron occupied. It is the energy spent removing the electron, so it should measure how deeply the electron was bound. Strictly, it is not, and several features of an XPS spectrum come from that difference.

The simplest picture is Koopmans' theorem. Under the frozen orbital approximation, which assumes that the remaining electrons' orbitals do not change at all after one electron is removed, the binding energy is the negative of the Hartree-Fock orbital energy $\varepsilon$.

$$ E_B^{K} = -\varepsilon_{\text{orbital}} $$

In other words, while one electron leaves, the others are assumed to hold still, like people in a photograph.

The binding energy, however, is defined not as an orbital energy but as a difference between the total energies of two states. With $E_i(N)$ the energy of the $N$-electron system before photoionisation and $E_f(N-1)$ the energy of the $(N-1)$-electron system afterwards,

$$ E_B = E_f(N-1) - E_i(N) $$

When a core electron leaves a real atom, it leaves behind a core hole, and one electron that was screening the nuclear charge is gone. The remaining electrons respond. Lift a ball off a mattress and the springs around it settle into the dip; in the same way, the remaining electrons feel the stronger pull of the nucleus and draw in towards the core hole to screen it. This is relaxation. It happens on essentially the same timescale as photoionisation and is reflected in the energy of the outgoing electron. Relaxation lowers the final-state energy $E_f(N-1)$ below the value computed with frozen orbitals. The measured binding energy is therefore **smaller** than the Koopmans prediction.

$$ E_B = E_B^{K} - E_R, \qquad E_R > 0 $$

$E_R$ is the relaxation energy. The sign is easy to get backwards, so it is worth restating. The more stable the final state, the less it costs to ionise the atom, and the surplus goes to the photoelectron. Kinetic energy goes up and binding energy goes down. Some introductory materials state this in the opposite direction, so care is needed. The Morgan primer in the references is one example.

The size of the relaxation energy is worth noting too. For the 1s level of a neon atom it is about 23 eV, an order of magnitude larger than the chemical shifts common in solids (0.5–4 eV, Section 4). A 10% change in relaxation energy with the surroundings is therefore enough to produce a shift the size of a chemical shift. Strictly, electron correlation and relativistic effects also enter as corrections. For light and medium elements, relaxation is the largest term: for Ne 1s the relativistic correction is about 1.1 eV and the correlation correction below 1 eV. The relativistic correction grows roughly as the fourth power of the atomic number, however, so for deep core levels of heavy elements it can exceed relaxation.

The conclusion from this equation is that two kinds of factors move the binding energy. One is the state of the atom before photoionisation: how much electron density it has lost or gained. This is the initial-state effect. The other is how effectively the surrounding electrons screen the core hole once it exists. This is the final-state effect. Final-state effects do not stop at relaxation. If a valence electron is additionally excited during relaxation, a small peak appears next to the main peak with some energy missing. That is the origin of the satellite structure covered in Part 3.

## 4. Why Does the Same Element Show Binding Energies a Few eV Apart?

Section 2 said that core levels barely participate in bonding. Yet the same core level of the same element appears at different binding energies in different chemical environments. "Barely" is doing real work in that sentence. Core electrons are not directly involved in bonds, but the valence charge redistributed by bonding changes the potential they feel.

When an atom bonds to a more electronegative neighbour, valence density is pulled towards that neighbour. The atom becomes partially positive, and the valence electrons screen less of its nuclear charge. The core electrons are bound more tightly, and the binding energy rises. An atom that gains electrons shows the opposite shift. This dependence of binding energy on chemical environment is the chemical shift. The simplest quantitative model is the charge potential model. It writes the shift for atom $i$ as the sum of a term from its own charge $q_i$ and the potential created by the surrounding atoms.

$$ \Delta E_{B,i} = k\,q_i + \sum_{j \neq i} \frac{q_j}{r_{ij}} $$

The first term is the effect of the atom's own charge. The second is a Madelung-type potential from neighbouring charges $q_j$ at distances $r_{ij}$. The second term often has the opposite sign to the first, so the two partly cancel. The model covers only the initial state; relaxation is not in it.

The textbook illustration of the chemical shift is ethyl trifluoroacetate ($\mathrm{CF_3COOCH_2CH_3}$). The molecule has four carbon atoms, each with very different neighbours. The methyl (CH₃) carbon is surrounded only by carbon and hydrogen. The CH₂ carbon of the ethyl group bonds to one oxygen, the carboxyl carbon to two, and the CF₃ carbon to three fluorine atoms, the most electronegative element. The C 1s spectrum of this molecule splits into four peaks, and their order matches the electronegativity of the neighbours exactly. Siegbahn's group chose it to demonstrate the chemical shift in their first ESCA monograph (1967), and it has been known as "the ESCA molecule" ever since.

<img src="/assets/img/posts/xps-photoemission-binding-energy/en/fig3-chemical-shift.png" alt="Component decomposition of an ethyl-trifluoroacetate-like C 1s spectrum" width="620">
_Figure 4. Simulated C 1s spectrum of the ethyl trifluoroacetate type. Component positions follow the values measured by Travnikova et al. (2012) (1.72, 4.33 and 7.46 eV relative to CH₃), and each component is a Gaussian with a FWHM of 1.0 eV. Grey dots are simulated data with noise; the black line is the sum of the four components. The more electronegative the neighbours, the further left, that is, the higher in binding energy, the peak sits._

In gas-phase measurements the four peaks spread over about 7.5 eV (about 8.3 eV in the 1973 measurement by Gelius et al.), with 1.7–3.1 eV between neighbouring peaks. Figure 4 is a spectrum synthesised with those spacings. The shift range is exceptionally large here because three fluorine atoms sit on one carbon. Chemical shifts commonly met in solid samples are typically 0.5–4 eV, and even large ones stay within about 10 eV. Silicon shows oxidation states most clearly. The Si 2p level sits at about 99.4 eV for elemental Si(0) and about 103.3 eV for Si(4+) in SiO₂, a shift of roughly 4 eV on oxidation. The intermediate oxidation states (Si¹⁺, Si²⁺, Si³⁺) fall in order between those two values.

Think of a binding energy as an address: the element is the city and the chemical state is the street number. In Figure 5, O 1s and C 1s alone are about 248 eV apart, while chemical shifts are 0.5–4 eV. A broad scan is therefore enough to tell which elements are present, but telling chemical states apart takes a second, high-resolution look at a narrow window. The two-stage measurement in Section 5 follows from this distinction.

Put a shift of 0.5–4 eV next to the linewidths in Section 2 and the problem is clear. Even with a monochromated source, the instrumental resolution (about 0.25 eV) combines with the natural width of the level to give peak widths of around 0.4–0.5 eV, the same order of magnitude as the shifts. Peaks from neighbouring chemical states often do not separate cleanly; they overlap into what looks like a single peak. Splitting such an envelope into components is peak fitting, and its difficulties and pitfalls are covered in Part 4.

One qualification is needed here. The chemical shift is mostly explained by initial-state effects, but not always. The relaxation energy $E_R$ from Section 3 also depends on chemical environment. Where free electrons are abundant, as in a metal, the core hole is screened well; in an insulator, less so. A measured shift $\Delta E_B$ contains both contributions, and binding energy alone cannot separate them. Part 3 shows how the Auger parameter, which combines a photoelectron peak with an Auger peak, pulls the two apart.

## 5. What Does a Real Spectrum Look Like, and What Should You Read First?

The horizontal axis is the first surprise. The instrument measures kinetic energy, yet the axis is labelled binding energy, and the numbers increase from right to left. Both conventions follow from the equation in Section 1. Since $E_K = h\nu - E_B - \phi_{sp}$, binding energy moves in the opposite direction to kinetic energy. Plotting in binding energy is convenient because, as the C 1s example in Section 1 showed, a given element's peaks land on the same numbers regardless of the source. Reversing the axis keeps the original measurement direction, with kinetic energy increasing to the right.

<img src="/assets/img/posts/xps-photoemission-binding-energy/en/fig2-survey-schematic.png" alt="Synthetic survey spectrum of a hypothetical SiO2/Si sample, showing photoelectron peaks, Auger peaks and a stepped background" width="680">
_Figure 5. Synthetic survey spectrum of a hypothetical SiO₂/Si sample (Al K$\alpha$). Blue labels mark photoelectron peaks and red labels Auger peaks. The background steps up on the high-binding-energy (left) side of each peak, and the small shoulder just right of Si 2p is Si(0) from the substrate under the oxide. Peak positions are approximate._

Measurements are usually done in two stages. A survey scan first covers a wide range, from 0 eV up to near the source energy. The analyser pass energy is set high, trading resolution for signal. The aim is a list of the elements present on the surface. Then high-resolution (narrow) scans cover a few tens of eV around the peaks of interest. The pass energy is lowered to improve resolution, and the lost signal is made up by averaging many repeated sweeps. Chemical-state analysis is done on these high-resolution scans. How pass energy trades resolution against sensitivity is covered with the analyser design in Part 4.

Peaks are named by element symbol, principal quantum number, orbital letter and total angular momentum $j$. For s orbitals such as O 1s and C 1s there is only one $j$, so it is omitted. p, d and f orbitals split into two levels through spin-orbit coupling, so $j$ is added, as in Ti 2p$_{3/2}$ or Cd 3d$_{5/2}$. Why the spacing and area ratio of the two levels take fixed values is covered in Part 3.

A survey spectrum contains more than photoelectron peaks. The O KLL and C KLL features marked in red in Figure 5 are Auger peaks. After a core hole forms, an electron from an outer level drops into it, and the energy released ejects another electron. The kinetic energy of that electron is set by the level spacings within the atom and does not depend on the source energy. Converting it to binding energy puts it on the same axis, but it is a different kind of peak from a photoelectron line. The broad Cd MNN series in a CdS sample arises the same way. Also notice that the background steps up to the left of each peak, on the high-binding-energy side. Electrons that lose energy through inelastic scattering on their way out of the sample pile up at lower kinetic energy, which is higher binding energy. Figure 5 was synthesised by attaching one such step to each peak. The small shoulder about 3.9 eV below Si 2p, just to its right, is Si(0) from the silicon substrate beneath the oxide. The oxidation shift of Section 4 shows up even in a survey.

```python
def peak_with_step(be, center, fwhm, area, step_frac):
    """Photoelectron peak plus a step background on the high-BE side."""
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    step = 0.5 * (1 + np.vectorize(erf)((be - center) / (np.sqrt(2) * sigma)))
    return gaussian(be, center, fwhm, area) + step_frac * area * step
```

The step height here is an arbitrary fraction. In a real spectrum, the shape chosen for this background when subtracting it changes the quantitative result by several percent. Auger peaks are covered in Part 3 and background models in Part 4. The full code is in [`_code/xps-photoemission-binding-energy/generate_figures.py`](https://github.com/eddykim/eddykim.github.io/blob/main/_code/xps-photoemission-binding-energy/generate_figures.py).

## 6. What Else Is Needed to Read Elements, States, Amounts and Depths from One Curve?

What is left if everything so far is reduced to one line? The instrument gives a single curve of counts versus kinetic energy $E_K$. We want four things. Which elements are present comes from peak positions; which chemical states, from chemical shifts; how much, from peak areas; and at what depth, from how much the signal is attenuated.

None of these four answers can be read directly off the curve. Converting kinetic energy to binding energy required the energy-conservation equation and a work-function calibration, plus a charge correction for insulators. Separating chemical states requires a line-shape model to split overlapping peaks, and a background model. Turning areas into composition requires sensitivity factors that differ by element and level. Talking about depth requires an attenuation model for how far electrons travel in a solid. A number in a report such as "x% SiO₂ component, oxide thickness y nm" is the output of all these models in sequence. If one model's assumption is wrong, the number still comes out with plausible decimal places, but it no longer means what it claims. The remaining four parts follow each model to the point where it breaks.

This is the same concern as in [Ellipsometry Foundations 1](/en/posts/ellipsometry-electromagnetic-fresnel/). There, the detector measured only intensity, so the phase of the electric field was lost, and comparing polarisation states recovered that information. Thickness and refractive index still had to be extracted by model fitting. In XPS too, the instrument measures only electron kinetic energy and models supply everything else. The physics differs, but the structure is the same: a model sits between what is measured and what we want to know.

## Summary and Next Part

This post traced how the XPS energy scale is built. The energy-conservation equation of the photoelectric effect, $E_K = h\nu - E_B - \phi_{sp}$, gives the binding energy, and the sample work function drops out because the sample and spectrometer Fermi levels are aligned. In an insulator that premise fails and charging appears. X-rays are used to reach the core levels that fingerprint each element, and a monochromator narrows the linewidth at the cost of intensity. The measured binding energy is not an orbital energy; because of relaxation, it is smaller than the Koopmans value. The factors that move the binding energy divide into initial-state and final-state effects. The chemical shift comes mainly from the former, with a contribution from the latter.

The [next part](/en/posts/xps-surface-sensitivity-imfp/) questions a premise this post took for granted. Measured by the distance over which their intensity falls to 1/e, Al K$\alpha$ X-rays penetrate about 8 µm into silicon, and more than 0.2 µm even into a heavy element such as gold. Why, then, is XPS called a technique that sees only the top few nanometres? The answer lies with the electrons, not the X-rays. Part 2 covers the distance a photoelectron can travel through a solid without losing energy, the inelastic mean free path (IMFP), and the information depth defined from it.

## References

- D. J. Morgan, "X-Ray Photoelectron Spectroscopy (XPS): An Introduction," Cardiff Catalysis Institute. <https://sites.cardiff.ac.uk/xpsaccess/files/2014/07/AccessXPS_Primer_Paper.pdf> (Page 1 states that the observed binding energy is larger than the Koopmans prediction, the opposite of the relaxation direction in Section 3 of this post.)
- P. S. Bagus, C. J. Nelin, and C. R. Brundle, "Chemical significance of x-ray photoelectron spectroscopy binding energy shifts: A Perspective," *J. Vac. Sci. Technol. A* 41, 068501 (2023). <https://pubs.aip.org/avs/jva/article/41/6/068501/2915551/Chemical-significance-of-x-ray-photoelectron>
- O. Travnikova et al., "The ESCA molecule—Historical remarks and new results," *J. Electron Spectrosc. Relat. Phenom.* 185, 191–197 (2012). <https://www.sciencedirect.com/science/article/abs/pii/S0368204812000552>
- L. Keller et al., *J. Chem. Phys.* 153, 114110 (2020), Table I (relaxation, relativistic and correlation corrections for Ne 1s). <https://arxiv.org/abs/2005.13465>
- M. A. Kelly and C. E. Tyler, *Hewlett-Packard Journal*, July 1973, p. 4 (linewidth and intensity of a quartz-crystal monochromator). <https://worldradiohistory.com/Archive-Company-Publications/HP-Journal/70s/HPJ-1973-07.pdf>
- J. F. Moulder et al., *Handbook of X-ray Photoelectron Spectroscopy*, Perkin-Elmer, 1992.
- S. Hüfner, *Photoelectron Spectroscopy: Principles and Applications*, 3rd ed., Springer, 2003, ch. 1–2.
- T. Koopmans, "Über die Zuordnung von Wellenfunktionen und Eigenwerten zu den einzelnen Elektronen eines Atoms," *Physica* 1, 104–113 (1934).

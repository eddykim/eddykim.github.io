---
title: "Ellipsometer Instrumentation 1 — Oblique Incidence and Rotating-Element Modulation"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-oblique-rotating-element/
page_id: ellipsometer-oblique-rotating-element
date: 2026-09-24 20:00:00 +0900
categories: [Optics, Instrumentation]
tags: [ellipsometry, instrumentation, rotating-analyzer, rotating-compensator, phase-modulation, stokes]
description: What separates ellipsometer designs is not what they measure but what they cannot. Null, RAE, RCE and PME sorted by the Stokes parameters each can reach.
math: true
---

The three posts of "Ellipsometry Foundations" were all forward calculations. [Post 1](/en/posts/ellipsometry-electromagnetic-fresnel/) set up the Fresnel reflection coefficients as complex numbers, [post 2](/en/posts/ellipsometry-polarization-mueller-matrix/) wrote polarization states as Stokes vectors and Mueller matrices, and [post 3](/en/posts/ellipsometry-thin-film-multilayer-reflectance/) computed $r_p$ and $r_s$ for a multilayer stack to obtain $\Psi$ and $\Delta$. Give it a thickness and a refractive index, and out come $\Psi$ and $\Delta$.

A real instrument runs the other way. You put a sample in and you have to measure $\Psi$ and $\Delta$. But a photodetector hands you a single real number, an intensity proportional to a current. How do you extract two complex numbers' worth of information from intensity measurements alone? The answers branched over the course of a century, and those branches are today's ellipsometer designs.

There are several ways to survey them, but the sharpest framing is a sentence Fujiwara puts at the top of section 4.2.1. Ellipsometry fundamentally measures Stokes parameters; the parameters that each design can reach differ; and the differences in measurement error follow directly from that. So this post does not simply list the methods. It takes what each one *cannot* measure as its axis, and traces how that gap comes back as a limit on the instrument.

## 1. Why Oblique Incidence

What happens if you send the beam in along the surface normal? On an isotropic sample there is no longer any direction that distinguishes p from s. The plane of incidence is undefined, so $r_p = r_s$ and $\rho = r_p/r_s = 1$. That pins $\Psi = 45°$ and $\Delta = 0°$ regardless of the sample, which tells you nothing about it.

This is why ellipsometry insists on oblique incidence. An angle of incidence $\theta$ defines the plane of incidence, gives p and s different Fresnel coefficients, and loads that difference onto $\rho$. The difference is largest near the Brewster angle.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig1-oblique-configuration.png" alt="Oblique-incidence geometry and the PSG-sample-PSA chain" width="720">
_Figure 1. Oblique-incidence geometry. A known polarization state goes in on the left; on the right a single intensity comes out. The $M_{sample}$ computed in post 3 sits in the middle of this chain._

The instrument splits into three blocks: the polarization state generator (PSG) on the source side, the sample, and the polarization state analyzer (PSA) on the detector side. The PSG turns unpolarized light from the source into a known polarization state and sends it to the sample; the PSA translates the polarization state of the reflected beam into an intensity the detector can read. The sample Mueller matrix $M_{sample}(\Psi, \Delta)$ from post 3 is one factor in the middle of this chain. The detected intensity comes from a product of this form:

$$I = \left[ M_{PSA} \, M_{sample} \, M_{PSG} \, \mathbf{S}_{in} \right]_0$$

The subscript 0 denotes the first component of the Stokes vector, the total intensity. Every calculation in this post is that one product evaluated repeatedly; the implementation is in `_code/ellipsometer-oblique-rotating-element/modulation.py`.

## 2. The Detector Gives Only One Intensity

State the problem precisely. Four Stokes parameters describe the polarization of the reflected beam completely. Send $+45°$ linear polarization through the $M_{sample}$ of post 3 and normalize by the total intensity: as long as the sample does not depolarize, the remaining three components carry $\Psi$ and $\Delta$ directly.

$$S_1 = -\cos 2\Psi, \qquad S_2 = \sin 2\Psi \cos\Delta, \qquad S_3 = -\sin 2\Psi \sin\Delta$$

The minus sign on $S_3$ is a matter of convention. $\Delta$ is defined as $\phi_p - \phi_s$, whereas the Stokes parameters take the phase difference in the order $\phi_s - \phi_p$, following the definition in post 2. It is the same sign that $M_{sample}$ of post 3 produces.

Knowing any two of these three yields $\Psi$ and $\Delta$. But the detector reads only $S_0$. It has no direct access to $S_1$, $S_2$ or $S_3$.

This is where modulation enters. Vary one component of the PSA in time and it will bring different Stokes components of the reflected beam into the $S_0$ slot in turn. The detected intensity then oscillates, and $S_1$, $S_2$ and $S_3$ take up residence in the Fourier coefficients of that oscillation. For a detector that reads a single intensity, this is the only route to polarization information: spend time to separate the unknowns one at a time.

What you vary, and how, is what separates the designs. A common notation writes the initials of the optical elements in order and attaches a subscript $R$ to the one that rotates: polarizer $P$, sample $S$, analyzer $A$, compensator $C$, photoelastic modulator $M$. A compensator is the same component as the retarder of post 2; here its job is to insert a deliberate phase difference between p and s.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig2-instrument-configurations.png" alt="Four optical configurations" width="720">
_Figure 2. The four configurations. RAE with a compensator and RCE have identical optics; only the choice of which element turns differs._

The second and third rows deserve attention. $PSCA_R$ and $PSC_RA$ contain exactly the same components. The only difference is whether the compensator is fixed and the analyzer turns, or the analyzer is fixed and the compensator turns — and yet, as shown below, what they can measure is not the same.

## 3. Null Ellipsometry

The earliest method varies nothing at all. In the technique Drude used at the end of the nineteenth century, the polarizer and analyzer are turned by hand until the detected intensity drops to zero, and the angles are read off at that point. The name comes from the search for that zero, the null. In a $PCSA$ configuration with the compensator at $C = 45°$ and a retardation of $\delta = 90°$, the extinction condition gives

$$\Psi = -A, \qquad \Delta = -2P + 90°$$

$\Psi$ and $\Delta$ come straight off the dials of the polarizer and analyzer. No Fourier analysis, no intensity calibration, no detector linearity required. Judging darkness is something the human eye can do unaided, which is how Drude was able to build an ellipsometer a century ago.

Four $(P, A)$ combinations satisfy the null condition, two for each of $C = \pm 45°$. Measuring all four and averaging cancels alignment errors in the components. This is four-zone averaging, and it is why null ellipsometry is still regarded as the most accurate approach. Because it never uses an absolute intensity — only the judgment of zero or not zero — it is free of detector nonlinearity in principle.

The price is time. A new $(P, A)$ pair must be found at every wavelength, which makes spectroscopic measurement extremely slow. The trade-off that runs through this whole post surfaces here: null ellipsometry refuses to trust the detector and pays in time, while the photometric methods that follow buy speed by trusting the intensity reading. Every error discussion in rotating-element instruments descends from that exchange.

## 4. RAE (Rotating-Analyzer Ellipsometry)

What happens if you simply spin the analyzer at constant speed? That is RAE. Only two polarizing components sit on the optical axis, the polarizer and the analyzer, which makes it the simplest of the four configurations. That simplicity is both its strength and its limit. What it gains and what it loses are worth taking in order.

### The Modulation Principle

With the polarizer at $P = 45°$ and the analyzer azimuth written as $A$, the detected intensity is

$$I = I_0 \left[ 1 + S_1 \cos 2A + S_2 \sin 2A \right]$$

There is a reason $2A$ appears rather than $A$. A transmission axis has no distinguishable up and down, so rotating it by $180°$ returns the original state. A $180°$ rotation of the analyzer is therefore one optical cycle. Spinning the analyzer at angular frequency $\omega$ so that $A = \omega t$ gives

$$I(t) = I_0 \left[ 1 + \alpha \cos 2\omega t + \beta \sin 2\omega t \right]$$

Extracting $\alpha$ and $\beta$ from the measured waveform and converting back to $\Psi$ and $\Delta$ uses

$$\tan\Psi = \sqrt{\frac{1+\alpha}{1-\alpha}} \, \lvert \tan P \rvert, \qquad \cos\Delta = \frac{\beta}{\sqrt{1-\alpha^2}}$$

One component on a motor, and that is the whole instrument. With no compensator there is no part whose behavior drifts with wavelength. This simplicity is RAE's greatest asset.

### Ambiguity in the Sign of the Phase Difference

Look again at the equation above and the problem appears: $\Delta$ enters only as $\cos\Delta$.

Suppose a measurement returns $\cos\Delta = 0.707$. Is $\Delta = +45°$ or $-45°$? The cosine is even, so both answers fit equally well, and nothing distinguishes them. The measurement range of RAE therefore shrinks to $0° \le \Delta \le 180°$, half of the full $-180° \le \Delta \le 180°$.

The root cause lies in the definitions of section 2. The sign of $\Delta$ is held by $S_3 = -\sin 2\Psi \sin\Delta$, and only $S_1$ and $S_2$ appear in the RAE intensity. Failing to measure $S_3$, which carries the handedness of circular polarization, is what erases the sign of the phase. This is the same statement as "an intensity measurement loses phase information" from post 2, repeating itself one level up, at the instrument.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig3-rae-delta-ambiguity.png" alt="Sign ambiguity of Delta in RAE" width="760">
_Figure 3. Two samples with identical $\Psi$ and opposite signs of $\Delta$. (a) In RAE the two waveforms differ by exactly zero. (b) In RCE the $2\omega$ component flips sign and separates them immediately._

The left panel of figure 3 confirms the ambiguity numerically. The RAE waveforms were computed for two samples, one at $\Delta = +60.4°$ and one at $\Delta = -60.4°$, and the maximum difference between the two curves is $0.0$ — identical at the level of floating-point representation, not merely to within rounding. No detector and no amount of averaging will separate these two samples with RAE.

The ambiguity is not the only problem. Differentiating $\beta$ with respect to $\Delta$ gives $\partial \beta / \partial \Delta \propto \sin\Delta$. The amplification factor with which measurement noise propagates into $\Delta$ is therefore $1/\lvert \sin\Delta \rvert$, which diverges for samples whose $\Delta$ lies near $0°$ or $180°$.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig4-rae-sensitivity.png" alt="Error amplification in Delta" width="760">
_Figure 4. (a) The RAE error amplification diverges at $\Delta \cong 0°, 180°$. (b) Combining four measurements whose operating points have been shifted by a compensator keeps the amplification near unity everywhere._

### Adding a Compensator

What changes if a compensator is inserted into RAE? A compensator leaves the amplitude ratio alone and shifts only the p–s phase difference, by $\delta$. So $\Psi$ is untouched and $\Delta$ merely moves to $\Delta' = \Delta - \delta$. The intensity becomes

$$I = I_0 \left[ 1 + S_1 \cos 2A + \left( S_2 \cos\delta - S_3 \sin\delta \right) \sin 2A \right]$$

$S_3$ has finally appeared. It is bundled together with $S_2$, however, so a single measurement cannot separate them. At least two measurements at different retardations are needed — $\delta = 0°$ for $S_2$, $\delta = 90°$ for $S_3$ — which is why this takes longer than plain RAE.

But the ability to move $\Delta$ is useful in its own right. RAE is accurate near $\Delta \cong 90°$ and poor at $0°$ and $180°$. If the operating point can be shifted, then measuring at several values of $\delta$ and stitching together only the portions where $\Delta' \sim \pm 90°$ solves the problem. The right panel of figure 4 shows this computed. Four measurements at $\delta = 0°, 45°, 90°, 135°$ each have their own divergence, but taking the best measurement at every point keeps the amplification near unity across the whole range. One extra component changes the measurement strategy itself.

## 5. RCE (Rotating-Compensator Ellipsometry)

What happens if the analyzer is held fixed and the compensator turns instead? That is RCE. The components are identical to the $PSCA_R$ of the previous section; only the choice of what rotates has changed. With the analyzer at $A = 0°$, a compensator retardation of $\delta = 90°$ and the compensator azimuth written as $C$:

$$I = I_0 \left[ (2 + S_1) - 2 S_3 \sin 2C + S_1 \cos 4C + S_2 \sin 4C \right]$$

Terms in $2C$ and $4C$ appear together, and — more importantly — $S_3$ rides on the $2C$ term while $S_1$ and $S_2$ ride on $4C$. Because the frequencies split, extracting the Fourier coefficients separates all three Stokes parameters at once. No second measurement at a different retardation is needed.

Why do the frequencies split? If the reflected light is linearly polarized, the beam leaving the rotating compensator returns to the same state every $90°$; if it is circularly polarized, every $180°$. The linear components $S_1$ and $S_2$ modulate the intensity with a $4C$ period, the circular component $S_3$ with a $2C$ period. A property of the polarization state is translated directly into a modulation frequency.

Obtaining $S_1$ through $S_3$ in a single measurement means $\Delta$ can be measured over the full $-180°$ to $180°$, and that the sensitivity is uniform across that range. This is why RCE became the standard in commercial instruments. There is a bonus as well: the degree of polarization can be extracted alongside the Fourier coefficients, which makes the depolarization discussed in post 2 measurable as a function of wavelength.

## 6. PME (Phase-Modulation Ellipsometry)

Every method so far turned a component mechanically. Can one modulate without turning anything? PME drives the retardation itself electrically, using a photoelastic modulator (PEM). Driving a quartz crystal into resonance with a piezoelectric transducer produces stress birefringence, and the retardation oscillates at the resonant frequency:

$$\delta(t) = F \sin \omega t$$

In a $PSMA$ configuration with $45°$ between the analyzer and the modulator, the intensity is

$$I = I_0 \left[ 1 - S_3 \sin\delta + \left( -S_1 \sin 2M + S_2 \cos 2M \right) \cos\delta \right]$$

Since $\sin\omega t$ sits inside $\sin\delta$ and $\cos\delta$, this has to be unfolded with the Jacobi-Anger expansion:

$$\sin\delta = 2\sum_{m=0}^{\infty} J_{2m+1}(F) \sin\left[(2m+1)\omega t\right], \qquad \cos\delta = J_0(F) + 2\sum_{m=1}^{\infty} J_{2m}(F) \cos 2m\omega t$$

The $J_k$ are Bessel functions. A practical trick attaches itself to this expansion. Adjusting the voltage on the PEM so that $F \cong 138°$ makes $J_0(F)$ vanish; the first zero of $J_0$ is at $137.79°$, and $138°$ is the rounded value. With that setting the sample information drops out of the DC term and the analysis simplifies:

$$I(t) = I_0 \left[ 1 + \sin 2\Psi \sin\Delta \cdot 2J_1(F) \sin\omega t + \sin 2\Psi \cos\Delta \cdot 2J_2(F) \cos 2\omega t \right]$$

With $2J_1 = 1.04$ and $2J_2 = 0.86$, $\Psi$ and $\Delta$ follow immediately from the coefficients at $\omega$ and $2\omega$. A fast Fourier transform of the signal from a real instrument shows two sharp peaks, at 50 kHz and 100 kHz. Note, though, that this expression keeps only the low-order terms; higher orders such as $2J_3 = 0.40$ are in fact not negligible.

Having no moving parts is PME's decisive advantage. The modulation rate of a rotating-element instrument is tied to how fast the element spins, which tops out at 10–100 Hz and sets a minimum measurement time of roughly 10 ms. A PEM resonates at 50 kHz, so one measurement takes 20 μs — three orders of magnitude apart. That is what makes it possible to follow processes on sub-millisecond timescales, such as the response of liquid crystal molecules.

## 7. One Structure, Different Modulation Frequencies

Taken separately the three methods look like different instruments, but set their equations side by side and the structure is one. All three modulate something periodically to make the detected intensity oscillate, then read $\Psi$ and $\Delta$ out of the Fourier coefficients of that waveform. The only difference is where the modulation frequency lands.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig5-waveform-comparison.png" alt="Detected intensity and harmonic content for three methods" width="780">
_Figure 5. One sample (SiO$_2$ 100 nm / Si at $75°$ incidence) measured three ways. Left, the detected intensity over one modulation cycle; right, the normalized Fourier coefficient magnitudes of that waveform._

Figure 5 shows this. For a single sample ($\Psi = 41.2°$, $\Delta = 60.4°$) the waveforms of all three methods were computed and their harmonic content placed side by side. RAE raises a single bar at $2\omega$; RCE raises two, at $2\omega$ and $4\omega$; PME sits mainly at $\omega$ and $2\omega$ but retains a non-negligible bar at $3\omega$. The higher-order Bessel terms mentioned in section 6 show up in the figure exactly as expected.

| Method | Modulation | Frequencies | What it can measure |
| --- | --- | --- | --- |
| RAE | rotating analyzer | $2\omega$ | $S_0, S_1, S_2$ |
| RAE + compensator | rotating analyzer | $2\omega$ | $S_0$–$S_3$ (two measurements) |
| RCE | rotating compensator | $2\omega, 4\omega$ | $S_0$–$S_3$ (single measurement) |
| PME | electrically driven retardation | $\omega, 2\omega$ | $S_0, S_1, S_3$ or $S_0, S_2, S_3$ |

This unified view has a practical consequence: how many times the detector must sample per revolution is set by the number of unknown coefficients. RAE has to determine $I_0$, $\alpha$ and $\beta$, so it needs at least three intervals; RCE has five coefficients and needs at least five. The more coefficients, the more samples per revolution, and the slower the measurement. How badly this scales for instruments that measure a full Mueller matrix is a subject for later posts.

## 8. Comparing the Methods

| Method | $\Delta$ range | Minimum measurement time | Wavelengths | Character |
| --- | --- | --- | --- | --- |
| RAE | $0°$ to $180°$ (half) | about 10 ms | about 200 | achromatic |
| RAE + compensator | $-180°$ to $180°$ | 10 ms or more | about 200 | chromatic |
| RCE | $-180°$ to $180°$ | about 10 ms | about 200 | chromatic |
| PME | errors in specific ranges | 20 μs | about 10 | chromatic |

The last column is the hook into the next post. RAE uses only a polarizer and an analyzer, and those components barely change their behavior across a wide spectral band. Such an instrument is called achromatic, and that is why it can measure as many wavelengths simultaneously as a photodiode array has pixels. RCE, by contrast, has a compensator whose retardation depends on wavelength, and PME requires a different drive voltage at every wavelength to hold the retardation constant. Both are chromatic. This is the fundamental reason the number of wavelengths PME can handle in real time is an order of magnitude smaller.

The wavelength dependence of a single component determines the spectral capability of the entire instrument — which is precisely the subject of the next post.

## Summary

Three things came out of this post.

First, ellipsometer designs are separated by what they cannot measure rather than by what they can. Because RAE misses $S_3$, the sign of $\Delta$ disappears, the measurement range is halved, and the error diverges at $\Delta \cong 0°, 180°$. Adding a compensator brings $S_3$ into the equation but bundled with $S_2$, so two measurements are required; rotating the compensator instead drops $S_3$ onto its own $2\omega$ term and settles the matter in one. The lineage of these instruments is a history of filling that gap.

Second, the exchange between accuracy and time determines the configuration. Null ellipsometry declines to trust the detector and pays by hunting for angles at every wavelength; the photometric methods trust the intensity reading and are fast. Why detector nonlinearity and calibration matter so much in rotating-element instruments is the other face of that exchange.

Third, the methods differ in appearance but not in structure. Modulate, extract Fourier coefficients, read $\Psi$ and $\Delta$ — the procedure is the same for all three, and only the modulation frequency differs, $2\omega$ or $2\omega$ and $4\omega$ or $\omega$ and $2\omega$. Figure 5 confirms that this difference is visible in the signal itself for one and the same sample.

The next post examines the components that make up this chain one by one: how polarizers, retarders and photoelastic modulators work, what limits each of them carries, and how the achromatic/chromatic split seen at the end of this post turns into concrete constraints on component choice.

## Reproducing the Calculations

Every waveform and Fourier coefficient in this post was computed directly as a product of component Mueller matrices. The code is in `_code/ellipsometer-oblique-rotating-element/`.

- `modulation.py` — computes the detected intensity $I(t)$ for all four configurations as a Mueller product. It uses no closed-form expressions.
- `verify_modulation.py` — checks those numerical results against Fujiwara's closed forms (4.18, 4.28, 4.32, 4.41, 4.44). The two paths have to be independent for a sign-convention error to show itself. This cross-check did in fact catch a misreading introduced while transcribing one of the reference equations.
- `generate_figures.py` — produces the five figures.
- `smm_tensor.py` — the scattering matrix implementation from post 3, which supplies the $\Psi$ and $\Delta$ of the reference sample in figure 5.

## References

- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, sections 4.2–4.3 (comparison of measurement methods and Fourier analysis).
- R. W. Collins, "Automatic rotating element ellipsometers: Calibration, operation, and real-time applications," *Rev. Sci. Instrum.* **61**, 2029 (1990).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1987, sections 5.2–5.3 and 5.7 (null ellipsometers and rotating-element methods).
- H. G. Tompkins and E. A. Irene (eds.), *Handbook of Ellipsometry*, William Andrew, 2005, chapters 5–6 (rotating-element and phase-modulation instruments).
- Youngjoon Kim, "Development of Snapshot Angle-Resolved Spectroscopic Ellipsometry Using a Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD thesis, Seoul National University, 2025, section 1.2.

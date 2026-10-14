---
title: "Ellipsometer Instrumentation 1 — Why Shine the Light at an Angle and Rotate the Polarizer?"
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

A real instrument runs the other way. You put a sample in and you have to measure $\Psi$ and $\Delta$. But a photodetector hands you a single real number, an intensity proportional to a current. How do you extract one complex number's worth of information, $\rho$, or two real numbers, $\Psi$ and $\Delta$, from intensity measurements alone? The answers branched over more than a century, and those branches are today's ellipsometer designs.

This series follows the questions an instrument user is likely to ask.

- Why does an ellipsometer send light in obliquely, and why near 70°? (section 1 of this post)
- Why can an instrument with a rotating analyzer not tell $\Delta = +45°$ from $-45°$? (section 4)
- Why do the latest commercial instruments rotate a compensator instead of the analyzer? (section 5)
- Why does a single component change the number of wavelengths measured at once by an order of magnitude? ([post 2](/en/posts/ellipsometer-components/))
- What has to rotate, and how, to measure all 16 elements of a Mueller matrix? ([post 3](/en/posts/ellipsometer-dual-rotating-compensator-mme/))
- Is calibration really necessary? (post 4)
- Can one measure in a single shot with nothing rotating? (posts 5 and 6)

There are several ways to survey the designs. Fujiwara sums up the question in one sentence in section 4.2.1. Ellipsometry fundamentally measures Stokes parameters. The parameters each design can reach differ, and so do the measurement errors, for exactly that reason. This post therefore does not simply list the methods. It takes what each one *cannot* measure as its axis and traces how that gap comes back as a limit on the instrument.

## 1. Why Is the Light Sent in Obliquely?

What happens if you send the beam in along the surface normal? On an isotropic sample there is no longer any direction that distinguishes p from s. The plane of incidence is undefined, so both polarizations reflect identically. In the Fresnel sign convention of Ellipsometry Foundations post 1 this appears as $r_p = -r_s$, giving $\rho = r_p/r_s = -1$. That pins $\Psi = 45°$ and $\Delta = 180°$ regardless of the sample, which tells you nothing about it.

This is why ellipsometry insists on oblique incidence. An angle of incidence $\theta$ defines the plane of incidence, gives p and s different Fresnel coefficients, and loads that difference onto $\rho$. The difference is largest near the Brewster angle. For an absorbing sample such as Si, the relevant angle is the pseudo-Brewster angle of [Ellipsometry Foundations post 1](/en/posts/ellipsometry-electromagnetic-fresnel/), near 75.6° for Si at 632.8 nm.

How large is the effect? Take 2 nm of SiO$_2$ on Si and ask how much $\Delta$ changes when the film grows by 1 nm (632.8 nm, SiO$_2$ index 1.457, Si $3.882 - 0.019j$). The change is $-0.017°$ at 10°, $-0.42°$ at 45° and $-2.84°$ at 70°. The sensitivity at 70° is about 170 times that at 10°. Near normal incidence the sample tells you almost nothing. This is also why Fujiwara gives 70–80° as the usual angle of incidence for semiconductors (section 2.3.4).

That does not mean one should insist on the pseudo-Brewster angle itself. On bare Si, $\Psi$ falls to 0.15° at 75.55°. With $\Psi$ that close to zero, almost no p-polarized light is reflected and $\Delta$ becomes hard to read. The thickness sensitivity above also changes sign, from $-13.1°$/nm at 75° to $+11.7°$/nm at 76°. Right next to the angle of greatest sensitivity, a small angle error moves the result a lot.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig1-oblique-configuration.png" alt="Oblique-incidence geometry and the PSG-sample-PSA chain" width="720">
_Figure 1. Oblique-incidence geometry. A known polarization state goes in on the left; on the right a single intensity comes out. The sample Mueller matrix $$M_{sample}$$ computed in Ellipsometry Foundations post 3 is the middle factor of this chain._

The instrument splits into three blocks: the polarization state generator (PSG) on the source side, the sample, and the polarization state analyzer (PSA) on the detector side. The PSG turns unpolarized light from the source into a known polarization state and sends it to the sample. The PSA converts the polarization state of the reflected beam into an intensity the detector can read. The sample Mueller matrix $M_{sample}(\Psi, \Delta)$ is one factor in the middle of this chain. The detected intensity comes from the product

$$I = \left[ M_{PSA} \, M_{sample} \, M_{PSG} \, \mathbf{S}_{in} \right]_0$$

The subscript 0 denotes the first component of the Stokes vector, the total intensity. Every waveform and Fourier coefficient in this post was computed directly as this product of component Mueller matrices; only figure 6 is drawn from a closed-form expression. The results were checked separately against Fujiwara's closed forms (4.18, 4.28, 4.32, 4.41, 4.44), sign conventions included.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/ext-ellipsometer-laas.jpg" alt="A variable-angle ellipsometer with two arms on a circular rail" width="720">
_Figure 2. A variable-angle ellipsometer. Two arms face each other across the sample stage, mounted on a circular rail. One arm is the polarization state generator and the other the polarization state analyzer; moving them along the rail changes the angle of incidence. (Source: Guillaume Paumier (user:guillom), [Ellipsometer at LAAS](https://commons.wikimedia.org/wiki/File:Ellipsometer_at_LAAS.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

Figure 2 shows a real variable-angle ellipsometer. The PSG and PSA of figure 1 appear as the two arms, and the sample sits on the stage between them.

## 2. How Can a Detector That Reads One Intensity Measure Polarization?

State the problem precisely. Four Stokes parameters describe the polarization of the reflected beam completely. Send $+45°$ linear polarization through the sample Mueller matrix and normalize by the total intensity. As long as the sample does not depolarize, the remaining three components carry $\Psi$ and $\Delta$ directly.

$$S_1 = -\cos 2\Psi, \qquad S_2 = \sin 2\Psi \cos\Delta, \qquad S_3 = -\sin 2\Psi \sin\Delta$$

The minus sign on $S_3$ is a matter of convention. $\Delta$ is defined as $\phi_p - \phi_s$, whereas the Stokes parameters take the phase difference in the order $\phi_s - \phi_p$, following the definition in Ellipsometry Foundations post 2. It is the same sign that $M_{sample}$ produces.

Knowing all three fixes $\Psi$ and $\Delta$ over their full ranges. With only two, the range shrinks; as we will see, $S_1$ and $S_2$ alone leave the sign of $\Delta$ undetermined. But the detector reads only $S_0$. It has no direct access to $S_1$, $S_2$ or $S_3$.

This is where modulation enters. Vary one component of the PSA in time and it brings different Stokes components of the reflected beam into the $S_0$ slot in turn. The detected intensity then oscillates, and $S_1$, $S_2$ and $S_3$ appear in the Fourier coefficients of that oscillation. With one detector at one wavelength, spending time to separate the unknowns one by one is the basic route. Spending wavelength instead of time is also possible; post 5 of this series (channeled spectroscopy) covers it.

What you vary, and how, is what separates the designs. A common notation writes the initials of the optical elements in order and attaches a subscript $R$ to the one that rotates: polarizer $P$, sample $S$, analyzer $A$, compensator $C$, photoelastic modulator $M$. A compensator is the same component as the retarder of Ellipsometry Foundations post 2. Here its job is to insert a deliberate phase difference between p and s.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig2-instrument-configurations.png" alt="Four optical configurations" width="720">
_Figure 3. The four configurations. RAE with a compensator and RCE have identical optics; only the choice of which element turns differs._

The second and third rows deserve attention. $PSCA_R$ and $PSC_RA$ contain exactly the same components. The only difference is whether the compensator is fixed and the analyzer turns, or the analyzer is fixed and the compensator turns. Yet, as shown below, they do not measure the same things.

## 3. Null Ellipsometry

Can $\Psi$ and $\Delta$ be measured with nothing rotating at all? The earliest method does exactly that. In the technique Drude used at the end of the nineteenth century, the polarizer and analyzer are turned by hand until the detected intensity drops to zero, and the angles are read off at that point. The name comes from the search for that zero, the null. In a $PCSA$ configuration with the compensator's fast axis at $C = 45°$ and a retardance of $\delta = 90°$, the extinction condition gives

$$\Psi = -A, \qquad \Delta = -2P + 90°$$

Here $P$ and $A$ are the polarizer and analyzer azimuths read at extinction. Several $(P, A)$ combinations give extinction, and they are sorted into zones; the equation above holds in the zone with $-A > 0$. In other words, $\Psi$ and $\Delta$ come straight off the dials. No Fourier analysis, no intensity calibration, no detector linearity is required. Judging darkness is something the human eye can do unaided, which is how Drude was able to build an ellipsometer some 130 years ago.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/ext-nbs1963-null-ellipsometer.jpg" alt="Layout of a PCSA null ellipsometer from a 1963 NBS paper" width="720">
_Figure 4. Layout of a PCSA null ellipsometer from a 1963 paper of the US National Bureau of Standards (NBS). Light from the source passes a monochromatic filter, and a polarizer and compensator make it elliptically polarized before it reaches the sample. When the reflected light is linearly polarized, the analyzer extinguishes it, and a photometer connected to a photocell shows whether the intensity is essentially zero. (Source: F. L. McCrackin et al., J. Res. NBS 67A, 363 (1963), Figure 5, public domain, [NIST original](https://nvlpubs.nist.gov/nistpubs/jres/67A/jresv67An4p363_A1b.pdf), cropped and border trimmed)_

Figure 4 shows the actual parts of a hand-read null ellipsometer. The polarizer and compensator sit on the source side, the analyzer and photocell on the detector side, and the polarizer and analyzer angles are adjusted in turn until the meter needle approaches zero.

Four $(P, A)$ combinations satisfy the null condition, two for each of $C = \pm 45°$. Measuring all four zones and averaging cancels alignment errors in the components. This is four-zone averaging. Because the method never uses an absolute intensity, only the judgment of zero or not zero, null ellipsometry is still regarded as the most accurate approach.

The price is time. A new $(P, A)$ pair must be found at every wavelength, which makes spectroscopic measurement extremely slow. The trade-off that runs through this whole post surfaces here. Null ellipsometry refuses to trust the detector and pays in time; the photometric methods that follow buy speed by trusting the intensity reading. Every error discussion in rotating-element instruments branches from that exchange.

## 4. RAE (Rotating-Analyzer Ellipsometry)

What happens if you simply spin the analyzer at constant speed? That is RAE. Only two polarizing components sit on the optical axis, the polarizer and the analyzer, which makes it the simplest of the four configurations. That simplicity is both its strength and its limit. We take what it gains and what it loses in turn.

### How Does the Intensity Oscillate When the Analyzer Turns?

With the polarizer at $P = 45°$ and the analyzer azimuth written as $A$, the detected intensity is

$$I = I_0 \left[ 1 + S_1 \cos 2A + S_2 \sin 2A \right]$$

There is a reason $2A$ appears rather than $A$. A transmission axis has no distinguishable up and down, so rotating it by $180°$ returns the original state. A $180°$ rotation of the analyzer is therefore one optical cycle. Spinning the analyzer at angular frequency $\omega$ so that $A = \omega t$ gives

$$I(t) = I_0 \left[ 1 + \alpha \cos 2\omega t + \beta \sin 2\omega t \right]$$

Extracting $\alpha$ and $\beta$ from the measured waveform and converting back to $\Psi$ and $\Delta$ uses

$$\tan\Psi = \sqrt{\frac{1+\alpha}{1-\alpha}} \, \lvert \tan P \rvert, \qquad \cos\Delta = \frac{\beta}{\sqrt{1-\alpha^2}}$$

So $\alpha$ mainly carries $\Psi$, and $\beta$ carries $\Delta$. One component on a motor is the whole instrument. With no compensator there is no part whose behavior drifts with wavelength. As section 8 shows, this becomes RAE's main asset in spectroscopic measurement.

### Why Does the Sign of the Phase Difference Disappear?

Look again at the equation above and the problem appears: $\Delta$ enters only as $\cos\Delta$.

Suppose a measurement returns $\cos\Delta = 0.707$. Is $\Delta = +45°$ or $-45°$? The cosine is even, so both answers fit equally well, and nothing distinguishes them. The measurement range of RAE therefore shrinks to $0° \le \Delta \le 180°$, half of the full $-180° \le \Delta \le 180°$.

The root cause lies in the definitions of section 2. The sign of $\Delta$ is held by $S_3 = -\sin 2\Psi \sin\Delta$, and only $S_1$ and $S_2$ appear in the RAE intensity. Failing to measure $S_3$, which carries the handedness of circular polarization, is what erases the sign of the phase. This is the statement "an intensity measurement loses phase information" from Ellipsometry Foundations post 1, repeating itself at the level of the instrument.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig3-rae-delta-ambiguity.png" alt="Sign ambiguity of Delta in RAE" width="760">
_Figure 5. The reference sample (SiO$$_2$$ 100 nm / Si, $75°$ incidence, $\Psi = 41.2°$, $\Delta = +60.4°$) and a sample with only the sign of $\Delta$ reversed. (a) In RAE the two waveforms differ by exactly zero. (b) In RCE the $2\omega$ component flips sign and separates them immediately._

The left panel of figure 5 confirms the ambiguity numerically. The RAE waveforms were computed for $\Delta = +60.4°$ and $\Delta = -60.4°$, and the maximum difference between the two curves is $0.0$. They are identical at the level of floating-point representation, not merely to within rounding. No detector and no amount of averaging will separate these two samples with RAE.

The ambiguity is not the only problem. Differentiating $\beta$ with respect to $\Delta$ gives $\partial \beta / \partial \Delta \propto \sin\Delta$. The factor by which measurement noise is amplified into $\Delta$ is therefore $1/(\sin 2\Psi \, \lvert \sin\Delta \rvert)$, which as a function of $\Delta$ alone goes as $1/\lvert \sin\Delta \rvert$. It diverges for samples whose $\Delta$ lies near $0°$ or $180°$.

Such samples are common. At 70° incidence, bare Si has $\Delta = 179.2°$, so $1/\lvert \sin\Delta \rvert$ is 74; Si under 2 nm of SiO$_2$ has $\Delta = 173.5°$ and a factor of 8.9. The 100 nm oxide sample of figure 7, by contrast, has $\Delta = 79.8°$ at the same angle and a factor of only 1.02. Thin films and bare substrates sit squarely in RAE's weak region. That is the practical reason for adding a compensator to RAE, or for rotating the compensator instead.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig4-rae-sensitivity.png" alt="Error amplification in Delta" width="760">
_Figure 6. (a) Without a compensator, the RAE error amplification diverges at $\Delta \cong 0°, 180°$. (b) Taking, at each point, the best of four measurements whose operating points are shifted by $0°, 45°, 90°$ and $135°$ keeps the amplification at most 1.08 across the whole range._

### What Changes When a Compensator Is Added?

A compensator leaves the amplitude ratio alone and shifts only the p–s phase difference, by $\delta$. So $\Psi$ is untouched and $\Delta$ merely moves to $$\Delta' = \Delta - \delta$$. The intensity becomes

$$I = I_0 \left[ 1 + S_1 \cos 2A + \left( S_2 \cos\delta - S_3 \sin\delta \right) \sin 2A \right]$$

$S_3$ has finally appeared. It is bundled together with $S_2$, however, so a single measurement cannot separate them. At least two measurements at different retardances are needed, for example $\delta = 0°$ for $S_2$ and $\delta = 90°$ for $S_3$. This is why the measurement takes longer than plain RAE.

But the ability to move $\Delta$ is useful in its own right. RAE is accurate near $\Delta \cong 90°$ and poor at $0°$ and $180°$. If the operating point can be shifted, one can measure at several values of $\delta$ and stitch together only the portions where $$\Delta' \sim \pm 90°$$. The right panel of figure 6 shows this computed. Four measurements at $\delta = 0°, 45°, 90°, 135°$ each have their own divergence, but taking the best measurement at every point keeps the amplification at most 1.08 ($= 1/\cos 22.5°$) across the whole range. One extra component changes the measurement strategy itself.

## 5. RCE (Rotating-Compensator Ellipsometry)

What happens if the analyzer is held fixed and the compensator turns instead? That is RCE. The components are identical to the $PSCA_R$ of the previous section; only the choice of what rotates has changed. With the analyzer at $A = 0°$, a compensator retardance of $\delta = 90°$ and the compensator azimuth written as $C$:

$$I = I_0 \left[ (2 + S_1) - 2 S_3 \sin 2C + S_1 \cos 4C + S_2 \sin 4C \right]$$

Here $C$ is the azimuth of the compensator's fast axis. The $$\Delta' = \Delta - \delta$$ of section 4 assumes the slow axis lies along p at $C = 0°$, so the two sections' $C$ differ by $90°$ (Fujiwara, end of 4.2.3). Mixing the conventions flips the sign of the $\sin 2C$ term and with it the sign of $\Delta$.

Terms in $2C$ and $4C$ appear together. More importantly, $S_3$ rides on the $2C$ term while $S_1$ and $S_2$ ride on $4C$. Because the frequencies split, extracting the Fourier coefficients separates all three Stokes parameters at once. No second measurement at a different retardance is needed.

Why do the frequencies split? If the reflected light is linearly polarized, the detected intensity behind the rotating compensator repeats every $90°$; if it is circularly polarized, every $180°$. For linear input, a $90°$ turn yields the same ellipse with opposite handedness, which the analyzer cannot tell apart. The linear components $S_1$ and $S_2$ therefore modulate the intensity with a $4C$ period, and the circular component $S_3$ with a $2C$ period. A property of the polarization state is translated directly into a modulation frequency.

Obtaining $S_1$ through $S_3$ in a single measurement means $\Delta$ can be measured over the full $-180°$ to $180°$, with uniform sensitivity across that range. This is why the latest commercial instruments use RCE. There is a bonus as well: the degree of polarization can be extracted from the same Fourier coefficients, which makes the depolarization discussed in Ellipsometry Foundations post 2 measurable as a function of wavelength.

## 6. PME (Phase-Modulation Ellipsometry)

Every method so far turned a component mechanically. Can one modulate without turning anything? PME drives the retardance itself electrically, using a photoelastic modulator (PEM). Driving a fused-quartz bar into resonance with a piezoelectric transducer makes the otherwise isotropic bar birefringent under stress, and the retardance oscillates at the resonant frequency:

$$\delta(t) = F \sin \omega t$$

In a $PSMA$ configuration with $45°$ between the analyzer and the modulator, the intensity is

$$I = I_0 \left[ 1 - S_3 \sin\delta + \left( -S_1 \sin 2M + S_2 \cos 2M \right) \cos\delta \right]$$

Since $\sin\omega t$ sits inside $\sin\delta$ and $\cos\delta$, this has to be unfolded with the Jacobi-Anger expansion:

$$\sin\delta = 2\sum_{m=0}^{\infty} J_{2m+1}(F) \sin\left[(2m+1)\omega t\right], \qquad \cos\delta = J_0(F) + 2\sum_{m=1}^{\infty} J_{2m}(F) \cos 2m\omega t$$

The $J_k$ are Bessel functions. A practical trick comes with this expansion. Adjusting the voltage on the PEM so that $F \cong 138°$ makes $J_0(F)$ vanish; the first zero of $J_0$ is at $137.79°$, and $138°$ is the rounded value. With that setting the sample information drops out of the DC term and the analysis simplifies. With the modulator at $M = 0°$ and the analyzer at $A = 45°$, the intensity becomes

$$I(t) = I_0 \left[ 1 + \sin 2\Psi \sin\Delta \cdot 2J_1(F) \sin\omega t + \sin 2\Psi \cos\Delta \cdot 2J_2(F) \cos 2\omega t \right]$$

With $2J_1 = 1.04$ and $2J_2 = 0.86$, $\Psi$ and $\Delta$ follow immediately from the coefficients at $\omega$ and $2\omega$. A fast Fourier transform of the signal from a real instrument shows two sharp peaks, at 50 kHz and 100 kHz. This expression keeps only the low-order terms, though; higher orders such as $2J_3 = 0.40$ are not negligible in practice.

Having no moving parts is PME's decisive advantage. The modulation rate of a rotating-element instrument is tied to how fast the element spins, which tops out at 10–100 Hz and sets a minimum measurement time of roughly 10 ms. A PEM resonates at 50 kHz, so one measurement takes 20 μs, about 500 times faster. That is what makes it possible to follow processes on sub-millisecond timescales, such as the response of liquid crystal molecules.

## 7. What Do the Three Methods Share, and Where Do They Differ?

Taken separately the three methods look like different instruments, but set their equations side by side and the structure is one. All three modulate something periodically to make the detected intensity oscillate, then read $\Psi$ and $\Delta$ out of the Fourier coefficients of that waveform. The only difference is where the modulation frequency lands.

<img src="/assets/img/posts/ellipsometer-oblique-rotating-element/en/fig5-waveform-comparison.png" alt="Detected intensity and harmonic content for three methods" width="780">
_Figure 7. One sample (SiO$$_2$$ 100 nm / Si at $75°$ incidence) measured three ways. Left, the detected intensity over one modulation cycle; right, the normalized Fourier coefficient magnitudes of that waveform. RAE shows a bar at $2\omega$ only and RCE at $2\omega$ and $4\omega$, while PME shows $\omega$ and $2\omega$ plus a residual $3\omega$ component._

Figure 7 shows this. For a single sample ($\Psi = 41.2°$, $\Delta = 60.4°$) the waveforms of all three methods were computed and their harmonic content placed side by side. The sample's $\Psi$ and $\Delta$ come from the scattering matrix method of Ellipsometry Foundations post 3. RAE raises a single bar at $2\omega$. RCE raises two, at $2\omega$ and $4\omega$. PME sits mainly at $\omega$ and $2\omega$ but keeps a non-negligible bar at $3\omega$; the higher-order Bessel terms of section 6 show up in the figure as expected.

| Method | Modulation | Frequencies | What it can measure |
| --- | --- | --- | --- |
| RAE | rotating analyzer | $2\omega$ | $S_0, S_1, S_2$ |
| RAE + compensator | rotating analyzer | $2\omega$ | $S_0$–$S_3$ (two measurements) |
| RCE | rotating compensator | $2\omega, 4\omega$ | $S_0$–$S_3$ (single measurement) |
| PME | electrically driven retardance | $\omega, 2\omega$ | $S_0, S_1, S_3$ or $S_0, S_2, S_3$ |

This view has a practical consequence. The number of detector samples needed per revolution is set by the number of unknown coefficients. RAE has to determine $I_0$, $\alpha$ and $\beta$, so it needs at least three intervals. In the configuration above the $\cos 2C$ term of RCE vanishes, but in general five coefficients are needed (DC plus the cosine and sine terms at $2\omega$ and $4\omega$), so RCE needs at least five. The more coefficients, the more samples per revolution, and the slower the measurement. How this grows for instruments that measure a full Mueller matrix is the subject of [post 3](/en/posts/ellipsometer-dual-rotating-compensator-mme/).

## 8. Why Do the Methods Differ in How Many Wavelengths They Measure?

Why, then, can PME handle only about ten wavelengths at once? Collecting the properties of each method in one table points to the answer.

| Method | $\Delta$ range | Minimum measurement time | Wavelengths | Character |
| --- | --- | --- | --- | --- |
| RAE | $0°$ to $180°$ (half) | about 10 ms | about 200 | achromatic |
| RAE + compensator | $-180°$ to $180°$ | 10 ms or more | about 200 | chromatic |
| RCE | $-180°$ to $180°$ | about 10 ms | about 200 | chromatic |
| PME | restricted, depending on configuration (errors grow in specific ranges) | 20 μs | about 10 | chromatic |

The answer is in the last column. RAE uses only a polarizer and an analyzer, and those components barely change their behavior across a wide spectral band. Such an instrument is called achromatic, and that is why it can measure as many wavelengths simultaneously as a photodiode array has pixels. RCE, by contrast, has a compensator whose retardance depends on wavelength, and PME requires a different drive voltage at every wavelength to hold the retardance constant. Both are chromatic. This is the fundamental reason the number of wavelengths PME can handle in real time is an order of magnitude smaller.

The wavelength dependence of a single component thus sets the spectral capability of the whole instrument. That is the subject of [post 2](/en/posts/ellipsometer-components/).

## Summary

Section 1 showed in numbers why the light goes in obliquely. Near normal incidence p and s do not separate, and the thickness sensitivity of a 2 nm oxide at 70° is about 170 times that at 10°. Beyond that, three things came out of this post.

First, ellipsometer designs are separated by what they cannot measure rather than by what they can. Because RAE misses $S_3$, the sign of $\Delta$ disappears, the measurement range is halved, and the error diverges at $\Delta \cong 0°, 180°$. Adding a compensator brings $S_3$ into the equation but bundled with $S_2$, so two measurements are required. Rotating the compensator instead puts $S_3$ on its own $2\omega$ term and settles the matter in one. The lineage of these instruments is a history of filling that gap.

Second, the exchange between accuracy and time determines the configuration. Null ellipsometry declines to trust the detector and pays by hunting for angles at every wavelength; the photometric methods trust the intensity reading and are fast. Why detector nonlinearity and calibration matter so much in rotating-element instruments is the other face of that exchange.

Third, the methods differ in appearance but not in structure. All three modulate, extract Fourier coefficients and read $\Psi$ and $\Delta$. Only the modulation frequencies differ: $2\omega$ for RAE, $2\omega$ and $4\omega$ for RCE, and $\omega$ and $2\omega$ for PME. Figure 7 shows that for one and the same sample the harmonic bars stand in different places.

[Post 2](/en/posts/ellipsometer-components/) examines the components that make up this chain one by one: how polarizers, retarders and photoelastic modulators work, what limits each of them carries, and how the achromatic/chromatic split seen at the end of this post turns into constraints on component choice.

## References

- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, section 2.3.4 (Brewster angle and choice of incidence angle), sections 4.2–4.3 (comparison of measurement methods and Fourier analysis).
- R. W. Collins, "Automatic rotating element ellipsometers: Calibration, operation, and real-time applications," *Rev. Sci. Instrum.* **61**, 2029 (1990).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, sections 5.2–5.3 and 5.7 (null ellipsometers and rotating-element methods).
- H. G. Tompkins and E. A. Irene (eds.), *Handbook of Ellipsometry*, William Andrew, 2005, chapters 5–6 (rotating-element and phase-modulation instruments).
- F. L. McCrackin, E. Passaglia, R. R. Stromberg, and H. L. Steinberg, "Measurement of the thickness and refractive index of very thin films and the optical properties of surfaces by ellipsometry," *J. Res. Natl. Bur. Stand.* **67A**, 363 (1963). [NIST](https://nvlpubs.nist.gov/nistpubs/jres/67A/jresv67An4p363_A1b.pdf)
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025, section 1.2.

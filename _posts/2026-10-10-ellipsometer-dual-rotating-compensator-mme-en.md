---
date: 2026-10-10 20:00:00 +0900
title: "Ellipsometer Instrumentation 3 — Dual Rotating Compensators and the Mueller Matrix"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-dual-rotating-compensator-mme/
page_id: ellipsometer-dual-rotating-compensator-mme
categories: [Optics, Instrumentation]
tags: [ellipsometry, mueller-matrix, dual-rotating-compensator, fourier-analysis, instrumentation]
description: Modulating only the generator or only the analyzer cannot deliver all sixteen elements at once. This post follows why two compensators must turn at different speeds.
math: true
---

[Post 1](/en/posts/ellipsometer-oblique-rotating-element/) sorted the methods by what they cannot measure. RAE misses $S_3$ and loses the sign of $\Delta$; a compensator recovers it at the cost of a second run; rotating that compensator settles it in one. [Post 2](/en/posts/ellipsometer-components/) followed where each of those components stops working across a spectral band.

This post asks something one level deeper. Every instrument in post 1 compressed the sample into two numbers, $\Psi$ and $\Delta$. What happens when that compression does not hold?

The short answer is that the whole Mueller matrix has to be measured, all sixteen elements of it, and that one compensator is not enough to do it. Azzam stated the obstruction outright in the opening of a 1978 paper.

> The simultaneous determination of all sixteen elements of the Mueller matrix by Fourier analysis of one detected signal is not possible if either the polarizing or the analyzing optics alone are modulated.

Post 1 asked what a method cannot measure, post 2 asked where a component stops working, and this one asks why one element will not do. All three start from something missing.

## 1. Why All Sixteen Are Needed

The number of independent parameters the sample demands comes first.

An isotropic sample leaves only three independent parameters in its Mueller matrix. They are the $S_1$, $S_2$ and $S_3$ defined in [post 2 of the foundations series](/en/posts/ellipsometry-polarization-mueller-matrix/), and the RCE of post 1 obtains all three in a single measurement, so nothing is lacking. An anisotropic sample raises that count to six, which existing instruments can still fill by measuring at several component angles.

The difficulty arrives when an anisotropic sample also depolarizes. At least seven parameters are then required, and the Jones matrix cannot describe the sample at all. The limitation of the Jones formalism to fully polarized light, discussed in foundations post 2, becomes a hardware requirement here. The only remaining route is to measure the Mueller matrix directly.

## 2. Why One Element Will Not Do

Attempts to measure a Mueller matrix with rotating elements began in the late 1970s. Hauge's 1980 survey of what each configuration reaches is the skeleton of this post.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig1-mme-configurations.png" alt="Mueller elements reachable by each rotating-element configuration" width="780">
_Figure 1. Four rotating-element configurations and the Mueller elements each one measures. Only the last, with two compensators turning at different speeds, fills all sixteen._

Rotating the polarizer and analyzer together misses the fourth row and the fourth column. With no component to produce circular polarization, that is only to be expected. Adding a single compensator improves matters, and where it goes decides the outcome: placed after the sample it yields the first three columns, placed before it the first three rows.

What deserves attention is that these two configurations have exactly the optical layout of the RAE-with-compensator and the RCE of post 1. Section 2 of that post argued that identical components give different reach depending on which one turns. The same argument reappears one level up, at the Mueller matrix. Where a component sits decides whether rows or columns come back.

Three rows or three columns do not make sixteen. Among rotating-element configurations, only two compensators turning at different frequencies deliver the whole matrix in a single measurement. The design was proposed in the late 1970s, built in the early 1990s, and an instrument capable of real-time spectroscopic measurement appeared in 2000.

Photoelastic modulators can also measure a Mueller matrix. Two modulators at different resonant frequencies give nine elements, and four measurements at different component angles fill the rest. The ordinary PME of post 1 reaches only the first three columns, and needs eight measurements to do it. A single measurement of all sixteen is available on the rotating side alone.

## 3. Turning Two of Them at Different Speeds

Why will one not do? Azzam's argument is compact. The detected signal comes down to one line.

$$I = c\,\mathbf{A}\,M\,\mathbf{P}$$

Here $\mathbf{P}$ is the Stokes vector leaving the polarizing optics and $\mathbf{A}$ is the first row of the Mueller matrix of the analyzing optics. Expanded, it reads

$$I = c\sum_{i,j} \mu_{ij}\,m_{ij}, \qquad \mu_{ij} = a_i\,p_j$$

The weight $\mu_{ij}$ decides how much each Mueller element contributes to the signal, and the fact that this weight is a product of $a_i$ and $p_j$ explains everything. Modulate only the generator and only $p_j$ varies in time, so $\mu_{ij}$ sweeps in the $j$ direction alone. Whole rows never get sampled. Both sides must be driven for $\mu_{ij}$ to sweep all sixteen positions.

The number of frequencies required falls out here too. Solving for sixteen elements needs sixteen independent non-zero Fourier amplitudes, and since each frequency contributes one cosine and one sine amplitude, at least eight distinct frequencies are necessary.

The specific instrument Azzam proposed rotates two quarter-wave plates together at $\omega$ and $5\omega$. The result is a periodic signal whose fundamental is $2\omega$, carrying harmonics up to twelfth order on top of a dc term — twenty-five amplitudes in all.

Not all twenty-five carry information. Four of them vanish identically, and five identities hold among the rest. Subtract four and then five from twenty-five and exactly sixteen remain, matching the number of unknowns.

The identities look like this: $a_7 = -a_3$, $b_7 = b_3$, $a_6 = -a_4$, $a_{11} = -a_9$, $b_{11} = -b_9$. Values read at different harmonics are tied together this way.

Recovering the elements proceeds in three layers. The five in the fourth row and column each follow from one amplitude directly. The four in the central $2\times2$ come from linear combinations of two amplitudes each. The remainder, in the first row and column, require substituting the values already obtained. The calculation is deliberately staged: rather than inverting one matrix whole, Azzam breaks it into small pieces.

All of that gives the normalized matrix. Absolute values need the proportionality constant separately, obtained by removing the sample and running the instrument straight through, where the Mueller matrix is the identity. Additional amplitudes vanish in that configuration, so simply checking whether they are zero amounts to an instrument check.

That Azzam did not leave the five identities idle matters. Departures from them are symptoms of azimuth errors, non-ideal components or detector nonlinearity, and he noted that they can be used to characterize exactly those faults. The surplus coefficients are not waste but a diagnostic built into the instrument.

Several remarks Azzam appended also touch this series directly. One is that a fixed analyzer at the end of the optical train avoids errors from polarization-dependent detector sensitivity. Section 5 of post 2 discussed why a depolarizer goes in front of the detector; here the layout itself solves the same problem.

The second remark matters more. Scanning wavelength requires the quarter-wave plates to be achromatic or tunable at every wavelength. Arbitrary and differing retardances still permit the Mueller matrix to be determined, he added, but the mathematics become complicated. The whole of post 2 hangs on that sentence, and the commercial design in the next section takes precisely the complicated route.

The last remark points ahead in the series. Replacing mechanical rotation with the optical rotation of a pair of Faraday cells would give an instrument with no moving parts at all. A 1978 paper was already looking toward doing away with rotation.

## 4. Which Element Rides on Which Harmonic

Collins and Koh set out the design of the commercial instrument in 1999. In the notation of post 1 the configuration is $PC_RSC_RA$; since the two compensators must be told apart, it is written here with numbers as $PC_{1r}(\omega_1)\,S\,C_{2r}(\omega_2)\,A$. A fixed polarizer and analyzer sit symmetrically about the sample. The rotation rates are $\omega_1 = 5\omega$ and $\omega_2 = 3\omega$, and $\pi/\omega$ is the fundamental optical period.

The detected intensity reads as follows, with $C$ the base rotation angle.

$$I(C) = I_0\Big\{1 + \sum_{n=1}^{16}\big[\alpha_{2n}\cos(2nC-\phi_{2n}) + \beta_{2n}\sin(2nC-\phi_{2n})\big]\Big\}$$

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig2-waveform-spectrum.png" alt="Waveform and harmonics for RCE and DRC-MME" width="780">
_Figure 2. The same sample seen with one compensator and with two. Four bars are missing in the lower spectrum, and the calculation puts the coefficients at exactly those positions to zero._

Figure 2 places this beside the RCE of post 1. One compensator gives two harmonics; two compensators raise twenty-four, from $2C$ up to $32C$. The waveform itself grows visibly more intricate.

Four of the sixteen positions are empty, however. The eight coefficients at $n$ equal to 9, 12, 14 and 15 vanish. Multiplying the Mueller matrices directly puts those positions at zero to within $10^{-17}$. The twenty-four that remain determine fifteen normalized Mueller elements, which leaves more equations than unknowns. The system is overdetermined.

There is a clear pattern to which element rides on which harmonic. A few entries give the flavour. Only $M_{44}$ appears at $4C$ and $16C$; the three elements of the fourth row gather at $6C$, the three of the fourth column at $10C$. The central $2\times2$ is spread over $8C$, $12C$, $20C$ and $32C$. That the same element appears at several harmonics is what overdetermination looks like in practice.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig3-coefficient-map.png" alt="Mueller elements carried by each harmonic" width="780">
_Figure 3. The elements carried at each harmonic and the factor that sets their strength. The red rows are the ones riding on $\sin\delta$._

Write $\delta_1$ and $\delta_2$ for the two retardances, with $s_j = \sin^2(\delta_j/2)$ and $c_j = \cos^2(\delta_j/2)$. Expanding what the generator contributes to the sample exposes the structure.

$$K_j = M_{j1} + \big[c_1\cos 2P' + s_1\cos(4C'_1-2P')\big]M_{j2} + \big[c_1\sin 2P' + s_1\sin(4C'_1-2P')\big]M_{j3} + \big[\sin\delta_1\sin(2C'_1-2P')\big]M_{j4}$$

The second and third columns are carried partly by a dc share and partly by a share modulated at $4C'_1$. Since $c_1 + s_1 = 1$, a larger retardance shifts weight from the dc term to the modulated one. The fourth column behaves differently. It rides on $\sin\delta_1$ alone, and only at $2C'_1$.

Section 3 of post 2 said that a retardance of 0° or 180° destroys the information. This is what that statement means precisely. At either value $\sin\delta_1$ is zero and the fourth column disappears wholesale, taking the sensitivity to circular polarization with it.

The calculation makes it plain. Set the first retardance to 0° or 180° and a sample with its fourth column zeroed out produces exactly the same signal as the original. Do the same to the second compensator and the fourth row vanishes in the same way. A constraint discussed in post 2 as a matter of choosing parts here trims the list of measurable quantities directly.

## 5. Why Five to Three

Why is the rotation ratio five to three? The answer is unexpected: hardware decided it, not mathematics.

The paper gives two considerations. There must be enough independent information to extract the fifteen elements. And the fundamental optical period must be long enough to permit the detector scans needed to reach the highest-order coefficient, yet short enough that continuously rotating motors do not become unstable at low frequency.

Then it states plainly that the combination giving all Mueller elements with the lowest highest-order Fourier coefficient is one to five. Azzam's original choice is optimal by that criterion, and the highest order is then $24C$.

Reality intervenes. The minimum scan time of a commercial photodiode array is 5.1 ms, and running the two motors at 25 Hz and 5 Hz in a one-to-five scheme leaves too little time to integrate. Dropping them to 17.5 Hz and 3.5 Hz gives 5.71 ms and satisfies that condition, but now the motor becomes unstable at 3.5 Hz.

The way out is to raise the slower ratio. The highest-order coefficient goes up, but the slow motor can run faster. Five to three is where the penalty in harmonic order meets the gain in motor stability.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig4-ratio-tradeoff.png" alt="What the rotation ratio decides" width="780">
_Figure 4. Left: the highest harmonic for each ratio and the integrations needed to extract it. Right: the condition numbers for those same ratios, where the two curves overlap into one._

The left panel of figure 4 is that trade. Computed directly, the highest harmonic comes out at $24C$ for one to five, $28C$ for five to two, $32C$ for five to three and $36C$ for five to four, matching the paper's table exactly. Extracting every coefficient takes one integration more than the highest order, so five to three needs 33 per optical period.

The design values follow. A fundamental frequency of 2.5 Hz gives a 200 ms optical period, and the two motors turn at 12.5 Hz and 7.5 Hz. Reading 36 times per period allows 5.56 ms each, which satisfies the detector.

How the reads are timed is part of the design too. The detector fires every 5° of base rotation, which is every 25° of the first compensator and every 15° of the second. A single 54 kHz clock divided by 3 drives the first motor and divided by 5 the second, while division by 300 and 7200 produces the encoder pulses and the reference trigger. Both motors and the detector branch off one clock, so nothing can drift apart.

Turning the compensators in the same direction or in opposite directions makes no difference. Counter-rotation only flips the sign of the terms belonging to the third and fourth columns, and reversing those signs recovers the same information.

Section 9 of post 1 said the number of samples needed is set by the number of unknown coefficients: three intervals for RAE, five for RCE. Here it is 33. The same principle, an order of magnitude larger.

## 6. The Ratio Does Not Set the Precision

If five to three was chosen for hardware reasons, was precision given up? The question is worth checking.

How far measurement noise is amplified into the Mueller elements is what the condition number of the data reduction matrix reports. Since the detected intensity is linear in the Mueller elements, feeding in unit matrices one at a time builds that matrix column by column.

The answer is unambiguous. Five to three and one to five have identical condition numbers across the whole range of retardance. That is why the two curves in the right panel of figure 4 merge into one.

The ratio plays no part in noise amplification. What it changes is the highest harmonic — how many times per period the detector must be read.

The result was unexpected at first. Given that five to three is a compromise forced by the detector and the motors, as section 5 showed, it is natural to assume that abandoning the mathematically optimal one to five costs something somewhere. The only cost turns out to be reading 33 times instead of 25. Not a digit of precision is lost. Collins and Koh were right to treat the ratio purely as a timing question.

This also fits what Smith obtained in 2002 by minimizing the condition number. With exactly sixteen measurements the relationship between condition number and angular increments is intricate, but once the measurement is overdetermined a broad range of increments performs essentially identically. Overdetermination relieves the burden of choosing a ratio.

So what does set the precision?

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig5-condition-number.png" alt="Condition number against retardance" width="780">
_Figure 5. The condition number of the data reduction matrix against compensator retardance. A quarter-wave plate is not the best choice._

The retardance does. Figure 5 sweeps it, and the minimum sits at 128°, not at 90°. A quarter-wave plate gives a condition number of 9.45; a 128° plate gives 2.49. For the same detector noise the error in the Mueller elements is 3.8 times smaller. Smith's optimum of 127° is effectively the same answer.

Why not 90° can be read off the structure in section 4. The fourth row and column ride on $\sin\delta$, which peaks at 90°. The modulated share of the second and third columns rides on $s = \sin^2(\delta/2)$, which keeps growing to 180°. The compromise between the two is pushed above 90°. Maximizing $\sin\delta\cdot\sin^2(\delta/2)$ alone gives 120°; accounting for the condition number of the full matrix lands near 128°.

Section 3 of post 2 said the retardance need not be exactly 90°. Here that sentence acquires a number. It is not merely unnecessary — 90° is a loss.

## 7. The Real Cost Is Speed

Where is the price paid for getting sixteen at once? In time.

The optical period at five to three is 200 ms, so that is how long one Mueller matrix takes. The minimum measurement time for RAE and RCE in the table of post 1 was about 10 ms, making this twenty times slower.

The chain is easy to follow. Two compensators push the highest harmonic from $24C$ to $32C$; a higher harmonic demands at least 33 integrations per period; each detector read costs over 5 ms; the period therefore stretches to 200 ms. Touch any link and the rest move with it.

What 200 ms means depends on the application. It is fast enough to follow film growth in real time, which was the original motivation for the instrument. But against the 20 μs with which the photoelastic modulator of post 1 tracked liquid crystal molecules, it is four orders of magnitude away. The price of all sixteen is paid in time.

The exchange from post 1 changes shape once more here. Null ellipsometry spent time to buy accuracy; photometric methods trusted the detector to save time. Mueller matrix measurement spends time again, this time to buy information. What is sought changes; what is surrendered sits in the same place.

## 8. Two Routes Back to Speed

Having established that it is slow, the next question is how to make it fast. Two approaches diverge.

One pushes on the hardware. The high-speed instrument Collins's group reported in 2001 took that path, refining multichannel detection and rotation drive together to gain speed on the same principle. The design table of section 5 already points that way: the ratio, the highest harmonic, the integration time and the motor frequency are linked in one chain, so a faster detector lets the whole chain be tightened. The ceiling of this approach is set by the detector.

The other redesigns when to read. The approach published in 2019 starts from the observation that there is no reason for sampling points in a continuous measurement to be evenly spaced. Building an objective that tolerates both systematic error and detection noise, and searching for optimal unevenly spaced points with a multi-objective genetic algorithm, reaches a globally minimal number of samples while improving error immunity rather than degrading it.

The contrast is sharp. The first wants to read faster; the second wants to read less.

The second is possible because of what section 6 established. Twenty-four coefficients for fifteen elements left room to cut samples from the start, and since the ratio does not touch the condition number, where the samples sit can be rearranged freely. Even spacing was never a requirement of the physics — it was a convenience imposed to keep the Fourier analysis simple.

This also connects to the later parts of the series. Raising the rotation speed eventually runs into the motor and the detector. What remains is to change how the reading is done, or to remove the rotation altogether. Azzam's 1978 mention of Faraday cells is one end of that thought; the instruments of posts 5 and 6 are the other.

## 9. What the Surplus Measurements Are For

Overdetermination has come up several times. Seeing what the surplus information actually does reveals the character of the instrument.

Thirty-six quantities are read per optical period. The signal contains fewer harmonics than that, so eleven of them ought to be zero from the start. If they are not, something in the instrument has shifted. The diagnostic Azzam performed with identities becomes, here, a matter of counting the coefficients that should have vanished.

The real use is calibration. Four things need calibrating: the two retardances, the two compensator phase offsets, and the angular offsets of the polarizer and analyzer. The order in which they are obtained is fixed.

Retardance should not depend on how the sample was mounted, so it is determined separately in a dedicated high-accuracy measurement. The phase offsets and angular offsets, by contrast, must be obtained in the actual configuration with the sample in place, because they change when the alignment does.

The method falls out of the coefficient structure. Squaring and summing the cosine and sine coefficients of a harmonic cancels the phase angle, which yields the polarizer and analyzer offsets; taking their ratio instead leaves the phase angle intact, which gives the compensator phase offsets.

Those phase angles are not spread evenly over the two compensators. Only the first compensator's phase appears at $10C$, only the second's at $12C$. A window onto each compensator separately is already open in the signal. At the other harmonics the two enter mixed as sums and differences.

Several routes to the same quantity therefore exist. For one combination of polarizer offset and first-compensator phase alone there are five different paths through different harmonics. Agreement among the five says the instrument is sound; divergence narrows down what has shifted.

Overdetermination is not waste. Twenty-four coefficients solving for fifteen elements leave nine over, and those nine watch the instrument. They were not kept to gain information but to make the instrument trustworthy.

## Summary

Three things came out of this post.

First, two compensators are needed because of the structure of the weight $\mu_{ij} = a_i p_j$. Driving one side moves the weight in one direction only, which never sweeps all sixteen positions. When a sample is anisotropic and depolarizing the Jones formalism breaks down, and from that point all sixteen elements are required.

Second, the rotation ratio sets timing, not precision. One to five gives the lowest highest-order harmonic, but the minimum detector scan time and the low-frequency instability of the motors pushed the design to five to three — and the condition number is identical either way. What sets the precision is the retardance, and its optimum is near 128° rather than a quarter wave.

Third, the instrument carries the material for checking itself inside its own signal. Azzam's five identities and the five calibration routes of Collins and Koh come from the same property. As section 9 showed, it is the surplus measurements that do that work.

The next post follows from there. Going beyond using surplus measurements for instrument checks leads to assuming the components are not ideal and extracting their non-ideality from the measurement as well. Post 4 takes up the calibration of rotating-element ellipsometers — the substance of what post 2 called using cheap, non-ideal parts and making up the difference in calibration.

## Reproducing the Calculations

Every calculation here comes from products of component Mueller matrices. The code is in `_code/ellipsometer-dual-rotating-compensator-mme/`.

- `mme.py` — detected intensity for the $PC_{1r}SC_{2r}A$ configuration, Fourier coefficients, the data reduction matrix and its condition number
- `verify_mme.py` — the cross-check. Nothing is computed from a closed form; the Mueller product alone produces the result, which is then compared against the coefficient table of Collins and Koh
- `generate_figures.py` — produces the five figures

Four things were confirmed. The vanishing harmonics fall exactly where the paper says; the highest harmonic for four different rotation ratios matches the paper's table; the fourth row and column contribute nothing to the signal when the retardance is 0° or 180°; and setting the first retardance to zero reduces the configuration to the rotating-compensator layout of post 1. Comparison against the coefficient table agreed to within $10^{-16}$ for three different combinations of retardance and component angle.

That comparison failed completely at first. The coefficients in the paper's table are normalized by the dc component, and only the numerators were being compared. Matching the normalization brought the two into agreement at machine precision.

## References

- R. M. A. Azzam, "Photopolarimetric measurement of the Mueller matrix by Fourier analysis of a single detected signal," *Optics Letters* **2**, 148 (1978) — the origin of the dual rotating retarder.
- R. W. Collins and J. Koh, "Dual rotating-compensator multichannel ellipsometer: instrument design for real-time Mueller matrix spectroscopy of surfaces and films," *J. Opt. Soc. Am. A* **16**, 1997 (1999) — the 5:3 design and the coefficient table.
- M. H. Smith, "Optimization of a dual-rotating-retarder Mueller matrix polarimeter," *Applied Optics* **41**, 2488 (2002) — condition-number minimization and the optimal retardance.
- J. Lee, J. Koh, R. W. Collins, "Dual rotating-compensator multichannel ellipsometer: instrument development for high-speed Mueller matrix spectroscopy of surfaces and thin films," *Rev. Sci. Instrum.* **72**, 1742 (2001).
- K. Meng, B. Jiang, C. D. Samolis, M. Alrished, K. Youcef-Toumi, "Unevenly spaced continuous measurement approach for dual rotating-retarder Mueller matrix ellipsometry," *Optics Express* **27**, 14736 (2019).
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, section 4.2.7 (configurations for Mueller matrix ellipsometry).
- Youngjoon Kim, "Development of Snapshot Angle-Resolved Spectroscopic Ellipsometry Using a Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD thesis, Seoul National University, 2025, section 1.2.

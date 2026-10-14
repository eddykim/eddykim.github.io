---
date: 2026-10-10 20:00:00 +0900
title: "Ellipsometer Instrumentation 3 — Why Spin Two Compensators at Different Speeds?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-dual-rotating-compensator-mme/
page_id: ellipsometer-dual-rotating-compensator-mme
categories: [Optics, Instrumentation]
tags: [ellipsometry, mueller-matrix, dual-rotating-compensator, fourier-analysis, instrumentation]
description: Modulating only the polarization state generator or only the analyzer cannot deliver all sixteen elements at once. This post follows why two compensators must turn at different speeds.
math: true
---

[Post 1](/en/posts/ellipsometer-oblique-rotating-element/) sorted the methods by what they cannot measure. RAE misses $S_3$ and loses the sign of $\Delta$; a compensator recovers it at the cost of a second run; rotating that compensator settles it in one. [Post 2](/en/posts/ellipsometer-components/) followed where each of those components stops working across a spectral band.

This post asks something one level deeper. Every instrument in post 1 compressed the sample into two numbers, $\Psi$ and $\Delta$. What happens when that compression does not hold?

The short answer is that the whole Mueller matrix has to be measured, all sixteen elements of it, and that one compensator is not enough to do it. Azzam stated the obstruction outright in the opening of a 1978 paper.

> The simultaneous determination of all 16 elements of M by Fourier analysis of one signal is not possible if either the polarizing or the analyzing optics alone are modulated.

## 1. Why Must All Sixteen Be Measured?

When can a sample not be summarized by $\Psi$ and $\Delta$? Counting the independent parameters the sample demands gives the answer.

The Mueller matrix of an isotropic sample is fixed by three values, $S_1$, $S_2$ and $S_3$: the same three that [section 2 of post 1](/en/posts/ellipsometer-oblique-rotating-element/) wrote in terms of $\Psi$ and $\Delta$. Without depolarization their squares sum to one, so only two of them, $\Psi$ and $\Delta$, are actually free. The RCE of post 1 obtains all three in a single measurement, so nothing is lacking.

An anisotropic sample needs the four complex elements of a Jones matrix. That is eight real numbers; the absolute phase cannot be measured and the intensity is normalized, which leaves six. Existing instruments can still fill these by measuring at several component angles.

The difficulty arrives when an anisotropic sample also depolarizes. At least seven parameters are then required, and the Jones matrix cannot describe the sample at all. The limitation of the Jones formalism to fully polarized light, discussed in [post 2 of the ellipsometry foundations series](/en/posts/ellipsometry-polarization-mueller-matrix/), becomes a hardware requirement here. Faced with an anisotropic, depolarizing sample, the two-number summary breaks down, and the only remaining route is to measure the Mueller matrix directly.

## 2. Which Elements Can Each Rotating Configuration Reach?

Attempts to measure a Mueller matrix with rotating elements began in the late 1970s. Hauge's 1980 survey of what each configuration reaches is the skeleton of this post.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig1-mme-configurations.png" alt="Mueller elements reachable by each rotating-element configuration" width="780">
_Figure 1. Four rotating-element configurations and the Mueller elements each one measures. In each 4×4 grid, red circles are measured elements and small grey dots are elements that are not. Only the last configuration, with two compensators turning at different speeds, fills all sixteen._

Rotating the polarizer and analyzer together misses the fourth row and the fourth column. With no component to produce or to read circular polarization, that is only to be expected. Adding a single compensator improves matters, and where it goes decides the outcome. Placed after the sample it yields the first three columns; placed before it, the first three rows. A compensator after the sample lets the analyzing side read circular polarization, which opens the fourth row. Nothing before the sample can produce circular polarization yet, so the fourth column stays closed. Placing it before the sample does the reverse.

What deserves attention is that the configuration with the compensator after the sample has exactly the optical layout of post 1's RAE-with-compensator ($PSCA_R$) and RCE ($PSC_RA$). Section 2 of that post argued that identical components give different reach depending on which one turns. The same argument reappears one level up, at the Mueller matrix. Where a component sits decides whether rows or columns come back.

Three rows or three columns do not make sixteen. Among rotating-element configurations, only two compensators turning at different frequencies deliver the whole matrix in a single measurement. The design was proposed in the late 1970s and built in the early 1990s. Collins's group reported an instrument capable of real-time spectroscopic measurement in 2000 (Lee et al. 2000).

Photoelastic modulators can also measure a Mueller matrix. Two modulators at different resonant frequencies give nine elements, and four measurements at different component angles fill the rest. The ordinary PME of post 1 reaches only the first three columns, and needs eight measurements to do it. Among the configurations Fujiwara (2007) surveys, only the rotating type gives all sixteen in one measurement. Doing it without rotation takes four modulators (electro-optic in Thompson et al. 1980, photoelastic in Arteaga et al. 2012), or moving the modulation onto the wavelength axis as in the channeled approach of post 5.

## 3. Why Must the Two Compensators Turn at Different Speeds?

The limit that figure 1 showed is no accident. Azzam's argument explains it compactly. The detected signal comes down to one line.

$$I = c\,\mathbf{A}\,M\,\mathbf{P}$$

Here $\mathbf{P}$ is the Stokes vector leaving the polarization state generator and $\mathbf{A}$ is the first row of the Mueller matrix of the polarization state analyzer. Expanded, it reads

$$I = c\sum_{i,j} \mu_{ij}\,m_{ij}, \qquad \mu_{ij} = a_i\,p_j$$

The weight $\mu_{ij}$ decides how much each Mueller element contributes to the signal, and the fact that this weight is a product of $a_i$ and $p_j$ explains everything. Modulate only the generator and only $p_j$ varies in time. The signal can then tell columns $j$ apart, but the four elements of a column always come out mixed in the fixed proportions $a_i$. Four bundles come back, not sixteen elements. Both sides must be driven for $\mu_{ij}$ to sweep all sixteen positions.

The number of frequencies required falls out here too. Solving for sixteen elements needs sixteen independent non-zero Fourier amplitudes, and since each frequency contributes one cosine and one sine amplitude, at least eight distinct frequencies are necessary.

Driving both sides is not enough on its own. If the two compensators turn at the same speed, the frequencies they generate pile onto the same positions, and different elements fall into the same harmonic. Computed, equal speeds leave only eight of the sixteen combinations independent. Section 4 returns to how the speeds must differ, in terms of harmonic "addresses".

The specific instrument Azzam proposed rotates two quarter-wave plates together, at $\omega$ in the generator and $5\omega$ in the analyzer. This post writes rotation ratios as fast:slow, so this is 5:1. The result is a periodic signal whose fundamental is $2\omega$, carrying harmonics up to twelfth order on top of a dc term, twenty-five amplitudes in all. The twelfth order comes from the fact that a rotating compensator modulates the signal up to four times its own rotation frequency. The highest frequency is $4\omega + 20\omega = 24\omega$, twelve times the $2\omega$ fundamental.

Not all twenty-five carry information. Four of them vanish identically, and five identities hold among the rest. Subtract four and then five from twenty-five and exactly sixteen remain, matching the number of unknowns.

Writing $a_n$ and $b_n$ for the cosine and sine amplitudes of the $n$th harmonic (a different symbol from the analyzer components $a_i$ above), the identities look like this: $a_7 = -a_3$, $b_7 = b_3$, $a_6 = -a_4$, $a_{11} = -a_9$, $b_{11} = -b_9$. Values read at different harmonics are tied together, so if one drifts its partner must drift with it. Azzam noted that departures from these identities are symptoms of azimuth errors, non-ideal components or detector nonlinearity, and can be used to characterize those faults. The surplus coefficients are a diagnostic built into the instrument.

Recovering the elements proceeds in three layers. The five in the fourth row and column each follow from one amplitude directly. The four in the central $2\times2$ come from linear combinations of two amplitudes each. The remainder, in the first row and column, require substituting the values already obtained. Rather than inverting one large matrix at once, Azzam breaks it into small pieces.

All of that gives the normalized matrix. Absolute values need the proportionality constant separately, obtained by removing the sample and running the instrument straight through, where the Mueller matrix is the identity. Additional amplitudes vanish in that configuration, so simply checking whether they are zero amounts to an instrument check.

Several remarks Azzam appended also touch this series directly. One is that a fixed analyzer at the end of the optical train avoids errors from polarization-dependent detector sensitivity. Section 5 of post 2 discussed why a depolarizer goes in front of the detector; here the layout itself solves the same problem.

The second remark matters more. Scanning wavelength requires the quarter-wave plates to be achromatic or tunable at every wavelength. Arbitrary and differing retardances still permit the Mueller matrix to be determined, he added, but the mathematics become complicated. The whole of post 2 hangs on that sentence, and the Collins–Koh design in the next section takes precisely the complicated route.

The last remark points ahead in the series. Replacing mechanical rotation with the optical rotation of a pair of Faraday cells would give an instrument with no moving parts at all. A 1978 paper was already looking toward doing away with rotation.

## 4. Which Element Rides on Which Harmonic?

Collins and Koh proposed the design of a real-time instrument in 1999. This post calls the method dual rotating compensator (DRC) Mueller matrix ellipsometry (MME), or DRC-MME. In the notation of post 1 the configuration is $PC_RSC_RA$; since the two compensators must be told apart, it is written here with numbers as $PC_{1r}(\omega_1)\,S\,C_{2r}(\omega_2)\,A$. A fixed polarizer and analyzer sit symmetrically about the sample. The rotation rates are $\omega_1 = 5\omega$ and $\omega_2 = 3\omega$ (ratio 5:3), and $\pi/\omega$ is the fundamental optical period. Here $\omega$ is the fundamental mechanical frequency, unlike post 1, where $\omega$ was the compensator's own rotation frequency.

The detected intensity reads as follows.

$$I(C) = I_0\Big\{1 + \sum_{n=1}^{16}\big[\alpha_{2n}\cos(2nC-\phi_{2n}) + \beta_{2n}\sin(2nC-\phi_{2n})\big]\Big\}$$

Here $C = \omega t$ is the base rotation angle and $I_0$ the dc intensity. $\alpha_{2n}$ and $\beta_{2n}$ are the normalized cosine and sine coefficients of the $2nC$ harmonic, and $\phi_{2n}$ is a phase set by the offset between each compensator's scale and its true azimuth, which section 9 determines.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig2-waveform-spectrum.png" alt="Waveform and harmonics for RCE and DRC-MME" width="780">
_Figure 2. Waveform and harmonics for the same sample matrix (a test matrix with randomly chosen elements) seen with one compensator (top) and with two (bottom). The top case sets the first retardance to zero, reducing the layout to the RCE of post 1; the compensator turning at 3C produces only 6C and 12C, twice and four times its rate. The bottom case turns two 90° compensators at 5:3, and the coefficients are zero at 18C, 24C, 28C and 30C (n = 9, 12, 14, 15)._

Figure 2 places this beside the RCE of post 1. One compensator gives two harmonics; two compensators raise twelve, from $2C$ up to $32C$, or twenty-four coefficients counting cosine and sine. The waveform itself grows visibly more intricate.

Four of the sixteen positions are empty, however. The eight coefficients at $n$ equal to 9, 12, 14 and 15 vanish. Multiplying the Mueller matrices directly puts those positions at zero to within $10^{-17}$. The same calculation agrees with the coefficient table of Collins and Koh to within $10^{-16}$ for three combinations of retardance and polarizer angle, once both are divided by the dc term as the paper's table is. The twenty-four that remain determine fifteen normalized Mueller elements, which leaves more equations than unknowns. The system is overdetermined.

Why those four? There is a rule for which element rides on which harmonic, and the rule answers the question. The first compensator, turning at 5C, modulates the light entering the sample and so carries the columns of the Mueller matrix: the fourth column at 10C, the second and third at dc and 20C. The second compensator, turning at 3C, modulates the light leaving the sample and so selects the rows: the fourth row at 6C, the second and third at dc and 12C. The detector receives the product of the two modulations, so each element rides on sums and differences of its column frequency and its row frequency.

| | Column 1 (dc) | Columns 2–3 (dc, 20C) | Column 4 (10C) |
|---|---|---|---|
| Row 1 (dc) | $M_{11}$: dc | $M_{12}, M_{13}$: 20C | $M_{14}$: 10C |
| Rows 2–3 (dc, 12C) | $M_{21}, M_{31}$: 12C | central $2\times2$: 8C, 12C, 20C, 32C | $M_{24}, M_{34}$: 2C, 10C, 22C |
| Row 4 (6C) | $M_{41}$: 6C | $M_{42}, M_{43}$: 6C, 14C, 26C | $M_{44}$: 4C, 16C |

Each cell is the sum and difference of its row and column frequencies. The central $2\times2$ rides on $20 \pm 12$, that is 8C and 32C, plus 12C and 20C when one side is dc. $M_{44}$ sits at $10 \pm 6$, 4C and 16C. Some entries drop to zero depending on the polarizer and analyzer azimuths, and the dc shares are listed only for $M_{11}$. Collecting every frequency that can be formed gives twelve: 2, 4, 6, 8, 10, 12, 14, 16, 20, 22, 26 and 32C. No sum or difference makes 18C, 24C, 28C or 30C. That is why $n$ = 9, 12, 14 and 15 are empty. The highest order follows from the same rule: 20C plus 12C is 32C, and for a ratio $m_1:m_2$ the highest order is $4(m_1+m_2)C$.

Each harmonic, in other words, is a broadcast frequency assigned to particular Mueller elements. That the same element is carried on several frequencies is what overdetermination looks like in practice.

The rule also shows why the speeds must differ. Two stations on one frequency are heard mixed. If the sums and differences of the column frequencies $$\{0, 2m_1, 4m_1\}$$ and the row frequencies $$\{0, 2m_2, 4m_2\}$$ coincide, different elements fall into the same harmonic and cannot be separated. Computed, 5:1, 5:2, 5:3 and 5:4 have no coinciding positions. Ratios 2:1, 3:1, 4:1 and 3:2 do, and cannot resolve all sixteen. Ratio 4:3 resolves them, but $16-12 = 12-8 = 4$ puts two addresses at 4C, and the condition number of section 6 rises to 12.8 against 9.45 for the collision-free ratios (at 90° retardance). With the fast side at 5, any slow side from 1 to 4 keeps the information unmixed.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig3-coefficient-map.png" alt="Mueller elements carried by each harmonic" width="780">
_Figure 3. For the 5:3 layout, the elements carried at each harmonic and the representative factor that sets their strength. Red rows ride on $\sin\delta$ and disappear when a retardance is 0° or 180°; blue rows ride on $s = \sin^2(\delta/2)$. Terms multiplied by the dc share $$c_j$$ are omitted._

Figure 3 attaches strength factors to the same map. Write $\delta_1$ and $\delta_2$ for the two retardances, with $s_j = \sin^2(\delta_j/2)$ and $c_j = \cos^2(\delta_j/2)$. Expanding what the generator contributes to the sample shows where the factors come from.

$$K_j = M_{j1} + \big[c_1\cos 2P' + s_1\cos(4C'_1-2P')\big]M_{j2} + \big[c_1\sin 2P' + s_1\sin(4C'_1-2P')\big]M_{j3} + \big[\sin\delta_1\sin(2C'_1-2P')\big]M_{j4}$$

Here $j$ is the row index of the Mueller matrix. $$P'$$ and $$C'_1$$ are the true angles after subtracting the offsets found by calibration, with $$C'_1 = 5(C - C_{S1})$$.

The second and third columns are carried partly by a dc share and partly by a share modulated at $$4C'_1$$. Since $c_1 + s_1 = 1$, a larger retardance shifts weight from the dc term to the modulated one. The fourth column behaves differently. It rides on $\sin\delta_1$ alone, and only at $$2C'_1$$.

Section 3 of post 2 said that a retardance of 0° or 180° destroys the information. This is what that statement means precisely. At either value $\sin\delta_1$ is zero and the fourth column disappears wholesale, taking the sensitivity to circular polarization with it. The calculation makes it plain. Set the first retardance to 0° or 180° and a sample with its fourth column zeroed out produces exactly the same signal as the original. Do the same to the second compensator and the fourth row vanishes in the same way. A constraint discussed in post 2 as a matter of choosing parts here trims the list of measurable quantities directly.

## 5. Why 5:3?

Why is the rotation ratio 5:3? As section 4 showed, with the fast side at 5 any slow side from 1 to 4 works. Hardware, not mathematics, picked one of the four.

The paper gives two considerations. There must be enough independent information to extract the fifteen elements. And the fundamental optical period must be long enough to permit the detector scans needed to reach the highest-order coefficient, yet short enough that continuously rotating motors do not become unstable at low frequency.

The paper states plainly that the combination giving all Mueller elements with the lowest highest-order Fourier coefficient is 5:1. Azzam's original choice is optimal by that criterion, and the highest order is then $24C$.

Reality intervenes. The minimum scan time of a commercial photodiode array is 5.1 ms, and running the two motors at 25 Hz and 5 Hz at 5:1 leaves too little time to integrate. Dropping them to 17.5 Hz and 3.5 Hz gives 5.71 ms and satisfies that condition, but the paper expected the motor to become unstable at 3.5 Hz.

A photodiode array is a detector that lines up its sensing elements in a row and reads them one after another (figure 4). Behind a spectrograph each element takes one wavelength, so many wavelengths are measured at once. Reading the whole row, though, takes a fixed time. The 5.1 ms is that time.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/ext-linear-photodiode-array.jpg" alt="Linear photodiode array chip for multichannel spectroscopy" width="700">
_Figure 4. A self-scanning linear photodiode array made for multichannel spectroscopy (Hamamatsu S3901-256Q). Photodiodes sit in a single row at 50 μm spacing along the long purple window, and one scan ends only when the whole row has been read. (Source: Mister rf, [S3901-256Q](https://commons.wikimedia.org/wiki/File:S3901-256Q.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

The way out is to raise the slow side. Keeping the fast side at 5 and raising the slow side from 1 to 4 lifts the highest harmonic to $4(5+m_2)C$. How many reads that takes follows from counting unknowns. At 5:3 the harmonics run to sixteenth order ($32C$), so the unknowns are one dc term plus sixteen cosine–sine pairs, 33 in all, and the detector must be read at least 33 times per optical period. In general $N = 4(m_1+m_2)+1$. More reads lengthen the period, but under the same detector limit they let the slow motor run faster.

| Ratio | Highest harmonic | Minimum reads | Minimum optical period (reads × 5.1 ms) | Slow-motor ceiling |
|---|---|---|---|---|
| 5:1 | 24C | 25 | 128 ms | 3.9 Hz |
| 5:2 | 28C | 29 | 148 ms | 6.8 Hz |
| 5:3 | 32C | 33 | 168 ms | 8.9 Hz |
| 5:4 | 36C | 37 | 189 ms | 10.6 Hz |

The last column is the trade. At 5:1 the slow motor cannot exceed 3.9 Hz; at 5:3 it can reach 8.9 Hz. The price is a minimum optical period stretched from 128 ms to 168 ms.

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig4-ratio-tradeoff.png" alt="What the rotation ratio decides" width="780">
_Figure 5. Left: the highest harmonic for each ratio and the integrations needed to extract it. Right: the condition number of the data reduction matrix against retardance for 5:3 and 5:1; the two curves overlap into one._

The left panel of figure 5 is that trade. Computed directly, the highest harmonic comes out at $24C$ for 5:1, $28C$ for 5:2, $32C$ for 5:3 and $36C$ for 5:4, matching the paper's table exactly.

The design values follow. A fundamental mechanical frequency $\omega/2\pi$ of 2.5 Hz gives a 200 ms optical period, and the two motors turn at 12.5 Hz and 7.5 Hz. Reading 36 times per period allows 5.56 ms each, which satisfies the detector.

Read against those 36 reads, the table shows why 5:3 is the end of the line. 5:4 needs 37 and does not fit; 5:3 needs 33 and does. So 5:3 is the largest slow side that fits 36 reads per period. This is arithmetic on the paper's numbers, though. Collins and Koh do not say that they fixed 36 reads first and chose the ratio from it.

How the reads are timed is part of the design too. The detector fires every 5° of base rotation, which is every 25° of the first compensator and every 15° of the second. A single 54 kHz clock divided by 3 drives the first motor and divided by 5 the second, while division by 300 and 7200 produces the encoder pulses and the reference trigger. Both motors and the detector branch off one clock, so nothing can drift apart.

Turning the compensators in the same direction or in opposite directions yields the same information. Counter-rotation flips the sign of the terms belonging to the third and fourth columns. If the polarizer is not at 0°, its angle must be read with the opposite sign as well. With those signs reversed, the same information comes back.

Section 7 of post 1 said the number of samples needed is set by the number of unknown coefficients: three intervals for RAE, five for RCE. Here it is 33. The same principle, an order of magnitude larger.

## 6. Does Choosing 5:3 Cost Precision?

If 5:3 was chosen for hardware reasons, was precision given up? The question is worth checking.

How far measurement noise is amplified into the Mueller elements is what the condition number of the data reduction matrix reports. The data reduction matrix maps the sixteen Mueller elements onto the Fourier coefficients. Since the detected intensity is linear in the Mueller elements, feeding in unit matrices one at a time builds it column by column. The condition numbers in this post are computed for this Fourier-coefficient matrix. A sample-domain matrix with equal noise on each sample gives different values, but the comparisons below come out the same.

The answer is unambiguous. 5:3 and 5:1 have identical condition numbers across the whole range of retardance. That is why the two curves in the right panel of figure 5 merge into one. 5:2, 5:4 and 7:3 also give 9.45 at 90°; only 4:3, with its one colliding address from section 4, is worse at 12.8.

Among ratios whose harmonic addresses do not collide, the ratio plays no part in noise amplification. What it changes is the highest harmonic, that is, how many times per period the detector must be read.

Given that 5:3 is a compromise forced by the detector and the motors, it is natural to assume that abandoning the mathematically optimal 5:1 costs something somewhere. The only cost turns out to be reading 33 times instead of 25. Not a digit of precision is lost. As long as addresses do not collide, the ratio can be treated purely as a timing question.

This also fits what Smith obtained in 2002 by minimizing the condition number. With exactly sixteen measurements the relationship between condition number and angular increments is intricate, but once the measurement is overdetermined a broad range of increments performs essentially identically. Overdetermination relieves the burden of choosing a ratio.

So what does set the precision?

<img src="/assets/img/posts/ellipsometer-dual-rotating-compensator-mme/en/fig5-condition-number.png" alt="Condition number against retardance" width="780">
_Figure 6. Condition number of the data reduction matrix (5:3) as both retardances are swept together. It is 9.45 at 90° and reaches its minimum of 2.49 at 128°, so a quarter-wave plate is not the best choice._

The retardance does. Figure 6 sweeps it, and the minimum sits at 128°, not at 90°. A quarter-wave plate gives a condition number of 9.45; a 128° plate gives 2.49. For the same detector noise, the error summed over the sixteen elements drops by about a factor of two. The condition-number ratio bounds the worst case; it is not the ratio of typical errors. Smith likewise said the noise could be cut by more than half. Smith's optimum of 127° is effectively the same answer.

128° does not help every element, though. Computing per-element errors for 36 reads per period, moving from 90° to 128° shrinks the error of fifteen elements by 1.3 to 2.7 times, while $M_{44}$ alone grows 1.6 times. $M_{44}$ rides only on $\sin\delta_1\sin\delta_2$ (figure 3), and at 128° that product is 0.62 of its value at 90°. If one circular-polarization element is the goal, 90° is better. The optimum the condition number points to is the optimum for measuring all sixteen evenly.

Why not 90° can be read off the structure in section 4. The fourth row and column ride on $\sin\delta$, which peaks at 90°. The modulated share of the second and third columns rides on $s = \sin^2(\delta/2)$, which keeps growing to 180°. The compromise between the two is pushed above 90°. Maximizing $\sin\delta\cdot\sin^2(\delta/2)$ alone gives 120°; accounting for the condition number of the full matrix lands near 128°.

Section 3 of post 2 said the retardance need not be exactly 90°. Here that sentence acquires a number. It is not merely unnecessary: 90° is a loss.

## 7. What Does Measuring Sixteen at Once Cost?

Time.

The optical period at 5:3 is 200 ms, so that is how long one Mueller matrix takes. The minimum measurement time for RAE and RCE in the table of post 1 was about 10 ms, making this twenty times slower.

The chain is easy to follow. Two compensators lift the highest harmonic from four times the RCE's compensator frequency to $24C$ and beyond. The motor-driven choice of 5:3 takes it to $32C$, and a higher harmonic demands at least 33 integrations per period. Each detector read costs over 5 ms, so the period stretches to 200 ms. Touch any link and the rest move with it. The author's dissertation (section 1.2) also counts this long measurement time among the drawbacks of the dual-rotating-compensator approach.

What 200 ms means depends on the application. It is fast enough to follow film growth in real time, which was the original motivation for the instrument. But against the 20 μs with which the photoelastic modulator of post 1 tracked liquid crystal molecules, it is four orders of magnitude away. The price of all sixteen is paid in time.

The exchange from post 1 changes shape once more here. Null ellipsometry spent time to buy accuracy; photometric methods trusted the detector to save time. Mueller matrix measurement spends time again, this time to buy information. What is sought changes; what is surrendered sits in the same place.

## 8. How Can the Slow Measurement Be Made Faster?

Two approaches diverge.

One pushes on the hardware. Collins's group described the development of this instrument in detail in 2001. With the fundamental mechanical frequency lowered to 2 Hz it took at least 0.25 s per Mueller matrix, a little slower than the 0.2 s of the design. The design table of section 5 already points that way: the ratio, the highest harmonic, the integration time and the motor frequency are linked in one chain, so a faster detector lets the whole chain be tightened. The ceiling of this approach is set by the detector.

The other redesigns when to read. The approach published in 2019 starts from the observation that there is no reason for sampling points in a continuous measurement to be evenly spaced. It builds an objective that tolerates both systematic error and detection noise and searches for optimal unevenly spaced points with a multi-objective genetic algorithm. The result reaches a globally optimal measurement matrix with the minimum number of samples (sixteen), while improving error immunity rather than degrading it.

The contrast is sharp. The first wants to read faster; the second wants to read less.

The second is possible because twenty-four coefficients for fifteen elements left room to cut samples from the start. Cut them close to sixteen, though, and the conclusion of section 6 no longer holds. With exactly sixteen measurements Smith found an intricate dependence on the angular increments, and Meng et al. report that the optimum is highly sensitive to the speed ratio. That is why where the samples sit had to be solved as an optimization problem. Even spacing was never a requirement of the physics; it was a convenience imposed to keep the Fourier analysis simple.

This also connects to the later parts of the series. Raising the rotation speed eventually runs into the motor and the detector. What remains is to change how the reading is done, or to remove the rotation altogether. Azzam's 1978 mention of Faraday cells is one end of that thought; the instruments of posts 5 and 6 are the other.

## 9. What Are the Surplus Measurements For?

Overdetermination has come up several times. Seeing what the surplus information actually does reveals the character of the instrument.

Thirty-six quantities are read per optical period. The signal contains fewer harmonics than that, so eleven of them ought to be zero from the start. If they are not, something in the instrument has shifted. The diagnostic Azzam performed with identities becomes, here, a matter of counting the coefficients that should have vanished.

The real use is calibration. Four things need calibrating: the two retardances, the two compensator phase offsets, and the angular offsets of the polarizer and analyzer. The order in which they are obtained is fixed.

Retardance should not depend on how the sample was mounted, so it is determined separately in a dedicated high-accuracy measurement. The phase offsets and angular offsets, by contrast, must be obtained in the actual configuration with the sample in place, because they change when the alignment does.

The method falls out of the coefficient structure. Squaring and summing the cosine and sine coefficients of a harmonic cancels the phase angle, which yields the polarizer and analyzer offsets. Taking their ratio instead leaves the phase angle intact, which gives the compensator phase offsets.

Those phase angles are not spread evenly over the two compensators. Only the first compensator's phase appears at $10C$ and $20C$, only the second's at $6C$ and $12C$. In the table of section 4 these are exactly the harmonics where only one compensator is involved. A window onto each compensator separately is already open in the signal. At the other harmonics the two enter mixed as sums and differences.

Several routes to the same quantity therefore exist. For one combination of polarizer offset and first-compensator phase alone there are five different paths through different harmonics. Agreement among the five says the instrument is sound; divergence narrows down what has shifted.

Overdetermination is not waste. Twenty-four coefficients solving for fifteen elements leave nine over, and those nine watch the instrument. They were not kept to gain information but to make the instrument trustworthy.

## Summary

Three things came out of this post.

First, two compensators are needed because of the structure of the weight $\mu_{ij} = a_i p_j$. Driving one side moves the weight in one direction only, which never sweeps all sixteen positions. The two must turn at different speeds because of harmonic addresses: the sums and differences of column and row frequencies must not coincide, or elements mix. When a sample is anisotropic and depolarizing the Jones formalism breaks down, and from that point all sixteen elements are required.

Second, the rotation ratio sets timing, not precision. 5:1 gives the lowest highest-order harmonic, but the minimum detector scan time and the low-frequency instability of the motors pushed the design to 5:3. Its addresses do not collide, so its condition number is identical to that of 5:1. What sets the precision is the retardance, and its optimum is near 128° rather than a quarter wave, although for $M_{44}$ alone 90° is better.

Third, the instrument carries the material for checking itself inside its own signal. Azzam's five identities and the five calibration routes of Collins and Koh come from the same property. As section 9 showed, it is the surplus measurements that do that work.

The next post follows from there. Going beyond using surplus measurements for instrument checks leads to assuming the components are not ideal and extracting their non-ideality from the measurement as well. Post 4 takes up the calibration of rotating-element ellipsometers, the substance of what post 2 called using cheap, non-ideal parts and making up the difference in calibration.

## References

- R. M. A. Azzam, "Photopolarimetric measurement of the Mueller matrix by Fourier analysis of a single detected signal," *Optics Letters* **2**, 148 (1978) — the origin of the dual rotating retarder.
- P. S. Hauge, "Recent developments in instrumentation in ellipsometry," *Surface Science* **96**, 108 (1980) — what each rotating configuration reaches.
- R. C. Thompson, J. R. Bottiger, E. S. Fry, "Measurement of polarized light interactions via the Mueller matrix," *Applied Optics* **19**, 1323 (1980) — simultaneous measurement with four electro-optic modulators.
- R. W. Collins and J. Koh, "Dual rotating-compensator multichannel ellipsometer: instrument design for real-time Mueller matrix spectroscopy of surfaces and films," *J. Opt. Soc. Am. A* **16**, 1997 (1999) — the 5:3 design and the coefficient table.
- J. Lee, J. Koh, R. W. Collins, "Multichannel Mueller matrix ellipsometer for real-time spectroscopy of anisotropic surfaces and films," *Optics Letters* **25**, 1573 (2000).
- J. Lee, J. Koh, R. W. Collins, "Dual rotating-compensator multichannel ellipsometer: instrument development for high-speed Mueller matrix spectroscopy of surfaces and thin films," *Rev. Sci. Instrum.* **72**, 1742 (2001).
- M. H. Smith, "Optimization of a dual-rotating-retarder Mueller matrix polarimeter," *Applied Optics* **41**, 2488 (2002) — condition-number minimization and the optimal retardance.
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, section 4.2.7 (configurations for Mueller matrix ellipsometry).
- O. Arteaga, J. Freudenthal, B. Wang, B. Kahr, "Mueller matrix polarimetry with four photoelastic modulators: theory and calibration," *Applied Optics* **51**, 6805 (2012).
- K. Meng, B. Jiang, C. D. Samolis, M. Alrished, K. Youcef-Toumi, "Unevenly spaced continuous measurement approach for dual rotating-retarder Mueller matrix ellipsometry," *Optics Express* **27**, 14736 (2019).
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD dissertation, Seoul National University, 2025, section 1.2.

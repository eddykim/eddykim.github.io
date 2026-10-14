---
title: "Computational Imaging and Ptychography 5 — Does Adding Acceleration Rescue a Stuck Reconstruction?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-momentum-difference-map/
page_id: ptycho-momentum-difference-map
date: 2026-12-01 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, mpie, difference-map, momentum, phase-retrieval]
description: "Revisiting the diffuser probe that rPIE could not solve with momentum (mPIE) and the difference map (DM), and sorting out what actually worked."
math: true
---

[Post 4](/en/posts/ptycho-probe-ambiguity-raster-rpie/) left one problem open. With a diffuser probe, a pinhole carrying a random phase plate, neither ePIE nor rPIE moved from an error of 0.99. Maiden, Johnson and Li, who proposed rPIE, offered a remedy in the same paper: mPIE, which adds momentum to rPIE. In their simulations mPIE halved rPIE's convergence time and was less sensitive to the pattern order.

This post checks that remedy directly, with three questions. Does mPIE really solve the diffuser? If so, what does the work? And how does the method that applies the difference map (DM) of Post 2 to ptychography handle the same problem? To give the answer up front: what solved the diffuser was not momentum but the small step that came with it. Momentum sped things up for some probes and made the reconstruction diverge for others. All calculations use the same specimen and 144-position Fermat spiral as Post 4, repeated with three to five pattern orders.

## 1. What Momentum Adds

A classic place for gradient descent to slow down is a long, narrow valley, steep in one direction and shallow in the other. A large step size bounces between the steep walls and diverges, so the step has to be sized for the steep direction. Then it barely moves along the shallow one. This is the limit of gradient descent treated in [Optimization Methods 1](/en/posts/optimization-gradient-descent/), and the same reason the Landweber iteration of [Post 1](/en/posts/ptycho-computational-imaging-inverse-problem/) was slow on ill-conditioned problems.

Momentum is usually pictured as a ball rolling downhill. The ball does not look only at the slope under it; it keeps the velocity it has built up. Motion back and forth between the walls keeps changing direction and cancels out, while motion along the floor in one direction accumulates. Maiden and colleagues adapted this to the PIE family. rPIE runs as usual, updating at one position at a time, and every $T$ updates the following two steps are added.

$$ v \leftarrow \eta\, v + \left(O - O_{\text{prev}}\right), \qquad O \leftarrow O + \eta\, v $$

Here $O_{\text{prev}}$ is the specimen $T$ updates ago, that is, right after the previous momentum step, and $v$ is the velocity. The value $\eta$, between 0 and 1, sets how much velocity to keep in the first equation (friction) and how much of it to add in the second (feedback). In short, it pushes the specimen once more in the direction the last $T$ updates moved it. With $\eta = 0$ it is rPIE. The probe has its own velocity, updated the same way. Maiden and colleagues used $T = 30$ and $\eta = 0.9$ and, in Table 3 of their paper, suggested a specimen step $\gamma_{\text{obj}}$ of 0.1–0.5. Since momentum enlarges the effective step, this can be read as a reduction. They used 0.2.

<img src="/assets/img/posts/ptycho-momentum-difference-map/en/fig1-momentum-valley.png" alt="Trajectories and distance to the minimum for gradient descent, moderate momentum, and excessive momentum on a narrow quadratic valley" width="900">
_Fig 1. On a narrow valley with the same step, trajectories (left) and distance to the minimum (right) for gradient descent, moderate momentum ($\eta$ = 0.6) and excessive momentum ($\eta$ = 0.9). Gradient descent bounces between the walls, moderate momentum builds speed along the floor, and too much overshoots the minimum._

Fig 1 uses $f = (x^2 + 25y^2)/2$ with the same step and only $\eta$ changed. Reaching a distance of 0.01 from the minimum takes 92 steps for gradient descent and 28 for $\eta = 0.6$. With $\eta = 0.9$ it takes 98, longer than without momentum. The accumulated velocity is so large that the trajectory overshoots the minimum, climbs the far wall, and comes back, again and again. So momentum is not better the more of it there is; each problem has a suitable amount. Fig 1 shows momentum in general, as a toy; it is not the mPIE update transcribed. Keeping this figure in mind makes Section 3 easier to read.

## 2. What Solved the Diffuser Probe?

The first step was to follow the prescription as given. With Maiden's values ($\alpha = 0.05$, $\beta = 1$, $T = 30$, $\eta = 0.9$, $\gamma_{\text{obj}} = 0.2$), mPIE reaches an error of 0.018–0.023 on the diffuser in all five runs. Against rPIE's 0.99, the remedy works.

But mPIE changes two things relative to rPIE. It adds momentum, and it reduces the specimen step from 1 to 0.2. Which one did the work? Removing momentum and changing only the steps answers this. The specimen step $\gamma_{\text{obj}}$, the probe step $\gamma_{\text{prb}}$, and the form of the probe update $\beta$ (1 for the ePIE form, 0.1 for the rPIE form) were each set to two values, and the eight combinations were run with three pattern orders.

| Specimen step $\gamma_{\text{obj}}$ | Probe step $\gamma_{\text{prb}}$ | Probe update $\beta$ | Specimen error after 300 it. |
|---|---|---|---|
| 1 | 1 | 1 | 0.997–1.000 |
| 1 | 1 | 0.1 | 0.997–1.000 |
| 1 | 0.25 | 1 | 0.994–0.999 |
| 1 | 0.25 | 0.1 | 0.999–1.000 |
| 0.25 | 1 | 1 | 0.020–0.026 |
| 0.25 | 1 | 0.1 | 0.036–0.082 |
| 0.25 | 0.25 | 1 | 0.019–0.021 |
| 0.25 | 0.25 | 0.1 | 0.019–0.067 |

The split fits in one line. With a specimen step of 1, all twelve runs fail whatever else is changed; with 0.25, all twelve go below 0.1. No momentum is involved. What solved the diffuser was the specimen step, not momentum.

Why does the step need to be small? Narrowing it down further, $\gamma_{\text{obj}} = 0.25$ was used only for the first $K$ iterations and then set back to 1. $K = 3$ was enough: all three runs reached 0.018–0.023, and $K = 10$ and $K = 30$ gave the same. The striking part is that after the third iteration the error was still 0.79–0.88. The reconstruction had barely started, yet the outcome was already decided.

In short, the outcome is decided in the first three iterations. One way to read this is as follows. During the first few iterations the probe is still close to its initial guess, a plain disc, and badly wrong. With a specimen step of 1, rPIE divides the corrected exit wave almost directly by this wrong probe. The specimen sets at once into a shape that absorbs the probe's error, and the probe then adapts to that specimen and stays wrong. The two hold each other in place and stop. With a small step the specimen moves slowly, and the probe finds its shape first. Three iterations are enough to establish that order.

<img src="/assets/img/posts/ptycho-momentum-difference-map/en/fig2-diffuser-step.png" alt="Specimen error and data residual for rPIE, rPIE with a small step for the first 10 iterations, and mPIE on the diffuser probe" width="900">
_Fig 2. Specimen error (left) and data residual (right) on the diffuser probe for rPIE, rPIE with the specimen step reduced to 0.25 for the first 10 iterations, and mPIE (lines: median of 3 orders, bands: range). Shrinking the step for just the first few iterations lets rPIE succeed, and momentum is only a little faster early._

Did momentum then do nothing? Fig 2 compares rPIE with a small step for the first 10 iterations (blue) against mPIE (green). At 30 iterations the errors are 0.056 and 0.036, so mPIE is somewhat faster. At 100 they are nearly equal, 0.025 and 0.023, and at 300 they are 0.022 and 0.021. The data residual on the right is, if anything, slightly lower for rPIE. On this problem, what momentum contributed was a modest early speed-up.

## 3. Why Does Momentum Wreck the Pinhole Reconstruction?

Running the same mPIE on the pinhole probe used since Post 3 tells a different story. With Maiden's values, mPIE reaches 0.024 in only one of five runs; the others fail. rPIE without momentum gives 0.00055–0.00072 on the same problem in all five. What goes wrong?

The cause is the scale exchange that Post 4 called harmless. Multiplying the specimen by $c$ and dividing the probe by $c$ leaves the measurement unchanged. Along that direction the cost does not change at all, so if the updates drift slightly that way, nothing pulls them back. In terms of the valley in Section 1, it is a perfectly flat floor with no wall, and a ball on it has no reason to stop. Gradient descent leaves the small steps accumulated along this direction as they are, but momentum turns that small drift into velocity and keeps it rolling. The left of Fig 3 shows the total probe power for mPIE run without holding the power fixed. After 300 iterations, three of five runs had grown by more than $10^{21}$ times the target, and one had shrunk to $10^{-84}$ times. The specimen moves the opposite way by the same factor.

It looks as if fixing the probe power should stop the drift. By Parseval's theorem, the total intensity of a pattern recorded where the transmission is 1 equals $\sum \lvert P \rvert^2$. So at every momentum step the probe power was set to the total intensity of the brightest pattern, and the specimen was multiplied by the reciprocal factor. The product $PO$ is unchanged, so the measurement is unaffected. Maiden and colleagues also mention a power-correction step in their supplement (the supplement was not consulted for this post). All rPIE and mPIE runs in this post use this control.

<img src="/assets/img/posts/ptycho-momentum-difference-map/en/fig3-drift-divergence.png" alt="Drift of the probe power without power control, and the amplitude of an unilluminated pixel with power control, for a pinhole probe" width="900">
_Fig 3. On the pinhole probe, without power control the probe power drifts by many orders of magnitude (left); with power control the same drift piles up in the specimen amplitude at unilluminated pixels (right). On the right, green is momentum on both, red on the specimen only, blue on the probe only, and the grey dashed line is rPIE. Only runs with specimen momentum blow up._

Yet with the power fixed, mPIE on the pinhole did not improve. In three of five runs the data residual blew up to $5\times10^4$–$6\times10^{11}$. What is growing becomes clear from a pixel that is never illuminated. rPIE never touches such a pixel, so its value is the accumulated product of the factors that power control has applied to the specimen. As the right of Fig 3 shows, it stays near 1 (1.04–1.37) for rPIE. For mPIE it grew to $10^{16}$ or shrank to $10^{-43}$. In other words, even with the power fixed, momentum keeps pushing along the scale direction at every step, and each time power control pulls the probe back, the factor accumulates on the specimen side. Bright pixels are pulled back into place by the data at the next update, but pixels that receive less than 1 % of the illumination have little to pull them back and simply carry the factor. That is why, in all three diverging runs, the largest specimen amplitude grew at pixels receiving less than 0.2 % of the peak illumination, and in two of them the maximum ended up at pixels that receive no illumination at all.

Splitting where the momentum is applied narrows the cause down. With momentum on the specimen only, four of five runs fail and the same accumulation appears. With momentum on the probe only, all five converge (0.004 four times, 0.022 once) and the accumulated factor stays at 1.2–1.4. On the pinhole, the culprit is the specimen momentum.

Then would a smaller $\eta$ help? Fig 4 shows five values each of $\eta$ and $\gamma_{\text{obj}}$ on all three probes.

<img src="/assets/img/posts/ptycho-momentum-difference-map/en/fig4-parameter-map.png" alt="Maps of the final error over momentum and specimen step for the pinhole, defocused focused, and diffuser probes" width="950">
_Fig 4. Final error (300 iterations, median of 3 orders) on three probes as momentum $\eta$ and specimen step $$\gamma_{obj}$$ vary. The number in each cell counts runs below 0.05 and the red cell marks the values used by Maiden and colleagues. The pinhole fails at large momentum, and the diffuser fails at large steps._

The three maps look different. For the pinhole, seven of the ten cells with $\eta$ of 0.8 or more never succeed. For the defocused probe from Post 4, every cell has a median below 0.003, and the error falls as momentum grows, below $10^{-4}$ at $\eta = 0.9$ (though two cells at $\eta = 0.9$ still fail in one run of three). For the diffuser, specimen steps above 0.35 fail regardless of momentum, and even within that range the successful cells are scattered. The $\eta = 0.7$ row fails in all five cells, while the $\eta = 0.8$ and 0.9 rows have successful cells. Maiden's values (red cell) succeed on the defocused and diffuser probes but fail on the pinhole.

The open-source toolbox PtyLab uses friction 0.7, feedback 0.3, and a step of 0.25 as its mPIE defaults. PtyLab, however, applies momentum per position with 5 % probability, whereas here it was applied every 20 iterations. These do the opposite: they work on the pinhole (0.0026) and the defocused probe (0.0006) and fail on the diffuser in all five runs. Only three cells, in the lower-left corner, succeed in all three runs on all three probes: no momentum with steps of 0.1 and 0.2, and momentum 0.5 with a step of 0.1. Momentum is a knob that has to be retuned for each probe. Maiden and colleagues themselves name the large set of parameters as mPIE's weakness in their conclusion.

## 4. How Is the Difference Map (DM) Different?

Every method so far has gone through the positions one at a time, correcting the specimen and probe as it goes. The DM that Thibault and colleagues applied to X-ray ptychography in 2008 starts differently. Its state is not the specimen but the whole set of exit waves at the 144 positions, $\Psi = (\psi_1, \dots, \psi_{144})$. The problem then becomes the intersection of two sets, as in [Post 2](/en/posts/ptycho-phase-problem-iterative-retrieval/). One set contains the bundles whose Fourier amplitudes match the measurement at every position (the Fourier set). The other contains the bundles that can be written as $\psi_j = P(r - R_j)O(r)$ for some specimen $O$ and probe $P$ (the overlap set).

The projection onto the Fourier set, $P_F$, is the same as in Post 3: at each position only the amplitude is replaced by the measurement. The projection onto the overlap set, $P_O$, is a least-squares problem, finding the $P \cdot O$ closest to the bundle of exit waves. With the probe fixed, the answer follows pixel by pixel.

$$ O(r) = \frac{\sum_j P^*(r - R_j)\,\psi_j(r)}{\sum_j \lvert P(r - R_j) \rvert^2} $$

Here the numerator takes every exit wave that illuminated the pixel, undoes the probe, and sums; the denominator is the total illumination the pixel received. In other words, it is an average of what the overlapping positions say about that pixel, weighted by illumination. The probe is found from the same formula with the roles exchanged, and the two are alternated a few times. One DM step is written with these two projections.

$$ \Psi \leftarrow \Psi + P_F\!\left(2P_O\Psi - \Psi\right) - P_O\Psi $$

In terms of reflections $R = 2P - I$, this equals $\tfrac{1}{2}(R_F R_O + I)\Psi$; computing one step both ways agrees to $10^{-16}$. It has the same form as RAAR with $\beta = 1$, written $\tfrac{1}{2}(R_s R_m + I)$ in Post 2, except that the order of the two sets is swapped, so the measurement (Fourier) reflection comes last. DM is Post 2's HIO-family algorithm applied to the two sets of ptychography.

A few things follow from this structure. There is no step size, no momentum (with $\beta$ fixed at 1), and no position order, so the same input always gives the same answer. The price is holding the whole bundle of exit waves. In these calculations the specimen is 0.5 MB while the exit-wave bundle is 9.4 MB. With 10,000 positions and 256×256-pixel patterns, the bundle alone is 5.2 GB in complex64 and 10.5 GB in complex128. In real measurements that becomes a large share of memory. Time per iteration was similar to rPIE.

<img src="/assets/img/posts/ptycho-momentum-difference-map/en/fig5-three-probes.png" alt="Convergence of rPIE, rPIE with a small-step start, two mPIE settings, and DM for the pinhole, defocused focused, and diffuser probes" width="950">
_Fig 5. Convergence of five methods on three probes (lines: median of 5 orders, bands: range; DM has no order). On the pinhole DM is lowest at 300 iterations (4×10⁻⁵); on the diffuser rPIE and mPIE (PtyLab settings) fail and DM oscillates._

In Fig 5, DM (grey) is slow at first. On the pinhole its error at 30 iterations is 0.075, larger than rPIE's 0.009. Once it settles, however, it falls quickly and reaches $4 \times 10^{-5}$ at 300 iterations, the lowest of the five methods. On the defocused probe it reaches 0.00026, the same as rPIE. On the diffuser it drops to 0.04 around 50 iterations, then rises again and oscillates near 0.12 at 300. Post 2 noted that HIO-family trajectories are sensitive to small perturbations, and that property seems to show up on a problem with a complicated probe like the diffuser. These calculations alone, however, do not establish the cause.

## 5. Why Does the Diffuser Reconstruction Stall at an Error of 0.02?

Even the methods that solve the diffuser stop near an error of 0.02. rPIE with a small step for the first 10 iterations is still at 0.018–0.023 after 500 iterations, with a data residual that stays at 3–4 %. With the probe given, the same data reach 0.0014, so the data are not the limit. Where does the remaining error come from?

<img src="/assets/img/posts/ptycho-momentum-difference-map/en/fig6-diffuser-stall.png" alt="True and recovered diffuser probes, and specimen error images with the probe recovered and with the probe known" width="950">
_Fig 6. The stalled diffuser reconstruction. The recovered probe has the right shape but is shifted by 2–4 pixels (red + : centre of mass of the true probe, cyan + : of the recovered probe, red line: half-maximum contour of the true probe). The two error images share one colour range (0–0.1), and a known probe gives an error more than ten times smaller._

Looking at the recovered probe, its shape is almost the same as the true probe, but it sits 2–4 pixels away. This is the same picture as the stagnation on the jittered grid in Post 4. Shifting the specimen and probe together is an ambiguity that leaves the measurement unchanged, but the edges of the specimen and the region outside the illumination are tied to their initial values and cannot move with the rest. That mismatch appears to remain as residual. Moving the probe's centre of mass back to its starting position at every iteration, and shifting the specimen with it, reduced the error to 0.012, but a 3 % residual still remained. Because the true probe's centre of mass differs from the centre of the initial guess, this alone cannot close the gap. This post did not find a way to remove the stagnation completely.

## Summary

| Probe | rPIE | rPIE, small step for first 10 | mPIE (Maiden values) | mPIE (PtyLab settings) | DM |
|---|---|---|---|---|---|
| Pinhole | 0.00065 | 0.00065 | 4 of 5 fail | 0.0026 | 0.00004 |
| Defocused focused | 0.00025 | 0.00024 | 0.00015 | 0.00061 | 0.00026 |
| Diffuser | fails (1.0) | 0.022 | 0.021 | fails (1.0) | 0.12 (oscillating) |

The values are specimen errors after 300 iterations, medians over five orders. What solved the diffuser probe was a small specimen step, not momentum, and it only needed to be small for the first three iterations. The key is to keep the specimen from setting too quickly while the probe is still wrong. Momentum lowers the error further and speeds up the start on problems that already go well, such as the defocused probe, but it drifts along the scale exchange. Even with the probe power held fixed, the drift from specimen momentum can pile up in weakly illuminated pixels and ruin the reconstruction, and the suitable values differ from probe to probe. DM has no step size or momentum and goes deepest on the pinhole, but on the diffuser it oscillated and did not settle.

In practice this suggests an order. Start with rPIE and a reduced specimen step for the first few iterations. Watch the data residual alongside to confirm convergence. If a faster or deeper answer is needed, add momentum or DM at that point, with the probe power fixed. Whatever the method, trust a result only after running it with several orders.

## Next

All the data so far have been ideal: no noise, exact scan positions, and fully coherent light. Post 6 breaks the first of these assumptions, noise. It covers a cost function suited to the Poisson noise of photon counting (maximum likelihood), and the noise-model discussion deferred from Post 1 returns there. Correcting the scan positions along with the specimen when they are wrong, and representing partial coherence with several probe modes, follow in Post 7.

## References

- A. Maiden, D. Johnson, P. Li, "[Further improvements to the ptychographical iterative engine](https://doi.org/10.1364/OPTICA.4.000736)," *Optica* 4(7), 736-745 (2017). The momentum equations of mPIE (Eqs. 19–22), the recommended parameter ranges (Table 3), the diffuse-probe simulations, and the limitation of many parameters.
- P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, "[High-resolution scanning X-ray diffraction microscopy](https://doi.org/10.1126/science.1158573)," *Science* 321, 379-382 (2008). DM with the exit-wave bundle as its state, the Fourier and overlap constraints, and an overlap projection that yields both specimen and probe.
- V. Elser, "[Phase retrieval by iterated projections](https://doi.org/10.1364/JOSAA.20.000040)," *Journal of the Optical Society of America A* 20(1), 40 (2003). The original difference map.
- L. Loetgering et al., "[PtyLab.m/py/jl: a cross-platform, open-source inverse modeling toolbox for conventional and Fourier ptychography](https://doi.org/10.1364/OE.485370)," *Optics Express* 31(9), 13763 (2023). The mPIE defaults (friction 0.7, feedback 0.3, step 0.25) were checked in the implementation in the [PtyLab.py repository](https://github.com/PtyLab/PtyLab.py).
- S. Ruder, "[An overview of gradient descent optimization algorithms](https://arxiv.org/abs/1609.04747)," arXiv:1609.04747 (2016). The overview Maiden and colleagues cite for momentum.

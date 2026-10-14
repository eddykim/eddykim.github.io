---
title: "Computational Imaging and Ptychography 7 — Can the Image Still Be Recovered When the Scan Positions Are Off and the Light Is Unsteady?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-position-correction-mixed-states/
page_id: ptycho-position-correction-mixed-states
date: 2026-12-07 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [ptychography, position-correction, partial-coherence, mixed-states, phase-retrieval]
description: "How a reconstruction breaks when scan positions are slightly wrong or the light is not fully coherent, and how the positions and coherent modes can be recovered from the data, by calculation."
math: true
---

[Post 6](/en/posts/ptycho-noise-maximum-likelihood/) dealt with measurements that carry noise. Real experiments break two more assumptions. One is the scan position. So far the post assumed we know exactly where the probe illuminated, but a real stage misses its commanded position slightly and drifts slowly with thermal expansion or vibration. The other is the coherence of the light, meaning how well its waves keep a fixed phase relation. So far the probe was a single wave, but a real source has a finite size, so several mutually incoherent waves arrive together.

The two problems are of the same kind: the model that explains the measurement is wrong. So the questions are the same too. How wrong does the model have to be before the reconstruction breaks? And can the wrong part be found from the data and corrected? The calculations use the pinhole probe and the 144-position Fermat spiral used since Post 3, without noise. The size of a position error is quoted as a root-mean-square (rms) value. Position errors were repeated with three error draws and coherence with three pattern orders.

## 1. Does a Slightly Wrong Position Ruin the Reconstruction?

With a scan step of 8 pixels, Gaussian errors were added to each position, and the reconstruction trusted those wrong positions as they were. At what error does trouble start?

<img src="/assets/img/posts/ptycho-position-correction-mixed-states/en/fig1-position-error-recon.png" alt="True specimen phase and phases reconstructed with positions wrong by 1.4 px rms, without correction, with annealing, and with the gradient method" width="950">
_Fig 1. Phase reconstructed from the same data with positions wrong by 1.4 px rms. Uncorrected it collapses; correcting the positions recovers it._

It breaks at far smaller errors than one might expect. With exact positions the specimen error is 0.0023. With positions wrong by only 0.35 px rms, one twenty-third of the scan step, it becomes 0.060, 27 times larger. At 1.4 px it is 0.23, and the image becomes blotchy, as in the second panel of Fig 1.

Why so sensitive? Recall how overlap provides information. A specimen pixel is illuminated by different parts of the probe from several positions, and the answer is fixed when those patterns point to the same specimen. With wrong positions, the algorithm believes the pixel was illuminated by the wrong part of the probe. Where the intensity changes strongly within one pixel, as at the edge and rings of a pinhole probe, a shift of 0.35 px changes the exit wave substantially. Each pattern then demands a different specimen, and the algorithm settles on a blotchy compromise between them.

## 2. Can Wrong Positions Be Found and Corrected from the Data?

Then the positions can be treated as unknowns and found together with the rest. In their handbook chapter, Rodenburg and Maiden describe two widely used methods. Both start from a reconstruction that ignores position errors and then move the positions slightly one at a time. This post compares annealing with a gradient method that follows directly from the cost function of Post 6.

The first is the annealing method of Maiden and colleagues. For each position, a few random candidates are generated around the current position, the diffraction pattern of each candidate is computed from the current specimen and probe, and the position moves to the candidate that best matches the measurement. The radius for drawing candidates shrinks as the iterations proceed: at first positions can wander far, and later they settle on a single point. Here the radius started at 1 pixel, comparing four candidates with the current position.

The second is a gradient method. Differentiating the amplitude cost of Post 6 with respect to the position tells which way to move the position to lower the cost.

$$ \frac{\partial E_j}{\partial \mathbf{d}} = 2\,\mathrm{Re}\left\langle \chi_j,\ \mathcal{F}\left[ -(\nabla P)(r - \mathbf{d})\, O_j \right] \right\rangle $$

Here $\mathbf{d}$ is the sub-pixel shift of position $j$, and $\chi_j = \left(1 - \sqrt{n_j}/\lvert\Psi_j\rvert\right)\Psi_j$ is the amplitude-cost gradient from Post 6. The term in brackets on the right is the direction in which the exit wave changes when the probe moves. In other words, it measures how much the change of the exit wave overlaps with the correction the data demand. The larger the probe gradient $\nabla P$, the larger this quantity, so the same property that made the reconstruction sensitive to position in Section 1 becomes the clue for finding it. Positions were moved by a Gauss–Newton step, this gradient divided by its squared size, limited to 0.5 pixel per step.

One point deserves attention. Moving all positions together in the same direction is the shift ambiguity from [Post 4](/en/posts/ptycho-probe-ambiguity-raster-rpie/). If the specimen moves along with them, the measurement is unchanged, so the data cannot fix the mean of the positions. The mean of the corrections was therefore reset to zero at every iteration.

<img src="/assets/img/posts/ptycho-position-correction-mixed-states/en/fig2-position-correction.png" alt="Final specimen error versus position error, and the position error reduced by annealing and the gradient method at 1.4 px rms" width="950">
_Fig 2. Left, final specimen error versus position error (3 error draws); right, convergence of the remaining position error at 1.4 px rms. Both methods cut the error sharply._

| Position error (rms) | No correction | Annealing | Gradient |
|---|---|---|---|
| 0.35 px | 0.060 | 0.0024 | 0.0032 |
| 0.71 px | 0.14 | 0.0026 | 0.0038 |
| 1.4 px | 0.23 | 0.011 | 0.0065 |
| 2.8 px | 0.38 | 0.010 | 0.036 |

The values are median specimen errors after 300 iterations. Both methods reduced the position errors sharply. An error of 0.71 px fell to 0.005 px with annealing and 0.003 px with the gradient method, and the specimen error came close to the 0.0023 of exact positions. The positions were recovered from the data to a few hundredths of a pixel or better.

<img src="/assets/img/posts/ptycho-position-correction-mixed-states/en/fig3-position-arrows.png" alt="Errors of the nominal positions and the errors left after annealing, drawn as arrows" width="900">
_Fig 3. Errors of the 144 positions drawn as arrows (magnified 5×). After annealing they shrink almost to dots._

The two methods part ways when the errors are large. Annealing searched only a 1-pixel radius at a time yet corrected an error of 2.8 px (0.03 px left). Because new candidates are drawn around the position it has just moved to, a position can walk a long way over many iterations. The gradient method left positions about 0.5 px off in the same case. A gradient uses only information near the current position, so it slows down when started outside the range where the cost decreases smoothly. Conversely, for errors of 0.7 px or less, the gradient method located the positions more precisely (0.002–0.003 px, against 0.005–0.007 px for annealing). The handbook notes that annealing and the cross-correlation method increase the total computation by about a factor of three. If the positions are badly wrong, start with annealing; if they are nearly right, refine with the gradient.

## 3. How Does the Measurement Change When the Light Is Not Fully Coherent?

Now the second assumption. If the source is a point, all the light from it is coherent. If the source has a size, light from different points of the source is mutually incoherent. Light from each point illuminates the pinhole at a slightly different angle, and the diffraction patterns from these points add at the detector as intensities, without interfering. That the degree of coherence is set by the source size and distance was covered in the spatial-coherence section of [Wave Optics Background 2](/en/posts/wave-optics-coherence/).

How is the measurement written then? Incoherent waves add as intensities, so with a probe $P_s$ and an intensity share $w_s$ for each source point,

$$ I_j = \sum_s w_s \left\lvert \mathcal{F}\{ P_s(r)\, O(r + R_j) \} \right\rvert^2 . $$

There can be tens or hundreds of source points. Fortunately this can be rewritten with far fewer orthogonal modes, the coherent-mode decomposition proposed by Wolf. A singular value decomposition (SVD) of the source-point probes gives mutually orthogonal modes $P_k$ and the intensity share each carries, its occupancy. The formula keeps the same form with $s$ replaced by $k$, and most of the intensity is concentrated in the first few modes.

<img src="/assets/img/posts/ptycho-position-correction-mixed-states/en/fig4-coherent-modes.png" alt="Coherent-mode amplitudes and occupancies of the probe for a coherent source and two extended sources" width="850">
_Fig 4. Coherent modes (amplitude) and occupancies of the probe made by an extended source. The wider the source, the smaller the first mode's share._

Fig 4 uses two widths for the source illuminating the pinhole. The first mode's occupancy is 72 % for the narrower one and 47 % for the wider one. The second and third modes are two-lobed and have equal occupancies. Thibault and Menzel reported that even in an X-ray experiment thought to be highly coherent the first mode carried about 80 %, and in an experiment with reduced coherence less than 40 %. Mixtures of this size are common in practice.

## 4. How Is Light That Is Not Fully Coherent Reconstructed?

What happens when data from mixed modes are reconstructed with a single coherent wave? The top middle of Fig 6 shows the result. With 72 % occupancy the specimen error is 0.86–0.98: effectively nothing is recovered. The same holds at 47 %. In Thibault and Menzel's experiment, the reconstruction assuming coherence also failed.

The cure is to fix the model. Thibault and Menzel represented the probe with several modes and reconstructed them together. This is called mixed-state ptychography. Only the Fourier projection changes. Since the measurement is the sum of the mode intensities, every mode is multiplied by the same factor so that the sum matches the measurement. Here $n$ is the measured photon count and $\Psi_k$ is the exit wave of mode $k$.

$$ \Psi'_k = \Psi_k \, \frac{\sqrt{n}}{\sqrt{\sum_{l} \lvert \Psi_l \rvert^2}} $$

With a single mode this is the Fourier projection of Post 3. The handbook gives an intuitive picture. When two incoherent waves arrive at one pixel, what is lost is two phases and two amplitudes, four unknowns. But the condition that the two intensities sum to the measured value remains, so there are really three. The formula above is exactly that condition. The specimen is corrected by summing the corrections of all modes, and each mode is corrected separately. Mixing the modes among themselves does not change the summed intensity, so the modes were orthogonalized by SVD at every iteration to fix their order.

<img src="/assets/img/posts/ptycho-position-correction-mixed-states/en/fig5-modes-vs-k.png" alt="Specimen error and data residual of mixed-state rPIE and ML refinement versus number of modes, with the residual of the truncated true state" width="950">
_Fig 5. Specimen error and data residual versus number of modes (3 orders). The dashed line is the residual of the true state truncated to K modes; at 72 % occupancy the reconstruction's residual is of the same size as this line._

| Modes | 72 % occupancy | 47 % occupancy |
|---|---|---|
| 1 | 0.97 (fails) | 0.99 (fails) |
| 2 | 0.12 | 0.89 (fails) |
| 3 | 0.028 | 0.92 (fails) |
| 5 | 0.040 | 0.076 (fails in one of three) |
| 8 | 0.025 | 0.034 |

At 72 % occupancy the reconstruction works from three modes. At 47 % it needs at least five, and eight are stable. The stronger the mixing, the more modes are needed. The recovered occupancies match the true ones well. With eight modes, the 72 % case gives 71/12/12/1.5 % (true 72/12/12/1.5 %) and the 47 % case gives 47/18/18/5.9 % (true 47/18/18/5.1 %).

<img src="/assets/img/posts/ptycho-position-correction-mixed-states/en/fig6-mode-reconstruction.png" alt="Phase reconstructed with one and three modes from data with 72 % occupancy, and amplitudes of the true and recovered modes" width="950">
_Fig 6. Phase reconstructed with one and three modes from data with a 72 % dominant mode, and the true and recovered modes. The first mode agrees; the second and third differ in orientation._

Comparing the modes at the bottom of Fig 6, the first mode is nearly identical, but the second and third have their lobes pointing in different directions. Their occupancies are equal, so any combination mixing the two gives the same data. Modes with equal occupancy have no fixed orientation. This is not a failure of the reconstruction but a freedom the measurement cannot distinguish, of the same kind as the ambiguities of Post 4.

## 5. Why Does Error Remain Even with the Right Number of Modes?

Even with the right number of modes, the specimen error is 0.025–0.04, far above the 0.0006 for coherent data. Where does it come from? As in Posts 4 and 6, the data residual separates the possibilities.

The dashed lines on the right of Fig 5 are the residuals of the true state truncated to its first $K$ modes. At 72 % occupancy, truncating to three modes gives a residual of $3.7\times10^{-2}$, and to five modes $1.1\times10^{-2}$. The reconstruction's residual is lower at three modes, $2.3\times10^{-2}$, and slightly higher at five, $1.6\times10^{-2}$; both are of the same size. With three or five modes, the remaining residual is a limit of the model. The intensity carried by the discarded modes cannot be explained by the remaining ones, so the specimen absorbs part of the difference.

With eight modes, the true state could bring the residual down to $1.1\times10^{-3}$. The reconstruction's residual, however, stopped at about $1\times10^{-2}$. This time it is not the model but the solution that stagnated. Continuing with the maximum-likelihood (ML) refinement of Post 6, adapted to modes, helps a little (0.028 → 0.020 for 72 % with three modes, 0.034 → 0.028 for 47 % with eight). It did not remove the stagnation. Thibault and Menzel noted that the data must grow as the number of modes grows. Whether 144 positions are too few to resolve eight modes, or something else is at work, was not settled here.

## Summary

| Situation | Without correction | With correction |
|---|---|---|
| Position error 0.35 px rms | 0.060 | 0.0024 (annealing) |
| Position error 2.8 px rms | 0.38 | 0.010 (annealing) |
| First mode 72 % | 0.97 (one mode) | 0.025–0.040 (three or more modes) |
| First mode 47 % | 0.99 (one mode) | 0.034 (eight modes) |

Positions wrong by only a twentieth of the scan step badly damage the reconstruction. The same data, however, can recover the positions to a few hundredths of a pixel or better. Annealing is robust to large errors, and the gradient method corrects small errors more precisely. The data cannot set the mean of the positions, so it is held fixed.

Treating incoherent light as a single wave makes the reconstruction fail outright. With the probe represented by several coherent modes, the data reveal the modes and their occupancies. The stronger the mixing, the more modes are needed. Too few modes leave the model wrong by the discarded modes; too many make the problem harder to solve. In these calculations, the smaller the first mode's share, the larger the remaining error.

Both problems teach the same lesson. When the data do not fit the model, add the variables that explain the mismatch and solve for them with the data. The generous information from ptychographic overlap is what makes room for that.

## Next

So far the probe has been moved over the specimen. Post 8 turns this around. With the specimen under a lens, several low-resolution images are taken while the illumination angle changes, giving a resolution beyond the numerical aperture of the lens. This is Fourier ptychography. It overlaps angles instead of positions and Fourier space instead of real space, but the mathematics has the same form as before.

## References

- J. M. Rodenburg, A. M. Maiden, "[Ptychography](https://eprints.whiterose.ac.uk/id/eprint/127795/)," in *Springer Handbook of Microscopy*, P. W. Hawkes, J. C. H. Spence (eds.), Springer, 2019, pp. 819-904, [doi:10.1007/978-3-030-00069-1_17](https://doi.org/10.1007/978-3-030-00069-1_17). Position correction (annealing, cross-correlation, conjugate gradients) and its cost, mixed states and coherent-mode decomposition (Sec. 8).
- P. Thibault, A. Menzel, "[Reconstructing state mixtures from diffraction measurements](https://doi.org/10.1038/nature11806)," *Nature* 494, 68-71 (2013). Mixed-state ptychography, the data requirement that grows with the number of modes, and the dominant-mode occupancy in X-ray experiments.
- A. M. Maiden, M. J. Humphry, M. C. Sarahan, B. Kraus, J. M. Rodenburg, "[An annealing algorithm to correct positioning errors in ptychography](https://doi.org/10.1016/j.ultramic.2012.06.001)," *Ultramicroscopy* 120, 64-72 (2012). The original annealing position correction. The original was not consulted for this post; the handbook's description was followed.
- M. Guizar-Sicairos, J. R. Fienup, "[Phase retrieval with transverse translation diversity: a nonlinear optimization approach](https://doi.org/10.1364/OE.16.007264)," *Optics Express* 16(10), 7264 (2008). An early study finding positions by gradients. The original was not consulted; cited via the handbook.

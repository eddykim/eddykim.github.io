---
title: "Computational Imaging and Ptychography 2 — Can Phase Be Recovered from Intensity Alone? Phase Retrieval and Iterative Methods"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ptycho-phase-problem-iterative-retrieval/
page_id: ptycho-phase-problem-iterative-retrieval
date: 2026-11-19 20:00:00 +0900
categories: [Computation, Computational Imaging]
tags: [phase-retrieval, coherent-diffraction-imaging, gerchberg-saxton, hybrid-input-output, oversampling]
description: "Why recovering phase from diffraction intensity is possible at all, and where it gets stuck, from Gerchberg–Saxton to HIO."
math: true
---

[Post 1](/en/posts/ptycho-computational-imaging-inverse-problem/) dealt with measurements that are linear in the unknowns, $y = Ax$. Noise wrecked the naive inversion, but every component of the unknown was at least present in the measurement at some strength. This post is about what happens when that premise fails.

Light and electron waves carry both an amplitude and a phase. Yet every camera and detector records only intensity, the squared amplitude. Visible light oscillates more than $10^{14}$ times per second, and no detector can follow the phase of that oscillation. The measurement becomes $y = \lvert Ax \rvert^2$, and the phase disappears from it entirely.

How much is lost with the phase? Fourier transform two photos, swap their magnitudes and phases, and the answer becomes clear.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig1-phase-swap.png" alt="Result of swapping the Fourier magnitude and phase of two images" width="850">
_Fig 1. Swap Fourier magnitude and phase, and the image follows the phase_

The magnitude of A with the phase of B looks like B, and the reverse looks like A. Most of an image's outlines and layout live in the phase. A diffraction experiment, however, throws exactly that phase away and keeps only the magnitude. X-ray crystallography, coherent diffraction imaging (CDI), and ptychography, the destination of this series, all have to solve this phase problem. This post follows the classical methods for recovering phase from a single pattern: the iterative method of Gerchberg and Saxton (1972), how Fienup improved it, and why it can still get stuck.

## 1. What Intensity Leaves Behind

What exactly survives an intensity measurement? Inverse Fourier transforming the diffraction intensity $\lvert F(k) \rvert^2$ gives the autocorrelation of the object.

$$ \mathcal{F}^{-1}\{\lvert F \rvert^2\}(r) = \sum_{r'} o(r' + r)\, o^*(r') $$

It is the object overlaid on itself, multiplied and summed for every displacement $r$. Measuring intensity means knowing this autocorrelation, and nothing more. Objects that share an autocorrelation therefore cannot be told apart by the measurement.

Three such pairs always exist: the object multiplied by a constant phase (global phase), the object shifted as a whole (translation), and the object inverted through the origin and complex conjugated (conjugate inversion). The last is called the twin image. None of the three changes the autocorrelation. Shechtman and colleagues call them the trivial ambiguities. No algorithm can distinguish them, so they are counted as the same answer from the outset.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig2-trivial-ambiguities.png" alt="The original object, a translated copy, and its twin image, with their identical diffraction intensity" width="850">
_Fig 2. Three different objects produce exactly the same diffraction intensity_

This matters when results are evaluated. A reconstruction that comes out flipped or shifted is not wrong. Every error in this post is measured after aligning translation, twin image, and global phase. In the first computations for this post the alignment was left out, and runs that had found the twin image exactly were tallied as failures with an error of 0.2.

## 2. One Dimension Does Not Work

Once the trivial ambiguities are set aside, is the answer unique? Not in one dimension.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig3-1d-counterexample.png" alt="Two 1D signals with the same support and the same autocorrelation" width="800">
_Fig 3. In 1D, two non-trivially different signals give the same intensity_

The example is from Shechtman and colleagues. $u = [1, 0, -2, 0, -2]$ and $v = [1-\sqrt{3}, 0, 1, 0, 1+\sqrt{3}]$ have the same length and the same nonzero positions. Neither can be obtained from the other by shifting, flipping, or changing sign. Yet both have the autocorrelation $[-2, 0, 2, 0, 9, 0, 2, 0, -2]$, identical to the last decimal place.

Polynomials explain why. The $z$-transform of a 1D signal is a polynomial in $z$, and a polynomial in one variable always factors into linear terms. Pick one factor, flip its root to the other side of the unit circle, and a different signal with the same magnitude spectrum appears. More factors mean more such partners.

From two dimensions on, things change. Most polynomials in two variables cannot be factored. On this basis Bruck and Sodin, Hayes, and Bates showed that, apart from exceptional factorable cases, an object in two or more dimensions is determined by its intensity up to the trivial ambiguities. The exceptions form a set of measure zero, so a randomly chosen object essentially never falls into it. This is why phase retrieval works in practice for images.

## 3. How Finely Must We Sample?

A unique answer requires enough measurements. How finely must the diffraction pattern be sampled?

Start by counting. If the object fits inside $s \times s$ pixels and is complex, there are $2s^2$ unknowns. The diffraction pattern supplies $N \times N$ intensities, that is, $N^2$ real numbers. So at least $N^2 \ge 2s^2$ is needed. Defining the oversampling ratio $\sigma$ as the area of the computational box divided by the area of the support the object occupies, $\sigma \ge 2$ is a necessary condition. A real object has half as many unknowns, but its intensity is also centrosymmetric, halving the independent measurements, so the conclusion is the same.

Sampling theory gives a stricter criterion. Rodenburg and Maiden put it this way: intensity is amplitude squared, so the fastest fringe in the amplitude becomes a fringe twice as fast in the intensity. In real space it is the same statement: as the left of Fig. 4 shows, the autocorrelation is twice as wide as the object in each direction. To sample the intensity at the [Nyquist rate](/en/posts/imgproc-sampling-noise-model/), the object must fit within half the computational box along each axis, which is $\sigma \ge 4$ in area. The uniqueness result summarized by Shechtman and colleagues likewise says that $2N-1$ samples per axis suffice.

What happens in between, for $2 < \sigma < 4$? The object size was varied, with eight random starts for each size.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig4-oversampling.png" alt="Supports of the object and its autocorrelation, and reconstruction successes versus oversampling ratio" width="850">
_Fig 4. The autocorrelation is twice the object's size; near σ = 2 a solution exists but is hard to find_

| σ | Real, 600 it. | Real, 3000 it. | Complex, 600 it. | Complex, 3000 it. |
|---|---|---|---|---|
| 2.12 | 0/8 | 0/8 | 0/8 | 0/8 |
| 2.56 | 0/8 | 7/8 | 0/8 | 4/8 |
| 3.16 | 8/8 | 8/8 | 8/8 | 8/8 |
| 4.00 | 8/8 | 8/8 | 8/8 | 8/8 |

At $\sigma = 2.12$ not a single run succeeded, even with five times as many iterations. At 2.56 every run failed at 600 iterations, but at 3000 the real object was found 7 times out of 8 and the complex object 4 times. Beyond 3 the problem is easy. The counting criterion $\sigma \ge 2$ only says whether an answer can exist, not whether an algorithm can find it. Close to the boundary, the path to the answer narrows sharply, which is why twofold sampling per axis is recommended in practice. Note also that an FFT treats the box as periodically repeated, the same assumption as the circular convolution of Post 1. If the autocorrelation overflows the box, it wraps around and overlaps the opposite edge.

## 4. Alternating Between Constraints

Suppose there are enough measurements. How is the phase found? Two things are known: in Fourier space the magnitude must equal the measurement, and in real space the object must lie within the support. Rodenburg and Maiden explain why the support constraint is so powerful. Each pixel of the diffraction pattern corresponds to one plane-wave component of the object, and changing its phase shifts that component sideways. For all components to cancel to zero outside the support, their phases must interlock exactly. Get a single phase wrong and something leaks outside the support.

So the two constraints can be enforced in turn. Gerchberg and Saxton did this when the intensities in both the image and the diffraction plane were known, and Fienup generalized it into what he called the error-reduction (ER) algorithm.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig5-er-loop.png" alt="Loop diagram of the ER algorithm" width="700">
_Fig 5. Gerchberg–Saxton and ER: enforce the constraints of each space in turn_

Fourier transform the current estimate $g$, replace the magnitude with the measurement while keeping the phase, and transform back. Something will now sit outside the support; set it to zero. For a real object, set negative values to zero too. Repeat. Both steps are projections: each moves the current point to the nearest point satisfying one constraint. Writing the modulus and support projections as $P_m$ and $P_s$, ER is one line.

$$ g_{k+1} = P_s \, P_m \, g_k $$

Fienup's 1982 paper showed two things. First, the mismatch with the measurement never increases under this iteration. Second, ER is equivalent to steepest descent on that mismatch, $E = \lVert\, \lvert F g \rvert - m \,\rVert^2$. Indeed, the gradient of $E$ is $2(g - P_m g)$, and stepping half the gradient and then projecting onto the support reproduces one ER step exactly. Checked against a finite-difference gradient, the two agree to about $10^{-8}$. Like the Landweber iteration of [Post 1](/en/posts/ptycho-computational-imaging-inverse-problem/), ER is a form of [gradient descent](/en/posts/optimization-gradient-descent/).

And it shares the weakness of gradient descent: on a non-convex objective it stops at a local minimum.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig6-er-vs-hio.png" alt="Error versus iteration for ER and HIO, and reconstructions after 600 iterations" width="900">
_Fig 6. ER stalls near 0.19 almost at once; HIO wanders before reaching the answer (here, the twin)_

For the real object, the ER error falls quickly from 0.61 after one iteration to 0.28 after 10 and 0.19 after 50, then stops. After 600 iterations it is still 0.185. Out of eight random starts, none reached the answer. This is stagnation. The guarantee that the mismatch never increases does not mean the iteration never stops short.

## 5. Deliberately Not Satisfying the Constraint

How can the iteration escape stagnation? Fienup's 1982 answer was a change of viewpoint. Group the Fourier transform, the magnitude replacement, and the inverse transform into one box. Whatever goes in, what comes out is consistent with the measurement. The input then need not be the current best estimate; it can be treated as a driving signal that pushes the output in a desired direction.

Where the output violates the support constraint, we want to push it towards zero. The box's output roughly follows changes in its input, so the input at those points is nudged a little in the direction opposite to the output. This is the hybrid input-output (HIO) algorithm.

$$ g_{k+1}(r) = \begin{cases} g'_k(r) & r \in S \\ g_k(r) - \beta\, g'_k(r) & r \notin S \end{cases}, \qquad g'_k = P_m g_k $$

Inside the support it matches ER; the difference is outside. ER erases the outside every time and forgets it, whereas HIO remembers the previous input and subtracts the output from it. If the output keeps appearing outside, the input at that point keeps accumulating in the opposite direction until it pushes the output out. In Fienup's comparison, the algorithm that reduced the error most was HIO with β close to 1.

The blue curve in Fig. 6 is HIO. For the first 50 iterations its error is actually larger than ER's and fluctuates. It then falls to 0.014 at 100 iterations and 0.0018 at 200, reaching $4 \times 10^{-4}$ at 600. The reconstruction on the right is flipped: HIO found the twin image, which, as Section 1 explained, is a correct answer indistinguishable from the original by measurement. HIO succeeded from all eight random starts, for the real object and for a complex object with both amplitude and phase alike, while ER succeeded in none.

Since the HIO input is not an object estimate, its error does not fall smoothly. Fienup therefore compared algorithms by appending a few ER iterations after the HIO iterations before reading the error, and alternating a few dozen HIO iterations with one ER iteration is still common practice. ER cannot escape stagnation, but it is good at polishing an answer once one is nearby.

## 6. Finding the Intersection of Two Sets

Why does HIO work where ER fails? Recasting phase retrieval as geometry makes it visible. There is the set $S$ of objects satisfying the support constraint and the set $M$ of objects with the measured magnitude. The answer is their intersection. $S$ is convex: mixing two objects inside the support leaves the result inside the support. $M$ is not convex: mixing two complex numbers of equal magnitude produces a smaller magnitude, just as the chord between two points on a circle passes inside it.

When both sets are convex and they intersect, simply alternating projections is known to reach the intersection. When one is not convex, that guarantee disappears. A 2D toy shows how. The straight line is $S$, the curve is $M$, and they meet only at the star.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig7-sets-toy.png" alt="Trajectories of ER and HIO searching for the intersection of a line and a curve" width="900">
_Fig 7. ER is trapped where the sets come closest; HIO is pushed along the gap and finds the intersection_

ER (red) slides into the valley where the curve comes closest to the line and stops there. At that point the two projections point at each other, leaving no reason to move. The remaining gap is 0.077. HIO (blue) behaves differently. As long as a gap remains, it keeps pushing its input in the gap's direction. It is pushed downwards for 200 iterations; the moment the nearest point on the curve jumps to the other branch, it moves sideways, circles the intersection, and settles exactly on it by about iteration 500. Marchesini summarized this as HIO continuing to move in the direction of the gap.

In this framework many algorithms are combinations of projections $P$ and reflections $R = 2P - I$. ER is $P_s P_m$, and the RAAR (relaxed averaged alternating reflections) algorithm proposed by Luke is $\tfrac{\beta}{2}(R_s R_m + I) + (1-\beta) P_m$. Without a positivity constraint and with $\beta = 1$, HIO and RAAR become the same formula; computing one step both ways agrees to $10^{-16}$. Running both from the same start for 50 iterations, though, the difference grows to $10^{-6}$ because rounding errors are amplified at each iteration. The trajectories of HIO-type algorithms are that sensitive to small perturbations. It ties in with Marchesini's toy examples, where adding a tiny random number at every iteration freed stagnated runs. Elser's difference map belongs to the same family, and Thibault and colleagues used it for ptychography; it returns in Post 4.

Different combinations behave differently. RAAR proved sensitive to β. For the complex object after 1000 iterations, the number of successes was as follows.

| β | 0.75 | 0.85 | 0.90 | 0.95 | 0.99 |
|---|---|---|---|---|---|
| RAAR successes (out of 8) | 0 | 0 | 0 | 2 | 8 |

For the real object too, RAAR with β = 0.9 never succeeded in 500 iterations, while HIO with the same β succeeded every time. Among these closely related algorithms, it is hard to name one as better in general; the ranking changes with the problem and its parameters.

More recently, the same problem has been attacked directly as optimization. The Wirtinger flow of Candès and colleagues takes the intensity mismatch as its objective and descends its gradient, after carefully choosing the starting point with a spectral method. For random measurement vectors it provably converges to the answer from a nearly minimal number of measurements. That guarantee, however, concerns random measurements or several diffraction patterns coded with different masks. For a single diffraction pattern, as in this post, no such proof exists. The limits of a single pattern are not the kind that a different algorithm overcomes.

## 7. Noise and the Twin

So far the measurement has been perfect. What about real measurements that count photons? Poisson noise was added to the diffraction intensity, and the reconstruction repeated for different total photon counts.

<img src="/assets/img/posts/ptycho-phase-problem-iterative-retrieval/en/fig8-noise.png" alt="Reconstruction error versus total photon count, and amplitude and phase of the reconstructions" width="900">
_Fig 8. A single pattern has no redundancy, so the reconstruction collapses quickly as photons drop_

With $10^9$ photons the error is 0.008. It grows to 0.032 at $10^8$, 0.16 at $10^7$, and 0.31 at $10^6$. A diffraction pattern is bright at the centre and dark towards the edges. As photons decrease, the high spatial frequencies at the edges are the first to sink into noise, and no other measurement exists to restore them. In the language of Post 1, the measurement has no redundancy. Methods that hold on by adding assumptions, like the regularization of Post 1, do exist, but fundamentally the information in a single pattern is fixed.

Ambiguity also remains. If the support is centrosymmetric, the original and the twin satisfy the support constraint equally well. The reconstruction may fail to commit to either and hover in a mixture of the two. Rodenburg and Maiden describe the two solutions as "fighting". In their example the remedy was simple. Recording and reconstructing diffraction patterns from four separate regions of an object independently produced reconstructions contaminated by the twin. Making the four regions overlap slightly, and requiring the overlapping parts to agree, removed the ambiguity, because no solution can be the original in one region and the twin in its neighbour at the same time.

## Summary

Measuring intensity loses the phase, and what remains is the autocorrelation. Objects with the same autocorrelation cannot be distinguished: translation, the twin image, and global phase always, and in 1D other partners as well. In two or more dimensions most objects are unique up to the trivial ambiguities, but only if the intensity is sampled finely enough. The counting criterion is $\sigma \ge 2$ and the sampling criterion is twofold per axis. In between, an answer may exist yet be hard to find.

Finding the answer started from enforcing the two constraints in turn. ER is gradient descent and never increases the mismatch with the measurement, but the non-convex modulus set makes it stagnate. HIO treats the input as a driving signal and keeps pushing along the gap until it escapes. Seen as combinations of projections and reflections, these algorithms form one family whose success depends on their parameters. And a single pattern has no slack against noise or the twin image.

## Next

The last example of Section 7 is the starting point of the next post. Moving the illumination and recording diffraction patterns from overlapping regions turns a problem that was marginal with one pattern into a comfortable one. Measurements far outnumber unknowns, overlapping parts are guarded jointly by neighbouring patterns, and the support no longer needs to be known in advance. In the handbook's words, oversampling is not a fundamental constraint in ptychography. Post 3 builds from scratch the PIE (ptychographical iterative engine), the first iterative implementation of this idea, and its refinement ePIE.

## References

- J. R. Fienup, "[Phase retrieval algorithms: a comparison](https://doi.org/10.1364/AO.21.002758)," *Applied Optics* 21(15), 2758-2769 (1982). The description of the Gerchberg–Saxton algorithm (R. W. Gerchberg, W. O. Saxton, *Optik* 35, 237, 1972) also follows this paper.
- J. M. Rodenburg, A. M. Maiden, "[Ptychography](https://eprints.whiterose.ac.uk/id/eprint/127795/)," in *Springer Handbook of Microscopy*, P. W. Hawkes, J. C. H. Spence (eds.), Springer, 2019, pp. 819-904, [doi:10.1007/978-3-030-00069-1_17](https://doi.org/10.1007/978-3-030-00069-1_17), Secs. 3.1–3.3 (intuition for the support constraint, intensity sampling, twin image and overlap).
- Y. Shechtman, Y. C. Eldar, O. Cohen, H. N. Chapman, J. Miao, M. Segev, "[Phase retrieval with application to optical imaging: a contemporary overview](https://doi.org/10.1109/MSP.2014.2352673)," *IEEE Signal Processing Magazine* 32(3), 87-109 (2015). Source of the trivial ambiguities, the 1D counterexample, 2D uniqueness (Bruck–Sodin, Hayes, Bates), and the phase-swap experiment.
- S. Marchesini, "[A unified evaluation of iterative projection algorithms for phase retrieval](https://doi.org/10.1063/1.2403783)," *Review of Scientific Instruments* 78, 011301 (2007).
- D. R. Luke, "[Relaxed averaged alternating reflections for diffraction imaging](https://doi.org/10.1088/0266-5611/21/1/004)," *Inverse Problems* 21, 37-50 (2005).
- E. J. Candès, X. Li, M. Soltanolkotabi, "[Phase retrieval via Wirtinger flow: theory and algorithms](https://doi.org/10.1109/TIT.2015.2399924)," *IEEE Transactions on Information Theory* 61(4), 1985-2007 (2015).
- P. Thibault, M. Dierolf, A. Menzel, O. Bunk, C. David, F. Pfeiffer, "[High-resolution scanning X-ray diffraction microscopy](https://doi.org/10.1126/science.1158573)," *Science* 321, 379-382 (2008) (difference map applied to ptychography).

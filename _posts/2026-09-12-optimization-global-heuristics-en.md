---
title: "Optimization 4 — Global Heuristics: Basin-Hopping and Simulated Annealing"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-global-heuristics/
page_id: optimization-global-heuristics
date: 2026-09-12 20:00:00 +0900
categories: [Computation, Optimization]
tags: [optimization, simulated-annealing, basin-hopping, global-optimization, metropolis, thin-film]
description: How Simulated Annealing and Basin-Hopping cross the basin boundaries that local optimization cannot, and how to choose the temperature — tested experimentally.
math: true
---

[Post 3](/en/posts/optimization-levenberg-marquardt/) showed that Levenberg-Marquardt (LM), by damping the step through the damping parameter $\mu$, prevents the divergence that comes from overly large Gauss-Newton steps; it reproduced this with a Jacobian close to zero. The slow convergence of large-residual problems, however, remains. And as the experiments at the end of that post established, LM has a further limit that $\mu$ cannot repair: starting from $d_0=1300$ nm it always converged near $d=1293$ nm, and from $d_0=1690$ nm always near $d=1689$ nm. Changing $\tau$, the initial damping scale, made no difference. Local methods all descend only the valley they happen to stand in; whether a deeper floor lies in the next valley is invisible to them. What, then, has to be added for a local optimizer to see the next valley?

This post takes on that basin problem directly; a basin is the valley that drains into one minimum. Two classical heuristics leave the local optimizer itself alone and instead layer randomness and a probabilistic acceptance rule on top of it, so the search can leave a bad valley on its own: simulated annealing (SA) and basin-hopping. Both are implemented here and tested to see how well they actually work. The landscape has a single unknown, the thickness, so it is one-dimensional; the many-parameter case gets a brief word in Section 6. Where posts 1 through 3 were a series about local methods, this one is the global search layer wrapped around them.

## 1. The three basins of the objective landscape

Take the objective $J(d)$ of the running example again (single SiO2/Si film, true thickness 1490 nm, seed=0 noise). Replotting Figure 2 of [post 1](/en/posts/optimization-gradient-descent/) over 1200 to 1800 nm on a log scale, with the minima and boundaries marked, shows how many basins there are and where their boundaries lie.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig1-objective-landscape.png" alt="The objective landscape and its three basins" width="600">
_Fig 1. The objective landscape (1200–1800 nm, log scale). The red dots are the three local minima and the dashed lines the basin boundaries. Only the central floor at 1490 nm is more than two orders of magnitude deeper._

There are three local minima in this range: $d\approx1293$ nm ($J\approx0.867$), $d\approx1490$ nm ($J\approx0.0024$, the global minimum), and $d\approx1689$ nm ($J\approx1.059$). The basin boundaries, the peaks of the landscape, sit at $d\approx1391$ nm and $d\approx1589$ nm, matching the "roughly 1390 and 1590" estimated in post 1. The boundary peaks reach $J\approx5.12$ and $J\approx4.88$, about 4.25 above the 1293 nm floor. It is also worth noting that the global minimum's $J$ is more than two orders of magnitude below the other two. If the initial position lies outside that 198 nm-wide basin, no amount of careful step adjustment lets LM leave the valley it is standing in.

## 2. Simulated Annealing (SA)

### 2.1 The principle

A local optimizer cannot leave its valley because it never takes a single step uphill. So why not let it go the wrong way now and then, on purpose? Simulated annealing turns that idea into one probabilistic rule.

Proposed by Kirkpatrick, Gelatt and Vecchi (1983), it takes its name from annealing in metallurgy: cool a metal slowly and its atoms settle into the lowest-energy crystal structure. The algorithm is simple. From the current position $d$, a randomly proposed $d'$ is accepted or rejected by the Metropolis criterion.

$$ P(\text{accept}) = \min\!\left(1,\ \exp(-\Delta J / T)\right), \qquad \Delta J = J(d') - J(d) $$

Here $\Delta J$ is the difference in the objective between the proposed and the current position, and the temperature $T$ is the handle on the acceptance probability. Proposals that decrease $J$ are always accepted; proposals that increase it are accepted with some probability. That willingness to occasionally move the wrong way is what lets the search escape a shallow valley instead of being trapped in it. As $T\to0$ only improving proposals are accepted and the method becomes pure greedy descent; as $T\to\infty$ it becomes effectively a random walk. $T$ is cooled geometrically, multiplied by `cooling` at every step ($T_{k+1}=\text{cooling}\cdot T_k$), so the search is broad at first and increasingly settles into one valley.

The shape of the formula tells you how to set the temperature. Since $\Delta J/T$ sits in the exponent, $T$ is in the units of $J$ and has to be judged against the $\Delta J$ that must be overcome. On this landscape the boundary peak stands about 4.25 above the 1293 nm floor, while floors differ from one another by 0.19 (1293 nm against 1689 nm) to 0.86 (1490 nm against 1293 nm). At $T=0.02$ the chance of climbing all the way to the peak is $e^{-212}$, effectively zero. At $T=5$ it is 0.43, and a proposal that leaves the global floor for the 1293 nm floor is accepted with probability 0.84. In other words, a $T$ well below the floor-to-floor differences all but closes off uphill moves, and a $T$ approaching the barrier height turns the search into little more than a random walk.

### 2.2 Applied to this problem

With thickness $d$ as the only parameter, a plain Gaussian random walk suffices as the proposal distribution. Its standard deviation (`step_sigma`) is 80 nm.

```python
# core of simulated_annealing.py (full code: _code/optimization-global-heuristics/)
T = T0
for _ in range(n_iter):
    d_prop = d + rng.normal(0, step_sigma)
    J_prop = objective(d_prop, wavelength_nm, measured_R)
    dJ = J_prop - J
    if dJ < 0 or rng.random() < np.exp(-dJ / T):
        d, J = d_prop, J_prop
    if J < best_J:
        best_d, best_J = d, J
    T *= cooling
```

The LM code of post 3 is not used here: SA crosses basins through its acceptance rule alone, with no local optimizer. What it does track separately is the lowest point visited so far (`best_d`, `best_J`) rather than the current position. Since the Metropolis criterion probabilistically accepts bad steps, there is no guarantee that where the walk currently stands is the best place it has been.

The 80 nm proposal width already predicts the next section's result. From the 1293 nm floor, the only place where $J$ drops below that floor's value (0.867) is the interval 1463–1517 nm. A single proposal lands there with probability about 1.4%. Such a proposal lowers $J$, so it is accepted whatever the temperature. The barrier is thus jumped in one go, with no uphill step at all, at a 1.4% chance per step, and the median number of steps to the first success is about 49.

### 2.3 Does changing the initial temperature change the outcome?

From $d_0=1300$ nm, the very starting point at which LM failed in post 3, the cooling rate was fixed (`cooling=0.97`) and only the initial temperature $T_0$ varied across three values, for 300 iterations. The short answer: on this problem the initial temperature does not decide the outcome, as long as it is not too high.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig2-sa-temperature-schedule.png" alt="SA search trajectories by temperature schedule" width="600">
_Fig 2. SA search trajectories for three initial temperatures T0 (d0 = 1300 nm, cooling 0.97, proposal standard deviation 80 nm, 300 iterations). All three reach the global minimum (dashed line, 1490 nm), but when the boundary is crossed depends on the random stream more than on the temperature._

With a very low $T_0$ (0.0005) the walk stays essentially fixed at $d\approx1293$ for the first 240-odd steps, since non-improving proposals are almost never accepted even probabilistically. But even at low probability a proposal that crosses the basin boundary eventually appears: at step 245 it crosses, at the next step it is within 5 nm of the global minimum, and it settles there immediately. With worse luck it might not have crossed within 300 steps at all.

With a moderate $T_0$ (0.02) the walk crosses the boundary at step 64 and is within 5 nm of the global minimum by step 74, settling there stably. Rerun over 200 seeds, the median step at which the boundary is crossed is 47 for both $T_0=0.0005$ and $T_0=0.02$, in line with the 49 predicted in Section 2.2. The gap between 64 and 245 is the luck of the random stream, not the temperature.

With a very high $T_0$ (5.0) the estimate wanders broadly between 1240 and 1722 nm for the first 110 or so steps. It does end up near the global minimum (lowest visited $J\approx0.0025$), but the trajectory is untidy and repeatedly re-enters several basins, with a real risk of never settling into one valley before the cooling schedule runs out.

All three found the global minimum within this 300-step budget, but widening the test to 200 seeds shows where they differ. A $T_0$ that is too high leaves the search unable to stay where it has already arrived, and the success rate falls to 141/200. $T_0=0.0005$ and $0.02$, by contrast, both reached 199/200: as Section 2.2 showed, a boundary on this landscape is crossed by one jump that lowers $J$, so temperature has nothing to contribute. A large temperature is not automatically a good one, and on this problem nothing was lost by choosing it low.

## 3. Basin-Hopping

### 3.1 The principle

SA compares arbitrary points at every step, so it keeps taking small steps even near a valley floor. What if the comparison were made between valley floors from the start? Basin-hopping has the same aim as SA, and this is exactly where it differs. Where SA accepts or rejects each coordinate $d$ by the Metropolis criterion, basin-hopping accepts or rejects the answer to "which basin's floor is the search on?"

1. Apply a large random perturbation to the current solution: $d_0' = d + \mathcal{N}(0,\sigma^2)$
2. Run a local optimization from $d_0'$ to convergence, finding that basin's floor $d'$
3. Compare the two floors, $J(d)$ against $J(d')$, by the Metropolis criterion, and accept or reject

Why this moves so cleanly becomes clear from the landscape basin-hopping effectively sees. Wales and Doye described the method as a transformation of the landscape: each point $d$ is replaced by the value of the floor that local optimization reaches from it. Applied to the landscape of Figure 1, every point between 1391 nm and 1589 nm takes the same value, 0.0024, because LM started anywhere in that range descends to the 1490 nm floor. The left basin becomes 0.867 and the right 1.059. The landscape turns into flat steps, one per basin, and the barriers between them vanish. The Metropolis criterion in basin-hopping therefore sees only the height of the steps (floor differences of 0.19–0.86), not the height of the barriers (about 4.25). It is unmoved by fine ripples in the landscape, such as those the noise creates, and moves far more cleanly in units of whole basins.

The name basin-hopping was given by Wales and Doye (1997). The method itself is essentially the Monte Carlo-minimization of Li and Scheraga (1987), as Wales and Doye themselves note.

### 3.2 An implementation with the LM of post 3 as the local step

The local optimization of step 2 is precisely the `levenberg_marquardt()` of post 3. The whole series converges here: post 3 shored up, with damping, the conditions under which the methods of posts 1 and 2 break down, and this post takes that same function wholesale and reuses it as the tool for descending one basin quickly and accurately.

```python
# core of basin_hopping.py (full code: _code/optimization-global-heuristics/)
d_hist_lm, J_hist_lm, _ = levenberg_marquardt(d0, wavelength_nm, measured_R, **lm_kwargs)
d, J = d_hist_lm[-1], J_hist_lm[-1]
T = T0
for _ in range(n_hops):
    d_trial0 = d + rng.normal(0, perturb_sigma)
    d_hist_trial, J_hist_trial, _ = levenberg_marquardt(d_trial0, wavelength_nm, measured_R, **lm_kwargs)
    d_trial, J_trial = d_hist_trial[-1], J_hist_trial[-1]
    dJ = J_trial - J
    if dJ < 0 or rng.random() < np.exp(-dJ / T):
        d, J = d_trial, J_trial
    if J < best_J:
        best_d, best_J = d, J
```

The perturbation size (`perturb_sigma`) has to be comparable to the basin spacing (about 200 nm, from Figure 1) for any chance of landing in a different basin. It is set to 150 nm here; since the neighbouring boundary is about 100 nm from the 1293 nm floor, each hop (one perturbation plus one local optimization) then has roughly a 23% chance of landing in the global basin.

### 3.3 Why does a high temperature fail?

Twenty-five hops were attempted from the same $d_0=1300$ nm, at a fixed temperature (no cooling).

<img src="/assets/img/posts/optimization-global-heuristics/en/fig3-basinhopping-temperature-schedule.png" alt="Basin-hopping trajectories by temperature" width="600">
_Fig 3. Basin-hopping trajectories at fixed temperatures (d0 = 1300 nm, perturbation standard deviation 150 nm, 25 hops, no cooling). T0 = 0.0005 and 0.05 give the same trajectory and overlap; T0 = 5 is pushed into basins worse than the starting one._

The very low temperature (0.0005) and the moderate one (0.05) give identical trajectories (the blue line is hidden behind the green). Both accept 17 of 25 hops, and both jump from $d=1293$ to $d=1490$ at the fourteenth hop and stay there. Within the 1200–1800 nm search range there are only three basins, widely spaced, so the low-temperature strategy of accepting only improving hops works perfectly well on this problem.

At a very high temperature (5.0), 24 of the 25 hops are accepted, effectively all of them. The result is a sequence $d=1293 \to 1022 \to 696 \to 530 \to 441$, thrown into basins all worse than the starting one ($J\approx2$–2.6). It never even passes near the global minimum, so `best_J` remains at the value of the starting basin, 0.867. The acceptance probability shows why. A hop from the 1293 nm floor to the 1022 nm floor has $\Delta J\approx1.10$; at $T=5$ it is accepted with probability 0.80, at $T=0.05$ with $3\times10^{-10}$.

I had expected that too low a temperature would reduce basin-hopping to multi-start. On this problem that was wrong; the lower temperature in fact worked better. The real obstacle was the opposite: at too high a temperature nearly every move into a worse basin is accepted, and the search ends up wandering among basins worse than the one it started in. This run went as far as 441 nm, outside the search range. The acceptance rule may allow such a hop, but physically it is a thickness outside the prior range. Bounding the search to the prior range prevents this kind of escape.

This difference cannot, however, be put down to the structure of SA versus basin-hopping. The basin-hopping here runs at a fixed temperature (`cooling=1.0`; Wales and Doye's original paper and scipy also keep T fixed), whereas SA cools to below one ten-thousandth of $T_0$ over its 300 steps. With the conditions matched over 200 seeds, holding both at $T=5$ leaves the final position within 5 nm of the global minimum in 1/200 SA runs and 20/200 basin-hopping runs, and cooling basin-hopping by 0.7 per hop raises it to 174/200. What prevents the loss at high temperature is cooling, not the choice of algorithm.

## 4. Over 50 random starts, how often and in how many steps?

Rather than looking at a single starting point, $d_0$ was drawn at random from $[1200, 1800]$ fifty times and the success rates of the three methods compared. "Success" is defined as reaching within 5 nm of the global minimum at 1490 nm. SA used $T_0=0.02$ with cooling 0.97; basin-hopping used a fixed $T_0=0.05$.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig4-success-rate-comparison.png" alt="Success rate over 50 random starts" width="600">
_Fig 4. Fraction of 50 random starts that reach within 5 nm of the global minimum (SA: T0 = 0.02, 300 iterations; basin-hopping: fixed T0 = 0.05, 15 hops). The 36% of a single LM run is close to the share of the search range taken up by the correct basin (33%)._

A single LM, with no retries, managed 18/50 = 36%, close to the 33% share that the basin width (198 nm) occupies of the whole search range (600 nm). That is exactly the limit of a local method that succeeds only when the initial value already lies in the right basin. SA, on a budget of 300 iterations, reached 50/50 = 100%, taking on average 61.8 iterations to first come within 5 nm. Basin-hopping, on a budget of 15 hops, also reached 50/50 = 100%, taking on average 2.2 hops.

Neither global method failed within these budgets. Iteration counts and hop counts must not be compared directly, though: one SA step is a single evaluation of the objective, whereas one basin-hopping hop contains an entire LM run (several iterations, Jacobians included) descending to a floor. Basin-hopping therefore needs far fewer outer steps, while SA keeps taking small refining steps even near a valley floor. A common yardstick is the number of reflectance (forward-model) evaluations. In the 50 runs above, the mean per starting point was about 121 for a single LM, 301 for SA and about 1,790 for basin-hopping. Basin-hopping takes fewer outer steps but costs about six times as much computation as SA. The counts that [post 5](/en/posts/optimization-nn-regression/) measured under different conditions (100 random starts; SA 301, basin-hopping 1,806) show the same ratio.

## 5. Basin-hopping seen on the landscape

Plotting on the landscape of Figure 1 the trial floor that LM reached at every hop of the $T_0=0.05$ run in Figure 3 makes it plain that what this algorithm does is compare floors. Some trial floors fall outside the search range, so the landscape is drawn over 1000–1900 nm.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig5-trajectory-on-landscape.png" alt="Basin-hopping trial floors on the landscape" width="600">
_Fig 5. Trial floors that LM reached at each hop of the T0 = 0.05 run in Figure 3 (grey band: the 1200–1800 nm search range). Filled dots mark floors with accepted tries, open dots floors whose tries were all rejected; the numbers give the tries and acceptances at each floor. Only one hop, the fourteenth, actually moved the search (arrow)._

The blue square is the first floor reached by LM from $d_0=1300$ ($d=1293$); the green star is the global minimum, moved to at the fourteenth hop. The 25 tries landed on six floors. Thirteen came back to the starting 1293 nm floor: eleven of them before the move, accepted with $\Delta J\approx0$, and two after the search had moved to 1490 nm, rejected with $\Delta J\approx0.86$. The 1689 nm floor (two tries) and three floors outside the search range (1022 nm, 1130 nm and near 1860 nm, four tries in all) were higher than the current floor and were all rejected. The only hop that actually moved the search was the fourteenth, the first to land on the 1490 nm floor. Basin-hopping does not sweep the landscape; it probes several floors and moves only when a lower one turns up. LM alone would have stayed at the blue square forever.

## 6. Summary: when is a local method enough?

Everything refined across posts 1 to 3 (gradient descent, Newton, Gauss-Newton, LM) rested on the premise of already being in the right basin. As this post shows, the moment that premise breaks, no amount of careful damping helps a local method: it has no way of seeing another valley in the first place.

In practice, four questions usually settle which to use.

| Question | When | Use | Why | Numbers in this post |
| --- | --- | --- | --- | --- |
| Prior information on the initial value | A previous measurement, process spec or design value narrows the initial guess to within a basin width | A single LM run | It starts inside the correct basin | Basin width 198 nm; single LM from random starts 36% |
| Cost of the objective | Milliseconds per evaluation (transfer matrix of a multilayer stack) or seconds (FDTD) | Global search only when needed | Global search makes several to over ten times as many calls | Calls per start: LM 121, SA 301, basin-hopping about 1,790 |
| Number and spacing of basins | Many, closely spaced basins (especially with several parameters) | Set temperature and perturbation with care; add cooling if needed | The high-temperature failure of wandering among worse basins appears readily | Three basins about 200 nm apart: a low temperature is enough |
| Local structure and local optimizer | A smooth landscape and a trustworthy local optimizer (LM) | Basin-hopping | It compares floors and moves basin by basin | Reached in 2.2 hops on average |
| Local structure and local optimizer | A rough landscape, or no good local optimizer | SA | It needs no local optimization stage and makes fewer calls | 61.8 iterations on average, 301 calls |

The first row is the most common answer. The basin width in this example was 198 nm, and having no prior information even at that precision is rarer than one might think. When a global search really is needed and the success rates are equal, SA is the cheaper of the two. Basin-hopping gives a cleaner trajectory and fewer hops but pays for a local optimization at every hop.

From gradient descent in post 1 to basin-hopping here, this series has worked step by step through a trade-off between two axes: how cleverly local information (gradient, curvature) is used, and how widely global structure is surveyed. Real metrology software commonly uses both: a global search such as a grid search to find the approximate location, and a local method such as LM for the final precision.

## References

- S. Kirkpatrick, C. D. Gelatt, M. P. Vecchi, "Optimization by Simulated Annealing," _Science_, 220(4598), 671-680, 1983.
- D. J. Wales, J. P. K. Doye, "Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms," _J. Phys. Chem. A_, 101(28), 5111-5116, 1997.
- C. P. Chang, Y. H. Lee, S. Y. Wu, "Optimization of a thin-film multilayer design by use of the generalized simulated-annealing method," _Opt. Lett._, 15(11), 595-597, 1990.
- Z. Li, H. A. Scheraga, "[Monte Carlo-minimization approach to the multiple-minima problem in protein folding](https://doi.org/10.1073/pnas.84.19.6611)," _Proc. Natl. Acad. Sci. USA_, 84(19), 6611-6615, 1987.
- [scipy.optimize.basinhopping documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html)

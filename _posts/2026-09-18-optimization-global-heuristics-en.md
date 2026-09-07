---
title: "Optimization 4 — Global Heuristics: Basin-Hopping and Simulated Annealing"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-global-heuristics/
page_id: optimization-global-heuristics
date: 2026-09-18 20:00:00 +0900
categories: [Computation, Optimization]
tags: [optimization, simulated-annealing, basin-hopping, global-optimization, metropolis, thin-film]
description: How Simulated Annealing and Basin-Hopping cross the basin boundaries that local optimization cannot, using a single temperature parameter — tested experimentally.
math: true
---

[Post 3](/en/posts/optimization-levenberg-marquardt/) showed that Levenberg-Marquardt (LM) rescues Gauss-Newton's failure conditions — large residuals, or a near-singular Jacobian — through the damping parameter $\mu$. But as the experiments at the end of that post established, LM has a limit that $\mu$ cannot repair: starting from $d_0=1300$ nm it always converged near $d=1293$ nm, and from $d_0=1690$ nm always near $d=1688$ nm. Changing $\tau$, the initial damping scale, made no difference. Local methods all descend only the valley they happen to stand in; whether a deeper floor lies in the next valley is invisible to them.

This post takes on that basin problem directly. Rather than fixing the local optimizer itself, it layers randomness and a probabilistic acceptance rule on top of it, so the search can leave a bad valley on its own. Two classical heuristics — Simulated Annealing (SA) and Basin-Hopping (BH) — are implemented and tested to see how well they actually work. Where posts 1 through 3 were a series about local methods, this one is the global search layer wrapped around them.

## 1. Looking at the landscape again — three basins

Plotting the objective $J(d)$ of the running example (single SiO2/Si film, true thickness 1490 nm, seed=0 noise) over a wider range than Figure 2 of [post 1](/en/posts/optimization-gradient-descent/) — 1200 to 1800 nm — makes it possible to count the basins exactly and see where the boundaries lie.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig1-objective-landscape.png" alt="The objective landscape and its three basins" width="600">
_Fig 1. The objective landscape — three basins (dashed lines: basin boundaries)_

There are exactly three local minima in this range: $d\approx1293$ nm ($J\approx0.867$), $d\approx1490$ nm ($J\approx0.0024$, the global minimum), and $d\approx1689$ nm ($J\approx1.059$). The basin boundaries — the peaks of the landscape — sit at $d\approx1391$ nm and $d\approx1589$ nm, matching exactly the "roughly 1390 and 1590" estimated in post 1. It is also worth noting that the global minimum's $J$ is more than two orders of magnitude below the other two. If the initial position lies outside that 198 nm-wide window, no amount of careful step adjustment lets LM leave the valley it is standing in.

## 2. Simulated Annealing — annealing the random walk itself

### 2.1 The principle

SA, proposed by Kirkpatrick, Gelatt and Vecchi (1983), takes its name from physical annealing: cool a metal slowly and its atoms settle into the globally most stable — lowest-energy — crystal structure. The algorithm is simple. From the current position $d$, a randomly proposed $d'$ is accepted or rejected by the Metropolis criterion.

$$ P(\text{accept}) = \min\!\left(1,\ \exp(-\Delta J / T)\right), \qquad \Delta J = J(d') - J(d) $$

Proposals that decrease $J$ are always accepted; proposals that increase it are accepted with some probability. That willingness to occasionally move the wrong way is what lets the search escape a shallow valley instead of being trapped in it. The temperature $T$ is the handle on that probability. As $T\to0$ only improving proposals are accepted and the method becomes pure greedy descent; as $T\to\infty$ it becomes effectively a random walk. Cooling $T$ geometrically at every step ($T_{k+1}=\text{cooling}\cdot T_k$) makes the search broad at first and increasingly settled into a single valley.

### 2.2 Applied to this problem

With thickness $d$ as the only parameter, a plain Gaussian random walk suffices as the proposal distribution.

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

The LM code of post 3 is not used here — SA crosses basins by annealing alone, with no local optimizer. What it does track separately is the lowest point visited so far (`best_d`, `best_J`) rather than the current position. Since the Metropolis criterion probabilistically accepts bad steps, there is no guarantee that where the walk currently stands is the best place it has been.

### 2.3 The experiment — is the temperature schedule really the issue?

From $d_0=1300$ nm, the very starting point at which LM failed in post 3, the cooling rate was fixed (`cooling=0.97`) and only the initial temperature $T_0$ varied across three values, for 300 iterations.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig2-sa-temperature-schedule.png" alt="SA search trajectories by temperature schedule" width="600">
_Fig 2. SA search trajectories by temperature schedule (d0=1300nm)_

With a very low $T_0$ (0.0005) the walk stays essentially fixed at $d\approx1293$ for the first 240-odd steps, since non-improving proposals are almost never accepted even probabilistically. But even at low probability a proposal that crosses the basin boundary eventually appears: at step 245 it crosses, at the next step it is within 5 nm of the global minimum, and it settles there immediately. With worse luck it might not have crossed within 300 steps at all.

With a moderate $T_0$ (0.02) the walk attempts moves into the neighbouring basin far more often from the start. It crosses the boundary at step 64 and is within 5 nm of the global minimum by step 74, settling there stably.

With a very high $T_0$ (5.0) the estimate wanders broadly between 1240 and 1722 nm. It does end up near the global minimum (final $J\approx0.0025$), but the trajectory is untidy and repeatedly re-enters several basins — with a real risk of never settling into one valley before the cooling schedule runs out.

All three found the global minimum within this 300-step budget, but a $T_0$ that is too low leaves the search waiting on luck to escape, and one that is too high leaves it unable to stay where it has already arrived. A good $T_0$ is a balance between exploration and convergence, not a matter of simply picking something large.

## 3. Basin-Hopping — annealing the local optimization

### 3.1 The principle

Basin-hopping, proposed by Wales and Doye (1997), has the same aim as SA but a different structure. Instead of annealing the coordinate $d$ itself at every step, it anneals the question of which basin's floor the search currently occupies.

1. Apply a large random perturbation to the current solution: $d_0' = d + \mathcal{N}(0,\sigma^2)$
2. Run a local optimization from $d_0'$ to convergence, finding that basin's floor $d'$
3. Compare the two floors, $J(d)$ against $J(d')$, by the Metropolis criterion, and accept or reject

Where SA compares arbitrary positions — not necessarily floors — at every step, basin-hopping always compares floor against floor. It is therefore unmoved by fine ripples in the landscape, such as those the noise creates, and moves far more cleanly in units of whole basins.

### 3.2 Implementation — reusing the LM of post 3 as the local step

The "local optimization" of step 2 is precisely the `levenberg_marquardt()` of post 3. The whole series converges here: posts 1 and 2 failed at a point that post 3 rescued with damping, and this post takes that same function wholesale and reuses it as the tool for descending one basin quickly and accurately.

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

The perturbation size (`perturb_sigma`) has to exceed the basin spacing — about 200 nm, from Figure 1 — for any chance of landing in a different basin. It is set to 150 nm here.

### 3.3 The experiment — too high a temperature makes it fail

Twenty-five hops were attempted from the same $d_0=1300$ nm.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig3-basinhopping-temperature-schedule.png" alt="Basin-hopping trajectories by temperature schedule" width="600">
_Fig 3. Basin-hopping trajectories by temperature schedule (d0=1300nm)_

The very low temperature (0.0005) and the moderate one (0.05) give identical trajectories — the blue line is hidden behind the green. Both accept 17 of 25 hops, and both jump from $d=1293$ to $d=1490$ at the fourteenth hop and stay there. With only three basins, widely spaced, the low-temperature strategy of accepting only improving hops works perfectly well on this problem.

At a very high temperature (5.0), 24 of the 25 hops are accepted — effectively all of them. The result is a sequence $d=1293 \to 1022 \to 696 \to 530 \to 441$, thrown into worse and worse basins. It never even passes near the global minimum, so `best_J` remains at the value of the starting basin, 0.867.

This runs contrary to the initial expectation. The guess that "too low a temperature will reduce it to multi-start" turned out to be wrong here — the lower temperature in fact worked better. The real obstacle was the opposite: at too high a temperature, basin-hopping accepts nearly every move into a worse basin, kicking away a good solution it had already found. SA and basin-hopping share the word "temperature," but where high temperature in SA merely widens the search, in basin-hopping it translates into the far more direct loss of an answer already in hand, because what Metropolis compares are basin floors and each step therefore drops much further.

## 4. Quantitative comparison — how often, and in how many steps

Rather than looking at a single starting point, $d_0$ was drawn at random from $[1200, 1800]$ fifty times and the success rates of the three methods compared. "Success" is defined as reaching within 5 nm of the global minimum at 1490 nm.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig4-success-rate-comparison.png" alt="Success rate over 50 random starts" width="600">
_Fig 4. Success rate over 50 random starts_

A single LM, with no retries, managed 18/50 = 36%. That is roughly the fraction the basin width (198 nm) occupies of the whole search range (600 nm) — exactly the limit of a local method that succeeds only when the initial value already lies in the right basin. SA, on a budget of 300 iterations, reached 50/50 = 100%, taking on average 61.8 iterations to first come within 5 nm. Basin-hopping, on a budget of 15 hops, also reached 50/50 = 100%, taking on average 2.2 hops.

Neither global method failed within these budgets. Iteration counts and hop counts must not be compared directly, though: one SA step is a single evaluation of the objective, whereas one basin-hopping hop contains an entire LM run — several iterations, Jacobians included — descending to a floor. So on a problem like this one, where the local structure is clear (smooth valleys, an accurate local optimizer), basin-hopping finishes in far fewer outer steps, while SA needs many more steps because it keeps taking small refining steps even near a valley floor.

## 5. Basin-hopping seen on the landscape

Overlaying the $T_0=0.05$ trajectory of Figure 3 onto the landscape of Figure 1 makes it obvious at a glance that what this algorithm does is hop directly between local minima.

<img src="/assets/img/posts/optimization-global-heuristics/en/fig5-trajectory-on-landscape.png" alt="A basin-hopping trajectory over the landscape" width="600">
_Fig 5. A basin-hopping trajectory overlaid on the landscape_

The blue square is the first floor reached by LM from $d_0=1300$ ($d=1293$); the green star is the global minimum, reached after fourteen hops and held from then on. LM alone would have stayed at the blue square forever. Basin-hopping keeps trying, through random perturbation, to get over another peak, and moves there whenever Metropolis judges the far side better — it does not sweep the landscape, it hops over a few peaks and compares floors.

## 6. Summary — when is a local method enough?

Everything refined across posts 1 to 3 — gradient descent, Newton, Gauss-Newton, LM — rested on the premise of already being in the right basin. As this post shows, the moment that premise breaks, no amount of careful damping helps a local method: it has no way of seeing another valley in the first place.

In practice, a few questions usually settle which to use. The first is whether prior information about the initial value exists. If a previous measurement, a process specification, or an approximate hardware design value can narrow the initial guess to within a basin width, a single LM run suffices. The basin width in this example was 198 nm — and having no prior information even at that precision is rarer than one might think.

The second is how expensive the objective is to evaluate. SA and basin-hopping both evaluate the objective, and more, far more often than a local method. Reflectance was effectively free here, so the cost did not matter; for a transfer-matrix calculation over a multilayer stack, or an FDTD model taking seconds per evaluation, the cost of a global search is not negligible.

Knowing how many basins there are and how far apart they lie matters too. This experiment had only three, widely spaced, which is why even a low temperature worked. With many closely spaced basins — particularly as the number of parameters grows — the temperature and perturbation size need more care, and the failure mode seen in Figure 3, where high temperature loses the answer, appears more readily.

Between SA and basin-hopping: where the local structure is clear and a trustworthy local optimizer such as LM is available, basin-hopping finishes in fewer outer steps. Where the landscape is not smooth, or a good local optimizer is hard to design, SA — which works without a local optimization stage at all — is the simpler and safer choice.

From gradient descent in post 1 to basin-hopping here, this series has worked step by step through a trade-off between two axes: how cleverly local information (gradient, curvature) is used, and how widely global structure is surveyed. Real metrology software commonly uses both — a global search to find the approximate location, a local method such as LM for the final precision.

## References

- S. Kirkpatrick, C. D. Gelatt, M. P. Vecchi, "Optimization by Simulated Annealing," _Science_, 220(4598), 671-680, 1983.
- D. J. Wales, J. P. K. Doye, "Global Optimization by Basin-Hopping and the Lowest Energy Structures of Lennard-Jones Clusters Containing up to 110 Atoms," _J. Phys. Chem. A_, 101(28), 5111-5116, 1997.
- C. P. Chang, Y. H. Lee, S. Y. Wu, "Optimization of a thin-film multilayer design by use of the generalized simulated-annealing method," _Opt. Lett._, 15(11), 595-597, 1990.
- [scipy.optimize.basinhopping documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.basinhopping.html)

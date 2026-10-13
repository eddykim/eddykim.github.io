---
title: "Optimization 3 — Why Does Most Metrology Software Use the Levenberg-Marquardt Method?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-levenberg-marquardt/
page_id: optimization-levenberg-marquardt
date: 2026-09-11 20:00:00 +0900
categories: [Computation, Optimization]
tags: [optimization, levenberg-marquardt, gauss-newton, least-squares, thin-film, python]
description: Levenberg-Marquardt moves between Gauss-Newton and steepest descent through a damping parameter, tested on synthetic data fitted with a wrong refractive index.
math: true
---

[Post 2](/en/posts/optimization-newton-gauss-newton/) applied Newton's method and Gauss-Newton (GN) to the same thickness-fitting problem. Newton diverged outright from the initial value $d_0=1540$ nm, where the sign of the curvature is inverted ($J''(d)<0$); Gauss-Newton, approximating the Hessian by $J^TJ$, avoided that failure and converged in three steps. The $J$ in $J^TJ$ is the Jacobian of the residuals, a different symbol from the objective $J(d)$. The end of that post noted that Gauss-Newton has a theoretical limit of its own: if the discarded term $\sum_i r_i r_i''$ is not negligible — large residuals, or a Jacobian close to singular — convergence slows or fails altogether.

This post covers Levenberg-Marquardt (LM), which uses a damping parameter $\mu$ to prevent the failure in which the Gauss-Newton step grows too large and diverges. The cause may be large residuals, as in section 7 of post 2, or a Jacobian close to zero, as reproduced later in this post with a wrong refractive index. Near the solution, however, LM uses the same approximate Hessian as Gauss-Newton, so it does not cure the slow convergence of large-residual problems. LM is the method most metrology software actually standardizes on. The experiments here construct a situation in which Gauss-Newton genuinely fails, to see how LM rescues it.

## 1. How can the Gauss-Newton step be tamed?

In post 2 the Gauss-Newton step came from solving the normal equation $(J^TJ)\,\Delta d = -J^Tr$ — here a scalar, with $J^TJ = \sum_i r_i'^2$ and $J^Tr=\sum_i r_i' r_i$. The trouble starts when the denominator $J^TJ$ approaches zero. The step is $-J^Tr$ divided by $J^TJ$, so as the denominator shrinks the step grows without bound. How can it be held in check? LM's answer is to add a positive $\mu$ to the diagonal of the normal equation.

$$ (J^TJ + \mu I)\, h_{lm} = -J^Tr $$

With thickness as the only parameter, $J^TJ$ is the scalar $A=\sum_i r_i'^2$ and the step reduces to:

$$ h_{lm} = -\frac{g}{A+\mu}, \qquad g = \sum_i r_i' r_i $$

Here $g$ is the gradient of the objective, $J'(d)$, and $\mu$ is the damping parameter. Its role is clearest at the extremes. As $\mu \to 0$, $h_{lm}$ becomes exactly the Gauss-Newton step. As $\mu \to \infty$, $A$ becomes negligible and $h_{lm} \approx -g/\mu$ — steepest descent, a step of $1/\mu$ times the gradient in the opposite direction. This is post 1's gradient descent with learning rate $\alpha=1/\mu$. A single parameter serves as the handle that slides between "safe but slow" and "fast but vulnerable to curvature" — a bridge between the two methods.

The form of the step has a practical consequence. Its size, $\lvert h_{lm}\rvert = \lvert g\rvert/(A+\mu)$, never exceeds either the Gauss-Newton step $\lvert g\rvert/A$ or the steepest-descent step $\lvert g\rvert/\mu$. Even if $A$ collapses to zero, the step stops at $\lvert g\rvert/\mu$. The two behaviours part ways at $\mu\approx A$: well below it LM acts like Gauss-Newton, well above it like steepest descent.

One point is worth making explicit. With a single parameter, Gauss-Newton and steepest descent point the same way, since both follow the sign of $-g$. In this example $\mu$ therefore changes only the length of the step. With several parameters, $\mu$ also turns the direction from the Gauss-Newton direction toward the gradient.

The initial damping is set as $\mu_0 = \tau \cdot A$. Scaling by $A$ gives $\mu$ the same units and magnitude as $A$, so $\tau$ is the ratio of damping to Gauss-Newton curvature: $\tau\ll1$ starts close to Gauss-Newton, $\tau\gg1$ close to steepest descent. The user chooses $\tau$: small ($10^{-6}$) if the initial value is believed to be close to the solution, larger ($10^{-3}$, or 1) if it is not — the recommendation of Madsen and co-authors.

```python
# core of levenberg_marquardt.py (full code: _code/optimization-levenberg-marquardt/)
mu = tau * A   # A = sum(Jr**2), g = sum(Jr*r)
nu = 2.0
for _ in range(n_iter):
    h_lm = -g / (A + mu)
    d_new = d + h_lm
    F_old, F_new = objective(d, ...), objective(d_new, ...)
    L0_minus_Lh = 0.5 * h_lm * (mu * h_lm - g)
    rho = (F_old - F_new) / L0_minus_Lh if L0_minus_Lh > 0 else -1.0

    if rho > 0:                    # accept the step
        d = d_new
        Jr, r = numerical_jacobian(d, ...), residual(d, ...)
        A, g = np.sum(Jr ** 2), np.sum(Jr * r)
        mu = mu * max(1/3, 1 - (2 * rho - 1) ** 3)
        nu = 2.0
    else:                           # reject the step, become more conservative
        mu = mu * nu
        nu = 2.0 * nu
```

## 2. A large initial damping against a small one

How large should $\mu_0$ be? Trying the two extremes shows directly what $\mu$ does. Under the conditions of posts 1 and 2 ($d_0=1540$ nm, true thickness 1490 nm), two cases were run: $\tau=10^{-6}$, close to Gauss-Newton, and $\tau=10^{6}$, close to steepest descent. The second value lies deliberately far outside Madsen's recommended range ($10^{-6}$ to 1). At this point $A \approx 0.002644$, giving $\mu_0 \approx 2.644\times10^{-9}$ and $\mu_0 \approx 2644$ respectively.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig1-mu0-comparison.png" alt="LM convergence trajectories by size of mu0" width="600">
_Fig 1. LM thickness estimates for different initial damping μ0 (d0=1540nm, true thickness 1490nm). τ=10⁻⁶ (blue) reaches the noise floor in 3 steps, τ=10⁶ (red) in about 17._

With $\tau=10^{-6}$ (blue) the trajectory is effectively Gauss-Newton's, dropping to the noise floor in three steps. The blue line ends at iteration 11 because the gradient $g$ fell below the stopping threshold ($10^{-12}$). With $\tau=10^{6}$ (red) the estimate barely moves from 1540 nm for nearly ten steps. $\mu$ is large enough that the step is crushed down to $-g/\mu$. Then, from around the eleventh step, it accelerates and reaches the same point within twenty.

The moment of acceleration can be predicted. Each accepted step multiplies $\mu$ by 1/3 (section 3), so bringing $\mu$ down from $10^6A$ to the level of $A$ takes $\log_3 10^6 \approx 12.6$ steps. In the run, $\mu/A$ is 1.87 at step 13 and 0.62 at step 14, and the largest steps (−10.2, −16.1, −12.9 nm) cluster right there. The transition from "slow and rarely overshooting" to "fast Gauss-Newton" is visible in a single plot.

## 3. When should the damping shrink, and when should it grow?

Nobody steers the transition in Fig 1. A criterion decides at every step whether $\mu$ goes down or up, and that criterion is the gain ratio $\rho$.

$$ \rho = \frac{F(d) - F(d+h_{lm})}{L(0) - L(h_{lm})}, \qquad L(0)-L(h_{lm}) = \frac{1}{2}h_{lm}(\mu h_{lm} - g) $$

Here $F(d)=\frac{1}{2}\sum_i r_i^2$ is the objective $J(d)$, and $L(h)=F(d)+gh+\frac{1}{2}Ah^2$ is the local quadratic model obtained by linearizing the residuals. The denominator is the reduction predicted by the local quadratic model; the numerator is the reduction actually achieved. In other words, $\rho$ scores the forecast against the outcome. Where the weather forecast has been reliable, you trust it and travel far; when it misses, you stay close.

If $\rho>0$ — the objective actually decreased — the step is accepted and $\mu \leftarrow \mu \cdot \max\{1/3,\ 1-(2\rho-1)^3\}$. If $\rho \le 0$, the step is rejected, $\mu \leftarrow \nu\mu$ enlarges the damping, and $\nu$ is doubled, so consecutive rejections grow $\mu$ by 2, 4, 8 times.

Reading the factor makes the rule plain. For $\rho$ above about 0.94 the factor sits at its floor of 1/3; at $\rho=0.5$ it is 1 and $\mu$ is unchanged; as $\rho$ approaches 0 it approaches 2. A good forecast cuts $\mu$ to a third and lengthens the stride; a poor one raises $\mu$ even when the step is accepted. Marquardt's original rule simply multiplied and divided $\mu$ by fixed factors, which made $\mu$ flutter up and down; the rule used here, due to Nielsen, damps that oscillation.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig2-mu-adaptation.png" alt="Adaptation of mu driven by the gain ratio" width="600">
_Fig 2. The damping parameter μ as set by the gain ratio (log scale, d0=1540nm). While steps are accepted μ falls by a factor of 3 each time, a straight line; the blue run grows again after reaching the noise floor, as rejections pile up._

In both cases $\mu$ shrinks by exactly a factor of 1/3 at every accepted step — a straight line on the log scale — because the local quadratic approximation fits well enough that $\rho$ stays above about 0.94 (values above 1 included), where the factor hits its floor. The $\tau=10^{-6}$ run (blue) behaves differently once it reaches the noise floor. When the objective $J$ stops decreasing, at the floor set by the noise, subsequent steps give $\rho \le 0$ and are rejected, and $\mu$ grows again by $\times2, \times4, \times8, \dots$ — the algorithm detecting that there is nothing left to reduce and turning conservative on its own.

## 4. Does damping cost anything on a problem that already converges?

Adding damping might be expected to slow down a problem that was already going well. It does not. Under the same conditions — $d_0=1540$ nm, ten steps — LM ($\tau=10^{-6}$) was run alongside the three methods from posts 1 and 2.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig3-four-methods-comparison.png" alt="The four methods compared" width="600">
_Fig 3. Thickness estimates of the four methods (d0=1540nm, 10 steps). GN (green) and LM (purple dashed) overlap and converge in 3 steps, gradient descent (blue) reaches the same value in 5, and Newton (red) diverges._

In Fig 3 the green Gauss-Newton line and the purple dashed LM line overlap from the start, the blue gradient-descent line reaches the same floor about two steps later, and the red Newton line wanders more than 100 nm from the true thickness.

| Method | Result at $d_0=1540$ nm (10 steps) |
|---|---|
| Gradient Descent ($\alpha=300$) | 1490.12 nm, converged in 5 steps |
| Newton | no convergence, over 100 nm from the minimum |
| Gauss-Newton | 1490.12 nm, converged in 3 steps |
| Levenberg-Marquardt ($\tau=10^{-6}$) | 1490.12 nm, converged in 3 steps (effectively identical to GN) |

No specific thickness is given for the Newton row, for the reason established in post 2: the diverging trajectory passes through regions where $J''$ is near zero, amplifying rounding error, so the individual iterates are not reproducible.

With $\tau$ small enough, LM traces a trajectory indistinguishable from Gauss-Newton's. That is as it must be: $\mu_0$ several orders of magnitude below $A$ leaves the normal equation essentially unchanged. The damping does its real work in the wrong-refractive-index demonstration.

## 5. Can stronger damping keep the fit out of a neighbouring basin?

Stepping cautiously with a large $\mu$ might seem a way to avoid sliding into a neighbouring basin. It is not. Damping changes the conditioning of the normal equation, not the shape of the objective, so the local-minimum basin problem seen in post 2 survives in LM untouched. Restarting from around 1300 nm and 1690 nm, and varying $\tau$ across five orders of magnitude ($10^{-3}$ to $10^{2}$), LM converges to exactly the same neighbouring minima as GN.

| Initial value $d_0$ | Gauss-Newton | LM ($\tau=10^{-3}$) | LM ($\tau=1$) | LM ($\tau=100$) |
|---|---|---|---|---|
| 1300 nm | 1293.12 nm | 1293.12 nm | 1293.12 nm | 1293.12 nm |
| 1690 nm | 1688.55 nm | 1688.55 nm | 1688.55 nm | 1688.55 nm |

By rejecting large steps, damping usually keeps the fit within one basin and carries it safely to that basin's minimum; which basin is entered is still decided by the initial value.

## 6. What happens when the refractive index is wrong?

Let us build both failure conditions from post 2 at once: large residuals and a near-singular Jacobian. This section is a demonstration constructed on synthetic data, not a measurement. The data are those of posts 1 and 2 (true thickness 1490 nm, SiO2 index 1.46); only the fitting model is given the wrong SiO2 index, 1.02. Being so close to air (1.0), that value all but removes the reflection at the air/SiO2 interface (reflectance 0.035 → about $10^{-4}$), and with it the thin-film interference contrast. Evaluating $\sum_i r_i'^2$ at $d=1490$ nm gives $4.34\times10^{-6}$ against $2.58\times10^{-3}$ for the correct model — smaller by a factor of about 595. The Jacobian has collapsed toward zero.

Running Gauss-Newton from $d_0=1540$ nm in this state, the estimate simply cycles among 1665.5 nm, 1543.4 nm and 1416.8 nm, showing no convergence at all over twenty steps. The denominator of the normal equation, $\sum_i r_i'^2$, is tiny while the residuals stay large. The discarded term $\sum_i r_i r_i''$ is therefore about three times the denominator at this model's minimum, it dominates the approximate Hessian, and every step overshoots badly. Near the minimum Gauss-Newton multiplies the error by about −3 at every step: each miss throws the estimate three times as far to the other side. LM ($\tau=10^{-3}$) under identical conditions converges steadily to the local minimum of this (wrong) model nearest the start, $d\approx1486.6$ nm (the lowest point between 1400 and 1700 nm).

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig4-gn-failure-lm-rescue.png" alt="GN cycles, LM converges" width="600">
_Fig 4. Thickness estimates with the SiO2 index wrongly set to 1.02 (synthetic data, d0=1540nm). GN (red) cycles among three values; LM (purple, τ=10⁻³) rejects its first four steps and then converges to this model's local minimum at 1486.6nm (dashed)._

Follow the purple line in Fig 4. The first four steps are all rejected ($\rho\approx-0.7$), while $\mu$ grows from 0.001 to 0.064 times $A$ and the estimate stays at 1540 nm. At the fifth step $\mu$ is about equal to $A$, and a step of −61 nm, half the Gauss-Newton step of −124 nm, is accepted. From then on $\rho$ stays between 0.35 and 0.44, the update factor is close to 1, and $\mu$ settles at about $2A$. The step becomes $A/(A+2A)=1/3$ of the Gauss-Newton step. The actual steps, +11.3 → −4.4 → +1.4 → −0.46 → +0.15 nm, alternate in sign and shrink by about a third each time. Where Gauss-Newton multiplied the error by −3, LM raises $\mu$ on its own until the error is multiplied by −1/3. Adding $\mu$ to the denominator keeps the normal equation away from singularity and suppresses the step size automatically. The failure condition Madsen and co-authors warned about in theory is thus confirmed on an example built to satisfy exactly that condition.

What LM prevents, however, is divergence, not slow convergence. An error that shrinks by a third per step is first-order (linear) convergence. With the correct model three steps brought the estimate within 0.01 nm; here it takes 13, including the five rejected steps. Near the solution LM uses the same approximate Hessian as Gauss-Newton, so the slow convergence of large-residual problems remains.

Nor does convergence make the answer right. The 1486.6 nm that LM finds is the minimum of a model with the wrong refractive index. There the RMS residual is 0.142, about 35 times the noise put into the data (standard deviation 0.004). At the correct model's solution (1490.12 nm) the RMS is 0.0040, equal to the noise. In terms of the objective that is 3.03 against 0.0024, and the latter matches the value expected from noise alone, $\frac{1}{2}\times300\times0.004^2=0.0024$. After convergence, then, the RMS residual should be compared with the noise level. If it is tens of times larger, suspect the model — refractive index, dispersion, layer structure — rather than the algorithm. LM fixed the convergence, not the model.

## Summary

From post 1 (gradient descent) through this one, four methods have been applied to the same thickness-fitting problem.

| Method | Step | Strength | Weakness |
|---|---|---|---|
| Gradient Descent | $-\alpha J'(d)$ | simplest to implement | the step size must be chosen by hand |
| Newton | $-J'(d)/J''(d)$ | quadratic convergence near the solution | diverges where the curvature sign flips in a non-convex region |
| Gauss-Newton | $-g/A$ (the $J^TJ$ approximation) | $J^TJ$ is always semidefinite, so the curvature problem disappears | slow or failing convergence with large residuals; overly long steps with a near-singular Jacobian |
| Levenberg-Marquardt | $-g/(A+\mu)$ | $\mu$ interpolates GN and steepest descent, preventing divergence from overly long steps, whether caused by large residuals or a near-zero Jacobian | convergence stays slow with large residuals; the basin problem remains; $\tau$ must be chosen |

Why Levenberg-Marquardt is the de facto standard in metrology software condenses into one row of that table: it keeps almost all of Gauss-Newton's speed while damping steps in as a safeguard at precisely the point where the Gauss-Newton step grows too long. What none of the four methods settles is which local minimum the fit ends up in — a problem all of them share.

In practice: start with $\tau=10^{-3}$, and lower it to $10^{-6}$ if the initial value is trustworthy. After convergence, compare the RMS residual with the noise level and suspect the model if it is tens of times larger. If the answer changes whenever the initial value changes, the problem is the basin, not the local method. [Post 4](/en/posts/optimization-global-heuristics/) takes on that basin problem directly with simulated annealing and with Basin-Hopping, which repeats local optimization while hopping between basins.

## References

- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

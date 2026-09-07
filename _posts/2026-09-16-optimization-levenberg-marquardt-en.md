---
title: "Optimization 3 — The Levenberg-Marquardt Method"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-levenberg-marquardt/
page_id: optimization-levenberg-marquardt
date: 2026-09-16 20:00:00 +0900
categories: [Computation, Optimization]
tags: [optimization, levenberg-marquardt, gauss-newton, least-squares, thin-film, python]
description: Levenberg-Marquardt moves smoothly between Gauss-Newton and steepest descent through a damping parameter, tested here against a case constructed to make Gauss-Newton fail.
math: true
---

[Post 2](/en/posts/optimization-newton-gauss-newton/) applied Newton's method and Gauss-Newton to the same thickness-fitting problem. Newton diverged outright from the initial value $d_0=1540$ nm, where the sign of the curvature is inverted ($J''(d)<0$); Gauss-Newton, approximating the Hessian by $J^TJ$, was free of that failure and converged in three steps. The end of that post noted that Gauss-Newton has a theoretical limit of its own: if the discarded term $\sum_i r_i r_i''$ is not negligible — large residuals, or a Jacobian close to singular — convergence slows or fails altogether.

This post covers Levenberg-Marquardt (LM), which compensates for that limit with a damping parameter $\mu$. It is the method most metrology software actually standardizes on, and the experiments here construct a situation in which Gauss-Newton genuinely fails, to see how LM rescues it.

## 1. Levenberg-Marquardt — a bridge between GN and steepest descent

In post 2 the Gauss-Newton step came from solving the normal equation $(J^TJ)\,\Delta d = -J^Tr$ — here a scalar, with $J^TJ = \sum_i r_i'^2$ and $J^Tr=\sum_i r_i' r_i$. Levenberg-Marquardt adds $\mu$ to the diagonal of that equation.

$$ (J^TJ + \mu I)\, h_{lm} = -J^Tr $$

With thickness as the only parameter, $J^TJ$ is the scalar $A=\sum_i r_i'^2$ and the step reduces to:

$$ h_{lm} = -\frac{g}{A+\mu}, \qquad g = \sum_i r_i' r_i $$

The role of $\mu$ is clearest at the extremes. As $\mu \to 0$, $h_{lm}$ becomes exactly the Gauss-Newton step. As $\mu \to \infty$, $A$ becomes negligible and $h_{lm} \approx -g/\mu$ — steepest descent, moving opposite the gradient by $1/\mu$. A single parameter serves as the handle that slides between "safe but slow" and "fast but vulnerable to curvature."

The initial damping is set as $\mu_0 = \tau \cdot A$. The scale factor $\tau$ is chosen by the user: small ($10^{-6}$) if the initial value is believed to be close to the solution, larger ($10^{-3}$, or 1) if it is not — the recommendation of Madsen and co-authors.

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

## 2. The experiment — a large $\mu_0$ against a small one

Under the conditions of posts 1 and 2 ($d_0=1540$ nm, true thickness 1490 nm), two cases were run: $\tau=10^{-6}$, close to Gauss-Newton, and $\tau=10^{6}$, close to steepest descent. At this point $A \approx 0.002644$, giving $\mu_0 \approx 2.644\times10^{-9}$ and $\mu_0 \approx 2644$ respectively.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig1-mu0-comparison.png" alt="LM convergence trajectories by size of mu0" width="600">
_Fig 1. LM convergence trajectories by the size of μ0 (d0=1540nm)_

With $\tau=10^{-6}$ (blue) the trajectory is effectively Gauss-Newton's, dropping to the noise floor in three steps. With $\tau=10^{6}$ (red) the estimate barely moves from 1540 nm for nearly ten steps — $\mu$ is large enough that the step is crushed down to $-g/\mu$ — and then accelerates from around the eleventh step, reaching the same point within twenty. The transition from "slow but never overshooting" to "fast Gauss-Newton" is visible in a single plot.

## 3. How the gain ratio adjusts $\mu$

That transition happens automatically, driven by a gain ratio $\rho$ computed at every step.

$$ \rho = \frac{F(d) - F(d+h_{lm})}{L(0) - L(h_{lm})}, \qquad L(0)-L(h_{lm}) = \frac{1}{2}h_{lm}(\mu h_{lm} - g) $$

The denominator is the reduction predicted by the local quadratic model; the numerator is the reduction actually achieved. If $\rho>0$, the model is holding up: the step is accepted and $\mu$ is reduced by $\mu \leftarrow \mu \cdot \max\{1/3,\ 1-(2\rho-1)^3\}$. If $\rho \le 0$, the step is rejected, $\mu \leftarrow \nu\mu$ enlarges the damping, and $\nu$ is doubled. Marquardt's original rule simply multiplied and divided, which produced flutter; the rule used here, due to Nielsen, damps that oscillation.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig2-mu-adaptation.png" alt="Adaptation of mu driven by the gain ratio" width="600">
_Fig 2. Adaptation of μ driven by the gain ratio_

In both cases $\mu$ shrinks by exactly a factor of 1/3 at every accepted step — a straight line on the log scale — because the local quadratic approximation fits well enough that $\rho$ stays close to 1. What happens after the $\tau=10^{-6}$ run (blue) reaches the noise floor is the interesting part. Once the objective $J$ stops decreasing, at the floor set by the noise, subsequent steps give $\rho \le 0$ and are rejected, and $\mu$ grows again by $\times2, \times4, \times8, \dots$ — the algorithm detecting that there is nothing left to reduce and turning conservative on its own.

## 4. The four methods side by side

Under the same conditions — $d_0=1540$ nm, ten steps — Levenberg-Marquardt ($\tau=10^{-6}$) was run alongside the three methods from posts 1 and 2.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig3-four-methods-comparison.png" alt="The four methods compared" width="600">
_Fig 3. The four methods compared (d0=1540nm)_

| Method | Result at $d_0=1540$ nm (10 steps) |
|---|---|
| Gradient Descent (alpha=300) | 1490.12 nm, converged in 5 steps |
| Newton | no convergence, over 100 nm from the minimum |
| Gauss-Newton | 1490.12 nm, converged in 3 steps |
| Levenberg-Marquardt ($\tau=10^{-6}$) | 1490.12 nm, converged in 3 steps (effectively identical to GN) |

No specific thickness is given for the Newton row, for the reason established in post 2: the diverging trajectory passes through regions where $J''$ is near zero, amplifying rounding error, so the individual iterates are not reproducible.

With $\tau$ small enough, LM traces a trajectory indistinguishable from Gauss-Newton's — as it must, since $\mu_0$ several orders of magnitude below $A$ leaves the normal equation essentially unchanged. Where LM earns its place comes next.

## 5. LM is not a cure-all — the basin problem remains

Damping changes the conditioning of the normal equation, not the shape of the objective. The local-minimum basin problem seen in post 2 therefore survives in LM untouched. Restarting from around 1300 nm and 1690 nm, and varying $\tau$ across three orders of magnitude, LM converges to exactly the same neighbouring minima as GN.

| Initial value $d_0$ | Gauss-Newton | LM ($\tau=10^{-3}$) | LM ($\tau=1$) | LM ($\tau=100$) |
|---|---|---|---|---|
| 1300 nm | 1293.12 nm | 1293.12 nm | 1293.12 nm | 1293.12 nm |
| 1690 nm | 1688.55 nm | 1688.55 nm | 1688.55 nm | 1688.55 nm |

Damping guarantees safe arrival at the minimum *within a basin*; which basin is entered is still decided by the initial value.

## 6. A real failure case — assuming the wrong refractive index

The theoretical limit of GN flagged in post 2 — failure when the residuals are large or the Jacobian near-singular — was reproduced deliberately. Suppose the fitting model assumes a refractive index of 1.02 for SiO2 instead of the actual 1.46. Being so close to air (1.0), that value all but erases the interference contrast at the SiO2/Si interface. Evaluating $\sum_i r_i'^2$ at $d=1490$ nm gives $4.34\times10^{-6}$ against $2.58\times10^{-3}$ for the correct model — smaller by a factor of about 595. The Jacobian has collapsed toward zero.

Running Gauss-Newton from $d_0=1540$ nm in this state, the estimate simply cycles among 1665.5 nm, 1543.4 nm and 1416.8 nm, showing no convergence at all over twenty steps. The denominator of the normal equation, $\sum_i r_i'^2$, is so small that every step overshoots badly. Levenberg-Marquardt under identical conditions converges steadily to the true minimum of this (wrong) model, $d\approx1486.6$ nm.

<img src="/assets/img/posts/optimization-levenberg-marquardt/en/fig4-gn-failure-lm-rescue.png" alt="GN fails, LM succeeds" width="600">
_Fig 4. With a wrongly assumed refractive index: GN fails, LM rescues_

Adding $\mu$ to the denominator keeps the normal equation away from singularity and suppresses the step size automatically. The failure condition Madsen warned about in theory has been confirmed directly, on an example built to satisfy exactly that condition.

## Summary

From post 1 (gradient descent) through this one, four methods have been applied to the same thickness-fitting problem.

| Method | Step | Strength | Weakness |
|---|---|---|---|
| Gradient Descent | $-\alpha J'(d)$ | simplest to implement | the step size must be chosen by hand |
| Newton | $-J'(d)/J''(d)$ | quadratic convergence near the solution | diverges where the curvature sign flips in a non-convex region |
| Gauss-Newton | $-g/A$ (the $J^TJ$ approximation) | $J^TJ$ is always semidefinite, so the curvature problem disappears | fails on large residuals or a near-singular Jacobian |
| Levenberg-Marquardt | $-g/(A+\mu)$ | $\mu$ interpolates GN and steepest descent, rescuing GN's failure mode | the basin problem remains; $\mu_0$ needs tuning |

Why Levenberg-Marquardt is the de facto standard in metrology software condenses into one row of that table: it keeps almost all of Gauss-Newton's speed while damping steps in as a safeguard at precisely the point where Gauss-Newton collapses. What none of the four methods settles is which local minimum the fit ends up in — a problem all of them share.

## References

- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

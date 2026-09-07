---
title: "Optimization 2 — Newton's Method and Gauss-Newton"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-newton-gauss-newton/
page_id: optimization-newton-gauss-newton
date: 2026-09-10 20:00:00 +0900
categories: [Computation, Optimization]
tags: [optimization, newton-method, gauss-newton, least-squares, thin-film, python]
description: Solving the same thin-film thickness fit with Newton's method and with Gauss-Newton, and seeing why Gauss-Newton is the standard in metrology practice.
math: true
---

[Post 1](/en/posts/optimization-gradient-descent/) solved the problem of recovering the thickness of a single SiO2/Si film from its reflectance spectrum using gradient descent. With a well-chosen step size `alpha` (300, say) it converged without trouble, but only slightly larger (1500) and it bounced across the valley and diverged. This post starts with Newton's method, which promises to settle that step-size dilemma on its own, and goes on to Gauss-Newton, which is used far more widely in metrology.

To state the conclusion first: the expectation that "Newton's method uses second derivatives, so it is theoretically faster and therefore always better" breaks down on this problem. Examining why shows precisely why metrology software mostly standardizes on Gauss-Newton — or its variant, Levenberg-Marquardt — rather than Newton.

## 1. Newton's method — letting the second derivative set the step

In section 4 of post 1, gradient descent came from the first-order Taylor expansion $J(d+\Delta d) \approx J(d) + J'(d)\Delta d$, taking $\Delta d$ opposite to the gradient. The open question was *how far* to move: the step size $\alpha$ had to be chosen by hand.

Newton's method goes one order further and uses the second-order expansion.

$$ J(d+\Delta d) \approx J(d) + J'(d)\Delta d + \frac{1}{2}J''(d)\Delta d^2 $$

Treating the right-hand side as a function of $\Delta d$ and minimizing it — differentiating with respect to $\Delta d$ and setting the result to zero — gives $\Delta d$ directly.

$$ J'(d) + J''(d)\Delta d = 0 \quad\Longrightarrow\quad \Delta d = -\frac{J'(d)}{J''(d)} $$

$$ d_{k+1} = d_k - \frac{J'(d_k)}{J''(d_k)} $$

Compared with gradient descent, $1/J''(d_k)$ has taken the place of $\alpha$. Where the local curvature is large — a narrow, steep valley — the step is automatically small; where the curvature is gentle, the step is large. Figure 1 shows the parabola that locally approximates the actual $J(d)$ at a thickness of $d_k=1500$ nm.

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig1-quadratic-approximation.png" alt="Local quadratic approximation of the objective and the Newton step" width="600">
_Fig 1. Local quadratic approximation of the objective and the Newton step_

At this point the parabola tracks the true curve closely enough to overlap it, and a single step ($d=1500 \to d\approx1489.8$) reaches the neighbourhood of the minimum. Both $J'(d)$ and $J''(d)$ were obtained by central difference, as in post 1.

```python
# core of newton.py (full code: _code/optimization-newton-gauss-newton/)
def numerical_grad_hess(thickness_nm, wavelength_nm, measured_R, h=1e-3):
    J0 = objective(thickness_nm, wavelength_nm, measured_R)
    J_plus = objective(thickness_nm + h, wavelength_nm, measured_R)
    J_minus = objective(thickness_nm - h, wavelength_nm, measured_R)
    grad = (J_plus - J_minus) / (2 * h)
    hess = (J_plus - 2 * J0 + J_minus) / (h ** 2)
    return grad, hess


def newton(d0, wavelength_nm, measured_R, n_iter=100, h=1e-3):
    d = d0
    d_hist, J_hist = [d], [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        grad, hess = numerical_grad_hess(d, wavelength_nm, measured_R, h)
        d = d - grad / hess
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)
```

In theory Newton's method converges quadratically near the solution — the error is squared at every step — whereas gradient descent manages only linear convergence. Applying Newton's method at the initial value of post 1, $d_0=1540$ nm, should therefore converge far faster than gradient descent did.

## 2. The experiment — running from the same initial value

Under exactly the conditions of post 1 (seed=0, true thickness 1490 nm, the same noise, $d_0=1540$ nm), gradient descent (alpha=300), Newton, and the Gauss-Newton of section 5 were run side by side.

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig2-thickness-vs-iteration.png" alt="Thickness estimate versus iteration" width="600">
_Fig 2. Thickness estimate vs iteration (d0=1540nm)_

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig3-objective-vs-iteration.png" alt="Objective function J versus iteration" width="600">
_Fig 3. Objective J vs iteration (d0=1540nm, log scale)_

The outcome was not the expected one. Gradient descent (blue) converges smoothly, as in post 1. Newton (red), however, jumps to 1876 nm on the very first step and then wanders between 1824 and 1921 nm without ever settling. For all that quadratic convergence is meant to be faster, at this initial value Newton's method simply does not work.

One caveat belongs here. The individual iterates along this diverging trajectory are not numbers worth quoting. The trajectory passes through regions where $J''$ is close to zero, which amplifies rounding error. Perturbing the initial value by as little as $10^{-12}$ nm moves the tenth iterate by hundreds of nanometres.

```text
sensitivity of the diverging trajectory to the initial value (10th-step thickness):
  d0=1540.000000000000 -> 1636.55 nm
  d0=1540.000000000001 -> 1855.29 nm
  d0=1540.000000000100 -> 1854.98 nm
```

What reproduces is the fact that the method fails to converge and wanders hundreds of nanometres from the minimum — not any particular iterate.

## 3. Why it failed — where the quadratic approximation collapses

The objective $J(d)$ plotted in section 5 of post 1 (Fig. 2 there) makes the reason visible. The reflectance of a thin film oscillates with thickness as interference fringes, and $J(d)$ was correspondingly non-convex: a global minimum with a row of local minima around it. For the quadratic approximation Newton relies on to hold, the curvature at that point — $J''(d)$ — must be positive, since the parabola has to open upward for the step to point toward a minimum. Evaluating $J''(d)$ at two points gives:

| Point | $J''(d)$ | Meaning |
|---|---|---|
| $d=1500$ | $+0.002487$ | convex — a well-behaved Newton step |
| $d=1540$ | $-0.000233$ | concave — the step is thrown the other way |

At $d=1540$ the sign of the curvature is itself inverted. The point sits near a local "hill" thrown up by the interference fringes, and the Newton step $-J'/J''$ moves in the direction of *increasing* objective whenever $J''$ is negative. Figure 4 shows the same method meeting entirely different fates from starting points only 40 nm apart.

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig4-newton-success-vs-failure.png" alt="Two faces of Newton's method" width="600">
_Fig 4. Two faces of Newton's method — success and failure according to the initial value_

From $d_0=1500$ (green) the curvature is well behaved and the method drops to the noise floor in two steps, textbook quadratic convergence. From $d_0=1540$ (red) it is still further from the minimum after ten steps than it started. This is not an implementation bug but a structural weakness of Newton's method. The guarantee that a Newton step is a descent direction holds only where the Hessian at that point is positive definite. For an objective as folded as this one, the fate of the method depends on which "valley" the initial value happens to sit in.

## 4. The limits of Newton's method

Here there was a single parameter, the thickness, so $J''(d)$ was a scalar. Real metrology often fits several parameters at once — thickness, refractive index, extinction coefficient. With $n$ parameters, $J''$ becomes an $n \times n$ Hessian that must be formed and inverted at every step. Beyond the cost, which grows quickly with $n$, the risk seen in this experiment remains: pass through a point where the Hessian is not positive definite and the step is thrown in a direction unrelated to any minimum. That is why Newton's method is not used bare in practice but always paired with a safeguard such as a line search or a trust region.

## 5. Gauss-Newton — approximating the Hessian for safety

Nonlinear least squares has a somewhat special structure. The objective is fixed as a sum of squared residuals $r_i(d) = R_{model}(d,\lambda_i) - R_{meas}(\lambda_i)$.

$$ J(d) = \frac{1}{2}\sum_i r_i(d)^2 $$

Differentiating this form gives the first and second derivatives directly.

$$ J'(d) = \sum_i r_i(d)\, r_i'(d), \qquad J''(d) = \sum_i \Big[ r_i'(d)^2 + r_i(d)\, r_i''(d) \Big] $$

Newton uses all of $J''(d)$. Gauss-Newton discards the second term $\sum_i r_i(d) r_i''(d)$ entirely and keeps only the first as an approximate Hessian.

$$ J''(d) \;\approx\; \sum_i r_i'(d)^2 $$

(In the general multi-parameter notation this approximate Hessian is written $J^TJ$, with $J$ the Jacobian of the residual vector.) The discarded term is a product of the residual $r_i$ and the residual's curvature $r_i''$. If the model fits the data well and the residuals are small — or if the model is locally near-linear, making $r_i''$ small — the term was small to begin with and costs little to drop. With this approximation the step becomes:

$$ d_{k+1} = d_k - \frac{\sum_i r_i'(d_k)\, r_i(d_k)}{\sum_i r_i'(d_k)^2} $$

In the general multi-parameter case this is the normal equation $(J^TJ)\,\Delta\mathbf{d} = -J^T\mathbf{r}$. In code:

```python
# core of gauss_newton.py (full code: _code/optimization-newton-gauss-newton/)
def gauss_newton(d0, wavelength_nm, measured_R, n_iter=100, h=1e-3):
    d = d0
    d_hist, J_hist = [d], [objective(d, wavelength_nm, measured_R)]
    for _ in range(n_iter):
        r = residual(d, wavelength_nm, measured_R)
        Jr = numerical_jacobian(d, wavelength_nm, measured_R, h)
        step = -np.sum(Jr * r) / np.sum(Jr ** 2)  # (J^T J) h = -J^T r
        d = d + step
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)
```

This approximate Hessian $J^TJ$ — here $\sum_i r_i'^2$ — has a property that Newton's true Hessian lacks. Being a sum of squares, $J^TJ$ is always positive semidefinite, and positive definite whenever the Jacobian has full rank. The Gauss-Newton step therefore has no point at which the sign of the curvature can flip, which frees it structurally from the very cause of Newton's divergence in section 3.

## 6. The experiment — Gauss-Newton at the same failure point

The green curve in Figures 2 and 3 is Gauss-Newton. From $d_0=1540$ nm, the exact starting point at which Newton diverged, it reaches the noise floor by the third step — faster even than gradient descent (alpha=300), which took five under the same conditions.

| Method | Result at $d_0=1540$ nm |
|---|---|
| Gradient Descent (alpha=300) | 1490.12 nm, converged in 5 steps |
| Newton | no convergence in 10 steps, over 100 nm from the minimum |
| Gauss-Newton | 1490.12 nm, converged in 3 steps |

The whole difference comes from one fact: whatever the sign of the true curvature, $\sum_i r_i'^2$ is never negative. That is exactly why metrology software standardizes on Gauss-Newton — or on Levenberg-Marquardt, the subject of the next post — rather than plain Newton. Saving the cost of second derivatives matters, but this structural stability matters more.

## 7. Gauss-Newton is not a cure-all either

Gauss-Newton is not safe in every situation. The discarded term $\sum_i r_i r_i''$ is not always negligible. Two failure conditions are known.

First, if the residuals themselves are large — the model fits the data poorly — the discarded term grows, and far from converging quadratically the method may manage only linear convergence, or diverge. On a zero-residual problem, where the residual vanishes at the solution, Gauss-Newton recovers the same quadratic convergence as Newton.

Second, if the starting point lies in the basin of a different local minimum, Gauss-Newton converges to that minimum, as it must. Restarting from around 1300 nm and around 1690 nm gives 1293 nm and 1689 nm respectively — precisely the neighbouring local minima visible in Figure 2 of post 1. The stability of Gauss-Newton does not mean "it finds the global minimum"; it means something closer to "once it has a direction, it goes there without wandering."

## Summary

Newton's method promises the fastest convergence in theory, but the promise holds only where the objective is locally convex; on an objective folded by interference fringes, as here, it can fail outright depending on the initial value. Gauss-Newton exploits the structure of nonlinear least squares to approximate the Hessian by $J^TJ$, removing that failure mode altogether — at the price of slower convergence when the residuals are large.

The next post covers Levenberg-Marquardt, which compensates for this weakness of Gauss-Newton — instability at large residuals — with a damping parameter. It is the de facto standard in metrology, and its structure moves smoothly between gradient descent and Gauss-Newton, taking the strengths of both.

## References

- D. Bindel, "Nonlinear Least Squares (Newton and Gauss-Newton)," Numerical Analysis lecture notes, Cornell University, 2023. [cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf](https://www.cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf)
- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

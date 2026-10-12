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

[Post 1](/en/posts/optimization-gradient-descent/) solved the problem of recovering the thickness of a single SiO2/Si film from its reflectance spectrum using gradient descent. With a well-chosen step size `alpha` (300, say) it converged without trouble, but at five times that (1500) it bounced across the valley and oscillated without converging. This post starts with Newton's method, which promises to settle that step-size dilemma on its own, and goes on to Gauss-Newton, which is used far more widely in metrology.

To state the conclusion first: the expectation that "Newton's method uses second derivatives, so it is theoretically faster and therefore always better" breaks down on this problem. The reason also explains why metrology software mostly standardizes on Gauss-Newton — or its variant, Levenberg-Marquardt (LM) — rather than Newton.

## 1. How does Newton's method choose its own step size?

In section 4 of post 1, gradient descent came from the first-order Taylor expansion $J(d+\Delta d) \approx J(d) + J'(d)\Delta d$, taking $\Delta d$ opposite to the gradient. The open question was *how far* to move: the step size $\alpha$ had to be chosen by hand. If we look not only at the slope of the curve but also at how it bends, could the method pick that value itself? That is the idea behind Newton's method.

Newton's method goes one order further and uses the second-order expansion. It approximates the curve by a parabola at the current point and lets that parabola set the step.

$$ J(d+\Delta d) \approx J(d) + J'(d)\Delta d + \frac{1}{2}J''(d)\Delta d^2 $$

Here $J'$ and $J''$ are the first and second derivatives of $J$ with respect to thickness, and $\Delta d$ is the distance to move on this step. Treating the right-hand side as a function of $\Delta d$ and minimizing it — differentiating with respect to $\Delta d$ and setting the result to zero — gives $\Delta d$ directly.

$$ J'(d) + J''(d)\Delta d = 0 \quad\Longrightarrow\quad \Delta d = -\frac{J'(d)}{J''(d)} $$

$$ d_{k+1} = d_k - \frac{J'(d_k)}{J''(d_k)} $$

In other words, a Newton step jumps straight to the vertex of the parabola fitted at the current point. Compared with gradient descent, $1/J''(d_k)$ has taken the place of $\alpha$. Where the local curvature is large — a narrow, steep valley — the step is automatically small; where the curvature is gentle, the step is large. Figure 1 shows the parabola that locally approximates the actual $J(d)$ at a thickness of $d_k=1500$ nm.

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig1-quadratic-approximation.png" alt="Local quadratic approximation of the objective and the Newton step" width="600">
_Fig 1. Quadratic approximation of the objective fitted at d=1500 nm (red dashed) and the Newton step. The parabola nearly overlaps the actual J(d) (blue), and a single step lands at 1489.8 nm, next to the minimum._

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

In theory Newton's method converges quadratically near the solution: each error is proportional to the square of the previous one. This problem shows it too. Starting from $d_0=1500$ nm, the distance to the minimum (1490.12 nm) shrinks 9.88 → 0.288 → 0.0001 nm, after which it hits the precision limit of the numerical derivatives. The closer the iterate gets, the larger the factor by which each step cuts the error. Gradient descent, by contrast, manages only linear convergence. Applying Newton's method at the initial value of post 1, $d_0=1540$ nm, should therefore converge far faster than gradient descent did.

## 2. Running Newton's method from the same initial value

Under exactly the conditions of post 1 (seed=0, true thickness 1490 nm, the same noise, $d_0=1540$ nm), gradient descent ($\alpha=300$), Newton, and the Gauss-Newton of section 5 were run side by side.

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig2-thickness-vs-iteration.png" alt="Thickness estimate versus iteration" width="600">
_Fig 2. Thickness estimate versus iteration for the three methods from an initial value of 1540 nm. Only Newton (red) fails to converge, circling hundreds of nanometres away from the minimum._

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig3-objective-vs-iteration.png" alt="Objective function J versus iteration" width="600">
_Fig 3. The objective J for the same runs (log scale). Gradient descent reaches the noise floor in five steps and Gauss-Newton in three, while Newton stays near its starting value._

The outcome was not the expected one. Gradient descent (blue) converges smoothly, as in post 1. Newton (red), however, jumps to 1876 nm on the very first step, wanders between 1824 and 1921 nm for steps 1 to 9, and on the tenth leaves even that range, dropping to 1637 nm. Either way it is more than 100 nm from the minimum (1490 nm). These individual values are sensitive to rounding error and are not worth quoting; section 3 closes with the reason.

The objective (Fig. 3) makes the difference sharper. Gradient descent reaches the noise floor ($J\approx0.0024$) in five steps — the value below which no thickness can go, because of the measurement noise. Newton's $J$, on the other hand, stays between 2.4 and 3.1 for all ten steps and never gets away from its starting value (2.57). The green curve (Gauss-Newton) is taken up in section 6. For all that quadratic convergence is meant to be faster, at this initial value Newton's method simply does not work.

## 3. Why did Newton's method fail?

The objective $J(d)$ plotted in section 2 of post 1 (Fig. 2 there) makes the reason visible: there are points where the quadratic approximation collapses. The reflectance of a thin film oscillates with thickness as interference fringes, and $J(d)$ was correspondingly non-convex: a global minimum with a row of local minima around it. A quadratic model can be built at any point. But it has a bottom only if the curvature there — $J''(d)$ — is positive, since the parabola has to open upward for the step to point toward a minimum. Evaluating $J''(d)$ at two points gives:

| Point | $J''(d)$ | Meaning |
|---|---|---|
| $d=1500$ | $+0.002487$ | convex, a well-behaved Newton step |
| $d=1540$ | $-0.000233$ | concave, the step is thrown the other way |

At $d=1540$ the sign of the curvature is itself inverted. The curvature changes sign at about 1537 nm, and the top of the local "hill" thrown up by the interference fringes is at 1589 nm, so 1540 nm sits on the slope just starting up that hill. The Newton step $-J'/J''$ moves in the direction of *increasing* objective whenever $J''$ is negative. Figure 4 shows the same method meeting entirely different outcomes from starting points only 40 nm apart.

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig4-newton-success-vs-failure.png" alt="Two faces of Newton's method" width="600">
_Fig 4. Two faces of Newton's method. From 1500 nm (green) it reaches the noise floor in two steps; from 1540 nm (red) it is still near its starting objective after ten._

From $d_0=1500$ (green) the curvature is well behaved and the method drops to the noise floor in two steps, textbook quadratic convergence. From $d_0=1540$ (red) it is still further from the minimum after ten steps than it started.

Drawing the failed first step the way Figure 1 was drawn shows what happened (Fig. 5). The parabola fitted at $d=1540$ opens downward. With no bottom to find, the Newton step, which looks for the point of zero slope, goes to the top of the parabola. That top lies 336 nm away at 1876 nm, and the first step lands there. The parabola predicts $J=15.7$ at that point; the actual $J$ is 2.83, higher than at the start (2.57).

<img src="/assets/img/posts/optimization-newton-gauss-newton/en/fig5-concave-quadratic-approximation.png" alt="Quadratic approximation at a point of negative curvature" width="600">
_Fig 5. Quadratic approximation fitted where the curvature is negative (d=1540 nm). The parabola (red dashed) opens downward, so the Newton step goes to its top at 1876 nm, where the parabola no longer resembles the actual J(d) (blue)._

So Newton's method does not care whether the vertex of its parabola is a bottom or a top. This is not an implementation bug but a structural weakness of Newton's method. The guarantee that a Newton step is a descent direction holds only where the Hessian at that point is positive definite. For an objective as folded as this one, the fate of the method depends on which "valley" the initial value happens to sit in.

One caveat belongs here. The individual iterates along the diverging trajectory in Figures 2 and 3 are not numbers worth quoting. After the first step the trajectory passes through regions where $J''$ is close to zero, which strongly amplifies the rounding error of the central-difference $J''$ (about $10^{-9}$). Perturbing the initial value by as little as $10^{-12}$ nm changes that rounding pattern and moves the tenth iterate by more than 200 nm.

```text
sensitivity of the diverging trajectory to the initial value (10th-step thickness):
  d0=1540.000000000000 -> 1636.55 nm
  d0=1540.000000000001 -> 1855.29 nm
  d0=1540.000000000100 -> 1854.98 nm
```

What does reproduce is that the method never converges and wanders hundreds of nanometres from the minimum — not any particular iterate.

## 4. The limits of Newton's method

Here there was a single parameter, the thickness, so $J''(d)$ was a scalar. What gets harder when there are several? Real metrology often fits several parameters at once — thickness, refractive index, extinction coefficient. With $n$ parameters, $J''$ becomes an $n \times n$ Hessian that must be formed and inverted at every step. Even though it is symmetric, it holds $n(n+1)/2$ distinct second derivatives: with two thicknesses and three dispersion coefficients, five parameters, that is 15 numerical second derivatives per step. The Gauss-Newton method of section 5 needs only first derivatives of the residuals.

Cost is not the only problem. The negative-curvature risk seen in section 3 remains. That is why Newton's method is not used bare in practice but always paired with a safeguard. A line search keeps the Newton direction and shrinks the step until $J$ actually decreases. A trust region sets a radius within which the quadratic model is believed and moves only inside it. The LM method of the next post belongs to the trust-region family.

## 5. How does Gauss-Newton approximate the Hessian?

Is there a way to avoid computing the full Hessian and dodge negative curvature at the same time? There is, if we use the fact that the objective is a sum of squares. Nonlinear least squares has a somewhat special structure. The objective is fixed as a sum of squared residuals $r_i(d) = R_{model}(d,\lambda_i) - R_{meas}(\lambda_i)$.

$$ J(d) = \frac{1}{2}\sum_i r_i(d)^2 $$

Differentiating this form gives the first and second derivatives directly.

$$ J'(d) = \sum_i r_i(d)\, r_i'(d), \qquad J''(d) = \sum_i \Big[ r_i'(d)^2 + r_i(d)\, r_i''(d) \Big] $$

Here $r_i'(d)$ is how sensitive the reflectance at wavelength $\lambda_i$ is to thickness — its sensitivity — and $r_i''(d)$ is how that sensitivity changes with thickness. Newton uses all of $J''(d)$. Gauss-Newton discards the second term $\sum_i r_i(d) r_i''(d)$ entirely and keeps only the first as an approximate Hessian.

$$ J''(d) \;\approx\; \sum_i r_i'(d)^2 $$

(In the general multi-parameter notation this approximate Hessian is written $J^TJ$, with $J$ the Jacobian of the residual vector. The letter clashes with the objective $J(d)$, but we follow convention: in $J^TJ$ and $J^T\mathbf{r}$, $J$ is the Jacobian; in $J(d)$, $J'$ and $J''$ it is the objective.) The discarded term is a product of the residual $r_i$ and the residual's curvature $r_i''$. If the model fits the data well and the residuals are small — or if the model is locally near-linear, making $r_i''$ small — the term was small to begin with and costs little to drop.

The same approximation can be seen from another side: at every step, straighten each residual into a line at the current point.

$$ r_i(d+\Delta d) \approx r_i(d) + r_i'(d)\,\Delta d $$

Choosing the $\Delta d$ that minimizes the sum of squares of these lines, $\frac{1}{2}\sum_i \big[ r_i(d) + r_i'(d)\,\Delta d \big]^2$, gives exactly the step formula below. Gauss-Newton, in other words, solves one straight-line fit (a linear least-squares problem) per step. It is like replacing a winding mountain path with a flat map of where you stand, choosing where to go, taking one step, and redrawing the map. The sum of squares of a linear model is quadratic in $\Delta d$ with coefficient $\frac{1}{2}\sum_i r_i'^2 \ge 0$ on $\Delta d^2$, so it is always a bowl that opens upward. The good property of the approximate Hessian discussed below comes from here. The discarded term $r_i r_i''$ reads off the same picture: how far the residual bends away from its line ($r_i''$), multiplied by the size of the residual ($r_i$).

With this approximation the step becomes:

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
        step = -np.sum(Jr * r) / np.sum(Jr ** 2)  # (J^T J) Δd = -J^T r
        d = d + step
        d_hist.append(d)
        J_hist.append(objective(d, wavelength_nm, measured_R))
    return np.array(d_hist), np.array(J_hist)
```

This approximate Hessian $J^TJ$ — here $\sum_i r_i'^2$ — has a property that Newton's true Hessian lacks. Being a sum of squares, $J^TJ$ is always positive semidefinite, and positive definite whenever the Jacobian has full rank. In this problem $\sum_i r_i'^2$ is a sum of squared sensitivities, so it is positive as long as a single wavelength responds to thickness. The Gauss-Newton step therefore has no point at which the sign of the curvature can flip, which frees it structurally from the very cause of Newton's divergence in section 3.

## 6. Running Gauss-Newton from the same failure point

The green curve in Figures 2 and 3 is Gauss-Newton. From $d_0=1540$ nm, the exact starting point at which Newton diverged, it reaches the noise floor by the third step — faster even than gradient descent ($\alpha=300$), which took five under the same conditions.

| Method | Result at $d_0=1540$ nm |
|---|---|
| Gradient descent ($\alpha=300$) | 1490.12 nm, converged in 5 steps |
| Newton | no convergence in 10 steps, over 100 nm from the minimum |
| Gauss-Newton | 1490.12 nm, converged in 3 steps |

The shape of the convergence is worth a look. The distance to the minimum falls 49.9 → 20.3 → 1.77 → 0.0036 nm, and from then on shrinks by a factor of about $7.3\times10^{-4}$ per step. A constant ratio means that, strictly speaking, this is linear convergence. The ratio equals the ratio of the discarded term to the kept one at the solution, $\lvert \sum_i r_i r_i'' \rvert / \sum_i r_i'^2$. Because the residual is nothing but noise, the discarded term is small, and the convergence, though linear, is too fast to tell apart from quadratic.

The whole difference comes from one fact: whatever the sign of the true curvature, $\sum_i r_i'^2$ is never negative. That is exactly why metrology software standardizes on Gauss-Newton — or on LM, the subject of the next post — rather than plain Newton. Saving the cost of second derivatives matters, but this structural stability matters more.

## 7. Gauss-Newton is not a cure-all either

Gauss-Newton is not safe in every situation. The discarded term $\sum_i r_i r_i''$ is not always negligible. Three failure conditions are known.

First, if the residuals themselves are large — the model fits the data poorly — the discarded term grows, and far from converging quadratically the method may manage only linear convergence, or diverge. On a zero-residual problem, where the residual vanishes at the solution, Gauss-Newton recovers the same quadratic convergence as Newton. Whether the residual is large can be judged by comparing the RMS of the residuals after the fit with the measurement noise. At the minimum of this problem the residual RMS is 0.0040, the same as the noise standard deviation (0.004). When the residual is at the noise level like this, the problem can safely be treated as close to zero-residual.

Second, if the Jacobian loses rank or comes close to singular — here, if $\sum_i r_i'^2$ approaches zero — the approximate Hessian $J^TJ$ nearly vanishes. The discarded term is then no longer negligible by comparison, and the solution of the normal equation, the step itself, becomes far too large. Section 6 of [post 3](/en/posts/optimization-levenberg-marquardt/) reproduces this situation with a model whose refractive index is wrongly set to 1.02.

Third, if the starting point lies in the basin of a different local minimum, Gauss-Newton converges to that minimum, as it must. Restarting from around 1300 nm and around 1690 nm gives 1293 nm and 1689 nm respectively — precisely the neighbouring local minima visible in Figure 2 of post 1. The stability of Gauss-Newton does not mean "it finds the global minimum"; it means something closer to "once it has a direction, it goes there without wandering."

## Summary

Newton's method promises quadratic convergence in theory. But the promise holds only where the objective is locally convex. On this problem, folded by interference fringes, it fails to converge at all from some initial values. Gauss-Newton exploits the structure of nonlinear least squares to approximate the Hessian by $J^TJ$, removing that failure mode altogether, at the price of slower convergence when the residuals are large. For a least-squares problem, then, there is little reason to use plain Newton. If the objective is a sum of squares and the residuals are small, use Gauss-Newton; if the steps jump around or the sensitivity is weak, use the LM method of the next post. Neither chooses which basin to end up in; [post 4](/en/posts/optimization-global-heuristics/) takes up that problem.

The next post covers LM, which uses a damping parameter to tame the cases where the Gauss-Newton step grows too large or the Jacobian is nearly singular. It is the de facto standard in metrology, and its structure moves smoothly between gradient descent and Gauss-Newton. Measured against the first condition of section 7, what LM prevents is the divergence that large residuals can cause; what remains is the slow convergence, because near the solution it uses the same approximate Hessian as Gauss-Newton.

## References

- D. Bindel, "Gauss-Newton and nonlinear least squares" (lecture notes, 2023-04-10), CS 4220 Numerical Analysis, Cornell University, Spring 2023. [cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf](https://www.cs.cornell.edu/courses/cs4220/2023sp/lec/2023-04-10.pdf)
- K. Madsen, H.B. Nielsen, O. Tingleff, "Methods for Non-Linear Least Squares Problems," 2nd ed., IMM, Technical University of Denmark, 2004. [imm.dtu.dk](https://www2.imm.dtu.dk/pubdb/edoc/imm3215.pdf)

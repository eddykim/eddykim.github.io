---
title: "Optimization 1 — When Thickness Cannot Be Solved Backwards from Reflectance, How Do We Search for It?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-gradient-descent/
page_id: optimization-gradient-descent
date: 2026-09-06 20:30:00 +0900
categories: [Computation, Optimization]
tags: [optimization, gradient-descent, least-squares, thin-film, python]
description: What numerical analysis and optimization theory address, and how gradient descent behaves — and fails — on a thin-film reflectance fitting problem.
math: true
---

Extracting a physical quantity from measurement data often reduces to solving an optimization problem. This series covers the iterative methods used for that purpose, beginning with gradient descent. Before the main discussion, the scope of numerical analysis and optimization theory is defined.

## 1. Numerical analysis and optimization theory

Numerical analysis is the field concerned with obtaining approximate solutions, through a finite number of computations, to mathematical problems whose exact analytic solutions are difficult or impossible to obtain. Root finding, numerical integration, the solution of differential equations, and interpolation between data points all belong to this field. Engineering at large — from predicting semiconductor device characteristics to simulating fluid flow — relies on numerical analysis wherever an algebraically exact solution is unavailable.

Optimization theory is the branch concerned with finding the parameter $x$ that minimizes (or maximizes) an objective function $J(x)$. In general form:

$$ x^* = \arg\min_x J(x) $$

When $J$ is linear or quadratic, $x^*$ can sometimes be obtained algebraically. When $J$ is nonlinear or the dimension of $x$ grows, however, closed-form solutions are the exception rather than the rule. This is where iterative methods are used. Starting from an initial estimate $x_0$, the parameter is updated at each step in a direction that decreases $J(x)$, approaching the minimum. This contrasts with direct methods, which reach the solution in a fixed number of operations.

The series has five parts. This one covers the simplest method, gradient descent. [Part 2](/posts/optimization-newton-gauss-newton/) treats Newton's method and the Gauss-Newton method, and [Part 3](/posts/optimization-levenberg-marquardt/) treats LM (Levenberg-Marquardt), the de facto standard in metrology. [Part 4](/posts/optimization-global-heuristics/) turns to global optimization heuristics that get past local minima, and [Part 5](/posts/optimization-nn-regression/) asks whether neural-network regression can replace iterative optimization.

## 2. Indirect measurement and the need for iteration

There is one reason numerical optimization is needed: when a model is nonlinear, the inverse function that recovers the underlying parameter from the measured value generally has no closed form. This is not specific to optics. Models that describe a phenomenon, in physics or statistics alike, are frequently nonlinear, and the same obstacle appears each time.

A representative example is determining film thickness from a thin-film reflectance spectrum. Since the model equation is given, it is tempting to assume the inverse can be derived, and that expectation is not entirely unfounded. For a film with no or weak absorption, thickness can be obtained in closed form from the extrema of the interference fringes in the transmission or reflectance spectrum; the envelope method family does exactly this. That case, however, is a special one. When absorption is strong enough that the extinction coefficient cannot be neglected, the refractive index becomes complex, $N = n - ik$. More importantly, the measured reflectance $R$ is a nonlinear function of the thickness $d$, so a closed-form inverse yielding $d$ from $R$ does not exist in general.

This is an instance of indirect measurement: when the quantity of interest cannot be measured directly, a model relating it to some measurable quantity is constructed, and the model is inverted to estimate the parameter. Instead of measuring thickness with a ruler, one measures the reflectance spectrum $R(\lambda)$ and inverts the theoretical model $R_{model}(d,\lambda)$ to obtain $d$. When the model cannot be inverted in closed form, as here, the iterative methods of Section 1 are required.

As a concrete structure, light arrives from air (the ambient), passes into the SiO2 layer of interest whose thickness is sought, and the stack is supported by an Si substrate. A structure with a single layer of this kind is called a single-layer sample. Although the structure is simple, the thickness enters the reflectance inside an exponential interference term, so the expression is nonlinear. At normal incidence it takes the following closed form (the Airy formula):

$$ r = \frac{r_{01} + r_{12} e^{-2i\beta}}{1 + r_{01} r_{12} e^{-2i\beta}}, \qquad \beta = \frac{2\pi n_1 d}{\lambda} $$

Here the subscripts 0, 1 and 2 denote air, SiO2 and Si. $r_{01}$ and $r_{12}$ are the Fresnel reflection coefficients at the two interfaces (air/layer, layer/substrate), and $R = \lvert r \rvert^2$ is the quantity actually measured (implementation in `_code/optimization-gradient-descent/reflectance_model.py`). $\beta$ is the phase the light accumulates in one pass through the film. In other words, $e^{-2i\beta}$ is the phase of one round trip, and each time the thickness changes by $\lambda/(2n_1)$ the phase turns through $2\pi$ and the same fringe pattern repeats. With $n_1=1.46$ for SiO2 (a constant that ignores dispersion, used for demonstration), this period is 154–257 nm over the 450–750 nm band used here. That is also why the local minima in Section 5 line up about 200 nm apart. This expression cannot be inverted for $d$ in closed form. Furthermore, the measurement contains noise, so no $d$ matches the model exactly; the problem becomes one of finding the $d$ that fits best.

This is where the least-squares objective function is needed.

$$ J(d) = \frac{1}{2}\sum_i \left( R_{model}(d, \lambda_i) - R_{meas}(\lambda_i) \right)^2 $$

Here $i$ indexes the wavelength points, and the factor $\tfrac{1}{2}$ is a convention that cancels the 2 produced by differentiation.

Why square the residuals and sum them? Other loss functions, including the sum of absolute values, could be used. But if the measurement noise is assumed independent and Gaussian with equal variance, the least-squares solution coincides with the maximum likelihood estimate, which provides a statistical justification. Measurement noise is the sum of several independent noise sources and is therefore often close to Gaussian by the central limit theorem, and for this reason least squares serves as the starting point for the entire series.

The mismatch between the initial model and the measurement is shown below.

<img src="/assets/img/posts/optimization-gradient-descent/en/fig1-model-vs-measurement.png" alt="Initial model assuming 1540 nm versus the measurement" width="600">
_Fig 1. The initial model assuming 1540 nm (blue line) against the synthetic measurement from a 1490 nm sample (grey dots). The fringe positions of the two spectra are misaligned._

The measurement is not raw spectrometer data; it is synthetic. It assumes the 1490 nm sample from an SiO2/Si thickness series (10–190 nm, 1490 nm) the author measured in the lab: the same physical model gives the reflectance at 300 points over 450–750 nm, and Gaussian noise with a standard deviation of 0.004 is added. At 1490 nm the film is optically thick enough that about four fringes ($2n_1d(1/\lambda_1-1/\lambda_2)\approx3.9$) fall within this band, so a change in thickness shows up clearly in the reflectance spectrum. In Fig 1 the fringe positions of the initial model (blue line, $d=1540$ nm) and the measurement (grey dots) are visibly misaligned, which shows how strongly the spectrum responds to a 50 nm change.

Plotting this mismatch as a function of thickness $d$ gives the shape of the objective function $J(d)$.

<img src="/assets/img/posts/optimization-gradient-descent/en/fig2-objective-landscape.png" alt="Objective function against thickness" width="600">
_Fig 2. Objective function J(d) against thickness. The grey dashed line marks the global minimum at 1490 nm. The two shallow valleys on either side are local minima._

The shallow valleys bottom out at 1293 nm and 1689 nm, where $J$ (0.87 and 1.06) is more than two orders of magnitude larger than at the global minimum (0.0024). Their spacing of about 200 nm from 1490 nm is the fringe period from the Airy formula, imprinted on the landscape. Sweeping the full range as in Fig 2 and locating the minimum by eye is not the actual procedure, however (Fig 2 is a wide scan computed in advance to illustrate the principle). Evaluating a single $d$ requires computing the reflectance across the entire wavelength band, so such a grid search is computationally expensive. The real goal is to find the $d$ that minimizes $J(d)$ using as few evaluations as possible.

## 3. What does the gradient tell us?

Section 2 established why an iterative method is needed. For such a method to proceed, each step must determine in which direction, and by how much, to move from the current estimate $d_k$ to the next estimate $d_{k+1}$. Moving in an arbitrary direction may increase $J$, so this decision requires a basis. That basis is provided by the gradient.

The derivative $dJ/dd$ indicates the direction through its sign and, through its magnitude, how steeply $J$ changes at that point. $J$ is large at the initial value in Fig 1 ($d=1540$ nm), so that point is not the optimum, and the gradient there indicates the direction of motion that decreases $J$. In this post the only parameter is thickness, so the gradient is a scalar derivative. In later posts, as parameters multiply (multilayer thicknesses, refractive indices, and so on), this scalar extends to a gradient vector and the second derivative to a Hessian matrix.

Repeatedly moving along the gradient eventually reaches a point with no direction left to move, that is, a point where the gradient is zero. The necessary condition at the minimum of $J(d)$ is therefore:

$$ \frac{dJ}{dd} = 0 $$

This is a necessary but not a sufficient condition. A point where $dJ/dd=0$ may be a global minimum, a local minimum, or even a local maximum such as the peaks in Fig 2. This issue is revisited in Section 5 and treated in detail in Part 4 (global optimization heuristics).

## 4. Gradient descent

The simplest approach is to move a small distance from the current position in the direction opposite the gradient. The justification for the opposite direction follows directly from a first-order Taylor expansion.

$$ J(d + \Delta d) \approx J(d) + \frac{dJ}{dd}\Delta d $$

Choosing $\Delta d$ with sign opposite to the gradient ($\Delta d = -\alpha \, dJ/dd$, $\alpha>0$) makes the second term on the right-hand side always negative, so $J$ decreases at least locally. Repeating this process is gradient descent. Because the direction opposite the gradient is the direction of steepest descent at that point, the method is also called steepest descent.

$$ d_{k+1} = d_k - \alpha \, \frac{dJ}{dd}(d_k) $$

Here $\alpha$ is the step size, called the learning rate in machine learning. Rather than deriving the analytic derivative, the gradient was computed numerically by central difference.

```python
# core of gradient_descent.py (full code: _code/optimization-gradient-descent/)
import numpy as np
from reflectance_model import reflectance

def gradient_descent(d0, wavelength_nm, measured_R, alpha, n_iter=100, h=1e-3):
    def J(d):
        return 0.5 * np.sum((reflectance(d, wavelength_nm) - measured_R) ** 2)

    d = d0
    d_hist, J_hist = [d], [J(d)]
    for _ in range(n_iter):
        grad = (J(d + h) - J(d - h)) / (2 * h)  # central difference for dJ/dd
        d = d - alpha * grad
        d_hist.append(d)
        J_hist.append(J(d))
    return np.array(d_hist), np.array(J_hist)
```

A deliberately large step size was tried first: `alpha=1500`, initial value `d0=1540` nm, 20 iterations. To compare that choice against an appropriate step size, `alpha=300` was run with the same initial value and iteration count and plotted together. Values of $\alpha$ in the hundreds or thousands look large only because of units: $\alpha$ converts the gradient into a distance in nanometres. At 1540 nm the gradient is about 0.078 per nm, so the first step is about 23 nm with `alpha=300` and about 117 nm with `alpha=1500`.

<img src="/assets/img/posts/optimization-gradient-descent/en/fig3-thickness-vs-iteration.png" alt="Thickness estimate versus iteration for two step sizes" width="600">
_Fig 3. Thickness estimate versus iteration. The initial value is 1540 nm; the blue line is $\alpha=300$ and the red line $\alpha=1500$. The blue line reaches 1490 nm in five iterations, while the red line wanders between 1415 and 1566 nm._

<img src="/assets/img/posts/optimization-gradient-descent/en/fig4-objective-vs-iteration.png" alt="Objective function versus iteration for two step sizes" width="600">
_Fig 4. Objective J versus iteration (log scale). Conditions and legend as in Fig 3. The blue line flattens at the noise floor of about 0.0024; the red line rises and falls and never settles at the floor._

The blue curve (`alpha=300`) reaches the true value (1490 nm, within 0.3 nm) after five iterations and holds there. The objective likewise decreases monotonically to roughly 0.0024, the noise floor, and then flattens (Fig 4).

That value is no accident. When the model sits at the true parameter, only noise remains, and the expected value of $J$ is then $\tfrac{1}{2}N\sigma^2$, where $N=300$ is the number of wavelength points and $\sigma=0.004$ the noise standard deviation. This gives $\tfrac{1}{2}\cdot300\cdot0.004^2=0.0024$. Once $J$ flattens near this value, all the mismatch the model can explain is gone, and further iterations will not improve the thickness. This serves as a stopping criterion: stop when $J$ reaches about $\tfrac{1}{2}N\sigma^2$, and suspect the model or the initial value if it flattens well above it.

The red curve (`alpha=1500`) behaves differently. It oscillates by tens to over a hundred nanometres per step — 1540 → 1423 → 1527 → 1416 → 1503 — and fails to converge through all twenty iterations (Fig 3). It passes near the minimum on occasion (at the seventh iteration $J$ falls to approximately 0.03, Fig 4) but departs again at the following step. The situation is that of descending a valley with a stride so long that it overshoots the floor and climbs the opposite slope.

Why the two curves part ways can be explained in one line by treating the floor as a parabola. Near the minimum $d^*$ the gradient is approximately $J''\,(d-d^*)$, and substituting this into the update gives

$$ e_{k+1} = \left(1 - \alpha J''\right) e_k , $$

where $e_k = d_k - d^*$ is the error of the $k$-th estimate and $J''$ the curvature of the floor. Each step multiplies the error by the same factor $1-\alpha J''$. The error shrinks only if that factor has magnitude below 1, which is the condition $\alpha < 2/J''$. With this problem's curvature, $J''\approx0.00258$, the limit is about 775.

For `alpha=300` the factor is 0.23. It keeps its sign, so the iteration never overshoots, and the error falls to nearly a quarter at each step. The error ratios of the actual trajectory settle at 0.23 after the first two steps (0.53, 0.31). For `alpha=1500` the factor is −2.87: each time the iterate comes near the floor, it is thrown nearly three times as far to the other side. Away from the floor, however, the slope flattens again, so the trajectory does not run off indefinitely but wanders between 1415 and 1566 nm. This is the characteristic failure mode of gradient descent under an excessive step size, and the same problem is noted in darkpgmr's article on optimization ([darkpgmr.tistory.com/133](https://darkpgmr.tistory.com/133), in Korean).

The same equation gives a rule for choosing the step size: take $\alpha$ on the order of the inverse curvature. At $\alpha = 1/J'' \approx 390$ the factor is zero and the iterate lands on the floor in one step. Somewhat larger values overshoot and swing back while converging, and beyond twice that (about 775) the iteration fails to converge. Newton's method in [Part 2](/posts/optimization-newton-gauss-newton/) computes exactly this $1/J''$ at every step and puts it in place of $\alpha$.

A single difference in step size separated convergence from endless oscillation.

## 5. Limitations

From this result alone, it might seem enough to choose $\alpha$ well. But the success here rests on a good initial value and a hand-picked $\alpha$. The gradient descent used here depends on two conditions.

First, the initial value must lie inside the basin that leads to the global minimum. As Fig 2 shows, the region near $d=1490$ is not the only minimum. Local minima also exist near $d \approx 1290$ and $d \approx 1690$, one interference fringe away, and starting outside the basin boundaries (approximately $d=1390$ and $d=1590$, the two peaks in Fig 2) causes the algorithm to converge to a different local minimum. This local-minimum problem is not addressed here; it reappears with Gauss-Newton in Part 2 and is taken up directly in Part 4.

Second, the convergence rate depends strongly on $\alpha$, and determining $\alpha$ by trial and error each time is inefficient. As the minimum is approached, $\lvert dJ/dd \rvert$ tends to zero and the steps become progressively smaller. As Section 4 showed, the error shrinks only by a fixed ratio per step (0.23 for `alpha=300`); this is called linear convergence. Set too large from the outset, the iteration oscillates without converging, as in Fig 3.

## Summary

Gradient descent moves downhill in small steps opposite the gradient. The gradient supplies only the direction; how far to go is set by the step size $\alpha$, and the right size for $\alpha$ is set by the curvature of the floor. In this problem convergence requires $\alpha<2/J''\approx775$, and `alpha=300` reached the noise floor in five iterations. The initial value must also lie inside the basin of the global minimum. In short, gradient descent works well only when someone already knows the curvature and the landscape on its behalf.

Newton's method, covered in the next post, determines the step size automatically from the Hessian, at the cost of computing the Hessian and inverting it.

## References

- Dark Programmer, "[Gradient Descent 탐색 방법 (Gradient descent search)](https://darkpgmr.tistory.com/133)" (in Korean). Describes, from the same angle, the oscillation under an excessive step size discussed in Section 4.
- J. Nocedal, S. J. Wright, _Numerical Optimization_, 2nd ed., Springer, 2006. Standard textbook on steepest descent, step sizes and linear convergence.
- S. Boyd, L. Vandenberghe, _Convex Optimization_, Cambridge University Press, 2004. Chapter 9 (unconstrained minimization) covers gradient descent and the step-size conditions under which the objective is guaranteed to decrease.
- H. Fujiwara, _Spectroscopic Ellipsometry: Principles and Applications_, Wiley, 2007. Source of the single-layer reflection coefficient in Section 2 and of the convention $N=n-ik$.
- J. C. Manifacier, J. Gasiot, J. P. Fillard, "A simple method for the determination of the optical constants n, k and the thickness of a weakly absorbing thin film," _J. Phys. E: Sci. Instrum._, 9, 1002, 1976. The original envelope method mentioned in Section 2.

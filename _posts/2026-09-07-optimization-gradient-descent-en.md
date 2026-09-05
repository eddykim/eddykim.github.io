---
title: "Optimization 1 — Least Squares and Gradient Descent"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-gradient-descent/
page_id: optimization-gradient-descent
date: 2026-09-07 20:00:00 +0900
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

When $J$ is linear or quadratic, $x^*$ can sometimes be obtained algebraically. When $J$ is nonlinear or the dimension of $x$ grows, however, closed-form solutions are the exception rather than the rule. This is where iterative methods are used. Starting from an initial estimate $x_0$, the parameter is updated at each step in a direction that decreases $J(x)$, approaching the minimum. This contrasts with direct methods, which produce an exact solution in a single computation.

This series treats gradient descent, Newton's method, the Gauss-Newton method, and Levenberg-Marquardt (LM), which is the de facto standard in metrology, in that order. The subject of this post is the simplest of them, gradient descent.

## 2. Indirect measurement and the need for iteration

The reason numerical optimization is required converges on a single point: when a model is nonlinear, an inverse function that recovers the underlying parameter from the measured value generally does not exist. This is not specific to optics. In physics or statistics alike, models that describe a phenomenon frequently carry nonlinearity, and the same obstacle appears each time. A linear model is solved by inverting a matrix; a nonlinear model is not.

A representative example is determining film thickness from a thin-film reflectance spectrum. Since the model equation is given, it is tempting to assume the inverse can be derived, and that expectation is not entirely unfounded. For a transparent, non-absorbing film with extinction coefficient $k=0$, methods do exist (the envelope method family) that obtain thickness in closed form from the extrema of the interference fringes in the reflectance spectrum. That case, however, is a special one. When the extinction coefficient is nonzero the refractive index becomes complex, $N = n - ik$, and more importantly the thickness $d$ enters the measured reflectance $R$ nonlinearly, so a closed-form inverse yielding $d$ from $R$ does not exist in general.

This is an instance of indirect measurement: when the quantity of interest cannot be measured directly, a model relating it to some measurable quantity is constructed, and the model is inverted to estimate the parameter. Instead of measuring thickness with a ruler, one measures the reflectance spectrum $R(\lambda)$ and inverts the theoretical model $R_{model}(d,\lambda)$ to obtain $d$. If the model is linear and a closed-form inverse exists, the solution follows algebraically (a direct method), but that is far less common. In the remaining cases the iterative methods defined in Section 1 are required.

As a concrete structure, consider air (the ambient layer, through which light is incident) above a single SiO2 layer of interest, supported by an Si substrate. A structure with a single layer of this kind is called a single-layer sample. Although the structure is simple, the reflectance is a nonlinear function of thickness containing an exponential interference term. At normal incidence it takes the following closed form (the Airy formula):

$$ r = \frac{r_{01} + r_{12} e^{-2i\beta}}{1 + r_{01} r_{12} e^{-2i\beta}}, \qquad \beta = \frac{2\pi n_1 d}{\lambda} $$

Here $r_{01}$ and $r_{12}$ are the Fresnel reflection coefficients at each interface (ambient/layer, layer/substrate), and $R = \lvert r \rvert^2$ is the quantity actually measured (implementation in `_code/optimization-gradient-descent/reflectance_model.py`). This expression cannot be inverted for $d$ in closed form. Furthermore, the measurement contains noise, so no $d$ matches the model exactly; the problem becomes one of finding the $d$ that fits best.

This is where the least-squares objective function is needed.

$$ J(d) = \frac{1}{2}\sum_i \left( R_{model}(d, \lambda_i) - R_{meas}(\lambda_i) \right)^2 $$

Why square the residuals and sum them? Other loss functions, including the sum of absolute values, could be used. But if the measurement noise is assumed Gaussian, the least-squares solution coincides with the maximum likelihood estimate, which provides a statistical justification. Measurement noise is the sum of several independent noise sources and is therefore often close to Gaussian by the central limit theorem, and for this reason least squares serves as the starting point for the entire series.

The goal is to find the $d$ that minimizes the sum of squared differences between model and measurement. The mismatch between the initial model and the measurement is shown below.

<img src="/assets/img/posts/optimization-gradient-descent/en/fig1-model-vs-measurement.png" alt="Initial model versus measurement" width="600">
_Fig 1. Initial model assuming 1540 nm against the measurement from a 1490 nm sample_

The measurement is not raw spectrometer data. It was synthesized by assuming a 1490 nm sample from a previously measured SiO2/Si thickness series (10–190 nm, 1490 nm) and adding Gaussian noise to the same physical model. At 1490 nm the film is optically thick enough that several interference fringes fall within the visible band, so a change in thickness is clearly reflected in the reflectance spectrum. In Fig 1 the fringe positions of the initial model (blue line, $d=1540$ nm) and the measurement (grey dots) are visibly misaligned, which demonstrates this.

Plotting this mismatch as a function of thickness $d$ gives the shape of the objective function $J(d)$.

<img src="/assets/img/posts/optimization-gradient-descent/en/fig2-objective-landscape.png" alt="Objective function against thickness" width="600">
_Fig 2. Objective function J(d) against thickness — the grey dashed line marks the global minimum_

Sweeping the full range as in Fig 2 and locating the minimum by eye is not the actual procedure (Fig 2 is a wide scan computed in advance to illustrate the principle). Evaluating a single $d$ requires computing the reflectance across the entire wavelength band, so such a grid search is computationally expensive. The real goal is to find the $d$ that minimizes $J(d)$ using as few evaluations as possible.

## 3. The objective function and the gradient condition

Section 2 established why an iterative method is needed. For such a method to proceed, each step must determine in which direction, and by how much, to move from the current estimate $d_k$ to the next estimate $d_{k+1}$. Moving in an arbitrary direction may increase $J$, so this decision requires a basis. That basis is provided by the gradient.

The derivative $dJ/dd$ indicates the direction through its sign and the severity of the current mismatch through its magnitude. That $J$ is large at the initial value in Fig 1 ($d=1540$ nm) means that point is not the optimum, and the gradient at that point indicates the direction of motion that decreases $J$. In this post the only parameter is thickness, so the gradient is a scalar derivative. In later posts, as parameters multiply (multilayer thicknesses, refractive indices, and so on), this scalar extends to a gradient vector and the second derivative to a Hessian matrix.

Repeatedly moving along the gradient eventually reaches a point with no direction left to move, that is, a point where the gradient is zero. The necessary condition at the minimum of $J(d)$ is therefore:

$$ \frac{dJ}{dd} = 0 $$

This is a necessary but not a sufficient condition. A point where $dJ/dd=0$ may be a global minimum or a local minimum. This issue is revisited in Section 5 and treated in detail in the next post (Newton and Gauss-Newton).

## 4. Gradient descent

The simplest approach is to move a small distance from the current position in the direction opposite the gradient. The justification for the opposite direction follows directly from a first-order Taylor expansion.

$$ J(d + \Delta d) \approx J(d) + \frac{dJ}{dd}\Delta d $$

Choosing $\Delta d$ with sign opposite to the gradient ($\Delta d = -\alpha \, dJ/dd$, $\alpha>0$) makes the second term on the right-hand side always negative, guaranteeing that $J$ decreases at least locally. Repeating this process is gradient descent.

$$ d_{k+1} = d_k - \alpha \, \frac{dJ}{dd}(d_k) $$

Here $\alpha$ is the step size (learning rate). Rather than deriving the analytic derivative, the gradient was computed numerically by central difference.

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

A deliberately large step size was tried first: `alpha=1500`, initial value `d0=1540` nm, 20 iterations. To compare that choice against an appropriate step size, `alpha=300` was run with the same initial value and iteration count and plotted together.

<img src="/assets/img/posts/optimization-gradient-descent/en/fig3-thickness-vs-iteration.png" alt="Thickness estimate versus iteration for two step sizes" width="600">
_Fig 3. Thickness estimate vs iteration — step size comparison_

<img src="/assets/img/posts/optimization-gradient-descent/en/fig4-objective-vs-iteration.png" alt="Objective function versus iteration for two step sizes" width="600">
_Fig 4. Objective J vs iteration — step size comparison_

The blue curve (`alpha=300`) reaches the true value (1490 nm, within 0.3 nm) after five iterations and holds there. The objective likewise decreases monotonically to roughly 0.0024, the noise floor, and then flattens (Fig 4).

The red curve (`alpha=1500`) behaves differently. It oscillates by tens to over a hundred nanometres per step — 1540 → 1423 → 1527 → 1416 → 1503 — and fails to converge through all twenty iterations (Fig 3). It passes near the minimum on occasion (at the seventh iteration $J$ falls to approximately 0.03, Fig 4) but departs again at the following step. The situation is that of descending a valley with a stride so long that it overshoots the floor and climbs the opposite slope. Since the gradient on the far side is steeper, the following step overshoots by even more. This is the characteristic failure mode of gradient descent under an excessive step size, and the same problem is noted in darkpgmr's article on optimization ([darkpgmr.tistory.com/133](https://darkpgmr.tistory.com/133), in Korean).

A single difference in step size separated convergence from divergence.

## 5. Limitations

The gradient descent used here depends on two conditions.

First, the initial value must lie inside the basin that leads to the global minimum. As Fig 2 shows, the region near $d=1490$ is not the only minimum. Local minima also exist near $d \approx 1290$ and $d \approx 1690$, one interference fringe away, and starting outside the basin boundaries (approximately $d=1390$ and $d=1590$, the two peaks in Fig 2) causes the algorithm to converge to a different local minimum. This local-minimum problem is not addressed here; it is discussed in the next post (Newton and Gauss-Newton).

Second, the convergence rate depends strongly on $\alpha$, and determining $\alpha$ by trial and error each time is inefficient. As the minimum is approached, $\lvert dJ/dd \rvert$ tends to zero and the steps become progressively smaller; set too large from the outset, the iteration diverges as in Fig 3. Newton's method, covered in the next post, determines this step size automatically from the Hessian, at the cost of computing the Hessian and inverting it.

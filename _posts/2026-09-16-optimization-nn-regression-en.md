---
title: "Optimization 5 — Neural Network Regression: Can Learning Replace Iterative Optimization?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/optimization-nn-regression/
page_id: optimization-nn-regression
date: 2026-09-16 20:00:00 +0900
categories: [Computation, Optimization]
tags: [optimization, neural-network, regression, inverse-problem, thin-film, python]
description: A neural network regressor is built in numpy and pitted against Levenberg-Marquardt on thin-film thickness fitting, to find where it wins, where it loses, and why losing points to its proper role.
math: true
---

[Post 4](/en/posts/optimization-global-heuristics/) ended on a practical note: real metrology software locates the approximate position with a global search, then hands the final precision to a local method such as LM. That sentence contains a gap. **Must the "approximate position" be found by searching at all?**

Every method in posts 1 through 4 was iterative. A measured spectrum arrives, and only then does the work begin — evaluate the objective, look at the gradient, take a step, evaluate again. When the next sample arrives, the whole process starts over. Yet the problem we solve is the same one every time: given a reflectance spectrum from an Air/SiO2/Si single layer, return a thickness. If that task repeats tens of thousands of times, could the input-output relation itself be learned once, so that later samples are answered directly?

This post takes up that question. A neural network that maps 300 spectral points to a single thickness is implemented from scratch in numpy, measured against LM on the same terms, and examined for where it wins, where it loses, and **why the losing side is what reveals the method's proper use.**

The example is unchanged from posts 1 through 4. Normal-incidence reflectance of a single SiO2 layer is sampled at 300 wavelengths from 450 to 750 nm (true thickness 1490 nm, Gaussian noise of standard deviation 0.004), and the objective is $J(d) = \frac{1}{2}\sum_i r_i^2$. Three local minima sit in the 1200-1800 nm range, at 1293, 1490 and 1689 nm, and the basin surrounding each is 198 nm wide.

## 1. Why the landscape argument changes

One thing deserves attention first. Local minima were the difficulty throughout posts 1 to 4. Yet the neural network literature does not treat local minima as nearly so serious a problem. Why?

Because the number of unknowns differs. Posts 1 to 4 had a single unknown, the thickness $d$. The network built here has tens of thousands of weights. Dauphin and co-workers argued that the real obstacle in high-dimensional non-convex minimization is not local minima but saddle points, and one strand of their argument is a counting argument over dimension. For a critical point to be a local minimum, the surface must curve upward along **every** direction. As the dimension $N$ grows, the probability of satisfying that condition in all directions falls exponentially, and saddle points — up along some directions, down along others — come to dominate.

Read the same argument backwards and our situation appears. **At $N=1$, an exact saddle point is a probability-zero event.** In a landscape with a single unknown, every critical point is an extremum, which is precisely why three basins caused so much trouble in posts 1 to 4. Weight space is a different world. That is why training the network below never raises the question of whether it is stuck in a local minimum.

One caveat matters, though. **The ambiguity of the problem does not disappear.** That 1293 nm and 1490 nm produce similar spectra is a fact of the measurement, not of the algorithm. In a neural network that ambiguity simply wears a different face, and section 8 shows which one.

## 2. Building the regressor in numpy

This series has transcribed algorithms directly rather than calling libraries. LM was written out as Algorithm 3.16 of Madsen, Nielsen and Tingleff; here the forward pass, backpropagation and Adam are likewise transcribed from the original papers into numpy. A fully connected network needs less machinery than one might expect.

Consider the training data first, where this problem shows its peculiarity — **the forward model is free.** The Airy formula is a single line, so a spectrum can be generated for any thickness, in any quantity. This inverts the usual situation in machine learning, where collecting data is the expensive step. Eight thousand thicknesses were drawn uniformly from $[1200, 1800]$ nm and their spectra synthesized.

```python
rng = np.random.default_rng(0)
d = rng.uniform(1200.0, 1800.0, size=8000)                # labels
X = np.stack([reflectance(x, wavelength_nm) for x in d])  # inputs, (8000, 300)
```

The network is plain: 300 inputs, 256 hidden units with tanh, one linear output. Weights follow the normalized initialization of Glorot and Bengio, which for a layer with $n_\text{in}$ inputs and $n_\text{out}$ outputs draws

$$W \sim U\left[-\sqrt{\frac{6}{n_\text{in}+n_\text{out}}},\ +\sqrt{\frac{6}{n_\text{in}+n_\text{out}}}\right]$$

The optimizer is Adam. Algorithm 1 of the paper is written as pseudocode and transcribes directly: accumulate first and second moments $m$ and $v$ as exponential moving averages, correct their bias, and divide.

```python
self._t += 1
m = b1 * m + (1 - b1) * grad
v = b2 * v + (1 - b2) * grad ** 2
W -= lr * (m / (1 - b1 ** self._t)) / (np.sqrt(v / (1 - b2 ** self._t)) + eps)
```

One practical point. Input reflectances span 0.09 to 0.35 while output thicknesses span 1200 to 1800, so without normalization the network does not train at all. Inputs were standardized per wavelength channel and outputs mapped linearly to $[-1,1]$, then inverted after prediction. The full implementation is in [`_code/optimization-nn-regression/`](https://github.com/eddykim/eddykim.github.io/tree/main/_code/optimization-nn-regression).

### A bigger network is not a better one

Three architectures were trained on identical data for 600 epochs each. The table reports both the error after running to the end and the error when the weights are restored to the epoch at which validation error was lowest.

| Architecture | After 600 epochs | Restored to best epoch | Best epoch |
| --- | ---: | ---: | ---: |
| 300 - 128 - 1 | 5.99 nm | 0.80 nm | 465 |
| **300 - 256 - 1** | 2.61 nm | **0.51 nm** | 505 |
| 300 - 256 - 128 - 1 | 2.64 nm | 0.63 nm | **55** |

Two things stand out.

First, **one hidden layer outperformed two** (0.51 against 0.63 nm). Lee and Jin reached the same conclusion when they validated neural thickness measurement against certified reference materials: configurations with more nodes and more hidden layers drifted further from the certified values, which the authors attributed to overtraining on idealized synthetic spectra.

Second, **overtraining is severe.** Run to the end, all three architectures fall back to three to seven times their best validation error. The deepest network reaches its minimum at epoch 55 and degrades from there — the more capacity available, the faster the theoretical spectra are memorized.

The remedy was written down thirty years ago. Fried and Masa, pre-evaluating spectroscopic ellipsometry data with a neural network in 1994, noted that training error keeps falling while test error begins to rise past some point, so the network must be tested periodically against unfamiliar data during training to locate the optimum. Every experiment below restores the weights from the epoch of lowest validation error.

## 3. Competing with LM — and losing

Feed the trained network the measured spectrum of the 1490 nm sample, then run LM on the same data.

| Method | Estimated thickness (nm) | Error (nm) | Final $J$ | Iterations |
| --- | ---: | ---: | ---: | ---: |
| Network | 1490.403 | +0.403 | 2.488e-03 | – |
| **LM ($d_0=1480$)** | **1490.122** | **+0.122** | **2.386e-03** | 4 |
| LM ($d_0=1300$) | 1293.125 | −196.875 | 8.670e-01 | 30 |

LM is more than three times as accurate. The last row is the failure already seen in post 3: with its starting point in a neighbouring basin, LM cannot leave.

Accuracy and precision must be **separated** here. Taking the same 1490 nm sample and varying only the noise realization across 50 measurements, the mean error (accuracy) and the scatter (precision) were recorded separately. Because a network also depends on its training seed, six were trained and the results reported as a range.

| Method | Error (nm) | Scatter, std (nm) |
| --- | ---: | ---: |
| Network (6 training seeds) | −0.98 to +0.33 | **0.10 to 2.68** |
| **LM** | **+0.001** | **0.081** |

LM is deterministic, so its value is fixed. The network is worse in both error and scatter. The striking part, though, is that **the scatter varies by a factor of 26 across training seeds.** Worse, the validation RMSE gives no warning of it: the network with a validation RMSE of 0.34 nm scattered by 2.68 nm, while the one at 0.53 nm scattered by 0.10 nm.

This is the most practical finding in the post. **A good aggregate metric does not guarantee the reliability of any individual sample.** Validation RMSE averages over thousands of cases; what we want to know is the thickness of the one wafer in front of us. The two come apart.

Stopping here, the conclusion would be that neural networks are useless. But only half the picture has been seen.

## 4. Except that LM is wrong two times out of three

The comparison in section 3 carried a hidden assumption. LM was handed $d_0=1480$, **a good starting value, for free.** Without knowing that the answer is 1490, how would anyone choose 1480?

Starting values were drawn at random from $[1200, 1800]$ over 100 trials. The conditions follow the structure of the table Fried and Masa published in 1994.

| Condition | Success rate | Median \|error\| (nm) | Mean iterations |
| --- | ---: | ---: | ---: |
| (a) Network alone | 100% | 0.930 | – |
| **(b) Network → LM** | **100%** | **0.060** | 6.6 |
| (c) Random start → LM | 37% | 196.873 | 22.8 |
| (d) $d_0=1300$ → LM | 0% | 196.945 | 29.8 |

"Success" is defined as in post 4: arriving within 5 nm of the global minimum.

The picture inverts. LM from a random start **returns the wrong answer roughly two times in three.** The figure of 37% matches the 36% measured in post 4, and corresponds roughly to the fraction of the search range (600 nm) occupied by the correct basin (198 nm). The errors on failure cluster consistently near **197 nm**, which is the signature of choosing an entire basin wrongly.

Seeded with the network's output, LM succeeds every time, and its median error of 0.060 nm is **15 times finer than the network alone** (0.930 nm). The network cannot deliver an accurate answer, but it **does pick the right basin**, and picking a basin demands a precision of only ±99 nm.

The weakness found in section 3 — an untrustworthy final nanometre, scattering by a factor of 26 with the training seed — therefore does not matter here. The precision required is two orders of magnitude coarser. An unreliable instrument has been handed a job where unreliability is affordable.

That was exactly the conclusion Fried and Masa drew in 1994. Neural pre-evaluation is neither more accurate nor more precise than ordinary multiparameter fitting, but it supplies such a good initial estimate for the subsequent regression that the time-consuming grid search becomes unnecessary. In their problem, an initial thickness wrong by as little as 50 nm was enough to land in a false minimum.

## 5. Placed beside the global searches of post 4

But post 4 already solved the basin problem. Simulated Annealing and Basin-Hopping both reached success rates near 100%. If the network also reaches 100%, nothing new has been gained.

By success rate, nothing has. **The difference lies in cost.**

Measuring fairly calls for a better yardstick than wall-clock time. The count of **forward-model evaluations** does not depend on implementation or hardware: how many times was the Airy formula called? That is a property of the algorithm.

<img src="/assets/img/posts/optimization-nn-regression/en/fig1-success-vs-cost.png" alt="Success rate and per-sample cost" width="600">
_Fig 1. Success rate and per-sample cost (100 random starts)_

| Method | Success rate | Forward-model calls / sample | Time / sample |
| --- | ---: | ---: | ---: |
| Single LM (random start) | 37% | 119 | 1.03 ms |
| Simulated Annealing | 93% | 301 | 2.96 ms |
| Basin-Hopping | 99% | 1,806 | 15.85 ms |
| **Network → LM** | **100%** | **35** | **0.34 ms** |

**The network does not beat basin-hopping on success rate.** It does the same job at 1/52 the forward-model calls and 1/47 the time. That is the real position of the method that "lost" in section 3.

It is not free, of course. Generating 8,000 training spectra costs 8,000 forward-model calls, and Adam training costs 31 seconds. How many samples it takes to recover that one-time cost depends on which yardstick is used.

| Yardstick | vs Simulated Annealing | vs Basin-Hopping |
| --- | ---: | ---: |
| Forward-model calls | 30.1 samples | **4.5 samples** |
| Wall clock (including 31 s of training) | 11,790 samples | 2,008 samples |

The two numbers diverge because the forward model in this example — one line of Airy formula — is far too cheap. Recovering 31 seconds of training through saved reflectance evaluations takes 2,000 samples. Put differently, **this example is the case least favourable to a neural network.** Where a single forward evaluation takes seconds, as in transfer-matrix calculations over 200-plus layers or rigorous coupled-wave analysis (RCWA), the call-count yardstick governs and recovery comes far sooner. This is the point Kwak and Kim make in their review of semiconductor multilayer metrology when they note that machine learning computation time does not scale with the number of layers.

## 6. Where it breaks — and what goes unreported

So far the training and measurement distributions have matched. What happens when they diverge? Three ways of breaking them were tried.

**Self-checking** was examined at the same time. Unlike LM, the network carries no internal residual signal. But the forward model is free, so the predicted thickness can be pushed back through the Airy formula and the residual recomputed. Whether that residual raises an alarm is the question.

<img src="/assets/img/posts/optimization-nn-regression/en/fig2-residual-alarm.png" alt="What recomputing the residual does and does not reveal" width="600">
_Fig 2. What recomputing the residual does and does not reveal_

| Case | Network error | Recomputed residual $J$ | LM error |
| --- | ---: | ---: | ---: |
| (a) Outside training range (2000 nm) | **−382.3 nm** | **1,130×** baseline | −469.5 nm |
| (b) Noise mismatch (std 0.02) | +1.5 nm | 23× | +0.6 nm |
| (c) Model error ($n_1$=1.50) | **+40.9 nm** | **4.4×** | **+40.9 nm** |

The three cases behave entirely differently.

**(a) Extrapolation fails badly, and the residual says so loudly.** Given a 2000 nm sample from outside the training range $[1200,1800]$, the network returns 1617.7 nm — not even pinned to the training boundary at 1800, but some arbitrary value inside the range. This is precisely what Fried and Masa warned of when they stated that a network extrapolates poorly and the training data must therefore cover the entire expected input space. Fortunately the residual rises to 1,130 times baseline, so the anomaly is immediately visible.

**(b) The noise mismatch was survived, contrary to expectation.** Five times the training noise left an error of only 1.5 nm. This had been listed as a failure mode while preparing the post, and it turned out to do almost no harm. Noise scatters randomly across the whole spectrum and so does not systematically corrupt the fringe positions that determine thickness.

**(c) Model error deceives both methods equally.** This is the most important result. Given a spectrum generated with the SiO2 index wrongly assumed to be 1.50 rather than 1.46, the network is off by +40.87 nm and LM by +40.91 nm. **The two answers agree to within 0.04 nm while both are wrong by 41 nm.** The residual, meanwhile, rises only to 4.4 times baseline — hardly an alarm next to the 1,130 of the extrapolation case.

Post 3 showed Gauss-Newton collapsing and LM rescuing it when the refractive index was wrongly assumed to be 1.02. There the difference between algorithms decided the outcome. Here it does not. **When the forward model itself is wrong, every method built on that model is wrong in the same direction by the same amount.** The network inherits the faulty physics through its training data, LM through its fitting model.

The contrast "the network has no self-check, LM does" is therefore only half true. Recomputing the residual catches distribution shift but largely misses model error, and that is **a blind spot the two methods share**.

## 7. What augmentation matched to real measurements repairs

Cases (a) and (c) differ in kind. The first calls for widening the training distribution; the second calls for correcting the physics. So how far does widening the training distribution toward real measurements actually get?

The augmentation axes follow the setup Kwak and co-workers used for 3D multilayer semiconductor metrology. Unable to obtain many normal samples from a production line, they expanded 125 into 5,000. What deserves attention is that the perturbations were not arbitrary noise but **the actual error modes of a spectrometer** — a vertical offset and a lateral (wavelength-axis) offset. The same axes were carried over here.

- **Noise**: a standard deviation drawn per sample from $[0, 0.02]$ and injected
- **Vertical**: $R \to (1+s)R + c$, with scale $s$ at ±5% and baseline $c$ at ±0.04
- **Lateral**: the wavelength axis shifted by ±6 nm

<img src="/assets/img/posts/optimization-nn-regression/en/fig3-augmentation-rmse.png" alt="Errors augmentation fixes, and errors it does not" width="700">
_Fig 3. Errors augmentation fixes, and errors it does not_

The test conditions must include the error modes the augmentation targets. RMSE was computed over 200 random thicknesses per condition.

| Training | Normal | 5× noise | Model error | **Wavelength +4 nm** | **Photometric ×1.03+0.02** |
| --- | ---: | ---: | ---: | ---: | ---: |
| 8,000 thicknesses, no augmentation | 2.72 | 8.50 | 40.50 | 14.65 | 16.13 |
| 8,000 thicknesses, **augmented** | 2.53 | 6.28 | 41.30 | **2.68** | **2.56** |
| 250 thicknesses, no augmentation | 9.69 | 13.48 | 39.48 | 21.22 | 22.24 |
| 250 thicknesses, **augmented** | **2.15** | 7.08 | 41.25 | **2.39** | **2.66** |

The result splits three ways.

Augmentation works **only on the perturbations it covers, but there it works dramatically.** On a spectrometer whose wavelength axis is off by 4 nm, RMSE falls from 14.65 to 2.68 nm, a factor of 5.5; on a photometric offset, from 16.13 to 2.56 nm, a factor of 6.3.

The price paid on clean data is **almost nothing.** With thicknesses plentiful, 2.72 nm becomes 2.53 nm, marginally better. In the sparse regime of only 250 training thicknesses, 9.69 nm becomes 2.15 nm, a large gain — here augmentation is also **acting as regularization**. That is exactly why Kwak and co-workers expanded 125 samples into 5,000.

On model error it does **nothing at all.** 40.50 nm becomes 41.30 nm, marginally worse. No amount of noise and offset variety repairs a forward model that is itself wrong. The conclusion matches case (c) of section 6.

In short, augmentation covers **the imperfections of the instrument**, not **the errors of the physics**. Robustness extends exactly as far as the error modes already known to the user.

## 8. Given ambiguity, the network returns an average

Now back to what section 1 left aside. What face does ambiguity wear in a neural network?

What a least-squares regressor converges to was settled by Bishop in 1994. In the limit of infinite data, minimizing squared error drives the network function to the **mean** of the target conditioned on the input.

$$f(\mathbf{x}; \mathbf{w}^*) = \langle t \mid \mathbf{x} \rangle = \int t\, p(t \mid \mathbf{x})\, dt$$

For classification this quantity is the posterior probability and therefore optimal, but inverse problems are another matter. When one input corresponds to several correct answers, **the average of several correct answers is not guaranteed to be correct itself.** Bishop noted that in such cases the conventional least-squares approach can give completely erroneous results.

Test it on our problem by deliberately amplifying the ambiguity. Normal-incidence reflectance repeats almost identically each time the thickness changes by $\lambda / 2n_1 \approx 205$ nm, so if the wavelength band is too narrow to contain even one fringe, thicknesses separated by multiples of that period become indistinguishable. About three periods fit inside the 600 nm training range.

<img src="/assets/img/posts/optimization-nn-regression/en/fig4-conditional-mean.png" alt="What the network returns on ambiguous input" width="600">
_Fig 4. On ambiguous input the network returns the conditional mean_

The band was narrowed to 595-605 nm (0.12 fringes) with the noise raised fivefold, and samples of various thicknesses were fed in.

| True thickness (nm) | Network prediction (nm) | Error (nm) |
| ---: | ---: | ---: |
| 1250 | 1415.9 | **+165.9** |
| 1350 | 1426.5 | +76.5 |
| 1450 | 1443.7 | −6.3 |
| **1500** | **1497.3** | **−2.7** |
| 1550 | 1489.4 | −60.6 |
| 1650 | 1456.2 | −193.8 |
| 1750 | 1574.7 | **−175.3** |

Bishop's prediction appears exactly. At 1500 nm, the centre of the training range, prediction nearly coincides with truth (−2.7 nm), because the candidates are symmetric about the centre and their average lands on the true value. Toward the edges, however, predictions are **pulled toward the centre by as much as 194 nm.** A 1750 nm sample reads as 1574.7 nm, a 1250 nm sample as 1415.9 nm.

The point worth attention is that **the value pulled to is not any basin at all.** 1574.7 nm is neither 1490 nor 1689, but some value between the candidates. No thickness with that solution exists physically; only the arithmetic mean sits there. LM would have converged to one of the three — wrong, perhaps, but at a physically meaningful value.

Repeating the test on the wide band (450-750 nm, 3.9 fringes) leaves every error below 2.5 nm. **The effect appears only when the ambiguity is real.** Nothing is wrong with the network; by the definition of the loss it has no choice but to average over what the data does not distinguish.

## 9. Summary — ambiguity is a property of the data, not the algorithm

Invert the conclusion of section 8 and the direction of the remedy follows. If the ambiguity comes from the data, the data is where it must be fixed. So does collecting more of it suffice?

The simplest approach was tried first: widen the wavelength range from 450-750 nm to 450-1000 nm, holding the number of points at 300.

| Range | Local minima | LM success, random start | Network validation RMSE |
| --- | ---: | ---: | ---: |
| 450-750 nm | 3 | 28% | 0.55 nm |
| 450-1000 nm | **3** | 37% | 0.45 nm |

**The count of local minima did not fall.** Because the order period scales with the centre wavelength (205 nm at 600 nm, 248 nm at 725 nm), widening the range leaves the number of basins unchanged.

Add an unknown to see it more sharply. What happens if the refractive index $n_1$ is fitted alongside the thickness? Now the **cross-correlation** of the two parameters becomes visible: from the covariance $C = (J^\top J)^{-1}$ it is defined as $\rho = C_{01}/\sqrt{C_{00}C_{11}}$.

| Range | $\rho(d, n_1)$ | cond($J^\top J$) | $d$ scatter (nm) |
| --- | ---: | ---: | ---: |
| 450-750 nm | −0.997866 | 2.47e+08 | 0.962 |
| 450-1000 nm | −0.996694 | 1.59e+08 | 0.967 |

Widening the range moves the correlation only from −0.9979 to −0.9967, and the scatter does not change at all. The reason lies in the physics. At normal incidence the phase is $\beta = 2\pi n_1 d / \lambda$, so reflectance depends chiefly on the **product of $n_1$ and $d$**. However much data of the same kind is collected, none of it constrains that degenerate direction.

### Adding a different kind of measurement

A different kind is needed, then. Consider oblique incidence. At an incidence angle $\theta_0$ the refraction angle $\theta_1$ follows from Snell's law, and $\cos\theta_1$ enters the phase.

$$\beta = \frac{2\pi n_1 d \cos\theta_1}{\lambda}, \qquad n_0 \sin\theta_0 = n_1 \sin\theta_1$$

Because $\cos\theta_1$ itself depends on $n_1$, the index now enters the phase differently from the thickness. The Fresnel coefficients also split into s and p polarizations, adding amplitude information. Changing the angle, in other words, creates **new information capable of separating $n_1$ from $d$.**

Holding the number of data points fixed at 300, only the kind of measurement was varied.

<img src="/assets/img/posts/optimization-nn-regression/en/fig5-correlation-ellipse.png" alt="(d, n1) estimate scatter by kind of measurement" width="600">
_Fig 5. (d, n1) estimate scatter — same 300 points, different kinds of measurement_

| Configuration (300 points total) | $\rho(d, n_1)$ | cond($J^\top J$) | $d$ scatter (nm) | $n_1$ scatter |
| --- | ---: | ---: | ---: | ---: |
| Normal incidence only, 450-750 nm | −0.997866 | 2.47e+08 | 1.303 | 0.00125 |
| Normal incidence only, 450-1000 nm | −0.996694 | 1.59e+08 | 1.176 | 0.00118 |
| **0°+60°, 150 points each** | **−0.976132** | **3.45e+07** | **0.458** | **0.00034** |
| **0°+45°+70°, 100 points each** | −0.975515 | 3.84e+07 | **0.340** | **0.00024** |

The same 300 points give different answers. Widening the wavelength range improves the thickness scatter by 9%, while **mixing angles improves it by a factor of 3.8 and the condition number by a factor of 6.4.**

This is the conclusion of the post. **What is needed is not "enough data" but "data that sees the degenerate direction."** Kwak and Kim make exactly this point for semiconductor multilayer metrology when they write that single-wavelength ellipsometry does not yield a unique multilayer thickness solution, so multiangle or multisample measurements are used alongside it, and that this additional information improves the uniqueness of the characterization.

And this is not a neural network story. Johnson, treating parameter cross-correlation in curve fitting generally, observed that once the correlation is nearly complete the estimation procedure cannot identify unique parameter values or their standard errors. The correlation coefficient in his two-exponential example was −0.999933 — the same kind of degeneracy as the −0.9979 seen here. More important is the qualification he adds: a correlation between two parameters does not imply that anything in the underlying mechanism is correlated, and **may simply be a consequence of estimating from a limited set of observations.**

### Closing the series

From gradient descent in post 1 to here, five posts have solved the same problem five ways. Posts 1 through 4 were about how well a given set of data can be mined: how to use the gradient, how to approximate the curvature, how to tune the damping, how to cross a valley.

Post 5 adds a layer on top. A learned function **cannot replace** iterative optimization — it loses on accuracy, its precision swings with the training seed, and it collapses outside its training range. What it can do is **pick the basin**, matching basin-hopping's success rate at 1/52 of the cost. Its place as an initial-value supplier rather than a replacement was already written down thirty years ago, and these experiments confirm it again.

There is one thing none of the five can do. **Create information that is not there.** In a measurement where index and thickness enter only as a product, neither LM nor the network can separate them. LM converges wherever the starting value sends it; the network quietly returns an average. What is needed then is not a better algorithm but one more angle of incidence.

The optimization story ends here, but the question of designing *what to measure* in the first place is a subject of its own. A later series will take it up.

## References

- M. Fried, P. Masa, "Backpropagation (neural) networks for fast pre-evaluation of spectroscopic ellipsometric measurements," _J. Appl. Phys._, 75(4), 2194-2201, 1994. DOI: 10.1063/1.356281
- J. Lee, J. Jin, "[A novel method to design and evaluate artificial neural network for thin film thickness measurement traceable to the length standard](https://doi.org/10.1038/s41598-022-06247-y)," _Scientific Reports_, 12, 2212, 2022.
- H. Kwak, J. Kim, "[Semiconductor Multilayer Nanometrology with Machine Learning](https://doi.org/10.1007/s41871-023-00193-7)," _Nanomanufacturing and Metrology_, 6, 15, 2023.
- H. Kwak et al., "[Non-destructive thickness characterisation of 3D multilayer semiconductor devices using optical spectral measurements and machine learning](https://doi.org/10.37188/lam.2021.001)," _Light: Advanced Manufacturing_, 2(1), 9-19, 2021.
- C. M. Bishop, "[Mixture Density Networks](https://publications.aston.ac.uk/id/eprint/373/1/NCRG_94_004.pdf)," Neural Computing Research Group Report NCRG/94/004, Aston University, 1994.
- M. L. Johnson, "[Parameter correlations while curve fitting](https://doi.org/10.1016/S0076-6879(00)21207-X)," _Methods in Enzymology_, 321, 424-446, 2000.
- Y. N. Dauphin et al., "[Identifying and attacking the saddle point problem in high-dimensional non-convex optimization](https://arxiv.org/abs/1406.2572)," _NIPS_, 2014.
- X. Glorot, Y. Bengio, "[Understanding the difficulty of training deep feedforward neural networks](https://proceedings.mlr.press/v9/glorot10a.html)," _AISTATS_, PMLR 9, 249-256, 2010.
- D. P. Kingma, J. Ba, "[Adam: A Method for Stochastic Optimization](https://arxiv.org/abs/1412.6980)," _ICLR_, 2015.

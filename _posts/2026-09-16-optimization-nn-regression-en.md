---
title: "Optimization 5 — Can a Neural Network Replace Iterative Fitting?"
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

[Post 4](/en/posts/optimization-global-heuristics/) ended on a practical note: real metrology software locates the approximate position with a global search, then hands the final precision to a local method such as LM (Levenberg-Marquardt). That sentence contains a gap. **Must the "approximate position" be found by searching at all?**

Every method in posts 1 through 4 was iterative. A measured spectrum arrives, and only then does the work begin: evaluate the objective, look at the gradient, take a step, evaluate again. When the next sample arrives, the whole process starts over. Yet the problem solved here is the same one every time: given a reflectance spectrum from an Air/SiO2/Si single layer, return a thickness. If that task repeats tens of thousands of times, could the input-output relation itself be learned once, so that later samples are answered directly?

This post takes up that question. A neural network that maps 300 spectral points to a single thickness is implemented from scratch in numpy, measured against LM on the same terms, and examined for where it wins, where it loses, and **why the losing side is what reveals the method's proper use.**

The short answer: on its own, the network loses to LM. Used as LM's starting value, however, it matches the success rate of the global searches in post 4 with 1/52 of the reflectance evaluations. And neither approach can fix model error or ambiguity; what has to change then is not the algorithm but the measurement (section 9).

The example is unchanged from posts 1 through 4. Normal-incidence reflectance of a single SiO2 layer is sampled at 300 wavelengths from 450 to 750 nm (true thickness 1490 nm, Gaussian noise of standard deviation 0.004), and the objective is $J(d) = \frac{1}{2}\sum_i r_i^2$. Three local minima sit in the 1200-1800 nm range, at 1293, 1490 and 1689 nm, and the basin surrounding each is 198 nm wide.

## 1. Why the landscape argument changes

One thing deserves attention first. Local minima were the difficulty throughout posts 1 to 4. Yet the neural network literature does not treat local minima as nearly so serious a problem. Why?

Because the number of unknowns differs. Posts 1 to 4 had a single unknown, the thickness $d$. The network built here has tens of thousands of weights. Dauphin and co-workers argued that the real obstacle in high-dimensional non-convex minimization is not local minima but saddle points, The first strand of their argument is a counting argument over dimension. For a critical point to be a local minimum, the surface must curve upward along **every** direction. As the dimension $N$ grows, the probability of satisfying that condition in all directions falls exponentially, and saddle points — up along some directions, down along others — come to dominate. Think of a mountain path. On a path that runs along a single line, any dip is a hollow. In a range that spreads in tens of thousands of directions, a hollow that rises on every side is far rarer than a pass that climbs one way and descends another.

Read the same argument backwards and the situation of this series appears. **At $N=1$, an exact saddle point is a probability-zero event.** In a landscape with a single unknown, every critical point is an extremum, which is precisely why three basins caused so much trouble in posts 1 to 4. Weight space is a different world. That is why training the network below never raises the question of whether it is stuck in a local minimum.

One caveat matters, though. **The ambiguity of the problem does not disappear.** That 1293 nm and 1490 nm produce similar spectra is a fact of the measurement, not of the algorithm. In a neural network that ambiguity simply wears a different face, and section 8 shows which one.

## 2. Building the regressor in numpy

This series has transcribed algorithms directly rather than calling libraries. LM was written out as Algorithm 3.16 of Madsen, Nielsen and Tingleff; here the forward pass, backpropagation and Adam are likewise transcribed from the original papers into numpy. A fully connected network needs less machinery than one might expect.

Consider the training data first, where this problem shows its peculiarity. **The forward model is free.** The Airy formula is a single line, so a spectrum can be generated for any thickness, in any quantity. This inverts the usual situation in machine learning, where collecting data is the expensive step. Eight thousand thicknesses were drawn uniformly from $[1200, 1800]$ nm and their spectra synthesized.

```python
rng = np.random.default_rng(0)
d = rng.uniform(1200.0, 1800.0, size=8000)                # labels
X = np.stack([reflectance(x, wavelength_nm) for x in d])  # inputs, (8000, 300)
```

The network is plain: 300 inputs, 256 hidden units with tanh, one linear output. Weights follow the normalized initialization of Glorot and Bengio, which for a layer with $n_\text{in}$ inputs and $n_\text{out}$ outputs draws

$$W \sim U\left[-\sqrt{\frac{6}{n_\text{in}+n_\text{out}}},\ +\sqrt{\frac{6}{n_\text{in}+n_\text{out}}}\right]$$

In other words, the width of the weights is chosen so that the variance of signals and gradients neither grows nor shrinks from layer to layer.

The optimizer is Adam. Algorithm 1 of the paper is written as pseudocode and transcribes directly: accumulate first and second moments $m$ and $v$ as exponential moving averages, correct their bias, and divide.

```python
self._t += 1
m = b1 * m + (1 - b1) * grad
v = b2 * v + (1 - b2) * grad ** 2
W -= lr * (m / (1 - b1 ** self._t)) / (np.sqrt(v / (1 - b2 ** self._t)) + eps)
```

One practical point. Input reflectances span 0.09 to 0.35 while output thicknesses span 1200 to 1800, so without normalization the network does not train at all. Inputs were standardized per wavelength channel and outputs mapped linearly to $[-1,1]$, then inverted after prediction. The full implementation is in [`_code/optimization-nn-regression/`](https://github.com/eddykim/eddykim.github.io/tree/main/_code/optimization-nn-regression).

### A bigger network is not a better one

Three architectures were trained on identical data for 600 epochs each. The table reports the error when the weights are restored to the epoch at which validation error was lowest.

| Architecture | Validation RMSE (restored to best epoch) |
| --- | ---: |
| 300 - 128 - 1 | 0.750 nm |
| **300 - 256 - 1** | **0.554 nm** |
| 300 - 256 - 128 - 1 | 0.755 nm |

Two things stand out.

First, **one hidden layer outperformed two** (0.554 against 0.755 nm). Lee and Jin reached the same conclusion when they validated neural thickness measurement against certified reference materials: configurations with more nodes and more hidden layers drifted further from the certified values, which the authors attributed to overfitting to idealized synthetic spectra.

Second, **training longer makes it worse.** The 300-256-1 network reaches its minimum validation error at epoch 217 and degrades from there; run to the full 600 epochs, its validation RMSE falls back to twice the best value (0.554 to 1.10 nm). Both are far below the ±99 nm needed to pick a basin (section 4), but where the last nanometre is at stake the difference cannot be ignored.

The remedy was written down thirty years ago. Fried and Masa, pre-evaluating spectroscopic ellipsometry data with a neural network in 1994, noted that training error keeps falling while test error begins to rise past some point, so the network must be tested periodically against unfamiliar data during training to locate the optimum. Every experiment below restores the weights from the epoch of lowest validation error.

## 3. Is the network more accurate than LM (Levenberg-Marquardt)?

The short answer is no. Feed the trained network the measured spectrum of the 1490 nm sample, then run LM on the same data.

| Method | Estimated thickness (nm) | Error (nm) | Final $J$ | Iterations |
| --- | ---: | ---: | ---: | ---: |
| Network | 1490.403 | +0.403 | 2.488e-03 | – |
| **LM ($d_0=1480$)** | **1490.122** | **+0.122** | **2.386e-03** | 4 |
| LM ($d_0=1300$) | 1293.125 | −196.875 | 8.670e-01 | 30 |

LM is more than three times as accurate. The last row is the failure already seen in post 3. With its starting point in a neighbouring basin, LM cannot leave.

Accuracy and precision must be **separated** here. Taking the same 1490 nm sample and varying only the noise realization across 50 measurements, the mean error (accuracy) and the scatter (precision) were recorded separately. Because a network also depends on its training seed and architecture, six were trained — two architectures (300-256-1 and 300-256-128-1), three seeds each — and the results reported as a range.

| Method | Error (nm) | Scatter, standard deviation (nm) |
| --- | ---: | ---: |
| Network (2 architectures × 3 seeds) | −0.98 to +0.33 | **0.10 to 2.68** |
| **LM** | **+0.001** | **0.081** |

LM has no training seed, so it occupies a single row. The network is worse in both error and scatter. The striking part, though, is that **the scatter varies by a factor of 26 across these networks** (up to 17 within one architecture, from the seed alone). Worse, the validation RMSE gives no warning of it: the network with a validation RMSE of 0.34 nm scattered by 1.72 nm, while the one at 0.53 nm scattered by 0.10 nm.

This is the point to watch in practice. **A good aggregate metric does not guarantee the reliability of any individual sample.** Validation RMSE averages over thousands of cases; what we want to know is the thickness of the one wafer in front of us. The two come apart. So when choosing a network, do not rely on validation RMSE alone; measure separately the scatter of one representative sample over repeated noise realizations.

Stopping here, the conclusion would be that neural networks are useless. But only half the picture has been seen.

## 4. Except that LM is wrong two times out of three

The comparison in section 3 carried a hidden assumption. LM was handed $d_0=1480$, **a good starting value, for free.** Without knowing that the answer is 1490, how would anyone choose 1480?

Starting values were drawn at random from $[1200, 1800]$ over 100 trials. The conditions follow the structure of the table Fried and Masa published.

| Condition | Success rate | Median \|error\| (nm) | Mean iterations |
| --- | ---: | ---: | ---: |
| (a) Network alone | 100% | 0.930 | – |
| **(b) Network → LM** | **100%** | **0.060** | 6.6 |
| (c) Random start → LM | 37% | 196.873 | 22.8 |
| (d) $d_0=1300$ → LM | 0% | 196.945 | 29.8 |

"Success" is defined as in post 4: arriving within 5 nm of the global minimum.

The picture inverts. LM from a random start **returns the wrong answer roughly two times in three.** The figure of 37% is close to the 36% measured in post 4, and corresponds roughly to the fraction of the search range (600 nm) occupied by the correct basin (198 nm). The errors on failure cluster consistently near **197 nm**, which is the signature of choosing an entire basin wrongly.

Seeded with the network's output, LM succeeds every time, and its median error of 0.060 nm is **15 times smaller than that of the network alone** (0.930 nm). The network cannot deliver an accurate answer, but it **does pick the right basin**, and picking a basin demands a precision of only ±99 nm.

The weakness found in section 3, an untrustworthy final nanometre and a scatter that swings by a factor of 26 from one network to the next, therefore does not matter here. The precision required is two orders of magnitude coarser. An unreliable instrument has been handed a job where unreliability is affordable.

That was the conclusion Fried and Masa drew. Neural pre-evaluation is neither more accurate nor more precise than ordinary multiparameter fitting, but it supplies such a good initial estimate for the subsequent regression that the time-consuming grid search becomes unnecessary. In their problem, an initial thickness wrong by as little as 50 nm was enough to land in a false minimum.

## 5. Compared with global search, what does it cost?

But post 4 already solved the basin problem. Simulated annealing and Basin-Hopping (a global search that repeats local optimization and hops between basins) never failed in post 4's 50 trials. Under the changed conditions of the 100 trials here they reach 93% and 99%. If the network also reaches 100%, nothing new has been gained.

By success rate, nothing has. **The difference lies in cost.**

Measuring fairly calls for a better yardstick than wall-clock time. The count of **forward-model evaluations** does not depend on implementation or hardware: how many times the Airy formula was called is a property of the algorithm, not of the implementation.

<img src="/assets/img/posts/optimization-nn-regression/en/fig1-success-vs-cost.png" alt="Success rate and per-sample cost" width="600">
_Fig 1. Success rate by method (bars, left axis) and forward-model calls per sample (dashed line, right log axis) over 100 random starts. Network → LM has the highest success rate and the fewest calls._

| Method | Success rate | Forward-model calls / sample | Time / sample |
| --- | ---: | ---: | ---: |
| Single LM (random start) | 37% | 119 | 1.03 ms |
| Simulated annealing | 93% | 301 | 2.96 ms |
| Basin-Hopping | 99% | 1,806 | 15.85 ms |
| **Network → LM** | **100%** | **35** | **0.34 ms** |

In Fig 1 the bars give the success rate and the dashed line the forward-model calls per sample. The bars for simulated annealing, Basin-Hopping and network → LM are similar in height, while the dashed line is lowest for network → LM.

**The network does not beat Basin-Hopping on success rate.** It does the same job at 1/52 the forward-model calls and 1/47 the time. That is the real position of the method that "lost" in section 3.

It is not free, of course. Generating 8,000 training spectra costs 8,000 forward-model calls, and Adam training costs 31 seconds. How many samples it takes to recover that one-time cost depends on which yardstick is used.

| Yardstick | vs Simulated annealing | vs Basin-Hopping |
| --- | ---: | ---: |
| Forward-model calls | 30.1 samples | **4.5 samples** |
| Wall clock (including 31 s of training) | 11,790 samples | 2,008 samples |

The two numbers diverge because the forward model in this example, one line of Airy formula, is far too cheap. Recovering 31 seconds of training through saved reflectance evaluations takes 2,000 samples. Put differently, **this example is the case least favourable to a neural network.** Where a single forward evaluation takes seconds, as in rigorous coupled-wave analysis (RCWA), or tens of milliseconds, as in transfer-matrix calculations over 200-plus layers, the call-count yardstick governs and recovery comes far sooner. Kwak and Kim make the same point in their review of semiconductor multilayer metrology, noting that machine learning computation time does not scale with the number of layers.

To check this against a problem of your own, one yardstick helps. Dividing the 31 seconds of training by the 8,000 training spectra gives about 4 ms. If a single forward evaluation costs more than that, training time is buried in the cost of generating the data. Assuming the call counts of this example carry over, recovery against Basin-Hopping takes 22 samples at 1 ms per forward evaluation, 9 at 4 ms, 6 at 10 ms and 4.5 at 1 s.

## 6. When training and measurement disagree, does the residual tell?

So far the training and measurement distributions have matched. What happens when they diverge? Three ways of breaking them were tried.

**Self-checking** was examined at the same time. Unlike LM, the network carries no internal residual signal. But the forward model is free, so the predicted thickness can be pushed back through the Airy formula and the residual recomputed. Whether that residual raises an alarm is the question.

<img src="/assets/img/posts/optimization-nn-regression/en/fig2-residual-alarm.png" alt="What recomputing the residual does and does not reveal" width="600">
_Fig 2. Residual $J$ recomputed by feeding the predicted thickness back through the Airy formula (multiple of the normal measurement, log axis), with the thickness error, for four cases. Extrapolation raises the residual 1,130-fold, model error only 4.4-fold._

| Case | Network error | Recomputed residual $J$ | Network → LM error |
| --- | ---: | ---: | ---: |
| (a) Outside training range (2000 nm) | **−382.3 nm** | **1,130×** baseline | −469.5 nm |
| (b) Noise mismatch (standard deviation 0.02) | +1.5 nm | 23× | +0.6 nm |
| (c) Model error ($n_1$=1.50) | **+40.9 nm** | **4.4×** | **+40.9 nm** |

The three cases behave entirely differently.

**(a) Extrapolation fails badly, and the residual says so loudly.** Given a 2000 nm sample from outside the training range $[1200,1800]$, the network returns 1617.7 nm. The value is not pinned to the training boundary at 1800 either; it is an arbitrary point inside the range. This is precisely what Fried and Masa warned of when they stated that a network extrapolates poorly and the training data must therefore cover the entire expected input space. Fortunately the residual rises to 1,130 times baseline, so the anomaly is immediately visible.

**(b) The noise mismatch was survived, contrary to expectation.** Five times the baseline measurement noise (0.004 to 0.02) left an error of only 1.5 nm. The author had listed this as a failure mode while preparing the post, and it turned out to do almost no harm. Noise scatters randomly across the whole spectrum and so does not systematically corrupt the fringe positions that determine thickness.

**(c) Model error deceives both methods equally.** This is the most important result. Given the spectrum of a film whose SiO2 index is actually 1.50, while both methods assume 1.46, the network is off by +40.87 nm and LM by +40.91 nm. **The two answers agree to within 0.04 nm while both are wrong by 41 nm.** The residual, meanwhile, rises only to 4.4 times baseline, hardly an alarm next to the 1,130 of the extrapolation case.

To see how small 4.4 is, look at how the residual relates to the noise. If the residual is noise alone, $J = \frac{1}{2}\sum_i r_i^2 \approx \frac{1}{2}N\sigma^2 = \frac{1}{2} \times 300 \times 0.004^2 = 0.0024$, which matches the baseline of the normal measurement (0.0025). The residual multiple, in other words, reads as a multiple of the noise variance. Fivefold noise should give 25 times; 23 was observed. A model error at 4.4 times is then indistinguishable from noise that has grown by $\sqrt{4.4} \approx 2.1$. On the log axis of Fig 2, only the extrapolation bar rises three orders of magnitude above baseline; fivefold noise and model error stay within one or two. The residual multiple alone cannot tell model error from a measurement that happens to be noisy today.

Post 3 showed Gauss-Newton collapsing and LM rescuing it when the refractive index was wrongly assumed to be 1.02. There the difference between algorithms decided the outcome. Here it does not. **When the forward model itself is wrong, every method built on that model is wrong in the same direction by the same amount.** The network inherits the faulty physics through its training data, LM through its fitting model.

The contrast "the network has no self-check, LM does" is therefore only half true. Recomputing the residual catches distribution shift but largely misses model error, and that is **a blind spot the two methods share**.

## 7. What augmentation matched to real measurements repairs

Cases (a) and (c) differ in kind. The first calls for widening the training distribution; the second calls for correcting the physics. So how far does widening the training distribution toward real measurements actually get?

The augmentation axes follow the setup Kwak and co-workers used for 3D multilayer semiconductor metrology. Unable to obtain many normal samples from a production line, they expanded 125 into 5,000. What deserves attention is that the perturbations were not arbitrary noise but **the actual error modes of a spectrometer**: a vertical offset and a lateral (wavelength-axis) offset. The same axes were carried over here.

- **Noise**: a standard deviation drawn per sample from $[0, 0.02]$ and injected
- **Vertical**: $R \to (1+s)R + c$, with scale $s$ at ±5% and baseline $c$ at ±0.04
- **Lateral**: the wavelength axis shifted by ±6 nm

<img src="/assets/img/posts/optimization-nn-regression/en/fig3-augmentation-rmse.png" alt="Errors augmentation fixes, and errors it does not" width="700">
_Fig 3. Prediction RMSE (log axis) without augmentation (grey) and with it (orange), for 8,000 (left) and 250 (right) training thicknesses. Of the five test conditions, augmentation fails to reduce the error only for model error._

The test conditions must include the error modes the augmentation targets. RMSE was computed over 200 random thicknesses per condition.

| Training | Normal | 5× noise | Model error | **Wavelength +4 nm** | **Photometric ×1.03+0.02** |
| --- | ---: | ---: | ---: | ---: | ---: |
| 8,000 thicknesses, no augmentation | 2.72 | 8.50 | 40.50 | 14.65 | 16.13 |
| 8,000 thicknesses, **augmented** | 2.53 | 6.28 | 41.30 | **2.68** | **2.56** |
| 250 thicknesses, no augmentation | 9.69 | 13.48 | 39.48 | 21.22 | 22.24 |
| 250 thicknesses, **augmented** | **2.15** | 7.08 | 41.25 | **2.39** | **2.66** |

In Fig 3, augmentation fails to reduce the error, and even raises it slightly, only for model error. The result splits three ways.

Augmentation works **only on the perturbations it covers, but there it works dramatically.** On a spectrometer whose wavelength axis is off by 4 nm, RMSE falls from 14.65 to 2.68 nm, a factor of 5.5; on a photometric offset, from 16.13 to 2.56 nm, a factor of 6.3.

The price paid on clean data is **almost nothing.** With thicknesses plentiful, 2.72 nm becomes 2.53 nm, marginally better. In the sparse regime of only 250 training thicknesses, 9.69 nm becomes 2.15 nm, a large gain; here augmentation is also **acting as regularization**. That is why Kwak and co-workers expanded 125 samples into 5,000.

On model error it does **nothing at all.** 40.50 nm becomes 41.30 nm, marginally worse. No amount of noise and offset variety repairs a forward model that is itself wrong. The conclusion matches case (c) of section 6.

In short, augmentation covers **the imperfections of the instrument**, not **the errors of the physics**. Robustness extends exactly as far as the error modes already known to the user.

## 8. Given ambiguity, the network returns an average

Now back to what section 1 left aside. What face does ambiguity wear in a neural network?

What a least-squares regressor converges to was settled by Bishop in 1994. In the limit of infinite data, minimizing squared error drives the network function to the **mean** of the target conditioned on the input.

$$f(\mathbf{x}; \mathbf{w}^*) = \langle t \mid \mathbf{x} \rangle = \int t\, p(t \mid \mathbf{x})\, dt$$

Here $t$ is the target (thickness), $\mathbf{x}$ the input spectrum, $\mathbf{w}^*$ the trained weights, and $p(t \mid \mathbf{x})$ the distribution of thicknesses that could produce that spectrum.

For classification this quantity is the posterior probability and therefore optimal, but inverse problems are another matter. When one input corresponds to several correct answers, **the average of several correct answers is not guaranteed to be correct itself.** At a fork where both the left and the right road are correct, driving down the middle puts you in a ditch. Bishop noted that in such cases the conventional least-squares approach can give completely erroneous results.

Test it on our problem by deliberately amplifying the ambiguity. Normal-incidence reflectance repeats almost identically each time the thickness changes by $\lambda / 2n_1 \approx 205$ nm, so if the wavelength band is too narrow to contain even one fringe, thicknesses separated by multiples of that period become indistinguishable. About three periods fit inside the 600 nm training range.

<img src="/assets/img/posts/optimization-nn-regression/en/fig4-conditional-mean.png" alt="What the network returns on ambiguous input" width="600">
_Fig 4. Network prediction against true thickness. Blue circles: wide band (450-750 nm, 3.9 fringes). Red squares: narrow band (595-605 nm) with fivefold noise. In the narrow band, predictions are pulled toward the centre of the training range (1500 nm)._

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

The pull toward the centre that Bishop predicts appears. At 1500 nm, the centre of the training range, prediction nearly coincides with truth (−2.7 nm), because the candidates are symmetric about the centre and their average lands on the true value. Toward the edges, however, predictions are **pulled toward the centre by as much as 194 nm.** A 1750 nm sample reads as 1574.7 nm, a 1250 nm sample as 1415.9 nm. The red curve in Fig 4 is not even monotonic: as the true thickness rises from 1550 to 1650 nm, the prediction falls from 1489.4 to 1456.2 nm.

The point worth attention is that **the value pulled to is not any basin at all.** 1574.7 nm matches none of the candidates (about 1341, 1545 and 1750 nm) that this band cannot distinguish from the 1750 nm sample; it lies between them. No thickness with that solution exists physically; only the arithmetic mean sits there. LM would have converged to one of the three. It might be wrong, but it would be a physically meaningful value.

Repeating the test on the wide band (450-750 nm, 3.9 fringes) leaves every error below 2.5 nm. **The effect appears only when the ambiguity is real.** Nothing is wrong with the network; by the definition of the loss it has no choice but to average over what the data does not distinguish.

## 9. Does collecting more of the same data resolve the ambiguity?

Invert the conclusion of section 8 and the direction of the remedy follows. Ambiguity is a property of the data, not of the algorithm. Then the data is where it must be fixed. Does collecting more of it suffice?

The simplest approach was tried first: widen the wavelength range from 450-750 nm to 450-1000 nm, holding the number of points at 300.

| Range | Local minima | LM success, random start | Network validation RMSE |
| --- | ---: | ---: | ---: |
| 450-750 nm | 3 | 28% | 0.55 nm |
| 450-1000 nm | **3** | 37% | 0.45 nm |

The 28% LM success rate is for the same 450-750 nm condition as the 37% of section 4, yet the values differ, because here a single noise realization was paired with a fresh set of 100 starting points. With 100 trials the standard error is about 5 percentage points, so neither the gap between 28% and 37% nor the change on widening the range is large enough to be statistically meaningful.

**The count of local minima did not fall.** Because the repeat period $\lambda/2n_1$ scales with the centre wavelength (205 nm at 600 nm, 248 nm at 725 nm), widening the range leaves the number of basins unchanged.

Add an unknown to see it more sharply. What happens if the refractive index $n_1$ is fitted alongside the thickness? Now the **cross-correlation** of the two parameters becomes visible: from the covariance $C = (J^\top J)^{-1}$ it is defined as $\rho = C_{01}/\sqrt{C_{00}C_{11}}$. Here $J$ is not the objective $J(d)$ but the Jacobian matrix of the residuals with respect to the two parameters. The cond column is the condition number of $J^\top J$; the larger it is, the harder the two parameters are to separate.

| Range | $\rho(d, n_1)$ | cond($J^\top J$) | $d$ scatter (nm, 30 noise trials) |
| --- | ---: | ---: | ---: |
| 450-750 nm | −0.997866 | 2.47e+08 | 0.962 |
| 450-1000 nm | −0.996694 | 1.59e+08 | 0.967 |

Widening the range moves the correlation only from −0.9979 to −0.9967, and the thickness scatter in the table above stays at 0.962 versus 0.967 nm. In the 300-point simulation below it shrinks by only about 10% (9.7%). The reason lies in the physics. At normal incidence the phase is $\beta = 2\pi n_1 d / \lambda$, so reflectance depends chiefly on the **product of $n_1$ and $d$**. However much data of the same kind is collected, none of it constrains that degenerate direction.

### Adding a different kind of measurement

A different kind is needed, then. Consider oblique incidence. At an incidence angle $\theta_0$ the refraction angle $\theta_1$ follows from Snell's law, and $\cos\theta_1$ enters the phase.

$$\beta = \frac{2\pi n_1 d \cos\theta_1}{\lambda}, \qquad n_0 \sin\theta_0 = n_1 \sin\theta_1$$

Because $\cos\theta_1$ itself depends on $n_1$, the index now enters the phase differently from the thickness. The Fresnel coefficients also split into s and p polarizations, adding amplitude information. Changing the angle, in other words, creates **new information capable of separating $n_1$ from $d$.**

Holding the number of data points fixed at 300, only the kind of measurement was varied.

<img src="/assets/img/posts/optimization-nn-regression/en/fig5-correlation-ellipse.png" alt="(d, n1) estimate scatter by kind of measurement" width="600">
_Fig 5. Scatter of (d, n1) estimates when the number of data points is fixed at 300 and only the kind of measurement changes (120 simulated measurements). Grey circles and blue triangles: normal incidence only. Orange squares: 0° and 60° combined. The points line up along the direction of constant n1·d, and mixing angles shortens their spread. The fourth configuration in the table (0°+45°+70°) is not plotted._

| Configuration (300 points total) | $\rho(d, n_1)$ | cond($J^\top J$) | $d$ scatter (nm) | $n_1$ scatter |
| --- | ---: | ---: | ---: | ---: |
| Normal incidence only, 450-750 nm | −0.997866 | 2.47e+08 | 1.303 | 0.00125 |
| Normal incidence only, 450-1000 nm | −0.996694 | 1.59e+08 | 1.176 | 0.00118 |
| **0°+60°, 150 points each** | **−0.976132** | **3.45e+07** | **0.458** | **0.00034** |
| **0°+45°+70°, 100 points each** | −0.975515 | 3.84e+07 | **0.340** | **0.00024** |

The same 300 points give different answers. Widening the wavelength range improves the thickness scatter by 9.7%, while **mixing angles improves it roughly threefold (3.1 from the linearized covariance, 3.8 in the 40-trial simulation) and the condition number by a factor of 6.4.**

The normal-incidence points in Fig 5 fall on a line with a slope of about −0.001 per nm. That line is the direction along which $n_1 d$ stays constant: keeping $n_1 d$ fixed requires $\Delta n_1 = -(n_1/d)\,\Delta d \approx -(1.46/1490)\,\Delta d = -0.00098\,\Delta d$, and lines fitted to the normal-incidence points have slopes of −0.00097 to −0.00098 per nm. The line the points follow is the degenerate direction. Mixing angles (orange) shortens the spread considerably.

This is the conclusion of section 9. **What is needed is not "enough data" but "data that sees the degenerate direction."** Kwak and Kim reach the same conclusion for semiconductor multilayer metrology. They write that single-wavelength ellipsometry does not yield a unique multilayer thickness solution, so multiangle or multisample measurements are used alongside it, and that this additional information improves the uniqueness of the characterization.

And this is not a neural network story. Johnson, treating parameter cross-correlation in curve fitting generally, observed that once the correlation is nearly complete the estimation procedure cannot identify unique parameter values or their standard errors. The correlation coefficient in his two-exponential example was −0.999933, the same kind of degeneracy as the −0.9979 seen here. More important is the qualification he adds: a correlation between two parameters does not imply that anything in the underlying mechanism is correlated, and **may simply be a consequence of estimating from a limited set of observations.**

## Closing the series

From gradient descent in post 1 to here, five posts have solved the same problem with seven methods: gradient descent, Newton's method, Gauss-Newton, LM, simulated annealing, Basin-Hopping and the network regressor. Posts 1 through 4 were about how well a given set of data can be mined: how to use the gradient ([post 1](/en/posts/optimization-gradient-descent/)), how to approximate the curvature ([post 2](/en/posts/optimization-newton-gauss-newton/)), how to tune the damping ([post 3](/en/posts/optimization-levenberg-marquardt/)), how to cross between basins ([post 4](/en/posts/optimization-global-heuristics/)).

Post 5 adds a layer on top. In this example, a learned function **could not replace** iterative optimization. It loses on accuracy, its precision swings with the training seed, and it collapses outside its training range. What it can do is **pick the basin**, matching Basin-Hopping's success rate with 1/52 of the forward-model calls. Its place as an initial-value supplier rather than a replacement was already written down thirty years ago, and these experiments confirm it again.

There is one thing none of the five can do. **Create information that is not there.** In a measurement where index and thickness enter only as a product, neither LM nor the network can separate them. LM converges wherever the starting value sends it; the network quietly returns an average. What is needed then is not a better algorithm but one more angle of incidence.

The optimization story ends here. Designing *what to measure* in the first place is a subject of its own. Why p and s light part ways at oblique incidence, and the instrument that measures it, are covered in [Ellipsometer Instrumentation 1](/en/posts/ellipsometer-oblique-rotating-element/).

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

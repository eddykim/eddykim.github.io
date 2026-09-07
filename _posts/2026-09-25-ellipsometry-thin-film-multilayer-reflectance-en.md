---
title: "Ellipsometry Foundations 3 — Multiple Reflection in Thin Films and Three Multilayer Methods"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometry-thin-film-multilayer-reflectance/
page_id: ellipsometry-thin-film-multilayer-reflectance
date: 2026-09-25 20:00:00 +0900
categories: [Optics, Thin-Film Modeling]
tags: [ellipsometry, thin-film, multilayer, transfer-matrix, scattering-matrix, interference]
description: Folding the multiple reflections inside a single film into one effective reflection coefficient, and how to choose among Rouard's method, the transfer matrix method and the scattering matrix method in practice.
math: true
---

[Post 1](/en/posts/ellipsometry-electromagnetic-fresnel/) treated reflection and refraction across a single interface. The Fresnel coefficient $r_{01}$ is a complete answer when there is exactly one interface, but a real thin-film sample has at least two — air/film and film/substrate. When the film thickness is comparable to the wavelength (tens to hundreds of nanometres), the beams reflected at those two interfaces overlap while still coherent, and the reflectance is no longer described by $r_{01}$ alone. The thickness information rides on precisely that interference pattern, which is what makes measuring film thickness by ellipsometry possible at all.

This post covers three things in order. First, the mathematical structure that folds the endlessly repeating internal reflections of a single film into one effective reflection coefficient. Second, the trade-offs among the three methods actually used when extending to a multilayer stack — Rouard's method, the Transfer Matrix Method (TMM), and the Scattering Matrix Method (SMM) — and why the situation decides which to use. Third, exactly where the resulting p and s reflection coefficients sit inside the Mueller formalism defined in [post 2](/en/posts/ellipsometry-polarization-mueller-matrix/) — the application of that formalism to a thin-film sample, promised at the end of that post.

## 1. Multiple reflection in a single film — folding a geometric series

Consider air ($N_0$), a film ($N_1$) and a substrate ($N_2$) in that order. Incident light $E_i$ reaching the upper interface (point O) is partly reflected ($E_{r1}$) and partly refracted into the film, reaching the lower interface (point A). There part of it transmits into the substrate ($E_{t1}$) and part reflects back to the upper interface (point B), where it transmits into the air ($E_{r2}$). In principle this repeats without end.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig1-concept-diagram.png" alt="Ray paths of multiple reflection in a single film and the optical path difference" width="640">
_Fig 1. Multiple reflection paths in a single film. The internal path O→A→B is longer than the reference path along the upper interface (OB), and phase accumulates by that difference._

What matters is how much extra optical path each successive reflection travels. $E_{r1}$ reflects straight off the upper interface, while $E_{r2}$ emerges only after one round trip through the film (O→A→B). Converting the excess of that round trip over the reference path — the dashed path along the upper interface in the figure — into an angle gives $\beta_{phase}$.

$$ \beta_{phase} = \frac{2\pi}{\lambda} N_1 d \cos\theta_1 $$

Here $d$ is the film thickness and $\theta_1$ the refraction angle inside the film. Each of $E_{r1}, E_{r2}, E_{r3}, \dots$ is retarded by a further $2\beta_{phase}$ per round trip and reduced in amplitude by one more factor of the interface coefficients $r_{12}$ (film to substrate) and $r_{10}$ (film to air).

$$ E_{r1} = r_{01} E_i, \qquad E_{r2} = t_{01} r_{12} t_{10}\, e^{-2j\beta_{phase}} E_i, \qquad E_{r3} = t_{01} r_{12} r_{10} r_{12} t_{10}\, e^{-4j\beta_{phase}} E_i,\ \dots $$

Summing every term gives a geometric series, and because the ratio satisfies $\lvert r_{10} r_{12} e^{-2j\beta_{phase}} \rvert < 1$ — always true physically — the infinite series converges to a closed form. Using the interface relations $r_{mn} = -r_{nm}$ and $t_{mn} t_{nm} = 1 - r_{mn}^2$, it collapses to a single effective reflection coefficient representing the whole film.

$$ r_{total} = \frac{r_{01} + r_{12}\, e^{-2j\beta_{phase}}}{1 + r_{01} r_{12}\, e^{-2j\beta_{phase}}} $$

This expression is the heart of the section. Physically the process is an unending sequence of reflections and transmissions; mathematically it folds into a closed expression in just two interface coefficients ($r_{01}, r_{12}$) and one phase ($\beta_{phase}$). The transmission coefficient folds the same way.

$$ t_{total} = \frac{t_{01} t_{12}\, e^{-j\beta_{phase}}}{1 + r_{01} r_{12}\, e^{-2j\beta_{phase}}} $$

Both carry the thickness $d$ inside $\beta_{phase}$, so measuring the reflectance $R = \lvert r_{total} \rvert^2$ while sweeping wavelength or angle of incidence records the thickness as an interference pattern. Section 8 computes how that pattern actually changes with thickness.

## 2. What goes wrong on extending to a multilayer

The expression above holds for exactly two interfaces — one film. Real samples often have several, such as deposited layers stacked on a native oxide. Each added layer multiplies the number of possible reflection paths, so enumerating every path and folding it into a geometric series as in section 1 becomes unmanageable as the count grows.

Practice therefore uses three methods, each answering differently the question of in what order and in what units to repeat the folding operation of section 1.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig2-methods-comparison.png" alt="Diagram comparing the computational structure of Rouard's method, TMM and SMM" width="720">
_Fig 2. How the three methods organize the repeated computation. Rouard folds values upward from the bottom; TMM binds interface and phase into one matrix and multiplies them in a chain; SMM keeps interfaces and layers in separate kinds of matrix._

## 3. Rouard's method — folding upward into equivalent interfaces

Rouard's method repeats the idea of section 1 recursively. Starting at the film layer closest to the substrate, everything below is bundled into an equivalent medium with an equivalent interface reflection coefficient $\rho$, and that value is updated layer by layer moving upward.

$$ \rho_{M+1} = r_{M,M+1}, \qquad \rho_m = \frac{r_{m-1,m} + \rho_{m+1}\, e^{-2j\beta_m}}{1 + r_{m-1,m}\, \rho_{m+1}\, e^{-2j\beta_m}}, \qquad r_{total} = \rho_1 $$

With $M$ layers, the recursion starts at the bottom interface coefficient $r_{M,M+1}$, proceeds upward for $m = M, M-1, \dots, 1$, and the surviving $\rho_1$ is the total reflection coefficient. The essential point is that what enters the numerator and denominator is not the bare interface coefficient $r_{m,m+1}$ of the layer below, but $\rho_{m+1}$, which has already folded in everything beneath it. Substituting $r_{m,m+1}$ instead breaks the recursion and accounts for only one layer down.

The structure is in fact identical in form to the expression of section 1 — everything from layer $m+1$ downward is treated as an "equivalent substrate" with a single layer $m$ on top, shrinking the problem step by step. As the left panel of Figure 2 shows, values are folded upward starting from the bottom interface. The implementation is in `_code/ellipsometry-thin-film-multilayer-reflectance/rouard.py`.

Because the computation updates a single scalar — the equivalent reflection coefficient — it is simple to implement and the fastest of the three for a single-wavelength, single-angle spectrum. It is optimized for reflection, though: obtaining the transmission coefficient separately is awkward, and p- and s-polarization must each be computed from scratch.

## 4. The transfer matrix method — phase and interface in one matrix

TMM takes a different approach entirely. Rather than folding the reflection coefficient through a recursion, it expresses as matrices the boundary condition that the tangential components of the electric and magnetic fields be continuous at each interface, and multiplies those matrices in order across the whole structure. The characteristic matrix $Q_m$, binding together the phase retardation through layer $m$ and that layer's index, is:

$$ Q_m = \begin{pmatrix} \cos\beta_m & \dfrac{j}{N_m}\sin\beta_m \\[4pt] jN_m\sin\beta_m & \cos\beta_m \end{pmatrix} $$

The total reflection and transmission coefficients are read directly off the matrix $M$ obtained by multiplying these in order between the two boundary matrices.

$$ M = \begin{pmatrix}1&1\\N_0&-N_0\end{pmatrix}^{-1} Q_1 Q_2 \cdots Q_m \begin{pmatrix}1&1\\N_t&-N_t\end{pmatrix}, \qquad r_{total} = \frac{M_{21}}{M_{11}}, \quad t_{total} = \frac{1}{M_{11}} $$

As the middle panel of Figure 2 shows, this is a chain: matrices multiplied end to end. Because each layer's phase change sits inside its own $Q_m$, the accumulation of phase can be followed in the order of the matrix product, which is intuitive. That intuitiveness is why TMM is widespread in optical design problems that ask how much phase to build into each layer to obtain a desired interference — anti-reflection coatings, dielectric mirrors, interference bandpass filters, Fabry-Pérot interferometers. The $N_m$ above is written for normal incidence; at oblique incidence the effective index (optical admittance) appropriate to s or p polarization is substituted and the formulation extends unchanged.

Below is the characteristic matrix implemented directly (full code in `_code/ellipsometry-thin-film-multilayer-reflectance/tmm.py`).

```python
# core of tmm.py — full code: _code/ellipsometry-thin-film-multilayer-reflectance/
def tmm_rt(n_list, d_list, wavelength, theta0=0.0, pol="s"):
    n0 = n_list[0]
    cos_list = [_cos_theta(n0, theta0, n) for n in n_list]
    eta = [_admittance(n, c, pol) for n, c in zip(n_list, cos_list)]

    D0 = np.array([[1, 1], [eta[0], -eta[0]]], dtype=complex)
    Dsub = np.array([[1, 1], [eta[-1], -eta[-1]]], dtype=complex)

    M = np.linalg.inv(D0)
    for n, d, c, e in zip(n_list[1:-1], d_list, cos_list[1:-1], eta[1:-1]):
        beta = 2 * np.pi / wavelength * n * d * c
        Q = np.array([[np.cos(beta), 1j / e * np.sin(beta)],
                      [1j * e * np.sin(beta), np.cos(beta)]])
        M = M @ Q
    M = M @ Dsub
    return M[1, 0] / M[0, 0], 1.0 / M[0, 0]   # r_total, t_total
```

## 5. The scattering matrix method — keeping interfaces and layers apart

SMM solves the same problem as TMM but keeps what happens at an interface — reflection and transmission — in a different kind of matrix from what happens inside a layer, the phase change. The interface matrix $I$, defined by the interface reflection and transmission coefficients:

$$ I_{m-1,m} = \frac{1}{t_{m-1,m}}\begin{pmatrix}1 & r_{m-1,m}\\ r_{m-1,m} & 1\end{pmatrix} $$

and the layer matrix $L$, expressing only the phase retardation inside a layer:

$$ L_m = \begin{pmatrix}e^{j\beta_m} & 0\\ 0 & e^{-j\beta_m}\end{pmatrix} $$

The whole structure is obtained by alternating the two kinds of matrix — interface, layer, interface, layer, and so on (right panel of Figure 2). The components of the resulting scattering matrix $S$ give the total reflection and transmission coefficients.

$$ \begin{pmatrix}E_0^{o+}\\E_0^{o-}\end{pmatrix} = I_{01}L_1 I_{12}L_2 \cdots L_m I_{m,m+1} \begin{pmatrix}E_{m+1}^{i+}\\E_{m+1}^{i-}\end{pmatrix} = \begin{pmatrix}S_{11}&S_{12}\\S_{21}&S_{22}\end{pmatrix}\begin{pmatrix}E_{m+1}^{i+}\\E_{m+1}^{i-}\end{pmatrix}, \qquad r_{total} = \frac{S_{21}}{S_{11}}, \quad t_{total} = \frac{1}{S_{11}} $$

In TMM the interface and phase information are mixed together inside a single $Q_m$; in SMM, $I$ and $L$ are independent blocks. That separation makes it flexible to fold in non-ideal conditions — reduced interface reflectance from surface roughness, absorption or scattering inside a film — by modifying just one kind of matrix. The cost is more matrices and more multiplications than TMM, and, as with the other two methods, separate computation per polarization.

## 6. Which method, and when

For an ideal sample all three give the same reflection coefficient. Cross-checking a single layer and a three-layer stack at normal and oblique (65°) incidence in both s and p polarization, the largest disagreement among the three is $3.6\times10^{-16}$ — floating-point noise (`verify_methods.py`). What differs is not the result but how the computation is organized, and the practical trade-offs that follow.

One trap deserves mention: the p-polarization sign convention. The $r_p$ read off a characteristic matrix built on the optical admittance $\eta_p = N/\cos\theta$ has the opposite sign to the Fresnel convention used in post 1. Reflectance $R=\lvert r\rvert^2$ is the same either way, so the difference goes unnoticed — but $\Delta$, the argument of $\rho = r_p/r_s$, is off by $180°$ the moment the conventions are mixed. All three implementations here are unified on post 1's Fresnel convention.

| Method | Speed | Transmission | Intuitiveness | Handling non-ideal conditions |
|---|---|---|---|---|
| Rouard | fastest for a single spectrum | awkward to obtain | equivalent-interface recursion, physically simple | difficult |
| TMM | reasonable | immediate | phase accumulation follows the matrix product | limited |
| SMM | relatively slow | immediate | interfaces and layers separated, easy to trace a cause | flexible (roughness, absorption applied locally) |

This table compares the situation of computing one spectrum at a single wavelength and single angle. Angle-resolved measurement, which iterates over a wavelength axis and an angle axis at once, changes the picture. Rouard's inner loop may be fast, but it still has to traverse every wavelength-angle combination, so its advantage in total operations disappears; TMM and SMM, whose repeated operation is a matrix product of identical structure, vectorize across both axes as tensors in one pass and end up faster. For reflection-based angle-resolved data collected to shorten measurement time, TMM and SMM beat Rouard on speed alone. Add the requirement to fold in non-ideal interface and film properties, and SMM becomes the final choice: the slowest of the three, but its slowness is offset by tensor operations while its flexibility is kept. Why angle resolution is needed and what that tensor computation looks like are treated directly in sections 9 and 10.

## 7. From reflection coefficients to the sample's Mueller matrix — the link to post 2

Computing $r_{total}$ independently for p- and s-polarized reflection from the film — running the same TMM or SMM calculation once per polarization with the appropriate effective index — the ratio of the two defines the ellipsometric parameters.

$$ \rho = \frac{r_p}{r_s} = \tan\Psi\, e^{j\Delta} $$

[Post 2](/en/posts/ellipsometry-polarization-mueller-matrix/) expressed the whole optical signal as a chain of Mueller matrices, $S_{out} = M_A M_C' M_{sample} M_C M_P S_{in}$. The sample's Mueller matrix $M_{sample}$, sitting at the centre of that chain, is exactly what the $r_p$ and $r_s$ just obtained determine.

$$ M_{sample} = \frac{r_p r_p^* + r_s r_s^*}{2} \begin{pmatrix} 1 & -\cos2\Psi & 0 & 0 \\ -\cos2\Psi & 1 & 0 & 0 \\ 0 & 0 & \sin2\Psi\cos\Delta & \sin2\Psi\sin\Delta \\ 0 & 0 & -\sin2\Psi\sin\Delta & \sin2\Psi\cos\Delta \end{pmatrix} $$

The leading factor $(r_p r_p^*+r_s r_s^*)/2$ sets the total intensity after reflection, the magnitude corresponding to $S_0$ in post 2. The two remaining $2\times2$ blocks split into a part mixing the linear polarization axes by $\cos2\Psi$ and a part rotating the $(S_2,S_3)$ plane by $\Delta$ — a structure resembling the polarizer and retarder Mueller matrices of section 6 of post 2. That is no coincidence: what a film does to reflected light is expressible in the same language of transmitting one axis preferentially and turning the phase. The expression above holds only for an ideal sample that preserves the polarization state entirely; a rough or multiply-scattering sample has a more general Mueller matrix with diagonal elements below 1, through the depolarization treated in section 3 of post 2.

At this point the three background posts join into one line. The Fresnel coefficients derived in post 1, extended to a multilayer by the methods of this post, give $r_p$ and $r_s$; their ratio gives $\Psi$ and $\Delta$; those build $M_{sample}$ in the formalism of post 2; and that matrix carries through the PSG-sample-PSA chain into an actual measured signal.

## 8. A computed experiment — interference fringes shifting with thickness

Section 1 claimed the thickness is written into the reflectance as an interference pattern; here it is directly. For an air/SiO$_2$/Si structure, the SiO$_2$ thickness was set to 100, 300, 600 and 1000 nm and the normal-incidence reflectance computed over 400–1000 nm with the TMM code above.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig3-thickness-fringes.png" alt="TMM reflectance spectra of SiO2/Si films, comparing interference fringes by thickness" width="640">
_Fig 3. TMM reflectance spectra of the same SiO2/Si structure with thickness alone varied. Thicker films give more fringes at narrower spacing._

At $d=100$ nm the fringe spacing is wide enough that the curve looks like a single gentle arc across this range, while $d=1000$ nm packs several narrow fringes into the same span. In $\beta_{phase} = 2\pi N_1 d\cos\theta_1/\lambda$, a larger $d$ makes $\beta_{phase}$ move further for the same small change in wavelength, so the constructive and destructive conditions are crossed more often over the same band. That is why the number and spacing of fringes alone already indicate the approximate thickness.

## 9. Why angle resolution — the need and the use

Everything so far assumed standard spectroscopic ellipsometry: one fixed angle of incidence, wavelength scanned. But looking at a single curve in Figure 3, there may be more than one $(d, N_1)$ pair producing a reflectance spectrum nearly indistinguishable from it — reduce the thickness slightly while raising the index slightly and the fringe positions barely move. This is essentially the problem seen in the [Newton/Gauss-Newton post](/en/posts/optimization-newton-gauss-newton/) and the [LM post](/en/posts/optimization-levenberg-marquardt/): too few observation equations from a wavelength axis alone cannot break the correlation between parameters, and the symptom appears in fitting as a Jacobian close to singular.

Varying the angle as well changes this. At the same wavelength, a different angle changes $\cos\theta_1$ inside $\beta_{phase}$, so each angle adds what amounts to an independent new observation. That sensitivity is also uneven across angle: for the same reason post 1 noted that p-polarized reflectance changes sharply near the Brewster angle, $\Psi$ and $\Delta$ swing strongly around certain angles for very small changes in thickness. Two models nearly indistinguishable on wavelength alone can therefore separate clearly once those sensitive angles are scanned. Whether that actually happens is checked at the end of section 10.

The obstacle is measurement time. Rotating the angle mechanically and repeating a wavelength scan at each setting multiplies the measurement time by the number of angles. The central idea of the back-focal-plane (BFP) configuration in the work this series draws on is to capture that entire angle scan in a single exposure, with no mechanical movement — at the back focal plane of an objective, light reflected at different angles of incidence lands at different spatial positions, so imaging that plane turns the angle axis into a spatial axis recorded all at once. Preserving that advantage on the computation side means computing every angle at once as a tensor, rather than one angle after another. That is the tensor computation promised in section 6.

## 10. Reflection coefficients and Mueller matrices over a spectral × angular tensor

The `tmm_rt` of section 4 is a scalar function: one wavelength, one angle, one value. Handling arrays of both could be done by wrapping it in Python loops over every combination, but the iteration count then grows as the number of angles times the number of wavelengths. Instead, keeping the SMM structure of section 5 and building the interface matrices $I$ and layer matrices $L$ as arrays carrying (angle, wavelength) axes from the start means only as many matrix products as there are layers — numpy's batched matrix multiply (`@`) broadcasts over the remaining two axes.

```python
# core of smm_tensor.py — full code: _code/ellipsometry-thin-film-multilayer-reflectance/
def smm_reflectance_tensor(n_list, d_list, theta0, wavelength, pol="s"):
    """theta0: (n_angle,), wavelength: (n_wav,) -> r_total: (n_angle, n_wav)"""
    cos_list = _cos_theta_grid(n_list, theta0)  # per-medium (n_angle, 1) cosines
    S = _interface_matrix(n_list[0], n_list[1], cos_list[0], cos_list[1], pol)
    for j in range(1, len(d_list) + 1):
        L = _layer_matrix(n_list[j], d_list[j - 1], cos_list[j], wavelength)
        I = _interface_matrix(n_list[j], n_list[j + 1], cos_list[j], cos_list[j + 1], pol)
        S = S @ L @ I          # (n_angle, n_wav, 2, 2) batched matrix product
    return S[..., 1, 0] / S[..., 0, 0]
```

`_interface_matrix` and `_layer_matrix` are simply the array versions of the interface and layer matrices defined in section 5 (included in the full code). From the two resulting tensors $r_p(\theta,\lambda)$ and $r_s(\theta,\lambda)$, computing $\rho = r_p/r_s$, $\Psi$, $\Delta$ and $M_{sample}$ from section 7 makes $M_{sample}$ itself a four-dimensional tensor carrying (angle, wavelength) axes — meaning the PSG-sample-PSA chain of post 2 must be evaluated not once per angle but simultaneously for every angle recorded in a single exposure.

Figure 4 shows $\Psi$ and $\Delta$ computed in one pass over the full grid of 0–85° and 400–1000 nm for the same SiO$_2$(300 nm)/Si structure.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig4-angle-wavelength-map.png" alt="Ellipsometric parameters Psi and Delta over angle and wavelength computed with an SMM tensor" width="760">
_Fig 4. $\Psi$ (left) and $\Delta$ (right) of a SiO2(300nm)/Si structure, computed with SMM in one pass over the full angle × wavelength grid._

At low angles $\Psi$ varies only gently with wavelength, but a bright band appears distinctly around 70–80° — near this structure's pseudo-Brewster angle, and precisely the sensitive angular region of section 9. The $\Delta$ map shows two sharp boundaries where the colour splits abruptly around that angle: the phase difference turns nearly a full cycle within that narrow angular span. Combinations of thickness and index hard to separate on a wavelength scan alone become more likely to separate once the $\Psi$ and $\Delta$ values across that sensitive band are compared.

### Comparing two models indistinguishable at a single angle

An impression that "the two-dimensional map looks richer" is not grounds for adding an angle axis. Establishing that information is genuinely gained requires picking two models indistinguishable at a single angle and comparing them along the angle axis. Assuming only a narrow band (630–670 nm) is available, the combination producing the spectrum most similar to a reference model ($d=300$ nm, $N_1=1.46$) at 65° incidence turns out to be $d=311$ nm, $N_1=1.428$ — 11 nm different in thickness and 0.03 in index.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig5-single-vs-angle-resolved.png" alt="Two film models that overlap at a single angle but separate under an angle scan" width="760">
_Fig 5. Two $(d, N_1)$ pairs indistinguishable over a narrow band. (a) At a single angle of 65°, the wavelength scans overlap; (b) at a single wavelength of 650 nm, the angle scan separates them clearly around 76°._

In (a) the two models differ by at most $0.06°$ in $\Psi$ and $0.10°$ in $\Delta$ — effectively the same curve. With this band and this angle alone there is no way to tell them apart. Fixing the wavelength at 650 nm and sweeping the angle instead, (b) shows $\Delta$ separating by more than $19°$ around 76°. Over that narrow span $\Delta$ drops steeply from near $150°$ to near $0°$, and the angle at which the drop occurs differs between the two models. A difference small enough to be buried in measurement noise on the gentle part of the curve is amplified by two orders of magnitude through this steep transition.

This advantage is conditional, however. Comparing the same two models over the wide 400–1000 nm band, the maximum difference reaches $9.36°$ even at a fixed angle of 65° — the $0.10°$ of the narrow band grows 93-fold from widening the band alone. A wide wavelength axis already contains several interference fringes, separating thickness from index before any angle axis is added. Angle resolution is decisive for information content when the wavelength axis is limited — a narrow band, monochromatic light, or many layers and therefore many unknowns. Under other conditions its value lies less in information than in the measurement time discussed in section 9.

## 11. What comes next — fitting this model to data

Everything so far has been the forward problem: assume thickness $d$ and index $N_1$ are known and compute reflectance. Real metrology needs the reverse — find the $d$ and $N_1$ best matching a measured reflectance (or $\Psi, \Delta$) spectrum. Fitting the nonlinear model $r_{total}(d, N_1, \theta, \lambda)$ to observations is exactly the nonlinear least-squares problem of the [Levenberg-Marquardt post](/en/posts/optimization-levenberg-marquardt/). The thickness fitted there is an input parameter of the $r_{total}$ model derived here, and updating $d$ and $N_1$ to reduce the residual between the TMM- or SMM-computed reflectance and the measurement is the core loop of real thin-film metrology software. When the observations form the full (angle, wavelength) grid as in section 10, the residual vector grows accordingly and the Jacobian carries that much more information, which helps break the correlation between parameters.

## Summary

This third post treated the multiple-reflection interference that arises once two or more interfaces overlap. Three points stand out. First, the mathematical structure by which even endlessly repeating physical reflections fold, as a geometric series, into a closed expression in a few interface coefficients and one phase. Second, that extending this folding to a multilayer yields three different strategies — Rouard, TMM, SMM — and that none is simply "best": which is favourable depends on the situation (single spectrum or angle-resolved, ideal sample or not), an engineering fact rather than a mathematical one. Third, that when the wavelength axis carries limited information — a narrow band, or many unknowns — scanning the angle as well separates combinations that wavelength alone leaves degenerate. Conversely, where a wide band is available, the added value of the angle axis lies in measurement time rather than information. Either way, realizing that value means treating the reflection-coefficient computation itself as a tensor over both angle and wavelength.

This concludes the "Ellipsometry Foundations" series. Post 1 established why measuring intensity alone loses the phase, and why a complex Fresnel reflection coefficient is therefore needed. Post 2 set out why that complex information must be handled with Stokes vectors and Mueller matrices rather than Jones vectors, and mapped the polarizer and retarder matrices onto real hardware. This third post connected how those reflection coefficients turn into an interference signal carrying thickness inside a film structure, and how the resulting $r_p$ and $r_s$ enter the $M_{sample}$ of post 2. The numerical side — fitting this reflection model to measured data — is already carried by the optimization series.

## References

- Youngjoon Kim, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발" [Development of snapshot angle-resolved ellipsometry using a line-scan spectrometer and back-focal-plane spectral interference], Ph.D. dissertation, Seoul National University, 2025 (in Korean), sections 2.4–2.5.
- O.S. Heavens, *Optical Properties of Thin Solid Films*, Dover, 1991 (the standard derivation of the characteristic matrix method).
- M. Born and E. Wolf, *Principles of Optics*, 7th ed., Cambridge University Press, 1999, ch. 1 (matrix formulation of multilayer reflection and transmission).
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 2, 5 (thin-film modelling and fitting of ellipsometric data).

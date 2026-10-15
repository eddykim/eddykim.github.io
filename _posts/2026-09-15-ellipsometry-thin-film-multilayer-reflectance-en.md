---
title: "Ellipsometry Foundations 3 — How Is Light Reflected Many Times Inside a Thin Film Added Up into One Value?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometry-thin-film-multilayer-reflectance/
page_id: ellipsometry-thin-film-multilayer-reflectance
date: 2026-09-15 20:00:00 +0900
categories: [Optics, Thin-Film Modeling]
tags: [ellipsometry, thin-film, multilayer, transfer-matrix, scattering-matrix, interference]
description: Folding the multiple reflections inside a single film into one effective reflection coefficient, and how to choose among Rouard's method, the transfer matrix method and the scattering matrix method in practice.
math: true
---

[Post 1](/en/posts/ellipsometry-electromagnetic-fresnel/) treated reflection and refraction at a single interface. With only one interface, the Fresnel coefficient $r_{01}$ is the whole answer. A real thin-film sample, however, has at least two: air/film and film/substrate. When the film is about as thick as a wavelength, tens to hundreds of nanometres, the beams reflected at the two interfaces overlap while still coherent. The reflectance is then no longer described by $r_{01}$ alone. The thickness information rides on exactly that interference, and it is the reason ellipsometry can measure film thickness at all.

This post covers three things in order. First, the mathematical structure that folds the endlessly repeated reflections inside a single film into one effective reflection coefficient. Second, the three methods used in practice for multilayer stacks: Rouard's method, the Transfer Matrix Method (TMM) and the Scattering Matrix Method (SMM). The point is what each one gains and what it gives up, and therefore why the situation decides which to use. Third, where the resulting p and s reflection coefficients sit in the Mueller formalism of [post 2](/en/posts/ellipsometry-polarization-mueller-matrix/). This is the step promised at the end of post 2: applying the Mueller formalism to reflection from a film sample.

## 1. How do we add up reflections that repeat without end inside a film?

Consider air ($N_0$), a film ($N_1$) and a substrate ($N_2$) in that order. Incident light $E_i$ reaching the upper interface (point O) is partly reflected ($E_{r1}$) and partly refracted into the film, reaching the lower interface (point A). There part of it transmits into the substrate ($E_{t1}$) and part reflects back to the upper interface (point B), where it transmits into the air ($E_{r2}$). In principle this repeats without end. How can infinitely many reflected beams be added into a single value?

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig1-concept-diagram.png" alt="Ray paths of multiple reflection in a single film and the optical path difference" width="640">
_Fig 1. Multiple reflection paths in a single film. The internal optical path O→A→B is longer than the reference path, and a phase difference $$2\beta_{phase}$$ accumulates by that excess. The reference path is the distance $$\overline{OB}\sin\theta_0$$ that the first reflection travels in air to the wavefront through B._

What matters is how much extra optical path each successive reflection travels. $E_{r1}$ reflects straight off the upper interface, while $E_{r2}$ emerges only after one round trip through the film (O→A→B). Converting the excess of that round trip over the reference path into an angle gives $2\beta_{phase}$. The reference path is the distance $\overline{OB}\sin\theta_0$ that the first reflection $E_{r1}$ travels in air to the wavefront through B (dashed in the figure). $\beta_{phase}$ is half of that, the phase thickness for a single crossing of the film.

$$ \beta_{phase} = \frac{2\pi}{\lambda} N_1 d \cos\theta_1 $$

Here $d$ is the film thickness and $\theta_1$ the refraction angle inside the film. Each of $E_{r1}, E_{r2}, E_{r3}, \dots$ is delayed by a further $2\beta_{phase}$ per round trip. Its amplitude also shrinks by one more factor of the interface coefficients $r_{12}$ (film to substrate) and $r_{10}$ (film to air).

$$ E_{r1} = r_{01} E_i, \qquad E_{r2} = t_{01} r_{12} t_{10}\, e^{-2j\beta_{phase}} E_i, \qquad E_{r3} = t_{01} r_{12} r_{10} r_{12} t_{10}\, e^{-4j\beta_{phase}} E_i,\ \dots $$

Summing every term gives a geometric series. Away from grazing incidence the magnitude of the ratio, $\lvert r_{10} r_{12} e^{-2j\beta_{phase}} \rvert$, is below 1, so the infinite series converges to a closed form. The interface coefficients obey $r_{mn} = -r_{nm}$ and $t_{mn} t_{nm} = 1 - r_{mn}^2$. Using these, the sum collapses to a single effective reflection coefficient for the whole film.

$$ r_{total} = \frac{r_{01} + r_{12}\, e^{-2j\beta_{phase}}}{1 + r_{01} r_{12}\, e^{-2j\beta_{phase}}} $$

This expression is the heart of the section. Physically, reflection and transmission repeat without end. Mathematically they fold into a geometric series that leaves only two interface coefficients ($r_{01}, r_{12}$) and one phase ($\beta_{phase}$).

Two limits can be read off its form. As the thickness goes to zero, $e^{-2j\beta_{phase}}\to1$ and $r_{total}$ becomes $(r_{01}+r_{12})/(1+r_{01}r_{12})$. That is none other than $r_{02}$, the coefficient of a bare air/substrate interface. With no film only the substrate remains, as it should. The other end is more interesting. For a non-absorbing film, whenever $2\beta_{phase}$ is a whole multiple of $2\pi$, $e^{-2j\beta_{phase}}$ returns to 1 and $r_{total} = r_{02}$ even though the film is there. So at one fixed wavelength and angle, the reflection coefficient repeats every time the thickness changes by $\lambda/(2N_1\cos\theta_1)$. At 633 nm, 70° and $N_1=1.46$ that period is 283.2 nm. Under these conditions, 100 nm and 383.2 nm of SiO$_2$ on Si ($3.88-0.02j$) give exactly the same $(\Psi, \Delta)$ of post 1: (41.21°, 79.53°). In other words, one wavelength cannot tell which period the thickness lies in. That is the first reason for measuring a spectrum.

The transmission coefficient folds the same way.

$$ t_{total} = \frac{t_{01} t_{12}\, e^{-j\beta_{phase}}}{1 + r_{01} r_{12}\, e^{-2j\beta_{phase}}} $$

Both $r_{total}$ and $t_{total}$ carry the thickness $d$ inside $\beta_{phase}$. Measuring the reflectance $R = \lvert r_{total} \rvert^2$ while sweeping wavelength or angle therefore records the thickness as interference fringes. Section 8 computes how those fringes change with thickness.

The same fringes are visible in everyday life. A soap film or an oil film on water shows different colours from place to place because its thickness varies. A different thickness means a different phase thickness $\beta_{phase}$, and so a different wavelength interferes constructively. The colour of the film is a clue to its thickness. Ellipsometry reads the same interference wavelength by wavelength, and reads the phase as well. Figure 2 shows the same effect at a single wavelength.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/ext-gasoline-film-589nm.jpg" alt="Bright and dark interference bands of a gasoline film on water under 589 nm light" width="640">
_Fig 2. A gasoline film on water lit by a 589 nm laser. Each bright or dark band is a contour of equal thickness, and between neighbouring bands the thickness changes by $$\lambda/(2N_1)$$ when viewed near normal incidence. (Source: Zaereth, [Thin film interference gasoline on water @ 589nm](https://commons.wikimedia.org/wiki/File:Thin_film_interference_gasoline_on_water_%40_589nm.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

## 2. What becomes hard once there are several layers?

The expression above holds for exactly two interfaces, one film. Real samples often have several, such as deposited layers stacked on a native oxide. Each added layer multiplies the number of reflection paths. Listing every path and folding it into a geometric series, as in section 1, soon becomes unmanageable.

Practice uses three methods instead. Each gives a different answer to the question of in what order, and in what units, to repeat the folding of section 1. Figure 3 shows the difference at a glance. Rouard's method (left) folds values upward from the bottom, TMM (centre) multiplies matrices in one chain, and SMM (right) alternates two kinds of matrix.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig2-methods-comparison.png" alt="Diagram comparing the computational structure of Rouard's method, TMM and SMM" width="720">
_Fig 3. How the three methods organize the repeated computation. Rouard's method starts from the reflection coefficient at the substrate and folds the equivalent coefficient $\rho$ upward. TMM binds interface and phase into one matrix $$Q_m$$ and multiplies them in order; SMM keeps interface matrices ($I$) and layer matrices ($L$) apart and alternates them._

## 3. How does Rouard's method fold upward from the bottom layer?

Could the expression of section 1 simply be reused layer by layer? Rouard's method does exactly that, recursively. First the notation. Layers are numbered $m = 1, 2, \dots, M$ from the top; air is medium 0 and the substrate is medium $M+1$. $\beta_m = 2\pi N_m d_m\cos\theta_m/\lambda$ is the phase thickness of layer $m$. Starting at the layer closest to the substrate, everything below is bundled into one equivalent medium with an equivalent reflection coefficient $\rho$. That value is then updated layer by layer moving upward.

$$ \rho_{M+1} = r_{M,M+1}, \qquad \rho_m = \frac{r_{m-1,m} + \rho_{m+1}\, e^{-2j\beta_m}}{1 + r_{m-1,m}\, \rho_{m+1}\, e^{-2j\beta_m}}, \qquad r_{total} = \rho_1 $$

The recursion starts at the bottom interface coefficient $r_{M,M+1}$, proceeds upward for $m = M, M-1, \dots, 1$, and the surviving $\rho_1$ is the total reflection coefficient. The essential point is what enters the numerator and denominator. It is not the bare coefficient $r_{m,m+1}$ of the interface below, but $\rho_{m+1}$, which has already folded in everything beneath. Substituting $r_{m,m+1}$ breaks the recursion and accounts for only one layer down.

The recursion has the same form as the expression of section 1. Everything from layer $m+1$ down is treated as an "equivalent substrate" with a single layer $m$ on top, and the problem shrinks step by step. The left panel of Figure 3 shows this. The implementation is in `_code/ellipsometry-thin-film-multilayer-reflectance/rouard.py`.

Because it updates a single scalar, the equivalent reflection coefficient, it is simple to implement. For one spectrum at one angle it is the fastest of the three. It is organized around reflection, though: the transmission coefficient needs a second recursion alongside the first (Fujiwara, eqs. 2.91 and 2.93). p- and s-polarization are each computed from scratch.

## 4. What does the transfer matrix method (TMM) multiply?

What changes if, instead of folding the reflection coefficient through a recursion, the boundary conditions themselves are written as matrices? TMM expresses as a matrix the condition that the tangential electric and magnetic fields be continuous at each interface, and multiplies those matrices in order across the whole structure. The characteristic matrix $Q_m$ combines the phase delay through layer $m$ with that layer's index:

$$ Q_m = \begin{pmatrix} \cos\beta_m & \dfrac{j}{N_m}\sin\beta_m \\[4pt] jN_m\sin\beta_m & \cos\beta_m \end{pmatrix} $$

The total reflection and transmission coefficients are read directly off the matrix $T$ obtained by multiplying these in order between the two boundary matrices. $N_t$ is the substrate index.

$$ T = \begin{pmatrix}1&1\\N_0&-N_0\end{pmatrix}^{-1} Q_1 Q_2 \cdots Q_M \begin{pmatrix}1&1\\N_t&-N_t\end{pmatrix}, \qquad r_{total} = \frac{T_{21}}{T_{11}}, \quad t_{total} = \frac{1}{T_{11}} $$

As the centre panel of Figure 3 shows, this is a chain of matrices multiplied end to end. Each layer's phase change sits inside its own $Q_m$, so the build-up of phase can be followed in the order of the product. That is why TMM dominates optical design problems that ask how much phase each layer needs for a desired interference: anti-reflection coatings, dielectric mirrors, interference bandpass filters, Fabry-Pérot interferometers. The $N_m$ above is written for normal incidence. At oblique incidence the effective index for s or p, the optical admittance, takes its place and the form carries over unchanged. For p-polarization, however, the $r_p$ read off this way has the opposite sign to post 1's convention and must be flipped once (section 6).

Figure 4 shows an anti-reflection coating of this kind.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/ext-antireflection-coating.jpg" alt="A window with an anti-reflection coating" width="640">
_Fig 4. A glass window with an anti-reflection coating. Its reflection changes with the viewing angle, because a larger angle shrinks $\cos\theta$ and moves each layer's phase thickness away from its design value. (Source: Zaereth, [Antireflection coating split pic](https://commons.wikimedia.org/wiki/File:Antireflection_coating_split_pic.jpg), [CC0](https://creativecommons.org/publicdomain/zero/1.0/deed.en))_

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
    r_total, t_total = M[1, 0] / M[0, 0], 1.0 / M[0, 0]
    if pol == "p":   # admittance convention -> post 1's Fresnel convention (section 6)
        r_total = -r_total
    return r_total, t_total
```

The variable `M` in the code is the matrix $T$ of the text.

## 5. Why does the scattering matrix method (SMM) separate interfaces from layers?

In TMM, interface information and phase information are mixed inside a single $Q_m$. What is gained by pulling them apart? SMM solves the same problem as TMM, but puts what happens at an interface (reflection and transmission) and what happens inside a layer (the phase change) into different kinds of matrix. There is the interface matrix $I$, defined by the interface reflection and transmission coefficients,

$$ I_{m-1,m} = \frac{1}{t_{m-1,m}}\begin{pmatrix}1 & r_{m-1,m}\\ r_{m-1,m} & 1\end{pmatrix} $$

and, separately, the layer matrix $L$, which carries only the phase delay inside a layer.

$$ L_m = \begin{pmatrix}e^{j\beta_m} & 0\\ 0 & e^{-j\beta_m}\end{pmatrix} $$

The whole structure is obtained by alternating the two kinds: interface, layer, interface, layer (right panel of Figure 3). The components of the resulting scattering matrix $S$ give the total reflection and transmission coefficients.

$$ \begin{pmatrix}E_0^{o+}\\E_0^{o-}\end{pmatrix} = I_{01}L_1 I_{12}L_2 \cdots L_M I_{M,M+1} \begin{pmatrix}E_{M+1}^{i+}\\E_{M+1}^{i-}\end{pmatrix} = \begin{pmatrix}S_{11}&S_{12}\\S_{21}&S_{22}\end{pmatrix}\begin{pmatrix}E_{M+1}^{i+}\\E_{M+1}^{i-}\end{pmatrix}, \qquad r_{total} = \frac{S_{21}}{S_{11}}, \quad t_{total} = \frac{1}{S_{11}} $$

The name follows Azzam and Bashara. This scattering matrix links the fields on the two sides of a layer. It is not the S-matrix of computational electromagnetics, which links incoming and outgoing waves.

Because $I$ and $L$ are independent blocks, non-ideal conditions can be added by changing only one of them. Reduced interface reflection from surface roughness goes into the interface matrix, and absorption or scattering inside a film goes into the layer matrix. The cost is more matrices and more multiplications than TMM. As with the other two methods, each polarization is computed separately.

## 6. Which method should be used, and when?

The conclusion first. For the spectrum of an ideal sample, Rouard's method is enough. When the transmission coefficient is needed too, or a coating is being designed, TMM is the convenient choice. When roughness or absorption must be added to a single interface or layer, use SMM.

The grounds are as follows. For an ideal sample all three give the same reflection coefficient. Cross-checking a single layer and a three-layer stack at normal and oblique (65°) incidence, in both s and p, the largest disagreement among the three is $3.6\times10^{-16}$, which is floating-point noise (`verify_methods.py`). What differs is not the result but how the computation is organized, and the practical gains and costs that follow.

One trap is the p-polarization sign convention. The $r_p$ read off a characteristic matrix built on the optical admittance $\eta_p = N/\cos\theta$ has the opposite sign to the Fresnel convention of post 1. Reflectance $R=\lvert r\rvert^2$ is the same either way, so the difference goes unnoticed. But $\Delta$, the argument of $\rho = r_p/r_s$, is off by $180°$ the moment the conventions are mixed. All three implementations here use post 1's Fresnel convention.

| Method | Speed | Transmission | Intuitiveness | Handling non-ideal conditions |
|---|---|---|---|---|
| Rouard | fastest for one spectrum at one angle | needs a second recursion | equivalent-interface recursion, physically simple | difficult |
| TMM | reasonable | immediate | phase build-up follows the matrix product | limited |
| SMM | relatively slow | immediate | interfaces and layers separated, easy to trace a cause | flexible (roughness, absorption applied locally) |

The table assumes one spectrum computed at one angle. Does the picture change when the wavelength and angle axes must be computed together, as in angle-resolved measurement? Then what decides the speed is not the method but the implementation: whether the (angle, wavelength) combinations are visited one by one in a Python loop or handled at once as array operations. On that point all three methods are in the same position.

Rouard's recursion is an independent scalar operation at each grid point, so it vectorizes just as directly. Done that way it needs the fewest operations and remains the fastest; with the code of this post, Rouard is 4–5 times faster than SMM on a 300×300 grid. The reason to choose SMM for angle-resolved data is therefore its flexibility: non-ideal interface and film properties can be added by changing one kind of matrix. Tensor operations merely keep the extra cost manageable. Nobody chooses SMM for speed. Why angle resolution is needed, and what the tensor computation looks like, are the subject of sections 9 and 10.

## 7. Where do the reflection coefficients enter the sample's Mueller matrix?

What we have so far are reflection coefficients. Where do they go in the Mueller chain of post 2? Computing $r_{total}$ separately for p- and s-polarized reflection from the film, by running the same TMM or SMM calculation once per polarization with the matching effective index, the ratio of the two defines the ellipsometric parameters.

$$ \rho = \frac{r_p}{r_s} = \tan\Psi\, e^{j\Delta} $$

Here $\rho$ without a subscript means this ratio, not the equivalent coefficient $\rho_m$ of section 3. [Post 2](/en/posts/ellipsometry-polarization-mueller-matrix/) expressed the whole optical signal as a chain of component Mueller matrices, $$S_{out} = M_A M_C' M_{sample} M_C M_P S_{in}$$. The sample's Mueller matrix $M_{sample}$, at the centre of that chain, is fixed by the $r_p$ and $r_s$ just obtained.

$$ M_{sample} = \frac{r_p r_p^* + r_s r_s^*}{2} \begin{pmatrix} 1 & -\cos2\Psi & 0 & 0 \\ -\cos2\Psi & 1 & 0 & 0 \\ 0 & 0 & \sin2\Psi\cos\Delta & \sin2\Psi\sin\Delta \\ 0 & 0 & -\sin2\Psi\sin\Delta & \sin2\Psi\cos\Delta \end{pmatrix} $$

The leading factor $$(r_p r_p^*+r_s r_s^*)/2$$ sets the total intensity after reflection, the quantity corresponding to $S_0$ in post 2. The Stokes components follow post 2's definition with $x$ taken as p and $y$ as s. Of the two remaining $2\times2$ blocks, the upper one mixes the linear polarization axes by $\cos2\Psi$ and the lower one rotates the $(S_2,S_3)$ plane by $\Delta$. The structure resembles the polarizer and retarder Mueller matrices of section 6 of post 2. That is no coincidence. What a film does to reflected light can be described in the same terms: pass one axis more strongly and turn the phase.

The expression holds only for an ideal sample that fully preserves the polarization state. A rough surface or multiple scattering pulls a sample away from this form through the depolarization discussed in section 3 of post 2. For an isotropic sample the form of the matrix stays the same, but the relation $N^2+C^2+S^2=1$ among $N=\cos2\Psi$, $C=\sin2\Psi\cos\Delta$ and $S=\sin2\Psi\sin\Delta$ breaks: the sum drops below 1, and so does the degree of polarization of the reflected light.

Here the three background posts join into one line. The Fresnel coefficients of post 1, extended to a multilayer by the methods of this post, give $r_p$ and $r_s$. Their ratio gives $\Psi$ and $\Delta$, from which $M_{sample}$ is built in the formalism of post 2. That matrix then carries through post 2's chain of polarization state generator (PSG), sample and polarization state analyzer (PSA) into an actual measured signal.

## 8. How do the interference fringes change with thickness?

Section 1 said the thickness is written into the reflectance as interference fringes. Here that is checked directly. For an air/SiO$_2$/Si structure, the SiO$_2$ thickness was set to 100, 300, 600 and 1000 nm and the normal-incidence reflectance computed over 400–1000 nm with the TMM code above. Refractive indices are held constant at their values near 633 nm (SiO$_2$ 1.46, Si $3.88-0.02j$).

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig3-thickness-fringes.png" alt="TMM reflectance spectra of SiO2/Si films, comparing interference fringes by thickness" width="640">
_Fig 5. Normal-incidence reflectance (TMM) of air/SiO2/Si with only the SiO2 thickness varied. Indices are fixed at their 633 nm values. Thicker films give denser fringes; at $d=1000$ nm five maxima fall within 400–1000 nm._

At $d=100$ nm the fringes are so widely spaced that the curve looks like one gentle arc across the band. At $d=1000$ nm several narrow fringes crowd into the same span. In $\beta_{phase} = 2\pi N_1 d\cos\theta_1/\lambda$, a larger $d$ makes $\beta_{phase}$ move further for the same small change in wavelength. The constructive and destructive conditions are therefore crossed more often over the same band.

The fringe spacing also gives a quick thickness estimate. At normal incidence, if two neighbouring maxima sit at $\lambda_1<\lambda_2$, their round-trip phases $2\beta_{phase}$ differ by exactly $2\pi$, so

$$ d \approx \frac{\lambda_1\lambda_2}{2N_1(\lambda_2-\lambda_1)} $$

The $d=1000$ nm curve of Figure 5 has maxima at 417.2, 486.7, 584.1, 730.1 and 973.6 nm, and every neighbouring pair returns 1000 nm. For $d=300$ nm the two maxima at 438.2 and 876.6 nm return 300 nm. At $d=100$ nm, by contrast, the band holds no maximum at all, only one minimum near 585 nm, and the formula cannot be used. The estimate works only when at least two maxima are visible, which in this band means from about 300 nm up. Thinner films need the fitting of section 11.

Thinner still, the reflectance becomes almost blind. A 1 nm SiO$_2$ layer on Si has a normal-incidence reflectance at 633 nm of 0.34825, against 0.34830 for bare Si: a relative difference of only 0.015%. Seen at 633 nm and 70°, the same 1 nm layer moves $\Delta$ from 179.19° to 176.32°, a shift of 2.87°, while $\Psi$ moves by only 0.014°. A film of a few ångströms is visible in $\Delta$ and nearly invisible in reflectance. That is why ellipsometry measures phase.

## 9. Why measure at more than one angle of incidence?

All the computations so far assumed standard spectroscopic ellipsometry: one fixed angle of incidence, wavelength scanned. Yet for any single curve in Figure 5, there may be several $(d, N_1)$ pairs that produce a nearly indistinguishable reflectance spectrum. Reduce the thickness slightly while raising the index slightly, and the fringe positions barely move. When thickness and index are fitted together, this correlation appears as two nearly parallel columns of the Jacobian, and the normal equations approach a singular matrix. The [LM post](/en/posts/optimization-levenberg-marquardt/) showed the normal equations approaching singularity because the Jacobian itself went to zero. That post fitted thickness alone, so the cause differs, but the symptom, an unstable fit, is the same.

Varying the angle as well changes this. At the same wavelength, a different angle changes $\cos\theta_1$ inside $\beta_{phase}$, so each angle adds what amounts to an independent observation. The sensitivity is also uneven across angle. Post 1 noted that the phase difference $\Delta$ changes sharply near the Brewster angle. With the Si of this post ($3.88-0.02j$ at 633 nm), the bare substrate's $\Delta$ falls from 179.19° at 70° through 88.2° at 75.55° to 0.77° at 80°. For the same reason, around certain angles $\Psi$ and $\Delta$ swing strongly for very small changes in thickness. Two models nearly indistinguishable on wavelength alone can therefore separate clearly once those sensitive angles are scanned. The end of section 10 checks whether that actually happens.

The obstacle is measurement time. Rotating the angle mechanically and repeating a wavelength scan at each setting multiplies the measurement time by the number of angles. The back-focal-plane (BFP) configuration adopted in the author's doctoral dissertation captures that whole angle scan in a single exposure, with no moving parts. At the back focal plane of an objective, light reflected at different angles of incidence lands at different positions, so imaging that plane turns the angle axis into a spatial axis recorded all at once. To keep that advantage in the computation, every angle must be computed at once as a tensor rather than one angle after another. This is the tensor computation announced in section 6.

## 10. How can every angle and wavelength be computed at once?

The `tmm_rt` of section 4 is a scalar function: one wavelength, one angle, one value. How should arrays of both be handled? Wrapping the function in loops over every combination works, but the iteration count grows as the number of angles times the number of wavelengths. Instead, keep the SMM structure of section 5 and build the interface matrices $I$ and layer matrices $L$ as arrays with (angle, wavelength) axes from the start. Then only as many matrix products as there are layers are needed, because numpy's batched matrix multiply (`@`) broadcasts over the remaining two axes.

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

`_interface_matrix` and `_layer_matrix` are the array versions of the interface and layer matrices of section 5 (included in the full code). From the two resulting tensors $r_p(\theta,\lambda)$ and $r_s(\theta,\lambda)$, computing $\rho$, $\Psi$, $\Delta$ and $M_{sample}$ as in section 7 makes $M_{sample}$ itself a four-dimensional tensor with (angle, wavelength) axes. The PSG-sample-PSA chain of post 2 must then be evaluated not once per angle but for every angle recorded in a single exposure at the same time.

Figure 6 shows $\Psi$ and $\Delta$ computed in one pass over the full grid of 0–85° and 400–1000 nm for the same SiO$_2$(300 nm)/Si structure.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig4-angle-wavelength-map.png" alt="Ellipsometric parameters Psi and Delta over angle and wavelength computed with an SMM tensor" width="760">
_Fig 6. $\Psi$ (left) and $\Delta$ (right) of a SiO2(300nm)/Si structure, computed with SMM in one pass over the full angle × wavelength grid. The maximum on the bright band at 446 nm and 70.5° is where $$r_s\approx0$$, the dark spot at 657 nm and 75.6° is where $$r_p\approx0$$, and in the $\Delta$ map the colours wind once around each spot._

At low angles $\Psi$ varies only gently with wavelength, but around 70–80° two spots stand out. The maximum on the bright band near 446 nm and 70.5° ($\Psi \approx 89.5°$) is where interference drives $r_s$ almost to zero. The dark spot near 657 nm and 75.6° ($\Psi \approx 0.3°$) is where $r_p$ almost vanishes, the pseudo-Brewster angle. At 657 nm and 75.6° the round-trip phase $2\beta_{phase}$ is almost exactly $2\pi$. As section 1 showed, the film then behaves as if absent, and the pseudo-Brewster angle of the Si substrate itself (about 75.6°) shows through. This is the sensitive angular region of section 9.

In the $\Delta$ map the colours wind once around each of the two spots. Each is a phase singularity: going once around it changes $\Delta$ by 360°. Fixing the wavelength and sweeping the angle close to one therefore flips $\Delta$ by nearly 150° within a few degrees (Figure 7(b)). Thickness and index combinations that a wavelength scan can hardly tell apart are likely to separate once the $\Psi$ and $\Delta$ values in this region are compared too.

### Do two models that look alike at one angle separate in an angle scan?

An impression that the two-dimensional map looks richer is not grounds for adding an angle axis. To establish that information is really gained, pick two models that one angle cannot distinguish and compare them along the angle axis. Suppose only a narrow band (630–670 nm) is available. The reference model is $d=300$ nm, $N_1=1.46$. Lowering the index to 1.428 and searching for the thickness whose spectrum at 65° is closest to the reference gives 311 nm. The pair differs by 11 nm in thickness and 0.03 in index.

<img src="/assets/img/posts/ellipsometry-thin-film-multilayer-reflectance/en/fig5-single-vs-angle-resolved.png" alt="Two film models that overlap at a single angle but separate under an angle scan" width="760">
_Fig 7. Two $$(d, N_1)$$ pairs indistinguishable over a narrow band. (a) In a wavelength scan at a single angle of 65°, the curves overlap. (b) In an angle scan at a single wavelength of 650 nm, $\Delta$ separates clearly around 76°. The dotted lines in (b) mark the angle of the largest difference between the two models, separately for $\Psi$ and $\Delta$, and the $\Psi$ axis of (a) is magnified._

In (a) the two models differ by at most $0.06°$ in $\Psi$ and $0.10°$ in $\Delta$, effectively the same curve. If noise masks differences of this size, there is no way to tell them apart with this band and this angle alone. Fixing the wavelength at 650 nm and sweeping the angle instead, (b) shows $\Delta$ separating by more than $19°$ around 76°. Over that narrow span $\Delta$ drops steeply from near $150°$ to near $0°$, and the angle of the drop differs between the two models. A difference small enough to be lost in measurement noise on the gentle part of the curve grows by two orders of magnitude across this steep transition.

The advantage is conditional. Comparing the same two models over the wide 400–1000 nm band, the maximum difference reaches $9.36°$ even at a fixed 65°. The $0.10°$ of the narrow band grows 93-fold from widening the band alone. A wide wavelength axis already contains several interference fringes, which separate thickness from index before any angle axis is added. For this one pair of thickness and index, angle resolution was decisive for information content when the wavelength axis was limited: a narrow band or monochromatic light. It may also matter with many layers and therefore many unknowns, but this post did not calculate that case. Nor did it compare the gain from adding an angle axis to a wide band, so whether the value of angle resolution there lies less in information than in the measurement time of section 9 is left open.

## 11. How do we work back from measurements to thickness?

Everything so far has been the forward problem: assume the thickness $d$ and index $N_1$ are known and compute the reflectance. Real metrology needs the reverse, finding the $d$ and $N_1$ that best match a measured reflectance (or $\Psi, \Delta$) spectrum. Fitting the nonlinear model $r_{total}(d, N_1, \theta, \lambda)$ to observations is exactly the nonlinear least-squares problem of the [Levenberg-Marquardt post](/en/posts/optimization-levenberg-marquardt/). The thickness fitted there is an input parameter of the $r_{total}$ model derived here. Updating $d$ and $N_1$ to reduce the residual between the TMM or SMM reflectance and the measured values is the core loop of real thin-film metrology software. When the observations cover a whole (angle, wavelength) grid, as in section 10, the residual vector grows accordingly and the Jacobian carries more information, which helps break the correlation between parameters.

## Summary

This third post treated the multiple-reflection interference that arises once two or more interfaces are stacked. Three points stand out.

First, even endlessly repeated reflections fold, as a geometric series, into a closed expression in a few interface coefficients and one phase. The same expression shows that at a single wavelength the reflection coefficient is periodic in thickness (283.2 nm for SiO$_2$ at 633 nm and 70°). For very thin films, $\Delta$ is far more sensitive than reflectance.

Second, extending this folding to a multilayer yields three strategies: Rouard, TMM and SMM. None is simply the best. Rouard for the spectrum of an ideal sample, TMM for transmission or coating design, SMM when non-ideal conditions must be added locally.

Third, for the single-thickness problem of this post, when the wavelength axis carries limited information because the band is narrow, scanning the angle as well separated a pair that wavelength alone could not. The effect with many unknown layers was not calculated. With a wide band, the wavelength axis alone already separated the pair. Either way, realizing that benefit means treating the reflection calculation as a tensor over the angle and wavelength axes.

This concludes the "Ellipsometry Foundations" series. Post 1 started from the fact that intensity measurement loses phase, and set up the Fresnel reflection coefficients as complex numbers. Post 2 compared two ways of writing a polarization state: the Jones formalism is usually enough, and the Mueller formalism takes over when there is depolarization. This post connected the $r_p$ and $r_s$ of a film to $M_{sample}$. Two paths remain. The numerical one, fitting this reflection model to measured data, is already covered by the [Levenberg-Marquardt post](/en/posts/optimization-levenberg-marquardt/) of the optimization series. The measurement one, putting this $M_{sample}$ into a real instrument and reading $\Psi$ and $\Delta$ from intensity signals, continues in [Ellipsometer Instrumentation 1](/en/posts/ellipsometer-oblique-rotating-element/). Its section 1 (Why Oblique Incidence) places the $M_{sample}$ of section 7 here between the polarization state generator and analyzer. Its section 2 (The Detector Gives Only One Intensity) starts from the $(S_1, S_2, S_3) = (-\cos2\Psi, \sin2\Psi\cos\Delta, -\sin2\Psi\sin\Delta)$ obtained for +45° linear input. The back-focal-plane configuration of section 9 is treated as an instrument in part 6 of that series (the back focal plane).

## References

- Youngjoon Kim, "라인 스캔 분광기와 후초점면 분광 간섭을 이용한 스냅샷 각도 분해 엘립소메트리 개발" [Development of snapshot angle-resolved ellipsometry using a line-scan spectrometer and back-focal-plane spectral interference], Ph.D. dissertation, Seoul National University, 2025 (in Korean), sections 2.4–2.5.
- O.S. Heavens, *Optical Properties of Thin Solid Films*, Dover, 1991 (the standard derivation of the characteristic matrix method).
- M. Born and E. Wolf, *Principles of Optics*, 7th ed., Cambridge University Press, 1999, ch. 1 (matrix formulation of multilayer reflection and transmission). It uses the time factor $e^{-i\omega t}$, so the sign of the imaginary unit is opposite to this series ($e^{j\omega t}$).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, section 4.6 (the scattering-matrix formalism with interface and layer matrices, eqs. 4.146–4.170).
- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, ch. 2, 5 (thin-film modelling and fitting of ellipsometric data).

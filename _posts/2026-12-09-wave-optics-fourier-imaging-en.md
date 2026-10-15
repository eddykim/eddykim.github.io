---
layout: post
date: 2026-12-09 20:00:00 +0900
title: "Wave Optics Foundations 4 — Why Can't a Microscope See Below Half a Wavelength?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/wave-optics-fourier-imaging/
page_id: wave-optics-fourier-imaging
categories: [Optics, Wave Optics]
tags: [wave-optics, fourier-optics, optical-transfer-function, mtf, phase-contrast]
description: "Treating the lens as a Fourier transformer, this post follows coherent and incoherent imaging, the MTF, the condenser diaphragm and the phase-contrast microscope through a single filter picture."
math: true
---

In front of a university building in Jena, Germany, stands a sphere carved from stone. The pedestal bears the name of Ernst Abbe (1840–1905), and the band around the sphere carries a single formula between diagonal-stripe decorations (Figure 1).

<img src="/assets/img/posts/wave-optics-fourier-imaging/ext-abbe-memorial.jpg" alt="The Ernst Abbe memorial in Jena and the resolution formula carved on its sphere" width="760">
_Figure 1. The Abbe memorial (1977) in front of the building at Fürstengraben 23, Jena. On the right, the part of the band bearing the formula $d = \lambda/(2n\sin\alpha)$ is cropped and enlarged. (Source: KaurJmeb at the English Wikipedia project, [Ernst Abbe memorial](https://commons.wikimedia.org/wiki/File:Ernst_Abbe_memorial.JPG), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/), right side cropped and enlarged)_

$d = \lambda/(2n\sin\alpha)$ is the smallest spacing a microscope can resolve. Here $n\sin\alpha$ is the numerical aperture (NA) of the objective. With 500 nm visible light and an NA 1.4 oil-immersion objective it gives about 0.18 µm. In air NA is below 1, so this value never drops below $\lambda/2$. The 0.18 µm is smaller than $\lambda/2$ = 0.25 µm because the wavelength inside the immersion oil is shorter. However good the lens, two lines closer than that look like one. Abbe derived this limit in 1873 while designing microscopes in Carl Zeiss's workshop [1, 2]. His idea fits in one sentence: the microscope image is the result of the diffraction pattern made by the specimen interfering again [3].

[Part 3](/en/posts/wave-optics-diffraction/) ended by saying that the pattern of light that has travelled far is the Fourier transform of the aperture, and that a convex lens pulls that pattern into its focal plane. This post starts there and answers the following questions.

- Why can't a microscope (in air) see things smaller than half a wavelength? Why is there a 2 in the formula above?
- Why does an image taken with coherent light, like a laser, differ from one taken with lamp light? Which is sharper?
- What does the modulation transfer function (MTF) on a lens data sheet measure? Why does the contrast of stripes reverse when the focus is slightly off?
- How far should a microscope's condenser diaphragm be opened?
- How do we make nearly transparent cells visible?

All the answers come from the same picture. Forming an image with a lens is a filtering operation: the object is split into spatial frequencies, the lens aperture filters out some of them, and the rest are gathered again. This view is called Fourier optics [4].

## 1. How does a lens split an object into spatial frequencies?

In Part 3, seeing the far-field pattern of a 2.5 cm aperture required going much further than about 1 km (at λ = 0.6 µm). How does a lens shorten that distance? One point of the far-field pattern is the light that left the aperture in a particular direction. A convex lens gathers parallel light arriving from one direction to one point of its focal plane. So the back focal plane of a lens holds the light sorted by direction, which is the far-field pattern.

A thin lens can be treated as a plate that multiplies the light by a quadratic phase $\exp[-ik(x^2+y^2)/2f]$. Place the object a distance $d$ in front of the lens and compute the field in the back focal plane, and the quadratic phase of the Fresnel integral cancels against that of the lens. In particular, when $d = f$ no quadratic phase remains at all, and the result is [4]

$$
U_f(u, v) = \frac{1}{i\lambda f}\iint U(\xi, \eta)\,\exp\!\left[-i\frac{2\pi}{\lambda f}(u\xi + v\eta)\right] d\xi\, d\eta
$$

Here $U$ is the field of the object placed in the front focal plane and $(u, v)$ are coordinates in the back focal plane. In other words, the field in the back focal plane is the exact Fourier transform of the object, and position $u$ corresponds to spatial frequency $\nu = u/\lambda f$. A pattern of period $\Lambda$ gathers its light at a point $\lambda f/\Lambda$ from the centre of the focal plane. With the object placed against the lens the intensity pattern is the same, but a quadratic phase remains, so it is not an exact transform once phase is included [4]. As a check in this post's calculation, propagating an object in the front focal plane by the angular spectrum method gives a back-focal-plane field that matches the directly computed Fourier transform to within 0.01 % in intensity and 0.3 % in phase.

The back focal plane is not an imaginary plane. The back focal plane of a microscope objective physically exists inside the tube, and it can be seen by removing the eyepiece or inserting a Bertrand lens [3]. [Ellipsometer part 6](/en/posts/ellipsometer-back-focal-plane-micro/) read reflections at all angles at once in this plane.

Applying the Fourier transform once more returns the original function (inverted). So if two lenses are placed $2f$ apart, with the object a distance $f$ in front of the first and the image a distance $f$ behind the second, the first lens transforms, the second reverses it, and an image of the object appears. The distance from object to image is $4f$, so this is called a 4f system. A mask placed in the Fourier plane in the middle is itself a filter on spatial frequencies [4, 5].

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig2-lens-fourier-4f.png" alt="Structure of the 4f system, the object, the back-focal-plane intensity, and the image passed through the lens aperture" width="780">
_Figure 2. (a) The 4f system. Light leaving the object plane at different heights in the same direction (red, green) gathers at one point of the Fourier plane. (b) An object made of vertical stripes (period 4 µm), horizontal stripes (period 0.8 µm) and a round hole. (c) Intensity in the Fourier plane (log scale). The vertical stripes make points on the horizontal axis and the horizontal stripes make points on the vertical axis. The cyan dashed circle is the range of spatial frequencies accepted by an NA 0.5 lens (wavelength 0.5 µm). (d) The coherent image made only from the light inside the dashed circle. The points of the 0.8 µm horizontal stripes lie outside the circle, so those stripes vanish from the image._

Let us read Figure 2. In (c) the vertical stripes repeat in the horizontal direction, so they make points on the horizontal axis; with a 4 µm period the points stand every 0.25 cycles/µm. The horizontal stripes have a 0.8 µm period and make points at 1.25 cycles/µm, outside the dashed circle. The lens does not accept that light, so in the image of (d) the horizontal stripes turn into a uniform grey. The lens aperture is a threshold in spatial frequency. The next section says the same thing in Abbe's language.

## 2. What must a lens accept for stripes to be visible?

What must a lens accept in order to see fine stripes in a specimen? Illuminate a grating of period $\Lambda$ head-on with a plane wave, and, by the grating equation of Part 3, Section 6, the 0th order goes straight, the ±1st orders go at $\sin\theta = \pm\lambda/\Lambda$, and the ±2nd orders at twice that angle. The objective accepts only the orders within its half-angle $\alpha$.

- If only the 0th order enters, the image has uniform brightness. All the information about the stripes is in the ±1st orders and higher.
- If the 0th and ±1st orders enter, the three interfere to make a cosine pattern of period $\Lambda$. The shape is blunt, but the period is right.
- The more high orders enter, the closer the shape comes to the original grating.

In other words, an image forms only if at least two orders enter the lens [2, 3]. With head-on illumination, the ±1st orders enter when $\lambda/\Lambda \le n\sin\alpha$, so the shortest resolvable period is

$$
\Lambda_{\min} = \frac{\lambda}{\mathrm{NA}}\qquad (\text{head-on illumination})
$$

Tilt the illumination and the 0th order can be sent to one edge of the lens aperture while the +1st order is caught at the opposite edge. The 0th order enters at one edge of the aperture and the 1st at the opposite edge, so the two orders can be up to $2\,\mathrm{NA}$ apart in angle, and the resolvable period halves:

$$
\Lambda_{\min} = \frac{\lambda}{2\,\mathrm{NA}}\qquad (\text{oblique illumination})
$$

This is the formula carved on the memorial in Figure 1, and the reason it contains a 2. In other words, the 2 appears because the 0th and 1st orders each use one edge of the lens aperture. Numerically, at NA 0.5 head-on illumination shows a grating of period 1.05 λ/NA and makes a 0.95 λ/NA grating perfectly uniform. Illumination tilted to the pupil edge shows gratings down to 0.52 λ/NA and becomes uniform at 0.48 λ/NA. In numbers, with an NA 1.4 oil-immersion objective and 500 nm light, the limit is 0.36 µm for head-on illumination and 0.18 µm for oblique illumination.

Abbe's theory was not readily accepted at first. In 1906 Alexander B. Porter showed how the image changes when masks of various shapes are placed in the back focal plane [2, 4]. Figure 3 repeats that experiment by calculation.

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig3-abbe-porter.png" alt="Images of a mesh object with six different masks in the Fourier plane" width="780">
_Figure 3. Calculation of the Abbe–Porter experiment. The object is a wire mesh (period 4 µm) with a few specks of dust on it, illuminated with coherent light. Top row: Fourier-plane intensity and the mask (blue is blocked); bottom row: image intensity. (a) All pass. (b) Only the central 0th order and its surroundings pass: the mesh vanishes and only blurry dust remains. (c) Horizontal slit: only the orders on the horizontal axis remain, so only the vertical wires are visible. (d) Only the 0th and ±1st orders pass: the mesh becomes a smooth dot pattern. (e) Only the 0th order blocked: the contrast reverses and the wires become bright. (f) Only the mesh orders blocked: the mesh vanishes and only the dust remains, bright._

A few things stand out in Figure 3. The horizontal slit of (c) passes only horizontal spatial frequencies, so it keeps only the vertical wires, which repeat in the horizontal direction. As in (e), blocking just the 0th order reverses the contrast. In a mesh with open fraction $m$, removing the 0th order makes the amplitude $1-m$ at the holes and $-m$ at the wires, so for $m > 1/2$ the wires are brighter [1]. (f) is the prototype of a technique used in practice. The orders of a periodic pattern gather at a few points in the Fourier plane, so blocking just those points erases the periodic pattern and leaves only the dust or defects that do not follow the period [1].

## 3. How does the image change under coherent light?

Can all of image formation be written with the single circle of Figure 2(d)? When the object is illuminated with coherent light, yes. The lens passes only those spatial-frequency components of the object field that lie inside the pupil and discards the rest. Without aberrations, the components that pass keep their magnitude and phase [4]. The amplitude spectrum of the image is the object's amplitude spectrum multiplied by

$$
H_c(\nu) = \begin{cases} 1, & \lvert\nu\rvert \le \mathrm{NA}/\lambda \\ 0, & \lvert\nu\rvert > \mathrm{NA}/\lambda \end{cases}
$$

Here $\nu$ is the spatial frequency in the object plane, and $H_c$ is called the amplitude transfer function. In other words, the transfer function of coherent imaging is nothing other than the shape of the pupil, with cutoff frequency $\mathrm{NA}/\lambda$ [2, 4]. In the spatial domain, the image amplitude is the object amplitude convolved with the amplitude point-spread function (PSF). This PSF is the Fourier transform of the pupil, the Airy amplitude $2J_1(v)/v$. Goodman's book writes the cutoff as $w/\lambda z_i$ in terms of the lens radius $w$ and image distance $z_i$; rewritten in terms of the object-side NA, it becomes the formula above.

This filter has an abrupt threshold. An abruptly cut filter leaves ripples in the image. The statement in [Metrology Image Processing 3](/en/posts/imgproc-frequency-domain-filtering/) that cutting frequencies abruptly leaves rings in space happens here in a lens. Figure 4(a) is the result.

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig4-coherent-incoherent-edge-points.png" alt="Images of an edge, two points and three bars under coherent and incoherent illumination" width="780">
_Figure 4. An imaging system with NA 0.5 at 0.5 µm. (a) Image of an edge. Coherent illumination (red) passes 1/4 intensity at the boundary, rises to 1.19 and ripples. Incoherent illumination (blue) passes 1/2 at the boundary and rises smoothly. (b) Two points separated by the Rayleigh spacing (0.61 λ/NA). Two incoherent points have a dip of 0.735. For two coherent points, depending on the phase difference, the dip disappears (0), stays the same (π/2) or drops all the way to 0 (π). (c), (d) Three bars spaced 0.7 λ/NA. Illuminated coherently in phase they nearly merge into one lump; illuminated incoherently they separate into three._

The image of the edge resembles the straight edge of Part 3, Section 3. The coherent image has exactly 1/4 intensity at the boundary [4]. On the bright side it overshoots to a maximum of 1.19 numerically. The incoherent image is at 1/2 at the boundary. So in an image taken with coherent illumination, taking the 50 % point as the boundary misreads the linewidth. Goodman points out that this difference matters when measuring linewidths on integrated-circuit masks [4].

The result for two points depends on the coherence of the illumination. Two incoherent points add in intensity, so the dip (the minimum intensity between the peaks divided by the maximum intensity of the peaks) is 0.735 as in Part 3's Rayleigh criterion. Two coherent points add in amplitude, so the phase difference between them changes the result (Figure 4(b)) [4]. In phase, the amplitudes reinforce in the middle and no dip forms. At a phase difference of π/2 the cross term vanishes and the result is exactly the incoherent one. At a phase difference of π the amplitudes cancel in the middle and the dip falls to 0. For two in-phase points to show a 0.735 dip, the spacing must widen to $0.82\lambda/\mathrm{NA}$ [2]. The calculation also gives 0.82. This is where the textbook values come from: microscopy textbooks quote a resolution of $0.82\lambda/\mathrm{NA}$ for coherent illumination and $0.61\lambda/\mathrm{NA}$ for incoherent illumination.

## 4. How much stripe contrast survives under incoherent light?

Now consider the case where each point of the object emits light incoherently with the others, as with lamps or fluorescence. The image of each point is the intensity of the Airy pattern $\lvert h\rvert^2$, and since they do not interfere, the intensities add. In other words, the image intensity is the object intensity convolved with $\lvert h\rvert^2$ [1, 4]. In the frequency domain, the spectrum of the object intensity is multiplied by the Fourier transform of $\lvert h\rvert^2$. That function, normalized to 1 at zero frequency, is the optical transfer function (OTF).

The Fourier transform of $\lvert h\rvert^2$ is the autocorrelation of the Fourier transform of $h$, that is, of the pupil. So the OTF can be written as [4]

$$
\mathrm{OTF}(\nu) = \frac{\iint P\!\left(f + \tfrac{\nu}{2}\right) P^*\!\left(f - \tfrac{\nu}{2}\right) df}{\iint \lvert P(f)\rvert^2\, df}
$$

Here $P$ is the pupil function (a circle of radius $\mathrm{NA}/\lambda$ in the spatial-frequency plane). Without aberrations, the numerator is the overlap area of two circles offset by $\nu$. In other words, the OTF is "overlap area / area of one circle" (Figure 5(b)). Why the overlap area? For a stripe of frequency $\nu$ to appear in the image intensity, light that passed through two places in the pupil separated by $\nu$ must interfere. The number of ways to choose such pairs inside the pupil is the overlap area [4]. For a circular pupil there is a closed form [4]:

$$
\mathrm{OTF}(\nu) = \frac{2}{\pi}\left[\arccos s - s\sqrt{1-s^2}\right],\qquad s = \frac{\nu}{2\,\mathrm{NA}/\lambda}
$$

It reaches 0 at $\nu = 2\mathrm{NA}/\lambda$, where the two circles separate completely. The cutoff frequency is twice the coherent one (Figure 5(a)) [2]. Numerically, the FFT of $\lvert h\rvert^2$, a direct calculation of the pupil overlap and the closed form above agree to within 0.1 %.

The OTF is in general complex. Its magnitude is called the MTF (modulation transfer function) and its phase the PTF (phase transfer function) [1]. The meaning of the MTF is simple. Measure the contrast of a pattern whose brightness varies sinusoidally by its modulation

$$
M = \frac{I_{\max} - I_{\min}}{I_{\max} + I_{\min}}
$$

and the modulation of the image is the modulation of the object multiplied by the MTF at that frequency [1]. In other words, the MTF is a table, frequency by frequency, of "how much of the contrast of stripes this fine survives". The MTF curve on a lens data sheet is exactly this. For a point on the axis of an axially symmetric system the PTF is 0; with asymmetric aberrations such as coma, the PTF shifts the stripes sideways [1].

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig5-transfer-functions.png" alt="The coherent amplitude transfer function and the OTF, the pupil-overlap picture, and images of an amplitude-cos object and an |cos| object" width="780">
_Figure 5. (a) For an aberration-free circular pupil, the amplitude transfer function of coherent illumination (red) and the OTF of incoherent illumination (blue). (b) The OTF is the overlap area of two pupils offset by spatial frequency ν, divided by the pupil area. (c), (d) Images of two objects whose intensity is the same cos². Object A has amplitude cos and object B has amplitude \|cos\|, with an amplitude frequency of 0.75 NA/λ. Under coherent illumination A is imaged perfectly and B as uniform grey. Under incoherent illumination both become the same image with contrast reduced by a factor of 0.14._

Is incoherent imaging, with twice the cutoff frequency, then twice as good? No [4, 5]. The coherent cutoff applies to the amplitude spectrum and the incoherent cutoff to the intensity spectrum. They are filters acting on different quantities, so comparing their cutoff frequencies is meaningless. Goodman shows this with one example (Figure 5(c), (d)) [4]. The intensity of both objects is $\cos^2$; only the sign of the amplitude differs. The amplitude spectrum of object A, whose amplitude is cos, lies at ±0.75 NA/λ and passes the pupil, so the coherent image is perfect. Object B, whose amplitude is \|cos\|, has a fundamental amplitude frequency of 1.5 NA/λ, outside the pupil, so the coherent image is uniform grey (intensity $4/\pi^2 \approx 0.41$). Incoherent illumination sees only intensity, so it cannot tell the two objects apart and makes the same image with contrast reduced by the OTF value of 0.14 at the intensity frequency 1.5 NA/λ. In other words, which illumination is better depends on the phase of the object.

The three bars of Figure 4(c), (d) tell the same story. At a spacing of 0.7 λ/NA the coherent image is nearly one lump (dip between bars 0.83), while the incoherent image separates them (dip 0.63). But widen the spacing to 0.9 λ/NA and the order reverses: the dip of the coherent image (0.14) is deeper than that of the incoherent image (0.32), because the ripples of the coherent image carve the space between the bars darker. Goodman's conclusion is practical. If you can, incoherent illumination is the better choice, because it has none of the side effects of coherence such as ripples and speckle. In holography and optical information processing, however, coherence cannot be avoided [4].

## 5. How far should the condenser diaphragm be opened?

Which of the two is a real microscope? In a transmission microscope the specimen does not emit light itself; it is illuminated by lamp light gathered through a condenser. Close the condenser diaphragm down tightly and the specimen receives nearly a single plane wave from one direction, which is close to coherent illumination. Open it fully and light from many directions mixes incoherently, which is close to incoherent illumination. In between is partial coherence [2]. In terms of the van Cittert–Zernike theorem of [Part 2](/en/posts/wave-optics-coherence/), the size of the condenser diaphragm sets the spatial coherence distance on the specimen.

The calculation can use Abbe's picture as it is. Each point of the condenser diaphragm illuminates the specimen with one plane wave in a particular direction. The image for one plane wave is the coherent image of Section 3. Light from different points of the diaphragm is mutually incoherent, so the intensity of the coherent image is computed for each point and the results are added. Harold H. Hopkins's theory of partially coherent imaging, formalized in 1951 and 1953, organizes this calculation into a function called the transmission cross-coefficient [2]. Let the NA ratio of condenser to objective be $$\sigma = \mathrm{NA}_{\text{condenser}}/\mathrm{NA}_{\text{objective}}$$. At the most tilted plane wave, the 0th order is pushed sideways within the pupil as in Section 2, so the highest frequency that can pass is

$$
\nu_{\max} = (1+\sigma)\,\frac{\mathrm{NA}}{\lambda}
$$

In other words, $\sigma = 0$ gives the coherent $\mathrm{NA}/\lambda$ and $\sigma = 1$ gives the incoherent $2\mathrm{NA}/\lambda$.

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig6-partial-coherence.png" alt="Images of a bar grating and an edge, and contrast transfer of a weak sinusoid, for condenser NA ratios σ" width="780">
_Figure 6. Partially coherent imaging with the condenser-diaphragm NA ratio σ set to 0, 0.3, 0.6 and 1 (objective NA 0.5, wavelength 0.5 µm). (a) Image of a bar grating with period 0.8 λ/NA. At σ = 0 the fundamental frequency 1.25 NA/λ lies outside the pupil, so the image is uniform; as σ grows the pattern comes back. (b) Image of an edge. As σ grows the ripples subside. (c) Contrast transfer of a weak sinusoid with 10 % contrast. σ = 0 cuts off abruptly at 1; as σ grows, the cutoff widens to (1 + σ)NA/λ, but contrast is shaved starting from the lower frequencies. At σ = 1 it coincides with the OTF of incoherent imaging (dashes)._

Figure 6 shows the price of closing and opening the condenser diaphragm. Opening it (larger σ) reveals finer patterns and calms the ripples at edges (a, b). In exchange, the contrast of medium-fine patterns drops (c). The σ = 0.6 curve passes contrast unchanged up to 0.4 NA/λ and then falls, while the σ = 0 curve keeps contrast 1 right up to the cutoff. This is the microscope user's experience that closing the diaphragm raises contrast and lowers resolution [3].

So how far should it be opened? The conclusion first: there is no single right answer. But there is a starting point. Murphy and Davidson advise setting the condenser diaphragm so that about 3/4 (70–75 %) of the diameter of the objective's back aperture appears lit, then adjusting for the specimen [3]. They also warn against closing the diaphragm to reduce brightness when the image is too bright. Brightness is reduced with the lamp voltage or ND filters; the diaphragm is used only to balance resolution and contrast [3]. Considering only two-point resolution, in the calculation of Hopkins and Barham (1950) the minimum resolvable distance falls from 0.82 λ/NA at σ = 0 to slightly below 0.61 λ/NA around σ ≈ 1.5 [2]. That result is for two pinholes of equal brightness, not a prescription for every specimen. Köhler illumination, introduced by August Köhler in 1893, images the lamp filament onto the plane of the condenser diaphragm. The specimen is therefore illuminated evenly, while the angular range of the illumination (σ) can be adjusted separately with the condenser diaphragm alone [2, 3]. It is the standard arrangement of transmission microscopes today.

## 6. How do aberrations erode the transfer function?

[Geometrical Optics Foundations 4](/en/posts/geometric-optics-aberrations/) wrote aberrations as the wavefront aberration $W$ over the pupil. In Fourier optics an aberration enters as a generalized pupil function, the pupil function multiplied by a phase [4]. Geometrical Optics Foundations 4 wrote that phase as $-kW$ for $W>0$, whereas Goodman writes $+kW$. This post measures $W$ in wavelengths and follows Goodman, writing the phase as $e^{i2\pi W}$. Only the sign convention differs; the results below are the same. In coherent imaging the cutoff frequency stays the same and only the phase of the passed components changes. In incoherent imaging the OTF is an autocorrelation, so when the phases of the two overlapping pupils fall out of step they partly cancel each other. That is why aberrations can never raise the MTF; they only lower it [4].

Take the most common aberration, defocus. Defocus is the pupil phase $W_{020}\rho^2$, and in the paraxial case it is related to the distance $z$ from the focal plane by $W_{020} = \tfrac{1}{2}\mathrm{NA}^2 z$ [2].

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig7-defocus-mtf-spoke.png" alt="OTF curves with defocus and images of a spoke target" width="780">
_Figure 7. (a) OTF of a circular pupil with defocus $$W_{020}$$ set to 0, 0.25, 0.5, 0.75 and 1 λ (incoherent imaging). At 0.75 λ and 1 λ there are ranges where the OTF becomes negative. (b) Image of a target with 36 spokes (in focus). Inside the cyan dashed circle the spokes are spaced more finely than the cutoff frequency and smear into grey. (c) The same target with $$W_{020}$$ = 1 λ. The spokes blur. The region between the two orange dotted rings (radius 3.6–14.7 µm) corresponds to the frequencies where the OTF is negative. The negative values are small, so the contrast reversal appears only weakly on the blurred spokes._

In Figure 7(a), the more the focus is off, the more the OTF collapses, starting from high frequencies. At 0.25 λ (Rayleigh's quarter-wave criterion) the curve is still gentle, but at 1 λ it already crosses 0 near 1/5 of the cutoff and goes negative. A negative OTF means reversed contrast: places that should be bright are dark and places that should be dark are bright. This can be followed in the spoke target (c). The spokes are thick on the outside and thin on the inside, so the spatial frequency changes along the radius, and in the rings corresponding to frequencies where the OTF is negative the contrast of the spokes flips [4]. For a circular pupil, this calculation shows the OTF first going negative when $W_{020}$ exceeds about 0.64 λ. For a square pupil the threshold is λ/2 [4]. Converting 0.64 λ with the relation above gives $z \approx 1.3\lambda/\mathrm{NA}^2$. In Part 3, Section 5, the distance at which the on-axis intensity first reaches 0 was $2\lambda/\mathrm{NA}^2$, so stripe contrast already begins to reverse somewhat closer than that. At NA 0.5 and 500 nm it is about 2.6 µm.

One point is worth noting when handling MTFs in practice. The MTFs of mutually independent elements, such as a lens, a sensor and vibration, multiply to give the overall MTF. But when the aberrations of one lens can be cancelled by another, as when two lenses are joined in tandem, they must not be multiplied [1]. Aberrations add at the amplitude (pupil) stage, and the MTF comes from the autocorrelation of the result.

## 7. How do we make nearly transparent cells visible?

Living cells are nearly transparent. They do not absorb light; their refractive index differs slightly from the surroundings and only changes the phase of the light passing through. Under a brightfield microscope they are almost invisible. Why? Expanding a transmittance with a slightly altered phase gives

$$
e^{i\varphi(x)} \approx 1 + i\varphi(x)
$$

Here $\varphi$ is the phase added by the specimen, assumed to be much smaller than 1. The 1 is the 0th-order light that went straight through the specimen, and $i\varphi$ is the diffracted light scattered by it. The $i$ attached to the diffracted light means it is 90° out of phase with the 0th order. Computing the intensity gives $\lvert 1 + i\varphi\rvert^2 = 1 + \varphi^2$: the cross term vanishes and only $\varphi^2$ remains [2, 4]. In other words, because the diffracted light is at right angles to the 0th order, it cannot change the brightness by interference. A phase of 0.1 rad changes the brightness by only 1 %.

Frits Zernike's idea was to change that 90°. In the back focal plane the 0th order gathers at one point (or one ring), while the diffracted light spreads widely. So place a plate a quarter-wave thick only where the 0th order gathers, shifting its phase by 90°. Advancing the 0th order by a quarter wave multiplies it by $-i$, giving

$$
I \approx \lvert -i + i\varphi\rvert^2 = (1-\varphi)^2 \approx 1 - 2\varphi
$$

In other words, the intensity now changes linearly with the phase, and 0.1 rad changes the brightness by 20 %. Advancing the 0th order makes regions of large phase (thick regions) look dark (positive phase contrast, following Hecht's convention; Goodman's is the opposite); retarding it gives $1 + 2\varphi$ and makes them look bright (negative phase contrast) [2, 3]. Adding a film that partly absorbs the 0th order raises the contrast further. Reducing the 0th-order amplitude by a factor $a$ makes the intensity $a^2 \mp 2a\varphi$, so the relative change in brightness grows to $2\varphi/a$. Cutting the 0th-order intensity to 1/9 ($a = 1/3$) triples the sensitivity [2]. The phase ring of a real phase-contrast objective cuts the 0th-order amplitude by 70–75 % [3]. Dark field, which blocks the 0th order completely, gives intensity of order $\varphi^2$ and is insensitive to weak phases [2].

Zernike published the method in 1934–35 and received the Nobel Prize in Physics in 1953 [2]. According to his Nobel lecture as related by Murphy and Davidson, the idea came not from microscopes but from diffraction gratings. To make visible the false spectral lines caused by periodic errors in ruling the grooves, he placed a phase plate in the grating spectrum, and through Abbe's grating analogy he carried this over to the microscope [3]. A real microscope puts an annular diaphragm in the condenser and a phase ring of the same shape where its image falls in the back focal plane. Zeiss introduced phase-contrast objectives from the early 1940s [3].

<img src="/assets/img/posts/wave-optics-fourier-imaging/ext-brightfield-phase-contrast.jpg" alt="Photographs of the same cheek epithelial cells in brightfield and phase contrast" width="600">
_Figure 8. The same cells seen in brightfield (left) and phase contrast (right). Cells whose outlines are barely visible in brightfield show their nuclei and internal structure in phase contrast. (Source: GcG, [Brightfield phase contrast cell image](https://commons.wikimedia.org/wiki/File:Brightfield_phase_contrast_cell_image.jpg), public domain)_

<img src="/assets/img/posts/wave-optics-fourier-imaging/en/fig9-phase-contrast.png" alt="Calculation of a weak phase object seen in brightfield, defocus, positive and negative phase contrast, and dark field" width="780">
_Figure 9. Calculation of a phase object imitating a cell (maximum 0.4 rad) seen with NA 0.5 at 0.5 µm. (a) The phase added by the specimen. (b) Brightfield, in focus. The modulation is 0.003, almost invisible. (c) Brightfield, defocused by 1.5 λ. Outlines and granules appear along with ripples. (d) Positive phase contrast (0th order advanced by λ/4 and its intensity reduced to 1/4). The modulation rises to 0.51; the thick nucleus and granules are dark, and a bright halo forms around the cell. (e) Negative phase contrast (0th order retarded by λ/4). The contrast is reversed. (f) Dark field. The edges glow brightly. (b) and (c) share the same brightness scale; the others each use their own._

Figure 8 shows this effect in a real photograph: cells whose outlines are barely visible in brightfield show their nuclei and internal structure in phase contrast. The calculation of Figure 9 shows two habits of phase-contrast images. First, a bright band appears around the cell: the halo. The phase plate covers not only the 0th order but also the very low frequencies around it. The slow phase variation across the whole cell is treated like the 0th order and has its phase shifted too, so a band of opposite contrast forms at the edges of large objects [3]. Second, shade-off, in which the middle of a large object drifts back toward the background brightness (the middle of the cell body in (d)). It arises for the same reason. So the brightness of a phase-contrast image is not a map of phase. It should be read as the output of a filter that emphasizes small structures.

Figure 9(c) is also worth noting. Deliberately defocusing makes a phase object visible. Defocus multiplies the pupil by a phase that differs from frequency to frequency, and at frequencies where that phase is close to 90° it does the same job as the phase plate. This is exactly how an electron microscope sees the phase of atoms. The contrast transfer function (CTF) and Scherzer defocus of [Electron Microscopy part 6](/en/posts/electron-microscopy-tem-phase-contrast-ctf/) are techniques for producing this 90° over a wide range of frequencies using aberrations and defocus. It is fair to regard them as Zernike phase contrast that uses the aberrations of an electron lens instead of an optical phase plate.

## Summary

- The exact Fourier transform of an object placed in the front focal plane appears in the back focal plane of the lens. Position $u$ corresponds to spatial frequency $u/\lambda f$. A mask in the middle plane of a 4f system, made of two lenses, becomes a spatial-frequency filter.
- Abbe's theory: the image is the interference of diffraction orders, and at least two orders must enter the lens. The resolvable period is $\lambda/\mathrm{NA}$ with head-on illumination and $\lambda/2\mathrm{NA}$ with oblique illumination. That is the 2 in the formula on the memorial.
- Coherent imaging is a filter on amplitude, and its transfer function is the pupil itself (cutoff $\mathrm{NA}/\lambda$). Edges ripple and the boundary sits at 1/4 height; the resolution of two points depends on their phase difference (0.82 λ/NA in phase).
- Incoherent imaging is a filter on intensity, and the OTF is the autocorrelation of the pupil (cutoff $2\mathrm{NA}/\lambda$). The MTF is the fraction of stripe contrast that survives. A doubled cutoff does not mean twice as good; which is better depends on the object.
- Under partial coherence the cutoff is $(1+\sigma)\mathrm{NA}/\lambda$. Closing the condenser diaphragm raises contrast and lowers resolution. The starting point is 70–75 % of the diameter of the objective's back aperture.
- Aberrations only lower the MTF. For a circular pupil, when defocus exceeds about 0.64 λ the OTF goes negative and the contrast of stripes reverses.
- In a weak phase object the diffracted light is 90° out of phase with the 0th order and is invisible in brightfield. Shifting the phase of the 0th order by λ/4 makes the intensity $1 \mp 2\varphi$, proportional to phase (Zernike phase contrast). Reducing the 0th order raises the sensitivity, with halos and shade-off as the price.

The light covered so far has been plane waves and spherical waves. Light from a laser is neither. It seems not to spread and then spreads after all, and its phase shifts by π through focus. Just as the size of the condenser diaphragm set the coherence in this post, a laser beam's own width sets its own spreading. The next part is the Gaussian beam.

## References

1. E. Hecht, *Optics*, 5th ed. (Pearson, 2017), Section 11.3 — linear systems and convolution, the point-spread function, transfer functions (MTF, PTF, OTF, modulation, when MTFs may be multiplied). Section 13.2 — the transform lens and Abbe's theory of imaging, spatial filtering, the Abbe–Porter experiments, phase contrast and dark field.
2. M. Born and E. Wolf, *Principles of Optics*, 7th ed. (Cambridge University Press, 1999), Section 8.6.3 — image formation in the microscope, Abbe's theory, resolution of two coherent points (0.82 λ/NA), Zernike phase contrast. Section 9.5 — imaging of extended objects and the frequency response, effect of defocus. Section 10.6 — partially coherent imaging, the condenser and resolution (Hopkins and Barham), Köhler illumination, the transmission cross-coefficient.
3. D. B. Murphy and M. W. Davidson, *Fundamentals of Light Microscopy and Electronic Imaging*, 2nd ed. (Wiley-Blackwell, 2013), Chapter 1 — Köhler illumination and setting the condenser diaphragm. Chapter 5 — diffraction and Abbe's theory, the diffraction pattern in the objective's back focal plane. Chapter 6 — resolution and contrast. Chapter 7 — phase-contrast microscopy and dark field, Zernike's idea.
4. J. W. Goodman, *Introduction to Fourier Optics*, 2nd ed. (McGraw-Hill, 1996), Chapter 5 — Fourier-transforming properties of lenses and image formation. Chapter 6 — transfer functions of coherent and incoherent imaging, the OTF of a circular pupil, aberrations and defocus, comparison of coherent and incoherent imaging. Section 8.1 — the Abbe–Porter experiments, Zernike phase contrast.
5. B. E. A. Saleh and M. C. Teich, *Fundamentals of Photonics*, 3rd ed. (Wiley, 2019), Section 4.4 — the 4f system and spatial filtering, the transfer function of single-lens imaging. Section 12.3 — comparison of the transfer functions of coherent and incoherent imaging.

---
title: "Ellipsometer Instrumentation 2 — Why Does an Ellipsometer Need Different Parts for Different Wavelengths?"
lang: en
lang-exclusive: ["en"]
permalink: /posts/ellipsometer-components/
page_id: ellipsometer-components
date: 2026-09-26 20:00:00 +0900
categories: [Optics, Instrumentation]
tags: [ellipsometry, instrumentation, polarizer, retarder, photoelastic-modulator, spectrometer]
description: Which wavelengths an ellipsometer can work at is decided by whether polarizing components exist that function there. This post takes the parts apart one by one.
math: true
---

[Post 1](/en/posts/ellipsometer-oblique-rotating-element/) sorted the methods by what they cannot measure. RAE (rotating-analyzer ellipsometry) misses $S_3$ and so loses the sign of $\Delta$; adding a compensator recovers it but costs a second measurement; rotating the compensator settles it in one. Near its end, that post claimed that the wavelength dependence of a single component decides the spectral capability of the whole instrument.

This post takes that claim from the component side. The idea is not new. Azzam and Bashara stated it as the first constraint back in 1977.

> The primary limiting factor in using ellipsometry in a given spectral interval is the availability of polarizing optical devices that function satisfactorily in that particular interval. Each spectral range will require its own polarizing optical elements, its own sources and its own detectors.

You do not pick an instrument and then choose a wavelength. You fix the wavelength, that narrows the components available, and those components decide the method. Where post 1 followed what each method cannot measure, this one follows where each component stops working.

## 1. Why Does the Same Wave Plate Become a Different Part at Every Wavelength?

Why should one crystal plate of fixed thickness behave differently from one wavelength to the next? Inside a birefringent crystal the speed of light depends on the polarization direction. While crossing the plate the fast component gradually gets ahead of the slow one, and a phase difference remains when the light leaves. A retarder is a part built on that phase difference, and a retarder made as a thin plate is a wave plate. The size of the phase difference, the retardance, is

$$\delta = \frac{2\pi d\,(n_e - n_o)}{\lambda}$$

where $d$ is the thickness and $n_e - n_o$ is the birefringence. The sign of that difference depends on the material, but only the magnitude of the retardance matters. In other words, the retardance is the path difference $d(n_e - n_o)$ counted in wavelengths. For the same path difference, a shorter wavelength fits more cycles into it, so the retardance is larger. Almost the whole of this post follows from this one line, because the retardance goes as $1/\lambda$.

For a single-wavelength instrument there is no difficulty. Pick a thickness that gives $\delta = 90°$ at that wavelength and you are done. Spectroscopic ellipsometry, however, has to cover a whole band at once. A part that is a quarter-wave plate at 550 nm becomes a half-wave plate when the wavelength halves, and an eighth-wave plate when it doubles. The same part becomes a different part at every wavelength.

<img src="/assets/img/posts/ellipsometer-components/en/fig1-wavelength-coverage.png" alt="Usable spectral range of each component" width="740">
_Figure 1. The spectral range each component can serve. The red band marks the 0.2–1 μm range a spectroscopic ellipsometer commonly covers, and an arrow means the source gives no upper limit. How far you push that band is what narrows the available materials._

Figure 1 is the summary of this post. Push toward the ultraviolet and calcite stops at 0.21 μm, so you switch to MgF$_2$. Push beyond the near-infrared edge (about 5–7 μm) and prism polarizers do not transmit, so you need a component built on an entirely different principle. The wire-grid bar at the bottom (2.5–100 μm) does not overlap the red band at all. Now take the parts one at a time.

## 2. How Far in Wavelength Can a Polarizer Go?

A polarizer extracts linear polarization from unpolarized light. How far it reaches depends on the transmission window of the crystal and on what sits between its two halves. Put it on the source side and it is a polarizer; put it on the detector side and it is an analyzer. It is the same object either way.

Performance is quoted as an extinction ratio.

$$\mathrm{ER} \equiv \frac{I_x}{I_y}$$

Here $I_x$ is the intensity along the transmission axis and $I_y$ the intensity perpendicular to it. An extinction ratio of $10^5$ means that one part in a hundred thousand of the unwanted component leaks through.

Most polarizers exploit the anisotropy of a birefringent crystal. Calcite at 2.1 eV has $n_o = 1.6584$ and $n_e = 1.4864$, a large split. Different refractive indices mean different critical angles, so choosing the incidence angle carefully removes one component by total internal reflection.

<img src="/assets/img/posts/ellipsometer-components/en/fig2-polarizer-types.png" alt="Structure of three polarizers" width="760">
_Figure 2. A Glan-type prism is a rectangular block cut along a diagonal with a thin layer between the halves. Whether that layer is air or glue decides whether the part survives in the ultraviolet. In the Rochon prism on the right the transmitted component is not deviated, so rotating the prism does not move the beam._

The Glan-Taylor prism is built that way. Separate two calcite prisms with an air gap and the condition for total internal reflection becomes $1/n_o < \sin\theta < 1/n_e$, which works out to $37.1° < \theta < 42.3°$; in that range the ordinary ray is removed by total internal reflection and only the extraordinary ray gets through. The extinction ratio reaches $10^5$ and the full calcite transmission window of 0.21–5 μm is available, which makes it the workhorse of spectroscopic ellipsometry.

The Glan-Thompson prism has almost the same structure but the halves are cemented. It accepts a much wider cone of incidence angles, a field of view of roughly 30° against about 10° for the air-spaced Glan types (Hecht), yet it is not used in spectroscopic ellipsometry, because the cement absorbs the ultraviolet. What cuts the band short is not the component but the glue holding it together.

The Rochon prism joins two prisms whose optic axes are perpendicular. The transmitted component sees the same refractive index in both halves, so the beam passes straight through without deviation. For an instrument that rotates a component, such as the $PSA_R$ of post 1, this matters a great deal, because turning the element does not shake the beam. Quartz, however, is optically active and needs a separate correction for it, while MgF$_2$ has no optical activity and transmits well enough in the ultraviolet to reach the deep UV around 6 eV.

Figure 3 shows a real Rochon prism splitting light. Rotate the prism and the deviated beam sweeps around the axis, while the straight beam stays put.

<img src="/assets/img/posts/ellipsometer-components/ext-rochon-prism-split.jpg" alt="A Rochon prism splitting a white laser beam into a straight beam and a deviated rainbow beam" width="680">
_Figure 3. A Rochon prism in a rotation mount with a white laser beam passing through it. One polarization goes straight through, while the orthogonal one leaves at an angle and is spread by wavelength into a rainbow. (Source: Peeter Piksarv, [Rochoni polarisaator](https://commons.wikimedia.org/wiki/File:Rochoni_polarisaator.jpg), [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0))_

Beyond about 5–7 μm, in the infrared, all of these prisms drop out. Over roughly 100–4000 cm$^{-1}$ a wire-grid polarizer is used instead: fine metal lines on a substrate. The field parallel to the wires drives electrons and is absorbed as heat or reflected, and only the perpendicular component gets through. This exploits anisotropy in absorption, which is dichroism rather than birefringence. The extinction ratio is about $10^3$, two orders below a Glan-Taylor, but in that band there is no alternative.

So the choice goes like this. A spectroscopic instrument that reaches the ultraviolet uses an air-spaced Glan-Taylor, with an extinction ratio of $10^5$ over 0.21–5 μm. An instrument that stays out of the ultraviolet and takes light in over a wide cone is better served by the wider field of a Glan-Thompson. A component that has to rotate calls for a Rochon, which does not move the beam, and beyond 2.5 μm in the infrared the wire grid is the answer, at an extinction ratio of $10^3$.

## 3. How Do You Build a Retarder for a Wide Band?

This is the heart of the post. If a polarizer produces linear polarization, a retarder turns it into circular or elliptical polarization. The compensator of post 1 is exactly this part.

### How Do Zero-Order and Multi-Order Plates Differ?

Write $\delta = 2\pi N$ in the equation of section 1 and $N$ counts the retardance in wavelengths. $N = 1/4$ gives a zero-order quarter-wave plate; $N = 1\frac{1}{4}$ or $2\frac{1}{4}$ gives a multi-order quarter-wave plate. At the design wavelength both deliver the same 90°.

A multi-order plate, though, amplifies both the wavelength dependence and the temperature dependence. It is several times thicker, so the same change in wavelength swings the retardance several times as far. Zero order is the better choice.

The difficulty is making one. Quartz has $n_e - n_o \approx 0.009$, and a zero-order quarter-wave plate at 633 nm works out to 17.6 μm thick. That nearly matches the roughly 18 μm quoted in the Handbook, and it is only a fraction of the thickness of a single sheet of paper. A plate that thin is hard to fabricate and hard to handle.

The practical answer is the pseudo zero-order plate. Bond two thick crystals with their fast axes perpendicular and their retardances subtract, leaving only the difference. Each piece is thick enough to handle, and the thicknesses are chosen so that what remains is a quarter wave.

### Why Does It Break Down Once the Band Gets Wide?

The most important fact in this post follows from this. The retardance does not have to be exactly 90°. It is enough that the analysis software knows the retardance at every wavelength.

There is an exception. If the retardance is zero or a multiple of 360°, no data can be obtained at all; if it is 180°, or 180° plus a multiple of 360°, $\Delta$ cannot be measured accurately. At 0° or 360° the part behaves like a sheet of plain glass, and at 180° it merely flips the polarization without producing any circular component.

Post 1 argued that RAE loses the sign of $\Delta$ because it cannot measure $S_3$. The same thing happens one level down, at the component: at particular retardances the part refuses to pass the information through.

<img src="/assets/img/posts/ellipsometer-components/en/fig3-retardance-vs-wavelength.png" alt="Wavelength dependence of retardance" width="760">
_Figure 4. Even with the thickness set for a quarter wave at 550 nm, the retardance crosses 180° toward the ultraviolet. Stacking two plates avoids that crossing, at the cost of an effective fast axis that swings with wavelength._

The blue curve in the upper panel of Figure 4 is that situation. A zero-order plate cut for 90° at 550 nm passes 180° at 275 nm. If the plan is visible light only, nothing is wrong; widen the band into the ultraviolet and the part becomes unusable.

This calculation holds the birefringence of quartz fixed at 0.009. The value is the same whether taken for red light (the Handbook) or for 400–700 nm (the author's dissertation), but in the ultraviolet the birefringence grows, reaching 0.013 at 200 nm (Ghosh 1999). With the real dispersion the 180° crossing moves up to about 306 nm. A real plate fails earlier than the approximation suggests.

### How Is a Wide-Band Retarder Designed?

So what is done over a wide band? There are three routes.

The first is to stack two plates with their axes offset. Overlapping plates designed for different wavelengths at an angle changes the shape of the combined retardance curve itself. Stacking quarter-wave plates cut for different wavelengths with their fast axes offset is an idea due to Johs and colleagues. The red curve in Figure 4 applies the same idea to a 300 nm plate and a 550 nm plate offset by 67.5°. These design values are not a reproduction of any published device; they were chosen here under the same approximation, and the curve was computed twice, once as a product of Jones matrices and once as a rotation on the Poincaré sphere, with the two agreeing. Across the whole 200–1000 nm range the retardance stays between 35.6° and 141.4°, reaching neither 0° nor 180°. With the real dispersion, though, the retardance near 200 nm climbs to 163–183° depending on how the thicknesses are chosen, so the margin all but vanishes at the ultraviolet end.

It is not free. The lower panel of Figure 4 shows the price. Stack two plates with offset axes and the effective fast axis of the combination is no longer fixed; it moves with wavelength, by 29° in this calculation. Both the retardance and the axis have to be known at every wavelength, so calibration becomes correspondingly harder. The idea of making up for imperfect components in calibration connects to Johs's regression calibration (1993), and that story belongs to post 4.

<img src="/assets/img/posts/ellipsometer-components/en/fig4-zero-vs-multi-order.png" alt="Zero order versus multi order" width="760">
_Figure 5. Both curves pass through 90° at 550 nm. Take the retardance modulo 360°, as in the lower panel, and the multi-order plate crosses 180° five times and 0° (360°) five times across the same band._

Figure 5 compares zero and multi order the same way. Retardance is $2\pi$-periodic, so what actually acts on the polarization is the value modulo 360°. In the lower panel the zero-order plate crosses 180° once across the band while the multi-order plate crosses five times. That both plates give exactly 90° at the 550 nm design wavelength only sharpens the contrast.

The second route abandons birefringence altogether. On total internal reflection the p and s components pick up different phase shifts, and the difference depends on the angle of incidence. A retarder built on that property alone is the Fresnel rhomb.

<img src="/assets/img/posts/ellipsometer-components/en/fig5-fresnel-rhomb.png" alt="Phase difference on total internal reflection" width="760">
_Figure 6. The p–s phase difference from a single total internal reflection. At an index of 1.51, entering at 42.3° or 74.7° gives 22.5° per bounce, and four bounces make a quarter wave._

Figure 6 is the calculation. In glass of index 1.51 the critical angle is 41.5°, and the phase difference rises to a maximum of 45.9° at 51.3°. Two angles give 22.5°: 42.3° and 74.7°. The double Fresnel rhomb uses the 74.7° branch, with four reflections summing to 90°, arranged so the exit beam lies on the same axis as the entrance beam. These angles were computed both from the phase difference of the Fresnel coefficients and from the closed-form expression, and the two agree.

What is gained? In Figure 6 the three curves for indices from 1.50 to 1.52 nearly coincide around 74.7° (four bounces give 89.5–90.5°); near the critical angle, at 42.3°, the same change swings the total from 70° to 106°, which is why that branch is not used. The refractive index depends only weakly on wavelength, so the retardance is nearly wavelength-independent as well, unlike a wave plate, which goes as $1/\lambda$. This is why an old method is still in use when a wide-band achromatic retarder is needed.

The third route is to make the retardance adjustable rather than fixed. The Babinet-Soleil compensator drives a wedge of birefringent material with a micrometer screw to change the thickness, so it can be set to 90° at whatever wavelength is wanted. It was used mainly in manual instruments before lasers.

The Berek compensator works differently. In an ordinary wave plate the optic axis lies in the plane of the crystal face; in a Berek compensator it is perpendicular to that face. Untilted, both polarizations see the same index and no phase shift occurs; tilt it and one component becomes partly extraordinary, producing a phase difference. Tilting the crystal tunes the retardance continuously.

Its use connects directly back to post 1. In rotating-analyzer instruments it serves to shift the operating point so that $\Delta$ does not sit near 0° or 180°. Section 4 of post 1 described measuring at several values of $\delta$ and stitching together only the stretches where $$\Delta' \sim \pm 90°$$. This is the hardware that carries out that strategy.

Choosing by material goes like this. Calcite, the workhorse of polarizers, is rarely used for compensators: its birefringence is so large that the required thickness becomes impractically small. A compensator that has to work in the ultraviolet is made of MgF$_2$ (transmitting above 0.12 μm); if nothing shorter than 0.29 μm is needed, mica will do. MgF$_2$ dominates now because of its ultraviolet transmission. Polarizers, too, must turn to MgF$_2$ Rochon prisms once they go below 0.21 μm.

## 4. Can the Retardance Change Without Moving Anything?

Every component so far has a retardance that is fixed once installed. Changing it means rotating or tilting. Can it be changed without moving anything mechanically?

An isotropic material under stress has its electron density pulled to one side and becomes anisotropic. The effect is photoelasticity, and it allows a time-varying birefringence in place of a stationary crystal. A device built on it is a photoelastic modulator (PEM).

The construction is simple. A quartz crystal cut in a particular orientation serves as a piezoelectric transducer, and a block of fused quartz is bonded to it with its length chosen to share the same resonant frequency. Driving the transducer at 50 kHz sets up a resonance whose stress is transmitted to the fused quartz, producing a periodic birefringence. The retardance then varies as

$$\delta(t) = F \sin \omega t$$

Here $F$ is the amplitude of the retardance and $\omega$ the drive angular frequency. Post 1 used this equation to describe PME (phase-modulation ellipsometry). What this post can add is what the amplitude $F$ depends on. It is proportional to the drive voltage $V$ and inversely proportional to the wavelength.

$$F \propto \frac{V}{\lambda}$$

At a fixed voltage, then, the retardance swings harder at shorter wavelengths.

Post 1 noted that PME handles an order of magnitude fewer wavelengths in real time, and this is the reason. Holding the retardance amplitude constant across wavelengths means re-setting the voltage at every wavelength. A rotating instrument is the opposite case: turning a single element serves all wavelengths at once. A photoelastic modulator is also highly sensitive to temperature, so accurate control of the retardance requires careful thermal regulation.

## 5. Why Would You Want a Part That Erases Polarization?

Why would an instrument built to create and analyze polarization contain a part whose job is to destroy it?

There are three reasons. The light leaving a source is not perfectly unpolarized; a slightly polarized component remains, known as source polarization. The diffraction efficiency of a grating spectrometer depends on polarization. And the sensitivity of a detector varies with the polarization state as well.

The last two cannot be distinguished from a polarization change caused by the sample. There is no way to tell whether the signal moved because of the sample or because the spectrometer responded differently to a different polarization. So a depolarizer is placed immediately before the spectrometer or detector, erasing the polarization and passing on the intensity alone.

The construction inverts the equation of section 1. In a wedge of birefringent crystal the thickness differs from point to point across the beam, and so does the retardance. Spread the retardance continuously over the beam cross section and the combined output is unpolarized. A Cornu prism, which exploits the optical activity of quartz, achieves the same end by another route.

## 6. How Are the Source, Detector and Spectrometer Chosen?

The classical source was a gas-discharge lamp. Mercury-arc lamps were used extensively, the strong 5461 Å mercury green line in particular. A discharge source also emits a continuum between its strong resonance lines, which makes it useful for spectroscopic work when paired with a monochromator. Lasers offer ideal monochromaticity and collimation, but the polarization state of the output must be stable, and unpolarized output is preferred, although linearly polarized output can easily be made circular with a suitable quarter-wave plate.

What do today's spectroscopic ellipsometers use? From the visible into the near ultraviolet, a high-pressure xenon lamp is the usual choice. Its output falls off steeply above 4 eV (310 nm), though, and the measurement noise on the ultraviolet side grows with it. Instruments that reach the deep ultraviolet around 6 eV (207 nm) therefore add a deuterium (D$_2$) lamp to the xenon lamp (Figure 7). In short: a laser for a single wavelength, a xenon lamp from the visible to the near ultraviolet, and a deuterium lamp on top of it to go deeper. An instrument that only covers 400–700 nm can get by with a white light-emitting diode (LED), as the reference ellipsometer in the author's dissertation does.

<img src="/assets/img/posts/ellipsometer-components/ext-deuterium-lamp.jpg" alt="A deuterium arc lamp with a metal shield box inside a glass bulb" width="380">
_Figure 7. A deuterium arc lamp. The discharge light comes out through a small hole in the metal shield inside the glass bulb, and it fills in the ultraviolet where a xenon lamp weakens. (Source: Latente Flickr, [Deuterium arc lamp](https://commons.wikimedia.org/wiki/File:Deuterium_arc_lamp.jpg), [CC BY-SA 2.0](https://creativecommons.org/licenses/by-sa/2.0))_

Detectors are chosen by spectral range. Photomultiplier tubes (PMTs) served visible ellipsometry for a long time and sometimes need cooling to stay stable; silicon photodiodes are widely used as well. But the spectral response curve is not the only thing to look at when choosing a detector. Its noise characteristic directly governs the precision of the measurement.

The previous sentence answers post 1 from the component side. Null ellipsometry never uses an absolute intensity, so in principle it is free of detector nonlinearity. How precisely the null can be located still depends on detector noise, which is why the Azzam–Bashara sentence above sits in their section on the null ellipsometer. Photometric methods trust the intensity reading and are fast in return, and the limit of that trust shows up as the noise characteristic of the detector.

Spectrometers moved from prisms to diffraction gratings. The $m$-th order of a grating interferes constructively at angles satisfying

$$\sin\phi + \sin\theta = m\lambda G$$

where $G$ is the number of grooves per unit length and $\phi$ and $\theta$ are the incidence and diffraction angles. What sets the angle is not the wavelength alone but the product $m\lambda$. A practical problem falls straight out of this: 800 nm in first order emerges at the same angle as 400 nm in second order. For a spectroscopic ellipsometer trying to cover the ultraviolet through the near infrared in one shot, this order overlap is a real constraint, and filters are needed to separate them.

Figure 8 shows the inside of a simple spectrometer built around a reflection grating.

<img src="/assets/img/posts/ellipsometer-components/ext-grating-spectrometer-inside.jpg" alt="White light split into colour bands inside an opened grating spectrometer" width="600">
_Figure 8. The inside of a simple grating spectrometer with its lid off. White light from the source on the left passes the mirrors and grating in the box and is split by wavelength, landing as colour bands on a part near the top and on the right-hand wall. The grating equation's message, that the same white light leaves in a different direction for each wavelength, is plain to see. (Source: Tigerzeng, [Simple grating spectrometer inside](https://commons.wikimedia.org/wiki/File:Simple_grating_spectrometer_inside.jpg), [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0))_

What choosing components actually involves is clearest in a concrete design. In the author's dissertation, a line-scan spectrometer was designed specifically to capture angle-resolved data from a back focal plane (BFP) in a single shot. With a 6.144 mm detector, a 400–700 nm band and a grating of 830 grooves per millimetre, the ideal focal lengths work out to 22.1 mm for the focusing lens and 24.6 mm for the collimating lens, and the spread of diffraction angles calls for a minimum numerical aperture (NA) of about 0.12–0.14. The nearest available lenses, 22.5 mm and 25 mm, were chosen as achromatic doublets to suppress chromatic and geometric aberration. The word achromatic, which sorted the methods in post 1, reappears here as the name of a lens.

At the slit the trade-off becomes numerical. A slit width of 40 μm gives a spectral resolution of about 1.9 nm. Narrowing it improves the resolution but admits less light. Meanwhile the slit must be longer than 3.8 mm to capture the full back focal plane of a 100× objective. Resolution, throughput and field of view are all settled by one slit.

## 7. How Do the Components Combine into the Methods of Post 1?

Lay out the components covered so far and the methods of post 1 fall out. Use only a polarizer and an analyzer and rotate the analyzer, and you have RAE. Insert a compensator between them and rotate it, and you have RCE (rotating-compensator ellipsometry). Put a photoelastic modulator in the compensator's place and you have PME.

The components also explain why null ellipsometry was, by and large, the standard for three quarters of the twentieth century. There was almost no need to quantify intensity, so the detector did not have to be good, and a single wavelength meant a simple quarter-wave plate cut for that one wavelength sufficed. It was a method whose components were allowed to be simple.

Moving to spectroscopy changes the terms. The spectroscopic condition can even reduce the number of components. Section 4 of post 1 covered the strategy of shifting the operating point with a compensator and combining several runs, yet many real spectroscopic ellipsometers carry no compensator at all. They simply discard the data at wavelengths where $\Delta$ falls near 0° or 180° and run the regression on what remains. With some two hundred wavelengths in hand, a few can be thrown away. For a single-wavelength instrument that choice does not exist.

There is a trick to the component angles as well. Setting the fixed polarizer or analyzer not at 45° but near $\Psi$ reduces the measurement uncertainty considerably. Moving it before each measurement based on the $\Psi$ from the previous one is called polarizer tracking. In post 1 the polarizer angle $P$ sat in the RAE relation as a free variable; this is what that freedom is for.

$$\tan\Psi = \sqrt{\frac{1+\alpha}{1-\alpha}} \, \lvert \tan P \rvert$$

Here $\alpha$ is the Fourier coefficient of the $\cos 2\omega t$ term in the RAE intensity of post 1, and $P$ is the azimuth of the fixed polarizer.

What links the components is a component too. An optical fibre lets the spectrometer sit somewhere convenient instead of riding on the goniometer arm. A fibre carries intensity through, however, but scrambles the polarization state. [Ellipsometry Foundations 2](/en/posts/ellipsometry-polarization-mueller-matrix/) said an intensity measurement loses phase, post 1 said the same thing repeats at the instrument level, and here it repeats once more in a single cable.

## 8. Components Do Not Arrive as Specified

One last point. Everything above assumed the retardance and axis orientation of each component are known. In practice they are not.

The birefringence quoted by a manufacturer carries an error. Precise work therefore starts by measuring the Mueller matrix of the component itself. In the author's dissertation, a separate dual-rotating-compensator ellipsometer was assembled to measure the Mueller matrices of the polarizer, the multi-order retarder and the beam splitter. Each matrix was decomposed by the Lu–Chipman (polar) decomposition into retarder, diattenuation and depolarization matrices, and the depolarization term, which is awkward to carry in the model, was dropped. What comes out is the actual retardance as a function of wavenumber, measured rather than assumed.

The component azimuths are not arbitrary either. There is an orientation that minimizes the condition number of the system of equations built from the measured quantities, and setting the components there makes the measurement less sensitive to noise. In the author's setup, with the polarizer at 0°, the compensator at 45° and the analyzer at 45°, the signals reduce to

$$\alpha' = \cos 2\Psi, \qquad \beta' = \sin\Delta \sin 2\Psi$$

Here $$\alpha'$$ and $$\beta'$$ are the coefficients of the $\cos\delta$ and $\sin\delta$ terms, in the retardance $\delta$ of the multi-order retarder, left after the spectral interference signal is divided by its DC part. The angles are measured from that setup's reference axis. In this arrangement $$\alpha'$$ carries $\Psi$ alone, and $\Delta$ appears only in $$\beta'$$.

Buying the component is not the end of it. Choosing it, orienting it, and measuring what it actually does before putting it in the model all belong together. That last step is calibration, and post 4 takes it up properly.

## Summary

This post makes one claim: the spectral range picks the components, and the components picked decide the method of the instrument.

That chain starts from one line, the retardance equation. Because $\delta$ goes as $1/\lambda$, a component that is perfect at one wavelength falls apart over a band. Where it falls apart is not arbitrary: it is at 0° and 180°, and there the information disappears. That is the component-level version of what post 1 called "what it cannot measure". So plates are stacked, total internal reflection is used, crystals are tilted, and stress is driven electrically. None of these is free: calibration grows harder, or the number of wavelengths shrinks.

Component by component the conclusions are simple. Going into the ultraviolet means no cement and MgF$_2$. Going into the infrared means abandoning prisms for a wire grid. Anything that will be rotated must not deviate the beam, which favours the Rochon prism. The photoelastic modulator, which needs its voltage re-set at every wavelength, is fast but gives up the number of wavelengths. For the source, going above 4 eV means adding a deuterium lamp to the xenon lamp.

In the language of post 1, then: an ellipsometer creates a variety of polarization states, sends them to the sample, and reads the states that come back. What those states can be made from is what fixes the wavelengths available.

[The next post](/en/posts/ellipsometer-dual-rotating-compensator-mme/) moves to using two of these parts at once. A compensator goes in both the polarization state generator and the polarization state analyzer, and the two are rotated at different speeds. Where post 1 obtained $S_1$ through $S_3$ in a single measurement, this becomes a question of obtaining the full Mueller matrix of the sample.

## References

- H. Fujiwara, *Spectroscopic Ellipsometry: Principles and Applications*, Wiley, 2007, section 3.2 (polarizers, compensators, photoelastic modulators and depolarizers).
- H. G. Tompkins and E. A. Irene (eds.), *Handbook of Ellipsometry*, William Andrew, 2005, chapter 4 (optical components and the simple PCSA ellipsometer).
- R. M. A. Azzam and N. M. Bashara, *Ellipsometry and Polarized Light*, North-Holland, 1977, section 5.2 (polarizing elements, light sources and detectors).
- B. Johs, J. Hale, N. J. Ianno, C. M. Herzinger, T. Tiwald, J. A. Woollam, "Recent developments in spectroscopic ellipsometry for in-situ applications," *Proc. SPIE* **4449**, 41–57 (2001), [doi:10.1117/12.450108](https://doi.org/10.1117/12.450108) — source of the two-element compensator.
- E. Hecht, *Optics*, 5th ed. (Pearson, 2017), section 8.4.3 (field of view of Glan-type polarizers).
- G. Ghosh, "Dispersion-equation coefficients for the refractive index and birefringence of calcite and quartz crystals," *Opt. Commun.* **163**, 95–102 (1999) — wavelength dispersion of the birefringence of quartz.
- B. Johs, "Regression calibration method for rotating element ellipsometers," *Thin Solid Films* **234**, 395–398 (1993) — making up for non-ideal components in calibration.
- Youngjoon Kim, "Snapshot Angle-Resolved Spectroscopic Ellipsometry Using Line-Scan Spectrometer and Back Focal Plane Spectral Interference," PhD thesis, Seoul National University, 2025, sections 4.4–4.5.

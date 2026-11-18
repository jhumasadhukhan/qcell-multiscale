### A. The reference state

All worked examples below are evaluated at the **apo (ligand-free) reference state in a polymer-free cytoplasm**, $\phi_P = 0$, unless stated. This is the state against which a conventional calibration is performed, so it is the natural origin for the error analysis. Values are taken from Supplementary Tables 4 and 5.

### B. Refractive index of the medium — equation (11)

The Lorentz–Lorenz (Clausius–Mossotti) mixing rule is volume-additive in the polarisability factor $f(n) = (n^2-1)/(n^2+2)$, not in $n$ itself:

$$f(n_{med}) = (1-\phi_P)\,f(n_c) + \phi_P\,f(n_P)$$

With $n_c = 1.36$ and $n_P = 1.49$,

$$f(n_c) = \frac{1.36^2-1}{1.36^2+2} = \frac{0.8496}{3.8496} = 0.22070,\qquad f(n_P) = 0.28912.$$

At $\phi_P = 0$ the mixture is pure cytosol, so $f = 0.22070$ and

$$n_{med} = \sqrt{\frac{1+2(0.22070)}{1-0.22070}} = 1.3600,$$

recovering $n_c$ exactly, as the rule requires. The mixing rule earns its place only once granules appear; §M below repeats the calculation in a granule-loaded cell.

### C. Orientation factor under partial averaging — equations (12)–(13)

**Step 1 — rotational correlation time (Stokes–Einstein–Debye).** At $\phi_P = 0$ the Krieger–Dougherty factor is unity, so $\eta = \eta_0 = 3.0$ mPa s. With an FP-barrel hydrodynamic volume $V_h = 8.000 \times 10^{-26}$ m³ at $T = 298.15$ K,

$$\theta_{rot} = \frac{\eta V_h}{k_B T} = \frac{(0.003)(8.000 \times 10^{-26})}{(1.381 \times 10^{-23})(298.15)} = 5.830 \times 10^{-8}\ \mathrm{s} = 58.3\ \mathrm{ns}.$$

**Step 2 — the dynamically averaged fraction.** The barrel tumbles on 58 ns while the donor lives only 4.0 ns, so most dipoles have *not* reorientated before transfer occurs:

$$f_{dyn} = \frac{\tau_D}{\theta_{rot}+\tau_D} = \frac{4.0}{58.3+4.0} = 0.0642.$$

Only 6.4 % of the ensemble is isotropically averaged.

**Step 3 — the effective orientation factor.** Mixing the averaged and frozen populations,

$$\kappa^2_{eff,apo} = \tfrac{2}{3}(0.0642) + 0.52(1-0.0642) = 0.5294,$$

$$\kappa^2_{eff,holo} = \tfrac{2}{3}(0.0642) + 1.15(1-0.0642) = 1.1190.$$

Both sit within 3 % of their static values and far from the isotropic 2/3. **This is the quantitative justification for not assuming $\kappa^2 = 2/3$**, an assumption that would inflate $R_0$ in the apo state by 3.9 % and distort the whole calibration.

### D. Förster radius — equation (5)

$$R_0 = 0.211\left[\kappa^2 n^{-4}\Phi_D J_{ov}\right]^{1/6}\ \mathrm{(Å,\ with\ }J_{ov}\mathrm{\ in\ M^{-1}cm^{-1}nm^4)}$$

Substituting the apo state:

$$\kappa^2 n^{-4}\Phi_D J_{ov} = (0.5294)(1.3600)^{-4}(0.93)(2.400 \times 10^{15}) = 3.454 \times 10^{14}$$

$$R_0 = 0.211\left[3.454 \times 10^{14}\right]^{1/6} = 0.211 \times 264.88 = 55.89\ \mathrm{Å}.$$

The reference value used to calibrate the overlap density (§E) is the isotropic one, $\kappa^2 = 2/3$ and $n = 1.33$, giving $R_0 = 58.95$ Å — consistent with published values of 5.8–5.9 nm for this donor–acceptor pair.

### E. Coulomb coupling and the golden-rule rate — equations (1)–(6)

**Step 1 — point-dipole coupling, equation (1).** In the units used here ($\mu$ in D, $R$ in Å) the prefactor $5034.1$ converts D²Å⁻³ to cm⁻¹:

$$J = 5034.1\,\frac{\kappa\,\mu_D\mu_A}{n^2 R^3}$$

With $\kappa = \sqrt{0.5294} = 0.7276$, $\mu_D = 4.37$ D, $\mu_A = 8.86$ D, $n = 1.3600$ and $R_{apo} = 68.0$ Å:

$$J = 5034.1\times\frac{(0.7276)(4.37)(8.86)}{(1.3600)^2(68.0)^3} = 5034.1\times\frac{28.172}{5.816 \times 10^{5}} = 0.2439\ \mathrm{cm^{-1}}.$$

Note that $\kappa$, not $\kappa^2$, enters here; the square appears only in rate expressions.

**Step 2 — the calibrated overlap density, equation (6).** $\mathcal{J}$ is fixed once, by requiring that the golden-rule rate equal the Förster rate at $R = R_0$. At the isotropic reference point ($R_0 = 58.95$ Å, $\kappa = \sqrt{2/3}$, $n = 1.33$) the coupling is

$$J(R_0) = 0.4392\ \mathrm{cm^{-1}},$$

and therefore

$$\mathcal{J} = \frac{1}{4\pi^2 c\,J(R_0)^2\,\tau_D} = \frac{1}{4\pi^2(2.998 \times 10^{10})(0.4392)^2(4.000 \times 10^{-9})} = 0.001095\ \mathrm{cm}.$$

The reciprocal, $1/\mathcal{J} \approx 913$ cm⁻¹, is an effective overlap bandwidth — a physically sensible width for two broad fluorescent-protein bands, and an internal consistency check on the layer.

**Step 3 — the transfer rate, equation (3).**

$$k_{EET} = 4\pi^2 c\,J^2\mathcal{J} = 4\pi^2(2.998 \times 10^{10})(0.2439)^2(0.001095) = 7.708 \times 10^{7}\ \mathrm{s^{-1}}.$$

**Cross-check against the Förster form.** The same rate must follow from $k = \tau_D^{-1}(R_0/R)^6$:

$$k = \frac{1}{4.000 \times 10^{-9}}\left(\frac{55.89}{68.0}\right)^6 = (2.500 \times 10^{8})(0.30830) = 7.708 \times 10^{7}\ \mathrm{s^{-1}},$$

agreeing with the golden-rule value to the last digit carried (relative difference 1.9 \times 10^{-16}). The two formalisms are numerically identical by construction, which is the point of the calibration in Step 2.

### F. Coherence correction — equation (8)

This term tests, rather than assumes, that transfer is incoherent hopping.

$$\omega_J = 2\pi c J = 2\pi(2.998 \times 10^{10})(0.2439) = 4.593 \times 10^{10}\ \mathrm{rad\ s^{-1}}$$

$$\gamma = 1/T_2^{\star} = 1/(6.000 \times 10^{-14}) = 1.667 \times 10^{13}\ \mathrm{s^{-1}}$$

$$\frac{2\omega_J}{\gamma} = 0.005512, \qquad k^{corr}_{EET} = \frac{k_{EET}}{1+(0.005512)^2} = k_{EET}\times\left(1 - 3.038 \times 10^{-5}\right).$$

The correction is 3.0e-05 in relative terms, about one part in 32,914. It is largest in the bound state, where the chromophores are closest and $J = 1.145$ cm⁻¹, giving 6.7e-04, one part in 1,492. Even that bound is three orders of magnitude below the quantities this paper reports. **The sensor is unambiguously in the Förster regime.** The term is retained in the code because it becomes material for strongly coupled pairs such as tandem-dimer fluorescent proteins, where $J$ reaches tens of cm⁻¹ and this ratio approaches unity.

### G. Transfer efficiency in both conformations — equation (7)

$$E = \frac{k_{EET}}{k_{EET}+\tau_D^{-1}} = \frac{1}{1+(R/R_0)^6}$$

**Apo state** ($R = 68.0$ Å):

$$E_{apo} = \frac{7.707 \times 10^{7}}{7.707 \times 10^{7}+2.500 \times 10^{8}} = 0.2356$$

**Holo state** ($R = 46.0$ Å, $\kappa^2_{eff} = 1.1190$, $R_0 = 63.32$ Å):

$$E_{holo} = \frac{1}{1+(46.0/63.32)^6} = 0.8717$$

The sensor therefore swings from $E = 0.236$ to $E = 0.872$ on binding — a dynamic range of 0.636 in efficiency. A 32 % contraction in distance produces a 270 % change in efficiency: this is the sixth-power amplification that makes FRET a usable ruler.

### H. Ligand occupancy and the conformational ensemble — equations (9)–(10)

At an intracellular 2-OG concentration $G = 1.0$ mM, with $K_D = 0.2$ mM and $h = 1.15$:

$$\theta = \frac{G^h}{G^h+K_D^h} = \frac{1.0^{1.15}}{1.0^{1.15}+0.2^{1.15}} = \frac{1.0000}{1.0000+0.1571} = 0.8642$$

Because conformational exchange (microseconds) is slow relative to the donor lifetime (nanoseconds), each molecule transfers from a *fixed* geometry and the ensemble averages **efficiencies, not geometries**:

$$E = (1-\theta)E_{apo} + \theta E_{holo} = (0.1358)(0.2356) + (0.8642)(0.8717) = 0.7854$$

Averaging the *distances* first and computing one efficiency from the mean geometry would give

$$\langle R\rangle = 48.99\ \mathrm{Å} \Rightarrow E(\langle R\rangle) = 0.8123,$$

an error of +3.4 %. The distinction is not pedantic: it is a direct consequence of the timescale separation, and getting it backwards biases every inferred concentration.

### I. The ratiometric observable — equation (19)

$$\mathcal{R} = \frac{I_A}{I_D} = \frac{\left[\Phi_A E + \alpha\Phi_A + b\,\Phi_D(1-E)\right]T_A + d}{\Phi_D(1-E)\,T_D + d}$$

In a cuvette with no turbidity ($T_D = T_A = 1$), at $E = 0.7854$:

$$I_D = \Phi_D(1-E) = (0.93)(1-0.7854) = 0.1996$$

$$I_A = \Phi_A E + \alpha\Phi_A + b\Phi_D(1-E)$$
$$I_A = (0.57)(0.7854) + (0.085)(0.57) + (0.055)(0.93)(1-0.7854)$$
$$I_A = 0.4477 + 0.0485 + 0.0110 = 0.5071$$

$$\mathcal{R} = \frac{0.5071+0.004}{0.1996+0.004} = 2.5100$$

Note that the spectroscopic cross-terms are not negligible: direct acceptor excitation and donor bleed-through together contribute 12 % of the acceptor channel. They are constants of the optical train, however, and a conventional calibration absorbs them correctly. The terms that follow do not behave so well.

### J. Band-resolved optical transport — equations (16)–(18)

**This is the step a conventional calibration omits, and it is where the paper's central result comes from.**

Take a mid-fermentation state: $C_X = 8.0$ g L⁻¹ total biomass of which $C_P = 3.2$ g L⁻¹ is polymer, so residual biomass is $C_R = C_X - C_P = 4.8$ g L⁻¹ and the polymer fraction is $f_P = 0.40$.

**Step 1 — attenuation coefficient at the reference wavelength.** Granules scatter *in addition to* the cell material around them:

$$\mu_{ref} = a_X C_R + (a_X + a_P)C_P$$
$$\mu_{ref} = (0.25)(4.8) + (0.25+0.12)(3.2) = 1.2000 + 1.1840 = 2.3840\ \mathrm{cm^{-1}}$$

**Step 2 — wavelength scaling.** Scattering falls with wavelength as $\lambda^{-p}$ with $p = 1.5$, referenced to $\lambda_{ref} = 500$ nm, over an effective path $L_{eff} = 0.3$ cm:

$$A_D = \mu_{ref}\left(\frac{500}{480}\right)^{1.5}L_{eff} = (2.3840)(1.0631)(0.3) = 0.7604$$

$$A_A = (2.3840)(0.9163)(0.3) = 0.6553$$

**Step 3 — transmissions.**

$$T_D = e^{-0.7604} = 0.4675, \qquad T_A = e^{-0.6553} = 0.5193$$

$$\boxed{\ \frac{T_A}{T_D} = \frac{0.5193}{0.4675} = 1.1107\ }$$

**The ratio is inflated by 11.1 % by transport alone, before the sensor has reported anything.** Three features of this result deserve emphasis.

1. **The sign is fixed.** Because $\lambda_D < \lambda_A$ always — the donor emits to the blue of the acceptor, by construction of any FRET pair — $A_D > A_A$ and hence $T_A/T_D > 1$ *always*. The bias can never be negative. It is a one-sided error, not noise.

2. **It is not removable by optical-density normalisation.** The quantity $\mu_{ref}$ depends on $C_P$ and $C_R$ separately, through different coefficients. An OD probe returns only their weighted sum. §K shows what that costs.

3. **It grows as the process succeeds.** $C_P$ rises monotonically through a production run, so the measurement degrades exactly when it matters most.

### K. Why correcting for cell density is not enough

A conventional in-situ calibration normalises by optical density, which in effect assumes all attenuation is proportional to total biomass:

$$\mu_{ref}^{OD} = a_X C_X = (0.25)(8.0) = 2.0000\ \mathrm{cm^{-1}} \quad\text{versus the true }2.3840\ \mathrm{cm^{-1}}.$$

That correction removes a drift of 9.21 %, leaving a residual of

$$\frac{T_A/T_D}{(T_A/T_D)^{OD}} = \frac{1.1107}{1.0921} = 1.0171,$$

i.e. 1.71 % of uncorrected, polymer-specific inflation. The missing term is $a_P C_P$ — and $C_P$ is precisely the quantity no non-destructive bulk measurement returns. **It has to come from the cellular model. This is why the two layers cannot be decoupled.**

### L. Propagation to the inferred metabolite concentration

Finally, the error that matters operationally. The true state is $G = 1.00$ mM at the biomass composition of §J ($C_X = 8.0$, $C_P = 3.2$ g L⁻¹, $\phi_P = 0.138$). The forward map returns an observed ratio $\mathcal{R} = 2.6666$, against $\mathcal{R} = 2.5100$ for the same $G$ in a clear cuvette. Inverting under the three correction levels gives:

| Correction applied | Inferred $G$ (mM) | Error |
|---|---|---|
| None (fixed calibration curve) | 1.1399 | +14.0 % |
| Optical density only | 0.9671 | -3.3 % |
| Full multiscale correction | 1.0000 | -0.00 % |

The full correction recovers the true value to within numerical tolerance, confirming that the inverse map is exact. The uncorrected error is **amplified relative to the optical bias itself**: a 11.1 % inflation of the ratio becomes a +14.0 % error in $G$, because the Hill function compresses the ratio axis near saturation and the inversion therefore expands it.

**How this relates to the figure in the main text.** The +14.0 % above is the error at one static operating point, chosen to make the arithmetic followable. The main text instead quotes +29 %, which is the same quantity averaged over a closed-loop run: the controller drives the culture to higher polymer content than this reference state, the bias grows with polymer content, and the average is taken over the production phase rather than at a single point. The two numbers are consistent, and Fig. 2c of the main text plots the whole curve from which both are read.

**A caution on the middle row.** The OD-only correction looks respectable at -3.3 %, but that figure is the near-cancellation of two larger errors of opposite sign. Isolating them — correcting one channel properly while leaving the other as OD-only does:

| Channel left uncorrected | Inferred $G$ (mM) | Error |
|---|---|---|
| Transport under-correction only (photophysics right) | 1.0302 | +3.02 % |
| Photophysics assumed at $\phi_P = 0$ only (transport right) | 0.9399 | -6.01 % |
| **Both, i.e. OD-only** | **0.9671** | **-3.29 %** |

The two channels are close to additive (+3.02 % and -6.01 % sum to -2.98 % against an actual -3.29 %), and they **pull in opposite directions**. Leaving the granule term out of the attenuation under-corrects the transport, biasing $G$ upward; assuming a polymer-free cytoplasm puts the $\kappa^2$ and refractive-index calibration in the wrong place, biasing it down about twice as hard. The net is a modest over-correction.

Two consequences follow, and they matter more than the headline number. First, **the apparent adequacy of an OD-only correction is an artefact of this particular parameter set**: the balance between the two channels depends on $a_P$, on the granule volume fraction and on the static orientation factors, none of which is independently measured (§3.1). A modest change in any of them breaks the cancellation. Second, the cancellation is **not conservative** — it hides two errors of roughly 6 % behind a reported one of 3 %, so an experimenter who validated an OD-corrected sensor at one operating point could find it badly wrong at another with no warning. This is the reasoning behind giving the dimensionless criterion $\Pi$ of §6.4 rather than a blanket recommendation.

### M. The same calculation in a granule-loaded cell

Repeating §B–§D at $\phi_P = 0.3$, a realistic late-stage granule volume fraction:

| Quantity | $\phi_P = 0$ | $\phi_P = 0.3$ | Change |
|---|---|---|---|
| $n_{med}$ | 1.3600 | 1.3978 | +2.78 % |
| $\eta$ (mPa s) | 3.00 | 8.25 | +175.1 % |
| $\theta_{rot}$ (ns) | 58.3 | 160.4 | +175.1 % |
| $f_{dyn}$ | 0.0642 | 0.0243 | -62.1 % |
| $\kappa^2_{eff,apo}$ | 0.5294 | 0.5236 | -1.10 % |
| $R_0$ (Å) | 55.89 | 54.78 | -1.99 % |

The refractive index rises by 2.78 %, and because $R_0 \propto n^{-2/3}$ this alone shortens the Förster radius. The viscosity rise further freezes the orientation factor. Neither effect is large individually, but both act in the same direction, and both are invisible to an instrument that sees only the emitted ratio.

**The key asymmetry:** these photophysical shifts are of order 2.0 %, whereas the optical-transport bias of §J is 11.1 %. The transport term dominates by roughly an order of magnitude. A reader who takes only one number from this Supplementary Note should take that one: **the dominant quantum-to-cellular coupling is not exotic photophysics inside the cell, it is ordinary scattering on the way out.**

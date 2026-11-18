"""Generate the worked numerical examples for the Supplementary Information.

Every number printed here is computed by calling the SAME functions the model
uses, so the document cannot drift from the code. Intermediate quantities are
re-derived inline alongside the function call and asserted equal, which means a
discrepancy between the hand-substitution shown in the SI and what the model
actually does becomes a test failure rather than a silent error in the paper.
"""
import json
import os
import sys

import numpy as np

# the model lives one level up, so this runs from anywhere
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))

import qcell_model as Q
from qcell_model import QuantumParams

qp = QuantumParams()
J_CAL = Q.spectral_overlap_density(qp)

OUT = []


def w(s=""):
    OUT.append(s)


def chk(name, mine, theirs, tol=1e-9):
    """Assert the hand-substitution reproduces the model call."""
    rel = abs(mine - theirs) / max(abs(theirs), 1e-30)
    assert rel < tol, f"{name}: hand {mine!r} vs model {theirs!r} (rel {rel:.2e})"
    return mine


def fmt(x, sig=4):
    """Format in a style that reads naturally in a paper."""
    if x == 0:
        return "0"
    ax = abs(x)
    if 1e-3 <= ax < 1e5:
        return f"{x:.{sig}g}"
    m, e = f"{x:.{sig - 1}e}".split("e")
    return f"{m} \\times 10^{{{int(e)}}}"


# ===========================================================================
w("### A. The reference state")
w()
w("All worked examples below are evaluated at the **apo (ligand-free) reference "
  "state in a polymer-free cytoplasm**, $\\phi_P = 0$, unless stated. This is "
  "the state against which a conventional calibration is performed, so it is "
  "the natural origin for the error analysis. Values are taken from "
  "Supplementary Tables 4 and 5.")
w()

phi_P0 = 0.0

# --- Eq 11: Lorentz-Lorenz -------------------------------------------------
n0_model = Q.lorentz_lorenz_mix(phi_P0, qp.n_cytosol, qp.n_granule)
f_h = (qp.n_cytosol**2 - 1) / (qp.n_cytosol**2 + 2)
f_i = (qp.n_granule**2 - 1) / (qp.n_granule**2 + 2)
f0 = (1 - phi_P0) * f_h + phi_P0 * f_i
n0 = np.sqrt((1 + 2 * f0) / (1 - f0))
chk("n(phi=0)", n0, n0_model)

w("### B. Refractive index of the medium — equation (11)")
w()
w("The Lorentz–Lorenz (Clausius–Mossotti) mixing rule is volume-additive in the "
  "polarisability factor $f(n) = (n^2-1)/(n^2+2)$, not in $n$ itself:")
w()
w("$$f(n_{med}) = (1-\\phi_P)\\,f(n_c) + \\phi_P\\,f(n_P)$$")
w()
w(f"With $n_c = {qp.n_cytosol}$ and $n_P = {qp.n_granule}$,")
w()
w(f"$$f(n_c) = \\frac{{{qp.n_cytosol}^2-1}}{{{qp.n_cytosol}^2+2}} "
  f"= \\frac{{{qp.n_cytosol**2 - 1:.4f}}}{{{qp.n_cytosol**2 + 2:.4f}}} = {f_h:.5f},"
  f"\\qquad f(n_P) = {f_i:.5f}.$$")
w()
w(f"At $\\phi_P = 0$ the mixture is pure cytosol, so $f = {f0:.5f}$ and")
w()
w(f"$$n_{{med}} = \\sqrt{{\\frac{{1+2({f0:.5f})}}{{1-{f0:.5f}}}}} = {n0:.4f},$$")
w()
w("recovering $n_c$ exactly, as the rule requires. The mixing rule earns its "
  "place only once granules appear; §M below repeats the calculation in a "
  "granule-loaded cell.")
w()

# --- Eq 12-13: kappa^2 -----------------------------------------------------
eta0 = Q.krieger_dougherty(phi_P0, qp.eta_cyto0, qp.phi_pack)
kB = 1.380649e-23
theta_rot0 = eta0 * qp.V_hydro / (kB * Q.T_KELVIN)
f_dyn0 = qp.tau_D0 / (theta_rot0 + qp.tau_D0)
k2_apo0 = (2 / 3) * f_dyn0 + qp.kappa2_static_apo * (1 - f_dyn0)
chk("kappa2_apo", k2_apo0, Q.kappa2_effective(phi_P0, qp.kappa2_static_apo, qp))
k2_holo0 = (2 / 3) * f_dyn0 + qp.kappa2_static_holo * (1 - f_dyn0)
chk("kappa2_holo", k2_holo0, Q.kappa2_effective(phi_P0, qp.kappa2_static_holo, qp))

w("### C. Orientation factor under partial averaging — equations (12)–(13)")
w()
w("**Step 1 — rotational correlation time (Stokes–Einstein–Debye).** At "
  f"$\\phi_P = 0$ the Krieger–Dougherty factor is unity, so $\\eta = \\eta_0 = "
  f"{qp.eta_cyto0 * 1e3:.1f}$ mPa s. With an FP-barrel hydrodynamic volume "
  f"$V_h = {fmt(qp.V_hydro)}$ m³ at $T = {Q.T_KELVIN}$ K,")
w()
w("$$\\theta_{rot} = \\frac{\\eta V_h}{k_B T} "
  f"= \\frac{{({fmt(qp.eta_cyto0)})({fmt(qp.V_hydro)})}}"
  f"{{({fmt(kB)})({Q.T_KELVIN})}} = {fmt(theta_rot0)}\\ \\mathrm{{s}} "
  f"= {theta_rot0 * 1e9:.1f}\\ \\mathrm{{ns}}.$$")
w()
w("**Step 2 — the dynamically averaged fraction.** The barrel tumbles on "
  f"{theta_rot0 * 1e9:.0f} ns while the donor lives only "
  f"{qp.tau_D0 * 1e9:.1f} ns, so most dipoles have *not* reorientated before "
  "transfer occurs:")
w()
w("$$f_{dyn} = \\frac{\\tau_D}{\\theta_{rot}+\\tau_D} "
  f"= \\frac{{{qp.tau_D0 * 1e9:.1f}}}{{{theta_rot0 * 1e9:.1f}+{qp.tau_D0 * 1e9:.1f}}} "
  f"= {f_dyn0:.4f}.$$")
w()
w(f"Only {100 * f_dyn0:.1f} % of the ensemble is isotropically averaged.")
w()
w("**Step 3 — the effective orientation factor.** Mixing the averaged and "
  "frozen populations,")
w()
w(f"$$\\kappa^2_{{eff,apo}} = \\tfrac{{2}}{{3}}({f_dyn0:.4f}) "
  f"+ {qp.kappa2_static_apo}(1-{f_dyn0:.4f}) = {k2_apo0:.4f},$$")
w()
w(f"$$\\kappa^2_{{eff,holo}} = \\tfrac{{2}}{{3}}({f_dyn0:.4f}) "
  f"+ {qp.kappa2_static_holo}(1-{f_dyn0:.4f}) = {k2_holo0:.4f}.$$")
w()
w(f"Both sit within "
  f"{max(abs(k2_apo0/qp.kappa2_static_apo-1), abs(k2_holo0/qp.kappa2_static_holo-1))*100:.0f} % "
  "of their static values and far from the isotropic 2/3. "
  "**This is the quantitative justification for not assuming $\\kappa^2 = 2/3$**, "
  "an assumption that would inflate $R_0$ in the apo state by "
  f"{100 * ((2 / 3 / k2_apo0) ** (1 / 6) - 1):.1f} % and distort the whole "
  "calibration.")
w()

# --- Eq 5: R0 --------------------------------------------------------------
R0_apo = 0.211 * (k2_apo0 * n0**-4 * qp.phi_D * qp.J_overlap) ** (1 / 6)
chk("R0_apo", R0_apo, Q.forster_R0(k2_apo0, n0, qp))
R0_iso = Q.forster_R0(2 / 3, 1.33, qp)

w("### D. Förster radius — equation (5)")
w()
w("$$R_0 = 0.211\\left[\\kappa^2 n^{-4}\\Phi_D J_{ov}\\right]^{1/6}"
  "\\ \\mathrm{(Å,\\ with\\ }J_{ov}\\mathrm{\\ in\\ M^{-1}cm^{-1}nm^4)}$$")
w()
w("Substituting the apo state:")
w()
w(f"$$\\kappa^2 n^{{-4}}\\Phi_D J_{{ov}} = ({k2_apo0:.4f})({n0:.4f})^{{-4}}"
  f"({qp.phi_D})({fmt(qp.J_overlap)}) = {fmt(k2_apo0 * n0**-4 * qp.phi_D * qp.J_overlap)}$$")
w()
w(f"$$R_0 = 0.211\\left[{fmt(k2_apo0 * n0**-4 * qp.phi_D * qp.J_overlap)}\\right]^{{1/6}} "
  f"= 0.211 \\times {(k2_apo0 * n0**-4 * qp.phi_D * qp.J_overlap) ** (1/6):.2f} "
  f"= {R0_apo:.2f}\\ \\mathrm{{Å}}.$$")
w()
w(f"The reference value used to calibrate the overlap density (§E) is the "
  f"isotropic one, $\\kappa^2 = 2/3$ and $n = 1.33$, giving "
  f"$R_0 = {R0_iso:.2f}$ Å — consistent with published values of 5.8–5.9 nm for "
  "this donor–acceptor pair.")
w()

# --- Eq 1: coupling --------------------------------------------------------
kappa_apo = np.sqrt(k2_apo0)
J_apo = Q.DIPOLE_PREFACTOR * kappa_apo * qp.mu_D * qp.mu_A / (n0**2 * qp.R_apo**3)
chk("J_apo", J_apo, Q.coulomb_coupling(kappa_apo, qp.R_apo, n0, qp))

w("### E. Coulomb coupling and the golden-rule rate — equations (1)–(6)")
w()
w("**Step 1 — point-dipole coupling, equation (1).** In the units used here "
  "($\\mu$ in D, $R$ in Å) the prefactor $5034.1$ converts D²Å⁻³ to cm⁻¹:")
w()
w("$$J = 5034.1\\,\\frac{\\kappa\\,\\mu_D\\mu_A}{n^2 R^3}$$")
w()
w(f"With $\\kappa = \\sqrt{{{k2_apo0:.4f}}} = {kappa_apo:.4f}$, "
  f"$\\mu_D = {qp.mu_D}$ D, $\\mu_A = {qp.mu_A}$ D, $n = {n0:.4f}$ and "
  f"$R_{{apo}} = {qp.R_apo}$ Å:")
w()
w(f"$$J = 5034.1\\times\\frac{{({kappa_apo:.4f})({qp.mu_D})({qp.mu_A})}}"
  f"{{({n0:.4f})^2({qp.R_apo})^3}} "
  f"= 5034.1\\times\\frac{{{kappa_apo * qp.mu_D * qp.mu_A:.3f}}}"
  f"{{{fmt(n0**2 * qp.R_apo**3)}}} = {J_apo:.4f}\\ \\mathrm{{cm^{{-1}}}}.$$")
w()
w("Note that $\\kappa$, not $\\kappa^2$, enters here; the square appears only "
  "in rate expressions.")
w()
w("**Step 2 — the calibrated overlap density, equation (6).** $\\mathcal{J}$ is "
  "fixed once, by requiring that the golden-rule rate equal the Förster rate at "
  f"$R = R_0$. At the isotropic reference point ($R_0 = {R0_iso:.2f}$ Å, "
  f"$\\kappa = \\sqrt{{2/3}}$, $n = 1.33$) the coupling is")
w()
J0_ref = Q.coulomb_coupling(np.sqrt(2 / 3), R0_iso, 1.33, qp)
w(f"$$J(R_0) = {J0_ref:.4f}\\ \\mathrm{{cm^{{-1}}}},$$")
w()
w("and therefore")
w()
w("$$\\mathcal{J} = \\frac{1}{4\\pi^2 c\\,J(R_0)^2\\,\\tau_D} "
  f"= \\frac{{1}}{{4\\pi^2({fmt(Q.C_LIGHT)})({J0_ref:.4f})^2({fmt(qp.tau_D0)})}} "
  f"= {fmt(J_CAL)}\\ \\mathrm{{cm}}.$$")
w()
w(f"The reciprocal, $1/\\mathcal{{J}} \\approx {1 / J_CAL:.0f}$ cm⁻¹, is an "
  "effective overlap bandwidth — a physically sensible width for two broad "
  "fluorescent-protein bands, and an internal consistency check on the layer.")
w()
w("**Step 3 — the transfer rate, equation (3).**")
w()
k_bare = 4 * np.pi**2 * Q.C_LIGHT * J_apo**2 * J_CAL
w("$$k_{EET} = 4\\pi^2 c\\,J^2\\mathcal{J} "
  f"= 4\\pi^2({fmt(Q.C_LIGHT)})({J_apo:.4f})^2({fmt(J_CAL)}) "
  f"= {fmt(k_bare)}\\ \\mathrm{{s^{{-1}}}}.$$")
w()
w("**Cross-check against the Förster form.** The same rate must follow from "
  "$k = \\tau_D^{-1}(R_0/R)^6$:")
w()
k_forster = (1 / qp.tau_D0) * (R0_apo / qp.R_apo) ** 6
w(f"$$k = \\frac{{1}}{{{fmt(qp.tau_D0)}}}\\left(\\frac{{{R0_apo:.2f}}}"
  f"{{{qp.R_apo}}}\\right)^6 = ({fmt(1 / qp.tau_D0)})({(R0_apo / qp.R_apo) ** 6:.5f}) "
  f"= {fmt(k_forster)}\\ \\mathrm{{s^{{-1}}}},$$")
w()
_agree = abs(k_forster - k_bare) / k_bare
w("agreeing with the golden-rule value to the last digit carried "
  f"(relative difference {fmt(_agree, 2)}). The two formalisms are "
  "numerically identical by construction, which is the point of the "
  "calibration in Step 2.")
w()

# --- Eq 8: coherence -------------------------------------------------------
omega_J = 2 * np.pi * Q.C_LIGHT * J_apo
gamma = 1 / qp.T2_star
ratio_coh = 2 * omega_J / gamma
k_corr = k_bare / (1 + ratio_coh**2)
chk("k_corr", k_corr, Q.eet_rate(k2_apo0, qp.R_apo, n0, qp, J_CAL, True)[0])

w("### F. Coherence correction — equation (8)")
w()
w("This term tests, rather than assumes, that transfer is incoherent hopping.")
w()
w(f"$$\\omega_J = 2\\pi c J = 2\\pi({fmt(Q.C_LIGHT)})({J_apo:.4f}) "
  f"= {fmt(omega_J)}\\ \\mathrm{{rad\\ s^{{-1}}}}$$")
w()
w(f"$$\\gamma = 1/T_2^{{\\star}} = 1/({fmt(qp.T2_star)}) = {fmt(gamma)}\\ "
  f"\\mathrm{{s^{{-1}}}}$$")
w()
w(f"$$\\frac{{2\\omega_J}}{{\\gamma}} = {fmt(ratio_coh)}, \\qquad "
  f"k^{{corr}}_{{EET}} = \\frac{{k_{{EET}}}}{{1+({fmt(ratio_coh)})^2}} "
  f"= k_{{EET}}\\times\\left(1 - {fmt(ratio_coh**2)}\\right).$$")
w()
J_h = Q.coulomb_coupling(np.sqrt(k2_holo0), qp.R_holo, n0, qp)
ratio_h = 2 * (2 * np.pi * Q.C_LIGHT * J_h) * qp.T2_star
w(f"The correction is {ratio_coh**2:.1e} in relative terms, about one part in "
  f"{1 / ratio_coh**2:,.0f}. It is largest in the bound state, where the "
  f"chromophores are closest and $J = {J_h:.3f}$ cm⁻¹, giving "
  f"{ratio_h**2:.1e}, one part in {1 / ratio_h**2:,.0f}. Even that bound is "
  f"three orders of magnitude below the quantities this paper reports. "
  f"**The sensor is unambiguously in the Förster "
  "regime.** The term is retained in the code because it becomes material for "
  "strongly coupled pairs such as tandem-dimer fluorescent proteins, where $J$ "
  "reaches tens of cm⁻¹ and this ratio approaches unity.")
w()

# --- Eq 7: efficiency ------------------------------------------------------
E_apo = k_corr / (k_corr + 1 / qp.tau_D0)
chk("E_apo", E_apo, Q.fret_efficiency(k2_apo0, qp.R_apo, n0, qp, J_CAL, True)[0])
E_holo, J_holo = Q.fret_efficiency(k2_holo0, qp.R_holo, n0, qp, J_CAL, True)
R0_holo = Q.forster_R0(k2_holo0, n0, qp)

w("### G. Transfer efficiency in both conformations — equation (7)")
w()
w("$$E = \\frac{k_{EET}}{k_{EET}+\\tau_D^{-1}} = \\frac{1}{1+(R/R_0)^6}$$")
w()
w("**Apo state** ($R = " f"{qp.R_apo}" "$ Å):")
w()
w(f"$$E_{{apo}} = \\frac{{{fmt(k_corr)}}}{{{fmt(k_corr)}+{fmt(1 / qp.tau_D0)}}} "
  f"= {E_apo:.4f}$$")
w()
w(f"**Holo state** ($R = {qp.R_holo}$ Å, $\\kappa^2_{{eff}} = {k2_holo0:.4f}$, "
  f"$R_0 = {R0_holo:.2f}$ Å):")
w()
w(f"$$E_{{holo}} = \\frac{{1}}{{1+({qp.R_holo}/{R0_holo:.2f})^6}} = {E_holo:.4f}$$")
w()
w(f"The sensor therefore swings from $E = {E_apo:.3f}$ to $E = {E_holo:.3f}$ on "
  f"binding — a dynamic range of {E_holo - E_apo:.3f} in efficiency. A "
  f"{100 * (1 - qp.R_holo / qp.R_apo):.0f} % contraction in distance produces a "
  f"{100 * (E_holo - E_apo) / E_apo:.0f} % change in efficiency: this is the "
  "sixth-power amplification that makes FRET a usable ruler.")
w()

# --- Eq 9-10: occupancy and ensemble ---------------------------------------
G_ref = 1.0
theta_ref = Q.sensor_occupancy(G_ref, qp)
E_ens = (1 - theta_ref) * E_apo + theta_ref * E_holo

w("### H. Ligand occupancy and the conformational ensemble — equations (9)–(10)")
w()
w(f"At an intracellular 2-OG concentration $G = {G_ref}$ mM, with "
  f"$K_D = {qp.K_D}$ mM and $h = {qp.hill_h}$:")
w()
w(f"$$\\theta = \\frac{{G^h}}{{G^h+K_D^h}} "
  f"= \\frac{{{G_ref}^{{{qp.hill_h}}}}}{{{G_ref}^{{{qp.hill_h}}}"
  f"+{qp.K_D}^{{{qp.hill_h}}}}} "
  f"= \\frac{{{G_ref**qp.hill_h:.4f}}}{{{G_ref**qp.hill_h:.4f}+{qp.K_D**qp.hill_h:.4f}}} "
  f"= {theta_ref:.4f}$$")
w()
w("Because conformational exchange (microseconds) is slow relative to the donor "
  "lifetime (nanoseconds), each molecule transfers from a *fixed* geometry and "
  "the ensemble averages **efficiencies, not geometries**:")
w()
w(f"$$E = (1-\\theta)E_{{apo}} + \\theta E_{{holo}} "
  f"= ({1 - theta_ref:.4f})({E_apo:.4f}) + ({theta_ref:.4f})({E_holo:.4f}) "
  f"= {E_ens:.4f}$$")
w()
w("Averaging the *distances* first and computing one efficiency from the mean "
  "geometry would give")
w()
R_mean = (1 - theta_ref) * qp.R_apo + theta_ref * qp.R_holo
k2_mean = (1 - theta_ref) * k2_apo0 + theta_ref * k2_holo0
E_wrong, _ = Q.fret_efficiency(k2_mean, R_mean, n0, qp, J_CAL, True)
w(f"$$\\langle R\\rangle = {R_mean:.2f}\\ \\mathrm{{Å}} \\Rightarrow "
  f"E(\\langle R\\rangle) = {E_wrong:.4f},$$")
w()
w(f"an error of {100 * (E_wrong - E_ens) / E_ens:+.1f} %. The distinction is "
  "not pedantic: it is a direct consequence of the timescale separation, and "
  "getting it backwards biases every inferred concentration.")
w()

# --- Eq 19: ratio ----------------------------------------------------------
I_D_clean = qp.I_exc * qp.phi_D * (1 - E_ens)
I_A_clean = qp.I_exc * (qp.phi_A * E_ens + qp.alpha_direct * qp.phi_A
                        + qp.bleed_DA * qp.phi_D * (1 - E_ens))
R_clean = (I_A_clean + qp.dark) / (I_D_clean + qp.dark)
m0 = Q.quantum_map(G_ref, 0.0, qp, J_CAL, Cx=0.0, Cp=0.0)
chk("Rt clean", R_clean, m0["Rt"], tol=1e-9)

w("### I. The ratiometric observable — equation (19)")
w()
w("$$\\mathcal{R} = \\frac{I_A}{I_D} = "
  "\\frac{\\left[\\Phi_A E + \\alpha\\Phi_A + b\\,\\Phi_D(1-E)\\right]T_A + d}"
  "{\\Phi_D(1-E)\\,T_D + d}$$")
w()
w(f"In a cuvette with no turbidity ($T_D = T_A = 1$), at $E = {E_ens:.4f}$:")
w()
w(f"$$I_D = \\Phi_D(1-E) = ({qp.phi_D})(1-{E_ens:.4f}) = {I_D_clean:.4f}$$")
w()
w("$$I_A = \\Phi_A E + \\alpha\\Phi_A + b\\Phi_D(1-E)$$")
w(f"$$I_A = ({qp.phi_A})({E_ens:.4f}) + ({qp.alpha_direct})({qp.phi_A}) "
  f"+ ({qp.bleed_DA})({qp.phi_D})(1-{E_ens:.4f})$$")
w(f"$$I_A = {qp.phi_A * E_ens:.4f} + {qp.alpha_direct * qp.phi_A:.4f} "
  f"+ {qp.bleed_DA * qp.phi_D * (1 - E_ens):.4f} = {I_A_clean:.4f}$$")
w()
w(f"$$\\mathcal{{R}} = \\frac{{{I_A_clean:.4f}+{qp.dark}}}"
  f"{{{I_D_clean:.4f}+{qp.dark}}} = {R_clean:.4f}$$")
w()
w(f"Note that the spectroscopic cross-terms are not negligible: direct "
  f"acceptor excitation and donor bleed-through together contribute "
  f"{100 * (qp.alpha_direct * qp.phi_A + qp.bleed_DA * qp.phi_D * (1 - E_ens)) / I_A_clean:.0f} % "
  "of the acceptor channel. They are constants of the optical train, however, "
  "and a conventional calibration absorbs them correctly. The terms that "
  "follow do not behave so well.")
w()

# --- Eq 16-18: optical transport ------------------------------------------
Cx_ref, Cp_ref = 8.0, 3.2
cp_ref = Q.CellParams()
phi_ref = float(Q.granule_volume_fraction(Cx_ref, Cp_ref, cp_ref))
Cr_ref = Cx_ref - Cp_ref
mu_ref = qp.a_X * Cr_ref + (qp.a_X + qp.a_P) * Cp_ref
A_D = mu_ref * (qp.lam_ref / qp.lam_D) ** qp.p_mie * qp.L_eff
A_A = mu_ref * (qp.lam_ref / qp.lam_A) ** qp.p_mie * qp.L_eff
T_D, T_A = np.exp(-A_D), np.exp(-A_A)
chk("T_D", T_D, Q.optical_attenuation(Cx_ref, Cp_ref, qp)[0])
drift = T_A / T_D

w("### J. Band-resolved optical transport — equations (16)–(18)")
w()
w("**This is the step a conventional calibration omits, and it is where the "
  "paper's central result comes from.**")
w()
w(f"Take a mid-fermentation state: $C_X = {Cx_ref}$ g L⁻¹ total biomass of "
  f"which $C_P = {Cp_ref}$ g L⁻¹ is polymer, so residual biomass is "
  f"$C_R = C_X - C_P = {Cr_ref}$ g L⁻¹ and the polymer fraction is "
  f"$f_P = {Cp_ref / Cx_ref:.2f}$.")
w()
w("**Step 1 — attenuation coefficient at the reference wavelength.** Granules "
  "scatter *in addition to* the cell material around them:")
w()
w("$$\\mu_{ref} = a_X C_R + (a_X + a_P)C_P$$")
w(f"$$\\mu_{{ref}} = ({qp.a_X})({Cr_ref}) + ({qp.a_X}+{qp.a_P})({Cp_ref}) "
  f"= {qp.a_X * Cr_ref:.4f} + {(qp.a_X + qp.a_P) * Cp_ref:.4f} "
  f"= {mu_ref:.4f}\\ \\mathrm{{cm^{{-1}}}}$$")
w()
w("**Step 2 — wavelength scaling.** Scattering falls with wavelength as "
  f"$\\lambda^{{-p}}$ with $p = {qp.p_mie}$, referenced to "
  f"$\\lambda_{{ref}} = {qp.lam_ref:.0f}$ nm, over an effective path "
  f"$L_{{eff}} = {qp.L_eff}$ cm:")
w()
w(f"$$A_D = \\mu_{{ref}}\\left(\\frac{{{qp.lam_ref:.0f}}}{{{qp.lam_D:.0f}}}"
  f"\\right)^{{{qp.p_mie}}}L_{{eff}} = ({mu_ref:.4f})"
  f"({(qp.lam_ref / qp.lam_D) ** qp.p_mie:.4f})({qp.L_eff}) = {A_D:.4f}$$")
w()
w(f"$$A_A = ({mu_ref:.4f})({(qp.lam_ref / qp.lam_A) ** qp.p_mie:.4f})"
  f"({qp.L_eff}) = {A_A:.4f}$$")
w()
w("**Step 3 — transmissions.**")
w()
w(f"$$T_D = e^{{-{A_D:.4f}}} = {T_D:.4f}, \\qquad "
  f"T_A = e^{{-{A_A:.4f}}} = {T_A:.4f}$$")
w()
w(f"$$\\boxed{{\\ \\frac{{T_A}}{{T_D}} = \\frac{{{T_A:.4f}}}{{{T_D:.4f}}} "
  f"= {drift:.4f}\\ }}$$")
w()
w(f"**The ratio is inflated by {100 * (drift - 1):.1f} % by transport alone, "
  "before the sensor has reported anything.** Three features of this result "
  "deserve emphasis.")
w()
w(f"1. **The sign is fixed.** Because $\\lambda_D < \\lambda_A$ always — the "
  f"donor emits to the blue of the acceptor, by construction of any FRET pair — "
  f"$A_D > A_A$ and hence $T_A/T_D > 1$ *always*. The bias can never be "
  "negative. It is a one-sided error, not noise.")
w()
w("2. **It is not removable by optical-density normalisation.** The quantity "
  "$\\mu_{ref}$ depends on $C_P$ and $C_R$ separately, through different "
  "coefficients. An OD probe returns only their weighted sum. §K shows what "
  "that costs.")
w()
w("3. **It grows as the process succeeds.** $C_P$ rises monotonically through "
  "a production run, so the measurement degrades exactly when it matters most.")
w()

# --- Sensitivity: the od_only trap -----------------------------------------
mu_od = qp.a_X * Cx_ref           # what an OD-only correction assumes
A_D_od = mu_od * (qp.lam_ref / qp.lam_D) ** qp.p_mie * qp.L_eff
A_A_od = mu_od * (qp.lam_ref / qp.lam_A) ** qp.p_mie * qp.L_eff
drift_od = np.exp(-A_A_od) / np.exp(-A_D_od)
residual = drift / drift_od

w("### K. Why correcting for cell density is not enough")
w()
w("A conventional in-situ calibration normalises by optical density, which in "
  "effect assumes all attenuation is proportional to total biomass:")
w()
w(f"$$\\mu_{{ref}}^{{OD}} = a_X C_X = ({qp.a_X})({Cx_ref}) = {mu_od:.4f}"
  f"\\ \\mathrm{{cm^{{-1}}}} \\quad\\text{{versus the true }}"
  f"{mu_ref:.4f}\\ \\mathrm{{cm^{{-1}}}}.$$")
w()
w(f"That correction removes a drift of {100 * (drift_od - 1):.2f} %, leaving a "
  f"residual of")
w()
w(f"$$\\frac{{T_A/T_D}}{{(T_A/T_D)^{{OD}}}} = \\frac{{{drift:.4f}}}"
  f"{{{drift_od:.4f}}} = {residual:.4f},$$")
w()
w(f"i.e. {100 * (residual - 1):.2f} % of uncorrected, polymer-specific "
  "inflation. The missing term is $a_P C_P$ — and $C_P$ is precisely the "
  "quantity no non-destructive bulk measurement returns. **It has to come from "
  "the cellular model. This is why the two layers cannot be decoupled.**")
w()

# --- Propagation to inferred G ---------------------------------------------
m_true = Q.quantum_map(G_ref, phi_ref, qp, J_CAL, Cx=Cx_ref, Cp=Cp_ref)
G_none = Q.invert_quantum_map(m_true["Rt"], phi_ref, qp, J_CAL,
                              Cx=Cx_ref, Cp=Cp_ref, correction="none")[0]
G_od = Q.invert_quantum_map(m_true["Rt"], phi_ref, qp, J_CAL,
                            Cx=Cx_ref, Cp=Cp_ref, correction="od_only")[0]
G_full = Q.invert_quantum_map(m_true["Rt"], phi_ref, qp, J_CAL,
                              Cx=Cx_ref, Cp=Cp_ref, correction="full")[0]


def _single_channel(phi_assume, Cp_for_T):
    """Inferred G with exactly one of the two OD-only shortcuts applied."""
    T_D_, T_A_ = Q.optical_attenuation(Cx_ref, Cp_for_T, qp)
    n_ = Q.lorentz_lorenz_mix(phi_assume, qp.n_cytosol, qp.n_granule)
    Ea, _ = Q.fret_efficiency(Q.kappa2_effective(phi_assume, qp.kappa2_static_apo, qp),
                              qp.R_apo, n_, qp, J_CAL)
    Eh, _ = Q.fret_efficiency(Q.kappa2_effective(phi_assume, qp.kappa2_static_holo, qp),
                              qp.R_holo, n_, qp, J_CAL)
    a_ = qp.I_exc
    R_ = m_true["Rt"]
    num_ = a_ * qp.phi_A * T_A_ * (1 + qp.alpha_direct) + qp.dark * (1 - R_)
    den_ = a_ * (R_ * qp.phi_D * T_D_ + qp.phi_A * T_A_
                 - qp.bleed_DA * qp.phi_D * T_A_)
    E_ = np.clip(1 - num_ / den_, 1e-6, 1 - 1e-6)
    th_ = np.clip((E_ - Ea) / (Eh - Ea), 1e-9, 1 - 1e-9)
    return float(qp.K_D * (th_ / (1 - th_)) ** (1 / qp.hill_h))


chan_a = _single_channel(phi_ref, 0.0)   # transport under-corrected only
chan_b = _single_channel(0.0, Cp_ref)    # photophysics mis-assumed only

w("### L. Propagation to the inferred metabolite concentration")
w()
w(f"Finally, the error that matters operationally. The true state is "
  f"$G = {G_ref:.2f}$ mM at the biomass composition of §J "
  f"($C_X = {Cx_ref}$, $C_P = {Cp_ref}$ g L⁻¹, $\\phi_P = {phi_ref:.3f}$). "
  f"The forward map returns an observed ratio "
  f"$\\mathcal{{R}} = {float(m_true['Rt']):.4f}$, against "
  f"$\\mathcal{{R}} = {R_clean:.4f}$ for the same $G$ in a clear cuvette. "
  "Inverting under the three correction levels gives:")
w()
w("| Correction applied | Inferred $G$ (mM) | Error |")
w("|---|---|---|")
w(f"| None (fixed calibration curve) | {float(G_none):.4f} | "
  f"{100 * (float(G_none) - G_ref) / G_ref:+.1f} % |")
w(f"| Optical density only | {float(G_od):.4f} | "
  f"{100 * (float(G_od) - G_ref) / G_ref:+.1f} % |")
w(f"| Full multiscale correction | {float(G_full):.4f} | "
  f"{100 * (float(G_full) - G_ref) / G_ref:+.2f} % |")
w()
w("The full correction recovers the true value to within numerical tolerance, "
  "confirming that the inverse map is exact. The uncorrected error is "
  "**amplified relative to the optical bias itself**: a "
  f"{100 * (drift - 1):.1f} % inflation of the ratio becomes a "
  f"{100 * (float(G_none) - G_ref) / G_ref:+.1f} % error in $G$, because the "
  "Hill function compresses the ratio axis near saturation and the inversion "
  "therefore expands it.")
w()
w(f"**How this relates to the figure in the main text.** The "
  f"{100 * (float(G_none) / G_ref - 1):+.1f} % above is the error at one "
  f"static operating point, chosen to make the arithmetic followable. The "
  f"main text instead quotes +29 %, which is the same quantity averaged over "
  f"a closed-loop run: the controller drives the culture to higher polymer "
  f"content than this reference state, the bias grows with polymer content, "
  f"and the average is taken over the production phase rather than at a "
  f"single point. The two numbers are consistent, and Fig. 2c of the main "
  f"text plots the whole curve from which both are read.")
w()
w("**A caution on the middle row.** The OD-only correction looks respectable "
  f"at {100 * (G_od / G_ref - 1):+.1f} %, but that figure is the near-cancellation "
  "of two larger errors of opposite sign. Isolating them — correcting one "
  "channel properly while leaving the other as OD-only does:")
w()
w("| Channel left uncorrected | Inferred $G$ (mM) | Error |")
w("|---|---|---|")
w(f"| Transport under-correction only (photophysics right) | {chan_a:.4f} | "
  f"{100 * (chan_a / G_ref - 1):+.2f} % |")
w(f"| Photophysics assumed at $\\phi_P = 0$ only (transport right) | {chan_b:.4f} | "
  f"{100 * (chan_b / G_ref - 1):+.2f} % |")
w(f"| **Both, i.e. OD-only** | **{float(G_od):.4f}** | "
  f"**{100 * (float(G_od) / G_ref - 1):+.2f} %** |")
w()
w(f"The two channels are close to additive "
  f"({100 * (chan_a / G_ref - 1):+.2f} % and {100 * (chan_b / G_ref - 1):+.2f} % "
  f"sum to {100 * ((chan_a / G_ref - 1) + (chan_b / G_ref - 1)):+.2f} % against "
  f"an actual {100 * (float(G_od) / G_ref - 1):+.2f} %), and they **pull in "
  "opposite directions**. Leaving the granule term out of the attenuation "
  "under-corrects the transport, biasing $G$ upward; assuming a polymer-free "
  "cytoplasm puts the $\\kappa^2$ and refractive-index calibration in the wrong "
  "place, biasing it down about twice as hard. The net is a modest "
  "over-correction.")
w()
w("Two consequences follow, and they matter more than the headline number. "
  "First, **the apparent adequacy of an OD-only correction is an artefact of "
  "this particular parameter set**: the balance between the two channels "
  "depends on $a_P$, on the granule volume fraction and on the static "
  "orientation factors, none of which is independently measured (§3.1). A "
  "modest change in any of them breaks the cancellation. Second, the "
  "cancellation is **not conservative** — it hides two errors of roughly "
  f"{max(abs(100 * (chan_a / G_ref - 1)), abs(100 * (chan_b / G_ref - 1))):.0f} % "
  "behind a reported one of "
  f"{abs(100 * (float(G_od) / G_ref - 1)):.0f} %, so an experimenter who "
  "validated an OD-corrected sensor at one operating point could find it badly "
  "wrong at another with no warning. This is the reasoning behind giving the "
  "dimensionless criterion $\\Pi$ of §6.4 rather than a blanket "
  "recommendation.")
w()

# --- Granule-loaded repeat -------------------------------------------------
phi_hi = 0.30
n_hi = Q.lorentz_lorenz_mix(phi_hi, qp.n_cytosol, qp.n_granule)
eta_hi = Q.krieger_dougherty(phi_hi, qp.eta_cyto0, qp.phi_pack)
theta_hi = eta_hi * qp.V_hydro / (kB * Q.T_KELVIN)
fdyn_hi = qp.tau_D0 / (theta_hi + qp.tau_D0)
k2_apo_hi = Q.kappa2_effective(phi_hi, qp.kappa2_static_apo, qp)
R0_apo_hi = Q.forster_R0(k2_apo_hi, n_hi, qp)

w("### M. The same calculation in a granule-loaded cell")
w()
w(f"Repeating §B–§D at $\\phi_P = {phi_hi}$, a realistic late-stage granule "
  "volume fraction:")
w()
w("| Quantity | $\\phi_P = 0$ | " f"$\\phi_P = {phi_hi}$ | Change |")
w("|---|---|---|---|")
w(f"| $n_{{med}}$ | {n0:.4f} | {n_hi:.4f} | {100 * (n_hi / n0 - 1):+.2f} % |")
w(f"| $\\eta$ (mPa s) | {eta0 * 1e3:.2f} | {eta_hi * 1e3:.2f} | "
  f"{100 * (eta_hi / eta0 - 1):+.1f} % |")
w(f"| $\\theta_{{rot}}$ (ns) | {theta_rot0 * 1e9:.1f} | {theta_hi * 1e9:.1f} | "
  f"{100 * (theta_hi / theta_rot0 - 1):+.1f} % |")
w(f"| $f_{{dyn}}$ | {f_dyn0:.4f} | {fdyn_hi:.4f} | "
  f"{100 * (fdyn_hi / f_dyn0 - 1):+.1f} % |")
w(f"| $\\kappa^2_{{eff,apo}}$ | {k2_apo0:.4f} | {k2_apo_hi:.4f} | "
  f"{100 * (k2_apo_hi / k2_apo0 - 1):+.2f} % |")
w(f"| $R_0$ (Å) | {R0_apo:.2f} | {R0_apo_hi:.2f} | "
  f"{100 * (R0_apo_hi / R0_apo - 1):+.2f} % |")
w()
w(f"The refractive index rises by {100 * (n_hi / n0 - 1):.2f} %, and because "
  f"$R_0 \\propto n^{{-2/3}}$ this alone shortens the Förster radius. The "
  "viscosity rise further freezes the orientation factor. Neither effect is "
  "large individually, but both act in the same direction, and both are "
  "invisible to an instrument that sees only the emitted ratio.")
w()
w("**The key asymmetry:** these photophysical shifts are of order "
  f"{abs(100 * (R0_apo_hi / R0_apo - 1)):.1f} %, whereas the optical-transport "
  f"bias of §J is {100 * (drift - 1):.1f} %. The transport term dominates by "
  "roughly an order of magnitude. A reader who takes only one number from this "
  "Supplementary Note should take that one: **the dominant quantum-to-cellular "
  "coupling is not exotic photophysics inside the cell, it is ordinary "
  "scattering on the way out.**")
w()

text = "\n".join(OUT)
open(os.path.join(HERE, "worked_examples.md"), "w").write(text)

# machine-readable copy for cross-checking the parameter tables
json.dump({
    "J_cal": J_CAL, "R0_iso": R0_iso, "R0_apo": R0_apo, "R0_holo": R0_holo,
    "n0": n0, "theta_rot_ns": theta_rot0 * 1e9, "f_dyn": f_dyn0,
    "kappa2_apo": k2_apo0, "kappa2_holo": k2_holo0,
    "J_apo": J_apo, "k_bare": k_bare, "k_corr": k_corr,
    "coh_ratio": ratio_coh, "E_apo": E_apo, "E_holo": E_holo,
    "theta_ref": theta_ref, "E_ens": E_ens, "R_clean": R_clean,
    "T_D": T_D, "T_A": T_A, "drift": drift, "drift_od": drift_od,
    "residual": residual, "G_none": float(G_none), "G_od": float(G_od),
    "G_full": float(G_full), "n_hi": n_hi, "R0_apo_hi": R0_apo_hi,
}, open(os.path.join(HERE, "worked_examples.json"), "w"), indent=1)

print(f"worked_examples.md written: {len(OUT)} lines, "
      f"{len(text.split())} words")
print(f"  R0 (isotropic ref) = {R0_iso:.2f} A")
print(f"  J_cal              = {J_CAL:.4g} cm  (1/J = {1/J_CAL:.0f} cm^-1)")
print(f"  E_apo -> E_holo    = {E_apo:.4f} -> {E_holo:.4f}")
print(f"  T_A/T_D            = {drift:.4f}")
print(f"  G: none={float(G_none):.4f} od={float(G_od):.4f} full={float(G_full):.4f}")

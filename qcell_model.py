"""
QCELL: a quantum-to-cellular multiscale framework for closed-loop PHA biosynthesis.

Couples
  (Q) a femtosecond-to-nanosecond excitonic / FRET photophysics layer describing a
      genetically encoded PII-based 2-oxoglutarate biosensor, with
  (C) a minute-to-hour cellular layer: a bioelectrochemical formate supply
      (Nernst-Monod, after Shemfe et al. 2018) feeding a dynamic physics-based
      PHA fermentation model (after Sadhukhan, Sen & Hu 2026).

The coupling is bidirectional:
  C -> Q : intracellular 2-OG G, PHA granule volume fraction phi_P (-> refractive
           index n_cyt and cytoplasmic viscosity eta -> kappa^2), sensor loading.
  Q -> C : the ratiometric observable Rt -> UKF state estimate -> nitrogen feed
           F_N and applied BES potential V_app; plus expression burden on mu_max.

Author: generated reference implementation.
Units: g/L, h, cm^-1, Angstrom, mM unless stated.
"""

import numpy as np
from dataclasses import dataclass, field, asdict

# ----------------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------------
C_LIGHT = 2.99792458e10      # cm/s
N_AVO = 6.02214076e23
F_FARADAY = 96485.0          # C/mol
R_GAS = 8.314                # J/mol/K
T_KELVIN = 298.15

# Debye^2 / Angstrom^3 -> cm^-1 for the point-dipole coupling
DIPOLE_PREFACTOR = 5034.1    # = 1e8 * 1e-18^2 / (4 pi eps0) / (h c) ... see docs


# ============================================================================
# QUANTUM LAYER
# ============================================================================

@dataclass
class QuantumParams:
    """Photophysics of the donor-acceptor FRET pair on the PII sensor scaffold."""
    # Donor / acceptor photophysics (mTurquoise2 -> mVenus like)
    phi_D: float = 0.93          # donor fluorescence quantum yield
    phi_A: float = 0.57          # acceptor fluorescence quantum yield
    tau_D0: float = 4.0e-9       # donor excited-state lifetime, s (no acceptor)
    J_overlap: float = 2.4e15    # spectral overlap, M^-1 cm^-1 nm^4
    mu_D: float = 4.37           # donor transition dipole, Debye
    mu_A: float = 8.86           # acceptor transition dipole, Debye
    nu_D: float = 23041.0        # donor emission centroid, cm^-1
    nu_A: float = 19417.0        # acceptor absorption centroid, cm^-1

    # Sensor geometry in the two conformational states
    R_apo: float = 68.0          # donor-acceptor separation, apo state, Angstrom
    R_holo: float = 46.0         # donor-acceptor separation, 2-OG bound, Angstrom
    kappa2_static_apo: float = 0.52
    kappa2_static_holo: float = 1.15

    # Ligand binding (PII trimer, 2-OG; negative cooperativity -> h slightly > 1)
    K_D: float = 0.20            # 2-OG dissociation constant, mM (PII sensor variant)
    hill_h: float = 1.15

    # Hydrodynamics of the fluorophore domains (for kappa^2 dynamic averaging)
    V_hydro: float = 8.0e-26     # m^3, effective rotational volume of an FP barrel
    eta_cyto0: float = 3.0e-3    # Pa s, cytoplasmic viscosity at phi_P = 0

    # Optical environment
    n_cytosol: float = 1.360     # refractive index of PHA-free cytosol
    n_granule: float = 1.490     # refractive index of a PHB granule
    # Random close packing of MONODISPERSE spheres is 0.64 (maximally random
    # jammed), not 0.68; see Torquato, Truskett & Debenedetti, Phys. Rev.
    # Lett. 84, 2064 (2000). 0.68 lies between random close packing and
    # crystalline close packing, where no random monodisperse packing sits.
    phi_pack: float = 0.64       # granule random-close-packing limit

    # Detection channel
    alpha_direct: float = 0.085  # direct acceptor excitation / donor excitation
    bleed_DA: float = 0.055      # donor bleed-through into acceptor channel
    I_exc: float = 1.0           # relative excitation irradiance
    dark: float = 0.004          # relative dark/background offset

    # Dephasing (links to the excitonic-coupling literature)
    T2_star: float = 60e-15      # pure dephasing time, s

    # Optical transport from the cell to an in-situ probe. At process scale the
    # ratio is read through a turbid suspension, so the donor and acceptor bands
    # are attenuated unequally by cell and granule scattering (Mie regime).
    lam_D: float = 480.0         # nm, donor emission band
    lam_A: float = 530.0         # nm, acceptor emission band
    lam_ref: float = 500.0       # nm, reference wavelength
    p_mie: float = 1.5           # scattering wavelength exponent
    a_X: float = 0.250           # L/g/cm, attenuation per unit residual biomass
    a_P: float = 0.120           # L/g/cm, EXTRA attenuation per unit PHA granule
    L_eff: float = 0.30          # cm, effective optical path of the probe


def lorentz_lorenz_mix(phi, n_host, n_incl):
    """Effective refractive index of granules dispersed in cytosol.

    Volume-additive polarisability (Lorentz-Lorenz / Clausius-Mossotti):
        (n^2-1)/(n^2+2) = (1-phi)(n_h^2-1)/(n_h^2+2) + phi (n_i^2-1)/(n_i^2+2)
    """
    f_h = (n_host**2 - 1.0) / (n_host**2 + 2.0)
    f_i = (n_incl**2 - 1.0) / (n_incl**2 + 2.0)
    f = (1.0 - phi) * f_h + phi * f_i
    return np.sqrt((1.0 + 2.0 * f) / (1.0 - f))


def krieger_dougherty(phi, eta0, phi_max=0.68, exponent=2.5):
    """Suspension viscosity as granules fill the cytoplasm."""
    x = np.clip(phi / phi_max, 0.0, 0.95)
    return eta0 * (1.0 - x) ** (-exponent * phi_max)


def kappa2_effective(phi_P, kappa2_static, qp: QuantumParams, T=T_KELVIN):
    """Orientation factor under partial rotational averaging.

    Steady-state anisotropy gives the FROZEN fraction, r_ss/r_0 = 1/(1+tau_D/theta);
    the complementary part has reoriented within the donor lifetime, so
        f_dyn = tau_D / (theta_rot + tau_D).
    Full averaging gives kappa^2 -> 2/3; a frozen geometry retains kappa2_static.
    Fluorescent-protein barrels tumble on tens of nanoseconds against a ~4 ns
    donor lifetime, so kappa^2 stays close to its static value throughout.
    """
    eta = krieger_dougherty(phi_P, qp.eta_cyto0, qp.phi_pack)
    theta_rot = eta * qp.V_hydro / (1.380649e-23 * T)      # s
    f_dyn = qp.tau_D0 / (theta_rot + qp.tau_D0)
    return (2.0 / 3.0) * f_dyn + kappa2_static * (1.0 - f_dyn)


def forster_R0(kappa2, n_med, qp: QuantumParams):
    """Forster radius in Angstrom. J_overlap in M^-1 cm^-1 nm^4."""
    return 0.211 * (kappa2 * n_med**-4 * qp.phi_D * qp.J_overlap) ** (1.0 / 6.0)


def coulomb_coupling(kappa, R_ang, n_med, qp: QuantumParams):
    """Point-dipole electronic coupling J in cm^-1, optically screened by n^2."""
    return DIPOLE_PREFACTOR * kappa * qp.mu_D * qp.mu_A / (n_med**2 * R_ang**3)


def spectral_overlap_density(qp: QuantumParams):
    """Effective density of accepting states (cm) calibrated so that the golden
    rule k = 4 pi^2 c J^2 Jcal reproduces the Forster rate at R = R0.

    This makes the quantum layer internally consistent: the same coupling J that
    appears in excitonic treatments reproduces the measured Forster radius.
    """
    kappa2_ref = 2.0 / 3.0
    n_ref = 1.33
    R0 = forster_R0(kappa2_ref, n_ref, qp)
    J0 = coulomb_coupling(np.sqrt(kappa2_ref), R0, n_ref, qp)
    return 1.0 / (4.0 * np.pi**2 * C_LIGHT * J0**2 * qp.tau_D0)


def eet_rate(kappa2, R_ang, n_med, qp: QuantumParams, Jcal, coherent=True):
    """Excitation energy transfer rate, s^-1.

    Golden-rule form k = 4 pi^2 c J^2 Jcal, which equals the Forster rate.
    With coherent=True a Haken-Strobl / Redfield correction is applied that
    attenuates the rate when the dephasing rate gamma is not large compared with
    the coupling: k_corr = k / (1 + (2 omega_J / gamma)^2).
    """
    kappa = np.sqrt(np.maximum(kappa2, 1e-12))
    J = coulomb_coupling(kappa, R_ang, n_med, qp)
    k = 4.0 * np.pi**2 * C_LIGHT * J**2 * Jcal
    if coherent:
        omega_J = 2.0 * np.pi * C_LIGHT * J          # rad/s
        gamma = 1.0 / qp.T2_star
        k = k / (1.0 + (2.0 * omega_J / gamma) ** 2)
    return k, J


def fret_efficiency(kappa2, R_ang, n_med, qp: QuantumParams, Jcal, coherent=True):
    k, J = eet_rate(kappa2, R_ang, n_med, qp, Jcal, coherent)
    E = k / (k + 1.0 / qp.tau_D0)
    return E, J


def optical_attenuation(Cx, Cp, qp: QuantumParams):
    """Band-resolved attenuation of the two emission channels, exp(-A).

    Granules raise the scattering coefficient beyond what biomass alone explains,
    and because scattering falls with wavelength the shorter donor band is
    attenuated more. The ratio is therefore inflated as polymer accumulates --
    a drift that normalising by optical density alone cannot remove, because it
    depends on the polymer content, not the cell count.
    """
    Cr = np.maximum(Cx - Cp, 0.0)
    mu_ref = qp.a_X * Cr + (qp.a_X + qp.a_P) * np.maximum(Cp, 0.0)
    A_D = mu_ref * (qp.lam_ref / qp.lam_D) ** qp.p_mie * qp.L_eff
    A_A = mu_ref * (qp.lam_ref / qp.lam_A) ** qp.p_mie * qp.L_eff
    return np.exp(-A_D), np.exp(-A_A)


def sensor_occupancy(G_mM, qp: QuantumParams):
    """Fractional 2-OG occupancy of the PII sensor (Hill)."""
    G = np.maximum(G_mM, 0.0)
    return G**qp.hill_h / (G**qp.hill_h + qp.K_D**qp.hill_h)


def quantum_map(G_mM, phi_P, qp: QuantumParams, Jcal, coherent=True,
                Cx=0.0, Cp=0.0):
    """The full C -> Q map.

    Returns the ratiometric observable and the intermediate quantities.
    Conformational exchange (microseconds) is slow relative to the donor
    lifetime (nanoseconds), so the ensemble averages efficiencies, not geometries.
    """
    n_med = lorentz_lorenz_mix(phi_P, qp.n_cytosol, qp.n_granule)
    k2_apo = kappa2_effective(phi_P, qp.kappa2_static_apo, qp)
    k2_holo = kappa2_effective(phi_P, qp.kappa2_static_holo, qp)

    E_apo, J_apo = fret_efficiency(k2_apo, qp.R_apo, n_med, qp, Jcal, coherent)
    E_holo, J_holo = fret_efficiency(k2_holo, qp.R_holo, n_med, qp, Jcal, coherent)

    theta = sensor_occupancy(G_mM, qp)
    E = (1.0 - theta) * E_apo + theta * E_holo

    # Ratiometric readout: acceptor channel over donor channel
    I_D = qp.I_exc * qp.phi_D * (1.0 - E)
    I_A = qp.I_exc * (qp.phi_A * E + qp.alpha_direct * qp.phi_A
                      + qp.bleed_DA * qp.phi_D * (1.0 - E))
    # Optical transport to the probe
    T_D, T_A = optical_attenuation(Cx, Cp, qp)
    I_D = I_D * T_D + qp.dark
    I_A = I_A * T_A + qp.dark
    Rt = I_A / I_D

    return dict(Rt=Rt, E=E, E_apo=E_apo, E_holo=E_holo, theta=theta,
                T_D=T_D, T_A=T_A, ratio_drift=T_A / T_D,
                n_med=n_med, kappa2_apo=k2_apo, kappa2_holo=k2_holo,
                R0_apo=forster_R0(k2_apo, n_med, qp),
                R0_holo=forster_R0(k2_holo, n_med, qp),
                J_apo=J_apo, J_holo=J_holo, I_D=I_D, I_A=I_A)


def invert_quantum_map(Rt, phi_P, qp: QuantumParams, Jcal, coherent=True,
                       phi_assumed=None, Cx=0.0, Cp=0.0, correction="full"):
    """The Q -> C inverse: recover 2-OG from the ratio.

    correction = 'full'    : the multiscale correction. Uses the model's own
                             estimates of granule loading and polymer content to
                             undo both the photophysical drift (n, kappa^2) and
                             the band-resolved optical transport.
    correction = 'od_only' : what a conventional calibration does. Normalises for
                             cell density, but assumes a polymer-free cytoplasm,
                             so it leaves the granule-specific terms uncorrected.
    correction = 'none'    : raw ratio against a fixed calibration curve.
    """
    if correction == "full":
        phi_use = phi_P if phi_assumed is None else phi_assumed
        T_D, T_A = optical_attenuation(Cx, Cp, qp)
    elif correction == "od_only":
        phi_use = 0.0
        T_D, T_A = optical_attenuation(Cx, 0.0, qp)   # polymer ignored
    else:
        phi_use = 0.0
        T_D = T_A = 1.0
    n_med = lorentz_lorenz_mix(phi_use, qp.n_cytosol, qp.n_granule)
    k2_apo = kappa2_effective(phi_use, qp.kappa2_static_apo, qp)
    k2_holo = kappa2_effective(phi_use, qp.kappa2_static_holo, qp)
    E_apo, _ = fret_efficiency(k2_apo, qp.R_apo, n_med, qp, Jcal, coherent)
    E_holo, _ = fret_efficiency(k2_holo, qp.R_holo, n_med, qp, Jcal, coherent)

    # Invert the detection model for E, then the mixing law for theta.
    # Forward:  I_D = a.phi_D.(1-E).T_D               + d
    #           I_A = a.[phi_A.E + alpha.phi_A + b.phi_D.(1-E)].T_A + d
    # Writing u = 1-E, the acceptor channel is a.T_A.[phi_A(1+alpha)
    # + u(b.phi_D - phi_A)] + d, so Rt.I_D = I_A is LINEAR in u:
    #   u[Rt.phi_D.T_D + phi_A.T_A - b.phi_D.T_A] = phi_A.T_A(1+alpha)
    #                                               + (d/a)(1 - Rt).
    # The transmissions must enter here rather than as a pre-scaling
    # Rt -> Rt.(T_D/T_A): that shortcut is exact only for a dark offset of
    # zero, because d is added AFTER attenuation and so does not scale with
    # it. With d = 0.004 the shortcut biases recovered 2-OG by -0.5 % at low
    # turbidity and by -10 % at Cx = 16 g/L, i.e. worst exactly where the
    # correction is most needed.
    a = qp.I_exc
    num = a * qp.phi_A * T_A * (1.0 + qp.alpha_direct) + qp.dark * (1.0 - Rt)
    den = a * (Rt * qp.phi_D * T_D + qp.phi_A * T_A - qp.bleed_DA * qp.phi_D * T_A)
    u = num / den
    E = np.clip(1.0 - u, 1e-6, 1 - 1e-6)

    theta = np.clip((E - E_apo) / (E_holo - E_apo), 1e-9, 1 - 1e-9)
    G = qp.K_D * (theta / (1.0 - theta)) ** (1.0 / qp.hill_h)
    return G, E, theta


# ============================================================================
# CELLULAR LAYER
# ============================================================================

@dataclass
class BESParams:
    """Lumped bioelectrochemical formate generator (after Shemfe et al. 2018)."""
    q_max: float = 7.92 / 24.0     # kgCOD/kgVS/h  (paper reports per day)
    K_SD: float = 0.174            # kgCOD/m3 == g/L
    E_KA: float = 0.0              # V, half-max rate potential (vs ref)
    V1: float = 0.228              # applied-potential linear coefficient
    V2: float = 0.0
    A_cat: float = 25000.0         # cm2 (2.5 m2 modular cathode stack)
    V_cat: float = 290.0           # cm3 per module
    FE: float = 0.949              # Faradaic efficiency to formate
    M_w: float = 46.03             # g/mol formate
    n_e: float = 2.0
    j_ref: float = 2.6             # A/m2 reference current density at eta_ref
    eta_ref: float = 0.30          # V
    S_cod_anode: float = 1.0       # g/L wastewater COD held at the anode
    M_res0: float = 60.0           # g formate initially in the buffer reservoir
    M_res_max: float = 400.0       # g, reservoir capacity


def bes_formate_rate(V_app, S_cod_anode, bp: BESParams, T=T_KELVIN):
    """Formate production rate from the BES cathode, g/h.

    Nernst-Monod oxidation of wastewater COD at the anode sets the current
    density (Shemfe Eqs. 1, 7-8); the current and the cathodic Faradaic
    efficiency set formate synthesis (Shemfe Eq. 15). S_cod_anode is the
    wastewater COD in the ANODE chamber, which is independent of the formate
    concentration in the downstream fermenter.
    """
    eta = bp.V1 * V_app - bp.V2 - bp.E_KA
    nernst = 1.0 / (1.0 + np.exp(-F_FARADAY * eta / (R_GAS * T)))
    monod = S_cod_anode / (bp.K_SD + S_cod_anode)
    j = bp.j_ref * (eta / bp.eta_ref) * nernst * monod          # A/m2
    j = np.maximum(j, 0.0)
    A_m2 = bp.A_cat * 1e-4
    I = j * A_m2                                                # A
    mol_per_s = I * bp.FE / (bp.n_e * F_FARADAY)
    return mol_per_s * bp.M_w * 3600.0                          # g/h


@dataclass
class CellParams:
    """Dynamic physics-based PHA model (structure after Sadhukhan, Sen & Hu 2026),
    parameterised for formate-grown Cupriavidus necator."""
    mu_xs_max: float = 0.18        # 1/h  (Gorret et al. 2015, formate)
    mu_xp_max: float = 0.040       # 1/h, growth on remobilised PHA
    mu_ps_max: float = 0.200       # g PHA / g cell / h

    Y_xs: float = 0.090            # g cell / g formate
    Y_pi: float = 0.100            # g PHA / g formate
    Y_xp: float = 0.60             # g cell / g PHA (remobilisation)
    Y_xN: float = 8.0              # g cell / g N
    m_s: float = 0.020             # g formate / g cell / h

    K_S: float = 0.10              # g/L
    K_IS: float = 2.00             # g/L, formate substrate inhibition
    K_N: float = 0.050             # g/L
    K_IN: float = 10.0             # g/L
    K_P: float = 0.10              # PHA fraction saturation for mu_xp
    K_PS: float = 0.20             # g/L
    K_PIS: float = 2.00            # g/L
    K_PIN: float = 0.020           # g/L, nitrogen conc. half-inhibiting PHA synth.

    # Active / inactive biomass split: total nitrogen starvation inactivates the
    # PHA synthase machinery, so the optimum is nitrogen-LIMITED, not nitrogen-free
    b_dec0: float = 0.0040         # 1/h, baseline inactivation
    b_decN: float = 0.0900         # 1/h, additional inactivation under N starvation
    K_Nv: float = 0.0100           # g/L, nitrogen conc. half-protecting viability

    Cx_max: float = 15.0           # g/L
    alpha: float = 1.0
    fP_max: float = 0.80
    beta: float = 2.0

    # 2-oxoglutarate node
    k_G: float = 1.35              # mM per (g formate / g cell / h)
    V_GS: float = 2.60             # mM/h, GS-GOGAT capacity
    K_G: float = 0.50              # mM
    K_NG: float = 0.020            # g/L
    G_basal: float = 0.030         # mM/h leak
    k_Gdeg: float = 0.08           # 1/h, 2-OG efflux / anaplerotic drain

    # Sensor expression and burden
    k_exp: float = 0.0115          # g sensor / g cell / h
    k_bleach: float = 0.004        # 1/h at unit irradiance
    burden: float = 0.55           # fractional mu loss per unit specific sensor

    # Feed composition
    C_SF: float = 600.0            # g/L formate (concentrated formic acid, also the pH titrant)
    C_NF: float = 40.0             # g/L nitrogen in concentrated nitrogen feed

    rho_P: float = 1.25            # g/cm3, PHB granule
    rho_R: float = 1.10            # g/cm3, residual biomass
    rho_cell_wet: float = 1.00     # g/cm3 approximate wet cell density
    dry_frac: float = 0.30         # dry weight fraction of wet cell volume


def granule_volume_fraction(Cx, Cp, cp: CellParams):
    """PHA granule volume fraction within the cytoplasm."""
    Cp = np.clip(Cp, 0.0, None)
    Cr = np.maximum(Cx - Cp, 1e-9)
    # Granules are essentially anhydrous; only the residual biomass is hydrated.
    v_P = Cp / cp.rho_P
    v_R_wet = (Cr / cp.dry_frac) / cp.rho_cell_wet
    return np.clip(v_P / (v_P + v_R_wet), 0.0, 0.95)


def inactivation_rate(Cn, cp: CellParams):
    """Active-biomass inactivation; rises sharply when nitrogen is fully absent."""
    return cp.b_dec0 + cp.b_decN * cp.K_Nv / (Cn + cp.K_Nv)


def specific_rates(Cs, Cn, Cx, Cp, cp: CellParams, burden_frac=0.0, phi_a=1.0):
    """mu_xs, mu_xp, mu_ps after Eqs (6)-(8) of the DPM."""
    Cx = np.maximum(Cx, 1e-9)
    fP = np.clip(Cp / Cx, 0.0, cp.fP_max * 0.999)

    gate_S = Cs / (cp.K_S + Cs + Cs**2 / cp.K_IS)
    gate_N = Cn / (cp.K_N + Cn + Cn**2 / cp.K_IN)
    gate_X = max(0.0, 1.0 - (Cx / cp.Cx_max) ** cp.alpha)

    scale = max(0.0, 1.0 - burden_frac)

    mu_xs = cp.mu_xs_max * scale * phi_a * gate_S * gate_N * gate_X
    mu_xp = cp.mu_xp_max * scale * (fP / (cp.K_P + fP)) * gate_N * gate_X

    gate_Nstarve = cp.K_PIN / (Cn + cp.K_PIN)
    gate_PS = Cs / (cp.K_PS + Cs + Cs**2 / cp.K_PIS)
    gate_sat = max(0.0, 1.0 - (fP / cp.fP_max) ** cp.beta)
    mu_ps = cp.mu_ps_max * scale * phi_a * gate_Nstarve * gate_PS * gate_sat

    return mu_xs, mu_xp, mu_ps, fP


def rhs(t, y, cp: CellParams, bp: BESParams, qp: QuantumParams, Jcal,
        control, log=None):
    """Right-hand side of the coupled cellular system.

    y = [Cs, Cn, Cx, Cp, G, Csens, V, M_res, phi_a]
    """
    Cs, Cn, Cx, Cp, G, Csens, V, M_res, phi_a = y
    phi_a = float(np.clip(phi_a, 0.0, 1.0))
    Cs = max(Cs, 0.0); Cn = max(Cn, 0.0); Cx = max(Cx, 1e-9)
    Cp = max(Cp, 0.0); G = max(G, 0.0); Csens = max(Csens, 0.0)
    V = max(V, 1e-6); M_res = max(M_res, 0.0)

    c_sens_specific = Csens / Cx
    burden_frac = cp.burden * c_sens_specific

    mu_xs, mu_xp, mu_ps, fP = specific_rates(Cs, Cn, Cx, Cp, cp, burden_frac, phi_a)
    mu_x = mu_xs + mu_xp

    # Control action: carbon feed, nitrogen feed, BES applied potential
    F_S, F_N, V_app = control(t, y)

    # The fermenter cannot draw carbon the BES has not yet made
    if M_res <= 0.0:
        F_S = min(F_S, bes_formate_rate(V_app, bp.S_cod_anode, bp) / cp.C_SF)
    F_S = max(F_S, 0.0)

    R_form = bes_formate_rate(V_app, bp.S_cod_anode, bp)      # g/h
    dM_res = R_form - F_S * cp.C_SF
    if M_res >= bp.M_res_max and dM_res > 0:
        dM_res = 0.0                                          # reservoir full, vent/idle

    F_tot = F_S + F_N
    D = F_tot / V

    # Substrate (Eq. 1): feed - dilution - growth/product/maintenance demand
    q_S = (mu_xs / cp.Y_xs + mu_ps / cp.Y_pi + cp.m_s)
    dCs = F_S * cp.C_SF / V - D * Cs - q_S * Cx

    # Nitrogen (Eq. 2)
    q_N = (mu_xs + mu_xp) / cp.Y_xN
    dCn = F_N * cp.C_NF / V - D * Cn - q_N * Cx

    # Biomass (Eq. 3) and PHA (Eq. 4)
    dCx = (mu_xs + mu_xp - D) * Cx
    dCp = (mu_ps - mu_xp / cp.Y_xp) * Cx - D * Cp

    # 2-oxoglutarate node: produced with carbon flux, drained by GS-GOGAT
    v_prod = cp.k_G * q_S + cp.G_basal
    v_gogat = cp.V_GS * (G / (cp.K_G + G)) * (Cn / (cp.K_NG + Cn))
    dG = v_prod - v_gogat - (cp.k_Gdeg + mu_x) * G

    # Sensor: constitutive expression, growth dilution, photobleaching
    dCsens = cp.k_exp * Cx - (mu_x + cp.k_bleach * qp.I_exc) * Csens

    dV = F_tot

    # Active fraction: lost to starvation-driven inactivation, restored by the
    # dilution of inactive material through new growth on nitrogen
    b_dec = inactivation_rate(Cn, cp)
    dphi_a = -b_dec * phi_a + mu_xs * (1.0 - phi_a)

    return np.array([dCs, dCn, dCx, dCp, dG, dCsens, dV, dM_res, dphi_a])


# ============================================================================
# SOLVER: co-simulation with adiabatic elimination of the quantum layer
# ============================================================================

def rk4_step(f, t, y, h, *args, **kw):
    k1 = f(t, y, *args, **kw)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1, *args, **kw)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2, *args, **kw)
    k4 = f(t + h, y + h * k3, *args, **kw)
    return y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


@dataclass
class Controller:
    """Closed-loop nitrogen feeding driven by the estimated 2-OG signal."""
    mode: str = "closed_multiscale"  # open | closed_raw | closed_od | closed_multiscale
    G_set: float = 0.122           # mM target 2-OG during accumulation
    Cn_set: float = 0.400          # g/L nitrogen setpoint during the growth phase
    Kc_N: float = 0.0060           # (L/h) per (g/L) nitrogen error
    tau_IN: float = 4.0            # h
    t_switch: float = 24.0         # h, open-loop switch time
    F_N_max: float = 0.00300       # L/h
    Kc: float = 0.00090            # proportional gain, L/h per mM
    tau_I: float = 3.0             # h, integral time
    V_app: float = 0.90            # V, nominal applied potential
    V_app_min: float = 0.70
    V_app_max: float = 1.05
    Cs_set: float = 0.45           # g/L formate setpoint (maximises the Haldane gate)
    F_S_max: float = 0.060         # L/h maximum carbon feed
    Kc_S: float = 0.030            # (L/h) per (g/L) formate error
    tau_IS: float = 1.0            # h
    M_res_set: float = 60.0        # g, formate buffer level setpoint
    Kc_V: float = 0.004            # V per g of reservoir error
    tau_IV: float = 10.0           # h
    Cx_target: float = 7.0         # g/L, biomass handover target
    noise_sd: float = 0.02         # fractional measurement noise on the ratio
    _integral: float = field(default=0.0, repr=False)
    _last_t: float = field(default=0.0, repr=False)
    _Ghat: float = field(default=0.0, repr=False)
    _phase: int = field(default=0, repr=False)
    _integral_S: float = field(default=0.0, repr=False)
    _integral_V: float = field(default=0.0, repr=False)
    _integral_N: float = field(default=0.0, repr=False)


def make_control(ctrl: Controller, cp: CellParams, qp: QuantumParams, Jcal, rng):
    """Build the control callback. This is where Q -> C closure happens."""

    def carbon_feed(Cs, dt):
        """Substrate-stat: hold residual formate below the inhibition knee."""
        err = ctrl.Cs_set - Cs
        ctrl._integral_S += err * dt
        ctrl._integral_S = np.clip(ctrl._integral_S, -0.5, 0.5)
        u = ctrl.Kc_S * (err + ctrl._integral_S / ctrl.tau_IS)
        return float(np.clip(u, 0.0, ctrl.F_S_max))

    def bes_potential(M_res, dt):
        """Hold the formate buffer near its level setpoint."""
        err = ctrl.M_res_set - M_res
        ctrl._integral_V += err * dt
        ctrl._integral_V = np.clip(ctrl._integral_V, -300.0, 300.0)
        v = ctrl.V_app + ctrl.Kc_V * (err + ctrl._integral_V / ctrl.tau_IV)
        return float(np.clip(v, ctrl.V_app_min, ctrl.V_app_max))

    def control(t, y):
        Cs, Cn, Cx, Cp, G, Csens, V, M_res, phi_a = y
        dt = max(t - ctrl._last_t, 0.0)
        ctrl._last_t = t

        F_S = carbon_feed(Cs, dt)
        V_app = bes_potential(M_res, dt)

        if ctrl.mode == "open":
            # Fair open-loop baseline: the growth phase is regulated exactly as in
            # the closed-loop cases (an N-stat, which is normal practice), and the
            # ONLY difference is that nitrogen is cut on a clock rather than on the
            # measured 2-oxoglutarate. This isolates the single decision compared.
            if t >= ctrl.t_switch:
                return F_S, 0.0, V_app
            errN = ctrl.Cn_set - Cn
            ctrl._integral_N += errN * dt
            ctrl._integral_N = np.clip(ctrl._integral_N, -5.0, 5.0)
            uN = ctrl.Kc_N * (errN + ctrl._integral_N / ctrl.tau_IN)
            return F_S, float(np.clip(uN, 0.0, ctrl.F_N_max)), V_app

        # --- measurement: the cells report through the quantum layer ---
        phi_P = granule_volume_fraction(Cx, Cp, cp)
        qm = quantum_map(max(G, 0.0), phi_P, qp, Jcal, Cx=Cx, Cp=Cp)
        Rt = qm["Rt"]
        if ctrl.noise_sd > 0:
            Rt = Rt * (1.0 + rng.normal(0.0, ctrl.noise_sd))

        # --- inversion at the chosen level of multiscale correction ---
        level = {"closed_raw": "none", "closed_od": "od_only",
                 "closed_multiscale": "full"}[ctrl.mode]
        Ghat, _, _ = invert_quantum_map(Rt, phi_P, qp, Jcal, Cx=Cx, Cp=Cp,
                                        correction=level)
        ctrl._Ghat = Ghat

        # --- phase logic + PI on the 2-OG setpoint ---
        if ctrl._phase == 0 and Cx >= ctrl.Cx_target:
            ctrl._phase = 1

        if ctrl._phase == 0:
            # Growth phase: N-stat holds nitrogen at a low, non-limiting level so
            # that control authority is available as soon as the phase changes.
            errN = ctrl.Cn_set - Cn
            ctrl._integral_N += errN * dt
            ctrl._integral_N = np.clip(ctrl._integral_N, -5.0, 5.0)
            uN = ctrl.Kc_N * (errN + ctrl._integral_N / ctrl.tau_IN)
            return F_S, float(np.clip(uN, 0.0, ctrl.F_N_max)), V_app

        err = ctrl.G_set - Ghat
        ctrl._integral += err * dt
        ctrl._integral = np.clip(ctrl._integral, -40.0, 40.0)
        u = -ctrl.Kc * (err + ctrl._integral / ctrl.tau_I)
        F_N = float(np.clip(u, 0.0, ctrl.F_N_max))
        return F_S, F_N, V_app

    return control


def simulate(cp: CellParams, bp: BESParams, qp: QuantumParams, ctrl: Controller,
             t_end=60.0, dt=0.01, y0=None, seed=0):
    Jcal = spectral_overlap_density(qp)
    rng = np.random.default_rng(seed)
    control = make_control(ctrl, cp, qp, Jcal, rng)

    if y0 is None:
        y0 = np.array([1.00, 0.60, 0.50, 0.005, 0.30, 0.0050, 1.0, bp.M_res0, 1.0])

    n = int(round(t_end / dt))
    ts = np.zeros(n + 1)
    ys = np.zeros((n + 1, len(y0)))
    aux = np.zeros((n + 1, 8))
    ts[0], ys[0] = 0.0, y0

    y = y0.copy()
    for i in range(n):
        t = i * dt
        y = rk4_step(rhs, t, y, dt, cp, bp, qp, Jcal, control)
        y = np.maximum(y, np.array([0, 0, 1e-9, 0, 0, 0, 1e-6, 0, 0]))
        y[8] = min(y[8], 1.0)
        ts[i + 1], ys[i + 1] = t + dt, y

        Cs, Cn, Cx, Cp, G, Csens, V, M_res, phi_a = y
        phi_P = granule_volume_fraction(Cx, Cp, cp)
        qm = quantum_map(G, phi_P, qp, Jcal)
        F_S, F_N, V_app = control(t + dt, y)
        aux[i + 1] = [phi_P, qm["Rt"], qm["E"], qm["n_med"],
                      V_app, ctrl._Ghat if ctrl.mode != 'open' else np.nan,
                      F_N, Cp / max(Cx, 1e-9)]

    aux[0] = aux[1]
    return dict(t=ts, y=ys, aux=aux, Jcal=Jcal,
                labels=["Cs", "Cn", "Cx", "Cp", "G", "Csens", "V", "M_res", "phi_a"],
                aux_labels=["phi_P", "Rt", "E", "n_med", "V_app", "Ghat",
                            "F_N", "fP"])


# ============================================================================
# DIAGNOSTICS
# ============================================================================

def metrics(res, cp: CellParams):
    t, y = res["t"], res["y"]
    Cs, Cn, Cx, Cp, G, Csens, V = [y[:, i] for i in range(7)]
    mass_P = Cp * V
    titre = Cp[-1]
    prod = mass_P[-1] / (V[-1] * t[-1])
    fP = Cp[-1] / max(Cx[-1], 1e-9)
    # carbon fed via the BES over the run
    phi_a = y[:, 8]
    viable = phi_a[-1] * Cx[-1]
    # Economic objective for a cell-recycle process: harvested polymer plus the
    # value of viable catalyst carried into the next cycle (w in g PHA-equivalent
    # per g viable residual biomass, i.e. the growth phase it saves).
    w_recycle = 0.35
    J = titre + w_recycle * viable
    return dict(titre_PHA=titre, fP_final=fP, Cx_final=Cx[-1],
                productivity=prod, mass_PHA=mass_P[-1], phi_a_final=phi_a[-1],
                viable_biomass=viable, objective_J=J,
                G_final=G[-1], residual_formate=Cs[-1])


def sensor_bias_scan(qp: QuantumParams, G_grid=None, phi_grid=None):
    """Systematic error in inferred 2-OG when granule loading is ignored."""
    Jcal = spectral_overlap_density(qp)
    if G_grid is None:
        G_grid = np.logspace(-2, 1.0, 120)
    if phi_grid is None:
        phi_grid = np.array([0.0, 0.10, 0.20, 0.30, 0.40, 0.50])
    out = np.zeros((len(phi_grid), len(G_grid)))
    for i, phi in enumerate(phi_grid):
        qm = quantum_map(G_grid, phi, qp, Jcal)
        Ghat, _, _ = invert_quantum_map(qm["Rt"], phi, qp, Jcal, phi_assumed=0.0)
        out[i] = (Ghat - G_grid) / G_grid * 100.0
    return G_grid, phi_grid, out


def sensor_design_map(qp: QuantumParams, KD_grid=None, Rholo_grid=None,
                      G_lo=0.25, G_hi=2.0, phi_op=0.35):
    """Design chart: sensitivity of the ratio to 2-OG across the operating window,
    as a function of binding affinity and bound-state separation."""
    Jcal = spectral_overlap_density(qp)
    if KD_grid is None:
        KD_grid = np.logspace(-1.3, 0.7, 60)
    if Rholo_grid is None:
        Rholo_grid = np.linspace(36.0, 62.0, 60)
    S = np.zeros((len(Rholo_grid), len(KD_grid)))
    for i, Rh in enumerate(Rholo_grid):
        for j, KD in enumerate(KD_grid):
            q = QuantumParams(**{**asdict(qp), "K_D": float(KD), "R_holo": float(Rh)})
            Jc = spectral_overlap_density(q)
            r_lo = quantum_map(G_lo, phi_op, q, Jc)["Rt"]
            r_hi = quantum_map(G_hi, phi_op, q, Jc)["Rt"]
            S[i, j] = abs(r_hi - r_lo) / (0.5 * (r_hi + r_lo))
    return KD_grid, Rholo_grid, S


if __name__ == "__main__":
    qp, cp, bp = QuantumParams(), CellParams(), BESParams()
    Jcal = spectral_overlap_density(qp)
    print(f"calibrated overlap density Jcal = {Jcal:.3e} cm")
    print(f"R0 (kappa2=2/3, n=1.33)       = {forster_R0(2/3, 1.33, qp):.2f} A")
    for mode, label in [("open", "open-loop"),
                        ("closed_naive", "closed, naive calibration"),
                        ("closed_quantum", "closed, multiscale corrected")]:
        ctrl = Controller(mode=mode)
        res = simulate(cp, bp, qp, ctrl)
        m = metrics(res, cp)
        print(f"{label:32s} PHA={m['titre_PHA']:6.2f} g/L  "
              f"fP={m['fP_final']:5.3f}  prod={m['productivity']:6.4f} g/L/h")

"""Tests for the quantum-to-cellular model.

These check the identities and invariants the paper's claims rest on, so that
a reader can confirm the framework is internally consistent without reading
every line, and so that a future change cannot quietly break a published
result. Two of them are regression tests for errors that were found and fixed
during the work; both are marked as such.

Run with:  python -m pytest tests -q
or simply: python tests/test_qcell.py
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import qcell_model as Q
from qcell_model import QuantumParams, CellParams, BESParams, Controller

QP = QuantumParams()
CP = CellParams()
BP = BESParams()
JCAL = Q.spectral_overlap_density(QP)


# ---------------------------------------------------------------- quantum layer
def test_golden_rule_equals_forster_at_R0():
    """The spectral overlap density is calibrated so the two formalisms agree.

    This is the identity that lets an excitonic coupling J be used in place of
    a Forster rate; if it fails, the quantum layer is not self-consistent.
    """
    R0 = Q.forster_R0(2 / 3, 1.33, QP)
    k_gr, _ = Q.eet_rate(2 / 3, R0, 1.33, QP, JCAL, coherent=False)
    k_forster = 1.0 / QP.tau_D0          # by definition, k = 1/tau_D at R = R0
    assert abs(k_gr - k_forster) / k_forster < 1e-10


def test_forster_sixth_power_law():
    """E = 1/(1+(R/R0)^6) must follow from the rate expression, not be imposed."""
    n = Q.lorentz_lorenz_mix(0.0, QP.n_cytosol, QP.n_granule)
    k2 = Q.kappa2_effective(0.0, QP.kappa2_static_apo, QP)
    R0 = Q.forster_R0(k2, n, QP)
    for R in (40.0, 55.0, 68.0, 90.0):
        E, _ = Q.fret_efficiency(k2, R, n, QP, JCAL, coherent=False)
        assert abs(E - 1.0 / (1.0 + (R / R0) ** 6)) < 1e-9


def test_efficiency_decreases_with_distance():
    n = Q.lorentz_lorenz_mix(0.0, QP.n_cytosol, QP.n_granule)
    k2 = Q.kappa2_effective(0.0, QP.kappa2_static_apo, QP)
    Es = [Q.fret_efficiency(k2, R, n, QP, JCAL)[0] for R in range(30, 100, 5)]
    assert all(a > b for a, b in zip(Es, Es[1:]))


def test_lorentz_lorenz_endpoints():
    assert abs(Q.lorentz_lorenz_mix(0.0, 1.36, 1.49) - 1.36) < 1e-12
    assert abs(Q.lorentz_lorenz_mix(1.0, 1.36, 1.49) - 1.49) < 1e-12
    # and it must be monotone in between
    xs = np.linspace(0, 1, 21)
    ns = [Q.lorentz_lorenz_mix(x, 1.36, 1.49) for x in xs]
    assert all(a < b for a, b in zip(ns, ns[1:]))


def test_kappa2_lies_between_static_and_isotropic():
    """Partial rotational averaging interpolates; it cannot overshoot either end."""
    for k2s in (QP.kappa2_static_apo, QP.kappa2_static_holo):
        for phi in (0.0, 0.1, 0.3, 0.5):
            k2 = Q.kappa2_effective(phi, k2s, QP)
            lo, hi = sorted((k2s, 2 / 3))
            assert lo - 1e-12 <= k2 <= hi + 1e-12


def test_crowding_freezes_the_orientation_factor():
    """More granules -> higher viscosity -> less averaging -> closer to static."""
    k2s = QP.kappa2_static_apo
    far = abs(Q.kappa2_effective(0.0, k2s, QP) - k2s)
    near = abs(Q.kappa2_effective(0.4, k2s, QP) - k2s)
    assert near < far


def test_hill_occupancy():
    assert Q.sensor_occupancy(0.0, QP) == 0.0
    assert Q.sensor_occupancy(1e6, QP) > 0.999
    half = Q.sensor_occupancy(QP.K_D, QP)
    assert abs(half - 0.5) < 1e-12        # exact at G = K_D for any Hill h


# ------------------------------------------------- the paper's central claim
def test_acceptor_band_always_better_transmitted():
    """T_A/T_D > 1 for ANY turbid sample. This is the paper's main result.

    It follows from lambda_D < lambda_A alone, so it must hold for every
    biomass and polymer loading and for every scattering exponent, not merely
    at the nominal parameters.
    """
    from dataclasses import replace
    for p in (0.5, 1.0, 1.5, 2.5, 4.0):
        qp = replace(QP, p_mie=p)
        for Cx in (0.5, 4.0, 11.0, 25.0):
            for fP in (0.0, 0.3, 0.75):
                T_D, T_A = Q.optical_attenuation(Cx, Cx * fP, qp)
                assert T_A > T_D, (p, Cx, fP)
                assert 0 < T_D <= 1 and 0 < T_A <= 1


def test_no_bias_in_a_clear_sample():
    """With nothing to scatter, the transport term must vanish identically."""
    T_D, T_A = Q.optical_attenuation(0.0, 0.0, QP)
    assert T_D == 1.0 and T_A == 1.0


def test_granule_contrast_sets_the_size_of_the_effect():
    """a_P = 0 means granules scatter no more than the biomass they displace."""
    from dataclasses import replace
    drift = []
    for aP in (0.0, 0.12, 0.5):
        qp = replace(QP, a_P=aP)
        T_D, T_A = Q.optical_attenuation(11.0, 5.0, qp)
        drift.append(T_A / T_D)
    assert all(a < b for a, b in zip(drift, drift[1:]))


# ----------------------------------------------------------- the inverse map
def test_inverse_map_is_exact():
    """REGRESSION. The forward map adds the dark offset AFTER attenuation, so
    the inverse cannot undo transport by pre-scaling the ratio by T_D/T_A --
    that shortcut is exact only for a zero dark offset and biased recovered
    2-OG by up to -10 % at high turbidity. The transmissions must enter the
    linear solve for 1-E directly.
    """
    worst = 0.0
    for Cx, Cp in [(2, 0.2), (5, 1.5), (8, 3.2), (12, 6.0), (16, 9.0), (20, 12.0)]:
        phi = float(Q.granule_volume_fraction(Cx, Cp, CP))
        for G in (0.05, 0.2, 1.0, 3.0, 10.0):
            m = Q.quantum_map(G, phi, QP, JCAL, Cx=Cx, Cp=Cp)
            G_hat, _, _ = Q.invert_quantum_map(m["Rt"], phi, QP, JCAL,
                                               Cx=Cx, Cp=Cp, correction="full")
            worst = max(worst, abs(float(G_hat) - G) / G)
    assert worst < 1e-9, f"round-trip error {worst:.2e}"


def test_uncorrected_estimate_always_reads_high():
    """The bias is one-sided: a raw ratio can over-read 2-OG but never under-read."""
    for Cx, Cp in [(4, 1.0), (11, 5.0), (16, 9.0)]:
        phi = float(Q.granule_volume_fraction(Cx, Cp, CP))
        for G in (0.1, 0.5, 2.0):
            m = Q.quantum_map(G, phi, QP, JCAL, Cx=Cx, Cp=Cp)
            G_raw, _, _ = Q.invert_quantum_map(m["Rt"], phi, QP, JCAL,
                                               Cx=Cx, Cp=Cp, correction="none")
            assert float(G_raw) > G


def test_density_only_correction_is_incomplete():
    """OD normalisation cannot remove the granule-specific part, by construction."""
    Cx, Cp, G = 11.0, 6.0, 1.0
    phi = float(Q.granule_volume_fraction(Cx, Cp, CP))
    m = Q.quantum_map(G, phi, QP, JCAL, Cx=Cx, Cp=Cp)
    g_od, _, _ = Q.invert_quantum_map(m["Rt"], phi, QP, JCAL, Cx=Cx, Cp=Cp,
                                      correction="od_only")
    g_full, _, _ = Q.invert_quantum_map(m["Rt"], phi, QP, JCAL, Cx=Cx, Cp=Cp,
                                        correction="full")
    assert abs(float(g_full) - G) < 1e-9
    assert abs(float(g_od) - G) > 1e-6


# -------------------------------------------------------------- cellular layer
def test_granule_volume_fraction_bounds():
    assert Q.granule_volume_fraction(10.0, 0.0, CP) == 0.0
    f = [Q.granule_volume_fraction(10.0, c, CP) for c in (0, 2, 5, 8)]
    assert all(a < b for a, b in zip(f, f[1:]))
    assert all(0.0 <= x <= 0.95 for x in f)


def test_packing_fraction_is_random_close_packing():
    """REGRESSION. An earlier version used 0.68, which lies between random and
    crystalline close packing, where no random monodisperse packing sits.
    """
    assert abs(QP.phi_pack - 0.64) < 1e-12


def test_clock_baseline_is_tuned_to_its_own_optimum():
    """The open-loop baseline must be run at ITS best switch time.

    REGRESSION. A baseline tuned anywhere else flatters feedback, which is the
    comparison the paper is careful to make fairly. run_cases.py asserts the
    constant still matches the scan; this checks the constant itself, cheaply,
    without re-running the scan.
    """
    import run_cases
    assert abs(run_cases.T_SWITCH_NOMINAL - 60.0) < 1e-9


def test_simulation_conserves_sanity():
    """A short run must keep every state physical."""
    ctrl = Controller(mode="open", t_switch=60.0)
    res = Q.simulate(CP, BP, QP, ctrl, t_end=24.0)
    y = res["y"]
    assert np.all(np.isfinite(y))
    assert np.all(y[:, :4] >= -1e-9), "negative concentration"
    Cx, Cp = y[:, 2], y[:, 3]
    assert np.all(Cp <= Cx + 1e-9), "more polymer than biomass"
    assert np.all((y[:, 8] >= -1e-9) & (y[:, 8] <= 1 + 1e-9)), "active fraction out of range"


def test_closed_loop_runs_and_estimates():
    ctrl = Controller(mode="closed_multiscale", G_set=0.16)
    ctrl.noise_sd = 0.0
    res = Q.simulate(CP, BP, QP, ctrl, t_end=24.0, seed=0)
    assert np.all(np.isfinite(res["y"]))
    assert res["aux"].shape[0] == res["y"].shape[0]


if __name__ == "__main__":
    fns = [(n, f) for n, f in sorted(globals().items())
           if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in fns:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"  FAIL  {name}: {e}")
    print(f"\n{len(fns) - failed}/{len(fns)} passed")
    sys.exit(1 if failed else 0)

"""Scenario runner: generates every number and figure datum used in the manuscript."""
import json
import os
import numpy as np

# outputs are written next to the code, so the scripts work from anywhere
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
from dataclasses import replace, asdict
import qcell_model as M
from qcell_model import (QuantumParams, CellParams, BESParams, Controller,
                         simulate, metrics, quantum_map, invert_quantum_map,
                         spectral_overlap_density, forster_R0, granule_volume_fraction,
                         lorentz_lorenz_mix, kappa2_effective, coulomb_coupling,
                         krieger_dougherty, optical_attenuation, eet_rate)

T_END = 96.0
G_STAR = 0.16
# The clock baseline must be run at ITS OWN optimum, or the comparison
# flatters feedback. This is the argmax of exp_switch_scan; main() asserts
# the two still agree, so the constant cannot silently drift out of date.
T_SWITCH_NOMINAL = 60.0
MODES = ["closed_raw", "closed_od", "closed_multiscale"]


# ------------------------------------------------------------------ runners
def run_clock(cp, bp, qp, t_switch=T_SWITCH_NOMINAL, t_end=T_END):
    ctrl = Controller(mode="open", t_switch=t_switch)
    return simulate(cp, bp, qp, ctrl, t_end=t_end)


def run_closed(cp, bp, qp, mode, t_end=T_END, seed=0, G_set=G_STAR, noise=0.02):
    ctrl = Controller(mode=mode, G_set=G_set)
    ctrl.noise_sd = noise
    return simulate(cp, bp, qp, ctrl, t_end=t_end, seed=seed)


def est_error(res):
    t, y, a = res["t"], res["y"], res["aux"]
    mask = (t > 45) & (y[:, 4] > 1e-6) & np.isfinite(a[:, 5])
    if mask.sum() == 0:
        return np.nan
    return float(np.mean((a[mask, 5] - y[mask, 4]) / y[mask, 4]) * 100)


# -------------------------------------------------------------- experiments
def exp_photophysics(qp):
    Jc = spectral_overlap_density(qp)
    R0 = forster_R0(2 / 3, 1.33, qp)
    J0 = coulomb_coupling(np.sqrt(2 / 3), R0, 1.33, qp)
    k_gr, _ = eet_rate(2 / 3, R0, 1.33, qp, Jc, coherent=False)
    theta0 = krieger_dougherty(0.0, qp.eta_cyto0) * qp.V_hydro / (1.380649e-23 * 298.15)
    return dict(Jcal_cm=Jc, R0_ref_A=R0, J_at_R0_cm1=J0,
                k_goldenrule_s=k_gr, k_forster_s=1.0 / qp.tau_D0,
                theta_rot_ns=theta0 * 1e9,
                kappa2_apo_eff=float(kappa2_effective(0.0, qp.kappa2_static_apo, qp)),
                kappa2_holo_eff=float(kappa2_effective(0.0, qp.kappa2_static_holo, qp)))


def exp_sensor_curve(qp, phis=(0.0, 0.15, 0.30)):
    Jc = spectral_overlap_density(qp)
    G = np.logspace(-2.3, 0.8, 160)
    out = {}
    for phi in phis:
        out[f"{phi:.2f}"] = quantum_map(G, phi, qp, Jc)["Rt"].tolist()
    return G.tolist(), out


def exp_drift(qp, cp, G_op=G_STAR):
    """Photophysical and transport drift along a representative trajectory."""
    Jc = spectral_overlap_density(qp)
    rows = []
    for fP in np.linspace(0.0, 0.75, 40):
        Cx, Cp = 11.0, 11.0 * fP
        phi = float(granule_volume_fraction(Cx, Cp, cp))
        qm = quantum_map(G_op, phi, qp, Jc, Cx=Cx, Cp=Cp)
        b = {}
        for c in ("none", "od_only", "full"):
            Gh, _, _ = invert_quantum_map(qm["Rt"], phi, qp, Jc, Cx=Cx, Cp=Cp,
                                          correction=c)
            b[c] = 100 * (Gh - G_op) / G_op
        rows.append([fP, phi, qm["n_med"], qm["kappa2_holo"], qm["R0_holo"],
                     qm["E"], qm["ratio_drift"], b["none"], b["od_only"], b["full"]])
    return np.array(rows)


def exp_switch_scan(cp, bp, qp):
    out = []
    for tsw in np.arange(16, 76, 4.0):
        r = run_clock(cp, bp, qp, t_switch=float(tsw))
        m = metrics(r, cp)
        out.append([tsw, m["titre_PHA"], m["fP_final"], m["phi_a_final"],
                    m["objective_J"], m["productivity"]])
    return np.array(out)


def exp_setpoint_scan(cp, bp, qp):
    out = []
    for gs in [0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.24, 0.28, 0.34]:
        r = run_closed(cp, bp, qp, "closed_multiscale", G_set=gs, noise=0.0)
        m = metrics(r, cp)
        out.append([gs, m["titre_PHA"], m["fP_final"], m["phi_a_final"],
                    m["objective_J"]])
    return np.array(out)


def exp_scattering_sensitivity(cp, bp, qp):
    """When is the full multiscale correction actually required?

    Sweeps the granule scattering contrast a_P. OD normalisation corrects the
    biomass term only, so its residual error grows with a_P.
    """
    out = []
    for aP in [0.0, 0.06, 0.12, 0.25, 0.40, 0.60, 0.90]:
        qx = replace(qp, a_P=aP)
        row = [aP]
        for mode in MODES:
            r = run_closed(cp, bp, qx, mode, noise=0.0)
            m = metrics(r, cp)
            row += [m["objective_J"], m["titre_PHA"], est_error(r)]
        out.append(row)
    return np.array(out)


def exp_robustness(cp, bp, qp):
    perturbations = [
        ("nominal", {}),
        ("mu_max -30%", {"mu_xs_max": cp.mu_xs_max * 0.70}),
        ("mu_max +30%", {"mu_xs_max": cp.mu_xs_max * 1.30}),
        ("Y_xs -20%", {"Y_xs": cp.Y_xs * 0.80}),
        ("Y_xs +20%", {"Y_xs": cp.Y_xs * 1.20}),
        ("mu_ps -25%", {"mu_ps_max": cp.mu_ps_max * 0.75}),
        ("K_IS -40%", {"K_IS": cp.K_IS * 0.60}),
    ]
    rows = []
    for name, kw in perturbations:
        cpx = replace(cp, **kw)
        e = dict(case=name)
        e["clock"] = metrics(run_clock(cpx, bp, qp), cpx)["objective_J"]
        for mode in MODES:
            e[mode] = metrics(run_closed(cpx, bp, qp, mode), cpx)["objective_J"]
        rows.append(e)
    return rows


def exp_mc(cp, bp, qp, n=40, seed=1):
    # n is the sample size reported throughout the paper (Fig. 4c,
    # Extended Data Table 3 and three claims in the text). Changing it
    # invalidates those; the strains are drawn sequentially from the seeded
    # generator, so a larger n extends this set rather than replacing it.
    rng = np.random.default_rng(seed)
    res = {k: [] for k in ["clock"] + MODES}
    for i in range(n):
        cpx = replace(
            cp,
            mu_xs_max=cp.mu_xs_max * float(rng.uniform(0.70, 1.30)),
            mu_ps_max=cp.mu_ps_max * float(rng.uniform(0.75, 1.25)),
            Y_xs=cp.Y_xs * float(rng.uniform(0.80, 1.20)),
            K_IS=cp.K_IS * float(rng.uniform(0.60, 1.40)),
            b_decN=cp.b_decN * float(rng.uniform(0.7, 1.3)),
        )
        res["clock"].append(metrics(run_clock(cpx, bp, qp), cpx)["objective_J"])
        for mode in MODES:
            res[mode].append(metrics(run_closed(cpx, bp, qp, mode, seed=i), cpx)["objective_J"])
    return {k: np.array(v) for k, v in res.items()}


def exp_design_map(qp, G_lo=0.06, G_hi=0.45, phi_op=0.20):
    KD = np.logspace(-2.0, 0.3, 40)
    Rh = np.linspace(38.0, 64.0, 40)
    S = np.zeros((len(Rh), len(KD)))
    for i, r in enumerate(Rh):
        for j, kd in enumerate(KD):
            q = replace(qp, K_D=float(kd), R_holo=float(r))
            Jc = spectral_overlap_density(q)
            a = quantum_map(G_lo, phi_op, q, Jc)["Rt"]
            b = quantum_map(G_hi, phi_op, q, Jc)["Rt"]
            S[i, j] = abs(b - a) / (0.5 * (a + b))
    return KD, Rh, S


# ------------------------------------------------------------------- driver
def main():
    qp, cp, bp = QuantumParams(), CellParams(), BESParams()
    store = {}

    pp = exp_photophysics(qp)
    store["photophysics"] = pp
    print("PHOTOPHYSICS")
    for k, v in pp.items():
        print(f"  {k:20s} {v:.5g}")

    print("\nSWITCH-TIME SCAN (clock control)")
    sc = exp_switch_scan(cp, bp, qp)
    t_best = float(sc[int(np.argmax(sc[:, 4])), 0])
    assert abs(t_best - T_SWITCH_NOMINAL) < 1e-9, (
        f"the clock baseline is tuned to t={T_SWITCH_NOMINAL} h but its optimum "
        f"is t={t_best} h; every clock-vs-feedback comparison below would be "
        f"unfair. Set T_SWITCH_NOMINAL = {t_best}.")
    store["switch_scan"] = sc.tolist()
    for r in sc:
        print("  t=%5.1f  PHA=%6.3f fP=%5.3f phi_a=%5.3f J=%6.3f" % tuple(r[:5]))
    print(f"  best clock: t={sc[np.argmax(sc[:,4]),0]:.1f} J={sc[:,4].max():.3f}")

    print("\nSETPOINT SCAN (closed loop)")
    sp = exp_setpoint_scan(cp, bp, qp)
    store["setpoint_scan"] = sp.tolist()
    for r in sp:
        print("  G*=%5.3f  PHA=%6.3f fP=%5.3f phi_a=%5.3f J=%6.3f" % tuple(r))

    print("\nTRAJECTORIES")
    traj = {}
    r = run_clock(cp, bp, qp)
    traj["clock"] = dict(t=r["t"][::50].tolist(), y=r["y"][::50].tolist(),
                         aux=r["aux"][::50].tolist())
    print(f"  clock              J={metrics(r,cp)['objective_J']:.3f} "
          f"PHA={metrics(r,cp)['titre_PHA']:.3f}")
    for mode in MODES:
        r = run_closed(cp, bp, qp, mode)
        m = metrics(r, cp)
        traj[mode] = dict(t=r["t"][::50].tolist(), y=r["y"][::50].tolist(),
                          aux=r["aux"][::50].tolist())
        print(f"  {mode:18s} J={m['objective_J']:.3f} PHA={m['titre_PHA']:.3f} "
              f"fP={m['fP_final']:.3f} est.err={est_error(r):+.2f}%")
    store["trajectories"] = traj

    print("\nDRIFT ALONG ACCUMULATION")
    dr = exp_drift(qp, cp)
    store["drift"] = dr.tolist()
    for r in dr[::8]:
        print("  fP=%4.2f phi=%5.3f n=%6.4f R0=%5.2f drift=%6.4f "
              "bias none/od/full = %+6.2f / %+5.2f / %+5.2f %%"
              % (r[0], r[1], r[2], r[4], r[6], r[7], r[8], r[9]))

    print("\nSCATTERING-CONTRAST SENSITIVITY (a_P)")
    ss = exp_scattering_sensitivity(cp, bp, qp)
    store["scattering"] = ss.tolist()
    print("   a_P    J_raw  J_od   J_ms   | err_raw err_od  err_ms")
    for r in ss:
        print(f"  {r[0]:5.2f}  {r[1]:6.3f} {r[4]:6.3f} {r[7]:6.3f}  | "
              f"{r[3]:+7.2f} {r[6]:+7.2f} {r[9]:+7.2f}")

    print("\nROBUSTNESS")
    rb = exp_robustness(cp, bp, qp)
    store["robustness"] = rb
    hdr = f"  {'case':14s} {'clock':>8s}" + "".join(f"{m.replace('closed_',''):>12s}" for m in MODES)
    print(hdr)
    for e in rb:
        print(f"  {e['case']:14s} {e['clock']:8.3f}" +
              "".join(f"{e[m]:12.3f}" for m in MODES))

    print("\nMONTE CARLO")
    mc = exp_mc(cp, bp, qp)
    store["monte_carlo"] = {k: v.tolist() for k, v in mc.items()}
    for k, v in mc.items():
        print(f"  {k:18s} mean={v.mean():6.3f} sd={v.std():5.3f} "
              f"p10={np.percentile(v,10):6.3f} min={v.min():6.3f}")

    print("\nDESIGN MAP")
    KD, Rh, S = exp_design_map(qp)
    store["design_map"] = dict(KD=KD.tolist(), Rh=Rh.tolist(), S=S.tolist())
    i, j = np.unravel_index(np.argmax(S), S.shape)
    print(f"  best: K_D={KD[j]:.4f} mM  R_holo={Rh[i]:.1f} A  sensitivity={S[i,j]:.3f}")

    gg, curves = exp_sensor_curve(qp)
    store["sensor_curve"] = dict(G=gg, curves=curves)

    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "results.json"), "w") as f:
        json.dump(store, f)
    print(f"\nwrote {os.path.join(OUT_DIR, 'results.json')}")


if __name__ == "__main__":
    main()

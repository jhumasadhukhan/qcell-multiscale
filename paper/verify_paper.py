"""Check every number quoted in the paper against results.json.

The claims below are transcribed from the published manuscript and
supplementary file. Each is paired with the expression that should reproduce
it, so running this script answers the question a referee actually has: does
the code still produce what the paper says it produces?

    python paper/verify_paper.py

Exits non-zero if any claim fails, so it can be used in CI.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))

RESULTS = os.path.join(HERE, os.pardir, "results", "results.json")
if not os.path.exists(RESULTS):
    sys.exit("results/results.json not found - run ./reproduce.sh first")

R = json.load(open(RESULTS))

# column layouts, from run_cases.py
#   drift      : fP phi n_med kappa2_holo R0_holo E ratio_drift bias_none bias_od bias_full
#   scattering : a_P then (J, titre, est_err) x (raw, od, multiscale)
#   setpoint   : G* titre fP phi_a J
DR = np.array(R["drift"])
SC = np.array(R["scattering"])
SP = np.array(R["setpoint_scan"])
RB = R["robustness"]
MC = {k: np.array(v) for k, v in R["monte_carlo"].items()}
nominal = next(r for r in RB if r["case"] == "nominal")
i_nom = int(np.argmin(np.abs(SC[:, 0] - 0.12)))

rows = []


def claim(where, text, published, actual, tol=0.05):
    """tol is relative, except against a published value of zero."""
    if published == 0:
        ok = abs(actual) <= tol
    else:
        ok = abs(actual - published) <= tol * abs(published)
    rows.append((where, text, published, actual, ok))


# ------------------------------------------------------------------ abstract
claim("Abstract", "bias in inferred 2-OG (%)", 29, SC[i_nom, 3])
claim("Abstract", "polymer titre lost (%)", 20,
      100 * (SC[i_nom, 8] - SC[i_nom, 2]) / SC[i_nom, 8])

# ------------------------------------------------------------ quantum layer
pp = R["photophysics"]
claim("Results", "Forster radius (Angstrom)", 58.95, pp["R0_ref_A"], tol=0.001)
claim("Results", "effective overlap width (cm-1)", 910, 1 / pp["Jcal_cm"], tol=0.01)
claim("Results", "rotational correlation time (ns)", 58, pp["theta_rot_ns"], tol=0.02)
claim("Results", "dynamically averaged fraction", 0.064,
      4.0 / (pp["theta_rot_ns"] + 4.0), tol=0.02)
claim("Results", "golden rule == Forster at R0", 0.0,
      abs(pp["k_goldenrule_s"] - pp["k_forster_s"]) / pp["k_forster_s"], tol=1e-9)

# --------------------------------------------------------------- sensor drift
claim("Results", "max cytoplasmic refractive index", 1.41, DR[:, 2].max(), tol=0.005)
claim("Results", "max transmission ratio T_A/T_D", 1.18, DR[:, 6].max(), tol=0.01)
claim("Results", "raw over-read, low end (%)", 15, DR[:, 7].min(), tol=0.08)
claim("Results", "raw over-read, high end (%)", 21, DR[:, 7].max(), tol=0.08)
claim("Results", "granule-specific residual at fP=0.62 (%)", -2.6,
      float(np.interp(0.62, DR[:, 0], DR[:, 8])), tol=0.10)
claim("Results", "objective lost to the raw ratio (%)", 5.6,
      100 * (nominal["closed_multiscale"] - nominal["closed_raw"])
      / nominal["closed_multiscale"])
claim("Results", "multiscale residual across the sweep (%)", 0.0,
      float(np.abs(SC[:, 9]).max()), tol=1e-3)

# ----------------------------------------------------- when the layer matters
claim("Results", "density-only residual at nominal a_P (%)", 0.4,
      abs(SC[i_nom, 6]), tol=0.30)
claim("Results", "density-only residual at a_P = 0 (%)", -4.4, SC[0, 6], tol=0.08)
claim("Results", "density-only worst error (%)", 30, SC[:, 6].max(), tol=0.08)
claim("Results", "objective loss, density-only, max contrast (%)", 6,
      100 * (SC[:, 7].max() - SC[:, 4].min()) / SC[:, 7].max(), tol=0.12)

# ------------------------------------------------------------------- control
best = SP[int(np.argmax(SP[:, 4]))]
claim("Results", "best 2-OG setpoint (mM)", 0.16, best[0], tol=0.02)
claim("Results", "best closed-loop objective", 7.456, float(best[4]), tol=0.002)
claim("Results", "best clock objective", 7.439, nominal["clock"], tol=0.002)
claim("Results", "closed-loop polymer (g/L)", 6.17, best[1], tol=0.01)

jc = np.array([r["clock"] for r in RB])
jm = np.array([r["closed_multiscale"] for r in RB])
claim("Results", "mean gain, 7 perturbations (%)", 3.5,
      100 * (jm.mean() - jc.mean()) / jc.mean(), tol=0.10)
claim("Results", "worst-case gain, 7 perturbations (%)", 17.9,
      100 * (jm.min() - jc.min()) / jc.min(), tol=0.10)
claim("Results", "mean gain, 40 strains (%)", 3.4,
      100 * (MC["closed_multiscale"].mean() - MC["clock"].mean()) / MC["clock"].mean(),
      tol=0.10)
claim("Results", "closed loop beats clock (of 40)", 24,
      int((MC["closed_multiscale"] > MC["clock"]).sum()), tol=1e-9)
# the paper quotes population s.d. (ddof=0)
claim("Results", "s.d. clock", 0.64, MC["clock"].std(ddof=0), tol=0.02)
claim("Results", "s.d. closed loop", 0.81, MC["closed_multiscale"].std(ddof=0), tol=0.02)
claim("Results", "multiscale beats raw (of 40)", 35,
      int((MC["closed_multiscale"] > MC["closed_raw"]).sum()), tol=1e-9)
claim("Results", "multiscale mean gain over raw (%)", 4.6,
      100 * (MC["closed_multiscale"].mean() - MC["closed_raw"].mean())
      / MC["closed_raw"].mean(), tol=0.08)

# --------------------------------------------------------------- sensor design
dm = R["design_map"]
S, KD, Rh = np.array(dm["S"]), np.array(dm["KD"]), np.array(dm["Rh"])
i_rh, i_kd = np.unravel_index(int(np.argmax(S)), S.shape)
claim("Results", "best K_D (mM)", 0.077, float(KD[i_kd]), tol=0.02)

# ------------------------------------------------------------------- report
w = max(len(t) for _, t, _, _, _ in rows)
print(f"{'where':<10} {'claim':<{w}} {'published':>10} {'computed':>12}   ")
print("-" * (w + 40))
bad = 0
for where, text, pub, act, ok in rows:
    if not ok:
        bad += 1
    print(f"{where:<10} {text:<{w}} {pub:>10g} {act:>12.4g}   "
          f"{'ok' if ok else '** MISMATCH **'}")
print("-" * (w + 40))
print(f"{len(rows) - bad}/{len(rows)} claims verified against {os.path.relpath(RESULTS)}")
sys.exit(1 if bad else 0)

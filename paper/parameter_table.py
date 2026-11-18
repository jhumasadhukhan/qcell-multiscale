"""Dump every model parameter with its value, as Supplementary Tables 4-6.

Supplementary Note 3 states that the parameter tables are generated directly
from the model's own objects rather than retyped. This is that generator: it
walks the three dataclasses and prints every field, so a reader can confirm
that the value printed in the paper is the value the code actually used.

What it does NOT reproduce is the provenance column, the literature citations
and the audit flags (A-dagger, L-double-dagger) of the published tables. Those
are editorial judgements about where each number came from, they are recorded
in Supplementary Tables 4-6 and in Supplementary Note 3.5, and they cannot be
derived from the code. Use the paper for provenance and this script to check
the numbers.

    python paper/parameter_table.py            human-readable
    python paper/parameter_table.py --csv      machine-readable
"""
import os
import sys
from dataclasses import fields

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir))

import qcell_model as Q

TABLES = [
    ("Supplementary Table 4 | Quantum-layer parameters", Q.QuantumParams()),
    ("Supplementary Table 5 | Cellular-layer parameters", Q.CellParams()),
    ("Supplementary Table 6 | Bioelectrochemical supply parameters", Q.BESParams()),
]

# the handful of quantities that are computed rather than set
DERIVED = [
    ("Jcal", "Calibrated spectral overlap density (cm)",
     lambda: Q.spectral_overlap_density(Q.QuantumParams())),
    ("R0_ref", "Forster radius at the isotropic reference point (Angstrom)",
     lambda: Q.forster_R0(2 / 3, 1.33, Q.QuantumParams())),
]


def fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}" if 1e-3 <= abs(v) < 1e5 or v == 0 else f"{v:.4e}"
    return str(v)


def main():
    as_csv = "--csv" in sys.argv
    if as_csv:
        print("table,parameter,value")
    for title, obj in TABLES:
        if not as_csv:
            print(f"\n{title}")
            print(f"{'parameter':<22} {'value':>14}")
            print("-" * 38)
        for f in fields(obj):
            v = fmt(getattr(obj, f.name))
            if as_csv:
                print(f'"{title.split("|")[0].strip()}",{f.name},{v}')
            else:
                print(f"{f.name:<22} {v:>14}")

    if not as_csv:
        print("\nDerived quantities (computed, not set)")
        print("-" * 38)
    for name, desc, fn in DERIVED:
        v = fmt(fn())
        if as_csv:
            print(f'"derived",{name},{v}')
        else:
            print(f"{name:<22} {v:>14}   {desc}")

    n = sum(len(fields(o)) for _, o in TABLES)
    if not as_csv:
        print(f"\n{n} parameters. Provenance codes and literature sources for "
              f"every one are in Supplementary Tables 4-6;\nthe independent "
              f"audit of those sources is Supplementary Note 3.5.")


if __name__ == "__main__":
    main()

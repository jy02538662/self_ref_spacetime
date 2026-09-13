"""Probe: does the DISCRETE crossed product (observer = add clock J, group Z_4)
give Type II?  -- the core blocker of Type III -> II, made concrete.

Core blocker (2026-09-12): Witten's crossed product works from Type III (continuous
field theory, no trace).  Quantum tide's D is Type I (finite-dim).  Direction 3
("observer = add a clock") would be the DISCRETE crossed product: couple the clock
J = [[0,1],[-1,0]] into D.

Key check: the discrete crossed product  A (x) Z_4  (A = M_2, Z_4 = <J>, J^4 = I)
is FINITE-dimensional (16-dim), so it is Type I (has minimal projections and a
trace) -- it does NOT give Type II.  Type II needs a CONTINUOUS clock (R-translation),
which is exactly the "discrete -> continuous" obstacle (Wielandt-Wintner: finite-dim
cannot give the continuous spectrum / [X,P]=i hbar).

So this probe SITS the blocker precisely: "observer = add clock" in the discrete
setting stays Type I; Type II requires the continuous clock, i.e. the same
"discrete -> continuous" step that the i-hbar litmus test already located.

Code: `py -m experiments.exp_crossed_product`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

J = np.array([[0, 1], [-1, 0]], complex)   # clock, J^2 = -I, J^4 = I  -> Z_4


def crossed_product_dim(A_dim, group_order):
    """Crossed product A (x) G has dim = dim(A) * |G| (as a vector space)."""
    return A_dim * group_order


def main():
    print("=== probe: discrete crossed product (observer = clock J, Z_4) -> Type II? ===")
    print()

    # A = M_2 (2x2 matrices, dim 4);  Z_4 = {I, J, J^2, J^3}
    print("A = M_2 (finite-dim, Type I);  clock J with J^4 = I -> group Z_4")
    J2 = J @ J
    J4 = J2 @ J2
    print(f"  J^2 = {J2.tolist()}  (=-I),   J^4 = {np.round(J4, 3).tolist()}  (=I) -> order 4")
    print()

    # crossed product A (x) Z_4: dimension = dim(A) * |Z_4| = 4 * 4 = 16
    dim = crossed_product_dim(4, 4)
    print(f"Part A. discrete crossed product A (x) Z_4:")
    print(f"  dim = dim(A) * |Z_4| = 4 * 4 = {dim}  -> FINITE-dimensional")
    print(f"  => Type I (has minimal projections, has a trace); entropy is finite but")
    print(f"     the algebra is the SAME kind as the starting finite-dim D (just bigger).")
    print()

    # contrast: continuous clock R
    print("Part B. contrast -- Witten's crossed product A (x) R (continuous clock):")
    print("  group R is CONTINUOUS -> the crossed product is INFINITE-dimensional,")
    print("  has no minimal projections, and acquires a (semi)finite trace -> Type II.")
    print()

    print("interpretation:")
    print("  - 'observer = add a clock' in the DISCRETE setting (J, Z_4) stays Type I:")
    print("    the crossed product is just a bigger finite-dim algebra (16-dim).")
    print("  - Type II needs the CONTINUOUS clock (R-translation), which is infinite-dim,")
    print("    i.e. the same 'discrete -> continuous' step the i-hbar litmus test located.")
    print("  - so Type III->II is NOT reachable by a finite clock; its blocker IS")
    print("    'discrete -> continuous' (observer internalization = infinite-dim, no shortcut).")

    summary = {
        "discrete_crossed_product_dim": dim,
        "discrete_crossed_product_type": "Type I (finite-dim)",
        "type_II_requires": "continuous clock (R), i.e. discrete->continuous",
        "note": "observer=add-clock (discrete J, Z_4) stays Type I; Type II needs the "
                "continuous clock = the discrete->continuous obstacle (Wielandt-Wintner).",
    }
    out = ROOT / "experiments" / "exp_crossed_product_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

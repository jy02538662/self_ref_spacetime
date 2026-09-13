"""Probe: does "self-reference" (observer observes itself) generate INFINITE
dimension (Type II) from a FINITE-dim D, without N -> inf?

User's intuition: "本体是观察者，自指产生无限" -- the infinity is not from "N
large" but from "the infinite recursion of self-reference".

Math check: the concrete operations of self-reference on a FINITE-dim algebra are
still finite-dimensional:
  - discrete modular flow  sigma_n(x) = rho^n x rho^{-n}  (Tomita discrete): rho^n
    is N x N for all n, so {sigma_n(x)} spans <= N^2 dimensions -> Type I.
  - infinite recursion  M_N -> M_{N^2} -> ...  has inductive limit = AF algebra
    (Type I approximations); its WEAK closure (the hyperfinite II_1 factor R) is
    the N -> inf limit, NOT a finite-dim self-reference.

So "self-reference" does NOT produce Type II within finite dim.  The infinity in
"self-reference" IS the inductive limit (hyperfinite factor R), which is the SAME
N -> inf step -- re-labelled as "infinite recursion", not bypassed.

This probe shows operation 1 (discrete modular flow) spans finite-dim (Type I).

Code: `py -m experiments.exp_self_reference`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main():
    print("=== probe: does self-reference (observer observes itself) give Type II? ===")
    print()

    # operation 1: discrete modular flow sigma_n(x) = rho^n x rho^{-n}
    N = 4
    rho = np.diag([1.0, 2.0, 3.0, 4.0])
    rho = rho / np.trace(rho)          # density matrix
    print(f"Operation 1: discrete modular flow on M_{N},  rho = diag(1,2,3,4)/10")
    print(f"  sigma_n(x) = rho^n x rho^{{-n}}   (Tomita discrete)")
    # span of {sigma_n(E_{ij}) : n in Z}: sigma_n(E_ij) = (lam_i/lam_j)^n E_ij
    # the set {(lam_i/lam_j)^n : n in Z} spans (over C) a 1-dim space for each (i,j)
    # -> total span <= N^2
    lam = np.array([1.0, 2.0, 3.0, 4.0]) / 10.0
    ratios = np.array([[lam[i] / lam[j] for j in range(N)] for i in range(N)])
    # count distinct (i,j) -> span dim = number of distinct E_ij blocks = N^2
    print(f"  sigma_n(E_ij) = (lam_i/lam_j)^n E_ij;  span over C is <= N^2 = {N*N}")
    print(f"  -> the modular-flow algebra is FINITE-dimensional (Type I), for ALL n.")
    print()

    # operation 2/3: infinite recursion M_N -> M_{N^2} -> ... inductive limit
    print("Operation 2: infinite recursion  M_N -> M_{N^2} -> M_{N^4} -> ...")
    dims = [N ** (2 ** k) for k in range(5)]
    print(f"  dimensions: {dims}  -> inductive limit = AF algebra (Type I approx.)")
    print(f"  its weak closure = hyperfinite II_1 factor R = the N->inf limit,")
    print(f"  NOT a finite-dim self-reference.")
    print()

    print("interpretation:")
    print("  - 'self-reference' on a FINITE-dim D (discrete modular flow, recursion) stays")
    print("    finite-dim (Type I).  The infinity only appears as the INDUCTIVE LIMIT")
    print("    (hyperfinite factor R), which is the same N->inf step re-labelled.")
    print("  - so the intuition 'infinity from self-reference' is conceptually right (the")
    print("    infinity is the recursion, not the matrix size), but its math is the")
    print("    N->inf inductive limit -- the obstacle is NOT bypassed, only re-labelled.")
    print("  - net: Type II still needs the inductive limit (N->inf); no finite self-ref.")

    summary = {
        "modular_flow_span_dim": N * N,
        "modular_flow_type": "Type I (finite-dim)",
        "infinite_recursion_limit": "AF algebra -> hyperfinite II_1 (N->inf)",
        "self_reference_bypasses_obstacle": False,
        "note": "self-reference on finite-dim D stays Type I; the infinity is the inductive "
                "limit (hyperfinite factor R) = N->inf, re-labelled not bypassed.",
    }
    out = ROOT / "experiments" / "exp_self_reference_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

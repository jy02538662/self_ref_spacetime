"""Layer 3, route 2, step 2: q-CG of V_{1/2} (x) V_{1/2} -- and a discovery.

The q-deformed Clebsch-Gordan of  V_{1/2} (x) V_{1/2} = V_0 (+) V_1  is the base of
the S^2_q fusion table.  This probe verifies the q-singlet and reveals a KEY fact:

  DISCOVERY: the q-deformation BREAKS standard-inner-product orthogonality.
  The q-singlet  |s_q> = |ud> - q|du>  (annihilated by Delta(E),Delta(F)) and the
  q-triplet m=0 state  |1,0> = Delta(F)|uu> / norm  are NOT orthogonal under the
  STANDARD inner product:  <s_q | 1,0> = 1 - q^{-2} != 0  for |q|=1, q != 1.

  Reason: Delta(F) = F(x)K^{-1} + 1(x)F is NOT Hermitian, so its image is not the
  standard-orthogonal complement of the kernel.  Orthogonality of the fusion
  decomposition requires the U_q(su2)-INVARIANT inner product (or the JW idempotent
  / quantum trace), NOT the standard one.  At q -> 1 the standard orthogonality is
  recovered.

This matches the earlier finding (exp_loop_lattice): the q-deformation lives in the
STRUCTURE (weights, coproduct, invariant inner product), not in a singlet vector.

Convention: K|up>=q|up>, K|down>=q^{-1}|down>;  Delta(E)=E(x)1+K(x)E,
Delta(F)=F(x)K^{-1}+1(x)F.  q = e^{i pi/(k+2)}.

Code: `py -m experiments.exp_s2q_cg`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

I2 = np.eye(2)


def rep12(q):
    K = np.diag([q, 1 / q])
    E = np.array([[0, 1], [0, 0]], complex)
    F = np.array([[0, 0], [1, 0]], complex)
    return K, E, F


def delta_E(q):
    K, E, F = rep12(q)
    return np.kron(E, I2) + np.kron(K, E)


def delta_F(q):
    K, E, F = rep12(q)
    return np.kron(F, np.linalg.inv(K)) + np.kron(I2, F)


def norm(v):
    return v / np.linalg.norm(v)


def main():
    print("=== layer 3 route 2 step 2: q-CG of V_{1/2}(x)V_{1/2} + orthogonality discovery ===")
    print()
    rows = []
    for k in (2, 3, 4, 5, 10):
        q = np.exp(1j * np.pi / (k + 2))
        dE = delta_E(q)
        dF = delta_F(q)

        s = norm(np.array([0, 1, -q, 0], complex))        # q-singlet
        annE = float(np.linalg.norm(dE @ s))
        annF = float(np.linalg.norm(dF @ s))

        t0 = norm(dF @ np.array([1, 0, 0, 0], complex))   # Delta(F)|uu> normalized
        overlap = abs(np.vdot(s, t0))                     # standard inner product
        expected = abs(1 - q ** -2) / 2.0                 # <s|1,0> = (1 - q^{-2})/2 (normalized)
        rows.append((k, round(float(overlap), 4), round(float(expected), 4)))

        print(f"  k={k}: q-singlet annihilated (dE={annE:.0e}, dF={annF:.0e})  "
              f"<s|1,0>_std = {overlap:.4f}  (=(1-q^-2)/2 = {expected:.4f})  "
              f"-> {'NOT orthogonal' if overlap > 1e-6 else 'orthogonal'}")

    print()
    print("classical limit q -> 1:")
    print("  <s|1,0>_std = (1-q^{-2})/2 -> 0 as q->1: standard orthogonality recovered.")
    print()
    print("interpretation:")
    print("  - V_{1/2}(x)V_{1/2} = V_0 (+) V_1: q-singlet is annihilated by Delta(E),Delta(F).")
    print("  - DISCOVERY: the fusion decomposition is NOT orthogonal under the STANDARD")
    print("    inner product (Delta(F) non-Hermitian).  Orthogonality needs the U_q(su2)")
    print("    invariant inner product / JW idempotent / quantum trace.")
    print("  - this is the real difficulty of the S^2_q fusion table (q-6j): the projection")
    print("    operators are the Jones-Wenzl idempotents, not standard orthogonal projectors.")
    print("  - next: build the JW idempotents f_j (already in exp_jones_wenzl) as the")
    print("    S^2_q fusion projections, then the V_1(x)V_1 noncommutative multiplication.")

    summary = {
        "fusion": "V_{1/2}(x)V_{1/2} = V_0 (+) V_1",
        "q_singlet_annihilated": True,
        "discovery": "standard inner product orthogonality is BROKEN by q-deformation "
                     "(<s|1,0> = (1-q^{-2})/2 != 0); recovered at q->1",
        "rows": [{"k": r[0], "overlap_std": r[1], "expected": r[2]} for r in rows],
        "note": "q-deformed fusion needs the invariant inner product / JW idempotents, not "
                "standard orthogonal projectors. Basis of the S^2_q (q-6j) fusion table.",
    }
    out = ROOT / "experiments" / "exp_s2q_cg_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

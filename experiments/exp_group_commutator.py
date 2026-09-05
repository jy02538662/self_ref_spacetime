"""Group commutator as a candidate linking invariant: L = 1 - (1/2) Tr([W1, W2]).

Candidate (correct as a "non-commutativity" measure, NOT a linking number):
    L(C1, C2) = 1 - (1/2) Tr( W(C1) W(C2) W(C1)^-1 W(C2)^-1 )
  - U(1) (abelian): W(C1), W(C2) commute  ->  L = 0  (always).
  - SU(2) (non-abelian): generically they do NOT commute  ->  L != 0.

Honest subtlety (must be stated): L measures whether the two Wilson loops COMMUTE, not
their LINKING NUMBER.  For generic (non-flat) SU(2) phases, two loop holonomies are
generic SU(2) matrices and generically do NOT commute, so L != 0 REGARDLESS of whether
the graph is "2D" or "3D".  The true "2D=0, 3D!=0" (Hopf linking) is a property of the
EMBEDDING and needs the Jones polynomial / quantum trace, NOT this commutator.

This script verifies:
  Part A: for random SU(2) matrices, L is generically in (0, 2) (non-zero), mean > 0.
  Part B: for U(1) (abelian) phases, L = 0 exactly.
  Part C: L vanishes only when the two holonomies commute (abelian / flat), so it does
          NOT distinguish 2D from 3D -- it distinguishes abelian from non-abelian.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from numpy.random import default_rng

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def random_su2(rng):
    q = rng.standard_normal(4)
    q = q / np.linalg.norm(q)
    a0, ax, ay, az = q
    return np.array(
        [[a0 + 1j * az, ay + 1j * ax],
         [-ay + 1j * ax, a0 - 1j * az]],
        dtype=complex,
    )


def commutator_L(A, B):
    """L = 1 - (1/2) Tr(A B A^-1 B^-1), with A^-1 = A^dag for SU(2)."""
    C = A @ B @ A.conj().T @ B.conj().T
    return 1.0 - 0.5 * np.real(np.trace(C))


if __name__ == "__main__":
    rng = default_rng(0)
    print("=== Group commutator L = 1 - (1/2) Tr([W1, W2]) ===")
    print()

    # Part A: SU(2) — L generically non-zero
    Ls = []
    for _ in range(5000):
        A = random_su2(rng)
        B = random_su2(rng)
        Ls.append(commutator_L(A, B))
    Ls = np.array(Ls)
    print("Part A. SU(2) (non-abelian), 5000 random pairs (W1, W2):")
    print(f"  L = 1 - (1/2)Tr([W1,W2]):  mean = {Ls.mean():.4f},  std = {Ls.std():.4f},  min = {Ls.min():.4f},  max = {Ls.max():.4f}")
    frac_nonzero = float(np.mean(Ls > 1e-6))
    print(f"  fraction with L > 0 = {frac_nonzero:.3f}   (generically NON-zero)")
    print()

    # Part B: U(1) — L = 0 exactly (abelian, commute)
    # represent U(1) as diagonal SU(2) matrices (phases), which commute
    max_L_u1 = 0.0
    for _ in range(5000):
        th1, th2 = rng.uniform(0, 2 * np.pi, 2)
        A = np.diag([np.exp(1j * th1), np.exp(-1j * th1)])
        B = np.diag([np.exp(1j * th2), np.exp(-1j * th2)])
        max_L_u1 = max(max_L_u1, commutator_L(A, B))
    print("Part B. U(1) (abelian), diagonal phases (commute):")
    print(f"  max L over 5000 pairs = {max_L_u1:.2e}   (L = 0 exactly, abelian)")
    print()

    # Part C: L vanishes iff W1, W2 commute
    print("Part C. L is a NON-commutativity measure, not a linking number:")
    print("  - L = 0  <=>  W1, W2 commute  (abelian / flat / U(1))")
    print("  - L != 0  <=>  W1, W2 non-abelian (generic SU(2) phases)")
    print("  => L distinguishes abelian vs non-abelian, NOT 2D vs 3D.")
    print("     (2D=0 / 3D!=0 for HOPF LINKING needs the Jones polynomial / embedding,")
    print("      not this commutator.)")
    print()

    summary = {
        "SU2_L_mean": round(float(Ls.mean()), 4),
        "SU2_L_std": round(float(Ls.std()), 4),
        "SU2_fraction_nonzero": round(float(frac_nonzero), 4),
        "U1_max_L": round(float(max_L_u1), 12),
        "note": "L = non-commutativity measure: 0 for abelian, non-zero for non-abelian; NOT a linking number (2D/3D needs Jones/embedding).",
    }
    out = ROOT / "experiments" / "exp_group_commutator_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

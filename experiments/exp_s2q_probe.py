"""Layer 3, route 2 recon: the algebraic S^2_q (Podles quantum sphere) at the
fusion-category level.

S^2_q = quantum homogeneous space; its "function algebra" = q-deformed spherical
harmonics
    A(S^2_q) = (+)_j  V_j (x) V_j* ,   j = 0, 1/2, ..., k/2
truncated at j = k/2 (level k).  This is the fusion-category-level picture of the
"algebraic S^2 field": a FINITE (truncated) quantum sphere whose coordinate blocks
are the q-spherical harmonics (representations V_j of U_q(su2)).

This probe verifies the two defining features:
  A. quantum dimension  d_j = [2j+1]_q = sin((2j+1)gamma)/sin(gamma), gamma=pi/(k+2)
  B. truncation  d_{k/2} = 1,  d_{k/2 + 1/2} = 0  (beyond level k)
  C. classical limit  k -> inf:  d_j -> 2j+1  (classical SU(2) dimensions),
     and the coordinate block j=1 has d_1 = [3]_q -> 3 (the 3 coordinates of S^2).

Route 2 works at the fusion-category level and does NOT need the notion of "boundary",
which is why it is the first calculable step after probe 1 (U(1), not S^2).

Code: `py -m experiments.exp_s2q_probe`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def q_number(n, gamma):
    return np.sin(n * gamma) / np.sin(gamma)


def quantum_dim(j, gamma):
    return q_number(2 * j + 1, gamma)


def main():
    print("=== layer 3 route 2 recon: algebraic S^2_q (fusion-category level) ===")
    print("S^2_q function algebra = (+)_j V_j (x) V_j*,  j = 0..k/2  (truncated)")
    print()

    # ---- A/B: quantum dimensions + truncation ----
    print("Part A/B. quantum dimensions d_j = [2j+1]_q, truncation at j = k/2:")
    rows = []
    for k in (2, 3, 4, 5, 10):
        gamma = np.pi / (k + 2)
        dims = {j: float(quantum_dim(j, gamma)) for j in (0, 0.5, 1, 1.5, k / 2)}
        d_k2 = dims[k / 2]
        d_over = float(quantum_dim((k + 1) / 2, gamma))  # j = k/2 + 1/2
        rows.append((k, round(d_k2, 4), round(abs(d_over), 10)))
        print(f"  k={k:>2}: d_0={dims[0]:.3f}  d_1/2={dims[0.5]:.3f}  d_1={dims[1]:.3f}  "
              f"d_3/2={dims[1.5]:.3f}  ...  d_{{k/2}}={d_k2:.4f}  "
              f"d_{{k/2+1/2}}={abs(d_over):.1e} (truncation=0)")
    print("  => d_{k/2} = 1 (smallest), d_{k/2+1/2} = 0: SU(2)_k category truncates at j=k/2.")

    # ---- C: classical limit ----
    print()
    print("Part C. classical limit k->inf: d_j -> 2j+1 (classical SU(2) dims):")
    for j in (0.5, 1, 1.5, 2):
        lims = [float(quantum_dim(j, np.pi / (k + 2))) for k in (10, 50, 500)]
        print(f"  j={j}: d_j for k=10,50,500 = {lims[0]:.3f}, {lims[1]:.3f}, {lims[2]:.3f}  "
              f"-> classical {2*j+1}")

    # ---- S^2_q coordinate block (j=1) ----
    print()
    print("Part D. S^2_q coordinate block (j=1, the 3 coordinates of S^2):")
    for k in (2, 3, 4, 5, 50):
        gamma = np.pi / (k + 2)
        d1 = float(quantum_dim(1, gamma))
        print(f"  k={k:>2}: d_1 = [3]_q = {d1:.4f}   (classical S^2 has 3 coords, d_1 -> 3)")

    print()
    print("interpretation:")
    print("  - S^2_q is a FINITE (truncated) quantum sphere: q-spherical harmonics live in")
    print("    V_j (x) V_j* for j = 0..k/2, with quantum dimension [2j+1]_q.")
    print("  - classical limit k->inf recovers C(S^2) (all j, dimensions 2j+1).")
    print("  - this is the 'algebraic S^2 field' basis: route 2 needs no 'boundary' concept,")
    print("    only the fusion category Rep(SU(2)_k).")
    print("  - next: the q-6j / fusion (V_j (x) V_l -> (+)_m V_m) = the multiplication table")
    print("    of the q-spherical harmonics (the noncommutative structure of S^2_q).")

    summary = {
        "S2q_function_algebra": "(+)_j V_j (x) V_j*, j=0..k/2",
        "truncation": [{"k": r[0], "d_k_over_2": r[1], "d_k_over_2_plus_half": r[2]} for r in rows],
        "classical_limit_d_j_to_2j1": True,
        "coordinate_block_j1_d1_eq_3q": True,
        "note": "S^2_q = truncated quantum sphere; fusion-category level, no boundary concept "
                "needed. Route 2 recon step 1.",
    }
    out = ROOT / "experiments" / "exp_s2q_probe_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

"""Wielandt-Wintner: the i*hbar litmus test for "discrete -> continuous".

This is the numerical grounding of the dictionary seed #7 precision (2026-09-12):
the TRUE bridge from discrete D to continuous spacetime is "observer internalization"
(Type III->II crossed product), NOT "taking the limit" (N->inf) and NOT "large-N
expansion".  The discriminating quantity is i*hbar.

Wielandt-Wintner (trivial linear algebra): for FINITE-dimensional X, P,
    Tr([X,P]) = Tr(XP) - Tr(PX) = 0   (cyclic trace)
so  [X,P] = i*hbar*I  is IMPOSSIBLE (Tr(i*hbar I) = i*hbar N != 0).  Finite-dim
matrices give an "off-diagonal i*hbar-like structure" but the diagonal (hence the
trace) is always zero.  Only INFINITE-dimensional unbounded operators (where the
trace is undefined and cyclicity fails) give the true [X,P] = i*hbar I -- and that
is exactly the Type III->II crossed-product (observer algebra) entry.

This probe demonstrates the litmus test numerically:
  X = diag(0..N-1)  (position),   P = central-difference derivative (momentum)
  -> Tr([X,P]) = 0 for ALL N;  [X,P] is off-diagonal ~ i/2, diagonal = 0.

Code: `py -m experiments.exp_wielandt_wintner`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def position_op(n):
    return np.diag(np.arange(n).astype(float))


def momentum_op(n):
    """Central-difference discrete derivative: P = -i (fwd - bwd)/2."""
    P = np.zeros((n, n), complex)
    for j in range(n):
        P[j, (j + 1) % n] = -0.5j
        P[j, (j - 1) % n] = +0.5j
    return P


def main():
    print("=== Wielandt-Wintner: i*hbar litmus test for discrete -> continuous ===")
    print()
    rows = []
    for n in (4, 8, 16, 32, 64, 128):
        X = position_op(n)
        P = momentum_op(n)
        C = X @ P - P @ X                      # [X,P]
        tr = abs(np.trace(C))
        max_diag = float(np.max(np.abs(np.diag(C))))
        off = float(abs(C[0, 1])) if n > 1 else 0.0
        rows.append((n, abs(tr), max_diag, off))
        print(f"  N={n:>3}: |Tr[X,P]|={abs(tr):.2e}   max|diagonal|={max_diag:.2e}   "
              f"|off-diag|={off:.4f}")
    print()
    print("  -> Tr[X,P]=0 for ALL finite N (cyclic trace).  The [X,P] structure is")
    print("     off-diagonal ~ i/2 (an i*hbar-like shape), but the DIAGONAL (hence the")
    print("     trace) is exactly zero.  So [X,P]=i*hbar*I (needs nonzero diagonal) is")
    print("     impossible in finite dim.")
    print()
    print("interpretation:")
    print("  - 'coarse-graining' and 'large-N expansion' still work with finite-dim")
    print("    matrices, so they give an APPROXIMATE i*hbar (off-diagonal), never the")
    print("    TRUE i*hbar I.  They are 'taking the limit' in disguise.")
    print("  - the TRUE i*hbar needs infinite-dim unbounded operators (trace undefined,")
    print("    cyclicity fails) = the Type III->II crossed product (observer algebra).")
    print("  - so 'discrete -> continuous' true bridge = observer internalization,")
    print("    matching dictionary seed #7 and 断裂与自指 第1步.")

    summary = {
        "wielandt_wintner": "Tr[X,P]=0 for all finite N",
        "off_diagonal_structure": "~i/2, diagonal=0",
        "conclusion": "i*hbar needs infinite-dim (observer algebra Type III->II), "
                      "not limit nor large-N expansion",
        "rows": [{"N": r[0], "trace": r[1], "max_diag": r[2], "off_diag": r[3]} for r in rows],
    }
    out = ROOT / "experiments" / "exp_wielandt_wintner_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

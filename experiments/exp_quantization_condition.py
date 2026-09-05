"""Where does the crossing amplitude A come from? classical vs quantization.

Honest status of "layer 3": the FORWARD chain (braid word -> Jones) is fully
closed, but the amplitude A in  sigma = A*1 + A^{-1}*e  is NOT fixed by classical
self-consistency -- it is a QUANTIZATION input (a root of unity).  This script
makes that boundary precise and verifies it numerically.

Convention (same as exp_jones / exp_pairing_to_braid / exp_braid_word_to_jones):
  A = q^{1/4},  q = e^{i pi/(k+2)}  (SU(2)_k root of unity),
  d = -A^2 - A^{-2} = -2 cos(pi/(2(k+2)))  (the Kauffman loop value).

Part A: classical self-consistency does NOT fix A.
        The skein inverse condition R R^{-1} = I forces d = -A^2 - A^{-2}, but
        this identity holds for ANY A != 0.  Verify with random A's.
        => A is a FREE classical parameter; self-consistency cannot grow it.

Part B: quantization = root of unity.  q^{k+2} = 1 fixes A = e^{i pi/(4(k+2))}.
        The quantum dimension of the fundamental rep is  -d = A^2 + A^{-2}
        = 2 cos(pi/(2(k+2))).  Tabulate k=1..5.

Part C: classical limit k -> inf:  A -> 1,  d -> -2 (classical dim 2).
        The "quantum" content = the finite-k deviation of -d from 2; at k=1 the
        quantum dimension collapses to sqrt(3) ~ 1.732, below the classical 2.

Part D: conclusion -- A is representation-theoretic data of SU(2)_k (a quantum
        dimension / root of unity), NOT a classical pairing weight.  Fixing A to
        a root of unity = choosing the SU(2)_k category = "quantum coherence".
        This is the one switch layer 3 still needs, and it is NOT derivable from
        classical pairing self-consistency (Part A shows A is free).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def skein_inverse_error(A: complex) -> float:
    """|R R^{-1} - I| for R = A*1 + A^{-1}*e, e = TL generator (e^2 = d e)."""
    d = -(A ** 2) - A ** (-2)
    eps = np.array([0.0, 1.0, -1.0, 0.0], dtype=complex)
    e = (d / 2.0) * np.outer(eps, eps.conj())
    R = A * np.eye(4) + (A ** -1) * e
    Rinv = (A ** -1) * np.eye(4) + A * e
    return float(np.max(np.abs(R @ Rinv - np.eye(4))))


if __name__ == "__main__":
    print("=== Where does A come from? classical self-consistency vs quantization ===")
    print()

    # Part A: classical self-consistency does NOT fix A
    print("Part A. skein inverse holds for ANY A != 0 (A is a free classical parameter):")
    rng = np.random.default_rng(0)
    rows = []
    for _ in range(5):
        # random A (arbitrary modulus & phase, not a root of unity)
        A = complex(rng.uniform(0.3, 2.0), rng.uniform(-1.0, 1.0))
        err = skein_inverse_error(A)
        rows.append((A, err))
        print(f"  A = {A.real:+.4f}{A.imag:+.4f}i   -> |R R^-1 - I| = {err:.2e}")
    all_free = all(err < 1e-9 for _, err in rows)
    print(f"  => R R^-1 = I for every A; d = -A^2 - A^{-2} is an identity, not a constraint.")
    print(f"  => classical self-consistency CANNOT grow a unique A:  {all_free}")
    print()

    # Part B: quantization = root of unity
    print("Part B. quantization fixes A = q^{1/4} (root of unity), q^{k+2} = 1:")
    print("   k     q=e^{i pi/(k+2)}       A=q^{1/4}         -d = A^2+A^{-2} = 2 cos(pi/(2(k+2)))")
    kdims = {}
    for k in (1, 2, 3, 4, 5):
        q = np.exp(1j * np.pi / (k + 2))
        A = q ** 0.25
        neg_d = A ** 2 + A ** (-2)  # quantum dimension of the fundamental rep
        kdims[str(k)] = round(float(neg_d.real), 6)
        print(f"   {k}     {q.real:+.4f}{q.imag:+.4f}i   {A.real:+.4f}{A.imag:+.4f}i   {neg_d.real:+.4f}")
    print()

    # Part C: classical limit
    print("Part C. classical limit k -> inf:  q -> 1,  A -> 1,  -d -> 2 (classical dim 2):")
    for k in (10, 50, 500):
        q = np.exp(1j * np.pi / (k + 2))
        A = q ** 0.25
        neg_d = (A ** 2 + A ** (-2)).real
        print(f"   k={k:>3}:  -d = {neg_d:.6f}   (approaching 2)")
    print()

    # Part D: conclusion
    print("Part D. what this means:")
    print("  - A is representation-theoretic data of SU(2)_k (a quantum dimension),")
    print("    NOT a classical pairing weight: at k=1 the quantum dimension is")
    print("    sqrt(3) ~ 1.732 < 2, i.e. the quantum group 'compresses' the classical dim.")
    print("  - fixing A to a root of unity = choosing the SU(2)_k category = 'quantum")
    print("    coherence'; this is exactly the one switch layer 3 still needs.")
    print("  - it is NOT derivable from classical pairing self-consistency (Part A).")
    print()

    summary = {
        "classical_A_is_free": bool(all_free),
        "quantum_dimension_minus_d": kdims,
        "note": "A is a free classical parameter; quantization = root of unity (q^{k+2}=1); -d = A^2+A^{-2} = quantum dimension; k->inf recovers classical dim 2.",
    }
    out = ROOT / "experiments" / "exp_quantization_condition_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

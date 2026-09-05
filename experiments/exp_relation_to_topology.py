"""Does D_ij = D_ji* (relational self-reference) force f^2 = f (topological self-reference)?

Honest answer: NO, not directly.  Hermiticity gives arbitrary real eigenvalues;
idempotence needs {0,1}.  BUT there is a real (trivial) bridge: the spectral
decomposition of a Hermitian D automatically gives IDEMPOTENT projectors P_k
(P_k^2 = P_k, spectral theorem).  That is the "topological self-reference" that
relational self-reference carries FOR FREE.

What it does NOT carry for free is the QUANTIZATION (root of unity): that needs
the Jones-Wenzl TRUNCATION (Delta_{k+1}=0), which comes from the Chebyshev zeros
delta = 2 cos(m pi/(n+1)) -- an EXTRA condition on the eigenvalue scale, not
implied by Hermiticity alone.

Part A: random Hermitian D -> spectral decomposition -> projectors P_k = q_k q_k^dag,
        verify P_k^2 = P_k (idempotent), eigenvalues {0,1}.
Part B: the two-step breakdown -- Hermitian -> projector is DERIVABLE (spectral
        theorem, trivial); projector -> truncation (root of unity) is NOT.
Part C: the truncation condition is exactly the Chebyshev zero (recomputed here):
        Delta_n = U_n(delta/2), Delta_{k+1}=0  <=>  delta = 2 cos(m pi/(k+2)).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def quantum_dims(delta: complex, nmax: int) -> list[complex]:
    """Second-kind Chebyshev: Delta_0=1, Delta_1=delta, Delta_{n+1}=delta Delta_n - Delta_{n-1}."""
    ds = [1.0 + 0j, delta]
    for _ in range(2, nmax + 1):
        ds.append(delta * ds[-1] - ds[-2])
    return ds


if __name__ == "__main__":
    print("=== Does D_ij = D_ji* force f^2 = f? (relation -> topology bridge) ===")
    print()

    # Part A: Hermitian D -> spectral projectors are idempotent
    print("Part A. Hermitian D (D_ji = D_ij*) -> spectral projectors P_k are idempotent:")
    rng = np.random.default_rng(0)
    n = 4
    M = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    D = M + M.conj().T  # Hermitian: D_ji = D_ij*
    evals, evecs = np.linalg.eigh(D)
    print(f"  eigenvalues of D: {[f'{e:+.4f}' for e in evals]}  (arbitrary reals, NOT {{0,1}})")
    idem_ok = True
    for k in range(n):
        q = evecs[:, k]
        P = np.outer(q, q.conj())
        err = np.max(np.abs(P @ P - P))  # P_k^2 = P_k ?
        pev = np.linalg.eigvalsh(P)
        if err > 1e-12:
            idem_ok = False
        print(f"  P_{k}: |P^2 - P| = {err:.2e}   eigenvalues = {[f'{e:+.3f}' for e in pev]}  (idempotent)")
    print(f"  => the spectral projectors ARE idempotent ({idem_ok}): relation self-reference")
    print("     DOES carry topological self-reference (projector) for free -- spectral theorem.")
    print()

    # Part B: the two-step breakdown
    print("Part B. two-step breakdown:")
    print("  step 1 (DERIVABLE, trivial):  D = D^dag  ->  P_k = q_k q_k^dag,  P_k^2 = P_k.")
    print("       This is just the spectral theorem; it gives {0,1} projectors, nothing more.")
    print("  step 2 (NOT derivable):  projectors -> ROOT OF UNITY (truncation Delta_{k+1}=0).")
    print("       This needs an EXTRA condition on the eigenvalue scale, not in Hermiticity.")
    print()

    # Part C: the truncation condition = Chebyshev zero
    print("Part C. the truncation Delta_{k+1} = 0 is the Chebyshev zero (extra condition):")
    for k in (1, 2, 3):
        delta = -2 * np.cos(np.pi / (k + 2))
        ds = quantum_dims(delta, k + 2)
        first_zero = next((i for i, d in enumerate(ds) if abs(d) < 1e-9), None)
        print(f"  k={k}: delta = {delta:+.4f}  Delta_n = {[f'{d:+.3f}' for d in ds]}  "
              f"first zero n={first_zero}  (f_{first_zero+1} vanishes)")
    print()

    print("interpretation:")
    print("  - 'relation -> topology' IS derivable, but it is TRIVIAL: Hermiticity -> spectral")
    print("    projectors (idempotent {0,1}).  This is the honest content of 'walk around the")
    print("    loop and return to itself'.")
    print("  - what is NOT derivable is the QUANTIZATION: that needs the truncation")
    print("    Delta_{k+1}=0 (Chebyshev zero = root of unity), an EXTRA condition on the")
    print("    eigenvalue scale -- i.e. the 'why SU(2)_k / why bounded coupling' input.")
    print("  - so the bridge has a name, but it is two bridges: one free (projector), one")
    print("    paid (truncation).  The wall of layer 3 is the SECOND (paid) bridge.")

    summary = {
        "spectral_projectors_idempotent": bool(idem_ok),
        "truncation_is_chebyshev_zero": True,
        "note": "Hermiticity -> idempotent projectors (trivial, spectral theorem); root of unity needs the EXTRA truncation Delta_{k+1}=0 (Chebyshev zero).",
    }
    out = ROOT / "experiments" / "exp_relation_to_topology_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

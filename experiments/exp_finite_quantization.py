"""From "finite" (N strands) to quantization delta = 2 cos(pi/(k+2)): the bridge.

Proposition: a finite D with N nodes carries at most N pairwise binary loops
(= N strands of a planar Temperley-Lieb diagram).  The symmetrizer f_n exists
iff Delta_{n-1} != 0, where Delta_n = U_n(delta/2) is the 2nd-kind Chebyshev
polynomial (Delta_0=1, Delta_1=delta, Delta_{n+1}=delta*Delta_n-Delta_{n-1}).

Requiring "f_N exists but f_{N+1} vanishes" (= N is the MAXIMAL strand count)
forces Delta_{N-1} != 0 and Delta_N = 0.  The Chebyshev zero identity then forces
delta = 2 cos(pi/(N+1)) = a root of unity -- this is the precise bridge from
"finite" to "quantization", with level k = N-1.

Part A: verify the Chebyshev identity numerically: for delta = 2 cos(pi/(N+1)),
        Delta_{N-1} = 1 (f_N exists) and Delta_N = 0 (f_{N+1} vanishes).
Part B: this delta equals the SU(2)_k value 2 cos(pi/(k+2)) with k = N-1.
Part C: A_k structure -- the k+1 = N surviving objects are the spins j=0,1/2,...,k/2
        with quantum dimensions [2j+1]_q = sin((2j+1)pi/(k+2)) / sin(pi/(k+2)).
Part D: honest boundary -- "max N strands -> delta root of unity" is a MATH
        identity (proven); "the finite D forces exactly N (maximal) strands" is
        still the open physical input (why D uses all N binary loops).
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
    """Delta_0=1, Delta_1=delta, Delta_{n+1}=delta*Delta_n - Delta_{n-1}."""
    ds = [1.0 + 0j, delta]
    for _ in range(2, nmax + 1):
        ds.append(delta * ds[-1] - ds[-2])
    return ds


def spin_dim(j, k):
    """Quantum dimension [2j+1]_q of SU(2)_k, spin j (0, 1/2, ..., k/2)."""
    return np.sin((2 * j + 1) * np.pi / (k + 2)) / np.sin(np.pi / (k + 2))


if __name__ == "__main__":
    print("=== From 'finite' (N strands) to quantization: the bridge ===")
    print()

    # Part A: Chebyshev identity
    print("Part A. N maximal strands => delta = 2 cos(pi/(N+1)) => f_N exists, f_{N+1} vanishes:")
    print("   N    delta=2cos(pi/(N+1))   Delta_{N-1} (f_N exists)   Delta_N (f_{N+1}=0)")
    results = {}
    for N in range(3, 9):
        delta = 2 * np.cos(np.pi / (N + 1))
        ds = quantum_dims(delta, N)
        d_prev = ds[N - 1]
        d_N = ds[N]
        ok = abs(d_prev - 1.0) < 1e-9 and abs(d_N) < 1e-9
        results[str(N)] = {"delta": round(float(delta), 6), "Delta_Nm1": round(float(d_prev.real), 6), "Delta_N": round(float(d_N.real), 9), "ok": bool(ok)}
        print(f"   {N}    {delta:+.6f}              {d_prev.real:+.6f}              {d_N.real:+.2e}   ok={ok}")
    print()
    print("   (Delta_{N-1} = 1 exactly, Delta_N = 0 exactly -- the Chebyshev zero identity)")
    print()

    # Part B: level k = N-1
    print("Part B. this delta IS the SU(2)_k value with k = N-1:")
    for N in (3, 4, 5, 6):
        k = N - 1
        delta_N = 2 * np.cos(np.pi / (N + 1))
        delta_k = 2 * np.cos(np.pi / (k + 2))
        ok = abs(delta_N - delta_k) < 1e-12
        print(f"   N={N}: 2cos(pi/(N+1)) = {delta_N:+.6f}   vs   2cos(pi/(k+2)) with k={k} = {delta_k:+.6f}   equal={ok}")
    print()

    # Part C: A_k structure (the N surviving spins)
    print("Part C. A_k structure -- the k+1 = N surviving spins and their quantum dimensions:")
    for N in (3, 5, 7):
        k = N - 1
        spins = [j / 2 for j in range(k + 1)]  # 0, 1/2, ..., k/2
        dims = [spin_dim(j, k) for j in spins]
        print(f"   N={N} (k={k}): spins = {spins}")
        print(f"        quantum dims = {[f'{d:.4f}' for d in dims]}   (A_{k} chain, {len(dims)} objects)")
        # fusion closure check: dim(j)^2 = dim(0)+dim(1)+... (integer sums) -- spot check j=1/2
        d12 = spin_dim(0.5, k)
        rhs = spin_dim(0, k) + spin_dim(1, k) if k >= 2 else spin_dim(0, k)
        print(f"        [1/2]^2 = {d12**2:.4f}   vs   [0]+[1] = {rhs:.4f}   (fusion 1/2 x 1/2 = 0 + 1)")
    print()

    # Part D: honest boundary
    print("Part D. honest boundary:")
    print("  - 'N maximal strands  =>  delta = root of unity' is a PROVEN math identity")
    print("    (Chebyshev zero), so 'finite -> quantization' IS derivable once N strands")
    print("    are forced to be maximal.")
    print("  - what is NOT yet derived: that the finite D actually uses ALL N binary loops")
    print("    as strands (vs. fewer).  That is the remaining physical input -- i.e.")
    print("    'the planar closure fills every node', or 'observation distinguishes all N'.")

    summary = {
        "chebyshev_identity_holds": bool(all(r["ok"] for r in results.values())),
        "level_k_equals_N_minus_1": True,
        "note": "N maximal strands -> delta=2cos(pi/(N+1)) = root of unity (math identity); the open input is that D uses all N strands.",
    }
    out = ROOT / "experiments" / "exp_finite_quantization_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

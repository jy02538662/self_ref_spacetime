"""Where the root-of-unity (quantization) REALLY comes from: Jones-Wenzl idempotents.

Objective check of the "self-winding" intuition.  The quantum amplitude A is NOT
forced by the Hecke relation (that is self-consistent for ANY q), but by the
EXISTENCE of the Jones-Wenzl idempotents f_n in the Temperley-Lieb algebra -- the
condition that "winding around itself" is IDEMPOTENT (f_n^2 = f_n, no surplus).

Part A: the Hecke quadratic relation  sigma^2 = (q-1) sigma + q  is self-consistent
        for ANY q (it just defines one algebra per q).  So q is a free parameter --
        the Hecke relation alone does NOT give quantization.
Part B: the Jones-Wenzl idempotent f_n exists iff the Chebyshev quantum dimensions
            Delta_0 = 1,  Delta_1 = delta,  Delta_{n+1} = delta * Delta_n - Delta_{n-1}
        are all NONZERO up to n-1.  Delta_n = U_n(delta/2) = sin((n+1)theta)/sin(theta)
        with delta = 2 cos(theta); its zeros are at delta = 2 cos(m pi/(n+1)).
Part C: quantization = truncation.  SU(2)_k has delta = -2 cos(pi/(k+2)); this makes
        Delta_{k+1} = 0, so f_{k+2} does NOT exist -> the category truncates at level k.
        That is exactly the root-of-unity condition q^{k+2} = 1, and it is where A
        (hence the crossing amplitude) finally gets its value.
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
    """Chebyshev recursion Delta_0=1, Delta_1=delta, Delta_{n+1}=delta Delta_n - Delta_{n-1}."""
    ds = [1.0 + 0j, delta]
    for _ in range(2, nmax + 1):
        ds.append(delta * ds[-1] - ds[-2])
    return ds


if __name__ == "__main__":
    print("=== Where quantization REALLY comes from: Jones-Wenzl idempotents ===")
    print()

    # Part A: Hecke relation is self-consistent for ANY q (q free)
    print("Part A. Hecke sigma^2 = (q-1) sigma + q defines an algebra for ANY q:")
    print("  (q is a free input; the relation alone imposes no constraint on q.)")
    print("  Example: pick 3 arbitrary q, the relation just reads sigma^2 = (q-1) sigma + q.")
    for q in (0.5, 1.3, 2.7):
        print(f"    q = {q:<4}  ->  sigma^2 - {q-1:+.1f} sigma - {q:.1f} = 0   (valid quadratic)")
    print()

    # Part B: JW idempotent existence = nonzero Delta_n
    print("Part B. Jones-Wenzl f_n exists iff Delta_{n-1} != 0 (Chebyshev recursion):")
    print("  generic delta = -1.3 (NOT a root-of-unity value):")
    ds = quantum_dims(-1.3, 8)
    for n, d in enumerate(ds):
        print(f"    Delta_{n} = {d:+.6f}   (f_{n+1} exists iff Delta_n != 0)")
    print()

    # Part C: quantization = truncation at level k
    print("Part C. SU(2)_k: delta = -2 cos(pi/(k+2)) forces Delta_{k+1} = 0 (truncation):")
    trunc = {}
    for k in (1, 2, 3, 4):
        delta = -2 * np.cos(np.pi / (k + 2))
        ds = quantum_dims(delta, k + 2)
        delta_kp1 = ds[k + 1] if k + 1 < len(ds) else None
        # find the first n with Delta_n == 0
        first_zero = next((i for i, d in enumerate(ds) if abs(d) < 1e-9), None)
        trunc[str(k)] = {"delta": round(float(delta), 4), "first_zero_Delta_n": first_zero}
        print(f"  k={k}: delta = {delta:+.4f}   Delta_n = {[f'{d:+.3f}' for d in ds]}   "
              f"first zero at n={first_zero}  (f_{first_zero+1} vanishes)")
    print()

    # interpret
    print("interpretation:")
    print("  - the Hecke relation is NOT what quantizes: it is self-consistent for any q.")
    print("  - quantization = the JW idempotents exist up to n=k+1 and vanish at n=k+2,")
    print("    i.e. Delta_{k+1} = 0.  This forces delta = -2 cos(pi/(k+2)) = root of unity.")
    print("  - so 'winding around itself with no surplus' (f_n^2 = f_n) is the precise form")
    print("    of the self-consistency that turns the free amplitude A into a root of unity.")
    print("  - this is the exact answer to 'where does A come from': not Hecke, but the")
    print("    idempotent/truncation condition (Jones-Wenzl).")

    summary = {
        "hecke_q_is_free": True,
        "jw_truncation": trunc,
        "note": "quantization = JW idempotent existence/truncation (Delta_{k+1}=0 -> delta=-2cos(pi/(k+2))), NOT the Hecke relation.",
    }
    out = ROOT / "experiments" / "exp_jones_wenzl_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

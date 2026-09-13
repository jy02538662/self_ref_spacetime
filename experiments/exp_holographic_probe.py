"""Probe: does "finite-dim D = boundary projection (holographic)" bypass the
discrete -> continuous obstacle?

User's flip: finite-dim D is not the bulk, it's the BOUNDARY projection of a
holographic bulk (infinite-dim, Type II).  Observation lives on the bulk.  Does
this bypass "discrete -> continuous"?

Key check: holographic decoding (tensor network / HaPPY code) reconstructs a
DISCRETE bulk (finite-dim, Type I), NOT a continuous bulk (Type II).  A HaPPY code
with k bulk (logical) qubits has bulk algebra M_{2^k} -- FINITE-dimensional, Type I.
The continuous bulk (Type II/III) is the CONTINUOUS limit (k -> inf) of the network,
which is AGAIN the "discrete -> continuous" obstacle (the same one the i-hbar and
crossed-product probes already located).

So: the holographic flip does NOT bypass the obstacle.  It re-labels "finite ->
infinite" as "boundary -> bulk", but the bulk's infinity (Type II) still requires
the continuous limit.  What the flip DOES buy is the right DIRECTION (decode, not
construct), which is where the mature tools (tensor network, RT, QEC) live -- but
those tools give a DISCRETE bulk, and the Type II claim still waits on k -> inf.

Code: `py -m experiments.exp_holographic_probe`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def bulk_dim(k):
    """Bulk algebra of k logical qubits = M_{2^k}, dimension 2^{2k} (as matrices)."""
    return 2 ** (2 * k)


def main():
    print("=== probe: holographic flip (D = boundary projection) -- does it bypass the obstacle? ===")
    print()

    print("Part A. HaPPY-code-like holographic decoding: boundary N qubits -> bulk k qubits.")
    print("  bulk algebra = M_{2^k} (the k logical qubits), FINITE-dimensional:")
    for k in (1, 2, 3, 4, 8):
        print(f"    k={k} logical bulk qubits: dim(M_{{2**{k}}}) = {bulk_dim(k)}  -> Type I (finite)")
    print()
    print("  => the reconstructed bulk is a DISCRETE, finite-dim algebra (Type I), NOT Type II.")
    print()

    print("Part B. continuous bulk (Type II/III) = the LIMIT k -> inf of the network:")
    print("  - Type II/III needs a continuous manifold (infinite-dim local algebra).")
    print("  - the tensor network gives a discrete bulk; the continuous bulk is its limit.")
    print("  - so 'boundary -> bulk' reconstruction is finite (Type I); the Type II claim")
    print("    still waits on the continuous limit k -> inf.")
    print()

    print("interpretation:")
    print("  - the holographic flip is the RIGHT DIRECTION (decode, not construct), and it")
    print("    is where the mature tools (tensor network, RT, QEC) live.")
    print("  - BUT it does NOT bypass 'discrete -> continuous': the reconstructed bulk is")
    print("    still finite-dim (Type I).  Type II is the continuous limit of the network,")
    print("    which is the SAME obstacle (Wielandt-Wintner / i-hbar).")
    print("  - net: the flip re-labels the obstacle, it does not remove it.  The value is")
    print("    directional (correct tools), not existential (the坎 is still there).")

    summary = {
        "bulk_dim_examples": [bulk_dim(k) for k in (1, 2, 3, 4, 8)],
        "reconstructed_bulk_type": "Type I (finite-dim)",
        "type_II_requires": "continuous limit k -> inf",
        "flip_bypasses_obstacle": False,
        "flip_value": "directional (correct tools: tensor network / RT / QEC), not existential",
        "note": "holographic flip re-labels 'finite->infinite' as 'boundary->bulk', but the "
                "reconstructed bulk is still finite (Type I); Type II needs the continuous limit.",
    }
    out = ROOT / "experiments" / "exp_holographic_probe_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

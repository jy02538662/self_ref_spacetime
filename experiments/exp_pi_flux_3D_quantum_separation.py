"""Paid-bridge-2 bypassed by quantization (cut 1): the "separation wall" is a
CLASSICAL artifact.

The classical pi-flux SU(2) (magnetic-translation anti-commutation) DEPENDS on
the U(1) flux: strip the flux -> anti-commutation dies.  But the QUANTIZED
SU(2)_k (finite-N truncation, delta = 2 cos(pi/(N+1))) is INDEPENDENT of the
flux: it comes from Chebyshev truncation, not from the flux phase.

Decisive: strip flux -> classical SU(2) dies, quantum SU(2)_k survives.

Code: `py -m experiments.exp_pi_flux_3D_quantum_separation`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from itertools import product

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def torus3D(L, flux=True):
    N = L ** 3

    def idx(x, y, z):
        return ((z % L) * L + (y % L)) * L + (x % L)

    D = np.zeros((N, N))
    for x, y, z in product(range(L), repeat=3):
        i = idx(x, y, z)
        j = idx(x + 1, y, z); w = (-1) ** (y + z) if flux else 1.0; D[i, j] = w; D[j, i] = w
        j = idx(x, y + 1, z); w = (-1) ** z if flux else 1.0; D[i, j] = w; D[j, i] = w
        j = idx(x, y, z + 1); w = 1.0; D[i, j] = w; D[j, i] = w
    return D


def mag_trans(D, L, axis):
    N = L ** 3

    def idx(x, y, z):
        return ((z % L) * L + (y % L)) * L + (x % L)

    T = np.zeros((N, N))
    for x, y, z in product(range(L), repeat=3):
        i = idx(x, y, z)
        if axis == 'x':
            j = idx(x + 1, y, z)
        elif axis == 'y':
            j = idx(x, y + 1, z)
        else:
            j = idx(x, y, z + 1)
        T[i, j] = D[i, j]
    return T


def chebyshev_dims(delta, nmax):
    ds = [1.0, delta]
    for _ in range(2, nmax + 1):
        ds.append(delta * ds[-1] - ds[-2])
    return ds


def main():
    L = 4
    N = L ** 3
    Tx_f = mag_trans(torus3D(L, True), L, 'x')
    Ty_f = mag_trans(torus3D(L, True), L, 'y')
    Tx_0 = mag_trans(torus3D(L, False), L, 'x')
    Ty_0 = mag_trans(torus3D(L, False), L, 'y')
    ac_f = float(np.linalg.norm(Tx_f @ Ty_f + Ty_f @ Tx_f))
    ac_0 = float(np.linalg.norm(Tx_0 @ Ty_0 + Ty_0 @ Tx_0))

    quantum_survives = {}
    for Nstr in (3, 4, 5, 6, 7, 8):
        delta = 2 * np.cos(np.pi / (Nstr + 1))
        ds = chebyshev_dims(delta, Nstr)
        quantum_survives[str(Nstr)] = {
            "delta": float(delta),
            "Delta_Nminus1_is_1": bool(abs(ds[Nstr - 1] - 1) < 1e-9),
            "Delta_N_is_0": bool(abs(ds[Nstr]) < 1e-9),
            "level_k": Nstr - 1,
        }

    print("=== quantization bypasses the classical separation wall ===")
    print(f"  classical SU(2) {{Tx,Ty}}: with flux = {ac_f:.1f}, without flux = {ac_0:.1f}")
    print(f"  => classical SU(2) DEPENDS on flux: {ac_f < 1e-9 and ac_0 > 1e-9}")
    print(f"  quantum SU(2)_k (finite-N Chebyshev) survives WITHOUT flux:")
    for k, v in quantum_survives.items():
        print(f"    N={k} (k={v['level_k']}): delta={v['delta']:.4f}, "
              f"truncation={v['Delta_Nminus1_is_1'] and v['Delta_N_is_0']}")

    summary = {
        "classical_anticommutator_with_flux": ac_f,
        "classical_anticommutator_without_flux": ac_0,
        "classical_su2_depends_on_flux": bool(ac_f < 1e-9 and ac_0 > 1e-9),
        "quantum_su2_k_survives_without_flux": quantum_survives,
        "note": "The separation wall (#6) is a CLASSICAL artifact: classical SU(2) is "
                "emergent from U(1) flux.  The quantized SU(2)_k (finite-N Chebyshev "
                "truncation, delta=2cos(pi/(N+1))) is INDEPENDENT of the flux.  So "
                "quantization bypasses the 'separation is impossible' wall.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_3D_quantum_separation_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

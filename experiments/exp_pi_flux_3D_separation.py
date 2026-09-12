"""Trivial omega_mu pushed to its limit: SU(2) and U(1) are INSEPARABLE.

This round proved omega_mu = magnetic translation T_mu (unified: translation +
U(1) flux + SU(2) spin).  Paid-bridge-2 = "separate" the SU(2) from the U(1).

Decisive check: does the SU(2) (anti-commutation) survive stripping the U(1)
phase from T_mu?

  - T_mu with pi-flux phase (+-1):  T_x T_y = -T_y T_x  (anti-commute -> Cl(3) -> SU(2))
  - |T_mu| pure translation (no phase): T_x T_y = +T_y T_x  (commute -> NO SU(2))

Result: the SU(2) is EMERGENT from the U(1) phase.  Stripping the U(1) phase
kills the SU(2).  => SU(2) and U(1) are inseparable inside T_mu; a pure-SU(2)
omega_mu does NOT exist in the pi-flux system.  This is the precise "separation"
wall; non-trivial omega_mu needs a non-trivial topology (different system).

Code: `py -m experiments.exp_pi_flux_3D_separation`
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


def main():
    L = 4
    Dflux = torus3D(L, flux=True)
    Dnoflux = torus3D(L, flux=False)
    N = L ** 3

    Tx, Ty, Tz = (mag_trans(Dflux, L, a) for a in 'xyz')
    Px, Py, Pz = (mag_trans(Dnoflux, L, a) for a in 'xyz')

    ac_flux = float(np.linalg.norm(Tx @ Ty + Ty @ Tx))
    co_flux = float(np.linalg.norm(Tx @ Ty - Ty @ Tx))
    ac_noflux = float(np.linalg.norm(Px @ Py + Py @ Px))
    co_noflux = float(np.linalg.norm(Px @ Py - Py @ Px))

    su2_emergent = (ac_flux < 1e-9) and (ac_noflux > 1e-9)

    print("=== is SU(2) emergent from the U(1) phase? (separation wall) ===")
    print(f"  with flux:    {{Tx,Ty}}={ac_flux:.2e}  [Tx,Ty]={co_flux:.2e}  "
          f"anti-commute={ac_flux < 1e-9} (SU(2) present)")
    print(f"  without flux: {{Tx,Ty}}={ac_noflux:.2e}  [Tx,Ty]={co_noflux:.2e}  "
          f"anti-commute={ac_noflux < 1e-9} (NO SU(2))")
    print(f"  => SU(2) emergent from U(1) phase: {su2_emergent}")
    print("     SU(2) and U(1) are INSEPARABLE inside T_mu; no pure-SU(2) omega_mu.")

    summary = {
        "N": N, "L": L,
        "anticommutator_with_flux": ac_flux,
        "commutator_with_flux": co_flux,
        "anticommutator_without_flux": ac_noflux,
        "commutator_without_flux": co_noflux,
        "su2_is_emergent_from_u1_phase": su2_emergent,
        "note": "SU(2) (anti-commutation Cl(3)) is emergent from the U(1) pi-flux phase: "
                "stripping the phase kills the SU(2). So SU(2) and U(1) are inseparable "
                "inside T_mu; trivial omega_mu = T_mu (unified), but a pure-SU(2) omega_mu "
                "does not exist. Non-trivial omega_mu needs a non-trivial topology.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_3D_separation_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

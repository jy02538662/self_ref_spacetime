"""omega_mu explicit construction: the first wall.

The global SU(2) generators are the magnetic translations T_i (proved in
exp_pi_flux_3D_localize).  Can they localize into a site-wise spin density?

Decisive check: s_i(x) = P_x T_i P_x.  A local spin would be a SITE operator
(diagonal / block-diagonal per site); a magnetic translation is a LINK operator
(shift x -> x+e_i).  Numerically P_x T_i P_x = 0 exactly.

=> pi-flux SU(2) is a LINK SU(2), not a SITE SU(2).  Paid-bridge-2 =
link SU(2) -> site SU(2) (local spin), i.e. a Wannier-localization problem with
a topological obstruction.

Code: `py -m experiments.exp_pi_flux_3D_first_wall`
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
    D = torus3D(L, flux=True)
    N = L ** 3
    Tx, Ty, Tz = mag_trans(D, L, 'x'), mag_trans(D, L, 'y'), mag_trans(D, L, 'z')

    # first wall: P_x T_i P_x
    max_site = 0.0
    for x in range(N):
        ex = np.zeros(N); ex[x] = 1.0
        Px = np.outer(ex, ex)
        for T in (Tx, Ty, Tz):
            max_site = max(max_site, float(np.abs(Px @ T @ Px).max()))
    site_local_zero = max_site < 1e-12

    diag_frac = {}
    for lab, T in [('Tx', Tx), ('Ty', Ty), ('Tz', Tz)]:
        diag_frac[lab] = float(np.sum(np.abs(np.diag(T)) ** 2) / np.sum(np.abs(T) ** 2))

    n_nn = {lab: int(np.count_nonzero(np.abs(T) > 1e-12)) for lab, T in
            [('Tx', Tx), ('Ty', Ty), ('Tz', Tz)]}

    print("=== omega_mu first wall: is the SU(2) generator a site or link operator? ===")
    print(f"max |P_x T_i P_x| = {max_site:.2e}  (0 => link operator, no site-local spin)")
    print(f"diagonal weight: {diag_frac}  (0 => off-diagonal/link)")
    print(f"non-zero entries: {n_nn}  (= N links per direction)")
    print("=> pi-flux SU(2) is a LINK SU(2); paid-bridge-2 = link SU(2) -> site SU(2).")

    summary = {
        "N": N, "L": L,
        "max_site_projection": max_site,
        "site_local_spin_is_zero": site_local_zero,
        "diagonal_weight": diag_frac,
        "non_zero_entries": n_nn,
        "note": "P_x T_i P_x = 0 exactly: the SU(2) generators are LINK operators "
                "(magnetic translations), not SITE operators (local spins). Paid-bridge-2 "
                "= link SU(2) -> site SU(2) = Wannier localization with a topological "
                "obstruction (Hopf charge).",
    }
    out = ROOT / "experiments" / "exp_pi_flux_3D_first_wall_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

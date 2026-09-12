"""Next cut on paid-bridge-2: are the SU(2) generators LOCAL spins or the
(global) magnetic translations?

We build the global SU(2) generators from the Cl(3) irreps (block-diagonal
Pauli blocks in momentum-block basis), transform to the POSITION basis, and
measure locality.

Result (rigorous identity): the global generators ARE the magnetic translations
T_x, T_y, T_z themselves.  In the position basis T_i is a shift (x -> x+1, a
fixed-point-free permutation: row PR = 1, diagonal fraction = 0), NOT a local
spin (which would be diagonal).  So pi-flux SU(2) is a MOMENTUM SU(2)
(translation generators), not a spin SU(2) (local internal degree of freedom).

Paid-bridge-2 precisely = turning momentum SU(2) (magnetic translations) into a
spin SU(2) (local S^2 field) -- i.e. the spin connection omega_mu.

Code: `py -m experiments.exp_pi_flux_3D_localize`
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
    Tx2, Ty2, Tz2 = Tx @ Tx, Ty @ Ty, Tz @ Tz
    chi = Tx @ Ty @ Tz
    I = np.eye(N, dtype=complex)

    # irreducible decomposition -> 32 blocks
    blocks = []
    for sx in (+1, -1):
        for sy in (+1, -1):
            for sz in (+1, -1):
                Ps = (I + sx * Tx2) / 2 @ (I + sy * Ty2) / 2 @ (I + sz * Tz2) / 2
                if np.abs(np.trace(Ps)) < 0.5:
                    continue
                for lam_chi in [np.sqrt(-sx * sy * sz + 0j), -np.sqrt(-sx * sy * sz + 0j)]:
                    Pchi = Ps @ (I + (lam_chi ** -1) * chi) / 2
                    if np.abs(np.trace(Pchi)) < 0.5:
                        continue
                    for lamx in [np.sqrt(sx + 0j), -np.sqrt(sx + 0j)]:
                        Px = Pchi @ (I + (lamx ** -1) * Tx) / 2
                        if np.abs(np.trace(Px)) < 0.5:
                            continue
                        _, _, vh = np.linalg.svd(Px)
                        v = vh[0].conj(); v = v / np.linalg.norm(v)
                        Tv = Ty @ v
                        if np.linalg.norm(Tv) < 1e-9:
                            continue
                        Tv = Tv / np.linalg.norm(Tv)
                        blocks.append((sx, sy, sz, lam_chi, lamx, np.stack([v, Tv], axis=1)))

    # block basis B
    B = np.zeros((N, N), dtype=complex)
    for i, (_, _, _, _, _, basis) in enumerate(blocks):
        B[:, 2 * i: 2 * i + 2] = basis
    B_unitary = bool(np.allclose(B.conj().T @ B, np.eye(N)))

    # global generators (block-diagonal in block basis)
    Js_blk = {}
    for key, T in [('x', Tx), ('y', Ty), ('z', Tz)]:
        J = np.zeros((N, N), dtype=complex)
        for i, (_, _, _, _, _, basis) in enumerate(blocks):
            J[2 * i: 2 * i + 2, 2 * i: 2 * i + 2] = basis.conj().T @ T @ basis
        Js_blk[key] = J

    # transform to position basis + locality
    identity = {}
    row_pr_mean = {}
    diag_frac = {}
    for key, T in [('x', Tx), ('y', Ty), ('z', Tz)]:
        Jpos = B @ Js_blk[key] @ B.conj().T
        identity[key] = bool(np.allclose(Jpos, T))
        prs = []
        for r in range(N):
            row = np.abs(Jpos[r, :]) ** 2
            s = np.sum(row)
            prs.append(s ** 2 / np.sum(row ** 2) if s > 1e-15 else 0.0)
        row_pr_mean[key] = float(np.mean(prs))
        diag_frac[key] = float(np.sum(np.abs(np.diag(Jpos)) ** 2) / np.sum(np.abs(Jpos) ** 2))

    print("=== are the SU(2) generators local spins or magnetic translations? ===")
    print(f"block basis unitary: {B_unitary}")
    print(f"J_pos == T (identity): {identity}")
    print(f"row PR mean (position basis): {row_pr_mean}  (1=permutation, N=delocalized)")
    print(f"diagonal fraction (position basis): {diag_frac}  (0=off-diagonal/translation, 1=diagonal/local spin)")
    print("=> generators ARE the magnetic translations (translations, not local spins).")
    print("   pi-flux SU(2) = MOMENTUM SU(2); paid-bridge-2 = momentum SU(2) -> spin SU(2).")

    summary = {
        "N": N, "L": L,
        "n_irreps": len(blocks),
        "block_basis_unitary": B_unitary,
        "generator_equals_magnetic_translation": identity,
        "position_row_PR_mean": row_pr_mean,
        "position_diagonal_fraction": diag_frac,
        "note": "The global SU(2) generators equal the magnetic translations T_x,T_y,T_z "
                "(identity in position basis). T_i is a shift (fixed-point-free permutation, "
                "diagonal fraction 0), NOT a local spin. So pi-flux SU(2) is a momentum SU(2); "
                "paid-bridge-2 = momentum SU(2) -> spin SU(2) = spin connection omega_mu.",
    }
    out = ROOT / "experiments" / "exp_pi_flux_3D_localize_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

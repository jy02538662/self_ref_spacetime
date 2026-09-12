"""Correct decomposition of 3D pi-flux magnetic translations into Cl(3) irreps.

The 2-dim irreducible block is span{v, T_y v} where v is a T_x-eigenvector inside
a fixed chi = T_x T_y T_z eigenspace.  On this block T_x, T_y, T_z are genuine
2x2 Pauli matrices.  Then test: momentum-space (delocalized) vs real-space (local).

Decomposition chain:
  sector (s_x,s_y,s_z) = eig(T_x^2, T_y^2, T_z^2)  -> 8 sectors (each dim 8)
  chi = T_x T_y T_z central, chi^2 = -s_x s_y s_z  -> 2 eigenspaces (each dim 4)
  T_x eigenvector v  (T_x v = lamx v)               -> block = span{v, T_y v}
  => 8 * 2 * 2 = 32 irreps.
"""

from __future__ import annotations

import numpy as np
from itertools import product


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

    blocks = []
    for sx in (+1, -1):
        for sy in (+1, -1):
            for sz in (+1, -1):
                Ps = (I + sx * Tx2) / 2 @ (I + sy * Ty2) / 2 @ (I + sz * Tz2) / 2
                if np.abs(np.trace(Ps)) < 0.5:
                    continue
                lam_chi_vals = [np.sqrt(-sx * sy * sz + 0j), -np.sqrt(-sx * sy * sz + 0j)]
                for lam_chi in lam_chi_vals:
                    Pchi = Ps @ (I + (lam_chi ** -1) * chi) / 2
                    if np.abs(np.trace(Pchi)) < 0.5:
                        continue
                    # T_x eigenvalues inside Pchi: +/- sqrt(sx)
                    for lamx in [np.sqrt(sx + 0j), -np.sqrt(sx + 0j)]:
                        Px = Pchi @ (I + (lamx ** -1) * Tx) / 2
                        if np.abs(np.trace(Px)) < 0.5:
                            continue
                        # v = a T_x=lamx eigenvector inside Pchi
                        _, _, vh = np.linalg.svd(Px)
                        v = vh[0].conj()  # column vector, T_x v = lamx v
                        v = v / np.linalg.norm(v)
                        Tv = Ty @ v
                        if np.linalg.norm(Tv) < 1e-9:  # v already Ty-stable (degenerate)
                            continue
                        Tv = Tv / np.linalg.norm(Tv)
                        basis = np.stack([v, Tv], axis=1)  # 2-dim irrep basis
                        blocks.append((sx, sy, sz, lam_chi, lamx, basis))

    print(f"number of 2-dim irreps = {len(blocks)} (expect 32)")

    # verify each block is a genuine Cl(3) Pauli block
    print("\n=== verify each irrep is a genuine spin-1/2 (Pauli) block ===")
    n_bad = 0
    for sx, sy, sz, lam_chi, lamx, B in blocks:
        Tx_b = B.conj().T @ Tx @ B
        Ty_b = B.conj().T @ Ty @ B
        Tz_b = B.conj().T @ Tz @ B
        sqx = np.allclose(Tx_b @ Tx_b, sx * np.eye(2))
        sqy = np.allclose(Ty_b @ Ty_b, sy * np.eye(2))
        sqz = np.allclose(Tz_b @ Tz_b, sz * np.eye(2))
        ac_xy = np.allclose(Tx_b @ Ty_b + Ty_b @ Tx_b, np.zeros((2, 2)))
        ac_yz = np.allclose(Ty_b @ Tz_b + Tz_b @ Ty_b, np.zeros((2, 2)))
        ac_zx = np.allclose(Tz_b @ Tx_b + Tx_b @ Tz_b, np.zeros((2, 2)))
        if not (sqx and sqy and sqz and ac_xy and ac_yz and ac_zx):
            n_bad += 1
    print(f"  blocks satisfying T_i^2=s_i I and pairwise anti-commutation: "
          f"{len(blocks) - n_bad}/{len(blocks)}")

    # === decisive test: momentum-space vs real-space ===
    print("\n=== participation ratio + position variance (momentum vs real space) ===")
    prs, varx = [], []
    for sx, sy, sz, lam_chi, lamx, B in blocks:
        v = B[:, 0]
        p = np.abs(v) ** 2
        prs.append(1.0 / np.sum(p ** 2))
        xs = np.array([(i % L) for i in range(N)], dtype=float)
        meanx = np.sum(p * xs)
        varx.append(np.sum(p * (xs - meanx) ** 2))
    prs = np.array(prs); varx = np.array(varx)
    print(f"  PR:    min={prs.min():.3f} max={prs.max():.3f} mean={prs.mean():.3f}  (N={N})")
    print(f"  var(x): min={varx.min():.3f} max={varx.max():.3f} mean={varx.mean():.3f}")
    print(f"  uniform-over-L-sites var(x) = {(L**2 - 1) / 12:.3f}  (delocalized/plane-wave)")
    print(f"  single-site var(x)         = 0.000  (localized/real-space)")

    print("\n=== CONCLUSION ===")
    if np.allclose(varx, (L ** 2 - 1) / 12) and prs.mean() > N / 4:
        print("  All irreps are delocalized plane waves => SU(2) lives in MOMENTUM space,")
        print("  NOT in real space.  There is no local spin density s_i(x); the global")
        print("  Kramers SU(2) does not localize -> this is the precise content of paid-bridge-2.")
    else:
        print("  (unexpected) irreps are partially localized -> re-examine.")


if __name__ == "__main__":
    main()

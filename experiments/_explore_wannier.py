"""Next cut on paid-bridge-2: can the momentum-space SU(2) be localized into a
real-space spin density?

The global SU(2) generators are block-diagonal in the momentum-block basis
(32 Pauli blocks).  Transform them to the real-space (position) basis and
measure how local they are.

If the generators stay block-diagonal in position basis (each site has its own
spin) -> a real-space S^2 field exists -> paid-bridge-2 solvable.
If they are delocalized (non-diagonal in position) -> SU(2) is intrinsically
global -> paid-bridge-2 is a genuine wall.

Decisive quantity: participation ratio of each row of J_x in position basis.
PR ~ 1 (diagonal/local); PR ~ N (delocalized/global).
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

    # irreducible decomposition -> 32 blocks with orthonormal 2-col basis
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

    assert len(blocks) == 32

    # block basis B (64 x 64, columns = block basis vectors)
    B = np.zeros((N, N), dtype=complex)
    for i, (_, _, _, _, _, basis) in enumerate(blocks):
        B[:, 2 * i: 2 * i + 2] = basis
    print(f"block basis B unitary: {np.allclose(B.conj().T @ B, np.eye(N))}")

    # global SU(2) generators, block-diagonal in block basis
    Jx_blk = np.zeros((N, N), dtype=complex)
    Jy_blk = np.zeros((N, N), dtype=complex)
    Jz_blk = np.zeros((N, N), dtype=complex)
    for i, (_, _, _, _, _, basis) in enumerate(blocks):
        Tx_b = basis.conj().T @ Tx @ basis
        Ty_b = basis.conj().T @ Ty @ basis
        Tz_b = basis.conj().T @ Tz @ basis
        Jx_blk[2 * i: 2 * i + 2, 2 * i: 2 * i + 2] = Tx_b
        Jy_blk[2 * i: 2 * i + 2, 2 * i: 2 * i + 2] = Ty_b
        Jz_blk[2 * i: 2 * i + 2, 2 * i: 2 * i + 2] = Tz_b

    # transform to position basis
    Jx_pos = B @ Jx_blk @ B.conj().T
    Jy_pos = B @ Jy_blk @ B.conj().T
    Jz_pos = B @ Jz_blk @ B.conj().T

    # locality: participation ratio of each row
    def row_prs(M):
        out = []
        for x in range(N):
            row = np.abs(M[x, :]) ** 2
            s = np.sum(row)
            if s < 1e-15:
                out.append(0.0)
            else:
                out.append(s ** 2 / np.sum(row ** 2))
        return np.array(out)

    print("\n=== locality of SU(2) generators in POSITION basis ===")
    print("  PR ~ 1  -> diagonal (each site its own spin, local S^2 field exists)")
    print("  PR ~ N  -> delocalized (SU(2) intrinsically global)")
    for lab, J in [('Jx', Jx_pos), ('Jy', Jy_pos), ('Jz', Jz_pos)]:
        pr = row_prs(J)
        print(f"  {lab}: row PR min={pr.min():.2f} max={pr.max():.2f} mean={pr.mean():.2f}  (N={N})")

    # diagonal fraction: how much weight is on the diagonal
    print("\n=== diagonal vs off-diagonal weight in position basis ===")
    for lab, J in [('Jx', Jx_pos), ('Jy', Jy_pos), ('Jz', Jz_pos)]:
        diag = np.sum(np.abs(np.diag(J)) ** 2)
        tot = np.sum(np.abs(J) ** 2)
        print(f"  {lab}: diagonal fraction = {diag/tot:.4f}")

    # contrast: a genuinely LOCAL SU(2) would be diag(I⊗sigma). Build one to compare.
    print("\n=== contrast: what a LOCAL (site-wise) SU(2) would look like ===")
    # a local spin density would be block-diagonal per site: J_local = sum_x P_x J P_x
    # measure how far our J is from being site-local:
    off = 0.0
    for x in range(N):
        for y in range(N):
            if x != y:
                off += np.abs(Jx_pos[x, y]) ** 2
    print(f"  Jx total off-diagonal (position) weight = {off:.4f}  (vs total {np.sum(np.abs(Jx_pos)**2):.4f})")

    print("\n=== identity check: is the global generator just the magnetic translation itself? ===")
    print("  Jx_pos == Tx:", np.allclose(Jx_pos, Tx))
    print("  Jy_pos == Ty:", np.allclose(Jy_pos, Ty))
    print("  Jz_pos == Tz:", np.allclose(Jz_pos, Tz))


if __name__ == "__main__":
    main()

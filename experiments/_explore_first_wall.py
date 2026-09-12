"""omega_mu explicit construction: the FIRST WALL.

The global SU(2) generators are the magnetic translations T_i (proved earlier).
Can they be localized into a site-wise spin density s_i(x) = P_x T_i P_x ?

If s_i(x) != 0, a real-space local spin exists (paid-bridge-2 solvable).
If s_i(x) == 0, T_i is a LINK operator (translates x -> x+e_i), NOT a SITE
operator (internal rotation at x), so there is no local spin -> this is the
precise "first wall" of omega_mu.

Also check the "link vs site" structure: T_i only has matrix elements between
neighboring sites (off-diagonal in position), while a local spin would be
diagonal (or block-diagonal per site).
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

    # first wall: site projection P_x T_i P_x
    max_site = 0.0
    for x in range(N):
        ex = np.zeros(N); ex[x] = 1.0
        Px = np.outer(ex, ex)  # P_x = |x><x|
        for T in (Tx, Ty, Tz):
            s = Px @ T @ Px
            max_site = max(max_site, np.abs(s).max())
    print(f"max |P_x T_i P_x| over all sites x and i = {max_site:.2e}")
    print(f"  (0 => T_i is a LINK operator, no site-local spin exists)")

    # link structure: T_i has elements only between neighboring sites
    # diagonal of T_i (in position basis) should be zero
    for lab, T in [('Tx', Tx), ('Ty', Ty), ('Tz', Tz)]:
        diag = np.sum(np.abs(np.diag(T)) ** 2)
        tot = np.sum(np.abs(T) ** 2)
        print(f"  {lab}: diagonal weight = {diag:.6f}, total = {tot:.1f}  "
              f"-> link (off-diagonal) operator: {diag < 1e-12}")

    # what DOES T_i look like? nearest-neighbor shift
    print("\n=== T_i is a nearest-neighbor shift (link), not a site rotation ===")
    print(f"  T_x non-zero entries: {np.count_nonzero(np.abs(Tx) > 1e-12)} (= N links x-dir)")
    print(f"  a site spin would have N diagonal blocks (each 2x2 for spin-1/2)")
    print("=> paid-bridge-2 = link SU(2) (magnetic translations) -> site SU(2) (local spin)")


if __name__ == "__main__":
    main()

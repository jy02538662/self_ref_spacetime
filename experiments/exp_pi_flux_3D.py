"""Spin connection / third layer: 3D pi-flux (思路 3).

The 2D pi-flux gives ONE anti-commuting pair of magnetic translations (上下 = T_xT_y = -T_yT_x).
The 3D pi-flux (each face holonomy -1) gives THREE anti-commuting pairs — a richer
non-abelian structure, and it is Kramers (T^2=-1, all-even spectrum) for even size.

Gauge: x-edge phase (-1)^(y+z), y-edge (-1)^z, z-edge +1  => every square face has holonomy -1.
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
    ev = np.linalg.eigvalsh(D)
    uniq, cnt = np.unique(np.round(ev, 8), return_counts=True)
    print(f"3D torus C_{L}^3 (N={N}), pi-flux: multiplicities={list(cnt)}  all-even(Kramers)={all(c%2==0 for c in cnt)}")

    Tx, Ty, Tz = mag_trans(D, L, 'x'), mag_trans(D, L, 'y'), mag_trans(D, L, 'z')
    print("\nanti-commutation of magnetic translations:")
    for A, B, lab in [(Tx, Ty, 'Tx,Ty'), (Ty, Tz, 'Ty,Tz'), (Tz, Tx, 'Tz,Tx')]:
        ac = np.linalg.norm(A @ B + B @ A)  # anticommutator {A,B}
        co = np.linalg.norm(A @ B - B @ A)  # commutator [A,B]
        print(f"  {{{lab}}}={ac:.2e}  [{lab}]={co:.2e}   => anti-commute: {ac < 1e-9 and co > 1}")

    print(f"\nT_i^2 = translation-by-2 (not I):  Tx^2==I? {np.allclose(Tx@Tx, np.eye(N))}")
    print(f"T_i^L = I (periodic):  Tx^4==I? {np.allclose(np.linalg.matrix_power(Tx, L), np.eye(N))}")

    print("\n=> 3D pi-flux 给三对反对易磁平移（3D「上下」，vs 2D 的一对）。")
    print("   这是「扭结/投影表示」的非交换结构（T_i^L=I、两两反交换），不是标准 Clifford（T_i^2=I）。")
    print("   物理内容 = 3D 里 π 磁通的非交换性更丰富（三对反对易 vs 一对）。")


if __name__ == "__main__":
    main()

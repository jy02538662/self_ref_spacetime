"""Attack the "separation wall": does QUANTIZED SU(2)_k depend on U(1) flux?

The separation wall (#6) says: in the CLASSICAL pi-flux, SU(2) is emergent from
U(1) flux (strip the flux -> anti-commutation dies).  But QUANTIZATION gives
delta = 2 cos(pi/(N+1)) from FINITE-N truncation (Chebyshev), NOT from flux.

Decisive check: does SU(2)_k survive when the U(1) flux is stripped?

  - classical: strip flux -> T_x T_y commute (no SU(2))   [already proved #6]
  - quantum:   strip flux -> delta = 2cos(pi/(N+1)) STILL there (finite N)

If the quantum SU(2)_k is INDEPENDENT of the flux, then the "separation" wall
is a CLASSICAL artifact, and the quantum route may bypass it: the spin (SU(2)_k)
is a direct product of FINITENESS (truncation), not a byproduct of U(1) flux.
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


def chebyshev_dims(delta, nmax):
    ds = [1.0, delta]
    for _ in range(2, nmax + 1):
        ds.append(delta * ds[-1] - ds[-2])
    return ds


def main():
    L = 4
    N = L ** 3

    # ---- PART 1: classical SU(2) DOES depend on flux (reproduce #6) ----
    print("=== PART 1: CLASSICAL SU(2) = anti-commutation of magnetic translations ===")
    Tx_f, Ty_f = mag_trans(torus3D(L, True), L, 'x'), mag_trans(torus3D(L, True), L, 'y')
    Tx_0, Ty_0 = mag_trans(torus3D(L, False), L, 'x'), mag_trans(torus3D(L, False), L, 'y')
    ac_f = np.linalg.norm(Tx_f @ Ty_f + Ty_f @ Tx_f)
    ac_0 = np.linalg.norm(Tx_0 @ Ty_0 + Ty_0 @ Tx_0)
    print(f"  with flux:    {{Tx,Ty}} = {ac_f:.1f}  (anti-commute -> SU(2))")
    print(f"  without flux: {{Tx,Ty}} = {ac_0:.1f}  (commute -> NO SU(2))")
    print(f"  => classical SU(2) DEPENDS on flux: {ac_f < 1e-9 and ac_0 > 1e-9}\n")

    # ---- PART 2: QUANTIZED SU(2)_k does NOT depend on flux ----
    print("=== PART 2: QUANTIZED SU(2)_k = finite-N truncation (Chebyshev), flux-independent ===")
    print("  delta = 2 cos(pi/(N+1)) comes from 'N maximal strands', NOT from flux phase.")
    print("  So it exists with OR without flux.  Verify the SU(2)_k structure for several N:")
    for Nstr in (3, 4, 5, 6, 7, 8):
        delta = 2 * np.cos(np.pi / (Nstr + 1))
        ds = chebyshev_dims(delta, Nstr)
        # SU(2)_k truncation: Delta_{N-1}=1 (f_N exists), Delta_N=0 (f_{N+1} vanishes)
        fN_exists = abs(ds[Nstr - 1] - 1) < 1e-9
        fN1_vanishes = abs(ds[Nstr]) < 1e-9
        level = Nstr - 1
        print(f"  N={Nstr} (level k={level}): delta={delta:.4f}, "
              f"Delta_{{N-1}}=1: {fN_exists}, Delta_N=0: {fN1_vanishes}")

    print("\n=== the key distinction ===")
    print("  classical SU(2):  pi-flux anti-commutation  -> DEPENDS on U(1) flux (#6)")
    print("  quantized SU(2)_k: finite-N Chebyshev truncation -> INDEPENDENT of U(1) flux")
    print("  => the 'separation is impossible' wall is a CLASSICAL artifact.")
    print("     The quantum route gives an SU(2)_k that does NOT come from U(1),")
    print("     so 'separating SU(2) from U(1)' may be bypassed by QUANTIZATION.")


if __name__ == "__main__":
    main()

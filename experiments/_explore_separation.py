"""Trivial omega_mu: how far can we push the "separation" of SU(2) from U(1)?

This round proved: omega_mu = magnetic translation T_mu (a UNIFIED object =
translation + U(1) flux + SU(2) spin).  The "separation" (paid-bridge-2) means
splitting off a PURE SU(2) part from the U(1) part.

Decisive check: does the anti-commutation (SU(2)) survive if we strip the U(1)
phase from T_mu?

  - T_mu (with phase +-1): anti-commute T_x T_y = -T_y T_x (Cl(3) -> SU(2))
  - |T_mu| (pure translation, no phase): COMMUTE (no Cl(3) -> no SU(2))

If anti-commutation vanishes when phase is stripped, then SU(2) is EMERGENT from
the U(1) phase (Aharonov-Bohm), i.e. SU(2) and U(1) are INSEPARABLE inside T_mu.
This is the precise content of the "separation" wall.
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
    Dflux = torus3D(L, flux=True)
    Dnoflux = torus3D(L, flux=False)
    N = L ** 3

    # with flux (phase +-1)
    Tx, Ty, Tz = (mag_trans(Dflux, L, a) for a in 'xyz')
    # no flux (pure translation, phase all +1)
    Px, Py, Pz = (mag_trans(Dnoflux, L, a) for a in 'xyz')

    print("=== does anti-commutation (SU(2)) survive stripping the U(1) phase? ===")
    for lab, A, B in [('flux', Tx, Ty), ('no-flux', Px, Py)]:
        ac = np.linalg.norm(A @ B + B @ A)   # anticommutator
        co = np.linalg.norm(A @ B - B @ A)   # commutator
        print(f"  {lab:8s}: {{Tx,Ty}}={ac:.2e}  [Tx,Ty]={co:.2e}  "
              f"anti-commute={ac < 1e-9}  (SU(2)/Cl(3) present={ac < 1e-9})")

    print("\n=== what the phase does ===")
    print("  flux:    T_x T_y = -T_y T_x  (anti-commute -> Cl(3) -> SU(2))")
    print("  no-flux: T_x T_y = +T_y T_x  (commute -> abelian, NO SU(2))")
    print("  => the SU(2) (anti-commutation) is EMERGENT from the U(1) phase (pi-flux).")
    print("     SU(2) and U(1) are INSEPARABLE inside T_mu: you cannot strip the U(1)")
    print("     phase without killing the SU(2).")

    # the "pure SU(2)" object people want for omega_mu: is there one?
    print("\n=== is there a PURE SU(2) (no U(1) phase) omega_mu in this system? ===")
    # the anti-commuting pair is T_x, T_y (magnetic translations), which CARRY the
    # U(1) phase.  There is no SU(2) without the phase.
    print("  The only SU(2) generators ARE the magnetic translations T_i,")
    print("  which carry the U(1) phase.  A pure-SU(2) omega_mu (phase-stripped)")
    print("  does NOT exist in the pi-flux system: SU(2) is the phase's non-commutativity.")
    print("  => trivial omega_mu = T_mu (unified U(1)xSU(2)); SEPARATION is the wall.")


if __name__ == "__main__":
    main()

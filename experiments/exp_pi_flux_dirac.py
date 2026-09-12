"""pi-flux (Z2 flux) toroidal spectrum -- Dirac point / zero-mode scout.

The bridge to test: pi-flux = Z2 = spin-1/2 (fractionalization).
Chain: pi-flux (Z2) -> Dirac points (massless fermion) -> topology (Chern) -> Majorana.

Step 1: does pi-flux toroidal D (even N) develop Dirac points (E=0 degeneracy)
vs no-flux (which has a Fermi surface, not isolated Dirac points)?
"""

from __future__ import annotations

import numpy as np

from experiments.exp_pure_spectral_anneal import toroidal_D


def spectrum_stats(D, tol=1e-9):
    eig = np.linalg.eigvalsh(D)
    zero_modes = int(np.sum(np.abs(eig) < tol))
    nz = eig[np.abs(eig) >= tol]
    gap = float(np.min(np.abs(nz))) if nz.size else 0.0
    # chiral symmetry: are eigenvalues +- paired (bipartite)?
    pos = eig[eig > tol]
    neg = eig[eig < -tol]
    chiral = float(np.max(np.abs(pos + np.sort(neg)[::-1]))) if pos.size == neg.size else float('nan')
    return eig, zero_modes, gap, chiral


def main():
    print(f"{'N':>4} {'N2':>4} {'flux':>7} {'n_zero':>7} {'gap':>8} {'chiral_err':>12} {'eig_min':>8} {'eig_max':>8}")
    for N in [4, 6, 8, 10, 12]:
        for flux, label in [(True, 'pi'), (False, '0')]:
            D = toroidal_D(N, pi_flux=flux)
            eig, nz, gap, chiral = spectrum_stats(D)
            print(f"{N:>4} {N*N:>4} {label:>7} {nz:>7} {gap:>8.4f} {chiral:>12.3e} {eig.min():>8.3f} {eig.max():>8.3f}")

    # detail: print the full spectrum for N=8 pi-flux vs no-flux
    print("\nFull spectrum N=8 (sorted):")
    for flux, label in [(True, 'pi-flux'), (False, 'no-flux')]:
        D = toroidal_D(8, pi_flux=flux)
        eig = np.sort(np.linalg.eigvalsh(D))
        print(f"  {label}: {np.round(eig, 3).tolist()}")


if __name__ == "__main__":
    main()

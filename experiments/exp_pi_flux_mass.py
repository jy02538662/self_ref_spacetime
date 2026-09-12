"""pi-flux toroidal + staggered mass -> open gap? (confirm Dirac points, not accidental zero modes).

If the 4 zero modes are genuine Dirac points, a staggered (sublattice) mass m
should open a gap ~ |m| (massive Dirac fermion). Then Chern number becomes +-2
(Hofstadter / Haldane-type), the topological origin of spin-1/2.

H = D + m * diag(+1,-1,+1,-1,...)  (sublattice mass on the bipartite square lattice)
"""

from __future__ import annotations

import numpy as np

from experiments.exp_pure_spectral_anneal import toroidal_D


def main():
    N = 8
    D = toroidal_D(N, pi_flux=True)
    # staggered mass: A/B sublattice. node idx = N*i + j; sublattice sign = (-1)^(i+j)
    signs = np.array([(-1) ** (i + j) for i in range(N) for j in range(N)])
    M = np.diag(signs.astype(float))

    print(f"N={N} pi-flux + staggered mass m:  gap and spectrum near E=0")
    for m in [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]:
        H = D + m * M
        eig = np.linalg.eigvalsh(H)
        nz = eig[np.abs(eig) > 1e-9]
        gap = float(np.min(np.abs(nz))) if nz.size else 0.0
        n_zero = int(np.sum(np.abs(eig) < 1e-9))
        # lowest few eigenvalues
        low = np.sort(np.abs(eig))[:6]
        print(f"  m={m:.1f}: n_zero={n_zero}  gap={gap:.4f}  low|eig|={np.round(low, 4).tolist()}")


if __name__ == "__main__":
    main()

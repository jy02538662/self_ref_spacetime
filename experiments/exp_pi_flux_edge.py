"""pi-flux + staggered mass on a strip (open x, periodic y) -> edge states?

If Chern number != 0, opening the gap with mass m creates in-gap edge states
localized on the open boundaries. This pins down "topological protection"
(the third link: pi-flux(Z2) -> Dirac -> Chern -> protected spin-1/2).

Hamiltonian: hopping (pi-flux: vertical phase = pi*i) + staggered mass m*(-1)^(i+j).
"""

from __future__ import annotations

import numpy as np


def build_strip(Nx, Ny, m):
    n = Nx * Ny
    H = np.zeros((n, n), complex)

    def idx(i, j):
        return i * Ny + j

    for i in range(Nx):
        for j in range(Ny):
            k = idx(i, j)
            # staggered mass
            H[k, k] = m * ((-1) ** (i + j))
            # horizontal: (i,j) <-> (i+1,j), phase 0
            if i + 1 < Nx:
                k2 = idx(i + 1, j)
                H[k, k2] += 1.0
                H[k2, k] += 1.0
            # vertical: (i,j) <-> (i,j+1), phase pi*i (pi-flux)
            j2 = (j + 1) % Ny
            k3 = idx(i, j2)
            ph = np.pi * i
            H[k, k3] += np.exp(1j * ph)
            H[k3, k] += np.exp(-1j * ph)
    return H


def main():
    Ny = 12
    for Nx in [8, 16]:
        for m in [0.5, 1.0]:
            H = build_strip(Nx, Ny, m)
            eig = np.linalg.eigvalsh(H)
            # in-gap states: |E| < m (inside the Dirac gap)
            in_gap = eig[np.abs(eig) < m - 1e-6]
            # edge-localized check for in-gap states
            print(f"Nx={Nx} (open), Ny={Ny}, m={m}: gap={m:.1f}  in-gap states={in_gap.size}  "
                  f"their |E|={np.round(np.abs(in_gap), 4).tolist() if in_gap.size else []}")
    # detailed: Nx=16, m=1.0, print in-gap eigenvalues + localization
    print("\nDetail Nx=16, Ny=12, m=1.0:")
    H = build_strip(16, 12, 1.0)
    eig, vec = np.linalg.eigh(H)
    in_gap_idx = np.where(np.abs(eig) < 1.0 - 1e-6)[0]
    for k in in_gap_idx:
        v = vec[:, k]
        # localization: weight near left edge (i=0,1) vs right edge (i=Nx-2,Nx-1)
        p = np.abs(v) ** 2
        p = p.reshape(16, 12).sum(axis=1)  # per-row (x) weight
        left = p[0] + p[1]
        right = p[-2] + p[-1]
        edge_frac = (left + right) / p.sum()
        print(f"  E={eig[k]:.6f}  edge_frac={edge_frac:.3f}  (left {left:.3f}, right {right:.3f})")


if __name__ == "__main__":
    main()

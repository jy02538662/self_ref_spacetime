"""pi-flux + staggered mass + Haldane next-nearest hopping -> Chern != 0?

Scheme A (Haldane-type): add a time-reversal-breaking next-nearest-neighbor
hopping t2*e^{i phi} (complex phase) to the pi-flux + staggered-mass Dirac fermion.
This un-balances the two valleys, so their +-1/2 Chern contributions no longer
cancel -> Chern = +-1 or +-2 -> in-gap edge states appear.

H = nearest-neighbor hopping (pi-flux) + staggered mass m*(-1)^(i+j)
  + next-nearest hopping t2*e^{i phi} (diagonal, TR-breaking).

Strip: open x (Nx), periodic y (Ny). Diagnose in-gap states (edge states).
"""

from __future__ import annotations

import numpy as np


def build_strip(Nx, Ny, m, t2, phi):
    n = Nx * Ny
    H = np.zeros((n, n), complex)

    def idx(i, j):
        return i * Ny + j

    for i in range(Nx):
        for j in range(Ny):
            k = idx(i, j)
            # staggered mass
            H[k, k] = m * ((-1) ** (i + j))
            # nearest-neighbor horizontal: (i,j) <-> (i+1,j), phase 0
            if i + 1 < Nx:
                k2 = idx(i + 1, j)
                H[k, k2] += 1.0
                H[k2, k] += 1.0
            # nearest-neighbor vertical: (i,j) <-> (i,j+1), phase pi*i (pi-flux)
            j2 = (j + 1) % Ny
            k3 = idx(i, j2)
            ph = np.pi * i
            H[k, k3] += np.exp(1j * ph)
            H[k3, k] += np.exp(-1j * ph)
            # next-nearest (diagonal): (i,j) <-> (i+1, j+1) and (i+1, j-1), phase +phi
            if i + 1 < Nx:
                jp = (j + 1) % Ny
                k4 = idx(i + 1, jp)
                H[k, k4] += t2 * np.exp(1j * phi)
                H[k4, k] += t2 * np.exp(-1j * phi)
                jm = (j - 1) % Ny
                k5 = idx(i + 1, jm)
                H[k, k5] += t2 * np.exp(-1j * phi)
                H[k5, k] += t2 * np.exp(1j * phi)
    return H


def count_in_gap(H, m, tol=1e-6):
    eig = np.linalg.eigvalsh(H)
    in_gap = eig[np.abs(eig) < m - tol]
    return in_gap


def main():
    Nx, Ny = 16, 12
    m = 1.0
    print(f"pi-flux + staggered mass m={m} + Haldane NNN t2*e^{{i phi}}  (strip Nx={Nx}, Ny={Ny})")
    print(f"{'t2':>5} {'phi/pi':>7} {'n_edge':>7} {'edge |E|':>20}")
    for t2 in [0.0, 0.2, 0.4, 0.6, 0.8]:
        for phi in [np.pi / 2]:
            H = build_strip(Nx, Ny, m, t2, phi)
            in_gap = count_in_gap(H, m)
            print(f"{t2:>5.2f} {phi/np.pi:>7.3f} {in_gap.size:>7} "
                  f"{np.round(np.abs(np.sort(in_gap)), 4).tolist() if in_gap.size else []}")

    # scan phi at fixed t2
    print("\nscan phi at t2=0.5:")
    t2 = 0.5
    for phi in [0.0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]:
        H = build_strip(Nx, Ny, m, t2, phi)
        in_gap = count_in_gap(H, m)
        print(f"  phi/pi={phi/np.pi:.3f}: n_edge={in_gap.size}  |E|={np.round(np.abs(np.sort(in_gap)),4).tolist() if in_gap.size else []}")

    # detailed edge localization for a nontrivial case
    print("\nDetail t2=0.5, phi=pi/2: edge-state localization")
    H = build_strip(Nx, Ny, m, 0.5, np.pi / 2)
    eig, vec = np.linalg.eigh(H)
    for k in np.where(np.abs(eig) < m - 1e-6)[0]:
        v = vec[:, k]
        p = (np.abs(v) ** 2).reshape(Nx, Ny).sum(axis=1)
        left = p[0] + p[1]
        right = p[-2] + p[-1]
        print(f"  E={eig[k]:.5f}  edge_frac={(left+right)/p.sum():.3f}")


if __name__ == "__main__":
    main()

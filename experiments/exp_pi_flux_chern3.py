"""Chern number of pi-flux square lattice (alpha=1/2) + staggered mass + TR-breaking term.

Bloch Hamiltonian (magnetic unit cell 2x1, Landau gauge vertical phase pi*i):
    H(k) = [ [2t cos ky + m,  2t cos kx e^{-i kx}],
             [2t cos kx e^{i kx},  -2t cos ky - m] ]

m = staggered (sublattice) mass.
TR-breaking term candidates (to add):
  - Haldane-type NNN:  + g * sin(2kx) * sin(ky) * sigma_z   (complex, TR-odd)
  - "valley-Same" mass: + g * sin(kx) * sigma_y

Compute Chern via discrete Berry curvature (Fukui-Hatsugai-Suzuki).
"""

from __future__ import annotations

import numpy as np


def H_of_k(kx, ky, m, g, kind):
    H = np.zeros((2, 2), complex)
    H[0, 0] = 2 * np.cos(ky) + m
    H[1, 1] = -2 * np.cos(ky) - m
    H[0, 1] = 2 * np.cos(kx) * np.exp(-1j * kx)
    H[1, 0] = np.conjugate(H[0, 1])
    if kind == "none":
        pass
    elif kind == "tr_sz":
        # sigma_z * sin(ky): odd in ky -> TR-odd; opposite mass at the two valleys
        H[0, 0] += g * np.sin(ky)
        H[1, 1] -= g * np.sin(ky)
    elif kind == "haldane_sz":
        # sigma_z, complex TR-odd term (vanishes at Dirac points, sin(2kx)=0 there)
        H[0, 0] += g * np.sin(2 * kx) * np.sin(ky)
        H[1, 1] -= g * np.sin(2 * kx) * np.sin(ky)
    elif kind == "tr_sy":
        # sigma_y, imaginary TR-odd
        H[0, 1] += 1j * g * np.sin(kx)
        H[1, 0] += -1j * g * np.sin(kx)
    return H


def chern(kind, m, g, Nx=80, Ny=80):
    # magnetic BZ: kx in [0, pi], ky in [0, 2pi]
    kxs = np.linspace(0, np.pi, Nx, endpoint=False)
    kys = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
    # occupied band projector
    occ = np.zeros((Nx, Ny, 2), complex)
    for i, kx in enumerate(kxs):
        for j, ky in enumerate(kys):
            H = H_of_k(kx, ky, m, g, kind)
            w, v = np.linalg.eigh(H)
            occ[i, j] = v[:, 0]  # lower band
    # link variables
    total = 0.0
    for i in range(Nx):
        for j in range(Ny):
            ip = (i + 1) % Nx
            jp = (j + 1) % Ny
            Ux = np.vdot(occ[i, j], occ[ip, j])
            Uy = np.vdot(occ[i, j], occ[i, jp])
            Uxp = np.vdot(occ[ip, j], occ[ip, jp])  # U_y at (i+1,j)
            Uyp = np.vdot(occ[i, jp], occ[ip, jp])  # U_x at (i,j+1)
            total += np.angle(Ux * Uxp / (Uyp * Uy))
    return total / (2 * np.pi)


def main():
    print("Chern number of pi-flux square lattice (alpha=1/2):")
    print(f"{'kind':>12} {'m':>4} {'g':>4} {'Chern':>8}")
    for m in [0.5]:
        for kind in ["none", "tr_sz", "haldane_sz", "tr_sy"]:
            for g in ([0.0] if kind == "none" else [0.3, 0.6, 0.8, 1.2]):
                C = chern(kind, m, g)
                print(f"{kind:>12} {m:>4} {g:>4} {C:>8.3f}")


if __name__ == "__main__":
    main()

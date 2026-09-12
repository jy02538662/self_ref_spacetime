"""Chern via d-hat winding for pi-flux square lattice + real NNN.

Key finding: for the plain pi-flux square lattice (Landau gauge),
    d_x = Re(H01) = 2 cos^2(kx) >= 0  ALWAYS,
so d-hat lives in the x>=0 hemisphere -> winding = 0 -> Chern = 0 identically.
This is the chiral structure of the pi-flux Dirac points.

To get Chern != 0 we must let d_x change sign: add a REAL next-nearest hopping
t2*cos(ky) (which contributes a real part to H01 that can be negative).

H01 = 2 cos(kx) e^{-i kx} + t2 cos(ky)
d = (Re H01, Im H01, 2 cos ky + m + g sin ky)

Chern = (1/4pi) int d^hat . (d_x d^hat x d_y d^hat) d^2 k
"""

from __future__ import annotations

import numpy as np


def dvec(kx, ky, m, g, t2):
    H01 = 2 * np.cos(kx) * np.exp(-1j * kx) + t2 * np.cos(ky)
    return np.array([np.real(H01), np.imag(H01), 2 * np.cos(ky) + m + g * np.sin(ky)])


def winding(m, g, t2, N=200):
    kxs = np.linspace(0, np.pi, N)
    kys = np.linspace(0, 2 * np.pi, N)
    dkx_step = np.pi / (N - 1)
    dky_step = 2 * np.pi / (N - 1)
    C = 0.0
    for i in range(N - 1):
        for j in range(N - 1):
            kx, ky = kxs[i], kys[j]
            d = dvec(kx, ky, m, g, t2)
            n = np.linalg.norm(d)
            if n < 1e-10:
                continue
            dx = (dvec(kx + dkx_step, ky, m, g, t2) - d) / dkx_step
            dy = (dvec(kx, ky + dky_step, m, g, t2) - d) / dky_step
            C += float(np.dot(d / n, np.cross(dx, dy)))
    return C / (4 * np.pi) * dkx_step * dky_step


def main():
    print("Chern via d-hat winding (pi-flux square lattice + real NNN t2):")
    print(f"{'m':>4} {'g':>4} {'t2':>5} {'Chern':>8}")
    for t2 in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5]:
        C = winding(0.5, 1.2, t2, 150)
        print(f"0.5  1.2  {t2:>5.1f}  {C:>8.3f}")
    # also t2 with g=0 (pure staggered mass + NNN, should give Chern too)
    print("\nwith g=0 (staggered mass m=0.5 + NNN t2):")
    for t2 in [1.0, 2.0, 2.5]:
        C = winding(0.5, 0.0, t2, 150)
        print(f"  0.5  0.0  {t2:>5.1f}  {C:>8.3f}")


if __name__ == "__main__":
    main()

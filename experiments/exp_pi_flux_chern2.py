"""pi-flux + staggered mass, tune flux alpha = 1/2 -> 1/2 + delta.

At alpha=1/2 the flux is Z2 (holonomy = pi = real), time-reversal symmetric,
so Chern = 0 (valley Hall). Tuning alpha away from 1/2 by delta makes the flux
complex -> breaks time-reversal -> the two valleys no longer cancel -> Chern != 0.

This is the cleanest TR-breaking knob. Diagnose in-gap EDGE states (edge_frac ~ 1).

H = nearest hopping (flux alpha: vertical phase = 2*pi*alpha*i) + staggered mass m.
"""

from __future__ import annotations

import numpy as np


def build_strip(Nx, Ny, m, alpha):
    n = Nx * Ny
    H = np.zeros((n, n), complex)

    def idx(i, j):
        return i * Ny + j

    for i in range(Nx):
        for j in range(Ny):
            k = idx(i, j)
            H[k, k] = m * ((-1) ** (i + j))
            if i + 1 < Nx:
                k2 = idx(i + 1, j)
                H[k, k2] += 1.0
                H[k2, k] += 1.0
            j2 = (j + 1) % Ny
            k3 = idx(i, j2)
            ph = 2 * np.pi * alpha * i
            H[k, k3] += np.exp(1j * ph)
            H[k3, k] += np.exp(-1j * ph)
    return H


def edge_states(H, m, tol=1e-6):
    Nx = int(round(np.sqrt(H.shape[0] / 1)))  # placeholder
    eig, vec = np.linalg.eigh(H)
    # infer Nx, Ny: total = Nx*Ny. We'll pass explicitly instead.
    return eig, vec


def main():
    Nx, Ny = 16, 12
    m = 1.0
    print(f"pi-flux alpha=1/2+delta + staggered mass m={m}  (strip Nx={Nx}, Ny={Ny})")
    print(f"{'delta':>7} {'n_edge':>7} {'edge |E| (top 6)':>28}")
    for delta in [0.0, 0.02, 0.05, 0.1, 0.15, 0.2]:
        H = build_strip(Nx, Ny, m, 0.5 + delta)
        eig, vec = np.linalg.eigh(H)
        in_gap = np.where(np.abs(eig) < m - 1e-6)[0]
        # localization of in-gap states
        edge_fracs = []
        for k in in_gap:
            v = vec[:, k]
            p = (np.abs(v) ** 2).reshape(Nx, Ny).sum(axis=1)
            edge_fracs.append((p[0] + p[1] + p[-2] + p[-1]) / p.sum())
        # count only TRUE edge states (edge_frac > 0.6)
        true_edge = sum(1 for f in edge_fracs if f > 0.6)
        topE = np.round(np.abs(np.sort(eig[in_gap]))[:6], 3).tolist() if len(in_gap) else []
        print(f"{delta:>7.2f} {len(in_gap):>7} {topE}")
        print(f"         true_edge_states (edge_frac>0.6): {true_edge}")

    # detail: delta=0.1, print edge states localization
    print("\nDetail delta=0.1: in-gap states and their edge_frac")
    H = build_strip(Nx, Ny, m, 0.6)
    eig, vec = np.linalg.eigh(H)
    for k in np.where(np.abs(eig) < m - 1e-6)[0]:
        v = vec[:, k]
        p = (np.abs(v) ** 2).reshape(Nx, Ny).sum(axis=1)
        edge_frac = (p[0] + p[1] + p[-2] + p[-1]) / p.sum()
        print(f"  E={eig[k]:.5f}  edge_frac={edge_frac:.3f}")


if __name__ == "__main__":
    main()

"""Task 3: spectral dimension vs system size N for toroidal grids.

Check whether specdim -> 2 as N grows (N=4,9,16,25,36 = 2x2..6x6 torus).
Compare D^2-based heat kernel vs standard graph Laplacian L = deg*I - D.
"""

from __future__ import annotations

import numpy as np

from experiments.exp_pure_spectral_anneal import toroidal_D


def specdim_from_eig(eig, t_lo=0.5, t_hi=50.0, n_pts=40):
    eig = np.clip(eig, 1e-14, None)
    ts = np.geomspace(t_lo, t_hi, n_pts)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    slope = np.polyfit(np.log(ts), logK, 1)[0]
    return float(-2.0 * slope)


def main():
    print(f"{'grid':>8} {'N':>4} {'specdim(D^2)':>14} {'specdim(Laplacian)':>20}")
    for dim in [2, 3, 4, 5, 6, 8, 10]:
        n = dim ** 2
        D = toroidal_D(dim, pi_flux=False)
        r = np.abs(D)
        deg = np.sum(r, axis=1)  # = 4 for toroidal
        L = np.diag(deg) - D  # graph Laplacian (weighted, no flux)
        D2 = D @ D
        eig_D2 = np.linalg.eigvalsh(D2)
        eig_L = np.linalg.eigvalsh(L)
        sd_D2 = specdim_from_eig(eig_D2)
        sd_L = specdim_from_eig(eig_L)
        print(f"{dim}x{dim}    {n:>4}  {sd_D2:>14.3f}  {sd_L:>20.3f}")


if __name__ == "__main__":
    main()

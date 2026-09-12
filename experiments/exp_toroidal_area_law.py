"""Area-law check for the spontaneously-emerged toroidal + pi-flux.

The pure-spectral action grew a toroidal 3x3 + pi-flux. Is that "space" in the
entanglement sense? Criterion (preprint 2): area law S ~ L^{d-1} <=> Fermi surface
<=> spatiality; volume law S ~ L^d <=> no Fermi surface <=> non-spatial.

Here we build toroidal NxN + pi-flux (the emerged structure), treat D as a free
fermion hopping, take the Fermi-sea ground state (lowest half of eigenstates),
and measure entanglement entropy of a l x l corner subregion via Peschel.

Contrast: pi-flux vs no-flux vs random graph (random -> volume law).
"""

from __future__ import annotations

import numpy as np

from experiments.exp_pure_spectral_anneal import toroidal_D


def entropy_from_corr(C_A):
    lam = np.linalg.eigvalsh(C_A)
    lam = np.clip(lam, 1e-15, 1.0 - 1e-15)
    return float(-np.sum(lam * np.log(lam) + (1.0 - lam) * np.log(1.0 - lam)))


def fermi_sea_corr(D, fill=0.5):
    N = D.shape[0]
    _, vec = np.linalg.eigh(D)
    n_occ = int(round(fill * N))
    occ = vec[:, :n_occ]
    return occ @ occ.conj().T


def corner_indices(N, l):
    return [i * N + j for i in range(l) for j in range(l)]


def random_graph(N, p=0.5, seed=0):
    rng = np.random.default_rng(seed)
    A = np.triu(rng.random((N, N)) < p, k=1).astype(float)
    A = A + A.T
    return A


def main():
    print("=== Area law for emerged toroidal + pi-flux (Fermi sea, half filling) ===")
    for N in [10, 14, 18]:
        D_pi = toroidal_D(N, pi_flux=True)
        D_0 = toroidal_D(N, pi_flux=False)
        C_pi = fermi_sea_corr(D_pi)
        C_0 = fermi_sea_corr(D_0)
        C_r = fermi_sea_corr(random_graph(N * N, p=0.5))
        print(f"\nN={N} (N^2={N*N}):")
        print(f"  {'l':>3} {'S_pi':>8} {'S_0':>8} {'S_rand':>8}  |  S_pi/l  S_pi/l^2")
        lp = []
        Sp = []
        for l in range(2, N // 2 + 1):
            idx = corner_indices(N, l)
            Spi = entropy_from_corr(C_pi[np.ix_(idx, idx)])
            S0 = entropy_from_corr(C_0[np.ix_(idx, idx)])
            Sr = entropy_from_corr(C_r[np.ix_(idx, idx)])
            lp.append(l)
            Sp.append(Spi)
            print(f"  {l:>3} {Spi:>8.3f} {S0:>8.3f} {Sr:>8.3f}  |  {Spi/l:7.3f}  {Spi/(l*l):8.3f}")
        # fit slope of S_pi vs l (log-log)
        lp = np.array(lp, float)
        Sp = np.array(Sp)
        slope = float(np.polyfit(np.log(lp), np.log(Sp), 1)[0])
        print(f"  slope d(S_pi)/d(log l) = {slope:.2f}  (area law ~1, volume ~2)")


if __name__ == "__main__":
    main()

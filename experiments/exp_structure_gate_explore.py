"""Structure-selection gate: explore Tr(D^4) for the 4x4 torus.

Goal: reproduce the note's numbers (pi-flux Tr=384, no-flux Tr=640), list the
24 four-cycles, and find the TRUE minimum over phase (gauge) assignments, not
just the standard pi-flux gauge.

Key identity (from the note, section 2 / 16), for Hermitian D, zero diagonal:

    Tr(D^4) = [2 sum_h d_h^2 - sum_{i!=j} r_ij^4]  +  sum_{distinct i,j,k,l} D_ij D_jk D_kl D_li

For a d-regular graph with all |D_ij|=r=1 on edges the modulus term is
Nd(2d-1) = 448 (N=16,d=4) constant across all 4-regular graphs.  The last term
(holonomy) = 8 * sum_{undirected 4-cycles C} cos(Phi_C).
"""

from __future__ import annotations

import numpy as np
from itertools import combinations


def torus_D(Lx, Ly, mode="pi"):
    N = Lx * Ly
    D = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (y % Ly) * Lx + (x % Lx)

    for y in range(Ly):
        for x in range(Lx):
            i = idx(x, y)
            j = idx(x, y + 1)          # vertical
            D[i, j] = 1.0
            D[j, i] = 1.0
            k = idx(x + 1, y)          # horizontal
            w = (-1.0) ** y if mode == "pi" else 1.0
            D[i, k] = w
            D[k, i] = w
    return D


def edge_set_from_D(D, thr=1e-9):
    n = D.shape[0]
    es = set()
    for i in range(n):
        for j in range(i + 1, n):
            if abs(D[i, j]) > thr:
                es.add(frozenset((i, j)))
    return es


def four_cycles(n, edge_set):
    """All distinct undirected 4-cycles, as cyclic-order tuples (a,b,c,d).

    For each 4-subset, the 3 possible cyclic orderings are the 3 distinct
    Hamiltonian cycles; keep those whose 4 consecutive edges all exist.
    """
    cycles = []
    for S in combinations(range(n), 4):
        a, b, c, d = S
        for cyc in [(a, b, c, d), (a, b, d, c), (a, c, b, d)]:
            if all(frozenset((cyc[i], cyc[(i + 1) % 4])) in edge_set for i in range(4)):
                cycles.append(cyc)
    return cycles


def holonomy_from_D(D, cycles):
    """8 * sum_C Re(prod_C), the exact holonomy term (works for any r)."""
    s = 0.0
    cosvals = []
    for (a, b, c, d) in cycles:
        prod = D[a, b] * D[b, c] * D[c, d] * D[d, a]
        s += 8.0 * float(np.real(prod))
        cosvals.append(float(np.real(prod)))
    return s, cosvals


def analyze(D, label):
    n = D.shape[0]
    tr4 = float(np.real(np.trace(np.linalg.matrix_power(D, 4))))
    es = edge_set_from_D(D)
    d_h = np.real(np.diag(D @ D))
    modulus_term = float(2 * np.sum(d_h ** 2) - np.sum(np.abs(D) ** 4))
    cyc = four_cycles(n, es)
    holo, cosvals = holonomy_from_D(D, cyc)
    print(f"\n=== {label} ===")
    print(f"Tr(D^4) direct      = {tr4:.6f}")
    print(f"modulus term        = {modulus_term:.6f}")
    print(f"n_4cycles           = {len(cyc)}")
    print(f"8*sum Re(prod)      = {holo:.6f}  -> modulus+holo = {modulus_term+holo:.6f}  err={abs(tr4-(modulus_term+holo)):.2e}")
    print(f"cos values: {sorted(round(c, 3) for c in cosvals)}")
    return tr4, modulus_term, cyc, cosvals, holo


def main():
    print("### 4x4 torus, r=1, N=16, d=4, modulus term should be 448")
    for mode in ["pi", "zero"]:
        analyze(torus_D(4, 4, mode=mode), f"{mode}-flux")

    # sanity: random Hermitian 7x7 identity check
    rng = np.random.default_rng(0)
    D = rng.normal(0, 1, (7, 7)) + 1j * rng.normal(0, 1, (7, 7))
    D = (D + D.conj().T) / 2
    np.fill_diagonal(D, 0.0)
    analyze(D, "random Hermitian 7x7")


if __name__ == "__main__":
    main()

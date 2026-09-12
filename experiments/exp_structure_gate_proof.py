"""Structure-selection gate: PROOF that pi-flux toroidal is the global minimum.

Key result (rigorous, no enumeration needed):

    For any N-vertex d-regular graph with uniform edge modulus |D_ij| = r = 1:
        Tr(D^2) = Nd                          (each of Nd/2 edges contributes 2)
        Tr(D^4) >= Nd^2                        (Cauchy-Schwarz on eigenvalues)
    Equality iff all eigenvalues of D are +-d  (equivalently D^2 = d^2 I).

    The 4x4 torus with FULL pi-flux frustration (all 24 four-cycles = pi,
    including the 8 non-contractible loops) has eigenvalues all +-2 (D^2 = 4 I),
    hence Tr(D^4) = 16 * 16 = 256 = Nd^2, the global minimum.

This supersedes the note's "pi-flux toroidal Tr(D^4) = 384" (that value used the
standard pi-flux gauge, where the 8 non-contractible loops sit at flux 0; the
global holonomy is an independent gauge-invariant d.o.f. and can be flipped to pi,
lowering Tr from 384 to 256).
"""

from __future__ import annotations

import numpy as np
from itertools import combinations


# ---------------- torus + cycles + GF(2) ----------------

def torus_D(Lx, Ly, mode="pi"):
    N = Lx * Ly
    D = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (y % Ly) * Lx + (x % Lx)

    for y in range(Ly):
        for x in range(Lx):
            i = idx(x, y)
            j = idx(x, y + 1)
            D[i, j] = 1.0
            D[j, i] = 1.0
            k = idx(x + 1, y)
            w = (-1.0) ** y if mode == "pi" else 1.0
            D[i, k] = w
            D[k, i] = w
    return D


def edges_from_D(D, thr=1e-9):
    n = D.shape[0]
    return [(i, j) for i in range(n) for j in range(i + 1, n) if abs(D[i, j]) > thr]


def four_cycles(n, edge_set):
    cycles = []
    for S in combinations(range(n), 4):
        a, b, c, d = S
        for cyc in [(a, b, c, d), (a, b, d, c), (a, c, b, d)]:
            if all(frozenset((cyc[i], cyc[(i + 1) % 4])) in edge_set for i in range(4)):
                cycles.append(cyc)
    return cycles


def gf2_solve(A, b):
    A = A.astype(np.uint8).copy()
    b = b.astype(np.uint8).copy()
    nr, nc = A.shape
    row = 0
    piv = []
    for col in range(nc):
        p = None
        for r in range(row, nr):
            if A[r, col]:
                p = r
                break
        if p is None:
            continue
        A[[row, p]] = A[[p, row]]
        b[[row, p]] = b[[p, row]]
        for r in range(nr):
            if r != row and A[r, col]:
                A[r] ^= A[row]
                b[r] ^= b[row]
        piv.append(col)
        row += 1
        if row == nr:
            break
    if np.any(b[row:]):
        return None
    x = np.zeros(nc, dtype=np.uint8)
    for r, col in enumerate(piv):
        if b[r]:
            x[col] = 1
    return x


def fully_frustrated_torus(Lx, Ly):
    """Torus with ALL four-cycles (plaquettes + non-contractible) at flux pi."""
    D = torus_D(Lx, Ly, "pi")
    n = Lx * Ly
    edges = edges_from_D(D)
    fs = set(frozenset(e) for e in edges)
    cycles = four_cycles(n, fs)
    idx = {frozenset(e): k for k, e in enumerate(edges)}
    A = np.zeros((len(cycles), len(edges)), dtype=np.uint8)
    for r, cyc in enumerate(cycles):
        for i in range(4):
            e = frozenset((cyc[i], cyc[(i + 1) % 4]))
            A[r, idx[e]] ^= 1
    signs = gf2_solve(A, np.ones(len(cycles), dtype=np.uint8))
    assert signs is not None, "all-pi not achievable"
    Df = np.zeros((n, n), dtype=complex)
    for k, (i, j) in enumerate(edges):
        s = 1.0 if signs[k] == 0 else -1.0
        Df[i, j] = s
        Df[j, i] = s
    return Df, len(cycles), edges


def main():
    print("### 1) universal bound check on several graphs\n")
    # torus 4x4, standard pi-flux gauge (non-contractible loops at flux 0)
    for label, D in [("torus 4x4 (std pi-gauge)", torus_D(4, 4, "pi")),
                     ("torus 4x4 (no-flux)", torus_D(4, 4, "zero"))]:
        tr2 = float(np.real(np.trace(D @ D)))
        tr4 = float(np.real(np.trace(np.linalg.matrix_power(D, 4))))
        N = D.shape[0]
        print(f"{label:28s}  Tr(D^2)={tr2:6.1f}  Tr(D^4)={tr4:6.1f}  "
              f"Nd^2={N * 16}  bound-holds={tr4 >= tr2 ** 2 / N}")

    # random 4-regular graph, random phases
    import networkx as nx
    rng = np.random.default_rng(1)
    G = nx.random_regular_graph(4, 16, seed=3)
    adj = nx.to_numpy_array(G)
    ph = rng.uniform(-np.pi, np.pi, (16, 16))
    ph = (ph - ph.T) / 2
    D = np.exp(1j * ph) * adj
    tr2 = float(np.real(np.trace(D @ D)))
    tr4 = float(np.real(np.trace(np.linalg.matrix_power(D, 4))))
    print(f"{'random 4-reg + phases':28s}  Tr(D^2)={tr2:6.1f}  Tr(D^4)={tr4:6.1f}  "
          f"Nd^2=256  bound-holds={tr4 >= tr2 ** 2 / 16}")

    print("\n### 2) fully-frustrated torus saturates the bound\n")
    Df, n4, edges = fully_frustrated_torus(4, 4)
    ev = np.linalg.eigvalsh(Df)
    tr2 = float(np.real(np.trace(Df @ Df)))
    tr4 = float(np.real(np.trace(np.linalg.matrix_power(Df, 4))))
    print(f"N=16, n_4cycles={n4}")
    print(f"eigenvalues: {sorted(np.round(ev, 3))}")
    print(f"D^2 == 4 I : {bool(np.allclose(Df @ Df, 4 * np.eye(16)))}")
    print(f"Tr(D^2) = {tr2}  (should be Nd = 64)")
    print(f"Tr(D^4) = {tr4}  (should be Nd^2 = 256)")
    print(f"saturates Cauchy-Schwarz: {np.isclose(tr4, tr2 ** 2 / 16)}")

    print("\n### 3) identity check: Tr(D^4) = modulus + 8 sum cos(Phi)\n")
    for label, D in [("std pi-gauge", torus_D(4, 4, "pi")),
                     ("full frustration", Df)]:
        d_h = np.real(np.diag(D @ D))
        modulus = float(2 * np.sum(d_h ** 2) - np.sum(np.abs(D) ** 4))
        es = set(frozenset(e) for e in edges_from_D(D))
        cyc = four_cycles(16, es)
        holo = sum(8.0 * float(np.real(D[a, b] * D[b, c] * D[c, d] * D[d, a]))
                   for (a, b, c, d) in cyc)
        tr4 = float(np.real(np.trace(np.linalg.matrix_power(D, 4))))
        print(f"{label:16s}  modulus={modulus:6.1f}  holonomy={holo:7.1f}  "
              f"sum={modulus + holo:6.1f}  Tr={tr4:6.1f}  err={abs(tr4 - (modulus + holo)):.1e}")


if __name__ == "__main__":
    main()

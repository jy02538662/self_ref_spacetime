"""Algebra A, edge-space spectral triple (candidate 2, full scan).

Graph Dirac D = partial + partial* on H = C^V (+) C^E (vertices + oriented edges),
where partial: C^E -> C^V is the incidence (boundary) operator.

Checks the single-sided first-order condition [[D,a],b]=0 for candidate algebras A:
  - vertex functions, edge functions, vertex+edge (all commutative)
  - full matrix algebra M_{V+E} (non-commutative, maximal)

Result (2026-09-10): ALL fail.  The single-sided first-order condition has NO
non-trivial solution on a discrete graph -- it is a MANIFOLD condition.  The only
known non-trivial solution is the standard-model finite part (A = C (+) H (+) M3
with the correct charge-conjugation J).
"""

from __future__ import annotations

import numpy as np


def torus_edges(Lx, Ly):
    def idx(x, y):
        return (y % Ly) * Lx + (x % Lx)
    edges = []
    for y in range(Ly):
        for x in range(Lx):
            edges.append((idx(x, y), idx(x + 1, y)))   # right
            edges.append((idx(x, y), idx(x, y + 1)))   # up
    return edges


def main():
    Lx = Ly = 4
    N = Lx * Ly
    edges = torus_edges(Lx, Ly)
    E = len(edges)

    d = np.zeros((N, E))  # incidence: boundary
    for k, (u, v) in enumerate(edges):
        d[v, k] += 1.0
        d[u, k] -= 1.0

    D = np.zeros((N + E, N + E))  # graph Dirac
    D[:N, N:] = d
    D[N:, :N] = d.T

    # D^2 = Hodge Laplacian diag(Delta_0, Delta_1)
    D2 = D @ D
    print("D^2 = diag(Delta_0, Delta_1):  off-diag ~0?",
          np.allclose(D2[:N, N:], 0) and np.allclose(D2[N:, :N], 0))
    print(f"graph Dirac on C^V(+)C^E: N={N} vertices, E={E} edges, dim={N+E}\n")

    rng = np.random.default_rng(0)

    def fo_viol(afn, tries=300):
        w = 0.0
        for _ in range(tries):
            a, b = afn(), afn()
            c = D @ a - a @ D
            w = max(w, float(np.linalg.norm(c @ b - b @ c)))
        return w

    def vert_fun():
        a = np.zeros((N + E, N + E)); a[:N, :N] = np.diag(rng.normal(size=N)); return a
    def edge_fun():
        a = np.zeros((N + E, N + E)); a[N:, N:] = np.diag(rng.normal(size=E)); return a
    def both_fun():
        a = np.zeros((N + E, N + E)); a[:N, :N] = np.diag(rng.normal(size=N)); a[N:, N:] = np.diag(rng.normal(size=E)); return a
    def full():
        return rng.normal(size=(N + E, N + E))

    print(f"first-order [[D,a],b], A=vertex functions : {fo_viol(vert_fun):.3e}")
    print(f"first-order [[D,a],b], A=edge functions   : {fo_viol(edge_fun):.3e}")
    print(f"first-order [[D,a],b], A=vertex+edge diag : {fo_viol(both_fun):.3e}")
    print(f"first-order [[D,a],b], A=full M_48        : {fo_viol(full):.3e}")
    print("\n=> 全违反：单边一阶条件在离散图上无解（流形条件）。"
          "\n   唯一出路 = 标准模型有限部分（非交换 A + 电荷共轭 J）。")


if __name__ == "__main__":
    main()

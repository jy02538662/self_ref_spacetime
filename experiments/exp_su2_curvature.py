"""Step 2B: SU(2) curvature term — algebraic validation.

U(1) curvature term rewards pi-flux (holonomy = e^{i pi} = -1).
SU(2) version: rewards plaquette holonomy W_p = -I (tr = -2), the SU(2) analog of
pi-flux (Z2 center {-I, I}).

Action curvature: F = sum_p (2 + Re tr W_p),  W_p = U_ij U_jk U_kl U_li.
F = 0  <=>  every plaquette W_p = -I  (non-Abelian SU(2) pi-flux).

Two parts:
1. Construct a NON-ABELIAN SU(2) pi-flux on toroidal 3x3: horizontal edges = i*sigma_x,
   vertical edges = i*sigma_y  ->  every plaquette W_p = -I, F = 0.
2. Cold-start optimize F over random SU(2) edges -> does it converge to W_p = -I?
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)


def su2_from_params(a):
    """U = exp(i a.sigma) = cos|a| I + i (sin|a|/|a|) a.sigma."""
    norm = np.linalg.norm(a)
    if norm < 1e-12:
        return I2.copy()
    u = np.cos(norm) * I2 + 1j * (np.sin(norm) / norm) * (a[0] * SX + a[1] * SY + a[2] * SZ)
    return u


def build_edges(N, edge_params):
    """edge_params: dict (i,j) -> 3-vector. Build SU(2) matrix per directed edge (i<j)."""
    U = {}
    for (i, j), a in edge_params.items():
        U[(i, j)] = su2_from_params(a)
        U[(j, i)] = U[(i, j)].conj().T
    return U


def plaquette_holonomies(N, U):
    """W_p for each plaquette (i,j)->(i,j+1)->(i+1,j+1)->(i+1,j)->(i,j)."""
    Ws = []
    for i in range(N):
        for j in range(N):
            jp = (j + 1) % N
            ip = (i + 1) % N
            a = i * N + j       # (i,j)
            b = i * N + jp      # (i,j+1)
            c = ip * N + jp     # (i+1,j+1)
            d = ip * N + j      # (i+1,j)
            W = U[(a, b)] @ U[(b, c)] @ U[(c, d)] @ U[(d, a)]
            Ws.append(W)
    return Ws


def curvature(Ws):
    return float(sum(2.0 + np.real(np.trace(W)) for W in Ws))


def main():
    N = 3
    # ---- Part 1: non-Abelian SU(2) pi-flux ----
    # horizontal edge (i,j)->(i,j+1) = i*sx ; vertical edge (i,j)->(i+1,j) = i*sy
    edge_params = {}
    for i in range(N):
        for j in range(N):
            jp = (j + 1) % N
            ip = (i + 1) % N
            # horizontal edge (i,j)-(i,jp)
            edge_params[tuple(sorted([i * N + j, i * N + jp]))] = np.array([np.pi / 2, 0.0, 0.0])  # i sx
            # vertical edge (i,j)-(ip,j)
            edge_params[tuple(sorted([i * N + j, ip * N + j]))] = np.array([0.0, np.pi / 2, 0.0])  # i sy
    U = build_edges(N, edge_params)
    Ws = plaquette_holonomies(N, U)
    F = curvature(Ws)
    trs = [np.real(np.trace(W)) for W in Ws]
    print(f"Part 1: non-Abelian SU(2) pi-flux (h=i*sx, v=i*sy)")
    print(f"  plaquette tr(W_p) = {np.round(trs, 3).tolist()}")
    print(f"  curvature F = {F:.6f}  (expect 0)")
    print()

    # ---- Part 2: cold-start optimize ----
    print("Part 2: cold-start optimize F over random SU(2) edges")
    # all edges (i<j) on toroidal 3x3
    all_edges = []
    for i in range(N):
        for j in range(N):
            jp = (j + 1) % N
            ip = (i + 1) % N
            e1 = tuple(sorted([i * N + j, i * N + jp]))
            e2 = tuple(sorted([i * N + j, ip * N + j]))
            for e in [e1, e2]:
                if e not in all_edges:
                    all_edges.append(e)
    n_edges = len(all_edges)
    print(f"  n_edges = {n_edges}, n_params = {n_edges * 3}")

    def make_fun(edges):
        def fun(x):
            ep = {edges[k]: x[3 * k:3 * k + 3] for k in range(len(edges))}
            U = build_edges(N, ep)
            Ws = plaquette_holonomies(N, U)
            return curvature(Ws)
        return fun

    fun = make_fun(all_edges)
    for seed in range(3):
        rng = np.random.default_rng(seed)
        x0 = rng.normal(0, 0.5, n_edges * 3)
        res = minimize(fun, x0, method="L-BFGS-B", options={"maxiter": 4000, "ftol": 1e-12, "gtol": 1e-9})
        ep = {all_edges[k]: res.x[3 * k:3 * k + 3] for k in range(len(all_edges))}
        U = build_edges(N, ep)
        Ws = plaquette_holonomies(N, U)
        trs = [np.real(np.trace(W)) for W in Ws]
        print(f"  seed={seed}: F={res.fun:.6f}  tr(W_p)={np.round(trs, 3).tolist()}")


if __name__ == "__main__":
    main()

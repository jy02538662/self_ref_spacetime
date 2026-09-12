"""Find the TRUE minimum of Tr(D^4) over phase (gauge) assignments for a graph.

For a 4-regular graph with |D_ij| = r = 1 on edges:
    Tr(D^4) = Nd(2d-1)  +  8 * sum_{4-cycles C} cos(Phi_C)
so minimizing Tr(D^4) = minimizing sum_C cos(Phi_C) over edge phases.

Two methods:
  (1) Z2 exact: phases in {0, pi} -> solve A eps = 1 over GF(2) (all cycles pi)
      when possible; else enumerate to maximize frustrated count.
  (2) Continuous: scipy optimize over edge phases (catches near-pi compromises).
"""

from __future__ import annotations

import numpy as np
from itertools import combinations
from scipy.optimize import minimize


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


def edge_list_from_D(D, thr=1e-9):
    n = D.shape[0]
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            if abs(D[i, j]) > thr:
                edges.append((i, j))
    return edges


def four_cycles(n, edge_set):
    cycles = []
    for S in combinations(range(n), 4):
        a, b, c, d = S
        for cyc in [(a, b, c, d), (a, b, d, c), (a, c, b, d)]:
            if all(frozenset((cyc[i], cyc[(i + 1) % 4])) in edge_set for i in range(4)):
                cycles.append(cyc)
    return cycles


def build_incidence(cycles, edges):
    """A[C,e] = 1 if edge e in cycle C (over GF(2)).  Returns integer arrays."""
    idx = {frozenset(e): k for k, e in enumerate(edges)}
    A = np.zeros((len(cycles), len(edges)), dtype=np.uint8)
    for r, cyc in enumerate(cycles):
        for i in range(4):
            e = frozenset((cyc[i], cyc[(i + 1) % 4]))
            A[r, idx[e]] ^= 1
    return A


def gf2_solve(A, b):
    """Solve A x = b over GF(2).  Returns x (int array) or None if inconsistent."""
    A = A.astype(np.uint8).copy()
    b = b.astype(np.uint8).copy()
    nrows, ncols = A.shape
    row = 0
    pivots = []
    for col in range(ncols):
        pivot = None
        for r in range(row, nrows):
            if A[r, col]:
                pivot = r
                break
        if pivot is None:
            continue
        A[[row, pivot]] = A[[pivot, row]]
        b[[row, pivot]] = b[[pivot, row]]
        for r in range(nrows):
            if r != row and A[r, col]:
                A[r] ^= A[row]
                b[r] ^= b[row]
        pivots.append(col)
        row += 1
        if row == nrows:
            break
    if np.any(b[row:]):
        return None  # inconsistent
    x = np.zeros(ncols, dtype=np.uint8)
    for r, col in enumerate(pivots):
        if b[r]:
            x[col] = 1
    return x


def apply_phases(n, edges, phase_signs):
    """Build D from edge signs (+1 / -1)."""
    D = np.zeros((n, n), dtype=complex)
    for k, (i, j) in enumerate(edges):
        s = 1.0 if phase_signs[k] == 0 else -1.0
        D[i, j] = s
        D[j, i] = s
    return D


def holonomy_sum(D, cycles):
    s = 0.0
    for (a, b, c, d) in cycles:
        s += float(np.real(D[a, b] * D[b, c] * D[c, d] * D[d, a]))
    return s


def continuous_opt(n, edges, cycles, restarts=40):
    """Minimize sum_C cos(Phi_C) over continuous edge phases."""
    m = len(edges)
    # orientation sign for each cycle/edge: +1 if edge traversed forward in cyc tuple
    edge_pos = {e: k for k, e in enumerate(edges)}
    # cycle -> list of (edge index, sign)
    coeffs = []
    for cyc in cycles:
        row = np.zeros(m)
        for i in range(4):
            u, v = cyc[i], cyc[(i + 1) % 4]
            e = (min(u, v), max(u, v))
            k = edge_pos[e]
            # sign = +1 if traversal u->v matches ascending order, else -1
            sgn = 1 if u < v else -1
            row[k] += sgn
        coeffs.append(row)
    C = np.array(coeffs)  # (nc, m)

    def f(theta):
        ph = C @ theta
        return float(np.sum(np.cos(ph)))

    best = 1e9
    best_theta = None
    rng = np.random.default_rng(1)
    for _ in range(restarts):
        theta0 = rng.uniform(-np.pi, np.pi, m)
        # also seed the pi-flux-like config (signs +/-)
        res = minimize(f, theta0, method="L-BFGS-B",
                       options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-10})
        if res.fun < best:
            best = res.fun
            best_theta = res.x
    return best, best_theta


def analyze_graph(n, edges, label):
    es = {frozenset(e) for e in edges}
    cycles = four_cycles(n, es)
    n4 = len(cycles)

    # --- Z2 exact: can all cycles be pi? ---
    A = build_incidence(cycles, edges)
    ones = np.ones(len(cycles), dtype=np.uint8)
    signs = gf2_solve(A, ones)
    if signs is not None:
        D = apply_phases(n, edges, signs)
        holo = holonomy_sum(D, cycles)
        # verify
        tr4 = float(np.real(np.trace(np.linalg.matrix_power(D, 4))))
        print(f"[{label}] N={n} n_4cycles={n4}  ALL-PI achievable: sum cos = {holo:.4f}, "
              f"Tr(D^4)={tr4:.4f}")
    else:
        print(f"[{label}] N={n} n_4cycles={n4}  all-pi NOT achievable (frustrated)")

    # --- continuous opt ---
    best, _ = continuous_opt(n, edges, cycles)
    print(f"    continuous min sum cos = {best:.4f}  -> Tr(D^4) = {n * 4 * 7 + 8 * best:.4f}")
    return n4, best


def k44_edges():
    """K_{4,4}: left 0..3, right 4..7."""
    edges = []
    for i in range(4):
        for j in range(4, 8):
            edges.append((i, j))
    return edges


def torus_edges(Lx, Ly):
    edges = []
    for y in range(Ly):
        for x in range(Lx):
            i = y * Lx + x
            j = y * Lx + ((x + 1) % Lx)          # horizontal
            k = ((y + 1) % Ly) * Lx + x          # vertical
            edges.append((min(i, j), max(i, j)))
            edges.append((min(i, k), max(i, k)))
    return list(dict.fromkeys(edges))


def two_k44_ladder():
    """Two K_{4,4} copies joined by a 2-edge swap -> connected 4-regular N=16."""
    # copy 1: L=0..3, R=4..7 ; copy 2: L=8..11, R=12..15
    edges = []
    for i in range(4):
        for j in range(4, 8):
            edges.append((i, j))
    for i in range(8, 12):
        for j in range(12, 16):
            edges.append((i, j))
    # swap: remove (0,4) and (8,12); add (0,8) and (4,12)
    edges = [e for e in edges if e not in [(0, 4), (8, 12)]]
    edges.append((0, 8))
    edges.append((4, 12))
    return edges


def main():
    print("### torus 4x4 (N=16, d=4): modulus term = 448")
    analyze_graph(16, torus_edges(4, 4), "4x4 torus")

    print("\n### torus 3x4 (N=12, d=4): modulus term = 336")
    analyze_graph(12, torus_edges(3, 4), "3x4 torus")

    print("\n### K_{4,4} (N=8, d=4): modulus term = 224")
    analyze_graph(8, k44_edges(), "K44")

    print("\n### two K_{4,4} ladder (N=16, d=4): modulus term = 448")
    analyze_graph(16, two_k44_ladder(), "2xK44 ladder")

    print("\n### circulant C16(1,2) (N=16, d=4)")
    n = 16
    edges = []
    for i in range(n):
        for step in (1, 2):
            j = (i + step) % n
            edges.append((min(i, j), max(i, j)))
    edges = list(dict.fromkeys(edges))
    analyze_graph(n, edges, "C16(1,2)")


if __name__ == "__main__":
    main()

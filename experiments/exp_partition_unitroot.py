"""阶段 2a 步骤 1：真正的路径积分选单位根（量子涨落）。

Loop-gas partition function on a toroidal 4-regular graph:
    Z(A) = sum_{FPL configs} d(A)^{#loops},   d(A) = -A^2 - A^{-2}
This is the Kauffman-bracket partition function (quantum dimension d as loop weight),
NOT the weightless approximation sum A^c.

Claim to test: |Z(A)| has peaks/non-analyticity at unit roots A = q^{1/4} (q^{k+2}=1),
i.e. the quantum dimension d(A) = 2cos(pi/(2(k+2))) is "selected" by the sum.

Enumerate all FPL (fully-packed-loop) configs on toroidal 3x3 (9 nodes, 4-degree),
count loops via union-find, sum d(A)^loops, scan A on the unit circle.
"""

from __future__ import annotations

import numpy as np


def toroidal_edges(N):
    """4-regular toroidal NxN. Return list of (node, port) pairs for each edge.
    port: 0=up, 1=right, 2=down, 3=left. Edge connects (i,j,port) to neighbor's opposite port."""
    edges = []
    for i in range(N):
        for j in range(N):
            idx = i * N + j
            # up: (i,j,0) <-> (i-1,j,2)
            up = ((i - 1) % N) * N + j
            edges.append((idx, 0, up, 2))
            # right: (i,j,1) <-> (i,j+1,3)
            r = i * N + ((j + 1) % N)
            edges.append((idx, 1, r, 3))
            # down: (i,j,2) <-> (i+1,j,0)
            dn = ((i + 1) % N) * N + j
            edges.append((idx, 2, dn, 0))
            # left: (i,j,3) <-> (i,j-1,1)
            lf = i * N + ((j - 1) % N)
            edges.append((idx, 3, lf, 1))
    return edges


# 3 pairings of the 4 ports per node: pair (0,1)&(2,3), (0,2)&(1,3), (0,3)&(1,2)
PAIRINGS = [
    [(0, 1), (2, 3)],
    [(0, 2), (1, 3)],
    [(0, 3), (1, 2)],
]


def count_loops(N, config):
    """config: per-node pairing index (0,1,2). Count closed loops via union-find.
    Ports are "half-edges": union each node's paired ports, and each edge's two ends."""
    n_ports = N * N * 4
    parent = list(range(n_ports))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    def pid(node, port):
        return node * 4 + port

    # union node's paired ports
    for node in range(N * N):
        p = PAIRINGS[config[node]]
        union(pid(node, p[0][0]), pid(node, p[0][1]))
        union(pid(node, p[1][0]), pid(node, p[1][1]))
    # union edge ends
    for (a, pa, b, pb) in toroidal_edges(N):
        union(pid(a, pa), pid(b, pb))

    # count connected components = number of closed loops
    roots = set()
    for x in range(n_ports):
        roots.add(find(x))
    return len(roots)


def partition_Z(N, A):
    """Z(A) = sum over FPL configs of d(A)^loops."""
    d = -A ** 2 - A ** (-2)
    n_nodes = N * N
    total = 0.0 + 0j
    # enumerate 3^n configs (ternary)
    for code in range(3 ** n_nodes):
        config = []
        c = code
        for _ in range(n_nodes):
            config.append(c % 3)
            c //= 3
        loops = count_loops(N, config)
        total += d ** loops
    return total


def main():
    N = 3  # toroidal 3x3, 9 nodes, 3^9 = 19683 configs
    print(f"Loop-gas partition function Z(A) on toroidal {N}x{N}, Z = sum d(A)^loops, d=-A^2-A^{-2}")
    print()
    # scan A on the unit circle
    print("scan A = e^{i theta}:  |Z(A)| vs theta")
    thetas = np.linspace(0, np.pi / 2, 13)
    for theta in thetas:
        A = np.exp(1j * theta)
        Z = partition_Z(N, A)
        print(f"  theta={theta/np.pi:.3f}pi  |Z|={abs(Z):.4f}")
    # unit roots: A = q^{1/4}, q = e^{i pi/(k+2)}, i.e. A^4 = e^{i pi/(k+2)}
    # => theta = pi/(4(k+2))
    print()
    print("unit roots A=q^{1/4}, q^{k+2}=1 -> theta_k = pi/(4(k+2)):")
    for k in [1, 2, 3, 4, 5]:
        theta = np.pi / (4 * (k + 2))
        A = np.exp(1j * theta)
        Z = partition_Z(N, A)
        print(f"  k={k}: theta={theta/np.pi:.4f}pi  |Z(A)|={abs(Z):.4f}")


if __name__ == "__main__":
    main()

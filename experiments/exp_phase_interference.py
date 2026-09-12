"""阶段 2a：相位积分 Z = int e^{iS[D]} 选出单位根？（量子干涉，非实权重）

User's core judgment: A is NOT an external parameter in the action; A emerges from
D's phase structure. The action stays pure-trace S[D] = -alpha Tr(D^2) + beta Tr(D^4)
whose phase-sensitive part is the 4-cycle holonomy term (beta * sum_plaquette cos Phi_p,
the "pi-flux" term). Quantumization = phase INTERFERENCE e^{iS} selects unit roots.

Experiment (faithful to user's steps):
1. toroidal 3x3, modulus r=1 fixed, random edge phases theta_ij;
2. per config: S[theta] = beta * sum_p cos(Phi_p)   (pure-trace phase part);
3. weight e^{iS[theta]};
4. after interference, is the plaquette holonomy Phi_p distribution peaked at unit-root
   values (Phi_p = pi, i.e. pi-flux), or uniform?

Diagnostic: weighted distribution of Phi_p (and of A_eff = e^{i Phi_bar}).
"""

from __future__ import annotations

import numpy as np


def toroidal_plaquettes(N):
    """plaquettes of toroidal NxN: each is 4 edges (i,j)->(i,j+1)->(i+1,j+1)->(i+1,j)->(i,j)."""
    plaquettes = []
    for i in range(N):
        for j in range(N):
            a = (i, j)
            b = (i, (j + 1) % N)
            c = ((i + 1) % N, (j + 1) % N)
            d = ((i + 1) % N, j)
            plaquettes.append([a, b, c, d])
    return plaquettes


def edges_of(N):
    """undirected edges of toroidal NxN, with a fixed orientation for the phase."""
    edges = {}
    idx = 0
    for i in range(N):
        for j in range(N):
            # horizontal (i,j)-(i,j+1), vertical (i,j)-(i+1,j)
            e1 = tuple(sorted([(i, j), (i, (j + 1) % N)]))
            e2 = tuple(sorted([(i, j), ((i + 1) % N, j)]))
            for e in [e1, e2]:
                if e not in edges:
                    edges[e] = idx
                    idx += 1
    return edges  # edge -> index


def main():
    N = 3
    plaquettes = toroidal_plaquettes(N)
    edges = edges_of(N)
    n_edges = len(edges)
    # orientation: assign each undirected edge a direction (first -> second) for phase sign
    edge_dir = {e: e for e in edges}  # phase theta_e; Phi_p = sum of signed theta around plaquette

    # build plaquette -> list of (edge, sign)
    def edge_sign(a, b):
        # phase accumulates as theta_{(a,b)} going from a to b
        e = tuple(sorted([a, b]))
        sign = +1 if (a, b) == e else -1
        return e, sign

    plaq_terms = []
    for plaq in plaquettes:
        terms = []
        for k in range(4):
            a, b = plaq[k], plaq[(k + 1) % 4]
            terms.append(edge_sign(a, b))
        plaq_terms.append(terms)

    beta = 1.0
    print(f"Phase interference: Z = int e^{{i beta sum_p cos(Phi_p)}} d theta, toroidal {N}x{N}, "
          f"{n_edges} edges, {len(plaquettes)} plaquettes")
    print()

    # Monte Carlo over phases with reweighting (importance sample uniform, weight e^{iS})
    n_samples = 200000
    rng = np.random.default_rng(0)
    theta = rng.uniform(0, 2 * np.pi, (n_samples, n_edges))

    # compute Phi_p for each sample
    Phi = np.zeros((n_samples, len(plaquettes)))
    for p, terms in enumerate(plaq_terms):
        ph = np.zeros(n_samples)
        for (e, sign) in terms:
            ph += sign * theta[:, edges[e]]
        Phi[:, p] = ph

    S = beta * np.cos(Phi).sum(axis=1)
    weight = np.exp(1j * S)  # complex phase weight

    # weighted distribution of Phi_p (collapse over plaquettes and samples)
    Phi_flat = Phi.ravel()
    w_flat = np.tile(weight, len(plaquettes))

    # histogram of Phi_p weighted by e^{iS}: take real part of weighted histogram
    bins = np.linspace(-np.pi, np.pi, 49)
    hist = np.zeros(len(bins) - 1)
    for k in range(len(bins) - 1):
        mask = (Phi_flat >= bins[k]) & (Phi_flat < bins[k + 1])
        hist[k] = np.real(np.sum(w_flat[mask]))

    print("weighted (e^{iS}) distribution of plaquette holonomy Phi_p:")
    print("  (positive = constructive interference at that Phi)")
    centers = 0.5 * (bins[:-1] + bins[1:])
    for c, h in zip(centers, hist):
        bar = '+' * max(0, int(h / max(1.0, abs(hist).max()) * 30))
        print(f"  Phi={c/np.pi:+.2f}pi  {h:+10.2f}  {bar}")

    # also: average |A_eff|, A_eff = e^{i Phi_bar}
    Phi_bar = Phi.mean(axis=1)
    A_eff = np.exp(1j * Phi_bar / 4.0)
    print()
    print(f"mean A_eff (weighted) = {np.sum(A_eff * weight) / np.sum(weight):.4f}")
    print(f"|mean A_eff| = {abs(np.sum(A_eff * weight) / np.sum(weight)):.4f}")


if __name__ == "__main__":
    main()

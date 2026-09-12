"""Pure Spectral Action -- exact analysis.

1. Verify the exact identity:
      Tr(D^4) = 2*sum_h d_h^2 - sum_{i!=j} r_ij^4 + sum_{distinct} D_ij D_jk D_kl D_li
   where d_h = (D^2)_hh. The 4-cycle term enumerates all 3 Hamilton cycles per 4-vertex subset.

2. Analytic comparison of candidate structures (dimer / ring / complete graph) S values,
   proving "under degree constraint d_h=c, the modulus part prefers fewest edges (dimer), not ring".
"""

from __future__ import annotations

import numpy as np


def decompose_TrD4(D):
    n = D.shape[0]
    r = np.abs(D)
    D2 = D @ D
    tr_D4 = float(np.real(np.trace(D2 @ D2)))
    d = np.real(np.diag(D2))
    term_deg = 2.0 * float(np.sum(d**2))
    term_r4 = float(np.sum(r**4))  # diag=0, so = sum_{i!=j} r^4
    term_quad = 0.0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for dd in range(c + 1, n):
                    for cyc in [(a, b, c, dd), (a, b, dd, c), (a, c, b, dd)]:
                        p, q, s, t = cyc
                        ph = D[p, q] * D[q, s] * D[s, t] * np.conjugate(D[p, t])
                        # 3 undirected Hamilton cycles per 4-set; each cycle = 8 directed
                        # walks (4 starts x 2 orientations) -> 8 Re(ph)
                        term_quad += 8.0 * float(np.real(ph))
    rec = term_deg - term_r4 + term_quad
    return tr_D4, term_deg, term_r4, term_quad, rec, abs(tr_D4 - rec)


def analytic_S(n, alpha, beta, gamma, c, r_vec, phases_kind):
    """r_vec: edge moduli r (upper triangle). phases_kind: all real positive or specified."""
    D = np.zeros((n, n), complex)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = r_vec[idx]
            D[j, i] = r_vec[idx]
            idx += 1
    D2 = D @ D
    tr2 = float(np.real(np.trace(D2)))
    tr4, *_ = decompose_TrD4(D)
    d = np.real(np.diag(D2))
    deg = float(np.sum((d - c) ** 2))
    return -alpha * tr2 + beta * tr4 + gamma * deg, tr2, tr4, d


def candidate_structures(n, c):
    """Build r_vec for each candidate structure (uniform edge length to match d_h=c)."""
    # dimer: N/2 edges r=sqrt(c)
    r_dimer = np.full(n * (n - 1) // 2, 0.0)
    # pairs (0,1),(2,3),...
    idx = 0
    edge_index = {}
    for i in range(n):
        for j in range(i + 1, n):
            edge_index[(i, j)] = idx
            idx += 1
    for k in range(0, n, 2):
        r_dimer[edge_index[(k, k + 1)]] = np.sqrt(c)
    # ring: N edges r=sqrt(c/2)
    r_ring = np.full(n * (n - 1) // 2, 0.0)
    for i in range(n):
        j = (i + 1) % n
        a, b = min(i, j), max(i, j)
        r_ring[edge_index[(a, b)]] = np.sqrt(c / 2)
    # complete: all edges r=sqrt(c/(n-1))
    r_comp = np.full(n * (n - 1) // 2, np.sqrt(c / (n - 1)))
    return r_dimer, r_ring, r_comp


def main():
    print("=" * 72)
    print("1. Verify identity  Tr(D^4) = 2 sum d^2 - sum r^4 + 2 sum_quad r cos(Phi)")
    print("=" * 72)
    rng = np.random.default_rng(0)
    for label, D in []:
        pass
    # (a) real positive ring C6
    n = 6
    D_ring = np.zeros((n, n), complex)
    for i in range(n):
        j = (i + 1) % n
        D_ring[i, j] = 1.0
        D_ring[j, i] = 1.0
    # (b) complete K6 real positive r=1
    D_comp = np.ones((n, n), complex) - np.eye(n, dtype=complex)
    # (c) random Hermitian with phases
    D_rand = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(i + 1, n):
            amp = rng.uniform(0.5, 1.5)
            ph = rng.uniform(0, 2 * np.pi)
            D_rand[i, j] = amp * np.exp(1j * ph)
            D_rand[j, i] = np.conjugate(D_rand[i, j])

    for label, D in [("ring C6 (real)", D_ring), ("K6 (real)", D_comp), ("random Hermitian", D_rand)]:
        tr4, tdeg, tr4_, tquad, rec, err = decompose_TrD4(D)
        print(f"  {label:22s}: TrD4={tr4:.6f}  2sum d^2={tdeg:.4f}  -sum r^4={tr4_:.4f}  quad={tquad:+.4f}  err={err:.2e}")

    print()
    print("=" * 72)
    print("2. Analytic S comparison of candidate structures")
    print("=" * 72)
    for (n, c) in [(6, 1.0), (6, 2.0), (8, 2.0)]:
        alpha, beta, gamma = 2.0, 0.5, 10.0
        rd, rr, rc = candidate_structures(n, c)
        print(f"\n  N={n}, c={c}, alpha={alpha}, beta={beta}, gamma={gamma}")
        for name, rv in [("dimer", rd), ("ring", rr), ("complete", rc)]:
            S, tr2, tr4, d = analytic_S(n, alpha, beta, gamma, c, rv, None)
            d_uniform = np.allclose(d, c, atol=1e-6)
            print(f"    {name:10s}: S={S:9.4f}  TrD2={tr2:6.3f}  TrD4={tr4:8.3f}  d_h~={np.round(d[0],4)}  (uniform deg {d_uniform})")

    print()
    print("=" * 72)
    print("3. Key decomposition numbers for ring vs complete (real positive)")
    print("=" * 72)
    for n, c in [(6, 2.0)]:
        rd, rr, rc = candidate_structures(n, c)
        for name, rv in [("ring", rr), ("complete", rc)]:
            idx = 0
            D = np.zeros((n, n), complex)
            for i in range(n):
                for j in range(i + 1, n):
                    D[i, j] = rv[idx]
                    D[j, i] = rv[idx]
                    idx += 1
            tr4, tdeg, tr4_, tquad, rec, err = decompose_TrD4(D)
            print(f"  {name}: TrD4={tr4:.4f} = {tdeg:.4f}(2sum d^2) {tr4_:+.4f}(-sum r^4) {tquad:+.4f}(quad)")


if __name__ == "__main__":
    main()

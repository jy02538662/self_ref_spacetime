"""Pure Spectral Action cold-start experiment.

User-proposed action:
    S[D] = -alpha*Tr(D^2) + beta*Tr(D^4) + gamma*sum_i [(D^2)_ii - c]^2

This script:
1. Reproduces the user script (N=6, alpha=2, beta=0.5, gamma=10, c=1) and reports what cold start finds.
2. Verifies the exact decomposition identity:
      Tr(D^4) = 2*sum_h d_h^2 - sum_{i!=j} r_ij^4 + sum_{distinct} D_ij D_jk D_kl D_li
   where d_h = (D^2)_hh = sum_j r_hj^2 is the node degree.
   This reveals: the first two terms depend only on moduli (degree), the third on phase (4-cycle flux).
3. Parameter sweeps: whether the modulus part of beta*Tr(D^4) prefers dimer or ring.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]


def vec_to_D(x, n):
    """Real vector -> zero-diagonal Hermitian matrix D."""
    num_edges = n * (n - 1) // 2
    real_parts = x[:num_edges]
    imag_parts = x[num_edges:]
    D = np.zeros((n, n), dtype=complex)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = real_parts[idx] + 1j * imag_parts[idx]
            D[j, i] = real_parts[idx] - 1j * imag_parts[idx]
            idx += 1
    return D


def spectral_action(x, n, alpha, beta, gamma, c):
    D = vec_to_D(x, n)
    D2 = D @ D
    D4 = D2 @ D2
    tr_D2 = np.real(np.trace(D2))
    tr_D4 = np.real(np.trace(D4))
    deg_constraint = 0.0
    for i in range(n):
        deg_constraint += (np.real(D2[i, i]) - c) ** 2
    return -alpha * tr_D2 + beta * tr_D4 + gamma * deg_constraint


def decompose_TrD4(D):
    """Verify identity Tr(D^4) = 2*sum d_h^2 - sum r^4 + sum_{distinct} D_ij D_jk D_kl D_li."""
    n = D.shape[0]
    D2 = D @ D
    tr_D4 = float(np.real(np.trace(D2 @ D2)))

    r = np.abs(D)
    d = np.real(np.diag(D2))  # d_h = sum_j r_hj^2

    term_deg = 2.0 * float(np.sum(d**2))
    term_r4 = float(np.sum(r**4))  # diag=0, so = sum_{i,j} r_ij^4 = sum_{i!=j} r_ij^4
    term_quad = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                for l in range(k + 1, n):
                    # 3 undirected Hamilton cycles per 4-set; each = 8 directed walks -> 8 Re(ph)
                    for cyc in [(i, j, k, l), (i, j, l, k), (i, k, j, l)]:
                        p, q, s, t = cyc
                        ph = D[p, q] * D[q, s] * D[s, t] * np.conjugate(D[p, t])
                        term_quad += 8.0 * float(np.real(ph))
    reconstructed = term_deg - term_r4 + term_quad
    return {
        "TrD4_exact": tr_D4,
        "term_2sum_d2": term_deg,
        "term_sum_r4": term_r4,
        "term_quad_flux": term_quad,
        "reconstructed": reconstructed,
        "identity_error": abs(tr_D4 - reconstructed),
    }


def classify(D, thr=0.2):
    """Threshold the modulus matrix; classify ring / dimer / other."""
    n = D.shape[0]
    r = np.abs(D)
    adj = r > thr
    np.fill_diagonal(adj, False)
    degrees = adj.sum(axis=1)
    # degrees
    deg_counts = {int(k): int(np.sum(degrees == k)) for k in set(degrees)}
    # ring check: each node degree=2 and connected into a single cycle
    is_ring = bool(np.all(degrees == 2))
    is_dimer = bool(np.all(degrees == 1))
    return {
        "degrees": degrees.astype(int).tolist(),
        "deg_counts": deg_counts,
        "is_ring": is_ring,
        "is_dimer": is_dimer,
    }


def run_once(n, alpha, beta, gamma, c, seed, maxiter=4000):
    rng = np.random.default_rng(seed)
    num_edges = n * (n - 1) // 2
    x0 = rng.normal(0, 0.1, num_edges * 2)
    res = minimize(
        spectral_action, x0, args=(n, alpha, beta, gamma, c),
        method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-9},
    )
    D = vec_to_D(res.x, n)
    D2 = D @ D
    d = np.real(np.diag(D2))
    return {
        "seed": seed,
        "success": bool(res.success),
        "S": float(res.fun),
        "TrD2": float(np.real(np.trace(D2))),
        "TrD4": float(np.real(np.trace(D2 @ D2))),
        "degrees_d": [round(float(v), 4) for v in d],
        "classify": classify(D),
        "decomp": decompose_TrD4(D),
        "D_modulus": np.round(np.abs(D), 3).tolist(),
    }


def main():
    # original params (user's)
    params = dict(n=6, alpha=2.0, beta=0.5, gamma=10.0, c=1.0)
    print("=" * 70)
    print(f"original params {params}")
    print("=" * 70)
    rows = []
    for seed in range(4):
        r = run_once(**params, seed=seed)
        rows.append(r)
        cl = r["classify"]
        print(f"\n[seed={seed}] S={r['S']:.5f}  TrD2={r['TrD2']:.4f}  TrD4={r['TrD4']:.4f}")
        print(f"  degrees d_h = {r['degrees_d']}")
        print(f"  classify: ring={cl['is_ring']} dimer={cl['is_dimer']}  deg_counts={cl['deg_counts']}")
        print(f"  modulus matrix:")
        for row in r["D_modulus"]:
            print("    " + " ".join(f"{v:.2f}" for v in row))
        dc = r["decomp"]
        print(f"  TrD4 decomp: 2sum d^2={dc['term_2sum_d2']:.4f}  -sum r^4={dc['term_sum_r4']:.4f}  quad_flux={dc['term_quad_flux']:.4f}  (identity error {dc['identity_error']:.2e})")

    # verify identity (on a graph with 4-cycles: complete graph + random phases)
    print("\n" + "=" * 70)
    print("Verify Tr(D^4) identity (complete graph K6 + random phases)")
    print("=" * 70)
    rng = np.random.default_rng(0)
    n = 6
    D = np.zeros((n, n), complex)
    for i in range(n):
        for j in range(i + 1, n):
            amp = rng.uniform(0.5, 1.5)
            ph = rng.uniform(0, 2 * np.pi)
            D[i, j] = amp * np.exp(1j * ph)
            D[j, i] = np.conjugate(D[i, j])
    dc = decompose_TrD4(D)
    print(f"  TrD4 = {dc['TrD4_exact']:.6f}  (reconstructed {dc['reconstructed']:.6f}, err {dc['identity_error']:.2e})")
    print(f"  2sum d^2={dc['term_2sum_d2']:.4f}  -sum r^4={dc['term_sum_r4']:.4f}  quad_flux(phase term)={dc['term_quad_flux']:.4f}")

    # sweep: c=2 (natural degree of a ring), vary beta/gamma
    print("\n" + "=" * 70)
    print("Sweep c=2 (2 neighbors per node = ring), vary beta/gamma")
    print("=" * 70)
    for (beta, gamma) in [(0.5, 10.0), (1.0, 10.0), (0.5, 50.0)]:
        for seed in range(2):
            r = run_once(n=6, alpha=2.0, beta=beta, gamma=gamma, c=2.0, seed=seed)
            cl = r["classify"]
            print(f"  beta={beta} gamma={gamma} seed={seed}: S={r['S']:.4f} ring={cl['is_ring']} dimer={cl['is_dimer']} deg={cl['deg_counts']}")

    # sweep: N=8 ring
    print("\n" + "=" * 70)
    print("N=8, c=2 (8-node ring)")
    print("=" * 70)
    for seed in range(3):
        r = run_once(n=8, alpha=2.0, beta=0.5, gamma=10.0, c=2.0, seed=seed)
        cl = r["classify"]
        print(f"  seed={seed}: S={r['S']:.4f} ring={cl['is_ring']} dimer={cl['is_dimer']} deg={cl['deg_counts']}")

    out = ROOT / "experiments" / "exp_pure_spectral_last_run.json"
    out.write_text(json.dumps(rows, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

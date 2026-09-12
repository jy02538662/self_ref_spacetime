"""Verify whether curvature term (nu) stably drives toroidal 3x3 + pi-flux.

Print nu=1 full modulus matrix + flux distribution + specdim per seed, count hit rate.
"""

from __future__ import annotations

import json
from pathlib import Path
from collections import Counter

import numpy as np
from scipy.optimize import minimize

from experiments.exp_pure_spectral_4reg import (
    vec_to_D, precompute_cycles, make_action, curvature_term,
    spectral_dimension_from_D,
)

ROOT = Path(__file__).resolve().parents[1]


def degree_distribution(D, thr=0.5):
    n = D.shape[0]
    r = np.abs(D)
    adj = r > thr
    np.fill_diagonal(adj, False)
    return adj.sum(axis=1).astype(int).tolist()


def plaquette_flux_stats(D, thr=0.5):
    n = D.shape[0]
    r = np.abs(D)
    adj = r > thr
    np.fill_diagonal(adj, False)
    cos_list = []
    n_cycles = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for dd in range(c + 1, n):
                    for cyc in [(a, b, c, dd), (a, b, dd, c), (a, c, b, dd)]:
                        p, q, s, t = cyc
                        if adj[p, q] and adj[q, s] and adj[s, t] and adj[t, p]:
                            ph = D[p, q] * D[q, s] * D[s, t] * np.conjugate(D[p, t])
                            denom = abs(D[p, q]) * abs(D[q, s]) * abs(D[s, t]) * abs(D[p, t])
                            if denom > 1e-12:
                                cos_list.append(float(np.real(ph) / denom))
                            n_cycles += 1
    cos = np.asarray(cos_list) if cos_list else np.array([])
    return n_cycles, cos


def main():
    n = 9
    c, m, c4 = 4.0, 4.0, 4.0
    alpha, gamma, delta, nu = 2.0, 10.0, 10.0, 1.0
    cycles = precompute_cycles(n)
    n_seeds = 8
    print(f"=== verify: nu={nu}, N={n}, c={c}, m={m} ===  (threshold 0.5)\n")
    toroidal_hits = 0
    for seed in range(n_seeds):
        rng = np.random.default_rng(seed)
        num_edges = n * (n - 1) // 2
        x0 = rng.normal(0, 0.3, num_edges * 2)
        fun = make_action(n, alpha, gamma, delta, nu, c, c4, cycles)
        res = minimize(fun, x0, method="L-BFGS-B", options={"maxiter": 4000, "ftol": 1e-10, "gtol": 1e-7})
        D = vec_to_D(res.x, n)
        deg = degree_distribution(D)
        n_cycles, cos = plaquette_flux_stats(D)
        pi_frac = float(np.mean(cos < -0.5)) if cos.size else 0.0
        is_4reg = all(dd == 4 for dd in deg)
        is_toroidal = is_4reg and n_cycles == 9
        if is_toroidal:
            toroidal_hits += 1
        print(f"seed={seed}: degree={deg}  4reg={is_4reg}  n_cycles={n_cycles}  pi_frac={pi_frac:.2f}  "
              f"specdim={spectral_dimension_from_D(D):.2f}  {'<<< TOROIDAL' if is_toroidal else ''}")
        if is_toroidal and seed < 3:
            print("  D_modulus:")
            for row in np.round(np.abs(D), 1).tolist():
                print("    " + " ".join(f"{v:.1f}" for v in row))
    print(f"\ntoroidal hits: {toroidal_hits}/{n_seeds}")


if __name__ == "__main__":
    main()

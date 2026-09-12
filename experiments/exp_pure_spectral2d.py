"""Pure Spectral Action -- kill-dimer comparison (step 1).

Action:
    S[D] = -alpha Tr(D^2) + beta Tr(D^4) + gamma sum_i (d_i - c)^2 + mu * sum_{i,j} r_ij^4

where d_i = (D^2)_ii = sum_j r_ij^2 is the weighted degree, r_ij = |D_ij|.

mu term penalizes "energy concentration" (large sum r^4 = dimer, by convexity).
Under degree constraint d_i=c:
    S = const + (mu - beta)*sum r^4 + beta*[4-cycle holonomy term]
- mu < beta  -> prefers concentration (dimer)
- mu > beta  -> prefers dispersion (complete graph direction)
- beta*[4-cycle]  -> prefers pi-flux 2D plaquettes

This script compares mu=0 vs mu>beta to see what degree structure cold start lands on.
Target: N=9, c=4 (2D toroidal grid = 4 neighbors per node).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]


def vec_to_D(x, n):
    num_edges = n * (n - 1) // 2
    rp = x[:num_edges]
    ip = x[num_edges:]
    D = np.zeros((n, n), complex)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = rp[idx] + 1j * ip[idx]
            D[j, i] = rp[idx] - 1j * ip[idx]
            idx += 1
    return D


def action(x, n, alpha, beta, gamma, mu, c):
    D = vec_to_D(x, n)
    D2 = D @ D
    tr2 = float(np.real(np.trace(D2)))
    tr4 = float(np.real(np.trace(D2 @ D2)))
    d = np.real(np.diag(D2))
    deg = float(np.sum((d - c) ** 2))
    r = np.abs(D)
    r4 = float(np.sum(r ** 4))  # sum_{i,j} r_ij^4 (diag 0) = 2*sum_{i<j} r^4
    return -alpha * tr2 + beta * tr4 + gamma * deg + mu * r4


def spectral_dimension_from_D(D, t_lo=0.5, t_hi=50.0, n_pts=40):
    D2 = D @ D
    eig = np.linalg.eigvalsh(D2)
    eig = np.clip(eig, 1e-14, None)
    ts = np.geomspace(t_lo, t_hi, n_pts)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    slope = np.polyfit(np.log(ts), logK, 1)[0]
    return float(-2.0 * slope)


def degree_distribution(D, thr=0.3):
    n = D.shape[0]
    r = np.abs(D)
    adj = r > thr
    np.fill_diagonal(adj, False)
    deg = adj.sum(axis=1).astype(int)
    return deg.tolist()


def plaquette_flux_stats(D, thr=0.3):
    """Find all true 4-cycles (4 distinct vertices + 3 Hamilton cycles whose 4 edges are all in adjacency), compute holonomy cos(Phi)."""
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
    return {
        "n_cycles": n_cycles,
        "n_fluxed": int(len(cos)),
        "pi_frac": float(np.mean(cos < -0.5)) if cos.size else 0.0,
        "mean_cos": float(np.mean(cos)) if cos.size else 0.0,
    }


def run(n, alpha, beta, gamma, mu, c, seed, maxiter=6000):
    rng = np.random.default_rng(seed)
    num_edges = n * (n - 1) // 2
    x0 = rng.normal(0, 0.3, num_edges * 2)
    res = minimize(action, x0, args=(n, alpha, beta, gamma, mu, c),
                   method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-9})
    D = vec_to_D(res.x, n)
    D2 = D @ D
    d = np.real(np.diag(D2))
    r = np.abs(D)
    return {
        "seed": seed,
        "S": float(res.fun),
        "TrD2": float(np.real(np.trace(D2))),
        "TrD4": float(np.real(np.trace(D2 @ D2))),
        "sum_r4": float(np.sum(r ** 4)),
        "d_mean": float(np.mean(d)),
        "d_std": float(np.std(d)),
        "deg_dist": degree_distribution(D),
        "spectral_dim": spectral_dimension_from_D(D),
        "flux": plaquette_flux_stats(D),
    }


def main():
    n, c, alpha, beta, gamma = 9, 4.0, 2.0, 0.5, 10.0
    print(f"=== Pure Spectral Action: kill-dimer comparison  (N={n}, c={c}, alpha={alpha}, beta={beta}, gamma={gamma}) ===")
    print(f"target: 2D toroidal 3x3 (4-regular, 9 plaquettes, pi-flux)\n")
    results = {}
    for mu in [0.0, 1.0, 2.0, 4.0]:
        print(f"--- mu = {mu} (beta={beta}) ---")
        rows = []
        for seed in range(3):
            r = run(n, alpha, beta, gamma, mu, c, seed)
            rows.append(r)
            deg = r["deg_dist"]
            # degree distribution
            from collections import Counter
            degc = dict(Counter(deg))
            print(f"  seed={seed}: S={r['S']:.3f}  d_mean={r['d_mean']:.3f}  d_std={r['d_std']:.3f}  "
                  f"specdim={r['spectral_dim']:.2f}  degree={degc}  pi-flux={r['flux']['pi_frac']:.2f}  "
                  f"(n_cycles={r['flux']['n_cycles']})")
        results[mu] = rows
    out = ROOT / "experiments" / "exp_pure_spectral2d_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

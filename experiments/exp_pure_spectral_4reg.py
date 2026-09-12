"""Pure Spectral Action -- lock 4-regular + curvature to select 2D (step 2, vectorized).

Action:
    S = -alpha Tr(D^2) + gamma sum(d_h - c)^2 + delta sum(sum_j r_hj^4 - c4)^2
        + nu * F_curv
    F_curv = sum_{4-cycles} sum_{3 Hamilton cycles} r_prod (1 + cos Phi)   (rewards pi-flux)

Curvature term uses precomputed cycle indices + numpy vectorization.
"""

from __future__ import annotations

import json
from pathlib import Path
from collections import Counter

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]


def precompute_cycles(n):
    cycles = []
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                for dd in range(c + 1, n):
                    for cyc in [(a, b, c, dd), (a, b, dd, c), (a, c, b, dd)]:
                        cycles.append(cyc)
    return np.asarray(cycles, dtype=int)  # (M, 4)


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


def curvature_term(D, cycles):
    r = np.abs(D)
    p, q, s, t = cycles[:, 0], cycles[:, 1], cycles[:, 2], cycles[:, 3]
    rprod = r[p, q] * r[q, s] * r[s, t] * r[t, p]
    re = np.real(D[p, q] * D[q, s] * D[s, t] * np.conjugate(D[p, t]))
    return float(np.sum(rprod + re))


def make_action(n, alpha, gamma, delta, nu, c, c4, cycles):
    def action(x):
        D = vec_to_D(x, n)
        D2 = D @ D
        tr2 = float(np.real(np.trace(D2)))
        d = np.real(np.diag(D2))
        deg = float(np.sum((d - c) ** 2))
        r = np.abs(D)
        qd = np.sum(r ** 4, axis=1)
        qc = float(np.sum((qd - c4) ** 2))
        curv = curvature_term(D, cycles)
        return -alpha * tr2 + gamma * deg + delta * qc + nu * curv
    return action


def spectral_dimension_from_D(D, t_lo=0.5, t_hi=50.0, n_pts=40):
    # Graph Laplacian L = diag(deg) - |D| (modulus weights, no phase).
    # NOTE: earlier versions used D^2 (adjacency squared) which is WRONG as a Laplacian.
    r = np.abs(D)
    deg = np.sum(r, axis=1)
    L = np.diag(deg) - r
    eig = np.linalg.eigvalsh(L)
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
    return adj.sum(axis=1).astype(int).tolist()


def plaquette_flux_stats(D, thr=0.3):
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
        "pi_frac": float(np.mean(cos < -0.5)) if cos.size else 0.0,
        "mean_cos": float(np.mean(cos)) if cos.size else 0.0,
    }


def run(n, alpha, gamma, delta, nu, c, c4, cycles, seed, maxiter=4000):
    rng = np.random.default_rng(seed)
    num_edges = n * (n - 1) // 2
    x0 = rng.normal(0, 0.3, num_edges * 2)
    fun = make_action(n, alpha, gamma, delta, nu, c, c4, cycles)
    res = minimize(fun, x0, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-7})
    D = vec_to_D(res.x, n)
    D2 = D @ D
    d = np.real(np.diag(D2))
    r = np.abs(D)
    qd = np.sum(r ** 4, axis=1)
    return {
        "seed": seed,
        "S": float(res.fun),
        "d_mean": float(np.mean(d)),
        "d_std": float(np.std(d)),
        "qd_mean": float(np.mean(qd)),
        "qd_std": float(np.std(qd)),
        "deg_dist": dict(Counter(degree_distribution(D))),
        "spectral_dim": spectral_dimension_from_D(D),
        "flux": plaquette_flux_stats(D),
        "D_modulus": np.round(np.abs(D), 2).tolist(),
    }


def main():
    n = 9
    c = 4.0
    m = 4.0
    c4 = c ** 2 / m  # 4
    alpha, gamma, delta = 2.0, 10.0, 10.0
    cycles = precompute_cycles(n)
    print(f"=== Lock 4-regular + curvature  (N={n}, c={c}, m={m}, c4={c4}) ===")
    print(f"alpha={alpha} gamma={gamma} delta={delta}; target: toroidal 3x3 (4-regular, pi-flux)\n")
    results = {}
    for nu in [0.0, 1.0, 2.0]:
        print(f"--- nu = {nu} ---")
        rows = []
        for seed in range(3):
            r = run(n, alpha, gamma, delta, nu, c, c4, cycles, seed)
            rows.append(r)
            print(f"  seed={seed}: S={r['S']:.3f}  d=({r['d_mean']:.2f}+/-{r['d_std']:.3f})  "
                  f"qd=({r['qd_mean']:.2f}+/-{r['qd_std']:.3f})  degree={r['deg_dist']}  "
                  f"specdim={r['spectral_dim']:.2f}  pi-flux={r['flux']['pi_frac']:.2f}  "
                  f"(n_cycles={r['flux']['n_cycles']})")
        results[nu] = rows
        if nu == 2.0:
            print("  D_modulus (seed 0):")
            for row in rows[0]["D_modulus"]:
                print("    " + " ".join(f"{v:.1f}" for v in row))
    out = ROOT / "experiments" / "exp_pure_spectral_4reg_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

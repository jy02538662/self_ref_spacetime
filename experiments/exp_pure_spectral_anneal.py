"""Task 1: pin down whether toroidal 3x3 + pi-flux is the GLOBAL optimum.

Three diagnostics:
1. Construct exact toroidal 3x3 + pi-flux D, compute its S (baseline).
2. Warm-start stability: toroidal + perturbation -> does L-BFGS return to toroidal?
3. Basin-hopping global search: what is the lowest S found, and is it toroidal?

Action: S = -alpha Tr(D^2) + gamma sum(d-c)^2 + delta sum(sum_j r_hj^4 - c4)^2 + nu * F_curv
F_curv = sum_{4-cycles} r_prod (1 + cos Phi)  (rewards pi-flux plaquettes)
"""

from __future__ import annotations

import json
from pathlib import Path
from collections import Counter

import numpy as np
from scipy.optimize import minimize, basinhopping

from experiments.exp_pure_spectral_4reg import (
    vec_to_D, precompute_cycles, make_action, curvature_term,
    spectral_dimension_from_D,
)

ROOT = Path(__file__).resolve().parents[1]


def toroidal_D(n_per_dim=3, pi_flux=True):
    """3x3 torus, nodes (i,j), idx = n_per_dim*i + j. Edges: right (phase 0), down (phase pi*i)."""
    n = n_per_dim ** 2
    D = np.zeros((n, n), complex)
    for i in range(n_per_dim):
        for j in range(n_per_dim):
            idx = n_per_dim * i + j
            # right: (i,j) -> (i, (j+1)%dim), phase 0
            jr = (j + 1) % n_per_dim
            idr = n_per_dim * i + jr
            D[idx, idr] += 1.0
            D[idr, idx] += 1.0
            # down: (i,j) -> ((i+1)%dim, j), phase pi*j (pi-flux per plaquette: holonomy = pi(j+1)-pi*j = pi)
            id_ = (i + 1) % n_per_dim
            idd = n_per_dim * id_ + j
            ph = np.pi * j if pi_flux else 0.0
            D[idx, idd] += np.exp(1j * ph)
            D[idd, idx] += np.exp(-1j * ph)
    return D


def D_to_x(D):
    n = D.shape[0]
    num_edges = n * (n - 1) // 2
    x = np.zeros(num_edges * 2)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            x[idx] = np.real(D[i, j])
            x[num_edges + idx] = np.imag(D[i, j])
            idx += 1
    return x


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


def summarize(D, fun):
    n = D.shape[0]
    D2 = D @ D
    d = np.real(np.diag(D2))
    r = np.abs(D)
    qd = np.sum(r ** 4, axis=1)
    n_cycles, cos = plaquette_flux_stats(D)
    pi_frac = float(np.mean(cos < -0.5)) if cos.size else 0.0
    deg = degree_distribution(D)
    return {
        "S": float(fun(D_to_x(D))),
        "d_mean": float(np.mean(d)),
        "qd_mean": float(np.mean(qd)),
        "deg": deg,
        "n_cycles": n_cycles,
        "pi_frac": pi_frac,
        "specdim": spectral_dimension_from_D(D),
    }


def main():
    n = 9
    c, m, c4 = 4.0, 4.0, 4.0
    alpha, gamma, delta, nu = 2.0, 10.0, 10.0, 1.0
    cycles = precompute_cycles(n)
    fun = make_action(n, alpha, gamma, delta, nu, c, c4, cycles)

    print("=" * 72)
    print("1. Exact toroidal 3x3 + pi-flux  (baseline S)")
    print("=" * 72)
    Dt = toroidal_D(3, pi_flux=True)
    st = summarize(Dt, fun)
    print(f"  toroidal(pi-flux): S={st['S']:.4f}  d_mean={st['d_mean']:.2f}  qd_mean={st['qd_mean']:.2f}  "
          f"deg={Counter(st['deg'])}  n_cycles={st['n_cycles']}  pi_frac={st['pi_frac']:.2f}  specdim={st['specdim']:.2f}")

    Dt0 = toroidal_D(3, pi_flux=False)
    st0 = summarize(Dt0, fun)
    print(f"  toroidal(no-flux): S={st0['S']:.4f}  deg={Counter(st0['deg'])}  n_cycles={st0['n_cycles']}  "
          f"pi_frac={st0['pi_frac']:.2f}  specdim={st0['specdim']:.2f}")

    print("\n" + "=" * 72)
    print("2. Warm-start stability: toroidal + perturbation -> L-BFGS")
    print("=" * 72)
    x_t = D_to_x(Dt)
    n_return = 0
    for seed in range(5):
        rng = np.random.default_rng(seed)
        x0 = x_t + rng.normal(0, 0.3, x_t.size)
        res = minimize(fun, x0, method="L-BFGS-B", options={"maxiter": 4000, "ftol": 1e-10, "gtol": 1e-7})
        D = vec_to_D(res.x, n)
        deg = degree_distribution(D)
        nc, cos = plaquette_flux_stats(D)
        pi_frac = float(np.mean(cos < -0.5)) if cos.size else 0.0
        back = all(dd == 4 for dd in deg) and nc == 9
        if back:
            n_return += 1
        print(f"  seed={seed}: S={res.fun:.4f}  deg={Counter(deg)}  n_cycles={nc}  pi_frac={pi_frac:.2f}  "
              f"{'<- RETURNED' if back else ''}")
    print(f"  warm-start return rate: {n_return}/5")

    print("\n" + "=" * 72)
    print("3. Basin-hopping global search (10 iterations, cold start)")
    print("=" * 72)
    rng = np.random.default_rng(123)
    num_edges = n * (n - 1) // 2
    hits = 0
    best_S = 1e9
    best_desc = ""
    for trial in range(3):
        x0 = rng.normal(0, 0.3, num_edges * 2)
        res = basinhopping(fun, x0, niter=12, T=0.5, stepsize=0.5,
                           minimizer_kwargs={"method": "L-BFGS-B", "options": {"maxiter": 2000, "ftol": 1e-9}})
        D = vec_to_D(res.x, n)
        deg = degree_distribution(D)
        nc, cos = plaquette_flux_stats(D)
        pi_frac = float(np.mean(cos < -0.5)) if cos.size else 0.0
        is_tor = all(dd == 4 for dd in deg) and nc == 9
        if is_tor:
            hits += 1
        if res.fun < best_S:
            best_S = res.fun
            best_desc = f"deg={Counter(deg)} n_cycles={nc} pi_frac={pi_frac:.2f} toroidal={is_tor}"
        print(f"  trial={trial}: S={res.fun:.4f}  deg={Counter(deg)}  n_cycles={nc}  pi_frac={pi_frac:.2f}  "
              f"{'<<< TOROIDAL' if is_tor else ''}")
    print(f"\n  basinhopping toroidal hits: {hits}/3")
    print(f"  best S found: {best_S:.4f}  ({best_desc})")
    print(f"  toroidal baseline S: {st['S']:.4f}")


if __name__ == "__main__":
    main()

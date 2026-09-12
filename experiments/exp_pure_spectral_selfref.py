"""Pure-trace Mexican-hat action (self-reference only), cold start.

User's insight: curvature and pi-flux come from SELF-REFERENCE (D_ij = D_ji^*
gives phases -> loop holonomy = curvature), NOT from a hand-designed curvature
reward term. Tr(D^4) already contains the 4-cycle holonomy (curvature).

So drop the hand-added degree constraint and curvature reward; keep only:
    S[D] = -alpha Tr(D^2) + beta Tr(D^4)      (Mexican hat: break + contain)

Cold-start and see what structure naturally emerges (pi-flux? strings? 3D?).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


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


def action(x, n, alpha, beta):
    D = vec_to_D(x, n)
    D2 = D @ D
    tr2 = float(np.real(np.trace(D2)))
    tr4 = float(np.real(np.trace(D2 @ D2)))
    return -alpha * tr2 + beta * tr4


def spectral_dim(D):
    r = np.abs(D)
    deg = np.sum(r, axis=1)
    L = np.diag(deg) - r
    eig = np.clip(np.linalg.eigvalsh(L), 1e-14, None)
    ts = np.geomspace(0.5, 50.0, 40)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    slope = np.polyfit(np.log(ts), logK, 1)[0]
    return float(-2.0 * slope)


def plaquette_flux(D, thr=0.3):
    n = D.shape[0]
    r = np.abs(D)
    adj = r > thr
    np.fill_diagonal(adj, False)
    cos_list = []
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
    cos = np.asarray(cos_list) if cos_list else np.array([])
    return float(np.mean(cos < -0.5)) if cos.size else 0.0, cos.size


def main():
    n = 9
    print(f"Pure-trace Mexican hat S = -alpha Tr(D^2) + beta Tr(D^4), N={n}, cold start\n")
    for (alpha, beta) in [(2.0, 0.5), (1.0, 0.2), (3.0, 1.0)]:
        print(f"--- alpha={alpha}, beta={beta} ---")
        for seed in range(3):
            rng = np.random.default_rng(seed)
            num_edges = n * (n - 1) // 2
            x0 = rng.normal(0, 0.3, num_edges * 2)
            res = minimize(action, x0, args=(n, alpha, beta), method="L-BFGS-B",
                           options={"maxiter": 6000, "ftol": 1e-12, "gtol": 1e-9})
            D = vec_to_D(res.x, n)
            r = np.abs(D)
            pi_frac, n_cyc = plaquette_flux(D)
            deg = (r > 0.3).sum(axis=1)
            print(f"  seed={seed}: S={res.fun:.3f}  max|z|={r.max():.3f}  specdim={spectral_dim(D):.2f}  "
                  f"pi-flux={pi_frac:.2f} (n_cycles={n_cyc})  degrees={sorted(deg.tolist())}")


if __name__ == "__main__":
    main()

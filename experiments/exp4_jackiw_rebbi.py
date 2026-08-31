"""Exp4: topological charge (winding) as zero modes of D — Jackiw-Rebbi / index theorem.

1D Dirac operator with a domain-wall mass m(x):
    D = kinetic (sigma_x hopping) + m(x) sigma_z   (2 spin components per site).

The winding number of m(x) across the wall (from -m0 to +m0) is a topological charge.
Jackiw-Rebbi: number of zero modes of D equals the winding number.

This is the minimal "topology emerges from the spectrum of D" statement: the same kind
of D that generates distance (Exp1) and signature (Exp3) now also carries topology.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def build_D(N, m0=1.0, t=1.0, wall=None, xi=10.0):
    """1D Dirac operator, dim = 2N. Basis: |x,up> = 2x, |x,down> = 2x+1.

    Kinetic: antisymmetric Hermitian hopping (-i sigma_x d/dx), sigma_x connects up<->down.
    Mass: m(x) sigma_z (up gets +m, down gets -m).
    wall=None -> uniform mass m0 (winding 0); wall=k -> smooth domain wall (winding 1).
    """
    dim = 2 * N
    D = np.zeros((dim, dim), dtype=np.complex128)
    for x in range(N):
        if wall is None:
            m = m0
        else:
            m = m0 * np.tanh((x - wall) / xi)  # smooth domain wall
        D[2 * x, 2 * x] = m
        D[2 * x + 1, 2 * x + 1] = -m
    for x in range(N - 1):
        # antisymmetric Hermitian kinetic = discretized -i sigma_x d/dx
        D[2 * x, 2 * (x + 1) + 1] = -1j * t / 2
        D[2 * (x + 1) + 1, 2 * x] = 1j * t / 2
        D[2 * x + 1, 2 * (x + 1)] = -1j * t / 2
        D[2 * (x + 1), 2 * x + 1] = 1j * t / 2
    return D


def analyze(D, label, winding, tol=1e-8):
    ev = np.linalg.eigvalsh(D)  # D Hermitian -> real eigenvalues
    zero = ev[np.abs(ev) < tol]
    n_zero = int(zero.size)
    n_pos = int(np.sum(ev > tol))
    n_neg = int(np.sum(ev < -tol))
    return {
        "label": label,
        "winding": winding,
        "n_sites": D.shape[0] // 2,
        "eig_min_abs": float(np.min(np.abs(ev))),
        "n_zero": n_zero,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "match": n_zero == winding,
    }


def zero_mode_locations(D, tol=1e-8):
    """Site positions of zero modes (spin-summed |psi|^2 centroid)."""
    ev, evec = np.linalg.eigh(D)
    n = D.shape[0] // 2
    locs = []
    for k in range(ev.size):
        if abs(ev[k]) < tol:
            psi = evec[:, k]
            rho = np.array([abs(psi[2 * x]) ** 2 + abs(psi[2 * x + 1]) ** 2 for x in range(n)])
            rho /= rho.sum()
            locs.append(float(np.sum(rho * np.arange(n))))
    return locs


def main(N=400, m0=1.0, t=1.0, xi=10.0, tol=1e-8):
    print(f"=== Exp4 Jackiw-Rebbi  N={N} m0={m0} t={t} xi={xi} ===")
    rows = [
        analyze(build_D(N, m0=m0, t=t, wall=None), "uniform mass", 0, tol=tol),
        analyze(build_D(N, m0=m0, t=t, wall=N // 2, xi=xi), "domain wall", 1, tol=tol),
    ]
    for r in rows:
        print(
            f"{r['label']:14s}  winding={r['winding']}  "
            f"zero_modes={r['n_zero']:2d}  eig_min_abs={r['eig_min_abs']:.2e}  "
            f"match={r['match']}"
        )
    wall_D = build_D(N, m0=m0, t=t, wall=N // 2, xi=xi)
    locs = zero_mode_locations(wall_D, tol=tol)
    print(f"domain wall zero-mode sites (x): {locs}")
    # naive lattice Dirac doubles fermions (k=0 and k=pi): winding 1 -> 2 zero modes, both at the wall
    n_wall = int(np.sum(np.abs(np.array(locs) - N / 2) < 2 * xi)) if locs else 0
    summary = {
        "uniform_zero": rows[0]["n_zero"],
        "wall_zero_total": rows[1]["n_zero"],
        "zero_mode_sites": locs,
        "fermion_doubling": 2,
        "pass": rows[0]["n_zero"] == 0 and n_wall == 2,
        "criterion": "uniform mass -> 0 zero modes; domain wall -> 2 zero modes at the wall "
        "(= winding 1 x fermion-doubling 2 of naive lattice Dirac)",
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    out = ROOT / "experiments" / "exp4_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp4 Jackiw-Rebbi zero modes / winding")
    p.add_argument("--n", type=int, default=200)
    p.add_argument("--m0", type=float, default=1.0)
    p.add_argument("--t", type=float, default=1.0)
    p.add_argument("--xi", type=float, default=10.0)
    p.add_argument("--tol", type=float, default=1e-8)
    args = p.parse_args()
    main(N=args.n, m0=args.m0, t=args.t, xi=args.xi, tol=args.tol)

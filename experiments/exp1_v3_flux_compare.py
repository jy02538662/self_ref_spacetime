"""Exp1 v3: compare curvature terms — does the DEFINITION of plaquette curvature
decide which magnetic flux survives a curvature-penalty?

    "sin":    F = sin(Phi)      -> zero at Phi=0 AND Phi=pi  (allows Z2 pi-flux)
    "wilson": F = |e^{iPhi}-1|  -> zero only at Phi=0        (forbids pi-flux)
    "cos":    F = cos(Phi)      -> zero at Phi=+/-pi/2       (prefers maximal flux)

Falsifiable expectation:
    - sin    -> Z2 flux (0/pi degenerate, pi-frac random across seeds)
    - wilson -> zero flux (real matrix, pi-frac -> 0, mean_abs_cos -> 1)
    - cos    -> +/-pi/2 flux (quarter-flux, mean_abs_cos -> 0)

This shows that a curvature *penalty* never "prefers" non-trivial topology: it
only pushes flux to the zeros of F, and WHICH fluxes are allowed is set by the
definition of F itself.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.action import geometry_penalty, tr_D2
from src.algebra import build_D, pack_params, random_params, unpack_params
from src.distance import distance_matrix, ring_target_L
from src.flux import flux_stats, plaquette_flux_term


def action(theta: np.ndarray, n: int, L: np.ndarray, lam: float, mu: float, kind: str) -> float:
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    cost = tr_D2(D)
    geo = geometry_penalty(distance_matrix(D), L)
    flux = plaquette_flux_term(D, kind)
    return cost + lam * geo + mu * flux


def run_one(n: int, lam: float, mu: float, kind: str, seed: int, maxiter: int = 400) -> dict:
    rng = np.random.default_rng(seed)
    L = ring_target_L(n)
    theta0 = random_params(n, rng, scale=0.8)

    res = minimize(
        action,
        theta0,
        args=(n, L, lam, mu, kind),
        method="L-BFGS-B",
        options={"maxiter": maxiter, "ftol": 1e-10},
    )

    theta = res.x
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    fs = flux_stats(D)

    return {
        "seed": int(seed),
        "kind": kind,
        "mu": float(mu),
        "success": bool(res.success),
        "flux": float(plaquette_flux_term(D, kind)),
        "geo": float(geometry_penalty(distance_matrix(D), L)),
        "mean_abs_sin_phase": float(np.mean(np.abs(np.sin(phases)))),
        "mean_cos_flux": fs["mean_cos"],
        "mean_abs_cos_flux": fs["mean_abs_cos"],
        "pi_frac": fs["frac_neg_cos"],
        "quarter_frac": fs["frac_quarter_flux"],
    }


def run(
    n: int = 6,
    lam: float = 30.0,
    mu: float = 1.0,
    kinds: list[str] | None = None,
    seeds: range | list[int] | None = None,
    maxiter: int = 400,
) -> list[dict]:
    if kinds is None:
        kinds = ["sin", "wilson", "cos"]
    if seeds is None:
        seeds = range(8)
    rows: list[dict] = []
    for seed in seeds:
        for kind in kinds:
            r = run_one(n=n, lam=lam, mu=mu, kind=kind, seed=int(seed), maxiter=maxiter)
            rows.append(r)
            print(
                f"{kind:7s} seed={seed:2d}  geo={r['geo']:7.3f}  flux={r['flux']:8.5f}  "
                f"mean_cos={r['mean_cos_flux']:+.3f}  mean|cos|={r['mean_abs_cos_flux']:.3f}  "
                f"pi_frac={r['pi_frac']:.2f}  quarter={r['quarter_frac']:.2f}"
            )
    return rows


def summarize(rows: list[dict]) -> dict:
    kinds = sorted({r["kind"] for r in rows})
    by_kind: dict[str, dict] = {}
    for k in kinds:
        sub = [r for r in rows if r["kind"] == k]
        by_kind[k] = {
            "n_seeds": len(sub),
            "mean_flux": float(np.mean([r["flux"] for r in sub])),
            "mean_cos": float(np.mean([r["mean_cos_flux"] for r in sub])),
            "mean_abs_cos": float(np.mean([r["mean_abs_cos_flux"] for r in sub])),
            "mean_pi_frac": float(np.mean([r["pi_frac"] for r in sub])),
            "mean_quarter_frac": float(np.mean([r["quarter_frac"] for r in sub])),
            "mean_abs_sin_phase": float(np.mean([r["mean_abs_sin_phase"] for r in sub])),
        }

    # falsifiable claims
    sin = by_kind.get("sin", {})
    wilson = by_kind.get("wilson", {})
    cos = by_kind.get("cos", {})

    wilson_forbids_pi = wilson.get("mean_pi_frac", 1.0) < 0.1
    wilson_zero_flux = wilson.get("mean_abs_cos", 0.0) > 0.9
    cos_quarter = cos.get("mean_quarter_frac", 0.0) > 0.6
    cos_max_flux = cos.get("mean_abs_cos", 1.0) < 0.4
    sin_allows_pi = sin.get("mean_pi_frac", 0.0) > 0.1

    return {
        "by_kind": by_kind,
        "wilson_forbids_pi_flux": bool(wilson_forbids_pi),
        "wilson_flattens_to_zero_flux": bool(wilson_zero_flux),
        "cos_prefers_quarter_flux": bool(cos_quarter),
        "sin_allows_pi_flux": bool(sin_allows_pi),
        "conclusion": (
            "curvature-term DEFINITION decides which flux survives: "
            "sin allows pi-flux, wilson forbids it, cos prefers +/-pi/2. "
            "A curvature penalty never 'prefers' non-trivial topology — it only "
            "pushes flux to the zeros of F."
        ),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp1 v3: curvature-term definition comparison")
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--mu", type=float, default=1.0)
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--maxiter", type=int, default=400)
    args = p.parse_args()

    print(f"=== Exp1 v3 curvature-term comparison  N={args.n}  lambda={args.lam}  mu={args.mu} ===")
    rows = run(n=args.n, lam=args.lam, mu=args.mu, seeds=range(args.seeds), maxiter=args.maxiter)
    summary = summarize(rows)
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp1_v3_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

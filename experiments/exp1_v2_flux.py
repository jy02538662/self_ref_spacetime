"""Exp1 v2: does a closed-loop magnetic-flux term let phases order?

S[D] = Tr(D^2) + lam * sum_{i<j}(d_ij - L_ij)^2 + mu * sum_{i<j<k} sin^2(Phi_ijk)

The flux term depends ONLY on phases (z/|z| = e^{i phi}), so it decouples from
the moduli.  Question: from random phases, does mu > 0 drive phases to a trivial
(zero-flux / real) state, or a non-trivial (pi-flux / textured) state?

Honest expectation: mu > 0 is a *penalty* on flux, so it flattens phases to
zero flux (real matrix) — NOT a non-trivial texture.  This is the concrete,
falsifiable version of "penalty != competition": a flux penalty alone cannot
spontaneously grow topology.
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
from src.flux import flux_stats, flux_term


def action_with_flux(theta: np.ndarray, n: int, L: np.ndarray, lam: float, mu: float) -> float:
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    cost = tr_D2(D)
    geo = geometry_penalty(distance_matrix(D), L)
    flux = flux_term(D)
    return cost + lam * geo + mu * flux


def run_one(n: int, lam: float, mu: float, seed: int, maxiter: int = 400, a: float = 1.0) -> dict:
    rng = np.random.default_rng(seed)
    L = ring_target_L(n, a=a)
    theta0 = random_params(n, rng, scale=0.8)

    moduli0, phases0 = unpack_params(theta0, n)
    D0 = build_D(moduli0, phases0, n)
    flux0 = flux_term(D0)

    res = minimize(
        action_with_flux,
        theta0,
        args=(n, L, lam, mu),
        method="L-BFGS-B",
        options={"maxiter": maxiter, "ftol": 1e-10},
    )

    theta = res.x
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    cost = tr_D2(D)
    geo = geometry_penalty(distance_matrix(D), L)
    flux = flux_term(D)
    fs = flux_stats(D)

    return {
        "seed": int(seed),
        "mu": float(mu),
        "success": bool(res.success),
        "S": float(cost + lam * geo + mu * flux),
        "TrD2": float(cost),
        "geo": float(geo),
        "flux": float(flux),
        "flux0": float(flux0),
        "mean_abs_sin_phase": float(np.mean(np.abs(np.sin(phases)))),
        "mean_cos_flux": fs["mean_cos"],
        "mean_sin2_flux": fs["mean_sin2"],
        "frac_pi_flux": fs["frac_neg_cos"],
        "n_triangles": fs["n_triangles"],
    }


def run(
    n: int = 6,
    lam: float = 30.0,
    mus: list[float] | None = None,
    seeds: range | list[int] | None = None,
    maxiter: int = 400,
) -> list[dict]:
    if mus is None:
        mus = [0.0, 0.1, 1.0, 10.0]
    if seeds is None:
        seeds = range(8)
    rows: list[dict] = []
    for seed in seeds:
        for mu in mus:
            r = run_one(n=n, lam=lam, mu=mu, seed=int(seed), maxiter=maxiter)
            rows.append(r)
            print(
                f"mu={mu:5.1f} seed={seed:2d}  S={r['S']:8.3f}  geo={r['geo']:7.3f}  "
                f"flux={r['flux']:7.4f} (init {r['flux0']:6.3f})  "
                f"mean|sin phi|={r['mean_abs_sin_phase']:.3f}  "
                f"mean_cos={r['mean_cos_flux']:+.3f}  pi-frac={r['frac_pi_flux']:.2f}"
            )
    return rows


def summarize(rows: list[dict]) -> dict:
    mus = sorted({r["mu"] for r in rows})
    by_mu: dict[float, dict] = {}
    for mu in mus:
        sub = [r for r in rows if r["mu"] == mu]
        by_mu[mu] = {
            "n_seeds": len(sub),
            "mean_flux": float(np.mean([r["flux"] for r in sub])),
            "mean_geo": float(np.mean([r["geo"] for r in sub])),
            "mean_abs_sin_phase": float(np.mean([r["mean_abs_sin_phase"] for r in sub])),
            "mean_cos_flux": float(np.mean([r["mean_cos_flux"] for r in sub])),
            "mean_pi_frac": float(np.mean([r["frac_pi_flux"] for r in sub])),
        }

    # Key claim to test: does the flux penalty flatten phases to zero flux?
    base = by_mu.get(0.0, {})
    strong = by_mu.get(10.0, {})
    flux_suppressed = (
        base.get("mean_flux", 1.0) > 0.1
        and strong.get("mean_flux", 1.0) < 0.1 * base.get("mean_flux", 1.0)
    )
    phases_flattened = strong.get("mean_abs_sin_phase", 1.0) < 0.1

    return {
        "by_mu": by_mu,
        "flux_penalty_suppresses_flux": bool(flux_suppressed),
        "phases_flattened_to_real": bool(phases_flattened),
        "conclusion": (
            "flux penalty drives phases to zero flux (real matrix) — "
            "penalty != competition; no non-trivial texture emerges."
        )
        if flux_suppressed and phases_flattened
        else "flux penalty did NOT simply flatten phases — inspect by_mu.",
        "criterion": (
            "mu=0 flux is substantial (>0.1) and mu=10 flux drops to <10% of it, "
            "and mean |sin phi| < 0.1 (real matrix)."
        ),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp1 v2: closed-loop magnetic flux term")
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--maxiter", type=int, default=400)
    args = p.parse_args()

    print(f"=== Exp1 v2 flux  N={args.n}  lambda={args.lam} ===")
    rows = run(n=args.n, lam=args.lam, seeds=range(args.seeds), maxiter=args.maxiter)
    summary = summarize(rows)
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp1_v2_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

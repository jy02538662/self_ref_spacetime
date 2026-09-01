"""Exp1 v4: Metropolis-Hastings thermal sampling — does finite temperature break
the Z2 (0/pi) flux degeneracy?

Samples exp(-beta S) over phases (moduli fixed to a ring).  This is CLASSICAL
thermal sampling (beta = 1/T_classical), the first step off the T=0 L-BFGS
corpse — NOT quantum (no e^{iS}, no coherence, no path integral).

Falsifiable question: at finite beta, is 0-flux vs pi-flux occupancy 1:1
(degeneracy intact, because sin^2(Phi) is symmetric around 0 and pi) or biased
(degeneracy broken by some hidden entropy asymmetry)?
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.algebra import build_D, pack_params, unpack_params
from src.flux import flux_stats, plaquette_flux_term
from src.optimize import minimize_action


def run_mcmc(
    n: int,
    lam: float,
    mu: float,
    kind: str,
    beta: float,
    seed: int,
    warm_seed: int,
    n_steps: int = 30000,
    step_scale: float = 0.3,
    record_every: int = 100,
) -> tuple[list[dict], float]:
    rng = np.random.default_rng(seed)
    # Warm-start: find ring moduli via L-BFGS (mu=0), keep only the moduli.
    warm = minimize_action(n=n, lam=lam, lam4=0.0, seed=warm_seed, maxiter=300, phase_mode="free")
    moduli, _ = unpack_params(warm.theta, n)
    m = n * (n - 1) // 2
    phases = rng.uniform(-np.pi, np.pi, size=m)

    def flux_of(ph: np.ndarray) -> float:
        return mu * plaquette_flux_term(build_D(moduli, ph, n), kind)

    S_cur = flux_of(phases)
    n_accept = 0
    records: list[dict] = []
    for step in range(n_steps):
        ph_new = phases + rng.normal(0.0, step_scale, size=m)
        S_new = flux_of(ph_new)
        dS = S_new - S_cur
        if dS <= 0.0 or rng.random() < np.exp(-beta * dS):
            phases = ph_new
            S_cur = S_new
            n_accept += 1
        if (step + 1) % record_every == 0:
            D = build_D(moduli, phases, n)
            fs = flux_stats(D)
            records.append(
                {
                    "step": step + 1,
                    "flux": float(plaquette_flux_term(D, kind)),
                    "pi_frac": fs["frac_neg_cos"],
                    "mean_cos": fs["mean_cos"],
                    "quarter_frac": fs["frac_quarter_flux"],
                }
            )
    accept_rate = n_accept / n_steps
    return records, accept_rate


def summarize_beta(rows: list[dict]) -> dict:
    betas = sorted({r["beta"] for r in rows})
    by_beta: dict[float, dict] = {}
    for b in betas:
        sub = [r for r in rows if r["beta"] == b]
        means = np.array([r["pi_frac_mean"] for r in sub])
        by_beta[b] = {
            "n_seeds": len(sub),
            "pi_frac_grand_mean": float(np.mean(means)),
            "pi_frac_spread_across_seeds": float(np.std(means)),
            "mean_flux": float(np.mean([r["flux_mean"] for r in sub])),
            "mean_accept_rate": float(np.mean([r["accept_rate"] for r in sub])),
        }
    # degeneracy intact <=> pi_frac grand mean ~ 0.5 (no bias toward 0 or pi)
    return {
        "by_beta": by_beta,
        "degeneracy_intact": all(
            abs(v["pi_frac_grand_mean"] - 0.5) < 0.15 for v in by_beta.values()
        ),
        "criterion": "pi_frac grand mean stays ~0.5 across betas => 0/pi flux occupancy symmetric",
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp1 v4: thermal sampling, Z2 degeneracy")
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--mu", type=float, default=1.0)
    p.add_argument("--kind", type=str, default="sin")
    p.add_argument("--seeds", type=int, default=4)
    p.add_argument("--steps", type=int, default=30000)
    args = p.parse_args()

    betas = [0.01, 0.1, 1.0, 10.0]
    print(f"=== Exp1 v4 thermal sampling  N={args.n}  mu={args.mu}  kind={args.kind} ===")
    rows: list[dict] = []
    for beta in betas:
        for seed in range(args.seeds):
            records, accept = run_mcmc(
                n=args.n,
                lam=args.lam,
                mu=args.mu,
                kind=args.kind,
                beta=beta,
                seed=seed,
                warm_seed=seed,
                n_steps=args.steps,
            )
            tail = records[len(records) // 2 :]
            pi = np.array([r["pi_frac"] for r in tail])
            flux = np.array([r["flux"] for r in tail])
            row = {
                "beta": beta,
                "seed": seed,
                "accept_rate": accept,
                "pi_frac_mean": float(pi.mean()),
                "pi_frac_std": float(pi.std()),
                "flux_mean": float(flux.mean()),
                "pi_frac_min": float(pi.min()),
                "pi_frac_max": float(pi.max()),
            }
            rows.append(row)
            print(
                f"beta={beta:5.2f} seed={seed}  accept={accept:.2f}  "
                f"pi_frac={row['pi_frac_mean']:.3f} +/- {row['pi_frac_std']:.3f}  "
                f"[{row['pi_frac_min']:.2f},{row['pi_frac_max']:.2f}]  flux={row['flux_mean']:.3f}"
            )
    summary = summarize_beta(rows)
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp1_v4_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

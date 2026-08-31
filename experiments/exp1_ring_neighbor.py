"""Exp1: ring target L_ij — does S-minimization yield a neighbor cycle?"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.optimize import minimize_action


def run(
    n: int = 6,
    lam: float = 30.0,
    seeds: range | list[int] | None = None,
    maxiter: int = 500,
) -> list[dict]:
    if seeds is None:
        seeds = range(8)
    rows: list[dict] = []
    for seed in seeds:
        r = minimize_action(n=n, lam=lam, lam4=0.0, seed=int(seed), maxiter=maxiter)
        row = {
            "n": n,
            "lam": lam,
            "seed": int(seed),
            "opt_success": r.success,
            "S": r.S,
            "TrD2": r.stats["TrD2"],
            "geo": r.stats["geo"],
            "mean_neighbor_|z|": r.stats["mean_neighbor_|z|"],
            "mean_other_|z|": r.stats["mean_other_|z|"],
            "strength_ratio": r.stats["strength_ratio_other_over_neigh"],
            "neighbor_hit_rate": r.stats["neighbor_hit_rate"],
            "non_neighbor_strong_rate": r.stats["non_neighbor_strong_rate"],
            "message": r.message,
        }
        rows.append(row)
        print(
            f"seed={seed:2d}  S={r.S:.4f}  geo={r.stats['geo']:.4f}  "
            f"hit={r.stats['neighbor_hit_rate']:.2f}  "
            f"|z|_other/|z|_ring={r.stats['strength_ratio_other_over_neigh']:.3f}  "
            f"ok={r.success}"
        )
    return rows


def seed_is_ring(row: dict, ratio_pass: float = 0.2, hit_pass: float = 0.99, geo_pass: float = 0.1) -> bool:
    return (
        row["strength_ratio"] <= ratio_pass
        and row["neighbor_hit_rate"] >= hit_pass
        and row["geo"] <= geo_pass
    )


def summarize(rows: list[dict], fraction_pass: float = 0.25) -> dict:
    """
    Landscape is multimodal: report fraction of seeds that reach a clean ring,
    plus best-S seed stats. Pass if enough seeds find the ring basin.
    """
    flags = [seed_is_ring(r) for r in rows]
    frac = float(np.mean(flags)) if rows else 0.0
    best = min(rows, key=lambda r: r["S"])
    return {
        "n_seeds": len(rows),
        "fraction_ring_basins": frac,
        "n_ring_basins": int(sum(flags)),
        "best_S": best["S"],
        "best_geo": best["geo"],
        "best_strength_ratio": best["strength_ratio"],
        "best_neighbor_hit_rate": best["neighbor_hit_rate"],
        "pass": frac >= fraction_pass and seed_is_ring(best),
        "criterion": (
            f"fraction of seeds reaching ring >= {fraction_pass}, "
            "and best-S seed is a clean ring (geo<=0.1, hit>=0.99, |z| ratio<=0.2)"
        ),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp1 ring neighbor emergence")
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--seeds", type=int, default=8)
    p.add_argument("--maxiter", type=int, default=500)
    args = p.parse_args()

    print(f"=== Exp1 ring neighbor  N={args.n}  lambda={args.lam} ===")
    rows = run(n=args.n, lam=args.lam, seeds=range(args.seeds), maxiter=args.maxiter)
    summary = summarize(rows)
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp1_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

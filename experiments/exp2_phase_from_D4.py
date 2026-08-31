"""Exp2: does Tr(D^4) favor imaginary ring phases?

Protocol A — fixed exact ring moduli, compare Tr(D^4) across phase patterns.
Protocol B — warm-start from Exp1 ring, optimize with lam4>0:
              phases free vs phases clamped real.

Honest note: for Hermitian D, spectrum of D^2 is always >=0 (squares of eig(D)).
Lorentzian signature is NOT visible as negative eig(D^2); Exp2 only tests the
weaker claim that Tr(D^4) prefers imaginary phases on the ring.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.action import action_from_D, tr_D4
from src.algebra import build_D, pack_params, ring_neighbor_mask, unpack_params
from src.distance import ring_target_L
from src.optimize import exact_ring_theta, minimize_action


def _ring_phases_pattern(n: int, pattern: str) -> np.ndarray:
    m = n * (n - 1) // 2
    phases = np.zeros(m)
    mask = ring_neighbor_mask(n)
    ring_idx = np.where(mask)[0]
    if pattern == "all_real":
        pass
    elif pattern == "one_imag":
        phases[ring_idx[0]] = np.pi / 2
    elif pattern == "two_imag":
        phases[ring_idx[0]] = np.pi / 2
        if len(ring_idx) > 1:
            phases[ring_idx[1]] = np.pi / 2
    elif pattern == "all_imag":
        phases[mask] = np.pi / 2
    elif pattern == "alternating":
        for k, idx in enumerate(ring_idx):
            phases[idx] = (np.pi / 2) if (k % 2 == 0) else 0.0
    else:
        raise ValueError(pattern)
    return phases


def protocol_A(n: int = 6, a: float = 1.0, lam: float = 30.0, lam4: float = 1.0) -> list[dict]:
    L = ring_target_L(n, a=a)
    rows = []
    for pattern in ("all_real", "one_imag", "two_imag", "all_imag", "alternating"):
        phases = _ring_phases_pattern(n, pattern)
        theta = exact_ring_theta(n, a=a, phases=phases)
        moduli, ph = unpack_params(theta, n)
        D = build_D(moduli, ph, n)
        parts = action_from_D(D, L, lam=lam, lam4=lam4)
        rows.append(
            {
                "pattern": pattern,
                "TrD2": parts["TrD2"],
                "TrD4": parts["TrD4"],
                "geo": parts["geo"],
                "S": parts["S"],
            }
        )
    # Sort by TrD4 ascending (claimed preferred)
    rows_sorted = sorted(rows, key=lambda r: r["TrD4"])
    print("=== Protocol A: fixed ring moduli, vary phases ===")
    for r in rows_sorted:
        print(
            f"  {r['pattern']:12s}  TrD4={r['TrD4']:.6f}  TrD2={r['TrD2']:.6f}  "
            f"geo={r['geo']:.4f}  S={r['S']:.4f}"
        )
    return rows


def protocol_B(
    n: int = 6,
    lam: float = 30.0,
    lam4: float = 1.0,
    seeds: range | list[int] | None = None,
    maxiter: int = 400,
) -> list[dict]:
    if seeds is None:
        seeds = range(6)
    rows: list[dict] = []
    print(f"=== Protocol B: warm-start ring, lam4={lam4}, free vs real ===")
    for seed in seeds:
        # First find a ring with lam4=0 (reuse Exp1 basin via multi-start)
        warm = minimize_action(n=n, lam=lam, lam4=0.0, seed=int(seed), maxiter=maxiter, phase_mode="free")
        # Keep moduli from warm, reset phases for the two branches
        moduli, _ = unpack_params(warm.theta, n)
        theta_real0 = pack_params(moduli, np.zeros_like(moduli))
        rng = np.random.default_rng(1000 + int(seed))
        theta_free0 = pack_params(moduli, rng.uniform(-np.pi, np.pi, size=moduli.shape))

        real = minimize_action(
            n=n,
            lam=lam,
            lam4=lam4,
            seed=int(seed),
            maxiter=maxiter,
            phase_mode="phases_only",
            theta0=theta_real0,
        )
        # real phases_only with zeros stays zero; also allow moduli+real clamp
        real_mod = minimize_action(
            n=n,
            lam=lam,
            lam4=lam4,
            seed=int(seed),
            maxiter=maxiter,
            phase_mode="real",
            theta0=theta_real0,
        )
        free = minimize_action(
            n=n,
            lam=lam,
            lam4=lam4,
            seed=int(seed),
            maxiter=maxiter,
            phase_mode="phases_only",
            theta0=theta_free0,
        )
        row = {
            "seed": int(seed),
            "warm_S": warm.S,
            "warm_geo": warm.stats["geo"],
            "real_S": real_mod.S,
            "real_TrD4": real_mod.stats["TrD4"],
            "real_sin2_ring": real_mod.stats["mean_sin2_phase_ring"],
            "free_S": free.S,
            "free_TrD4": free.stats["TrD4"],
            "free_sin2_ring": free.stats["mean_sin2_phase_ring"],
            "free_imag_frac": free.stats["ring_imag_like_fraction"],
            "delta_S_free_minus_real": free.S - real_mod.S,
            "delta_TrD4_free_minus_real": free.stats["TrD4"] - real_mod.stats["TrD4"],
        }
        rows.append(row)
        print(
            f"seed={seed:2d}  warm_geo={warm.stats['geo']:.3f}  "
            f"real_S={real_mod.S:.4f} TrD4={real_mod.stats['TrD4']:.4f}  "
            f"free_S={free.S:.4f} TrD4={free.stats['TrD4']:.4f} "
            f"sin2={free.stats['mean_sin2_phase_ring']:.3f}  "
            f"dS={row['delta_S_free_minus_real']:+.4f}"
        )
    return rows


def summarize(A: list[dict], B: list[dict]) -> dict:
    best_A = min(A, key=lambda r: r["TrD4"])
    real_A = next(r for r in A if r["pattern"] == "all_real")
    # Claim: some imaginary pattern has strictly lower TrD4 than all_real
    claim_A = best_A["TrD4"] < real_A["TrD4"] - 1e-9 and best_A["pattern"] != "all_real"

    # Among seeds with decent warm ring (geo small), does free beat real?
    good = [r for r in B if r["warm_geo"] < 0.5]
    if not good:
        good = B
    frac_free_lower_S = float(np.mean([r["delta_S_free_minus_real"] < -1e-6 for r in good]))
    frac_free_lower_TrD4 = float(np.mean([r["delta_TrD4_free_minus_real"] < -1e-6 for r in good]))
    mean_free_sin2 = float(np.mean([r["free_sin2_ring"] for r in good]))

    # Pass Exp2 only if BOTH: imag patterns can lower TrD4, AND optimization prefers free phases
    passed = claim_A and frac_free_lower_TrD4 >= 0.5 and mean_free_sin2 > 0.2
    return {
        "protocol_A_lowest_TrD4_pattern": best_A["pattern"],
        "protocol_A_imag_beats_real": claim_A,
        "protocol_A_best_TrD4": best_A["TrD4"],
        "protocol_A_real_TrD4": real_A["TrD4"],
        "protocol_B_n_good_seeds": len(good),
        "protocol_B_frac_free_lower_S": frac_free_lower_S,
        "protocol_B_frac_free_lower_TrD4": frac_free_lower_TrD4,
        "protocol_B_mean_free_sin2_ring": mean_free_sin2,
        "pass": passed,
        "note": (
            "Hermitian D => eig(D^2)>=0 always; Exp2 does not test Lorentz signature, "
            "only whether Tr(D^4) favors imaginary phases."
        ),
        "criterion": (
            "A: some imag pattern has TrD4 < all_real; "
            "B: >=50% good seeds free TrD4 < real and mean sin2_ring > 0.2"
        ),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp2 Tr(D^4) / imaginary phase")
    p.add_argument("--n", type=int, default=6)
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--lam4", type=float, default=1.0)
    p.add_argument("--seeds", type=int, default=6)
    p.add_argument("--maxiter", type=int, default=400)
    args = p.parse_args()

    A = protocol_A(n=args.n, lam=args.lam, lam4=args.lam4)
    B = protocol_B(n=args.n, lam=args.lam, lam4=args.lam4, seeds=range(args.seeds), maxiter=args.maxiter)
    summary = summarize(A, B)
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp2_last_run.json"
    out.write_text(json.dumps({"A": A, "B": B, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

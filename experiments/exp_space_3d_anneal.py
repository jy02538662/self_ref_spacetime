"""空间涌现第一步：D 自发复现三维（退火 + L-BFGS 精修）。

纯随机 + L-BFGS 会卡在局部极小（hit 0.62）。这里用 Metropolis 退火
（温度从高到低）帮它跳出局部极小，再 L-BFGS 精修，看能否自发收敛到三维
（hit_top6 -> 1.0）。

诚实标注：仍是"复现三维"（给三维目标距离），不是"自发选维"。
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
from src.distance import distance_matrix


def cubic_coords(n_per_dim: int) -> list[tuple[int, int, int]]:
    return [(x, y, z) for x in range(n_per_dim) for y in range(n_per_dim) for z in range(n_per_dim)]


def manhattan(a, b) -> int:
    return abs(a[0] - b[0]) + abs(a[1] - b[1]) + abs(a[2] - b[2])


def cubic_target_L(n_per_dim: int) -> np.ndarray:
    coords = cubic_coords(n_per_dim)
    n = len(coords)
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            L[i, j] = manhattan(coords[i], coords[j])
    return L


def anneal_and_refine(
    n_per_dim: int,
    lam: float,
    seed: int,
    T0: float = 150.0,
    T_final: float = 0.001,
    n_steps: int = 15000,
    maxiter: int = 600,
) -> dict:
    n = n_per_dim**3
    L = cubic_target_L(n_per_dim)
    coords = cubic_coords(n_per_dim)

    neighbors = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if manhattan(coords[i], coords[j]) == 1:
                neighbors[i].append(j)

    rng = np.random.default_rng(seed)
    theta = random_params(n, rng, scale=0.8)

    def S_of(theta):
        moduli, phases = unpack_params(theta, n)
        D = build_D(moduli, phases, n)
        return tr_D2(D) + lam * geometry_penalty(distance_matrix(D), L)

    # 退火
    S_cur = S_of(theta)
    for step in range(n_steps):
        T = T0 * (T_final / T0) ** (step / n_steps)
        step_scale = 0.025 * np.sqrt(T) + 1e-5
        theta_new = theta + rng.normal(0.0, step_scale, size=len(theta))
        S_new = S_of(theta_new)
        dS = S_new - S_cur
        if dS <= 0.0 or rng.random() < np.exp(-dS / T):
            theta = theta_new
            S_cur = S_new

    # L-BFGS 精修
    res = minimize(S_of, theta, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-10})
    theta = res.x
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    absD = np.abs(D)

    hit_list = []
    for i in range(n):
        row = absD[i].copy()
        row[i] = -1.0
        top6 = np.argsort(row)[-6:]
        true_neigh = set(neighbors[i])
        n_hit = sum(1 for j in top6 if j in true_neigh)
        hit_list.append(n_hit / len(true_neigh) if true_neigh else 0.0)

    hit = float(np.mean(hit_list))
    geo = float(geometry_penalty(distance_matrix(D), L))
    return {"seed": seed, "lam": lam, "hit_top6": hit, "geo": geo, "success": bool(res.success)}


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="空间涌现第一步：退火+精修 自发复现三维")
    p.add_argument("--n-per-dim", type=int, default=3)
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--seeds", type=int, default=5)
    args = p.parse_args()

    print(f"=== 退火+精修 自发复现三维  N={args.n_per_dim}^3={args.n_per_dim**3}  lam={args.lam} ===")
    rows = []
    for seed in range(args.seeds):
        r = anneal_and_refine(args.n_per_dim, args.lam, seed)
        rows.append(r)
        print(f"seed={seed}  hit_top6={r['hit_top6']:.3f}  geo={r['geo']:.3f}")

    mean_hit = float(np.mean([r["hit_top6"] for r in rows]))
    best = max(rows, key=lambda r: r["hit_top6"])
    summary = {
        "mean_hit_top6": mean_hit,
        "best_hit_top6": best["hit_top6"],
        "best_geo": best["geo"],
        "n_clean": int(sum(1 for r in rows if r["hit_top6"] > 0.9)),
        "n_seeds": len(rows),
        "conclusion": (
            "退火让纯随机也自发收敛到三维" if mean_hit > 0.9 else f"退火后 mean_hit={mean_hit:.3f}，还需调"
        ),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp_space_3d_anneal_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

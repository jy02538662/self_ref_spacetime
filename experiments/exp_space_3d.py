"""空间涌现第一步：D 复现三维近邻图（Exp1 的三维推广）——推进干净版。

N=27（3³），内部态有 6 个近邻、边界效应小。多种子 + 扫 λ。
判据：每个态取 |z| 最大的 6 条边，看命中三维近邻（曼哈顿距离=1）的比例。

诚实标注：验证"D 复现三维几何"，不是"D 自发选维"。
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


def run(n_per_dim: int, lam: float, seed: int, maxiter: int = 300, warm_start: bool = False) -> dict:
    n = n_per_dim**3
    L = cubic_target_L(n_per_dim)
    coords = cubic_coords(n_per_dim)

    # 三维近邻集合（每个态的真实三维邻居）
    neighbors = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if manhattan(coords[i], coords[j]) == 1:
                neighbors[i].append(j)

    rng = np.random.default_rng(seed)
    if warm_start:
        # 暖启动：三维近邻边 |z|=1，非近邻 |z|=0.05，phases 全 0
        m = n * (n - 1) // 2
        moduli = np.full(m, 0.05)
        idx = 0
        for i in range(n):
            for j in range(i + 1, n):
                if manhattan(coords[i], coords[j]) == 1:
                    moduli[idx] = 1.0
                idx += 1
        theta0 = pack_params(moduli, np.zeros(m))
    else:
        theta0 = random_params(n, rng, scale=0.8)

    def fun(theta):
        moduli, phases = unpack_params(theta, n)
        D = build_D(moduli, phases, n)
        return tr_D2(D) + lam * geometry_penalty(distance_matrix(D), L)

    res = minimize(fun, theta0, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-10})
    theta = res.x
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    absD = np.abs(D)

    # top-6 判据：每个态取 |z| 最大的 6 条边，看命中真实三维近邻的比例
    hit_list = []
    for i in range(n):
        # 排除对角，取 |z| 最大的 6 条边（内部态有 6 个近邻）
        row = absD[i].copy()
        row[i] = -1.0  # 排除自边
        top6 = np.argsort(row)[-6:]
        true_neigh = set(neighbors[i])
        n_hit = sum(1 for j in top6 if j in true_neigh)
        hit_list.append(n_hit / len(true_neigh) if true_neigh else 0.0)

    hit = float(np.mean(hit_list))
    geo = float(geometry_penalty(distance_matrix(D), L))

    # 每个态的实际近邻数（|z| > 0.3*max 的强边数）
    thr = 0.3 * np.max(absD)
    degree = (absD > thr).sum(axis=1)  # 对角 |z|=0 不是强边，无需 -1
    mean_degree = float(np.mean(degree))

    return {
        "n": n,
        "n_per_dim": n_per_dim,
        "lam": lam,
        "seed": seed,
        "success": bool(res.success),
        "geo": geo,
        "hit_top6": hit,
        "mean_degree": mean_degree,
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="空间涌现第一步：D 复现三维近邻（干净版）")
    p.add_argument("--n-per-dim", type=int, default=3)
    p.add_argument("--lams", type=str, default="30")
    p.add_argument("--seeds", type=int, default=4)
    p.add_argument("--warm-start", action="store_true")
    args = p.parse_args()

    lams = [float(x) for x in args.lams.split(",")]
    print(f"=== 空间涌现第一步  N={args.n_per_dim}^3={args.n_per_dim**3}  lams={lams}  warm={args.warm_start} ===")
    rows = []
    for lam in lams:
        for seed in range(args.seeds):
            r = run(args.n_per_dim, lam, seed, warm_start=args.warm_start)
            rows.append(r)
            print(
                f"lam={lam:5.1f} seed={seed}  geo={r['geo']:8.3f}  hit_top6={r['hit_top6']:.3f}  "
                f"mean_degree={r['mean_degree']:.2f}"
            )

    by_lam = {}
    for lam in lams:
        sub = [r for r in rows if r["lam"] == lam]
        by_lam[lam] = {
            "mean_hit_top6": float(np.mean([r["hit_top6"] for r in sub])),
            "mean_geo": float(np.mean([r["geo"] for r in sub])),
            "mean_degree": float(np.mean([r["mean_degree"] for r in sub])),
        }

    best = max(by_lam.items(), key=lambda kv: kv[1]["mean_hit_top6"])
    summary = {
        "by_lam": by_lam,
        "best_lam": best[0],
        "best_hit_top6": best[1]["mean_hit_top6"],
        "conclusion": (
            f"最优 lam={best[0]}，hit_top6={best[1]['mean_hit_top6']:.3f}"
            + ("（干净，接近 1.0）" if best[1]["mean_hit_top6"] > 0.9 else "（待继续调）")
        ),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp_space_3d_clean_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

"""空间涌现 · 档 2：维度对比实验（2D vs 3D vs 4D）—— 只优化模长版。

关键：作用量 S = Tr(D²) + λ·geo_penalty(d[D], L) 只依赖 |z|（模长），不依赖相位
（Tr(D²)=Σ|z|²，距离用 1/|z|）。所以相位固定为 0，只优化 m = n(n-1)/2 个模长，参数砍半、提速。

问题：D 复现三维近邻结构，是「三维特殊」，还是「任何维度都能同样复现」？
判据 hit：每个态取 |z| 最大的 2d 条边，看命中真实 d 维近邻（曼哈顿距离=1）的比例。
baseline = 2d/(N-1)；excess = (hit - baseline)/(1 - baseline) 归一化到 [0,1]。

诚实边界：这是「复现维度对比」，不是「自发选维」。只能回答「三维对 D 是否（相对）特殊」。
"""

from __future__ import annotations

import json
import sys
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.action import geometry_penalty, tr_D2
from src.algebra import build_D
from src.distance import distance_matrix


def lattice_coords(d: int, n_per_dim: int) -> list[tuple[int, ...]]:
    return list(product(range(n_per_dim), repeat=d))


def manhattan(a, b) -> int:
    return sum(abs(x - y) for x, y in zip(a, b))


def target_L(coords) -> np.ndarray:
    n = len(coords)
    L = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            L[i, j] = manhattan(coords[i], coords[j])
    return L


def neighbor_lists(coords) -> list[list[int]]:
    n = len(coords)
    neighbors = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if manhattan(coords[i], coords[j]) == 1:
                neighbors[i].append(j)
    return neighbors


def warm_moduli(coords, n: int) -> np.ndarray:
    """暖启动模长：d 维近邻边 =1，非近邻 =0.05。"""
    m = n * (n - 1) // 2
    moduli = np.full(m, 0.05)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if manhattan(coords[i], coords[j]) == 1:
                moduli[idx] = 1.0
            idx += 1
    return moduli


def hit_top_k(absD: np.ndarray, neighbors, coord: int) -> float:
    n = absD.shape[0]
    hit_list = []
    for i in range(n):
        row = absD[i].copy()
        row[i] = -1.0
        top = np.argsort(row)[-coord:]
        true = set(neighbors[i])
        nh = sum(1 for j in top if j in true)
        hit_list.append(nh / len(true) if true else 0.0)
    return float(np.mean(hit_list))


def run(d: int, n_per_dim: int, lam: float, seed: int, warm: bool, maxiter: int = 150) -> dict:
    coords = lattice_coords(d, n_per_dim)
    n = len(coords)
    m = n * (n - 1) // 2
    L = target_L(coords)
    neighbors = neighbor_lists(coords)
    coord = 2 * d

    rng = np.random.default_rng(seed)
    mod0 = warm_moduli(coords, n) if warm else rng.uniform(0.05, 0.8, size=m)

    def fun(moduli):
        D = build_D(np.abs(moduli), np.zeros(m), n)  # 相位固定 0
        return tr_D2(D) + lam * geometry_penalty(distance_matrix(D), L)

    res = minimize(
        fun, mod0, method="L-BFGS-B", bounds=[(0.0, None)] * m,
        options={"maxiter": maxiter, "ftol": 1e-10},
    )
    D = build_D(np.abs(res.x), np.zeros(m), n)
    absD = np.abs(D)
    hit = hit_top_k(absD, neighbors, coord)
    geo = float(geometry_penalty(distance_matrix(D), L))
    baseline = coord / (n - 1)
    excess = (hit - baseline) / (1.0 - baseline) if baseline < 1.0 else 0.0
    return {
        "d": d, "n_per_dim": n_per_dim, "N": n, "coord": coord,
        "lam": lam, "seed": seed, "warm": warm,
        "hit": hit, "geo": geo, "baseline": baseline, "excess": excess,
        "success": bool(res.success), "niter": int(res.nit),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="档 2：维度对比（2D/3D/4D 复现难易，模长优化版）")
    p.add_argument("--lam", type=float, default=30.0)
    p.add_argument("--seeds", type=int, default=2)
    p.add_argument("--maxiter", type=int, default=150)
    args = p.parse_args()

    # (dim, n_per_dim)：参数≤351（N≤27）。4D n=2 是 16 角点超立方（无内部态），仅作对照。
    configs = [(2, 5), (3, 3), (4, 2)]  # N = 25, 27, 16

    print(f"=== 档 2 维度对比（模长版） lam={args.lam} seeds={args.seeds} maxiter={args.maxiter} ===", flush=True)
    rows = []
    for d, npd in configs:
        for warm in (True, False):
            for seed in range(args.seeds):
                r = run(d, npd, args.lam, seed, warm=warm, maxiter=args.maxiter)
                rows.append(r)
                tag = "warm" if warm else "cold"
                print(
                    f"d={d} N={r['N']:3d} {tag:4s} seed={seed}  "
                    f"hit={r['hit']:.3f}  excess={r['excess']:.3f}  geo={r['geo']:.2f}  niter={r['niter']}",
                    flush=True,
                )

    by_key = {}
    for r in rows:
        by_key.setdefault((r["d"], r["warm"]), []).append(r)

    summary = {}
    for (d, warm), sub in sorted(by_key.items()):
        tag = "warm" if warm else "cold"
        summary[f"d={d}_{tag}"] = {
            "N": sub[0]["N"], "coord": sub[0]["coord"],
            "mean_hit": float(np.mean([x["hit"] for x in sub])),
            "mean_excess": float(np.mean([x["excess"] for x in sub])),
            "baseline": sub[0]["baseline"],
        }

    warm_hits = {d: summary[f"d={d}_warm"]["mean_hit"] for d in (2, 3, 4)}
    cold_hits = {d: summary[f"d={d}_cold"]["mean_hit"] for d in (2, 3, 4)}
    cold_excess = {d: summary[f"d={d}_cold"]["mean_excess"] for d in (2, 3, 4)}
    conclusion = {
        "warm_stability": {
            d: ("稳定(hit≈1)" if warm_hits[d] > 0.9 else f"不稳(hit={warm_hits[d]:.2f})") for d in (2, 3, 4)
        },
        "cold_findability": {d: round(cold_hits[d], 3) for d in (2, 3, 4)},
        "cold_excess": {d: round(cold_excess[d], 3) for d in (2, 3, 4)},
        "cold_rank": sorted((2, 3, 4), key=lambda d: -cold_excess[d]),
        "note": "4D n=2 为 16 角点超立方（无内部态），坐标=8 但每角点真实近邻仅 4，仅作弱对照",
    }
    print("--- conclusion ---", flush=True)
    print(json.dumps(conclusion, indent=2, ensure_ascii=False), flush=True)

    out = ROOT / "experiments" / "exp_space_dim_compare_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary, "conclusion": conclusion}, indent=2), encoding="utf-8")
    print(f"wrote {out}", flush=True)

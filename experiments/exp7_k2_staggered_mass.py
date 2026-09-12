"""Exp7 后续：交错质量分离 —— 把 Dirac Kramers 对（T²=−1）从格点平移简并里分出来。

背景（v4 第十节逼三维）：π 磁通偶数尺寸的谱是「自对偶」（存在反幺 T、T²=−1），但
k=2 的 Kramers 对（二重简并）被格点平移的 4 重简并（±2、0）污染。

目标：证明 T²=−1 是 Dirac 点的内禀性质，而不是格点平移的偶然产物。

方案（用户建议，且物理正确）：在 π 磁通背景上加**交错（Semenoff）质量 m·(−1)^(i+j)**。
  - Semenoff 质量是**时间反演不变**质量（区别于 Haldane 质量破 T），所以保持 T²=−1，
    → Kramers 对（二重简并）应**存活**；
  - 交错质量是 sublattice 调制，破坏格点平移 → 抬掉平移造成的 4 重简并。

解析预期（π 磁通能带 E(k)=±2√(sin²kx+sin²ky)，加质量后 E=±√(m²+4(sin²kx+sin²ky))）：
  - 零模（Dirac 点，E=0）→ ±m
  - ±2 → ±√(m²+4)
  - ±2√2（带边）→ ±√(m²+8)

判据：扫 m，Dirac 相关能级按上述公式移动；若 Kramers 对（精确二重简并）在 m>0 后
仍然精确二重、且 4 重格点简并抬成 2+2（或单态），则「纯 k=2」是真的 → T²=−1 内禀。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
from experiments.exp_pure_spectral_anneal import toroidal_D


def staggered_H(n_per_dim: int, m: float) -> np.ndarray:
    """π 磁通 toroidal + 交错质量 m·(−1)^(i+j)。返回厄米 H。"""
    D = toroidal_D(n_per_dim, pi_flux=True)
    n = n_per_dim
    signs = np.array([(-1) ** (i + j) for i in range(n) for j in range(n)], dtype=float)
    return D + m * np.diag(signs)


def multiplicities(evals: np.ndarray, tol: float) -> np.ndarray:
    evals = np.sort(evals)
    mults = []
    start = 0
    for i in range(1, len(evals)):
        if evals[i] - evals[i - 1] > tol:
            mults.append(i - start)
            start = i
    mults.append(len(evals) - start)
    return np.asarray(mults, dtype=int)


def main():
    print("=" * 80)
    print("交错质量分离：π 磁通 + m·(−1)^(i+j)，看 Kramers 对能否从格点简并里分出来")
    print("判据：m>0 后 Kramers 对保持精确二重简并，4 重格点简并被抬掉")
    print("=" * 80)

    results = {}
    ms = [0.0, 0.2, 0.5, 1.0, 1.5, 2.0]

    for n_per_dim in (4, 6, 8):
        n = n_per_dim ** 2
        print(f"\n{'#' * 80}\n# N = {n}  ({n_per_dim}×{n_per_dim} π-flux)\n{'#' * 80}")
        for m in ms:
            H = staggered_H(n_per_dim, m)
            evals, evecs = np.linalg.eigh(H)
            evals = np.sort(evals)
            scale = float(evals.max() - evals.min())
            mults = multiplicities(evals, 1e-8 * scale if scale > 1e-12 else 1e-8)
            all_even = bool(np.all(mults % 2 == 0))
            # 统计 k=2（能级占比）
            total = n
            k2 = float(2 * np.sum(mults == 2)) / total
            k4 = float(4 * np.sum(mults == 4)) / total
            k1 = float(np.sum(mults == 1)) / total
            results[f"N{n}_m{m}"] = {"evals": np.round(evals, 4).tolist(),
                                     "mults": mults.tolist(), "all_even": all_even,
                                     "k1": round(k1, 3), "k2": round(k2, 3), "k4": round(k4, 3)}
            print(f"  m={m:.1f}: mults={mults.tolist()}  自对偶={all_even}  k1={k1:.2f} k2={k2:.2f} k4={k4:.2f}")
            print(f"         evals={np.round(evals, 3).tolist()}")

    # 单独验证解析式：N=16，取 0 模 / ±2 / ±2√2 三条轨迹 vs √(m²+4 sin²)
    print("\n" + "=" * 80)
    print("轨迹验证（N=16）：0→±m、±2→±√(m²+4)、±2√2→±√(m²+8)")
    print("=" * 80)
    for m in [0.0, 0.5, 1.0, 2.0]:
        H = staggered_H(4, m)
        ev = np.sort(np.linalg.eigvalsh(H))
        # 找最接近 0 / 2 / 2.828 的能级（正半轴）
        pos = ev[ev > 0]
        closest = {}
        for target, name in [(0.0, "0模"), (2.0, "±2"), (2.828, "±2√2")]:
            closest[name] = pos[np.argmin(np.abs(pos - target))]
        pred0 = m
        pred2 = np.sqrt(m ** 2 + 4)
        predE = np.sqrt(m ** 2 + 8)
        print(f"  m={m}: 实测 0模={closest['0模']:.4f} (预测±{pred0:.4f})  "
              f"±2={closest['±2']:.4f} (预测±{pred2:.4f})  "
              f"±2√2={closest['±2√2']:.4f} (预测±{predE:.4f})")

    out = ROOT / "experiments" / "exp7_k2_staggered_mass_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

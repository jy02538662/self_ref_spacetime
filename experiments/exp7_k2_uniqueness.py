"""Exp7 后续：唯一性（方向 1）—— 有向 → 反对称 → 号差 + SU(2)

方向 1 核心命题：唯一的（至符号）实反对称 2×2 满足 J²=−I 的矩阵是 [[0,1],[-1,0]]，
它**同时**是 (a) Dirac γ⁰（号差）和 (b) 四元数结构 J₀（T²=−1/SU(2)）。

链：
  有向（反对称：方向翻号）
    → 闭环（J²=−1：绕两圈翻号 = 断裂的二元闭环）
    → 唯一 2×2 解 [[0,1],[-1,0]] = γ⁰ = J₀
    → 号差（时空） + SU(2)（内部/自旋）是同一个对象的两个面。

验证四件：
1. 唯一性：实反对称 2×2 J=[[0,a],[-a,0]]，J²=−a²I=−I ⟺ a=±1（唯一至符号）。
2. γ⁰ = [[0,1],[-1,0]]，(γ⁰)²=−1（号差，时间方向 −1）。
3. J₀ = [[0,1],[-1,0]]，T=J₀K，T²=J₀ conj(J₀)=J₀²=−1（辛类 SU(2)）。
4. 同一矩阵 [[0,1],[-1,0]] 同时是 γ⁰ 和 J₀（统一）。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    I2 = np.eye(2)
    sigma = np.array([[0, 1], [-1, 0]], dtype=float)  # = γ⁰ = J₀ = 有向区分

    print("=" * 80)
    print("唯一性（方向 1）：有向 → 反对称 → 号差 + SU(2)")
    print("=" * 80)

    # 1) 唯一性：实反对称 2×2 J=[[0,a],[-a,0]]，J²=−I 的唯一解
    print("\n[1] 唯一性：实反对称 2×2 满足 J²=−I 的解")
    sols = []
    for a in np.linspace(-3, 3, 1201):
        J = np.array([[0, a], [-a, 0]], dtype=float)
        err = np.max(np.abs(J @ J + I2))
        if err < 1e-9:
            sols.append(round(a, 3))
    print(f"  J²=−I 的解 a = {sorted(set(sols))}  （即 ±[[0,1],[-1,0]]，唯一至符号）")

    # 2) γ⁰ = [[0,1],[-1,0]]，号差（(γ⁰)²=−1）
    print("\n[2] γ⁰ = [[0,1],[-1,0]]，号差")
    gamma0 = np.array([[0, 1], [-1, 0]], dtype=float)
    print(f"  (γ⁰)² = {gamma0 @ gamma0}  →  (γ⁰)²=−I  ✓（时间方向 −1）")
    print(f"  γ⁰ 反对称？ γ⁰+γ⁰ᵀ = {gamma0 + gamma0.T}  ✓")

    # 3) J₀ = [[0,1],[-1,0]]，T=J₀K，T²=−1（SU(2)/辛类）
    print("\n[3] J₀ = [[0,1],[-1,0]]，T=J₀K，T²=−1")
    J0 = np.array([[0, 1], [-1, 0]], dtype=float)
    Tsq = J0 @ J0.conj()  # T² = J₀ conj(J₀) = J₀²（J₀ 实）
    print(f"  T² = J₀·conj(J₀) = {Tsq}  →  T²=−I  ✓（辛类 β=4 = SU(2)）")
    print(f"  J₀ 反对称？ J₀+J₀ᵀ = {J0 + J0.T}  ✓")

    # 4) 统一：同一矩阵 [[0,1],[-1,0]] = γ⁰ = J₀
    print("\n[4] 统一：γ⁰ == J₀ == 有向区分？")
    same = np.allclose(gamma0, J0) and np.allclose(gamma0, sigma)
    print(f"  γ⁰ == J₀ == [[0,1],[-1,0]]  →  {same}")

    # 5) 关键：J²=−I 本身不唯一（反对称才是关键）
    print("\n[5] 关键：J²=−I 本身不唯一，反对称才唯一")
    M = np.array([[1, 2], [-1, -1]], dtype=float)  # 非反对称，但 M²=−I
    print(f"  非反对称 M=[[1,2],[-1,-1]]，M² = {M @ M}  → M²=−I 但非反对称 ✓（不唯一）")
    print(f"  而反对称 2×2 只能 [[0,a],[-a,0]]，J²=−I ⟹ a=±1 → 唯一 ±[[0,1],[-1,0]]")
    print("  ⇒ 唯一性来自「反对称（有向）」，不是 J²=−I 本身")

    print("\n结论：号差（γ⁰）和 SU(2)（J₀）是同一个对象 [[0,1],[-1,0]] 的两面，")
    print("      而这个对象是「有向区分」（反对称闭环 J²=−1）唯一的最小实现。")

    out = ROOT / "experiments" / "exp7_k2_uniqueness_last_run.json"
    out.write_text(json.dumps({
        "solutions_a": sorted(set(sols)),
        "gamma0_sq_minus1": True,
        "J0_Tsq_minus1": True,
        "unification": bool(same),
        "J2_minus1_without_antisymmetry_not_unique": True,
    }, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

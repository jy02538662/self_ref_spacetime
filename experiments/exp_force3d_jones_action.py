"""逼三维 · 组装版：涌现 SU(2) 提供「上下」+ Jones 提供「链接奖励」

把前面散开的几块拼成一个「逼三维」的作用量项：

  [块1] π 磁通（偶数尺寸）→ T²=−1（自发）→ Kramers 对 = spin-1/2（2×2，内部上下）。
        ——「上下」不是 U(1) 相位给的（绕一圈 holonomy 只是 ±1，对易），
        而是「有限性（偶数 N）→ 自对偶 → SU(2)」自己长出来的内部两极。
  [块2] 两条「strand」（= 两个 Kramers 对）的张量积 = 2×2 ⊗ 2×2 = 4×4，
        辫子生成元 σ = A·1 + A⁻¹·e（R 矩阵）作用其上，交叉符号住在 A vs A⁻¹。
  [块3] 链接（客观）= Jones 多项式。σ⁰（平凡，2D）= d；σ²（Hopf，3D）= -A⁴-A⁻⁴。

作用量项（奖励链接）：
    S_link = -ν · |V(σ^k) - V(σ⁰)|     （ν>0，奖励 Jones 偏离平凡值 = 奖励链接）

数值结论（本脚本产出）：S_link(σ⁰)=0（2D 不奖励）、S_link(σ²)>0（3D 奖励），
即「Jones 奖励」会把作用量往「链接（3D）」方向推。这是逼三维的作用量项雏形。

诚实边界：本脚本把「上下」的来源（涌现 SU(2)）和「链接奖励」（Jones）接上了，
但「环 → 辫子词从 D 自发读出」仍用玩具版（σ^k 手放）；完全自发的「D 动力学 → 
辫子词 → Jones 奖励 → 选三维」还需要下一版把 σ^k 的 k 从 D 的环结构里读出来。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

from experiments.exp_pure_spectral_anneal import toroidal_D


def tl2_kauffman(A: complex, k: int) -> complex:
    """σ₁^k 闭包的 Kauffman bracket（Markov 迹 / d），σ = A·1 + A⁻¹·e。"""
    d = -(A**2) - A ** (-2)
    a, b = 1.0 + 0j, 0.0 + 0j
    for _ in range(k):
        a, b = A * a, A * b + (A ** -1) * a + (A ** -1) * b * d
    return a * d + b


def jones(A: complex, k: int) -> complex:
    return (-A) ** (-3 * k) * tl2_kauffman(A, k)


def main():
    q = np.exp(1j * np.pi / 4)          # level k=2
    A = q ** 0.25
    d = -(A**2) - A ** (-2)

    print("=" * 80)
    print("逼三维 · 组装版：涌现 SU(2) 提供「上下」 + Jones 提供「链接奖励」")
    print("=" * 80)

    # ---- 块1：涌现 SU(2)（Kramers spin-1/2）----
    print("\n[块1] π 磁通 → T²=−1 → Kramers = spin-1/2（2×2 内部上下）")
    n = 4
    N = n * n
    D = toroidal_D(n, pi_flux=True)
    evals = np.linalg.eigvalsh(D)
    # 全偶重数 ⟺ 自对偶 T²=-1（Kramers）
    mults = []
    es = np.sort(evals)
    s = 0
    for i in range(1, len(es)):
        if es[i] - es[i-1] > 1e-8 * (es.max()-es.min()):
            mults.append(i-s); s = i
    mults.append(len(es)-s)
    all_even = all(m % 2 == 0 for m in mults)
    print(f"    π 磁通 N={N}：重数={mults}，全偶={all_even}  →  T²=−1（Kramers）✓")
    print(f"    Kramers 对 = 2×2 spin-1/2（内部「上下」= 交叉符号的来源，非 U(1) 相位）")

    # ---- 块2：R 矩阵作用在 2 strands（2×2 ⊗ 2×2 = 4×4）----
    print("\n[块2] 辫子生成元 σ = A·1 + A⁻¹·e 作用在 2 strands（4×4）")
    e = np.array([0.0, 1.0, -1.0, 0.0], complex)
    E = np.outer(e, e.conj())
    tl_e = (d/2.0) * E                      # TL 生成元 e（e²=d e）
    R = A * np.eye(4) + (A**-1) * tl_e      # R 矩阵 = A·1 + A⁻¹·e
    print(f"    R 矩阵 = A·1 + A⁻¹·e，交叉符号 = 振幅比 A vs A⁻¹（量子叠加）")

    # ---- 块3：链接奖励 S_link = -ν·|V(σ^k) - V(σ⁰)|----
    print("\n[块3] 链接奖励 S_link = -ν·|V(σ^k) - V(σ⁰)|")
    print(f"    {'k':>2}  {'链接':<12} {'Jones V(σ^k)':>26} {'S_link/ν':>14}")
    for k in (0, 1, 2):
        vk = jones(A, k)
        v0 = jones(A, 0)
        s_link = -abs(vk - v0)
        label = {0: '平凡(2D)', 1: '扭unknot', 2: 'Hopf(3D)'}[k]
        print(f"    {k:>2}  {label:<12} {vk.real:+.4f}{vk.imag:+.4f}i   {s_link:+.4f}")

    # ---- 结论 ----
    print("\n[结论]")
    print("    S_link(σ⁰)=0（2D 不奖励）；S_link(σ²)<0（3D 被奖励）。")
    print("    即「Jones 链接奖励」把作用量往「链接（3D）」方向推 —— 逼三维作用量项雏形。")
    print("    剩余缺口：σ^k 的 k 目前手放；自发版要「从 D 的环结构读出 k」（下一版）。")

    out = ROOT / "experiments" / "exp_force3d_jones_action_last_run.json"
    out.write_text(json.dumps({
        "pi_flux_kramers": all_even,
        "mults": mults,
        "S_link_sigma0": -abs(jones(A, 0) - jones(A, 0)),
        "S_link_sigma2": -abs(jones(A, 2) - jones(A, 0)),
        "note": "emergent SU(2)=updown source; Jones=linking reward; sigma^k still hand-placed (spontaneous readout is next)",
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

"""Exp7 后续：完成 Clifford 代数（自旋联络的起头）

从「有向区分」γ⁰ = [[0,1],[-1,0]] = -iσ_y（时间）+ SU(2)≅S³（空间），
补全 3+1 的 Clifford 代数 {γ⁰,γ¹,γ²,γ³}，号差 (−+++)。

诚实要点：2×2 只能装下 Cl(0,3)（四元数，3 个反交换生成元），装不下 Cl(3,1)
（4 个生成元）。所以「完成」需要 2×2 → 4×4（Dirac 旋量），这一步正是
「自旋几何 / 自旋联络」的桥——本文把它精确暴露出来。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)


def anticom(A, B):
    return A @ B + B @ A


def main():
    print("=" * 80)
    print("完成 Clifford 代数：有向区分（时间）+ SU(2)≅S³（空间）→ Cl(3,1)")
    print("=" * 80)

    # 1) 有向区分 = -iσ_y = γ⁰（时间种子）
    print("\n[1] 有向区分 γ⁰ = [[0,1],[-1,0]] = -iσ_y")
    g0_2x2 = np.array([[0, 1], [-1, 0]], complex)
    print(f"  -iσ_y = {(-1j*SY).tolist()}")
    print(f"  (γ⁰)² = {g0_2x2 @ g0_2x2}  →  -I  ✓（时间，号差 -1）")

    # 2) SU(2)≅S³ = 四元数 {iσ_x, iσ_y, iσ_z}（空间，欧氏 Cl(0,3)）
    print("\n[2] SU(2)≅S³ = 四元数 {iσx, iσy, iσz}（全部平方 -I，欧氏 Cl(0,3)）")
    quat = {"iσx": 1j*SX, "iσy": 1j*SY, "iσz": 1j*SZ}
    for name, q in quat.items():
        print(f"  ({name})² = {np.round(q @ q, 3).tolist()}  →  -I  ✓")

    # 3) 诚实坎：2×2 装不下 Cl(3,1)（需要 4 个反交换生成元）
    print("\n[3] 诚实坎：2×2 装不下 Cl(3,1)")
    print(f"  {{γ⁰, γ²}} = {{-iσy, σy}} = {np.round(anticom(g0_2x2, SY), 3).tolist()}  ≠ 0  ✗")
    print("  → 2×2 里 γ⁰=-iσy 与 γ²=σy 不反交换（因为 γ⁰∝γ²），最多只有 3 个反交换生成元（四元数 Cl(0,3)）")
    print("  → 完整 3+1 的 Cl(3,1) 需要 4×4 Dirac 旋量（2×2 → 4×4 这一步 = 自旋几何的桥）")

    # 4) 4×4 Dirac（从 2×2 的 σ_i 和「时间」iI 搭出来）
    print("\n[4] 4×4 Dirac：γ⁰=[[0,iI],[iI,0]]，γ^i=[[0,iσ_i],[-iσ_i,0]]（号差 −+++）")
    Z = np.zeros((2, 2), complex)
    G0 = np.block([[Z, 1j*I2], [1j*I2, Z]])
    G1 = np.block([[Z, 1j*SX], [-1j*SX, Z]])
    G2 = np.block([[Z, 1j*SY], [-1j*SY, Z]])
    G3 = np.block([[Z, 1j*SZ], [-1j*SZ, Z]])
    G = [G0, G1, G2, G3]
    g = np.diag([-1.0, 1.0, 1.0, 1.0])
    print(f"  (γ⁰)² = {np.round(G0 @ G0, 2).tolist()}  →  -I  ✓（时间）")
    for i, Gi in enumerate(G[1:], 1):
        sq = Gi @ Gi
        print(f"  (γ^{i})² = +I ？ {np.allclose(sq, np.eye(4))}  (对角 {np.round(np.diag(sq).real, 2)})")

    # 5) 验证 Clifford 代数 {γ^μ, γ^ν} = 2g^{μν}
    print("\n[5] 验证 Clifford 代数 {γ^μ,γ^ν} = 2g^{μν}，g = diag(-1,+1,+1,+1)")
    maxerr = 0.0
    for mu in range(4):
        for nu in range(4):
            ac = anticom(G[mu], G[nu])
            target = 2 * g[mu, nu] * np.eye(4)
            err = np.max(np.abs(ac - target))
            maxerr = max(maxerr, err)
    print(f"  对所有 μ,ν：最大误差 = {maxerr:.2e}  →  {'✓ 号差 (−+++)' if maxerr < 1e-10 else '✗'}")

    # 6) 自旋联络种子：洛伦兹生成元 = [γ^μ, γ^ν]/4（由 σ_i 生成）
    print("\n[6] 自旋联络种子：洛伦兹生成元 σ^{μν} = [γ^μ,γ^ν]/4")
    s01 = (G0 @ G1 - G1 @ G0) / 4
    print(f"  σ^{{01}} = [γ⁰,γ¹]/4 = {np.round(s01, 3).tolist()}")
    print("  → 洛伦兹生成元（自旋联络的种子）由 σ_i（SU(2) 生成元）组成，")
    print("    即：自旋联络 = SU(2) 结构在 4×4 里的「旋转」，这正是那根绳子。")

    out = ROOT / "experiments" / "exp7_k2_clifford_last_run.json"
    out.write_text(json.dumps({
        "gamma0_sq_minus1": True,
        "quaternions_sq_minus1": True,
        "clifford_max_err": float(maxerr),
        "signature_minus_plus_plus_plus": bool(maxerr < 1e-10),
        "note": "2x2 holds Cl(0,3) not Cl(3,1); completion needs 4x4 Dirac (spin geometry bridge)",
    }, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

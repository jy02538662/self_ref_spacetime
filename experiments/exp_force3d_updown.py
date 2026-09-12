"""逼三维 · 上下 = π 磁通反对易（最终结论的数值验证）

核心结论（2026-09-09 收口）：「上下」（交叉符号）不在外部 3D 嵌入里，而在 D 自己的
π 磁通反对易里。π 磁通的磁平移 T_x（绕 x）、T_y（绕 y）反对易：

    T_x T_y = -T_y T_x     （非阿贝尔 / SU(2) / 内部三维 SU(2)≅S³）

每个 plaquette 的 holonomy = -1（π 磁通）。

这个「-1 符号 / 反对易」就是「上下」——有限性（π 磁通的 Z₂）产生反对易（上下），
正是「自指的有限性认识产生上下」的精确落点。它是非阿贝尔（SU(2）结构 = 内部三维，
而 SU(2) 已自发涌现（π 磁通 → T²=−1 → Kramers，见 D谱k2 记录）。

本脚本只验证这一个事实；完整推导见 vault《逼三维：涌现SU(2)的上下与Jones链接奖励》。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

from experiments.exp_pure_spectral_anneal import toroidal_D


def magnetic_translations(n):
    """π 磁通 toroidal 的磁平移 T_x（绕 x，相位 0）、T_y（绕 y，相位 πj）。"""
    N = n * n
    Tx = np.zeros((N, N), complex)
    Ty = np.zeros((N, N), complex)
    for i in range(n):
        for j in range(n):
            idx = n * i + j
            Tx[idx, n * i + ((j + 1) % n)] = 1.0
            Ty[idx, n * ((i + 1) % n) + j] = np.exp(1j * np.pi * j)
    return Tx, Ty


def main():
    print("=" * 72)
    print("上下 = π 磁通反对易：T_x T_y = -T_y T_x")
    print("=" * 72)
    results = {}
    for n in (4, 6, 8):
        Tx, Ty = magnetic_translations(n)
        comm = float(np.max(np.abs(Tx @ Ty - Ty @ Tx)))
        anti = float(np.max(np.abs(Tx @ Ty + Ty @ Tx)))
        hol = Tx @ Ty @ Tx.conj().T @ Ty.conj().T
        hol_diag = float(np.real(np.trace(hol))) / (n * n)  # 每个 plaquette 的 holonomy（平均）
        anticommute = anti < 1e-10 and comm > 1.0
        results[f"n{n}"] = {"comm_norm": comm, "anti_norm": anti, "holonomy_mean": hol_diag, "anticommute": anticommute}
        print(f"n={n}: |T_xT_y - T_yT_x|={comm:.2e}  |T_xT_y + T_yT_x|={anti:.2e}  "
              f"holonomy均值={hol_diag:+.2f}  →  {'反对易（非阿贝尔=SU(2)=内部三维）' if anticommute else '对易（阿贝尔）'}")

    print("\n结论：π 磁通给磁平移反对易 T_xT_y = -T_yT_x（holonomy=-1），")
    print("      这个「-1 符号」就是「上下」，在 D 自己结构里，不需要外部 3D 嵌入。")
    print("      反对易 = 非阿贝尔 = SU(2) = 内部三维（SU(2)≅S³）——逼三维（内部三维）闭环。")

    out = ROOT / "experiments" / "exp_force3d_updown_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

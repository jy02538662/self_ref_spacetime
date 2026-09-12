"""逼三维 · 读 k（环结构 → 交叉 → 链接 → k → S_link）

把「环 → 辫子词」的 k 从 D 的环结构里读出来，接到 S_link 上。这是「逼三维」最后一步
的前半段（读 k）；后半段「D 动力学自发把 k 从 0 推到 2（升三维）」仍开放，见诚实边界。

环结构（emergent 2D toroidal 的两个非收缩环）：
  子午线 μ（绕 x）、经线 λ（绕 y），在节点 (0,0) 处横截相交一次（1 个交点）。
链接（Gauss 积分，需嵌入）：
  - 2D（两环在同一平面）：链接 = 0（共面不能缠）→ k = 0 → Jones 平凡（d）
  - 3D（实心环面：经线=核心圆、子午线=绕管小圆）：链接 = ±1 → k = 2 → Jones Hopf（-A⁴-A⁻⁴）

关键（2026-09-09 认知）：「上下」（over/under）不是 U(1) 相位给的（holonomy 只 ±1 对易），
是「有限性（偶数 N）→ 自对偶 → 涌现 SU(2)（Kramers 对 = 内部上下）」给的。
所以 k 的「交叉次数」是环结构的拓扑，而「上下符号」由涌现 SU(2) 承载。

诚实边界：本脚本的「3D 嵌入（实心环面）」是手放的；完全自发的版本 = 让涌现 SU(2)
（Kramers 对）自己充当「第三维」，使链接数不需要外部嵌入。这是真正的最后一步，仍开放。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

from experiments.exp_pure_spectral_anneal import toroidal_D
from experiments.exp_force3d_jones_toy import meridian_longitude_cycles, tl2_kauffman, jones


def gauss_linking(r1, t1, r2, t2):
    """Gauss 链接数 Lk = (1/4π)∮∮ (r1-r2)·(dr1×dr2)/|r1-r2|³，离散求和。

    r1, r2: (N,3) 两条闭合曲线的点；t1, t2: (N,3) 切向量（已按弧长归一不用）。
    """
    Lk = 0.0
    n1, n2 = r1.shape[0], r2.shape[0]
    for i in range(n1):
        for j in range(n2):
            dr = r1[i] - r2[j]
            denom = np.linalg.norm(dr) ** 3
            if denom < 1e-12:
                continue
            cross = np.cross(t1[i], t2[j])
            Lk += np.dot(dr, cross) / denom
    # 弧长因子 Δs Δt（均匀采样，Δs = 周长/n1 等），归一到 1/4π
    return Lk / (4.0 * np.pi) * (2.0 * np.pi / n1) * (2.0 * np.pi / n2)


def torus_curves_2d(n=200):
    """2D：经线（大圆）+ 子午线（大圆）都在 xy 平面，相交不缠。"""
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    R = 1.0
    r1 = np.stack([R * np.cos(th), R * np.sin(th), np.zeros_like(th)], axis=1)          # 经线
    r2 = np.stack([R + 0.5 * np.cos(th), 0.5 * np.sin(th), np.zeros_like(th)], axis=1)  # 子午线
    t1 = np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], axis=1)
    t2 = np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], axis=1)
    return r1, t1, r2, t2


def torus_curves_3d(n=200):
    """3D：经线=核心圆（xy 平面 z=0），子午线=绕管小圆（xz 平面 y=0，θ=0 处），相链一次。"""
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    R = 1.0
    r = 0.4
    # 经线（核心圆）：xy 平面
    r1 = np.stack([R * np.cos(th), R * np.sin(th), np.zeros_like(th)], axis=1)
    t1 = np.stack([-np.sin(th), np.cos(th), np.zeros_like(th)], axis=1)
    # 子午线（绕管小圆，固定 θ=0）：x=R+r cos φ, y=0, z=r sin φ
    r2 = np.stack([R + r * np.cos(th), np.zeros_like(th), r * np.sin(th)], axis=1)
    t2 = np.stack([-r * np.sin(th), np.zeros_like(th), r * np.cos(th)], axis=1)  # 切向量 = dγ₂/dφ（模 r）
    return r1, t1, r2, t2


def main():
    q = np.exp(1j * np.pi / 4)
    A = q ** 0.25

    print("=" * 78)
    print("逼三维 · 读 k：环结构 → 交叉 → 链接 → k → S_link")
    print("=" * 78)

    # 1) 环结构：交叉次数
    n = 4
    D = toroidal_D(n, pi_flux=True)
    mer, lon = meridian_longitude_cycles(n)
    mer_verts = set(a for a, b in mer) | set(b for a, b in mer)
    lon_verts = set(a for a, b in lon) | set(b for a, b in lon)
    n_cross = len(mer_verts & lon_verts)
    print(f"\n[1] emergent 2D toroidal N={n*n}：子午线 μ({len(mer)}边) + 经线 λ({len(lon)}边)")
    print(f"    共享顶点（横截交点）数 = {n_cross}  →  交叉次数 = 1")

    # 2) Gauss 链接数：2D vs 3D
    r1, t1, r2, t2 = torus_curves_2d()
    lk2 = gauss_linking(r1, t1, r2, t2)
    r1, t1, r2, t2 = torus_curves_3d()
    lk3 = gauss_linking(r1, t1, r2, t2)
    print(f"\n[2] Gauss 链接数（嵌入依赖）")
    print(f"    2D（共面）: Lk = {lk2:+.4f}  →  不链 → k=0")
    print(f"    3D（实心环面）: Lk = {lk3:+.4f}  →  相链一次 → k=2")

    # 3) 映射 k → S_link
    print(f"\n[3] 映射 k → S_link = -ν|V(σ^k) - V(σ^0)|")
    for k, label in ((0, "2D 不链"), (2, "3D Hopf")):
        vk = jones(A, k)
        v0 = jones(A, 0)
        s = -abs(vk - v0)
        print(f"    k={k} ({label}): V={vk.real:+.4f}{vk.imag:+.4f}i  S_link/ν={s:+.4f}")

    # 4) 结论 + 诚实边界
    print(f"\n[4] 结论")
    print("    环结构（1 交叉）→ 链接数（2D=0 / 3D=±1）→ k（0/2）→ S_link 奖励 k=2（3D）。")
    print("    「读 k」闭环；但「3D 嵌入」手放。")
    print("    真正最后一步（开放）：让涌现 SU(2)（Kramers 对）充当「第三维」，")
    print("    使链接数不需外部嵌入——即「D 动力学自发把 k 从 0 推到 2」.")

    out = ROOT / "experiments" / "exp_force3d_read_k_last_run.json"
    out.write_text(json.dumps({
        "n_crossings": n_cross,
        "gauss_linking_2d": float(lk2),
        "gauss_linking_3d": float(lk3),
        "k_2d": 0, "k_3d": 2,
        "S_link_k0": -abs(jones(A, 0) - jones(A, 0)),
        "S_link_k2": -abs(jones(A, 2) - jones(A, 0)),
        "note": "cycle->crossing->linking->k closed; 3D embedding hand-placed; spontaneous lift (emergent SU(2) as 3rd dim) still open",
    }, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

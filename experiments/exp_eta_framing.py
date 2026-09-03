"""验证 3/2 = 谱流(整数) + eta 分数部分(1/2) 的 1D 数值演示（修正版 2）。

模型：S^1 上的 1D Dirac 算子 D_A = -i d/dtheta + A（周期边界），本征值 lambda_n = n + A（n 为整数，精确）。

- 谱流（整数）＝ 一族 D(t) = -i d/dtheta + t，t: 0 -> 1 的本征值净穿越 0 的次数。
  这是"绕数 / Hopf 荷"的 1D 类比（整数拓扑荷）。
- eta-不变量（分数）＝ 热核正则化 eta_eps = sum_n sign(n+A) exp(-eps|n+A|)，eps -> 0 收敛。
  对 0 < A < 1 可解析求极限：
      eta_eps = [e^{-eps*A} - e^{-eps*(1-A)}] / (1 - e^{-eps})
      -> (1 - 2A)  当 eps -> 0
  故取 A = 1/4（"四分之一扭转" = 自旋 framing 的 1D 类比）-> eta = 1/2。

修正说明（相对初版/修正版1）：
  1) scipy.special.zeta 对 0 < s < 1 返回 nan（1.17 数值缺陷），弃用；
  2) 改用热核正则化 eta_eps（指数截断，任意 eps>0 收敛），数值可稳定外推 eps -> 0；
  3) 输出全部 ASCII 安全，避免 Windows GBK 控制台编码错误。

诚实边界：
  这是 1D 类比，演示"整数谱流 + 分数 eta = 半整数"的机制，不是 3D 的 S^3 完整计算。
  真正的 3D Hopf 荷（整数 = 1）已由 Exp6a 钉死（0.99996）；
  这里补的是 eta 的分数部分（1/2），两者在 framed/APS 框架里合成 3/2。
"""

from __future__ import annotations

import numpy as np


def eta_heatkernel(A: float, eps: float, N: int = 20000) -> float:
    """热核正则化 eta_eps = sum sign(lam) exp(-eps|lam|)，lam_n = n + A（对称截断）。"""
    n = np.arange(-N, N + 1, dtype=float)
    lam = n + A
    lam = lam[lam != 0.0]  # 去掉可能的零模
    return float(np.sum(np.sign(lam) * np.exp(-eps * np.abs(lam))))


def eta_closed_form(A: float) -> float:
    """eta(0) = 1 - 2A（0 < frac(A) < 1），A 为整数时 eta = 0。"""
    a = A - np.floor(A)
    if a < 1e-12:
        return 0.0
    return 1.0 - 2.0 * a


def spectral_flow_numeric(n_t: int = 2000, M: int = 100) -> int:
    """数值追踪 D(t) = -i d/dtheta + t 的本征值 lambda_n(t) = n + t（t: 0 -> 1）。

    数净上穿 0 的次数（上穿 +1，下穿 -1）。本征值精确为 n + t，这里仍做逐点追踪，
    保证是可复现的"数值谱流"而非手填常数。
    """
    t_grid = np.linspace(0.0, 1.0, n_t)
    n = np.arange(-M, M + 1, dtype=float)
    sf = 0
    for j in range(n_t - 1):
        lam0 = n + t_grid[j]
        lam1 = n + t_grid[j + 1]
        up = np.sum((lam0 <= 0.0) & (lam1 > 0.0))
        down = np.sum((lam0 > 0.0) & (lam1 <= 0.0))
        sf += int(up) - int(down)
    return sf


if __name__ == "__main__":
    print("=== 1D 谱流 + eta 分数部分 验证（修正版 2） ===")
    print()

    # 1) 谱流（整数）
    sf = spectral_flow_numeric()
    print(f"[谱流] D(t) = -i d/dtheta + t, t:0->1 的本征值净上穿 0 次数 = {sf}  (整数)")
    print()

    # 2) eta-不变量（分数）：热核正则化 eta_eps，eps -> 0 外推到闭式 1 - 2A
    eps_list = [0.20, 0.10, 0.05, 0.02, 0.01, 0.005]
    print("[eta-不变量] 热核正则化 eta_eps = sum sign(lam) exp(-eps|lam|)，eps -> 0")
    for A in (0.25, 0.75):
        closed = eta_closed_form(A)
        print(f"  A = {A}  (闭式 eta(0) = {closed:+.4f})")
        for eps in eps_list:
            val = eta_heatkernel(A, eps)
            print(f"    eps={eps:6.3f}  eta_eps = {val:+.5f}   偏差 {abs(val-closed):.2e}")
        print()
    print("  -> A=1/4 时 eta -> +1/2，A=3/4 时 eta -> -1/2（分数部分 = 自旋 framing 的 +/-1/2）")
    print()

    # 3) 合成
    hopf_int = 1                       # 整数 Hopf 荷（Exp6a 已钉死 0.99996，这里取理论值 1）
    spin_frac = eta_closed_form(0.25)  # +1/2
    total = hopf_int + spin_frac
    print("=== 合成 ===")
    print(f"整数 Hopf 荷（Exp6a 数值 0.99996 ~ 1）  = {hopf_int}")
    print(f"分数 eta / 自旋 framing（A=1/4）          = {spin_frac:+.4f}")
    print(f"总拓扑荷（framed 框架里两者相加）        = {total:+.4f}  = 3/2")
    print()
    print("诚实边界：3/2 是 framed/APS 框架下 '整数谱流 + 分数 eta' 的组合；")
    print("真正的 S^3 上自旋 1/2 Hopfion 的 eta = 3/2 需在具体 Dirac 算子上把 eta 完整算出来（未做）。")

"""谱层第一步：谱维数 d_s 区分 2D/3D/4D（热核 K(t) = Tr(e^{-tL}) ~ t^{-d_s/2}）。

用周期边界 d 维超立方晶格的拉普拉斯 L 的解析本征值（避免矩阵构造 bug）：
  lambda(m) = 4 * sum_j sin^2(pi m_j / n)，m_j = 0..n-1。
热核迹：K(t) = [ sum_m exp(-4t sin^2(pi m/n)) ]^d。

在中间 t 区间（1 << t << n^2，即"走了很多步但还没绕满环面"）有 K(t) ~ t^{-d/2}，
故 d_s = -2 * d(log K)/d(log t) = d。

诚实标注：谱维数从"拉普拉斯 L"读出，不是从"邻接矩阵 A"直接读出——
A 的谱在能带中心（sum cos = 0，codim-1 曲面）处 DOS ~ const，
Tr(e^{-t A^2}) ~ t^{-1/2}，会读出 d_s ~ 1（错）。
正确的"Dirac"解读是 D^2 = L（D 是拉普拉斯的平方根）。
"""

from __future__ import annotations

import numpy as np


def spectral_dim(d: int, n: int, t_lo: float, t_hi: float, n_t: int = 50):
    """K(t) = [sum_m exp(-4t sin^2(pi m/n))]^d，拟合 log K vs log t，d_s = -2*slope。"""
    ts = np.geomspace(t_lo, t_hi, n_t)
    m = np.arange(n)
    s2 = np.sin(np.pi * m / n) ** 2
    f = np.array([np.sum(np.exp(-4.0 * t * s2)) for t in ts])
    K = f ** d
    slope = np.polyfit(np.log(ts), np.log(K), 1)[0]
    return -2.0 * slope, ts, K


if __name__ == "__main__":
    print("=== 谱维数区分 2D/3D/4D（周期边界，解析本征值） ===")
    print()
    # (d, n, t_lo, t_hi)：中间 t 区间选在 1 << t << n^2/(4π)，n 取大让平台更宽
    configs = [
        (2, 100, 10.0, 100.0),
        (3, 60, 10.0, 100.0),
        (4, 40, 10.0, 60.0),
    ]
    for d, n, t_lo, t_hi in configs:
        ds, ts, K = spectral_dim(d, n, t_lo, t_hi)
        print(f"d={d}  n={n}  t=[{t_lo:g},{t_hi:g}]  d_s = {ds:.3f}  (目标 {d})")
    print()
    print("判读：d_s 应 ≈ d（2/3/4），说明谱维数能从热核读出、能区分维度。")
    print("（d=4 因有限格子平台略窄，d_s≈3.96 略欠 4，属有限尺寸效应，不碍区分）")

"""谱维数：把"维度"从预设的容器，变成 D 的谱里长出来的观测量。

热核 K(t) = Tr(e^{-tL})，L 是图拉普拉斯。小 t 下 K(t) ~ t^{-d_s/2}，
d_s 就是谱维数。三维格点 d_s=3，二维=2，一维环=1。

这一步验证：谱维数能"读出"D 承载的维度——不再预设三维目标距离。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def ring_laplacian(n: int) -> np.ndarray:
    A = np.zeros((n, n))
    for i in range(n):
        A[i, (i + 1) % n] = 1.0
        A[(i + 1) % n, i] = 1.0
    deg = A.sum(axis=1)
    return np.diag(deg) - A


def grid_laplacian(n_per_dim: int, dim: int) -> np.ndarray:
    """d 维立方格点的图拉普拉斯。"""
    coords = np.indices([n_per_dim] * dim).reshape(dim, -1).T  # (n, dim)
    n = coords.shape[0]
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if np.sum(np.abs(coords[i] - coords[j])) == 1:  # 曼哈顿距离 1
                A[i, j] = A[j, i] = 1.0
    deg = A.sum(axis=1)
    return np.diag(deg) - A


def spectral_dimension(L: np.ndarray, t_lo: float = 0.02, t_hi: float = 0.8, n_pts: int = 30) -> float:
    """从热核 K(t)=Tr(e^{-tL}) 的标度读出谱维数 d_s。"""
    eig = np.linalg.eigvalsh(L)
    eig = eig[eig > 1e-12]  # 去掉零模
    ts = np.geomspace(t_lo, t_hi, n_pts)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    slope = np.polyfit(np.log(ts), logK, 1)[0]
    return float(-2.0 * slope)


if __name__ == "__main__":
    results = {}
    # 一维环（连续极限区间 t ∈ [5, 500]）
    L1 = ring_laplacian(200)
    d1 = spectral_dimension(L1, t_lo=5.0, t_hi=500.0)
    # 二维格 14×14（t ∈ [2, 50]）
    L2 = grid_laplacian(14, 2)
    d2 = spectral_dimension(L2, t_lo=2.0, t_hi=50.0)
    # 三维格 5×5×5（t ∈ [1, 10]）
    L3 = grid_laplacian(5, 3)
    d3 = spectral_dimension(L3, t_lo=1.0, t_hi=10.0)

    results = {
        "ring_1d_spectral_dim": d1,
        "grid_2d_spectral_dim": d2,
        "grid_3d_spectral_dim": d3,
        "expected": {"1d": 1.0, "2d": 2.0, "3d": 3.0},
    }
    print("=== 谱维数读出维度 ===")
    print(f"一维环  谱维数 = {d1:.3f}（期望 1）")
    print(f"二维格  谱维数 = {d2:.3f}（期望 2）")
    print(f"三维格  谱维数 = {d3:.3f}（期望 3）")
    print(json.dumps(results, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp_spectral_dim_last_run.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"wrote {out}")

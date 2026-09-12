"""Exp7 后续：自发 vs 手放 —— 纯迹作用量是否「自发」破 C₄ / 自发加质量？

问题（v4 逼三维）：纯 k=2 需要 破 C₄（各向异性）+ 交错质量（gap Dirac）。这两个是
自发的（作用量自己想要的），还是手放的（必须外力）？

判据 = 作用量在对称点（ε=0 各向同性、m=0 无质量）的**曲率（二阶导）**：
  - 曲率 < 0（对称点是不稳定极大）→ 自发破缺（墨西哥帽式自发倾斜）；
  - 曲率 > 0（对称点稳定极小）→ 手放（需外力倾斜）。

解析预期（墨西哥帽 S = −α Tr(D²) + β Tr(D⁴)）：
  - Tr(D²) = N m² + 2N(t_x²+t_y²)：对 ε 与 m 都是凸的 → −α Tr(D²) 罚「各向异性 + 质量」；
  - π 磁通奖励（plaquette 项 ∝ t_x² t_y²，AM-GM 在 t_x=t_y 最大）→ 各向同性最被奖励；
  → 预期：纯迹作用量偏好 各向同性 + 无质量 → 破 C₄ 与加质量都是「手放」（对称点是稳定极小）。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def anisotropic_piflux_D(n_per_dim: int, tx: float, ty: float, m: float = 0.0) -> np.ndarray:
    """π 磁通 toroidal，各向异性跳变 + 交错质量。"""
    n = n_per_dim
    D = np.zeros((n * n, n * n), complex)
    for i in range(n):
        for j in range(n):
            idx = n * i + j
            jr = (j + 1) % n
            D[idx, n * i + jr] += tx
            D[n * i + jr, idx] += tx
            id_ = (i + 1) % n
            ph = np.exp(1j * np.pi * j)
            D[idx, n * id_ + j] += ty * ph
            D[n * id_ + j, idx] += ty * np.conj(ph)
    if m != 0.0:
        signs = np.array([(-1) ** (i + j) for i in range(n) for j in range(n)], dtype=float)
        D += m * np.diag(signs)
    return D


def action(D: np.ndarray, alpha: float, beta: float) -> dict:
    D2 = D @ D
    tr2 = float(np.real(np.trace(D2)))
    tr4 = float(np.real(np.trace(D2 @ D2)))
    return {"S": -alpha * tr2 + beta * tr4, "TrD2": tr2, "TrD4": tr4}


def curvature(vals: np.ndarray, xs: np.ndarray) -> float:
    """二阶导（中心差分）在对称点附近。"""
    # 用最小二乘二次拟合 a x² + b x + c，返回 2a（曲率）
    coeff = np.polyfit(xs, vals, 2)
    return float(2.0 * coeff[0])


def main():
    n_per_dim = 4
    alpha, beta = 2.0, 0.5  # 墨西哥帽（exp_pure_spectral_selfref 同款）

    print("=" * 80)
    print(f"自发 vs 手放：纯迹作用量 S = -{alpha} Tr(D²) + {beta} Tr(D⁴)（N=16 π-flux）")
    print("判据：对称点(ε=0, m=0)曲率 >0 = 稳定极小 = 手放；<0 = 不稳定 = 自发")
    print("=" * 80)

    results = {}

    # 1) 各向异性 ε（t_y = 1+ε），**固定标度 Tr(D²) 不变**，只看方向偏好
    #    修正：S = -α Tr(D²) + β Tr(D⁴) 在固定 Tr(D²) 下 ⟺ 最小化 Tr(D⁴)。
    #    所以归一化 D 使 Tr(D²)=Tr(D²_iso)，再看 Tr(D⁴) 在 ε=0（各向同性）还是 ε≠0 最小。
    print("\n[各向异性 ε：t_x=1, t_y=1+ε，固定标度 Tr(D²)]")
    D_iso = anisotropic_piflux_D(n_per_dim, 1.0, 1.0, 0.0)
    tr2_iso = float(np.real(np.trace(D_iso @ D_iso)))
    eps = np.linspace(-0.5, 0.5, 11)
    TrD4s = []
    for e in eps:
        D = anisotropic_piflux_D(n_per_dim, 1.0, 1.0 + e, 0.0)
        D = D * np.sqrt(tr2_iso / float(np.real(np.trace(D @ D))))  # 固定 Tr(D²)
        TrD4s.append(float(np.real(np.trace((D @ D) @ (D @ D)))))
    TrD4s = np.array(TrD4s)
    curv_e = curvature(TrD4s, eps)
    print(f"  ε: {np.round(eps, 2).tolist()}")
    print(f"  Tr(D⁴)@固定标度: {np.round(TrD4s, 2).tolist()}")
    print(f"  Tr(D⁴)(ε=0)={TrD4s[len(eps)//2]:.2f}  最小={TrD4s.min():.2f}  在 ε={eps[TrD4s.argmin()]:.2f}")
    print(f"  曲率 = {curv_e:.2f}  → {'各向同性是稳定极小（手放破 C₄）' if curv_e > 0 else '各向同性不稳定（自发破 C₄）'}")
    results["anisotropy"] = {"eps": eps.tolist(), "TrD4_fixed_scale": TrD4s.tolist(), "curvature": curv_e,
                             "verdict": "手放" if curv_e > 0 else "自发"}

    # 2) 交错质量 m，各向同性
    print("\n[交错质量 m：各向同性 t_x=t_y=1]")
    ms = np.linspace(0.0, 1.0, 11)
    Ss = []
    for m in ms:
        D = anisotropic_piflux_D(n_per_dim, 1.0, 1.0, m)
        Ss.append(action(D, alpha, beta)["S"])
    Ss = np.array(Ss)
    curv_m = curvature(Ss, ms)
    print(f"  m: {np.round(ms, 2).tolist()}")
    print(f"  S: {np.round(Ss, 3).tolist()}")
    print(f"  S(m=0)={Ss[0]:.3f}  最小 S={Ss.min():.3f}  在 m={ms[Ss.argmin()]:.2f}")
    print(f"  曲率(二阶导) = {curv_m:.3f}  → {'稳定极小（手放加质量）' if curv_m > 0 else '不稳定（自发加质量）'}")
    results["mass"] = {"m": ms.tolist(), "S": Ss.tolist(), "curvature": curv_m,
                       "verdict": "手放" if curv_m > 0 else "自发"}

    # 3) 各向异性 + 质量联合（m=0.5 固定，扫 ε）
    print("\n[联合：m=0.5 固定，扫 ε]")
    eps = np.linspace(-0.4, 0.4, 9)
    Ss = []
    for e in eps:
        D = anisotropic_piflux_D(n_per_dim, 1.0, 1.0 + e, 0.5)
        Ss.append(action(D, alpha, beta)["S"])
    Ss = np.array(Ss)
    curv_e_m = curvature(Ss, eps)
    print(f"  S: {np.round(Ss, 3).tolist()}")
    print(f"  曲率 = {curv_e_m:.3f}  → {'稳定极小（手放）' if curv_e_m > 0 else '不稳定（自发）'}")

    out = ROOT / "experiments" / "exp7_k2_spontaneous_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

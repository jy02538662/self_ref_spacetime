"""Exp7 后续：破 C₄（各向异性跳变）—— 看 ±2 的 4 重能否劈成 2+2，逼近「纯 k=2」。

背景：π 磁通 N=16 谱 = [±2√2(2)、±2(4)、0(4)]。其中 ±2(4) 是 C₄ 旋转（k_x↔k_y）
关联的偶然简并，0(4) 是 Dirac 点。用各向异性跳变（t_x ≠ t_y）破 C₄。

能带（各向异性，t_x=水平、t_y=垂直）：E(k) = ±√((2t_y cos k_x)² + (2t_x cos k_y)²)。
  - Dirac 点（E=0，cos kx=cos ky=0）在 (±π/2,±π/2)：两项都 0 → **不随 t 变，0 模不劈**；
  - ±2(4) 来自 (π/2,0)/(0,π/2) 等：各向异性下 → ±2t_x(2) + ±2t_y(2)，**4 重劈成 2+2**；
  - 带边 ±2√2 → ±2√(t_x²+t_y²)(2)，仍 2 重。

判据：
  - 只破 C₄：±2(4) → 2+2，但 0(4) 不劈（Dirac 点保护）→ 不是「纯 k=2」；
  - 破 C₄ + 交错质量 m：0(4) → ±m(2)，于是全谱 → 8 个二重 Kramers 对 = **纯 k=2**。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def anisotropic_piflux_D(n_per_dim: int, tx: float, ty: float) -> np.ndarray:
    """π 磁通 toroidal，各向异性跳变：水平(j)=tx，垂直(i)=ty·e^{iπj}。"""
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
    return D


def staggered_mass(n_per_dim: int) -> np.ndarray:
    n = n_per_dim
    signs = np.array([(-1) ** (i + j) for i in range(n) for j in range(n)], dtype=float)
    return np.diag(signs)


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
    n_per_dim = 4
    n = n_per_dim ** 2
    print("=" * 80)
    print(f"破 C₄：各向异性跳变 t_x=1, t_y=1+ε  （N={n} π-flux）")
    print("=" * 80)

    # 1) 只破 C₄（无质量）
    print("\n[只破 C₄，无交错质量]")
    for eps in [0.0, 0.1, 0.3, 0.5, 1.0]:
        D = anisotropic_piflux_D(n_per_dim, 1.0, 1.0 + eps)
        ev = np.sort(np.linalg.eigvalsh(D))
        mults = multiplicities(ev, 1e-8 * (ev.max() - ev.min()))
        # 正半轴能级（用于核对公式）
        pos = ev[ev > 1e-9]
        print(f"  eps={eps}: mults={mults.tolist()}")
        print(f"         evals={np.round(ev, 3).tolist()}")

    # 2) 破 C₄ + 交错质量 → 逼近纯 k=2
    print("\n[破 C₄ + 交错质量 m → 逼近纯 k=2]")
    M = staggered_mass(n_per_dim)
    for eps, m in [(0.3, 0.5), (0.3, 1.0), (1.0, 1.0)]:
        D = anisotropic_piflux_D(n_per_dim, 1.0, 1.0 + eps) + m * M
        ev = np.sort(np.linalg.eigvalsh(D))
        mults = multiplicities(ev, 1e-8 * (ev.max() - ev.min()))
        all_even = bool(np.all(mults % 2 == 0))
        pure_k2 = bool(np.all(mults == 2))
        print(f"  eps={eps}, m={m}: mults={mults.tolist()}  全偶={all_even}  纯k2={pure_k2}")
        print(f"         evals={np.round(ev, 3).tolist()}")

    # 3) 核对公式：±2t_x / ±2t_y / ±2√(t_x²+t_y²)（无质量、各向异性）
    print("\n[公式核对，无质量各向异性 eps=0.3]")
    tx, ty = 1.0, 1.3
    D = anisotropic_piflux_D(n_per_dim, tx, ty)
    ev = np.sort(np.linalg.eigvalsh(D))
    pos = ev[ev > 1e-9]
    print(f"  预测 2t_x={2*tx:.3f}, 2t_y={2*ty:.3f}, 2√(t_x²+t_y²)={2*np.sqrt(tx**2+ty**2):.3f}")
    print(f"  实测正半轴能级: {np.round(pos, 3).tolist()}")

    out = ROOT / "experiments" / "exp7_k2_break_c4_last_run.json"
    out.write_text(json.dumps({"note": "break C4 anisotropic hopping"}, indent=2), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

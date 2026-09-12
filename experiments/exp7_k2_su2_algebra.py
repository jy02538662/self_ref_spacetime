"""Exp7 后续：完整 SU(2) 代数（三算符存在性）

问题：π 磁通「纯 k=2」系统（8 个 Kramers 对）除了 T²=−1 反幺对称，是否还有完整的
SU(2) 代数（三个厄米 J_x, J_y, J_z，[J_i,J_j]=iε_ijk J_k，[J_i,H]=0）？

定理：Kramers（全偶重数）⇒ 完整 SU(2) 代数（自动）。构造 = 在每个 2 重简并子空间放
σ_i/2，J_i = V (⨁σ_i/2) V†。

验证四条：厄米、[J_i,H]=0、[J_x,J_y]=iJ_z（及循环）、Casimir J²=(3/4)I（自旋 1/2）。

诚实边界：这个 SU(2) 是「形式/涌现」的——对任何 Kramers 系统都成立（是对简并对的重新标号、
非局域的），不是「实质/物理」的（π 磁通是 spinless 费米子，无物理自旋；这里的「自旋 1/2」
是 Kramers 结构给的 pseudospin = 子格/谷自由度）。「实质」在 T²=−1（涌现时间反演）这一步。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

SX = np.array([[0, 1], [1, 0]], complex) / 2.0
SY = np.array([[0, -1j], [1j, 0]], complex) / 2.0
SZ = np.array([[1, 0], [0, -1]], complex) / 2.0


def anisotropic_piflux_D(n_per_dim: int, tx: float, ty: float, m: float = 0.0) -> np.ndarray:
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


def construct_su2(H: np.ndarray) -> dict:
    """构造 J_i = V (⨁σ_i/2) V† 并验证 SU(2) 代数。"""
    evals, V = np.linalg.eigh(H)
    n = H.shape[0]
    mults = multiplicities(evals, 1e-8 * (evals.max() - evals.min()))
    Jx = np.zeros((n, n), complex)
    Jy = np.zeros((n, n), complex)
    Jz = np.zeros((n, n), complex)
    pos = 0
    for m in mults:
        if m % 2 == 0:
            for _ in range(m // 2):
                Jx[pos:pos + 2, pos:pos + 2] = SX
                Jy[pos:pos + 2, pos:pos + 2] = SY
                Jz[pos:pos + 2, pos:pos + 2] = SZ
                pos += 2
        else:
            pos += m
    Jx = V @ Jx @ V.conj().T
    Jy = V @ Jy @ V.conj().T
    Jz = V @ Jz @ V.conj().T

    I = np.eye(n)
    hermiticity = max(np.max(np.abs(Jx - Jx.conj().T)),
                      np.max(np.abs(Jy - Jy.conj().T)),
                      np.max(np.abs(Jz - Jz.conj().T)))
    comm_H = max(np.max(np.abs(Jx @ H - H @ Jx)),
                 np.max(np.abs(Jy @ H - H @ Jy)),
                 np.max(np.abs(Jz @ H - H @ Jz)))
    alg_xy = np.max(np.abs((Jx @ Jy - Jy @ Jx) - 1j * Jz))
    alg_yz = np.max(np.abs((Jy @ Jz - Jz @ Jy) - 1j * Jx))
    alg_zx = np.max(np.abs((Jz @ Jx - Jx @ Jz) - 1j * Jy))
    casimir = np.max(np.abs(Jx @ Jx + Jy @ Jy + Jz @ Jz - 0.75 * I))
    # 非局域性：J 在「本征基」外是否稀疏（形式 SU(2) 应是稠密的、混所有 site）
    frac_nonsparse = float(np.mean(np.abs(Jx) > 1e-8))
    return {
        "mults": mults.tolist(),
        "hermiticity_err": float(hermiticity),
        "comm_H_err": float(comm_H),
        "su2_alg_err": float(max(alg_xy, alg_yz, alg_zx)),
        "casimir_err": float(casimir),
        "spin12": bool(casimir < 1e-8),
        "frac_nonsparse": frac_nonsparse,
        "Jz_evals": np.round(np.linalg.eigvalsh(Jz), 4).tolist(),
    }


def main():
    n_per_dim = 4
    print("=" * 80)
    print("完整 SU(2) 代数：构造 J_i=V(⨁σ_i/2)V†，验证 [J_x,J_y]=iJ_z、[J_i,H]=0、J²=3/4")
    print("=" * 80)

    results = {}
    for label, (tx, ty, m) in [("纯k2（破C4+质量）", (1.0, 1.3, 0.5)),
                                ("只交错质量", (1.0, 1.0, 0.5)),
                                ("只破C4", (1.0, 1.3, 0.0)),
                                ("纯π磁通（各向同性无质量）", (1.0, 1.0, 0.0))]:
        H = anisotropic_piflux_D(n_per_dim, tx, ty, m)
        r = construct_su2(H)
        results[label] = r
        print(f"\n[{label}]  重数={r['mults']}")
        print(f"  厄米误差={r['hermiticity_err']:.1e}  [J_i,H]误差={r['comm_H_err']:.1e}")
        print(f"  SU(2)代数误差={r['su2_alg_err']:.1e}  Casimir误差={r['casimir_err']:.1e}")
        print(f"  自旋1/2(J²=3/4)={r['spin12']}  J_z本征值={r['Jz_evals']}  非稀疏占比={r['frac_nonsparse']:.2f}")

    out = ROOT / "experiments" / "exp7_k2_su2_algebra_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

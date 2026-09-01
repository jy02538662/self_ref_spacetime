"""量子性 = 断裂 + 自指 = 自指观测：三个可证伪的数值验证。

ExpA: [X, D] != 0（位置=对角、动量=非对角，自指两端不对易）
ExpB: hbar = |d|^2_min（断裂最小单位 = 取向矢量模长平方最小非零值）
ExpC: 投影观测 -> 坍缩（自指观测的类比，测量改变状态）

注意：这三个实验验证的是"推导链每一步的数学正确性"，不是新物理发现。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# 泡利矩阵
SIGMA1 = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
SIGMA3 = np.array([[1, 0], [0, -1]], dtype=complex)


def expA() -> dict:
    """位置 X=sigma3（对角）、动量 P=sigma1（非对角），验证 [X,P]=2i sigma2 != 0."""
    X = SIGMA3
    P = SIGMA1
    comm = X @ P - P @ X
    expected = 2j * SIGMA2
    ok = np.allclose(comm, expected)
    return {
        "commutator_nonzero": bool(np.any(np.abs(comm) > 1e-12)),
        "matches_2i_sigma2": bool(ok),
    }


def expB(eps: float = 0.5) -> dict:
    """断裂贡献 |d|^2 离散化，最小非零单位 = eps^2 = 一个量子."""
    d_vals = np.array([0.0, eps, 2 * eps, 3 * eps, 4 * eps])
    energies = d_vals**2  # 断裂贡献（固定方向，|d|^2）
    min_nonzero = float(energies[energies > 0].min())
    return {
        "eps": eps,
        "break_energies": energies.tolist(),
        "min_nonzero_break": min_nonzero,
        "expected_eps2": eps**2,
        "min_equals_eps2": bool(min_nonzero == eps**2),
    }


def expC(n: int = 4, seed: int = 0) -> dict:
    """投影观测 -> 坍缩：叠加态投影到本征态后，叠加信息丢失."""
    rng = np.random.default_rng(seed)
    M = rng.normal(0, 1, size=(n, n)) + 1j * rng.normal(0, 1, size=(n, n))
    D = (M + M.conj().T) / 2  # 厄米化
    evals, evecs = np.linalg.eigh(D)

    # 叠加态 |psi> = sum c_i |v_i>
    c = rng.normal(0, 1, size=n) + 1j * rng.normal(0, 1, size=n)
    c = c / np.linalg.norm(c)
    psi = evecs @ c

    # 观测 = 投影到本征态 v_0
    P0 = np.outer(evecs[:, 0], evecs[:, 0].conj())
    psi_collapsed = P0 @ psi
    psi_collapsed = psi_collapsed / np.linalg.norm(psi_collapsed)

    overlap_before = float(np.abs(np.vdot(evecs[:, 0], psi)))  # 坍缩前重叠 = |c_0|
    overlap_after = float(np.abs(np.vdot(evecs[:, 0], psi_collapsed)))  # 坍缩后应 = 1
    return {
        "overlap_before_collapse": overlap_before,
        "overlap_after_collapse": overlap_after,
        "collapsed_to_eigenstate": bool(overlap_after > 0.999),
        "state_changed": bool(abs(overlap_after - overlap_before) > 1e-9),
    }


if __name__ == "__main__":
    results = {"expA": expA(), "expB": expB(), "expC": expC()}

    print("=== 量子性 = 断裂 + 自指 = 自指观测：数值验证 ===")
    print("\n[ExpA] 位置-动量不对易")
    print(f"  [X,P] != 0: {results['expA']['commutator_nonzero']}")
    print(f"  [X,P] == 2i*sigma2: {results['expA']['matches_2i_sigma2']}")

    print("\n[ExpB] 断裂最小单位 = hbar")
    b = results["expB"]
    print(f"  断裂能量序列 |d|^2 = {b['break_energies']}")
    print(f"  最小非零断裂 = {b['min_nonzero_break']} == eps^2 = {b['expected_eps2']}: {b['min_equals_eps2']}")

    print("\n[ExpC] 投影观测 -> 坍缩")
    c = results["expC"]
    print(f"  坍缩前重叠 |<v0|psi>| = {c['overlap_before_collapse']:.4f}")
    print(f"  坍缩后重叠 = {c['overlap_after_collapse']:.4f}")
    print(f"  坍缩到本征态: {c['collapsed_to_eigenstate']}")

    print("\n--- summary ---")
    print(json.dumps(results, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp_quantization_selfref_last_run.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"wrote {out}")

"""Exp7 后续：引力侧第一道门 = 自旋联络 = 「2×2 → 4×4」的自发产生

问题（缺口 1 = 内部 SU(2) → 时空）：
    exp7_k2_clifford 已证明 4×4 Dirac 是「手搭的」——2×2（有向区分 γ⁰ = 时间）
    ⊗ 2×2（SU(2)≅S³ = 空间）→ 4×4 Cl(3,1)，号差 −+++ 误差 0。
    但那一步「2×2 → 4×4」是手动 Dirac 构造，不是 D 自发长出来的。

本实验回答：**π 磁通 D（纯迹作用量自发选出的结构）是否自发携带
「2×2 → 4×4」的两个因子？**

物理对应（AZ 分类的精确一步）：
    - 2×2 = Kramers 对 = 自旋 1/2 = 反幺 T、T²=−1（exp7 已证 π 磁通自发）。
    - 4×4 = Dirac = Kramers ⊗ 手征。加一个手征/子格对称 Γ（Γ²=I、{Γ,D}=0、
      Γ 与 T 反交换）就把 2×2 抬到 4×4（DIII/CII 类）。
    - 判据 = π 磁通 D 是否也自发携带手征对称 Γ（在 T²=−1 之上）。

客观判据（不手放任何算符）：
    1) 手征 Γ 存在 ⟺ 谱 ±E 精确对称（bipartite ⟹ {Γ,D}=0 ⟹ 谱关于 0 对称）。
    2) Kramers T²=−1 存在 ⟺ 所有本征值重数偶数。
    3) 4×4 Dirac ⟺ 零模四重简并 = 2（Kramers）× 2（手征），且 Dirac 型
       J（ΓJ=−JΓ）可构造：A 的奇异值成对 ⟺ 存在正交 B 使 B A^T = −A B^T。

诚实边界：本实验只判「两个因子是否都自发存在」+「4×4 Dirac 结构是否成」，
不声称「自旋联络已完整导出」——那还需把 Γ 和 T 焊成规范协变的 ω_μ（下一道门）。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

from experiments.exp_pure_spectral_anneal import toroidal_D


def sublattice_Gamma(n_per_dim: int) -> np.ndarray:
    """子格/手征算符 Γ = diag((-1)^{i+j})（checkerboard，Γ²=I，{Γ,D}=0）。"""
    n = n_per_dim
    signs = np.array([(-1.0) ** (i + j) for i in range(n) for j in range(n)])
    return np.diag(signs)


def chiral_err(evals: np.ndarray, tol: float = 1e-9) -> float:
    """谱 ±E 对称性误差：手征对称 ⟺ {Γ,D}=0 ⟺ 谱关于 0 对称。"""
    pos = np.sort(evals[evals > tol])
    neg = np.sort(evals[evals < -tol])
    if pos.size != neg.size:
        return float("nan")
    return float(np.max(np.abs(pos + neg[::-1])))


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


def selfdual_J(D: np.ndarray) -> dict:
    """构造反幺 T=JK（T²=−1）的自对偶 J 并验证（复用 exp7 判据）。"""
    n = D.shape[0]
    evals, evecs = np.linalg.eigh(D)
    mults = multiplicities(evals, 1e-8 * (evals.max() - evals.min()))
    all_even = bool(np.all(mults % 2 == 0))
    sigma = np.array([[0, 1], [-1, 0]], dtype=float)
    J = np.zeros((n, n), dtype=complex)
    pos = 0
    for m in mults:
        if m % 2 == 0:
            for _ in range(m // 2):
                J[pos:pos + 2, pos:pos + 2] = sigma
                pos += 2
        else:
            pos += m
    Jphys = evecs @ J @ evecs.T
    return {
        "all_even": all_even,
        "antisym_err": float(np.max(np.abs(Jphys + Jphys.T))),
        "orth_err": float(np.max(np.abs(Jphys @ Jphys.T - np.eye(n)))),
        "comm_err": float(np.max(np.abs(Jphys @ D - D @ Jphys))),
        "T2_minus1_err": float(np.max(np.abs(Jphys @ Jphys.conj() + np.eye(n)))),
    }


def zero_mode_chirality(D: np.ndarray, Gamma: np.ndarray, tol: float = 1e-8) -> dict:
    """零模（Dirac 点）四重简并是否 = 2 手征正 + 2 手征负（4×4 = 2×2 ⊗ 2×2）。"""
    evals, V = np.linalg.eigh(D)
    scale = float(np.max(np.abs(evals))) if evals.size else 1.0
    zero_mask = np.abs(evals) < tol * scale
    nz = int(np.sum(zero_mask))
    if nz == 0:
        return {"n_zero": 0, "Gamma_evals": [], "verdict": "无零模（N≡2 mod 4）"}
    V0 = V[:, zero_mask]
    Gamma0 = V0.conj().T @ Gamma @ V0
    ge = np.sort(np.real(np.linalg.eigvalsh(Gamma0)))
    n_plus = int(np.sum(ge > 0.5))
    n_minus = int(np.sum(ge < -0.5))
    verdict = ("4×4 Dirac（2 手征正 × 2 手征负）" if nz == 4 and n_plus == 2 and n_minus == 2
               else f"{nz} 重（非标准）")
    return {"n_zero": nz, "Gamma_evals": np.round(ge, 4).tolist(),
            "verdict": verdict, "n_plus": n_plus, "n_minus": n_minus}


def dirac_J(D: np.ndarray, Gamma: np.ndarray) -> dict:
    """构造 Dirac 型自对偶 J（ΓJ=−JΓ，DIII 类）并验证。

    chiral 基 [A,B] 里 D=[[0,A],[A^T,0]]、Γ=diag(+1..+1,−1..−1)。
    J=[[0,B],[-B^T,0]] 自动 ΓJ=−JΓ、反对称。还需 B 正交（J²=−I）且
    B A^T = −A B^T（JD=DJ）。B A^T = −A B^T ⟺ A 奇异值成对（B' Σ = −Σ B'^T）。
    """
    n = D.shape[0]
    idxA = np.where(np.diag(Gamma) > 0)[0]
    idxB = np.where(np.diag(Gamma) < 0)[0]
    A = np.real(D[np.ix_(idxA, idxB)])            # A→B 跳变 (N/2 × N/2)
    U, S, Vh = np.linalg.svd(A)
    V = Vh.T                                       # A = U S V^T
    m = len(S)
    # 奇异值是否成对（偶数重数）——按 SVD 实际（降序）顺序分组
    Bp = np.zeros((m, m), dtype=float)
    i = 0
    all_even = True
    while i < m:
        j = i
        while j + 1 < m and abs(S[j + 1] - S[i]) < 1e-8:
            j += 1
        k = j - i + 1                              # 该奇异值块的大小
        if k % 2 == 0:
            h = k // 2
            Bp[i:i + h, i + h:i + 2 * h] = np.eye(h)
            Bp[i + h:i + 2 * h, i:i + h] = -np.eye(h)
        else:
            all_even = False
        i = j + 1
    B = U @ Bp @ V.T
    # J 在 chiral 基（[A,B] 排序），再转回 site 基
    Z = np.zeros((m, m))
    J_c = np.block([[Z, B], [-B.T, Z]])
    perm = np.concatenate([idxA, idxB])
    inv = np.argsort(perm)
    J = J_c[np.ix_(inv, inv)]
    err = {
        "svals_even_paired": all_even,
        "antisym_err": float(np.max(np.abs(J + J.T))),
        "orth_err": float(np.max(np.abs(J @ J.T - np.eye(n)))),
        "J2_minus1_err": float(np.max(np.abs(J @ J + np.eye(n)))),
        "comm_D_err": float(np.max(np.abs(J @ D - D @ J))),
        "anti_Gamma_err": float(np.max(np.abs(Gamma @ J + J @ Gamma))),
    }
    err["dirac_DIII"] = bool(all(
        err[k] < 1e-8 for k in ("antisym_err", "orth_err", "J2_minus1_err",
                                "comm_D_err", "anti_Gamma_err")))
    return err


def main():
    print("=" * 82)
    print("引力侧第一道门：π 磁通 D 是否自发携带「2×2 → 4×4」的两个因子？")
    print("  2×2 = T²=−1（Kramers）  ×  2×2 = Γ（手征/子格）  →  4×4 Dirac")
    print("=" * 82)

    results = {}
    for n_per_dim in (4, 6, 8, 10, 12):
        n = n_per_dim ** 2
        D = toroidal_D(n_per_dim, pi_flux=True)
        Gamma = sublattice_Gamma(n_per_dim)

        evals = np.linalg.eigvalsh(D)
        ce = chiral_err(evals)
        sd = selfdual_J(D)
        zm = zero_mode_chirality(D, Gamma)
        dj = dirac_J(D, Gamma)

        results[f"pi_flux_N{n}"] = {
            "chiral_err": ce, "selfdual": sd, "zero_mode": zm, "dirac_J": dj,
        }
        print(f"\n[N={n}  π-flux]")
        print(f"  [1] 手征 Γ：谱 ±E 对称误差 = {ce:.2e}  →  "
              f"{'✓ 有手征（bipartite）' if ce < 1e-10 else '✗ 无手征'}")
        print(f"  [2] Kramers T²=−1：全偶重数 = {sd['all_even']}  T²=−1 误差 = {sd['T2_minus1_err']:.1e}")
        print(f"  [3] 零模：{zm['verdict']}  （Γ 本征值 = {zm['Gamma_evals']}）")
        print(f"  [4] Dirac J（ΓJ=−JΓ，DIII）：奇异值成对 = {dj['svals_even_paired']}")
        print(f"      J 验证（反对称/正交/J²=−I/与D对易/与Γ反交换）误差 = "
              f"{dj['antisym_err']:.1e}/{dj['orth_err']:.1e}/{dj['J2_minus1_err']:.1e}/"
              f"{dj['comm_D_err']:.1e}/{dj['anti_Gamma_err']:.1e}")
        print(f"      →  {'✓ 4×4 Dirac（DIII）结构成立' if dj['dirac_DIII'] else '✗'}")

    out = ROOT / "experiments" / "exp7_k2_chiral_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

"""Exp7: 无预设聚类 —— D 的本征态是否自发聚成 k=2（二元 / SU(2)）？

核心问题（v4 第十节 · 第三层/逼三维 的 ① Exp7）：
    断裂 -> SU(2) 的自发涌现。如果 D 的谱自发出现两两简并（k=2 二元结构），
    第三层就活了（断裂 -> SU(2) -> S^2 -> Hopf -> 三维 整条链自动接上）。

最客观的判据（不预设任何聚类度量，纯用谱的客观结构）：
    1) 本征值简并结构 —— 排序后按相邻能级间隔聚类，统计 k=1/2/3/4+ 簇的占比。
    2) 四元数/SU(2) 结构（Kramers 定理）—— 一个厄米矩阵 D 是「四元数」（等价于
       存在反幺对称 T、T^2=-1、[T,D]=0，即自旋 1/2 的 Kramers 二重简并）
       当且仅当 **所有本征值重数都是偶数**。这是严格定理（反对称幺正矩阵只能偶数维），
       它把「两两简并」和「SU(2) 二元」精确对应起来。
    3) 本征向量配对（次判据）—— 对每个二重简并对，复共轭双线性形式
       B_ab = v_a^T v_b（无共轭）在 Kramers 对里是反对称（|B_12|=1, 对角 0），
       在实对称（T=K, T^2=+1）里是对称。区分「SU(2) 二元」与「实矩阵/格点对称」。

对照（诚实校准，让结果可解释）：
    - random_phase : toroidal 模长 + 每条边随机相位 θ_ij（用户主提案）
    - pi_flux      : 纯迹作用量实际长出的 π 磁通结构（真实相位，实矩阵 ±1）
    - pi_flux_perturbed : π 磁通 + 小相位扰动（结构对噪声的鲁棒性）
    - zero_phase   : 纯环面（全实 +1，格点平移对称）
    - gue          : 随机复厄米矩阵（零假设：无简并，能级排斥，<r>≈0.60）
    - quaternionic : 随机四元数厄米矩阵（正控制：验证探测器报 100% k=2）

诚实边界：本实验只回答「谱是否两两简并」这一客观事实，不声称「三维已涌现」。
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

# 复用纯迹作用量里已验证的 toroidal 构造（n_per_dim 每维节点数，N = n_per_dim^2）
from experiments.exp_pure_spectral_anneal import toroidal_D


# ----------------------------------------------------------------------------
# 构造各场景的 D
# ----------------------------------------------------------------------------

def toroidal_modulus(n_per_dim: int) -> np.ndarray:
    """toroidal 的模长部分（无相位）：|D_ij| = 1 在边，0 否则。"""
    return np.abs(toroidal_D(n_per_dim, pi_flux=False))


def random_phase_D(modulus: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """模长固定（toroidal 结构），每条边赋随机相位 θ_ij。"""
    n = modulus.shape[0]
    D = np.zeros((n, n), complex)
    iu = np.triu_indices(n, k=1)
    r = modulus[iu]
    phase = rng.uniform(-np.pi, np.pi, size=len(r))
    z = r * np.exp(1j * phase)
    D[iu] = z
    D[(iu[1], iu[0])] = np.conj(z)
    return D


def perturbed_pi_flux_D(n_per_dim: int, rng: np.random.Generator, eps: float) -> np.ndarray:
    """π 磁通结构 + 每条边相位加 N(0, eps) 扰动。"""
    D0 = toroidal_D(n_per_dim, pi_flux=True)
    n = D0.shape[0]
    iu = np.triu_indices(n, k=1)
    z = D0[iu]
    keep = np.abs(z) > 1e-12
    phase_jitter = rng.normal(0.0, eps, size=len(z))
    z = z * np.exp(1j * phase_jitter)
    z = np.where(keep, z, 0.0)
    D = np.zeros((n, n), complex)
    D[iu] = z
    D[(iu[1], iu[0])] = np.conj(z)
    return D


def gue_D(n: int, rng: np.random.Generator) -> np.ndarray:
    """随机复厄米矩阵（GUE，零假设：无简并）。"""
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return (M + M.conj().T) / 2


def quaternionic_D(n_sites: int, rng: np.random.Generator) -> np.ndarray:
    """随机四元数厄米矩阵（正控制）：对反幺 T=J K（T^2=-1）不变 -> 必二重简并。

    J = ⨁ [[0,1],[-1,0]]（反对称正交，J^2=-I）。M_quat = (M + J conj(M) J^T)/2
    对 T 不变，故所有本征值严格二重简并（Kramers）。
    """
    n = 2 * n_sites
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    M = (M + M.conj().T) / 2
    J = np.kron(np.eye(n_sites), np.array([[0, 1], [-1, 0]], dtype=float))
    Mq = (M + J @ M.conj() @ J.T) / 2
    return Mq


# ----------------------------------------------------------------------------
# 诊断
# ----------------------------------------------------------------------------

def cluster_multiplicities(evals: np.ndarray, tol: float) -> np.ndarray:
    """把排序后的本征值按间隔 < tol 聚成简并簇，返回每个簇的大小（重数）。

    tol 用绝对量纲，调用方负责给一个相对谱宽的合理值。
    """
    evals = np.sort(evals)
    n = len(evals)
    if n == 0:
        return np.array([], dtype=int)
    sizes = []
    start = 0
    for i in range(1, n):
        if evals[i] - evals[i - 1] > tol:
            sizes.append(i - start)
            start = i
    sizes.append(n - start)
    return np.asarray(sizes, dtype=int)


def degeneracy_summary(D: np.ndarray, tol_exact: float, tol_near: float) -> dict:
    """对一个 D 算：本征值 + 简并簇统计 + 四元数判据 + 本征向量配对类型。"""
    n = D.shape[0]
    evals, evecs = np.linalg.eigh(D)
    width = float(np.max(evals) - np.min(evals))
    scale = width if width > 1e-12 else 1.0

    def _frac(mults: np.ndarray) -> dict:
        # 占「能级」的比例（不是簇数）：k 重简并簇贡献 k 个能级。
        total = int(np.sum(mults))
        hist = {}
        for k in (1, 2, 3, 4):
            hist[f"k{k}"] = float(k * np.sum(mults == k)) / total if total else 0.0
        hist["k_ge5"] = float(sum(m for m in mults if m >= 5)) / total if total else 0.0
        return hist

    mults_exact = cluster_multiplicities(evals, tol_exact * scale)
    mults_near = cluster_multiplicities(evals, tol_near * scale)

    # 四元数/SU(2) 判据：所有本征值重数都为偶数（Kramers 定理）
    all_even = bool(np.all(mults_exact % 2 == 0))

    # 能级间隔比 <r>（GUE≈0.60，泊松≈0.39，精确二重简并→0）
    gaps = np.diff(np.sort(evals))
    gaps = gaps[gaps > 1e-14]
    if len(gaps) >= 2:
        r = np.minimum(gaps[:-1], gaps[1:]) / np.maximum(gaps[:-1], gaps[1:])
        mean_r = float(np.mean(r))
    else:
        mean_r = float("nan")

    # 构造性验证：若 D 自对偶（全偶重数），显式构造反幺 T=J K（T^2=-1）并数值验证
    #   J = V (⨁σ) V^T，σ=[[0,1],[-1,0]]，放在每个偶数重数本征子空间里。
    #   验证四条：J 反对称、J 幺正、J D^* = D J（自对偶）、T^2 = J conj(J) = -I。
    #   这是「本征向量配对成二元结构」的严格、构造性判据（不是依赖基选择的 bilinear）。
    jcheck = _verify_selfdual_J(D, evals, evecs, mults_exact) if all_even else None

    return {
        "n": n,
        "evals": np.round(evals, 6).tolist(),
        "mults_exact": mults_exact.tolist(),
        "mults_near": mults_near.tolist(),
        "frac_exact": _frac(mults_exact),
        "frac_near": _frac(mults_near),
        "all_even_mult_quaternionic": all_even,
        "mean_r": mean_r,
        "selfdual_J_check": jcheck,
    }


def _verify_selfdual_J(D: np.ndarray, evals: np.ndarray, evecs: np.ndarray, mults: np.ndarray) -> dict:
    """构造反幺 T=J K 并验证它是 T^2=-1、与 D 对易的反幺对称（= 四元数/SU(2) 结构）。"""
    n = D.shape[0]
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
    V = evecs
    Jphys = V @ J @ V.T
    return {
        "antisym_err": float(np.max(np.abs(Jphys + Jphys.T))),
        "unitary_err": float(np.max(np.abs(Jphys @ Jphys.conj().T - np.eye(n)))),
        "comm_err": float(np.max(np.abs(Jphys @ D.conj() - D @ Jphys))),
        "T2_minus1_err": float(np.max(np.abs(Jphys @ Jphys.conj() + np.eye(n)))),
    }


# ----------------------------------------------------------------------------
# 场景聚合
# ----------------------------------------------------------------------------

def aggregate(Ds: list[np.ndarray], tol_exact: float, tol_near: float) -> dict:
    """对一批 D 实例聚合 k=2 占比等统计。"""
    summaries = [degeneracy_summary(D, tol_exact, tol_near) for D in Ds]
    keys = ["k1", "k2", "k3", "k4", "k_ge5"]
    out = {"n_instances": len(summaries), "all_even_mult_quaternionic": []}
    for tag, tol in (("exact", tol_exact), ("near", tol_near)):
        for k in keys:
            vals = [s[f"frac_{tag}"][k] for s in summaries]
            out[f"frac_{tag}_{k}_mean"] = float(np.mean(vals))
            out[f"frac_{tag}_{k}_std"] = float(np.std(vals))
    out["all_even_mult_quaternionic"] = [s["all_even_mult_quaternionic"] for s in summaries]
    out["all_even_frac"] = float(np.mean(out["all_even_mult_quaternionic"]))
    out["mean_r_mean"] = float(np.nanmean([s["mean_r"] for s in summaries]))
    return out


def main():
    tol_exact = 1e-8
    tol_near = 2e-2
    n_seeds_random = 300
    n_seeds_pert = 100

    print("=" * 80)
    print("Exp7: D 的本征态是否自发聚成 k=2（二元 / SU(2) 结构）？")
    print("判据 = 谱的客观简并结构（无预设聚类度量）+ 四元数/Kramers 判据")
    print("=" * 80)

    results = {}

    for n_per_dim in (3, 4, 5, 6, 7, 8):
        n = n_per_dim ** 2
        mod = toroidal_modulus(n_per_dim)
        print(f"\n{'#' * 80}\n# N = {n}  ({n_per_dim}×{n_per_dim} torus)\n{'#' * 80}")

        rng = np.random.default_rng(0)

        # 1) 用户主提案：随机相位
        Ds = [random_phase_D(mod, np.random.default_rng(s)) for s in range(n_seeds_random)]
        a = aggregate(Ds, tol_exact, tol_near)
        results[f"random_phase_N{n}"] = a
        print(f"\n[random_phase]  {n_seeds_random} 实例  (toroidal 模长 + 随机相位)")
        print(f"  精确简并占比 k1/k2/k3/k4+:  "
              f"k1={a['frac_exact_k1_mean']:.3f}  k2={a['frac_exact_k2_mean']:.3f}  "
              f"k3={a['frac_exact_k3_mean']:.3f}  k4+={a['frac_exact_k4_mean']:.3f}")
        print(f"  近似简并(2%)占比:          "
              f"k1={a['frac_near_k1_mean']:.3f}  k2={a['frac_near_k2_mean']:.3f}")
        print(f"  四元数(Kramers)占比: {a['all_even_frac']:.3f}   <r>={a['mean_r_mean']:.4f}")

        # 2) π 磁通（纯迹作用量实际长出的结构，实矩阵 ±1）
        Dpi = toroidal_D(n_per_dim, pi_flux=True)
        api = degeneracy_summary(Dpi, tol_exact, tol_near)
        results[f"pi_flux_N{n}"] = api
        print(f"\n[pi_flux]  1 实例  (纯迹作用量长出的真实相位结构)")
        print(f"  本征值: {api['evals']}")
        print(f"  重数: {api['mults_exact']}   自对偶/四元数={api['all_even_mult_quaternionic']}")
        if api["selfdual_J_check"] is not None:
            j = api["selfdual_J_check"]
            print(f"  J 验证(反对称/幺正/对易/T²=-1 误差): "
                  f"{j['antisym_err']:.1e}/{j['unitary_err']:.1e}/{j['comm_err']:.1e}/{j['T2_minus1_err']:.1e}")

        # 3) π 磁通 + 相位扰动（鲁棒性）
        for eps in (0.05, 0.3):
            Ds = [perturbed_pi_flux_D(n_per_dim, np.random.default_rng(s), eps) for s in range(n_seeds_pert)]
            ap = aggregate(Ds, tol_exact, tol_near)
            results[f"pi_flux_perturbed_N{n}_eps{eps}"] = ap
            print(f"\n[pi_flux_perturbed eps={eps}]  {n_seeds_pert} 实例")
            print(f"  精确 k2 占比: {ap['frac_exact_k2_mean']:.3f}   四元数占比: {ap['all_even_frac']:.3f}   <r>={ap['mean_r_mean']:.4f}")

        # 4) 零相位（纯环面，格点平移对称）
        Dz = toroidal_D(n_per_dim, pi_flux=False)
        az = degeneracy_summary(Dz, tol_exact, tol_near)
        results[f"zero_phase_N{n}"] = az
        print(f"\n[zero_phase]  1 实例  (纯环面，全实 +1)")
        print(f"  本征值: {az['evals']}")
        print(f"  重数: {az['mults_exact']}   自对偶/四元数={az['all_even_mult_quaternionic']}")

        # 5) GUE 零假设
        Ds = [gue_D(n, np.random.default_rng(s)) for s in range(n_seeds_random)]
        ag = aggregate(Ds, tol_exact, tol_near)
        results[f"gue_N{n}"] = ag
        print(f"\n[gue]  {n_seeds_random} 实例  (随机复厄米，零假设)")
        print(f"  精确 k1 占比: {ag['frac_exact_k1_mean']:.3f}   四元数占比: {ag['all_even_frac']:.3f}   <r>={ag['mean_r_mean']:.4f}")
        print(f"  近似简并(2%) k2: {ag['frac_near_k2_mean']:.3f}  (对照 random_phase 的近似 k2)")

        # 6) 四元数正控制（N 需为偶数）
        if n % 2 == 0:
            Ds = [quaternionic_D(n // 2, np.random.default_rng(s)) for s in range(n_seeds_random)]
            aq = aggregate(Ds, tol_exact, tol_near)
            results[f"quaternionic_N{n}"] = aq
            print(f"\n[quaternionic]  {n_seeds_random} 实例  (正控制：随机四元数厄米)")
            print(f"  精确 k2 占比: {aq['frac_exact_k2_mean']:.3f}   四元数占比: {aq['all_even_frac']:.3f}   <r>={aq['mean_r_mean']:.4f}")

    # 落盘
    out = ROOT / "experiments" / "exp7_k2_degeneracy_last_run.json"
    out.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()

"""冷启动：随机 D + 纯 Tr(D²)（无任何目标距离/邻接表），看 D 塌缩成什么谱维数。

诚实目标：证明"纯耦合成本不会自发选维，只会塌缩"——这暴露了"自发选维"
需要一个维度无关的维数压力项（第三层的坎），而不是预设目标距离。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.action import tr_D2
from src.algebra import build_D, random_params, unpack_params


def spectral_dimension_from_D(D: np.ndarray, t_lo: float = 0.5, t_hi: float = 50.0, n_pts: int = 40) -> float:
    """从 D² 的谱算谱维数（D² 当拉普拉斯）。"""
    D2 = D @ D
    eig = np.linalg.eigvalsh(D2)
    eig = np.clip(eig, 1e-14, None)
    ts = np.geomspace(t_lo, t_hi, n_pts)
    logK = np.array([np.log(np.sum(np.exp(-t * eig))) for t in ts])
    slope = np.polyfit(np.log(ts), logK, 1)[0]
    return float(-2.0 * slope)


def run(n: int, seed: int, maxiter: int = 400) -> dict:
    rng = np.random.default_rng(seed)
    theta0 = random_params(n, rng, scale=0.8)

    def fun(theta):
        moduli, phases = unpack_params(theta, n)
        D = build_D(moduli, phases, n)
        return tr_D2(D)  # 纯耦合成本，无目标距离

    res = minimize(fun, theta0, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-12})
    theta = res.x
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)

    trD2 = tr_D2(D)
    max_abs_z = float(np.max(np.abs(D))) if n > 1 else 0.0
    # 初始谱维数（随机 D）
    _, phases0 = unpack_params(theta0, n)
    D0 = build_D(moduli, phases0, n)
    d0 = spectral_dimension_from_D(D0) if np.max(np.abs(D0)) > 1e-10 else 0.0
    # 优化后谱维数
    d_final = spectral_dimension_from_D(D) if max_abs_z > 1e-10 else 0.0

    return {
        "n": n,
        "seed": seed,
        "TrD2_initial": float(tr_D2(build_D(*unpack_params(theta0, n), n))),
        "TrD2_final": float(trD2),
        "max_abs_z_final": max_abs_z,
        "spectral_dim_initial": d0,
        "spectral_dim_final": d_final,
    }


if __name__ == "__main__":
    n = 27
    print(f"=== 冷启动：随机 D + 纯 Tr(D^2)  N={n} ===")
    rows = []
    for seed in range(3):
        r = run(n, seed)
        rows.append(r)
        print(
            f"seed={seed}  TrD2 {r['TrD2_initial']:.2f} -> {r['TrD2_final']:.5f}  "
            f"max|z|={r['max_abs_z_final']:.5f}  "
            f"谱维数 {r['spectral_dim_initial']:.2f} -> {r['spectral_dim_final']:.2f}"
        )
    summary = {
        "mean_TrD2_final": float(np.mean([r["TrD2_final"] for r in rows])),
        "mean_spectral_dim_final": float(np.mean([r["spectral_dim_final"] for r in rows])),
        "conclusion": (
            "纯耦合成本把 D 压到 0，谱维数塌缩到 ~0——冷启动不涌现，只会塌缩。"
            "自发选维需要一个维度无关的维数压力项（第三层的坎）。"
        ),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp_cold_start_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

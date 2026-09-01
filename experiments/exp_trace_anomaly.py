"""Trace anomaly: does the one-loop fluctuation Tr ln M acquire a non-trivial
log-lambda term under D -> lambda D?

S = (1/2) Tr(D^2) + (g/4) Tr(D^4),  D real symmetric.

One-loop effective action:  Gamma = S + (1/2) Tr ln M,   M = Hessian of S.

Under D -> lambda D:
    Tr ln M[lambda D] = Tr ln(M2 + g lambda^2 M4) + 2m log lambda,
    m = N(N+1)/2,  M2 = Hessian of (1/2)Tr(D^2) (constant),  M4 = Hessian of (1/4)Tr(D^4).

The trivial measure term is 2m log lambda.  The TRACE ANOMALY is any residual
log-lambda dependence of Tr ln(M2 + g lambda^2 M4).

Falsifiable: g=0 -> no anomaly (residual flat);  g!=0 -> anomaly (residual has log lambda).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def unpack_real_sym(x: np.ndarray, n: int) -> np.ndarray:
    D = np.zeros((n, n))
    idx = 0
    for i in range(n):
        for j in range(i, n):
            D[i, j] = x[idx]
            D[j, i] = x[idx]
            idx += 1
    return D


def pack_real_sym(D: np.ndarray, n: int) -> np.ndarray:
    return np.array([D[i, j] for i in range(n) for j in range(i, n)])


def S_val(D: np.ndarray, g: float) -> float:
    D2 = D @ D
    return 0.5 * np.trace(D2) + (g / 4.0) * np.trace(D2 @ D2)


def hessian(fun, x: np.ndarray, h: float = 1e-3) -> np.ndarray:
    m = len(x)
    H = np.zeros((m, m))
    for a in range(m):
        for b in range(a, m):
            ea = np.zeros(m)
            ea[a] = h
            eb = np.zeros(m)
            eb[b] = h
            fpp = fun(x + ea + eb)
            fpm = fun(x + ea - eb)
            fmp = fun(x - ea + eb)
            fmm = fun(x - ea - eb)
            H[a, b] = (fpp - fpm - fmp + fmm) / (4 * h * h)
            H[b, a] = H[a, b]
    return H


def tr_log_M(x: np.ndarray, n: int, g: float, h: float = 1e-3) -> tuple[float, float]:
    fun = lambda xx: S_val(unpack_real_sym(xx, n), g)
    H = hessian(fun, x, h)
    eig = np.linalg.eigvalsh(H)
    return float(np.sum(np.log(np.maximum(eig, 1e-300)))), float(np.min(eig))


def run(n: int, gs: list[float], lam: np.ndarray, seed: int = 0) -> list[dict]:
    m = n * (n + 1) // 2
    rng = np.random.default_rng(seed)
    xbar = pack_real_sym(unpack_real_sym(rng.normal(0.0, 1.0, size=m), n), n)

    # Sanity: g=0 Hessian of (1/2)Tr(D^2) is diagonal with entries {1 (diag), 2 (offdiag)}
    # so Tr ln M2 = (m - n) log 2.
    analytic_g0 = (m - n) * np.log(2.0)

    rows: list[dict] = []
    for g in gs:
        for l in lam:
            tl, mineig = tr_log_M(l * xbar, n, g)
            resid = tl - 2 * m * np.log(l)
            rows.append(
                {
                    "g": g,
                    "lambda": float(l),
                    "ln_lambda": float(np.log(l)),
                    "Tr_ln_M": tl,
                    "resid_minus_2m_lnl": resid,
                    "min_eig": mineig,
                }
            )
            print(
                f"g={g:4.1f}  λ={l:5.2f}  lnλ={np.log(l):+.3f}  "
                f"TrlnM={tl:8.3f}  resid={resid:8.3f}  min_eig={mineig:.2e}"
            )
    return rows, analytic_g0


def summarize(rows: list[dict], analytic_g0: float) -> dict:
    # residual slope vs ln(lambda) per g: flat => no anomaly, non-zero => anomaly
    by_g: dict[float, dict] = {}
    for g in sorted({r["g"] for r in rows}):
        sub = sorted([r for r in rows if r["g"] == g], key=lambda r: r["lambda"])
        x = np.array([r["ln_lambda"] for r in sub])
        y = np.array([r["resid_minus_2m_lnl"] for r in sub])
        # slope of residual vs ln lambda (linear fit)
        slope = float(np.polyfit(x, y, 1)[0])
        by_g[g] = {
            "resid_spread": float(np.max(y) - np.min(y)),
            "resid_slope_vs_lnl": slope,
        }
    return {
        "analytic_g0_TrlnM2": analytic_g0,
        "by_g": by_g,
        "anomaly_present_for_g_nonzero": any(
            abs(v["resid_slope_vs_lnl"]) > 0.5 for g, v in by_g.items() if g > 0
        ),
        "g0_flat": abs(by_g.get(0.0, {}).get("resid_slope_vs_lnl", 99.0)) < 0.1,
        "criterion": (
            "g=0 residual flat (no anomaly); g!=0 residual has non-zero log-lambda slope (anomaly)"
        ),
    }


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Trace anomaly of D under scale")
    p.add_argument("--n", type=int, default=4)
    args = p.parse_args()

    n = args.n
    lam = np.array([0.5, 0.7, 1.0, 1.4, 2.0])
    gs = [0.0, 0.3, 1.0]
    print(f"=== Trace anomaly  N={n} (m={n*(n+1)//2} params) ===")
    rows, analytic_g0 = run(n, gs, lam)
    summary = summarize(rows, analytic_g0)
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp_trace_anomaly_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")

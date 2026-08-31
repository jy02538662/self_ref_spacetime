"""Exp3: Dirac structure gives Lorentzian signature (indefinite D²) on a 1+1 lattice.

D = γ⁰ ⊗ D_t + γ¹ ⊗ D_x, with γ⁰ antisymmetric ((γ⁰)²=−1) and γ¹ symmetric ((γ¹)²=+1).
Then D² = (γ⁰)²⊗D_t² + (γ¹)²⊗D_x² = −D_t² + D_x²  (cross terms cancel by γ anti-commutation).

Compare against the Euclidean version (both γ symmetric), where D² = +D_t² + D_x² ≥ 0.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def hopping_1d(n):
    """Symmetric nearest-neighbor hopping on a 1D chain (open boundary)."""
    H = np.zeros((n, n))
    for i in range(n - 1):
        H[i, i + 1] = 1.0
        H[i + 1, i] = 1.0
    return H


def build_D(Nt, Nx, lorentzian=True):
    """D = γ⁰⊗D_t⊗I_x + γ¹⊗I_t⊗D_x, 2 Dirac components per site.

    lorentzian=True : γ⁰ antisymmetric → (γ⁰)²=−I → time signature −
    lorentzian=False: γ⁰ symmetric     → (γ⁰)²=+I → Euclidean (all +)
    """
    if lorentzian:
        g0 = np.array([[0.0, 1.0], [-1.0, 0.0]])  # (γ⁰)² = −I
    else:
        g0 = np.array([[0.0, 1.0], [1.0, 0.0]])   # (γ⁰)² = +I (Euclidean)
    g1 = np.array([[0.0, 1.0], [1.0, 0.0]])        # (γ¹)² = +I
    Dt = hopping_1d(Nt)
    Dx = hopping_1d(Nx)
    Ix = np.eye(Nx)
    It = np.eye(Nt)
    return np.kron(g0, np.kron(Dt, Ix)) + np.kron(g1, np.kron(It, Dx))


def report(D, label):
    D2 = D @ D
    ev = np.real(np.linalg.eigvals(D2))
    pos = ev[ev > 1e-9]
    neg = ev[ev < -1e-9]
    zero = ev[np.abs(ev) <= 1e-9]
    r = {
        "label": label,
        "D_hermitian": bool(np.allclose(D, D.conj().T)),
        "eig_D2_min": float(np.min(ev)),
        "eig_D2_max": float(np.max(ev)),
        "n_pos": int(pos.size),
        "n_neg": int(neg.size),
        "n_zero": int(zero.size),
        "indefinite": bool(neg.size > 0 and pos.size > 0),
    }
    print(
        f"{label:12s}  hermitian={str(r['D_hermitian']):5s}  "
        f"eig(D^2) in [{r['eig_D2_min']:.3f}, {r['eig_D2_max']:.3f}]  "
        f"pos={r['n_pos']:3d} neg={r['n_neg']:3d} zero={r['n_zero']:3d}  "
        f"indefinite={r['indefinite']}"
    )
    return r


def main(Nt=8, Nx=8):
    print(f"=== Exp3 directed signature  1+1 lattice  Nt={Nt} Nx={Nx} ===")
    rows = [
        report(build_D(Nt, Nx, lorentzian=False), "Euclidean"),
        report(build_D(Nt, Nx, lorentzian=True), "Lorentzian"),
    ]
    summary = {
        "lorentzian_indefinite": rows[1]["indefinite"],
        "euclidean_indefinite": rows[0]["indefinite"],
        "pass": rows[1]["indefinite"] and not rows[0]["indefinite"],
        "criterion": "Lorentzian D^2 is indefinite (pos+neg) while Euclidean D^2 is not",
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    out = ROOT / "experiments" / "exp3_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp3 Dirac structure / Lorentzian signature")
    p.add_argument("--nt", type=int, default=8)
    p.add_argument("--nx", type=int, default=8)
    args = p.parse_args()
    main(Nt=args.nt, Nx=args.nx)

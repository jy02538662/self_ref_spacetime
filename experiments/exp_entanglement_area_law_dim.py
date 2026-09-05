"""Entanglement entropy area law in 2D and 3D: S ~ boundary = L^(d-1), not volume L^d.

The 1D version (exp_entanglement_area_law.py) shows local (Fermi sea) -> sub-volume,
random -> volume, but 1D's "boundary" is only 2 points, so the area law degenerates to a
logarithm.  This script goes to 2D and 3D, where the area law is clean:

    ground state (Fermi sea, = Gibbs at beta -> infinity) of the hopping Hamiltonian
    on a d-dim hypercubic lattice; subregion A = a corner cube of side l;
    S(rho_A) ~ l^(d-1)   (perimeter in 2D, area in 3D)  --  NOT l^d (volume).

The scaling exponent (d-1) is the DIMENSION, so this is where the quantum leg
(entanglement) and the geometric leg (dimension) meet in one computation.

Honest boundary:
  - Kinematic entanglement: the lattice (2D grid / 3D cube) is a preset "space", and the
    state is the Fermi sea (ground state), which is D's natural state but not derived from
    a self-reference D dynamics (that dynamics is still missing).
  - Demonstrates the AREA LAW and how the exponent gives dimension; does NOT show that
    self-reference D spontaneously chooses 3D (the lattice is preset).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def entropy_from_corr(C_A: np.ndarray) -> float:
    """Peschel: S = -sum [lam ln lam + (1-lam) ln(1-lam)] over eigenvalues of C_A."""
    lam = np.linalg.eigvalsh(C_A)
    lam = np.clip(lam, 1e-15, 1.0 - 1e-15)
    return float(-np.sum(lam * np.log(lam) + (1.0 - lam) * np.log(1.0 - lam)))


def hopping_2d(L: int) -> np.ndarray:
    """Nearest-neighbour hopping on an L x L grid (open BC)."""
    N = L * L
    H = np.zeros((N, N))
    for i in range(L):
        for j in range(L):
            idx = i * L + j
            if i + 1 < L:
                H[idx, (i + 1) * L + j] = H[(i + 1) * L + j, idx] = -1.0
            if j + 1 < L:
                H[idx, i * L + (j + 1)] = H[i * L + (j + 1), idx] = -1.0
    return H


def hopping_3d(L: int) -> np.ndarray:
    """Nearest-neighbour hopping on an L x L x L cube (open BC)."""
    N = L * L * L
    H = np.zeros((N, N))
    for i in range(L):
        for j in range(L):
            for k in range(L):
                idx = (i * L + j) * L + k
                if i + 1 < L:
                    H[idx, ((i + 1) * L + j) * L + k] = H[((i + 1) * L + j) * L + k, idx] = -1.0
                if j + 1 < L:
                    H[idx, (i * L + (j + 1)) * L + k] = H[(i * L + (j + 1)) * L + k, idx] = -1.0
                if k + 1 < L:
                    H[idx, (i * L + j) * L + (k + 1)] = H[(i * L + j) * L + (k + 1), idx] = -1.0
    return H


def fermi_sea_corr(H: np.ndarray, fill: float = 0.5) -> np.ndarray:
    """Correlation matrix of the Fermi sea: occupy the lowest fill*N eigenstates."""
    N = H.shape[0]
    _, eigvecs = np.linalg.eigh(H)
    n_occ = int(round(fill * N))
    occ = eigvecs[:, :n_occ]
    return occ @ occ.T


def subregion_indices_2d(L: int, l: int):
    """Indices of the l x l corner subregion on an L x L grid."""
    return [i * L + j for i in range(l) for j in range(l)]


def subregion_indices_3d(L: int, l: int):
    """Indices of the l x l x l corner subregion on an L x L x L cube."""
    return [(i * L + j) * L + k for i in range(l) for j in range(l) for k in range(l)]


if __name__ == "__main__":
    print("=== Entanglement area law: 2D and 3D Fermi sea (ground state) ===")
    print()

    # ---- 2D ----
    L2 = 40
    H2 = hopping_2d(L2)
    C2 = fermi_sea_corr(H2)
    print(f"2D grid {L2}x{L2}, Fermi sea (half filling), subregion = l x l corner")
    print("  l      S      S/l (perimeter-law if ~const)      S/l^2 (area if ~const)")
    rows2 = []
    for l in [2, 4, 6, 8, 10, 12]:
        idx = subregion_indices_2d(L2, l)
        S = entropy_from_corr(C2[np.ix_(idx, idx)])
        rows2.append((l, S, S / l, S / (l * l)))
        print(f"  {l:2d}   {S:7.3f}      {S / l:7.3f}                    {S / (l * l):7.3f}")
    print()

    # ---- 3D ----
    L3 = 12
    H3 = hopping_3d(L3)
    C3 = fermi_sea_corr(H3)
    print(f"3D cube {L3}^3, Fermi sea (half filling), subregion = l x l x l corner")
    print("  l      S      S/l^2 (area-law if ~const)      S/l^3 (volume if ~const)")
    rows3 = []
    for l in [2, 3, 4, 5, 6]:
        idx = subregion_indices_3d(L3, l)
        S = entropy_from_corr(C3[np.ix_(idx, idx)])
        rows3.append((l, S, S / (l * l), S / (l * l * l)))
        print(f"  {l:2d}   {S:7.3f}      {S / (l * l):7.3f}                    {S / (l * l * l):7.3f}")
    print()

    # ---- verdict ----
    l2 = np.array([r[0] for r in rows2], dtype=float)
    S2 = np.array([r[1] for r in rows2])
    l3 = np.array([r[0] for r in rows3], dtype=float)
    S3 = np.array([r[1] for r in rows3])

    # 2D: is S linear in l (perimeter) rather than l^2 (area)?  fit log S vs log l
    slope2 = float(np.polyfit(np.log(l2), np.log(S2), 1)[0])
    slope3 = float(np.polyfit(np.log(l3), np.log(S3), 1)[0])
    print(f"scaling exponent: 2D  d(S) ~ l^{slope2:.2f}   (perimeter=1, area=2)")
    print(f"                  3D  d(S) ~ l^{slope3:.2f}   (area=2, volume=3)")
    print()
    area_law_2d = 0.7 < slope2 < 1.6  # between perimeter and (a bit above) -> boundary scaling
    area_law_3d = 1.7 < slope3 < 2.6
    print(f"verdict: 2D boundary-scaling = {area_law_2d}, 3D boundary-scaling = {area_law_3d}")
    print()

    summary = {
        "2D": {"L": L2, "rows": [[int(r[0]), round(r[1], 4), round(r[2], 4), round(r[3], 4)] for r in rows2], "slope": round(slope2, 3)},
        "3D": {"L": L3, "rows": [[int(r[0]), round(r[1], 4), round(r[2], 4), round(r[3], 4)] for r in rows3], "slope": round(slope3, 3)},
        "area_law_2D": bool(area_law_2d),
        "area_law_3D": bool(area_law_3d),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2))
    out = ROOT / "experiments" / "exp_entanglement_area_law_dim_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    print()
    print("Honest boundary:")
    print("  Preset lattice (2D grid / 3D cube) + Fermi-sea ground state; NOT spontaneous 3D from D.")
    print("  Shows the area law and that the exponent gives dimension; the lattice is preset.")

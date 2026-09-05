"""Entanglement entropy: area law (local ground state) vs volume law (random state).

The next spine after the Born weak proposition: "from self-reference D define entanglement,
and show entanglement gives geometry in the low-energy limit."  The key quantity is
    S(rho_A) = -Tr(rho_A ln rho_A)
where rho_A is the reduced density matrix of a subregion A.

The AREA law (S ~ boundary) is the signature of a LOCAL ground state; a generic/random
state gives the VOLUME law (S ~ |A|).  This script demonstrates the distinction on a 1D
chain using free-fermion correlation matrices (Peschel's method):

  Part A (local): ground state of the hopping Hamiltonian (Fermi sea) -> S ~ (1/3) log L.
  Part B (random): random rank-(N/2) projector (random Slater determinant) -> S ~ L.
  Part C (contrast): local grows logarithmically (sub-volume), random grows linearly (volume).

Honest boundary:
  - This is a KINEMATIC entanglement computation: the bipartition uses the chain's
    graph structure (a preset "space"), and the "state" is an ansatz (Fermi sea / random
    projector), not derived from D's dynamics (which is still missing).
  - It demonstrates the area-law/volume-law MECHANISM, not that self-reference D
    spontaneously produces an area law -- that needs a local Hamiltonian + its ground
    state, i.e. the (currently missing) dynamics.
  - In 1D the TRUE area law is S -> const (gapped system); the Fermi sea is gapless, so
    it gives the logarithmic "area-law violation" S ~ (1/3) log L.  Both are sub-volume
    (the signature of locality), in contrast to the volume law S ~ L.
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
    """Peschel formula: S = -sum [lam ln lam + (1-lam) ln(1-lam)] over eigenvalues of C_A."""
    lam = np.linalg.eigvalsh(C_A)
    lam = np.clip(lam, 1e-15, 1.0 - 1e-15)
    return float(-np.sum(lam * np.log(lam) + (1.0 - lam) * np.log(1.0 - lam)))


def fermi_sea_corr(N: int) -> np.ndarray:
    """Correlation matrix of the Fermi sea: ground state of 1D hopping (half filling)."""
    H = np.zeros((N, N))
    for i in range(N - 1):
        H[i, i + 1] = H[i + 1, i] = -1.0
    _, eigvecs = np.linalg.eigh(H)
    occ = eigvecs[:, : N // 2]  # occupy the N/2 lowest single-particle states
    return occ @ occ.T


def random_projector_corr(N: int, seed: int = 0) -> np.ndarray:
    """Correlation matrix of a random rank-(N/2) projector (random Slater determinant)."""
    rng = np.random.default_rng(seed)
    X = rng.standard_normal((N, N // 2))
    Q, _ = np.linalg.qr(X)
    return Q @ Q.T


if __name__ == "__main__":
    N = 400
    Ls = [10, 20, 40, 80, 160, 200]  # up to N/2 (symmetry S(L)=S(N-L))

    C_local = fermi_sea_corr(N)
    C_random = random_projector_corr(N)

    print(f"=== Entanglement entropy: local ground state vs random state (1D chain, N={N}) ===")
    print()
    print("  L      S_local (Fermi sea)      S_random (volume law)")
    rows = []
    for L in Ls:
        s_local = entropy_from_corr(C_local[:L, :L])
        s_random = entropy_from_corr(C_random[:L, :L])
        rows.append((L, s_local, s_random))
        print(f"  {L:4d}      {s_local:8.4f}              {s_random:8.4f}")
    print()

    L_arr = np.array([r[0] for r in rows], dtype=float)
    s_local_arr = np.array([r[1] for r in rows])
    s_random_arr = np.array([r[2] for r in rows])

    logL = np.log(L_arr)
    slope_local_logL = float(np.polyfit(logL, s_local_arr, 1)[0])
    slope_random_L = float(np.polyfit(L_arr, s_random_arr, 1)[0])

    print("growth check:")
    print(f"  S_local  ~ {slope_local_logL:.3f} * log L   (sub-volume: locality signature)")
    print(f"  S_random ~ {slope_random_L:.3f} * L        (volume law: non-local)")
    print()

    # contrast: local grows much slower than random (log vs linear)
    ratio = slope_local_logL / (s_local_arr[-1] / L_arr[-1])
    sub_volume = slope_local_logL < 0.5 * slope_random_L
    print(f"contrast: local = log-growth, random = linear-growth.  local is sub-volume = {sub_volume}")
    print()

    summary = {
        "N": N,
        "Ls": Ls,
        "S_local": [round(x, 4) for x in s_local_arr],
        "S_random": [round(x, 4) for x in s_random_arr],
        "slope_local_vs_logL": round(slope_local_logL, 4),
        "slope_random_vs_L": round(slope_random_L, 4),
        "local_is_sub_volume": bool(sub_volume),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2))
    out = ROOT / "experiments" / "exp_entanglement_area_law_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    print()
    print("Honest boundary:")
    print("  Kinematic entanglement (graph bipartition + ansatz state), NOT background-independent;")
    print("  demonstrates area-law/volume-law mechanism, NOT that self-reference D yields an area law.")

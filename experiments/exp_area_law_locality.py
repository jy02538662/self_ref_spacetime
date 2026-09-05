"""Area law as a locality signature: lattice (Fermi surface) vs random graph.

Key idea (corrected): the entanglement-entropy scaling (area vs volume) is a signature
of whether the ground state has a FERMI SURFACE, i.e. whether the graph is
translation-invariant / "local".

    - 1D chain (translation-invariant, plane-wave eigenstates, Fermi surface)
        -> Fermi sea -> area law  (S ~ log L)
    - Erdos-Renyi random graph (breaks translation invariance, delocalized random
        eigenstates, NO Fermi surface)
        -> "random Fermi sea" -> volume law  (S ~ L)

NOTE (honest correction to an earlier attempt): long-range hopping 1/r^alpha is STILL
translation-invariant, so it still has a Fermi surface and gives the area law regardless
of alpha.  "Locality" in the entanglement sense is NOT the hopping range; it is whether
the ground state has a Fermi surface (translation-invariant) or not (random graph).

Honest boundary:
  - Preset graph (chain / random graph) + Fermi-sea ground state.
  - Shows area law = Fermi-surface/locality signature, volume law = no-Fermi-surface.
  - Does NOT show self-reference D spontaneously becomes local (needs dynamics).
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


def chain_hopping(N: int) -> np.ndarray:
    """1D nearest-neighbour hopping (translation-invariant, has a Fermi surface)."""
    H = np.zeros((N, N))
    for i in range(N - 1):
        H[i, i + 1] = H[i + 1, i] = -1.0
    return H


def random_graph_hopping(N: int, p: float, seed: int = 0) -> np.ndarray:
    """Erdos-Renyi random graph G(N,p) hopping (breaks translation invariance)."""
    rng = np.random.default_rng(seed)
    H = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            if rng.random() < p:
                H[i, j] = H[j, i] = -1.0
    return H


def fermi_sea_corr(H: np.ndarray, fill: float = 0.5) -> np.ndarray:
    """Correlation matrix of the Fermi sea: occupy the lowest fill*N eigenstates."""
    N = H.shape[0]
    _, eigvecs = np.linalg.eigh(H)
    n_occ = int(round(fill * N))
    occ = eigvecs[:, :n_occ]
    return occ @ occ.T


def classify(S_arr: np.ndarray, L_arr: np.ndarray) -> str:
    """Area (log) vs volume (linear): compare growth ratio S(L_max)/S(L_min).

    log growth: S ~ (1/6) ln L -> over 10x L, ratio ~ 1.4;  linear: over 10x L, ratio ~ 10.
    Threshold 2.5 cleanly separates the two (avoids the saturation near L = N/2).
    """
    ratio = float(S_arr[-1] / S_arr[0])
    return "area(log)" if ratio < 2.5 else "volume(linear)"


if __name__ == "__main__":
    N = 200
    Ls = [10, 20, 40, 80, 100]
    L_arr = np.array(Ls, dtype=float)

    print("=== Area law as locality signature: chain (Fermi surface) vs random graph ===")
    print(f"N={N}, Fermi sea (half filling), subregion = first L sites")
    print()

    results = {}

    # 1D chain: area law
    C_chain = fermi_sea_corr(chain_hopping(N))
    S_chain = [entropy_from_corr(C_chain[:L, :L]) for L in Ls]
    kind_chain = classify(np.array(S_chain), L_arr)
    results["chain"] = {"S": {str(L): round(float(s), 4) for L, s in zip(Ls, S_chain)}, "law": kind_chain}
    print(f"  1D chain (Fermi surface):     S(L=10)={S_chain[0]:.3f}  S(L=100)={S_chain[-1]:.3f}  -> {kind_chain}")

    # random graph: volume law
    for p in [0.3, 0.5]:
        C_rg = fermi_sea_corr(random_graph_hopping(N, p, seed=1))
        S_rg = [entropy_from_corr(C_rg[:L, :L]) for L in Ls]
        kind_rg = classify(np.array(S_rg), L_arr)
        results[f"random_graph_p{p}"] = {"S": {str(L): round(float(s), 4) for L, s in zip(Ls, S_rg)}, "law": kind_rg}
        print(f"  random graph G(N,{p}) (no FS): S(L=10)={S_rg[0]:.3f}  S(L=100)={S_rg[-1]:.3f}  -> {kind_rg}")
    print()

    verdict = kind_chain == "area(log)" and all(
        results[k]["law"] == "volume(linear)" for k in results if k.startswith("random")
    )
    print(f"verdict: chain -> {kind_chain}; random graph -> volume(linear).  area law = Fermi-surface/locality signature = {verdict}")
    print()

    summary = {
        "N": N,
        "Ls": Ls,
        "results": results,
        "area_law_is_locality_signature": bool(verdict),
        "note": "long-range hopping 1/r^alpha is still translation-invariant (Fermi surface) -> area law; non-locality (volume law) needs a random graph (no Fermi surface).",
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2))
    out = ROOT / "experiments" / "exp_area_law_locality_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    print()
    print("Honest boundary:")
    print("  Preset graph + Fermi-sea ground state; shows area law = Fermi-surface signature,")
    print("  NOT that self-reference D spontaneously becomes local (needs dynamics).")

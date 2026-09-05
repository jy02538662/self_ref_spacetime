"""String-net condensation sanity check (toy) -- the two "if" conditions, tested in isolation.

Background (theory route, see vault notes):
    self-referential closed strings -> SU(2)_k string-net phase -> spectral triple -> GR,
    and separately self-reference -> Born rule -> entanglement -> spacetime.
    The step "closed strings -> string-net phase" is an if-then theorem with TWO if's:

        if-1  string tension T -> 0  (strings proliferate and condense)   [checkable numerically]
        if-2  the fusion structure is SU(2)_k (gives the topological data) [= "why SU(2)", the postulate]

    This script sanity-checks the two if's SEPARATELY, and is explicit that it does
    NOT assemble them into a proof of "self-referential D -> string-net phase".

Part A (analytic): SU(2)_k string-net topological data.
    quantum dimensions  d_j = [2j+1]_q = sin((2j+1) pi/(k+2)) / sin(pi/(k+2)),
    total quantum dimension  D = sum_j d_j^2,
    topological entanglement entropy  gamma = log D.
    This is the TEE TARGET a SU(2)_k string-net phase must produce (if-2's fingerprint).

Part B (toy loop condensation): a closed-loop gas on a small torus.
    Hilbert space = loop configs (even edge degree at every vertex).
    H = T * (total length) - lam * sum_p (plaquette flip).
    Sweep T from large to 0 and show the ground state goes vacuum -> loop condensate.
    This demonstrates if-1's MECHANISM (tension -> 0 proliferates loops), but in the
    EQUAL-weight (single loop type, D = 2, toric-code-like) limit, NOT the SU(2)_k
    weighted condensate.

Honest boundary:
    - Part A gives the TEE target for SU(2)_k; Part B gives the condensation mechanism
      for the trivial (Z_2-like) category.  Neither proves "closed strings -> SU(2)_k
      string-net phase".  That still needs (i) the actual self-referential Hamiltonian's
      low-energy behaviour and (ii) the full Levin-Wen machinery (F-symbols, B_p).
"""

from __future__ import annotations

import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# -----------------------------------------------------------------------------
# Part A: SU(2)_k topological data (analytic)
# -----------------------------------------------------------------------------

def su2k_quantum_dims(k: int):
    """SU(2)_k quantum dimensions d_j = [2j+1]_q with q = exp(i pi/(k+2)).

    Returns (spins, d) where spins = [0, 1/2, ..., k/2] and d_j = [2j+1]_q.
    """
    spins = np.arange(0, k / 2 + 1e-9, 0.5)
    n = 2.0 * spins + 1.0
    d = np.sin(n * np.pi / (k + 2.0)) / np.sin(np.pi / (k + 2.0))
    return spins, d


def su2k_total_dim(k: int) -> float:
    """Total quantum dimension D = sum_j d_j^2, and TEE gamma = log D."""
    _, d = su2k_quantum_dims(k)
    return float(np.sum(d**2))


# -----------------------------------------------------------------------------
# Part B: toy loop gas on a torus
# -----------------------------------------------------------------------------

def build_lattice(Lx: int, Ly: int):
    """Lx x Ly torus. Returns (E, V, edge_verts, plaquettes).

    edge_verts[e] = (v1, v2) endpoints of edge e.
    plaquettes[p] = tuple of 4 edge indices around plaquette p.
    Edge index: horizontal e = i*Ly+j, vertical e = Lx*Ly + i*Ly+j.
    """
    E = 2 * Lx * Ly
    V = Lx * Ly
    edge_verts = []
    for i in range(Lx):
        for j in range(Ly):
            v = i * Ly + j
            vright = i * Ly + (j + 1) % Ly
            edge_verts.append((v, vright))
    for i in range(Lx):
        for j in range(Ly):
            v = i * Ly + j
            vdown = ((i + 1) % Lx) * Ly + j
            edge_verts.append((v, vdown))
    plaquettes = []
    for i in range(Lx):
        for j in range(Ly):
            h_bot = i * Ly + j
            v_left = Lx * Ly + i * Ly + j
            h_top = ((i + 1) % Lx) * Ly + j
            v_right = Lx * Ly + i * Ly + (j + 1) % Ly
            plaquettes.append((h_bot, v_left, h_top, v_right))
    return E, V, edge_verts, plaquettes


def enumerate_loops(E: int, V: int, edge_verts):
    """All closed-loop configs: binary vectors over E edges with even degree at every vertex."""
    loops = []
    for bits in itertools.product((0, 1), repeat=E):
        degree = [0] * V
        for e, occ in enumerate(bits):
            if occ:
                v1, v2 = edge_verts[e]
                degree[v1] += 1
                degree[v2] += 1
        if all(d % 2 == 0 for d in degree):
            loops.append(bits)
    return loops


def build_hamiltonian(loops, E: int, plaquettes, T: float, lam: float):
    """H = T * length - lam * sum_p (plaquette flip), in the loop-config basis (dense)."""
    n = len(loops)
    idx = {c: k for k, c in enumerate(loops)}
    H = np.zeros((n, n))
    for k, c in enumerate(loops):
        H[k, k] = T * float(sum(c))
    for p in plaquettes:
        for k, c in enumerate(loops):
            c2 = list(c)
            for e in p:
                c2[e] ^= 1
            c2 = tuple(c2)
            if c2 in idx:
                H[k, idx[c2]] -= lam
    return H


def sweep_tension(loops, E: int, plaquettes, lam: float, Ts):
    """Ground state observables vs tension T: loop density and participation ratio."""
    occs = np.array([float(sum(c)) for c in loops])
    n = len(loops)
    out = []
    for T in Ts:
        H = build_hamiltonian(loops, E, plaquettes, T, lam)
        w, v = np.linalg.eigh(H)
        gs = v[:, 0]
        p = gs**2
        density = float(p @ occs) / E
        # participation ratio: ~1 (localized) -> ~n (delocalized / condensed)
        pr = float(1.0 / np.sum(p**2))
        out.append((T, w[0], density, pr))
    return out


# -----------------------------------------------------------------------------
# Part C: weighted string-net condensate (B_p operator, toy, single plaquette)
# -----------------------------------------------------------------------------

def b_p_operator(k: int, weighted: bool = True):
    """Single-plaquette string-net kinetic operator B_p (one-insertion toy).

    States: |vacuum> (index 0) + |loop_j> for j in {1/2, ..., k/2}, N = k/2 + 1.
    Weighted (correct string-net): B_p = (1/D) |v><v|,  |v> = sum_j d_j |loop_j>
        (d_0 = 1, D = sum_j d_j^2).  Ground state = sum_j d_j |loop_j> / sqrt(D).
    Equal (naive):                B_p = (1/N) |w><w|,  |w> = sum_j |loop_j>.
        Ground state = sum_j |loop_j> / sqrt(N).
    The difference is exactly the quantum-dimension distribution d_j^2 / D,
    which is the fingerprint of the string-net condensate (and the source of log D).
    """
    spins, d = su2k_quantum_dims(k)
    N = len(spins)
    if weighted:
        D = float(np.sum(d**2))
        Bp = np.outer(d, d) / D
    else:
        Bp = np.outer(np.ones(N), np.ones(N)) / N
    return spins, Bp


def condensate_distribution(k: int):
    """Loop-type distribution of the single-plaquette condensate (weighted vs equal).

    Returns (spins, p_weighted, p_equal, D, N) where p = |gs|^2 (ground-state weights).
    """
    spins, Bp_w = b_p_operator(k, weighted=True)
    _, Bp_e = b_p_operator(k, weighted=False)
    _, v_w = np.linalg.eigh(Bp_w)
    _, v_e = np.linalg.eigh(Bp_e)
    p_w = v_w[:, -1] ** 2  # largest-eigenvalue eigenvector
    p_e = v_e[:, -1] ** 2
    _, d = su2k_quantum_dims(k)
    D = float(np.sum(d**2))
    N = len(spins)
    return spins, p_w, p_e, D, N


if __name__ == "__main__":
    print("=== String-net condensation sanity check (toy) ===")
    print()

    # ---- Part A: SU(2)_k topological data / TEE target ----
    print("Part A. SU(2)_k string-net topological data (TEE target = log D)")
    print("  (D = sum_j d_j^2,  d_j = [2j+1]_q,  q = exp(i pi/(k+2)))")
    part_a = {}
    for k in (1, 2, 3, 4, 5):
        spins, d = su2k_quantum_dims(k)
        D = float(np.sum(d**2))
        gamma = float(np.log(D))
        part_a[k] = {"D": D, "TEE_logD": gamma, "spins": list(spins), "d": list(np.round(d, 4))}
        ds = ", ".join(f"{x:.3f}" for x in d)
        print(f"  k={k}:  spins=[{', '.join(f'{s:g}' for s in spins)}]  "
              f"d=[{ds}]  D={D:.4f}  TEE=logD={gamma:.4f}")
    print()

    # ---- Part B: toy loop condensation (tension -> 0) ----
    Lx = Ly = 3
    E, V, edge_verts, plaquettes = build_lattice(Lx, Ly)
    loops = enumerate_loops(E, V, edge_verts)
    print(f"Part B. loop gas on {Lx}x{Ly} torus: E={E} edges, {len(loops)} closed-loop configs")
    lam = 1.0
    Ts = [8.0, 4.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.0]
    print(f"  H = T*(length) - lam*sum(plaquette flip),  lam={lam}")
    print("  sweep T large -> 0, ground-state loop density + participation ratio:")
    print(f"  {'T':>5}  {'E0':>9}  {'density':>8}  {'participation':>13}")
    part_b = []
    for T, E0, density, pr in sweep_tension(loops, E, plaquettes, lam, Ts):
        part_b.append({"T": T, "E0": E0, "density": density, "participation": pr})
        print(f"  {T:5.2f}  {E0:9.4f}  {density:8.4f}  {pr:13.2f}")
    print()

    # condensation verdict: does density rise monotonically as T -> 0?
    dens = [r["density"] for r in part_b]
    condensed = dens[0] < dens[-1] and dens[-1] > 0.5 * max(dens)
    print(f"  condensation (density large-T < small-T, and small-T substantial): {condensed}")
    print()

    # ---- Part C: weighted vs equal condensate (single plaquette, B_p operator) ----
    print("Part C. single-plaquette string-net condensate: quantum-dimension weighted vs equal")
    print("  (B_p ground state -> loop-type distribution p_j;  vacuum = j=0)")
    part_c = {}
    for k in (2, 3, 4):
        spins, p_w, p_e, D, N = condensate_distribution(k)
        part_c[k] = {
            "p_weighted": list(np.round(p_w, 4)),
            "p_equal": list(np.round(p_e, 4)),
            "D": D,
            "N": N,
            "vacuum_weighted": float(p_w[0]),
            "vacuum_equal": float(p_e[0]),
        }
        print(f"  k={k}: spins=[{', '.join(f'{s:g}' for s in spins)}]  D={D:.4f}  N={N}")
        print(f"    weighted p_j = [{', '.join(f'{x:.3f}' for x in p_w)}]  vacuum p0={p_w[0]:.3f}")
        print(f"    equal    p_j = [{', '.join(f'{x:.3f}' for x in p_e)}]  vacuum p0={p_e[0]:.3f}")
    print()
    print("  => weighted condensate suppresses vacuum (p0 = 1/D) and favours high-d_j spins;")
    print("     equal condensate is uniform (p0 = 1/N).  The d_j^2/D distribution is the")
    print("     string-net fingerprint, and D (not N) is the correct normalization (log D = TEE).")
    print()

    summary = {
        "criterion": "if-1 (tension->0 condenses loops): density 0 -> substantial; if-2 (SU(2)_k TEE): log D; part-C (weighted B_p -> d_j^2/D distribution)",
        "part_A_SU2k_TEE": {str(k): round(part_a[k]["TEE_logD"], 4) for k in part_a},
        "part_B_loop_gas": {
            "lattice": [Lx, Ly],
            "n_loops": len(loops),
            "density_vs_T": [(r["T"], round(r["density"], 4)) for r in part_b],
            "condensed": bool(condensed),
        },
        "part_C_weighted_condensate": {
            str(k): {
                "p_weighted": part_c[k]["p_weighted"],
                "p_equal": part_c[k]["p_equal"],
                "vacuum_weighted": part_c[k]["vacuum_weighted"],
                "vacuum_equal": part_c[k]["vacuum_equal"],
            }
            for k in part_c
        },
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2))
    out = ROOT / "experiments" / "exp_string_net_condensation_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    print()
    print("Honest boundary:")
    print("  Part A = TEE target for SU(2)_k (if-2 fingerprint, analytic).")
    print("  Part B = condensation MECHANISM for the trivial (equal-weight, D=2) category (if-1).")
    print("  Part C = the string-net kinetic term B_p weights loops by d_j (not equal),")
    print("           giving p_j = d_j^2/D -- the first concrete step toward Levin-Wen,")
    print("           but still single-plaquette / one-insertion (no F-symbols, no fusion at vertices).")
    print("  This does NOT prove 'self-referential closed strings -> SU(2)_k string-net phase'.")

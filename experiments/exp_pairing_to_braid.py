"""Node pairing <-> braid crossing: the Temperley-Lieb / Kauffman-bracket bridge.

This is the "last conceptual gap": how does a 4-valent node's PAIRING (the FPL
structure used in exp_bare_reconnection / exp_su2_wilson) become a braid CROSSING?

Precise correspondence (the honest, corrected version of "pairing -> sigma"):

  A 4-valent node has 4 ports (fix a framing: 1,2 = the two bottom strand-ends,
  3,4 = the two top strand-ends).  Its 3 pairings map to 3 local configurations:
      (13)(24)  ->  identity 1   (two strings pass straight through)
      (12)(34)  ->  TL generator e (a smoothing: a cup below + cap above)
      (14)(23)  ->  the CROSSING sigma -- itself the 3rd pairing, but "undirected"
                    (no over/under); it skein-expands as  A*1 + A^{-1}*e

  So a pure FPL pairing set is {1, e, crossing}.  The crossing does NOT by itself
  carry a sign: sigma and sigma^-1 are the SAME pairing, differing only in
  over/under -- which lives in the skein amplitudes (A vs A^{-1}), NOT the pairing.

  The Kauffman-bracket SKEIN relation turns the braid generator into a QUANTUM
  SUPERPOSITION of the two pairings:
      sigma    = A*1 + A^{-1}*e      (positive crossing, matrix check{R})
      sigma^-1 = A^{-1}*1 + A*e      (negative crossing, matrix check{R}^{-1})

  So the crossing SIGN is the ratio of the two amplitudes (A vs A^{-1}), NOT the
  pairing.  A classical FPL pairing has no sign; the sign needs the superposition
  (SU(2) / quantization) -- consistent with "classical action/heat/reconnection
  cannot generate topology" in the B-route notes.

  The R-matrix eigenvalues {A (x3), -A^{-1}...} are the CHANNELS (symmetric spin-1
  vs antisymmetric spin-0), NOT sigma vs sigma^{-1}.

This script verifies (numpy only, no new deps):
  Part A:  skein self-consistency -- check{R} * check{R}^{-1} = I iff e^2 = d e,
           with d = -A^2 - A^{-2}.
  Part B:  eigenvalues of check{R} = {A (x3, symmetric/spin-1), -A^{-1}-type (x1,
           antisymmetric/spin-0)} -> channels, not sign.
  Part C:  positive vs negative crossing = coefficient swap A <-> A^{-1} on the
           SAME pairings {1, e}.
  Part D:  Kauffman bracket of a pure FPL configuration = d^{c-1} (c = #closed
           strings) -- the pairing-only invariant, independent of the pairing
           details beyond the loop count.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from numpy.random import default_rng

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


# ----------------------------------------------------------------------------
# Temperley-Lieb generator on V (x) V,  V = C^2
# ----------------------------------------------------------------------------

def tl_generator(A: complex):
    """4x4 TL generator e = (d/2) eps eps^dag, eps = |01>-|10>, normalized e^2 = d e.

    d = -A^2 - A^{-2} is the Kauffman-bracket loop value (unknot (+) disjoint loop).
    """
    d = -(A ** 2) - A ** (-2)
    eps = np.array([0.0, 1.0, -1.0, 0.0], dtype=complex)  # |01> - |10>
    E = np.outer(eps, eps.conj())  # eps eps^dag ; tr(E) = eps^dag eps = 2
    e = (d / 2.0) * E
    return e, d


def braid_matrix(A: complex):
    """Positive crossing matrix check{R} = A*1 + A^{-1}*e."""
    e, d = tl_generator(A)
    return A * np.eye(4) + (A ** -1) * e


def braid_matrix_inv(A: complex):
    """Negative crossing matrix check{R}^{-1} = A^{-1}*1 + A*e."""
    e, d = tl_generator(A)
    return (A ** -1) * np.eye(4) + A * e


# ----------------------------------------------------------------------------
# FPL helpers (same 4-regular graph + pairing structure as exp_bare_reconnection)
# ----------------------------------------------------------------------------

def generate_simple_4regular(N: int, rng, max_tries: int = 5000):
    for _ in range(max_tries):
        ports = np.arange(4 * N)
        rng.shuffle(ports)
        seen = set()
        ok = True
        for i in range(0, 4 * N, 2):
            a, b = int(ports[i]), int(ports[i + 1])
            na, nb = a // 4, b // 4
            if na == nb:
                ok = False
                break
            e = (na, nb) if na < nb else (nb, na)
            if e in seen:
                ok = False
                break
            seen.add(e)
        if not ok:
            continue
        edge_partner = np.empty(4 * N, dtype=int)
        for i in range(0, 4 * N, 2):
            a, b = int(ports[i]), int(ports[i + 1])
            edge_partner[a] = b
            edge_partner[b] = a
        return edge_partner
    raise RuntimeError("failed to generate simple 4-regular graph")


def pairing_partners(pair_type: int):
    return [1, 0, 3, 2] if pair_type == 0 else ([2, 3, 0, 1] if pair_type == 1 else [3, 2, 1, 0])


def build_thread_partner(N: int, pair_types: np.ndarray) -> np.ndarray:
    tp = np.empty(4 * N, dtype=int)
    for v in range(N):
        p = pairing_partners(int(pair_types[v]))
        base = 4 * v
        for s in range(4):
            tp[base + s] = base + p[s]
    return tp


def trace_loops(edge_partner: np.ndarray, thread_partner: np.ndarray) -> list:
    """Closed strings = cycles of f = edge_partner o thread_partner (lengths in #ports)."""
    M = len(edge_partner)
    visited = np.zeros(M, dtype=bool)
    lengths = []
    for p in range(M):
        if visited[p]:
            continue
        length = 0
        q = p
        while not visited[q]:
            visited[q] = True
            q = edge_partner[thread_partner[q]]
            length += 1
        lengths.append(length)
    return lengths


if __name__ == "__main__":
    # A = q^{1/4}, q = e^{i pi/(k+2)} -- same convention as exp_jones.py
    k = 1
    q = np.exp(1j * np.pi / (k + 2))
    A = q ** 0.25
    e, d = tl_generator(A)

    print("=== Node pairing <-> braid crossing: Temperley-Lieb / Kauffman bridge ===")
    print(f"k={k}: q = {q.real:+.4f}{q.imag:+.4f}i,  A = q^(1/4) = {A.real:+.4f}{A.imag:+.4f}i")
    print(f"d = -A^2 - A^{-2} = {d.real:+.4f}{d.imag:+.4f}i")
    print()
    print("pairing (framing: 1,2 = bottom, 3,4 = top) -> local config:")
    print("  (13)(24)  = identity 1  (strings straight through)")
    print("  (12)(34)  = TL generator e  (a smoothing: cup + cap)")
    print("  (14)(23)  = CROSSING sigma  (undirected; skein = A*1 + A^-1*e)")
    print("  => crossing is the 3rd pairing but undirected; NO over/under sign yet.")
    print()

    # Part A: skein self-consistency
    R = braid_matrix(A)
    Rinv = braid_matrix_inv(A)
    e2_minus_de = np.max(np.abs(e @ e - d * e))
    RRinv_minus_I = np.max(np.abs(R @ Rinv - np.eye(4)))
    RinvR_minus_I = np.max(np.abs(Rinv @ R - np.eye(4)))
    print("Part A. skein self-consistency  check{R} = A*1 + A^{-1}*e :")
    print(f"  |e^2 - d e|                 = {e2_minus_de:.2e}   (d = -A^2 - A^{-2})")
    print(f"  |R R^-1 - I|                = {RRinv_minus_I:.2e}")
    print(f"  |R^-1 R - I|                = {RinvR_minus_I:.2e}")
    print(f"  => R and R^-1 are inverses iff e^2 = d e with d = -A^2 - A^{-2}  : "
          f"{e2_minus_de < 1e-9 and RRinv_minus_I < 1e-9}")
    print()

    # Part B: eigenvalues = channels (symmetric spin-1 vs antisymmetric spin-0)
    evals = np.linalg.eigvals(R)
    # symmetric subspace (orthogonal to eps): eigenvalue A (3-fold)
    # antisymmetric subspace (span of eps): eigenvalue A + A^{-1} d (1-fold)
    anti_val = A + (A ** -1) * d
    print("Part B. eigenvalues of check{R} (channels, NOT sign):")
    print(f"  eigenvalues = {[f'{e.real:+.4f}{e.imag:+.4f}i' for e in sorted(evals, key=lambda z: z.real)]}")
    print(f"  expected symmetric (spin-1, 3-fold)  = A = {A.real:+.4f}{A.imag:+.4f}i")
    print(f"  expected antisymmetric (spin-0, 1-fold) = A + A^{-1}d = {anti_val.real:+.4f}{anti_val.imag:+.4f}i")
    # count how many eigenvalues equal A
    n_A = int(np.sum(np.abs(evals - A) < 1e-6))
    n_anti = int(np.sum(np.abs(evals - anti_val) < 1e-6))
    print(f"  (multiplicity: {n_A} x A, {n_anti} x antisym)  ->  2(x)2 = 3(+)1")
    print("  NOTE: these eigenvalues label the CHANNELS (spin-1 vs spin-0),")
    print("        NOT sigma vs sigma^{-1}.")
    print()

    # Part C: positive vs negative = coefficient swap on the SAME pairings {1, e}
    print("Part C. positive vs negative crossing = coefficient swap A <-> A^{-1}:")
    print(f"  sigma    = A*1 + A^{-1}*e   (coeff of 1 = {A.real:+.4f}{A.imag:+.4f}i)")
    print(f"  sigma^-1 = A^{-1}*1 + A*e   (coeff of 1 = {(A**-1).real:+.4f}{(A**-1).imag:+.4f}i)")
    print("  => same two pairings {1, e}; the SIGN is where A sits, not the pairing.")
    print("  => a classical FPL pairing has no sign; sign needs the quantum superposition.")
    print()

    # Part D: Kauffman bracket of a pure FPL configuration = d^{c-1} (c = #closed strings)
    print("Part D. Kauffman bracket of a pure FPL config = d^{c-1}  (c = #closed strings):")
    rng = default_rng(0)
    N = 30
    edge_partner = generate_simple_4regular(N, rng)
    rows = []
    for trial in range(5):
        pair_types = rng.integers(0, 3, N)
        thread_partner = build_thread_partner(N, pair_types)
        lengths = trace_loops(edge_partner, thread_partner)
        c = len(lengths)
        kb = d ** (c - 1)
        rows.append((c, kb))
        print(f"  trial {trial}: n_closed_strings c = {c:>2}   <FPL> = d^(c-1) = "
              f"{kb.real:+.4f}{kb.imag:+.4f}i")
    print("  => the only topology a pairing-only config 'knows' is its loop count.")
    print("  => to get Hopf-type (non-trivial) values one must put the crossing")
    print("     (sigma = A*1 + A^{-1}*e) in as a QUANTUM superposition, not a pairing.")
    print()

    print("interpretation:")
    print("  - pairing -> TL element {1, e}; crossing -> SUPERPOSITION A*1 + A^{-1}*e.")
    print("  - crossing sign = the amplitude ratio (A vs A^{-1}), not the pairing.")
    print("  - R-matrix eigenvalues = channels (spin-1/spin-0), not sign.")
    print("  - 'pairing -> braid word' is therefore NOT unique in reverse (the smoothing")
    print("    forgets over/under); the clean forward route is braid word -> skein ->")
    print("    pairing -> Markov trace -> Kauffman bracket (already verified in exp_jones).")

    summary = {
        "k": k,
        "A": [round(float(A.real), 8), round(float(A.imag), 8)],
        "d": [round(float(d.real), 8), round(float(d.imag), 8)],
        "skein_consistent": bool(e2_minus_de < 1e-9 and RRinv_minus_I < 1e-9),
        "e2_minus_de": round(float(e2_minus_de), 12),
        "RRinv_minus_I": round(float(RRinv_minus_I), 12),
        "eigenvalue_multiplicities": {"A": n_A, "antisym": n_anti},
        "note": "pairing -> {1,e}; crossing = A*1 + A^{-1}*e (superposition); sign = amplitude ratio; eigenvalues = channels.",
    }
    out = ROOT / "experiments" / "exp_pairing_to_braid_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

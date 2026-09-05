"""SU(2) Wilson loops on a 4-regular graph: relative linking (skein relation).

Upgrade the edge "phase" from U(1) (e^{i theta}) to SU(2) (2x2 matrices).  The Wilson
loop is now  W(C) = tr( prod_{edges in C} U_ij ), and the ORDER of multiplication matters.

Key physics (correct): U(1) is abelian -> only the overall phase (the seed Q).  SU(2) is
non-abelian -> the order matters, so two loops can "braid", giving the RELATIVE linking
(Hopf / Jones).  This is "why SU(2)": only non-abelian structure gives relative linking.

This script:
  Part A: the individual SU(2) Wilson-loop traces are NON-trivial (in (-2,2)) and CHANGE
          under node reconnection, in contrast to U(1) where Q = prod W(C) was a strict
          invariant (single number).
  Part B: verifies the SU(2) Fierz / skein identity  tr(A)tr(B) = tr(AB) + tr(A B^dag),
          which is the "relative linking" structure (the cross term A B^dag is the braid).

Honest note: the SU(2) invariant under surgery is the Jones polynomial (a polynomial in a
variable), NOT a single number like the U(1) Q.  So this gives the skein relation, not a
clean integer.
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


def set_node_pairing(tp: np.ndarray, v: int, pair_type: int):
    p = pairing_partners(pair_type)
    base = 4 * v
    for s in range(4):
        tp[base + s] = base + p[s]


def random_su2(rng):
    """Random SU(2) matrix from a uniform point on S^3."""
    q = rng.standard_normal(4)
    q = q / np.linalg.norm(q)
    a0, ax, ay, az = q
    return np.array(
        [[a0 + 1j * az, ay + 1j * ax],
         [-ay + 1j * ax, a0 - 1j * az]],
        dtype=complex,
    )


def assign_edge_su2(edge_partner, rng):
    """edge_U[t] = SU(2) matrix when crossing edge from port t to edge_partner[t]."""
    M = len(edge_partner)
    edge_U = np.empty((M, 2, 2), dtype=complex)
    for a in range(M):
        b = int(edge_partner[a])
        if a < b:
            U = random_su2(rng)
            edge_U[a] = U
            edge_U[b] = U.conj().T  # U^dag (crossing b -> a)
    return edge_U


def trace_loops_holonomy(edge_partner, thread_partner, edge_U):
    """Return list of SU(2) holonomies (one per closed string)."""
    M = len(edge_partner)
    visited = np.zeros(M, dtype=bool)
    holonomies = []
    for p in range(M):
        if visited[p]:
            continue
        U = np.eye(2, dtype=complex)
        q = p
        while not visited[q]:
            visited[q] = True
            t = int(thread_partner[q])
            if t < int(edge_partner[t]):  # canonical direction: count each edge once
                U = U @ edge_U[t]
            q = int(edge_partner[t])
        holonomies.append(U)
    return holonomies


def traces_of(holonomies):
    return np.array([np.real(np.trace(U)) for U in holonomies])


if __name__ == "__main__":
    rng = default_rng(0)
    N = 30

    edge_partner = generate_simple_4regular(N, rng)
    edge_U = assign_edge_su2(edge_partner, rng)
    pair_types = rng.integers(0, 3, N)
    thread_partner = build_thread_partner(N, pair_types)

    hol = trace_loops_holonomy(edge_partner, thread_partner, edge_U)
    tr0 = traces_of(hol)
    print("=== SU(2) Wilson loops under reconnection (4-regular, N=30) ===")
    print(f"n_loops initial = {len(hol)}")
    print(f"initial traces tr(W(C)) = {[f'{t:.3f}' for t in tr0]}")
    print(f"  (non-trivial: in (-2,2), NOT all = 2 or = 1 like U(1) seed)")
    print()

    # Part A: traces change under reconnection (non-invariant)
    steps = 10_000
    changed = False
    for t in range(1, steps + 1):
        v = int(rng.integers(0, N))
        old = int(pair_types[v])
        new = int(rng.integers(0, 3))
        while new == old:
            new = int(rng.integers(0, 3))
        pair_types[v] = new
        set_node_pairing(thread_partner, v, new)
        if t % 2000 == 0:
            hol = trace_loops_holonomy(edge_partner, thread_partner, edge_U)
            tr = traces_of(hol)
            # multiset comparison: different if loop count differs, else if sorted traces differ
            is_diff = len(tr) != len(tr0) or np.max(np.abs(np.sort(tr) - np.sort(tr0))) > 1e-9
            if is_diff:
                changed = True
            print(f"t={t:>5}: n_loops={len(tr):>2}  traces = {[f'{x:.2f}' for x in np.sort(tr)]}")

    print()
    print(f"SU(2) traces CHANGE under reconnection = {changed}  (vs U(1) Q was strictly invariant)")
    print()

    # Part B: Fierz / skein identity  tr(A)tr(B) = tr(AB) + tr(A B^dag)
    print("Part B: SU(2) Fierz / skein identity  tr(A)tr(B) = tr(AB) + tr(A B^dag)")
    max_err = 0.0
    for _ in range(200):
        A = random_su2(rng)
        B = random_su2(rng)
        lhs = np.real(np.trace(A)) * np.real(np.trace(B))
        rhs = np.real(np.trace(A @ B)) + np.real(np.trace(A @ B.conj().T))
        max_err = max(max_err, abs(lhs - rhs))
    print(f"  max |tr(A)tr(B) - [tr(AB)+tr(AB^dag)]| over 200 random pairs = {max_err:.2e}")
    print(f"  skein identity holds = {max_err < 1e-9}")
    print()

    summary = {
        "N": N,
        "n_loops_initial": len(hol),
        "initial_traces": [round(float(t), 4) for t in tr0],
        "traces_change_under_reconnection": bool(changed),
        "skein_identity_max_err": round(float(max_err), 12),
        "skein_identity_holds": bool(max_err < 1e-9),
        "note": "SU(2) non-abelian -> relative linking (skein/Jones), not the U(1) single-number seed Q.",
    }
    out = ROOT / "experiments" / "exp_su2_wilson_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

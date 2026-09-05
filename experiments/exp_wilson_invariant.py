"""Global invariant under node reconnection: product of all Wilson loops.

Claim (correct): on a 4-regular graph with edge phases U_ij = e^{i theta_ij}, the product
over ALL closed strings of their Wilson loops
    Q = prod_{closed strings C} W(C),   W(C) = prod_{edges in C} U_ij
equals the product over all edges (each edge appears exactly once in the loop
decomposition), so it is INDEPENDENT of the pairing and hence a STRICT invariant under
node reconnection (merge/split of loops preserves it: W(C3)=W(C1)W(C2) on merge).

This is the "seed" of the topological charge: the overall phase conservation; the
relative linking (Hopf) is the next step.  It comes straight from the self-reflexivity
axiom  D_ij = D_ji*  =>  U_ij = e^{i theta_ij}  =>  W(C)  =>  Q.

This script verifies Q is exactly invariant under random reconnection, and that
Q_wilson = Q_direct = exp(i sum_{edges} theta).
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


def assign_edge_phases(edge_partner: np.ndarray, rng) -> np.ndarray:
    """edge_phase[a] = phase when crossing the edge from port a to edge_partner[a]."""
    M = len(edge_partner)
    phase = np.zeros(M)
    for a in range(M):
        b = int(edge_partner[a])
        if a < b:
            theta = rng.uniform(-np.pi, np.pi)
            phase[a] = theta
            phase[b] = -theta
    return phase


def trace_wilson(edge_partner, thread_partner, edge_phase):
    """Return (wilson_loops, total_Q) = (list of W(C), prod W(C)).

    Each edge is crossed twice in the f-cycle (once from each end); count it ONCE,
    in the canonical direction t -> edge_partner[t] (t < edge_partner[t]).
    """
    M = len(edge_partner)
    visited = np.zeros(M, dtype=bool)
    wilson = []
    for p in range(M):
        if visited[p]:
            continue
        ph = 0.0
        q = p
        while not visited[q]:
            visited[q] = True
            t = int(thread_partner[q])   # thread mate (same node)
            if t < int(edge_partner[t]):  # canonical direction: count each edge once
                ph += edge_phase[t]
            q = int(edge_partner[t])     # cross the edge
        wilson.append(np.exp(1j * ph))
    Q = np.prod(np.array(wilson))
    return wilson, Q


if __name__ == "__main__":
    rng = default_rng(0)
    N = 30

    edge_partner = generate_simple_4regular(N, rng)
    edge_phase = assign_edge_phases(edge_partner, rng)

    # direct invariant: Q_direct = exp(i * sum over edges a<b of theta)
    Q_direct = np.exp(1j * sum(edge_phase[a] for a in range(4 * N) if a < int(edge_partner[a])))

    pair_types = rng.integers(0, 3, N)
    thread_partner = build_thread_partner(N, pair_types)

    wilson, Q0 = trace_wilson(edge_partner, thread_partner, edge_phase)
    print("=== Global Wilson-loop invariant under reconnection (4-regular, N=30) ===")
    print(f"n_loops initial = {len(wilson)}")
    print(f"Q_wilson (prod W(C)) = {Q0:.6f}  arg = {np.angle(Q0):+.6f}")
    print(f"Q_direct (prod edges) = {Q_direct:.6f}  arg = {np.angle(Q_direct):+.6f}")
    print(f"match Q_wilson ~= Q_direct: {abs(Q0 - Q_direct) < 1e-6}")
    print()

    # random reconnection, track Q_wilson invariance
    steps = 50_000
    record_every = 5000
    max_dev = 0.0
    for t in range(1, steps + 1):
        v = int(rng.integers(0, N))
        old = int(pair_types[v])
        new = int(rng.integers(0, 3))
        while new == old:
            new = int(rng.integers(0, 3))
        pair_types[v] = new
        set_node_pairing(thread_partner, v, new)

        if t % record_every == 0:
            _, Q = trace_wilson(edge_partner, thread_partner, edge_phase)
            dev = abs(Q - Q_direct)
            max_dev = max(max_dev, dev)
            print(f"t={t:>6}: |Q_wilson - Q_direct| = {dev:.2e}   arg(Q_wilson) = {np.angle(Q):+.6f}")

    print()
    print(f"max |Q_wilson - Q_direct| over reconnection = {max_dev:.2e}")
    print(f"INVARIANT (dev < 1e-9): {max_dev < 1e-9}")
    print()
    print("interpretation:")
    print("  Q = prod W(C) = prod edges e^{i theta} is a strict global invariant under reconnection.")
    print("  It is the SEED of the topological charge (overall phase conservation);")
    print("  the relative linking (Hopf) is the next step, built on top of this seed.")

    summary = {
        "N": N,
        "n_loops_initial": len(wilson),
        "Q_wilson": [round(float(Q0.real), 8), round(float(Q0.imag), 8)],
        "Q_direct": [round(float(Q_direct.real), 8), round(float(Q_direct.imag), 8)],
        "max_deviation": round(float(max_dev), 12),
        "invariant": bool(max_dev < 1e-9),
        "note": "Q = prod W(C) = prod edges e^{i theta} is invariant under reconnection (seed of topological charge).",
    }
    out = ROOT / "experiments" / "exp_wilson_invariant_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

"""Bare reconnection dynamics on a 4-regular graph (toy, "fully packed loop" + surgery).

Model (per design):
  - N nodes, each degree 4 (graph FIXED; reconnection does NOT change edges).
  - At each node, its 4 ports are paired into 2 "threads" (3 possible pairings).
  - A closed string = cycle of the composite permutation  f = edge_partner o thread_partner.
  - One reconnection = pick a random node, switch its pairing to a different one, accept all.

Question: does this bare (no energy / no time / no space) dynamics self-organize?

Honest expectation: since every move is accepted, this is a SYMMETRIC random walk on the
3^N pairing configurations, so it converges to the UNIFORM measure over FPL configs.
Under uniform measure there is no preference for large loops -> outcome B (random, no
self-organization), NOT C.  That is itself a meaningful negative result: bare dynamics
without a preference does not organise (same lesson as "classical action can't generate
topology" and "classical heat can't break Z2").

Note on d_s: the underlying graph is fixed (only pairings change), so the graph
Laplacian / spectral dimension is constant in time; for a random regular graph it is not
even a lattice-like power law.  So the meaningful observables are the loop statistics.
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
    """Random simple 4-regular graph via configuration model + rejection (no self/multi-edge)."""
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
    """Thread pairing of a node's 4 ports: 0->(0,1)(2,3); 1->(0,2)(1,3); 2->(0,3)(1,2)."""
    return [1, 0, 3, 2] if pair_type == 0 else ([2, 3, 0, 1] if pair_type == 1 else [3, 2, 1, 0])


def build_thread_partner(N: int, pair_types: np.ndarray) -> np.ndarray:
    """thread_partner[port] = its thread mate (same node)."""
    tp = np.empty(4 * N, dtype=int)
    for v in range(N):
        p = pairing_partners(int(pair_types[v]))
        base = 4 * v
        for s in range(4):
            tp[base + s] = base + p[s]
    return tp


def set_node_pairing(tp: np.ndarray, v: int, pair_type: int):
    """Update thread_partner for node v to a new pairing."""
    p = pairing_partners(pair_type)
    base = 4 * v
    for s in range(4):
        tp[base + s] = base + p[s]


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


def loop_stats(lengths: list):
    """(n_loops, max_len_in_nodes, avg_len_in_nodes). length in nodes = length/2."""
    n = len(lengths)
    if n == 0:
        return 0, 0, 0.0
    nodes = [L // 2 for L in lengths]
    return n, max(nodes), sum(nodes) / n


if __name__ == "__main__":
    rng = default_rng(0)
    N = 30

    edge_partner = generate_simple_4regular(N, rng)
    pair_types = rng.integers(0, 3, N)
    thread_partner = build_thread_partner(N, pair_types)

    # sanity: total ports = 4N, sum of loop lengths = 4N
    lengths = trace_loops(edge_partner, thread_partner)
    assert sum(lengths) == 4 * N, f"loop partition wrong: {sum(lengths)} != {4*N}"

    def report(t):
        n, mx, avg = loop_stats(lengths)
        return f"t={t:>6}: n_loops={n:>3}  max_len={mx:>3}  avg_len={avg:6.3f}"

    print("=== Bare reconnection dynamics (4-regular graph, N=30) ===")
    print("observables: n_loops, max_len, avg_len (in nodes).  No energy, all moves accepted.")
    print()
    print(report(0))

    steps = 100_000
    record_every = 1000
    history = []
    for t in range(1, steps + 1):
        v = int(rng.integers(0, N))
        old = int(pair_types[v])
        new = int(rng.integers(0, 3))
        while new == old:
            new = int(rng.integers(0, 3))
        pair_types[v] = new
        set_node_pairing(thread_partner, v, new)

        if t % record_every == 0:
            lengths = trace_loops(edge_partner, thread_partner)
            n, mx, avg = loop_stats(lengths)
            history.append([t, n, mx, avg])
            print(report(t))
    print()

    # verdict: did it self-organise (few large loops) or stay random?
    n0, mx0, avg0 = loop_stats(trace_loops(edge_partner, thread_partner))
    # "self-organise" would mean n_loops -> small and max_len -> large, consistently.
    # Compare: uniform FPL on a 4-regular graph has many loops of moderate size.
    print(f"final: n_loops={n0}, max_len={mx0}, avg_len={avg0:.3f}")
    print()
    print("interpretation:")
    print("  all moves accepted => symmetric random walk on 3^N configs => uniform FPL measure.")
    print("  uniform measure has no large-loop preference => outcome B (random, no self-organisation).")
    print()

    summary = {
        "N": N,
        "steps": steps,
        "record_every": record_every,
        "initial": [history[0][1], history[0][2], history[0][3]] if history else [n0, mx0, avg0],
        "final": [n0, mx0, avg0],
        "history": [[int(t), int(n), int(mx), round(avg, 4)] for t, n, mx, avg in history],
        "note": "graph fixed => d_s constant; uniform random walk => no self-organisation (outcome B).",
    }
    out = ROOT / "experiments" / "exp_bare_reconnection_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

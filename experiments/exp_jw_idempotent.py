"""DEPRECATED (2026-09-08): this experiment was judged unnecessary and is NOT run.

The premise was wrong: f_n^2 = f_n is an ALGEBRAIC IDENTITY (holds for ALL delta != 0,
verified by hand for f_2 = 1 - e_1/delta). Unit roots act via TRUNCATION (f_{k+2} vanishes,
Delta_{k+1}=0), NOT via idempotency breaking. So "scan delta, look at ||f_n^2-f_n||" has no
information -- it is zero everywhere. Stage 2a is closed (quantization = finiteness + truncation).

The TL string-diagram basis below is also buggy (used single-row pairings instead of the
two-row TL convention). Left as-is for the record; do not fix, do not run.

--- original docstring (premise now rejected) ---
阶段 2a：Jones-Wenzl 幂等元 f_n 的存在性 = 单位根条件（弦图空间上的量子选择）。

User's insight: "quantumization selects unit roots" lives NOT in plaquette phase integrals
(those are separable -> uniform), but in the TEMPERLEY-LIEB (string diagram) space, where
the idempotent f_n (Jones-Wenzl projector) satisfies f_n^2 = f_n exactly at unit roots.

Experiment: scan delta continuously; for each n, construct f_n explicitly in the TL_n
string-diagram basis; measure the idempotency residual ||f_n^2 - f_n||; see whether it
is EXACTLY zero at delta = 2cos(pi/(N+1)) (Chebyshev zero) and how it breaks away.

Recursive definition: f_1 = 1;  f_n = f_{n-1} - (Delta_{n-2}/Delta_{n-1}) f_{n-1} e_{n-1} f_{n-1},
Delta_0=1, Delta_1=delta, Delta_{n+1}=delta*Delta_n - Delta_{n-1}.
"""

from __future__ import annotations

import numpy as np


def catalan_pairings(n):
    """All non-crossing pairings of 2n points, as list of pairs (a<b)."""
    def rec(pts):
        if not pts:
            return [[]]
        out = []
        a = pts[0]
        for k in range(1, len(pts), 2):  # pair a with pts[k] (k odd keeps segments even)
            b = pts[k]
            left = pts[1:k]
            right = pts[k + 1:]
            for L in rec(left):
                for R in rec(right):
                    out.append([(a, b)] + L + R)
        return out
    return rec(list(range(2 * n)))


def pair_map(n):
    """(pairing, pair list) -> index dictionary."""
    return [dict((x, y) for (x, y) in pairing) | dict((y, x) for (x, y) in pairing)
            for pairing in catalan_pairings(n)]


def e_matrix(n, i, delta):
    """e_i acting on TL_n string-diagram basis (i in 0..n-2, cup-cap at positions i,i+1)."""
    pairings = catalan_pairings(n)
    pmaps = pair_map(n)
    dim = len(pairings)
    M = np.zeros((dim, dim), float)
    for a in range(dim):
        # e_i * diagram_a: add cup (i,i+1) top and cap (n+i, n+i+1) bottom, then resolve
        # top cup connects top point i and i+1; bottom cap connects bottom point n+i, n+i+1
        # build the combined loop system and count closed loops
        # diagram a: pairs among 2n points (top 0..n-1, bottom n..2n-1)
        pa = pmaps[a]
        # new pairing = pa restricted to points not in {i, i+1, n+i, n+i+1}, plus loops
        # cup: (i, i+1); cap: (n+i, n+i+1)
        # follow connections and count closed loops
        # simpler: build union-find over the 2n points with edges: pa's pairs, cup, cap
        parent = list(range(2 * n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x, y):
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        for (x, y) in pairings[a]:
            union(x, y)
        union(i, i + 1)
        union(n + i, n + i + 1)
        # count loops = number of components
        roots = set(find(x) for x in range(2 * n))
        loops = len(roots)
        # the result is delta^{loops} times the "straight" diagram (with cup/cap removed)
        # Here we only need the *algebraic* structure: e_i^2 = delta e_i, etc.
        # For idempotency check we instead build f_n via the recursive FORMULA below.
        M[a, a] = delta ** (loops - n)  # normalized so e_i has correct TL relations
    return M


def f_matrices(n, delta):
    """Recursive Jones-Wenzl projectors f_1..f_n as matrices in TL_n (via nested bases).
    For simplicity we build f_n directly from the recurrence in the TL_n basis using
    explicit e_i matrices and the embedding f_{n-1} -> TL_n (straight last strand)."""
    # Chebyshev
    def Delta(k):
        if k == 0:
            return 1.0
        if k == 1:
            return delta
        d0, d1 = 1.0, delta
        for _ in range(2, k + 1):
            d0, d1 = d1, delta * d1 - d0
        return d1

    pairings = catalan_pairings(n)
    dim = len(pairings)
    E = [e_matrix(n, i, delta) for i in range(n - 1)]

    # embedding of TL_{n-1} into TL_n: f_{n-1} acts on first n-1 strands, last strand straight
    # We build f_n step by step: f_1 = id in TL_1, then embed and apply recurrence.
    # --- TL_1: dim C_1 = 1, f_1 = 1 ---
    f = np.array([[1.0]])
    for m in range(2, n + 1):
        # embed f (dim C_{m-1}) into TL_m (dim C_m): last strand straight
        pairings_m = catalan_pairings(m)
        dim_m = len(pairings_m)
        # map a TL_{m-1} diagram to TL_m diagram by adding a straight pair (n-1, n) ... 
        # (use the "connect top m-1 to bottom m-1" = straight last strand)
        # Build embedding matrix P: (dim_m, dim_prev)
        P = np.zeros((dim_m, f.shape[0]))
        prev_pairs = catalan_pairings(m - 1)
        for b, pp in enumerate(prev_pairs):
            # add straight strand: pair top point (m-1) with bottom point (2m-2)
            new_pairing = [(x, y) for (x, y) in pp] + [(m - 1, 2 * m - 2)]
            # find its index in pairings_m
            # normalize: sort pairs
            key = tuple(sorted(tuple(sorted(p)) for p in new_pairing))
            keys = {tuple(sorted(tuple(sorted(p)) for p in pl)): k for k, pl in enumerate(pairings_m)}
            P[keys[key], b] = 1.0
        f_emb = P @ f @ P.T  # embedded projector
        # e_{m-1} in TL_m
        e = E[m - 2]
        f = f_emb - (Delta(m - 2) / Delta(m - 1)) * (f_emb @ e @ f_emb)
    return f, E, Delta


def main():
    print("Jones-Wenzl idempotency residual ||f_n^2 - f_n|| vs delta")
    print()
    n = 4
    print(f"n={n}: scan delta in [0, 2], check ||f_n^2 - f_n||")
    print(f"  unit root: delta = 2cos(pi/(N+1)). For f_n (n={n}), Chebyshev zero Delta_n=0 at")
    print(f"  delta = 2cos(pi m/(n+1)): m=1 -> {2*np.cos(np.pi/(n+1)):.4f}")
    print()
    for delta in np.linspace(0.0, 2.0, 17):
        f, E, Delta = f_matrices(n, delta)
        resid = np.linalg.norm(f @ f - f)
        # also print Delta_{n-1} (denominator of recurrence)
        print(f"  delta={delta:.3f}  ||f_n^2-f_n||={resid:.3e}  Delta_{n-1}={Delta(n-1):+.3e}  Delta_n={Delta(n):+.3e}")


if __name__ == "__main__":
    main()

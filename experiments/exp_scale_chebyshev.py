"""阶段 2b 步骤 1：Chebyshev 零点间隔 = 尺度（有限 N 生成尺度）。

Insight (v4 9.2): the gap between adjacent unit-root eigenvalues 2cos(pi m/(N+1))
is O(1/N) in the bulk but O(1/N^2) at the edge (near ground state). As N -> infinity
the gap -> 0 (wall-3: no intrinsic scale); for finite N the gap != 0 (scale generation).

Verify: edge gap ~ N^{-2}, bulk gap ~ N^{-1}, both -> 0 as N -> infinity.
"""

from __future__ import annotations

import numpy as np


def chebyshev_zeros(N):
    return np.array([2.0 * np.cos(np.pi * m / (N + 1)) for m in range(1, N + 1)])


def main():
    print("Chebyshev zeros 2cos(pi m/(N+1)): adjacent gaps vs N")
    print(f"{'N':>5} {'edge_gap':>12} {'edge*N^2':>12} {'bulk_gap':>12} {'bulk*N':>12}")
    for N in [5, 10, 20, 40, 80, 160, 320]:
        z = np.sort(chebyshev_zeros(N))  # ascending: -2 ... +2
        # edge gap: between the two largest zeros (near +2, i.e. m=1,2 -> ground state edge)
        edge_gap = z[-1] - z[-2]
        # bulk gap: near zero
        mid = int(np.argmin(np.abs(z)))
        bulk_gap = z[mid + 1] - z[mid] if mid + 1 < len(z) else float('nan')
        print(f"{N:>5} {edge_gap:>12.6f} {edge_gap * N**2:>12.4f} {bulk_gap:>12.6f} {bulk_gap * N:>12.4f}")

    print()
    print("interpretation:")
    print("  - edge_gap * N^2 -> const  => edge gap ~ N^{-2} (scale, survives longest)")
    print("  - bulk_gap * N  -> const  => bulk gap ~ N^{-1}")
    print("  - both -> 0 as N -> infinity  => wall-3 (no intrinsic scale in the continuum limit)")
    print("  - finite N keeps a non-zero edge gap => finite N GENERATES a scale")


if __name__ == "__main__":
    main()

"""Two closed strings on D, trying to link (twist into a helix) -- can it succeed?

Core intuition (spiral = double-strand winding = link = 3D): two closed curves
can only truly link (non-zero linking number) in 3D; in 2D the linking number is
identically 0 (they can only overlap, never thread through each other).

This verifies that fact with the Gauss linking integral, discretized:
    Lk = (1/4pi) sum_i sum_j (r1_i - r2_j) . (dr1_i x dr2_j) / |r1_i - r2_j|^3

Cases:
  - 2D embedding: both curves in the z=0 plane, separated -> Lk = 0.
  - 3D embedding (Hopf link): curve1 = unit circle in xy-plane,
                              curve2 = unit circle in yz-plane threading curve1 -> Lk = 1.
"""

from __future__ import annotations

import numpy as np


def gauss_linking(r1, r2):
    """Discretized Gauss linking integral. r1, r2: (n,3) arrays of curve points."""
    dr1 = np.roll(r1, -1, axis=0) - r1  # segment r1[i] -> r1[i+1 mod n]
    dr2 = np.roll(r2, -1, axis=0) - r2
    total = 0.0
    for a in range(len(r1)):
        for b in range(len(r2)):
            delta = r1[a] - r2[b]
            dist3 = np.dot(delta, delta) ** 1.5
            if dist3 < 1e-18:
                continue
            total += np.dot(delta, np.cross(dr1[a], dr2[b])) / dist3
    return total / (4.0 * np.pi)


def circle_xy(theta, R=1.0):
    return np.column_stack([R * np.cos(theta), R * np.sin(theta), np.zeros_like(theta)])


def circle_xy_shifted(theta, cx=2.0):
    return np.column_stack([cx + np.cos(theta), np.sin(theta), np.zeros_like(theta)])


def circle_threading(theta, R=2.0):
    # circle in the xz-plane, center (R,0,0) (ON circle1), radius 1 -> threads circle1's disk once
    return np.column_stack([R + np.cos(theta), np.zeros_like(theta), np.sin(theta)])


def main():
    n = 300
    theta = np.linspace(0, 2 * np.pi, n, endpoint=False)

    c1 = circle_xy(theta, R=2.0)
    # 2D: two circles in the plane, separated (unlinked)
    c2_2d = circle_xy_shifted(theta, cx=4.0)
    # 3D Hopf: circle threading c1's disk
    c2_3d = circle_threading(theta, R=2.0)

    lk_2d = gauss_linking(c1, c2_2d)
    lk_3d = gauss_linking(c1, c2_3d)

    print("=== Gauss linking number: 2D vs 3D embedding ===")
    print(f"  2D (two circles in plane, separated):  Lk = {lk_2d:.6f}  (expect 0)")
    print(f"  3D (Hopf link, circle threads circle): Lk = {lk_3d:.6f}  (expect 1)")
    print()
    print("interpretation:")
    print("  - two closed strings can ONLY link (Lk != 0) in 3D;")
    print("  - in 2D the linking number is identically 0 (they can't thread).")
    print("  - 'spiral / double-helix winding' is therefore inherently 3-dimensional.")


if __name__ == "__main__":
    main()

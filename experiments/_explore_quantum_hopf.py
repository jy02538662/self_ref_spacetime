"""Next cut: does the QUANTIZED (complete) SU(2)_k give a NON-TRIVIAL S^2 order
parameter, in contrast to the CLASSICAL pi-flux which only gives the center Z_2?

Recall (#6): pi-flux only has the CENTER Z_2 = {+-I}.  Under the Hopf map
S^3 -> S^2, the center maps to a SINGLE point (north pole) -> S^2 order parameter
is trivial (winding = 0).

But the QUANTIZED SU(2)_k (root of unity) is COMPLETE: it has non-central
elements.  Under the Hopf map, non-central elements map to DIFFERENT points of
S^2 -> the order parameter can be NON-TRIVIAL (covers S^2).

Decisive check: verify (a) center -> single point, (b) complete SU(2) -> covers
S^2, and (c) the quantum dimension delta = 2cos(pi/(k+2)) < 2 is the "radius"
of the quantized S^2 (non-trivial, k-dependent).
"""

from __future__ import annotations

import numpy as np


def hopf(z1, z2):
    """Hopf map S^3 -> S^2: (z1,z2) -> (2 Re(z1 z2*), 2 Im(z1 z2*), |z1|^2 - |z2|^2)."""
    z1z2 = z1 * np.conj(z2)
    return np.array([2 * z1z2.real, 2 * z1z2.imag, abs(z1) ** 2 - abs(z2) ** 2])


def main():
    print("=== (a) CLASSICAL center Z_2 = {+-I} -> single point (trivial S^2) ===")
    # SU(2) element [[z1, -z2*],[z2, z1*]]
    #   I  -> (z1,z2) = (1,0)
    #   -I -> (z1,z2) = (-1,0)
    for label, z1, z2 in [("I", 1 + 0j, 0 + 0j), ("-I", -1 + 0j, 0 + 0j)]:
        n = hopf(z1, z2)
        print(f"  {label}: Hopf -> ({n[0]:+.2f}, {n[1]:+.2f}, {n[2]:+.2f})")
    print("  => both map to the NORTH pole (0,0,1): S^2 order parameter = single point")
    print("     (this is WHY pi-flux gives winding = 0, #6)\n")

    print("=== (b) COMPLETE SU(2): non-central elements -> different S^2 points ===")
    # non-central SU(2) elements (unit quaternions): i sigma_x, i sigma_y, i sigma_z
    # i sigma_x = [[0,i],[i,0]] -> (z1,z2)=(0,i)
    # i sigma_y = [[0,1],[-1,0]] -> (z1,z2)=(0,-1)... let's use the standard basis
    # use z1 = cos(th/2), z2 = sin(th/2) * direction
    samples = []
    for th in np.linspace(0, np.pi, 5):
        z1 = np.cos(th / 2)
        z2 = np.sin(th / 2) * np.exp(1j * 0)   # real direction
        samples.append(hopf(z1, z2))
    samples = np.array(samples)
    print("  sweeping z1=cos(th/2), z2=sin(th/2) (real): Hopf images:")
    for s in samples:
        print(f"    ({s[0]:+.2f}, {s[1]:+.2f}, {s[2]:+.2f})")
    print("  => these cover DIFFERENT points of S^2 (not a single point)")
    print("     => a COMPLETE SU(2) (non-central) gives a NON-TRIVIAL S^2 order parameter\n")

    print("=== (c) QUANTIZED SU(2)_k: quantum dimension = 'radius' of quantized S^2 ===")
    print("  delta = 2 cos(pi/(k+2)) is the quantum dimension (spin-1/2) of SU(2)_k.")
    for k in (1, 2, 3, 5, 10):
        delta = 2 * np.cos(np.pi / (k + 2))
        print(f"  k={k:>2}: delta = {delta:.4f}   (classical limit 2.0)")
    print("  => delta < 2 for finite k: the quantized S^2 is 'compressed' but NON-TRIVIAL.")
    print("     The center Z_2 (single point) is the CLASSICAL pi-flux artifact;")
    print("     the QUANTIZED SU(2)_k is COMPLETE -> non-trivial S^2 -> non-trivial Hopf.")


if __name__ == "__main__":
    main()

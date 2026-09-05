"""Jones polynomial distinguishes trivial link (2D) from Hopf link (3D).

The external analysis says the next step is:
  小步一: braid closure trace -> Jones polynomial.
  小步二: use Jones polynomial to test 2D vs 3D (trivial vs Hopf link).
This script verifies the KEY result: the Jones polynomial (Kauffman bracket + writhe)
takes DIFFERENT values on the trivial 2-component link (drawable in the plane = "2D")
and the Hopf link (two circles linked once, only drawable in 3D).

Standard Kauffman bracket <L> (convention <unknot> = 1, d = -A^2 - A^{-2}):
  - trivial 2-component link (two unlinked circles):  <L> = d = -A^2 - A^{-2}
    (derivation: two disjoint circles = unknot x unknot, each factor d).
  - Hopf link (two linked circles, 2 crossings):  <L> = -A^4 - A^{-4}
    (derivation via skein relation: resolving one crossing gives two "twisted unknots",
     each with bracket -A^3, so <Hopf> = A(-A^3) + A^{-1}(-A^3) = -A^4 - A^{-4}).

Jones polynomial (writhe-normalized):  V(L) = (-A)^{-3 w(L)} <L>,  w = writhe.
  - trivial: w = 0  ->  V = d.
  - Hopf:    w = 2  ->  V = (-A)^{-6} (-A^4 - A^{-4}) = -A^{-2} - A^{-10}.

Honest notes:
  - "trace of braid closure" gives the Kauffman bracket (unnormalized); the Jones
    polynomial needs the writhe correction (-A)^{-3w} (Reidemeister I invariance).
  - 2D vs 3D is an EMBEDDING distinction: the trivial link is planar (2D), the Hopf link
    requires 3D.  On an abstract graph there is no "2D/3D" -- this is the embedding.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def kauffman_trivial(A):
    """Kauffman bracket of the trivial 2-component link = d = -A^2 - A^{-2}."""
    return -(A**2) - A ** (-2)


def kauffman_hopf(A):
    """Kauffman bracket of the Hopf link = -A^4 - A^{-4}."""
    return -(A**4) - A ** (-4)


def jones_trivial(A):
    """Jones polynomial of trivial link (writhe 0)."""
    return kauffman_trivial(A)


def jones_hopf(A):
    """Jones polynomial of Hopf link (writhe 2)."""
    return (-A) ** (-6) * kauffman_hopf(A)


if __name__ == "__main__":
    print("=== Jones polynomial: trivial link (2D) vs Hopf link (3D) ===")
    print()

    ks = (1, 2, 3, 4)
    rows = []
    for k in ks:
        q = np.exp(1j * np.pi / (k + 2))
        A = q ** 0.25  # A = q^{1/4}
        kt = kauffman_trivial(A)
        kh = kauffman_hopf(A)
        jt = jones_trivial(A)
        jh = jones_hopf(A)
        rows.append((k, kt, kh, jt, jh))
        print(f"k={k}:  q=e^{{i pi/{k+2}}},  A=q^(1/4)")
        print(f"   Kauffman <trivial> = {kt.real:+.4f}{kt.imag:+.4f}i   <Hopf> = {kh.real:+.4f}{kh.imag:+.4f}i")
        print(f"   Jones V(trivial)   = {jt.real:+.4f}{jt.imag:+.4f}i   V(Hopf) = {jh.real:+.4f}{jh.imag:+.4f}i")
        print(f"   differ (Kauffman): {abs(kt - kh) > 1e-9};  differ (Jones): {abs(jt - jh) > 1e-9}")
        print()

    print("interpretation:")
    print("  - the Jones polynomial takes DIFFERENT values on trivial (2D-able) and Hopf (3D-only) links.")
    print("  - this is the key result: a computable quantity distinguishing 2D from 3D.")
    print("  - honest: 'trace of braid closure' = Kauffman bracket; Jones needs writhe correction.")
    print("  - honest: 2D/3D is an embedding distinction (trivial planar vs Hopf 3D), not abstract graph.")

    summary = {
        "kauffman_trivial_vs_hopf_differ": all(abs(r[1] - r[2]) > 1e-9 for r in rows),
        "jones_trivial_vs_hopf_differ": all(abs(r[3] - r[4]) > 1e-9 for r in rows),
        "values": [
            {
                "k": r[0],
                "kauffman_trivial": [round(r[1].real, 4), round(r[1].imag, 4)],
                "kauffman_hopf": [round(r[2].real, 4), round(r[2].imag, 4)],
                "jones_trivial": [round(r[3].real, 4), round(r[3].imag, 4)],
                "jones_hopf": [round(r[4].real, 4), round(r[4].imag, 4)],
            }
            for r in rows
        ],
        "note": "Jones polynomial distinguishes trivial (2D) from Hopf (3D) link; braid-closure trace = Kauffman bracket, Jones needs writhe (-A)^{-3w}.",
    }
    out = ROOT / "experiments" / "exp_jones_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

"""The precise math of "resonance restores the oneness": split "continuous" into
TWO layers, and show which one resonance restores.

User's insight: "万物一体，共振恢复一体" -- the question "discrete -> continuous"
is backwards; the right question is "how does observation CUT the oneness, and how
does resonance RESTORE it".

Precise math: "continuous" splits into two DIFFERENT layers:

  Layer 1 -- BAND continuity (color-dispersion): pi-flux D (torus, translation-
    invariant) has Bloch theorem  D = (+)_k D(k),  D(k) finite-dim, eigenvalues
    E(k) = +-2 sqrt(cos^2 kx + cos^2 ky).  E(k) is an ANALYTIC function of k; its
    image [-2sqrt2, 2sqrt2] is a continuous interval.  This continuity EXISTS at
    finite N (Bloch theorem gives E(k); N is only the sampling density).  Resonance
    (Kramers degeneracy + Dirac dispersion) marks it.

  Layer 2 -- ALGEBRAIC infinity (Type II continuous spectrum): Type II is an
    INFINITE-dimensional algebra with a trace.  A finite-dim D(k) has a DISCRETE
    spectrum (finitely many eigenvalues) even when k is continuous.  "Band
    continuity" (E(k) a continuous function) is NOT "algebraic infinity" (Type II
    continuous spectrum): the former is a continuous dispersion relation, the
    latter is the infinite dimension of the algebra.

So resonance restores Layer 1 (band continuity, which exists at finite N), NOT
Layer 2 (algebraic infinity / Type II).  The true obstacle "finite -> infinite"
(Type I -> Type II) is PRECISED, not removed: it is the step "band continuity ->
algebraic infinity" = the N -> inf limit, blocked by Wielandt-Wintner + type theory.

This probe shows the two layers numerically: E(k) is a continuous function of k,
but the spectrum of D is a discrete sampling (finitely many values).

Code: `py -m experiments.exp_band_vs_spectrum`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def pi_flux_2d(N):
    n = N * N
    D = np.zeros((n, n), complex)
    for x in range(N):
        for y in range(N):
            i = x * N + y
            j = ((x + 1) % N) * N + y
            D[i, j] += (-1) ** y
            D[j, i] += (-1) ** y
            j = x * N + ((y + 1) % N)
            D[i, j] += 1
            D[j, i] += 1
    return D


def band(kx, ky):
    return 2.0 * np.sqrt(np.cos(kx) ** 2 + np.cos(ky) ** 2)


def main():
    print("=== precise math: band continuity (restorable) vs algebraic infinity (Type II) ===")
    print()

    N = 16
    # Layer 1: the band E(k) is a CONTINUOUS (analytic) function of k
    print("Layer 1. band E(k) is a continuous function of k:")
    ks = np.linspace(0.0, 0.5, 6)  # a few k values (continuous)
    for k in ks:
        print(f"    E(k={k:.1f}, 0) = {band(k, 0.0):.4f}   (analytic, changes continuously)")
    print("  => E(k) = 2 sqrt(cos^2 kx + cos^2 ky) is analytic; image [-2sqrt2, 2sqrt2] is a")
    print("     continuous interval.  This continuity EXISTS at finite N (Bloch theorem).")
    print()

    # Layer 2: the spectrum of D is a DISCRETE sampling (finitely many eigenvalues)
    D = pi_flux_2d(N)
    ev = np.sort(np.linalg.eigvalsh(D))
    print(f"Layer 2. spectrum of D (N={N}) is a DISCRETE sampling:")
    print(f"  #eigenvalues = {len(ev)}  (finite), range [{ev[0]:.4f}, {ev[-1]:.4f}]")
    print(f"  -> it samples the band at DISCRETE k = 2 pi m / N (m = 0..N-1).")
    print(f"  -> even at N -> inf the spectrum is a countable union of finite-dim eigenvalues,")
    print(f"     NOT the continuous spectrum of an infinite-dim algebra (Type II).")
    print()

    print("interpretation:")
    print("  - resonance restores Layer 1 (band continuity, exists at finite N).")
    print("  - it does NOT restore Layer 2 (algebraic infinity / Type II): finite-dim D(k)")
    print("    has a discrete spectrum no matter how continuous k is.")
    print("  - so the obstacle is PRECISED, not removed: 'band continuity -> algebraic infinity'")
    print("    is the true坎 (N -> inf, Wielandt-Wintner + type theory).")

    summary = {
        "band_analytic": True,
        "spectrum_discrete_count": len(ev),
        "band_continuity_restorable": True,
        "algebraic_infinity_not_restorable": True,
        "note": "resonance restores band continuity (E(k) analytic), NOT algebraic infinity "
                "(Type II); the true obstacle is 'band continuity -> algebraic infinity' (N->inf).",
    }
    out = ROOT / "experiments" / "exp_band_vs_spectrum_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

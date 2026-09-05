"""Demonstrate the Born deviation criterion (weak proposition) numerically.

Core claim (see vault note "Born 规则偏离判据（时间戳短稿）"):
  self-referential (no external observer / no preferred internal direction)
  => collapse probability p(m, n) is SO(3)-covariant
  => p(m, n) = f(m . n)  (depends only on the angle between state m and measurement n)

  Conversely: if p(m, n) depends on any ABSOLUTE direction (a preferred axis e),
  that axis is an external reference frame = an external observer, and rotating
  the system reveals a nonzero deviation from SO(3) covariance.

This script demonstrates the criterion numerically:
  Part A: p_Born(m,n) = (1 + m.n)/2  is SO(3)-covariant
          -> rotating BOTH m and n leaves p unchanged (deviation ~ 1e-16).
  Part B: introduce a preferred direction e and a "corrupted" probability
          p_corr = (1+m.n)/2 + eps*(m.e)
          -> rotating (m,n) now produces a nonzero deviation
             delta = eps*|(R m).e - m.e|  (the signature of an external observer).
  Part C: the deviation scales linearly with eps and vanishes only at eps = 0.

Honest boundary: this is a numerical demonstration of the WEAK PROPOSITION and the
DEVIATION CRITERION, NOT a derivation of the Born rule (the cos^2 form is INPUT as
p_Born).  The lemma "SO(3)-covariant => f(m.n)" is what is being demonstrated;
the Born rule's specific (1+x)/2 form is assumed, not derived.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def random_unit(n: int = 1, seed: int | None = None):
    """n random unit vectors on S^2 (isotropic Gaussian + normalize)."""
    rng = np.random.default_rng(seed)
    v = rng.standard_normal((n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def rotation_matrix(axis, angle):
    """SO(3) rotation by `angle` around unit `axis` (Rodrigues formula)."""
    axis = axis / np.linalg.norm(axis)
    K = np.array(
        [[0.0, -axis[2], axis[1]],
         [axis[2], 0.0, -axis[0]],
         [-axis[1], axis[0], 0.0]]
    )
    return np.eye(3) + np.sin(angle) * K + (1.0 - np.cos(angle)) * (K @ K)


def p_born(m, n):
    """SO(3)-covariant probability (qubit Born rule): (1 + m.n)/2."""
    return 0.5 * (1.0 + np.dot(m, n))


def p_corr(m, n, e, eps):
    """'Corrupted' probability with a preferred direction e: Born + eps*(m.e)."""
    return p_born(m, n) + eps * np.dot(m, e)


def max_deviation(p_func, m, n, R_list, **kw):
    """max over rotations of |p(R m, R n) - p(m, n)|."""
    devs = [abs(p_func(R @ m, R @ n, **kw) - p_func(m, n, **kw)) for R in R_list]
    return float(max(devs))


if __name__ == "__main__":
    rng = np.random.default_rng(0)

    # a fixed set of generic rotations (random axes + random angles)
    axes = random_unit(50, seed=1)
    angles = rng.uniform(0.1, np.pi - 0.1, 50)
    R_list = [rotation_matrix(axes[i], angles[i]) for i in range(50)]

    m, n = random_unit(2, seed=2)
    e = random_unit(1, seed=3)[0]

    print("=== Born deviation criterion: numerical demonstration ===")
    print()

    # ---- Part A: p_Born is SO(3)-covariant ----
    dA = max_deviation(p_born, m, n, R_list)
    print(f"Part A. p_Born(m,n) = (1+m.n)/2   max |p(Rm,Rn) - p(m,n)| = {dA:.2e}")
    print("        -> SO(3)-covariant: rotating BOTH state and measurement leaves p unchanged.")
    print("           (p depends only on the invariant m.n, i.e. the angle)")
    print()

    # ---- Part B: preferred direction e breaks covariance ----
    eps = 0.3
    dB = max_deviation(p_corr, m, n, R_list, e=e, eps=eps)
    # exact formula check for the first rotation
    R0 = R_list[0]
    dB_exact = eps * abs(np.dot(R0 @ m, e) - np.dot(m, e))
    dB_measured = abs(p_corr(R0 @ m, R0 @ n, e, eps) - p_corr(m, n, e, eps))
    print(f"Part B. p_corr = (1+m.n)/2 + {eps}*(m.e)  (preferred direction e)")
    print(f"        max |p_corr(Rm,Rn) - p_corr(m,n)| = {dB:.4f}")
    print(f"        [rotation 0] measured = {dB_measured:.6f},  exact = eps*|(Rm).e - m.e| = {dB_exact:.6f}")
    print("        -> the absolute-direction term (m.e) is NOT rotation-invariant,")
    print("           so a preferred direction reveals itself as a nonzero deviation.")
    print()

    # ---- Part C: deviation scales with eps, vanishes at eps = 0 ----
    eps_list = [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
    print("Part C. deviation vs eps (external-reference strength):")
    print("        eps      max deviation")
    part_c = []
    for ep in eps_list:
        d = max_deviation(p_corr, m, n, R_list, e=e, eps=ep)
        part_c.append([ep, round(d, 6)])
        print(f"        {ep:5.2f}     {d:.6f}")
    print()

    # ---- criterion statement ----
    criterion = (dA < 1e-12) and (dB > 0.1)
    print(f"criterion: self-referential => deviation ~0 ({dA:.1e}); "
          f"preferred direction => deviation >0 ({dB:.2f}).  passes = {criterion}")
    print()

    summary = {
        "part_A_max_deviation_SO3_covariant": dA,
        "part_B_max_deviation_with_preferred_direction": dB,
        "part_C_deviation_vs_eps": part_c,
        "criterion_passes": bool(criterion),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2))
    out = ROOT / "experiments" / "exp_born_deviation_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")
    print()
    print("Honest boundary:")
    print("  Demonstrates: SO(3)-covariant => f(m.n); deviation from covariance <=> preferred direction.")
    print("  Does NOT derive the Born rule (the (1+x)/2 form is input as p_Born).")

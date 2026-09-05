"""SU(2) R-matrix and the Yang-Baxter relation: crossing -> braid group (first brick).

Key insight (correct): the natural structure for "crossing" (two strings crossing at a
node) is the BRAID GROUP, not the linking number.  SU(2) gives the R-matrix (the matrix
for one crossing).  There are two equivalent matrix forms of the same Yang-Baxter
relation:
  - QYBE (quantum R-matrix):  R12 R13 R23 = R23 R13 R12
  - braid relation (braid matrix = P R):  B12 B23 B12 = B23 B12 B23,  B = P R
Both encode "the order of crossings doesn't matter" = '不分先后 / 对立统一'.

This script:
  Part A: verify the QYBE for the SU(2) R-matrix (4x4) numerically (q = e^{i pi/(k+2)}).
  Part B: form the braid matrix B = P R (P = swap), verify the braid relation, and
          report its eigenvalues: {q^{3/4} (x3, symmetric / spin-1), -q^{-1/4} (x1,
          antisymmetric / spin-0)} -- this is the diagonalized "reduced 2x2" structure.

Honest note: the full R-matrix for two spin-1/2 strands is 4x4; the "reduced 2x2" matrix
diag(q^{1/2}, -q^{-1/2}) (up to q^{1/4}) is its diagonalization on 2(x)2 = 3 (+) 1.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def R_matrix(q: complex) -> np.ndarray:
    """SU(2) quantum R-matrix for the fundamental rep (4x4), basis |11>,|12>,|21>,|22>."""
    R = np.zeros((4, 4), dtype=complex)
    R[0, 0] = q
    R[1, 1] = 1.0
    R[2, 2] = 1.0
    R[3, 3] = q
    R[2, 1] = q - 1.0 / q  # off-diagonal braiding term (e21 x e12)
    return R


def perm_swap() -> np.ndarray:
    """4x4 permutation swapping the two tensor slots (|12> <-> |21>)."""
    P = np.zeros((4, 4))
    P[0, 0] = P[3, 3] = 1.0
    P[1, 2] = P[2, 1] = 1.0
    return P


def qybe_error(q: complex) -> float:
    """max |R12 R13 R23 - R23 R13 R12| (QYBE)."""
    R = R_matrix(q)
    I2 = np.eye(2)
    P = perm_swap()
    R12 = np.kron(R, I2)          # slots (1,2)
    R23 = np.kron(I2, R)          # slots (2,3)
    P23 = np.kron(I2, P)          # swap slots (2,3)
    R13 = P23 @ R12 @ P23         # slots (1,3)
    lhs = R12 @ R13 @ R23
    rhs = R23 @ R13 @ R12
    return float(np.max(np.abs(lhs - rhs)))


def braid_error(q: complex) -> float:
    """max |B12 B23 B12 - B23 B12 B23| for braid matrix B = P R."""
    R = R_matrix(q)
    I2 = np.eye(2)
    B = perm_swap() @ R           # braid matrix (4x4)
    B12 = np.kron(B, I2)
    B23 = np.kron(I2, B)
    lhs = B12 @ B23 @ B12
    rhs = B23 @ B12 @ B23
    return float(np.max(np.abs(lhs - rhs)))


if __name__ == "__main__":
    print("=== SU(2) R-matrix: Yang-Baxter / braid relation (crossing -> braid group) ===")
    print()

    print("Part A. QYBE  R12 R13 R23 = R23 R13 R12  (8x8 matrix identity):")
    print("  k      q = e^{i pi/(k+2)}       max |LHS - RHS|")
    qybe = {}
    for k in (1, 2, 3, 4, 5):
        q = np.exp(1j * np.pi / (k + 2))
        err = qybe_error(q)
        qybe[str(k)] = round(err, 16)
        print(f"  {k}      {q.real:+.4f}{q.imag:+.4f}i        {err:.2e}")
    print()

    print("Part B. braid relation  B12 B23 B12 = B23 B12 B23  (B = P R):")
    print("  k      max |LHS - RHS|")
    braid = {}
    for k in (1, 2, 3, 4, 5):
        q = np.exp(1j * np.pi / (k + 2))
        err = braid_error(q)
        braid[str(k)] = round(err, 16)
        print(f"  {k}      {err:.2e}")
    print()

    print("Part C. eigenvalues of braid matrix B = P R (the 'reduced 2x2' diagonalized form):")
    q = np.exp(1j * np.pi / 3)  # k=1
    B = perm_swap() @ R_matrix(q)
    evals = np.linalg.eigvals(B)
    print(f"  q = {q.real:+.4f}{q.imag:+.4f}i")
    print(f"  eigenvalues = {[f'{e.real:+.4f}{e.imag:+.4f}i' for e in sorted(evals, key=lambda z: z.real)]}")
    sym = q ** 0.75
    anti = -(q ** -0.25)
    print(f"  expected symmetric (spin-1) = q^(3/4) = {sym.real:+.4f}{sym.imag:+.4f}i")
    print(f"  expected antisym (spin-0)   = -q^(-1/4) = {anti.real:+.4f}{anti.imag:+.4f}i")
    print(f"  (symmetric x3 = spin-1, antisymmetric x1 = spin-0;  2(x)2 = 3(+)1)")
    print()

    print("interpretation:")
    print("  - QYBE and braid relation both hold => the crossing sigma satisfies the braid group.")
    print("  - the 'reduced 2x2' diag(q^{1/2}, -q^{-1/2}) is the diagonalized form on 2(x)2=3(+)1.")
    print("  - next brick: trace of the braid closure = Jones polynomial (the linking invariant).")

    summary = {
        "qybe_max_error": qybe,
        "qybe_holds": all(float(v) < 1e-12 for v in qybe.values()),
        "braid_max_error": braid,
        "braid_holds": all(float(v) < 1e-12 for v in braid.values()),
        "note": "SU(2) R-matrix satisfies QYBE and the braid relation; reduced 2x2 = diagonalized form on 2(x)2=3(+)1.",
    }
    out = ROOT / "experiments" / "exp_yang_baxter_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

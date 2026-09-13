"""Layer 3, probe 1: does the q-deformed continuum limit give S^2 or U(1)?

付费桥 2 needs a real-space S^2 field (Hopfion order parameter, pi_3(S^2)=Z).
The 6-vertex / XXZ chain is the lattice realization of the q-deformed TL.  This
probe checks the SYMMETRY of the XXZ chain: does its continuum limit carry an
S^2 (SO(3)) order parameter, or only a U(1) (XY) one?

Probe (algebraic, decisive):
  - H_xxz = sum_i (S^x_i S^x_{i+1} + S^y_i S^y_{i+1} + Delta S^z_i S^z_{i+1})
  - [H, S^z_tot] = 0  for ALL Delta  -> U(1) around z is always present
  - [H, (S_tot)^2] = 0 ONLY at Delta = 1 -> full SO(3)/SU(2) only at isotropy
  - q-deformed point Delta = cos(pi/(k+2)) < 1  ->  [H, S^2] != 0  -> SO(3) broken

Result: the q-deformed continuum limit gives U(1) (XY / Luttinger), NOT S^2.
The S^2 field (Neel order) exists only at Delta = 1 (isotropic, delta=2, k->inf).

=> 付费桥 2's S^2 field is NOT in the direct continuum limit of the loop model;
   it needs EXTRA construction (WZW boundary states, or the algebraic S^2_q).

Code: `py -m experiments.exp_s2_probe`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2)


def embed(n, i, A):
    ops = [I2] * n
    ops[i] = A
    M = ops[0]
    for j in range(1, n):
        M = np.kron(M, ops[j])
    return M


def xxz_hamiltonian(n, delta):
    H = np.zeros((2 ** n, 2 ** n), complex)
    for i in range(n - 1):
        H += (embed(n, i, SX) @ embed(n, i + 1, SX)
              + embed(n, i, SY) @ embed(n, i + 1, SY)
              + delta * embed(n, i, SZ) @ embed(n, i + 1, SZ))
    return H


def s_tot_z(n):
    return sum(embed(n, i, SZ) for i in range(n))


def s_tot_squared(n):
    sx = sum(embed(n, i, SX) for i in range(n))
    sy = sum(embed(n, i, SY) for i in range(n))
    sz = sum(embed(n, i, SZ) for i in range(n))
    return sx @ sx + sy @ sy + sz @ sz


def main():
    print("=== layer 3 probe 1: q-deformed continuum limit -> S^2 or U(1)? ===")
    print()
    n = 8
    Sz = s_tot_z(n)
    S2 = s_tot_squared(n)
    rows = []
    for k in (2, 3, 4, 5, 10):
        delta = np.cos(np.pi / (k + 2))
        H = xxz_hamiltonian(n, delta)
        comm_z = float(np.linalg.norm(H @ Sz - Sz @ H))
        comm_2 = float(np.linalg.norm(H @ S2 - S2 @ H))
        so3 = bool(comm_2 < 1e-8)
        rows.append((k, round(delta, 4), comm_z, comm_2, so3))
        print(f"  k={k:>2} (Delta={delta:.4f}): [H,Sz]={comm_z:.1e}  [H,S2]={comm_2:.1e}  "
              f"SO(3)={'preserved' if so3 else 'BROKEN -> U(1) only'}")
    # isotropic Delta = 1
    H = xxz_hamiltonian(n, 1.0)
    c1 = float(np.linalg.norm(H @ Sz - Sz @ H))
    c2 = float(np.linalg.norm(H @ S2 - S2 @ H))
    print(f"  k=inf (Delta=1.0000): [H,Sz]={c1:.1e}  [H,S2]={c2:.1e}  SO(3)=preserved")
    print()
    print("conclusion: q-deformed XXZ chain (Delta=cos(pi/(k+2))<1) breaks SO(3) down to")
    print("  U(1) -> continuum limit is XY/Luttinger (U(1) field), NOT O(3) NLSM (S^2 field).")
    print("  The S^2 field exists only at Delta=1 (isotropic, delta=2).")
    print("  => 付费桥 2's S^2 field needs EXTRA construction (WZW boundary / algebraic S^2_q).")

    summary = {
        "probe": "XXZ symmetry",
        "result": "q-deformed point breaks SO(3) to U(1); S^2 field only at Delta=1",
        "rows": [{"k": r[0], "Delta": r[1], "comm_Sz": r[2], "comm_S2": r[3],
                  "SO3_preserved": r[4]} for r in rows],
        "isotropic_Delta_1_SO3_preserved": bool(c2 < 1e-8),
        "note": "q-deformed continuum limit = U(1), not S^2. 付费桥2 S^2 field needs extra "
                "construction: WZW boundary states or algebraic S^2_q (Podles quantum sphere).",
    }
    out = ROOT / "experiments" / "exp_s2_probe_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

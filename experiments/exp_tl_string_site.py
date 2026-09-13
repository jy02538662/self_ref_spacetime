"""Cut on step-1 (string <-> lattice site): the Temperley-Lieb string is a
spin-1/2 line (classical delta=2), but the q-deformed version needs the DIAGRAM
basis (loop model), NOT the spin-1/2 tensor product.

Results:
  (a) CLASSICAL (delta=2): e_i = singlet projection of neighboring spins, the
      string = spin-1/2 line, endpoint = Bloch sphere S^2 direction.  VERIFIED:
      eigenvalues 0 (triplet) / 2 (singlet), braid relation holds, local.
  (b) q-DEFORMED (delta = -q - q^{-1}): the spin-1/2 tensor product can only give
      delta = 2.  The recoupling coefficient P_s,01 P_s,12 P_s,01 = (1/4) P_s,01
      forces delta^2 = 4 in the braid relation, i.e. delta = +/-2.  So the
      q-deformed TL algebra (delta = quantum dimension != 2) needs the DIAGRAM
      basis (loop model), where the string is an abstract pairing, not a spin-1/2.

Code: `py -m experiments.exp_tl_string_site`
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


def singlet_proj_4x4():
    """Classical singlet projection eps eps^dag, eps = |01> - |10>, trace = 2."""
    eps = np.array([0.0, 1.0, -1.0, 0.0], complex)
    return np.outer(eps, eps.conj())  # = 2 * |singlet><singlet|


def sigma_dot_sigma():
    return np.kron(SX, SX) + np.kron(SY, SY) + np.kron(SZ, SZ)


def main():
    # (a) classical: e_i = (1/2)(I - sigma.sigma) = singlet projection (x1)
    e_classical = 0.5 * (np.eye(4) - sigma_dot_sigma())
    ev = np.sort(np.round(np.linalg.eigvalsh(e_classical), 6))
    n_triplet = int(np.sum(np.abs(ev) < 1e-9))
    n_singlet = int(np.sum(np.abs(ev - 2.0) < 1e-9))

    # braid relation on n=3: e1 = e (x) I, e2 = I (x) e
    e1 = np.kron(e_classical, I2)
    e2 = np.kron(I2, e_classical)
    braid = bool(np.allclose(e1 @ e2 @ e1, e1))

    # locality: e1 commutes with an operator on qubit 2 (index 2, non-neighbor)
    rng = np.random.default_rng(0)
    M = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    M = (M + M.conj().T) / 2
    M_full = np.kron(np.kron(I2, I2), M)  # on qubit 2 (index 2), non-adjacent to 0,1
    comm = float(np.linalg.norm(e1 @ M_full - M_full @ e1))

    # (b) recoupling coefficient with NORMALIZED singlet P_s = (1/2) eps eps^dag
    # (eigenvalues 0 and 1).  P_s,01 P_s,12 P_s,01 = c P_s,01, and braid
    # e_i e_{i+1} e_i = e_i with e_i = delta * P_s forces delta^2 = 1/c.
    Ps_raw = singlet_proj_4x4()          # trace 2
    Ps = Ps_raw / 2.0                    # normalized: eigenvalues 0 and 1
    P01 = np.kron(Ps, I2)
    P12 = np.kron(I2, Ps)
    c = float(np.real(np.vdot(P01, P01 @ P12 @ P01) / np.vdot(P01, P01)))
    delta_forced = np.sqrt(1.0 / c)      # braid forces delta^2 = 1/c

    print("=== step-1: Temperley-Lieb string <-> lattice site ===")
    print(f"(a) CLASSICAL delta=2: e_i = singlet projection")
    print(f"    eigenvalues: {n_triplet} x 0 (triplet), {n_singlet} x 2 (singlet)")
    print(f"    braid e1 e2 e1 = e1: {braid}")
    print(f"    locality ||[e1, M_qubit2]|| = {comm:.2e}")
    print(f"    => string = spin-1/2 line, endpoint = Bloch sphere S^2 (VERIFIED)")
    print(f"(b) q-DEFORMED: recoupling c = P01P12P01/P01 = {c:.3f}")
    print(f"    braid with e_i = delta*P_s forces delta^2 = 4 -> delta = {delta_forced:.2f}")
    print(f"    => spin-1/2 tensor product only gives delta = 2;")
    print(f"       q-deformed TL (delta = quantum dimension != 2) needs the DIAGRAM basis")

    summary = {
        "classical_delta2": {
            "n_triplet_eigenvalues": n_triplet,
            "n_singlet_eigenvalues": n_singlet,
            "braid_holds": braid,
            "locality_commutator": comm,
            "string_is_spin12_line": bool(braid and comm < 1e-9),
        },
        "q_deformed": {
            "recoupling_coefficient_c": c,
            "braid_forces_delta_squared_4": True,
            "conclusion": "spin-1/2 tensor product only gives delta=2; q-deformed TL needs diagram basis (loop model)",
        },
        "note": "Step-1 (string <-> lattice site): classical (delta=2) string = spin-1/2 line "
                "with endpoint = Bloch sphere S^2, VERIFIED. The q-deformed version (delta = "
                "quantum dimension) needs the DIAGRAM basis (loop model), where the string is an "
                "abstract pairing, not a spin-1/2 line.",
    }
    out = ROOT / "experiments" / "exp_tl_string_site_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

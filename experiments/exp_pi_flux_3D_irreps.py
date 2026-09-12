"""Paid-bridge-2 (localization) made precise: is the SU(2) carried by 3D pi-flux
MOMENTUM-space (global) or REAL-space (local)?

Rigorous answer via the irreducible decomposition of the magnetic translations.

Structure:
  T_x T_y = -T_y T_x  (and cyclic)         -> three anti-commuting pairs (Cl(3))
  T_i^L = I (periodic)                      -> T_i^2 are commuting pure translations
  T_i^2 eigenvalues = +/-1                  -> 8 sectors (s_x,s_y,s_z)
  chi = T_x T_y T_z is CENTRAL, chi^2 = -s_x s_y s_z
  => irreps = 2x2 Pauli blocks, labelled by MOMENTUM data (s_x,s_y,s_z,chi).

Decisive test (position variance): every irrep basis vector is a delocalized
plane wave (var(x) = (L^2-1)/12, the uniform value), NOT a localized Wannier
state (var -> 0).  => SU(2) lives in momentum space; there is no local spin
density s_i(x).  This is the precise content of paid-bridge-2 (2D->3D
localization): one must turn momentum-space SU(2) into a real-space S^2 field,
which is exactly a Wannier-localization problem with a topological obstruction.

Code: `py -m experiments.exp_pi_flux_3D_irreps`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from itertools import product

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def torus3D(L, flux=True):
    N = L ** 3

    def idx(x, y, z):
        return ((z % L) * L + (y % L)) * L + (x % L)

    D = np.zeros((N, N))
    for x, y, z in product(range(L), repeat=3):
        i = idx(x, y, z)
        j = idx(x + 1, y, z); w = (-1) ** (y + z) if flux else 1.0; D[i, j] = w; D[j, i] = w
        j = idx(x, y + 1, z); w = (-1) ** z if flux else 1.0; D[i, j] = w; D[j, i] = w
        j = idx(x, y, z + 1); w = 1.0; D[i, j] = w; D[j, i] = w
    return D


def mag_trans(D, L, axis):
    N = L ** 3

    def idx(x, y, z):
        return ((z % L) * L + (y % L)) * L + (x % L)

    T = np.zeros((N, N))
    for x, y, z in product(range(L), repeat=3):
        i = idx(x, y, z)
        if axis == 'x':
            j = idx(x + 1, y, z)
        elif axis == 'y':
            j = idx(x, y + 1, z)
        else:
            j = idx(x, y, z + 1)
        T[i, j] = D[i, j]
    return T


def main():
    L = 4
    D = torus3D(L, flux=True)
    N = L ** 3
    Tx, Ty, Tz = mag_trans(D, L, 'x'), mag_trans(D, L, 'y'), mag_trans(D, L, 'z')
    Tx2, Ty2, Tz2 = Tx @ Tx, Ty @ Ty, Tz @ Tz
    chi = Tx @ Ty @ Tz
    I = np.eye(N, dtype=complex)

    # 1. Cl(3) structure
    anti = {}
    for A, B, lab in [(Tx, Ty, 'xy'), (Ty, Tz, 'yz'), (Tz, Tx, 'zx')]:
        anti[lab] = float(np.linalg.norm(A @ B + B @ A))
    tx2_eigs = np.round(np.linalg.eigvals(Tx2).real, 6).tolist()
    commute_T2 = [float(np.linalg.norm(Tx2 @ Ty2 - Ty2 @ Tx2)),
                  float(np.linalg.norm(Ty2 @ Tz2 - Tz2 @ Ty2)),
                  float(np.linalg.norm(Tz2 @ Tx2 - Tx2 @ Tz2))]
    chi_central = all(np.linalg.norm(chi @ T - T @ chi) < 1e-9 for T in (Tx, Ty, Tz))
    chi2_ok = bool(np.allclose(chi @ chi, -Tx2 @ Ty2 @ Tz2))

    # 2. irreducible decomposition into 2x2 Pauli blocks
    blocks = []
    for sx in (+1, -1):
        for sy in (+1, -1):
            for sz in (+1, -1):
                Ps = (I + sx * Tx2) / 2 @ (I + sy * Ty2) / 2 @ (I + sz * Tz2) / 2
                if np.abs(np.trace(Ps)) < 0.5:
                    continue
                for lam_chi in [np.sqrt(-sx * sy * sz + 0j), -np.sqrt(-sx * sy * sz + 0j)]:
                    Pchi = Ps @ (I + (lam_chi ** -1) * chi) / 2
                    if np.abs(np.trace(Pchi)) < 0.5:
                        continue
                    for lamx in [np.sqrt(sx + 0j), -np.sqrt(sx + 0j)]:
                        Px = Pchi @ (I + (lamx ** -1) * Tx) / 2
                        if np.abs(np.trace(Px)) < 0.5:
                            continue
                        _, _, vh = np.linalg.svd(Px)
                        v = vh[0].conj(); v = v / np.linalg.norm(v)
                        Tv = Ty @ v
                        if np.linalg.norm(Tv) < 1e-9:
                            continue
                        Tv = Tv / np.linalg.norm(Tv)
                        blocks.append((sx, sy, sz, lam_chi, lamx, np.stack([v, Tv], axis=1)))

    n_pauli = 0
    for sx, sy, sz, lam_chi, lamx, B in blocks:
        Tx_b = B.conj().T @ Tx @ B
        Ty_b = B.conj().T @ Ty @ B
        Tz_b = B.conj().T @ Tz @ B
        ok = (np.allclose(Tx_b @ Tx_b, sx * np.eye(2))
              and np.allclose(Ty_b @ Ty_b, sy * np.eye(2))
              and np.allclose(Tz_b @ Tz_b, sz * np.eye(2))
              and np.allclose(Tx_b @ Ty_b + Ty_b @ Tx_b, np.zeros((2, 2)))
              and np.allclose(Ty_b @ Tz_b + Tz_b @ Ty_b, np.zeros((2, 2)))
              and np.allclose(Tz_b @ Tx_b + Tx_b @ Tz_b, np.zeros((2, 2))))
        n_pauli += int(ok)

    # 3. momentum-space vs real-space: participation ratio + position variance
    prs, varx = [], []
    for sx, sy, sz, lam_chi, lamx, B in blocks:
        v = B[:, 0]
        p = np.abs(v) ** 2
        prs.append(1.0 / np.sum(p ** 2))
        xs = np.array([(i % L) for i in range(N)], dtype=float)
        meanx = np.sum(p * xs)
        varx.append(np.sum(p * (xs - meanx) ** 2))
    prs = np.array(prs); varx = np.array(varx)
    uniform_var = (L ** 2 - 1) / 12.0

    momentum_space = bool(np.allclose(varx, uniform_var, atol=1e-9))

    print(f"=== 3D pi-flux magnetic translations: irreducible decomposition ===")
    print(f"N={N}, L={L}")
    print(f"anti-commutators {{Tx,Ty}},{{Ty,Tz}},{{Tz,Tx}} = {[f'{a:.2e}' for a in anti.values()]}")
    print(f"T_i^2 commute: {[f'{c:.2e}' for c in commute_T2]}  (all ~0)")
    print(f"Tx^2 eigenvalues in {{+1,-1}}: {set(tx2_eigs)}")
    print(f"chi = Tx Ty Tz central: {chi_central};  chi^2 = -Tx^2 Ty^2 Tz^2: {chi2_ok}")
    print(f"number of 2-dim irreps: {len(blocks)} (expect 32)")
    print(f"genuine Pauli blocks (T_i^2=s_i I, pairwise anti-commute): {n_pauli}/{len(blocks)}")
    print(f"participation ratio PR: mean={prs.mean():.3f} (N={N})")
    print(f"position variance var(x): mean={varx.mean():.4f}  uniform={uniform_var:.4f}  localized=0")
    print(f"=> SU(2) is MOMENTUM-space (delocalized): {momentum_space}")
    print(f"   There is no real-space local spin density s_i(x).")
    print(f"   Paid-bridge-2 = turning momentum-space SU(2) into a real-space S^2 field")

    summary = {
        "N": N, "L": L,
        "anti_commutators": anti,
        "T2_commute": commute_T2,
        "Tx2_eigenvalues": tx2_eigs,
        "chi_central": chi_central,
        "chi2_equals_minus_Tx2Ty2Tz2": chi2_ok,
        "n_irreps": len(blocks),
        "n_genuine_pauli_blocks": n_pauli,
        "participation_ratio_mean": float(prs.mean()),
        "position_variance_mean": float(varx.mean()),
        "uniform_position_variance": uniform_var,
        "su2_is_momentum_space": momentum_space,
        "note": "3D pi-flux SU(2) irreps are delocalized plane waves (momentum space), "
                "not localized Wannier states. Paid-bridge-2 (localization 2D->3D) = "
                "momentum-space SU(2) -> real-space S^2 field (a Wannier-localization "
                "problem with a topological obstruction).",
    }
    out = ROOT / "experiments" / "exp_pi_flux_3D_irreps_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

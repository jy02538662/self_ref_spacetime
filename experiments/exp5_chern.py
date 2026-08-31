"""Exp5: 2D Dirac + magnetic field -> zero modes = Chern number (= total flux).

2D lattice, 2 Dirac components per site. Dirac kinetic (antisymmetric Hermitian)
plus Peierls phases for a uniform perpendicular magnetic field.
Atiyah-Singer (2D): # zero modes = Chern number = total flux (in units of Phi0 = 2pi).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def build_D_2d(Nx, Ny, Phi, m0=0.0):
    """2D Dirac operator, dim = 2*Nx*Ny. Basis: |x,y,up> = 2*idx, |x,y,down> = 2*idx+1.

    Kinetic: (i/2) sigma_x (x-hopping) + (i/2) sigma_y (y-hopping with Peierls phase).
    Landau gauge A_x=0, A_y = 2*pi*Phi*x.  Phi = flux per plaquette (units of 2pi).
    """
    N = Nx * Ny
    dim = 2 * N
    D = np.zeros((dim, dim), dtype=np.complex128)

    def idx(x, y):
        return y * Nx + x

    # on-site mass m sigma_z (optional)
    if m0 != 0:
        for y in range(Ny):
            for x in range(Nx):
                i = idx(x, y)
                D[2 * i, 2 * i] = m0
                D[2 * i + 1, 2 * i + 1] = -m0

    # x-direction kinetic: (i/2) sigma_x (|x><x+1| - |x+1><x|)
    for y in range(Ny):
        for x in range(Nx - 1):
            i = idx(x, y)
            j = idx(x + 1, y)
            D[2 * i, 2 * j + 1] = 1j / 2
            D[2 * i + 1, 2 * j] = 1j / 2
            D[2 * j + 1, 2 * i] = -1j / 2
            D[2 * j, 2 * i + 1] = -1j / 2

    # y-direction kinetic with Peierls phase p = exp(i 2 pi Phi x)
    for x in range(Nx):
        for y in range(Ny - 1):
            i = idx(x, y)
            j = idx(x, y + 1)
            p = np.exp(1j * 2 * np.pi * Phi * x)
            pc = np.conj(p)
            D[2 * i, 2 * j + 1] = 0.5 * p
            D[2 * i + 1, 2 * j] = -0.5 * p
            D[2 * j + 1, 2 * i] = -0.5 * pc
            D[2 * j, 2 * i + 1] = 0.5 * pc

    return D


def analyze(D, Nx, Ny, Phi, tol=1e-8):
    ev = np.linalg.eigvalsh(D)
    zero = ev[np.abs(ev) < tol]
    n_zero = int(zero.size)
    n_plaquettes = (Nx - 1) * (Ny - 1)
    chern = n_plaquettes * Phi  # total flux = Chern number (uniform field)
    return {
        "Nx": Nx,
        "Ny": Ny,
        "Phi": Phi,
        "n_plaquettes": n_plaquettes,
        "chern_flux": chern,
        "eig_min_abs": float(np.min(np.abs(ev))),
        "n_zero": n_zero,
        "n_pos": int(np.sum(ev > tol)),
        "n_neg": int(np.sum(ev < -tol)),
    }


def main(Nx=20, Ny=20, p=1, m0=0.0, tol=1e-8):
    # p magnetic-flux quanta total: Phi = p / n_plaquettes
    n_plaq = (Nx - 1) * (Ny - 1)
    Phi = p / n_plaq
    print(f"=== Exp5 2D Dirac + B  Nx={Nx} Ny={Ny}  flux quanta p={p}  Phi={Phi:.4f} ===")
    rows = [
        analyze(build_D_2d(Nx, Ny, 0.0, m0=m0), Nx, Ny, 0.0, tol=tol),
        analyze(build_D_2d(Nx, Ny, Phi, m0=m0), Nx, Ny, Phi, tol=tol),
    ]
    for r in rows:
        print(
            f"Phi={r['Phi']:.4f}  chern_flux={r['chern_flux']:.3f}  "
            f"zero_modes={r['n_zero']:2d}  eig_min_abs={r['eig_min_abs']:.2e}"
        )
    chern = int(round(rows[1]["chern_flux"]))
    doubling = 2  # naive lattice Dirac fermion doubling / chiral pairing
    summary = {
        "no_field_zero": rows[0]["n_zero"],
        "field_zero": rows[1]["n_zero"],
        "chern_flux": rows[1]["chern_flux"],
        "expected_zero": doubling * chern,
        "pass": rows[0]["n_zero"] == 0 and rows[1]["n_zero"] == doubling * chern,
        "criterion": f"B=0 -> 0 zero modes; B!=0 -> zero modes = {doubling} x Chern number (= total flux)",
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    out = ROOT / "experiments" / "exp5_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp5 2D Chern number / zero modes")
    p.add_argument("--nx", type=int, default=20)
    p.add_argument("--ny", type=int, default=20)
    p.add_argument("--p", type=int, default=1)
    p.add_argument("--m0", type=float, default=0.0)
    p.add_argument("--tol", type=float, default=1e-8)
    args = p.parse_args()
    main(Nx=args.nx, Ny=args.ny, p=args.p, m0=args.m0, tol=args.tol)

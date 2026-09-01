"""Closed-loop (plaquette) magnetic flux — the first phase-sensitive term.

Phase = gauge connection.  Flux around a triangle (i,j,k) is
    Phi_ijk = phi_ij + phi_jk - phi_ik   (Hermitian => phi_ki = -phi_ik).

sin(Phi) is computed from complex products, NOT arg(), so there is no branch
cut and the term stays smooth for autodiff / L-BFGS.

Important: sin(Phi_ijk) depends ONLY on phases (z/|z| = e^{i phi}), so the
flux term decouples from the moduli.
"""

from __future__ import annotations

import numpy as np

EPS = 1e-12


def triangle_unit_phase(D: np.ndarray, i: int, j: int, k: int) -> complex:
    """e^{i Phi_ijk} = (z_ij * z_jk * conj(z_ik)) / (|z_ij| |z_jk| |z_ik|)."""
    denom = abs(D[i, j]) * abs(D[j, k]) * abs(D[i, k])
    if denom < EPS:
        return 1.0 + 0j
    return (D[i, j] * D[j, k] * np.conjugate(D[i, k])) / denom


def flux_term(D: np.ndarray) -> float:
    """S_flux = sum_{i<j<k} sin^2(Phi_ijk).  (Z2-permissive: Phi=pi allowed.)"""
    return plaquette_flux_term(D, kind="sin")


def plaquette_flux_term(D: np.ndarray, kind: str = "sin") -> float:
    """Sum over triangles of F^2, where the plaquette curvature F is defined by `kind`.

    - "sin":    F = sin(Phi)          -> zero at Phi=0 AND Phi=pi  (allows Z2 pi-flux)
    - "wilson": F = |e^{i Phi} - 1|   -> zero only at Phi=0        (forbids pi-flux)
    - "cos":    F = cos(Phi)          -> zero at Phi=+/-pi/2       (prefers maximal flux)
    """
    n = D.shape[0]
    total = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                u = triangle_unit_phase(D, i, j, k)
                if kind == "sin":
                    F = u.imag
                elif kind == "wilson":
                    F = abs(u - 1.0)
                elif kind == "cos":
                    F = u.real
                else:
                    raise ValueError(f"unknown kind {kind!r}")
                total += F * F
    return total


def flux_stats(D: np.ndarray) -> dict:
    """Distribution of triangle fluxes (unit phases e^{i Phi})."""
    n = D.shape[0]
    sin2_list: list[float] = []
    cos_list: list[float] = []
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                u = triangle_unit_phase(D, i, j, k)
                sin2_list.append(u.imag**2)
                cos_list.append(u.real)
    sin2 = np.asarray(sin2_list)
    cos = np.asarray(cos_list)
    return {
        "n_triangles": int(len(sin2)),
        "mean_sin2": float(np.mean(sin2)) if sin2.size else 0.0,
        "max_abs_sin": float(np.max(np.abs(np.sqrt(sin2)))) if sin2.size else 0.0,
        "mean_cos": float(np.mean(cos)) if cos.size else 0.0,
        "mean_abs_cos": float(np.mean(np.abs(cos))) if cos.size else 0.0,
        # fraction of triangles with e^{i Phi} ~ -1 (pi-flux)
        "frac_neg_cos": float(np.mean(cos < -0.5)) if cos.size else 0.0,
        # fraction with |cos Phi| < 0.5, i.e. Phi near +/-pi/2 (quarter flux)
        "frac_quarter_flux": float(np.mean(np.abs(cos) < 0.5)) if cos.size else 0.0,
    }

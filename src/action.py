"""Self-referential action S[D]."""

from __future__ import annotations

import numpy as np

from .algebra import build_D, unpack_params
from .distance import distance_matrix


def tr_D2(D: np.ndarray) -> float:
    # Tr(D^2) = sum_{i,j} |D_ij|^2 for Hermitian D (incl. diag 0)
    return float(np.real(np.trace(D @ D)))


def tr_D4(D: np.ndarray) -> float:
    D2 = D @ D
    return float(np.real(np.trace(D2 @ D2)))


def geometry_penalty(d: np.ndarray, L: np.ndarray) -> float:
    # Upper triangle only, exclude diagonal
    n = d.shape[0]
    iu = np.triu_indices(n, k=1)
    return float(np.sum((d[iu] - L[iu]) ** 2))


def action_from_D(
    D: np.ndarray,
    L: np.ndarray,
    lam: float = 10.0,
    lam4: float = 0.0,
) -> dict[str, float]:
    cost = tr_D2(D)
    # For Hermitian zero-diag, Tr(D^2) = 2 sum_{i<j} |z|^2; keep raw Tr for theory match
    geo = geometry_penalty(distance_matrix(D), L)
    s4 = tr_D4(D)
    total = cost + lam * geo + lam4 * s4
    return {"S": total, "TrD2": cost, "geo": geo, "TrD4": s4}


def action_from_theta(
    theta: np.ndarray,
    n: int,
    L: np.ndarray,
    lam: float = 10.0,
    lam4: float = 0.0,
) -> float:
    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    return action_from_D(D, L, lam=lam, lam4=lam4)["S"]

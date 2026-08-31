"""Self-referential distances from coupling moduli via shortest paths."""

from __future__ import annotations

import numpy as np
from scipy.sparse.csgraph import shortest_path


EPS = 1e-12


def edge_lengths_from_D(D: np.ndarray, floor: float = EPS) -> np.ndarray:
    """Length matrix: 1/|z_ij| on edges with |z|>floor, else inf (no edge)."""
    abs_z = np.abs(D)
    n = D.shape[0]
    L = np.full((n, n), np.inf, dtype=np.float64)
    np.fill_diagonal(L, 0.0)
    mask = abs_z > floor
    L[mask] = 1.0 / abs_z[mask]
    return L


def distance_matrix(D: np.ndarray, floor: float = EPS) -> np.ndarray:
    """Graph shortest-path distances d_ij[D]. Disconnected pairs -> large penalty value."""
    lengths = edge_lengths_from_D(D, floor=floor)
    d = shortest_path(lengths, directed=False, method="D")
    # Replace unreachable with a large finite value so the action stays smooth-ish
    finite = np.isfinite(d)
    if not finite.all():
        max_finite = float(np.max(d[finite])) if finite.any() else 1.0
        d = np.where(finite, d, max_finite * 10.0 + 100.0)
    return d


def ring_target_L(n: int, a: float = 1.0) -> np.ndarray:
    """Target distances on a cycle graph with unit spacing a."""
    L = np.zeros((n, n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            L[i, j] = a * min(abs(i - j), n - abs(i - j))
    return L

"""Hermitian coupling operator D parameterized by upper-triangle moduli and phases."""

from __future__ import annotations

import numpy as np


def n_params(n: int) -> int:
    """Number of real parameters for off-diagonal Hermitian D (diag fixed to 0)."""
    return n * (n - 1)  # (n choose 2) moduli + (n choose 2) phases


def pack_params(moduli: np.ndarray, phases: np.ndarray) -> np.ndarray:
    return np.concatenate([moduli.ravel(), phases.ravel()])


def unpack_params(theta: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    m = n * (n - 1) // 2
    moduli = np.abs(theta[:m])
    phases = theta[m : 2 * m]
    return moduli, phases


def build_D(moduli: np.ndarray, phases: np.ndarray, n: int) -> np.ndarray:
    """Assemble Hermitian D with zero diagonal; z_ij = r * exp(i phi)."""
    D = np.zeros((n, n), dtype=np.complex128)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            z = moduli[idx] * np.exp(1j * phases[idx])
            D[i, j] = z
            D[j, i] = np.conjugate(z)
            idx += 1
    return D


def edge_moduli_matrix(D: np.ndarray) -> np.ndarray:
    return np.abs(D)


def random_params(n: int, rng: np.random.Generator, scale: float = 0.5) -> np.ndarray:
    m = n * (n - 1) // 2
    moduli = rng.uniform(0.05, scale, size=m)
    phases = rng.uniform(-np.pi, np.pi, size=m)
    return pack_params(moduli, phases)


def ring_neighbor_mask(n: int) -> np.ndarray:
    """Boolean upper-triangle mask: True on ring edges (i, i+1) and (0, n-1)."""
    m = n * (n - 1) // 2
    mask = np.zeros(m, dtype=bool)
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            if j == i + 1 or (i == 0 and j == n - 1):
                mask[idx] = True
            idx += 1
    return mask

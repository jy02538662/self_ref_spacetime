"""Diagnostics: sparsity, neighbor hit-rate, phase order parameters."""

from __future__ import annotations

import numpy as np

from .algebra import ring_neighbor_mask, unpack_params


def upper_moduli(theta: np.ndarray, n: int) -> np.ndarray:
    moduli, _ = unpack_params(theta, n)
    return moduli


def neighbor_stats(theta: np.ndarray, n: int, ratio: float = 0.2) -> dict:
    """
    Classify edges by |z| relative to mean neighbor |z|.
    An edge is 'strong' if |z| >= ratio * mean(|z|_neighbor).
    """
    moduli = upper_moduli(theta, n)
    mask = ring_neighbor_mask(n)
    neigh = moduli[mask]
    other = moduli[~mask]
    mean_n = float(np.mean(neigh)) if neigh.size else 0.0
    thr = ratio * mean_n if mean_n > 0 else 0.0
    strong = moduli >= thr
    # hit: fraction of true ring edges that are strong
    hit = float(np.mean(strong[mask])) if mask.any() else 0.0
    # false positive: fraction of non-ring edges that are strong
    fp = float(np.mean(strong[~mask])) if (~mask).any() else 0.0
    return {
        "mean_neighbor_|z|": mean_n,
        "mean_other_|z|": float(np.mean(other)) if other.size else 0.0,
        "strength_ratio_other_over_neigh": (
            float(np.mean(other) / mean_n) if mean_n > 1e-12 and other.size else float("inf")
        ),
        "neighbor_hit_rate": hit,
        "non_neighbor_strong_rate": fp,
        "n_strong_edges": int(np.sum(strong)),
        "n_ring_edges": int(np.sum(mask)),
    }


def phase_stats(theta: np.ndarray, n: int) -> dict:
    moduli, phases = unpack_params(theta, n)
    mask = ring_neighbor_mask(n)
    # sin^2(phi): ~1 for pure imaginary, ~0 for real
    s2 = np.sin(phases) ** 2
    ring_phases = phases[mask]
    # fraction of ring edges with |sin phi| > 0.7 (~ within ~45° of ±π/2)
    imag_like = float(np.mean(np.abs(np.sin(ring_phases)) > 0.7)) if mask.any() else 0.0
    return {
        "mean_sin2_phase_all": float(np.mean(s2)),
        "mean_sin2_phase_ring": float(np.mean(s2[mask])) if mask.any() else 0.0,
        "ring_imag_like_fraction": imag_like,
        "max_abs_sin_ring": float(np.max(np.abs(np.sin(ring_phases)))) if mask.any() else 0.0,
    }


def spectrum_D2_stats(D: np.ndarray) -> dict:
    """Sign structure of D^2 eigenvalues (causal-ish diagnostic)."""
    D2 = D @ D
    # D Hermitian => D2 Hermitian positive-semidefinite if D real-symmetric? 
    # Actually for Hermitian D, D^2 is PSD: eigenvalues of D^2 are λ(D)^2 >= 0.
    # So "signature of D^2" is never Lorentzian for plain Hermitian D^2!
    # We still report eig(D) and eig(D^2) for transparency.
    ev = np.linalg.eigvalsh(D)
    ev2 = np.linalg.eigvalsh(D2)
    return {
        "eig_D_min": float(np.min(ev)),
        "eig_D_max": float(np.max(ev)),
        "eig_D2_min": float(np.min(ev2)),
        "eig_D2_max": float(np.max(ev2)),
        "TrD4_over_TrD2sq": float(
            np.real(np.trace(D2 @ D2)) / (np.real(np.trace(D2)) ** 2 + 1e-15)
        ),
    }

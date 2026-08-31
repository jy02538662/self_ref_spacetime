"""Sanity: known ring D reproduces analytic cycle distances."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.algebra import build_D, pack_params, ring_neighbor_mask
from src.distance import distance_matrix, ring_target_L


def test_ring_distance_analytic(n: int = 6, a: float = 1.0) -> None:
    mask = ring_neighbor_mask(n)
    m = n * (n - 1) // 2
    moduli = np.zeros(m)
    moduli[mask] = 1.0 / a  # edge length = a
    phases = np.zeros(m)
    D = build_D(moduli, phases, n)
    d = distance_matrix(D)
    L = ring_target_L(n, a=a)
    err = np.max(np.abs(d - L))
    assert err < 1e-9, f"max |d-L|={err}"
    print(f"OK ring distance N={n} max|d-L|={err:.2e}")


def test_pack_roundtrip(n: int = 5) -> None:
    rng = np.random.default_rng(1)
    m = n * (n - 1) // 2
    moduli = rng.random(m)
    phases = rng.uniform(-np.pi, np.pi, m)
    theta = pack_params(moduli, phases)
    from src.algebra import unpack_params

    m2, p2 = unpack_params(theta, n)
    assert np.allclose(m2, moduli)
    assert np.allclose(p2, phases)
    print("OK pack/unpack roundtrip")


if __name__ == "__main__":
    test_pack_roundtrip()
    test_ring_distance_analytic(6)
    test_ring_distance_analytic(8)
    print("all tests passed")

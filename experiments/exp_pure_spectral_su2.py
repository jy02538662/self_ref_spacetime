"""Self-reference -> SU(2) naturally: each edge is a 2x2 Hermitian matrix.

Home base (vault note "break & self-reference", step 0): the minimal unit of a BREAK is NOT a
complex number (U(1) phase) but a 2x2 Hermitian matrix D = d0*I + d.sigma,
because a break creates a BINARY loop (i<->j conjugate pair), and the minimal
self-referential structure of a binary is 2x2; d.sigma = SU(2) (Pauli).

So the natural "U(1) -> SU(2)" upgrade is NOT hand-swapping complex -> 2x2, but
making each edge a 2x2 Hermitian matrix from the start (the minimal unit of the break).

S = -alpha Tr(D^2) + beta Tr(D^4),  Tr = spatial trace + internal 2x2 trace.
Cold start; see what structure (SU(2) holonomy / linking / 3D) naturally emerges.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)


def build_D(x, n):
    """x: per-edge 4 params (d0, dx, dy, dz). D_ij = d0 I + d.sigma (2x2 Hermitian)."""
    D = np.zeros((n, n, 2, 2), complex)
    num_edges = n * (n - 1) // 2
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            d0, dx, dy, dz = x[idx:idx + 4]
            M = d0 * I2 + dx * SX + dy * SY + dz * SZ
            D[i, j] = M
            D[j, i] = M.conj().T
            idx += 4
    return D


def matmul2(A, B):
    # A,B: (n,n,2,2); (AB)[i,j] = sum_k A[i,k] @ B[k,j]
    return np.einsum('ikac,kjcb->ijab', A, B)


def trace_spatial_internal(D):
    """Tr = sum over nodes of internal 2x2 trace of diagonal blocks."""
    return float(np.real(np.trace(np.trace(D, axis1=2, axis2=3))))


def action(x, n, alpha, beta):
    D = build_D(x, n)
    D2 = matmul2(D, D)
    D4 = matmul2(D2, D2)
    tr2 = trace_spatial_internal(D2)
    tr4 = trace_spatial_internal(D4)
    return -alpha * tr2 + beta * tr4


def main():
    n = 9
    alpha = 2.0
    print(f"Self-ref SU(2): each edge = d0 I + d.sigma (2x2 Hermitian), S = -{alpha} Tr(D^2) + b Tr(D^4)\n")
    for beta in [0.5, 0.2, 1.0]:
        print(f"--- beta={beta} ---")
        for seed in range(3):
            rng = np.random.default_rng(seed)
            num_edges = n * (n - 1) // 2
            x0 = rng.normal(0, 0.3, num_edges * 4)
            res = minimize(action, x0, args=(n, alpha, beta), method="L-BFGS-B",
                           options={"maxiter": 6000, "ftol": 1e-12, "gtol": 1e-9})
            D = build_D(res.x, n)
            # edge modulus (Frobenius norm of each 2x2 edge)
            r = np.sqrt(np.einsum('ijab,ijab->ij', D.conj(), D).real)
            deg = (r > 0.3).sum(axis=1)
            print(f"  seed={seed}: S={res.fun:.3f}  max|edge|={r.max():.3f}  degrees={sorted(deg.tolist())}")


if __name__ == "__main__":
    main()

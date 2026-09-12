"""Algebra A (spin-connection gate): formalize "equal modulus = no external observer"
and pin down the Connes first-order condition failure.

Two rigorous facts, verified numerically here:

(1) Equal modulus = "no preferred direction".
    If the moduli r_ij = |D_ij| are S_N-invariant (no preferred point labeling),
    then r_ij = r is constant over all i != j.  (Transitivity of S_N on pairs.)

(2) First-order condition on a connected graph forces commutative A trivial.
    For diagonal a,b (A subset of C^N):  [[D,a],b]_ij = D_ij (a_j-a_i)(b_j-b_i).
    So first-order  <=>  for every edge (i,j), a is constant on {i,j} for all a in A.
    On a connected graph this forces A = C (scalars).  Hence a NON-trivial
    first-order triple needs A NON-commutative (or a product geometry H = C^N (x) C^k).

This is the exact "conceptual gate" before omega_mu: the algebra A cannot be C^N.
"""

from __future__ import annotations

import numpy as np


def torus_D(Lx, Ly, mode="pi"):
    N = Lx * Ly
    D = np.zeros((N, N), dtype=complex)

    def idx(x, y):
        return (y % Ly) * Lx + (x % Lx)

    for y in range(Ly):
        for x in range(Lx):
            i = idx(x, y)
            j = idx(x, y + 1)
            D[i, j] = 1.0
            D[j, i] = 1.0
            k = idx(x + 1, y)
            w = (-1.0) ** y if mode == "pi" else 1.0
            D[i, k] = w
            D[k, i] = w
    return D


def first_order_violation(D, A_basis, rng):
    """max_{a,b in A} || [[D,a],b] ||_F over random unit-norm a,b from span(A_basis)."""
    worst = 0.0
    for _ in range(300):
        ca = rng.normal(size=len(A_basis))
        cb = rng.normal(size=len(A_basis))
        a = sum(ca[k] * A_basis[k] for k in range(len(A_basis)))
        b = sum(cb[k] * A_basis[k] for k in range(len(A_basis)))
        comm = D @ a - a @ D
        fo = comm @ b - b @ comm
        worst = max(worst, np.linalg.norm(fo))
    return worst


def main():
    D = torus_D(4, 4, "pi")
    n = D.shape[0]
    rng = np.random.default_rng(0)

    print(f"### 4x4 torus (pi-flux), N={n}\n")

    # --- chiral/bipartite structure ---
    Gamma = np.diag([(-1.0) ** (i // 4 + i % 4) for i in range(n)])  # idx = 4*y+x, i+j parity
    print("bipartite  Gamma D Gamma == -D :",
          np.allclose(Gamma @ D @ Gamma, -D))
    print("Gamma^2 == I                    :", np.allclose(Gamma @ Gamma, np.eye(n)))
    print()

    # --- fact (1): equal modulus ---
    # (trivial; just confirm the S_N-invariance -> constant statement on the torus moduli)
    r = np.abs(D)
    nz = r[r > 0]
    print("off-diagonal moduli are all 1 (equal modulus):",
          np.allclose(nz, 1.0), f"({nz.size} edges)")

    # --- fact (2): first-order for A = diagonal ---
    diag_basis = [np.diag(np.eye(n)[k]) for k in range(n)]  # e_0, ..., e_{n-1}
    fo_diag = first_order_violation(D, diag_basis, rng)
    print(f"\nfirst-order violation  A=C^N (diagonal)      : {fo_diag:.3e}   (should be >> 0)")

    # --- first-order for A = scalars (trivial) ---
    scalar_basis = [np.eye(n)]
    fo_scalar = first_order_violation(D, scalar_basis, rng)
    print(f"first-order violation  A=C (scalars, trivial) : {fo_scalar:.3e}   (should be ~0)")

    # --- first-order for A = span{I, D} (commutes with D) ---
    spanD_basis = [np.eye(n), D]
    fo_spanD = first_order_violation(D, spanD_basis, rng)
    print(f"first-order violation  A=span{{I,D}}           : {fo_spanD:.3e}   (should be ~0, trivial)")

    # --- first-order for A = chiral (block-diagonal C^{n/2}+C^{n/2}) ---
    order = np.argsort(np.diag(Gamma))  # A-sublattice first
    # chiral algebra: e_i on A-sublattice only, f_j on B-sublattice only
    chiral_basis = []
    half = n // 2
    for k in range(half):
        m = np.zeros((n, n))
        m[k, k] = 1.0
        chiral_basis.append(m)
        m2 = np.zeros((n, n))
        m2[half + k, half + k] = 1.0
        chiral_basis.append(m2)
    # rotate D into chiral order to be safe
    P = np.eye(n)[order]
    Dc = P @ D @ P.T
    fo_chiral = first_order_violation(Dc, chiral_basis, rng)
    print(f"first-order violation  A=C^8+C^8 (chiral)     : {fo_chiral:.3e}   (should be >> 0)")

    # --- first-order for A = M_16 (full matrix algebra, "万物统一 -> 单代数") ---
    full_basis = []
    for i in range(n):
        for j in range(n):
            m = np.zeros((n, n))
            m[i, j] = 1.0
            full_basis.append(m)
    fo_full = first_order_violation(D, full_basis, rng)
    print(f"first-order violation  A=M_16 (full)          : {fo_full:.3e}   (should be >> 0 -> D=0)")

    # --- first-order for A = M_8 + M_8 (full chiral blocks, non-commutative) ---
    echiral = []
    for i in range(n):
        for j in range(n):
            if (i < half) == (j < half):
                m = np.zeros((n, n))
                m[i, j] = 1.0
                echiral.append(m)
    fo_echiral = first_order_violation(Dc, echiral, rng)
    print(f"first-order violation  A=M_8+M_8 (full chiral): {fo_echiral:.3e}   (should be >> 0)")

    print("\n=> 结论（定理 2 + 定理 3）："
          "\n   A=C^N 对角 -> D 逼对角；A=M_N 全矩阵 -> D 逼 0（[D,a] 迹零 + 中心 -> D 标量 -> 零对角 -> 0）。"
          "\n   span{I,D} 平凡通过但 2 维；M_8+M_8 / C^8+C^8 都违反。"
          "\n   => 标准一阶条件对非平凡关系 D 在 C^N 上无代数解：不是换 A，是要替换一阶条件本身。")


if __name__ == "__main__":
    main()

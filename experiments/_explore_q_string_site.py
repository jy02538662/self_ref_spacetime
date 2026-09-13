"""Cut-3 / problem-2, q-deformed version (step A) -- CORRECTED construction.

⚠️ WARNING (2026-09-12): the "key discovery" below is WRONG.  This script's own
output shows `e e e = e: False` for ALL k, because the classical singlet forces
delta^2 = 4 in the braid relation (recoupling c = 1/4).  The q-deformed TL
(delta = quantum dimension != 2) needs the DIAGRAM basis (loop model), NOT the
spin-1/2 tensor product.  Superseded by exp_loop_model.py.

Key discovery (from exp_pairing_to_braid.py, #1) [FALSE, see above]: the TL generator is
    e = (d/2) eps eps^dag,  eps = |01> - |10>  (CLASSICAL singlet, no q-deform)
    d = -A^2 - A^{-2}  (q-dependent loop value),  A = q^{1/4}, q = e^{i pi/(k+2)}

i.e. the PAIRING part (singlet) is classical, only the LOOP VALUE d is quantized.
This means: the "string = spin-1/2 line, endpoint = S^2 direction" holds with the
CLASSICAL S^2 (Bloch sphere) even at finite k -- quantization only truncates d.
[FALSE: the braid relation e e e = e FAILS for d != 2, so the pairing cannot be
classical; the loop value AND the representation both change in the q-deformed case.]

Verify the TL relations (e^2 = d e, e_i e_{i+1} e_i = e_i, far commute) for finite k,
and confirm the endpoints are still classical spin-1/2.
"""

import numpy as np

I2 = np.eye(2)


def tl_4x4(A):
    """4x4 TL generator e = (d/2) eps eps^dag, eps = |01>-|10> (classical singlet)."""
    d = -(A ** 2) - A ** (-2)
    eps = np.array([0.0, 1.0, -1.0, 0.0], complex)  # |01> - |10>
    e = (d / 2.0) * np.outer(eps, eps.conj())
    return e, d


def embed(n, i, P4):
    """Embed a 4x4 matrix on qubits i,i+1 into the n-qubit space (16 dim for n=4)."""
    right = 2 ** i              # qubits to the LEFT of the pair (lower index)
    left = 2 ** (n - i - 2)     # qubits to the RIGHT of the pair
    return np.kron(np.kron(np.eye(right), P4), np.eye(left))


def check_TL(n, A):
    d = -(A ** 2) - A ** (-2)
    P4, _ = tl_4x4(A)
    e = [embed(n, i, P4) for i in range(n - 1)]
    ok_sq = all(np.allclose(e[i] @ e[i], d * e[i]) for i in range(n - 1))
    ok_braid = all(np.allclose(e[i] @ e[i + 1] @ e[i], e[i]) for i in range(n - 2))
    ok_far = all(np.allclose(e[i] @ e[j], e[j] @ e[i])
                 for i in range(n - 1) for j in range(i + 2, n - 1))
    return ok_sq, ok_braid, ok_far, d


def main():
    n = 4
    print("=== q-deformed TL via CLASSICAL singlet + q-dependent d ===")
    print("  e = (d/2)|s><s|, |s>=|01>-|10> (classical), d = -A^2 - A^{-2}, A=q^{1/4}")
    print()
    for k in (1, 2, 3, 4, 5, 10, 50):
        q = np.exp(1j * np.pi / (k + 2))
        A = q ** 0.25
        ok_sq, ok_braid, ok_far, d = check_TL(n, A)
        print(f"  k={k:>2}: d={d.real:+.4f}  e^2=d e:{ok_sq}  e e e=e:{ok_braid}  "
              f"far-commute:{ok_far}")

    print("\n=== classical limit k->inf: d -> -2, recovers singlet projector x(-2) ===")
    A1 = 1.0 + 0j
    P4, d1 = tl_4x4(A1)
    eps = np.array([0.0, 1.0, -1.0, 0.0], complex)
    P_cl = np.outer(eps, eps.conj())  # |s><s|
    print(f"  d(k->inf) = {d1.real:+.2f}")
    print(f"  e == (d/2)|s><s| == -|s><s| (since d=-2): {np.allclose(P4, -P_cl)}")

    print("\n=== endpoints are still CLASSICAL spin-1/2 (Bloch sphere) ===")
    print("  representation space = 2^n = tensor product of n spin-1/2,")
    print("  pairing = classical singlet |01>-|10|, so endpoint = classical S^2.")
    print("  => 'string = spin-1/2 line, endpoint = S^2 direction' holds at finite k;")
    print("     quantization only truncates the loop value d = -2cos(pi/(k+2)),")
    print("     it does NOT deform the S^2 (no quantum-sphere complication).")

    print("\n=== locality at finite k ===")
    q = np.exp(1j * np.pi / 5)  # k=3
    A = q ** 0.25
    P4, _ = tl_4x4(A)
    e1 = embed(n, 0, P4)
    # random operator on qubit 2 (index 1), which is NOT adjacent to qubits 0,1
    rng = np.random.default_rng(1)
    M = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
    M = (M + M.conj().T) / 2
    M_full = embed(n, 1, np.kron(M, I2))  # wait, embed expects 4x4; do it manually
    M_full = np.kron(np.kron(I2, M), np.eye(2 ** (n - 2)))
    comm = np.linalg.norm(e1 @ M_full - M_full @ e1)
    print(f"  ||[e_1, M_on_qubit2]|| = {comm:.2e}  (0 => e_1 acts only on qubits 0,1)")


if __name__ == "__main__":
    main()

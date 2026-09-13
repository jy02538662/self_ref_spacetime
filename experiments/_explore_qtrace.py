"""Step A continued: does the QUANTUM TRACE tr_q = tr(K .) give the q-dependent
loop value (quantum dimension) instead of the fixed 2?

Classical singlet projection has loop value 2 (ordinary trace), independent of q.
The q-deformed Temperley-Lieb needs loop value = -q - q^{-1} (or quantum dimension
q + q^{-1}).  Check: with the q-trace tr_q(A) = tr(K A), does the singlet
projection acquire the q-dependent loop value?

U_q(su_2) convention: K|up>=q|up>, K|down>=q^{-1}|down>;  Delta(E)=E(x)1 + K(x)E,
Delta(F)=F(x)K^{-1} + 1(x)F.
"""

import numpy as np


def q_singlet(q):
    """Solve Delta(E)|s>=0 and Delta(F)|s>=0 for |s> = a|ud> + b|du>."""
    # E|up>=0, E|down>=|up>;  F|up>=|down>, F|down>=0
    # Delta(E)|ud> = E|u>(x)|d> + K|u>(x)E|d> = 0 + q|u>(x)|u> = q |uu>
    # Delta(E)|du> = E|d>(x)|u> + K|d>(x)E|u> = |u>(x)|u> + 0 = |uu>
    # => Delta(E)|s> = (a q + b) |uu> = 0  =>  b = -a q
    # Delta(F)|ud> = F|u>(x)K^{-1}|d> + |u>(x)F|d> = |d>(x)q^{-1}|d> + 0 = q^{-1}|dd>
    # Delta(F)|du> = F|d>(x)K^{-1}|u> + |d>(x)F|u> = 0 + |d>(x)|d> = |dd>
    # => Delta(F)|s> = (a q^{-1} + b) |dd> = 0  =>  b = -a q^{-1}
    # consistent iff q = q^{-1}??  -> this signals a convention subtlety; test numerically.
    # We'll just test the two natural candidates numerically below.
    return None


def q_trace(A, K2):
    """Quantum trace on the 2-qubit space: tr_q(A) = tr((K x K) A)."""
    return np.trace(K2 @ A)


def main():
    q = np.exp(1j * np.pi / 5)  # k=3
    K = np.diag([q, 1 / q])
    K2 = np.kron(K, K)
    print(f"q = e^{{i pi/5}},  quantum dimension tr_q(1) = tr(K x K) = "
          f"{np.trace(K2).real:+.4f}")
    print(f"  (expected q+q^{-1} = {q + 1/q:+.4f}... wait, tr(K x K) = (tr K)^2)")
    print(f"  tr(K) = {np.trace(K).real:+.4f},  so tr_q(1_2qubit) = (q+q^{-1})^2")
    print()

    # candidates for the q-singlet |s> = a|ud> + b|du>
    print("=== q-trace of candidate singlets (does it pick up q-dependence?) ===")
    cands = {
        "|ud>-|du> (classical)": np.array([0.0, 1.0, -1.0, 0.0], complex),
        "|ud>-q|du>": np.array([0.0, 1.0, -q, 0.0], complex),
        "|ud>-q^{-1}|du>": np.array([0.0, 1.0, -1/q, 0.0], complex),
    }
    for name, s in cands.items():
        # normalize so <s|s>=1 (ordinary norm) to compare
        s = s / np.linalg.norm(s)
        Ps = np.outer(s, s.conj())
        trq = q_trace(Ps, K2)
        print(f"  {name:24s}: tr_q(|s><s|) = {trq.real:+.4f}")

    print()
    print("=== key question: is there a singlet whose q-trace = quantum dimension? ===")
    # The quantum trace on the FULL n-qubit space is tr_q = tr(K^{x n} .).
    # For a SINGLE cup-cap (e_i), the "loop" in the diagram basis closes ONE loop,
    # whose q-trace value should be tr(K) = q + q^{-1} (the quantum dimension).
    # On the 2-qubit block, the relevant trace is over ONE qubit, not both.
    # So the correct q-trace of a single e_i is tr over the "looped" qubit only.
    print("  Subtlety: the loop of a single e_i closes ONE qubit, so its q-trace")
    print("  should be tr(K) = q+q^{-1} (quantum dimension), not tr(K x K).")
    print("  => the q-deformation enters via the QUANTUM TRACE tr(K .), which is the")
    print("     quantum dimension, exactly the object from quantization (Chebyshev).")


if __name__ == "__main__":
    main()

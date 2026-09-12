"""Algebra A, candidate 1: bimodule (real) first-order condition.

Checks whether the Connes real-spectral-triple condition [[D,a], J b* J^-1] = 0
relaxes the single-sided condition for the torus D.

Three cases (see vault note 代数A section 4):
  (1) J = K (complex conjugation, J^2=+1, trivial real structure), A = C^N
      -> order-one = single-sided (no relaxation; D still forced diagonal).
  (2) J = T = J_Kramers * K (the theory's natural anti-unitary, T^2=-1), A = C^N
      -> order-zero [a, T b* T^-1] FAILS (T does not commute with the point algebra).
  (3) J = T, A = M_{N/2}(H) (quaternionic algebra = operators commuting with T)
      -> order-zero FAILS (T b* T^-1 = -b^T, needs commutativity, A is non-comm.).

Conclusion: candidate 1 is NOT a "quick relax".  The theory's natural J=T (Kramers)
is incompatible with the POINT algebra C^N; the real resolution is the product
geometry A = A_space (x) H with the charge-conjugation J (not Kramers T).
"""

from __future__ import annotations

import numpy as np


def torus_D(Lx, Ly):
    N = Lx * Ly
    D = np.zeros((N, N))

    def idx(x, y):
        return (y % Ly) * Lx + (x % Lx)

    for y in range(Ly):
        for x in range(Lx):
            i = idx(x, y)
            j = idx(x, y + 1)
            D[i, j] = 1.0
            D[j, i] = 1.0
            k = idx(x + 1, y)
            D[i, k] = (-1.0) ** y
            D[k, i] = (-1.0) ** y
    return D


def kramers_J(D):
    """Real antisymmetric orthogonal J with J^2=-I, JD=DJ (Kramers structure).

    Exists because D (real symmetric) has all-even eigenvalue multiplicities.
    """
    n = D.shape[0]
    ev, V = np.linalg.eigh(D)
    J = np.zeros((n, n))
    used = np.zeros(n, bool)
    for i in range(n):
        if used[i]:
            continue
        for j in range(i + 1, n):
            if not used[j] and abs(ev[i] - ev[j]) < 1e-9:
                a, b = V[:, i], V[:, j]
                J += np.outer(a, b) - np.outer(b, a)
                used[i] = used[j] = True
                break
        else:
            used[i] = True
    return J


def violation(fn, tries=300, rng=None):
    rng = rng or np.random.default_rng(0)
    w = 0.0
    for _ in range(tries):
        w = max(w, float(np.linalg.norm(fn(rng))))
    return w


def main():
    D = torus_D(4, 4)
    n = D.shape[0]
    J = kramers_J(D)
    print("Kramers J: J^2=-I", np.allclose(J @ J, -np.eye(n)),
          "antisym", np.allclose(J.T, -J), "JD=DJ", np.allclose(J @ D, D @ J))
    rng = np.random.default_rng(1)

    def rdiag():
        return np.diag(rng.normal(size=n) + 1j * rng.normal(size=n))

    def proj(X):
        return (X + J @ X.conj() @ J.T) / 2  # quaternionic algebra M_{N/2}(H)

    def rand_quat():
        return proj(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))

    print("D in M_{N/2}(H)?", np.allclose(proj(D), D))

    # (1) J = K (conjugation): [[D,a], conj(b)] == single-sided
    fo1 = violation(lambda r: (lambda a, b: (D @ a - a @ D) @ b.conj() - b.conj() @ (D @ a - a @ D))(rdiag(), rdiag()), rng=rng)
    print(f"\n(1) J=K  order-one [[D,a],conj(b)]  (A=C^N): {fo1:.3e}  (= single-sided)")

    # (2) J = T (Kramers), A = C^N: order-zero [a, T b* T^-1]
    def TbT(b):  # T b* T^-1 = J conj(b^*) J^T = J b^T J^T
        return J @ b.conj().T @ J.T
    oz2 = violation(lambda r: (lambda a, b: a @ TbT(b) - TbT(b) @ a)(rdiag(), rdiag()), rng=rng)
    print(f"(2) J=T  order-zero [a,Tb*T^-1] (A=C^N): {oz2:.3e}  (>>0 => FAILS)")

    # (3) J = T, A = M_{N/2}(H): order-zero
    oz3 = violation(lambda r: (lambda a, b: a @ TbT(b) - TbT(b) @ a)(rand_quat(), rand_quat()), rng=rng)
    print(f"(3) J=T  order-zero [a,Tb*T^-1] (A=M_N/2(H)): {oz3:.3e}  (>>0 => FAILS)")

    print("\n=> 候选 1 不算快修：J=K 退化成单边；J=T(Kramers) 与点代数/四元数代数都不兼容(order-zero)。"
          "\n   真解 = 标准模型有限部分（A=C(+)H(+)M3 + 电荷共轭 J，非 Kramers T）。")


if __name__ == "__main__":
    main()

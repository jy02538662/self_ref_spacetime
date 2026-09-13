"""Find the correct q-deformed singlet so that BOTH the braid relation
(e_i e_{i+1} e_i = e_i) AND the q-dependent loop value d = -q - q^{-1} hold.

From U_q(su(2)) (K|up>=q^{1/2}|up>, Delta(E)=E(x)K^{-1}+1(x)E):
  Delta(E)|up,down> = |up,up>
  Delta(E)|down,up> = q^{-1/2} |up,up>
  => singlet  |s_q> = |up,down> - q^{1/2} |down,up>   (annihilated by E)

Try several candidate singlet forms and check which one gives the braid relation
with the q-dependent loop value.
"""

import numpy as np


def check(n, svec, d):
    """e_i = (d/2) |s><s| embedded on spins (i,i+1); check braid + loop value."""
    Ps = np.outer(svec, svec.conj())
    I2 = np.eye(2)
    # n=3: e1 on (0,1), e2 on (1,2)
    e1 = np.kron(Ps, I2)
    e2 = np.kron(I2, Ps)
    # loop value of e1: e1^2 = d' e1 with d' = trace(Ps)/... let's compute
    d_actual = np.vdot(e1, e1 @ e1) / np.vdot(e1, e1)
    braid = np.allclose(e1 @ e2 @ e1, e1)
    return braid, d_actual


def main():
    print("=== candidate q-deformed singlets |s> = a|up,down> + b|down,up> ===")
    for k in (2, 3, 4):
        q = np.exp(1j * np.pi / (k + 2))
        d = -(q + 1 / q)  # loop value
        # basis order of the 4-dim block: |uu>,|ud>,|du>,|dd>
        cands = {
            "|ud> - q^{1/2}|du>": np.array([0.0, 1.0, -(q ** 0.5), 0.0], complex),
            "q^{-1/2}|ud> - q^{1/2}|du>": np.array([0.0, q ** (-0.5), -(q ** 0.5), 0.0], complex),
            "|ud> - q|du>": np.array([0.0, 1.0, -q, 0.0], complex),
            "q^{1/2}|ud> - q^{-1/2}|du>": np.array([0.0, q ** 0.5, -(q ** (-0.5)), 0.0], complex),
        }
        print(f"\n  k={k}, q=e^{{i pi/{k+2}}}, loop value d = {d.real:+.4f}")
        for name, svec in cands.items():
            braid, d_actual = check(3, svec, d)
            # note: check() computes d_actual from e1^2 = d_actual e1, independent of d
            print(f"    {name:28s}: braid={braid},  actual loop value={d_actual.real:+.4f}")

    print("\n=== conclusion ===")
    print("  The singlet that gives braid=True AND loop value = -q-q^{-1} is the correct")
    print("  q-deformed singlet.  Then e_i = (d/2)|s_q><s_q| is the q-deformed TL generator.")


if __name__ == "__main__":
    main()

"""Diagnose the correct normalization of the Temperley-Lieb generator so that
e_i e_{i+1} e_i = e_i holds (the braid relation), not just e_i^2 = d e_i."""

import numpy as np

eps = np.array([0.0, 1.0, -1.0, 0.0], complex)  # |01> - |10>
Ps = np.outer(eps, eps.conj())  # = 2 * |singlet><singlet|
I2 = np.eye(2)

# n=3 spins: P on (0,1) and on (1,2)
P01 = np.kron(Ps, I2)
P12 = np.kron(I2, Ps)

prod = P01 @ P12 @ P01
c = np.vdot(P01, prod) / np.vdot(P01, P01)

print("Ps = eps eps^dag,  trace(Ps) =", np.trace(Ps).real)
print("P01 @ P12 @ P01 = c * P01,  c =", c.real)
print()

# candidate normalizations
for name, e1, e2, d0 in [
    ("e_i = Ps (raw)", np.kron(Ps, I2), np.kron(I2, Ps), None),
    ("e_i = (d/2) Ps, d=-2cos(pi/5)", (None,), None, None),
]:
    pass

# correct normalization: e_i = (1/c) * Ps, so that e_i e_{i+1} e_i = e_i
e1 = (1.0 / c) * np.kron(Ps, I2)
e2 = (1.0 / c) * np.kron(I2, Ps)
print("with e_i = (1/c) Ps:")
print("  braid e1 e2 e1 == e1 :", np.allclose(e1 @ e2 @ e1, e1))
# what loop value d does this give?  e_i^2 = d e_i
d_prime = np.vdot(e1, e1 @ e1) / np.vdot(e1, e1)
print("  loop value d' =", d_prime.real, " (should relate to trace(Ps)/c)")

# Now relate to the q-deformed loop value d = -A^2 - A^{-2}.
# The TL generator with loop value d is e_i = d * |singlet><singlet| = (d/2) * Ps.
# For this to satisfy the braid relation, we found e_i = (1/c) Ps.
# Hence (d/2) Ps must equal (1/c) Ps, i.e. d/2 = 1/c, d = 2/c.
# But the q-deformed d = -A^2 - A^{-2} is NOT 2/c in general.
print()
print("=== reconciliation ===")
print("raw Ps: braid requires e_i = (1/c) Ps, i.e. loop value d = 2/c =", (2.0 / c).real)
print("q-deformed loop value d = -A^2 - A^{-2} = -2cos(pi/(k+2)) (for A=q^{1/4})")
print()
print("=> the CLASSICAL singlet eps=|01>-|10> gives loop value 2/c (independent of q).")
print("   To get the q-dependent loop value d=-2cos(pi/(k+2)), the singlet itself")
print("   MUST be q-deformed (eps_q), so that the braid relation + loop value both hold.")
print("   This is the key subtlety: the q-deformation is in the SINGLET, not just d.")

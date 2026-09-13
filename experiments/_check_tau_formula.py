"""Test the tau-mass formula m = (1/2pi) exp(2pi H / sqrt3) against the hard
constraint: if H is the Hopf charge (pi_3(S^2)=Z, an INTEGER), then the mass
ratio between adjacent H is FIXED to exp(2pi/sqrt3) ~ 37.6.  Compare with the
lepton ratios (m_mu/m_e ~ 206.8, m_tau/m_mu ~ 16.8).

This is the decisive test: does ANY integer (or half-integer) assignment of H
reproduce the lepton masses?
"""

import numpy as np

# the formula
def m_of_H(H):
    return (1 / (2 * np.pi)) * np.exp((2 * np.pi / np.sqrt(3)) * H)

# fixed ratio between adjacent integer H
ratio_int = np.exp(2 * np.pi / np.sqrt(3))
ratio_half = np.exp(np.pi / np.sqrt(3))  # for half-integer steps

print("=== fixed mass ratios predicted by the formula ===")
print(f"  adjacent INTEGER H: exp(2pi/sqrt3)  = {ratio_int:.4f}")
print(f"  adjacent HALF-integer H: exp(pi/sqrt3) = {ratio_half:.4f}")
print()
print("=== lepton ratios (experiment) ===")
print(f"  m_mu / m_e    = 206.77")
print(f"  m_tau / m_mu  = 16.82")
print()
print("=== does it match? ===")
print(f"  integer-H ratio {ratio_int:.1f} vs experiment (206.8, 16.8) -> "
      f"{'match' if abs(ratio_int-206.77)<1 or abs(ratio_int-16.82)<1 else 'NO MATCH'}")
print(f"  half-integer-H ratio {ratio_half:.1f} vs experiment -> "
      f"{'match' if abs(ratio_half-206.77)<1 or abs(ratio_half-16.82)<1 else 'NO MATCH'}")
print()
print("=== masses for integer H (if H = Hopf charge, must be integer) ===")
for H in range(0, 5):
    print(f"  H={H}: m = {m_of_H(H):.2f} MeV")
print("  (leptons for reference: e=0.511, mu=105.66, tau=1776.86)")
print()
print("=== conclusion ===")
print("  If H is an INTEGER (Hopf charge), the formula's mass ratios are FIXED to")
print("  37.6 / 6.13, which do NOT match the lepton ratios 206.8 / 16.8.")
print("  => this specific formula (m = exp(2pi H/sqrt3), integer H) is EXCLUDED by data.")

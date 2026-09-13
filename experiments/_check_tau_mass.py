"""Minimal numerical check of the tau-mass formula:
    m_tau = (1/2pi) * exp((2pi/sqrt(3)) * H)

Invert it: H = (sqrt(3)/(2pi)) * ln(2pi * m), for the three lepton masses.
Decisive question: are the inferred H values "structured" (integers / simple
rationals / theory constants), or just arbitrary numbers (=> the formula is a
re-parametrization with no predictive power)?
"""

import numpy as np

m_tau = 1776.86   # MeV
m_mu = 105.66
m_e = 0.511


def H_of(m):
    return (np.sqrt(3) / (2 * np.pi)) * np.log(2 * np.pi * m)


H_tau = H_of(m_tau)
H_mu = H_of(m_mu)
H_e = H_of(m_e)

print("=== invert the formula to get H for each lepton ===")
print(f"  H_tau = {H_tau:.4f}")
print(f"  H_mu  = {H_mu:.4f}")
print(f"  H_e   = {H_e:.4f}")
print()
print("=== differences (these would be the 'topological charge gaps') ===")
print(f"  H_tau - H_mu = {H_tau - H_mu:.4f}")
print(f"  H_mu  - H_e  = {H_mu - H_e:.4f}")
print(f"  H_tau - H_e  = {H_tau - H_e:.4f}")
print()
print("=== theory constants for comparison ===")
consts = {
    "1": 1.0, "2": 2.0, "3": 3.0, "5/2": 2.5, "3/2": 1.5, "1/2": 0.5,
    "1/3": 1 / 3, "2/3": 2 / 3,
    "sqrt2": np.sqrt(2), "sqrt3": np.sqrt(3), "sqrt6": np.sqrt(6),
    "phi": (1 + np.sqrt(5)) / 2, "pi/2": np.pi / 2, "pi": np.pi,
    "2cos(pi/6)": 2 * np.cos(np.pi / 6),  # = sqrt3, k=4
    "2cos(pi/5)": 2 * np.cos(np.pi / 5),  # = phi, k=3
    "2cos(pi/7)": 2 * np.cos(np.pi / 7),  # k=5
}
for name, v in consts.items():
    print(f"    {name:14s} = {v:.4f}")

print()
print("=== mass ratios (for reference, should match experiment) ===")
print(f"  m_mu/m_e  = {m_mu / m_e:.2f}  (exp ~206.77)")
print(f"  m_tau/m_mu = {m_tau / m_mu:.2f}  (exp ~16.82)")
print(f"  m_tau/m_e  = {m_tau / m_e:.2f}  (exp ~3477)")

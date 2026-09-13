"""Bridge A step 1: time-inclusion -- the consistency of "time direction J" and
"time evolution H" from the SAME D.

Precise proposition (2026-09-12): can we go from "time direction J = gamma^0" to
"time evolution H = gamma^0 (gamma.p + m)" via the action S[D]?  #11 already fixed
the two pieces separately; here we probe their CONSISTENCY from one D, and the
relation between J (direction) and the spatial part K (which is the hidden obstacle).

Part A: signature.  D = gamma^0 D_t + gamma^1 D_x (1+1, signature -+).
    D^2 = -D_t^2 + D_x^2  -> time term is NEGATIVE = "J enters S[D] via gamma^0^2=-1".
Part B: evolution.  D psi = 0  =>  i d_t psi = H psi,  H = gamma^0 gamma^1 p (massless).
    H is Hermitian (unitary evolution).
Part C: the hidden obstacle.  H = J . K with J = gamma^0 (direction), K = gamma^1 p
    (spatial).  [J, K] = [gamma^0, gamma^1 p] = 2 gamma^0 gamma^1 p != 0:
    "direction" and "space" do NOT commute -> time is not a direct product.

Code: `py -m experiments.exp_time_action`
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# 2x2 signature -+ (matching 家底 #11: gamma^0 anti-Hermitian, gamma^1 Hermitian)
G0 = np.array([[0, 1], [-1, 0]], complex)   # gamma^0,  G0^2 = -I,  G0^dag = -G0
G1 = np.array([[0, 1], [1, 0]], complex)    # gamma^1,  G1^2 = +I,  G1^dag = +G1
I2 = np.eye(2)


def main():
    print("=== bridge A step 1: time direction J vs time evolution H (consistency) ===")
    print()

    # ---- Part A: signature (J enters S[D] via gamma^0^2 = -1) ----
    print("Part A. signature: D = g0 D_t + g1 D_x  ->  D^2 = -D_t^2 + D_x^2 (time NEGATIVE):")
    Dt, Dx = 1.0, 1.0  # symbolic scalars
    D2 = (G0 * Dt + G1 * Dx) @ (G0 * Dt + G1 * Dx)
    # expected: (-Dt^2 + Dx^2) * I
    expected = (-Dt**2 + Dx**2) * I2
    sig_ok = bool(np.allclose(D2, expected))
    print(f"  g0^2 = {G0 @ G0.tolist()}  (=-I),  g1^2 = {G1 @ G1.tolist()}  (=+I)")
    print(f"  D^2 = {D2.tolist()}  vs  (-D_t^2+D_x^2)I = {expected.tolist()}  -> {sig_ok}")
    print(f"  => 'time direction' J = gamma^0 enters S[D] = Tr(D^2) via gamma^0^2 = -1")
    print()

    # ---- Part B: evolution (H from D psi = 0) ----
    print("Part B. evolution: D psi = 0  =>  i d_t psi = H psi,  H = g0 g1 p:")
    for p in (0.5, 1.0, 2.0):
        H = G0 @ G1 * p                    # H = gamma^0 gamma^1 p (massless)
        H_dag = H.conj().T
        herm = bool(np.allclose(H, H_dag))
        print(f"  p={p}: H = {np.round(H, 3).tolist()}   Hermitian: {herm}")
    print("  => H is Hermitian (unitary time evolution), H comes from the SAME D via D psi=0.")
    print()

    # ---- Part C: hidden obstacle [J, K] != 0 ----
    print("Part C. hidden obstacle: H = J.K,  [J, K] = [g0, g1 p] != 0:")
    J = G0                                  # direction
    for p in (0.5, 1.0):
        K = G1 * p                          # spatial part
        comm = J @ K - K @ J                # [J, K]
        print(f"  p={p}: [J, K] = {np.round(comm, 3).tolist()}   (nonzero -> direction/space entangled)")
    print()
    print("interpretation:")
    print("  - A/B: 'time direction' (J, signature) and 'time evolution' (H) BOTH come from")
    print("    the same D, consistently: J enters S[D] via gamma^0^2=-1 (negative time term),")
    print("    H emerges from D psi = 0 and is Hermitian.")
    print("  - C: BUT J (direction) and K (space) do NOT commute ([J,K]=2 g0 g1 p != 0).")
    print("    So 'time = direction x space' is NOT a direct product -- this is the hidden")
    print("    obstacle of 'time inclusion' (bridge A): the time-direction and spatial parts")
    print("    are entangled through the Clifford anti-commutation {g0, g1} = 0.")
    print("  - next: whether this non-commutativity is a FEATURE (spin/Dirac structure) or a")
    print("    BLOCKER (needs the observer algebra to resolve) is the next question.")

    summary = {
        "signature_time_negative": sig_ok,
        "H_hermitian": True,
        "J_K_do_not_commute": True,
        "note": "time direction (J, signature via g0^2=-1) and evolution (H, via D psi=0) are "
                "consistent from one D; but [J,K]!=0 so time != direction x space (entangled).",
    }
    out = ROOT / "experiments" / "exp_time_action_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

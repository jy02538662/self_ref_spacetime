"""Phase as discrete connection: Wilson loop / face flux, and the Maxwell term goes flat.

Key insight (CORRECT): the phase theta_ij = arg(D_ij) is a discrete U(1) connection
(D_ji = D_ij*  =>  theta_ji = -theta_ij).  The topological object is the Wilson loop
W(C) = prod e^{i theta_ij} = e^{i Phi_C}, and the face flux Phi_p = sum_{edges in p} theta
is the curvature.  This is the right answer to "topological charge from phase" (layer 3),
and connects to Exp6a (where FFT solves the emergent gauge field A; theta = integral of A).

This script:
  1. builds a 2D grid with natural square faces;
  2. assigns random phases;
  3. computes face flux Phi_p for each square;
  4. minimizes the compact U(1) Maxwell term  S_top = -sum_p cos(Phi_p)  over phases;
  5. shows Phi_p -> 0  (flat / pure gauge / trivial).

Honest note: the face flux is the ABELIAN (U(1)) curvature = Maxwell term, which prefers
FLAT (trivial).  The "3D" / Hopf / linking is the NON-abelian Chern-Simons object, a
separate higher structure -- not produced by this U(1) Maxwell term.  So this experiment
gives the U(1) flat solution, not "3D".
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from numpy.random import default_rng

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def grid_graph(Lx: int, Ly: int):
    """2D grid. vertices (i,j) with index i*Ly+j. Returns (V, edges, faces).

    edges: list of (u, v) vertex index pairs (horizontal then vertical).
    faces: list of (e_left, e_bottom, e_right, e_top) edge indices, oriented.
    """
    V = Lx * Ly
    # horizontal edges: (i,j)-(i,j+1); vertical edges: (i,j)-(i+1,j)
    h_edges = []
    for i in range(Lx):
        for j in range(Ly - 1):
            h_edges.append((i * Ly + j, i * Ly + (j + 1)))
    v_edges = []
    for i in range(Lx - 1):
        for j in range(Ly):
            v_edges.append((i * Ly + j, (i + 1) * Ly + j))
    edges = h_edges + v_edges
    Eh = len(h_edges)
    # faces (squares): (i,j) bottom-left
    faces = []
    for i in range(Lx - 1):
        for j in range(Ly - 1):
            # bottom h edge (i,j)-(i,j+1), right v edge (i,j+1)-(i+1,j+1),
            # top h edge (i+1,j)-(i+1,j+1), left v edge (i,j)-(i+1,j)
            bottom = i * (Ly - 1) + j            # h edge index
            top = (i + 1) * (Ly - 1) + j
            right = Eh + i * Ly + (j + 1)        # v edge index (i,j+1)
            left = Eh + i * Ly + j               # v edge index (i,j)
            faces.append((bottom, right, top, left))
    return V, edges, faces


def random_phases(E: int, rng):
    """Random antisymmetric phase theta per edge: theta[edge] in [-pi, pi]."""
    return rng.uniform(-np.pi, np.pi, E)


def maxwell_action(phases, faces):
    """S_top = -sum_p cos(Phi_p), and the face fluxes Phi_p.

    Phi_p = sum of phases around each square face (the U(1) face curvature).
    """
    fluxes = []
    S = 0.0
    for f in faces:
        Phi = sum(phases[e] for e in f)
        fluxes.append(Phi)
        S -= np.cos(Phi)
    return S, fluxes


def minimize_maxwell(phases, faces, lr=0.1, steps=2000):
    """Gradient descent on S_top = -sum cos(Phi_p) over edge phases."""
    phases = phases.copy()
    history = []
    for t in range(steps):
        # gradient: dS/dtheta_e = sum_{faces containing e} sin(Phi_face) * (sign)
        # for simplicity, use numerical gradient via the face structure
        grad = np.zeros_like(phases)
        for f in faces:
            Phi = sum(phases[e] for e in f)
            for e in f:
                grad[e] += np.sin(Phi)  # d(-cos Phi)/d theta_e = sin(Phi) (for + sign)
        phases -= lr * grad
        if t % 200 == 0:
            S, fluxes = maxwell_action(phases, faces)
            # max angular distance of face fluxes from 0 (mod 2*pi) -> 0 when flat
            dev = float(max(abs(np.angle(np.exp(1j * f))) for f in fluxes)) if fluxes else 0.0
            history.append((t, S, dev))
    return phases, history


if __name__ == "__main__":
    rng = default_rng(0)
    Lx, Ly = 3, 3  # 3x3 grid: 9 vertices, 12 edges, 4 faces
    V, edges, faces = grid_graph(Lx, Ly)
    E = len(edges)

    phases = random_phases(E, rng)

    S0, fluxes0 = maxwell_action(phases, faces)
    print("=== Phase as U(1) connection: Maxwell term goes flat ===")
    print(f"{Lx}x{Ly} grid: V={V}, E={E}, faces={len(faces)}")
    print(f"initial fluxes (mod 2pi): {[f'{f % (2*np.pi):.2f}' for f in fluxes0]}")
    print(f"initial S_top = {S0:.3f}   (minimum = -#faces = {-len(faces)})")
    print()

    phases_opt, history = minimize_maxwell(phases, faces)
    Sf, fluxesf = maxwell_action(phases_opt, faces)
    print("after gradient descent on S_top:")
    print(f"  final fluxes (mod 2pi): {[f'{f % (2*np.pi):.2f}' for f in fluxesf]}")
    print(f"  final S_top = {Sf:.3f}   (target minimum = {-len(faces)})")
    print()
    # flat = every face flux sits at a multiple of 2*pi (cos=1, pure gauge).
    # Use angular distance (np.angle of e^{i f}) so that 0, +/-2pi, +/-4pi, ... all count as flat.
    fluxes_flat = all(abs(np.angle(np.exp(1j * f))) < 0.05 for f in fluxesf)
    print(f"  flux -> 0 (flat, mod 2pi) = {fluxes_flat}")
    print()

    print("interpretation:")
    print("  phase-as-connection is CORRECT (theta = arg D is a U(1) connection, Wilson loop")
    print("  is the topological object).  But the face flux is the ABELIAN U(1) curvature =")
    print("  Maxwell term, which minimizes to FLAT (Phi -> 0, pure gauge) -- trivial, NOT")
    print("  topological.  '3D'/Hopf/linking is the NON-abelian Chern-Simons object, separate.")
    print()

    summary = {
        "grid": [Lx, Ly],
        "V": V,
        "E": E,
        "faces": len(faces),
        "initial_S_top": round(S0, 4),
        "final_S_top": round(Sf, 4),
        "min_target": -len(faces),
        "flux_flat": bool(fluxes_flat),
        "note": "phase=U(1) connection (correct); Maxwell term -> flat (trivial); Hopf/3D = non-abelian Chern-Simons (separate).",
    }
    out = ROOT / "experiments" / "exp_phase_connection_last_run.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"wrote {out}")

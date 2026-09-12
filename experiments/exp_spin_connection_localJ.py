"""Spin connection omega_mu: frame source (local Kramers J_p) + Z2 -> S^2 -> winding.

Follows the "framework source" attack (方向 1 continued):

1. Local subgraph Kramers: the STAR (vertex+4 neighbors) is NOT Kramers (eigenvalue 0
   multiplicity 3, odd); the PLAQUETTE (4-cycle) IS Kramers (eigenvalues +-sqrt2 each x2).
   => the minimal Kramers object is the FACE (plaquette), not the vertex.

2. Z2 -> S^2 -> winding: the pi-flux is a Z2 monopole.  Abelian winding = Chern number = 8
   (= 8 merons, each pi-flux plaquette = Chern 1/2).  But the S^2 (SU(2)) promotion is
   BLOCKED: the global Kramers J does NOT restrict to plaquettes (J_l^2 != -I).

3. KEY positive result: the LOCAL J_p computed FRESH from the local D_p (NOT by restricting
   the global J) IS valid for all 16 plaquettes (J_p^2 = -I, J_p D_p = D_p J_p).
   => the frame exists at the plaquette level.

4. But the twist is GAUGE: J_p varies with y (because the pi-flux gauge puts -1 on
   horizontal edges (-1)^y).  The gauge-invariant content = uniform holonomy -1 (pi-flux).

Conclusion: the frame (local J_p) exists and is a gauge choice; the physical content is the
pi-flux (Z2, abelian).  The non-abelian SU(2) (Hopf) still needs the S^2 order parameter
(SU(2) ~ S^3 -> S^2) = the third layer.
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


def kramers_J(Dp):
    """Kramers J (J^2=-I, JD=DJ) for a matrix with all-even eigenvalue multiplicities."""
    n = Dp.shape[0]
    ev, V = np.linalg.eigh(Dp)
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


def main():
    D = torus_D(4, 4)
    N = 16

    def idx(x, y):
        return (y % 4) * 4 + (x % 4)

    def plaquette(x, y):
        return [idx(x, y), idx(x + 1, y), idx(x + 1, y + 1), idx(x, y + 1)]

    # 1. local subgraph: star vs plaquette
    print("=== 1) minimal Kramers object: star vs plaquette ===")
    nb = [j for j in range(N) if abs(D[0, j]) > 1e-9]
    Dl = D[np.ix_([0] + nb, [0] + nb)]
    ev = np.linalg.eigvalsh(Dl)
    uniq, cnt = np.unique(np.round(ev, 8), return_counts=True)
    print(f"  star (vertex+4 neighbors): mults={list(cnt)}  all-even={all(c%2==0 for c in cnt)}")
    p = plaquette(0, 0)
    Dp = D[np.ix_(p, p)]
    ev = np.linalg.eigvalsh(Dp)
    uniq, cnt = np.unique(np.round(ev, 8), return_counts=True)
    print(f"  plaquette (4-cycle):       mults={list(cnt)}  all-even={all(c%2==0 for c in cnt)}")

    # 2. Z2 field + abelian winding (Chern)
    print("\n=== 2) Z2 field + abelian winding ===")
    flux = 0.0
    signs = []
    for y in range(4):
        for x in range(4):
            i = idx(x, y)
            h = D[i, idx(x + 1, y)] * D[idx(x + 1, y), idx(x + 1, y + 1)] * D[idx(x + 1, y + 1), idx(x, y + 1)] * D[idx(x, y + 1), i]
            signs.append(int(h))
            flux += np.angle(h)
    print(f"  16 plaquette holonomies = {set(signs)} (pi-flux); Chern number = {flux/(2*np.pi):.1f} (= 8 merons)")

    # 3. global J does NOT localize (restriction fails)
    print("\n=== 3) global J restrict -> plaquette ===")
    Jg = kramers_J(D)
    bad = 0
    for y in range(4):
        for x in range(4):
            inds = plaquette(x, y)
            Jl = Jg[np.ix_(inds, inds)]
            if not np.allclose(Jl @ Jl, -np.eye(4)):
                bad += 1
    print(f"  global J restriction invalid on {bad}/16 plaquettes (J_l^2 != -I)")

    # 4. LOCAL J_p from local D_p: valid + gauge-dependent
    print("\n=== 4) LOCAL J_p from local D_p (fresh Kramers) ===")
    ok = True
    for y in range(4):
        for x in range(4):
            inds = plaquette(x, y)
            Dp = D[np.ix_(inds, inds)]
            Jp = kramers_J(Dp)
            if not (np.allclose(Jp @ Jp, -np.eye(4)) and np.allclose(Jp @ Dp, Dp @ Jp)):
                ok = False
    print(f"  all 16 local J_p valid (J^2=-I, JD=DJ): {ok}")
    # gauge dependence: J_p varies with y (pi-flux gauge puts -1 on horizontal edges)
    J0 = kramers_J(D[np.ix_(plaquette(0, 0), plaquette(0, 0))])
    J1 = kramers_J(D[np.ix_(plaquette(0, 1), plaquette(0, 1))])
    print(f"  J_p(0,0) == J_p(0,1)? {np.allclose(J0, J1)}  (differs with y => gauge, not physical)")

    print("\n=> 框架存在（局部 J_p），但是 gauge；物理内容 = 均匀 π 磁通（Z2，abelian）。"
          "\n   非交换 SU(2)（Hopf）仍需 S^2 序参量 = 第三层。")


if __name__ == "__main__":
    main()

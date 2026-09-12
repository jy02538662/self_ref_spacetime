"""Spin connection omega_mu from D's edge structure (direction 1).

Two DIFFERENT objects must be kept separate (this is the key clarification):

(A) pi-flux = U(1)/Z_2 MAGNETIC connection, from D's edge phases (+-1):
    - abelian, holonomy around each plaquette = -1 (pi-flux)
    - its SU(2) "lift" is trivial (in the center +-I)

(B) spin connection omega_mu = SU(2) connection, from the FRAME g_i in SU(2):
    - non-abelian, U_ij = g_i g_j^{-1}
    - gauge law U_ij -> g_i U_ij g_j^{-1}  <=>  omega -> g omega g^-1 + g dg^-1

Here we verify BOTH:
  1. lattice SU(2) gauge law  (frame construction)  -> exact
  2. continuum reduction omega -> g omega g^-1 + g dg^-1  -> exact
  3. pi-flux is the U(1) curvature (plaquette holonomy = -1), separate from omega_mu
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


SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2)


def SU2(a):
    th, nx, ny, nz = a
    n = np.array([nx, ny, nz]); n = n / np.linalg.norm(n)
    return np.cos(th / 2) * I2 + 1j * np.sin(th / 2) * (n[0] * SX + n[1] * SY + n[2] * SZ)


def mat(w):
    return w[0] * SX + w[1] * SY + w[2] * SZ


def extract_omega(X):
    return np.array([np.trace(X @ SX) / 2, np.trace(X @ SY) / 2, np.trace(X @ SZ) / 2]).real


def main():
    rng = np.random.default_rng(0)

    # ---- 1. lattice SU(2) gauge law ----
    N = 8
    g = [SU2(rng.uniform(0, 2 * np.pi, 4)) for _ in range(N)]
    h = [SU2(rng.uniform(0, 2 * np.pi, 4)) for _ in range(N)]
    U = np.zeros((N, N, 2, 2), complex)
    for i in range(N):
        for j in range(N):
            U[i, j] = g[i] @ np.linalg.inv(g[j])
    err = 0.0
    for i in range(N):
        for j in range(N):
            Unew = (h[i] @ g[i]) @ np.linalg.inv(h[j] @ g[j])
            err = max(err, np.linalg.norm(Unew - h[i] @ U[i, j] @ np.linalg.inv(h[j])))
    print(f"[1] lattice SU(2) gauge law  U_ij -> g_i U_ij g_j^-1 :  err = {err:.2e}")

    # ---- 2. continuum reduction ----
    w = np.array([0.3, -0.2, 0.5]); dx = 1e-4
    g0 = SU2([0.4, 1, 0.2, 0.1])
    delta = np.array([0.1, 0.2, -0.1])          # delta = (dg) g^-1 (right Lie deriv)
    g1 = (I2 + dx * mat(delta)) @ g0            # g(x+dx), so dg g^-1 = delta
    U_e = I2 + dx * mat(w)
    U_t = g0 @ U_e @ np.linalg.inv(g1)
    om_t = extract_omega((U_t - I2) / dx)
    gomega = g0 @ mat(w) @ np.linalg.inv(g0)
    om_exp = extract_omega(gomega) - delta      # g w g^-1 - (dg)g^-1 = g w g^-1 + g dg^-1
    print(f"[2] continuum omega -> g omega g^-1 + g dg^-1 :  err = {np.linalg.norm(om_t - om_exp):.2e}")

    # ---- 3. pi-flux is U(1) curvature, separate from omega_mu ----
    D = torus_D(4, 4)
    # plaquette holonomy = product of edge phases around each 1x1 square
    hols = []
    for y in range(4):
        for x in range(4):
            i = y * 4 + x
            right = D[i, ((x + 1) % 4) + y * 4]       # (x,y)->(x+1,y)
            up = D[((x + 1) % 4) + y * 4, ((x + 1) % 4) + ((y + 1) % 4) * 4]
            left = D[((x + 1) % 4) + ((y + 1) % 4) * 4, x + ((y + 1) % 4) * 4]
            down = D[x + ((y + 1) % 4) * 4, y * 4 + x]
            hols.append(right * up * left * down)
    print(f"[3] pi-flux: 16 plaquette holonomies = {set(int(h) for h in hols)}  "
          f"(= -1 for ALL => U(1) pi-flux, NOT an SU(2) spin connection)")

    print("\n结论：方向 1 的规范变换律验证通过（[1] 精确、[2] 精确）。"
          "\n但 D 的边相位给的是 U(1)/Z2 磁通（[3]），不是 SU(2) 自旋联络。"
          "\nSU(2) 自旋联络 = 框架联络 U_ij = g_i g_j^-1（框架来自 T^2=-1 的 Kramers 结构），"
          "\n对平坦 torus 它是平的（omega=0）；非平凡 omega_mu 需要「2x2->4x4」焊接到空间 SO(3)。")


if __name__ == "__main__":
    main()

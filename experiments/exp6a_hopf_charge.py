"""Exp6a: Hopf 荷 = odd Chern-Simons 绕数（谱不变量）的数值验证。

给定一个单位矢量场 n(x): R^3 -> S^2（Hopfion），它的 Hopf 荷（Hopf 不变量）
是三维拓扑不变量 Q_H in Z = pi_3(S^2)。数值算法（标准，见磁性 skyrmion/Hopfion 文献）：

1. 算 emergent 磁场（skyrmion 密度）：
       F_ij = n . (∂_i n × ∂_j n),   B_i = (1/2) ε_ijk F_jk = n . (∂_j n × ∂_k n)
   即 B = 曲率 2-形式 F 的 Hodge 对偶。

2. 在库仑规范 ∇·A = 0 下解 ∇×A = B（用 FFT）：
       傅里叶空间：ik × A_hat = B_hat,  ik · A_hat = 0
       =>  A_hat = -i (k × B_hat) / |k|^2

3. Hopf 荷 = (1/16π²) ∫ A·B d³x  （因为 A∧F = A·B d³x，见推导）

可证伪判据：
  - 标准 Hopf 映射（Q=1） ->  Q_H ≈ +1
  - Q 阶 Hopfion（z_1 -> z_1^Q 相位缠绕） ->  Q_H ≈ Q
  - 平凡场 n = 常数（Q=0） ->  Q_H = 0

Hopfion 场用 S^3 参数化 + 立体投影构造（处处光滑，无奇点）：
  立体投影 R^3 -> S^3:
      X1=2x/(1+r²), X2=2y/(1+r²), X3=2z/(1+r²), X4=(1-r²)/(1+r²)
  Q 阶 Hopf 映射（复坐标 z1=X1+iX2, z2=X3+iX4）:
      n = (2Re(z1^Q z̄2^Q), 2Im(z1^Q z̄2^Q), |z1^Q|² - |z2^Q|²) / (|z1^Q|²+|z2^Q|²)

注意：这里直接算的是 odd Chern-Simons 绕数（Hopf 荷的积分形式）。
"把 Q_H 表成一族 1D Dirac 算子的谱流"是 Exp6b（族指标/Bott，开放），
本脚本只做 (a) odd Chern-Simons 直接算。
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def hopfion_field(X: np.ndarray, Y: np.ndarray, Z: np.ndarray, Q: int) -> np.ndarray:
    """Q 阶 Hopfion 场 n: R^3 -> S^2（S^3 立体投影 + Hopf 纤维化的 Q 次幂）。

    逆立体投影 R^3 -> S^3:
        X1=2x/(1+r²), X2=2y/(1+r²), X3=2z/(1+r²), X4=(1-r²)/(1+r²)
    复坐标 z1=X1+iX2, z2=X3+iX4（|z1|²+|z2|²=1），Hopf 映射：
        n = (2Re(z1^Q z̄2^Q), 2Im(z1^Q z̄2^Q), |z1^Q|²-|z2^Q|²) / (|z1^Q|²+|z2^Q|²)

    Q=1 是标准 Hopf 映射（Hopf 荷 = 1）。
    Q≥2 时 n = inv.stereo((z1/z2)^Q) 是在目标 S^2 上复合了度 Q 的映射，
    由复合律 H(g∘f)=(deg g)²H(f) 得 Hopf 荷 = Q²（不是 Q）。
    这恰是复合律的一个数值验证，也证明算法能分辨不同拓扑荷。
    """
    r2 = X * X + Y * Y + Z * Z
    inv = 1.0 / (1.0 + r2)
    X1 = 2.0 * X * inv
    X2 = 2.0 * Y * inv
    X3 = 2.0 * Z * inv
    X4 = (1.0 - r2) * inv
    z1 = X1 + 1j * X2
    z2 = X3 + 1j * X4
    z1Q = z1**Q
    z2Q = z2**Q
    m2 = np.abs(z1Q) ** 2 + np.abs(z2Q) ** 2
    w = z1Q * np.conj(z2Q)
    n1 = 2.0 * np.real(w) / m2
    n2 = 2.0 * np.imag(w) / m2
    n3 = (np.abs(z1Q) ** 2 - np.abs(z2Q) ** 2) / m2
    return np.stack([n1, n2, n3], axis=0)


def _spectral_grad(f: np.ndarray, k_axis: np.ndarray, axis: int) -> np.ndarray:
    """谱微分 ∂_axis f = ifft(i k_axis * fft(f))，谱精度（远超中心差分）。"""
    shape = [1, 1, 1]
    shape[axis] = f.shape[axis]
    kk = k_axis.reshape(shape)
    return np.fft.ifftn(1j * kk * np.fft.fftn(f)).real


def _cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """3 分量矢量场的叉积，输入 (3,N,N,N)。"""
    return np.stack(
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ],
        axis=0,
    )


def emergent_B(n: np.ndarray, k1d: np.ndarray) -> np.ndarray:
    """emergent 磁场 B_i = n·(∂_j n × ∂_k n)，(i,j,k) 循环。返回 (3,N,N,N)。

    k1d 是一维频率（三个方向格距相同，共用同一 fftfreq）。
    """
    dx = np.stack([_spectral_grad(n[i], k1d, 0) for i in range(3)])
    dy = np.stack([_spectral_grad(n[i], k1d, 1) for i in range(3)])
    dz = np.stack([_spectral_grad(n[i], k1d, 2) for i in range(3)])
    B1 = np.sum(n * _cross(dy, dz), axis=0)  # n·(∂_y n × ∂_z n)
    B2 = np.sum(n * _cross(dz, dx), axis=0)  # n·(∂_z n × ∂_x n)
    B3 = np.sum(n * _cross(dx, dy), axis=0)  # n·(∂_x n × ∂_y n)
    return np.stack([B1, B2, B3], axis=0)


def solve_A_coulomb(B: np.ndarray, kx: np.ndarray, ky: np.ndarray, kz: np.ndarray, k2: np.ndarray) -> np.ndarray:
    """库仑规范 ∇·A=0 下解 ∇×A=B：A_hat = -i(k×B_hat)/|k|²。"""
    k2s = k2.copy()
    k2s[k2s < 1e-12] = 1.0
    B1h = np.fft.fftn(B[0])
    B2h = np.fft.fftn(B[1])
    B3h = np.fft.fftn(B[2])
    A1h = -1j * (ky * B3h - kz * B2h) / k2s
    A2h = -1j * (kz * B1h - kx * B3h) / k2s
    A3h = -1j * (kx * B2h - ky * B1h) / k2s
    A1h[0, 0, 0] = A2h[0, 0, 0] = A3h[0, 0, 0] = 0.0  # DC 分量置零
    A1 = np.fft.ifftn(A1h).real
    A2 = np.fft.ifftn(A2h).real
    A3 = np.fft.ifftn(A3h).real
    return np.stack([A1, A2, A3], axis=0)


def hopf_charge(
    n: np.ndarray,
    h: float,
    k1d: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray,
    kz: np.ndarray,
    k2: np.ndarray,
) -> dict:
    """算 Hopf 荷 Q_H = (1/16π²)∫A·B d³x，并附带散度健康检查。"""
    B = emergent_B(n, k1d)
    A = solve_A_coulomb(B, kx, ky, kz, k2)
    dot = A[0] * B[0] + A[1] * B[1] + A[2] * B[2]
    QH = float(np.sum(dot) * h**3 / (16.0 * np.pi**2))

    # 健康检查：B 应无散（dF=0 的数值体现）
    divB = (
        _spectral_grad(B[0], k1d, 0) + _spectral_grad(B[1], k1d, 1) + _spectral_grad(B[2], k1d, 2)
    )
    divB_rms = float(np.sqrt(np.mean(divB**2)))
    B_rms = float(np.sqrt(np.mean(B[0] ** 2 + B[1] ** 2 + B[2] ** 2)))
    return {
        "QH": QH,
        "divB_rms_over_B_rms": (divB_rms / B_rms) if B_rms > 1e-15 else 0.0,
        "B_rms": B_rms,
    }


def main(N: int = 64, L: float = 5.0, Qs: tuple = (0, 1, 2)) -> None:
    h = 2.0 * L / N
    xs = np.linspace(-L, L - h, N)  # 周期盒子，去掉重复端点
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing="ij")
    k1d = 2.0 * np.pi * np.fft.fftfreq(N, d=h)
    kx, ky, kz = np.meshgrid(k1d, k1d, k1d, indexing="ij")
    k2 = kx**2 + ky**2 + kz**2

    print(f"=== Exp6a Hopf 荷 = odd Chern-Simons  N={N}  L={L} ===")
    rows = []
    for Q in Qs:
        if Q == 0:
            n = np.stack([np.zeros_like(X), np.zeros_like(X), np.ones_like(X)], axis=0)
        else:
            n = hopfion_field(X, Y, Z, Q)
        r = hopf_charge(n, h, k1d, kx, ky, kz, k2)
        rows.append({"Q_target": Q, **r})
        print(
            f"Q_target={Q}  Q_H={r['QH']:+.4f}  "
            f"divB_rms/B_rms={r['divB_rms_over_B_rms']:.2e}"
        )

    summary = {
        "criterion": "Q_H ~= 1 (Hopf map); Q=0 -> 0; Q>=2 -> Q^2 (composition law H(g.f)=(deg g)^2 H(f))",
        "Q_H_values": [r["QH"] for r in rows],
        "Q1_close_to_1": abs(rows[1]["QH"] - 1.0) < 0.05,
        "Q0_is_zero": abs(rows[0]["QH"]) < 1e-6,
        "Q2_close_to_Q2": abs(rows[2]["QH"] - 2.0**2) < 0.05,
        "pass": (
            abs(rows[1]["QH"] - 1.0) < 0.05
            and abs(rows[0]["QH"]) < 1e-6
            and abs(rows[2]["QH"] - 2.0**2) < 0.05
        ),
    }
    print("--- summary ---")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    out = ROOT / "experiments" / "exp6a_last_run.json"
    out.write_text(json.dumps({"rows": rows, "summary": summary}, indent=2), encoding="utf-8")
    print(f"wrote {out}")


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Exp6a Hopf charge via odd Chern-Simons")
    p.add_argument("--N", type=int, default=64)
    p.add_argument("--L", type=float, default=5.0)
    args = p.parse_args()
    main(N=args.N, L=args.L)

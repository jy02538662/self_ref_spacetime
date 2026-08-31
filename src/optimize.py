"""Ring moduli helpers and phase-mode optimization for Exp2."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

from .action import action_from_D, action_from_theta
from .algebra import build_D, pack_params, random_params, ring_neighbor_mask, unpack_params
from .distance import ring_target_L
from .metrics import neighbor_stats, phase_stats, spectrum_D2_stats


@dataclass
class OptResult:
    success: bool
    S: float
    theta: np.ndarray
    D: np.ndarray
    stats: dict
    message: str
    phase_mode: str


def exact_ring_theta(n: int, a: float = 1.0, phases: np.ndarray | None = None) -> np.ndarray:
    """Moduli = 1/a on ring edges, 0 elsewhere; optional phases on upper triangle."""
    m = n * (n - 1) // 2
    moduli = np.zeros(m)
    mask = ring_neighbor_mask(n)
    moduli[mask] = 1.0 / a
    if phases is None:
        phases = np.zeros(m)
    return pack_params(moduli, phases)


def minimize_action(
    n: int,
    lam: float = 20.0,
    lam4: float = 0.0,
    seed: int = 0,
    maxiter: int = 400,
    a: float = 1.0,
    phase_mode: str = "free",
    theta0: np.ndarray | None = None,
) -> OptResult:
    """
    phase_mode:
      - free: optimize moduli + phases
      - real: phases clamped to 0
      - phases_only: moduli fixed from theta0, optimize phases only
    """
    rng = np.random.default_rng(seed)
    L = ring_target_L(n, a=a)
    if theta0 is None:
        theta0 = random_params(n, rng, scale=0.8)
        if phase_mode == "real":
            moduli, _ = unpack_params(theta0, n)
            theta0 = pack_params(moduli, np.zeros_like(moduli))

    m = n * (n - 1) // 2

    if phase_mode == "real":
        moduli0, _ = unpack_params(theta0, n)

        def fun(moduli: np.ndarray) -> float:
            theta = pack_params(np.abs(moduli), np.zeros(m))
            return action_from_theta(theta, n, L, lam=lam, lam4=lam4)

        res = minimize(fun, moduli0, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-10})
        theta = pack_params(np.abs(res.x), np.zeros(m))
        message = str(res.message)
        success = bool(res.success)
    elif phase_mode == "phases_only":
        moduli0, phases0 = unpack_params(theta0, n)

        def fun(phases: np.ndarray) -> float:
            theta = pack_params(moduli0, phases)
            return action_from_theta(theta, n, L, lam=lam, lam4=lam4)

        res = minimize(fun, phases0, method="L-BFGS-B", options={"maxiter": maxiter, "ftol": 1e-10})
        theta = pack_params(moduli0, res.x)
        message = str(res.message)
        success = bool(res.success)
    else:
        def fun(theta: np.ndarray) -> float:
            return action_from_theta(theta, n, L, lam=lam, lam4=lam4)

        res = minimize(
            fun,
            theta0,
            method="L-BFGS-B",
            options={"maxiter": maxiter, "ftol": 1e-10},
        )
        theta = res.x
        message = str(res.message)
        success = bool(res.success)

    moduli, phases = unpack_params(theta, n)
    D = build_D(moduli, phases, n)
    parts = action_from_D(D, L, lam=lam, lam4=lam4)
    stats = {
        **parts,
        **neighbor_stats(theta, n),
        **phase_stats(theta, n),
        **spectrum_D2_stats(D),
    }
    return OptResult(
        success=success,
        S=float(parts["S"]),
        theta=theta,
        D=D,
        stats=stats,
        message=message,
        phase_mode=phase_mode,
    )

"""Dormand-Prince 5(4), the same stepper as cpp/include/dhj/core/ode.hpp.

atol = rtol = 1e-12. On each accepted step the infinity-norm tolerance budget
is updated

    error <- error * exp(L * h) + δ,

where δ = max_i (atol + rtol * max(|y_i|, |ŷ_i|)). That budget bounds the
embedded local-error estimate and is grown with the vector-field Lipschitz
constant. It is not a validated enclosure of the true truncation error.
"""
from __future__ import annotations

import math

import numpy as np

ATOL = 1e-12
RTOL = 1e-12

_A21 = 1.0 / 5
_A31, _A32 = 3.0 / 40, 9.0 / 40
_A41, _A42, _A43 = 44.0 / 45, -56.0 / 15, 32.0 / 9
_A51, _A52, _A53, _A54 = 19372.0 / 6561, -25360.0 / 2187, 64448.0 / 6561, -212.0 / 729
_A61, _A62, _A63, _A64, _A65 = 9017.0 / 3168, -355.0 / 33, 46732.0 / 5247, 49.0 / 176, -5103.0 / 18656
_A71, _A73, _A74, _A75, _A76 = 35.0 / 384, 500.0 / 1113, 125.0 / 192, -2187.0 / 6784, 11.0 / 84
_E1, _E3, _E4, _E5, _E6, _E7 = 71.0 / 57600, -71.0 / 16695, 71.0 / 1920, -17253.0 / 339200, 22.0 / 525, -1.0 / 40


def integrate(f, y, h_total: float, lipschitz: float, error: float) -> tuple[np.ndarray, float]:
    """Advance y over h_total. Returns the new state and the updated error budget."""
    y = np.array(y, dtype=float, copy=True)
    if h_total <= 0.0:
        return y, error
    k1 = np.asarray(f(y), dtype=float)
    sc = ATOL + RTOL * np.abs(y)
    d0 = float(np.max(np.abs(y) / sc))
    d1 = float(np.max(np.abs(k1) / sc))
    h = 1e-6 if (d0 < 1e-5 or d1 < 1e-5) else 0.01 * (d0 / d1)
    h = min(h, h_total)
    t = 0.0
    steps = 0
    while t < h_total:
        h = min(h, h_total - t)
        y2 = y + h * _A21 * k1
        k2 = np.asarray(f(y2), dtype=float)
        y3 = y + h * (_A31 * k1 + _A32 * k2)
        k3 = np.asarray(f(y3), dtype=float)
        y4 = y + h * (_A41 * k1 + _A42 * k2 + _A43 * k3)
        k4 = np.asarray(f(y4), dtype=float)
        y5 = y + h * (_A51 * k1 + _A52 * k2 + _A53 * k3 + _A54 * k4)
        k5 = np.asarray(f(y5), dtype=float)
        y6 = y + h * (_A61 * k1 + _A62 * k2 + _A63 * k3 + _A64 * k4 + _A65 * k5)
        k6 = np.asarray(f(y6), dtype=float)
        y7 = y + h * (_A71 * k1 + _A73 * k3 + _A74 * k4 + _A75 * k5 + _A76 * k6)
        k7 = np.asarray(f(y7), dtype=float)
        e = h * (_E1 * k1 + _E3 * k3 + _E4 * k4 + _E5 * k5 + _E6 * k6 + _E7 * k7)
        sc = ATOL + RTOL * np.maximum(np.abs(y), np.abs(y7))
        err = float(np.max(np.abs(e) / sc))
        steps += 1
        if err <= 1.0:
            delta = float(np.max(ATOL + RTOL * np.maximum(np.abs(y), np.abs(y7))))
            error = error * math.exp(lipschitz * h) + delta
            t += h
            y = y7
            k1 = k7
        factor = 5.0 if err == 0.0 else min(5.0, max(0.2, 0.9 * err ** -0.2))
        h *= factor
        if h < 1e-14:
            h = 1e-14
        if steps > 100000:
            raise RuntimeError("ODE stepper exceeded 100000 steps")
    return y, error

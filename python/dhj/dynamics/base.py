"""Base class for a dynamics model.

A dynamics file (see dubins.py, evasion.py, double_integrator.py) bundles
EVERYTHING the solver needs to know about one system:

  * the state box, the finite action set, the ODE right-hand side,
  * the failure function l(x) (> 0 outside the obstacle) and, for reach-avoid,
    the reward function r(x) (> 0 inside the target),
  * the Lipschitz constants L_f (of the vector field, infinity norm), L_l, L_r,
  * which state dimensions are periodic (angles) - these wrap instead of
    leaving the domain, both in the successor queries and in the plots,
  * how to draw the result (which 2-D slices, what overlays).

The solver core never looks at the concrete system; it only calls the methods
below. To add a system, subclass Dynamics in a new file and register it in
dynamics/__init__.py (two lines).

Numerical conventions are kept identical to the paper scripts: trajectories
are integrated with scipy.integrate.odeint at atol = rtol = 1e-12 evaluated
at the checkpoint times dt, 2dt, ..., tau, and periodic coordinates are
re-wrapped with arctan2(sin, cos).
"""
from abc import ABC, abstractmethod
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.integrate import odeint

ODE_ATOL = 1e-12
ODE_RTOL = 1e-12


def wrap_angle(a: float) -> float:
    """Wrap to (-pi, pi], exactly as the paper scripts do."""
    return float(np.arctan2(np.sin(a), np.cos(a)))


class Overlay:
    """A shape drawn on top of every (x, y) panel (obstacle, target, ...)."""

    def __init__(self, kind: str, center, radius: float, color: str, linestyle: str = "-",
                 linewidth: float = 2.0, when: str = "always"):
        self.kind = kind            # currently only "circle"
        self.center = tuple(float(c) for c in center)
        self.radius = float(radius)
        self.color = color
        self.linestyle = linestyle
        self.linewidth = linewidth
        self.when = when            # "always" | "reach_avoid" (hidden in avoid-only mode)


class Slice:
    """A 2-D slice of the state space for plotting.

    `fixed` maps a state dimension to the value at which it is pinned; the
    remaining two dimensions (`Dynamics.plot_dims`) span the panel.
    """

    def __init__(self, fixed: Dict[int, float], label: str):
        self.fixed = dict(fixed)
        self.label = label


class Dynamics(ABC):
    # ---- class-level declarations every subclass must provide -------------
    name: str = "base"                 # CLI name, used in results paths
    dim: int = 0                       # state dimension
    periodic_dims: Tuple[int, ...] = ()  # indices of angular (wrapping) states
    plot_dims: Tuple[int, int] = (0, 1)  # the two dims spanning each panel
    state_names: Tuple[str, ...] = ()  # for CSV headers / plot labels

    def __init__(self, dt: float, tau: float):
        ratio = tau / dt
        n_steps = int(np.round(ratio))
        if not np.isclose(ratio, n_steps, rtol=1e-9, atol=1e-9):
            raise ValueError(
                f"tau ({tau}) must be evenly divisible by dt ({dt}). "
                f"Current tau/dt = {ratio:.10f} (should be ~{n_steps}). "
                f"Try tau={dt * n_steps} or dt={tau / n_steps}")
        self.dt = dt
        self.tau = tau
        self.n_steps = n_steps
        print(f"  Time discretization: {self.n_steps} steps of dt={dt}s over tau={tau}s")
        # Subclasses set these in __init__:
        self.L_f: float = float("nan")
        self.L_l: float = float("nan")
        self.L_r: float = float("nan")

    # ---- required interface ----------------------------------------------
    @abstractmethod
    def get_state_bounds(self) -> np.ndarray:
        """(dim, 2) array of [lo, hi] per state."""

    @abstractmethod
    def get_action_space(self) -> List:
        """Finite list of actions (floats or tuples)."""

    @abstractmethod
    def ode(self, state: np.ndarray, t: float, action) -> Sequence[float]:
        """Right-hand side f(x, u), odeint signature."""

    @abstractmethod
    def failure_function(self, state: np.ndarray) -> float:
        """l(x): > 0 safe, <= 0 in the failure set."""

    @abstractmethod
    def reward_function(self, state: np.ndarray) -> float:
        """r(x): > 0 inside the target set (reach-avoid only)."""

    # ---- optional overrides ------------------------------------------------
    def overlays(self) -> List[Overlay]:
        """Shapes drawn on each panel (obstacle, target)."""
        return []

    def slices(self) -> List[Slice]:
        """Default: pin every non-plotted dimension at the middle of its box.
        Subclasses normally override to choose meaningful slices."""
        b = self.get_state_bounds()
        fixed = {d: 0.5 * (b[d, 0] + b[d, 1]) for d in range(self.dim) if d not in self.plot_dims}
        label = ", ".join(f"{self.state_name(d)}={v:.2f}" for d, v in fixed.items())
        return [Slice(fixed, label)]

    def describe(self) -> dict:
        """Parameters recorded in params.json. Subclasses extend this."""
        return {
            "name": self.name, "dim": self.dim, "dt": self.dt, "tau": self.tau,
            "n_steps": self.n_steps, "L_f": self.L_f, "L_l": self.L_l, "L_r": self.L_r,
            "periodic_dims": list(self.periodic_dims),
            "state_bounds": self.get_state_bounds().tolist(),
            "actions": [list(a) if isinstance(a, (tuple, list, np.ndarray)) else a
                        for a in self.get_action_space()],
        }

    # ---- shared helpers ----------------------------------------------------
    def get_lipschitz_constants(self) -> Tuple[float, float, float]:
        return self.L_f, self.L_l, self.L_r

    def get_state_dim(self) -> int:
        return self.dim

    def state_name(self, d: int) -> str:
        return self.state_names[d] if d < len(self.state_names) else f"x{d}"

    def normalize_state(self, state: np.ndarray) -> np.ndarray:
        """Wrap periodic coordinates in place (arctan2(sin, cos), as the paper
        scripts). Periodic dimensions must be the box [-pi, pi]."""
        for d in self.periodic_dims:
            state[d] = wrap_angle(state[d])
        return state

    def dynamics(self, state: np.ndarray, action) -> np.ndarray:
        """Single step of length tau (not used by the solver, kept for parity)."""
        sol = odeint(self.ode, state, [0, self.tau], args=(action,), atol=ODE_ATOL, rtol=ODE_RTOL)
        return self.normalize_state(sol[-1])

    def dynamics_multi_step(self, state: np.ndarray, action, duration: float, dt: float) -> List[np.ndarray]:
        """States at the checkpoint times dt, 2dt, ..., duration (periodic
        coordinates wrapped), from one odeint call over the whole horizon."""
        n_steps = int(np.ceil(duration / dt))
        t_eval = np.linspace(0, duration, n_steps + 1)
        solution = odeint(self.ode, state, t_eval, args=(action,), atol=ODE_ATOL, rtol=ODE_RTOL)
        states = []
        for s in solution[1:]:
            self.normalize_state(s)
            states.append(s.copy())
        return states

    def validate(self):
        """Sanity checks run once at start-up."""
        b = self.get_state_bounds()
        assert b.shape == (self.dim, 2), f"state bounds must be ({self.dim}, 2)"
        for d in self.periodic_dims:
            assert np.isclose(b[d, 0], -np.pi) and np.isclose(b[d, 1], np.pi), \
                f"periodic dimension {d} must span [-pi, pi]"
        assert all(0 <= d < self.dim for d in self.plot_dims)
        assert np.isfinite(self.L_f) and np.isfinite(self.L_l) and np.isfinite(self.L_r), \
            "set L_f, L_l, L_r in __init__"

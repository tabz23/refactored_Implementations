"""Autonomous van der Pol oscillator, 2-D state (x, y) with y = dx/dt:

    x' = y
    y' = mu * (1 - x^2) * y - x

No control: the action list is the singleton {0}, so the solver's max over
actions is a no-op. There is no periodic coordinate.

For mu = 1 the stable limit cycle (period about 6.67) lies in
x in [-2.01, 2.01], y in [-2.68, 2.68], and its closest point to the origin
is at distance about 1.53. The state box is [-2.5, 2.5] x [-2.9, 2.9]:
wide enough to hold the cycle (about 0.5 in x and 0.2 in y of margin) and
no larger, because L_f is set by the corners of this box.

Reach-avoid sets, chosen so the cycle is the object the solver has to follow:

  * obstacle: disc of radius 0.7 at the origin. This is inside the cycle
    (the unstable focus). The orbit itself never enters it.
  * target: disc of radius 0.4 at (2, 0), which the cycle passes through
    (its rightmost point is about (2.009, 0)).

Every trajectory except the origin approaches the cycle, and the cycle hits
the target, so the converged positive set is essentially the whole box
outside the obstacle. With a short action duration the positive set grows
along the orbit one step per value-iteration sweep, which is the picture
of the periodic cycle. Because L_f = 20.75, keep tau small (0.02 or 0.05);
at tau = 0.3 the Gronwall tube is thousands of times a cell and the cycle
is no longer visible.

Lipschitz constants on [-2.5, 2.5] x [-2.9, 2.9], infinity norm. Df has rows
(0, 1) and (-1 - 2 mu x y, mu (1 - x^2)), so

    ||Df||_inf = max( 1, |1 + 2 mu x y| + |mu (1 - x^2)| ).

On this box the maximum is 1 + 2*mu*2.5*2.9 + mu*(2.5^2 - 1) = 20.75,
attained at the corners where x y > 0. l and r are Euclidean signed distances,
hence 1-Lipschitz in the Euclidean norm and sqrt(2)-Lipschitz in the
infinity norm.
"""
from typing import List

import numpy as np

from .base import Dynamics, Overlay, Slice


class VanDerPol(Dynamics):
    name = "van_der_pol"
    dim = 2
    periodic_dims = ()
    plot_dims = (0, 1)
    state_names = ("x", "y")

    def __init__(self, dt: float, tau: float, mu: float = 1.0,
                 obstacle_radius: float = 0.7, target_radius: float = 0.4,
                 obstacle_position=(0.0, 0.0), target_position=(2.0, 0.0)):
        super().__init__(dt, tau)
        self.mu = mu
        self.state_bounds = np.array([[-2.5, 2.5], [-2.9, 2.9]])
        self.obstacle_position = np.array(obstacle_position, dtype=float)
        self.obstacle_radius = obstacle_radius
        self.target_position = np.array(target_position, dtype=float)
        self.target_radius = target_radius
        self.actions = [0.0]
        a = 0.5 * (self.state_bounds[0, 1] - self.state_bounds[0, 0])  # 3
        b = 0.5 * (self.state_bounds[1, 1] - self.state_bounds[1, 0])
        # Sharp ||Df||_inf on the symmetric box; valid for a >= sqrt(2).
        self.L_f = max(1.0, 1.0 + 2.0 * mu * a * b + mu * (a * a - 1.0))
        self.L_l = np.sqrt(2.0)
        self.L_r = np.sqrt(2.0)

    def get_state_bounds(self) -> np.ndarray:
        return self.state_bounds

    def get_action_space(self) -> List:
        return self.actions

    def ode(self, state, t, action):
        x, y = state
        return [y, self.mu * (1.0 - x * x) * y - x]

    def failure_function(self, state) -> float:
        pos = state[:2]
        return np.linalg.norm(pos - self.obstacle_position) - self.obstacle_radius

    def reward_function(self, state) -> float:
        pos = state[:2]
        return -(np.linalg.norm(pos - self.target_position) - self.target_radius)

    def overlays(self):
        return [Overlay("circle", self.obstacle_position, self.obstacle_radius, "darkblue"),
                Overlay("circle", self.target_position, self.target_radius, "orange",
                        linestyle="--", when="reach_avoid")]

    def slices(self):
        return [Slice({}, "phase plane")]

    def describe(self) -> dict:
        d = super().describe()
        d.update(mu=self.mu, obstacle_position=self.obstacle_position.tolist(),
                 obstacle_radius=self.obstacle_radius,
                 target_position=self.target_position.tolist(),
                 target_radius=self.target_radius)
        return d


class VanDerPolVelocityAvoid(VanDerPol):
    """JuliaReach / Althoff safety spec: autonomous van der Pol, unsafe set y >= 2.75.

    No target. Use with --mode avoid_nodiscount. The published check is that the
    box x in [1.25, 1.55], y in [2.35, 2.45] never reaches y = 2.75.
    l(x, y) = 2.75 - y, so L_l = 1 in the infinity norm.
    """

    name = "van_der_pol_avoid"

    def __init__(self, dt: float, tau: float, mu: float = 1.0, y_unsafe: float = 2.75):
        super().__init__(dt, tau, mu=mu)
        self.y_unsafe = y_unsafe
        self.L_l = 1.0
        self.L_r = 1.0

    def failure_function(self, state) -> float:
        return self.y_unsafe - state[1]

    def overlays(self):
        return []

    def describe(self) -> dict:
        d = super().describe()
        d["y_unsafe"] = self.y_unsafe
        d["spec"] = "juliareach_velocity_avoid"
        return d

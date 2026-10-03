"""Planar double integrator, 4-D state (x, y, vx, vy), control a = (ax, ay):

    x' = vx,  y' = vy,  vx' = ax,  vy' = ay,     a in {-a_max, 0, a_max}^2

Same reach-avoid sets as the Dubins car: obstacle disc radius 1.3 at the
origin, target disc radius 0.5 at (2.5, 0). Position box [-3, 3]^2, velocity
box [-v_max, v_max]^2. No periodic dimension. A trajectory whose reach box
leaves the velocity box is treated like one leaving the position box: the
out-of-bounds successor (value -1) is added, so the velocity bounds act as a
state constraint.

Lipschitz constants (infinity norm):
  L_f = 1       |f(x)-f(y)|_inf = max(|dvx|, |dvy|, 0, 0) <= |x - y|_inf
  L_l = L_r = sqrt(2)   l and r depend on (x, y) only, exactly as for Dubins.
"""
from itertools import product
from typing import List

import numpy as np

from .base import Dynamics, Overlay, Slice


class DoubleIntegrator4D(Dynamics):
    name = "double_integrator"
    dim = 4
    periodic_dims = ()
    plot_dims = (0, 1)
    state_names = ("x", "y", "vx", "vy")

    def __init__(self, dt: float, tau: float, a_max: float = 1.0, v_max: float = 1.0,
                 obstacle_radius: float = 1.3, target_radius: float = 0.5,
                 obstacle_position=(0.0, 0.0), target_position=(2.5, 0.0)):
        super().__init__(dt, tau)
        self.a_max = a_max
        self.v_max = v_max
        self.state_bounds = np.array([[-3.0, 3.0], [-3.0, 3.0], [-v_max, v_max], [-v_max, v_max]])
        self.obstacle_position = np.array(obstacle_position, dtype=float)
        self.obstacle_radius = obstacle_radius
        self.target_position = np.array(target_position, dtype=float)
        self.target_radius = target_radius
        levels = [-a_max, 0.0, a_max]
        self.actions = [tuple(a) for a in product(levels, levels)]  # 9 actions
        self.L_f = 1.0
        self.L_l = np.sqrt(2)
        self.L_r = np.sqrt(2)

    def get_state_bounds(self) -> np.ndarray:
        return self.state_bounds

    def get_action_space(self) -> List:
        return self.actions

    def ode(self, state, t, action):
        ax, ay = action
        return [state[2], state[3], ax, ay]

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
        v = 0.5 * self.v_max
        pins = [(0.0, 0.0), (v, 0.0), (-v, 0.0), (0.0, v), (0.0, -v), (v, v)]
        return [Slice({2: vx, 3: vy}, f"vx={vx:.2f}, vy={vy:.2f}") for vx, vy in pins]

    def describe(self) -> dict:
        d = super().describe()
        d.update(a_max=self.a_max, v_max=self.v_max,
                 obstacle_position=self.obstacle_position.tolist(), obstacle_radius=self.obstacle_radius,
                 target_position=self.target_position.tolist(), target_radius=self.target_radius)
        return d

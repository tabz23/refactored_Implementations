"""3-D evasion (pursuer-relative Dubins) with state (x1, x2, x3) and
control u in linspace(-1, 1, 5):

    x1' = -v + v cos(x3) + u x2
    x2' =  v sin(x3) - u x1
    x3' = -u

Obstacle: disc of radius 1.0 at the origin.  Target: disc of radius 0.5 at (2.5, 0).

Lipschitz constants (infinity norm): with |u| <= 1 the Jacobian rows are
bounded by v + |u| <= 1 + v, so L_f = 1 + v. L_l = L_r = sqrt(2) as for Dubins.
"""
from typing import List

import numpy as np

from .base import Dynamics, Overlay, Slice


class Evasion(Dynamics):
    name = "evasion"
    dim = 3
    periodic_dims = (2,)
    plot_dims = (0, 1)
    state_names = ("x", "y", "θ")

    def __init__(self, dt: float, tau: float, v_const: float = 1.0, obstacle_radius: float = 1.0,
                 target_radius: float = 0.5, obstacle_position=(0.0, 0.0), target_position=(2.5, 0.0)):
        super().__init__(dt, tau)
        self.v_const = v_const
        self.state_bounds = np.array([[-3.0, 3.0], [-3.0, 3.0], [-np.pi, np.pi]])
        self.obstacle_position = np.array(obstacle_position, dtype=float)
        self.obstacle_radius = obstacle_radius
        self.target_position = np.array(target_position, dtype=float)
        self.target_radius = target_radius
        self.actions = np.linspace(-1.0, 1.0, 5)
        self.L_f = 1 + self.v_const
        self.L_l = np.sqrt(2)
        self.L_r = np.sqrt(2)

    def get_state_bounds(self) -> np.ndarray:
        return self.state_bounds

    def get_action_space(self) -> List:
        return self.actions.tolist()

    def ode(self, state, t, u):
        x1, x2, x3 = state
        v = self.v_const
        return [-v + v * np.cos(x3) + u * x2, v * np.sin(x3) - u * x1, -u]

    def failure_function(self, state) -> float:
        x1, x2 = state[0], state[1]
        dx, dy = x1 - self.obstacle_position[0], x2 - self.obstacle_position[1]
        return np.sqrt(dx ** 2 + dy ** 2) - self.obstacle_radius

    def reward_function(self, state) -> float:
        pos = state[:2]
        return -(np.linalg.norm(pos - self.target_position) - self.target_radius)

    def overlays(self):
        return [Overlay("circle", self.obstacle_position, self.obstacle_radius, "darkblue"),
                Overlay("circle", self.target_position, self.target_radius, "orange",
                        linestyle="--", when="reach_avoid")]

    def slices(self):
        thetas = [0, np.pi, np.pi / 4, -np.pi / 4, np.pi / 2, -np.pi / 2]
        return [Slice({2: th}, f"θ={th:.2f} rad") for th in thetas]

    def describe(self) -> dict:
        d = super().describe()
        d.update(v_const=self.v_const, obstacle_position=self.obstacle_position.tolist(),
                 obstacle_radius=self.obstacle_radius, target_position=self.target_position.tolist(),
                 target_radius=self.target_radius)
        return d

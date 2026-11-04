"""Dubins car, state (x, y, theta), control u = turn rate in {-1, 0, 1}.

    x' = v cos(theta),  y' = v sin(theta),  theta' = u

Obstacle: disc of radius 1.3 at the origin (failure set).
Target:   disc of radius 0.5 at (2.5, 0) (reach set; ignored in avoid mode).

Lipschitz constants (infinity norm, as used by the Gronwall reach sets):
  L_f = v     |f(x)-f(y)|_inf = v * max(|cos t1 - cos t2|, |sin t1 - sin t2|) <= v |t1 - t2|
  L_l = L_r = sqrt(2)   l, r are 1-Lipschitz in the Euclidean position; a cell of
              half-width r (inf-norm) has positions at most sqrt(2) r away.
"""
from typing import List

import numpy as np

from .base import Dynamics, Overlay, Slice


class DubinsCar(Dynamics):
    name = "dubins"
    dim = 3
    periodic_dims = (2,)
    plot_dims = (0, 1)
    state_names = ("x", "y", "θ")

    def __init__(self, dt: float, tau: float, v_const: float = 1.0, obstacle_radius: float = 1.3,
                 target_radius: float = 0.5, obstacle_position=(0.0, 0.0), target_position=(2.5, 0.0)):
        super().__init__(dt, tau)
        self.v_const = v_const
        self.state_bounds = np.array([[-3.0, 3.0], [-3.0, 3.0], [-np.pi, np.pi]])
        self.obstacle_position = np.array(obstacle_position, dtype=float)
        self.obstacle_radius = obstacle_radius
        self.target_position = np.array(target_position, dtype=float)
        self.target_radius = target_radius
        self.L_f = v_const
        self.L_l = np.sqrt(2)
        self.L_r = np.sqrt(2)
        self.actions = [-1.0, 0.0, 1.0]

    def get_state_bounds(self) -> np.ndarray:
        return self.state_bounds

    def get_action_space(self) -> List:
        return self.actions

    def ode(self, state, t, action):
        x, y, theta = state
        return np.array([self.v_const * np.cos(theta), self.v_const * np.sin(theta), action])

    def failure_function(self, state) -> float:
        pos = state[:2]
        return np.linalg.norm(pos - self.obstacle_position) - self.obstacle_radius

    def reward_function(self, state) -> float:
        pos = state[:2]
        return -(np.linalg.norm(pos - self.target_position) - self.target_radius)

    # ---- plotting ----------------------------------------------------------
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

"""5-D kinematic bicycle, state (x, y, θ, v, δ), wheelbase L = 1:

    x' = v cos θ,  y' = v sin θ,  θ' = v tan δ,  v' = a,  δ' = ω

Controls: a and ω each in {-a_max, 0, a_max} and {-ω_max, 0, ω_max} (9 actions).
Same reach-avoid sets as the Dubins car: obstacle disc radius 1.3 at the
origin, target disc radius 0.5 at (2.5, 0). Position box [-3, 3]^2.
θ is periodic on [-π, π]. Speed box [-v_max, v_max], steering box
[-δ_max, δ_max]. Leaving either box is out of bounds (value -1).
δ_max = 0.6 keeps tan δ Lipschitz.

Lipschitz constants (infinity norm) on this box:
  ||Df||_∞ = max( sqrt(v_max^2 + 1), |tan δ_max| + v_max sec^2(δ_max) )
  L_l = L_r = sqrt(2)   l and r depend on (x, y) only, exactly as for Dubins.
"""
from itertools import product
from typing import List

import numpy as np

from .base import Dynamics, Overlay, Slice


class Bicycle5D(Dynamics):
    name = "bicycle"
    dim = 5
    periodic_dims = (2,)
    plot_dims = (0, 1)
    state_names = ("x", "y", "θ", "v", "δ")

    def __init__(self, dt: float, tau: float, a_max: float = 1.0, omega_max: float = 1.0,
                 v_max: float = 1.0, delta_max: float = 0.6, wheelbase: float = 1.0,
                 obstacle_radius: float = 1.3, target_radius: float = 0.5,
                 obstacle_position=(0.0, 0.0), target_position=(2.5, 0.0)):
        super().__init__(dt, tau)
        self.a_max = a_max
        self.omega_max = omega_max
        self.v_max = v_max
        self.delta_max = delta_max
        self.wheelbase = wheelbase
        self.state_bounds = np.array([
            [-3.0, 3.0], [-3.0, 3.0], [-np.pi, np.pi], [-v_max, v_max], [-delta_max, delta_max]])
        self.obstacle_position = np.array(obstacle_position, dtype=float)
        self.obstacle_radius = obstacle_radius
        self.target_position = np.array(target_position, dtype=float)
        self.target_radius = target_radius
        td = np.tan(delta_max)
        row_xy = np.sqrt(v_max * v_max + 1.0)
        row_th = (abs(td) + v_max * (1.0 + td * td)) / wheelbase
        self.L_f = float(max(row_xy, row_th))
        self.L_l = np.sqrt(2)
        self.L_r = np.sqrt(2)
        self.actions = [tuple(u) for u in product([-a_max, 0.0, a_max], [-omega_max, 0.0, omega_max])]

    def get_state_bounds(self) -> np.ndarray:
        return self.state_bounds

    def get_action_space(self) -> List:
        return self.actions

    def ode(self, state, t, action):
        a, w = action
        th, v, delta = state[2], state[3], state[4]
        return [v * np.cos(th), v * np.sin(th), (v / self.wheelbase) * np.tan(delta), a, w]

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
        # All at full speed. Stopped (v = 0) is omitted: the car can then
        # accelerate either way and the slice does not keep a tail.
        v, dmax = self.v_max, self.delta_max
        pins = [(0.0, v, 0.0),
                (np.pi, v, 0.0),
                (0.5 * np.pi, v, 0.0),
                (-0.5 * np.pi, v, 0.0),
                (0.0, v, dmax),
                (0.0, v, -dmax)]
        return [Slice({2: th, 3: v, 4: d}, f"θ={th:.2f}, v={v:.2f}, δ={d:.2f}") for th, v, d in pins]

    def describe(self) -> dict:
        d = super().describe()
        d.update(wheelbase=self.wheelbase, a_max=self.a_max, omega_max=self.omega_max,
                 v_max=self.v_max, delta_max=self.delta_max,
                 obstacle_position=self.obstacle_position.tolist(), obstacle_radius=self.obstacle_radius,
                 target_position=self.target_position.tolist(), target_radius=self.target_radius)
        return d

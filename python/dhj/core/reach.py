"""Gronwall reachability bounds at intermediate checkpoints.

A cell of half-width r (infinity norm) around its centre reaches, after time
t under a fixed action, a set contained in the box of half-width
r * exp(L_f t) around the centre trajectory. Boxes are taken at the checkpoint
times dt, 2 dt, ..., tau and every leaf meeting any of them is a successor.
"""
import numpy as np


class GronwallReachabilityAnalyzer:
    def __init__(self, dyn, use_infinity_norm: bool = True):
        self.dyn = dyn
        self.use_infinity_norm = use_infinity_norm
        self.L = dyn.L_f
        self.dt = dyn.dt
        self.tau = dyn.tau
        self.n_checkpoints = int(np.round(self.tau / self.dt))
        self.checkpoint_times = np.linspace(self.dt, self.tau, self.n_checkpoints)
        self.growth_factors = [np.exp(self.L * t) for t in self.checkpoint_times]

        print("\nGrönwall Reachability Initialized:")
        print(f"  Lipschitz constant L = {self.L}")
        print(f"  Checkpoint interval dt = {self.dt}")
        print(f"  Control horizon τ = {self.tau}")
        print(f"  Number of checkpoints = {self.n_checkpoints}")
        print(f"  Checkpoint times: {[f'{t:.3f}' for t in self.checkpoint_times]}")
        print(f"  Growth factors: {[f'{gf:.6f}' for gf in self.growth_factors]}")

    def describe(self) -> dict:
        return {"L_f": self.L, "n_checkpoints": self.n_checkpoints,
                "checkpoint_times": [float(t) for t in self.checkpoint_times],
                "growth_factors": [float(g) for g in self.growth_factors],
                "use_infinity_norm": self.use_infinity_norm}

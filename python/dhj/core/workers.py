"""multiprocessing worker functions.

The dynamics object is installed once per worker process (Pool initializer)
instead of being pickled into every task as the paper scripts did; the
arithmetic per task is unchanged.
"""
import numpy as np

_DYN = None


def init_worker(dyn):
    global _DYN  # noqa: PLW0603
    _DYN = dyn


def cell_init_worker(task):
    """(cell_id, bounds, center, L_l, L_r) -> (cell_id, l_lower, l_upper, r_lower, r_upper)."""
    cell_id, bounds, center, L_l, L_r = task
    l_center = _DYN.failure_function(center)
    r_center = _DYN.reward_function(center)
    ranges = bounds[:, 1] - bounds[:, 0]
    r = 0.5 * np.max(ranges)
    return (cell_id, l_center - L_l * r, l_center + L_l * r, r_center - L_r * r, r_center + L_r * r)


def trajectory_worker(task):
    """(cell_id, center, action_idx, action, tau, dt) -> (cell_id, action_idx, checkpoint states)."""
    cell_id, center, action_idx, action, tau, dt = task
    return (cell_id, action_idx, _DYN.dynamics_multi_step(center, action, tau, dt))

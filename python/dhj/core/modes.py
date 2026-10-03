"""Solver modes. The three paper scripts differ only in the Bellman backup:

    RA_NODISCOUNT     V <- min(l, max(r, best))            (RA_nodiscount.py)
    RA_DISCOUNT       V <- min(l, max(r, gamma * best))    (RA_discount.py)
    AVOID_NODISCOUNT  V <- min(l, best)                    (avoid_nodiscount.py)

where `best` is the max over actions of the min (lower bound) / max (upper
bound) of V over the successor cells of that action.
"""
from enum import Enum


class Mode(str, Enum):
    RA_NODISCOUNT = "ra_nodiscount"
    RA_DISCOUNT = "ra_discount"
    AVOID_NODISCOUNT = "avoid_nodiscount"

    @property
    def discounted(self) -> bool:
        return self is Mode.RA_DISCOUNT

    @property
    def has_target(self) -> bool:
        """Reach-avoid modes use the reward r(x); the avoid mode has no target."""
        return self is not Mode.AVOID_NODISCOUNT

    @property
    def script_name(self) -> str:
        """Name of the original paper script this mode reproduces."""
        return {
            Mode.RA_NODISCOUNT: "RA_nodiscount.py",
            Mode.RA_DISCOUNT: "RA_discount.py",
            Mode.AVOID_NODISCOUNT: "avoid_nodiscount.py",
        }[self]

    @property
    def default_vi_iterations_alg1(self) -> int:
        """`--iterations` default of the original script (Algorithm 1 only)."""
        return {Mode.RA_NODISCOUNT: 200, Mode.RA_DISCOUNT: 2000, Mode.AVOID_NODISCOUNT: 20000}[self]

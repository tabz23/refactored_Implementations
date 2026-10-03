"""Registry of dynamics models.

To add a new system:
  1. create dynamics/<my_system>.py with a subclass of Dynamics (see base.py),
  2. import it here and add it to REGISTRY under the name used on the CLI.
"""
from .base import Dynamics, Overlay, Slice  # noqa: F401
from .double_integrator import DoubleIntegrator4D
from .dubins import DubinsCar
from .evasion import Evasion
from .van_der_pol import VanDerPol, VanDerPolVelocityAvoid

REGISTRY = {
    DubinsCar.name: DubinsCar,
    Evasion.name: Evasion,
    DoubleIntegrator4D.name: DoubleIntegrator4D,
    VanDerPol.name: VanDerPol,
    VanDerPolVelocityAvoid.name: VanDerPolVelocityAvoid,
}


def make_dynamics(name: str, args) -> Dynamics:
    """Instantiate a registered dynamics from the parsed CLI arguments."""
    if name not in REGISTRY:
        raise SystemExit(f"unknown dynamics '{name}'; available: {', '.join(REGISTRY)}")
    cls = REGISTRY[name]
    if cls is DubinsCar:
        return DubinsCar(dt=args.dt, tau=args.tau, v_const=args.velocity, obstacle_radius=1.3)
    if cls is Evasion:
        return Evasion(dt=args.dt, tau=args.tau, v_const=args.velocity, obstacle_radius=1.0)
    if cls is DoubleIntegrator4D:
        return DoubleIntegrator4D(dt=args.dt, tau=args.tau, a_max=args.a_max, v_max=args.v_max)
    if cls is VanDerPol:
        return VanDerPol(dt=args.dt, tau=args.tau)
    if cls is VanDerPolVelocityAvoid:
        return VanDerPolVelocityAvoid(dt=args.dt, tau=args.tau)
    return cls(dt=args.dt, tau=args.tau)

# Python implementation

Reference implementation: it reproduces `Discrete_HJ_submission/*.py` **bit for
bit** (same `odeint`, same `rtree`, same order of operations) while being split
into a small package.

```
python/
├── run.py                 CLI (algorithm 1 = fixed grid VI, algorithm 2 = adaptive refinement)
├── plot_from_csv.py       render any value_function_phase_N.csv with matplotlib
├── requirements.txt       numpy scipy matplotlib rtree
└── dhj/
    ├── core/
    │   ├── modes.py           Mode enum (ra_nodiscount / ra_discount / avoid_nodiscount)
    │   ├── cell.py            Cell, CellTree (rtree index, periodic wrap-around queries, bisection)
    │   ├── reach.py           Gronwall growth factors at the checkpoint times
    │   ├── workers.py         multiprocessing worker functions (cell init, ODE trajectories)
    │   ├── vi.py              SafetyValueIterator: successor sets, vectorised Bellman sweeps, conservative stopping
    │   ├── refine.py          AdaptiveRefinement: phase loop, checkpoints, phases.csv
    │   ├── checkpoint.py      .npz checkpoints, --resume auto
    │   ├── params.py          params.json, phases.csv
    │   └── logging_utils.py   RunLogger: tee stdout/stderr -> run.log
    └── dynamics/
        ├── base.py            Dynamics ABC + Overlay/Slice (plotting hooks) + dynamics_multi_step
        ├── plotting.py        generic 3 x n_slices figure (periodic-aware slice matching)
        ├── dubins.py          Dubins car (3-D, theta periodic)
        ├── evasion.py         pursuit/evasion (3-D, theta periodic)
        ├── double_integrator.py  planar double integrator (4-D)
        └── __init__.py        REGISTRY / make_dynamics
```

## Install and run

```bash
pip install -r requirements.txt          # rtree needs libspatialindex (conda: `conda install rtree`)
python run.py --help
```

Reproducing the paper's dubins run (parameters from `results_all/.../logs359053.txt`):

```bash
python run.py --mode ra_nodiscount --dynamics dubins --algorithm 2 \
    --dt 0.3 --tau 0.3 --gamma 0.9 --delta-min 1e-8 --delta-max 1e-8 \
    --epsilon 0.05 --initial-resolution 20 --vi-iterations 100033 --conservative
```

Other modes / systems: `--mode ra_discount --gamma 0.96`, `--mode avoid_nodiscount`,
`--dynamics evasion`, `--dynamics double_integrator --a-max 1 --v-max 1`.

Useful flags: `--workers N` (default: CPUs − 1), `--no-plot`, `--dump-csv`,
`--out-dir`, `--results-root` (default: `refactored_Implementations/results`), `--tag`, `--checkpoint-every`, `--keep-checkpoints`,
`--resume auto|PATH`, `--phase0-tol` (see the top-level README), `--precompute`.

### Resume

```bash
python run.py <same arguments as before> --resume auto
```

The latest `checkpoints/phase_NNNN.npz` in the derived output directory is
loaded; mode and dynamics must match (checked), the initial resolution is taken
from the checkpoint (a warning is printed if the flag differs). The run continues
with phase NNNN+1 and appends to `run.log` / `phases.csv`.

### Plotting a CSV (e.g. produced by the C++)

```bash
python plot_from_csv.py --dynamics dubins --mode ra_nodiscount \
    /path/to/value_function_phase_5.csv --out phase5.png
```

## Adding a new dynamics

Create `dhj/dynamics/my_system.py`:

```python
import numpy as np
from .base import Dynamics, Overlay, Slice

class MySystem(Dynamics):
    name = "my_system"            # --dynamics my_system ; used in results paths
    dim = 3
    periodic_dims = (2,)          # indices of angle states (their box MUST be [-pi, pi])
    plot_dims = (0, 1)            # the two state dims spanning each figure panel
    state_names = ("x", "y", "θ")

    def __init__(self, dt, tau, speed=1.0):
        super().__init__(dt, tau)
        self.speed = speed
        self.state_bounds = np.array([[-3, 3], [-3, 3], [-np.pi, np.pi]], dtype=float)
        self.actions = [-1.0, 0.0, 1.0]       # floats or tuples, anything ode() understands
        self.L_f = speed                      # Lipschitz constant of f in the infinity norm
        self.L_l = np.sqrt(2)                 # Lipschitz constant of l w.r.t. the inf-norm on the state
        self.L_r = np.sqrt(2)                 # same for r

    def get_state_bounds(self): return self.state_bounds
    def get_action_space(self): return self.actions

    def ode(self, state, t, action):          # odeint signature f(x, t, u)
        x, y, th = state
        return [self.speed * np.cos(th), self.speed * np.sin(th), action]

    def failure_function(self, state):       # l(x) > 0 outside the obstacle
        return np.linalg.norm(state[:2]) - 1.3

    def reward_function(self, state):        # r(x) > 0 inside the target (ignored in avoid mode)
        return -(np.linalg.norm(state[:2] - np.array([2.5, 0.0])) - 0.5)

    # ---- plotting -------------------------------------------------------------
    def overlays(self):                      # shapes drawn on every panel
        return [Overlay("circle", (0, 0), 1.3, "darkblue"),
                Overlay("circle", (2.5, 0), 0.5, "orange", linestyle="--", when="reach_avoid")]

    def slices(self):                        # which 2-D slices to draw (one column each)
        return [Slice({2: th}, f"θ={th:.2f} rad") for th in (0, np.pi, np.pi / 2, -np.pi / 2)]

    def describe(self):                      # extra fields for params.json (optional)
        d = super().describe(); d.update(speed=self.speed); return d
```

Then register it in `dhj/dynamics/__init__.py`:

```python
from .my_system import MySystem
REGISTRY = {..., MySystem.name: MySystem}
# and, in make_dynamics(), construct it from the CLI args:
if name == "my_system": return MySystem(args.dt, args.tau, speed=args.velocity)
```

Notes:

* **Lipschitz constants.** `L_f` must satisfy ‖f(x,u) − f(y,u)‖_∞ ≤ L_f ‖x − y‖_∞
  for every action; it drives the Grönwall inflation r·exp(L_f t) of the reach
  boxes. `L_l`/`L_r` must satisfy |l(x) − l(y)| ≤ L_l ‖x − y‖_∞; for a signed
  distance in the (x, y) plane this is √2 (a cell of half-width r in the
  ∞-norm has points at Euclidean distance up to √2 r from its centre).
* **Periodic dims** need nothing beyond `periodic_dims`: wrapping after
  integration, ±2π index queries, exclusion from the out-of-bounds test and
  modulo-2π slice matching are done by the core and by `plotting.py`.
* **Plotting** is generic (`dhj/dynamics/plotting.py`): every leaf whose box
  contains the slice's pinned coordinates is drawn as a coloured rectangle on
  the `plot_dims` plane (upper bound, lower bound, classification), with the
  overlays on top. For a 4-D system pin two coordinates per slice
  (`Slice({2: vx, 3: vy}, label)`), see `double_integrator.py`. If you need a
  completely different figure, override `plot(self, cells, filename, iteration)`
  — `run.py` will use it if present (`make_plotter`).
* **Workers.** The dynamics object is pickled once to each worker process
  (`workers.init_worker`), so keep it picklable (no lambdas/open files).

## Implementation notes (what is identical to the paper scripts)

* grid: `np.linspace` edges, `itertools.product` order, cell ids in creation order;
* bisection along the widest dimension (ties → lowest index, 1e-10 tolerance);
* successor sets: `odeint` at the checkpoint times `linspace(0, τ, n+1)[1:]`,
  inflation `0.5·max_range·exp(L_f t_i)`, inclusive rtree intersection, θ queries
  shifted by 0/±2π, out-of-bounds sink if a box leaves a non-periodic bound;
* Bellman sweeps: vectorised with `np.minimum.reduceat` / `np.maximum.reduceat`
  over a CSR successor graph — numerically identical to the per-cell `min`/`max`
  of the scripts (no reordering of floating-point sums is involved);
* conservative stopping and correction, max-iteration fallback, boundary
  identification (`V_upper > 0 and V_lower <= 0`), refinable test
  (`max_range > ε / (2 L_l)`), phase-0 tolerance quirk.

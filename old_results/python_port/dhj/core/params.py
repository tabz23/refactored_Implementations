"""params.json / phases.csv: a machine-readable record of every run.

params.json  - all CLI arguments, the dynamics description (bounds, actions,
               Lipschitz constants, ...), derived quantities (eta_min, growth
               factors), and the environment (host, versions, start time).
phases.csv   - one row per completed refinement phase: |cells|, boundary
               cells, VI sweeps, timings. Appended as the run goes.
"""
import csv
import datetime as _dt
import json
import os
import platform
import subprocess
import sys


def _git_hash(path: str):
    try:
        out = subprocess.run(["git", "-C", path, "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None
    except Exception:  # noqa: BLE001 - best effort only
        return None


def _jsonable(x):
    try:
        json.dumps(x)
        return x
    except TypeError:
        return str(x)


def write_params(out_dir: str, args, dynamics, extra: dict, filename: str = "params.json") -> str:
    """Write params.json. Called at every (re)start; previous copies are kept as
    params_<timestamp>.json so a resumed run never loses its history."""
    import numpy  # noqa: PLC0415
    import scipy  # noqa: PLC0415

    path = os.path.join(out_dir, filename)
    if os.path.exists(path):
        stamp = _dt.datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y%m%d_%H%M%S")
        os.replace(path, os.path.join(out_dir, f"params_{stamp}.json"))

    record = {
        "command": " ".join(sys.argv),
        "started": _dt.datetime.now().isoformat(timespec="seconds"),
        "args": {k: _jsonable(v) for k, v in sorted(vars(args).items())},
        "dynamics": dynamics.describe(),
        "derived": extra,
        "environment": {
            "host": platform.node(),
            "platform": platform.platform(),
            "python": sys.version.split()[0],
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
            "git": _git_hash(os.path.dirname(os.path.abspath(__file__))),
        },
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(record, f, indent=2)
    return path


PHASE_COLUMNS = [
    "phase", "n_leaves", "boundary_cells", "refinable_cells", "refined_parents",
    "vi_iterations", "converged", "time_vi_s", "time_successors_s", "time_phase_s",
    "wall_since_start_s", "timestamp",
]


class PhaseLog:
    """Appends one CSV row per phase to <out_dir>/phases.csv."""

    def __init__(self, out_dir: str, filename: str = "phases.csv"):
        self.path = os.path.join(out_dir, filename)
        if not os.path.exists(self.path):
            with open(self.path, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(PHASE_COLUMNS)

    def append(self, **row):
        row.setdefault("timestamp", _dt.datetime.now().isoformat(timespec="seconds"))
        with open(self.path, "a", newline="", encoding="utf-8") as f:
            csv.writer(f).writerow([row.get(c, "") for c in PHASE_COLUMNS])

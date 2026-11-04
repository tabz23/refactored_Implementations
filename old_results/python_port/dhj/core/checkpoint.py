"""Checkpoint / resume for the adaptive refinement loop.

A checkpoint is written after every completed phase (phase 0 = the initial
value iteration) and holds exactly what is needed to continue: the leaf boxes
with their ids and interval values, the id counter, and the loop counters.
The successor cache is deliberately not stored: it is a pure function of the
leaf set and is rebuilt in full on the first local VI after a resume, which is
exact (and far cheaper than replaying the phases).

Format: <out_dir>/checkpoints/phase_NNNN.npz (numpy arrays + a JSON header).
"""
import glob
import json
import os
import re
from typing import List, Optional, Tuple

import numpy as np

from .cell import Cell, CellTree

FORMAT_VERSION = 1


def checkpoint_dir(out_dir: str) -> str:
    return os.path.join(out_dir, "checkpoints")


def checkpoint_path(out_dir: str, phase: int) -> str:
    return os.path.join(checkpoint_dir(out_dir), f"phase_{phase:04d}.npz")


def save_checkpoint(out_dir: str, tree: CellTree, phase: int, total_refined: int, meta: dict) -> str:
    """Atomically write the checkpoint for the state at the top of `phase`."""
    os.makedirs(checkpoint_dir(out_dir), exist_ok=True)
    leaves = tree.leaves
    n, dim = len(leaves), tree.dim
    bounds = np.empty((n, dim, 2))
    cell_id = np.empty(n, dtype=np.int64)
    parent_id = np.full(n, -1, dtype=np.int64)
    vals = np.empty((n, 6))
    for i, c in enumerate(leaves):
        bounds[i] = c.bounds
        cell_id[i] = c.cell_id
        if c.parent is not None:
            parent_id[i] = c.parent.cell_id
        vals[i] = (c.V_lower, c.V_upper, c.l_lower, c.l_upper, c.r_lower, c.r_upper)
    header = dict(meta)
    header.update(format_version=FORMAT_VERSION, phase=int(phase), total_refined=int(total_refined),
                  next_id=int(tree.next_id), dim=int(dim), n_leaves=int(n),
                  root_bounds=tree.root_bounds.tolist(), periodic_dims=list(tree.periodic_dims))
    path = checkpoint_path(out_dir, phase)
    tmp = path + ".tmp.npz"
    np.savez(tmp, header=np.array(json.dumps(header)), bounds=bounds, cell_id=cell_id,
             parent_id=parent_id, values=vals)
    os.replace(tmp, path)
    return path


def load_checkpoint(path: str) -> Tuple[CellTree, dict]:
    """Rebuild a CellTree (leaves only) and return it with the header dict."""
    with np.load(path, allow_pickle=False) as z:
        header = json.loads(str(z["header"]))
        bounds, cell_id, parent_id, vals = z["bounds"], z["cell_id"], z["parent_id"], z["values"]
    if header.get("format_version") != FORMAT_VERSION:
        raise RuntimeError(f"{path}: unsupported checkpoint format {header.get('format_version')}")
    leaves: List[Cell] = []
    parents = {}
    for i in range(len(cell_id)):
        c = Cell(bounds[i].copy(), int(cell_id[i]))
        c.V_lower, c.V_upper, c.l_lower, c.l_upper, c.r_lower, c.r_upper = (float(v) for v in vals[i])
        pid = int(parent_id[i])
        if pid >= 0:
            # a stub parent so cache invalidation keeps working for these leaves
            if pid not in parents:
                stub = Cell(bounds[i].copy(), pid)
                stub.is_leaf = False
                stub.is_refined = True
                parents[pid] = stub
            c.parent = parents[pid]
        leaves.append(c)
    tree = CellTree.from_leaves(np.array(header["root_bounds"]), leaves, int(header["next_id"]),
                                periodic_dims=tuple(header.get("periodic_dims", ())))
    return tree, header


def latest_checkpoint(out_dir: str) -> Optional[str]:
    best, best_phase = None, -1
    for p in glob.glob(os.path.join(checkpoint_dir(out_dir), "phase_*.npz")):
        m = re.fullmatch(r"phase_(\d+)\.npz", os.path.basename(p))
        if m and int(m.group(1)) > best_phase:
            best, best_phase = p, int(m.group(1))
    return best


def prune_checkpoints(out_dir: str, keep: int):
    """Keep only the newest `keep` files (0 keeps everything)."""
    if keep <= 0:
        return
    files = []
    for p in glob.glob(os.path.join(checkpoint_dir(out_dir), "phase_*.npz")):
        m = re.fullmatch(r"phase_(\d+)\.npz", os.path.basename(p))
        if m:
            files.append((int(m.group(1)), p))
    files.sort()
    for _, p in files[:-keep]:
        try:
            os.remove(p)
        except OSError:
            pass

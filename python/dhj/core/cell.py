"""Hyper-rectangular cells and the adaptively refined cell tree with an R-tree
spatial index. Dimension-agnostic port of Cell / CellTree from the paper
scripts; the only behavioural extension is that the wrap-around queries are
done for every periodic dimension (the scripts hard-coded theta = dim 2).
"""
import time
from itertools import product
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import rtree


class Cell:
    """A leaf (or refined) box of the state space with its interval values."""

    __slots__ = ("bounds", "cell_id", "center", "V_upper", "V_lower", "l_upper", "l_lower",
                 "r_upper", "r_lower", "children", "is_leaf", "is_refined", "parent")

    def __init__(self, bounds: np.ndarray, cell_id: int = 0):
        self.bounds = bounds
        self.cell_id = cell_id
        self.center = np.mean(bounds, axis=1)
        self.V_upper = None
        self.V_lower = None
        self.l_upper = None
        self.l_lower = None
        self.r_upper = None
        self.r_lower = None
        self.children = []
        self.is_leaf = True
        self.is_refined = False
        self.parent = None

    def get_range(self, dim: int) -> float:
        return self.bounds[dim, 1] - self.bounds[dim, 0]

    def get_max_range_dim(self) -> int:
        """Dimension of maximum range; ties (within 1e-10) go to the lowest index."""
        ranges = [self.get_range(j) for j in range(len(self.bounds))]
        max_range = max(ranges)
        tolerance = 1e-10
        candidates = [i for i, r in enumerate(ranges) if abs(r - max_range) < tolerance]
        return candidates[0]

    def get_max_range(self) -> float:
        return self.get_range(self.get_max_range_dim())

    def contains_point(self, point: np.ndarray) -> bool:
        for j in range(len(point)):
            if point[j] < self.bounds[j, 0] or point[j] > self.bounds[j, 1]:
                return False
        return True

    def intersects(self, other_bounds: np.ndarray) -> bool:
        for j in range(len(self.bounds)):
            if self.bounds[j, 1] < other_bounds[j, 0] or self.bounds[j, 0] > other_bounds[j, 1]:
                return False
        return True

    def split(self, next_id: int) -> Tuple["Cell", "Cell"]:
        """Bisect along the widest dimension; children get ids next_id, next_id+1."""
        dim = self.get_max_range_dim()
        mid = (self.bounds[dim, 0] + self.bounds[dim, 1]) / 2.0
        bounds1 = self.bounds.copy()
        bounds1[dim, 1] = mid
        bounds2 = self.bounds.copy()
        bounds2[dim, 0] = mid
        child1 = Cell(bounds1, next_id)
        child2 = Cell(bounds2, next_id + 1)
        child1.parent = self
        child2.parent = self
        self.children = [child1, child2]
        self.is_leaf = False
        return child1, child2


class CellTree:
    """Leaf list + R-tree over the leaves. Leaf order is creation order:
    surviving leaves keep their order, new children are appended."""

    def __init__(self, initial_bounds: np.ndarray, initial_resolution: int = 10,
                 periodic_dims: Sequence[int] = (), build: bool = True):
        self.root_bounds = np.asarray(initial_bounds, dtype=float)
        self.dim = len(self.root_bounds)
        self.periodic_dims = tuple(periodic_dims)
        self.next_id = 0
        self.leaves: List[Cell] = []
        self.spatial_index = None
        self.index_build_seconds = 0.0
        if build:
            self._create_initial_grid(initial_resolution)
            self._build_spatial_index()

    # ---- construction -------------------------------------------------------
    @classmethod
    def from_leaves(cls, root_bounds: np.ndarray, leaves: List[Cell], next_id: int,
                    periodic_dims: Sequence[int] = ()) -> "CellTree":
        """Rebuild a tree from a saved leaf list (checkpoint resume)."""
        tree = cls(root_bounds, 0, periodic_dims, build=False)
        tree.leaves = list(leaves)
        tree.next_id = next_id
        tree._build_spatial_index()
        return tree

    def _create_initial_grid(self, resolution: int):
        ranges = [np.linspace(a, b, resolution + 1) for a, b in self.root_bounds]
        for idx in product(*[range(resolution) for _ in range(self.dim)]):
            bounds = np.zeros((self.dim, 2))
            for j in range(self.dim):
                bounds[j, 0] = ranges[j][idx[j]]
                bounds[j, 1] = ranges[j][idx[j] + 1]
            self.leaves.append(Cell(bounds, self.next_id))
            self.next_id += 1

    def _build_spatial_index(self):
        print(f"  Building spatial index for {len(self.leaves)} cells...")
        t0 = time.time()
        p = rtree.index.Property()
        p.dimension = self.dim
        dim = self.dim

        def stream():
            for i, cell in enumerate(self.leaves):
                b = cell.bounds
                # interleaved=True: (min_0, ..., min_{d-1}, max_0, ..., max_{d-1})
                yield i, tuple(b[j, 0] for j in range(dim)) + tuple(b[j, 1] for j in range(dim)), None

        self.spatial_index = rtree.index.Index(stream(), properties=p)
        self.index_build_seconds = time.time() - t0
        print(f"   Spatial index built in {self.index_build_seconds:.2f}s")

    def rebuild_spatial_index(self):
        self._build_spatial_index()

    # ---- refinement ----------------------------------------------------------
    def refine_cells(self, cells: Iterable[Cell]) -> List[Cell]:
        """Split every cell in `cells` (leaf order preserved, children appended
        in refinement order). Returns the new children. Equivalent to calling
        the scripts' refine_cell in a loop, without the O(n) list.remove."""
        to_refine = [c for c in cells if c.is_leaf]
        refining = set(id(c) for c in to_refine)
        new_cells: List[Cell] = []
        for cell in to_refine:
            child1, child2 = cell.split(self.next_id)
            self.next_id += 2
            cell.is_refined = True
            new_cells.extend((child1, child2))
        kept = [c for c in self.leaves if id(c) not in refining]
        kept.extend(new_cells)
        self.leaves = kept
        return new_cells

    # ---- queries -------------------------------------------------------------
    def get_intersecting_cells(self, bounds: np.ndarray) -> List[Cell]:
        """Leaves intersecting `bounds` (inclusive of touching faces), with
        periodic dimensions compared modulo 2*pi."""
        if self.spatial_index is None:
            return [cell for cell in self.leaves if cell.intersects(bounds)]
        hits = set()
        for shifted in self._wrapped_boxes(bounds):
            hits.update(self.spatial_index.intersection(shifted))
        return [self.leaves[i] for i in hits]

    def get_intersecting_indices(self, bounds: np.ndarray) -> set:
        """Same as above but returns leaf positions (faster for the solver)."""
        hits = set()
        for shifted in self._wrapped_boxes(bounds):
            hits.update(self.spatial_index.intersection(shifted))
        return hits

    def _wrapped_boxes(self, bounds: np.ndarray):
        lo = [bounds[j, 0] for j in range(self.dim)]
        hi = [bounds[j, 1] for j in range(self.dim)]
        if not self.periodic_dims:
            yield tuple(lo) + tuple(hi)
            return
        shifts = (0.0, -2 * np.pi, 2 * np.pi)
        for combo in product(shifts, repeat=len(self.periodic_dims)):
            lo2, hi2 = list(lo), list(hi)
            for d, s in zip(self.periodic_dims, combo):
                lo2[d] = lo[d] + s
                hi2[d] = hi[d] + s
            yield tuple(lo2) + tuple(hi2)

    def find_leaf_containing(self, point: np.ndarray):
        pt = tuple(float(v) for v in point)
        best = None
        for i in self.spatial_index.intersection(pt + pt):
            c = self.leaves[i]
            if best is None or c.cell_id < best.cell_id:
                best = c
        return best

    def get_leaves(self) -> List[Cell]:
        return self.leaves

    def get_num_leaves(self) -> int:
        return len(self.leaves)

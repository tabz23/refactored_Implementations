// cell.hpp - hyperrectangular cells and the adaptively refined cell tree.
//
// Cells live in one arena (`cells_`) addressed by id; `leaves_` holds the ids
// of the current leaves in the order the Python list keeps them (surviving
// leaves first, then the children appended in refinement order).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "rtree.hpp"
#include "timing.hpp"
#include "types.hpp"

namespace dhj {

template <std::size_t N>
struct Cell {
    Bounds<N> bounds{};
    State<N> center{};
    std::uint32_t cell_id = 0;
    std::int64_t parent = -1;
    std::int64_t child0 = -1, child1 = -1;
    bool is_leaf = true;

    double V_upper = 0.0, V_lower = 0.0;
    double l_upper = 0.0, l_lower = 0.0;
    double r_upper = 0.0, r_lower = 0.0;

    double range(std::size_t d) const { return bounds[d][1] - bounds[d][0]; }

    // Dimension of maximum range; ties (within 1e-10) go to the lowest index,
    // matching Cell.get_max_range_dim in the paper scripts.
    std::size_t max_range_dim() const {
        double mx = range(0);
        for (std::size_t d = 1; d < N; ++d) mx = std::max(mx, range(d));
        constexpr double tol = 1e-10;
        for (std::size_t d = 0; d + 1 < N; ++d)
            if (std::fabs(range(d) - mx) < tol) return d;
        return N - 1;
    }

    double max_range() const { return range(max_range_dim()); }
};

template <std::size_t N>
inline State<N> center_of(const Bounds<N>& b) {
    State<N> c{};
    for (std::size_t d = 0; d < N; ++d) c[d] = 0.5 * (b[d][0] + b[d][1]);
    return c;
}

template <std::size_t N>
class CellTree {
public:
    using Cell = dhj::Cell<N>;

    CellTree(const Bounds<N>& root_bounds, int initial_resolution, const std::array<bool, N>& periodic)
        : root_bounds_(root_bounds), periodic_(periodic) {
        create_initial_grid(initial_resolution);
        build_spatial_index();
    }

    // Restore constructor (checkpoint resume): arena and leaf list are already
    // consistent; only the spatial index has to be rebuilt.
    CellTree(const Bounds<N>& root_bounds, std::vector<Cell>&& cells, std::vector<std::uint32_t>&& leaves,
             const std::array<bool, N>& periodic)
        : root_bounds_(root_bounds), periodic_(periodic), cells_(std::move(cells)), leaves_(std::move(leaves)) {
        build_spatial_index();
    }

    // The index points into leaf_boxes_, so copying/moving would dangle.
    CellTree(const CellTree&) = delete;
    CellTree& operator=(const CellTree&) = delete;
    CellTree(CellTree&&) = delete;
    CellTree& operator=(CellTree&&) = delete;

    const Bounds<N>& root_bounds() const { return root_bounds_; }
    const std::array<bool, N>& periodic() const { return periodic_; }
    const std::vector<std::uint32_t>& leaves() const { return leaves_; }
    std::size_t num_leaves() const { return leaves_.size(); }
    Cell& cell(std::uint32_t id) { return cells_[id]; }
    const Cell& cell(std::uint32_t id) const { return cells_[id]; }
    std::size_t num_cells() const { return cells_.size(); }

    // Bisects every cell in `to_refine` (leaf order) at the midpoint of its
    // widest dimension. Returns the new child ids (in creation order).
    std::vector<std::uint32_t> refine_cells(const std::vector<std::uint32_t>& to_refine) {
        std::vector<bool> refining(cells_.size(), false);
        for (std::uint32_t id : to_refine) refining[id] = true;

        std::vector<std::uint32_t> new_cells;
        new_cells.reserve(2 * to_refine.size());
        for (std::uint32_t id : to_refine) {
            Cell& parent = cells_[id];
            if (!parent.is_leaf) continue;
            const std::size_t d = parent.max_range_dim();
            const double mid = 0.5 * (parent.bounds[d][0] + parent.bounds[d][1]);

            Cell c0{}, c1{};
            c0.bounds = parent.bounds; c0.bounds[d][1] = mid;
            c1.bounds = parent.bounds; c1.bounds[d][0] = mid;
            c0.center = center_of<N>(c0.bounds);
            c1.center = center_of<N>(c1.bounds);
            c0.cell_id = static_cast<std::uint32_t>(cells_.size());
            c1.cell_id = c0.cell_id + 1;
            c0.parent = c1.parent = static_cast<std::int64_t>(id);

            cells_[id].is_leaf = false;
            cells_[id].child0 = c0.cell_id;
            cells_[id].child1 = c1.cell_id;
            cells_.push_back(c0);
            cells_.push_back(c1);
            new_cells.push_back(c0.cell_id);
            new_cells.push_back(c1.cell_id);
        }

        std::vector<std::uint32_t> kept;
        kept.reserve(leaves_.size() + new_cells.size());
        for (std::uint32_t id : leaves_)
            if (id >= refining.size() || !refining[id]) kept.push_back(id);
        kept.insert(kept.end(), new_cells.begin(), new_cells.end());
        leaves_.swap(kept);
        return new_cells;
    }

    void rebuild_spatial_index() { build_spatial_index(); }

    // Leaves intersecting `b`, with every periodic dimension compared modulo
    // 2*pi (the query is repeated shifted by 0, -2pi, +2pi along each periodic
    // axis, as CellTree.get_intersecting_cells does). `fn` is called once per
    // distinct leaf; `stamp`/`epoch` implement the de-duplication (the caller
    // bumps `epoch` once per de-duplication group).
    template <typename Fn>
    void for_each_intersecting(const Bounds<N>& b, std::vector<std::uint32_t>& stamp, std::uint32_t epoch,
                               Fn&& fn) const {
        static const double shifts[3] = {0.0, -2.0 * kPi, 2.0 * kPi};
        std::size_t n_periodic = 0;
        std::array<std::size_t, N> pdims{};
        for (std::size_t d = 0; d < N; ++d) if (periodic_[d]) pdims[n_periodic++] = d;

        std::size_t combos = 1;
        for (std::size_t k = 0; k < n_periodic; ++k) combos *= 3;

        double lo[N], hi[N];
        for (std::size_t c = 0; c < combos; ++c) {
            for (std::size_t d = 0; d < N; ++d) { lo[d] = b[d][0]; hi[d] = b[d][1]; }
            std::size_t rem = c;
            for (std::size_t k = 0; k < n_periodic; ++k) {
                const double sh = shifts[rem % 3];
                rem /= 3;
                lo[pdims[k]] += sh;
                hi[pdims[k]] += sh;
            }
            index_.query(lo, hi, [&](std::size_t entry) {
                const std::uint32_t leaf_pos = static_cast<std::uint32_t>(entry);
                if (stamp[leaf_pos] == epoch) return;
                stamp[leaf_pos] = epoch;
                fn(leaves_[leaf_pos]);
            });
        }
    }

    double index_build_seconds() const { return index_build_seconds_; }

private:
    void create_initial_grid(int resolution) {
        std::vector<std::vector<double>> edges(N);
        for (std::size_t d = 0; d < N; ++d) {
            edges[d].resize(static_cast<std::size_t>(resolution) + 1);
            const double a = root_bounds_[d][0], b = root_bounds_[d][1];
            // np.linspace(a, b, resolution + 1)
            for (int i = 0; i <= resolution; ++i)
                edges[d][static_cast<std::size_t>(i)] =
                    a + (b - a) * static_cast<double>(i) / static_cast<double>(resolution);
            edges[d][static_cast<std::size_t>(resolution)] = b;
        }
        std::size_t total = 1;
        for (std::size_t d = 0; d < N; ++d) total *= static_cast<std::size_t>(resolution);
        cells_.reserve(total);
        leaves_.reserve(total);

        // itertools.product ordering: the last dimension varies fastest.
        std::array<std::size_t, N> idx{};
        for (std::size_t n = 0; n < total; ++n) {
            Cell c{};
            for (std::size_t d = 0; d < N; ++d) c.bounds[d] = {{edges[d][idx[d]], edges[d][idx[d] + 1]}};
            c.center = center_of<N>(c.bounds);
            c.cell_id = static_cast<std::uint32_t>(cells_.size());
            leaves_.push_back(c.cell_id);
            cells_.push_back(c);
            for (std::size_t d = N; d-- > 0;) {
                if (++idx[d] < static_cast<std::size_t>(resolution)) break;
                idx[d] = 0;
            }
        }
    }

    void build_spatial_index() {
        std::printf("  Building spatial index for %zu cells...\n", leaves_.size());
        const double wall0 = now_seconds();
        leaf_boxes_.resize(leaves_.size());
        for (std::size_t i = 0; i < leaves_.size(); ++i) leaf_boxes_[i] = cells_[leaves_[i]].bounds;
        index_.build(leaf_boxes_);
        index_build_seconds_ = now_seconds() - wall0;
        std::printf("   Spatial index built in %.2fs\n", index_build_seconds_);
    }

    Bounds<N> root_bounds_;
    std::array<bool, N> periodic_{};
    std::vector<Cell> cells_;
    std::vector<std::uint32_t> leaves_;
    std::vector<Bounds<N>> leaf_boxes_;
    BoxIndex<N> index_;
    double index_build_seconds_ = 0.0;
};

}  // namespace dhj

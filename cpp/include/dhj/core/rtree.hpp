// rtree.hpp - static bounding-volume hierarchy over axis-aligned N-D boxes.
//
// Replaces the `rtree` (libspatialindex) index used by the Python CellTree.
// Bulk-built by recursive median splits on the axis with the widest centroid
// spread. Intersection is inclusive of touching faces, matching both
// libspatialindex and Cell.intersects in the paper scripts.
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "types.hpp"

namespace dhj {

template <std::size_t N>
class BoxIndex {
public:
    static constexpr std::size_t kLeafCapacity = 8;
    static constexpr int kMaxTraversalDepth = 64;

    void build(const std::vector<Bounds<N>>& boxes) {
        boxes_ = &boxes;
        const std::size_t n = boxes.size();
        order_.resize(n);
        std::iota(order_.begin(), order_.end(), std::size_t{0});
        nodes_.clear();
        if (n == 0) return;
        nodes_.reserve(2 * (n / kLeafCapacity + 1) + 1);
        build_range(0, n);
    }

    bool empty() const { return nodes_.empty(); }

    // Calls fn(entry_index) for every stored box intersecting [lo, hi].
    template <typename Fn>
    void query(const double lo[N], const double hi[N], Fn&& fn) const {
        if (nodes_.empty()) return;
        std::size_t stack[kMaxTraversalDepth];
        int sp = 0;
        stack[sp++] = 0;
        while (sp > 0) {
            const Node& nd = nodes_[stack[--sp]];
            bool miss = false;
            for (std::size_t d = 0; d < N; ++d)
                if (nd.lo[d] > hi[d] || nd.hi[d] < lo[d]) { miss = true; break; }
            if (miss) continue;
            if (nd.count > 0) {  // leaf
                for (std::size_t k = 0; k < nd.count; ++k) {
                    const std::size_t e = order_[nd.first + k];
                    const Bounds<N>& b = (*boxes_)[e];
                    bool box_miss = false;
                    for (std::size_t d = 0; d < N; ++d)
                        if (b[d][0] > hi[d] || b[d][1] < lo[d]) { box_miss = true; break; }
                    if (box_miss) continue;
                    fn(e);
                }
            } else {
                if (sp + 2 > kMaxTraversalDepth) throw std::runtime_error("BoxIndex: traversal stack overflow");
                stack[sp++] = nd.left;
                stack[sp++] = nd.right;
            }
        }
    }

private:
    struct Node {
        double lo[N], hi[N];
        std::size_t first = 0, count = 0;  // count > 0 marks a leaf
        std::size_t left = 0, right = 0;
    };

    std::size_t build_range(std::size_t first, std::size_t last) {
        const std::size_t self = nodes_.size();
        nodes_.push_back(Node{});
        Node nd{};
        for (std::size_t d = 0; d < N; ++d) { nd.lo[d] = 1e300; nd.hi[d] = -1e300; }
        for (std::size_t i = first; i < last; ++i) {
            const Bounds<N>& b = (*boxes_)[order_[i]];
            for (std::size_t d = 0; d < N; ++d) {
                nd.lo[d] = std::min(nd.lo[d], b[d][0]);
                nd.hi[d] = std::max(nd.hi[d], b[d][1]);
            }
        }
        const std::size_t n = last - first;
        if (n <= kLeafCapacity) {
            nd.first = first;
            nd.count = n;
            nodes_[self] = nd;
            return self;
        }

        std::size_t axis = 0;
        double best_spread = -1.0;
        for (std::size_t d = 0; d < N; ++d) {
            double lo = 1e300, hi = -1e300;
            for (std::size_t i = first; i < last; ++i) {
                const Bounds<N>& b = (*boxes_)[order_[i]];
                const double c = 0.5 * (b[d][0] + b[d][1]);
                lo = std::min(lo, c);
                hi = std::max(hi, c);
            }
            if (hi - lo > best_spread) { best_spread = hi - lo; axis = d; }
        }

        const std::size_t mid = first + n / 2;
        std::nth_element(order_.begin() + static_cast<std::ptrdiff_t>(first),
                         order_.begin() + static_cast<std::ptrdiff_t>(mid),
                         order_.begin() + static_cast<std::ptrdiff_t>(last),
                         [&](std::size_t a, std::size_t b) {
                             const Bounds<N>& ba = (*boxes_)[a];
                             const Bounds<N>& bb = (*boxes_)[b];
                             const double ca = ba[axis][0] + ba[axis][1];
                             const double cb = bb[axis][0] + bb[axis][1];
                             return ca != cb ? ca < cb : a < b;
                         });
        nd.count = 0;
        nd.left = build_range(first, mid);
        nd.right = build_range(mid, last);
        nodes_[self] = nd;
        return self;
    }

    const std::vector<Bounds<N>>* boxes_ = nullptr;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;
};

}  // namespace dhj

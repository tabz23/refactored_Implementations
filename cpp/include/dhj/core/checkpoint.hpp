// checkpoint.hpp - binary checkpoint/resume for the adaptive refinement loop.
//
// A checkpoint holds what is needed to restart Algorithm 2 at the top of a
// refinement phase: the cell arena (geometry, tree links, the six interval
// values per cell), the leaf list in order, and the loop counters. The BVH is
// rebuilt on load and the successor cache is deliberately NOT stored: it is a
// pure function of the leaf set and the dynamics, and the first local sweep
// after a resume re-derives it in full (exact, not an approximation).
//
// On-disk format (little-endian), identical to the earlier ComputingHJ/cpp
// port so its checkpoints remain readable:
//
//   magic "DHJCKPT1" | version u32 | dim u32 | mode u32 | phase u32 |
//   total_refined u64 | initial_resolution u32 | pad u32 | epsilon f64 |
//   root_bounds f64[2*dim] | num_cells u64 | num_leaves u64 |
//   cells: { bounds f64[2*dim], parent i32, child0 i32,
//            V_upper, V_lower, l_upper, l_lower, r_upper, r_lower (f64) } x num_cells |
//   leaves u32[num_leaves]
#pragma once

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "cell.hpp"
#include "types.hpp"

namespace dhj {

constexpr char kCkptMagic[8] = {'D', 'H', 'J', 'C', 'K', 'P', 'T', '1'};
constexpr std::uint32_t kCkptVersion = 1;

template <std::size_t N>
struct CheckpointMeta {
    std::uint32_t phase = 0;
    std::uint64_t total_refined = 0;
    std::uint32_t initial_resolution = 0;
    double epsilon = 0.0;
    Bounds<N> root_bounds{};
};

namespace detail {
inline void put_u32(unsigned char*& p, std::uint32_t v) { std::memcpy(p, &v, 4); p += 4; }
inline void put_i32(unsigned char*& p, std::int32_t v) { std::memcpy(p, &v, 4); p += 4; }
inline void put_u64(unsigned char*& p, std::uint64_t v) { std::memcpy(p, &v, 8); p += 8; }
inline void put_f64(unsigned char*& p, double v) { std::memcpy(p, &v, 8); p += 8; }
inline std::uint32_t get_u32(const unsigned char*& p) { std::uint32_t v; std::memcpy(&v, p, 4); p += 4; return v; }
inline std::int32_t get_i32(const unsigned char*& p) { std::int32_t v; std::memcpy(&v, p, 4); p += 4; return v; }
inline std::uint64_t get_u64(const unsigned char*& p) { std::uint64_t v; std::memcpy(&v, p, 8); p += 8; return v; }
inline double get_f64(const unsigned char*& p) { double v; std::memcpy(&v, p, 8); p += 8; return v; }
inline void must_write(std::FILE* f, const void* p, std::size_t n, const std::string& path) {
    if (std::fwrite(p, 1, n, f) != n) throw std::runtime_error("checkpoint: short write to " + path);
}
inline void must_read(std::FILE* f, void* p, std::size_t n, const std::string& path) {
    if (std::fread(p, 1, n, f) != n) throw std::runtime_error("checkpoint: truncated file " + path);
}
}  // namespace detail

template <std::size_t N>
constexpr std::size_t cell_record_bytes() { return 16 * N + 56; }

// Written to <path>.tmp and renamed, so a kill mid-write leaves the previous
// checkpoint intact.
template <std::size_t N>
inline void save_checkpoint(const std::string& path, Mode mode, const CellTree<N>& tree, const CheckpointMeta<N>& meta) {
    const std::string tmp = path + ".tmp";
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(path).parent_path(), ec);
    std::FILE* f = std::fopen(tmp.c_str(), "wb");
    if (!f) throw std::runtime_error("checkpoint: cannot open " + tmp + " for writing");
    try {
        unsigned char head[8 + 4 * 4 + 8 + 4 + 4 + 8 + 16 * N + 16];
        unsigned char* p = head;
        std::memcpy(p, kCkptMagic, 8); p += 8;
        detail::put_u32(p, kCkptVersion);
        detail::put_u32(p, static_cast<std::uint32_t>(N));
        detail::put_u32(p, static_cast<std::uint32_t>(mode));
        detail::put_u32(p, meta.phase);
        detail::put_u64(p, meta.total_refined);
        detail::put_u32(p, meta.initial_resolution);
        detail::put_u32(p, 0);
        detail::put_f64(p, meta.epsilon);
        for (std::size_t d = 0; d < N; ++d) { detail::put_f64(p, meta.root_bounds[d][0]); detail::put_f64(p, meta.root_bounds[d][1]); }
        detail::put_u64(p, static_cast<std::uint64_t>(tree.num_cells()));
        detail::put_u64(p, static_cast<std::uint64_t>(tree.num_leaves()));
        detail::must_write(f, head, static_cast<std::size_t>(p - head), tmp);

        constexpr std::size_t kBlock = 65536;
        constexpr std::size_t rec = cell_record_bytes<N>();
        std::vector<unsigned char> buf(kBlock * rec);
        std::size_t written = 0;
        while (written < tree.num_cells()) {
            const std::size_t n = std::min(kBlock, tree.num_cells() - written);
            unsigned char* q = buf.data();
            for (std::size_t i = 0; i < n; ++i) {
                const auto& c = tree.cell(static_cast<std::uint32_t>(written + i));
                for (std::size_t d = 0; d < N; ++d) { detail::put_f64(q, c.bounds[d][0]); detail::put_f64(q, c.bounds[d][1]); }
                detail::put_i32(q, static_cast<std::int32_t>(c.parent));
                detail::put_i32(q, static_cast<std::int32_t>(c.child0));
                detail::put_f64(q, c.V_upper); detail::put_f64(q, c.V_lower);
                detail::put_f64(q, c.l_upper); detail::put_f64(q, c.l_lower);
                detail::put_f64(q, c.r_upper); detail::put_f64(q, c.r_lower);
            }
            detail::must_write(f, buf.data(), n * rec, tmp);
            written += n;
        }
        detail::must_write(f, tree.leaves().data(), tree.num_leaves() * sizeof(std::uint32_t), tmp);
    } catch (...) {
        std::fclose(f);
        std::filesystem::remove(tmp, ec);
        throw;
    }
    if (std::fclose(f) != 0) { std::filesystem::remove(tmp, ec); throw std::runtime_error("checkpoint: failed to close " + tmp); }
    std::filesystem::rename(tmp, path, ec);
    if (ec) throw std::runtime_error("checkpoint: cannot rename " + tmp + " -> " + path);
}

// Loads `path` into a new CellTree (out-parameter) and returns the loop state.
// Throws if the file has a different dimension or mode.
template <std::size_t N>
inline CheckpointMeta<N> load_checkpoint(const std::string& path, Mode mode, const std::array<bool, N>& periodic,
                                         std::unique_ptr<CellTree<N>>& tree_out) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) throw std::runtime_error("checkpoint: cannot open " + path + " for reading");
    CheckpointMeta<N> meta;
    std::size_t num_cells = 0, num_leaves = 0;
    std::vector<Cell<N>> cells;
    std::vector<std::uint32_t> leaves;
    try {
        unsigned char head[8 + 4 * 4 + 8 + 4 + 4 + 8 + 16 * N + 16];
        const std::size_t head_bytes = 8 + 4 + 4 + 4 + 4 + 8 + 4 + 4 + 8 + 16 * N + 8 + 8;
        detail::must_read(f, head, head_bytes, path);
        const unsigned char* p = head;
        if (std::memcmp(p, kCkptMagic, 8) != 0) throw std::runtime_error("checkpoint: " + path + " is not a DHJ checkpoint");
        p += 8;
        const std::uint32_t version = detail::get_u32(p);
        if (version != kCkptVersion) throw std::runtime_error("checkpoint: unsupported version " + std::to_string(version));
        const std::uint32_t dim = detail::get_u32(p);
        if (dim != N) throw std::runtime_error("checkpoint: " + path + " is " + std::to_string(dim) + "-D but the dynamics is " + std::to_string(N) + "-D");
        const std::uint32_t file_mode = detail::get_u32(p);
        if (file_mode != static_cast<std::uint32_t>(mode))
            throw std::runtime_error("checkpoint: " + path + " was written in mode '" + mode_name(static_cast<Mode>(file_mode)) + "', not '" + mode_name(mode) + "'");
        meta.phase = detail::get_u32(p);
        meta.total_refined = detail::get_u64(p);
        meta.initial_resolution = detail::get_u32(p);
        (void)detail::get_u32(p);
        meta.epsilon = detail::get_f64(p);
        for (std::size_t d = 0; d < N; ++d) { meta.root_bounds[d][0] = detail::get_f64(p); meta.root_bounds[d][1] = detail::get_f64(p); }
        num_cells = static_cast<std::size_t>(detail::get_u64(p));
        num_leaves = static_cast<std::size_t>(detail::get_u64(p));

        cells.resize(num_cells);
        constexpr std::size_t kBlock = 65536;
        constexpr std::size_t rec = cell_record_bytes<N>();
        std::vector<unsigned char> buf(kBlock * rec);
        std::size_t done = 0;
        while (done < num_cells) {
            const std::size_t n = std::min(kBlock, num_cells - done);
            detail::must_read(f, buf.data(), n * rec, path);
            const unsigned char* q = buf.data();
            for (std::size_t i = 0; i < n; ++i) {
                Cell<N>& c = cells[done + i];
                for (std::size_t d = 0; d < N; ++d) { c.bounds[d][0] = detail::get_f64(q); c.bounds[d][1] = detail::get_f64(q); }
                const std::int32_t parent = detail::get_i32(q);
                const std::int32_t child0 = detail::get_i32(q);
                c.V_upper = detail::get_f64(q); c.V_lower = detail::get_f64(q);
                c.l_upper = detail::get_f64(q); c.l_lower = detail::get_f64(q);
                c.r_upper = detail::get_f64(q); c.r_lower = detail::get_f64(q);
                c.center = center_of<N>(c.bounds);
                c.cell_id = static_cast<std::uint32_t>(done + i);
                c.parent = parent;
                c.child0 = child0;
                c.child1 = child0 < 0 ? -1 : child0 + 1;
                c.is_leaf = child0 < 0;
            }
            done += n;
        }
        leaves.resize(num_leaves);
        detail::must_read(f, leaves.data(), num_leaves * sizeof(std::uint32_t), path);
    } catch (...) {
        std::fclose(f);
        throw;
    }
    std::fclose(f);

    const std::int64_t n = static_cast<std::int64_t>(num_cells);
    for (std::uint32_t id : leaves)
        if (id >= num_cells || !cells[id].is_leaf) throw std::runtime_error("checkpoint: inconsistent leaf list in " + path);
    for (std::size_t i = 0; i < num_cells; ++i) {
        const Cell<N>& c = cells[i];
        if (c.parent < -1 || c.parent >= n || c.child0 < -1 || c.child0 >= n || c.child1 < -1 || c.child1 >= n)
            throw std::runtime_error("checkpoint: out-of-range parent/child index in " + path);
    }
    tree_out.reset(new CellTree<N>(meta.root_bounds, std::move(cells), std::move(leaves), periodic));
    return meta;
}

inline std::string checkpoint_path(const std::string& dir, int phase) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "/checkpoint_phase_%04d.bin", phase);
    return dir + buf;
}

// Highest-numbered checkpoint in `dir`, or "" if none (for --resume auto).
inline std::string latest_checkpoint(const std::string& dir) {
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return {};
    std::string best;
    int best_phase = -1;
    const std::string prefix = "checkpoint_phase_", suffix = ".bin";
    for (const auto& e : std::filesystem::directory_iterator(dir, ec)) {
        const std::string name = e.path().filename().string();
        if (name.size() <= prefix.size() + suffix.size()) continue;
        if (name.compare(0, prefix.size(), prefix) != 0) continue;
        if (name.compare(name.size() - suffix.size(), suffix.size(), suffix) != 0) continue;
        const std::string digits = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
        if (digits.empty() || digits.find_first_not_of("0123456789") != std::string::npos) continue;
        const int phase = std::atoi(digits.c_str());
        if (phase > best_phase) { best_phase = phase; best = e.path().string(); }
    }
    return best;
}

}  // namespace dhj

// compare_checkpoints.cpp - cell-by-cell comparison of two DHJ checkpoints.
//
//   dhj_compare_checkpoints A.bin B.bin [--tol 1e-12]
//
// Both files must have the same dimension. Leaves are matched by cell id and
// their bounds and the six interval values compared; the classification
// (safe / unsafe / boundary) is compared as well. Works on checkpoints from
// this package and from the earlier ComputingHJ/cpp port (same format).
// Exit code 0 = match within tolerance, 1 = mismatch.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

#include "dhj/core/checkpoint.hpp"

namespace {

std::uint32_t peek_dim_and_mode(const std::string& path, std::uint32_t& mode) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", path.c_str()); std::exit(2); }
    unsigned char head[8 + 4 + 4 + 4];
    if (std::fread(head, 1, sizeof head, f) != sizeof head) { std::fclose(f); std::fprintf(stderr, "short file %s\n", path.c_str()); std::exit(2); }
    std::fclose(f);
    std::uint32_t dim;
    std::memcpy(&dim, head + 12, 4);
    std::memcpy(&mode, head + 16, 4);
    return dim;
}

template <std::size_t N>
int compare(const std::string& a, const std::string& b, dhj::Mode mode, double tol) {
    using namespace dhj;
    std::array<bool, N> periodic{};
    std::unique_ptr<CellTree<N>> ta, tb;
    const CheckpointMeta<N> ma = load_checkpoint<N>(a, mode, periodic, ta);
    const CheckpointMeta<N> mb = load_checkpoint<N>(b, mode, periodic, tb);

    std::printf("A: %s\n   phase %u, %zu leaves, %zu cells, total_refined %llu\n", a.c_str(), ma.phase, ta->num_leaves(),
                ta->num_cells(), static_cast<unsigned long long>(ma.total_refined));
    std::printf("B: %s\n   phase %u, %zu leaves, %zu cells, total_refined %llu\n", b.c_str(), mb.phase, tb->num_leaves(),
                tb->num_cells(), static_cast<unsigned long long>(mb.total_refined));

    bool ok = true;
    if (ta->num_leaves() != tb->num_leaves()) { std::printf("LEAF COUNT MISMATCH\n"); ok = false; }
    if (ta->num_cells() != tb->num_cells()) { std::printf("CELL COUNT MISMATCH\n"); ok = false; }
    if (ma.phase != mb.phase) std::printf("note: phases differ (%u vs %u)\n", ma.phase, mb.phase);

    const std::size_t n = std::min(ta->num_leaves(), tb->num_leaves());
    std::size_t id_mismatch = 0, class_mismatch = 0, bound_a = 0, bound_b = 0;
    double d_bounds = 0.0, d_vl = 0.0, d_vu = 0.0, d_lr = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const auto& ca = ta->cell(ta->leaves()[i]);
        const auto& cb = tb->cell(tb->leaves()[i]);
        if (ca.cell_id != cb.cell_id) ++id_mismatch;
        for (std::size_t d = 0; d < N; ++d) {
            d_bounds = std::max(d_bounds, std::fabs(ca.bounds[d][0] - cb.bounds[d][0]));
            d_bounds = std::max(d_bounds, std::fabs(ca.bounds[d][1] - cb.bounds[d][1]));
        }
        d_vl = std::max(d_vl, std::fabs(ca.V_lower - cb.V_lower));
        d_vu = std::max(d_vu, std::fabs(ca.V_upper - cb.V_upper));
        d_lr = std::max(d_lr, std::fabs(ca.l_lower - cb.l_lower));
        d_lr = std::max(d_lr, std::fabs(ca.l_upper - cb.l_upper));
        d_lr = std::max(d_lr, std::fabs(ca.r_lower - cb.r_lower));
        d_lr = std::max(d_lr, std::fabs(ca.r_upper - cb.r_upper));
        const int cls_a = ca.V_lower > 0 ? 0 : (ca.V_upper < 0 ? 1 : 2);
        const int cls_b = cb.V_lower > 0 ? 0 : (cb.V_upper < 0 ? 1 : 2);
        if (cls_a != cls_b) ++class_mismatch;
        if (ca.V_upper > 0 && ca.V_lower <= 0) ++bound_a;
        if (cb.V_upper > 0 && cb.V_lower <= 0) ++bound_b;
    }
    std::printf("leaf id mismatches : %zu\n", id_mismatch);
    std::printf("max |d bounds|     : %.3e\n", d_bounds);
    std::printf("max |d V_lower|    : %.3e\n", d_vl);
    std::printf("max |d V_upper|    : %.3e\n", d_vu);
    std::printf("max |d l,r|        : %.3e\n", d_lr);
    std::printf("class mismatches   : %zu\n", class_mismatch);
    std::printf("boundary cells A/B : %zu / %zu\n", bound_a, bound_b);
    ok = ok && id_mismatch == 0 && d_bounds <= tol && d_vl <= tol && d_vu <= tol && class_mismatch == 0;
    std::printf("RESULT: %s\n", ok ? "MATCH" : "MISMATCH");
    return ok ? 0 : 1;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) { std::fprintf(stderr, "usage: %s A.bin B.bin [--tol T]\n", argv[0]); return 2; }
    double tol = 1e-12;
    for (int i = 3; i + 1 < argc; ++i) if (std::strcmp(argv[i], "--tol") == 0) tol = std::atof(argv[i + 1]);
    std::uint32_t mode_a, mode_b;
    const std::uint32_t dim_a = peek_dim_and_mode(argv[1], mode_a);
    const std::uint32_t dim_b = peek_dim_and_mode(argv[2], mode_b);
    if (dim_a != dim_b) { std::fprintf(stderr, "dimension mismatch: %u vs %u\n", dim_a, dim_b); return 1; }
    if (mode_a != mode_b) { std::fprintf(stderr, "mode mismatch: %u vs %u\n", mode_a, mode_b); return 1; }
    const dhj::Mode mode = static_cast<dhj::Mode>(mode_a);
    try {
        if (dim_a == 3) return compare<3>(argv[1], argv[2], mode, tol);
        if (dim_a == 4) return compare<4>(argv[1], argv[2], mode, tol);
        if (dim_a == 2) return compare<2>(argv[1], argv[2], mode, tol);
        if (dim_a == 5) return compare<5>(argv[1], argv[2], mode, tol);
        if (dim_a == 6) return compare<6>(argv[1], argv[2], mode, tol);
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
    std::fprintf(stderr, "unsupported dimension %u\n", dim_a);
    return 2;
}

// svg.hpp - generic slice plotter (SVG) and per-leaf CSV dump.
//
// Reproduces the figure of the paper scripts: for each slice declared by the
// dynamics a column with three panels (upper bound, lower bound, cell
// classification), every leaf intersecting the slice drawn as a coloured
// rectangle on the plot_dims() plane, obstacle/target overlays, colour bars.
// Periodic dimensions are handled by also matching the pinned value shifted
// by +-2*pi, exactly like dhj/dynamics/plotting.py.
//
// The same data can be rendered with matplotlib from the CSV dump:
//   python/plot_from_csv.py --dynamics <name> <csv>
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include <zlib.h>

#include "../dynamics/base.hpp"
#include "cell.hpp"
#include "types.hpp"

namespace dhj {

struct PlotStyle {
    bool draw_value_cells = true;
    bool cell_labels = false;   // print the value inside large cells (discounted mode)
    bool draw_target = true;    // false in avoid-only mode
    // 800 keeps the vector SVG (and, in Python, the 800 dpi figure). A smaller
    // value writes a raster PNG at that dpi instead.
    int save_dpi = 800;
};

// matplotlib's RdYlGn (11-anchor ColorBrewer map, linear interpolation).
inline void rd_yl_gn(double t, int rgb[3]) {
    static const int anchors[11][3] = {
        {0xa5, 0x00, 0x26}, {0xd7, 0x30, 0x27}, {0xf4, 0x6d, 0x43}, {0xfd, 0xae, 0x61},
        {0xfe, 0xe0, 0x8b}, {0xff, 0xff, 0xbf}, {0xd9, 0xef, 0x8b}, {0xa6, 0xd9, 0x6a},
        {0x66, 0xbd, 0x63}, {0x1a, 0x98, 0x50}, {0x00, 0x68, 0x37}};
    t = std::min(1.0, std::max(0.0, t));
    const double x = t * 10.0;
    int i = static_cast<int>(x);
    if (i > 9) i = 9;
    const double f = x - i;
    for (int c = 0; c < 3; ++c)
        rgb[c] = static_cast<int>(std::lround(anchors[i][c] + f * (anchors[i + 1][c] - anchors[i][c])));
}

inline std::string hex_color(const int rgb[3]) {
    char buf[16];
    std::snprintf(buf, sizeof buf, "#%02x%02x%02x", rgb[0], rgb[1], rgb[2]);
    return buf;
}

inline std::string with_svg_extension(const std::string& path) {
    const std::size_t dot = path.rfind(".png");
    if (dot != std::string::npos && dot + 4 == path.size()) return path.substr(0, dot) + ".svg";
    return path + ".svg";
}

// True if the cell box contains every pinned coordinate of the slice
// (inclusive), with wrap-around on periodic dimensions.
template <std::size_t N>
inline bool cell_in_slice(const Bounds<N>& b, const Slice& sl, const std::array<bool, N>& periodic) {
    for (const auto& fv : sl.fixed) {
        const std::size_t d = fv.first;
        const double v = fv.second;
        const double lo = b[d][0], hi = b[d][1];
        if (periodic[d]) {
            const bool ok = (lo <= v && v <= hi) || (lo <= v - 2 * kPi && v - 2 * kPi <= hi) ||
                            (lo <= v + 2 * kPi && v + 2 * kPi <= hi);
            if (!ok) return false;
        } else if (!(lo <= v && v <= hi)) {
            return false;
        }
    }
    return true;
}

inline void png_chunk(std::ofstream& out, const char type[4], const unsigned char* data, std::uint32_t n) {
    unsigned char len[4] = {static_cast<unsigned char>(n >> 24), static_cast<unsigned char>(n >> 16),
                            static_cast<unsigned char>(n >> 8), static_cast<unsigned char>(n)};
    out.write(reinterpret_cast<char*>(len), 4);
    out.write(type, 4);
    if (n) out.write(reinterpret_cast<const char*>(data), n);
    unsigned crc = static_cast<unsigned>(::crc32(0, reinterpret_cast<const Bytef*>(type), 4));
    if (n) crc = static_cast<unsigned>(::crc32(crc, data, n));
    unsigned char cb[4] = {static_cast<unsigned char>(crc >> 24), static_cast<unsigned char>(crc >> 16),
                           static_cast<unsigned char>(crc >> 8), static_cast<unsigned char>(crc)};
    out.write(reinterpret_cast<char*>(cb), 4);
}

// RGB, row-major, 3 bytes per pixel.
inline bool write_png(const std::string& path, int w, int h, const std::vector<unsigned char>& rgb) {
    std::ofstream out(path, std::ios::binary);
    if (!out) return false;
    const unsigned char sig[8] = {137, 80, 78, 71, 13, 10, 26, 10};
    out.write(reinterpret_cast<const char*>(sig), 8);
    unsigned char ihdr[13] = {static_cast<unsigned char>(w >> 24), static_cast<unsigned char>(w >> 16),
                              static_cast<unsigned char>(w >> 8), static_cast<unsigned char>(w),
                              static_cast<unsigned char>(h >> 24), static_cast<unsigned char>(h >> 16),
                              static_cast<unsigned char>(h >> 8), static_cast<unsigned char>(h),
                              8, 2, 0, 0, 0};
    png_chunk(out, "IHDR", ihdr, 13);
    const std::size_t stride = static_cast<std::size_t>(w) * 3;
    std::vector<unsigned char> raw(static_cast<std::size_t>(h) * (stride + 1));
    for (int y = 0; y < h; ++y) {
        raw[static_cast<std::size_t>(y) * (stride + 1)] = 0;
        std::memcpy(&raw[static_cast<std::size_t>(y) * (stride + 1) + 1], &rgb[static_cast<std::size_t>(y) * stride], stride);
    }
    uLongf bound = compressBound(static_cast<uLong>(raw.size()));
    std::vector<unsigned char> comp(bound);
    if (::compress(comp.data(), &bound, raw.data(), static_cast<uLong>(raw.size())) != Z_OK) return false;
    png_chunk(out, "IDAT", comp.data(), static_cast<std::uint32_t>(bound));
    png_chunk(out, "IEND", nullptr, 0);
    return static_cast<bool>(out);
}

template <std::size_t N>
class SlicePlotter {
public:
    SlicePlotter(const Dynamics<N>& dyn, PlotStyle style) : dyn_(dyn), style_(style), periodic_(dyn.periodic_mask()) {}

    void plot(const CellTree<N>& tree, const std::string& png_filename, int iteration) const {
        // Default (800) is the full SVG. A smaller --plot-dpi rasterizes a PNG.
        if (style_.save_dpi > 0 && style_.save_dpi < 800) {
            plot_raster(tree, png_filename, iteration);
            return;
        }
        const std::vector<Slice> slices = dyn_.slices();
        const int ncol = static_cast<int>(slices.size());
        constexpr int kPanel = 460, kPadL = 60, kPadR = 96, kPadT = 44, kPadB = 52;
        const int cell_w = kPanel + kPadL + kPadR, cell_h = kPanel + kPadT + kPadB;
        const int width = cell_w * ncol, height = cell_h * 3 + 40;

        const std::string path = with_svg_extension(png_filename);
        std::ofstream out(path);
        if (!out) { std::fprintf(stderr, "  [plot] could not open %s\n", path.c_str()); return; }
        out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height
            << "\" viewBox=\"0 0 " << width << " " << height << "\">\n"
            << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n"
            << "<text x=\"" << width / 2 << "\" y=\"26\" font-family=\"sans-serif\" font-size=\"22\""
               " text-anchor=\"middle\">Safety Value Function - Iteration " << iteration << "</text>\n";
        for (int row = 0; row < 3; ++row)
            for (int col = 0; col < ncol; ++col) {
                out << "<g>\n";
                panel(out, tree, slices[static_cast<std::size_t>(col)], row, col * cell_w + kPadL,
                      40 + row * cell_h + kPadT, kPanel);
                out << "</g>\n";
            }
        out << "</svg>\n";
    }

    // Generic per-leaf dump: cell_id, d<k>_lo, d<k>_hi (k < N), values. Same
    // header as dhj.core.vi.write_cells_csv.
    static void dump_csv(const CellTree<N>& tree, const std::string& path) {
        std::ofstream out(path);
        if (!out) return;
        out << "cell_id";
        for (std::size_t k = 0; k < N; ++k) out << ",d" << k << "_lo,d" << k << "_hi";
        out << ",V_lower,V_upper,l_lower,l_upper,r_lower,r_upper\n";
        for (std::uint32_t id : tree.leaves()) {
            const auto& c = tree.cell(id);
            out << c.cell_id;
            for (std::size_t k = 0; k < N; ++k) out << ',' << fmt_g17(c.bounds[k][0]) << ',' << fmt_g17(c.bounds[k][1]);
            out << ',' << fmt_g17(c.V_lower) << ',' << fmt_g17(c.V_upper) << ',' << fmt_g17(c.l_lower) << ','
                << fmt_g17(c.l_upper) << ',' << fmt_g17(c.r_lower) << ',' << fmt_g17(c.r_upper) << '\n';
        }
    }

private:
    // Raster PNG at style_.save_dpi. Figure matches the Python figsize (5 in per
    // slice, 14 in tall). Cells are filled only: a 1 px edge hides the color
    // once cells are only a few pixels wide.
    void plot_raster(const CellTree<N>& tree, const std::string& png_filename, int iteration) const {
        const std::vector<Slice> slices = dyn_.slices();
        const int ncol = std::max(1, static_cast<int>(slices.size()));
        const int dpi = style_.save_dpi;
        const int width = std::max(1, 5 * ncol * dpi);
        const int height = std::max(1, 14 * dpi);
        const int margin_x = std::max(8, width / 40);
        const int margin_y = std::max(8, height / 28);
        const int panel = std::min((width - 2 * margin_x) / ncol, (height - 2 * margin_y) / 3);
        std::vector<unsigned char> img(static_cast<std::size_t>(width) * height * 3, 255);
        auto put = [&](int x, int y, int r, int g, int b) {
            if (x < 0 || y < 0 || x >= width || y >= height) return;
            unsigned char* p = &img[(static_cast<std::size_t>(y) * width + x) * 3];
            p[0] = static_cast<unsigned char>(r);
            p[1] = static_cast<unsigned char>(g);
            p[2] = static_cast<unsigned char>(b);
        };
        auto fill = [&](int x0, int y0, int x1, int y1, int r, int g, int b) {
            if (x1 < x0) std::swap(x0, x1);
            if (y1 < y0) std::swap(y0, y1);
            if (x1 == x0) ++x1;
            if (y1 == y0) ++y1;
            for (int y = y0; y < y1; ++y)
                for (int x = x0; x < x1; ++x) put(x, y, r, g, b);
        };
        const auto pd = dyn_.plot_dims();
        const auto& gb = dyn_.state_bounds();
        const double x_min = gb[pd[0]][0], x_max = gb[pd[0]][1];
        const double y_min = gb[pd[1]][0], y_max = gb[pd[1]][1];
        for (int col = 0; col < ncol; ++col) {
            const Slice& sl = slices[static_cast<std::size_t>(col)];
            std::vector<std::uint32_t> hits;
            std::vector<double> upper, lower;
            for (std::uint32_t id : tree.leaves()) {
                const auto& c = tree.cell(id);
                if (!cell_in_slice<N>(c.bounds, sl, periodic_)) continue;
                hits.push_back(id);
                upper.push_back(c.V_upper);
                lower.push_back(c.V_lower);
            }
            auto range_of = [](const std::vector<double>& v, double& vmin, double& vmax) {
                if (v.empty()) { vmin = 0; vmax = 1; return; }
                vmin = *std::min_element(v.begin(), v.end());
                vmax = *std::max_element(v.begin(), v.end());
                if (std::fabs(vmax - vmin) <= 1e-8 * std::max(std::fabs(vmin), std::fabs(vmax))) { vmin -= 1e-12; vmax += 1e-12; }
                if (vmax <= 0.0) vmax = 0.0;
                else if (vmin >= 0.0) vmin = 0.0;
            };
            double umin, umax, lmin, lmax;
            range_of(upper, umin, umax);
            range_of(lower, lmin, lmax);
            const int ox = margin_x + col * ((width - 2 * margin_x) / ncol);
            for (int row = 0; row < 3; ++row) {
                const int oy = margin_y + row * ((height - 2 * margin_y) / 3);
                const double sx = panel / (x_max - x_min), sy = panel / (y_max - y_min);
                if (row < 2 && style_.draw_value_cells) {
                    const double vmin = row == 0 ? umin : lmin, vmax = row == 0 ? umax : lmax;
                    for (std::size_t i = 0; i < hits.size(); ++i) {
                        const auto& c = tree.cell(hits[i]);
                        const double value = row == 0 ? upper[i] : lower[i];
                        const double t = (vmax > vmin) ? (value - vmin) / (vmax - vmin) : 0.5;
                        int rgb[3];
                        rd_yl_gn(t, rgb);
                        const int x0 = ox + static_cast<int>((c.bounds[pd[0]][0] - x_min) * sx);
                        const int x1 = ox + static_cast<int>((c.bounds[pd[0]][1] - x_min) * sx);
                        const int y0 = oy + static_cast<int>((y_max - c.bounds[pd[1]][1]) * sy);
                        const int y1 = oy + static_cast<int>((y_max - c.bounds[pd[1]][0]) * sy);
                        fill(x0, y0, x1, y1, rgb[0], rgb[1], rgb[2]);
                    }
                } else if (row == 2) {
                    for (std::uint32_t id : hits) {
                        const auto& c = tree.cell(id);
                        int r, g, b;
                        if (c.V_lower > 0.0) { r = 0x2c; g = 0xa0; b = 0x2c; }
                        else if (c.V_upper <= 0.0) { r = 0xd6; g = 0x27; b = 0x28; }
                        else { r = 0x7f; g = 0x7f; b = 0x7f; }
                        const int x0 = ox + static_cast<int>((c.bounds[pd[0]][0] - x_min) * sx);
                        const int x1 = ox + static_cast<int>((c.bounds[pd[0]][1] - x_min) * sx);
                        const int y0 = oy + static_cast<int>((y_max - c.bounds[pd[1]][1]) * sy);
                        const int y1 = oy + static_cast<int>((y_max - c.bounds[pd[1]][0]) * sy);
                        fill(x0, y0, x1, y1, r, g, b);
                    }
                }
                for (const Overlay& ov : dyn_.overlays()) {
                    if (ov.reach_avoid_only && !style_.draw_target) continue;
                    if (ov.kind != "circle") continue;
                    const int cx = ox + static_cast<int>((ov.cx - x_min) * sx);
                    const int cy = oy + static_cast<int>((y_max - ov.cy) * sy);
                    const int rx = std::max(1, static_cast<int>(ov.radius * sx));
                    const int ry = std::max(1, static_cast<int>(ov.radius * sy));
                    int cr = 0, cg = 0, cb = 139;
                    if (ov.color == "orange") { cr = 255; cg = 165; cb = 0; }
                    for (int a = 0; a < 720; ++a) {
                        const double th = a * 3.14159265358979323846 / 360.0;
                        put(cx + static_cast<int>(rx * std::cos(th)), cy + static_cast<int>(ry * std::sin(th)), cr, cg, cb);
                    }
                }
                for (int x = ox; x < ox + panel; ++x) { put(x, oy, 0, 0, 0); put(x, oy + panel - 1, 0, 0, 0); }
                for (int y = oy; y < oy + panel; ++y) { put(ox, y, 0, 0, 0); put(ox + panel - 1, y, 0, 0, 0); }
            }
        }
        (void)iteration;
        std::string path = png_filename;
        const std::size_t dot = path.rfind('.');
        if (dot != std::string::npos) path.resize(dot);
        path += ".png";
        if (!write_png(path, width, height, img))
            std::fprintf(stderr, "  [plot] could not write %s\n", path.c_str());
    }

    void panel(std::ofstream& out, const CellTree<N>& tree, const Slice& sl, int row, int ox, int oy, int size) const {
        const auto pd = dyn_.plot_dims();
        const std::size_t px_d = pd[0], py_d = pd[1];
        const auto& gb = dyn_.state_bounds();
        const double x_min = gb[px_d][0], x_max = gb[px_d][1], y_min = gb[py_d][0], y_max = gb[py_d][1];
        const double sx = size / (x_max - x_min), sy = size / (y_max - y_min);
        auto px = [&](double x) { return ox + (x - x_min) * sx; };
        auto py = [&](double y) { return oy + size - (y - y_min) * sy; };

        std::vector<std::uint32_t> hits;
        std::vector<double> values;
        for (std::uint32_t id : tree.leaves()) {
            const auto& c = tree.cell(id);
            if (!cell_in_slice<N>(c.bounds, sl, periodic_)) continue;
            hits.push_back(id);
            if (row < 2) values.push_back(row == 0 ? c.V_upper : c.V_lower);
        }

        static const char* titles[3] = {"Upper Bound V&#773;_&#947;", "Lower Bound V_&#947;", "Cell Classification"};
        out << "<text x=\"" << ox + size / 2 << "\" y=\"" << oy - 12
            << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">" << titles[row] << " ("
            << xml_escape(sl.label) << ")</text>\n";

        if (row < 2 && values.empty()) {
            out << "<text x=\"" << ox + size / 2 << "\" y=\"" << oy + size / 2
                << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">No data for this slice</text>\n";
            frame(out, ox, oy, size, x_min, x_max, y_min, y_max, px_d, py_d);
            return;
        }

        double vmin = 0.0, vmax = 0.0;
        if (row < 2) {
            vmin = *std::min_element(values.begin(), values.end());
            vmax = *std::max_element(values.begin(), values.end());
            if (std::fabs(vmax - vmin) <= 1e-8 * std::max(std::fabs(vmin), std::fabs(vmax))) { vmin -= 1e-12; vmax += 1e-12; }
            if (vmax <= 0.0) vmax = 0.0; else if (vmin >= 0.0) vmin = 0.0;
        }

        if (row == 2 || style_.draw_value_cells) {
            std::size_t vi = 0;
            for (std::uint32_t id : hits) {
                const auto& c = tree.cell(id);
                std::string fill;
                double value = 0.0, tnorm = 0.0;
                if (row < 2) {
                    value = values[vi++];
                    tnorm = (vmax > vmin) ? (value - vmin) / (vmax - vmin) : 0.5;
                    int rgb[3];
                    rd_yl_gn(tnorm, rgb);
                    fill = hex_color(rgb);
                } else {
                    fill = c.V_lower > 0.0 ? "#2ca02c" : (c.V_upper <= 0.0 ? "#d62728" : "#7f7f7f");
                }
                const double x0 = px(c.bounds[px_d][0]), x1 = px(c.bounds[px_d][1]);
                const double y0 = py(c.bounds[py_d][1]), y1 = py(c.bounds[py_d][0]);
                out << "<rect x=\"" << x0 << "\" y=\"" << y0 << "\" width=\"" << (x1 - x0) << "\" height=\""
                    << (y1 - y0) << "\" fill=\"" << fill
                    << "\" stroke=\"black\" stroke-width=\"0.05\" shape-rendering=\"crispEdges\"/>\n";
                if (row < 2 && style_.cell_labels) {
                    const double w = c.bounds[px_d][1] - c.bounds[px_d][0], h = c.bounds[py_d][1] - c.bounds[py_d][0];
                    if (w > 0.05 && h > 0.05) {
                        const double fs = std::min(std::min(w * 15.0, h * 15.0), 6.0);
                        char lbl[32];
                        std::snprintf(lbl, sizeof lbl, "%.1e", value);
                        out << "<text x=\"" << 0.5 * (x0 + x1) << "\" y=\"" << 0.5 * (y0 + y1)
                            << "\" font-family=\"sans-serif\" font-size=\"" << fs
                            << "\" text-anchor=\"middle\" dominant-baseline=\"middle\" fill=\""
                            << (tnorm > 0.5 ? "black" : "white") << "\">" << lbl << "</text>\n";
                    }
                }
            }
        }

        for (const Overlay& ov : dyn_.overlays()) {
            if (ov.reach_avoid_only && !style_.draw_target) continue;
            if (ov.kind != "circle") continue;
            out << "<ellipse cx=\"" << px(ov.cx) << "\" cy=\"" << py(ov.cy) << "\" rx=\"" << ov.radius * sx
                << "\" ry=\"" << ov.radius * sy << "\" fill=\"none\" stroke=\"" << ov.color << "\" stroke-width=\""
                << ov.linewidth << "\"" << (ov.dashed ? " stroke-dasharray=\"6 4\"" : "") << "/>\n";
        }
        frame(out, ox, oy, size, x_min, x_max, y_min, y_max, px_d, py_d);
        if (row < 2) colorbar(out, ox + size + 14, oy, size, vmin, vmax);
    }

    static std::string xml_escape(const std::string& s) {
        std::string o;
        for (char c : s) {
            if (c == '<') o += "&lt;"; else if (c == '>') o += "&gt;"; else if (c == '&') o += "&amp;"; else o += c;
        }
        return o;
    }

    void frame(std::ofstream& out, int ox, int oy, int size, double x_min, double x_max, double y_min, double y_max,
               std::size_t px_d, std::size_t py_d) const {
        out << "<rect x=\"" << ox << "\" y=\"" << oy << "\" width=\"" << size << "\" height=\"" << size
            << "\" fill=\"none\" stroke=\"black\" stroke-width=\"1\"/>\n";
        auto label = [&](double x, double y, const std::string& txt, const char* anchor) {
            out << "<text x=\"" << x << "\" y=\"" << y << "\" font-family=\"sans-serif\" font-size=\"11\" text-anchor=\""
                << anchor << "\">" << txt << "</text>\n";
        };
        label(ox, oy + size + 16, py_fixed(x_min, 1), "middle");
        label(ox + size, oy + size + 16, py_fixed(x_max, 1), "middle");
        label(ox - 6, oy + size, py_fixed(y_min, 1), "end");
        label(ox - 6, oy + 10, py_fixed(y_max, 1), "end");
        out << "<text x=\"" << ox + size / 2 << "\" y=\"" << oy + size + 34
            << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">" << xml_escape(dyn_.state_name(px_d))
            << "</text>\n";
        out << "<text x=\"" << ox - 40 << "\" y=\"" << oy + size / 2
            << "\" font-family=\"sans-serif\" font-size=\"13\" text-anchor=\"middle\">" << xml_escape(dyn_.state_name(py_d))
            << "</text>\n";
    }

    void colorbar(std::ofstream& out, int x, int y, int size, double vmin, double vmax) const {
        constexpr int kSteps = 64;
        const int w = 14;
        for (int i = 0; i < kSteps; ++i) {
            int rgb[3];
            rd_yl_gn(1.0 - static_cast<double>(i) / (kSteps - 1), rgb);
            out << "<rect x=\"" << x << "\" y=\"" << y + i * size / kSteps << "\" width=\"" << w << "\" height=\""
                << (size / kSteps + 1) << "\" fill=\"" << hex_color(rgb) << "\" shape-rendering=\"crispEdges\"/>\n";
        }
        out << "<rect x=\"" << x << "\" y=\"" << y << "\" width=\"" << w << "\" height=\"" << size
            << "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.7\"/>\n";
        char b[32];
        std::snprintf(b, sizeof b, "%.1e", vmax);
        out << "<text x=\"" << x + w + 4 << "\" y=\"" << y + 9 << "\" font-family=\"sans-serif\" font-size=\"10\">" << b << "</text>\n";
        std::snprintf(b, sizeof b, "%.1e", vmin);
        out << "<text x=\"" << x + w + 4 << "\" y=\"" << y + size << "\" font-family=\"sans-serif\" font-size=\"10\">" << b << "</text>\n";
        if (std::fabs(vmin) >= 0.1 && std::fabs(vmax) >= 0.1 && vmin < 0.0 && vmax > 0.0) {
            const double t = (0.0 - vmin) / (vmax - vmin);
            out << "<text x=\"" << x + w + 4 << "\" y=\"" << y + size - t * size
                << "\" font-family=\"sans-serif\" font-size=\"10\">0.0e+00</text>\n";
        }
    }

    const Dynamics<N>& dyn_;
    PlotStyle style_;
    std::array<bool, N> periodic_;
};

}  // namespace dhj

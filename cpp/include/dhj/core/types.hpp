// types.hpp - fixed-size state/box types shared by the whole solver, the
// solver mode enum, and small Python-compatible formatting helpers.
//
// Everything geometric is templated on the state dimension N so that one build
// serves 3-D (Dubins, evasion) and 4-D (double integrator) systems. The
// dimension is a property of the Dynamics class, not of the build.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

namespace dhj {

constexpr double kPi = 3.14159265358979323846;
constexpr double kInf = std::numeric_limits<double>::infinity();

// Successor id standing for "the reach set left the (non-periodic) domain".
// The sink carries the fixed value -1, exactly as in the paper scripts.
constexpr std::int32_t kOOB = -1;

template <std::size_t N> using State = std::array<double, N>;
template <std::size_t N> using Bounds = std::array<std::array<double, 2>, N>;  // [dim][lo/hi]

// ----------------------------------------------------------------- modes ----

// Which Bellman backup is used. Names match the Python package (dhj.core.modes)
// and the --mode CLI values.
enum class Mode : std::uint32_t {
    RANoDiscount = 0,     // V <- min(l, max(r, best))            (RA_nodiscount.py)
    RADiscount = 1,       // V <- min(l, max(r, gamma * best))    (RA_discount.py)
    AvoidNoDiscount = 2,  // V <- min(l, best)                    (avoid_nodiscount.py)
};

inline const char* mode_name(Mode m) {
    switch (m) {
        case Mode::RANoDiscount: return "ra_nodiscount";
        case Mode::RADiscount: return "ra_discount";
        default: return "avoid_nodiscount";
    }
}

inline const char* mode_script(Mode m) {
    switch (m) {
        case Mode::RANoDiscount: return "RA_nodiscount";
        case Mode::RADiscount: return "RA_discount";
        default: return "avoid_nodiscount";
    }
}

inline bool mode_discounted(Mode m) { return m == Mode::RADiscount; }
inline bool mode_has_target(Mode m) { return m != Mode::AvoidNoDiscount; }

// Default --iterations for Algorithm 1, per original script.
inline int mode_default_iterations(Mode m) {
    switch (m) {
        case Mode::RANoDiscount: return 200;
        case Mode::RADiscount: return 2000;
        default: return 20000;
    }
}

inline bool parse_mode(const std::string& s, Mode& out) {
    if (s == "ra_nodiscount") { out = Mode::RANoDiscount; return true; }
    if (s == "ra_discount") { out = Mode::RADiscount; return true; }
    if (s == "avoid_nodiscount") { out = Mode::AvoidNoDiscount; return true; }
    return false;
}

// ------------------------------------------------- Python-like formatting ----

// Equivalent of Python's repr(float) / str(float): shortest round-tripping
// representation, scientific when exp < -4 or exp >= 16. Used so that the
// results-directory names coincide with the Python package's.
inline std::string py_repr(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    if (v == 0.0) return std::signbit(v) ? "-0.0" : "0.0";

    char buf[64];
    int prec = 0;
    for (; prec <= 17; ++prec) {
        std::snprintf(buf, sizeof buf, "%.*e", prec, v);
        if (std::strtod(buf, nullptr) == v) break;
    }
    std::string s(buf);
    const std::size_t epos = s.find('e');
    std::string mant = s.substr(0, epos);
    const int exp = std::atoi(s.c_str() + epos + 1);

    bool neg = false;
    if (!mant.empty() && mant[0] == '-') { neg = true; mant.erase(0, 1); }
    std::string digits;
    for (char c : mant) if (c != '.') digits += c;
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    const int ndig = static_cast<int>(digits.size());

    std::string out;
    if (exp < -4 || exp >= 16) {
        out = digits.substr(0, 1);
        if (ndig > 1) out += "." + digits.substr(1);
        char eb[16];
        std::snprintf(eb, sizeof eb, "e%c%02d", exp < 0 ? '-' : '+', std::abs(exp));
        out += eb;
    } else if (exp >= 0) {
        if (ndig > exp + 1) out = digits.substr(0, exp + 1) + "." + digits.substr(exp + 1);
        else out = digits + std::string(static_cast<std::size_t>(exp + 1 - ndig), '0') + ".0";
    } else {
        out = "0." + std::string(static_cast<std::size_t>(-exp - 1), '0') + digits;
    }
    return neg ? "-" + out : out;
}

// Python's format(v, '.<prec>e'); C's %e already matches (>= 2 exponent digits).
inline std::string py_exp(double v, int prec) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.*e", prec, v);
    return buf;
}

inline std::string py_fixed(double v, int prec) {
    char buf[512];
    std::snprintf(buf, sizeof buf, "%.*f", prec, v);
    return buf;
}

inline std::string py_bool(bool b) { return b ? "True" : "False"; }

// %r-style output of a double for CSV / JSON (17 significant digits, always
// round-trips). Not the shortest form, but exact.
inline std::string fmt_g17(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

}  // namespace dhj

// base.hpp - the interface a dynamics model implements.
//
// A dynamics header (see dubins.hpp, evasion.hpp, double_integrator.hpp)
// bundles EVERYTHING the solver needs to know about one system:
//
//   * the state box, the finite action set, the ODE right-hand side,
//   * the failure function l(x) (> 0 outside the obstacle) and, for
//     reach-avoid, the reward function r(x) (> 0 inside the target),
//   * the Lipschitz constants L_f (vector field, infinity norm), L_l, L_r,
//   * which state dimensions are periodic (angles): these wrap instead of
//     leaving the domain, in the successor queries and in the plots,
//   * how to draw the result (which 2-D slices, which overlays).
//
// The solver core is templated on the state dimension N and only talks to
// the virtual interface below. To add a system: write a header with a class
// deriving from Dynamics<N>, then register its CLI name in registry.hpp.
//
// Actions are addressed by INDEX. ode() receives the index and looks the
// action up in its own table, so an action may be a scalar (Dubins turn rate)
// or a vector (double integrator acceleration pair) without the core caring.
#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "../core/ode.hpp"
#include "../core/types.hpp"

namespace dhj {

// A shape drawn on top of every panel (obstacle, target, ...).
struct Overlay {
    std::string kind = "circle";   // only circles are needed so far
    double cx = 0.0, cy = 0.0;     // in plot_dims coordinates
    double radius = 0.0;
    std::string color = "black";   // SVG colour name / hex
    bool dashed = false;
    double linewidth = 2.0;
    bool reach_avoid_only = false; // hidden in avoid-only mode (no target set)
};

// A 2-D slice of the state space for plotting: `fixed` pins every dimension
// that is not one of plot_dims().
struct Slice {
    std::vector<std::pair<std::size_t, double>> fixed;
    std::string label;
};

template <std::size_t N>
class Dynamics {
public:
    static constexpr std::size_t kDim = N;
    using State = dhj::State<N>;
    using Bounds = dhj::Bounds<N>;

    Dynamics(double dt, double tau) : dt_(dt), tau_(tau) {
        const double ratio = tau / dt;
        const long n = std::lround(ratio);
        if (std::fabs(ratio - static_cast<double>(n)) > 1e-9 * (1.0 + std::fabs(static_cast<double>(n))))
            throw std::invalid_argument("tau must be evenly divisible by dt");
        n_steps_ = static_cast<int>(n);
    }
    virtual ~Dynamics() = default;

    // ---- required interface -------------------------------------------------
    virtual const char* name() const = 0;                       // CLI name, used in paths
    virtual const Bounds& state_bounds() const = 0;             // [lo, hi] per state
    virtual std::size_t num_actions() const = 0;
    virtual std::vector<double> action(std::size_t i) const = 0;  // the i-th action, as numbers
    virtual State ode(const State& s, std::size_t action_idx) const = 0;  // f(x, u_i)
    virtual double failure_function(const State& s) const = 0;  // l(x)
    virtual double reward_function(const State& s) const = 0;   // r(x)

    // ---- optional overrides -------------------------------------------------
    // Indices of angular (2*pi-periodic) states. Their box must be [-pi, pi].
    virtual std::vector<std::size_t> periodic_dims() const { return {}; }
    // The two dimensions spanning each plot panel.
    virtual std::array<std::size_t, 2> plot_dims() const { return {0, 1}; }
    virtual std::string state_name(std::size_t d) const { return "x" + std::to_string(d); }
    virtual std::vector<Overlay> overlays() const { return {}; }
    // Default: one slice pinning every non-plotted dimension at the middle of
    // its box. Override to choose meaningful slices.
    virtual std::vector<Slice> slices() const {
        Slice s;
        const auto pd = plot_dims();
        for (std::size_t d = 0; d < N; ++d) {
            if (d == pd[0] || d == pd[1]) continue;
            const double v = 0.5 * (state_bounds()[d][0] + state_bounds()[d][1]);
            s.fixed.push_back({d, v});
            if (!s.label.empty()) s.label += ", ";
            s.label += state_name(d) + "=" + py_fixed(v, 2);
        }
        return {s};
    }
    // Extra "key": value pairs (already JSON-formatted) for params.json.
    virtual std::string describe_extra_json() const { return ""; }

    // ---- shared helpers -------------------------------------------------------
    double L_f() const { return L_f_; }
    double L_l() const { return L_l_; }
    double L_r() const { return L_r_; }
    double dt() const { return dt_; }
    double tau() const { return tau_; }
    int n_steps() const { return n_steps_; }

    std::array<bool, N> periodic_mask() const {
        std::array<bool, N> m{};
        for (std::size_t d : periodic_dims()) m[d] = true;
        return m;
    }

    // Wrap periodic coordinates (arctan2(sin, cos), as the paper scripts).
    // Applied to each reported state, never fed back into the integrator.
    void normalize(State& s) const {
        for (std::size_t d : periodic_dims()) s[d] = wrap_angle(s[d]);
    }

    // States at the checkpoint times dt, 2dt, ..., duration, periodic
    // coordinates wrapped. Mirrors Dynamics.dynamics_multi_step in the Python.
    void dynamics_multi_step(const State& s0, std::size_t action_idx, double duration, double dt,
                             std::vector<State>& out) const {
        const int n = static_cast<int>(std::ceil(duration / dt));
        out.clear();
        out.reserve(static_cast<std::size_t>(n));
        State y = s0;
        double t_prev = 0.0;
        const auto f = [this, action_idx](const State& s) { return this->ode(s, action_idx); };
        for (int i = 1; i <= n; ++i) {
            const double t = duration * static_cast<double>(i) / static_cast<double>(n);
            y = integrate<N>(f, y, t - t_prev);
            t_prev = t;
            State w = y;
            normalize(w);
            out.push_back(w);
        }
    }

    // Sanity checks run once at start-up.
    void validate() const {
        for (std::size_t d : periodic_dims()) {
            if (d >= N) throw std::invalid_argument("periodic dimension out of range");
            const auto& b = state_bounds()[d];
            if (std::fabs(b[0] + kPi) > 1e-12 || std::fabs(b[1] - kPi) > 1e-12)
                throw std::invalid_argument("periodic dimension " + std::to_string(d) + " must span [-pi, pi]");
        }
        const auto pd = plot_dims();
        if (pd[0] >= N || pd[1] >= N) throw std::invalid_argument("plot_dims out of range");
        if (!(std::isfinite(L_f_) && std::isfinite(L_l_) && std::isfinite(L_r_)))
            throw std::invalid_argument("set L_f_, L_l_, L_r_ in the dynamics constructor");
        if (num_actions() == 0) throw std::invalid_argument("dynamics has no actions");
    }

    // JSON description for params.json (base fields + describe_extra_json()).
    std::string describe_json(const std::string& indent) const {
        std::string j = "{\n";
        const std::string in = indent + "  ";
        j += in + "\"name\": \"" + name() + "\",\n";
        j += in + "\"dim\": " + std::to_string(N) + ",\n";
        j += in + "\"dt\": " + py_repr(dt_) + ",\n";
        j += in + "\"tau\": " + py_repr(tau_) + ",\n";
        j += in + "\"n_steps\": " + std::to_string(n_steps_) + ",\n";
        j += in + "\"L_f\": " + py_repr(L_f_) + ",\n";
        j += in + "\"L_l\": " + py_repr(L_l_) + ",\n";
        j += in + "\"L_r\": " + py_repr(L_r_) + ",\n";
        j += in + "\"periodic_dims\": [";
        bool first = true;
        for (std::size_t d : periodic_dims()) { j += (first ? "" : ", ") + std::to_string(d); first = false; }
        j += "],\n";
        j += in + "\"state_names\": [";
        for (std::size_t d = 0; d < N; ++d) j += std::string(d ? ", " : "") + "\"" + state_name(d) + "\"";
        j += "],\n";
        j += in + "\"state_bounds\": [";
        for (std::size_t d = 0; d < N; ++d)
            j += std::string(d ? ", " : "") + "[" + py_repr(state_bounds()[d][0]) + ", " +
                 py_repr(state_bounds()[d][1]) + "]";
        j += "],\n";
        j += in + "\"actions\": [";
        for (std::size_t i = 0; i < num_actions(); ++i) {
            const auto a = action(i);
            j += (i ? ", " : "");
            if (a.size() == 1) j += py_repr(a[0]);
            else {
                j += "[";
                for (std::size_t k = 0; k < a.size(); ++k) j += std::string(k ? ", " : "") + py_repr(a[k]);
                j += "]";
            }
        }
        j += "]";
        const std::string extra = describe_extra_json();
        if (!extra.empty()) j += ",\n" + in + extra;
        j += "\n" + indent + "}";
        return j;
    }

protected:
    double dt_, tau_;
    int n_steps_ = 1;
    // Subclasses MUST set these in their constructor.
    double L_f_ = std::nan(""), L_l_ = std::nan(""), L_r_ = std::nan("");
};

}  // namespace dhj

// double_integrator.hpp - planar double integrator, 4-D state (x, y, vx, vy),
// control a = (ax, ay):
//
//     x' = vx,  y' = vy,  vx' = ax,  vy' = ay,     a in {-a_max, 0, a_max}^2
//
// Same reach-avoid sets as the Dubins car: obstacle disc radius 1.3 at the
// origin, target disc radius 0.5 at (2.5, 0). Position box [-3, 3]^2,
// velocity box [-v_max, v_max]^2. No periodic dimension. A reach box leaving
// the velocity box is treated like one leaving the position box (the
// out-of-bounds successor, value -1, is added), so the velocity bounds act as
// a state constraint.
//
// Lipschitz constants (infinity norm):
//   L_f = 1          |f(x)-f(y)|_inf = max(|dvx|, |dvy|, 0, 0) <= |x - y|_inf
//   L_l = L_r = sqrt2  l and r depend on (x, y) only, exactly as for Dubins.
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "base.hpp"

namespace dhj {

class DoubleIntegrator4D : public Dynamics<4> {
public:
    DoubleIntegrator4D(double dt, double tau, double a_max = 1.0, double v_max = 1.0, double obstacle_radius = 1.3,
                       double target_radius = 0.5)
        : Dynamics<4>(dt, tau), a_max_(a_max), v_max_(v_max), obstacle_radius_(obstacle_radius),
          target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-v_max, v_max}}, {{-v_max, v_max}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};
        L_f_ = 1.0;
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        // itertools.product([-a, 0, a], [-a, 0, a]): ax outer, ay inner.
        const double levels[3] = {-a_max, 0.0, a_max};
        for (double ax : levels) for (double ay : levels) actions_.push_back({ax, ay});
    }

    const char* name() const override { return "double_integrator"; }
    const Bounds& state_bounds() const override { return bounds_; }
    std::size_t num_actions() const override { return actions_.size(); }
    std::vector<double> action(std::size_t i) const override { return {actions_[i][0], actions_[i][1]}; }
    std::string state_name(std::size_t d) const override {
        static const char* names[4] = {"x", "y", "vx", "vy"};
        return names[d];
    }

    State ode(const State& s, std::size_t a) const override {
        return {s[2], s[3], actions_[a][0], actions_[a][1]};
    }

    double failure_function(const State& s) const override {
        const double dx = s[0] - obstacle_position_[0], dy = s[1] - obstacle_position_[1];
        return std::sqrt(dx * dx + dy * dy) - obstacle_radius_;
    }

    double reward_function(const State& s) const override {
        const double dx = s[0] - target_position_[0], dy = s[1] - target_position_[1];
        return -(std::sqrt(dx * dx + dy * dy) - target_radius_);
    }

    std::vector<Overlay> overlays() const override {
        Overlay obs; obs.cx = obstacle_position_[0]; obs.cy = obstacle_position_[1]; obs.radius = obstacle_radius_; obs.color = "darkblue";
        Overlay tgt; tgt.cx = target_position_[0]; tgt.cy = target_position_[1]; tgt.radius = target_radius_; tgt.color = "orange";
        tgt.dashed = true; tgt.reach_avoid_only = true;
        return {obs, tgt};
    }

    // (x, y) panels at six velocity pins.
    std::vector<Slice> slices() const override {
        const double v = 0.5 * v_max_;
        const double pins[6][2] = {{0.0, 0.0}, {v, 0.0}, {-v, 0.0}, {0.0, v}, {0.0, -v}, {v, v}};
        std::vector<Slice> out;
        for (const auto& p : pins)
            out.push_back(Slice{{{2, p[0]}, {3, p[1]}}, "vx=" + py_fixed(p[0], 2) + ", vy=" + py_fixed(p[1], 2)});
        return out;
    }

    std::string describe_extra_json() const override {
        return "\"a_max\": " + py_repr(a_max_) + ", \"v_max\": " + py_repr(v_max_) + ", \"obstacle_position\": [" +
               py_repr(obstacle_position_[0]) + ", " + py_repr(obstacle_position_[1]) + "], \"obstacle_radius\": " +
               py_repr(obstacle_radius_) + ", \"target_position\": [" + py_repr(target_position_[0]) + ", " +
               py_repr(target_position_[1]) + "], \"target_radius\": " + py_repr(target_radius_);
    }

private:
    double a_max_, v_max_, obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds bounds_{};
    std::vector<std::array<double, 2>> actions_;
};

}  // namespace dhj

// bicycle.hpp - 5-D kinematic bicycle, state (x, y, θ, v, δ), wheelbase L = 1:
//
//     x' = v cos θ,   y' = v sin θ,   θ' = v tan δ,   v' = a,   δ' = ω
//
// Controls: a and ω each in {-a_max, 0, a_max} and {-ω_max, 0, ω_max} (9 actions).
// Same reach-avoid sets as the Dubins car: obstacle disc radius 1.3 at the
// origin, target disc radius 0.5 at (2.5, 0). Position box [-3, 3]^2.
// θ is periodic on [-π, π]. Speed box [-v_max, v_max], steering box
// [-δ_max, δ_max]. Leaving either box is out of bounds (value -1), as for
// the double integrator's velocity box. δ_max = 0.6 keeps tan δ Lipschitz.
//
// Lipschitz constants (infinity norm) on this box:
//   ||Df||_∞ = max( sqrt(v_max^2 + 1),  |tan δ_max| + v_max sec^2(δ_max) )
//   L_l = L_r = √2    l and r depend on (x, y) only, exactly as for Dubins.
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "base.hpp"

namespace dhj {

class Bicycle5D : public Dynamics<5> {
public:
    Bicycle5D(double dt, double tau, double a_max = 1.0, double omega_max = 1.0, double v_max = 1.0,
              double delta_max = 0.6, double wheelbase = 1.0, double obstacle_radius = 1.3,
              double target_radius = 0.5)
        : Dynamics<5>(dt, tau), a_max_(a_max), omega_max_(omega_max), v_max_(v_max), delta_max_(delta_max),
          wheelbase_(wheelbase), obstacle_radius_(obstacle_radius), target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}, {{-v_max, v_max}}, {{-delta_max, delta_max}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};
        // θ' = (v / L) tan δ, so the heading row of Df scales by 1/L.
        const double td = std::tan(delta_max_);
        const double row_xy = std::sqrt(v_max_ * v_max_ + 1.0);
        const double row_th = (std::fabs(td) + v_max_ * (1.0 + td * td)) / wheelbase_;
        L_f_ = std::max(row_xy, row_th);
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        const double a_levels[3] = {-a_max_, 0.0, a_max_};
        const double w_levels[3] = {-omega_max_, 0.0, omega_max_};
        for (double a : a_levels)
            for (double w : w_levels) actions_.push_back({a, w});
    }

    const char* name() const override { return "bicycle"; }
    const Bounds& state_bounds() const override { return bounds_; }
    std::size_t num_actions() const override { return actions_.size(); }
    std::vector<double> action(std::size_t i) const override { return {actions_[i][0], actions_[i][1]}; }
    std::vector<std::size_t> periodic_dims() const override { return {2}; }
    std::string state_name(std::size_t d) const override {
        static const char* names[5] = {"x", "y", "θ", "v", "δ"};
        return names[d];
    }

    State ode(const State& s, std::size_t a) const override {
        const double th = s[2], v = s[3], delta = s[4];
        return {v * std::cos(th), v * std::sin(th), (v / wheelbase_) * std::tan(delta), actions_[a][0], actions_[a][1]};
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
        Overlay obs;
        obs.cx = obstacle_position_[0];
        obs.cy = obstacle_position_[1];
        obs.radius = obstacle_radius_;
        obs.color = "darkblue";
        Overlay tgt;
        tgt.cx = target_position_[0];
        tgt.cy = target_position_[1];
        tgt.radius = target_radius_;
        tgt.color = "orange";
        tgt.dashed = true;
        tgt.reach_avoid_only = true;
        return {obs, tgt};
    }

    // (x, y) panels, all at full speed so the car is committed to a direction.
    // Straight at the target, straight away, both sideways headings, then both
    // steering stops while aimed at the target. No stopped slice: v = 0 can
    // accelerate either way and does not keep a tail.
    std::vector<Slice> slices() const override {
        const double pins[6][3] = {
            {0.0, v_max_, 0.0},
            {kPi, v_max_, 0.0},
            {0.5 * kPi, v_max_, 0.0},
            {-0.5 * kPi, v_max_, 0.0},
            {0.0, v_max_, delta_max_},
            {0.0, v_max_, -delta_max_},
        };
        std::vector<Slice> out;
        for (const auto& p : pins) {
            out.push_back(Slice{{{2, p[0]}, {3, p[1]}, {4, p[2]}},
                                 "θ=" + py_fixed(p[0], 2) + ", v=" + py_fixed(p[1], 2) + ", δ=" + py_fixed(p[2], 2)});
        }
        return out;
    }

    std::string describe_extra_json() const override {
        return "\"wheelbase\": " + py_repr(wheelbase_) + ", \"a_max\": " + py_repr(a_max_) + ", \"omega_max\": " +
               py_repr(omega_max_) + ", \"v_max\": " + py_repr(v_max_) + ", \"delta_max\": " + py_repr(delta_max_) +
               ", \"obstacle_position\": [" + py_repr(obstacle_position_[0]) + ", " + py_repr(obstacle_position_[1]) +
               "], \"obstacle_radius\": " + py_repr(obstacle_radius_) + ", \"target_position\": [" +
               py_repr(target_position_[0]) + ", " + py_repr(target_position_[1]) + "], \"target_radius\": " +
               py_repr(target_radius_);
    }

private:
    double a_max_, omega_max_, v_max_, delta_max_, wheelbase_, obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds bounds_{};
    std::vector<std::array<double, 2>> actions_;
};

}  // namespace dhj

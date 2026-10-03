// evasion.hpp - 3-D pursuit/evasion in relative coordinates (x, y, theta).
//
//     x' = -v + v cos(theta) + u y,  y' = v sin(theta) - u x,  theta' = -u,
//     u in linspace(-1, 1, 5)
//
// Obstacle: disc of radius 1.0 at the origin.  Target: disc of radius 0.5 at (2.5, 0).
//
// Lipschitz constants (infinity norm): L_f = 1 + v (|u| <= 1 contributes 1 via
// the rotation terms, the cos/sin terms contribute v); L_l = L_r = sqrt(2).
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "base.hpp"

namespace dhj {

class Evasion : public Dynamics<3> {
public:
    Evasion(double dt, double tau, double v_const = 1.0, double obstacle_radius = 1.0, double target_radius = 0.5)
        : Dynamics<3>(dt, tau), v_const_(v_const), obstacle_radius_(obstacle_radius), target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};
        L_f_ = 1.0 + v_const;
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        actions_.resize(5);
        for (int i = 0; i < 5; ++i) actions_[static_cast<std::size_t>(i)] = -1.0 + 2.0 * i / 4.0;  // np.linspace(-1, 1, 5)
    }

    const char* name() const override { return "evasion"; }
    const Bounds& state_bounds() const override { return bounds_; }
    std::size_t num_actions() const override { return actions_.size(); }
    std::vector<double> action(std::size_t i) const override { return {actions_[i]}; }
    std::vector<std::size_t> periodic_dims() const override { return {2}; }
    std::string state_name(std::size_t d) const override { return d == 0 ? "x" : d == 1 ? "y" : "θ"; }

    State ode(const State& s, std::size_t a) const override {
        const double u = actions_[a], v = v_const_;
        return {-v + v * std::cos(s[2]) + u * s[1], v * std::sin(s[2]) - u * s[0], -u};
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

    std::vector<Slice> slices() const override {
        const double thetas[6] = {0.0, kPi, kPi / 4, -kPi / 4, kPi / 2, -kPi / 2};
        std::vector<Slice> out;
        for (double th : thetas) out.push_back(Slice{{{2, th}}, "θ=" + py_fixed(th, 2) + " rad"});
        return out;
    }

    std::string describe_extra_json() const override {
        return "\"v_const\": " + py_repr(v_const_) + ", \"obstacle_position\": [" + py_repr(obstacle_position_[0]) + ", " +
               py_repr(obstacle_position_[1]) + "], \"obstacle_radius\": " + py_repr(obstacle_radius_) +
               ", \"target_position\": [" + py_repr(target_position_[0]) + ", " + py_repr(target_position_[1]) +
               "], \"target_radius\": " + py_repr(target_radius_);
    }

private:
    double v_const_, obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds bounds_{};
    std::vector<double> actions_;
};

}  // namespace dhj

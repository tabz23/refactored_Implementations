// dubins.hpp - Dubins car, state (x, y, theta), control u = turn rate.
//
//     x' = v cos(theta),  y' = v sin(theta),  theta' = u,   u in {-1, 0, 1}
//
// Obstacle: disc of radius 1.3 at the origin (failure set).
// Target:   disc of radius 0.5 at (2.5, 0) (reach set; ignored in avoid mode).
//
// Lipschitz constants (infinity norm, as used by the Gronwall reach sets):
//   L_f = v           |f(x)-f(y)|_inf = v max(|cos t1 - cos t2|, |sin t1 - sin t2|) <= v |t1 - t2|
//   L_l = L_r = sqrt2  l, r are 1-Lipschitz in the Euclidean position; a cell of
//                     half-width r (inf norm) has positions at most sqrt(2) r away.
#pragma once

#include <array>
#include <cmath>
#include <vector>

#include "base.hpp"

namespace dhj {

class DubinsCar : public Dynamics<3> {
public:
    DubinsCar(double dt, double tau, double v_const = 1.0, double obstacle_radius = 1.3, double target_radius = 0.5)
        : Dynamics<3>(dt, tau), v_const_(v_const), obstacle_radius_(obstacle_radius), target_radius_(target_radius) {
        bounds_ = {{{{-3.0, 3.0}}, {{-3.0, 3.0}}, {{-kPi, kPi}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.5, 0.0};
        L_f_ = v_const;
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        actions_ = {-1.0, 0.0, 1.0};
    }

    const char* name() const override { return "dubins"; }
    const Bounds& state_bounds() const override { return bounds_; }
    std::size_t num_actions() const override { return actions_.size(); }
    std::vector<double> action(std::size_t i) const override { return {actions_[i]}; }
    std::vector<std::size_t> periodic_dims() const override { return {2}; }
    std::string state_name(std::size_t d) const override { return d == 0 ? "x" : d == 1 ? "y" : "θ"; }

    State ode(const State& s, std::size_t a) const override {
        const double u = actions_[a];
        return {v_const_ * std::cos(s[2]), v_const_ * std::sin(s[2]), u};
    }

    double failure_function(const State& s) const override {
        const double dx = s[0] - obstacle_position_[0], dy = s[1] - obstacle_position_[1];
        return std::sqrt(dx * dx + dy * dy) - obstacle_radius_;
    }

    double reward_function(const State& s) const override {
        const double dx = s[0] - target_position_[0], dy = s[1] - target_position_[1];
        return -(std::sqrt(dx * dx + dy * dy) - target_radius_);
    }

    // ---- plotting: six theta slices, obstacle (solid) + target (dashed) --------
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

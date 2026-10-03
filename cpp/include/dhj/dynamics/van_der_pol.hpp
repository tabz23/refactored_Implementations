// van_der_pol.hpp - autonomous van der Pol oscillator, 2-D state (x, y), y = dx/dt:
//
//     x' = y,    y' = mu (1 - x^2) y - x
//
// No control: one action, the value 0, so the max over actions does nothing.
// No periodic coordinate. For mu = 1 the stable limit cycle (period about
// 6.67) lies in x in [-2.01, 2.01], y in [-2.68, 2.68], and stays at least
// 1.53 from the origin. The box [-2.5, 2.5] x [-2.9, 2.9] holds the cycle
// with a small margin and no more: L_f is fixed by the corners.
//
// Reach-avoid sets:
//   obstacle  disc radius 0.7 at the origin (inside the cycle; the unstable focus)
//   target    disc radius 0.4 at (2, 0), which the cycle passes through
//
// Lipschitz constants on [-2.5, 2.5] x [-2.9, 2.9], infinity norm:
//   Df = [ 0 , 1 ;  -1 - 2 mu x y , mu (1 - x^2) ]
//   ||Df||_inf = max(1, |1 + 2 mu x y| + |mu (1 - x^2)|) = 20.75
//   at the corners where x y > 0. l and r are Euclidean signed distances, so
//   L_l = L_r = sqrt(2).
//
// L_f = 20.75 still makes exp(L_f * tau) large at the paper's
// tau = 0.3. Use tau = dt = 0.02 or 0.05 if the cycle itself should stay
// visible; the positive set then grows along the orbit one step per sweep.
#pragma once

#include <cmath>
#include <vector>

#include "base.hpp"

namespace dhj {

class VanDerPol : public Dynamics<2> {
public:
    VanDerPol(double dt, double tau, double mu = 1.0, double obstacle_radius = 0.7,
              double target_radius = 0.4)
        : Dynamics<2>(dt, tau), mu_(mu), obstacle_radius_(obstacle_radius), target_radius_(target_radius) {
        bounds_ = {{{{-2.5, 2.5}}, {{-2.9, 2.9}}}};
        obstacle_position_ = {0.0, 0.0};
        target_position_ = {2.0, 0.0};
        const double a = 0.5 * (bounds_[0][1] - bounds_[0][0]);
        const double b = 0.5 * (bounds_[1][1] - bounds_[1][0]);
        // Sharp ||Df||_inf on this symmetric box (a >= sqrt(2)).
        L_f_ = std::max(1.0, 1.0 + 2.0 * mu_ * a * b + mu_ * (a * a - 1.0));
        L_l_ = std::sqrt(2.0);
        L_r_ = std::sqrt(2.0);
        actions_ = {0.0};
    }

    const char* name() const override { return "van_der_pol"; }
    const Bounds& state_bounds() const override { return bounds_; }
    std::size_t num_actions() const override { return actions_.size(); }
    std::vector<double> action(std::size_t i) const override { return {actions_[i]}; }
    std::string state_name(std::size_t d) const override { return d == 0 ? "x" : "y"; }

    State ode(const State& s, std::size_t) const override {
        const double x = s[0], y = s[1];
        return {y, mu_ * (1.0 - x * x) * y - x};
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

    std::vector<Slice> slices() const override { return {Slice{{}, "phase plane"}}; }

    std::string describe_extra_json() const override {
        return "\"mu\": " + py_repr(mu_) + ", \"obstacle_position\": [" + py_repr(obstacle_position_[0]) + ", " +
               py_repr(obstacle_position_[1]) + "], \"obstacle_radius\": " + py_repr(obstacle_radius_) +
               ", \"target_position\": [" + py_repr(target_position_[0]) + ", " + py_repr(target_position_[1]) +
               "], \"target_radius\": " + py_repr(target_radius_);
    }

protected:
    double mu_, obstacle_radius_, target_radius_;
    std::array<double, 2> obstacle_position_{}, target_position_{};
    Bounds bounds_{};
    std::vector<double> actions_;
};

// JuliaReach / Althoff spec: no target, unsafe set y >= 2.75.
// l = 2.75 - y is 1-Lipschitz in the infinity norm. Use with avoid_nodiscount.
class VanDerPolVelocityAvoid : public VanDerPol {
public:
    explicit VanDerPolVelocityAvoid(double dt, double tau, double y_unsafe = 2.75)
        : VanDerPol(dt, tau), y_unsafe_(y_unsafe) {
        L_l_ = 1.0;
        L_r_ = 1.0;
    }

    const char* name() const override { return "van_der_pol_avoid"; }

    double failure_function(const State& s) const override { return y_unsafe_ - s[1]; }

    std::vector<Overlay> overlays() const override { return {}; }

    std::string describe_extra_json() const override {
        return "\"mu\": " + py_repr(mu_) + ", \"y_unsafe\": " + py_repr(y_unsafe_) +
               ", \"spec\": \"juliareach_velocity_avoid\"";
    }

private:
    double y_unsafe_;
};

}  // namespace dhj

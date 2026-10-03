// reach.hpp - Gronwall reachability bounds at the intermediate checkpoints.
//
// A trajectory started anywhere in a cell of half-width r (infinity norm)
// stays within r * exp(L_f t) of the trajectory started at the centre; the
// reach box at checkpoint t_i is the centre trajectory inflated by that much.
#pragma once

#include <cmath>
#include <cstdio>
#include <vector>

namespace dhj {

class GronwallReachabilityAnalyzer {
public:
    GronwallReachabilityAnalyzer(double L_f, double dt, double tau, bool use_infinity_norm = true)
        : use_infinity_norm_(use_infinity_norm), L_(L_f), dt_(dt), tau_(tau) {
        n_checkpoints_ = static_cast<int>(std::lround(tau_ / dt_));
        checkpoint_times_.resize(static_cast<std::size_t>(n_checkpoints_));
        growth_factors_.resize(static_cast<std::size_t>(n_checkpoints_));
        for (int i = 0; i < n_checkpoints_; ++i) {
            // np.linspace(dt, tau, n_checkpoints)
            const double t = n_checkpoints_ == 1
                                 ? tau_
                                 : dt_ + (tau_ - dt_) * static_cast<double>(i) / static_cast<double>(n_checkpoints_ - 1);
            checkpoint_times_[static_cast<std::size_t>(i)] = t;
            growth_factors_[static_cast<std::size_t>(i)] = std::exp(L_ * t);
        }
        std::printf("\nGronwall Reachability Initialized:\n");
        std::printf("  Lipschitz constant L = %g\n", L_);
        std::printf("  Checkpoint interval dt = %g\n", dt_);
        std::printf("  Control horizon tau = %g\n", tau_);
        std::printf("  Number of checkpoints = %d\n", n_checkpoints_);
        std::printf("  Checkpoint times: [");
        for (std::size_t i = 0; i < checkpoint_times_.size(); ++i)
            std::printf("%s'%.3f'", i ? ", " : "", checkpoint_times_[i]);
        std::printf("]\n  Growth factors: [");
        for (std::size_t i = 0; i < growth_factors_.size(); ++i)
            std::printf("%s'%.6f'", i ? ", " : "", growth_factors_[i]);
        std::printf("]\n");
    }

    bool use_infinity_norm() const { return use_infinity_norm_; }
    double dt() const { return dt_; }
    double tau() const { return tau_; }
    int n_checkpoints() const { return n_checkpoints_; }
    const std::vector<double>& growth_factors() const { return growth_factors_; }
    const std::vector<double>& checkpoint_times() const { return checkpoint_times_; }

private:
    bool use_infinity_norm_;
    double L_, dt_, tau_;
    int n_checkpoints_ = 1;
    std::vector<double> checkpoint_times_, growth_factors_;
};

}  // namespace dhj

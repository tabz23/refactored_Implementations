// ode.hpp - adaptive Dormand-Prince 5(4) integrator.
//
// Stands in for scipy.integrate.odeint(..., atol=1e-12, rtol=1e-12) used by the
// Python. The systems treated here are smooth and non-stiff, so an explicit RK
// pair reaches the same 1e-12 accuracy the LSODA call asks for; the two agree
// to ~1e-15 on the checkpoint states.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "types.hpp"

namespace dhj {

constexpr double kAtol = 1e-12;
constexpr double kRtol = 1e-12;

namespace dopri5 {
constexpr double a21 = 1.0 / 5;
constexpr double a31 = 3.0 / 40, a32 = 9.0 / 40;
constexpr double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
constexpr double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561, a54 = -212.0 / 729;
constexpr double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247, a64 = 49.0 / 176, a65 = -5103.0 / 18656;
constexpr double a71 = 35.0 / 384, a73 = 500.0 / 1113, a74 = 125.0 / 192, a75 = -2187.0 / 6784, a76 = 11.0 / 84;
// b - b_hat, for the embedded 4th-order error estimate.
constexpr double e1 = 71.0 / 57600, e3 = -71.0 / 16695, e4 = 71.0 / 1920,
                 e5 = -17253.0 / 339200, e6 = 22.0 / 525, e7 = -1.0 / 40;
}  // namespace dopri5

// Integrates the autonomous system y' = f(y) from t = 0 to t = h_total and
// returns y(h_total). F is callable as `State<N> f(const State<N>&)`.
template <std::size_t N, typename F>
inline State<N> integrate(const F& f, State<N> y, double h_total) {
    if (h_total <= 0.0) return y;

    State<N> k1 = f(y);

    // Starting step: Hairer's heuristic.
    double d0 = 0.0, d1 = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        const double sc = kAtol + kRtol * std::fabs(y[i]);
        d0 = std::max(d0, std::fabs(y[i]) / sc);
        d1 = std::max(d1, std::fabs(k1[i]) / sc);
    }
    double h = (d0 < 1e-5 || d1 < 1e-5) ? 1e-6 : 0.01 * (d0 / d1);
    h = std::min(h, h_total);

    double t = 0.0;
    while (t < h_total) {
        h = std::min(h, h_total - t);

        State<N> y2, y3, y4, y5, y6, y7;
        for (std::size_t i = 0; i < N; ++i) y2[i] = y[i] + h * dopri5::a21 * k1[i];
        const State<N> k2 = f(y2);
        for (std::size_t i = 0; i < N; ++i) y3[i] = y[i] + h * (dopri5::a31 * k1[i] + dopri5::a32 * k2[i]);
        const State<N> k3 = f(y3);
        for (std::size_t i = 0; i < N; ++i)
            y4[i] = y[i] + h * (dopri5::a41 * k1[i] + dopri5::a42 * k2[i] + dopri5::a43 * k3[i]);
        const State<N> k4 = f(y4);
        for (std::size_t i = 0; i < N; ++i)
            y5[i] = y[i] + h * (dopri5::a51 * k1[i] + dopri5::a52 * k2[i] + dopri5::a53 * k3[i] + dopri5::a54 * k4[i]);
        const State<N> k5 = f(y5);
        for (std::size_t i = 0; i < N; ++i)
            y6[i] = y[i] + h * (dopri5::a61 * k1[i] + dopri5::a62 * k2[i] + dopri5::a63 * k3[i] +
                                dopri5::a64 * k4[i] + dopri5::a65 * k5[i]);
        const State<N> k6 = f(y6);
        for (std::size_t i = 0; i < N; ++i)
            y7[i] = y[i] + h * (dopri5::a71 * k1[i] + dopri5::a73 * k3[i] + dopri5::a74 * k4[i] +
                                dopri5::a75 * k5[i] + dopri5::a76 * k6[i]);
        const State<N> k7 = f(y7);

        double err = 0.0;
        for (std::size_t i = 0; i < N; ++i) {
            const double e = h * (dopri5::e1 * k1[i] + dopri5::e3 * k3[i] + dopri5::e4 * k4[i] +
                                  dopri5::e5 * k5[i] + dopri5::e6 * k6[i] + dopri5::e7 * k7[i]);
            const double sc = kAtol + kRtol * std::max(std::fabs(y[i]), std::fabs(y7[i]));
            err = std::max(err, std::fabs(e / sc));
        }

        if (err <= 1.0) {   // accept; FSAL: k7 == f(y7)
            t += h;
            y = y7;
            k1 = k7;
        }
        const double factor = (err == 0.0) ? 5.0 : std::min(5.0, std::max(0.2, 0.9 * std::pow(err, -0.2)));
        h *= factor;
        if (h < 1e-14) h = 1e-14;  // guard against stalling
    }
    return y;
}

// arctan2(sin a, cos a): the exact wrap the paper scripts use.
inline double wrap_angle(double a) { return std::atan2(std::sin(a), std::cos(a)); }

}  // namespace dhj

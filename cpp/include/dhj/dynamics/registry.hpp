// registry.hpp - maps the --dynamics name to a dynamics class.
//
// To add a system: include its header and add one `if` block below. The
// callback receives the constructed object, so the state dimension is
// deduced automatically from the class.
#pragma once

#include <string>
#include <vector>

#include "../core/args.hpp"
#include "double_integrator.hpp"
#include "dubins.hpp"
#include "evasion.hpp"
#include "van_der_pol.hpp"

namespace dhj {

inline std::vector<std::string> registered_dynamics() {
    return {"dubins", "evasion", "double_integrator", "van_der_pol", "van_der_pol_avoid"};
}

// Calls fn(dyn) with the dynamics selected by args.dynamics. Returns false if
// the name is unknown.
template <typename Fn>
bool with_dynamics(const Args& args, Fn&& fn) {
    if (args.dynamics == "dubins") {
        DubinsCar dyn(args.dt, args.tau, args.velocity, /*obstacle_radius=*/1.3);
        fn(dyn);
        return true;
    }
    if (args.dynamics == "evasion") {
        Evasion dyn(args.dt, args.tau, args.velocity, /*obstacle_radius=*/1.0);
        fn(dyn);
        return true;
    }
    if (args.dynamics == "double_integrator") {
        DoubleIntegrator4D dyn(args.dt, args.tau, args.a_max, args.v_max);
        fn(dyn);
        return true;
    }
    if (args.dynamics == "van_der_pol") {
        VanDerPol dyn(args.dt, args.tau);
        fn(dyn);
        return true;
    }
    if (args.dynamics == "van_der_pol_avoid") {
        VanDerPolVelocityAvoid dyn(args.dt, args.tau);
        fn(dyn);
        return true;
    }
    return false;
}

}  // namespace dhj

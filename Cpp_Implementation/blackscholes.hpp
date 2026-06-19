// blackscholes.hpp
// Closed-form Black-Scholes price for a European  option.

#ifndef BASELINE_BLACKSCHOLES_HPP
#define BASELINE_BLACKSCHOLES_HPP

#include <cmath>

namespace baseline {

struct OptionSpec {
    double S0    = 100.0;   // spot
    double K     = 100.0;   // strike
    double r     = 0.05;    // risk-free rate
    double sigma = 0.20;    // volatility
    double T     = 1.0;     // maturity (years)
    bool   call  = true;    // true = call, false = put
};

// Standard normal CDF via erfc.
inline double norm_cdf(double x) {
    return 0.5 * std::erfc(-x * 0.70710678118654752440); // 1/sqrt(2)
}

inline double black_scholes(const OptionSpec& o) {
    double vsqrt = o.sigma * std::sqrt(o.T);
    double d1 = (std::log(o.S0 / o.K) + (o.r + 0.5 * o.sigma * o.sigma) * o.T) / vsqrt;
    double d2 = d1 - vsqrt;
    double disc = std::exp(-o.r * o.T);
    if (o.call)
        return o.S0 * norm_cdf(d1) - o.K * disc * norm_cdf(d2);
    else
        return o.K * disc * norm_cdf(-d2) - o.S0 * norm_cdf(-d1);
}

} // namespace baseline
#endif // BASELINE_BLACKSCHOLES_HPP

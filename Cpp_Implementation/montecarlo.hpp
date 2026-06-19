// montecarlo.hpp
// European-option Monte Carlo pricer

//  S_T = A * exp(B * Z),  A = S0*exp((r - sigma^2/2)T),  B = sigma*sqrt(T)
// so A and B are precomputed once (on the ARM in hardware) and the per-path work
// is one exp, one multiply, a subtract, a max, and an accumulate. The design is
// therefore RNG-throughput-bound.

// The accumulator uses Kahan compensated summation: over tens of millions of
// payoffs naive double summation loses low-order bits, and on the FPGA this
// becomes an explicit accumulator-width question, so the reference models the
// "good" accumulation against which fixed-point width can be judged.

#ifndef BASELINE_MONTECARLO_HPP
#define BASELINE_MONTECARLO_HPP

#include <cmath>
#include <cstdint>
#include "rng.hpp"
#include "gaussian.hpp"
#include "blackscholes.hpp"

namespace baseline {

struct McResult {
    double price;      // discounted MC estimate
    double std_error;  // standard error of the estimate
    double ci_lo;      // 95% confidence interval
    double ci_hi;
    uint64_t paths;
};

// Kahan compensated accumulator.
class Kahan {
    double sum_ = 0.0, c_ = 0.0;
public:
    inline void add(double x) {
        double y = x - c_;
        double t = sum_ + y;
        c_ = (t - sum_) - y;
        sum_ = t;
    }
    double value() const { return sum_; }
};

// Pricer function
inline McResult price_european(Rng32& rng, const OptionSpec& o, uint64_t paths) {
    const double drift = (o.r - 0.5 * o.sigma * o.sigma) * o.T;
    const double A     = o.S0 * std::exp(drift);
    const double B     = o.sigma * std::sqrt(o.T);
    const double disc  = std::exp(-o.r * o.T);

    Kahan sum, sumsq;
    GaussianBoxMuller gauss;
    for (uint64_t n = 0; n < paths; ++n) {
        double z = gauss.next([&]{ return rng.next_unit(); });
        double ST = A * std::exp(B * z);
        double payoff = o.call ? (ST - o.K) : (o.K - ST);
        if (payoff < 0.0) payoff = 0.0;
        sum.add(payoff);
        sumsq.add(payoff * payoff);
    }

    double n   = static_cast<double>(paths);
    double mean = sum.value() / n;
    double var  = (sumsq.value() - sum.value() * sum.value() / n) / (n - 1.0);
    if (var < 0.0) var = 0.0; // For roundoff error

    McResult res;
    res.price     = disc * mean;
    res.std_error = disc * std::sqrt(var / n);
    res.ci_lo     = res.price - 1.959963985 * res.std_error;
    res.ci_hi     = res.price + 1.959963985 * res.std_error;
    res.paths     = paths;
    return res;
}

} // namespace baseline
#endif // BASELINE_MONTECARLO_HPP

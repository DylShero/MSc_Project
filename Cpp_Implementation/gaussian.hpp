// gaussian.hpp
// Uniform -> Gaussian transform via the Box-Muller method

#ifndef BASELINE_GAUSSIAN_HPP
#define BASELINE_GAUSSIAN_HPP

#include <cmath>

namespace baseline {

constexpr double TWO_PI = 6.28318530717958647692;

// Stateful one-Gaussian-at-a-time interface: generates a pair, returns one, and
// caches the spare for the next call. 
class GaussianBoxMuller {
    bool   has_spare_ = false;
    double spare_     = 0.0;
public:
    template <class UnitSource>
    double next(UnitSource&& unit) {
        if (has_spare_) { has_spare_ = false; return spare_; }
        double u1 = unit();
        double u2 = unit();
        double r     = std::sqrt(-2.0 * std::log(u1));
        double theta = TWO_PI * u2;
        spare_     = r * std::sin(theta);
        has_spare_ = true;
        return r * std::cos(theta);
    }
    void reset() { has_spare_ = false; }
};

// Stateless variant: two uniforms in, two Gaussians out. Will be same as FPGA implementation
inline void box_muller(double u1, double u2, double& z0, double& z1) {
    double r     = std::sqrt(-2.0 * std::log(u1));
    double theta = TWO_PI * u2;
    z0 = r * std::cos(theta);
    z1 = r * std::sin(theta);
}

} // namespace baseline
#endif // BASELINE_GAUSSIAN_HPP

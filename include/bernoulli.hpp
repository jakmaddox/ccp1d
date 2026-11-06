#pragma once

#include <cmath>

namespace ccp {

/// Bernoulli function B(x) = x / (exp(x) - 1), with B(0) = 1.
///
/// This is the weighting function at the heart of the Scharfetter-Gummel
/// exponential-fitting scheme. A naive evaluation fails in three places:
///
///   * x -> 0    catastrophic cancellation in (exp(x) - 1), then 0/0
///   * x >> 1    exp(x) overflows to inf, giving B = 0 instead of x*exp(-x)
///   * x << -1   exp(x) underflows to 0, but the limit -x is exact
///
/// Each branch below is the asymptotically correct form, so B is smooth,
/// finite and strictly positive over the whole real line.
inline double bernoulli(double x) noexcept {
    constexpr double kSmall = 1.0e-4;  // below this, use the Taylor series
    constexpr double kLarge = 5.0e2;   // above this, exp(x) is not representable

    if (x > kLarge) {
        // B(x) = x / (e^x - 1) ~ x e^{-x}; computed this way it underflows
        // smoothly to zero rather than dividing by inf.
        return x * std::exp(-x);
    }
    if (x < -kLarge) {
        // e^x -> 0, so B(x) -> x / (-1) = -x, exactly.
        return -x;
    }
    if (std::abs(x) < kSmall) {
        // B(x) = 1 - x/2 + x^2/12 - x^4/720 + O(x^6)
        return 1.0 - 0.5 * x + x * x / 12.0;
    }
    // expm1 computes e^x - 1 without the cancellation that (exp(x) - 1) suffers.
    return x / std::expm1(x);
}

}  // namespace ccp

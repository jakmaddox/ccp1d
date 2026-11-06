#pragma once

#include <cassert>
#include <cstddef>
#include <vector>

namespace ccp {

/// Tridiagonal system solved by the Thomas algorithm (unpivoted LU).
///
/// Every implicit update in this code -- Poisson, both continuity equations,
/// and the electron energy equation -- reduces to one of these in 1D, so the
/// scratch arrays are kept as members and reused across the ~10^5 timesteps
/// of a run rather than reallocated each call.
///
/// No pivoting is performed. That is safe here because Scharfetter-Gummel
/// discretisation of the transport operator produces an M-matrix: the
/// diagonal is positive, the off-diagonals are non-positive, and the system
/// is diagonally dominant. Those are exactly the conditions under which the
/// Thomas algorithm is stable, and they are also what guarantees the solution
/// stays non-negative -- densities cannot go negative by a numerical accident.
class TridiagonalSolver {
public:
    explicit TridiagonalSolver(std::size_t n)
        : lower_(n, 0.0), diag_(n, 0.0), upper_(n, 0.0), rhs_(n, 0.0),
          c_star_(n, 0.0), d_star_(n, 0.0), n_(n) {}

    std::size_t size() const noexcept { return n_; }

    /// Zero the system so coefficients can be accumulated for a new solve.
    void reset() noexcept {
        for (std::size_t i = 0; i < n_; ++i) {
            lower_[i] = 0.0;
            diag_[i] = 0.0;
            upper_[i] = 0.0;
            rhs_[i] = 0.0;
        }
    }

    double& lower(std::size_t i) noexcept { return lower_[i]; }
    double& diag(std::size_t i) noexcept { return diag_[i]; }
    double& upper(std::size_t i) noexcept { return upper_[i]; }
    double& rhs(std::size_t i) noexcept { return rhs_[i]; }

    /// Solve in place, writing the result into `out` (which must hold n values).
    void solve(std::vector<double>& out) {
        assert(out.size() == n_);
        if (n_ == 0) return;

        // Forward elimination.
        c_star_[0] = upper_[0] / diag_[0];
        d_star_[0] = rhs_[0] / diag_[0];
        for (std::size_t i = 1; i < n_; ++i) {
            const double m = diag_[i] - lower_[i] * c_star_[i - 1];
            c_star_[i] = upper_[i] / m;
            d_star_[i] = (rhs_[i] - lower_[i] * d_star_[i - 1]) / m;
        }

        // Back substitution.
        out[n_ - 1] = d_star_[n_ - 1];
        for (std::size_t i = n_ - 1; i-- > 0;) {
            out[i] = d_star_[i] - c_star_[i] * out[i + 1];
        }
    }

private:
    std::vector<double> lower_;
    std::vector<double> diag_;
    std::vector<double> upper_;
    std::vector<double> rhs_;
    std::vector<double> c_star_;
    std::vector<double> d_star_;
    std::size_t n_;
};

}  // namespace ccp

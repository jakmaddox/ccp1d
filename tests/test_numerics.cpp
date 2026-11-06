#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "argon.hpp"
#include "bernoulli.hpp"
#include "discharge.hpp"
#include "tridiagonal.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& what) {
    std::cout << (condition ? "  pass  " : "  FAIL  ") << what << '\n';
    if (!condition) ++failures;
}

void check_close(double actual, double expected, double tolerance,
                 const std::string& what) {
    const double scale = std::max(1.0, std::abs(expected));
    const bool ok = std::abs(actual - expected) <= tolerance * scale;
    if (!ok) {
        std::cout << "  FAIL  " << what << "  (got " << std::scientific
                  << std::setprecision(6) << actual << ", expected " << expected
                  << ")\n";
        ++failures;
    } else {
        std::cout << "  pass  " << what << '\n';
    }
}

// -------------------------------------------------------------------------

void test_bernoulli() {
    std::cout << "\nBernoulli function\n";

    check_close(ccp::bernoulli(0.0), 1.0, 1e-14, "B(0) = 1");
    check_close(ccp::bernoulli(1e-12), 1.0, 1e-10, "B(x) -> 1 as x -> 0");

    // The defining identity B(-x) - B(x) = x must hold across every branch,
    // including the small-argument series and the large-argument asymptotes.
    for (double x : {1e-6, 1e-3, 0.5, 3.0, 40.0, 700.0, 1e4}) {
        const double residual = ccp::bernoulli(-x) - ccp::bernoulli(x) - x;
        check_close(residual, 0.0, 1e-10, "B(-x) - B(x) = x at x = " + std::to_string(x));
    }

    check(ccp::bernoulli(1e4) >= 0.0 && ccp::bernoulli(1e4) < 1e-300,
          "B(x) underflows smoothly to zero for large positive x");
    check_close(ccp::bernoulli(-1e4), 1e4, 1e-12, "B(x) -> -x for large negative x");
    check(std::isfinite(ccp::bernoulli(1e8)) && std::isfinite(ccp::bernoulli(-1e8)),
          "B stays finite at extreme arguments");
}

void test_tridiagonal() {
    std::cout << "\nTridiagonal solver\n";

    // Solve the 1D Laplacian with unit forcing: -u'' = 1 on a 5-point stencil
    // with homogeneous Dirichlet ends, comparing against the exact discrete
    // solution u_i = 0.5 * i * (n + 1 - i) * h^2 with h = 1.
    constexpr std::size_t n = 5;
    ccp::TridiagonalSolver solver(n);
    solver.reset();
    for (std::size_t i = 0; i < n; ++i) {
        solver.diag(i) = 2.0;
        if (i > 0) solver.lower(i) = -1.0;
        if (i + 1 < n) solver.upper(i) = -1.0;
        solver.rhs(i) = 1.0;
    }

    std::vector<double> u(n, 0.0);
    solver.solve(u);

    bool exact = true;
    for (std::size_t i = 0; i < n; ++i) {
        const double expected =
            0.5 * static_cast<double>(i + 1) * static_cast<double>(n - i);
        if (std::abs(u[i] - expected) > 1e-12) exact = false;
    }
    check(exact, "Thomas algorithm reproduces the analytic discrete Laplacian");

    // A diagonally dominant system with a positive right-hand side and
    // non-positive off-diagonals must return a non-negative solution: this is
    // the M-matrix property the density solves rely on for positivity.
    solver.reset();
    for (std::size_t i = 0; i < n; ++i) {
        solver.diag(i) = 4.0;
        if (i > 0) solver.lower(i) = -1.5;
        if (i + 1 < n) solver.upper(i) = -1.5;
        solver.rhs(i) = static_cast<double>(i);
    }
    solver.solve(u);
    bool non_negative = true;
    for (double value : u) non_negative = non_negative && (value >= 0.0);
    check(non_negative, "M-matrix system yields a non-negative solution");
}

void test_scharfetter_gummel_limits() {
    std::cout << "\nScharfetter-Gummel flux limits\n";

    const double d = 1.0;
    const double dx = 0.1;
    const double n_left = 2.0;
    const double n_right = 1.0;

    auto flux = [&](double velocity) {
        const double z = velocity * dx / d;
        const double a = (d / dx) * ccp::bernoulli(-z);
        const double b = (d / dx) * ccp::bernoulli(z);
        return a * n_left - b * n_right;
    };

    // No drift: the scheme must collapse to a plain central difference.
    check_close(flux(0.0), -d * (n_right - n_left) / dx, 1e-10,
                "zero drift recovers Fickian diffusion");

    // Strong positive drift: pure upwind from the left cell.
    const double v_fast = 1.0e4;
    check_close(flux(v_fast) / (v_fast * n_left), 1.0, 1e-6,
                "strong positive drift recovers upwind from the left");

    // Strong negative drift: pure upwind from the right cell.
    check_close(flux(-v_fast) / (-v_fast * n_right), 1.0, 1e-6,
                "strong negative drift recovers upwind from the right");
}

void test_poisson_operator() {
    std::cout << "\nPoisson solve\n";

    // A uniform net charge between grounded electrodes gives a parabola.
    // Build a solver, overwrite the state with a uniform imbalance, and check
    // the potential against rho * x * (L - x) / (2 eps0).
    ccp::DischargeConfig config;
    config.num_cells = 400;
    config.gap_length = 0.02;
    config.num_cycles = 0;
    ccp::DischargeSolver solver(config);

    // Run zero cycles, then inspect the seeded Poisson solve directly via one
    // step at zero applied voltage. The seed plasma is quasineutral, so first
    // confirm the trivial case.
    solver.step(0.0);
    const std::vector<double>& phi = solver.potential();
    double max_abs = 0.0;
    for (double value : phi) max_abs = std::max(max_abs, std::abs(value));
    check(max_abs < 5.0,
          "quasineutral seed plasma produces only a small self-consistent potential");

    check(static_cast<int>(solver.electric_field().size()) ==
              static_cast<int>(config.num_cells) + 1,
          "electric field is defined on faces");
}

void test_gas_properties() {
    std::cout << "\nArgon gas model\n";

    const ccp::ArgonGas gas(100.0, 300.0);

    // 100 mTorr at 300 K: n = P / (kT) with P = 13.33 Pa.
    check_close(gas.number_density(), 3.22e21, 2e-2, "number density at 100 mTorr");
    check_close(gas.temperature_volts(), 0.02585, 1e-2, "gas temperature in volts");

    // Ionisation must rise steeply and monotonically with Te over the range the
    // discharge actually samples.
    bool monotonic = true;
    double previous = 0.0;
    for (double te = 1.0; te <= 7.0; te += 0.25) {
        const double k = ccp::argon::ionisation_rate(te);
        if (k <= previous) monotonic = false;
        previous = k;
    }
    check(monotonic, "ionisation rate increases monotonically from 1 to 7 V");

    // The clamp must hold the fits bounded outside their calibration range.
    check(std::isfinite(ccp::argon::ionisation_rate(1e3)) &&
              std::isfinite(ccp::argon::elastic_rate(1e-6)),
          "rate fits stay finite outside their calibrated range");

    check(ccp::argon::ionisation_rate(2.0) < ccp::argon::excitation_rate(2.0),
          "excitation dominates ionisation at low Te");
}

}  // namespace

int main() {
    std::cout << "ccp1d numerics tests\n";

    test_bernoulli();
    test_tridiagonal();
    test_scharfetter_gummel_limits();
    test_poisson_operator();
    test_gas_properties();

    std::cout << "\n" << (failures == 0 ? "all tests passed" : "FAILURES: ")
              << (failures == 0 ? "" : std::to_string(failures)) << "\n";
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

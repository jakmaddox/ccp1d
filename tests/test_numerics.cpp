#include <algorithm>
#include <cstdio>
#include <memory>
#include <stdexcept>
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

    // The identity above is blind to the even part of B: B(-x) - B(x) cancels
    // every even-order Taylor term, so a wrong x^2 coefficient in the series
    // branch would still pass it to 1e-17. Check absolute accuracy against a
    // long-double reference on both sides of the branch point instead.
    for (double x : {1e-6, -1e-6, 5e-5, -5e-5, 9.9e-5, -9.9e-5, 1.1e-4, 1e-2, 1.0, 10.0}) {
        const long double lx = static_cast<long double>(x);
        const double reference = static_cast<double>(lx / std::expm1l(lx));
        const double rel = std::abs(ccp::bernoulli(x) - reference) / reference;
        check(rel < 1e-14, "B(x) matches a long-double reference at x = " +
                               std::to_string(x));
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

    constexpr double kE = ccp::constants::kElementaryCharge;
    constexpr double kEps0 = ccp::constants::kVacuumPermittivity;
    constexpr double kLength = 0.02;
    constexpr double kIons = 1.0e14;

    // A uniform net charge between grounded electrodes gives the parabola
    //     phi(x) = rho x (L - x) / (2 eps0).
    // The interior stencil is exact for a quadratic; the half-cell ghost
    // elimination at the electrodes is not, so the error should be O(dx^2)
    // and fall by four on each halving.
    auto parabola_error = [&](std::size_t cells) {
        ccp::DischargeConfig config;
        config.num_cells = cells;
        config.gap_length = kLength;
        ccp::DischargeSolver solver(config);
        solver.seed_uniform(0.0, kIons, 1.0);
        solver.solve_field_explicit(0.0);

        const double rho = kE * kIons;
        double err = 0.0;
        double peak = 0.0;
        for (std::size_t i = 0; i < cells; ++i) {
            const double x = solver.grid().cell_centre(i);
            const double exact = rho * x * (kLength - x) / (2.0 * kEps0);
            err = std::max(err, std::abs(solver.potential()[i] - exact));
            peak = std::max(peak, exact);
        }
        return err / peak;
    };
    const double e100 = parabola_error(100);
    const double e200 = parabola_error(200);
    const double e400 = parabola_error(400);
    check(e400 < 1.0e-5, "uniform charge reproduces the analytic parabola");
    const double order = std::log2(e200 / e400);
    check(order > 1.9 && order < 2.1,
          "Poisson error is second order in dx (p = " + std::to_string(order) + ")");
    check(e100 > e200 && e200 > e400, "Poisson error falls monotonically with dx");

    // Charge-free gap with the powered electrode at V: the potential is
    // linear, which the ghost elimination represents exactly, and the field
    // must be uniform at -V/L on every face, the two wall faces included.
    // This pins the sign and the half-cell spacing of the wall-face field.
    {
        ccp::DischargeConfig config;
        config.num_cells = 64;
        config.gap_length = kLength;
        ccp::DischargeSolver solver(config);
        solver.seed_uniform(0.0, 0.0, 1.0);
        constexpr double kVolts = 100.0;
        solver.solve_field_explicit(kVolts);

        double phi_err = 0.0;
        for (std::size_t i = 0; i < config.num_cells; ++i) {
            const double x = solver.grid().cell_centre(i);
            phi_err = std::max(phi_err, std::abs(solver.potential()[i] -
                                                 kVolts * (1.0 - x / kLength)));
        }
        double e_err = 0.0;
        for (double e : solver.electric_field())
            e_err = std::max(e_err, std::abs(e - kVolts / kLength));
        check(phi_err < 1e-9 * kVolts, "vacuum gap gives a linear potential");
        check(e_err < 1e-9 * kVolts / kLength,
              "vacuum gap gives E = V/L on every face, walls included");
        check(static_cast<int>(solver.electric_field().size()) ==
                  static_cast<int>(config.num_cells) + 1,
              "electric field is defined on faces");
    }
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

void test_transport_table() {
    std::cout << "\nTransport table\n";

    // A power law is exact under log-log interpolation.
    const ccp::TransportTable coarse({{1.0, 1.0, 1.0, 1.0, 0.0, 0.0},
                                      {100.0, 1.0e-2, 1.0e2, 1.0, 1.0, 1.0}});
    check_close(coarse.mobility_n(10.0), 0.1, 1e-12, "power law interpolates exactly");
    check_close(coarse.diffusivity_n(0.1), 1.0, 1e-12, "low end value is held");
    check_close(coarse.diffusivity_n(1.0e4), 100.0, 1e-12, "high end value is held");
    check(coarse.k_ionisation(0.5) < 1.0e-250, "a zero rate below threshold stays negligible");

    const ccp::TransportTable table = ccp::TransportTable::from_analytic();
    const ccp::ArgonGas gas(1000.0, 300.0);
    const ccp::ElectronKinetics analytic(gas, nullptr);
    const ccp::ElectronKinetics tabulated(
        gas, std::make_shared<const ccp::TransportTable>(table));

    // Over the Te the discharge samples, including the steep onset of k_iz.
    double worst = 0.0;
    for (double te = 0.2; te <= 40.0; te *= 1.0173) {
        const auto error = [](double a, double b) { return std::abs(a / b - 1.0); };
        worst = std::max({worst, error(tabulated.mobility(te), analytic.mobility(te)),
                          error(tabulated.diffusivity(te), analytic.diffusivity(te)),
                          error(tabulated.ionisation_rate(te), analytic.ionisation_rate(te)),
                          error(tabulated.excitation_rate(te), analytic.excitation_rate(te)),
                          error(tabulated.elastic_rate(te), analytic.elastic_rate(te))});
    }
    check(worst < 2.0e-4, "the analytic table reproduces the built-in model");
    check_close(tabulated.ionisation_rate(1.0e3), ccp::argon::ionisation_rate(1.0e3), 1e-12,
                "end value matches the clamp on the fits");

    const std::string path = "ccp1d_transport_roundtrip.tmp";
    table.save(path);
    const ccp::TransportTable reread = ccp::TransportTable::load(path);
    check(reread.rows().size() == table.rows().size() &&
              reread.rows()[1000].k_ionisation == table.rows()[1000].k_ionisation,
          "a saved table reads back exactly");
    std::remove(path.c_str());

    const auto rejected = [](std::vector<ccp::TransportTable::Row> rows) {
        try {
            const ccp::TransportTable bad(std::move(rows));
        } catch (const std::invalid_argument&) {
            return true;
        }
        return false;
    };
    check(rejected({{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}}), "a single row is rejected");
    check(rejected({{2.0, 1.0, 1.0, 1.0, 1.0, 1.0}, {1.0, 1.0, 1.0, 1.0, 1.0, 1.0}}),
          "descending energies are rejected");
    check(rejected({{1.0, 0.0, 1.0, 1.0, 1.0, 1.0}, {2.0, 1.0, 1.0, 1.0, 1.0, 1.0}}),
          "a zero mobility is rejected");
    check(rejected({{1.0, 1.0, 1.0, -1.0, 1.0, 1.0}, {2.0, 1.0, 1.0, 1.0, 1.0, 1.0}}),
          "a negative rate is rejected");
}

void test_divergence_stops_the_run() {
    std::cout << "\nDivergence detection\n";

    // Explicit Poisson on a dense seed breaks the dielectric relaxation limit
    // at once. NaN compares false with everything, so without a check the run
    // finishes normally and reports a drift cap duty of zero.
    ccp::DischargeConfig config;
    config.num_cells = 50;
    config.num_cycles = 2;
    config.semi_implicit_poisson = false;
    ccp::DischargeSolver unstable(config);
    unstable.seed_uniform(1.0e17, 1.0e17, 4.0);
    std::string message;
    try {
        unstable.run(false);
    } catch (const std::runtime_error& error) {
        message = error.what();
    }
    check(message.find("diverged in cycle 0") != std::string::npos,
          "a diverged run throws and names the cycle");

    config.semi_implicit_poisson = true;
    ccp::DischargeSolver stable(config);
    stable.run(false);
    stable.require_finite(config.num_cycles);
    check(true, "a healthy run passes the finite check");
}

}  // namespace

int main() {
    std::cout << "ccp1d numerics tests\n";

    test_bernoulli();
    test_tridiagonal();
    test_scharfetter_gummel_limits();
    test_poisson_operator();
    test_gas_properties();
    test_transport_table();
    test_divergence_stops_the_run();

    std::cout << "\n" << (failures == 0 ? "all tests passed" : "FAILURES: ")
              << (failures == 0 ? "" : std::to_string(failures)) << "\n";
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

// Solution-verification tests: conservation, symmetry and convergence.
//
// test_numerics checks the building blocks in isolation. These tests check
// that the assembled solver obeys the laws the building blocks are supposed
// to add up to. Each is measured from the solver's public state by the
// independent reconstructions in verification.hpp, so a sign error or a
// missing factor in an update shows up as a broken balance instead of being
// reproduced by the check.
//
// Runs at 100 cells and 250-500 steps per cycle to keep the suite under a
// minute. Tolerances were set from measured values with margin; the numbers
// those measurements gave are quoted next to each check.

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "discharge.hpp"
#include "verification.hpp"

namespace {

int failures = 0;

void check(bool condition, const std::string& what, const std::string& detail = "") {
    std::cout << (condition ? "  pass  " : "  FAIL  ") << what;
    if (!detail.empty()) std::cout << "  (" << detail << ")";
    std::cout << '\n';
    if (!condition) ++failures;
}

std::string sci(double value) {
    std::ostringstream out;
    out << std::scientific << std::setprecision(3) << value;
    return out.str();
}

std::string fix(double value) {
    std::ostringstream out;
    out << std::fixed << std::setprecision(3) << value;
    return out.str();
}

ccp::DischargeConfig test_config() {
    ccp::DischargeConfig config;
    config.num_cells = 100;
    config.steps_per_cycle = 500;
    return config;
}

// -------------------------------------------------------------------------

void test_energy_balance() {
    std::cout << "\nGlobal electron energy balance\n";

    // Ohmic + secondary = inelastic + elastic + wall + d/dt(stored), measured
    // over one cycle after 30. The smallest physical term, secondary-electron
    // energy, is 1.4e-3 of input, so 1e-3 is tight enough that a sign error
    // even there fails the test.
    //
    // The residual is not zero because the solver freezes rate coefficients
    // for a Gummel sweep while this check evaluates losses at the end-of-step
    // state. That lag vanishes faster than dt: measured 5.5e-4, 9.3e-5 and
    // 1.4e-5 of input at 500, 1000 and 2000 steps per cycle. The test runs at
    // 1000 to keep a tenfold margin under the tolerance.
    ccp::DischargeConfig config = test_config();
    config.steps_per_cycle = 1000;
    ccp::DischargeSolver solver(config);
    double time = 0.0;
    ccp::verify::advance_cycles(solver, time, 30);
    const ccp::verify::EnergyBudget b =
        ccp::verify::energy_budget_over_cycle(solver, time);

    check(b.relative_residual() < 1.0e-3, "energy budget closes to < 1e-3 of input",
          "residual " + sci(b.relative_residual()) + ", Ohmic " + sci(b.ohmic) +
              " W/m^2");
    check(b.ohmic > 0.0 && b.inelastic > 0.0 && b.elastic > 0.0 && b.wall > 0.0 &&
              b.secondary > 0.0,
          "every budget term has its physical sign");
    check(b.inelastic > 0.8 * b.ohmic,
          "threshold losses dominate the sink at 1 Torr",
          "inelastic/Ohmic " + fix(b.inelastic / b.ohmic));
}

void test_charge_balance() {
    std::cout << "\nGlobal charge balance\n";

    // Ionisation makes electron-ion pairs, so the net charge in the gap can
    // change only through the net wall current, independent of k_iz. This is
    // the check on the species sources: a 1% electron/ion production mismatch
    // reads exactly 1.0e-2 here. Measured 5.0e-6 of gross production at 1000
    // steps per cycle; before the rate coefficients were frozen per Gummel
    // sweep it was 8.6e-5, which this tolerance rejects.
    ccp::DischargeConfig config = test_config();
    config.steps_per_cycle = 1000;
    ccp::DischargeSolver solver(config);
    double time = 0.0;
    ccp::verify::advance_cycles(solver, time, 30);
    const ccp::verify::ChargeBalance q = ccp::verify::charge_balance_over_cycle(solver, time);
    check(q.relative < 2.5e-5, "net charge changes only by the net wall current",
          "imbalance " + sci(q.relative) + " of ion production");
}

void test_undriven_symmetry() {
    std::cout << "\nUndriven decay\n";

    // With no applied voltage and a mirror-symmetric seed, every field must
    // stay mirror-symmetric to roundoff, and E antisymmetric. This is the
    // direct check on the hand-derived wall boundary conditions: any
    // left/right sign inconsistency breaks it at O(1). Measured 4.9e-14.
    ccp::DischargeConfig config = test_config();
    config.rf_voltage = 0.0;
    ccp::DischargeSolver solver(config);
    double time = 0.0;

    auto inventory = [&]() {
        double sum = 0.0;
        for (double value : solver.electron_density()) sum += value;
        return sum;
    };

    const double start = inventory();
    double previous = start;
    bool monotone = true;
    double worst_mirror = 0.0;
    for (int cycle = 0; cycle < 20; ++cycle) {
        for (int s = 0; s < config.steps_per_cycle; ++s) {
            solver.step(time);
            time += solver.timestep();
            using ccp::verify::mirror_distance;
            worst_mirror = std::max(
                {worst_mirror,
                 mirror_distance(solver.electron_density(), solver.electron_density()),
                 mirror_distance(solver.ion_density(), solver.ion_density()),
                 mirror_distance(solver.energy_density(), solver.energy_density())});
        }
        const double now = inventory();
        if (now > previous) monotone = false;
        previous = now;
    }

    const std::vector<double>& e = solver.electric_field();
    double antisym = 0.0;
    double scale = 0.0;
    for (std::size_t f = 0; f < e.size(); ++f) {
        antisym = std::max(antisym, std::abs(e[f] + e[e.size() - 1 - f]));
        scale = std::max(scale, std::abs(e[f]));
    }

    check(worst_mirror < 1.0e-10, "n_e, n_i, n_eps stay mirror-symmetric",
          "worst " + sci(worst_mirror));
    check(antisym < 1.0e-10 * scale, "E stays antisymmetric",
          "worst " + sci(antisym / scale) + " of max |E|");
    check(monotone && previous < start, "electron inventory decays monotonically",
          "N_end/N_start " + fix(previous / start));
}

void test_driven_symmetry() {
    std::cout << "\nDriven half-period symmetry\n";

    // In periodic steady state n(x,t) = n(L-x, t+T/2). Away from steady state
    // the violation should be set by the secular drift, so the test compares
    // the two rather than demanding zero. Measured ratio 1.40 at 40 cycles and
    // 1.02 at 160. A boundary-condition error gives an O(1) residual at any
    // drift.
    ccp::DischargeSolver solver(test_config());
    double time = 0.0;
    ccp::verify::advance_cycles(solver, time, 40);
    const double m1 = ccp::verify::mean_density_over_cycle(solver, time);
    const double m2 = ccp::verify::mean_density_over_cycle(solver, time);
    const double drift = std::abs(m2 - m1) / m1;
    const double mirror = ccp::verify::half_period_mirror_residual(solver, time);

    check(mirror < 3.0 * drift, "mirror residual is bounded by the cycle drift",
          "residual " + sci(mirror) + ", drift " + sci(drift));
    check(mirror < 0.05, "mirror residual is small in absolute terms",
          "residual " + sci(mirror));
}

double mean_density_at(std::size_t cells, int steps, int cycles) {
    ccp::DischargeConfig config = test_config();
    config.num_cells = cells;
    config.steps_per_cycle = steps;
    ccp::DischargeSolver solver(config);
    double time = 0.0;
    ccp::verify::advance_cycles(solver, time, cycles);
    return ccp::verify::mean_density_over_cycle(solver, time);
}

void test_convergence_order() {
    std::cout << "\nObserved order of accuracy\n";

    // Scharfetter-Gummel reduces to upwinding in a drift-dominated sheath and
    // backward Euler is first order, so the expected order is 1 in both. The
    // estimate is p = log2(|q1 - q0| / |q2 - q1|) on successive halvings.
    // Measured 0.84 in space (100/200/400 cells) and 1.15 in time.
    constexpr int kCycles = 20;
    const double c1 = mean_density_at(100, 250, kCycles);
    const double c2 = mean_density_at(200, 250, kCycles);
    const double c3 = mean_density_at(400, 250, kCycles);
    const double p_space = std::log2(std::abs(c2 - c1) / std::abs(c3 - c2));
    check(p_space > 0.6 && p_space < 1.4, "first-order convergence in space",
          "p = " + fix(p_space));

    const double t1 = mean_density_at(100, 125, kCycles);
    const double t2 = c1;  // 100 cells, 250 steps
    const double t3 = mean_density_at(100, 500, kCycles);
    const double p_time = std::log2(std::abs(t2 - t1) / std::abs(t3 - t2));
    check(p_time > 0.6 && p_time < 1.4, "first-order convergence in time",
          "p = " + fix(p_time));
}

void test_current_continuity_scaling() {
    std::cout << "\nTotal current continuity\n";

    // The semi-implicit Poisson predictor enforces Gauss's law against an
    // estimate of the next charge density, not the one the transport solves
    // produce, so J_cond + eps0 dE/dt is not exactly uniform in x. That error
    // must be O(dt): it is the price of removing the dielectric relaxation
    // limit, and this test pins down that it is that and nothing worse.
    auto spread_at = [](int steps) {
        ccp::DischargeConfig config = test_config();
        config.steps_per_cycle = steps;
        ccp::DischargeSolver solver(config);
        double time = 0.0;
        ccp::verify::advance_cycles(solver, time, 20);
        return ccp::verify::current_continuity_over_cycle(solver, time).mean_spread;
    };
    const double s1 = spread_at(250);
    const double s2 = spread_at(500);
    const double s3 = spread_at(1000);
    const double p = std::log2(s2 / s3);
    check(s1 > s2 && s2 > s3, "continuity error falls as dt is refined",
          sci(s1) + " / " + sci(s2) + " / " + sci(s3));
    check(p > 0.6 && p < 1.4, "continuity error is first order in dt",
          "p = " + fix(p));
}

}  // namespace

int main() {
    std::cout << "ccp1d conservation and symmetry tests\n";

    test_energy_balance();
    test_charge_balance();
    test_undriven_symmetry();
    test_driven_symmetry();
    test_convergence_order();
    test_current_continuity_scaling();

    std::cout << "\n" << (failures == 0 ? "all tests passed" : "FAILURES: ")
              << (failures == 0 ? "" : std::to_string(failures)) << "\n";
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

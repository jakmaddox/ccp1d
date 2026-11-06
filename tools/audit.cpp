// Long-form solution-verification report.
//
// The ctest suite checks that conservation, symmetry and convergence hold at
// a cheap resolution. This tool measures the same quantities at a resolution
// you choose, plus the things that are diagnostics rather than pass/fail:
// how much of the answer the numerical guards are setting, how far the rate
// fits are being extrapolated, and Richardson estimates of discretisation
// error.
//
//   ./build/ccp1d_audit                         # defaults, 400 cycles
//   ./build/ccp1d_audit --cycles 3000           # near periodic steady state
//   ./build/ccp1d_audit --convergence           # adds grid and dt studies

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "discharge.hpp"
#include "verification.hpp"

namespace {

using ccp::DischargeConfig;
using ccp::DischargeSolver;

struct Options {
    DischargeConfig config;
    bool convergence = false;
    int convergence_cycles = 100;
};

Options parse(int argc, char** argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        const bool v = (i + 1 < argc);
        if (a == "--cycles" && v) o.config.num_cycles = std::stoi(argv[++i]);
        else if (a == "--cells" && v) o.config.num_cells = std::stoul(argv[++i]);
        else if (a == "--steps" && v) o.config.steps_per_cycle = std::stoi(argv[++i]);
        else if (a == "--pressure" && v) o.config.pressure_mtorr = std::stod(argv[++i]);
        else if (a == "--voltage" && v) o.config.rf_voltage = std::stod(argv[++i]);
        else if (a == "--gummel" && v) o.config.gummel_iterations = std::stoi(argv[++i]);
        else if (a == "--drift-cap" && v) {
            const double f = std::stod(argv[++i]);
            o.config.drift_speed_limit_factor = (f > 0.0) ? f : 1.0e12;
        } else if (a == "--convergence") o.convergence = true;
        else if (a == "--convergence-cycles" && v) o.convergence_cycles = std::stoi(argv[++i]);
        else {
            std::fprintf(stderr,
                         "usage: %s [--cycles n] [--cells n] [--steps n] [--pressure mTorr]\n"
                         "          [--voltage V] [--gummel n] [--drift-cap f]\n"
                         "          [--convergence] [--convergence-cycles n]\n",
                         argv[0]);
            std::exit(EXIT_FAILURE);
        }
    }
    return o;
}

/// Fraction of the cycle's ionisation that comes from cells hotter than
/// te_limit, and the hottest Te sampled.
void rate_fit_extrapolation(DischargeSolver& solver, double& time, double te_limit,
                            double& hot_fraction, double& te_max) {
    const int steps = solver.config().steps_per_cycle;
    const double ngas = solver.gas().number_density();
    double total = 0.0;
    double hot = 0.0;
    te_max = 0.0;
    for (int s = 0; s < steps; ++s) {
        solver.step(time);
        time += solver.timestep();
        const auto& ne = solver.electron_density();
        for (std::size_t i = 0; i < ne.size(); ++i) {
            const double te = (2.0 / 3.0) * solver.mean_energy_at(i);
            const double src = ccp::argon::ionisation_rate(te) * ngas * ne[i];
            total += src;
            if (te > te_limit) hot += src;
            te_max = std::max(te_max, te);
        }
    }
    hot_fraction = total > 0.0 ? hot / total : 0.0;
}

void convergence_study(const DischargeConfig& base, int cycles) {
    auto mean_ne = [&](std::size_t cells, int steps) {
        DischargeConfig c = base;
        c.num_cells = cells;
        c.steps_per_cycle = steps;
        DischargeSolver s(c);
        double t = 0.0;
        ccp::verify::advance_cycles(s, t, cycles);
        return ccp::verify::mean_density_over_cycle(s, t);
    };

    auto report = [](const char* label, const std::vector<int>& n,
                     const std::vector<double>& q) {
        std::printf("\n  %-8s  mean n_e (m^-3)\n", label);
        for (std::size_t k = 0; k < n.size(); ++k)
            std::printf("  %8d  %.6e\n", n[k], q[k]);
        for (std::size_t k = 2; k < q.size(); ++k) {
            const double p = std::log2(std::abs(q[k - 1] - q[k - 2]) /
                                       std::abs(q[k] - q[k - 1]));
            // Richardson with refinement ratio 2: q_exact ~ q[k] + dq / (2^p - 1),
            // so the error in q[k-1] is dq + dq / (2^p - 1), with dq = q[k] - q[k-1].
            const double err = (q[k] - q[k - 1]) / (std::pow(2.0, p) - 1.0) +
                               (q[k] - q[k - 1]);
            std::printf("  observed order %d-%d-%d: p = %.2f, Richardson error at %d: %.2f %%\n",
                        n[k - 2], n[k - 1], n[k], p, n[k - 1], 100.0 * err / q[k - 1]);
        }
    };

    std::printf("\n=== CONVERGENCE (%d cycles each, final-cycle mean n_e) ===\n", cycles);
    const std::vector<int> cells = {100, 200, 400, 800};
    std::vector<double> qc;
    for (int c : cells)
        qc.push_back(mean_ne(static_cast<std::size_t>(c), base.steps_per_cycle));
    report("cells", cells, qc);

    const std::vector<int> steps = {250, 500, 1000, 2000};
    std::vector<double> qt;
    for (int s : steps) qt.push_back(mean_ne(base.num_cells, s));
    report("steps", steps, qt);
}

}  // namespace

int main(int argc, char** argv) {
    const Options opt = parse(argc, argv);
    const DischargeConfig& cfg = opt.config;
    std::printf("ccp1d audit: %.0f mTorr, %.0f V, %zu cells, %d steps/cycle, %d cycles\n",
                cfg.pressure_mtorr, cfg.rf_voltage, cfg.num_cells, cfg.steps_per_cycle,
                cfg.num_cycles);

    DischargeSolver solver(cfg);
    const ccp::CycleDiagnostics d = solver.run(false);
    double time = static_cast<double>(cfg.num_cycles) *
                  static_cast<double>(cfg.steps_per_cycle) * solver.timestep();

    std::printf("\n=== STATE AFTER %d CYCLES ===\n", d.cycles_completed);
    std::printf("  mean n_e (cycle avg)       %.4e m^-3\n", d.mean_electron_density);
    std::printf("  peak n_e (end of cycle)    %.4e m^-3\n", d.peak_electron_density);
    std::printf("  electron Joule power       %.4e W/m^3\n", d.electron_power_density);
    std::printf("  ion Joule power            %.4e W/m^3  (%.1f %% of total)\n",
                d.ion_power_density, 100.0 * d.ion_power_density / d.total_power_density);
    std::printf("  total absorbed power       %.4e W/m^3\n", d.total_power_density);
    std::printf("  drift over last cycle      %.3e\n", d.cycle_drift);

    std::printf("\n=== NUMERICAL GUARDS ===\n");
    std::printf("  peak mean energy in cycle  %.2f eV   (ceiling %.0f eV: %s)\n",
                d.peak_mean_energy_in_cycle, cfg.mean_energy_ceiling,
                d.peak_mean_energy_in_cycle >= cfg.mean_energy_ceiling ? "BINDS"
                                                                       : "inert");
    std::printf("  electron drift cap duty    %.1f %% of interior-face evaluations\n",
                100.0 * d.drift_cap_duty);
    std::printf("  max |mu_e E| / v_cap       %.2f\n", d.drift_cap_max_ratio);

    const ccp::verify::EnergyBudget b = ccp::verify::energy_budget_over_cycle(solver, time);
    std::printf("\n=== ELECTRON ENERGY BALANCE (one cycle, W/m^2 of electrode) ===\n");
    std::printf("  in   Ohmic heating         %+.4e\n", b.ohmic);
    std::printf("       secondary electrons   %+.4e\n", b.secondary);
    std::printf("  out  inelastic             %+.4e\n", b.inelastic);
    std::printf("       elastic               %+.4e\n", b.elastic);
    std::printf("       wall energy flux      %+.4e\n", b.wall);
    std::printf("       d/dt stored           %+.4e\n", b.stored);
    std::printf("  residual                   %+.4e  (%.2e of input)\n", b.residual(),
                b.relative_residual());

    const double m1 = ccp::verify::mean_density_over_cycle(solver, time);
    const double m2 = ccp::verify::mean_density_over_cycle(solver, time);
    const double drift = std::abs(m2 - m1) / m1;
    const double mirror = ccp::verify::half_period_mirror_residual(solver, time);
    std::printf("\n=== HALF-PERIOD MIRROR SYMMETRY ===\n");
    std::printf("  n_e(x,t) vs n_e(L-x,t+T/2)  %.3e   cycle drift %.3e   ratio %.2f\n",
                mirror, drift, mirror / drift);
    std::printf("  (a ratio of order 1 means the residual is the transient, not the BCs)\n");

    const ccp::verify::CurrentContinuity cc =
        ccp::verify::current_continuity_over_cycle(solver, time);
    std::printf("\n=== TOTAL CURRENT CONTINUITY ===\n");
    std::printf("  (max J - min J) / max |J|   mean %.3e   worst %.3e\n", cc.mean_spread,
                cc.worst_spread);
    std::printf("  (O(dt) by construction of the semi-implicit Poisson predictor)\n");

    double hot = 0.0;
    double te_max = 0.0;
    rate_fit_extrapolation(solver, time, 7.0, hot, te_max);
    std::printf("\n=== RATE-FIT EXTRAPOLATION (fits calibrated ~1-7 V) ===\n");
    std::printf("  hottest Te sampled          %.2f V\n", te_max);
    std::printf("  ionisation from Te > 7 V    %.1f %%\n", 100.0 * hot);

    if (opt.convergence) convergence_study(cfg, opt.convergence_cycles);
    return EXIT_SUCCESS;
}

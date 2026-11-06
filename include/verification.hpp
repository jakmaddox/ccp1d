#pragma once

// Solution-verification measurements built only on the solver's public
// interface. Nothing here reaches into the solver's internals, so each check
// is an independent reconstruction of a conservation law or symmetry from
// the state the solver exposes, not a replay of the solver's own arithmetic.
// That independence is what lets a check catch a sign error or a missing
// factor in the update it is checking.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "argon.hpp"
#include "discharge.hpp"

namespace ccp::verify {

/// Electron energy budget over one RF cycle, per unit electrode area [W/m^2].
///
/// Conservation requires
///     ohmic + secondary = inelastic + elastic + wall + stored
/// where `stored` is the change in total electron energy over the cycle
/// divided by the period. It holds whether or not the discharge has reached
/// periodic steady state, so a short run is enough to check it.
struct EnergyBudget {
    double ohmic = 0.0;       // -Gamma_e . E, integrated over the gap
    double secondary = 0.0;   // energy carried in by secondary electrons
    double inelastic = 0.0;   // ionisation + excitation threshold losses
    double elastic = 0.0;     // elastic recoil losses
    double wall = 0.0;        // energy flux out through both electrodes
    double stored = 0.0;      // d/dt of the total electron energy

    double sources() const { return ohmic + secondary; }
    double sinks() const { return inelastic + elastic + wall + stored; }
    double residual() const { return sources() - sinks(); }
    double relative_residual() const {
        return std::abs(residual()) / std::max(1.0e-300, std::abs(sources()));
    }
};

/// Advance one RF cycle from `time` and return the energy budget over it.
inline EnergyBudget energy_budget_over_cycle(DischargeSolver& solver, double& time) {
    constexpr double kE = constants::kElementaryCharge;
    const DischargeConfig& cfg = solver.config();
    const ArgonGas& gas = solver.gas();
    const std::size_t n = solver.grid().num_cells();
    const double dx = solver.grid().dx();
    const double dt = solver.timestep();
    const int steps = cfg.steps_per_cycle;
    const double ngas = gas.number_density();
    const double tgas = gas.temperature_volts();
    const double mass_ratio = constants::kElectronMass / argon::kMass;

    auto total_energy = [&]() {
        double u = 0.0;
        for (double value : solver.energy_density()) u += value * dx;
        return u;
    };

    EnergyBudget b;
    const double u0 = total_energy();
    for (int s = 0; s < steps; ++s) {
        solver.step(time);
        time += dt;

        const auto& ge = solver.electron_flux();
        const auto& gi = solver.ion_flux();
        const auto& ef = solver.electric_field();
        const auto& ne = solver.electron_density();
        const auto& neps = solver.energy_density();

        for (std::size_t i = 0; i < n; ++i) {
            b.ohmic += -0.5 * (ge[i] * ef[i] + ge[i + 1] * ef[i + 1]) * dx;

            const double te = (2.0 / 3.0) * solver.mean_energy_at(i);
            b.inelastic += (argon::ionisation_rate(te) * argon::kIonisationThreshold +
                            argon::excitation_rate(te) * argon::kExcitationThreshold) *
                           ngas * ne[i] * dx;
            b.elastic += 3.0 * mass_ratio * argon::elastic_rate(te) * ngas * ne[i] *
                         std::max(te - tgas, 0.0) * dx;
        }

        // Escaping electrons carry 2 Te each on a quarter-thermal flux, which
        // is v_th n_eps / 3 for a Maxwellian.
        const double te_left = (2.0 / 3.0) * solver.mean_energy_at(0);
        const double te_right = (2.0 / 3.0) * solver.mean_energy_at(n - 1);
        b.wall += ArgonGas::electron_thermal_speed(te_left) * neps[0] / 3.0 +
                  ArgonGas::electron_thermal_speed(te_right) * neps[n - 1] / 3.0;

        // Ion wall fluxes are signed outward-negative on the left.
        b.secondary += cfg.secondary_emission_yield * cfg.secondary_electron_energy *
                       (-gi[0] + gi[n]);
    }
    const double period = dt * static_cast<double>(steps);
    const double to_watts = kE / static_cast<double>(steps);  // eV/s -> W, cycle mean
    b.ohmic *= to_watts;
    b.inelastic *= to_watts;
    b.elastic *= to_watts;
    b.wall *= to_watts;
    b.secondary *= to_watts;
    b.stored = (total_energy() - u0) * kE / period;
    return b;
}

/// Relative L2 distance between a and the mirror image of b.
inline double mirror_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double num = 0.0;
    double den = 0.0;
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) {
        const double d = a[i] - b[n - 1 - i];
        num += d * d;
        den += a[i] * a[i];
    }
    return den > 0.0 ? std::sqrt(num / den) : 0.0;
}

/// Half-period mirror symmetry of a driven discharge, over one RF cycle.
///
/// With equal electrode areas and no blocking capacitor, driving the left
/// electrode at V sin(wt) differs from the mirrored problem only by a uniform
/// potential offset and a half-period shift, so in periodic steady state
///     n(x, t) = n(L - x, t + T/2).
/// Returns the mean over a half cycle of the relative L2 violation for the
/// electron density. Away from steady state the residual tracks the
/// cycle-to-cycle drift, so compare it with that, not with zero.
inline double half_period_mirror_residual(DischargeSolver& solver, double& time) {
    const int steps = solver.config().steps_per_cycle;
    const int half = steps / 2;
    std::vector<std::vector<double>> snapshots(static_cast<std::size_t>(steps));
    for (int s = 0; s < steps; ++s) {
        solver.step(time);
        time += solver.timestep();
        snapshots[static_cast<std::size_t>(s)] = solver.electron_density();
    }
    double sum = 0.0;
    for (int s = 0; s < half; ++s) {
        sum += mirror_distance(snapshots[static_cast<std::size_t>(s)],
                               snapshots[static_cast<std::size_t>(s + half)]);
    }
    return sum / static_cast<double>(half);
}

/// Spatial uniformity of the total current density over one RF cycle.
///
/// In 1D the conduction plus displacement current
///     J = e (Gamma_i - Gamma_e) + eps0 dE/dt
/// must be independent of x. Returns the cycle mean and worst value of
/// (max_x J - min_x J), normalised by the largest |J| seen anywhere in the
/// cycle. Normalising by the instantaneous |J| instead is misleading: near
/// the RF current zero crossing the denominator is ~0.5% of its peak and a
/// tiny absolute spread reads as an O(1) violation.
struct CurrentContinuity {
    double mean_spread = 0.0;
    double worst_spread = 0.0;
};

inline CurrentContinuity current_continuity_over_cycle(DischargeSolver& solver,
                                                      double& time) {
    constexpr double kE = constants::kElementaryCharge;
    constexpr double kEps0 = constants::kVacuumPermittivity;
    const int steps = solver.config().steps_per_cycle;
    const double dt = solver.timestep();

    std::vector<double> spreads;
    spreads.reserve(static_cast<std::size_t>(steps));
    double j_cycle_max = 0.0;
    std::vector<double> e_prev = solver.electric_field();
    for (int s = 0; s < steps; ++s) {
        solver.step(time);
        time += dt;
        const auto& ge = solver.electron_flux();
        const auto& gi = solver.ion_flux();
        const auto& ef = solver.electric_field();
        double jmin = 1.0e300;
        double jmax = -1.0e300;
        for (std::size_t f = 0; f < ef.size(); ++f) {
            const double j = kE * (gi[f] - ge[f]) + kEps0 * (ef[f] - e_prev[f]) / dt;
            jmin = std::min(jmin, j);
            jmax = std::max(jmax, j);
            j_cycle_max = std::max(j_cycle_max, std::abs(j));
        }
        spreads.push_back(jmax - jmin);
        e_prev = ef;
    }

    CurrentContinuity c;
    const double norm = std::max(1.0e-300, j_cycle_max);
    for (double spread : spreads) {
        c.mean_spread += spread / norm;
        c.worst_spread = std::max(c.worst_spread, spread / norm);
    }
    c.mean_spread /= static_cast<double>(steps);
    return c;
}

/// Global charge balance over one RF cycle.
///
/// Ionisation creates electrons and ions in pairs, so whatever the rate
/// coefficients are, the net charge in the gap can only change through the
/// net particle current into the electrodes:
///     d/dt sum (n_i - n_e) dx = -(Gamma_i[N] - Gamma_i[0]) + (Gamma_e[N] - Gamma_e[0])
/// Interior fluxes telescope, so this is independent of the transport
/// discretisation and isolates the sources and the wall fluxes. `relative` is
/// the imbalance divided by the gross ion production over the cycle, which is
/// the scale an electron/ion source mismatch would show up against.
struct ChargeBalance {
    double imbalance = 0.0;   // m^-2 over the cycle
    double production = 0.0;  // m^-2 over the cycle
    double relative = 0.0;
};

inline ChargeBalance charge_balance_over_cycle(DischargeSolver& solver, double& time) {
    const int steps = solver.config().steps_per_cycle;
    const double dt = solver.timestep();
    const double dx = solver.grid().dx();
    const double ngas = solver.gas().number_density();
    const std::size_t n = solver.grid().num_cells();

    auto net = [&]() {
        double q = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            q += (solver.ion_density()[i] - solver.electron_density()[i]) * dx;
        return q;
    };

    ChargeBalance c;
    const double q0 = net();
    double wall = 0.0;
    for (int s = 0; s < steps; ++s) {
        solver.step(time);
        time += dt;
        const auto& ge = solver.electron_flux();
        const auto& gi = solver.ion_flux();
        wall += (-(gi[n] - gi[0]) + (ge[n] - ge[0])) * dt;
        for (std::size_t i = 0; i < n; ++i) {
            const double te = (2.0 / 3.0) * solver.mean_energy_at(i);
            c.production += argon::ionisation_rate(te) * ngas *
                            solver.electron_density()[i] * dx * dt;
        }
    }
    c.imbalance = (net() - q0) - wall;
    c.relative = std::abs(c.imbalance) / std::max(1.0e-300, c.production);
    return c;
}

/// Advance a whole number of RF cycles.
inline void advance_cycles(DischargeSolver& solver, double& time, int cycles) {
    const int steps = solver.config().steps_per_cycle;
    for (int c = 0; c < cycles; ++c) {
        for (int s = 0; s < steps; ++s) {
            solver.step(time);
            time += solver.timestep();
        }
    }
}

/// Gap- and cycle-averaged electron density over one RF cycle.
inline double mean_density_over_cycle(DischargeSolver& solver, double& time) {
    const int steps = solver.config().steps_per_cycle;
    double sum = 0.0;
    for (int s = 0; s < steps; ++s) {
        solver.step(time);
        time += solver.timestep();
        double cell_sum = 0.0;
        for (double value : solver.electron_density()) cell_sum += value;
        sum += cell_sum / static_cast<double>(solver.grid().num_cells());
    }
    return sum / static_cast<double>(steps);
}

}  // namespace ccp::verify

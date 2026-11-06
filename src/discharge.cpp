#include "discharge.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#include "bernoulli.hpp"

namespace ccp {

namespace {
constexpr double kE = constants::kElementaryCharge;
constexpr double kEps0 = constants::kVacuumPermittivity;
}  // namespace

DischargeSolver::DischargeSolver(DischargeConfig config)
    : config_(config),
      grid_(config.num_cells, config.gap_length),
      gas_(config.pressure_mtorr, config.gas_temperature_k),
      kinetics_(gas_, config.transport_table),
      dt_(1.0 / (config.rf_frequency * static_cast<double>(config.steps_per_cycle))),
      omega_(2.0 * constants::kPi * config.rf_frequency),
      // Free fall through the largest potential difference in the gap. With
      // no drive that is the afterglow plasma potential, a few Te, so the
      // ionisation potential is a safe upper bound. Without the floor a zero
      // amplitude sets the cap to zero and silently freezes ion drift; without
      // the abs a negative amplitude makes it NaN.
      ion_speed_cap_(std::sqrt(2.0 * kE *
                               std::max(std::abs(config.rf_voltage),
                                        argon::kIonisationThreshold) /
                               argon::kMass)),
      ne_(config.num_cells, 0.0),
      ni_(config.num_cells, 0.0),
      neps_(config.num_cells, 0.0),
      ne_old_(config.num_cells, 0.0),
      ni_old_(config.num_cells, 0.0),
      neps_old_(config.num_cells, 0.0),
      phi_(config.num_cells, 0.0),
      efield_(config.num_cells + 1, 0.0),
      gamma_e_(config.num_cells + 1, 0.0),
      gamma_i_(config.num_cells + 1, 0.0),
      de_face_(config.num_cells + 1, 0.0),
      mue_face_(config.num_cells + 1, 0.0),
      sg_a_(config.num_cells + 1, 0.0),
      sg_b_(config.num_cells + 1, 0.0),
      eps_face_(config.num_cells + 1, constants::kVacuumPermittivity),
      jdiff_face_(config.num_cells + 1, 0.0),
      te_face_(config.num_cells + 1, 0.0),
      kiz_(config.num_cells, 0.0),
      kex_(config.num_cells, 0.0),
      kel_(config.num_cells, 0.0),
      matrix_(config.num_cells) {
    if (config.num_cells < 4) {
        throw std::invalid_argument("need at least 4 cells");
    }

    // Seed plasma: a cosine hump that tapers toward the electrodes, which is
    // closer to the eventual steady profile than a flat fill and so shortens
    // the transient by a few RF cycles.
    const std::size_t n = grid_.num_cells();
    for (std::size_t i = 0; i < n; ++i) {
        const double xi = grid_.cell_centre(i) / grid_.length();
        const double shape = std::sin(constants::kPi * xi);
        const double density =
            std::max(config_.density_floor, config_.initial_density * shape);
        ne_[i] = density;
        ni_[i] = density;
        neps_[i] = density * config_.initial_mean_energy;
    }
    ne_old_ = ne_;
    ni_old_ = ni_;
    neps_old_ = neps_;
}

double DischargeSolver::mean_energy(std::size_t i) const {
    const double ne = std::max(ne_[i], config_.density_floor);
    return std::clamp(neps_[i] / ne, config_.mean_energy_floor,
                      config_.mean_energy_ceiling);
}

double DischargeSolver::limited_drift(double mobility, double efield,
                                      double max_speed) {
    return std::clamp(mobility * efield, -max_speed, max_speed);
}

double DischargeSolver::temperature(std::size_t i) const {
    return (2.0 / 3.0) * mean_energy(i);
}

std::vector<double> DischargeSolver::electron_temperature() const {
    std::vector<double> te(grid_.num_cells());
    for (std::size_t i = 0; i < te.size(); ++i) te[i] = temperature(i);
    return te;
}

void DischargeSolver::seed_uniform(double ne, double ni, double mean_energy_ev) {
    std::fill(ne_.begin(), ne_.end(), ne);
    std::fill(ni_.begin(), ni_.end(), ni);
    std::fill(neps_.begin(), neps_.end(), ne * mean_energy_ev);
    ne_old_ = ne_;
    ni_old_ = ni_;
    neps_old_ = neps_;
}

void DischargeSolver::solve_field_explicit(double applied_voltage) {
    compute_electron_transport();
    std::fill(eps_face_.begin(), eps_face_.end(), kEps0);
    std::fill(jdiff_face_.begin(), jdiff_face_.end(), 0.0);

    const bool saved = config_.semi_implicit_poisson;
    config_.semi_implicit_poisson = false;
    solve_poisson(applied_voltage);
    compute_field(applied_voltage);
    config_.semi_implicit_poisson = saved;
}

// ---------------------------------------------------------------------------
// Poisson
// ---------------------------------------------------------------------------

/// Build the effective face permittivity for the semi-implicit Poisson solve.
///
/// Explicitly coupling Poisson to the continuity equations imposes the
/// dielectric relaxation limit dt < eps0 / (e mu_e n_e), which is about 45 ps
/// at 10^16 m^-3 and about 16 ps at the 2.8e16 m^-3 peak, below the 74 ps
/// default timestep. Following Ventzek, the charge density in
/// Poisson's equation is advanced one step first,
///
///     div(eps0 E^{n+1}) = rho^n - dt div(J^{n+1}),
///
/// and the drift part of J^{n+1} is written against the unknown field. That
/// moves the conductivity onto the left-hand side as an effective permittivity
///
///     eps_eff = eps0 + dt * e * (mu_i n_i + mu_e n_e),
///
/// leaving the same tridiagonal structure but removing the relaxation limit:
/// the scheme is now stable at the timestep the RF cycle actually dictates.
void DischargeSolver::compute_effective_permittivity() {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const double mu_i = gas_.ion_mobility();
    const double d_i = gas_.ion_diffusivity();

    for (std::size_t f = 0; f <= n; ++f) {
        const std::size_t left = (f == 0) ? 0 : f - 1;
        const std::size_t right = (f == n) ? n - 1 : f;

        const double ne_f = 0.5 * (ne_[left] + ne_[right]);
        const double ni_f = 0.5 * (ni_[left] + ni_[right]);

        // The conductivity must be the one the transport solve will actually
        // apply. Where the drift limiter binds, the electron flux responds to
        // the field with mobility cap/|E|, not mu_e; using the bare mobility
        // there overstates sigma, so Gauss's law is enforced against a current
        // the continuity equation never carries. The limiter state is taken
        // from the current field iterate, which is the linearisation point.
        const double e_mag = std::abs(efield_[f]);
        const double e_cap = config_.drift_speed_limit_factor *
                             ArgonGas::electron_thermal_speed(te_face_[f]);
        const double mu_e_eff = (mue_face_[f] * e_mag > e_cap) ? e_cap / e_mag
                                                               : mue_face_[f];
        const double mu_i_eff = (mu_i * e_mag > ion_speed_cap_) ? ion_speed_cap_ / e_mag
                                                                : mu_i;

        const double sigma = kE * (mu_i_eff * ni_f + mu_e_eff * ne_f);
        eps_face_[f] = config_.semi_implicit_poisson ? kEps0 + dt_ * sigma : kEps0;

        // Diffusive conduction current, evaluated at the current iterate like
        // the rest of J^{n+1}. Only the drift part was absorbed into eps_eff,
        // so this remainder stays on the right-hand side.
        if (f == 0 || f == n) {
            jdiff_face_[f] = 0.0;
        } else {
            const double grad_ne = (ne_[f] - ne_[f - 1]) / dx;
            const double grad_ni = (ni_[f] - ni_[f - 1]) / dx;
            const double gamma_e_diff = -de_face_[f] * grad_ne;
            const double gamma_i_diff = -d_i * grad_ni;
            jdiff_face_[f] = kE * (gamma_i_diff - gamma_e_diff);
        }
    }
}

void DischargeSolver::solve_poisson(double applied_voltage) {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const double dx2 = dx * dx;

    matrix_.reset();
    for (std::size_t i = 0; i < n; ++i) {
        const double eps_left = eps_face_[i] / dx2;
        const double eps_right = eps_face_[i + 1] / dx2;

        // The predictor advances rho^n, not the iterate: after the first sweep
        // the iterate already holds this step's transport, and starting from
        // it would count that transport twice.
        const bool predictor = config_.semi_implicit_poisson;
        double rho = predictor ? kE * (ni_old_[i] - ne_old_[i]) : kE * (ni_[i] - ne_[i]);
        if (predictor) {
            rho -= dt_ * (jdiff_face_[i + 1] - jdiff_face_[i]) / dx;
        }

        double diag = eps_left + eps_right;

        if (i > 0) {
            matrix_.lower(i) = -eps_left;
        } else {
            // Dirichlet electrodes sit on faces, half a cell outside the first
            // and last cell centres. Eliminating the ghost value
            // phi_{-1} = 2 V_wall - phi_0 doubles the boundary-face term and
            // moves the wall potential to the right-hand side.
            diag += eps_left;
            rho += 2.0 * eps_left * applied_voltage;
        }

        if (i + 1 < n) {
            matrix_.upper(i) = -eps_right;
        } else {
            diag += eps_right;
            rho += 2.0 * eps_right * 0.0;  // grounded electrode
        }

        matrix_.diag(i) = diag;
        matrix_.rhs(i) = rho;
    }

    matrix_.solve(phi_);
}

void DischargeSolver::compute_field(double applied_voltage) {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();

    // Interior faces: centred difference between neighbouring cell centres.
    for (std::size_t f = 1; f < n; ++f) {
        efield_[f] = -(phi_[f] - phi_[f - 1]) / dx;
    }
    // Wall faces sit half a cell from the adjacent centre.
    efield_[0] = -(phi_[0] - applied_voltage) / (0.5 * dx);
    efield_[n] = -(0.0 - phi_[n - 1]) / (0.5 * dx);
}

// ---------------------------------------------------------------------------
// Transport coefficients and fluxes
// ---------------------------------------------------------------------------

void DischargeSolver::compute_electron_transport() {
    const std::size_t n = grid_.num_cells();

    for (std::size_t f = 0; f <= n; ++f) {
        // Upwind-biased face temperature would be more robust in strong
        // gradients; the arithmetic mean is adequate here and keeps the
        // Jacobian symmetric between neighbours.
        double te;
        if (f == 0) {
            te = temperature(0);
        } else if (f == n) {
            te = temperature(n - 1);
        } else {
            te = 0.5 * (temperature(f - 1) + temperature(f));
        }
        te_face_[f] = te;
        mue_face_[f] = kinetics_.mobility(te);
        de_face_[f] = kinetics_.diffusivity(te);
    }
}

void DischargeSolver::compute_rates() {
    for (std::size_t i = 0; i < grid_.num_cells(); ++i) {
        const double te = temperature(i);
        kiz_[i] = kinetics_.ionisation_rate(te);
        kex_[i] = kinetics_.excitation_rate(te);
        kel_[i] = kinetics_.elastic_rate(te);
    }
}

void DischargeSolver::compute_electron_fluxes() {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();

    for (std::size_t f = 1; f < n; ++f) {
        // Electrons are negatively charged, so they drift against E.
        const double v = limited_drift(-mue_face_[f], efield_[f],
                                       config_.drift_speed_limit_factor *
                                           ArgonGas::electron_thermal_speed(te_face_[f]));
        const double d = de_face_[f];
        const double z = v * dx / d;
        const double a = (d / dx) * bernoulli(-z);
        const double b = (d / dx) * bernoulli(z);
        gamma_e_[f] = a * ne_[f - 1] - b * ne_[f];
    }

    const double vth0 = ArgonGas::electron_thermal_speed(temperature(0));
    const double vthn = ArgonGas::electron_thermal_speed(temperature(n - 1));
    gamma_e_[0] = -0.25 * vth0 * ne_[0] -
                  config_.secondary_emission_yield * wall_ion_flux_left();
    gamma_e_[n] = 0.25 * vthn * ne_[n - 1] -
                  config_.secondary_emission_yield * wall_ion_flux_right();
}

void DischargeSolver::compute_ion_fluxes() {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const double mu = gas_.ion_mobility();
    const double d = gas_.ion_diffusivity();

    for (std::size_t f = 1; f < n; ++f) {
        const double v = limited_drift(mu, efield_[f], ion_speed_cap_);
        const double z = v * dx / d;
        const double a = (d / dx) * bernoulli(-z);
        const double b = (d / dx) * bernoulli(z);
        gamma_i_[f] = a * ni_[f - 1] - b * ni_[f];
    }
    gamma_i_[0] = wall_ion_flux_left();
    gamma_i_[n] = wall_ion_flux_right();
}

/// Ions reaching an electrode are absorbed; none are emitted. The wall flux is
/// therefore one-way: thermal escape plus whatever drift is directed outward.
double DischargeSolver::wall_ion_flux_left() const {
    const double outward_drift = std::max(0.0, -gas_.ion_mobility() * efield_[0]);
    const double speed = 0.25 * gas_.ion_thermal_speed() + outward_drift;
    return -speed * ni_[0];
}

double DischargeSolver::wall_ion_flux_right() const {
    const std::size_t n = grid_.num_cells();
    const double outward_drift = std::max(0.0, gas_.ion_mobility() * efield_[n]);
    const double speed = 0.25 * gas_.ion_thermal_speed() + outward_drift;
    return speed * ni_[n - 1];
}

// ---------------------------------------------------------------------------
// Implicit species updates
// ---------------------------------------------------------------------------

void DischargeSolver::update_electrons() {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const double ngas = gas_.number_density();
    const double gamma_se = config_.secondary_emission_yield;

    // Scharfetter-Gummel weights on interior faces.
    for (std::size_t f = 1; f < n; ++f) {
        const double v = limited_drift(-mue_face_[f], efield_[f],
                                       config_.drift_speed_limit_factor *
                                           ArgonGas::electron_thermal_speed(te_face_[f]));
        const double d = de_face_[f];
        const double z = v * dx / d;
        sg_a_[f] = (d / dx) * bernoulli(-z);
        sg_b_[f] = (d / dx) * bernoulli(z);
    }

    matrix_.reset();
    for (std::size_t i = 0; i < n; ++i) {
        const double kiz = kiz_[i];

        double diag = 1.0 / dt_;
        double rhs = ne_old_[i] / dt_;

        // Ionisation is linear in n_e, so it goes on the diagonal. It lowers
        // the diagonal, but only by k_iz*N ~ 10^6 /s against 1/dt ~ 10^10 /s,
        // so diagonal dominance is never in danger at RF timesteps.
        diag -= kiz * ngas;

        if (i > 0) {
            matrix_.lower(i) = -sg_a_[i] / dx;
            diag += sg_b_[i] / dx;
        } else {
            // Powered electrode: thermal loss out, secondary electrons in.
            const double vth = ArgonGas::electron_thermal_speed(temperature(0));
            diag += 0.25 * vth / dx;
            rhs -= gamma_se * wall_ion_flux_left() / dx;
        }

        if (i + 1 < n) {
            matrix_.upper(i) = -sg_b_[i + 1] / dx;
            diag += sg_a_[i + 1] / dx;
        } else {
            const double vth = ArgonGas::electron_thermal_speed(temperature(n - 1));
            diag += 0.25 * vth / dx;
            rhs += gamma_se * wall_ion_flux_right() / dx;
        }

        matrix_.diag(i) = diag;
        matrix_.rhs(i) = rhs;
    }

    matrix_.solve(ne_);
    for (double& value : ne_) value = std::max(value, config_.density_floor);
}

void DischargeSolver::update_ions() {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const double ngas = gas_.number_density();
    const double mu = gas_.ion_mobility();
    const double d = gas_.ion_diffusivity();

    for (std::size_t f = 1; f < n; ++f) {
        const double v = limited_drift(mu, efield_[f], ion_speed_cap_);
        const double z = v * dx / d;
        sg_a_[f] = (d / dx) * bernoulli(-z);
        sg_b_[f] = (d / dx) * bernoulli(z);
    }

    matrix_.reset();
    for (std::size_t i = 0; i < n; ++i) {
        const double kiz = kiz_[i];

        double diag = 1.0 / dt_;
        // Ion production is proportional to n_e, not n_i, so it is an explicit
        // source here rather than a diagonal contribution.
        double rhs = ni_old_[i] / dt_ + kiz * ngas * ne_[i];

        if (i > 0) {
            matrix_.lower(i) = -sg_a_[i] / dx;
            diag += sg_b_[i] / dx;
        } else {
            const double outward = std::max(0.0, -mu * efield_[0]);
            diag += (0.25 * gas_.ion_thermal_speed() + outward) / dx;
        }

        if (i + 1 < n) {
            matrix_.upper(i) = -sg_b_[i + 1] / dx;
            diag += sg_a_[i + 1] / dx;
        } else {
            const double outward = std::max(0.0, mu * efield_[n]);
            diag += (0.25 * gas_.ion_thermal_speed() + outward) / dx;
        }

        matrix_.diag(i) = diag;
        matrix_.rhs(i) = rhs;
    }

    matrix_.solve(ni_);
    for (double& value : ni_) value = std::max(value, config_.density_floor);
}

void DischargeSolver::update_energy() {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const double ngas = gas_.number_density();
    const double gamma_se = config_.secondary_emission_yield;
    const double eps_se = config_.secondary_electron_energy;
    const double tgas = gas_.temperature_volts();
    const double mass_ratio = constants::kElectronMass / argon::kMass;

    // Energy is carried at 5/3 of the particle mobility and diffusivity.
    constexpr double kEnergyFactor = 5.0 / 3.0;
    for (std::size_t f = 1; f < n; ++f) {
        const double v = kEnergyFactor *
                         limited_drift(-mue_face_[f], efield_[f],
                                       config_.drift_speed_limit_factor *
                                           ArgonGas::electron_thermal_speed(te_face_[f]));
        const double d = kEnergyFactor * de_face_[f];
        const double z = v * dx / d;
        sg_a_[f] = (d / dx) * bernoulli(-z);
        sg_b_[f] = (d / dx) * bernoulli(z);
    }

    compute_electron_fluxes();

    matrix_.reset();
    for (std::size_t i = 0; i < n; ++i) {
        const double te = temperature(i);
        const double ne = ne_[i];

        // Ohmic heating, -Gamma_e . E, averaged from the two bounding faces.
        // Dividing by e is implicit: fluxes are in particles, energies in eV.
        const double heating =
            -0.5 * (gamma_e_[i] * efield_[i] + gamma_e_[i + 1] * efield_[i + 1]);

        // Collisional losses, all proportional to n_e.
        const double kiz = kiz_[i];
        const double kex = kex_[i];
        const double kel = kel_[i];
        const double inelastic = (kiz * argon::kIonisationThreshold +
                                  kex * argon::kExcitationThreshold) *
                                 ngas * ne;
        const double elastic =
            3.0 * mass_ratio * kel * ngas * ne * std::max(te - tgas, 0.0);
        const double loss = inelastic + elastic;

        double diag = 1.0 / dt_;
        double rhs = neps_old_[i] / dt_ + std::max(heating, 0.0);

        // Patankar linearisation: a sink written as (L / phi_old) * phi_new
        // raises the diagonal instead of subtracting from the source, which
        // keeps the energy density strictly positive however stiff the losses.
        const double neps_ref = std::max(neps_[i], config_.density_floor *
                                                       config_.mean_energy_floor);
        diag += loss / neps_ref;

        // Compression heating from a negative Ohmic term (possible transiently
        // in the sheath) is handled the same way, as a sink.
        if (heating < 0.0) diag += -heating / neps_ref;

        if (i > 0) {
            matrix_.lower(i) = -sg_a_[i] / dx;
            diag += sg_b_[i] / dx;
        } else {
            // Escaping electrons carry a mean energy of 2*Te, which for a
            // Maxwellian is (4/3) of the local mean energy: the one-quarter
            // thermal flux becomes a one-third coefficient on n_eps.
            const double vth = ArgonGas::electron_thermal_speed(temperature(0));
            diag += vth / (3.0 * dx);
            rhs -= gamma_se * wall_ion_flux_left() * eps_se / dx;
        }

        if (i + 1 < n) {
            matrix_.upper(i) = -sg_b_[i + 1] / dx;
            diag += sg_a_[i + 1] / dx;
        } else {
            const double vth = ArgonGas::electron_thermal_speed(temperature(n - 1));
            diag += vth / (3.0 * dx);
            rhs += gamma_se * wall_ion_flux_right() * eps_se / dx;
        }

        matrix_.diag(i) = diag;
        matrix_.rhs(i) = rhs;
    }

    matrix_.solve(neps_);
    for (std::size_t i = 0; i < n; ++i) {
        neps_[i] = std::clamp(neps_[i], ne_[i] * config_.mean_energy_floor,
                              ne_[i] * config_.mean_energy_ceiling);
    }
}

// ---------------------------------------------------------------------------
// Time integration
// ---------------------------------------------------------------------------

void DischargeSolver::step(double time) {
    ne_old_ = ne_;
    ni_old_ = ni_;
    neps_old_ = neps_;

    // Evaluate the electrode waveform at the end of the step, consistent with
    // the backward-Euler treatment of everything else.
    const double voltage = config_.rf_voltage * std::sin(omega_ * (time + dt_));

    for (int sweep = 0; sweep < config_.gummel_iterations; ++sweep) {
        compute_electron_transport();
        compute_rates();
        compute_effective_permittivity();
        solve_poisson(voltage);
        compute_field(voltage);
        update_electrons();
        update_ions();
        update_energy();
    }

    compute_electron_fluxes();
    compute_ion_fluxes();
}

CycleDiagnostics DischargeSolver::run(bool verbose) {
    const std::size_t n = grid_.num_cells();
    const double dx = grid_.dx();
    const int steps = config_.steps_per_cycle;
    const double ngas = gas_.number_density();
    const double inv_cells = 1.0 / static_cast<double>(n);
    const double inv_steps = 1.0 / static_cast<double>(std::max(1, steps));

    CycleDiagnostics diagnostics;
    double time = 0.0;
    double previous_mean_ne = 0.0;

    if (verbose) {
        std::cout << "  gas density   " << std::scientific << std::setprecision(3)
                  << ngas << " m^-3\n"
                  << "  electron mu   ";
        if (kinetics_.tabulated()) {
            std::cout << "from transport table\n";
        } else {
            std::cout << gas_.electron_mobility() << " m^2/(V s)\n";
        }
        std::cout
                  << "  timestep      " << dt_ << " s\n"
                  << "  cells         " << n << ", dx = " << dx << " m\n"
                  << "  Poisson       "
                  << (config_.semi_implicit_poisson ? "semi-implicit" : "explicit charge")
                  << "\n\n"
                  << "  cycle    mean n_e (m^-3)   peak eps (eV)   P_e (W/m^3)"
                  << "   P_tot (W/m^3)   cap duty   drift/cycle\n";
    }

    for (int cycle = 0; cycle < config_.num_cycles; ++cycle) {
        double pe_accum = 0.0;
        double pi_accum = 0.0;
        double ionisation_accum = 0.0;
        double ne_accum = 0.0;
        double peak_eps_cycle = 0.0;
        double cap_ratio_max = 0.0;
        long long cap_hits = 0;
        long long cap_evals = 0;

        for (int s = 0; s < steps; ++s) {
            step(time);
            time += dt_;

            double pe = 0.0;
            double pi = 0.0;
            double ionisation = 0.0;
            double ne_sum = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                // Cell-averaged Joule heating from the two bounding faces.
                // Electrons drift against E, so -Gamma_e.E is the work done
                // on them; ions drift along E, so +Gamma_i.E.
                pe += -0.5 * (gamma_e_[i] * efield_[i] + gamma_e_[i + 1] * efield_[i + 1]);
                pi += 0.5 * (gamma_i_[i] * efield_[i] + gamma_i_[i + 1] * efield_[i + 1]);
                ionisation += kinetics_.ionisation_rate(temperature(i)) * ngas * ne_[i];
                ne_sum += ne_[i];

                // Unclamped ratio: reaching the ceiling means the ceiling bound.
                const double raw = neps_[i] / std::max(ne_[i], config_.density_floor);
                peak_eps_cycle = std::max(peak_eps_cycle, raw);
            }
            pe_accum += pe * inv_cells;
            pi_accum += pi * inv_cells;
            ionisation_accum += ionisation * inv_cells;
            ne_accum += ne_sum * inv_cells;

            // Drift-cap activity on interior faces, against the end-of-step
            // field and the face temperatures the transport solve used.
            for (std::size_t f = 1; f < n; ++f) {
                const double cap = config_.drift_speed_limit_factor *
                                   ArgonGas::electron_thermal_speed(te_face_[f]);
                const double ratio = std::abs(mue_face_[f] * efield_[f]) / cap;
                cap_ratio_max = std::max(cap_ratio_max, ratio);
                if (ratio > 1.0) ++cap_hits;
                ++cap_evals;
            }
        }

        const double mean_ne = ne_accum * inv_steps;
        const double drift =
            (previous_mean_ne > 0.0)
                ? std::abs(mean_ne - previous_mean_ne) / previous_mean_ne
                : 1.0;
        previous_mean_ne = mean_ne;

        const double p_e = pe_accum * inv_steps * kE;
        const double p_i = pi_accum * inv_steps * kE;
        const double duty = cap_evals > 0 ? static_cast<double>(cap_hits) /
                                                static_cast<double>(cap_evals)
                                          : 0.0;

        double peak_eps = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            peak_eps = std::max(peak_eps, mean_energy(i));

        const bool final_cycle = (cycle == config_.num_cycles - 1);
        if (verbose && (cycle % 20 == 0 || final_cycle)) {
            std::cout << "  " << std::setw(5) << cycle << "    " << std::scientific
                      << std::setprecision(4) << mean_ne << "        " << std::fixed
                      << std::setprecision(2) << std::setw(6) << peak_eps << "      "
                      << std::scientific << std::setprecision(3) << p_e << "     "
                      << p_e + p_i << "      " << std::fixed << std::setprecision(3)
                      << duty << "      " << std::scientific << std::setprecision(3)
                      << drift << "\n";
        }

        if (final_cycle) {
            diagnostics.electron_power_density = p_e;
            diagnostics.ion_power_density = p_i;
            diagnostics.total_power_density = p_e + p_i;
            diagnostics.ionisation_rate = ionisation_accum * inv_steps;
            diagnostics.mean_electron_density = mean_ne;
            diagnostics.peak_electron_density =
                *std::max_element(ne_.begin(), ne_.end());
            diagnostics.peak_ion_density = *std::max_element(ni_.begin(), ni_.end());
            diagnostics.peak_mean_energy = peak_eps;
            diagnostics.peak_mean_energy_in_cycle = peak_eps_cycle;
            diagnostics.drift_cap_duty = duty;
            diagnostics.drift_cap_max_ratio = cap_ratio_max;
            diagnostics.cycle_drift = drift;
            diagnostics.cycles_completed = cycle + 1;
        }
    }

    return diagnostics;
}

void DischargeSolver::write_profiles(const std::string& path) const {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open " + path);

    out << "x_m,n_e_m3,n_i_m3,mean_energy_eV,Te_V,potential_V,charge_density_Cm3\n";
    out << std::scientific << std::setprecision(8);
    for (std::size_t i = 0; i < grid_.num_cells(); ++i) {
        out << grid_.cell_centre(i) << ',' << ne_[i] << ',' << ni_[i] << ','
            << mean_energy(i) << ',' << temperature(i) << ',' << phi_[i] << ','
            << kE * (ni_[i] - ne_[i]) << '\n';
    }
}

}  // namespace ccp

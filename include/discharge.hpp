#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "argon.hpp"
#include "grid.hpp"
#include "tridiagonal.hpp"

namespace ccp {

/// Everything a run needs, with defaults describing a fairly standard
/// capacitively-coupled argon discharge in the GEC reference cell regime.
struct DischargeConfig {
    std::size_t num_cells = 200;
    double gap_length = 0.02;            // m, electrode separation
    /// 1 Torr. Below roughly 300 mTorr in a 20 mm gap the electron mean free
    /// path approaches the gap width and the drift-diffusion closure stops
    /// being defensible -- see the validity note in the README.
    double pressure_mtorr = 1000.0;
    double gas_temperature_k = 300.0;

    double rf_voltage = 300.0;           // V, amplitude on the powered electrode
    double rf_frequency = 13.56e6;       // Hz
    double secondary_emission_yield = 0.05;
    double secondary_electron_energy = 5.0;  // eV

    double initial_density = 1.0e15;     // m^-3, uniform seed plasma
    double initial_mean_energy = 4.0;    // eV

    int steps_per_cycle = 1000;
    int num_cycles = 400;
    int gummel_iterations = 3;

    double density_floor = 1.0e6;        // m^-3
    double mean_energy_floor = 0.05;     // eV

    /// Ceiling on the local mean electron energy. In the sheath the electron
    /// density collapses toward the floor while the energy density decays more
    /// slowly, so the ratio n_eps/n_e can run away numerically even though no
    /// real electron population is at that energy. The ceiling bounds that
    /// ratio; it is a numerical guard, not a physical statement.
    double mean_energy_ceiling = 60.0;   // eV

    /// Constant-mobility drift-diffusion has no high-field saturation, so a
    /// 200 V sheath would drive electrons to a sizeable fraction of c. Real
    /// codes solve this with an E/N-dependent mobility table from a Boltzmann
    /// solver; lacking one, the drift speed is capped at this multiple of the
    /// local electron thermal speed, which is where the fluid closure fails
    /// anyway.
    double drift_speed_limit_factor = 1.0;

    /// Use the semi-implicit Poisson formulation, which folds the plasma
    /// conductivity into an effective permittivity and removes the dielectric
    /// relaxation timestep restriction.
    bool semi_implicit_poisson = true;
};

/// Cycle-averaged quantities, accumulated over the final RF period.
struct CycleDiagnostics {
    double peak_electron_density = 0.0;   // m^-3
    double peak_ion_density = 0.0;        // m^-3
    double mean_electron_density = 0.0;   // m^-3
    double peak_mean_energy = 0.0;        // eV
    double absorbed_power_density = 0.0;  // W/m^3, gap-averaged
    double ionisation_rate = 0.0;         // m^-3 s^-1, gap-averaged

    /// Relative change in gap-averaged electron density over the last RF
    /// cycle. The discharge is in periodic steady state when this approaches
    /// zero; it decays slowly because the equilibration timescale is set by
    /// ion transit, which is many RF periods.
    double cycle_drift = 0.0;
    int cycles_completed = 0;
};

/// 1D drift-diffusion model of a capacitively-coupled RF discharge.
///
/// State variables, all at cell centres:
///   n_e      electron density                 [m^-3]
///   n_i      singly-charged ion density       [m^-3]
///   n_eps    electron energy density          [eV m^-3]
///
/// The local mean energy approximation closes the system: the electron
/// temperature is Te = (2/3) * n_eps / n_e, and every rate coefficient and
/// diffusivity is evaluated at that local Te.
///
/// Each timestep runs a Gummel (sequential Picard) iteration: Poisson, then
/// electrons, then ions, then energy, repeated a few times so the strong
/// nonlinear coupling between space charge and field is converged before the
/// step is accepted.
class DischargeSolver {
public:
    explicit DischargeSolver(DischargeConfig config);

    /// Advance the discharge for the configured number of RF cycles.
    /// Returns diagnostics averaged over the final cycle.
    CycleDiagnostics run(bool verbose = true);

    /// Advance a single timestep. `time` is the time at the start of the step.
    void step(double time);

    double timestep() const noexcept { return dt_; }
    const Grid& grid() const noexcept { return grid_; }

    const std::vector<double>& electron_density() const noexcept { return ne_; }
    const std::vector<double>& ion_density() const noexcept { return ni_; }
    const std::vector<double>& potential() const noexcept { return phi_; }
    const std::vector<double>& electric_field() const noexcept { return efield_; }

    /// Local electron temperature in volts, cell by cell.
    std::vector<double> electron_temperature() const;

    /// Write cell-centred profiles to a CSV file.
    void write_profiles(const std::string& path) const;

private:
    // --- stages of a single Gummel sweep -------------------------------
    void solve_poisson(double applied_voltage);
    void compute_field(double applied_voltage);
    void update_electrons();
    void update_ions();
    void update_energy();

    // --- helpers -------------------------------------------------------
    void compute_electron_transport();
    void compute_electron_fluxes();
    void compute_ion_fluxes();
    void compute_effective_permittivity();
    double wall_ion_flux_left() const;
    double wall_ion_flux_right() const;
    double mean_energy(std::size_t i) const;
    double temperature(std::size_t i) const;

    /// Drift velocity mobility*E, clamped to +/- max_speed. Electrons are
    /// capped at their thermal speed, ions at free-fall through the full
    /// electrode voltage -- the two points where a mobility model stops being
    /// meaningful for each species.
    static double limited_drift(double mobility, double efield, double max_speed);

    DischargeConfig config_;
    Grid grid_;
    ArgonGas gas_;
    double dt_;
    double omega_;
    double ion_speed_cap_;

    // State at the new time level (the current Gummel iterate).
    std::vector<double> ne_, ni_, neps_;
    // State at the previous time level.
    std::vector<double> ne_old_, ni_old_, neps_old_;

    // Field quantities.
    std::vector<double> phi_;      // cells
    std::vector<double> efield_;   // faces
    std::vector<double> gamma_e_;  // faces
    std::vector<double> gamma_i_;  // faces

    // Face-centred electron transport coefficients, refreshed each sweep.
    std::vector<double> de_face_, mue_face_;

    // Scharfetter-Gummel face weights, reused by every implicit transport
    // solve. Held as members so the inner loop never allocates: a run is
    // ~10^5 steps x 3 Gummel sweeps x 3 species.
    std::vector<double> sg_a_, sg_b_;

    // Semi-implicit Poisson: face permittivity eps0 + dt*sigma, and the
    // diffusive part of the conduction current that its right-hand side needs.
    std::vector<double> eps_face_, jdiff_face_, te_face_;

    TridiagonalSolver matrix_;
};

}  // namespace ccp

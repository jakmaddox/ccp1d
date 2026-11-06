#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "argon.hpp"
#include "grid.hpp"
#include "transport.hpp"
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
    /// ratio; it is a numerical guard, not a physical statement. It has not
    /// been seen to bind at 100 mTorr or 1 Torr (peak ~20 eV); compare with
    /// CycleDiagnostics::peak_mean_energy_in_cycle.
    double mean_energy_ceiling = 60.0;   // eV

    /// Constant-mobility drift-diffusion has no high-field saturation, so in
    /// the sheath mu_e E reaches several times the electron thermal speed
    /// (about 5x at the defaults, ~1e7 m/s). A Boltzmann table in
    /// `transport_table` carries the mobility's energy dependence; either way
    /// the drift speed is capped at this multiple of the local thermal speed,
    /// which is where the fluid closure fails anyway.
    ///
    /// This is not a passive guard. At the defaults it binds on about 28% of
    /// interior faces over a cycle and removing it changes the absorbed power
    /// by roughly 25%, so the capped region is modelled as free-streaming at
    /// v_th. CycleDiagnostics::drift_cap_duty reports how much it is doing.
    double drift_speed_limit_factor = 1.0;

    /// Electron transport and rates against mean energy, usually from a
    /// Boltzmann solver. Null keeps the built-in constant mobility, Einstein
    /// relation and Lieberman fits. The drift cap and the 5/3 energy transport
    /// factor apply either way. Ion transport does not come from this table.
    std::shared_ptr<const TransportTable> transport_table;

    /// Use the semi-implicit Poisson formulation, which folds the plasma
    /// conductivity into an effective permittivity and removes the dielectric
    /// relaxation timestep restriction.
    bool semi_implicit_poisson = true;
};

/// Diagnostics from the final RF period.
///
/// Note which fields are cycle-averaged and which are instantaneous: the
/// `peak_*_density` values and `peak_mean_energy` are read off the state at
/// the end of the last step, while the power, ionisation rate and
/// `mean_electron_density` are averaged over the cycle. The three are
/// different numbers and were previously reported under one heading.
struct CycleDiagnostics {
    double peak_electron_density = 0.0;   // m^-3, end-of-cycle snapshot
    double peak_ion_density = 0.0;        // m^-3, end-of-cycle snapshot
    double mean_electron_density = 0.0;   // m^-3, cycle- and gap-averaged
    double peak_mean_energy = 0.0;        // eV, end-of-cycle snapshot

    /// Largest local mean energy seen at any point during the cycle. This is
    /// the number to compare against `mean_energy_ceiling`: the end-of-cycle
    /// snapshot can sit well below a ceiling that binds mid-cycle.
    double peak_mean_energy_in_cycle = 0.0;  // eV

    /// Joule heating, gap-averaged over the cycle. The electron term is the
    /// one that drives the energy equation; the ion term is comparable in
    /// magnitude at these pressures and belongs in any statement about the
    /// power the discharge absorbs.
    double electron_power_density = 0.0;  // W/m^3
    double ion_power_density = 0.0;       // W/m^3
    double total_power_density = 0.0;     // W/m^3

    double ionisation_rate = 0.0;         // m^-3 s^-1, gap-averaged

    /// Fraction of interior-face drift evaluations over the cycle at which
    /// |mu_e E| exceeded the cap and the limiter was active. A non-trivial
    /// value means the electron drift in those cells is set by the cap and
    /// not by the mobility model: see the validity note in the README.
    double drift_cap_duty = 0.0;

    /// Largest |mu_e E| / v_cap ratio seen over the cycle. 1.0 means the cap
    /// was just touched; larger means it was actively binding.
    double drift_cap_max_ratio = 0.0;

    /// Relative change in gap-averaged electron density over the last RF
    /// cycle. The discharge is in periodic steady state when this approaches
    /// zero; it decays slowly because the equilibration timescale is set by
    /// ambipolar ion transport across the gap, which is a few thousand RF
    /// periods at the default settings.
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

    /// Throw std::runtime_error naming `cycle` if the state is not finite. A
    /// NaN spreads through every field within a step, and run() calls this
    /// after each cycle; a caller stepping by hand should do the same.
    void require_finite(int cycle) const;

    double timestep() const noexcept { return dt_; }
    const Grid& grid() const noexcept { return grid_; }

    const std::vector<double>& electron_density() const noexcept { return ne_; }
    const std::vector<double>& ion_density() const noexcept { return ni_; }
    const std::vector<double>& potential() const noexcept { return phi_; }
    const std::vector<double>& electric_field() const noexcept { return efield_; }

    /// Electron energy density [eV m^-3], cell-centred.
    const std::vector<double>& energy_density() const noexcept { return neps_; }

    /// Face-centred particle fluxes [m^-2 s^-1] as of the end of the last
    /// step. Faces 0 and N carry the wall fluxes, secondary emission included.
    const std::vector<double>& electron_flux() const noexcept { return gamma_e_; }
    const std::vector<double>& ion_flux() const noexcept { return gamma_i_; }

    const DischargeConfig& config() const noexcept { return config_; }
    const ArgonGas& gas() const noexcept { return gas_; }
    const ElectronKinetics& kinetics() const noexcept { return kinetics_; }

    /// Local electron temperature in volts, cell by cell.
    std::vector<double> electron_temperature() const;

    /// Local mean electron energy in eV for cell i, floor and ceiling applied.
    double mean_energy_at(std::size_t i) const { return mean_energy(i); }

    /// Replace the state with a uniform fill. The constructor's cosine seed is
    /// one reasonable initial condition; this is another, and it is the one an
    /// analytic verification of the Poisson operator needs.
    void seed_uniform(double ne, double ni, double mean_energy_ev);

    /// Solve the field for the current state without advancing any species,
    /// with the plain (eps0) Poisson operator rather than the semi-implicit
    /// one. This is what lets the discrete Laplacian and the half-cell
    /// Dirichlet elimination be checked against a known charge distribution.
    void solve_field_explicit(double applied_voltage);

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
    void compute_rates();
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
    ElectronKinetics kinetics_;
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

    // Rate coefficients [m^3/s], frozen for the duration of one Gummel sweep.
    // Evaluating them once per sweep, before any species moves, means the
    // electron source, the ion source and the energy cost of ionisation all
    // see the same k_iz, so ionisation creates electrons and ions in equal
    // numbers and charges each event exactly once.
    std::vector<double> kiz_, kex_, kel_;

    TridiagonalSolver matrix_;
};

}  // namespace ccp

#pragma once

#include <algorithm>
#include <cmath>

namespace ccp {

/// Physical constants (CODATA 2018).
namespace constants {
inline constexpr double kElementaryCharge = 1.602176634e-19;   // C
inline constexpr double kVacuumPermittivity = 8.8541878128e-12; // F/m
inline constexpr double kElectronMass = 9.1093837015e-31;      // kg
inline constexpr double kBoltzmann = 1.380649e-23;             // J/K
inline constexpr double kAtomicMassUnit = 1.66053906660e-27;   // kg
inline constexpr double kPi = 3.14159265358979323846;
}  // namespace constants

/// Electron-impact reaction thresholds for argon, in eV.
namespace argon {
inline constexpr double kIonisationThreshold = 15.76;
inline constexpr double kExcitationThreshold = 11.55;
inline constexpr double kMass = 39.948 * constants::kAtomicMassUnit;  // kg

/// Reduced mobilities, quoted as mu*N so they scale correctly with pressure.
inline constexpr double kReducedElectronMobility = 4.0e24;   // 1/(m V s)
inline constexpr double kReducedIonMobility = 4.11e21;       // 1/(m V s)

/// Rate-coefficient fits from Lieberman & Lichtenberg, "Principles of Plasma
/// Discharges and Materials Processing". Electron temperature Te is in volts,
/// results in m^3/s. The fits are calibrated for roughly 1-7 V; the clamp below
/// keeps the extrapolation bounded when a cell transiently runs hot in the
/// sheath, which would otherwise feed an unphysical ionisation spike back into
/// the energy equation.
inline double clamp_temperature(double te_volts) noexcept {
    return std::clamp(te_volts, 0.1, 60.0);
}

inline double ionisation_rate(double te_volts) noexcept {
    const double te = clamp_temperature(te_volts);
    return 2.34e-14 * std::pow(te, 0.59) * std::exp(-17.44 / te);
}

inline double excitation_rate(double te_volts) noexcept {
    const double te = clamp_temperature(te_volts);
    return 2.48e-14 * std::pow(te, 0.33) * std::exp(-12.78 / te);
}

inline double elastic_rate(double te_volts) noexcept {
    const double te = clamp_temperature(te_volts);
    const double l = std::log(te);
    return 2.336e-14 * std::pow(te, 1.609) *
           std::exp(0.0618 * l * l - 0.1171 * l * l * l);
}
}  // namespace argon

/// Background gas state and the transport coefficients derived from it.
///
/// Everything that depends on pressure is resolved once, here, so the solver
/// only ever sees absolute mobilities and diffusivities.
class ArgonGas {
public:
    /// @param pressure_mtorr  neutral gas pressure in mTorr
    /// @param temperature_k   neutral gas temperature in kelvin
    ArgonGas(double pressure_mtorr, double temperature_k)
        : temperature_k_(temperature_k),
          temperature_volts_(constants::kBoltzmann * temperature_k /
                             constants::kElementaryCharge),
          number_density_(pressure_mtorr * 0.1333224 /
                          (constants::kBoltzmann * temperature_k)),
          electron_mobility_(argon::kReducedElectronMobility / number_density_),
          ion_mobility_(argon::kReducedIonMobility / number_density_) {}

    double number_density() const noexcept { return number_density_; }
    double temperature_volts() const noexcept { return temperature_volts_; }

    double electron_mobility() const noexcept { return electron_mobility_; }
    double ion_mobility() const noexcept { return ion_mobility_; }

    /// Einstein relation, using the local electron temperature in volts.
    double electron_diffusivity(double te_volts) const noexcept {
        return electron_mobility_ * te_volts;
    }

    /// Ions are assumed thermalised with the neutral background.
    double ion_diffusivity() const noexcept {
        return ion_mobility_ * temperature_volts_;
    }

    /// Mean thermal speed of electrons at temperature te_volts, in m/s.
    static double electron_thermal_speed(double te_volts) noexcept {
        const double te = std::max(te_volts, 1.0e-3);
        return std::sqrt(8.0 * constants::kElementaryCharge * te /
                         (constants::kPi * constants::kElectronMass));
    }

    /// Mean thermal speed of argon ions at the gas temperature, in m/s.
    double ion_thermal_speed() const noexcept {
        return std::sqrt(8.0 * constants::kElementaryCharge * temperature_volts_ /
                         (constants::kPi * argon::kMass));
    }

private:
    double temperature_k_;
    double temperature_volts_;
    double number_density_;
    double electron_mobility_;
    double ion_mobility_;
};

}  // namespace ccp

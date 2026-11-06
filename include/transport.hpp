#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "argon.hpp"

namespace ccp {

/// Electron swarm data tabulated against local mean energy.
///
/// The usual source is a Boltzmann solver run over a range of E/N, with mean
/// energy as the common axis so that the table can close a local mean energy
/// model. See tools/make_transport_table.py.
///
/// Mobility and diffusivity are stored reduced (times N), so one table holds
/// at every pressure. Queries interpolate linearly in log energy and log value
/// and hold the end values outside the table, which is the same behaviour as
/// the clamp on the analytic fits.
class TransportTable {
public:
    struct Row {
        double mean_energy;    // eV
        double mobility_n;     // mu*N, 1/(m V s)
        double diffusivity_n;  // D*N, 1/(m s)
        double k_elastic;      // momentum transfer, m^3/s
        double k_excitation;   // m^3/s, weighted to cost kExcitationThreshold
        double k_ionisation;   // m^3/s
    };

    /// @throws std::invalid_argument unless energies strictly increase and
    /// every coefficient is finite, with positive transport and rates >= 0.
    explicit TransportTable(std::vector<Row> rows);

    /// Read whitespace-separated columns in Row order. Blank lines and lines
    /// starting with '#' are skipped.
    static TransportTable load(const std::string& path);

    /// The built-in closure in table form: constant mu*N, the Einstein
    /// relation and the Lieberman fits. It exists so the table path can be
    /// checked against the analytic path.
    static TransportTable from_analytic(std::size_t points = 2000);

    void save(const std::string& path) const;

    const std::vector<Row>& rows() const noexcept { return rows_; }

    double mobility_n(double mean_energy) const { return value(0, mean_energy); }
    double diffusivity_n(double mean_energy) const { return value(1, mean_energy); }
    double k_elastic(double mean_energy) const { return value(2, mean_energy); }
    double k_excitation(double mean_energy) const { return value(3, mean_energy); }
    double k_ionisation(double mean_energy) const { return value(4, mean_energy); }

private:
    double value(std::size_t column, double mean_energy) const;

    std::vector<Row> rows_;
    std::vector<double> log_energy_;
    std::vector<std::array<double, 5>> log_values_;
};

/// The electron closure the solver and its checks read: every mobility,
/// diffusivity and rate coefficient comes from here, so a balance can never
/// be measured against different rates from the ones the solver used.
///
/// Without a table this is the built-in model: constant mu*N, the Einstein
/// relation and the Lieberman fits. Arguments are Te in volts, the variable
/// the solver carries; tables are read at the mean energy 3/2 Te.
class ElectronKinetics {
public:
    ElectronKinetics(const ArgonGas& gas, std::shared_ptr<const TransportTable> table)
        : table_(std::move(table)),
          number_density_(gas.number_density()),
          constant_mobility_(gas.electron_mobility()) {}

    bool tabulated() const noexcept { return table_ != nullptr; }

    double mobility(double te_volts) const {
        return table_ ? table_->mobility_n(1.5 * te_volts) / number_density_
                      : constant_mobility_;
    }
    double diffusivity(double te_volts) const {
        return table_ ? table_->diffusivity_n(1.5 * te_volts) / number_density_
                      : constant_mobility_ * te_volts;
    }
    double ionisation_rate(double te_volts) const {
        return table_ ? table_->k_ionisation(1.5 * te_volts)
                      : argon::ionisation_rate(te_volts);
    }
    double excitation_rate(double te_volts) const {
        return table_ ? table_->k_excitation(1.5 * te_volts)
                      : argon::excitation_rate(te_volts);
    }
    double elastic_rate(double te_volts) const {
        return table_ ? table_->k_elastic(1.5 * te_volts) : argon::elastic_rate(te_volts);
    }

private:
    std::shared_ptr<const TransportTable> table_;
    double number_density_;
    double constant_mobility_;
};

}  // namespace ccp

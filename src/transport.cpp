#include "transport.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ccp {

namespace {
// Rates are zero below threshold in a Boltzmann table. This stands in for zero
// on the log axis and is far below any rate that matters.
constexpr double kSmallestValue = 1.0e-300;

double log_value(double value) { return std::log(std::max(value, kSmallestValue)); }
}  // namespace

TransportTable::TransportTable(std::vector<Row> rows) : rows_(std::move(rows)) {
    if (rows_.size() < 2) {
        throw std::invalid_argument("transport table needs at least two rows");
    }
    log_energy_.reserve(rows_.size());
    log_values_.reserve(rows_.size());
    for (std::size_t i = 0; i < rows_.size(); ++i) {
        const Row& r = rows_[i];
        const bool finite = std::isfinite(r.mean_energy) && std::isfinite(r.mobility_n) &&
                            std::isfinite(r.diffusivity_n) && std::isfinite(r.k_elastic) &&
                            std::isfinite(r.k_excitation) && std::isfinite(r.k_ionisation);
        if (!finite || r.mean_energy <= 0.0 || r.mobility_n <= 0.0 ||
            r.diffusivity_n <= 0.0 || r.k_elastic < 0.0 || r.k_excitation < 0.0 ||
            r.k_ionisation < 0.0) {
            throw std::invalid_argument("transport table row " + std::to_string(i) +
                                        " has an invalid value");
        }
        if (i > 0 && r.mean_energy <= rows_[i - 1].mean_energy) {
            throw std::invalid_argument("transport table energies must strictly increase");
        }
        log_energy_.push_back(std::log(r.mean_energy));
        log_values_.push_back({log_value(r.mobility_n), log_value(r.diffusivity_n),
                               log_value(r.k_elastic), log_value(r.k_excitation),
                               log_value(r.k_ionisation)});
    }
}

TransportTable TransportTable::load(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open transport table: " + path);

    std::vector<Row> rows;
    std::string line;
    for (int number = 1; std::getline(in, line); ++number) {
        const auto first = line.find_first_not_of(" \t\r");
        if (first == std::string::npos || line[first] == '#') continue;

        std::istringstream fields(line);
        Row r{};
        std::string extra;
        if (!(fields >> r.mean_energy >> r.mobility_n >> r.diffusivity_n >> r.k_elastic >>
              r.k_excitation >> r.k_ionisation) ||
            (fields >> extra)) {
            throw std::runtime_error(path + ":" + std::to_string(number) +
                                     ": expected six numeric columns");
        }
        rows.push_back(r);
    }
    return TransportTable(std::move(rows));
}

TransportTable TransportTable::from_analytic(std::size_t points) {
    if (points < 2) throw std::invalid_argument("analytic table needs two points");
    // The span of argon::clamp_temperature, so holding the end values
    // reproduces the clamp on the fits.
    const double low = std::log(0.1);
    const double high = std::log(60.0);

    std::vector<Row> rows(points);
    for (std::size_t i = 0; i < points; ++i) {
        const double s = static_cast<double>(i) / static_cast<double>(points - 1);
        const double te = std::exp(low + s * (high - low));
        rows[i] = {1.5 * te,
                   argon::kReducedElectronMobility,
                   argon::kReducedElectronMobility * te,
                   argon::elastic_rate(te),
                   argon::excitation_rate(te),
                   argon::ionisation_rate(te)};
    }
    return TransportTable(std::move(rows));
}

void TransportTable::save(const std::string& path) const {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write transport table: " + path);
    out << "# mean_energy_eV mobility_N_per_m_V_s diffusivity_N_per_m_s "
           "k_elastic_m3_s k_excitation_m3_s k_ionisation_m3_s\n"
        << std::scientific << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const Row& r : rows_) {
        out << r.mean_energy << ' ' << r.mobility_n << ' ' << r.diffusivity_n << ' '
            << r.k_elastic << ' ' << r.k_excitation << ' ' << r.k_ionisation << '\n';
    }
    if (!out) throw std::runtime_error("failed writing transport table: " + path);
}

double TransportTable::value(std::size_t column, double mean_energy) const {
    if (mean_energy <= rows_.front().mean_energy) return std::exp(log_values_.front()[column]);
    if (mean_energy >= rows_.back().mean_energy) return std::exp(log_values_.back()[column]);

    const double x = std::log(mean_energy);
    const auto upper = std::upper_bound(log_energy_.begin(), log_energy_.end(), x);
    const auto i = static_cast<std::size_t>(upper - log_energy_.begin());
    const double w = (x - log_energy_[i - 1]) / (log_energy_[i] - log_energy_[i - 1]);
    const double a = log_values_[i - 1][column];
    const double b = log_values_[i][column];
    return std::exp(a + w * (b - a));
}

}  // namespace ccp

// Cycle-averaged mid-gap values, for comparison with published cases.
//
// The profiles CSV is a snapshot at the end of the last step. Published
// discharge results and probe measurements are cycle averages, so comparing
// the two needs this instead.
//
//   ccp1d_benchmark --voltage 75 --yield 0.1 --emission-energy 0.75
//                   --cycles 3000 --transport-table argon.table

#include <cstdio>
#include <cstdlib>
#include <exception>
#include <memory>
#include <string>

#include "discharge.hpp"

int main(int argc, char** argv) {
    ccp::DischargeConfig config;
    config.num_cycles = 3000;
    try {
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            const bool v = i + 1 < argc;
            if (a == "--pressure" && v) config.pressure_mtorr = std::stod(argv[++i]);
            else if (a == "--voltage" && v) config.rf_voltage = std::stod(argv[++i]);
            else if (a == "--gap" && v) config.gap_length = std::stod(argv[++i]);
            else if (a == "--cells" && v) config.num_cells = std::stoul(argv[++i]);
            else if (a == "--steps" && v) config.steps_per_cycle = std::stoi(argv[++i]);
            else if (a == "--cycles" && v) config.num_cycles = std::stoi(argv[++i]);
            else if (a == "--yield" && v) config.secondary_emission_yield = std::stod(argv[++i]);
            else if (a == "--emission-energy" && v)
                config.secondary_electron_energy = std::stod(argv[++i]);
            else if (a == "--transport-table" && v) {
                config.transport_table = std::make_shared<const ccp::TransportTable>(
                    ccp::TransportTable::load(argv[++i]));
            } else {
                std::fprintf(stderr,
                             "usage: %s [--pressure mTorr] [--voltage amplitude] [--gap m]\n"
                             "          [--cells n] [--steps n] [--cycles n] [--yield g]\n"
                             "          [--emission-energy eV] [--transport-table path]\n",
                             argv[0]);
                return EXIT_FAILURE;
            }
        }

        const int final_cycle = config.num_cycles;
        config.num_cycles -= 1;
        ccp::DischargeSolver solver(config);
        solver.run(false);

        double time = static_cast<double>(config.num_cycles) /
                      config.rf_frequency;
        const std::size_t right = solver.grid().num_cells() / 2;
        const std::size_t left = right - 1;
        const int steps = config.steps_per_cycle;
        double ne = 0.0;
        double ni = 0.0;
        double energy = 0.0;
        double potential = 0.0;
        for (int s = 0; s < steps; ++s) {
            solver.step(time);
            time += solver.timestep();
            ne += 0.5 * (solver.electron_density()[left] + solver.electron_density()[right]);
            ni += 0.5 * (solver.ion_density()[left] + solver.ion_density()[right]);
            energy += 0.5 * (solver.mean_energy_at(left) + solver.mean_energy_at(right));
            potential += 0.5 * (solver.potential()[left] + solver.potential()[right]);
        }
        solver.require_finite(config.num_cycles);
        const double n = static_cast<double>(steps);
        std::printf("after %d cycles, mid-gap cycle averages\n", final_cycle);
        std::printf("  electron density   %.4e m^-3\n", ne / n);
        std::printf("  ion density        %.4e m^-3\n", ni / n);
        std::printf("  mean energy        %.3f eV\n", energy / n);
        std::printf("  Te (2/3 energy)    %.3f V\n", (2.0 / 3.0) * energy / n);
        std::printf("  plasma potential   %.2f V\n", potential / n);
    } catch (const std::exception& error) {
        std::fprintf(stderr, "error: %s\n", error.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

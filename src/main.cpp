#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "discharge.hpp"

namespace {

void print_usage(const char* program) {
    std::cout
        << "usage: " << program << " [options]\n\n"
        << "  --pressure <mTorr>    neutral argon pressure       (default 100)\n"
        << "  --voltage <V>         RF amplitude                 (default 200)\n"
        << "  --frequency <Hz>      drive frequency              (default 13.56e6)\n"
        << "  --gap <m>             electrode separation         (default 0.02)\n"
        << "  --cells <n>           grid cells                   (default 200)\n"
        << "  --cycles <n>          RF cycles to simulate        (default 120)\n"
        << "  --steps <n>           timesteps per RF cycle       (default 2000)\n"
        << "  --gummel <n>          Gummel sweeps per step       (default 3)\n"
        << "  --output <path>       profile CSV                  (default profiles.csv)\n"
        << "  --explicit-poisson    disable the semi-implicit Poisson solve\n"
        << "                        (expect divergence at the default timestep)\n"
        << "  --quiet               suppress per-cycle output\n"
        << "  --help                this message\n";
}

bool matches(const std::string& arg, const char* name) { return arg == name; }

}  // namespace

int main(int argc, char** argv) {
    ccp::DischargeConfig config;
    std::string output_path = "profiles.csv";
    bool verbose = true;

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            const bool has_value = (i + 1 < argc);

            if (matches(arg, "--help")) {
                print_usage(argv[0]);
                return EXIT_SUCCESS;
            } else if (matches(arg, "--quiet")) {
                verbose = false;
            } else if (matches(arg, "--explicit-poisson")) {
                config.semi_implicit_poisson = false;
            } else if (matches(arg, "--pressure") && has_value) {
                config.pressure_mtorr = std::stod(argv[++i]);
            } else if (matches(arg, "--voltage") && has_value) {
                config.rf_voltage = std::stod(argv[++i]);
            } else if (matches(arg, "--frequency") && has_value) {
                config.rf_frequency = std::stod(argv[++i]);
            } else if (matches(arg, "--gap") && has_value) {
                config.gap_length = std::stod(argv[++i]);
            } else if (matches(arg, "--cells") && has_value) {
                config.num_cells = static_cast<std::size_t>(std::stoul(argv[++i]));
            } else if (matches(arg, "--cycles") && has_value) {
                config.num_cycles = std::stoi(argv[++i]);
            } else if (matches(arg, "--steps") && has_value) {
                config.steps_per_cycle = std::stoi(argv[++i]);
            } else if (matches(arg, "--gummel") && has_value) {
                config.gummel_iterations = std::stoi(argv[++i]);
            } else if (matches(arg, "--output") && has_value) {
                output_path = argv[++i];
            } else {
                std::cerr << "unrecognised argument: " << arg << "\n\n";
                print_usage(argv[0]);
                return EXIT_FAILURE;
            }
        }

        std::cout << "1D capacitively-coupled argon discharge\n"
                  << "  " << config.pressure_mtorr << " mTorr, "
                  << config.rf_voltage << " V at "
                  << config.rf_frequency / 1.0e6 << " MHz, "
                  << config.gap_length * 1000.0 << " mm gap\n\n";

        ccp::DischargeSolver solver(config);
        const ccp::CycleDiagnostics diagnostics = solver.run(verbose);

        std::cout << "\nfinal-cycle diagnostics\n"
                  << std::scientific << std::setprecision(4)
                  << "  peak electron density   " << diagnostics.peak_electron_density
                  << " m^-3\n"
                  << "  peak ion density        " << diagnostics.peak_ion_density
                  << " m^-3\n"
                  << "  mean electron density   " << diagnostics.mean_electron_density
                  << " m^-3\n"
                  << "  peak mean energy        " << std::fixed
                  << std::setprecision(2) << diagnostics.peak_mean_energy << " eV\n"
                  << std::scientific << std::setprecision(4)
                  << "  absorbed power density  " << diagnostics.absorbed_power_density
                  << " W/m^3\n"
                  << "  ionisation rate         " << diagnostics.ionisation_rate
                  << " m^-3 s^-1\n"
                  << "  drift over last cycle   " << diagnostics.cycle_drift
                  << "  (after " << diagnostics.cycles_completed << " cycles)\n";
        if (diagnostics.cycle_drift > 1.0e-3) {
            std::cout << "\n  note: the discharge has not reached periodic steady\n"
                      << "        state. Equilibration is set by ion transit, so\n"
                      << "        allow more cycles for a converged answer.\n";
        }

        solver.write_profiles(output_path);
        std::cout << "\nprofiles written to " << output_path << "\n";

    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

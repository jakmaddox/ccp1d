#include <cstdlib>
#include <exception>
#include <iomanip>
#include <iostream>
#include <string>

#include "discharge.hpp"

namespace {

void print_usage(const char* program) {
    // Defaults are read from DischargeConfig itself so this text cannot drift
    // out of step with the code again.
    const ccp::DischargeConfig d;
    std::cout
        << "usage: " << program << " [options]\n\n"
        << "  --pressure <mTorr>    neutral argon pressure       (default " << d.pressure_mtorr << ")\n"
        << "  --voltage <V>         RF amplitude                 (default " << d.rf_voltage << ")\n"
        << "  --frequency <Hz>      drive frequency              (default " << d.rf_frequency << ")\n"
        << "  --gap <m>             electrode separation         (default " << d.gap_length << ")\n"
        << "  --cells <n>           grid cells                   (default " << d.num_cells << ")\n"
        << "  --cycles <n>          RF cycles to simulate        (default " << d.num_cycles << ")\n"
        << "  --steps <n>           timesteps per RF cycle       (default " << d.steps_per_cycle << ")\n"
        << "  --gummel <n>          Gummel sweeps per step       (default " << d.gummel_iterations << ")\n"
        << "  --drift-cap <f>       electron drift cap, multiple of v_th (default "
        << d.drift_speed_limit_factor << "; 0 disables)\n"
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
            } else if (matches(arg, "--drift-cap") && has_value) {
                const double f = std::stod(argv[++i]);
                config.drift_speed_limit_factor = (f > 0.0) ? f : 1.0e12;
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
                  << "  mean electron density   " << diagnostics.mean_electron_density
                  << " m^-3   (cycle-averaged)\n"
                  << "  peak electron density   " << diagnostics.peak_electron_density
                  << " m^-3   (end of cycle)\n"
                  << "  peak ion density        " << diagnostics.peak_ion_density
                  << " m^-3   (end of cycle)\n"
                  << std::fixed << std::setprecision(2)
                  << "  peak mean energy        " << diagnostics.peak_mean_energy
                  << " eV     (end of cycle)\n"
                  << "  peak mean energy        " << diagnostics.peak_mean_energy_in_cycle
                  << " eV     (max over cycle; ceiling " << config.mean_energy_ceiling << ")\n"
                  << std::scientific << std::setprecision(4)
                  << "  electron Joule power    " << diagnostics.electron_power_density
                  << " W/m^3\n"
                  << "  ion Joule power         " << diagnostics.ion_power_density
                  << " W/m^3\n"
                  << "  total absorbed power    " << diagnostics.total_power_density
                  << " W/m^3\n"
                  << "  ionisation rate         " << diagnostics.ionisation_rate
                  << " m^-3 s^-1\n"
                  << std::fixed << std::setprecision(3)
                  << "  electron drift cap duty " << diagnostics.drift_cap_duty
                  << "          (max |mu_e E|/v_cap " << std::setprecision(2)
                  << diagnostics.drift_cap_max_ratio << ")\n"
                  << std::scientific << std::setprecision(4)
                  << "  drift over last cycle   " << diagnostics.cycle_drift
                  << "  (after " << diagnostics.cycles_completed << " cycles)\n";
        if (diagnostics.cycle_drift > 1.0e-4) {
            std::cout << "\n  note: the discharge has not reached periodic steady\n"
                      << "        state. Equilibration is set by ambipolar ion\n"
                      << "        transport, a few thousand RF periods at the\n"
                      << "        defaults, so allow more cycles for a converged answer.\n";
        }
        if (diagnostics.drift_cap_duty > 0.01) {
            std::cout << "\n  note: the electron drift cap was active on "
                      << std::fixed << std::setprecision(1)
                      << 100.0 * diagnostics.drift_cap_duty << " % of face\n"
                      << "        evaluations. Where it binds, electron transport is set\n"
                      << "        by the cap rather than the mobility model.\n";
        }

        solver.write_profiles(output_path);
        std::cout << "\nprofiles written to " << output_path << "\n";

    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

#!/usr/bin/env python3
"""Plot the cell-centred profiles written by ccp1d.

usage: python3 tools/plot_profiles.py profiles.csv [output.png]
"""

import csv
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_profiles(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle))
    columns = {key: [float(row[key]) for row in rows] for key in rows[0]}
    columns["x_mm"] = [x * 1000.0 for x in columns["x_m"]]
    return columns


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1

    source = sys.argv[1]
    destination = sys.argv[2] if len(sys.argv) > 2 else "profiles.png"
    data = read_profiles(source)
    x = data["x_mm"]

    figure, axes = plt.subplots(3, 1, figsize=(7.5, 9.0), sharex=True)

    axes[0].semilogy(x, data["n_e_m3"], label=r"$n_e$", linewidth=1.8)
    axes[0].semilogy(x, data["n_i_m3"], label=r"$n_i$", linewidth=1.8, linestyle="--")
    axes[0].set_ylabel(r"density  [m$^{-3}$]")
    axes[0].legend(frameon=False)
    axes[0].set_title("1D capacitively-coupled argon discharge")

    axes[1].plot(x, data["mean_energy_eV"], color="tab:red", linewidth=1.8)
    axes[1].set_ylabel("mean electron\nenergy  [eV]")

    axes[2].plot(x, data["potential_V"], color="tab:green", linewidth=1.8,
                 label="potential")
    axes[2].set_ylabel("potential  [V]")
    axes[2].set_xlabel("position  [mm]")

    twin = axes[2].twinx()
    twin.plot(x, data["charge_density_Cm3"], color="tab:purple", linewidth=1.2,
              alpha=0.7, label="net charge")
    twin.set_ylabel(r"net charge  [C m$^{-3}$]", color="tab:purple")
    twin.tick_params(axis="y", colors="tab:purple")

    for axis in axes:
        axis.grid(alpha=0.25)

    figure.tight_layout()
    figure.savefig(destination, dpi=150)
    print(f"wrote {destination}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

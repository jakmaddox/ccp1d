#!/usr/bin/env python3
"""Build a ccp1d transport table from electron cross sections, using BOLOS.

Solves the two-term Boltzmann equation over a range of E/N and writes mean
energy, mu*N, D*N and the elastic, excitation and ionisation rate coefficients
in the format that `ccp1d --transport-table` reads.

    pip install bolos numpy scipy
    python tools/make_transport_table.py argon_cross_sections.txt -o argon.table

The input is a BOLSIG+/LXCat format file. Cross sections are not distributed
with ccp1d: see "Transport tables" in the README for where to get them. Copy
the reference block from that file into anything you publish.

`--check` runs the solver on a built-in made-up gas instead. Its cross section
falls as 1/sqrt(energy), so the collision frequency is constant and the exact
mobility is e/(m k_m) at every E/N. No real data is involved.
"""

import argparse
import datetime
import io
import sys
import warnings

import numpy as np
import scipy.integrate

# BOLOS 0.2 calls this name, which SciPy has since removed.
if not hasattr(scipy.integrate, "simps"):
    scipy.integrate.simps = scipy.integrate.simpson

from bolos import grid, parser, solver  # noqa: E402

# The far tail of the distribution underflows inside BOLOS, which is harmless.
warnings.filterwarnings("ignore", category=RuntimeWarning, module="bolos")

E_OVER_M = 1.75882001076e11  # C/kg
EXCITATION_COST = 11.55  # eV, argon::kExcitationThreshold in include/argon.hpp


def synthetic_cross_sections(k_m):
    """A BOLSIG+ file for an elastic gas with constant momentum-transfer rate."""
    speed_factor = np.sqrt(2.0 * 1.602176634e-19 / 9.1093837015e-31)
    energies = np.logspace(-3, 3, 400)
    sigma = k_m / (speed_factor * np.sqrt(energies))
    rows = "\n".join(f" {e:.8e}\t{s:.8e}" for e, s in zip(energies, sigma))
    dashes = "-----------------------------"
    # BOLOS fails with no growth process, so add one too small to matter.
    negligible = " 300.0\t1.0e-40\n 1000.0\t1.0e-40"
    return (
        "ELASTIC\nAr\n 1.3734e-5\nCOMMENT: made-up gas with constant collision frequency\n"
        f"{dashes}\n{rows}\n{dashes}\n\n"
        f"IONIZATION\nAr -> Ar^+\n 300.0\n{dashes}\n{negligible}\n{dashes}\n"
    )


def build_solver(stream, gas, energy_max, cells):
    boltzmann = solver.BoltzmannSolver(grid.QuadraticGrid(0.0, energy_max, cells))
    boltzmann.load_collisions(parser.parse(stream))
    boltzmann.target[gas].density = 1.0
    boltzmann.kT = 300.0 * 8.617333262e-5
    return boltzmann


def rate_groups(boltzmann):
    """Split collisions into the three kinds of rate ccp1d carries."""
    elastic, excitation, ionisation = [], [], []
    for _, process in boltzmann.iter_all():
        kind = process.kind.upper()
        if kind in ("ELASTIC", "EFFECTIVE"):
            elastic.append(process)
        elif kind == "EXCITATION":
            excitation.append(process)
        elif kind == "IONIZATION":
            ionisation.append(process)
    return elastic, excitation, ionisation


def solve_point(boltzmann, en_townsend, seed):
    boltzmann.EN = en_townsend * solver.TOWNSEND
    boltzmann.init()
    return boltzmann.converge(seed, maxn=400, rtol=1e-5)


def sweep(boltzmann, fields):
    """Yield one table row per E/N.

    ccp1d charges every excitation the fixed cost EXCITATION_COST, so each
    level's rate is weighted by threshold/cost. That keeps the energy lost to
    excitation right when the data has many levels. k_excitation is therefore
    an energy-loss rate, not a count of excitation events.
    """
    elastic, excitation, ionisation = rate_groups(boltzmann)
    if not elastic:
        raise SystemExit("no ELASTIC or EFFECTIVE process in the cross-section file")

    seed = boltzmann.maxwell(1.0)
    for en in fields:
        f0 = solve_point(boltzmann, en, seed)
        seed = f0
        yield (
            float(boltzmann.mean_energy(f0)),
            float(boltzmann.mobility(f0)),
            float(boltzmann.diffusion(f0)),
            sum(float(boltzmann.rate(f0, p)) for p in elastic),
            sum(float(boltzmann.rate(f0, p)) * p.threshold / EXCITATION_COST
                for p in excitation),
            sum(float(boltzmann.rate(f0, p)) for p in ionisation),
        )


def monotonic(rows):
    """Keep rows whose mean energy rises, the axis ccp1d interpolates on.

    Below about 1 Td a gas barely heats, so neighbouring rows differ by solver
    noise and are dropped here.
    """
    kept = []
    for row in sorted(rows, key=lambda r: r[0]):
        if not kept or row[0] > kept[-1][0] * (1.0 + 1e-3):
            kept.append(row)
    return kept


def run_check():
    k_m = 1.0e-13
    exact = E_OVER_M / k_m
    stream = io.StringIO(synthetic_cross_sections(k_m))
    boltzmann = build_solver(stream, "Ar", 400.0, 1000)
    seed = boltzmann.maxwell(1.0)
    worst = 0.0
    for en in (0.3, 1.0, 3.0, 10.0):
        f0 = solve_point(boltzmann, en, seed)
        seed = f0
        error = float(boltzmann.mobility(f0)) / exact - 1.0
        worst = max(worst, abs(error))
        print(f"E/N {en:6.1f} Td  mean energy {float(boltzmann.mean_energy(f0)):10.4f} eV"
              f"  mu*N error {error:+.2e}")
    print(f"exact mu*N {exact:.5e}, worst error {worst:.2e}")
    return 0 if worst < 5.0e-3 else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("cross_sections", nargs="?", help="BOLSIG+/LXCat format file")
    ap.add_argument("-o", "--output", help="table to write")
    ap.add_argument("--gas", default="Ar", help="target name in the file (default Ar)")
    ap.add_argument("--en-min", type=float, default=0.0005, help="lowest E/N, Td")
    ap.add_argument("--en-max", type=float, default=3000.0, help="highest E/N, Td")
    ap.add_argument("--points", type=int, default=120, help="number of E/N values")
    ap.add_argument("--energy-max", type=float, default=400.0, help="grid top, eV")
    ap.add_argument("--cells", type=int, default=1000, help="energy grid cells")
    ap.add_argument("--check", action="store_true", help="run the exact test and exit")
    args = ap.parse_args()

    if args.check:
        return run_check()
    if not args.cross_sections or not args.output:
        ap.error("need a cross-section file and --output, or --check")

    with open(args.cross_sections, encoding="utf-8", errors="replace") as stream:
        boltzmann = build_solver(stream, args.gas, args.energy_max, args.cells)

    fields = np.logspace(np.log10(args.en_min), np.log10(args.en_max), args.points)
    rows = monotonic(list(sweep(boltzmann, fields)))
    with open(args.output, "w", encoding="utf-8") as out:
        out.write(f"# generated {datetime.date.today()} by tools/make_transport_table.py\n"
                  f"# from {args.cross_sections}, BOLOS two-term, "
                  f"E/N {args.en_min}-{args.en_max} Td\n"
                  "# mean_energy_eV mobility_N_per_m_V_s diffusivity_N_per_m_s "
                  "k_elastic_m3_s k_excitation_m3_s k_ionisation_m3_s\n")
        for row in rows:
            out.write(" ".join(f"{max(value, 0.0):.8e}" for value in row) + "\n")
    print(f"wrote {len(rows)} rows, mean energy {rows[0][0]:.3f} to {rows[-1][0]:.3f} eV")
    return 0


if __name__ == "__main__":
    sys.exit(main())

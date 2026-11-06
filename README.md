# ccp1d

A 1D drift-diffusion model of a capacitively-coupled RF discharge in argon,
written in C++17. Electron and ion continuity are coupled to Poisson's equation
and to an electron energy equation, closed with the local mean energy
approximation.

It reproduces the structure you expect from a CCP: a quasineutral bulk, two
space-charge sheaths, and an electron temperature that is low in the bulk and
rises sharply where the field does work. It also checks itself: the test suite
verifies global energy conservation, mirror symmetry, and the observed order of
accuracy, and the solver reports at runtime how much of the answer its
numerical guards are setting.

![Discharge profiles](docs/profiles.png)

*End-of-cycle snapshot after 3000 RF cycles at the default settings.*

## Results

1 Torr, 300 V, 13.56 MHz, 20 mm gap, 200 cells, 1000 steps per cycle, run for
3000 RF cycles to a cycle-to-cycle drift of 7.6e-6 (`./build/ccp1d --cycles 3000`):

| quantity | value |
| --- | --- |
| mean electron density (gap and cycle average) | 1.04e16 m^-3 |
| peak electron density | 2.80e16 m^-3 |
| bulk electron temperature | 1.57 V |
| peak mean electron energy (end of cycle / max over cycle) | 13.0 / 18.2 eV |
| electron Joule power | 8.19e3 W m^-3 |
| ion Joule power | 1.12e4 W m^-3 |
| total absorbed power | 1.94e4 W m^-3 |
| electron drift cap duty | 27.7% of face evaluations |

Two things about these numbers are easy to get wrong:

- **Convergence takes thousands of cycles, not hundreds.** The equilibration
  time is the ambipolar diffusion time, `(L/pi)^2 / (mu_i Te)`, about 1.9e-4 s
  or roughly 2600 RF periods at these settings. At the default of 400 cycles,
  which keeps a run near a minute, the mean density is still about 36% below
  its periodic steady state.
  The solver warns while the cycle drift is above 1e-4.
- **The ions absorb more than half the power** (58% here). The electron term
  `-Gamma_e . E` is what drives the energy equation, but ion Joule heating
  `e Gamma_i . E` in the sheaths is larger at this pressure, so both are
  reported. Earlier versions reported only the electron term as "absorbed
  power".

## Model

State variables, all at cell centres:

| symbol | quantity | units |
| --- | --- | --- |
| `n_e` | electron density | m^-3 |
| `n_i` | Ar+ density | m^-3 |
| `n_eps` | electron energy density | eV m^-3 |

The governing system is

```
  d(n_e)/dt   + div(Gamma_e)   =  k_iz N n_e
  d(n_i)/dt   + div(Gamma_i)   =  k_iz N n_e
  d(n_eps)/dt + div(Gamma_eps) = -Gamma_e . E  -  losses
  div(eps0 E) = e (n_i - n_e)
```

with drift-diffusion fluxes `Gamma = +-mu n E - D grad(n)`, energy transported at
5/3 of the particle mobility and diffusivity, and every rate coefficient and
diffusivity evaluated at the local electron temperature `Te = (2/3) n_eps / n_e`.

Reaction set: elastic scattering, excitation (11.55 eV), and ionisation
(15.76 eV), using the argon rate fits from Lieberman & Lichtenberg. Secondary
electron emission at the electrodes is included with a configurable yield.

## Numerics

**Scharfetter-Gummel flux discretisation.** Central differencing of a
drift-diffusion flux goes unstable once the cell Peclet number exceeds two;
naive upwinding is stable but smears the sheath edge across several cells.
Scharfetter-Gummel exponential fitting interpolates between the two, and only
the ions reach the upwind end. Their cell Peclet number peaks near 370 at the
defaults and exceeds two on 64% of face evaluations, so ion transport is first
order in the sheath. Electrons stay near the central end: the drift cap holds
their cell Peclet number below 0.47 at 200 cells, and below 0.92 at 100. The
first order in the convergence study below is consistent with the upwinded
ions, but no test isolates them as its cause.

The transport part of each density matrix is an M-matrix. The electron
equation also puts ionisation on the diagonal, `-k_iz N`, which erodes
diagonal dominance; at the defaults the largest `k_iz N dt` is about 0.06, so
the margin is real but not unlimited. The density floors after each solve are
a backstop and have not been seen to activate.

**Semi-implicit Poisson.** Coupling Poisson explicitly to the continuity
equations imposes the dielectric relaxation limit `dt < eps0 / (e mu_e n_e)`,
which is about 45 ps at 1e16 m^-3 and about 16 ps at the 2.8e16 m^-3 peak,
below the 74 ps default timestep. Advancing the charge density inside
Poisson's equation and writing the drift current against the unknown field moves
the conductivity onto the left-hand side as an effective permittivity,

```
  eps_eff = eps0 + dt e (mu_i n_i + mu_e n_e)
```

which preserves the tridiagonal structure and removes the restriction.

It is not free. Gauss's law is enforced against a predicted charge density
rather than the one the transport solves then produce, so the total current
`J_cond + eps0 dE/dt` is not exactly uniform in x. At the default timestep the
violation averages 0.40% of the peak current over a cycle after 100 cycles
(1.7% worst) and 0.46% at periodic steady state (1.7% worst). It is first order
in dt; `test_conservation` pins it to that scaling.

The predictor has to start from the charge density at the start of the step.
Each sweep then refines `J`: after 20 cycles at the defaults, three sweeps give
0.33% against 0.59% for one. An earlier version started from the current Gummel
iterate, which counted the step's transport twice from the second sweep on. Its
mean violation was 1.84% (3.95% worst), and sweeps made it worse. No balance or
order check saw this, so `test_conservation` now requires the sweeps to help.

The conductivity in `eps_eff` is built from the same limited drift the
transport solve applies. Where the electron drift cap binds, the flux responds
to the field with mobility `v_cap/|E|`, not `mu_e`; using the bare mobility
there overstated sigma and enforced Gauss's law against a current the
continuity equation never carried. Matching the two cut the continuity
violation by about 12x, from 4.7% to 0.40%.

`--explicit-poisson` restores the naive treatment, and the failure mode is worth
seeing. It tracks the semi-implicit solution to three significant figures while
the discharge is still filling (3.9268e15 against 3.9261e15 at cycle 120) and
then diverges at cycle 127, as the rising density pushes the dielectric
relaxation time below the timestep. The instability is invisible until the
plasma density crosses the threshold, which is exactly what makes it worth
designing out rather than tuning around. The run loop checks the state after
each cycle, so a diverged run stops with a non-zero exit and names the cycle.

**Gummel iteration.** Each timestep sweeps Poisson, electrons, ions, and energy
in sequence, repeated a few times so the strongly nonlinear coupling between
space charge and field is converged before the step is accepted. Three sweeps
is converged: going to 6 or 12 changes the absorbed power by less than 0.02%.
Rate coefficients are evaluated once per sweep, before any species moves, so
the electron source, the ion source and the energy cost of ionisation all see
the same `k_iz`. Re-evaluating them between the electron and ion updates, as an
earlier version did, created electrons and ions in slightly different numbers:
a spurious net charge source of 8.6e-5 of the ionisation rate, now 3.5e-6 and
vanishing with dt. Each sub-solve is a tridiagonal system handled by the
Thomas algorithm.

## Verification

`ctest` runs two suites and three exit-code checks. `test_numerics` checks the
building blocks in isolation; `test_conservation` checks that the assembled
solver obeys the laws those blocks are supposed to add up to; the exit-code
checks confirm that `ccp1d` and `ccp1d_benchmark` fail on a NaN run and name the
cycle. Every quantity in the second suite is
reconstructed from the solver's public state by `include/verification.hpp`,
independently of the solver's own arithmetic, so a sign error or a missing
factor shows up as a broken balance instead of being reproduced by the check.

| check | measured | tolerance |
| --- | --- | --- |
| Poisson operator vs analytic parabola, observed order | 2.000 | 1.9 to 2.1 |
| vacuum gap: linear potential, E = V/L on every face | roundoff | 1e-9 |
| Bernoulli function vs long-double reference | roundoff | 1e-14 |
| global electron energy balance over one cycle | 3.0e-6 of input | 1e-3 |
| global charge balance: net charge vs net wall current | 3.5e-6 of production | 2.5e-5 |
| undriven decay stays mirror-symmetric | 4.8e-14 | 1e-10 |
| driven: n(x,t) vs n(L-x,t+T/2), relative to cycle drift | 1.45 | < 3 |
| observed order in space / time | 0.84 / 1.05 | 0.6 to 1.4 |
| current-continuity error, order in dt | 0.97 | 0.6 to 1.4 |
| continuity error, three Gummel sweeps over one | 0.63 | < 0.8 |

The energy tolerance is set below the smallest physical term in the budget
(secondary-electron energy, about 1.3e-3 of input) so that an error even there
fails it. The two residuals are not zero because the solver freezes rate
coefficients for a Gummel sweep while the checks evaluate the end-of-step
state; both vanish faster than dt (energy: 1.0e-5, 3.0e-6, 4.7e-7 at 500,
1000, 2000 steps per cycle) and with more sweeps (energy: 9.0e-8 at 6 sweeps,
8.2e-11 at 12).

To confirm the suite can fail, bugs were planted one at a time:

| planted bug | caught by |
| --- | --- |
| wall energy coefficient 1/3 changed to 1/4 | energy balance |
| secondary-electron energy sign flipped at one wall | energy balance, undriven symmetry |
| Ohmic heating doubled | six checks |
| electron secondary-emission sign flipped at one wall | energy balance, both symmetry tests, convergence |
| elastic loss factor 3 changed to 2 | energy balance |
| electron and ion ionisation sources mismatched by 1% | charge balance (reads exactly 1.0e-2) |
| Poisson predictor started from the Gummel iterate, not the old charge | Gummel sweep check (ratio 2.03) |
| energy transport factor 5/3 changed to 1 | **not caught** |

The last row is a real gap. Interior fluxes telescope, so no global balance can
see a wrong transport coefficient, and the symmetry and order checks are
indifferent to it. Closing it needs a manufactured-solution test of the
transport operators.

`ccp1d_audit` measures the same things at any resolution, plus diagnostics
that are not pass/fail. With `--convergence` it adds grid and timestep studies
with Richardson error estimates. At the defaults after 100 cycles
(`./build/ccp1d_audit --cycles 100 --convergence`):

| refinement | observed order | Richardson error at the middle level |
| --- | --- | --- |
| cells 100 / 200 / 400 | 0.86 | 200 cells reads 2.8% low |
| cells 200 / 400 / 800 | 0.94 | 400 cells reads 1.4% low |
| steps 250 / 500 / 1000 | 0.91 | 500 steps reads 1.3% high |
| steps 500 / 1000 / 2000 | 0.93 | 1000 steps reads 0.67% high |

So the default resolution carries about 3% spatial and 0.7% temporal error in
the mean density, of opposite sign. The two partly cancel at the defaults; that
is coincidence, not design, and it will not hold at other settings.

## Validity

The drift-diffusion closure needs the electron mean free path to be short
compared with the gradient scale. In a 20 mm gap, taking the thermal speed at
the 1.57 V bulk temperature over the collision rate (0.59 mm at 1 Torr from the
built-in mobility, 0.54 mm from the elastic rate fit):

| pressure | mean free path | mfp / gap |
| --- | --- | --- |
| 100 mTorr | 6 mm | 0.3 |
| 300 mTorr | 2 mm | 0.1 |
| 1 Torr | 0.6 mm | 0.03 |

Three limits matter even at 1 Torr, and the code now measures each of them:

- **The electron drift cap is doing physics.** A constant mobility has no
  high-field saturation, so in the sheath `mu_e E` reaches about five times the
  thermal speed. The drift is capped at `v_th`, and at the defaults that cap is
  active on about 28% of interior-face evaluations, out to roughly 5 mm from
  each electrode. Removing it (`--drift-cap 0`) lowers the electron power by
  about 19% and the ion power by about 21%. Inside the capped region the model
  is free-streaming at `v_th`, not drift-diffusion. The duty is printed every
  run.
- **The rate fits are extrapolated.** The Lieberman fits are calibrated for
  about 1 to 7 V. At 1 Torr the hottest cells reach about 12 V, and 16% of all
  ionisation comes from above 7 V after 100 cycles, 22% at steady state.
- **The mean-energy ceiling is inert.** The 60 eV ceiling has not been seen to
  bind; the peak mean energy over a cycle is about 18 eV at 1 Torr and 21 eV at
  100 mTorr.

At 100 mTorr the model is outside its range and says so in its diagnostics:
the drift cap is active on 93% of face evaluations, so electron transport
across nearly the whole gap is set by the cap rather than the mobility, and
98.5% of ionisation comes from above the fits' calibrated range
(`./build/ccp1d_audit --pressure 100 --cells 100 --steps 500 --cycles 60`).
Below a few hundred mTorr, a kinetic treatment is required.

## Building

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Requires CMake 3.16+ and a C++17 compiler. No external dependencies. Every
target, tests and tools included, builds warning-clean under
`-Wall -Wextra -Wpedantic -Wshadow -Wconversion`. The test suite takes about
twenty seconds.

## Running

```
./build/ccp1d                                    # defaults: 1 Torr, 300 V, 13.56 MHz
./build/ccp1d --cycles 3000                      # to periodic steady state
./build/ccp1d --pressure 500 --voltage 200
./build/ccp1d --drift-cap 0                      # remove the electron drift cap
./build/ccp1d --transport-table argon.table      # Boltzmann transport and rates
./build/ccp1d --help

./build/ccp1d_audit                              # verification report
./build/ccp1d_audit --cycles 100 --convergence   # plus grid and dt studies
./build/ccp1d_benchmark --voltage 75 --yield 0.1 --emission-energy 0.75 --cycles 6000
```

Profiles are written as CSV and can be plotted with

```
python3 tools/plot_profiles.py profiles.csv docs/profiles.png
```

## Transport tables

By default the electron closure is a constant `mu N`, the Einstein relation and
the Lieberman rate fits. `--transport-table` replaces all of these with
coefficients tabulated against mean energy, usually from a Boltzmann solver. It
is optional, and ccp1d builds and runs without it.

Tables are made from electron cross sections, which ccp1d does not include.
LXCat asks third parties to send users to its site for the data, so:

1. Download a complete argon set from [LXCat](https://www.lxcat.net) in
   BOLSIG+ format. The site asks you to agree to its terms first, and the
   file's header gives the references to cite.
2. Install the generator's dependencies and build the table:

   ```
   pip install bolos numpy scipy
   python tools/make_transport_table.py argon_cross_sections.txt -o argon.table
   ```

3. Run with `./build/ccp1d --transport-table argon.table`.

LXCat's terms say its databases are not to be distributed by third parties or
used for commercial purposes, so neither the cross sections nor tables made from
them belong in this repository (`*.table` is ignored). Commercial use needs
permission from the data's contributor: see LXCat's redistribution policy.

`python tools/make_transport_table.py --check` tests the generator with no real
data, on a made-up gas whose exact mobility is known. It agrees to 0.2% from
0.3 to 10 Td.

The table is `mean_energy mu_N D_N k_elastic k_excitation k_ionisation` in SI
units and eV. Mean energy is the axis because the model is local mean energy,
and the drift cap and 5/3 energy transport factor still apply. The excitation
column is weighted by threshold so that a set with many levels loses the right
energy at the fixed 11.55 eV cost; it is not a count of excitation events.

The table path is checked against the built-in model: the built-in closure
written as a table reproduces a 20-cycle discharge to 8e-6 in `n_e`, and the
energy and charge balances close with a closure unrelated to the built-in one.
### Phelps argon

The default case at 3000 cycles, built-in closure against a table from the Phelps
argon set:

| quantity | built-in | Phelps table |
| --- | --- | --- |
| mean electron density | 1.04e16 m^-3 | 7.77e15 m^-3 |
| peak electron density | 2.80e16 m^-3 | 1.13e16 m^-3 |
| bulk mean energy | 2.4 eV (Te 1.57 V) | 5.65 eV at mid-gap |
| total absorbed power | 1.94e4 W m^-3 | 1.01e5 W m^-3 |
| electron drift cap duty | 27.7% | 2.9% |
| drift over last cycle | 7.6e-6 | 3.2e-12 |

Power rises about fivefold, mainly because the Boltzmann mobility near 5 eV is
about 1e24 1/(m V s), four times below the built-in 4.0e24. The drift cap almost
stops binding. Mesh convergence holds: 200 cells reads 3.2% low at 100 cycles,
with observed order 0.95.

These numbers are not validated, as the benchmark below shows. Do not read the
table column as more correct than the built-in one.

The energy budget closes less well with a table: 2.1% of input against under
0.1%. Measured directly, nearly all of that is negative Ohmic heating in a few
cold cells at the sheath edge. There the table's `D/mu` exceeds what the local
mean energy can supply, so the Patankar limiter drops the part of the sink that
would drive the energy negative. Keep the table's lower end near thermal energy:
stopping it at 0.8 eV gave a 10% residual, which is why the generator starts at
5e-4 Td.

Data: Phelps database, www.lxcat.net, retrieved on November 6, 2026, based on
Yamabe, Buckman and Phelps, *Phys. Rev. A* **27** (1983) 1345, revised 1997.

### Benchmark against a measured discharge

Tsagkaridis et al. (2025) report a symmetric 13.56 MHz argon discharge at 1 Torr
with a 2 cm gap, 150 V peak-to-peak, one electrode grounded. ccp1d ran this at
75 V amplitude, `gamma` = 0.1 and 0.75 eV emitted electrons, for 6000 cycles.
Values are cycle averages at mid-gap, from `ccp1d_benchmark`:

| quantity | measured | their fluid model | ccp1d built-in | ccp1d Phelps |
| --- | --- | --- | --- | --- |
| plasma potential | 46.3 V | 47.6 V | 40.5 V | 59.3 V |
| electron density | ~2e15 m^-3 | ~3.5e16 m^-3 | 4.4e15 m^-3 | 1.5e15 m^-3 |
| electron temperature | not reported | 2.92 V | 1.53 V | 3.88 V |

Neither closure reproduces the discharge. The plasma potential, the firmer
measurement, is 13% low with the built-in closure and 28% high with Phelps. The
densities are read from a figure, and the authors rate the probe density as an
order-of-magnitude estimate, so being within a factor of two of it shows little.
Their own model, with far more physics, overpredicts it by a factor of 16.

The pattern is informative, though. With a Maxwellian distribution their model
gives 1.24 V and 38.6 V, close to the built-in closure, which also assumes one.
With Boltzmann distributions it gives 2.92 V and 47.6 V, and the Phelps table
lands above both. A likely reason, not tested here, is that ccp1d has no
metastables: without two-step ionisation the discharge must run hotter to
ionise from the ground state alone. Lymberopoulos and Economou (1993) found
two-step ionisation to be the main electron source at 1 Torr.

## Next steps

In rough order of value:

1. Metastable argon as a fourth species, to test whether two-step ionisation
   closes the benchmark gap above. Boltzmann tables should not become the
   default before that.
2. A manufactured-solution test for the transport operators, the one class of
   error the current suite cannot see.
3. Kokkos or CUDA offload of the flux assembly and tridiagonal solves. The
   assembly is embarrassingly parallel across faces; the tridiagonal solves need
   cyclic reduction rather than Thomas.
4. MPI domain decomposition, for extension to 2D where the tridiagonal structure
   is lost and a Krylov solve with a multigrid preconditioner is needed instead.

## References

- Scharfetter & Gummel, *IEEE Trans. Electron Devices* **16** (1969) 64.
- Ventzek, Hoekstra & Kushner, *J. Vac. Sci. Technol. B* **12** (1994) 461.
- Lieberman & Lichtenberg, *Principles of Plasma Discharges and Materials
  Processing*, 2nd ed., Wiley (2005).
- Hagelaar & Pitchford, *Plasma Sources Sci. Technol.* **14** (2005) 722.
- BOLOS, a two-term Boltzmann solver: https://github.com/aluque/bolos
- Lymberopoulos & Economou, *J. Appl. Phys.* **73** (1993) 3668.
- Tsagkaridis et al., *Plasma Sources Sci. Technol.* **34** (2025) 105013.

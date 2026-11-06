# ccp1d

A 1D drift-diffusion model of a capacitively-coupled RF discharge in argon,
written in C++17. Electron and ion continuity are coupled to Poisson's equation
and to an electron energy equation, closed with the local mean energy
approximation.

It reproduces the structure you expect from a low-pressure CCP: a quasineutral
bulk, two space-charge sheaths, and an electron temperature that is low in the
bulk and rises sharply where the field does work.

![Discharge profiles](docs/profiles.png)

At 1 Torr, 300 V, 13.56 MHz across a 20 mm gap, after 200 RF cycles:

| quantity | value |
| --- | --- |
| peak electron density | 1.06e16 m^-3 |
| peak ion density | 1.06e16 m^-3 |
| bulk electron temperature | ~1.7 V |
| peak mean electron energy | 12.9 eV (in the sheath) |
| absorbed power density | 9.3e3 W m^-3 |
| sheath thickness | ~3 mm |

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

with drift-diffusion fluxes `Gamma = ±mu n E - D grad(n)`, energy transported at
5/3 of the particle mobility and diffusivity, and every rate coefficient and
diffusivity evaluated at the local electron temperature `Te = (2/3) n_eps / n_e`.

Reaction set: elastic scattering, excitation (11.55 eV), and ionisation
(15.76 eV), using the argon rate fits from Lieberman & Lichtenberg. Secondary
electron emission at the electrodes is included with a configurable yield.

## Numerics

Three choices do most of the work.

**Scharfetter-Gummel flux discretisation.** Central differencing of a
drift-diffusion flux goes unstable once the cell Peclet number exceeds two,
which in a sheath happens immediately; naive upwinding is stable but smears the
sheath edge across several cells. Scharfetter-Gummel exponential fitting
interpolates between the two automatically, and the resulting matrix is an
M-matrix — positive diagonal, non-positive off-diagonals, diagonally dominant —
so densities are guaranteed non-negative without clipping. `tests/` checks that
the scheme collapses to Fickian diffusion at zero field and to pure upwinding at
large drift.

**Semi-implicit Poisson.** Coupling Poisson explicitly to the continuity
equations imposes the dielectric relaxation limit `dt < eps0 / (e mu_e n_e)`,
which at 1e16 m^-3 is about 4 ps — roughly twenty times smaller than the
timestep needed to resolve 13.56 MHz. Advancing the charge density inside
Poisson's equation and writing the drift current against the unknown field moves
the conductivity onto the left-hand side as an effective permittivity,

```
  eps_eff = eps0 + dt e (mu_i n_i + mu_e n_e)
```

which preserves the tridiagonal structure and removes the restriction entirely.

`--explicit-poisson` restores the naive treatment, and the failure mode is worth
seeing. It tracks the semi-implicit solution to four significant figures while
the discharge is still filling — 3.9263e15 against 3.9287e15 at cycle 120 — and
then goes to NaN between cycles 120 and 140, as the rising density pushes the
dielectric relaxation time below the timestep. The instability is invisible
until the plasma density crosses the threshold, which is exactly what makes it
worth designing out rather than tuning around.

**Gummel iteration.** Each timestep sweeps Poisson, electrons, ions, and energy
in sequence, repeated a few times so the strongly nonlinear coupling between
space charge and field is converged before the step is accepted. Each sub-solve
is a tridiagonal system handled by the Thomas algorithm.

## Validity

The drift-diffusion closure needs the electron mean free path to be short
compared with the gradient scale. In a 20 mm gap:

| pressure | mean free path | mfp / gap |
| --- | --- | --- |
| 100 mTorr | 31 mm | 1.6 |
| 300 mTorr | 10 mm | 0.5 |
| 1 Torr | 3.1 mm | 0.16 |

At 100 mTorr the model is being asked to describe a nearly collisionless plasma
and it fails visibly — electrons free-stream to the walls, quasineutrality never
establishes, and the ion density exceeds the electron density across the entire
gap. That is the model reporting its own breakdown rather than a bug, and it is
why the default pressure is 1 Torr. Below a few hundred mTorr, a kinetic
treatment is required.

Two other limits are worth stating plainly:

- **Constant mobility.** Real codes take `mu(E/N)` from a Boltzmann solver such
  as BOLSIG+. A constant mobility has no high-field saturation, so a 300 V
  sheath would otherwise drive electrons to an appreciable fraction of *c*. The
  drift speed is capped at the local thermal speed as a numerical guard, which
  is roughly where the fluid closure stops being meaningful anyway.
- **Convergence is slow.** Equilibration is set by ion transit across the gap,
  which is hundreds of RF periods. The solver reports the cycle-to-cycle drift
  in gap-averaged density and warns when it has not settled.

## Building

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build
```

Requires CMake 3.16+ and a C++17 compiler. No external dependencies. The build
is warning-clean under `-Wall -Wextra -Wpedantic -Wshadow -Wconversion`.

## Running

```
./build/ccp1d                                    # defaults: 1 Torr, 300 V, 13.56 MHz
./build/ccp1d --pressure 500 --voltage 200
./build/ccp1d --cells 400 --steps 2000 --cycles 800
./build/ccp1d --help
```

Profiles are written as CSV and can be plotted with

```
python3 tools/plot_profiles.py profiles.csv docs/profiles.png
```

## Next steps

In rough order of value:

1. Tabulated `mu(E/N)`, `D(E/N)` and rate coefficients from a Boltzmann solver,
   replacing the constant-mobility closure and its speed cap.
2. Kokkos or CUDA offload of the flux assembly and tridiagonal solves. The
   assembly is embarrassingly parallel across faces; the tridiagonal solves need
   cyclic reduction rather than Thomas.
3. MPI domain decomposition, for extension to 2D where the tridiagonal structure
   is lost and a Krylov solve with a multigrid preconditioner is needed instead.
4. Surface chemistry beyond secondary emission — sticking, recombination, and a
   deposition or etch model at the electrodes.
5. Metastable argon as a fourth species, which matters for the ionisation
   balance well before anything else on this list does.

## References

- Scharfetter & Gummel, *IEEE Trans. Electron Devices* **16** (1969) 64.
- Ventzek, Hoekstra & Kushner, *J. Vac. Sci. Technol. B* **12** (1994) 461.
- Lieberman & Lichtenberg, *Principles of Plasma Discharges and Materials
  Processing*, 2nd ed., Wiley (2005).
- Hagelaar & Pitchford, *Plasma Sources Sci. Technol.* **14** (2005) 722.

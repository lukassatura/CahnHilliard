# Cahn-Hilliard
Dynamic simulation of density profile evolution in multi-phase liquid-vapor and liquid-liquid systems.

## Author
**Lukáš Šatura**

## Date
06/04/2021

This project simulates the evolution of density profiles in one dimension using the Peng-Robinson equation of state. It is a computational (*in-silico*) base for the research work done in [my thesis](https://repozitar.vscht.cz/theses/31278) (it explains the fundamental thermodynamic concepts quite exhaustively).

For more thorough studies, see my papers:
- [A Robust Physics-Based Calculation of Evolving Gas–Liquid Interfaces](https://doi.org/10.1515/jnet-2021-0080)
- [*Ab Initio* Prediction of Surface Tension from Fundamental Equations of State using Density Gradient Theory (DGT)](https://doi.org/10.1016/B978-0-443-28824-1.50090-9)

## Single-Component Simulation of Interfacial Density Profile Evolution in 1D with 'Periodic' Boundary Conditions (hexane)
![Single-Component Simulation of Interfacial Density Profile Evolution in 1D with 'Periodic' Boundary Conditions (hexane)](/n_alkanes_1D/graphs/hexane_rho_periodic_BC.png)

## Binary Simulation of Interfacial Density Profile Evolution in 1D with Dirichlet Boundary Conditions
![Binary Simulation of Interfacial Density Profile Evolution in 1D with Dirichlet Boundary Conditions](/n_hexane_N2_1D/hexane_N2_Dirichlet_1D.png)

## Files
The project includes the following files (yes, good old Fortran is still used!):

- `io.f90`
- `main.f90`
- `potentials.f90`
- `rho_profiles.plt`
- `rho_search.m`
- `potentials/mu_profiles.plt`

## Prerequisites
This project was developed in Fortran 90 (for historic reasons of the research group 🤓 and for computational efficiency). I used Intel® Fortran Compiler (`ifort`) on my Mac's Intel chip as it supported compiled code optimisation for this processor and support for OpenMP (parallelisation, in fact, turned out to be non-realisable for my numerical schema though). If you are an Apple Silicon user, consider using the GNU Fortran (GFortran) compiler.

## Compilation
To compile the `.f90` files, use the following commands in the terminal (again, `ifort` compiler for Intel processor is used here with the `-xHost` flag for automatic code optimisation):

```bash
ifort -xHost io.f90 main.f90
```
or
```bash
ifort -xHost io.f90 potentials.f90
```

## File Descriptions

### `io.f90`
- Contains variables for simulation parameter modification, such as:
  - Number of iterations (`M`)
  - Grid density (`N`)
  - Time step (`k`)
  - Periodic boundary conditions (`per_BC`)
- Critical parameters for each substance are sourced from `critical_parameters.xlsx` (source to be provided).
- Density values (`rho_L` and `rho_V`) were calculated in MATLAB using `rho_search.m` for a specific temperature (`T = 298.15 K`). These values significantly affect simulation stability.

### `main.f90`
- Implements the explicit Euler method with fixed or periodic boundary conditions.
- Recommended parameters for *n*-hexane:
  - `N = 100`
  - `M = 2000000`
  - `k = 1.0d-7`
  - `L = 1.0d-8`
  (satisfactory results even with periodic boundary conditions).
- Outputs surface tension values to the terminal (e.g., `16.58 mN/m` for hexane). These can be compared against reference values from [Surface Tension Reference](http://www.surface-tension.de/).
- Notes on stability:
  - Fourth-order accuracy formulas for boundary points require smaller time steps and denser grids. Neumann boundary conditions (`dy_dx = 0`) improve stability.

### `potentials.f90`
- Computes generalized chemical potential profiles from calculated density profiles.
- Fluxes `J` can be calculated by uncommenting the relevant section.
- Verifies steady-state simulation by checking if the resulting profile is constant across the domain. If not, increasing `M` may help.

## Plotting Profiles
To plot the profiles, use **gnuplot**:

```bash
gnuplot rho_profiles.plt
gnuplot mu_profiles.plt
```

## Final Remarks
Do you have any comments or questions to my Cahn-Hilliard model? Are you researching on a similar topic? Do let me know!

I put significant amount of research work into this project. I'll be grateful if you consider citing it properly according to scientific standards should you find the research findings useful.

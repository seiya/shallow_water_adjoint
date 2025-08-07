# shallow_water_adjoint
Shallow water equation model with adjoint model.

The testcases included here follow the benchmark proposed by Williamson
et al. (1992). Instead of the original spherical geometry, they solve a
channel flow in Cartesian coordinates with free-slip boundaries in the
y direction and periodic boundaries in the x direction.

## Test cases and executables

The project provides three test cases (1, 2 and 5). For each test case
three executables are built:

- `shallow_water_testX.out` – original program
- `shallow_water_testX_forward.out` – forward mode example
- `shallow_water_testX_reverse.out` – reverse mode example

where `X` is `1`, `2` or `5`.

## Command-line arguments

All executables accept an optional first command-line argument specifying
the snapshot output interval. The default is `48`; specifying `0` outputs
only the final-time snapshot, and a negative value disables snapshot
output entirely. Many programs also accept additional arguments to read
initial conditions or perturbations from binary files; see the source
code for details.

## Building

Build executables (required fautodiff modules are fetched automatically):

```bash
cd build
make
```

The build step will clone and copy `fautodiff` modules on demand.

# shallow_water_adjoint
Shallow water equation model with adjoint model

The test program `shallow_water_test1` accepts two optional command-line
arguments:

1. Rotation angle `alpha` in degrees.
2. Snapshot output interval. The default is 48; specifying `0` outputs only
   the final-time snapshot, and `-1` disables snapshot output entirely.

## Building and running

Build executables (required fautodiff modules are fetched automatically):

```bash
cd build
make
```

This produces three binaries under `build`:

- `shallow_water_test1` – original program
- `shallow_water_forward` – forward mode example
- `shallow_water_reverse` – reverse mode example

The build step will clone and copy `fautodiff` modules on demand.

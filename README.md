# shallow_water_adjoint
Shallow water equation model with adjoint model

The test program `shallow_water_test1` accepts two optional command-line
arguments:

1. Rotation angle `alpha` in degrees.
2. A flag to control snapshot output. Set this to `0` to disable writing
   `snapshot_*.bin` files; any non-zero value enables snapshots (default).

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

#!/usr/bin/env bash
set -euo pipefail

# Copy required Fortran modules from the fautodiff source tree
# into the current working directory (build directory).

FAUTODIFF_SRC="../scripts/fautodiff/fortran_modules"
if [ -d "$FAUTODIFF_SRC" ]; then
  cp -u "$FAUTODIFF_SRC"/*.f90 "$FAUTODIFF_SRC"/*.fadmod .
else
  tmp_dir=$(mktemp -d)
  git clone --depth 1 https://github.com/seiya/fautodiff "$tmp_dir" >/tmp/fautodiff_copy.log
  tail -n 20 /tmp/fautodiff_copy.log
  cp -u "$tmp_dir/fortran_modules"/*.f90 "$tmp_dir/fortran_modules"/*.fadmod .
  rm -rf "$tmp_dir"
fi

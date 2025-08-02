#!/usr/bin/env bash
set -euo pipefail

# Copy required Fortran modules from the installed fautodiff package
# into the current working directory (build directory).

src_dir=$(python3 - <<'PY'
from pathlib import Path
try:
    import fautodiff
    print(Path(fautodiff.__file__).resolve().parent.parent / "fortran_modules")
except ModuleNotFoundError:
    pass
PY
)

# If the installed package is unavailable or does not contain the expected
# Fortran sources, fall back to cloning the repository.
if [ -n "$src_dir" ] && [ -d "$src_dir" ] && ls "$src_dir"/*.f90 >/dev/null 2>&1; then
  cp "$src_dir"/*.f90 "$src_dir"/*.fadmod .
else
  tmp_dir=$(mktemp -d)
  git clone --depth 1 https://github.com/seiya/fautodiff "$tmp_dir" >/tmp/fautodiff_copy.log
  tail -n 20 /tmp/fautodiff_copy.log
  src_dir="$tmp_dir/fortran_modules"
  cp "$src_dir"/*.f90 "$src_dir"/*.fadmod .
  rm -rf "$tmp_dir"
fi

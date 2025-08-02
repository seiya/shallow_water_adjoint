#!/usr/bin/env bash
set -e
if ! command -v fautodiff >/dev/null 2>&1; then
  echo "Installing fautodiff..."
  # Install fautodiff into the global Python environment so that the
  # module and corresponding ``fautodiff`` executable are importable from
  # subsequent scripts.  ``--break-system-packages`` is required on some
  # systems (e.g. Debian/Ubuntu with PEP 668) to allow installation into
  # the system site-packages when not using a virtual environment.
  python3 -m pip install --break-system-packages \
    git+https://github.com/seiya/fautodiff >/tmp/fautodiff_install.log
  tail -n 20 /tmp/fautodiff_install.log
  # When using pyenv, newly installed executables are only available via
  # shims after ``pyenv rehash``.  Do this silently if pyenv exists.
  if command -v pyenv >/dev/null 2>&1; then
    pyenv rehash >/dev/null 2>&1 || true
  fi
else
  echo "fautodiff already installed."
fi

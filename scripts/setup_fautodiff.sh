#!/usr/bin/env bash
set -e
if ! command -v fautodiff >/dev/null 2>&1; then
  echo "Installing fautodiff..."
  pip install git+https://github.com/seiya/fautodiff >/tmp/fautodiff_install.log
  tail -n 20 /tmp/fautodiff_install.log
else
  echo "fautodiff already installed."
fi

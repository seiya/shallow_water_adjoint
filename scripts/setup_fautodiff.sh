#!/usr/bin/env bash
set -e
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FAUTODIFF_DIR="$SCRIPT_DIR/fautodiff"
if [ ! -d "$FAUTODIFF_DIR" ]; then
  echo "Cloning fautodiff..."
  git clone --depth 1 https://github.com/seiya/fautodiff "$FAUTODIFF_DIR" >/tmp/fautodiff_clone.log
  tail -n 20 /tmp/fautodiff_clone.log
else
  echo "fautodiff already cloned."
fi

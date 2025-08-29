#!/usr/bin/env bash
set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FAUTODIFF_DIR="$SCRIPT_DIR/fautodiff"

if [ ! -d "$FAUTODIFF_DIR/.git" ]; then
  echo "Cloning fautodiff..."
  # Try to clone, but do not fail build if network is unavailable.
  if git clone --depth 1 https://github.com/seiya/fautodiff "$FAUTODIFF_DIR" >/tmp/fautodiff_clone.log 2>&1; then
    tail -n 20 /tmp/fautodiff_clone.log || true
  else
    echo "Warning: network unavailable; skipping fautodiff clone (using existing files if present)."
  fi
else
  echo "Updating fautodiff..."
  # Try to update, but do not fail build if network is unavailable.
  if git -C "$FAUTODIFF_DIR" fetch --depth 1 origin main >/tmp/fautodiff_update.log 2>&1; then
    git -C "$FAUTODIFF_DIR" reset --hard origin/main >>/tmp/fautodiff_update.log 2>&1 || true
    tail -n 20 /tmp/fautodiff_update.log || true
  else
    echo "Warning: network unavailable; skipping fautodiff update."
  fi
fi

#!/usr/bin/env sh
set -eu

# Allow overriding the gcov binary via GCOV_BIN, otherwise use gcov from PATH.
GCOV_BIN="${GCOV_BIN:-gcov}"

if ! command -v "$GCOV_BIN" >/dev/null 2>&1; then
  echo "gcov wrapper: unable to find '$GCOV_BIN'" >&2
  exit 127
fi

# Resolve repository root from the location of this script.
REPO_DIR=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)

cd "$REPO_DIR/src"
exec "$GCOV_BIN" "$@"

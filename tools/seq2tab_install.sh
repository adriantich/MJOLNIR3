#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/bin/seq2tab"

mkdir -p "$ROOT/bin"

g++ -std=c++17 -O2 \
  -I"$ROOT/src" \
  "$ROOT/src/seq2tab_core.cpp" \
  "$ROOT/tools/seq2tab_cli.cpp" \
  -o "$OUT"

echo "Built CLI: $OUT"
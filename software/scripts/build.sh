#!/usr/bin/env bash
# Сборка всего суперпроекта software/ в каталог build/.
# Использование: scripts/build.sh [Release|Debug] [доп. аргументы cmake]
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
type="${1:-Release}"; shift || true
cmake -S "$here" -B "$here/build" -G Ninja -DCMAKE_BUILD_TYPE="$type" "$@" 2>/dev/null \
  || cmake -S "$here" -B "$here/build" -DCMAKE_BUILD_TYPE="$type" "$@"
cmake --build "$here/build" -j"$(nproc)"

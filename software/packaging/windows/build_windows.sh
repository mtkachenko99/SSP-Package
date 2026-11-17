#!/usr/bin/env bash
# Сборка Windows-дистрибутива в Docker (Fedora + MinGW-w64). Результат: software/packaging/out/ssp-windows-x64(.zip)
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src="$(cd "$here/../.." && pwd)"
out="$src/packaging/out"
mkdir -p "$out"
docker run --rm -e HOST_UID="$(id -u)" -e HOST_GID="$(id -g)" -v "$src":/src:ro -v "$out":/out -v ssp-mingw-cache:/opt fedora:41 bash /src/packaging/windows/build_in_container.sh

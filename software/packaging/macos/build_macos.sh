#!/usr/bin/env bash
# Сборка SSP Nano на macOS (Apple Silicon или Intel). Требуется Homebrew (https://brew.sh).
# Запуск из каталога software/:   bash packaging/macos/build_macos.sh
# Результат: packaging/out/ssp-macos-<arch>.zip — ssp_gui.app (с Qt внутри), ssp_aggregate, ssp_benchmark, examples/.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
src="$(cd "$here/../.." && pwd)"
arch="$(uname -m)"
brew install cmake ninja eigen fftw qt libomp pkg-config >/dev/null
build="$src/build-macos"
cmake -S "$src" -B "$build" -G Ninja -DCMAKE_BUILD_TYPE=Release -DSSP_BUILD_TESTS=OFF \
  -DCMAKE_PREFIX_PATH="$(brew --prefix qt);$(brew --prefix eigen);$(brew --prefix fftw)" \
  -DOpenMP_ROOT="$(brew --prefix libomp)"
cmake --build "$build"
out="$src/packaging/out/ssp-macos-$arch"
rm -rf "$out" && mkdir -p "$out"
cp -R "$build/gui/ssp_gui.app" "$out/"
cp "$build/excitons/ssp_aggregate" "$build/validation/ssp_benchmark" "$out/"
# Qt внутрь .app; консольные программы используют те же библиотеки из .app
"$(brew --prefix qt)/bin/macdeployqt" "$out/ssp_gui.app" -executable="$out/ssp_gui.app/Contents/MacOS/ssp_gui"
for b in ssp_aggregate ssp_benchmark; do
  install_name_tool -add_rpath "@executable_path/ssp_gui.app/Contents/Frameworks" "$out/$b" 2>/dev/null || true
done
cp "$(brew --prefix fftw)/lib/libfftw3.3.dylib" "$(brew --prefix libomp)/lib/libomp.dylib" "$out/ssp_gui.app/Contents/Frameworks/" 2>/dev/null || true
mkdir -p "$out/ssp_gui.app/Contents/MacOS/examples"
rsync -a --exclude results "$src/excitons/examples/" "$out/ssp_gui.app/Contents/MacOS/examples/"
rsync -a --exclude results "$src/excitons/examples/" "$out/examples/"
cp "$src/packaging/README_dist.txt" "$out/README.txt"
# подпись ad-hoc (без сертификата разработчика): иначе Gatekeeper на Apple Silicon не запустит бинарник
codesign --force --deep -s - "$out/ssp_gui.app" "$out/ssp_aggregate" "$out/ssp_benchmark"
(cd "$src/packaging/out" && rm -f "ssp-macos-$arch.zip" && ditto -c -k --keepParent "ssp-macos-$arch" "ssp-macos-$arch.zip")
echo "Done: $src/packaging/out/ssp-macos-$arch.zip"
echo "Colleagues: unzip, then right-click ssp_gui.app -> Open (first launch), or: xattr -dr com.apple.quarantine ssp-macos-$arch"

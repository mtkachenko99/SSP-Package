#!/usr/bin/env bash
# Кросс-сборка ssp_aggregate, ssp_gui и ssp_benchmark под Windows x86_64 (MinGW-w64, Qt6) внутри контейнера Fedora.
# Запуск (из каталога software/):  packaging/windows/build_windows.sh
# Внутри контейнера: /src — исходники (только чтение), /out — результат.
set -euo pipefail
dnf -y -q install mingw64-gcc-c++ mingw64-qt6-qtbase mingw64-qt6-qtsvg mingw64-winpthreads-static \
  cmake ninja-build make wget eigen3-devel pkgconf-pkg-config file zip >/dev/null
dnf -y -q install mingw64-libgomp >/dev/null 2>&1 || echo "note: mingw64-libgomp not available — serial build"
# FFTW3 для MinGW — из исходников (в Fedora нет пакета)
if [ ! -f /opt/fftw-mingw/lib/libfftw3.a ]; then
  cd /tmp && wget -q https://www.fftw.org/fftw-3.3.10.tar.gz && tar xzf fftw-3.3.10.tar.gz && cd fftw-3.3.10
  mingw64-configure --prefix=/opt/fftw-mingw --disable-shared --enable-static --disable-fortran \
    --with-our-malloc --enable-sse2 --enable-avx >/dev/null
  make -j"$(nproc)" >/dev/null && make install >/dev/null
fi
export PKG_CONFIG_PATH=/opt/fftw-mingw/lib/pkgconfig
export PKG_CONFIG_LIBDIR=/opt/fftw-mingw/lib/pkgconfig:/usr/x86_64-w64-mingw32/sys-root/mingw/lib/pkgconfig
mkdir -p /build && cd /build
mingw64-cmake -G Ninja /src -DCMAKE_BUILD_TYPE=Release -DSSP_BUILD_TESTS=OFF -DSSP_USE_OPENMP=ON -DBUILD_SHARED_LIBS=OFF \
  -DEigen3_DIR="$(dirname "$(find /usr/share -name Eigen3Config.cmake | head -1)")" -DCMAKE_PREFIX_PATH=/opt/fftw-mingw \
  -DCMAKE_FIND_ROOT_PATH_MODE_PACKAGE=BOTH
ninja -k 0 > /out/ninja_windows.log 2>&1 || { grep -E -B2 -A6 "error:|FAILED" /out/ninja_windows.log | head -80; exit 1; }
# Дистрибутив: программы + все нужные DLL + плагины Qt + примеры
D=/out/ssp-windows-x64
rm -rf "$D" && mkdir -p "$D/platforms" "$D/styles" "$D/imageformats"
for e in ssp_aggregate ssp_gui ssp_benchmark; do cp "$(find /build -name "$e.exe" | head -1)" "$D/"; done
SR=/usr/x86_64-w64-mingw32/sys-root/mingw
cp "$SR/lib/qt6/plugins/platforms/qwindows.dll" "$D/platforms/"
cp "$SR"/lib/qt6/plugins/styles/*.dll "$D/styles/" 2>/dev/null || true
cp "$SR"/lib/qt6/plugins/imageformats/qsvg.dll "$D/imageformats/" 2>/dev/null || true
# зависимости DLL — рекурсивно по таблице импорта
deps() { x86_64-w64-mingw32-objdump -p "$1" | awk '/DLL Name:/{print $3}'; }
queue=("$D"/*.exe "$D"/platforms/*.dll "$D"/styles/*.dll "$D"/imageformats/*.dll)
declare -A seen
while [ ${#queue[@]} -gt 0 ]; do
  f="${queue[0]}"; queue=("${queue[@]:1}")
  for dll in $(deps "$f"); do
    [ -n "${seen[$dll]:-}" ] && continue; seen[$dll]=1
    src=$(find "$SR/bin" -maxdepth 1 -iname "$dll" | head -1)
    if [ -n "$src" ]; then cp -n "$src" "$D/"; queue+=("$D/$(basename "$src")"); fi
  done
done
missing=""
for f in "$D"/*.exe "$D"/*.dll "$D"/platforms/*.dll; do
  for dll in $(deps "$f"); do
    [ -f "$D/$dll" ] && continue
    case "${dll,,}" in kernel32.dll|user32.dll|gdi32.dll|advapi32.dll|shell32.dll|ole32.dll|oleaut32.dll|msvcrt.dll|ws2_32.dll|\
      comdlg32.dll|imm32.dll|winmm.dll|version.dll|uxtheme.dll|dwmapi.dll|shcore.dll|shlwapi.dll|userenv.dll|setupapi.dll|\
      netapi32.dll|mpr.dll|authz.dll|wtsapi32.dll|d3d9.dll|d3d11.dll|d3d12.dll|dxgi.dll|dwrite.dll|api-ms-win-*|bcrypt.dll|crypt32.dll|secur32.dll|dnsapi.dll|iphlpapi.dll|ncrypt.dll|rpcrt4.dll|winspool.drv|opengl32.dll|comctl32.dll) ;;
      *) missing="$missing $dll($(basename "$f"))" ;;
    esac
  done
done
[ -n "$missing" ] && { echo "MISSING DLLs:$missing"; exit 1; }
echo "DLL check: all non-system imports present"
echo "[Paths]" > "$D/qt.conf"; echo "Plugins = ." >> "$D/qt.conf"
# примеры (без результатов) и описание
mkdir -p "$D/examples" && (cd /src/excitons/examples && find . -path ./results -prune -o -type f -print | while read -r f; do mkdir -p "$D/examples/$(dirname "$f")"; cp "$f" "$D/examples/$f"; done)
cp /src/packaging/README_dist.txt "$D/README.txt" 2>/dev/null || true
cd /out && rm -f ssp-windows-x64.zip && zip -qr ssp-windows-x64.zip ssp-windows-x64 && chown -R "${HOST_UID:-0}:${HOST_GID:-0}" /out
ls -la "$D" | head -50

# Общие зависимости расчётных модулей.
find_package(Eigen3 QUIET NO_MODULE)  # 3.4 или 5.x
if(NOT TARGET Eigen3::Eigen)
  message(FATAL_ERROR "Не найдена библиотека Eigen3 (3.4+ или 5.x). Установите: sudo pacman -S eigen")
endif()
find_package(PkgConfig REQUIRED)
pkg_check_modules(FFTW3 REQUIRED IMPORTED_TARGET fftw3)

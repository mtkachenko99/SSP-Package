cmake_minimum_required(VERSION 3.24)
project(SSPNano VERSION 0.1.0 LANGUAGES CXX)

# Общий суперпроект: каждый модуль — отдельный подкаталог со своим CMakeLists.txt.
# Модуль подключается автоматически, как только в нём появляется CMakeLists.txt.

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(SSP_BUILD_GUI   "Собирать Qt6 GUI"            ON)
option(SSP_BUILD_TESTS "Собирать тесты (ctest)"      ON)
option(SSP_USE_OPENMP  "Параллелизация через OpenMP" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(SSP_USE_OPENMP)
  find_package(OpenMP)
endif()
include(SSPDeps)

# Флаги предупреждений для собственного кода
set(SSP_WARNINGS -Wall -Wextra -Wpedantic)

if(SSP_BUILD_GUI)
  # Модуль для графиков (Qt Graphs / QCustomPlot) подключается в gui/CMakeLists.txt.
  find_package(Qt6 6.5 REQUIRED COMPONENTS Core Widgets)
  qt_standard_project_setup()
endif()

if(SSP_BUILD_TESTS)
  enable_testing()
endif()

set(SSP_MODULES
  common
  data_fetch
  plasmonics
  excitons
  spectroscopy/absorption
  spectroscopy/fluorescence
  spectroscopy/raman
  spectroscopy/brillouin
  validation
)
if(SSP_BUILD_GUI)
  list(APPEND SSP_MODULES gui)
endif()

foreach(mod IN LISTS SSP_MODULES)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${mod}/CMakeLists.txt)
    add_subdirectory(${mod})
  endif()
endforeach()

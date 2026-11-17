# ssp_excitons: агрегаты Френкеля произвольной формы, формы линий, Редфилд, HEOM, связь с плазмонной НЧ.
add_library(ssp_excitons
  src/geometry.cpp
  src/hamiltonian.cpp
  src/bath.cpp
  src/spectra.cpp
  src/relaxation.cpp
  src/modified_redfield.cpp
  src/vibronic.cpp
  src/monomer_fit.cpp
  src/huckel.cpp
  src/ensemble.cpp
  src/heom.cpp
  src/plexciton.cpp
)
target_include_directories(ssp_excitons PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ssp_excitons PUBLIC ssp_common ssp_plasmonics Eigen3::Eigen PRIVATE PkgConfig::FFTW3)
target_compile_options(ssp_excitons PRIVATE ${SSP_WARNINGS})
if(SSP_USE_OPENMP AND OpenMP_CXX_FOUND)
  target_link_libraries(ssp_excitons PUBLIC OpenMP::OpenMP_CXX)
endif()

# CLI: JSON-конфиг (Qt6 Core) → CSV/JSON-результаты
find_package(Qt6 6.5 REQUIRED COMPONENTS Core)
add_library(ssp_aggregate_core app/config.cpp app/output.cpp app/runner.cpp app/scan.cpp app/calibrate.cpp)
target_include_directories(ssp_aggregate_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/app)
target_link_libraries(ssp_aggregate_core PUBLIC ssp_excitons Qt6::Core)
target_compile_options(ssp_aggregate_core PRIVATE ${SSP_WARNINGS})
add_executable(ssp_aggregate app/main.cpp)
target_link_libraries(ssp_aggregate PRIVATE ssp_aggregate_core)
target_compile_options(ssp_aggregate PRIVATE ${SSP_WARNINGS})

if(SSP_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(test_excitons tests/test_excitons.cpp tests/test_heom.cpp tests/test_plexciton.cpp tests/test_scan.cpp tests/test_calibration.cpp)
  target_link_libraries(test_excitons PRIVATE ssp_aggregate_core GTest::gtest_main)
  target_compile_definitions(test_excitons PRIVATE
    SSP_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../data/reference/optical_constants/refractiveindex.info-database/database/data")
  add_test(NAME excitons COMMAND test_excitons)
endif()

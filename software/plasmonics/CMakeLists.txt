# ssp_plasmonics: материалы, теория Ми, квазистатический отклик сферы и сфероидов.
add_library(ssp_plasmonics
  src/material.cpp
  src/mie.cpp
  src/particle.cpp
  src/layered.cpp
  src/dipole_sphere.cpp
)
target_include_directories(ssp_plasmonics PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(ssp_plasmonics PUBLIC ssp_common Eigen3::Eigen)
target_compile_options(ssp_plasmonics PRIVATE ${SSP_WARNINGS})

if(SSP_BUILD_TESTS)
  find_package(GTest REQUIRED)
  add_executable(test_plasmonics tests/test_plasmonics.cpp)
  target_link_libraries(test_plasmonics PRIVATE ssp_plasmonics GTest::gtest_main)
  target_compile_definitions(test_plasmonics PRIVATE
    SSP_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../../data/reference/optical_constants/refractiveindex.info-database/database/data")
  add_test(NAME plasmonics COMMAND test_plasmonics)
endif()

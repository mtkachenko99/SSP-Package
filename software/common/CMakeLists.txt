# ssp_common: общие константы и единицы (header-only).
add_library(ssp_common INTERFACE)
target_include_directories(ssp_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(ssp_common INTERFACE cxx_std_20)

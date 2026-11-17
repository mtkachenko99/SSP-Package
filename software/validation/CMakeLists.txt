# ssp_benchmark: сравнение моделей программы с независимыми эталонами (точная электродинамика, литература).
add_executable(ssp_benchmark src/benchmark.cpp)
target_link_libraries(ssp_benchmark PRIVATE ssp_excitons ssp_plasmonics)
target_compile_options(ssp_benchmark PRIVATE ${SSP_WARNINGS})

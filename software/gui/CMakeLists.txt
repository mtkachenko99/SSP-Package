# ssp_gui: Qt6-оболочка над ssp_aggregate_core (параметры, запуск, графики, сравнение с экспериментом).
find_package(Qt6 6.5 REQUIRED COMPONENTS Widgets Concurrent Svg)
qt_add_executable(ssp_gui
  main.cpp
  mainwindow.cpp mainwindow.hpp
  plotwidget.cpp plotwidget.hpp
  geometryview.cpp geometryview.hpp
  formpanel.cpp formpanel.hpp
  csvtable.cpp csvtable.hpp
  experiment.cpp experiment.hpp
)
target_link_libraries(ssp_gui PRIVATE ssp_aggregate_core Qt6::Widgets Qt6::Concurrent Qt6::Svg)
target_compile_definitions(ssp_gui PRIVATE SSP_EXAMPLES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../excitons/examples")
target_compile_options(ssp_gui PRIVATE ${SSP_WARNINGS})
set_target_properties(ssp_gui PROPERTIES WIN32_EXECUTABLE ON MACOSX_BUNDLE ON)

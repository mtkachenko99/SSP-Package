#pragma once
// Калибровка по мономеру из JSON (раздел "monomer_fit"): чтение спектров, подгонка (ssp/excitons/monomer_fit),
// запись результатов и готового фрагмента конфига для расчёта агрегата.

#include "runner.hpp"

#include "ssp/excitons/monomer_fit.hpp"

#include <string>

namespace ssp::app {

struct CalibrationOutput {
  excitons::MonomerFitResult fit;
  excitons::Spectrum abs, em;      // экспериментальные данные в см⁻¹ (как подгонялись)
  std::string output_dir;
  std::string config_snippet;       // JSON: exciton/bath для конфига агрегата
};

bool is_calibration_config(const std::string& json);
std::string config_dir_of(const std::string& config_path);  // каталог файла (абсолютный)
CalibrationOutput run_calibration(const std::string& json, const std::string& base_dir, const Log& log);

}  // namespace ssp::app

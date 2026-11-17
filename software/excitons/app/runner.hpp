#pragma once
// Выполнение расчёта по конфигурации (общий код CLI и GUI).

#include "config.hpp"

#include <atomic>
#include <functional>
#include <map>
#include <vector>
#include <string>

namespace ssp::app {

using Log = std::function<void(const std::string& line)>;
using Progress = std::function<void(int index, int total, const std::string& name)>;

// Основные результаты одной формы (для сборки сканов в памяти).
struct ShapeResult {
  std::string name;
  int n_molecules = 0;
  std::vector<double> w, abs, em;               // спектры агрегата (средние по ансамблю)
  std::vector<double> vib_abs;                  // поглощение в модели Хольстейна (если включена)
  std::map<std::string, double> observables;    // средние значения
  bool has_np = false;
  std::vector<double> np_w, ext, abs_sys, sca;  // система агрегат+НЧ
  std::vector<double> np_ext, np_abs, np_sca;   // только НЧ
  double mode_omega_cm = 0, mode_gamma_cm = 0;  // дипольная мода НЧ (полюс α)
  double g_micro_cm = 0;                        // коллективная связь моды с агрегатом (√⟨g²⟩ по осям)
  std::vector<double> pl_hybrid, pl_free;       // ФЛ у НЧ и без неё (на сетке np_w)
  double qy_hybrid = 0, qy_free = 0;
};
struct RunResult {
  std::vector<ShapeResult> shapes;
};

// Выполняет все конфигурации, пишет результаты в rc.output_dir. cancel проверяется между формами.
RunResult run(const RunConfig& rc, const std::string& config_json, const Log& log,
              const std::atomic<bool>* cancel = nullptr, const Progress& progress = {});

}  // namespace ssp::app

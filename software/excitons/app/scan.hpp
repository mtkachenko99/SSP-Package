#pragma once
// Сканы по параметру конфигурации: для каждого значения — полный расчёт (run), затем сборка карт спектров,
// выделение гибридных ветвей и подгонка моделью двух связанных осцилляторов (плазмон + экситон):
//   H = [[ω_p − iγ_p/2, g], [g, ω_x − iγ_x/2]],  ω± = Re eig(H)  (γ — FWHM),
// с критериями сильной связи: 2g > (γ_p + γ_x)/2 (видимое расщепление), g > |γ_p − γ_x|/4 (антипересечение),
// C = 4g²/(γ_p γ_x) > 1 (кооперативность).

#include "runner.hpp"

#include <string>
#include <vector>

namespace ssp::app {

struct ScanSpec {
  bool enabled = false;
  std::string parameter;              // путь в JSON: "nanoparticle.semi_axis_xy_nm", "exciton.monomer_energy_cm",
                                      // "geometry_params.N" (все формы), "configurations[name].geometry_params.N"
  std::vector<double> values;
  std::vector<std::string> also_set;  // дополнительные пути, получающие то же значение
};

struct BranchFit {
  std::string spectrum;  // abs | ext | sca | micro
  double g_cm = 0.0, rms_cm = 0.0;
  double g_lo = 0.0, g_hi = 0.0;  // интервал Δχ² ≤ 1 (σ² = остаточная дисперсия)
  bool resolved = false;          // ≥ 2 точек скана с двумя пиками у разных ветвей
  int two_peak_points = 0;
  int n_peaks = 0;
  double gamma_p = 0.0, gamma_x = 0.0, detuning_min = 0.0, observed_splitting = 0.0;
  bool visible = false, anticrossing = false;
  double cooperativity = 0.0;
};

struct ScanShape {
  std::string name;
  std::vector<double> w_agg;                    // сетка спектров агрегата
  std::vector<std::vector<double>> abs, em;     // [значение][ω]
  std::vector<std::map<std::string, double>> observables;
  bool has_np = false;
  std::vector<double> w_np;
  std::vector<std::vector<double>> ext, abs_sys, sca, np_abs;
  std::vector<double> w_p, gamma_p, w_x, gamma_x;  // ω_p, γ_p — полюс дипольной моды; ω_x, γ_x — агрегат без НЧ
  std::vector<double> g_micro;
  std::vector<std::vector<double>> peaks_abs, peaks_ext, peaks_sca;
  std::vector<BranchFit> fits;
};

struct ScanResult {
  std::vector<double> values;
  std::vector<ScanShape> shapes;
};

// Разбор раздела "scan" (пустой spec.enabled = false, если раздела нет).
ScanSpec parse_scan(const std::string& config_json);

ScanResult run_scan(const std::string& config_json, const std::string& base_dir, const std::string& config_path,
                    const Log& log, const std::atomic<bool>* cancel = nullptr, const Progress& progress = {});

// Анализ ветвей (доступен отдельно для тестов): подгонка g по положениям пиков.
BranchFit fit_coupled_oscillators(const std::vector<double>& w_p, const std::vector<double>& gamma_p,
                                  const std::vector<double>& w_x, const std::vector<double>& gamma_x,
                                  const std::vector<std::vector<double>>& peaks, const std::string& name);
// Re ω± модели двух осцилляторов.
void coupled_branches(double wp, double gp, double wx, double gx, double g, double& lower, double& upper);

}  // namespace ssp::app

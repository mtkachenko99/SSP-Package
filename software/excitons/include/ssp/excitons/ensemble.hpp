#pragma once
// Усреднение по статическому беспорядку: спектры поглощения и флуоресценции, кинетика, скалярные наблюдаемые
// с оценкой стандартной ошибки; число реализаций — до достижения заданной точности.

#include "ssp/excitons/bath.hpp"
#include "ssp/excitons/geometry.hpp"
#include "ssp/excitons/hamiltonian.hpp"
#include "ssp/excitons/relaxation.hpp"
#include "ssp/excitons/spectra.hpp"

#include <Eigen/Dense>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ssp::excitons {

enum class PumpType { OscillatorStrength, Laser };
// Теория переноса между экситонными состояниями: секулярный Редфилд или модифицированный Редфилд
// (непертурбативно по диагональной части экситон-фононной связи; для N ≲ 300 — O(N²·n_t)).
enum class RelaxationTheory { Redfield, ModifiedRedfield };

struct ExcitonModel {
  CouplingParams coupling;
  double monomer_energy_cm = 17000.0;
  double disorder_sigma_cm = 60.0;
  Bath bath;
  DecayParams decay;
  PumpType pump = PumpType::OscillatorStrength;
  double excitation_cm = 17500.0;
  double n_medium = 1.0;  // для абсолютного сечения поглощения
  // Беспорядок сверх некоррелированного диагонального:
  double disorder_corr_length_nm = 0.0;   // корреляция энергий сайтов: ⟨δ_iδ_j⟩ = σ² exp(−r_ij/ξ)
  double orientation_jitter_deg = 0.0;    // структурный: разброс направлений диполей (в каждой реализации)
  double position_jitter_nm = 0.0;        // структурный: разброс положений
  RelaxationTheory relaxation = RelaxationTheory::Redfield;
  DispersionParams dispersion;           // зависящий от окружения сдвиг энергий сайтов
  Eigen::VectorXd site_shift_cm;         // постоянный сдвиг энергий сайтов (если пусто — нет); заполняется
                                         // run_ensemble из dispersion
};

struct SpectrumGrid {
  double wmin_cm = 15500.0;
  double wmax_cm = 18000.0;
  std::size_t nw = 1800;
  double t_max_ps = 100.0;  // временное окно кумулянтных форм линий
  std::vector<double> points() const;
};

struct EnsembleSettings {
  int min_realizations = 20;
  int max_realizations = 400;
  double rel_tolerance = 0.01;  // max SE(ω)/max mean(ω) для поглощения и флуоресценции
  std::uint64_t seed = 1234;
  int batch = 16;
  int bootstrap = 200;                 // ресэмплов для погрешности наблюдаемых усреднённого спектра
  std::vector<double> dyn_times_ps;    // если не пусто — кинетика населённостей по ансамблю
  std::string dyn_initial = "brightest";  // brightest | highest | pump
  double extra_wmin_cm = 0.0, extra_wmax_cm = 0.0;  // дополнительный диапазон для форм линий (спектры с НЧ)
};

// Среднее и стандартная ошибка (накопление по Уэлфорду).
struct RunningStat {
  long n = 0;
  double mean = 0.0, m2 = 0.0;
  double se_fixed = -1.0;  // если ≥ 0 — погрешность задана извне (бутстрэп)
  void add(double x);
  double se() const;
  double sd() const;
};

struct Realization {
  Eigen::MatrixXd H;            // см⁻¹
  ExcitonStates states;
  Eigen::MatrixXd K_ps;         // скорости Редфилда, пс⁻¹
  Eigen::VectorXd hwhm_cm;      // однородная полуширина по Редфилду
  Eigen::VectorXd gamma_rad_ps;
  Eigen::VectorXd db_energies;  // энергии детального баланса K: ε_n (Редфилд) или ε_n − λ s_n (модифицированный)
  Eigen::VectorXd pump;
  Eigen::VectorXd P;            // стационарные населённости
};

struct EnsembleResult {
  std::vector<double> w;
  std::vector<double> abs_mean, abs_se;  // σ_abs, нм² на агрегат
  std::vector<double> em_mean, em_se;    // фотонов/(пс·см⁻¹) на единицу накачки
  std::map<std::string, RunningStat> observables;
  int realizations = 0;
  bool converged = false;
  double max_lineshape_tail = 0.0;
  Realization representative;  // первая реализация (seed)
  // Кинетика по ансамблю: t, затем [total, lowest, brightest, mean_energy] × время (среднее и SE)
  std::vector<double> dyn_t;
  std::vector<std::vector<RunningStat>> dyn;  // 4 × n_t
  std::vector<Realization> stored;  // первые n_store реализаций (для плекситонного расчёта)
  std::shared_ptr<LineshapeEngine> engine;  // тот же движок форм линий (для согласованного отклика с НЧ)
};

// Один прогон реализации k (детерминированное зерно seed + k).
// chol — множитель Холецкого ковариации энергий (nullptr — некоррелированный беспорядок).
// engine — нужен для модифицированного Редфилда (g(t)); без него используется обычный Редфилд.
Realization make_realization(const Aggregate& agg, const Eigen::MatrixXd& J, const ExcitonModel& m,
                             std::uint64_t seed, const Eigen::MatrixXd* chol = nullptr,
                             const LineshapeEngine* engine = nullptr);

EnsembleResult run_ensemble(const Aggregate& agg, const ExcitonModel& model, const SpectrumGrid& grid,
                            const EnsembleSettings& settings, int n_store = 0);

// Спектральные наблюдаемые усреднённого по ансамблю спектра (как в эксперименте).
std::map<std::string, double> spectral_metrics(const std::vector<double>& w, const std::vector<double>& abs,
                                               const std::vector<double>& em, double monomer_cm);

// Описание скалярных наблюдаемых (имя → единицы и смысл) для отчёта.
const std::map<std::string, std::string>& observable_descriptions();

}  // namespace ssp::excitons

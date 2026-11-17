#pragma once
// Калибровка по мономеру (аудит П4): подгонка спектров поглощения и ФЛ одной молекулы кумулянтной моделью
//   A(ω) ∝ ω·[L_abs ⊗ G_σ](ω),  F(ω) ∝ ω³·[L_em ⊗ G_σ](ω),
// L — кумулянтные формы линии с баней «Друде (λ, ω_c) + броуновские моды (ω_k, S_k, γ_k)», G_σ — статическое
// неоднородное уширение. Масштаб каждого спектра — аналитически (линейный МНК). Левенберг–Марквардт с
// численным якобианом; погрешности — из (JᵀJ)⁻¹·χ²_ν.
// Из абсолютного спектра ε(ν̃) (М⁻¹см⁻¹): μ² = ∫σ(ν̃)/ν̃ dν̃ / P(n) — та же формула, что и в σ_abs модели агрегата;
// излучательное время — та же формула Γ = nω³μ²/(3πε₀ħc³) при ν̃ = ⟨ν̃_em⁻³⟩^{−1/3}, и для сравнения
// Стриклер–Берг: k_r = 2.880·10⁻⁹ n² ⟨ν̃_f⁻³⟩⁻¹ ∫ε dν̃/ν̃ (с⁻¹).

#include "ssp/excitons/bath.hpp"

#include <string>
#include <vector>

namespace ssp::excitons {

struct Spectrum {
  std::vector<double> w, y, sigma;  // см⁻¹; sigma пусто — равные веса (σ = 2 % от максимума)
};

struct MonomerFitParams {
  double E0_cm = 18000.0;        // вертикальная энергия (центр поглощения без мод)
  double lambda_cm = 100.0;      // реорганизация Друде
  double cutoff_cm = 100.0;      // ω_c Друде
  double sigma_inh_cm = 100.0;   // статическое неоднородное уширение (стандартное отклонение)
  std::vector<VibrationalMode> modes;
  // какие параметры подгоняются
  bool fit_E0 = true, fit_lambda = true, fit_cutoff = true, fit_sigma = true, fit_mode_S = true,
       fit_mode_omega = false;
  double temperature_K = 300.0;
  // Начальное приближение по моментам спектров (точные правила кумулянта): E0 = ⟨ν̃⟩_abs (по A/ν̃),
  // 2λ_tot = ⟨ν̃⟩_abs − ⟨ν̃⟩_em (по F/ν̃³), σ² = Var_abs − Var_одн(λ, ω_c, моды, T); S мод и ω_c — из start.
  bool moment_start = true;
  // Необязательное ограничение: измеренный сдвиг пиков поглощения и ФЛ (см⁻¹) ± погрешность — когда спектр ФЛ
  // недоступен, но сдвиг известен из литературы. 0 — не используется.
  double stokes_peak_cm = 0.0, stokes_sigma_cm = 0.0;
};

struct MonomerFitResult {
  MonomerFitParams p;
  std::vector<std::string> names;
  std::vector<double> values, errors;
  double chi2_red = 0.0;
  int iterations = 0;
  bool converged = false;
  // подогнанные кривые на сетках эксперимента (в масштабе эксперимента)
  std::vector<double> abs_model, em_model;
  double abs_scale = 0.0, em_scale = 0.0;
  // из абсолютного ε (если задано absolute_extinction)
  double mu_D = 0.0;            // переходный диполь, Д
  double tau_rad_ns = 0.0;      // излучательное время по формуле модели
  double tau_sb_ns = 0.0;       // Стриклер–Берг
  double nu_em_cm = 0.0;        // ⟨ν̃⁻³⟩^{−1/3} по спектру ФЛ (модельному)
  double stokes_peak_cm = 0.0;  // модельный сдвиг пиков поглощения и ФЛ
};

// abs — поглощение (если absolute_extinction, то ε в М⁻¹см⁻¹), em — ФЛ (может быть пуст). n_medium — для μ и τ.
MonomerFitResult fit_monomer(const Spectrum& abs, const Spectrum& em, const MonomerFitParams& start,
                             bool absolute_extinction, double n_medium, int max_iter = 100);

// Модельные спектры (без масштаба: площадь L = 1) на сетке w.
void monomer_model(const MonomerFitParams& p, const std::vector<double>& w, std::vector<double>& abs,
                   std::vector<double>& em);

// Сдвиг пиков поглощения и ФЛ модели (см⁻¹), поиск на сетке 2 см⁻¹ вокруг E0.
double monomer_stokes_peak(const MonomerFitParams& p);

// μ² (Д²) из ε(ν̃) (М⁻¹см⁻¹) в среде n, согласованно с absorption_cross_section_prefactor_nm2.
double mu2_from_extinction(const Spectrum& eps, double n_medium);

}  // namespace ssp::excitons

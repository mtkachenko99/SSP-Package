#pragma once
// Линейный отклик системы «молекулярный агрегат + плазмонная наночастица» (НЧ в начале координат).
//   H_eff(ω) = H + W(ω) − i·Γ_hom,  W_ij = +(5.034μ²/ε_m)(μ̂_i·∇)(μ̂_j·∇')K_ind(r_i, r_j; ω)
//   α_tot(ω) = α_NP(ω) + (5.034μ²/ε_m)·Mᵀ (H_eff − ω)⁻¹ M,  M_i = L_iᵀ μ̂_i,  L_i — тензор локального поля
//   σ_ext = 4πk Im Tr α/3,  σ_sca = (8π/3)k⁴ Σ_ab|α_ab|²/3,  σ_abs = σ_ext − σ_sca   (α в нм³, k в нм⁻¹)
// Дипольные (l = 1) компоненты отклика НЧ поправлены на запаздывание: сфера — точный a₁ Ми, сфероид — MLWA.

#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/geometry.hpp"
#include "ssp/excitons/spectra.hpp"
#include "ssp/plasmonics/material.hpp"
#include "ssp/plasmonics/layered.hpp"
#include "ssp/plasmonics/particle.hpp"

#include <Eigen/Dense>

#include <array>
#include <complex>
#include <memory>
#include <vector>

namespace ssp::excitons {

struct NanoparticleSpec {
  double a_xy_nm = 15.0;
  double a_z_nm = 15.0;
  int l_max = 40;
  double n_medium = 1.33;
  bool dynamic_dipole = true;      // сфера: динамика всех l ≤ 8 по Ми; сфероид: MLWA для l = 1
  double shell_thickness_nm = 0.0;  // диэлектрическая оболочка (лиганды) — только для сферы
  double shell_n = 1.45;
  std::complex<double> d_perp_nm = 0.0;  // Фейбельман d⊥ на поверхности металла — только для сферы
  // Внешний радиус сферы (с оболочкой) — его использует геометрия частицы и молекул
  // Быстрый решатель для больших N: мультиполи l ≤ fast_l — точно (низкий ранг, Вудбери), l > fast_l — только
  // диагональ в базисе экситонов (секулярно). 0 — точное решение N×N на каждой частоте; < 0 — автоматически.
  int fast_l = 0;
  // Распределение размеров (нормальное, относительная σ радиуса) — усреднение спектров по узлам Гаусса–Эрмита
  double size_sigma_rel = 0.0;
  int size_nodes = 5;
  double outer_xy() const { return a_xy_nm + shell_thickness_nm; }
  double outer_z() const { return a_z_nm + shell_thickness_nm; }
};

struct PlexcitonSpectra {
  std::vector<double> w;
  std::vector<double> ext, abs, sca;           // система целиком, нм²
  std::vector<double> np_ext, np_abs, np_sca;  // только НЧ
};

class PlexcitonSystem {
 public:
  PlexcitonSystem(const Aggregate& agg, const NanoparticleSpec& spec, const plasmonics::Material& mat, double mu_D);

  // Спектры для набора реализаций беспорядка; среднее по реализациям и стандартная ошибка.
  // Отклик каждого экситона — тот же, что в спектре агрегата без НЧ: кумулянтная резольвента i·X_n(ω) из engine
  // (если engine == nullptr или fixed_hwhm_cm ≥ 0 — лоренциан с полушириной fixed_hwhm_cm или Редфилда).
  void spectra(const std::vector<Realization>& reals, const LineshapeEngine* engine, double fixed_hwhm_cm,
               const std::vector<double>& w, PlexcitonSpectra& mean, PlexcitonSpectra& se) const;

  // Для каждой молекулы при частоте ν̃: усреднённое по поляризации |E_loc/E₀|², сдвиг Re W_ii (см⁻¹)
  // и скорость распада в НЧ −2 Im W_ii (пс⁻¹).
  struct MoleculeNearField {
    double field_enhancement = 1.0;
    double shift_cm = 0.0;
    double decay_ps = 0.0;            // Γ_tot − Γ₀ = Γ_abs + Γ₀(|μ_eff/μ|² − 1) (без Γ₀ — только Γ_abs)
    double gap_nm = 0.0;
    double rad_enhancement = 1.0;     // Γ_rad/Γ_rad,0 = |L_iᵀμ̂_i|²
    double nonrad_to_np_ps = 0.0;     // Γ_abs = −2 Im W_ii без радиационного затухания НЧ: потери в металле
    double quantum_yield = 0.0;       // мономера у НЧ (с собственным Γ_nr и Γ₀ из modelи)
    double quantum_yield_free = 0.0;  // мономера без НЧ
  };
  // gamma0_ps — излучательная скорость мономера без НЧ, gamma_nr_ps — собственная безызлучательная.
  std::vector<MoleculeNearField> near_field(double nu_cm, double gamma0_ps = 0.0, double gamma_nr_ps = 0.0) const;

  // Излучение гибрида: для каждой реализации — излучательные и безызлучательные (в металл) скорости экситонов,
  // стационарные населённости, квантовый выход и спектр ФЛ; среднее по реализациям. Скорости — на сетке узлов
  // по энергии (интерполяция по состояниям).
  struct HybridEmission {
    std::vector<double> w, pl, pl_free;  // фотонов/(пс·см⁻¹) на единицу накачки
    double qy = 0.0, qy_free = 0.0;
    double rad_enhancement_bright = 0.0;  // ⟨Γ_rad^NP/Γ_rad⟩ для самого яркого состояния
    double lifetime_ps = 0.0, lifetime_free_ps = 0.0;
  };
  HybridEmission emission(const std::vector<Realization>& reals, const LineshapeEngine& engine,
                          const std::vector<double>& w, double nonradiative_ps, bool laser_pump,
                          double excitation_cm, int nodes = 12) const;

  // Дипольная (l = 1) плазмонная мода: подгонка α_a(ω) ≈ A/(ω_p² − ω² − iγω) + B по окну вокруг резонанса
  // и коллективная константа связи моды со всеми молекулами g² = (A/2ω_p)·Tr K_a, где W^(l=1,a) = α_a(ω)·K_a.
  // Не зависит от беспорядка (след инвариантен), не требует разрешённых гибридных пиков.
  struct DipoleModeCoupling {
    double omega_p_cm = 0, gamma_p_cm = 0, residue = 0, fit_rel_rms = 0;  // для оси с наибольшим g
    double omega_axis[3] = {0, 0, 0}, gamma_axis[3] = {0, 0, 0};
    double g_cm[3] = {0, 0, 0};  // по осям x, y, z
    double g_rms_cm = 0;         // √(⟨g²⟩ по осям) — для изотропного возбуждения
  };
  DipoleModeCoupling dipole_mode_coupling(const std::vector<double>& w) const;

  // Сходимость мультипольного ряда: относительный вклад членов с l > L/2 в W (макс. по элементам) при ν̃.
  double truncation_estimate(double nu_cm) const;

  const plasmonics::QuasistaticParticle& particle() const { return *particle_; }
  int fast_l() const { return spec_.fast_l; }

 private:
  const Aggregate agg_;
  NanoparticleSpec spec_;
  const plasmonics::Material& mat_;
  double mu_D_;
  double eps_m_;
  std::unique_ptr<plasmonics::QuasistaticParticle> particle_;
  plasmonics::SphereLayers layers_;
  bool layered_ = false;  // сфера: отклик по слоям (оболочка, d⊥) с динамикой Ми
  Eigen::MatrixXcd Vc_, Vs_;  // N × T: μ̂_i·∇f_t(r_i)
  // Сфера: D̂ зависит только от l ⇒ W = pref·Σ_l D̂_l(ω)·G_l, G_l = Σ_m w_lm (V V^T) — предвычисленные N×N
  std::vector<Eigen::MatrixXd> Gl_;
  // ∇f для дипольных членов (l=1) в точках молекул: grad1_[i][t], t ∈ {(1,0),(1,1)}
  std::vector<Eigen::Vector3cd> g10c_, g11c_, g11s_;

  struct Freq {
    Eigen::MatrixXcd W;                 // N × N, см⁻¹
    std::vector<Eigen::Matrix3cd> L;    // тензоры локального поля
    Eigen::Vector3cd alpha_np;          // нм³
    Eigen::VectorXcd d;                 // w_t·D_t по членам (без pref)
    Eigen::VectorXcd d_abs;             // то же без радиационного затухания НЧ (поглощение в металле)
    double k = 0.0;                     // нм⁻¹
  };
  Freq at(double nu_cm, int l_cut = -1, bool with_W = true) const;
  // Проекции V на экситонный базис: Ṽc = UᵀVc, Ṽs = UᵀVs (N × T)
  void project(const Eigen::MatrixXd& U, Eigen::MatrixXcd& Pc, Eigen::MatrixXcd& Ps) const;
};

// Квадратура для нормального распределения N(0,1): узлы x_k и веса (Σw = 1), Гаусс–Эрмит (Голуб–Уэлш).
std::vector<std::pair<double, double>> normal_quadrature(int n);

// Агрегат при изменении полуосей ядра НЧ в s раз: каждая молекула сдвигается вдоль своего радиального направления
// на (s − 1)·ρ(r̂), ρ — радиус ядра в этом направлении, так что зазор до поверхности сохраняется (лиганды
// «следуют» за поверхностью). Межмолекулярные расстояния меняются в ~s раз — для σ ≲ 0.1 малая поправка к J,
// которая не учитывается (экситонные состояния берутся из исходного агрегата).
Aggregate follow_core_scale(const Aggregate& agg, const NanoparticleSpec& spec, double s);

}  // namespace ssp::excitons

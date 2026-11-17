#pragma once
// Вибронная модель Хольстейна в одночастичном приближении (OPA, Spano): одна внутримолекулярная мода ω₀ с
// фактором Хуанга–Риса S на молекулу. Базис |n, ν⟩ — молекула n возбуждена с ν квантами в возбуждённом потенциале,
// остальные — в колебательном основном состоянии.
//   ⟨nν|H|nν⟩ = H_nn − Sω₀ + νω₀   (H_nn — вертикальная энергия сайта: первый момент спектра мономера = H_nn)
//   ⟨nν|H|mν'⟩ = J_nm F_ν F_ν',   F_ν = ⟨χ_g0|χ_eν⟩ = e^{−S/2} S^{ν/2}/√ν!
//   μ_k = Σ_{n,ν} c_{nν,k} F_ν μ_n
// В отличие от кумулянта с модой в g(t) (вибронный повтор каждой экситонной линии), OPA описывает перенормировку
// связи для вибронных полос: отношение A₀₋₀/A₀₋₁ растёт для J- и падает для H-агрегатов.
// Мода ω₀ не должна одновременно входить в bath.modes (двойной учёт).

#include "ssp/excitons/ensemble.hpp"
#include "ssp/excitons/spectra.hpp"

#include <Eigen/Dense>

#include <vector>

namespace ssp::excitons {

struct HolsteinParams {
  double omega_cm = 1300.0;
  double huang_rhys = 0.6;
  int nu_max = 4;  // ν = 0..nu_max
};

struct VibronicStates {
  Eigen::VectorXd energies;  // см⁻¹, по возрастанию
  Eigen::VectorXd mu2;       // Д²
  Eigen::VectorXd ipr;       // Σ_n (Σ_ν c²_{nν,k})² — для формы линии по оставшейся бане
};

// H — электронный гамильтониан реализации (см⁻¹), mu_sites — переходные диполи сайтов (Д), строка n.
VibronicStates holstein_opa(const Eigen::MatrixXd& H, const Eigen::MatrixX3d& mu_sites, const HolsteinParams& p);

// Переходные диполи сайтов из экситонных состояний реализации: μ_n = Σ_k c_nk μ_k.
Eigen::MatrixX3d site_dipoles(const ExcitonStates& st);

// Поглощение (σ_abs, нм², как в EnsembleResult::abs_mean) по реализациям: среднее и SE. Форма линии вибронного
// состояния — кумулянт с IPR по сайтам и скоростью ухода электронного экситона, ближайшего по энергии.
void holstein_absorption(const std::vector<Realization>& reals, const HolsteinParams& p, const LineshapeEngine& engine,
                         double n_medium, const std::vector<double>& w, std::vector<double>& mean,
                         std::vector<double>& se);

}  // namespace ssp::excitons

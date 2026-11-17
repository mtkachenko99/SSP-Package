#pragma once
// Диэлектрическая функция материала наночастицы: таблицы refractiveindex.info ("tabulated nk") или модель Друде,
// с опциональной поправкой на размерное затухание (Крейбиг): γ(R) = γ_b + A·v_F/R.

#include <complex>
#include <string>
#include <vector>

namespace ssp::plasmonics {

struct DrudeParams {
  double eps_inf = 1.0;
  double omega_p_eV = 9.0;
  double gamma_eV = 0.07;
};

class Material {
 public:
  // Таблица n,k из YAML-файла refractiveindex.info (первый блок "tabulated nk").
  static Material from_refractiveindex_yaml(const std::string& path);
  static Material drude(const DrudeParams& p, const std::string& name = "drude");

  // ε(ν̃), ν̃ в см⁻¹. За пределами таблицы — исключение.
  std::complex<double> epsilon(double nu_cm) const;

  // Крейбиг: ε → ε + ω_p²/(ω² + iγ_bω) − ω_p²/(ω² + i(γ_b + A v_F/R)ω). Выключено при A = 0.
  void set_size_damping(double A, double vF_m_per_s, double R_eff_nm, double omega_p_eV, double gamma_bulk_eV);

  const std::string& name() const { return name_; }
  const std::string& reference() const { return reference_; }
  double lambda_min_um() const { return lam_.empty() ? 0.0 : lam_.front(); }
  double lambda_max_um() const { return lam_.empty() ? 0.0 : lam_.back(); }

 private:
  std::string name_, reference_;
  bool tabulated_ = false;
  std::vector<double> lam_, n_, k_;  // мкм
  DrudeParams drude_;
  double kreibig_extra_eV_ = 0.0, kreibig_wp_eV_ = 0.0, kreibig_gb_eV_ = 0.0;
};

}  // namespace ssp::plasmonics

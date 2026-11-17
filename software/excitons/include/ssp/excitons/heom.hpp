#pragma once
// Иерархические уравнения движения (HEOM) для одноэкситонного блока с независимыми банями Друде–Лоренца
// на каждом сайте (Tanimura; Ishizaki–Fleming, J. Chem. Phys. 130, 234111 (2009)).
// Численно точный эталон для проверки кумулянтных форм линий и Редфилда на малых агрегатах.

#include <Eigen/Dense>

#include <complex>
#include <vector>

namespace ssp::excitons {

struct HeomParams {
  double lambda_cm = 50.0;
  double gamma_cm = 100.0;
  double temperature_K = 300.0;
  int matsubara = 1;        // число мацубаровских членов K
  int depth = 4;            // глубина иерархии L
  bool terminator = true;   // учёт отброшенных мацубаровских членов (Ishizaki–Fleming)
  double dt_fs = 0.0;       // шаг RK4; 0 — автоматически
  double t_max_ps = 1.0;
};

struct HeomInfo {
  std::size_t n_ado = 0;
  double dt_fs = 0.0;
  std::size_t steps = 0;
};

struct HeomAbsorption {
  std::vector<double> t_ps;
  std::vector<std::complex<double>> response;  // R(t) в системе отсчёта, сдвинутой на e0_cm
  double e0_cm = 0.0;
  std::vector<double> w_cm;
  std::vector<double> lineshape;  // нормированная форма (площадь ≈ 1), без множителя ω
  HeomInfo info;
};

struct HeomPopulations {
  std::vector<double> t_ps;
  Eigen::MatrixXd site;     // N × n_t
  Eigen::MatrixXd exciton;  // N × n_t (в собственном базисе H)
  HeomInfo info;
};

std::size_t heom_ado_count(int n_sites, int matsubara, int depth);

// Линейное поглощение: R(t) = Σ_pol μ_polᵀ σ₀(t), σ(0) = μ_pol. site_mu — N×3 (Д).
HeomAbsorption heom_absorption(const Eigen::MatrixXd& H_cm, const Eigen::MatrixX3d& site_mu, const HeomParams& p,
                               const std::vector<double>& w_cm);

// Динамика населённостей из начальной матрицы плотности rho0 (N×N).
HeomPopulations heom_populations(const Eigen::MatrixXd& H_cm, const Eigen::MatrixXcd& rho0, const HeomParams& p,
                                 int n_out);

}  // namespace ssp::excitons
